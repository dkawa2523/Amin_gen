from __future__ import annotations

from pathlib import Path
from datetime import date, datetime
from numbers import Number
from xml.sax.saxutils import escape
import zipfile

import pandas as pd


def export_excel(sheets: dict[str, list[dict]], output_xlsx: str | Path) -> None:
    output = Path(output_xlsx)
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with pd.ExcelWriter(output, engine="openpyxl") as writer:
            for sheet_name, rows in sheets.items():
                df = pd.DataFrame(rows)
                if df.empty:
                    df = pd.DataFrame([{"message": "no rows"}])
                safe_name = sheet_name[:31]
                df.to_excel(writer, index=False, sheet_name=safe_name)
                ws = writer.book[safe_name]
                ws.freeze_panes = "A2"
                for column_cells in ws.columns:
                    length = max(len(str(cell.value)) if cell.value is not None else 0 for cell in column_cells)
                    ws.column_dimensions[column_cells[0].column_letter].width = min(max(length + 2, 10), 60)
    except ModuleNotFoundError as exc:
        if exc.name != "openpyxl":
            raise
        _export_minimal_xlsx(sheets, output)


def _export_minimal_xlsx(sheets: dict[str, list[dict]], output: Path) -> None:
    sheet_items = [(name[:31], _rows_with_headers(rows)) for name, rows in sheets.items()]
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", _content_types_xml(len(sheet_items)))
        zf.writestr("_rels/.rels", _root_rels_xml())
        zf.writestr("xl/workbook.xml", _workbook_xml([name for name, _rows in sheet_items]))
        zf.writestr("xl/_rels/workbook.xml.rels", _workbook_rels_xml(len(sheet_items)))
        zf.writestr("xl/styles.xml", _styles_xml())
        for idx, (_name, rows) in enumerate(sheet_items, start=1):
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", _worksheet_xml(rows))


def _rows_with_headers(rows: list[dict]) -> list[list[object]]:
    if not rows:
        return [["message"], ["no rows"]]
    headers: list[str] = []
    for row in rows:
        for key in row:
            if key not in headers:
                headers.append(key)
    return [headers] + [[row.get(header) for header in headers] for row in rows]


def _worksheet_xml(rows: list[list[object]]) -> str:
    body = []
    for r_idx, row in enumerate(rows, start=1):
        cells = []
        for c_idx, value in enumerate(row, start=1):
            cells.append(_cell_xml(_cell_ref(c_idx, r_idx), value))
        body.append(f'<row r="{r_idx}">{"".join(cells)}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
        'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
        f'<sheetData>{"".join(body)}</sheetData>'
        "</worksheet>"
    )


def _cell_xml(ref: str, value: object) -> str:
    if value is None:
        return f'<c r="{ref}"/>'
    if isinstance(value, bool):
        return f'<c r="{ref}" t="b"><v>{1 if value else 0}</v></c>'
    if isinstance(value, Number):
        return f'<c r="{ref}"><v>{value}</v></c>'
    if isinstance(value, (datetime, date)):
        value = value.isoformat()
    text = escape(str(value))
    return f'<c r="{ref}" t="inlineStr"><is><t>{text}</t></is></c>'


def _cell_ref(col_idx: int, row_idx: int) -> str:
    letters = ""
    while col_idx:
        col_idx, rem = divmod(col_idx - 1, 26)
        letters = chr(65 + rem) + letters
    return f"{letters}{row_idx}"


def _content_types_xml(sheet_count: int) -> str:
    sheets = "".join(
        f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for i in range(1, sheet_count + 1)
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        f"{sheets}</Types>"
    )


def _root_rels_xml() -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        "</Relationships>"
    )


def _workbook_xml(sheet_names: list[str]) -> str:
    sheets = "".join(
        f'<sheet name="{escape(name)}" sheetId="{idx}" r:id="rId{idx}"/>'
        for idx, name in enumerate(sheet_names, start=1)
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{sheets}</sheets></workbook>"
    )


def _workbook_rels_xml(sheet_count: int) -> str:
    sheets = "".join(
        f'<Relationship Id="rId{i}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        f'Target="worksheets/sheet{i}.xml"/>'
        for i in range(1, sheet_count + 1)
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f"{sheets}</Relationships>"
    )


def _styles_xml() -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
        '<borders count="1"><border/></borders>'
        '<cellStyleXfs count="1"><xf/></cellStyleXfs><cellXfs count="1"><xf/></cellXfs>'
        "</styleSheet>"
    )
