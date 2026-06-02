from gas_screening_mvp.classification.pfas_rules import detect_fully_fluorinated_cf2_or_cf3, has_c_f_bond


def test_pfas_motif_cf4_positive():
    hit, motifs = detect_fully_fluorinated_cf2_or_cf3("FC(F)(F)F")
    assert hit
    assert motifs


def test_difluoromethane_possible_not_positive():
    hit, motifs = detect_fully_fluorinated_cf2_or_cf3("FCF")
    assert not hit
    assert has_c_f_bond("FCF")
