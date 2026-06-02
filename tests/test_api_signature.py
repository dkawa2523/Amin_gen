from gas_screening_mvp.domain.api import ApiRequest


def test_api_signature_order_insensitive():
    a = ApiRequest("P", "/x", "GET", {"b": 2, "a": 1}, None).signature()
    b = ApiRequest("P", "/x", "GET", {"a": 1, "b": 2}, None).signature()
    assert a == b
