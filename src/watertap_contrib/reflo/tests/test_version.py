import watertap_contrib.seto


def test_dunder_version_module_attribute():
    assert hasattr(watertap_contrib.seto, "__version__")
