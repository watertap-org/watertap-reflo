import watertap_contrib.reflo


def test_dunder_version_module_attribute():
    assert hasattr(watertap_contrib.reflo, "__version__")
