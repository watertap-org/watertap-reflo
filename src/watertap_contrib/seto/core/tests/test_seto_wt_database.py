"""
Tests for SETO WaterTAP database wrapper
"""

import pytest
import os

from watertap_contrib.seto.core.seto_wt_database import SETODatabase


@pytest.mark.unit
def test_default_path():
    seto_db = SETODatabase()

    assert (
        os.path.normpath(seto_db._dbpath).casefold()
        == os.path.normpath(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..",
                "..",
                "data",
                "technoeconomic",
            )
        ).casefold()
    )


@pytest.mark.unit
def test_invalid_path():
    with pytest.raises(
        OSError,
        match="Could not find requested path foobar. Please "
        "check that this path exists.",
    ):
        SETODatabase(dbpath="foobar")


@pytest.mark.unit
def test_custom_path():
    # Pick a path we know will exist, even if it isn't a data folder
    seto_db = SETODatabase(
        dbpath=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"
        )
    )

    assert seto_db._dbpath == os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"
    )


class TestSETODatabase:
    @pytest.fixture(scope="class")
    def seto_db(self):
        return SETODatabase()

    @pytest.mark.unit
    def test_component_list(self, seto_db):
        assert seto_db._component_list is None

        assert isinstance(seto_db.component_list, dict)

    @pytest.mark.unit
    def test_get_technology(self, seto_db):
        assert seto_db._cached_files == {}

        data = seto_db._get_technology("solar_energy")

        assert "solar_energy" in seto_db._cached_files
        assert seto_db._cached_files["solar_energy"] is data

        assert "default" in data
        assert "pv" in data
        assert "capital_cost" in data["default"]
        assert "operating_cost" in data["pv"]
        assert "variable_op_by_generation" in data["pv"]["operating_cost"]

    @pytest.mark.unit
    def test_get_technology_invalid(self, seto_db):
        with pytest.raises(
            KeyError, match="Could not find entry for foobar in database."
        ):
            seto_db._get_technology("foobar")

        assert len(seto_db._cached_files) == 1
        assert "foobar" not in seto_db._cached_files

    @pytest.mark.unit
    def test_get_unit_operation_parameters_default(self, seto_db):
        data = seto_db.get_unit_operation_parameters("solar_energy")

        assert data == seto_db._cached_files["solar_energy"]["default"]

        # Check for a few expected keys to check what we got back looks right
        assert "removal_frac_mass_comp" not in data
        assert "capital_cost" in data
        assert "reference_state" in data["capital_cost"]

    @pytest.mark.unit
    def test_get_unit_operation_parameters_invalid_subtype(self, seto_db):
        with pytest.raises(
            KeyError,
            match="Received unrecognised subtype foobar for "
            "technology solar_energy.",
        ):
            seto_db.get_unit_operation_parameters("solar_energy", subtype="foobar")

    @pytest.mark.unit
    def test_get_unit_operation_parameters_single_subtype(self, seto_db):

        seto_db._cached_files["solar_energy"]["coal"] = {
            "solar_radiation": "sufficient",
            "new_param": True,
        }

        data = seto_db.get_unit_operation_parameters("solar_energy", subtype="coal")

        for k, v in data.items():
            if k == "solar_radiation":
                assert v == "sufficient"
            elif k == "new_param":
                assert v is True
            else:
                assert v == seto_db._cached_files["solar_energy"]["default"][k]

    @pytest.mark.unit
    def test_get_unit_operation_parameters_multi_subtype(self, seto_db):
        # First, insert some data for a 2nd subtype into nanofiltration entry
        seto_db._cached_files["solar_energy"]["natural_gas"] = {
            "solar_radiation": "still_sufficient",
            "new_param_2": False,
        }

        # Load data for subtype
        data = seto_db.get_unit_operation_parameters(
            "solar_energy", subtype="natural_gas"
        )

        # Check data
        for k, v in data.items():
            if k == "solar_radiation":
                assert v == "still_sufficient"
            elif k == "new_param":
                assert v is True
            elif k == "new_param_2":
                assert v is False
            else:
                # All other entries should match defaults
                assert v == seto_db._cached_files["solar_energy"]["default"][k]

    @pytest.mark.unit
    def test_get_unit_operation_parameters_subtype_argument(self, seto_db):
        with pytest.raises(
            TypeError,
            match="Unexpected type for subtype 400: must be " "string or list like.",
        ):
            seto_db.get_unit_operation_parameters("solar_energy", subtype=400)

    @pytest.mark.unit
    def test_get_solute_set_default(self, seto_db):
        comp_set = seto_db.get_solute_set()

        assert comp_set == [
            "boron",
            "bromide",
            "calcium",
            "chloride",
            "magnesium",
            "potassium",
            "sodium",
            "strontium",
            "sulfate",
            "tds",
            "tss",
        ]

    @pytest.mark.unit
    def test_get_solute_set_specified(self, seto_db):
        comp_set = seto_db.get_solute_set("seawater")

        assert comp_set == [
            "boron",
            "bromide",
            "calcium",
            "chloride",
            "magnesium",
            "potassium",
            "sodium",
            "strontium",
            "sulfate",
            "tds",
            "tss",
        ]

    @pytest.mark.unit
    def test_get_solute_set_no_default(self, seto_db):
        # First, delete default entry from database
        del seto_db._cached_files["water_sources"]["default"]

        with pytest.raises(
            KeyError,
            match="Database has not defined a default water "
            "source and none was provided.",
        ):
            seto_db.get_solute_set()

    @pytest.mark.unit
    def test_flush_cache(self, seto_db):
        seto_db.flush_cache()

        assert seto_db._cached_files == {}
