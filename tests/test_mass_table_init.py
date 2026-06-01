import pytest
import os
from pyzeo.netstorage import AtomNetwork
from pyzeo.extension import lookupMass

@pytest.mark.parametrize("filename, read_kwargs, expected_value", [
 
    (
        "EDI.cssr",
        {"mass_file": "SiO.mass"},
        15,
        ),

    (
        "EDI.cssr",
        {},
        15.999,
        ),

])
def test_mass_table_init(data_dir, filename, read_kwargs, expected_value):
    """
    Test for global state pollution of the atomic masses table
    when running sequential calculations with different
    atomic masses for the same atom types.
    """
    structure_path = os.path.join(data_dir, filename)

    current_kwargs = read_kwargs.copy()
    if "mass_file" in current_kwargs:
        current_kwargs["mass_file"] = os.path.join(data_dir, current_kwargs["mass_file"])

    atmnet = AtomNetwork.read_from_CSSR(structure_path, **current_kwargs)

    element = "O"
    encoded_str = element.encode("utf-8")

    element_mass = lookupMass(encoded_str)

    assert element_mass == expected_value, \
            f"Wrong element mass. Expected {expected_value}, got {element_mass}"

