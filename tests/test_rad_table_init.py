import pytest
import os
from pyzeo.netstorage import AtomNetwork
from pyzeo.extension import lookupRadius

@pytest.mark.parametrize("filename, read_kwargs, expected_value", [
 
    (
        "MgO_vac1.cssr",
        {"rad_file": "MgO.rad"},
        1.84,
        ),

    (
        "EDI.cssr",
        {},
        1.52,
        ),

])
def test_rad_table_init(data_dir, filename, read_kwargs, expected_value):
    """
    Test for global state pollution of the atomic radii table
    when running sequential calculations with different
    atomic radii for the same atom types.
    """

    if filename == "EDI.cssr":
        data_path = os.path.join(data_dir, filename)

    else:
        data_path = filename

    atmnet = AtomNetwork.read_from_CSSR(data_path, **read_kwargs)

    element = "O"
    encoded_str = element.encode("utf-8")

    element_radius = lookupRadius(encoded_str)

    assert element_radius == expected_value, \
            f"Wrong element radius. Expected {expected_value}, got {element_radius}"

