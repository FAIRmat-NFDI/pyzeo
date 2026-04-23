import pytest
from pyzeo.netstorage import AtomNetwork

@pytest.mark.parametrize("filename, read_kwargs, read_method", [
    # Custom radii scenario
    ("EDI.cif", {"rad_file": "SiO.mass"}, AtomNetwork.read_from_CIF),
    ("EDI.v1", {"rad_file": "SiO.mass"}, AtomNetwork.read_from_V1),
    
    # Default radii scenario
    ("EDI.cif", {}, AtomNetwork.read_from_CIF),
    ("EDI.v1", {}, AtomNetwork.read_from_V1),
])
def test_read_file(filename, read_kwargs, read_method):
    """
    Tests whether an AtomNetwork object can be successfully
    created from some file formats under different
    radius initialization conditions.
    """
    atmnet = read_method(filename, **read_kwargs)

    assert atmnet is not None
