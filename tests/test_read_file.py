import pytest
import os
from pyzeo.netstorage import AtomNetwork

@pytest.mark.parametrize("filename, read_kwargs, read_method", [
    # Custom radii scenario
    ("EDI.cif", {"rad_file": "SiO.rad"}, AtomNetwork.read_from_CIF),
    ("EDI.v1", {"rad_file": "SiO.rad"}, AtomNetwork.read_from_V1),
    
    # Default radii scenario
    ("EDI.cif", {}, AtomNetwork.read_from_CIF),
    ("EDI.v1", {}, AtomNetwork.read_from_V1),
])
def test_read_file(data_dir, filename, read_kwargs, read_method):
    """
    Tests whether an AtomNetwork object can be successfully
    created from some file formats under different
    radius initialization conditions.
    """

    structure_path = os.path.join(data_dir, filename)

    current_kwargs = read_kwargs.copy()
    if "rad_file" in current_kwargs:
        current_kwargs["rad_file"] = os.path.join(data_dir, current_kwargs["rad_file"])

    atmnet = read_method(structure_path, **current_kwargs)

    assert atmnet is not None
