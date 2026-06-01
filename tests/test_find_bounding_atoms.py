import pytest
import os
from pyzeo.netstorage import AtomNetwork  

def test_find_bounding_atoms(data_dir):  
    """  
    Test the channel finding and the identification of atoms that
    bound channels. 
    """  

    data_path = os.path.join(data_dir, "EDI.cssr")

    atmnet = AtomNetwork.read_from_CSSR(data_path)
    probe_radius = 1.5 

    vornet, _, _ = atmnet.perform_voronoi_decomposition()
    channels, access_info = vornet.find_channels(probe_radius)
    num_accessible = access_info.count(True)

    print(f"Found {len(channels)} channels in the structure")  
    print(f"Accessibility info for {len(access_info)} nodes")
    print(f"Number of accessible nodes {num_accessible}")

    assert len(channels) == 1, f"Expected 1 channel, got {len(channels)}"
    assert len(access_info) == 70, f"Expected 70 Voronoi nodes, got {len(access_info)}"
    assert num_accessible == 25, f"Expected 25 accessible nodes, got {num_accessible}"

    channel = channels[0]
    atom_ids = channel.find_bounding_atoms(atmnet, vornet)  

    print("\nChannel 0:")  
    print(f"  Number of bounding atoms: {len(atom_ids)}")  
    print(f"  Bounding atom IDs: {atom_ids}")  

    expected_atom_ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    assert len(atom_ids) == 14, f"Expected 14 bounding atoms, got {len(atom_ids)}"
    assert atom_ids == expected_atom_ids, "Bounding atom IDs do not match expectations"

    expected_coords = {
            0: ("O", (1.315, 0.000, 2.276)),
            1: ("O", (0.000, 5.611, 4.134)),
            2: ("O", (0.000, 1.315, 4.134))
    }

    for atom_id in atom_ids[:3]: 
        atom = atmnet.atoms[atom_id]
        coords = atom.coords  

        print(f"    Atom {atom_id}: type={atom.type}, "  
              f"coords=({coords[0]:.3f}, {coords[1]:.3f}, {coords[2]:.3f})")  
        
        expected_type, (ex_x, ex_y, ex_z) = expected_coords[atom_id]
        assert atom.type == expected_type, f"Atom {atom_id} expected type {expected_type}"
        assert coords[0] == pytest.approx(ex_x, abs=1e-3)
        assert coords[1] == pytest.approx(ex_y, abs=1e-3)
        assert coords[2] == pytest.approx(ex_z, abs=1e-3)
