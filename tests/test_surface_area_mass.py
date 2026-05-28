import pytest
import re
import math
import os
from pyzeo.netstorage import AtomNetwork
from pyzeo.area_volume import surface_area
from pyzeo.high_accuracy import high_accuracy_atomnet

@pytest.mark.parametrize("read_kwargs, expected_values", [
    # No mass scenario
    (
        {"mass_flag": False},
        {
            "density": 0.0,
            "asa_m2_g": math.inf,
            "nasa_m2_g": math.nan,
         }
        ),
    
    # Default mass scenario
    (
        {},
        {
            "density": 1.62239,
            "asa_m2_g": 1218.21,
            "nasa_m2_g": 0.0,
         }
        ),
    
    # Custom mass file scenario
    (
        {"mass_file": "SiO.mass"},
        {
            "density": 1.53911,
            "asa_m2_g": 1284.12,
            "nasa_m2_g": 0.0,
         }
        ), 
])
def test_surface_area_mass(data_dir, read_kwargs, expected_values):
    """
    Test that the surface_area calculation computes the correct
    mass-related results under different mass initialization conditions.
    """

    structure_path = os.path.join(data_dir, "EDI.cif")

    current_kwargs = read_kwargs.copy()
    if "mass_file" in current_kwargs:
        current_kwargs["mass_file"] = os.path.join(data_dir, current_kwargs["mass_file"])

    atmnet = AtomNetwork.read_from_CIF(structure_path, **current_kwargs)
    ha_atmnet = atmnet.copy()
    high_accuracy_atomnet(ha_atmnet, "DEF")

    sa_str = surface_area(
        atmnet, 
        channel_radius=1.2, 
        probe_radius=1.2, 
        mc_sampling_no=2000, 
        high_accuracy=True,
        high_accuracy_atmnet=ha_atmnet)

    decoded_str = sa_str.decode("utf-8")

    match_density = re.search(r"Density:\s*([^\s]+)", decoded_str)
    match_asa = re.search(r"ASA_m\^2/g:\s*([^\s]+)", decoded_str)
    match_nasa = re.search(r"NASA_m\^2/g:\s*([^\s]+)", decoded_str)

    assert match_density, "Could not find 'Density:' in output."
    assert match_asa, "Could not find 'ASA_m^2/g:' in output."
    assert match_nasa, "Could not find 'NASA_m^2/g:' in output."
    
    actual_density = float(match_density.group(1))
    actual_asa = float(match_asa.group(1))
    actual_nasa = float(match_nasa.group(1))

    assert actual_density == pytest.approx(expected_values["density"], rel=1e-3), \
        f"Density failed. Expected {expected_values['density']}, got {actual_density}"
        
    assert actual_asa == pytest.approx(expected_values["asa_m2_g"], rel=1e-3, abs=1e-6, nan_ok=True), \
        f"ASA_m^2/g failed. Expected {expected_values['asa_m2_g']}, got {actual_asa}"
        
    assert actual_nasa == pytest.approx(expected_values["nasa_m2_g"], rel=1e-3, abs=1e-6, nan_ok=True), \
        f"NASA_m^2/g failed. Expected {expected_values['nasa_m2_g']}, got {actual_nasa}"
