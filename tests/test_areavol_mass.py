import pytest
import re
import math
import os
from pyzeo.netstorage import AtomNetwork
from pyzeo.area_volume import volume, surface_area
from pyzeo.high_accuracy import high_accuracy_atomnet

@pytest.mark.parametrize("read_kwargs, expected_values", [
    # No mass scenario
    (
        {"mass_flag": False},
        {
            "density": 0.0,
            "av_cm3_g": math.inf,
            "nav_cm3_g": math.nan,
         }
        ),
    
    # Default mass scenario
    (
        {},
        {
            "density": 1.62239,
            "av_cm3_g": 0.0454022,
            "nav_cm3_g": 0.0,
         }
        ),
    
    # Custom mass file scenario
    (
        {"mass_file": "SiO.mass"},
        {
            "density": 1.53911,
            "av_cm3_g": 0.0478587,
            "nav_cm3_g": 0.0,
         }
        ),

])
def test_areavol_mass(data_dir, read_kwargs, expected_values):
    """
    Test that the volume calculations computes the correct
    mass-related results under different mass initialization conditions.
    """

    structure_path = os.path.join(data_dir, "EDI.cssr")

    current_kwargs = read_kwargs.copy()
    if "mass_file" in current_kwargs:
        current_kwargs["mass_file"] = os.path.join(data_dir, current_kwargs["mass_file"])

    atmnet = AtomNetwork.read_from_CSSR(structure_path, **current_kwargs)
    ha_atmnet = atmnet.copy()
    high_accuracy_atomnet(ha_atmnet, "DEF")

    vol_str = volume(
            atmnet,
            1.2,
            1.2,
            50000,
            high_accuracy=True,
            high_accuracy_atmnet=ha_atmnet
            )

    decoded_str = vol_str.decode("utf-8")

    match_density = re.search(r"Density:\s*([^\s]+)", decoded_str)
    match_av = re.search(r"AV_cm\^3/g:\s*([^\s]+)", decoded_str)
    match_nav = re.search(r"NAV_cm\^3/g:\s*([^\s]+)", decoded_str)

    assert match_density, "Could not find 'Density:' in output."
    assert match_av, "Could not find 'AV_cm^3/g:' in output."
    assert match_nav, "Could not find 'NAV_cm^3/g:' in output."
    
    actual_density = float(match_density.group(1))
    actual_av = float(match_av.group(1))
    actual_nav = float(match_nav.group(1))

    assert actual_density == pytest.approx(expected_values["density"], rel=1e-3), \
        f"Density failed. Expected {expected_values['density']}, got {actual_density}"
        
    assert actual_av == pytest.approx(expected_values["av_cm3_g"], rel=1e-3, abs=1e-6, nan_ok=True), \
        f"AV_cm^3/g failed. Expected {expected_values['av_cm3_g']}, got {actual_av}"
        
    assert actual_nav == pytest.approx(expected_values["nav_cm3_g"], rel=1e-3, abs=1e-6, nan_ok=True), \
        f"NAV_cm^3/g failed. Expected {expected_values['nav_cm3_g']}, got {actual_nav}"
