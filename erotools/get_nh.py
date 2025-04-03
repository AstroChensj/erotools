import re
import os
import subprocess

def get_nh(RA,DEC):
    """
    Get the Galactic NH from NASA"s HEASARC tool `NH` (https://heasarc.gsfc.nasa.gov/Tools/w3nh_help.html).
    Please ensure the HEASOFT env has been set up.
    
    Parameters
    ----------
    RA : float
    DEC : float

    Returns
    -------
    nh_val : float
        nh values in units of 1 cm^-2
    """
    nh_script = f"""
    (
    echo 2000
    echo {RA}
    echo {DEC}
    ) | nh
    """
    result = subprocess.run(nh_script,shell=True,executable="/bin/bash",stdout=subprocess.PIPE)
    lines = result.stdout.decode("utf-8").splitlines()

    pattern1 = r"Weighted average nH \(cm\*\*-2\)\s+([0-9.E+-]+)"
    nh_val = [re.search(pattern1,line).group(1) for line in lines if re.search(pattern1,line)]
    if len(nh_val)>0:
        return float(nh_val[0])
    else:
        pattern2 = r"h1_nh_HI4PI.fits >> nH \(cm\*\*-2\)\s+([0-9.E+-]+)"
        nh_val = [re.search(pattern2,line).group(1) for line in lines if re.search(pattern2,line)]
        if len(nh_val) > 0:
            return float(nh_val[0])
        else:
            raise Exception(f"Invalid RA ({RA:.4f}), DEC({DEC:.4f}) for nh!")