"""
This is an example script using run_apetool
"""
import numpy as np
import sys
import os
import glob
import subprocess
from joblib import Parallel, delayed
from tqdm import tqdm
from erotools.erocat import find_erode_skytile


filelst = sys.argv[1]
ero_archive_dir = '/data3/public/data/eROSITA/ero_archive'
usecpu = 10


def parse_sdss_name(sdss_name):
    """
    Extract RA and Dec from an SDSS object name.

    Parameters
    ----------
    sdss_name : str 
        SDSS name in the format 'Jhhmmss.ss±ddmmss.s'.
    
    Returns
    -------
    tuple: (RA in degrees, Dec in degrees)
    """
    # Remove the leading 'J'
    coords = sdss_name[1:]
    
    # Determine the split character (+ or -)
    split_char = '+' if '+' in coords else '-'
    
    # Split into RA and Dec parts
    ra_part, dec_part = coords.split(split_char, 1)
    dec_part = split_char + dec_part  # Add back the sign to Dec

    # Parse RA
    ra_h = int(ra_part[:2])
    ra_m = int(ra_part[2:4])
    ra_s = float(ra_part[4:])
    ra_deg = (ra_h + ra_m / 60 + ra_s / 3600) * 15

    # Parse Dec
    dec_sign = -1 if dec_part[0] == '-' else 1
    dec_d = int(dec_part[1:3])
    dec_m = int(dec_part[3:5])
    dec_s = float(dec_part[5:])
    dec_deg = dec_sign * (dec_d + dec_m / 60 + dec_s / 3600)

    return ra_deg, dec_deg


def process_entry(i):
    sdssname = sdssname_lst[i]
    ra = ra_lst[i]
    dec = dec_lst[i]
    z = z_lst[i]
    RAtile = RAtile_lst[i]
    DEtile = DEtile_lst[i]
    skytile = skytile_lst[i]

    os.system("mkdir -p %s"%(sdssname))
    img = glob.glob("%s/%s/%s/EXP_010/e[mnb]01_%s_020_EventList_c010.fits.gz"%(ero_archive_dir,DEtile,RAtile,skytile))[0]
    os.system("cp %s %s"%(img,sdssname))
    img_new = glob.glob("%s/e[mnb]01_%s_020_EventList_c010.fits.gz"%(sdssname,skytile))[0]

    print(img_new)

    output = subprocess.run([
        'run_apetool',
        '%s'%(img_new),
        '%f'%(ra),
        '%f'%(dec),
        '--emin', '0.5',
        '--emax', '8.0',
        '--target_z', '%f'%(z),
        '--outdir', '%s/ape_005_080'%(sdssname),
    ],capture_output=True,text=True)
    with open("%s/ape_005_080/ape_005_080.log"%(sdssname),"w") as log_file:
        log_file.write(output.stdout)

    output = subprocess.run([
        'run_apetool',
        '%s'%(img_new),
        '%f'%(ra),
        '%f'%(dec),
        '--emin', '2.0',
        '--emax', '10.0',
        '--target_z', '%f'%(z),
        '--outdir', '%s/ape_020_100'%(sdssname),
    ],capture_output=True,text=True)
    with open("%s/ape_020_100/ape_020_100.log"%(sdssname),"w") as log_file:
        log_file.write(output.stdout)

    return


sdssname_lst = []
ra_lst = []
dec_lst = []
z_lst = []
skytile_lst = []
RAtile_lst = []
DEtile_lst = []

with open(filelst,'r') as f:
    lines = f.readlines()

for line in lines:
    sdssname = line.strip().split('\t')[0]
    z = float(line.strip().split('\t')[1])
    ra,dec = parse_sdss_name(sdssname)
    skytile = find_erode_skytile(ra,dec)
    RAtile = skytile[:3]
    DEtile = skytile[3:]

    sdssname_lst.append(sdssname)
    ra_lst.append(ra)
    dec_lst.append(dec)
    z_lst.append(z)
    skytile_lst.append(skytile)
    RAtile_lst.append(RAtile)
    DEtile_lst.append(DEtile)


Parallel(n_jobs=usecpu)(delayed(process_entry)(i) for i in tqdm(range(len(sdssname_lst))))
