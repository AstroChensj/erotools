#!/usr/bin/env python3
"""
Script to look up sky tiles for eROSITA.

Credits: Jeremy Sanders (2022)
"""
import sys
import os.path
import argparse

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from erotools.radec2tile import tile_local



def degreesIfNumeric(val):
    """Add a d suffix for degrees if numeric."""
    try:
        float(val)
    except ValueError:
        return val
    return val + 'd'

errormessage=f'''Usage: {sys.argv[0]} [frame] lon lat

frame is optional, and includes: icrs, fk5, fk4, galactic, (bary|geo)centric(mean|true)ecliptic
lon is longitude (e.g. 1.23 [degrees assumed], 6.0d or 1h5m4.3s)
lat is latitude (e.g. -3.12 [degrees assumed], -3.12d or -3d5m2.3s)
'''

def main():
    """Call the function to lookup the tile."""

    if len(sys.argv) not in (3, 4):
        # error if number of arguments is not right
        print(errormessage, file=sys.stderr)
        sys.exit(1)

    # read arguments (not using argparse as it doesn't like -ve numbers)
    args = sys.argv[1:]
    if len(args) == 3:
        frame = args[0]
        args = args[1:]
    else:
        frame = 'icrs'
    lon = degreesIfNumeric(args[0])
    lat = degreesIfNumeric(args[1])
 
    # do frame conversion using astropy.SkyCoord
    coord = SkyCoord(lon, lat, frame=frame)
    coord_icrs = coord.transform_to(frame='icrs')
    coord_gal = coord.transform_to(frame='galactic')

    print('RA,Dec (ICRS): ', coord_icrs.to_string())
    print('l,b:           ', coord_gal.to_string())

    # get scalar decimal degrees of coordinate
    ra = coord_icrs.ra.to_value(u.deg)
    dec = coord_icrs.dec.to_value(u.deg)

    # do the work
    tile = tile_local(ra,dec)

    # print the results
    print('Sky tile:      ', tile)

if __name__ == '__main__':
    main()