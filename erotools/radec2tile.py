#!/usr/bin/env python3
"""
Get eROSITA:DE skytile. This can either be done locally (`tile_local`) or through eROSITA API (`tile_api`).

Credits to: J. Sanders (2022)
"""
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import requests
import os


skymapfilename = os.path.join(os.path.dirname(__file__),"SKYMAPS_052022_MPE.fits")


def tile_local(ra,dec,chatter=0):
    """
    Return the sky tile for a given ra, dec (decimal). Returns only MPE or MPE+IKI skytile.

    Parameters
    ----------
    ra : float
        Right Ascension in decimal degrees.
    dec : float 
        Declination in decimal degrees.
    chatter : int
        Level of output log. Print output if >0.

    Returns
    -------
    The resulting skytile (None if not found).
    """

    # get list of tiles and filter for DE ownership
    with fits.open(skymapfilename) as skyf:
        d = skyf["SMAPS"].data
        sky_nr = d["SRVMAP"]
        sky_owner = d["OWNER"]
        sky_ra_min = d["RA_MIN"]
        sky_ra_max = d["RA_MAX"]
        sky_de_min = d["DE_MIN"]
        sky_de_max = d["DE_MAX"]

    # get index in table for tile
    matches = np.where(
        (ra >= sky_ra_min) & (ra <= sky_ra_max) &
        (dec >= sky_de_min) & (dec <= sky_de_max))[0]
    if len(matches) != 1:
        for idx in matches: print(sky_nr[idx])
        # raise RuntimeError("No tile or multiple tiles found")
        if chatter>0:
            print("No tile or multiple tiles found")
        return None
    idx = matches[0]
    skytile = f"{sky_nr[idx]:06d}"
    owner_tile = {0: "Both", 1: "IKI", 2: "MPE"}[sky_owner[idx]]

    if owner_tile not in ["Both","MPE"]:
        if chatter>0:
            print("Target source not in eROSITA:DE sky!")
        return None
    
    else:
        return skytile



def tile_api(ra,dec,radius=0):
    """
    Query the eROSITA sky tiles for a given sky position and radius.

    Parameters
    ----------
    ra : float
        Right Ascension in decimal degrees.
    dec : float 
        Declination in decimal degrees.
    radius : float, optional
        Search radius in decimal degrees. Default is 0 (point).

    Returns
    -------
    The resulting skytile (None if not found).
    """
    # Construct the API URL
    base_url = "https://erosita.mpe.mpg.de/dr1/erodat/skyview/skytile_search_api/"
    url = f"{base_url}?RA={ra}&DEC={dec}&RAD={radius}"
    try:
        # Query the API
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        # Parse the JSON response
        data = response.json()
        # Check if the response contains tiles
        if "tiles" in data:
            de_sky = data["tiles"][0]["de_sky"]
            if de_sky:
                return f"{data["tiles"][0]["srvmap"]:06d}"
            else:
                print("Target source not in eROSITA:DE sky!")
                return None
        else:
            print("No tiles found for the given position and radius.")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Error querying the API: {e}")
        return None