#!/usr/bin/env python3
"""
Perform forced spectroscopy on target position.

Given the target position (RA&DEC), first look for it from the eRASS1 MAIN+SUPP catalog. 
* If detection, feed `srctool` with the skytile catalog, and generate SPEC, BKGSPEC, ARF, RMF under `AUTO` mode.
* If non-detection: first perform `apetool` to estimate `APE_CTS`, `APE_BKG`, convert to `ML_CTS`, `ML_BKG`, and add to existing skytile catalog; then feed `srctool` with the new skytile catalog.

In both cases, we extract only the spectral/response files for the target source. To do that, we set `AUTO_EXTRACT`=1, `AUTO_EXCLUDE`=0 for the target source, while `AUTO_EXTRACT`=0, `AUTO_EXCLUDE`=1 for the rest sources (we view them as potential contamination to our target area) in the skytile catalog. 

Output: SPEC, BKGSPEC, ARF, RMF.

"""
import numpy as np
from astropy.io import fits
# from astropy.coordinates import search_around_sky, SkyCoord
import astropy.units as u
import pandas as pd
import subprocess
import sys
import os
import argparse
import logging
from erotools.erocat import find_erode_skytile,fake_srclist,search_catalog,look_for_confusion,create_aperture
# from erotools.eroecf import get_eroecf
# from erotools.detlike import cal_detlike,counts_twoside_poisson,counts_ul_poisson



# define logger
logger = logging.getLogger("forcedspec")
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler("forcedspec.log")
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


# define argparser
class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write("error: %s\n" % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Shi-Jiang Chen, Johannes Buchner and Teng Liu (C) 2024 <JohnnyCsj666@gmail.com>""",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("science_img",type=str,help="the input image/events file")
parser.add_argument("target_ra",type=float,help="target source RA [degree, icrs]")
parser.add_argument("target_dec",type=float,help="target source DEC [degree, icrs]")
# parser.add_argument("--emin",type=float,default=0.2,help="rest-frame minimum energy [keV]")
# parser.add_argument("--emax",type=float,default=2.3,help="rest-frame maximum energy [keV]")
parser.add_argument("--target_z",type=float,default=0,help="target redshift (real photometry extracted from emin/z, emax/z)")
parser.add_argument("--outdir",type=str,default="outdir",help="output directory name")
parser.add_argument("--outname",type=str,default="result.fits",help="name of fits file storing the aperture photometry results, saved under outdir")
parser.add_argument("--prefix",type=str,default="",help="prefix for all products")
parser.add_argument("--suffix",type=str,default="",help="suffix for all products")
parser.add_argument("--R_match",type=float,default=15,help="matching radius (arcsec) between target and dr1 catalog (recommend: 15)")
parser.add_argument("--R_confusion",type=float,default=60,help="confusion radius (arcsec), where sources in annulus of R_match ~ R_confusion lead to confusion issues (recommend: 60)")
parser.add_argument("--eRASS_CAT_DIR",type=str,default=None,help="overwrite environmental variable eRASS_CAT_DIR")
args = parser.parse_args()


if args.eRASS_CAT_DIR is not None:
	eRASS_CAT_DIR = args.eRASS_CAT_DIR
else:
	eRASS_CAT_DIR = os.environ.get("eRASS_CAT_DIR") # remember to set the eRASS_CAT_DIR environmental variable!
main_cat = f"{eRASS_CAT_DIR}/eRASS1_Main.v1.1.fits"
supp_cat = f"{eRASS_CAT_DIR}/eRASS1_Supp.v1.1.fits"
if not (os.path.exists(main_cat) and os.path.exists(supp_cat)):
	logger.error(f"{main_cat} or {supp_cat} does not exist!")
	raise



skytile = find_erode_skytile(args.ra,args.dec,radius=0)	# point search
if skytile is None: # not in eRO-DE sky
	logger.error(f"Target source (RA={args.target_ra}, DEC={args.target_dec}) not in eROSITA:DE sky!")
	raise

# generate a fake box list (in catprep standards)
box_name = f"{args.outdir}/fake_srclist.fits"
fake_srclist(args.target_ra,args.target_dec,box_name,args.science_img,skytile=skytile,style="cat")



# update the boxlist: add `AUTO_EXTRACT`, `AUTO_EXCLUDE` columns
# find the closest match from eRASS1 MAIN+SUPP catalog
entry = search_catalog(args.target_ra,args.target_dec,args.R_match)
if entry is not None: # detection
	id_target = entry["ID_SRC"]
	with fits.open(box_name,mode="update") as hdu:
		box = hdu[1]
		ids = box.data["ID_SRC"]
		auto_extract = np.array([1 if ids[i]==id_target else 0 for i in range(len(ids))])
		auto_exclude = np.array([0 if ids[i]==id_target else 1 for i in range(len(ids))])
		if "AUTO_EXTRACT" in box.columns.names:
			box.columns.del_col("AUTO_EXTRACT")
		if "AUTO_EXCLUDE" in box.columns.names:
			box.columns.del_col("AUTO_EXCLUDE")
		AUTO_EXTRACT = fits.Column(name="AUTO_EXTRACT",format="I",array=auto_extract)
		AUTO_EXCLUDE = fits.Column(name="AUTO_EXCLUDE",format="I",array=auto_exclude)
		box.data = fits.BinTableHDU.from_columns(box.columns+AUTO_EXTRACT+AUTO_EXCLUDE).data
else:	# non-detection
	subprocess.run([
		"eroforcedphot",
		f"{args.science_img}",
        f"{args.target_ra}",
        f"{args.target_dec}",
        "--emin", "0.2",
        "--emax", "2.3",
        "--target_z", "0.",
        "--outdir", f"{args.outdir}/forcedphot",
		"--R_match", f"{args.R_match}",
		"--R_confusion", f"{args.R_confusion}",
	])
	with fits.open(f"{args.outdir}/forcedphot/result.fits") as hdu:
		forcedphot_data = hdu[1].data
	ape_eef = forcedphot_data["APE_EEF"]
	ape_radius = forcedphot_data["APE_RADIUS"] * 4 / 60	# 1 pixel = 4 arcsec; 1 arcsec = 1/60 arcmin
	scts = forcedphot_data["ME"] / ape_eef	# posterior median counts, corrected for psf loss
	bcts_arcmin = forcedphot_data["APE_BKG"] / (np.pi * ape_radius**2)	# bkg counts per arcmin^2
	
	with fits.open(box_name,mode="updatee") as hdu:
		box = hdu[1]
		box_data = box.data
		target_row = np.copy(box_data[0])
		target_row["ID_SRC"] = "99999"	# for undetected source
		target_row["RA"] = args.target_ra
		target_row["RA_LOWERR"] = 3.72	# average value for eRASS1 MAIN catalog; doesn't matter here
		target_row["RA_UPERR"] = 3.72
		target_row["DEC"] = args.target_dec
		target_row["DEC_LOWERR"] = 3.72
		target_row["DEC_UPERR"] = 3.72
		target_row["RADEC_ERR"] = 5.31
		# target_row["LII"] = 
		# target_row["BII"] = 
		target_row["EXT"] = 0
		target_row["EXT_ERR"] = 0
		target_row["EXT_LOWERR"] = 0
		target_row["EXT_UPERR"] = 0
		target_row["EXT_LIKE"] = 0
		target_row[""]

		

	pass

