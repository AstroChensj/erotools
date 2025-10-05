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
from astropy.coordinates import search_around_sky,SkyCoord
import astropy.units as u
import pandas as pd
import subprocess
import sys
import os
import argparse
import logging
import time
from erotools.erocat import fake_srclist,search_catalog
from erotools.radec2tile import tile_local


# define argparser
class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write("error: %s\n" % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Shi-Jiang Chen, Johannes Buchner and Teng Liu (C) 2025 <JohnnyCsj666@gmail.com>""",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("science_evt",type=str,help="the input image/events file, from eRASS archive")
parser.add_argument("target_ra",type=float,help="target source RA [degree, icrs]")
parser.add_argument("target_dec",type=float,help="target source DEC [degree, icrs]")
parser.add_argument("--detprefix",type=str,default=None,help="`prefix` to be passed to `erosrcdet`")
parser.add_argument("--detsuffix",type=str,default="",help="`suffix` to be passed to `erosrcdet`")
parser.add_argument("--prefix",type=str,default="./out/forcedspec_",help="prefix for all products, defaults to './out/forcedspec_', will create a directory if necessary")
parser.add_argument("--suffix",type=str,default="",help="suffix for all products, defaults to ''")
parser.add_argument("--R_match_optical",type=float,default=3,help="matching radius [arcsec] between target and CTP catalog, defaults to 3")
parser.add_argument("--R_match_xray",type=float,default=15,help="matching radius [arcsec] between target and dr1 catalog, defaults to 15 (recommended)")
parser.add_argument("--R_confusion",type=float,default=60,help="confusion radius [arcsec], where sources in annulus of R_match ~ R_confusion lead to confusion issues, defaults to 60 (recommended)")
parser.add_argument("--eRASS_CAT_DIR",type=str,default=None,help="overwrite environmental variable eRASS_CAT_DIR? Defaults to None")
parser.add_argument("--redshift",type=float,default=0,help="source redshift")
parser.add_argument("--record_redshift",action="store_true",help="make redshift file")
parser.add_argument("--galnh",type=float,default=0,help="source galactic nh [1 cm^-2]")
parser.add_argument("--record_galnh",action="store_true",help="make galnh file")
args = parser.parse_args()


# make output directory if necessary
outdir = os.path.dirname(args.prefix)
os.makedirs(outdir,exist_ok=True)
log_dir = f"{outdir}/log"
os.makedirs(log_dir,exist_ok=True)


# define logger
logger = logging.getLogger("forcedspec")
logger.setLevel(logging.DEBUG)
logname = f"{log_dir}/{os.path.basename(args.prefix)}forcedspec{args.suffix}.log"
os.system(f"rm -rf {logname}")
file_handler = logging.FileHandler(logname)
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def run_evtool(science_evt,science_img):
	"""
	Run evtool if necessary.

	Parameters
	----------
	science_evt : str
		Input science events file.
	science_img : str
		Output science image file.

	Returns
	-------
	None
	"""
	img_cmd = [
		"evtool",
		f"eventfiles={science_evt}",
		f"outfile={science_img}",
		"image=yes",
		"emin=0.2",
		"emax=2.3",
	]
	logger.info(" ".join(img_cmd))
	mkimg_log = f"{log_dir}/{os.path.basename(args.prefix)}mkimg{args.suffix}.log"
	with open(mkimg_log,"w") as log_file:
		subprocess.run(img_cmd,stdout=log_file)

	return


def main():
	t0 = time.time()
	# check if eRASS_CAT_DIR exists
	logger.info("Checking if eRASS_CAT_DIR exists ...")
	if args.eRASS_CAT_DIR is not None:
		eRASS_CAT_DIR = args.eRASS_CAT_DIR
	else:
		eRASS_CAT_DIR = os.environ.get("eRASS_CAT_DIR") # remember to set the eRASS_CAT_DIR environmental variable!
	main_cat = f"{eRASS_CAT_DIR}/eRASS1_Main.v1.1.fits"
	supp_cat = f"{eRASS_CAT_DIR}/eRASS1_Supp.v1.1.fits"
	if not (os.path.exists(main_cat) and os.path.exists(supp_cat)):
		logger.error(f"{main_cat} or {supp_cat} does not exist!")
		raise


	# check if in eRO:DE sky
	logger.info("Checking if in eROSITA:DE sky ...")
	skytile = tile_local(args.target_ra,args.target_dec)	# point search
	if skytile is None: # not in eRO-DE sky
		logger.error(f"Target source (RA={args.target_ra}, DEC={args.target_dec}) not in eROSITA:DE sky!")
		raise


	# generate a fake box list (in catprep standards)
	box_name = f"{outdir}/fake_srclist.fits"
	science_img = f"{args.detprefix}sciimg{args.detsuffix}.fits"
	if not os.path.exists(science_img):
		science_img = f"{args.prefix}sciimg{args.suffix}.fits"
		if not os.path.exists(science_img):
			run_evtool(args.science_evt,science_img)
	logger.info(f"Generating a fake box list ({box_name}) ...")
	fakesrc_log = f"{log_dir}/{os.path.basename(args.prefix)}fakesrc{args.suffix}.log"
	with open(fakesrc_log,"w") as log_file:
		sys.stdout = log_file
		sys.stderr = log_file
		fake_srclist(skytile,box_name,science_img,style="cat")	# NOTE: need to input image rather than events file!
		sys.stdout = sys.__stdout__
		sys.stderr = sys.__stderr__


	# update the boxlist: add `AUTO_EXTRACT`, `AUTO_EXCLUDE` columns
	logger.info("Updating the boxlist with `AUTO_EXTRACT`, `AUTO_EXCLUDE` ...")
	# find the closest match from eRASS1 MAIN+SUPP catalog
	logger.info("Looking for closest match from eRASS1 MAIN+SUPP catalog ...")
	entry = search_catalog(args.target_ra,args.target_dec,args.R_match_optical,args.R_match_xray)

	if entry is not None: # detection
		id_target = entry["ID_SRC"]
		logger.info(f"Found. ID_SRC is {id_target}. Updating the boxlist ...")
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
		logger.info("Non-detection. Running forced photometry first ...")
		# TODO: need to modify here!
		forcedphot_cmd = [
			"eroforcedphot",
			f"{args.science_evt}",
			f"{args.target_ra}",
			f"{args.target_dec}",
			"--emin", "0.2",
			"--emax", "2.3",
			"--target_z", "0.",
			"--detprefix",f"{args.detprefix}",
			"--detsuffix",f"{args.detsuffix}",
			"--skip_exist_srcdet",
			"--prefix",f"{outdir}/forcedphot/{os.path.basename(args.prefix)}",
			"--suffix",f"{args.suffix}",
			"--R_match", f"{args.R_match_xray}",
			"--R_confusion", f"{args.R_confusion}",
		]
		logger.info(" ".join(forcedphot_cmd))
		forcedphot_log = f"{log_dir}/{os.path.basename(args.prefix)}forcedphot{args.suffix}.log"
		with open(forcedphot_log,"w") as log_file:
			subprocess.run(forcedphot_cmd,stdout=log_file)
		logger.info("Done. ID_SRC will be 99999. Writing fake `ML_CTS_0`, `ML_BKG_0` ...")

		with fits.open(f"{outdir}/forcedphot/{os.path.basename(args.prefix)}apesummary{args.suffix}.fits") as hdu:
			forcedphot_data = hdu[1].data
		ape_eef = forcedphot_data["APE_EEF"]
		ape_radius = forcedphot_data["APE_RADIUS"] * 4 / 60	# 1 pixel = 4 arcsec; 1 arcsec = 1/60 arcmin
		scts = forcedphot_data["ME"] / ape_eef	# posterior median counts, bkg subtracted, corrected for psf loss
		scts_lo = forcedphot_data["LO"] / ape_eef
		scts_hi = forcedphot_data["HI"] / ape_eef
		bcts_arcmin = forcedphot_data["APE_BKG"] / (np.pi * ape_radius**2)	# bkg counts per arcmin^2
		ape_exp = forcedphot_data["APE_EXP"]	# aperture exposure [s], vignetting corrected
		cr_me = forcedphot_data["CR_ME"]	# source count rate, psf loss corrected, vignetting corrected
		cr_lo = forcedphot_data["CR_LO"]
		cr_hi = forcedphot_data["CR_HI"]
		flux_me = forcedphot_data["FLUX_ME"]	# source flux, psf loss corrected, vignetting corrected
		flux_lo = forcedphot_data["FLUX_LO"]
		flux_hi = forcedphot_data["FLUX_HI"]
		detlike = forcedphot_data["like"]
		# convert ra, dec to lii, bii
		coord = SkyCoord(ra=args.target_ra*u.degree,dec=args.target_dec*u.degree,frame="icrs")
		galactic = coord.galactic
		lii, bii = galactic.l.degree, galactic.b.degree
		# append our target source to existing box list
		with fits.open(box_name,mode="update") as hdu:
			box_data = hdu[1].data
			target_row = np.zeros(1,dtype=box_data.dtype)
			target_row["ID_SRC"] = "99999"	# for undetected source
			target_row["RA"] = args.target_ra
			target_row["RA_LOWERR"] = 3.72	# average value for eRASS1 MAIN catalog; doesn't matter here
			target_row["RA_UPERR"] = 3.72
			target_row["DEC"] = args.target_dec
			target_row["DEC_LOWERR"] = 3.72
			target_row["DEC_UPERR"] = 3.72
			target_row["RADEC_ERR"] = 5.31
			target_row["LII"] = lii
			target_row["BII"] = bii
			target_row["EXT"] = 0
			target_row["EXT_ERR"] = 0
			target_row["EXT_LOWERR"] = 0
			target_row["EXT_UPERR"] = 0
			target_row["EXT_LIKE"] = 0
			target_row["ML_CTS_0"] = scts
			target_row["ML_CTS_ERR_0"] = scts_lo + scts_hi
			target_row["ML_CTS_LOWERR_0"] = scts_lo
			target_row["ML_CTS_UPERR_0"] = scts_hi
			target_row["ML_RATE_0"] = cr_me		# actually not used; only to keep in catprep format
			target_row["ML_RATE_ERR_0"] = cr_lo + cr_hi
			target_row["ML_RATE_LOWERR_0"] = cr_lo
			target_row["ML_RATE_UPERR_0"] = cr_hi
			target_row["ML_FLUX_0"] = flux_me	# actually not used; only to keep in catprep format
			target_row["ML_FLUX_ERR_0"] = flux_lo + flux_hi
			target_row["ML_FLUX_LOWERR_0"] = flux_lo
			target_row["ML_FLUX_UPERR_0"] = flux_hi
			target_row["DET_LIKE_0"] = detlike
			target_row["ML_BKG_0"] = bcts_arcmin	# bkg counts per acmin^2
			target_row["ML_EXP_0"] = ape_exp
			target_row["ML_EEF_0"] = ape_eef
			box_data_updated = np.append(box_data,target_row)
			hdu[1].data = box_data_updated
		# and then add AUTO_EXTRACT, AUTO_EXCLUDE
		logger.info("Updating boxlist ...")
		with fits.open(box_name,mode="update") as hdu:
			box = hdu[1]
			ids = box.data["ID_SRC"]
			id_target = 99999	# ids are integers, not str
			auto_extract = np.array([1 if ids[i]==id_target else 0 for i in range(len(ids))])
			auto_exclude = np.array([0 if ids[i]==id_target else 1 for i in range(len(ids))])
			if "AUTO_EXTRACT" in box.columns.names:
				box.columns.del_col("AUTO_EXTRACT")
			if "AUTO_EXCLUDE" in box.columns.names:
				box.columns.del_col("AUTO_EXCLUDE")
			AUTO_EXTRACT = fits.Column(name="AUTO_EXTRACT",format="I",array=auto_extract)
			AUTO_EXCLUDE = fits.Column(name="AUTO_EXCLUDE",format="I",array=auto_exclude)
			box.data = fits.BinTableHDU.from_columns(box.columns+AUTO_EXTRACT+AUTO_EXCLUDE).data
			
	
	# run srctool
	logger.info("Running srctool ...")
	outsrcreg = f"{args.prefix}src{args.suffix}.reg"
	outbackreg = f"{args.prefix}bkg{args.suffix}.reg"
	srctool_cmd = [
		"srctool",
		f"eventfiles={args.science_evt}",
		f"srccoord={box_name}",
		f"prefix={args.prefix}",
		f"suffix={args.suffix}",
		f"todo=SPEC ARF RMF",	# do not add "" !
		f"insts=1 2 3 4 6",		# skip TM 5 and 7 for potential light leak
		f"srcreg=AUTO",
		f"backreg=AUTO",
		f"outsrcreg={outsrcreg}",
		f"outbackreg={outbackreg}",
		# f"tstep=0.05",
		# f"xgrid=1.0 2.0",
		f"clobber=yes",
	]
	logger.info(" ".join(srctool_cmd))
	srctool_log = f"{log_dir}/{os.path.basename(args.prefix)}srctool{args.suffix}.log"
	with open(srctool_log,"w") as log_file:
		subprocess.run(srctool_cmd,stdout=log_file)

		
	# make redshift file (optional)
	if args.record_redshift:
		zfile = f"{args.prefix}020_SourceSpec_{id_target:05d}{args.suffix}.fits.z"
		with open(zfile,"w") as f:
			f.writelines(f"{args.redshift}")


	# make nh file (optional)
	if args.record_galnh:
		nhfile = f"{args.prefix}020_SourceSpec_{id_target:05d}{args.suffix}.fits.nh"
		with open(nhfile,"w") as f:
			f.writelines(f"{args.galnh}")


	# make a link to science image (for later visual check)
	link_name = f"{args.prefix}sciimg{args.suffix}.fits"
	try:
		os.symlink(science_img,link_name)
	except FileExistsError:
		os.remove(link_name)
		os.symlink(science_img,link_name)


	t1 = time.time()
	logger.info(f"Total time use: {t1-t0} s")



if __name__ == "__main__":
	main()

