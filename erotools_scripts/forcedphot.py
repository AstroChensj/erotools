#!/usr/bin/env python3
"""
Perform aperture photometry based on eRASS1 MAIN+SUPP catalog.

Output: source detection likelihood, count rate (median and upper limit), flux (median and upper limit).

"""
import numpy as np
from astropy.io import fits
import astropy.units as u
import pandas as pd
import subprocess
import sys
import os
import argparse
import logging
from erotools.erocat import fake_srclist,look_for_confusion,create_aperture
from erotools.eroecf import get_eroecf
from erotools.detlike import cal_detlike,counts_twoside_poisson,counts_ul_poisson
from erotools.radec2tile import tile_local


class HelpfulParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)

parser = HelpfulParser(description=__doc__,
    epilog="""Shi-Jiang Chen, Johannes Buchner and Teng Liu (C) 2025 <JohnnyCsj666@gmail.com>""",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("science_evt",type=str,help="the input image/events file")
parser.add_argument("target_ra",type=float,help="target source RA (degree, icrs)")
parser.add_argument("target_dec",type=float,help="target source DEC (degree, icrs)")
parser.add_argument("--emin",type=float,default=0.2,help="rest-frame minimum energy [keV]")
parser.add_argument("--emax",type=float,default=2.3,help="rest-frame maximum energy [keV]")
parser.add_argument("--target_z",type=float,default=0.,help="target redshift (real photometry extracted from emin/z, emax/z)")
parser.add_argument("--detprefix",type=str,default=None,help="`prefix` to be passed to `erosrcdet`")
parser.add_argument("--detsuffix",type=str,default="",help="`suffix` to be passed to `erosrcdet`")
parser.add_argument("--skip_exist_srcdet",action="store_true",help="`skip_exist` parameter to be passed to `erosrcdet`")
parser.add_argument("--prefix",type=str,default="",help="prefix for all products")
parser.add_argument("--suffix",type=str,default="",help="suffix for all products")
parser.add_argument("--R_match",type=float,default=15,help="matching radius (arcsec) between target and dr1 catalog (recommend: 15)")
parser.add_argument("--R_confusion",type=float,default=60,help="confusion radius (arcsec), where sources in annulus of R_match ~ R_confusion lead to confusion issues (recommend: 60)")
parser.add_argument("--conf_limit",type=float,default=0.90,help="confidence limit for upper limit calculation, defaults to 0.90")
args = parser.parse_args()


# make output directory if necessary
outdir = os.path.dirname(args.prefix)
os.makedirs(outdir,exist_ok=True)
# make log_dir
log_dir = f"{outdir}/log"
os.makedirs(log_dir,exist_ok=True)


# define logger
logger = logging.getLogger("forcedphot")
logger.setLevel(logging.DEBUG)
logname = f"{log_dir}/{os.path.basename(args.prefix)}forcedphot{args.suffix}.log"
os.system(f"rm -rf {logname}")
file_handler = logging.FileHandler(logname)
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def main():

    skytile = tile_local(args.target_ra,args.target_dec)
    if skytile is None: # not in eRO-DE sky
        logger.error(f"Target source (RA={args.target_ra}, DEC={args.target_dec}) not in eROSITA:DE sky!")
        raise

    skip_exist_str = "--skip_exist" if args.skip_exist_srcdet else ""
    srcdet_cmd = [
        "erosrcdet",
        f"{args.science_evt}",
        "--skytile",f"{str(skytile)}",
        "--emin","0.2",
        "--emax","2.3",
        "--target_z","0",
        "--prefix",f"{args.detprefix}",
        "--suffix",f"{args.detsuffix}",
        skip_exist_str,
    ]
    logger.info(" ".join(srcdet_cmd))
    srcdet_log = f"{log_dir}/{os.path.basename(args.prefix)}srcdet{args.suffix}.log"
    with open(srcdet_log,"w") as log_file:
        subprocess.run(srcdet_cmd,stdout=log_file)

    # source detection files (SKYTILE-wise)
    science_img = f"{args.detprefix}sciimg{args.detsuffix}.fits"
    expmap = f"{args.detprefix}expmap{args.detsuffix}.fits"
    detmask = f"{args.detprefix}detmask{args.detsuffix}.fits"
    boxlist = f"{args.detprefix}fake_boxlist{args.detsuffix}.fits"
    bkgmap = f"{args.detprefix}bkgmap{args.detsuffix}.fits"
    cheesemask = f"{args.detprefix}cheesemask{args.detsuffix}.fits"
    psfmap = f"{args.detprefix}psfmap{args.detsuffix}.fits"
    mllist = f"{args.detprefix}fake_mllist{args.detsuffix}.fits"
    srcmap = f"{args.detprefix}srcmap{args.detsuffix}.fits"

    # aperture photometry files (SKYTILE&SRCID-wise)
    apelist = f"{args.prefix}apelist{args.suffix}.fits"
    apelistout = f"{args.prefix}apelistout{args.suffix}.fits"
    apesummary = f"{args.prefix}apesummary{args.suffix}.fits"

    obs_emin = args.emin / (1+args.target_z)
    obs_emax = args.emax / (1+args.target_z)

    # create an aperture file
    logger.info("###############################################################")
    logger.info("#################### Creating aperture list ###################")
    logger.info("###############################################################")
    create_aperture(args.target_ra,args.target_dec,0.75,0.75,apelist)


    # Finally, perform aperture photometry
    # det/non-det in MAIN+SUPP catalog does not matter; 
    # we use stackflag=yes for both cases
    logger.info("###############################################################")
    logger.info("############### Performing aperture photometry ################")
    logger.info("###############################################################")
    """
    apetool apelist="fake_apelist_2.fits" apelistout="fake_ape_out_2.fits" images="events_image_comb.fits" expimages="../output_expmap.fits" bkgimages="test_bkgmap.fits" psfmaps="psf_map.fits" detmasks="../detmask.fits" stackflag=yes emin="200" emax="2300" eefextract=0.75 eindex="1" srcimages="test_bkgmap.fits"
    """
    # check if there are source confusion problems, TODO: add log here
    confusion_log = f"{log_dir}/{os.path.basename(args.prefix)}confusion{args.suffix}.log"
    with open(confusion_log,"w") as log_file:
        sys.stdout = log_file
        sys.stderr = log_file
        n_confusion = look_for_confusion(
            args.target_ra,
            args.target_dec,
            R_match=args.R_match*u.arcsec,
            R_confusion=args.R_confusion*u.arcsec
        ) # find all sources from catalog within 60 arcsec
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
    

    # if there are no contamination sources nearby (15 ~ 60 arcsec), we can safely perform aperture tool without caring about SRCMAP
    if n_confusion == 0:
        os.system(f"rm -rf {apelistout}")                   # clobber=yes
        apetool_cmd = [
            "apetool",
            f"apelist={apelist}",                           # APELIST: Src list for count extraction
            f"apelistout={apelistout}",                     # Output src list with extracted counts
            f"images={science_img}",                        # Input images/event files
            f"psfmaps={psfmap}",                            # Input/Output PSF size maps
            f"expimages={expmap}",                          # Input exposure maps
            f"detmasks={detmask}",                          # Input detection masks
            f"bkgimages={bkgmap}",                          # Input background images
            f"srcimages={bkgmap}",                          # Input source maps (here we take bkgmap as no need to remove confusion sources)
            f"emin={obs_emin*1000}",                        # Minimum energies [PI channels]
            f"emax={obs_emax*1000}",                        # Maximum energies [PI channels]
            "eindex=1",                                     # Lists of indeces to map the input images to the ermldet band indeces
            "eefextract=0.75",                              # EEF for count extraction
            "pthresh=4e-6",                                 # Poisson false detection threshold
            "cutrad=15",                                    # Cut radius for source fitting
            "stackflag=yes",                                # Extract counts at generic srclist positions?  
        ]

    # if there are contamination sources nearby (15 ~ 60 arcsec), we have to care about SRCMAP and remove these contaminations
    elif n_confusion > 0:
        # TODO: need to check an example with n_src>1 if this really work
        os.system(f"rm -rf {apelistout}")                   # clobber=yes
        apetool_cmd = [
            "apetool",
            f"mllist={mllist}",                             # MLLIST from ERMLDET, used in combination with srcimages for source subtraction within RR from SRCMAP
            f"apelist={apelist}",                           # APELIST: Src list for count extraction
            f"apelistout={apelistout}",                     # Output src list with extracted counts
            f"images={science_img}",                        # Input images/event files
            f"psfmaps={psfmap}",                            # Input/Output PSF size maps
            f"expimages={expmap}",                          # Input exposure maps
            f"detmasks={detmask}",                          # Input detection masks
            f"bkgimages={bkgmap}",                          # Input background images
            f"srcimages={srcmap}",                          # (*) Input source maps
            f"emin={obs_emin*1000}",                        # Minimum energies [PI channels]
            f"emax={obs_emax*1000}",                        # Maximum energies [PI channels]
            "eindex=1",                                     # Lists of indeces to map the input images to the ermldet band indeces
            "eefextract=0.75",                              # EEF for count extraction
            "pthresh=4e-6",                                 # Poisson false detection threshold
            "cutrad=15",                                    # Cut radius for source fitting
            "stackflag=yes",                                # Extract counts at generic srclist positions?
        ]

    else:
        raise Exception(f"{n_confusion} matching source(s) is not supported.")
    
    logger.info(" ".join(apetool_cmd))
    apetool_log = f"{log_dir}/{os.path.basename(args.prefix)}apetool{args.suffix}.log"
    with open(apetool_log,"w") as log_file:
        subprocess.run(apetool_cmd,stdout=log_file)


    # read apelistout
    with fits.open(apelistout) as hdu:
        data = hdu[1].data
    ape_cts = data["APE_CTS"][0]
    ape_bkg = data["APE_BKG"][0]
    ape_exp = data["APE_EXP"][0]
    ape_eef = data["APE_EEF"][0]
    ape_radius = data["APE_RADIUS"][0]  # NOTE: the unit here is actually `pixels`, not `arcsec`!
    ape_pois = data["APE_POIS"][0]

    
    ecflike_log = f"{log_dir}/{os.path.basename(args.prefix)}ecflike{args.suffix}.log"
    with open(ecflike_log,"w") as log_file:
        sys.stdout = log_file
        sys.stderr = log_file
        # calculate ecf (credit: T.Liu)
        ecf = get_eroecf(obs_emin,obs_emax,obs_emin,obs_emax)
        # calculate detection likelihood
        detlike = cal_detlike(data["APE_CTS"],data["APE_BKG"])
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
    

    # calculate source count rate, flux
    ME,LO,HI = counts_twoside_poisson(data["APE_CTS"],data["APE_BKG"],0.68)
    CR_ME = ME / data["APE_EEF"] / data["APE_EXP"]
    CR_LO = LO / data["APE_EEF"] / data["APE_EXP"]
    CR_HI = HI / data["APE_EEF"] / data["APE_EXP"]
    FLUX_ME = CR_ME / ecf
    FLUX_LO = CR_LO / ecf
    FLUX_HI = CR_HI / ecf

    # calculate upper limit (credit: A.Ruiz)
    UL = counts_ul_poisson(data["APE_CTS"],data["APE_BKG"],conf_limit=args.conf_limit) # estimate 90% (or user-specified) upper limit in units of counts
    CR_UL =  UL / data["APE_EEF"] / data["APE_EXP"]     # estimate 90% (or user-specified) upper limit on the count rate inlcuding correction for the EEF
    FLUX_UL =  CR_UL / ecf  # estimate 90% (or user-specified) upper limit on the X-ray flux assuming an ECF correction

    # write data
    hdu_lst = fits.HDUList()
    hdu_primary = fits.PrimaryHDU()
    hdu_lst.append(hdu_primary)

    col_detlike = fits.Column(name="like",format="E",array=detlike)

    col_ME = fits.Column(name="ME",format="E",array=ME,unit="cts")
    col_LO = fits.Column(name="LO",format="E",array=LO,unit="cts")
    col_HI = fits.Column(name="HI",format="E",array=HI,unit="cts")
    col_UL = fits.Column(name="UL",format="E",array=UL,unit="cts")

    col_CR_ME = fits.Column(name="CR_ME",format="E",array=CR_ME,unit="cts/s")
    col_CR_LO = fits.Column(name="CR_LO",format="E",array=CR_LO,unit="cts/s")
    col_CR_HI = fits.Column(name="CR_HI",format="E",array=CR_HI,unit="cts/s")
    col_CR_UL = fits.Column(name="CR_UL",format="E",array=CR_UL,unit="cts/s")

    col_FLUX_ME = fits.Column(name="FLUX_ME",format="E",array=FLUX_ME,unit="erg/cm^2/s")
    col_FLUX_LO = fits.Column(name="FLUX_LO",format="E",array=FLUX_LO,unit="erg/cm^2/s")
    col_FLUX_HI = fits.Column(name="FLUX_HI",format="E",array=FLUX_HI,unit="erg/cm^2/s")
    col_FLUX_UL = fits.Column(name="FLUX_UL",format="E",array=FLUX_UL,unit="erg/cm^2/s")

    original_cols = fits.ColDefs(data)
    new_cols = fits.ColDefs([col_detlike,
                            col_ME,col_LO,col_HI,col_UL,
                            col_CR_ME,col_CR_LO,col_CR_HI,col_CR_UL,
                            col_FLUX_ME,col_FLUX_LO,col_FLUX_HI,col_FLUX_UL])
    combined_cols = original_cols + new_cols
    hdu_data = fits.BinTableHDU.from_columns(combined_cols,name="APERTURE")
    hdu_lst.append(hdu_data)

    hdu_lst.writeto(apesummary,overwrite=True)



    # print results
    logger.info("**************** Aperture photometry results ******************")
    logger.info(f"Rest-frame {args.emin} -- {args.emax} keV")
    logger.info(f"Obs-frame {obs_emin} -- {obs_emax} keV")
    logger.info("")
    logger.info(f"Aperture radius: {ape_radius} pixel or {ape_radius*4} arcsec (corresponding to {ape_eef} EEF)")
    logger.info(f"Aperture total counts: {ape_cts}")
    logger.info(f"Aperture bkg counts: {ape_bkg}")
    logger.info(f"")
    logger.info(f"Source detection likelihood: {detlike[0]}")
    logger.info(f"Source count rate (EEF corrected, 68%): {CR_ME[0]} (-{CR_LO[0]},+{CR_HI[0]})")
    logger.info(f"Upper limit for source count rate (EEF corrected, {args.conf_limit*100}%): {CR_UL[0]}")
    logger.info(f"Source flux (EEF corrected, 68%): {FLUX_ME[0]:.3e} (-{FLUX_LO[0]:.3e},+{FLUX_HI[0]:.3e})")
    logger.info(f"Upper limit for source flux (EEF corrected, {args.conf_limit*100}%): {FLUX_UL[0]:.3e}")
    logger.info(f"")
    logger.info(f"Results file saved to {apesummary}")


if __name__ == "__main__":
    main()
