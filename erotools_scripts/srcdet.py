#!/usr/bin/env python3
"""
Perform skytile-wise source detection.

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
import time
from erotools.erocat import fake_srclist,look_for_confusion,create_aperture
from erotools.eroecf import get_eroecf
from erotools.detlike import cal_detlike,counts_twoside_poisson,counts_ul_poisson


class HelpfulParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)

parser = HelpfulParser(description=__doc__,
    epilog="""Shi-Jiang Chen, Johannes Buchner and Teng Liu (C) 2025 <JohnnyCsj666@gmail.com>""",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("science_evt",type=str,help="the input image/events file")
parser.add_argument("--skytile",type=str,default=None,help="6-digit skytile; if not provided, will read from `science_evt`")
parser.add_argument("--emin",type=float,default=0.2,help="rest-frame minimum energy [keV]")
parser.add_argument("--emax",type=float,default=2.3,help="rest-frame maximum energy [keV]")
parser.add_argument("--target_z",type=float,default=0.,help="target redshift (real photometry extracted from emin/z, emax/z)")
parser.add_argument("--prefix",type=str,default="./out/sourcedet_",help="prefix for all products, defaults to './out/sourcedet_', will create a directory if necessary")
parser.add_argument("--suffix",type=str,default="",help="suffix for all products")
parser.add_argument("--skip_exist",action="store_true",help="skip source detection if the files already exist, which could save a lot of time for batch processing")
args = parser.parse_args()


# parse skytile
if args.skytile is None:
    skytile = str(os.path.basename(args.science_evt).split("_")[1])     # e.g., eb01_182057_020_EventList_c010.fits.gz
else:
    skytile = args.skytile


# make output directory if necessary
outdir = os.path.dirname(args.prefix)
os.makedirs(outdir,exist_ok=True)
# make log_dir
log_dir = f"{outdir}/log"
os.makedirs(log_dir,exist_ok=True)


# define logger
logger = logging.getLogger("srcdet")
logger.setLevel(logging.DEBUG)
logname = f"{log_dir}/{os.path.basename(args.prefix)}srcdet{args.suffix}.log"
os.system(f"rm -rf {logname}")
file_handler = logging.FileHandler(logname)
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def main():
    t0 = time.time()

    science_img = f"{args.prefix}sciimg{args.suffix}.fits"
    expmap = f"{args.prefix}expmap{args.suffix}.fits"
    detmask = f"{args.prefix}detmask{args.suffix}.fits"
    boxlist = f"{args.prefix}fake_boxlist{args.suffix}.fits"
    bkgmap = f"{args.prefix}bkgmap{args.suffix}.fits"
    cheesemask = f"{args.prefix}cheesemask{args.suffix}.fits"
    psfmap = f"{args.prefix}psfmap{args.suffix}.fits"
    mllist = f"{args.prefix}fake_mllist{args.suffix}.fits"
    srcmap = f"{args.prefix}srcmap{args.suffix}.fits"

    obs_emin = args.emin / (1+args.target_z)
    obs_emax = args.emax / (1+args.target_z)


    # select events within obs-frame emin -- emax
    if (not os.path.exists(science_img)) or (not args.skip_exist):
        logger.info("###############################################################")
        logger.info("###################### Filtering events #######################")
        logger.info("###############################################################")
        os.system(f"rm -rf {science_img}")           # clobber=yes
        evtool_cmd = [
            "evtool",
            f"eventfiles={args.science_evt}",
            f"outfile={science_img}",
            "image=yes",
            f"emin={obs_emin}",
            f"emax={obs_emax}",
        ]
        logger.info(" ".join(evtool_cmd))
        evtool_log = f"{log_dir}/{os.path.basename(args.prefix)}evtool{args.suffix}.log"
        with open(evtool_log,"w") as log_file:
            subprocess.run(evtool_cmd,stdout=log_file)


    # Catalog preparation:
    # first generate a fake boxlist for all detected sources in the target skytile
    if (not os.path.exists(boxlist)) or (not args.skip_exist):
        logger.info("###############################################################")
        logger.info("################### Preparing fake catalog ####################")
        logger.info("###############################################################")
        # outname = "%s/%s"%(args.outdir,boxlist)
        os.system(f"rm -rf {boxlist}")
        fksrc_log = f"{log_dir}/{os.path.basename(args.prefix)}fksrc{args.suffix}.log"
        with open(fksrc_log,"w") as log_file:
            sys.stdout = log_file
            sys.stderr = log_file
            fake_srclist(
                skytile,
                outname=boxlist,
                science_img=science_img,
                style="box",
            )
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__


    # generate EXPMAP
    if (not os.path.exists(expmap)) or (not args.skip_exist):
        logger.info("###############################################################")
        logger.info("################### Generating exposure map ###################")
        logger.info("###############################################################")
        os.system(f"rm -rf {expmap}")                           # clobber=yes
        expmap_cmd = [
            "expmap",
            f"inputdatasets={science_img}",                     # Event files
            f"templateimage={science_img}",                     # Template image
            "withmergedmaps=yes",                               # Create all-telescope merged exposure maps?
            f"mergedmaps={expmap}",                             # Name of all-telescope merged exposure maps
            f"emin={obs_emin}",                                 # Minimum energy [keV]
            f"emax={obs_emax}",                                 # Maximum energy [keV]
            "withvignetting=yes",                               # With vignetting?
            "withweights=yes",                                  # Apply camera area weight factors?
            "plindex=-1.7"                                      # Power-law index of energy for vignet weighting.
        ]
        logger.info(" ".join(expmap_cmd))
        expmap_log = f"{log_dir}/{os.path.basename(args.prefix)}expmap{args.suffix}.log"
        with open(expmap_log,"w") as log_file:
            subprocess.run(expmap_cmd,stdout=log_file)


    # generate DETMASK
    if (not os.path.exists(detmask)) or (not args.skip_exist):
        logger.info("###############################################################")
        logger.info("################## Generating detection map ###################")
        logger.info("###############################################################")
        os.system(f"rm -rf {detmask}")                          # clobber=yes
        ermask_cmd = [
            "ermask",
            f"expimage={expmap}",                               # Exposure map
            f"detmask={detmask}",                               # Output detection mask
            "threshold1=0.1",                                   # rel. exposure threshold
            "threshold2=1.0",                                   # exposure gradient threshold
        ]
        logger.info(" ".join(ermask_cmd))
        ermask_log = f"{log_dir}/{os.path.basename(args.prefix)}ermask{args.suffix}.log"
        with open(ermask_log,"w") as log_file:
            subprocess.run(ermask_cmd,stdout=log_file)


    # generate BKGMAP
    # this assumes all but DET_LIKE_0>5 events (sources) arising from BKG
    if (not os.path.exists(bkgmap)) or (not os.path.exists(cheesemask)) or (not args.skip_exist):
        logger.info("###############################################################")
        logger.info("################## Generating background map ##################")
        logger.info("###############################################################")
        erbackmap_cmd = [
            "erbackmap",
            f"image={science_img}",                             # Input image
            f"boxlist={boxlist}",                               # (Fake) boxdetect list
            f"bkgimage={bkgmap}",                               # Output background map
            "expima_flag=yes",                                  # Use exposure map?
            "expima2_flag=no",                                  # Use unvignetted exposure map? Not yet supported
            "detmask_flag=yes",                                 # Use detection mask?
            "cheesemask_flag=yes",                              # Write cheesed mask?
            f"expimage={expmap}",                               # Exposure map
            "expimage2=''",                                     # Exposure map (unvignetted), not yet supported
            f"detmask={detmask}",                               # Detection mask
            f"cheesemask={cheesemask}",                         # Output cheesed mask
            "idband=1",                                         # Energy band index in boxlist
            f"emin={obs_emin*1000}",                            # Minimum energy [PI channels]
            f"emax={obs_emax*1000}",                            # Maximum energy [PI channels]
            "scut=0.0001",                                      # Surface brightness limit [cts/pix]
            "mlmin=5",                                          # Minimum boxlist likelihood
            "maxcut=0.5",                                       # Maximum masked out area fraction
            "fitmethod=smooth",                                 # Fit method (spline / smooth)
            "snr=40",                                           # Signal to noise ratio for adaptive smoothing
            "smoothval=15",                                     # Smallest kernel size
            "nfitrun=3",                                        # Number of iterations to reject bins
            "excesssigma=1000",                                 # Chi square limit for rejected bins
            "nsplinenodes=36",                                  # Number of spline nodes / spatial bins per dimension
            "degree=2",                                         # Degree of the spline (1 <= k <= 5)
            "smoothflag=yes",                                   # Use smoothing spline?    
            "clobber=yes"                                       # Overwrite existing output file(s)
        ]
        logger.info(" ".join(erbackmap_cmd))
        erbackmap_log = f"{log_dir}/{os.path.basename(args.prefix)}erbackmap{args.suffix}.log"
        with open(erbackmap_log,"w") as log_file:
            subprocess.run(erbackmap_cmd,stdout=log_file)


    # generate PSFMAP
    if (not os.path.exists(psfmap)) or (not args.skip_exist):
        logger.info("###############################################################")
        logger.info("###################### Generating PSF map #####################")
        logger.info("###############################################################")
        os.system(f"rm -rf {psfmap}")                           # clobber=yes
        psfmap_cmd = [
            "apetool",
            f"images={science_img}",                            # Input images/event files
            f"psfmaps={psfmap}",                                # Output PSF size maps
            "psfmapflag=yes"                                    # Produce PSF map?
        ]
        logger.info(" ".join(psfmap_cmd))
        psfmap_log = f"{log_dir}/{os.path.basename(args.prefix)}psfmap{args.suffix}.log"
        with open(psfmap_log,"w") as log_file:
            subprocess.run(psfmap_cmd,stdout=log_file)


    # generate MLLIST
    if (not os.path.exists(mllist)) or (not os.path.exists(srcmap)) or (not args.skip_exist):
        logger.info("###############################################################")
        logger.info("###################### Generating MLlist ######################")
        logger.info("###############################################################")
        os.system(f"rm -rf {mllist} {srcmap}")
        ecf_log = f"{log_dir}/{os.path.basename(args.prefix)}ecf{args.suffix}.log"
        with open(ecf_log,"w") as log_file:
            sys.stdout = log_file
            sys.stderr = log_file
            ecf = get_eroecf(obs_emin,obs_emax,obs_emin,obs_emax)
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
        ermldet_cmd = [
            "ermldet",
            f"boxlist={boxlist}",                               # boxdetect input list (here taken from existing MAIN+SUPP catalog)
            f"images={science_img}",                            # Input images
            f"mllist={mllist}",                                 # ermldet output list
            f"expimages={expmap}",                              # Exposure maps
            f"detmasks={detmask}",                              # Detection masks
            f"bkgimages={bkgmap}",                              # Background maps
            f"srcimages={srcmap}",                              # Output source maps
            f"emin={obs_emin*1000}",                            # Minimum energies [PI channels]
            f"emax={obs_emax*1000}",                            # Maximum energies [PI channels]
            f"ecf={ecf}",                                       # energy conversion factor
            "likemin=5",                                        # Minimum ML likelihood (here we take 5 as this is the minimum like for SUPP catalog)
            "extentmodel=beta",                                 # Extent model (gaussian | beta)
            "nmaxfit=3",                                        # Max. number of sources for simultaneous fit
            "nmulsou=2",                                        # Maximum new sources for source splitting
            "thres_flag=no",                                    # Use threshold for source splitting
            "thres_col=like",                                   # Column to apply threshold (LIKE|RATE|SCTS)
            "multrad=15",                                       # Search radius for multiple source fitting
            "cutrad=15",                                        # Cut radius for source fitting
            "srcima_flag=yes",                                  # Write source model maps?
            "expima_flag=yes",                                  # Use exposure maps?
            "detmask_flag=yes",                                 # Use detection mask?
            "shapelet_flag=yes",                                # Use shapelet PSF?
        ]
        logger.info(" ".join(ermldet_cmd))
        ermldet_log = f"{log_dir}/{os.path.basename(args.prefix)}ermldet{args.suffix}.log"
        with open(ermldet_log,"w") as log_file:
            subprocess.run(ermldet_cmd,stdout=log_file)


    t1 = time.time()
    logger.info(f"Total time use: {t1-t0} s")



if __name__ == "__main__":
    main()
