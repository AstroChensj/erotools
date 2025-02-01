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
from ero_tools.erocat import fake_boxlist,look_for_confusion,create_aperture
from ero_tools.eroecf import get_eroecf
from ero_tools.detlike import cal_detlike,counts_twoside_poisson,counts_ul_poisson


class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Shi-Jiang Chen, Johannes Buchner and Teng Liu (C) 2024 <JohnnyCsj666@gmail.com>""",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('science_img',type=str,help='the input image/events file')
parser.add_argument('target_ra',type=float,help='target source RA (degree, icrs)')
parser.add_argument('target_dec',type=float,help='target source DEC (degree, icrs)')
parser.add_argument('--emin',type=float,default=0.2,help='rest-frame minimum energy [keV]')
parser.add_argument('--emax',type=float,default=2.3,help='rest-frame maximum energy [keV]')
parser.add_argument('--target_z',type=float,default=0,help='target redshift (real photometry extracted from emin/z, emax/z)')
parser.add_argument('--outdir',type=str,default='outdir',help='output directory name')
parser.add_argument('--outname',type=str,default='result.fits',help='name of fits file storing the aperture photometry results, saved under outdir')
parser.add_argument('--prefix',type=str,default='',help='prefix for all products')
parser.add_argument('--suffix',type=str,default='',help='suffix for all products')
parser.add_argument('--R_match',type=float,default=15,help='matching radius (arcsec) between target and dr1 catalog (recommend: 15)')
parser.add_argument('--R_confusion',type=float,default=60,help='confusion radius (arcsec), where sources in annulus of R_match ~ R_confusion lead to confusion issues (recommend: 60)')
# science_img = sys.argv[1]
# emin = float(sys.argv[2])                       # in units of eV, e.g., 200
# emax = float(sys.argv[3])                       # in units of eV, e.g., 2300
# target_ra = float(sys.argv[4])                  # target RA, icrs
# target_dec = float(sys.argv[5])                 # target DEC, icrs
# outdir = sys.argv[6]
# prefix = sys.argv[7]
# suffix = sys.argv[8]

args = parser.parse_args()


def main():

    events = '%sevents%s.fits'%(args.prefix,args.suffix)
    expmap = '%sexpmap%s.fits'%(args.prefix,args.suffix)
    detmask = '%sdetmask%s.fits'%(args.prefix,args.suffix)
    boxlist = '%sfake_boxlist%s.fits'%(args.prefix,args.suffix)
    bkgmap = '%sbkgmap%s.fits'%(args.prefix,args.suffix)
    cheesemask = '%scheesemask%s.fits'%(args.prefix,args.suffix)
    psfmap = '%spsfmap%s.fits'%(args.prefix,args.suffix)
    apelist = '%sfake_apelist%s.fits'%(args.prefix,args.suffix)
    apelistout = '%sapelistout%s.fits'%(args.prefix,args.suffix)
    mllist = '%sfake_mllist%s.fits'%(args.prefix,args.suffix)
    srcmap = '%ssrcmap%s.fits'%(args.prefix,args.suffix)

    obs_emin = args.emin / (1+args.target_z)
    obs_emax = args.emax / (1+args.target_z)

    os.system('mkdir -p %s'%(args.outdir))


    # select events within obs-frame emin -- emax
    print("###############################################################")
    print("###################### Filtering events #######################")
    print("###############################################################")
    os.system('rm -rf %s/%s'%(args.outdir,events))          # clobber=yes
    subprocess.run([
        'evtool',
        'eventfiles=%s'%(args.science_img),
        'outfile=%s/%s'%(args.outdir,events),
        'image=yes',
        'emin=%f'%(obs_emin),
        'emax=%f'%(obs_emax)
    ])


    # Catalog preparation:
    # first generate a fake boxlist for all detected sources within target's skytile
    print("###############################################################")
    print("################### Preparing fake catalog ####################")
    print("###############################################################")
    outname = '%s/%s'%(args.outdir,boxlist)
    fake_boxlist(
         args.target_ra,
         args.target_dec,
         outname,
         '%s/%s'%(args.outdir,events)
    )


    # generate EXPMAP
    print("###############################################################")
    print("################### Generating exposure map ###################")
    print("###############################################################")
    os.system('rm -rf %s/%s'%(args.outdir,expmap))          # clobber=yes
    subprocess.run([
        'expmap',
        'inputdatasets=%s/%s'%(args.outdir,events),         # Event files
        'templateimage=%s/%s'%(args.outdir,events),         # Template image
        'withmergedmaps=yes',                               # Create all-telescope merged exposure maps?
        'mergedmaps=%s/%s'%(args.outdir,expmap),            # Name of all-telescope merged exposure maps
        'emin=%f'%(obs_emin),                               # Minimum energy [keV]
        'emax=%f'%(obs_emax),                               # Maximum energy [keV]
        'withvignetting=yes',                               # With vignetting?
        'withweights=yes',                                  # Apply camera area weight factors?
        'plindex=-1.7'                                      # Power-law index of energy for vignet weighting.
    ])


    # generate DETMASK
    print("###############################################################")
    print("################## Generating detection map ###################")
    print("###############################################################")
    os.system('rm -rf %s/%s'%(args.outdir,detmask))         # clobber=yes
    subprocess.run([
        'ermask',
        'expimage=%s/%s'%(args.outdir,expmap),              # Exposure map
        'detmask=%s/%s'%(args.outdir,detmask),              # Output detection mask
        'threshold1=0.1',                                   # rel. exposure threshold
        'threshold2=1.0',                                   # exposure gradient threshold
    ])


    # generate BKGMAP
    print("###############################################################")
    print("################## Generating background map ##################")
    print("###############################################################")
    subprocess.run([
        'erbackmap',
        'image=%s/%s'%(args.outdir,events),                 # Input image
        'boxlist=%s/%s'%(args.outdir,boxlist),              # (Fake) boxdetect list
        'bkgimage=%s/%s'%(args.outdir,bkgmap),              # Output background map
        'expima_flag=yes',                                  # Use exposure map?
        'expima2_flag=no',                                  # Use unvignetted exposure map? Not yet supported
        'detmask_flag=yes',                                 # Use detection mask?
        'cheesemask_flag=yes',                              # Write cheesed mask?
        'expimage=%s/%s'%(args.outdir,expmap),              # Exposure map
        'expimage2=""',                                     # Exposure map (unvignetted), not yet supported
        'detmask=%s/%s'%(args.outdir,detmask),              # Detection mask
        'cheesemask=%s/%s'%(args.outdir,cheesemask),        # Output cheesed mask
        'idband=1',                                         # Energy band index in boxlist
        'emin=%d'%(obs_emin*1000),                          # Minimum energy [PI channels]
        'emax=%d'%(obs_emax*1000),                          # Maximum energy [PI channels]
        'scut=0.0001',                                      # Surface brightness limit [cts/pix]
        'mlmin=5',                                          # Minimum boxlist likelihood
        'maxcut=0.5',                                       # Maximum masked out area fraction
        'fitmethod=smooth',                                 # Fit method (spline / smooth)
        'snr=40',                                           # Signal to noise ratio for adaptive smoothing
        'smoothval=15',                                     # Smallest kernel size
        'nfitrun=3',                                        # Number of iterations to reject bins
        'excesssigma=1000',                                 # Chi square limit for rejected bins
        'nsplinenodes=36',                                  # Number of spline nodes / spatial bins per dimension
        'degree=2',                                         # Degree of the spline (1 <= k <= 5)
        'smoothflag=yes',                                   # Use smoothing spline?    
        'clobber=yes'                                       # Overwrite existing output file(s)
    ])


    # generate PSFMAP
    print("###############################################################")
    print("###################### Generating PSF map #####################")
    print("###############################################################")
    os.system('rm -rf %s/%s'%(args.outdir,psfmap))          # clobber=yes
    subprocess.run([
        'apetool',
        'images=%s/%s'%(args.outdir,events),                # Input images/event files
        'psfmaps=%s/%s'%(args.outdir,psfmap),               # Output PSF size maps
        'psfmapflag=yes'                                    # Produce PSF map?
    ])


    # create an aperture file
    print("###############################################################")
    print("#################### Creating aperture list ###################")
    print("###############################################################")
    outname = '%s/%s'%(args.outdir,apelist)
    create_aperture(args.target_ra,args.target_dec,0.75,0.75,outname)


    # Finally, perform aperture photometry
    # det/non-det in MAIN+SUPP catalog does not matter; 
    # we use stackflag=yes for both cases
    print("###############################################################")
    print("############### Performing aperture photometry ################")
    print("###############################################################")
    '''
    apetool apelist="fake_apelist_2.fits" apelistout="fake_ape_out_2.fits" images="events_image_comb.fits" expimages="../output_expmap.fits" bkgimages="test_bkgmap.fits" psfmaps="psf_map.fits" detmasks="../detmask.fits" stackflag=yes emin="200" emax="2300" eefextract=0.75 eindex="1" srcimages="test_bkgmap.fits"
    '''
    # check if there are source confusion problems
    n_confusion = look_for_confusion(
         args.target_ra,
         args.target_dec,
         R_match=args.R_match*u.arcsec,
         R_confusion=args.R_confusion*u.arcsec
    ) # find all sources from catalog within 60 arcsec

    # if there are no contamination sources nearby (15 ~ 60 arcsec), we can safely perform aperture tool without caring about SRCMAP
    if n_confusion == 0:
        os.system('rm -rf %s/%s'%(args.outdir,apelistout))  # clobber=yes
        subprocess.run([
            'apetool',
            'apelist=%s/%s'%(args.outdir,apelist),          # APELIST: Src list for count extraction
            'apelistout=%s/%s'%(args.outdir,apelistout),    # Output src list with extracted counts
            'images=%s/%s'%(args.outdir,events),            # Input images/event files
            'psfmaps=%s/%s'%(args.outdir,psfmap),           # Input/Output PSF size maps
            'expimages=%s/%s'%(args.outdir,expmap),         # Input exposure maps
            'detmasks=%s/%s'%(args.outdir,detmask),         # Input detection masks
            'bkgimages=%s/%s'%(args.outdir,bkgmap),         # Input background images
            'srcimages=%s/%s'%(args.outdir,bkgmap),         # Input source maps (here we take bkgmap as no need to remove confusion sources)
            'emin=%d'%(obs_emin*1000),                      # Minimum energies [PI channels]
            'emax=%d'%(obs_emax*1000),                      # Maximum energies [PI channels]
            'eindex=1',                                     # Lists of indeces to map the input images to the ermldet band indeces
            'eefextract=0.75',                              # EEF for count extraction
            'pthresh=4e-6',                                 # Poisson false detection threshold
            'cutrad=15',                                    # Cut radius for source fitting
            'stackflag=yes',                                # Extract counts at generic srclist positions?  
        ])

    # if there are contamination sources nearby (15 ~ 60 arcsec), we have to care about SRCMAP and remove these contaminations
    elif n_confusion > 0:
        # first generate SRCMAP based on boxlist (note, this is not standard routine)
        '''
        ermldet mllist="mllist.fits" boxlist="../src_catalog/test_src.fits" images="events_image_comb.fits" expimages="../output_expmap.fits" detmasks="../detmask.fits" bkgimages="test_bkgmap.fits" extentmodel=beta srcimages="sourceimage.fits" emin=200 emax=2300 likemin=5
        '''
        os.system('rm -rf %s/%s %s/%s'%(args.outdir,mllist,args.outdir,srcmap))
        subprocess.run([
            'ermldet',
            'boxlist=%s/%s'%(args.outdir,boxlist),          # boxdetect input list (here taken from existing MAIN+SUPP catalog)
            'images=%s/%s'%(args.outdir,events),            # Input images
            'mllist=%s/%s'%(args.outdir,mllist),            # ermldet output list
            'expimages=%s/%s'%(args.outdir,expmap),         # Exposure maps
            'detmasks=%s/%s'%(args.outdir,detmask),         # Detection masks
            'bkgimages=%s/%s'%(args.outdir,bkgmap),         # Background maps
            'srcimages=%s/%s'%(args.outdir,srcmap),         # Output source maps
            'emin=%d'%(obs_emin*1000),                      # Minimum energies [PI channels]
            'emax=%d'%(obs_emax*1000),                      # Maximum energies [PI channels]
            'ecf=1',                                        # TODO: merge get_eroecf.py and replace ecf here
            'likemin=5',                                    # Minimum ML likelihood (here we take 5 as this is the minimum like for SUPP catalog)
            'extentmodel=beta',                             # Extent model (gaussian | beta)
            'nmaxfit=3',                                    # Max. number of sources for simultaneous fit
            'nmulsou=2',                                    # Maximum new sources for source splitting
            'thres_flag=no',                                # Use threshold for source splitting
            'thres_col=like',                               # Column to apply threshold (LIKE|RATE|SCTS)
            'multrad=15',                                   # Search radius for multiple source fitting
            'cutrad=15',                                    # Cut radius for source fitting
            'srcima_flag=yes',                              # Write source model maps?
            'expima_flag=yes',                              # Use exposure maps?
            'detmask_flag=yes',                             # Use detection mask?
            'shapelet_flag=yes',                            # Use shapelet PSF?
        ])
        # then run aperture photometry, with SRCMAP (contamination now accounted for in the background)
        # TODO: need to check an example with n_src>1 if this really work
        os.system('rm -rf %s/%s'%(args.outdir,apelistout))  # clobber=yes
        subprocess.run([
            'apetool',
            'mllist=%s/%s'%(args.outdir,mllist),            # MLLIST from ERMLDET, used in combination with srcimages for source subtraction within RR from SRCMAP
            'apelist=%s/%s'%(args.outdir,apelist),          # APELIST: Src list for count extraction
            'apelistout=%s/%s'%(args.outdir,apelistout),    # Output src list with extracted counts
            'images=%s/%s'%(args.outdir,events),            # Input images/event files
            'psfmaps=%s/%s'%(args.outdir,psfmap),           # Input/Output PSF size maps
            'expimages=%s/%s'%(args.outdir,expmap),         # Input exposure maps
            'detmasks=%s/%s'%(args.outdir,detmask),         # Input detection masks
            'bkgimages=%s/%s'%(args.outdir,bkgmap),         # Input background images
            'srcimages=%s/%s'%(args.outdir,srcmap),         # (*) Input source maps
            'emin=%d'%(obs_emin*1000),                      # Minimum energies [PI channels]
            'emax=%d'%(obs_emax*1000),                      # Maximum energies [PI channels]
            'eindex=1',                                     # Lists of indeces to map the input images to the ermldet band indeces
            'eefextract=0.75',                              # EEF for count extraction
            'pthresh=4e-6',                                 # Poisson false detection threshold
            'cutrad=15',                                    # Cut radius for source fitting
            'stackflag=yes',                                # Extract counts at generic srclist positions?
        ])

    else:
        raise Exception('%d matching source(s) is not supported.'%(n_confusion))


    # read apelistout
    with fits.open('%s/%s'%(args.outdir,apelistout)) as hdu:
        data = hdu[1].data
    ape_cts = data['APE_CTS'][0]
    ape_bkg = data['APE_BKG'][0]
    ape_exp = data['APE_EXP'][0]
    ape_eef = data['APE_EEF'][0]
    ape_radius = data['APE_RADIUS'][0]
    ape_pois = data['APE_POIS'][0]

    # calculate ecf (credit: T.Liu)
    ecf = get_eroecf(obs_emin,obs_emax,obs_emin,obs_emax)

    # calculate detection likelihood
    detlike = cal_detlike(data['APE_CTS'],data['APE_BKG'])

    # calculate source count rate, flux
    ME,LO,HI = counts_twoside_poisson(data['APE_CTS'],data['APE_BKG'],0.68)
    CR_ME = ME / data['APE_EEF'] / data['APE_EXP']
    CR_LO = LO / data['APE_EEF'] / data['APE_EXP']
    CR_HI = HI / data['APE_EEF'] / data['APE_EXP']
    FLUX_ME = CR_ME / ecf
    FLUX_LO = CR_LO / ecf
    FLUX_HI = CR_HI / ecf

    # calculate upper limit (credit: A.Ruiz)
    UL = counts_ul_poisson(data['APE_CTS'],data['APE_BKG'],conf_limit=0.90) # estimate 90% upper limit in units of counts
    CR_UL =  UL / data['APE_EEF'] / data['APE_EXP']     # estimate 90% upper limit on the count rate inlcuding correction for the EEF
    FLUX_UL =  CR_UL / ecf  # estimate 90% upper limit on the X-ray flux assuming an ECF correction

    # write data
    hdu_lst = fits.HDUList()
    hdu_primary = fits.PrimaryHDU()
    hdu_lst.append(hdu_primary)

    col_detlike = fits.Column(name='like',format='E',array=detlike)

    col_ME = fits.Column(name='ME',format='E',array=ME,unit='cts')
    col_LO = fits.Column(name='LO',format='E',array=LO,unit='cts')
    col_HI = fits.Column(name='HI',format='E',array=HI,unit='cts')
    col_UL = fits.Column(name='UL',format='E',array=UL,unit='cts')

    col_CR_ME = fits.Column(name='CR_ME',format='E',array=CR_ME,unit='cts/s')
    col_CR_LO = fits.Column(name='CR_LO',format='E',array=CR_LO,unit='cts/s')
    col_CR_HI = fits.Column(name='CR_HI',format='E',array=CR_HI,unit='cts/s')
    col_CR_UL = fits.Column(name='CR_UL',format='E',array=CR_UL,unit='cts/s')

    col_FLUX_ME = fits.Column(name='FLUX_ME',format='E',array=FLUX_ME,unit='erg/cm^2/s')
    col_FLUX_LO = fits.Column(name='FLUX_LO',format='E',array=FLUX_LO,unit='erg/cm^2/s')
    col_FLUX_HI = fits.Column(name='FLUX_HI',format='E',array=FLUX_HI,unit='erg/cm^2/s')
    col_FLUX_UL = fits.Column(name='FLUX_UL',format='E',array=FLUX_UL,unit='erg/cm^2/s')

    original_cols = fits.ColDefs(data)
    new_cols = fits.ColDefs([col_detlike,
                            col_ME,col_LO,col_HI,col_UL,
                            col_CR_ME,col_CR_LO,col_CR_HI,col_CR_UL,
                            col_FLUX_ME,col_FLUX_LO,col_FLUX_HI,col_FLUX_UL])
    combined_cols = original_cols + new_cols
    hdu_data = fits.BinTableHDU.from_columns(combined_cols,name='APERTURE')
    hdu_lst.append(hdu_data)

    hdu_lst.writeto('%s/%s'%(args.outdir,args.outname),overwrite=True)



    # print results
    print("**************** Aperture photometry results ******************")
    print("Rest-frame %f -- %f keV"%(args.emin,args.emax))
    print("Obs-frame %f -- %f keV"%(obs_emin,obs_emax))
    print('')
    print("Aperture radius: %f pixel or %f arcsec (corresponding to %f EEF)"%(ape_radius,ape_radius*4,ape_eef))
    print("Aperture total counts: %f"%(ape_cts))
    print("Aperture bkg counts: %f"%(ape_bkg))
    print('')
    print("Source detection likelihood: %f"%(detlike[0]))
    print("Source count rate (EEF corrected, 68%%): %f (-%f,+%f)"%(CR_ME[0],CR_LO[0],CR_HI[0]))
    print("Upper limit for source count rate (EEF corrected, 90%%): %f"%(CR_UL[0]))
    print("Source flux (EEF corrected, 68%%): %.3e (-%.3e,+%.3e)"%(FLUX_ME[0],FLUX_LO[0],FLUX_HI[0]))
    print("Upper limit for source flux (EEF corrected, 90%%): %.3e"%(FLUX_UL[0]))
    print('')
    print("Results file saved to %s/%s"%(args.outdir,args.outname))


if __name__ == '__main__':
    main()