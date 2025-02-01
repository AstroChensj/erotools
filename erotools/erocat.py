import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord,search_around_sky
import astropy.units as u
from astropy.wcs import WCS
from astropy_healpix import HEALPix
import requests

#############################################
########## CATALOG & SCIENCE IMAGE ##########
#############################################
main_cat = '/data/chensj/eROSITA/catalogs/eRASS1_Main.v1.1.fits'
supp_cat = '/data/chensj/eROSITA/catalogs/eRASS1_Supp.v1.1.fits'
main_cat_5 = '/data/chensj/eROSITA/catalogs/all_s4_s5_SourceCat1B_c030_240905_poscorr_mpe_photom.fits'


#############################################
############## MAIN FUNCTIONS ###############
#############################################
def fake_boxlist(target_ra,target_dec,outname,science_img):
    """
    Look for a target in eRASS1 Main+Supp catalogs, and generte a `boxlist` for all detected sources in the skytile where it belongs to.

    Parameters
    ----------
    target_ra : float
    target_dec : float
    outname : str
    science_img : str

    Returns
    -------
    None
    """
    with fits.open(main_cat) as hdu:
        main_data = hdu[1].data
    with fits.open(supp_cat) as hdu:
        supp_data = hdu[1].data
    data = np.append(main_data,supp_data)

    skytile = find_erode_skytile(target_ra,target_dec)

    if skytile is not None:
        skytile = int(skytile)  # to match catalog standard
        data_inskytile = data[data['SKYTILE']==skytile]
        print("********************** Creating boxlist ***********************")
        # get wcs
        with fits.open(science_img) as hdu:
            img_head = hdu[0].header  # Use the primary HDU or the appropriate extension
            evt_head = hdu[1].header
            wcs = WCS(img_head)
        # get hdu[1].header
        excluded_keys = [
            'XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','TFIELDS','TTYPE','TFORM','TUNIT','EXTNAME',
            'TSCAL','TZERO' # this two header keyword will change column values! 
        ]
        science_head = fits.Header({
            key: value for key,value in evt_head.items()
            if not any(key.startswith(prefix) for prefix in excluded_keys)
        })
        # id
        id_src_inskytile = data_inskytile['ID_SRC']
        id_band_inskytile = np.ones(len(id_src_inskytile))
        id_inst_inskytile = np.zeros(len(id_src_inskytile))
        dist_nn_inskytile = np.zeros(len(id_src_inskytile)) # set to 0 for the moment
        # coordinates
        ra_inskytile = data_inskytile['RA']
        raerr_inskytile = (data_inskytile['RA_LOWERR'] + data_inskytile['RA_UPERR']) / 2 / 3600     # LOWERR, UPERR in units of arcsec
        dec_inskytile = data_inskytile['DEC']
        decerr_inskytile = (data_inskytile['DEC_LOWERR'] + data_inskytile['DEC_UPERR']) / 2 / 3600
        radec_err_inskytile = data_inskytile['RADEC_ERR']
        x_ima_inskytile,y_ima_inskytile = wcs.all_world2pix(ra_inskytile,dec_inskytile,1)  # image pixel starts with 1
        x1,y1 = wcs.all_world2pix(ra_inskytile+raerr_inskytile,dec_inskytile,1)
        x2,y2 = wcs.all_world2pix(ra_inskytile,dec_inskytile+decerr_inskytile,1)
        pxpra_dra = x1 - x_ima_inskytile    # partial(x)/partial(ra) * d ra
        pxpdec_ddec = x2 - x_ima_inskytile  # partial(x)/partial(dec) * d dec
        x_ima_err_inskytile = np.sqrt(pxpra_dra**2+pxpdec_ddec**2)
        pypra_dra = y1 - y_ima_inskytile    # partial(y)/partial(ra) * d ra
        pypdec_ddec = y2 - y_ima_inskytile  # partial(y)/partial(dec) * d dec
        y_ima_err_inskytile = np.sqrt(pypra_dra**2+pypdec_ddec**2)
        lII_inskytile = data_inskytile['LII']
        bII_inskytile = data_inskytile['BII']
        # counts, rate, flux
        like_inskytile = data_inskytile['DET_LIKE_0']
        scts_inskytile = data_inskytile['ML_CTS_1']
        scts_err_inskytile = data_inskytile['ML_CTS_ERR_1']
        box_cts_inskytile = scts_inskytile
        bg_map_inskytile = data_inskytile['ML_BKG_1']
        bg_raw_inskytile = data_inskytile['ML_BKG_1']
        flux_inskytile = data_inskytile['ML_FLUX_1']
        flux_err_inskytile = data_inskytile['ML_FLUX_ERR_1']
        rate_inskytile = data_inskytile['ML_RATE_1']
        rate_err_inskytile = data_inskytile['ML_RATE_ERR_1']
        exp_map_inskytile = data_inskytile['ML_EXP_1']
        box_size_inskytile = 4*np.ones(len(id_src_inskytile))   # assumes boxsize of 4 as recommended
        eef_inskytile = data_inskytile['ML_EEF_1']
        # make fits file
        arrays = [
            id_src_inskytile,id_inst_inskytile,id_band_inskytile,
            scts_inskytile,scts_err_inskytile,box_cts_inskytile,
            x_ima_inskytile,x_ima_err_inskytile,y_ima_inskytile,y_ima_err_inskytile,
            like_inskytile,bg_map_inskytile,bg_raw_inskytile,exp_map_inskytile,
            flux_inskytile,flux_err_inskytile,rate_inskytile,rate_err_inskytile,
            ra_inskytile,dec_inskytile,radec_err_inskytile,lII_inskytile,bII_inskytile,
            box_size_inskytile,eef_inskytile,dist_nn_inskytile
        ]
        colnames = [
            'id_src','id_inst','id_band',
            'scts','scts_err','box_cts',
            'x_ima','x_ima_err','y_ima','y_ima_err',
            'like','bg_map','bg_raw','exp_map',
            'flux','flux_err','rate','rate_err',
            'ra','dec','radec_err','lII','bII',
            'box_size','eef','dist_nn'
        ]
        formats = [
            '1J','1J','1J',
            '1E','1E','1E',
            '1E','1E','1E','1E',
            '1E','1E','1E','1E',
            '1E','1E','1E','1E',
            '1E','1E','1E','1E','1E',
            '1J','1E','1E'
        ]
        ## add a fake row with id_inst=1; without it ermldet cannot allocate memory to `tmplist` and likely raise an error (0x42d645,ermldet_in,line 630)
        arrays = [np.append(arr[0],arr) for arr in arrays]
        arrays[1][0] = 1
        
        hdu_lst = fits.HDUList()
        hdu_primary = fits.PrimaryHDU()
        hdu_lst.append(hdu_primary)

        columns = [fits.Column(name=colname_,format=format_,array=array_) for colname_,format_,array_ in zip(colnames,formats,arrays)]
        hdu_data = fits.BinTableHDU.from_columns(columns,name='SRCLIST')
        # copy the header from science_image here
        for key,value in science_head.items():
            hdu_data.header[key] = value
        hdu_lst.append(hdu_data)

        hdu_lst.writeto(outname,overwrite=True)

    else:
        raise Exception('Target source not in eROSITA:DE sky. Aperture photometry cannot be proceeded.')
    
    return


def look_for_confusion(target_ra,target_dec,R_match=15*u.arcsec,R_confusion=60*u.arcsec,cat_ra='RA',cat_dec='DEC'):
    """
    Check if there are source confusion issues.

    Around the target position, we are interested in two regions: 
    - (Region I) The circle of with radius of `R_match`, within which we find our target source in the catalog. For eROSITA `R_match` can be chosen as 15'' (following Merloni+24). 
    - (Region II) The annulus with inner radius of `R_match` and outer radius of `R_confusion`. Sources in this region are not associated with our target. But their PSF wings can extend to the source extraction radius (typically 75% EEF, or ~30''), leading to source confusion issues. eROSITA's HEW is ~30'', so a typical value for `R_confusion` can be 60''.

    Number of source in this Region I should ideally be 1. It doesn't matter if =0 (non detection). >1 means target source lying in very dense regions, which is problematic and current code cannot treat it.  

    Parameters
    ----------

    Returns
    -------
    """
    with fits.open(main_cat) as hdu:
        main_data = hdu[1].data
    with fits.open(supp_cat) as hdu:
        supp_data = hdu[1].data
    data = np.append(main_data,supp_data)
    cat = SkyCoord(data[cat_ra],data[cat_dec],unit='deg',frame='icrs')
    target = SkyCoord([target_ra],[target_dec],unit='deg',frame='icrs')

    print("****** Looking for nearby sources in MAIN+SUPP catalog ********")
    result_in = search_around_sky(target,cat,R_match)
    idx_in = result_in[1]
    sep_in = result_in[2].arcsec
    result_out = search_around_sky(target,cat,R_confusion)
    idx_out = result_out[1]
    sep_out = result_out[2].arcsec
    n_match = len(idx_in)
    n_confusion = len(idx_out) - len(idx_in)

    if n_match == 1:
        print("There is 1 matched source within a matching radius of %s."%(R_match))
    elif n_match == 0:
        print("Target source is not detected in eRASS1 Main+Supp catalog. It could be very faint (det_like_0<5).")
    else:
        print("There are multiple sources matched within a matching radius of %s. Does your target source lie within a dense region (e.g., cluster)?")
        # TODO: how to deal with this situation?

    if n_confusion > 0:
        print("But there is (are) %d confusion source(s) in the annulus of %s ~ %s around target source."%(n_confusion,R_match,R_confusion))
        print("eROSITA's FoV-averaged HEW is ~30 arcsec. Since there is (are) nearby source(s) around target position, the target source extraction region can be contaminated by it (them). Need to generate source map to exclude the contamination(s).")  
    else:
        print("There are no nearby sources in the annulus of %s ~ %s. The source extraction should be safe."%(R_match,R_confusion))

    return n_confusion


def create_aperture(RA,DEC,RE=0.75,RR=0.75,outname='fake_ape.fits'):
    """
    Create a apelist file as input for apetool.

    Parameters
    ----------
    RA : int, float, or numpy.ndarray
    DEC : int, float, or numpy.ndarray
    RE : int, float, or numpy.ndarray, optional
    RR : int, float, or numpy.ndarray, optional
    """
    hdu_lst = fits.HDUList()
    hdu_primary = fits.PrimaryHDU()
    hdu_lst.append(hdu_primary)

    if isinstance(RA,(int,float)):
        RA = [RA]
    if isinstance(DEC,(int,float)):
        DEC = [DEC]
    if isinstance(RE,(int,float)):
        RE = [RE]
    if isinstance(RR,(int,float)):
        RR = [RR]
    colnames = ['RA','DEC','RE','RR']
    arrays = [RA,DEC,RE,RR]
    formats = ['1E','1E','1E','1E']
    columns = [fits.Column(name=colname_,format=format_,array=array_) for colname_,format_,array_ in zip(colnames,formats,arrays)]
    hdu_data = fits.BinTableHDU.from_columns(columns,name='SRCLIST')
    hdu_lst.append(hdu_data)

    hdu_lst.writeto(outname,overwrite=True)
    return


#############################################
############### MISCELLANEOUS ###############
#############################################
def find_erode_skytile(ra,dec,radius=0):
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
    The resulting skytile.
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
                return '{:06d}'.format(data["tiles"][0]["srvmap"])
            else:
                print('Target source not in eROSITA:DE sky!')
                return None
        else:
            print("No tiles found for the given position and radius.")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Error querying the API: {e}")
        return None
    

def erass_flux(ra_target,dec_target,radius=10*u.arcsec,cat='erass1',band='1'):
    '''
    Get source flux from eRASS1 or eRASS:5 catalog.
    '''
    if cat == 'erass1':
        fits_name = main_cat
    elif cat == 'erass:5':
        fits_name = main_cat_5
    else:
        raise Exception('Available catalog: `erass1` or `erass:5` !')
    with fits.open(fits_name) as hdu:
        data = hdu[1].data
    cat = SkyCoord(ra=data['RA'],dec=data['DEC'],unit='deg',frame='icrs')
    target = SkyCoord(ra=[ra_target],dec=[dec_target],unit='deg',frame='icrs')
    result = search_around_sky(target,cat,radius)
    if len(result[1]) == 0:
        print('Target source NONDET in catalog.')
        return None,None,None
    idx = result[1]
    sep = result[2].arcsec # separation in arcsec
    flux_name = 'ML_FLUX_'+band
    expo_name = 'ML_EXP_'+band
    flux = data[idx][flux_name][0]
    expo = data[idx][expo_name][0]
    return flux,expo,sep


def erass1_upperlimit(ra,dec,band='024'):
    '''
    Get source flux upper limit from eRASS1 UpperLimit server.
    '''
    hpix = HEALPix(nside=2**16,order='nested',frame='icrs')
    coord = SkyCoord(ra,dec,unit='deg',frame='icrs')
    hpidx = hpix.skycoord_to_healpix(coord)
    url = f'https://sciserver.mpe.mpg.de/erosita-ul/ULbyHP/{band}/{hpidx}'
    req = requests.get(url)
    assert req.status_code == 200
    if len(req.json()[band]) == 0:
        print('Target source not in eROSITA:DE sky!')
        return None,None,None
    else:
        expo = req.json()[band][0]['Exposure']
        ul_b = req.json()[band][0]['UL_B']
        ul_s = req.json()[band][0]['UL_S']
        return expo,ul_b,ul_s


#############################################
################ DEPRECATED #################
#############################################

# def skytile2boxlist(target_ra,target_dec,outname,science_img,match_radius=15*u.arcsec):
#     """
#     Look for a target in eRASS1 Main+Supp catalogs, and generte a `boxlist` for all detected sources in the skytile where it belongs to.

#     Parameters
#     ----------
#     target_ra : float
#     target_dec : float
#     outname : str
#     science_img : str
#     match_radius : astropy.units.quantity.Quantity

#     Returns
#     -------
#     None
#     """
#     # first try main catalog (det_like_0>6)
#     print("********** Looking for target source in MAIN catalog **********")
#     idx = source_match(main_cat,target_ra,target_dec,match_radius=match_radius)
#     if idx is not None:
#         with fits.open(main_cat) as hdu:
#             data = hdu[1].data
#         skytile = data['SKYTILE'][idx]
#         srcid = data['ID_SRC'][idx]
#         data_inskytile = data[data['SKYTILE']==skytile]
#         print("Target source located in skytile %06d, consisting of %d sources. Target source id: %05d."%(skytile,len(data_inskytile),srcid))
#         # create sub-catalog in target skytile
#         print("********************** Creating boxlist ***********************")
#         fake_boxlist(data_inskytile,science_img,outname)
#     else: # this means det_like_0<6 --> try supp catalog
#         print("Target source undetected in MAIN catalog.")
#         print("********** Looking for target source in SUPP catalog **********")
#         idx = source_match(supp_cat,target_ra,target_dec,match_radius=match_radius)
#         if idx is not None:
#             with fits.open(supp_cat) as hdu:
#                 supp_data = hdu[1].data
#             skytile = supp_data['SKYTILE'][idx]
#             srcid = supp_data['SKYTILE'][idx]
#             supp_data_inskytile = supp_data[supp_data['SKYTILE']==skytile]
#             with fits.open(main_cat) as hdu:
#                 main_data = hdu[1].data
#             main_data_inskytile = main_data[main_data['SKYTILE']==skytile]
#             data_inskytile = np.append(main_data_inskytile,supp_data_inskytile) # concatenate main+supp catalog
#             print("Target source located in skytile %06d, consisting of %d sources. Target source id: %05d."%(skytile,len(data_inskytile),srcid))
#             print("********************** Creating boxlist ***********************")
#             fake_boxlist(data_inskytile,science_img,outname)
#         else:
#             # TODO: find skytile from online api
#             raise Exception("Target source undetected in neither MAIN (det_like_0>6) nor SUPP (5<det_like_0<=6) catalog. Please check if it is in eROSITA-DE sky, or too faint to be seen!")
#     return 


# def source_match(cat_name,target_ra,target_dec,target_unit='deg',target_frame='icrs',match_radius=15*u.arcsec,cat_ra='RA',cat_dec='DEC',cat_unit='deg',cat_frame='icrs'):
#     """
#     Match target source in a given catalog.

#     Parameters
#     ----------
#     """
#     with fits.open(cat_name) as hdu:
#         data = hdu[1].data
#     cat = SkyCoord(data[cat_ra],data[cat_dec],unit=cat_unit,frame=cat_frame)
#     target = SkyCoord([target_ra],[target_dec],unit=target_unit,frame=target_frame)
#     result = search_around_sky(target,cat,match_radius)
#     idx = result[1]
#     sep = result[2].arcsec
#     if len(idx) == 0:
#         print('Target source NONDET in catalog.')
#         return None
#     idx_closest = idx[np.argsort(sep)[0]]
#     return idx_closest