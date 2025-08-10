import numpy as np
from astropy.io import fits
import os
import sys
from xspec import *

# on-axis eROSITA response files
default_ARF = os.path.join(os.path.dirname(__file__),"onaxis_tm0_arf_filter_2023-01-17.fits.gz")
default_RMF = os.path.join(os.path.dirname(__file__),"onaxis_tm0_rmf_2023-01-17.fits.gz")


def get_eroecf(emin_rate,emax_rate,emin_flux,emax_flux,arf="onaxis",rmf="onaxis",restore_file=None,xmodel="tbabs*powerlaw",xmodel_par={1:0.01,2:2.0,3:1},chatter=0,abund="wilm",lmod=None,xset=None):
    """
    Energy conversion factor (ECF) calculation for eROSITA observations, where ECF defined as: rate (Photons/s) / flux (erg/cm^2/s).
    User should either supply an xcm file (`restore_file`, from Xset.save or Xspec command line), or a spectral model (`xmodel`) along with detailed settings (`xmodel_par`).
    Note that on-axis eROSITA response files are used, unless specified otherwise.

    Parameters
    ----------
    emin_rate : float
        emin for rate calculation.
    emax_rate : float
        emax for rate calculation.
    emin_flux : float
        emin for flux calculation.
    emax_flux : float
        emax for flux calculation.
    arf : str, optional
        Path to the ARF file. Defaults to 'onaxis', i.e., using the on-axis ARF.
    rmf : str, optional
        Path to the RMF file. Defaults to 'onaxis', i.e., using the on-axis RMF.
    restore_file : str, optional
        Path to an XSPEC xcm file to restore the model from. If provided, `xmodel` and `xmodel_par` are ignored.
    xmodel : str
        XSPEC model for the calculation. Defaults to "TBabs*powerlaw".
    xmodel_par : dict
        Dictionary of model parameters. Defaults to {1:0.01,2:2.0,3:1}. Note that the key must be an integer.
    chatter : int, optional
        Chatter level in XSPEC. Defaults to 0 (no output information).
    abund : str, optional
        Abund in XSPEC. Defaults to "wilm".
    lmod : dict, optional
        XSPEC lmod command, must be a dict with package name as key and package directory as value.
    xset : dict, optional
        XSPEC xset command, must be a dict with key as command name and value as command value.

    Returns
    -------
    ECF : float
        Energy conversion factor.

    """
    # input settings
    emin_rate = f"{emin_rate:.2f}"
    emax_rate = f"{emax_rate:.2f}"
    emin_flux = f"{emin_flux:.2f}"
    emax_flux = f"{emax_flux:.2f}"
    Xset.chatter = chatter
    Xset.abund = abund
    if lmod is not None:
        for package,package_dir in lmod.items():
            AllModels.lmod(package,package_dir)
    if xset is not None:
        for key,val in xset.items():
            Xset.addModelString(key,val)
    if arf == "onaxis":
        arf = default_ARF
    if rmf == "onaxis":
        rmf = default_RMF
    # set model
    AllModels.clear()
    if restore_file is not None:    # we can either load a xcm file
        # check for PyXspec indicator in first line
        xcm_content = open(restore_file,"r").read()
        if not xcm_content[0].startswith("#PyXspec"):
            print("Not a PyXspec log file. Append #PyXspec to first line.")
            with open(restore_file,"w") as f:
                f.write(f"#PyXspec: Output generated from Xset.save().  DO NOT MODIFY.\n{xcm_content}")
        Xset.restore(restore_file)
        m1 = AllModels(1)
        xmodel = m1.expression  # read model expression
    else:   # or specify the model manually
        m1 = Model(xmodel)
        m1.setPars(xmodel_par)

    AllData.clear()
    AllData.fakeit(1,FakeitSettings(response=rmf,arf=arf,exposure=10000),noWrite=True)
    AllData.ignore(f"0.-{emin_rate} {emax_rate}-**")
    # calculate model flux
    AllModels.calcFlux(f"{emin_flux} {emax_flux}") # to get emin -- emax flux, you must first run calcFlux
    flux = AllData(1).flux[0]     # (value, errLow, errHigh (in ergs/cm^2/s), value, errLow, errHigh (in photons/cm^2/s))
    # calculate count rate
    rate = AllData(1).rate[3]       # (current net rate, net rate uncertainty, total rate, predicted model rate)
    ECF = rate / flux

    print("################# ECF CALCULATION #################")
    print(f"Assuming RMF: {rmf}")
    print(f"Assuming ARF: {arf}")
    print(f"Assuming model: {xmodel}")
    for i in range(m1.nParameters):
        print(f"{m1(i+1).name}:\t\t{m1(i+1).values[0]}")
    print(f"{emin_rate} -- {emax_rate} keV count rate: {rate} Photons/s")
    print(f"{emin_flux} -- {emax_flux} keV model flux: {flux} erg/cm^2/s")
    print(f"ECF (rate/flux): {ECF:.3e}")

    return ECF