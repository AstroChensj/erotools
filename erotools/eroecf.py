import numpy as np
from astropy.io import fits
import os
import sys
from xspec import *

default_ARF = os.path.join(os.path.dirname(__file__),'onaxis_tm0_arf_filter_2023-01-17.fits.gz')
default_RMF = os.path.join(os.path.dirname(__file__),'onaxis_tm0_rmf_2023-01-17.fits.gz')


def get_eroecf(emin_rate,emax_rate,emin_flux,emax_flux,xmodel="tbabs*powerlaw",xmodel_par={1:0.01,2:2.0,3:1},chatter=0,abund="wilm",lmod=None,xset=None):
    """
    Energy conversion factor (ECF) calculation, where ECF defined as: rate (Photons/s) / flux (erg/cm^2/s).
    Assuming model for the conversion: TBabs*powerlaw

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
    Gamma : float, optional
        Powerlaw photon index. Defaults to 2.0.
    logNH : float, optional
        Logarithm of Galactic NH. Defaults to 0.
    chatter : int, optional
        Chatter level in XSPEC. Defaults to 0 (no output information).
    abund : str, optional
        Abund in XSPEC. Defaults to "wilm".

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
    # xmodel = "TBabs*powerlaw"
    # set model
    AllModels.clear()
    m1 = Model(xmodel)

    # m_vals = {}
    # # initialize the component model
    # for paridx,parval in xmodel_par.items():
    #     paridx = int(paridx)    # the key must be an integer
    #     # if isinstance(chain_idx,str) and len(chain_idx.split("__"))>1:   # e.g., PhoIndex__1
    #     #     chain_val = chain_data[chain_idx][i]
    #     # else:   # e.g., 0.01234 or "1 0.1 0.1 0.1 500 500"
    #     # chain_val = chain_idx
    #     m_vals[comp_idx] = chain_val
    m1.setPars(xmodel_par)

    # m1.TBabs.nH.values = 10**logNH*1e-22
    # m1.powerlaw.PhoIndex.values = Gamma
    AllData.clear()
    AllData.fakeit(1,FakeitSettings(response=default_RMF,arf=default_ARF,exposure=100),noWrite=True)
    AllData.ignore(f"0.-{emin_rate} {emax_rate}-**")
    # calculate model flux
    AllModels.calcFlux(f"{emin_flux} {emax_flux}") # to get emin -- emax flux, you must first run calcFlux
    flux = AllData(1).flux[0]     # (value, errLow, errHigh (in ergs/cm^2/s), value, errLow, errHigh (in photons/cm^2/s))
    # calculate count rate
    rate = AllData(1).rate[3]       # (current net rate, net rate uncertainty, total rate, predicted model rate)
    ECF = rate / flux

    print("################# ECF CALCULATION #################")
    print(f"Assuming model: {xmodel}")
    for i in range(m1.nParameters):
        print(f"{m1(i+1).name}:\t\t{m1(i+1).values[0]}")
    print(f"{emin_rate} -- {emax_rate} keV count rate: {rate} Photons/s")
    print(f"{emin_flux} -- {emax_flux} keV model flux: {flux} erg/cm^2/s")
    print(f"ECF (rate/flux): {ECF:.3e}")

    return ECF