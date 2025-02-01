import numpy as np
from astropy.io import fits
import os
import sys
from xspec import *

default_ARF = os.path.join(os.path.dirname(__file__),'onaxis_tm0_arf_filter_2023-01-17.fits.gz')
default_RMF = os.path.join(os.path.dirname(__file__),'onaxis_tm0_rmf_2023-01-17.fits.gz')


def get_eroecf(emin_rate,emax_rate,emin_flux,emax_flux,Gamma=2.0,logNH=0,chatter=0,abund="wilm"):
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
    emin_rate = '{:.2f}'.format(emin_rate)
    emax_rate = '{:.2f}'.format(emax_rate)
    emin_flux = '{:.2f}'.format(emin_flux)
    emax_flux = '{:.2f}'.format(emax_flux)
    Xset.chatter = chatter
    Xset.abund = abund
    xmodel = "TBabs*powerlaw"
    # set model
    AllModels.clear()
    m1 = Model(xmodel)
    m1.TBabs.nH.values = 10**logNH*1e-22
    m1.powerlaw.PhoIndex.values = Gamma
    AllData.clear()
    AllData.fakeit(1,FakeitSettings(response=default_RMF,arf=default_ARF,exposure=100),noWrite=True)
    AllData.ignore("0.-%s %s-**"%(emin_rate,emax_rate))
    # calculate model flux
    AllModels.calcFlux("%s %s"%(emin_flux,emax_flux)) # to get emin -- emax flux, you must first run calcFlux
    flux = AllData(1).flux[0]     # (value, errLow, errHigh (in ergs/cm^2/s), value, errLow, errHigh (in photons/cm^2/s))
    # calculate count rate
    rate = AllData(1).rate[3]       # (current net rate, net rate uncertainty, total rate, predicted model rate)
    ECF = rate / flux

    print("################# ECF CALCULATION #################")
    print("Assuming model: %s"%(xmodel))
    print("Gamma: %f"%(Gamma))
    print("logNH: %f"%(logNH))
    print("%s -- %s keV count rate: %f Photons/s"%(emin_rate,emax_rate,rate))
    print("%s -- %s keV model flux: %e erg/cm^2/s"%(emin_flux,emax_flux,flux))
    print("ECF (rate/flux): %.3e"%(ECF))

    return ECF