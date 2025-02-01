from astropy.table import Table
from scipy.stats import poisson
from scipy.special import gammainc, gammaincc, gammaincinv 
import scipy.optimize as optimize
import astropy.units as u
import numpy as np


def cal_detlike(tot,bkg):
    """
    Calculate detection likelihood, given the total counts (`tot`) and expected bkg counts (`bkg`) in a certain region. The bkg estimation is assumed to be infinitely precise (it is estimated from a very large area with sufficient photon counts).

    Our question is at what confidence level can we claim we detect sth. The null hypothesis is that the `tot` comes entirely from the fluctuation of background. The received bkg counts follows Poisson distribution with mean of `bkg`. What we really observe now is `tot`. The significance ($\alpha$) of such observation, assuming null hypothesis is correct, is calculated as:

    $\alpha = 1 - \sum_{k=0}^{tot} bkg^k * e^{-bkg} / k! = P_\Gamma(tot+1,bkg)$

    where $P_\Gamma$ is the regularized incomplete gamma function. $\alpha$ tells us the frequency that we observe counts as large or larger than `tot`. Our significance is detecting sth (detlike) is related to $\alpha$:

    $detlike = -1*\log\alpha$

    A rule of thumb for detection is detlike>=5 (corresponding to significance of ~0.1), as used by eRASS1 Main+Supp catalog. In this case, the source count rate and uncertainty can be calculated by the `counts_twoside_poisson` function below. 

    Otherwise if detlike<5, it should be considered as non-detection, and an upper limit for counts can be calculated via `counts_ul_poisson` below.


    Parameters
    ----------
    tot : float
        Total counts within an aperture.
    bkg : float
        Expected bkg counts in the same aperture.

    Returns
    -------
    detlike : float
        The detection likelihood. Examples:
        - detlike=3 --> significance~0.05
        - detlike=5 --> significance~7e-3
        - detlike=6 --> significance~2e-3
    """
    detlike = -1 * np.log(gammainc(tot+1,bkg))
    return detlike


def counts_twoside_poisson(tot,bkg,conf_limit=0.68):
    """
    Calculate the source net counts (median) and two side uncertainty, based on the source count cdf:

    $cdf(src) = (\gamma(tot+1,src+bkg) - \gamma(tot+1,bkg)) / \Gamma(tot+1,bkg)$

    Parameters
    ----------
    tot : float
        Total number of counts with the aperture.
    bkg : float
        Background level within the aperture.
    conf_limit : float
        Confidence level. Will be divided equally to find lower error and upper error.

    Returns
    -------
    src_me : float
        Source net counts.
    src_lo : float
        Lower error of source counts. 
    src_hi : float
        Upper error of source counts.
    
    """
    src_me = counts_ul_poisson(tot,bkg,0.5)
    src_left = counts_ul_poisson(tot,bkg,0.5-conf_limit/2)
    src_right = counts_ul_poisson(tot,bkg,0.5+conf_limit/2)
    src_lo = src_me - src_left
    src_hi = src_right - src_me
    
    return src_me,src_lo,src_hi


def counts_ul_poisson(tot,bkg,conf_limit=0.954):
    '''
    the upper limit function. Credit : A. Ruiz.
    Correspomds to Eq 8 of Ruiz et al. (2022)
    https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4265R/abstract
    
    Parameters
    ----------
    tot : float
        Total number of counts with the aperture.
    bkg : float
        Background level within the aperture.
    conf_limit : float
        Confidence limit of the upper limit. For example:
        - 0.8412 correspond to the one-sided 1-sigma confidence level or a normal distribution (0.50+0.34)
        - 0.9772 : one-sided 2-sigma
        - 0.9987 : one-sided 3-sigma
        
    Returns
    -------
    crul : float
        Lpper limit in units of counts at the confidence level conf_limit. 
        No correction is applied for the EEF of the extraction aperture. 
        This needs to be applied separately. 
    
    '''
    crul_integrand = conf_limit * gammaincc(tot + 1, bkg) + gammainc(tot + 1, bkg)
    crul = gammaincinv(tot + 1, crul_integrand) - bkg
    return crul 
