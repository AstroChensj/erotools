#!/usr/bin/env python3
"""
Wrapper function for erotools.eroecf.
Energy conversion factor (ECF) calculation for eROSITA observations, where ECF defined as: rate (Photons/s) / flux (erg/cm^2/s).
User should either supply an xcm file (`restore_file`, from Xset.save or Xspec command line), or a spectral model (`xmodel`) along with detailed settings (`xmodel_par`).
Note that on-axis eROSITA response files are used, unless specified otherwise.

"""
import numpy as np
from astropy.io import fits
import astropy.units as u
import pandas as pd
import argparse
import sys
from erotools.eroecf import get_eroecf


# define argparser
class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write("error: %s\n" % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Shi-Jiang Chen, Johannes Buchner and Teng Liu (C) 2024 <JohnnyCsj666@gmail.com>""",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("emin_rate",type=float,help="emin for rate calculation.")
parser.add_argument("emax_rate",type=float,help="emax for rate calculation.")
parser.add_argument("emin_flux",type=float,help="emin for flux calculation.")
parser.add_argument("emax_flux",type=float,help="emax for flux calculation.")
parser.add_argument("--arf",type=str,default="default_ARF",help="Path to the ARF file. Defaults to `default_ARF`.")
parser.add_argument("--rmf",type=str,default="default_RMF",help="Path to the RMF file. Defaults to `default_RMF`.")
parser.add_argument("--restore_file",type=str,default=None,help="Path to an XSPEC xcm file to restore the model from. If provided, `xmodel` and `xmodel_par` are ignored.")
parser.add_argument("--xmodel",type=str,default="TBabs*powerlaw",help="XSPEC model for the calculation. Defaults to 'TBabs*powerlaw'.")
parser.add_argument("--xmodel_par",type=str,default="1:0.01,2:2.0,3:1",help="Dictionary of model parameters. Defaults to '1:0.01,2:2.0,3:1'. Note that the key must be an integer.")
parser.add_argument("--chatter",type=int,default=0,help="Chatter level in XSPEC. Defaults to 0 (no output information).")
parser.add_argument("--abund",type=str,default="wilm",help="Abund in XSPEC. Defaults to 'wilm'.")
parser.add_argument("--lmod",type=str,default=None,help="XSPEC lmod command, must be a dict with package name as key and package directory as value.")
parser.add_argument("--xset",type=str,default=None,help="XSPEC xset command, must be a dict with key as command name and value as command value.")


def main():
    """
    Main function to parse arguments and call get_eroecf.
    """
    args = parser.parse_args()
    
    # Convert xmodel_par from string to dictionary
    xmodel_par = {int(k): float(v) for k, v in (pair.split(":") for pair in args.xmodel_par.split(","))}
    
    # Call get_eroecf with parsed arguments
    ECF = get_eroecf(
        emin_rate=args.emin_rate,
        emax_rate=args.emax_rate,
        emin_flux=args.emin_flux,
        emax_flux=args.emax_flux,
        arf=args.arf,
        rmf=args.rmf,
        restore_file=args.restore_file,
        xmodel=args.xmodel,
        xmodel_par=xmodel_par,
        chatter=args.chatter,
        abund=args.abund,
        lmod=args.lmod,
        xset=args.xset
    )
    
    print(f"Energy Conversion Factor (ECF): {ECF}")


if __name__ == "__main__":
    main()