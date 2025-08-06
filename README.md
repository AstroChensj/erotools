# Helpful tools for eROSITA

These include:
- `eroforcedphot`: performing forced photometry on eROSITA data.
  - This is a wrapper code for [apetool](https://erosita.mpe.mpg.de/dr1/eSASS4DR1/eSASS4DR1_tasks/apetool_doc.html).
- `eroforcedspec`: performing forced spectroscopy on eROSITA data.
  - This is a wrapper code for [srctool](https://erosita.mpe.mpg.de/dr1/eSASS4DR1/eSASS4DR1_tasks/srctool_doc.html).
- `erotile_api` and `erotile_local`: get the sky tile (eROSITA) for a given ra, dec.
  - This is copied from [Jeremy Sanders' code](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/apis.html), used for other functions in this code.

## Installation

```shell
git clone https://github.com/AstroChensj/erotools.git
cd erotools
python -m pip install .  # this should set alias needed
```

You would also need the eRASS1 source catalog (both MAIN and SUPP) for forced spectroscopy. These two catalogs can be downloaded from eRASS1 archive:
- [MAIN catalog](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Main.tar.gz)
- [SUPP catalog](https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Supp.tar.gz)

Move the two catalogs under your favorite folder, whose name should be stored as an environmental variable (`eRASS_CAT_DIR`):
```shell
mv eRASS1_Main.tar.gz eRASS1_Supp.tar.gz /path/to/your/catalogs
export eRASS_CAT_DIR=/path/to/your/catalogs
```


## Usage
More details to be filled ...

But for now, you could check the shell command doc for basic introduction, e.g., `eroforcedphot -h`.

### Forced photometry
```shell
eroforcedphot /path/to/your/science_events_file target_ra target_dec
```
- `eroforcedphot` is an alias for `python ./erotools_scripts/forcedphot.py`.
- Replace `target_ra` and `target_dec` with your target RA and DEC.

### Forced spectroscopy
```shell
eroforcedspec /path/to/your/science_events_file target_ra target_dec --redshift 0.1 --record_redshift --galnh 1e20 --record_galnh --prefix ./out/forcedspec_
```
- `eroforcedspec` is an alias for `python ./erotools_scripts/forcedspec.py`.
- `/path/to/your/science_events_file` is absolute path to the event file you download from eRASS1 archive. The event file is sth like `051/108/EXP_010/eb01_108051_020_EventList_c010.fits.gz`.
- Replace `target_ra` and `target_dec` with your target RA and DEC.
- Replace `0.1` with the real redshift, and `1e20` with the real Galactic NH (in units of 1 cm^-2). These two information is not necessary for forced spectroscopy, but could be useful for later spectral stacking, with e.g., [`Xstack`](https://github.com/AstroChensj/Xstack).
- Replace `./out/forcedspec_` with the desired output prefix.

The output will be SPEC, BKGSPEC, ARF, RMF, alongside the redshift and galnh files (if `--record_redshift` and `--record_galnh`).


## Credit
Credits to eSASS team, especially Jeremy Sanders, Antonis Georgakakis, Miriam Ramos-Ceja, and many others.

Part of the tools are available from [eSASS page](https://erosita.mpe.mpg.de/dr1/eSASS4DR1/eSASS4DR1_tasks/).
