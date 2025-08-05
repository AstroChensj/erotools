# Helpful tools for eROSITA

These include:
- `eroforcedphot`: performing forced photometry on eROSITA data;
- `eroforcedspec`: performing forced spectroscopy on eROSITA data;
- `erotile_api` and `erotile_local`: get the sky tile (eROSITA) for a given ra, dec;

## Installation

```shell
git clone
cd
python -m pip install .
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

### Forced spectroscopy
```shell
eroforcedspec /path/to/your/science_events_file target_ra target_dec --redshift 0.1 --record_redshift --galnh 1e20 --record_galnh
```
- `/path/to/your/science_events_file` is absolute path to the event file you download from eRASS1 archive. The event file is sth like `051/108/EXP_010/eb01_108051_020_EventList_c010.fits.gz`.
- replace `target_ra` and `target_dec` with your target RA and DEC.
- assuming the target source has redshift of 0.1, and assuming the Galactic NH of 1e20 cm^-2. These two information is not necessary for forced spectroscopy, but could be useful for later spectral stacking, with e.g., [`Xstack`](https://github.com/AstroChensj/Xstack).

The output will be SPEC, BKGSPEC, ARF, RMF, alongside the redshift and galnh files (if `--record_redshift` and `--record_galnh`).


## Credit
Credits to eSASS team, especially Jeremy Sanders, Antonis Georgakakis, Miriam Ramos-Ceja, etc.

Part of the tools are available from [eSASS page](https://erosita.mpe.mpg.de/dr1/eSASS4DR1/eSASS4DR1_tasks/).
