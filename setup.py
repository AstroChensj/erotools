from setuptools import setup, find_packages

setup(
    name="erotools",
    version="0.1.0",
    description="Python tools for eROSITA data analysis.",
    author="Shi-Jiang Chen",
    author_email="JohnnyCsj666@gmail.com",
    url="https://github.com/AstroChensj/erotools.git",
    packages=["erotools","erotools_scripts","examples"],
    install_requires=[
        "astropy",
        "numpy",
        "scipy",
        "pandas",
        "astropy_healpix",
        # "tqdm",
        # "numba",
        # "joblib",
    ],
    package_data={"erotools": ["*.fits.gz","*.fits"]},
    entry_points={
        "console_scripts": [
            "erosrcdet = erotools_scripts.srcdet:main",
            "eroforcedphot = erotools_scripts.forcedphot:main",
            "eroforcedspec = erotools_scripts.forcedspec:main",
            "erotile_api = erotools_scripts.gettile_api:main",
            "erotile_local = erotools_scripts.gettile_local:main",
            "eroecf = erotools_scripts.geteroecf:main",
        ]
    }
)
