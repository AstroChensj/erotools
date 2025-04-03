from setuptools import setup, find_packages

setup(
    name='erotools',
    version='0.1.0',
    description='Python tools for eROSITA data analysis.',
    author='Shi-Jiang Chen',
    author_email='JohnnyCsj666@gmail.com',
    #url='https://github.com/AstroChensj/Xstack.git',
    packages=['erotools','erotools_scripts','examples'],
    install_requires=[
        'astropy',
        'numpy',
        'scipy',
        'pandas',
        'astropy_healpix',
        # 'tqdm',
        # 'numba',
        # 'joblib',
    ],
    package_data={'erotools': ["*.fits.gz",'*.fits']},
    entry_points={
        'console_scripts': [
            'run_apetool = erotools_scripts.apetool_autoscript:main',
            'get_ero_skytile = erotools_scripts.get_ero_skytile:main'
        ]
    }
)
