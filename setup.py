#!/usr/bin/env python3

from setuptools import find_packages, setup
import crossmatch

description = ""
"Program for cross-match beteewn two catalogs: Pantheon sample of supernovae Ia and Meta-Catalog of X-ray Detected Cluster of galaxies. "
"Determinination SNe Ia in clusters of galaxies are needed for cosmological analysis and building Hubble diagram"

setup(
    name=crossmatch.__name__,
    version=crossmatch.__version__,
    url='https://github.com/alyonabalakinaSAI/Cross-match',
    author='Elena Balakina',
    author_email='balakina.ea15@physics.msu.ru',
    description=description,
    packages=find_packages(),
    test_suite='test',
    entry_points={
        'console_scripts': ['match = crossmatch.CrossMatch:main']
    },
    include_package_data=True,
    install_requires=['numpy', 'astropy', 'astroquery', 'scipy', 'argparse'],
    classifiers=[
        'Topic :: Astrophysics, Stellar evolution',
        'Programming Language :: Python :: 3.7',
    ],
)
