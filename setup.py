#!/usr/bin/env python3

from setuptools import find_packages, setup, Extension
import os
# import glob

packagename = 'balakina'

# Cython related
try:
    from Cython.Build import cythonize
    import numpy as np
    ext_modules = cythonize([
        Extension(
            'balakina.f_speed.cy_contfrac',
            ['src/cy_contfrac.pyx'],
            include_dirs=[np.get_include()]),
        Extension(
            'balakina.f_speed.erf',
            ['src/erf.pyx'],
            include_dirs=[np.get_include()]),
    ])
except ImportError:
    ext_modules = []

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'requirements.txt')) as f:
    install_requires = f.read().splitlines()

description = """
Program for cross-match beteewn two catalogs: Pantheon sample of supernovae Ia and Meta-Catalog of X-ray Detected Cluster of galaxies.
Determinination SNe Ia in clusters of galaxies are needed for cosmological analysis and building Hubble diagram
""".replace('\n', ' ')

setup(
    name=packagename,
    version='0.0.1',
    url='https://github.com/alyonabalakinaSAI/Cross-match',
    author='Elena Balakina',
    author_email='balakina.ea15@physics.msu.ru',
    description=description,
    packages=find_packages(),
    scripts=['/CrossMatch',
             'bin/sci_py_import_all',
             'bin/sci_py_test_style'],
    data_files=[
        ('FITRES', ['Pantheon.FITRES']),
    ],
    install_requires=install_requires,
    classifiers=[
        'Intended Audience :: Science',
        # 'License :: OSI Approved :: MIT License',
        'Topic :: Astrophysics, Stellar evolution',

        'Programming Language :: Python :: 3.7',
    ],
    keywords='sample education',
)
