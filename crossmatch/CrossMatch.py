#!/usr/bin/env python3

import numpy
import astropy.units as u
import astropy.coordinates as acoords
from astropy.io import ascii
from astroquery.vizier import Vizier
import astropy.table as astable
import os
import sys
from argparse import ArgumentParser


def crossmatch(catalog1, catalog2, colRA1='RA', colDec1='DEC', colZ1='z',
               colRA2='RA', colDec2='DEC', colZ2='z', max_sep=10, use_z=False):
    """
    NAME:
         Crossmatch
    PURPOSE:
         cross-match two catalogs (including proper motion in catalog2 if epochs are different)
    INPUT:
         catalog1 - First catalog
         name1 - objects name in catalog1
         colRA1= ('RA') name of the tag in catalog1 with the right ascension in degree in catalog1 (assumed to be ICRS)
         colDec1= ('DEC') name of the tag in catalog1 with the declination in degree in catalog1 (assumed to be ICRS)
         colZ1= ('z') name of the tag in catalog1 with the value of redshift

         catalog2 - Second catalog
         name2 - objects name in catalog2
         colRA2= ('RA') name of the tag in catalog2 with the right ascension in degree in catalog2 (assumed to be ICRS)
         colDec2= ('DEC') name of the tag in catalog2 with the declination in degree in catalog2 (assumed to be ICRS)
         colZ2= ('z') name of the tag in catalog2 with the value of redshift

         max_sep = (10) maximum value of angular separation in arcmin
    OUTPUT:
         (joined catalog of matching catalogs, additional column with redshift separation)

    """
    work_catalog1 = acoords.SkyCoord(ra=catalog1[colRA1], dec=catalog1[colDec1],
                                     unit=(u.degree, u.degree))
    work_catalog2 = acoords.SkyCoord(ra=catalog2[colRA2], dec=catalog2[colDec2],
                                     unit=(u.hourangle, u.degree))

    idx, d2d, d3d = work_catalog1.match_to_catalog_sky(work_catalog2)
    d2d = Column(d2d, name='Angle_Sep')
    catalog2_sort = catalog2[idx]

    joined_catalog = astable.hstack([catalog1, catalog2_sort])
    max_separ_all, d3d_all = [], []
    for z in joined_catalog['zCMB']:  # R200 + 5 sigma = 2.92 Mpc, H0=69.8
        max_separ = (2.92 * 69.8 * 206265) / (z * 299792.458 * math.pow((z + 1),2))
        max_separ_all.append(max_separ)
    max_separ_col = Column(max_separ_all * u.arcsec, unit=u.deg, name='Max_Sep')
    joined_catalog = astable.hstack([joined_catalog, max_separ_col, d2d])
    diff = joined_catalog[d2d < max_separ_col]

    if use_z is True:
        delta_z = numpy.abs(diff[colZ1] - diff[colZ2])
        delta_z.name = 'DELTA_Z'
        diff = astable.hstack([diff, delta_z])
        velocity = delta_z * 299792.458
        velocity.name = 'Vel'
        diff = astable.hstack([diff, velocity])
        diff = diff[diff['DELTA_Z'] < 0.02]
    else:
        delta_z = 0
    i = 0
    for z in diff['zCMB']:
        d3d = (z * 299792.458 / 69.8) * (diff['Angle_Sep'][i] / 57.296)
        d3d_all.append(d3d)
        i = i + 1
    d3d_col = Column(d3d_all * u.Mpc, unit=u.Mpc, name='Mpc_Sep')
    diff = astable.hstack([diff, d3d_col])


def parser():
    """Parse command line arguments"""
    parser = ArgumentParser()
    parser.add_argument(
        '-s', '--sep', default=15, type=float,
        help='Maximum value of angular separation in arcmin')
    parser.add_argument(
        '-z', default=True, type=bool,
        help='Set False to not use z separation criterion')

    return parser


def main():

    path = os.path.dirname(os.path.abspath(__file__))

    sne_table = ascii.read(os.path.join(path, 'Pantheon.FITRES'))
    sne_table = sne_table[sne_table['RA'] != 0]

    vizier = Vizier(row_limit=2000)
    mcxc_table = vizier.get_catalogs('J/A+A/534/A109/mcxc')[0]

    namespace = parser().parse_args(sys.argv[1:])
    sep = namespace.sep
    z = namespace.z

    table = crossmatch(catalog1=sne_table, catalog2=mcxc_table, colRA1='RA', colDec1='DECL', colZ1='zCMB',
                       colRA2='RAJ2000', colDec2='DEJ2000', colZ2='z', max_sep=sep, use_z=z)
    final_table = table['CID', 'MCXC', 'RAJ2000', 'RA', 'DEJ2000', 'DECL', 'DELTA_Z']
    final_table['CID'].name, final_table['MCXC'].name = 'SN name', 'Cluster'
    
    ascii.write(final_table, os.path.join(path, 'CrossTable.csv'), format='csv', overwrite=True)


if __name__ == '__main__':
    main()
