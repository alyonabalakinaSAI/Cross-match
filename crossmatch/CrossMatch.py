#!/usr/bin/env python3

import numpy
import astropy.units as u
import astropy.coordinates as acoords
from astropy.io import ascii
# from astroquery.vizier import Vizier
import astropy.table as astable
import os


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
    catalog2_sort = catalog2[idx]

    joined_catalog = astable.hstack([catalog1, catalog2_sort])
    joined_catalog = joined_catalog[d2d < (max_sep * u.arcmin)]

    if use_z is True:
        delta_z = numpy.abs(joined_catalog[colZ1] - joined_catalog[colZ2])
        delta_z.name = 'DELTA_Z'
        joined_catalog = astable.hstack([joined_catalog, delta_z])
        joined_catalog = joined_catalog[joined_catalog['DELTA_Z'] < 0.015]
    else:
        delta_z = 0

    return (joined_catalog)


def main():
    # path = os.path.dirname(os.path.abspath(__file__))

    # sne_table = ascii.read(os.path.join(path, 'Pantheon.FITRES'))
    # sne_table = sne_table[sne_table['RA'] != 0]

    # vizier = Vizier(row_limit=2000)
    # mcxc_table = vizier.get_catalogs('J/A+A/534/A109/mcxc')[0]

    # table = crossmatch(catalog1=sne_table, catalog2=mcxc_table, colRA1='RA', colDec1='DECL', colZ1='zCMB', colRA2='RAJ2000', colDec2='DEJ2000', colZ2='z', max_sep=20, use_z=True)
    # final_table = table['CID', 'MCXC', 'RAJ2000', 'RA', 'DEJ2000', 'DECL', 'DELTA_Z']
    # final_table['CID'].name, final_table['MCXC'].name = 'SN name', 'Cluster'
    # print(final_table)
    # ascii.write(final_table, os.path.join(path, 'CrossTable'), format='csv', overwrite=True)

    # testfile = os.path.dirname(os.path.abspath(__file__))
    # test_cat_path = os.path.join('../../../ ', testfile, 'crossmatch/testcatalog')
    cat1 = ascii.read(os.path.join('cat1.csv'), delimiter=',', guess=False, fast_reader=False)
    cat2 = ascii.read(os.path.join('cat2.csv'), delimiter=',', guess=False, fast_reader=False)
    table = crossmatch(cat1, cat2, 'RA1', 'DEC1', 'Z1', 'RA2', 'DEC2', 'Z2', 20, use_z=True)
    final_table = table['NAME1', 'NAME2', 'RA1', 'RA2', 'DEC1', 'DEC2', 'DELTA_Z']
    ascii.write(final_table, 'TEST')


if __name__ == '__main__':
    main()
