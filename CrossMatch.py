import warnings
import numpy
import astropy.units as u
import astropy.coordinates as acoords
from astropy.io import ascii
from astroquery.vizier import Vizier
from astropy.cosmology import FlatLambdaCDM
import astropy.table as astable

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


def check_epoch(catalog, epoch):
    """
    NAME:
         check_epoch
    PURPOSE:
         determination is there difference beteewn epochs in the first and the second catalog
    INPUT:
         catalog - catalog used in Crossmatch function
         epoch - epoch of observation in catalog
    """

    warn_about_epoch = False
    if 'ref_epoch' in catalog.dtype.fields:
        if 'designation' not in catalog.dtype.fields:  # Assume this is DR1
            if numpy.any(numpy.fabs(epoch - 2015.) > 0.01):
                warn_about_epoch = True
        elif 'Gaia DR2' in catalog['designation'][0].decode('utf-8'):
            if numpy.any(numpy.fabs(epoch - 2015.5) > 0.01):
                warn_about_epoch = True
    if warn_about_epoch:
        warnings.warn("You appear to be using a Gaia catalog, set the epoch (DR1) or 2015.5 (DR2)")
    return None


def Crossmatch(catalog1, name1, catalog2, name2, colRA1='RA', colDec1='DEC', colZ1='z', epoch1=None,
               colRA2='RA', colDec2='DEC', colZ2='z', epoch2=None, max_sep=10,
               colpmRA2='pmra', colpmDec2='pmdec',
               swap=False, use_z=False, use_pm=False):
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
         epoch1= (2000.) epoch of the coordinates in catalog1

         catalog2 - Second catalog
         name2 - objects name in catalog2
         colRA2= ('RA') name of the tag in catalog2 with the right ascension in degree in catalog2 (assumed to be ICRS)
         colDec2= ('DEC') name of the tag in catalog2 with the declination in degree in catalog2 (assumed to be ICRS)
         colZ2= ('z') name of the tag in catalog2 with the value of redshift
         epoch2= (2000.) epoch of the coordinates in catalog2

         colpmRA2= ('pmra') name of the tag in catalog2 with the proper motion in right ascension in degree in catalog2 (assumed to be ICRS; includes cos(Dec)) [only used when epochs are different]
         colpmDec2= ('pmdec') name of the tag in catalog2 with the proper motion in declination in degree in catalog2 (assumed to be ICRS) [only used when epochs are different]
         max_sep = (10) maximum value of angular separation in arcmin
         swap= (False) if False, find closest matches in catalog2 for each catalog1 source, if True do the opposite (important when one of the catalogs has duplicates)
    OUTPUT:
         (index into catalog1 of matching objects,
            index into catalog2 of matching objects,
            angular separation between matching objects,
            separation in redshift between matching objects)

    """

    if epoch1 is None:
        if 'ref_epoch' in catalog1.dtype.fields:
            epoch1 = catalog1['ref_epoch']
        else:
            epoch1 = 2000.0
    if epoch2 is None:
        if 'ref_epoch' in catalog2.dtype.fields:
            epoch2 = catalog2['ref_epoch']
        else:
            epoch2 = 2000.0
    check_epoch(catalog1, epoch1)
    check_epoch(catalog2, epoch2)
    depoch = epoch2 - epoch1
    if numpy.any(depoch != 0.):
        # Use proper motion to get both catalogs at the same time
        dra = catalog2[colpmRA2] / numpy.cos(catalog2[colDec2] / 180.0 * numpy.pi)\
            / 3600000.0 * depoch
        ddec = catalog2[colpmDec2] / 3600000.0 * depoch
    else:
        dra = 0.0
        ddec = 0.0
    work_catalog1 = acoords.SkyCoord(catalog1[colRA1], catalog1[colDec1],
                                     unit=(u.degree, u.degree), frame='icrs')
    if use_pm is True:
        work_catalog2 = acoords.SkyCoord(catalog2[colRA2] - dra, catalog2[colDec2] - ddec,
                                         unit=(u.degree, u.degree), frame='icrs')
    else:
        work_catalog2 = acoords.SkyCoord(catalog2[colRA2], catalog2[colDec2],
                                         unit=(u.degree, u.degree), frame='icrs')
    if swap:
        idx, d2d, d3d = work_catalog2.match_to_catalog_sky(work_catalog1)
    else:
        idx, d2d, d3d = work_catalog1.match_to_catalog_sky(work_catalog2)

    catalog2 = catalog2[idx]

    if use_z is True:
        delta_z = numpy.abs(catalog1[colZ1] - catalog2[colZ2])
    else:
        delta_z = 0

    catalog2 = catalog2[delta_z < 0.015]
    if swap:
        return (catalog2)
    else:
        return (catalog2)


sne_table = ascii.read('Pantheon.FITRES')
sne_table = sne_table[sne_table['RA'] != 0]
name_sne = sne_table['CID']
vizier = Vizier(row_limit=2000)
mcxc_table = vizier.get_catalogs('J/A+A/534/A109/mcxc')[0]
name_mcxc = mcxc_table['MCXC']

table = Crossmatch(sne_table, name_sne, mcxc_table, name_mcxc, 'RA', 'DECL', 'zCMB', 2000.0, 'RAJ2000', 'DEJ2000', 'z', 2000.0, 1000, use_z=True)

fin_table = table['MCXC', 'RAJ2000', 'DEJ2000']
full_table = astable.hstack([fin_table, name_sne])
full_table.write('TEST', format='ascii')
