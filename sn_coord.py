import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astroquery.vizier import Vizier
from astropy.cosmology import FlatLambdaCDM
import astropy.table as astable
import numpy as np

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

sne_table = ascii.read('Pantheon.FITRES')  # load pantheon data
sne_table = sne_table[sne_table['RA'] != 0]   # only RA != 0 data
sne_dist = cosmo.comoving_distance(sne_table['zCMB'])   # calculate distanse 
sne = SkyCoord(ra=sne_table['RA'], dec=sne_table['DECL'], distance=sne_dist, 
               unit=(u.deg, u.deg, u.Mpc))
 
vizier = Vizier(row_limit=2000)  # load mcxc from vizier
mcxc_table = vizier.get_catalogs('J/A+A/534/A109/mcxc')[0]   # mcxc data
mcxc_dist = cosmo.comoving_distance(mcxc_table['z'])       # calculate distanse 
mcxc = SkyCoord(ra=mcxc_table['RAJ2000'], dec=mcxc_table['DEJ2000'], distance=mcxc_dist,
                unit=(u.hourangle, u.deg, u.Mpc))

mcxc_idx, d2d, d3d = sne.match_to_catalog_3d(mcxc)  # cross-match sne.match_to_catalod_3d
mcxc_sorted=mcxc_table[mcxc_idx]       # sort mcxc like pantheon
joined_table = astable.hstack([sne_table, mcxc_sorted]) # mcxc + pantheon
res_sne = joined_table[d3d < 5*u.Mpc]   # chose good sne with clusters
delta_z = np.abs(res_sne['zCMB'] - res_sne['z'])
great_table = astable.hstack([res_sne, delta_z])     # full good data
fin_table = great_table[delta_z < 0.015]
short_table = fin_table['CID', ['MCXC'], 'delta_z']

# short_table.write('short_table.dat', format = 'ascii')
print(short_table)
