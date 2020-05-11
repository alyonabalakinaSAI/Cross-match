#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib import *
from matplotlib.colors import ListedColormap

# Plot Setup
rcParams['font.size'] = 16
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True
font = {'family': 'normal', 'size': 16}
rc('axes', linewidth=1)
rc("text", usetex=True)
rc('font', family='serif')
rc('font', serif='Times')
rc('legend', fontsize=14)
rc('xtick.major', size=5, width=1.5)
rc('ytick.major', size=5, width=1.5)
rc('xtick.minor', size=3, width=1)
rc('ytick.minor', size=3, width=1)


def sigma_velocity(redshift, sigma_magnitude, light_speed):
    """
    NAME:
        sigma_velocity
    PURPOSE:
        Determination of peculiar velocity dispersion knowing int. dispersion sigma_mag
        (it is inverse problem, in real calculation we do it vice versa)
    INPUT:
        redshift - cosmological redshift of the supernova [float]
        sigma_magnitude - intrinsic dispersion of the SN magnitude [float]
        light_speed - set value of speed of light [float]
    OUTPUT:
        value of pecular velocity dispersion of host galaxy in cluster [float]
    """
    result = (light_speed * redshift * sigma_magnitude * np.log(10)) / 5
    return result


def sigma_magnitude(redshift, sigma_velocity, light_speed):
    """
    NAME:
        sigma_magnitude
    PURPOSE:
        Determination of intrinsic dispersion of peculiar velocities of SNe Ia
    INPUT:
        redshift - cosmological redshift of the supernova [float]
        sigma_velocity - velocity dispersion of host galaxy in cluster [float]
        light_speed - set value of speed of light [float]
    OUTPUT:
        value of the dispersion in SN magnitude taking into account SN pecular velocity [float]
    """
    result = (5 * sigma_velocity) / (light_speed * redshift * np.log(10))
    return result


def cyan_map():
    """
    NAME:
        cyan_map
    PURPOSE:
        Create my cyan color map for the plot
    """
    n_color = 256
    vals = np.ones((n_color, 4))
    vals[:, 0] = np.linspace(53 / 256, 1, n_color)
    vals[:, 1] = np.linspace(81 / 256, 1, n_color)
    vals[:, 2] = np.linspace(80 / 256, 1, n_color)
    Cyan_map = ListedColormap(vals)
    return Cyan_map


# Set the constant values
N = 300
c = 299792
z_max = 0.20
sigma_m_max = 0.3

z = np.linspace(0.001, z_max, N)
sigma_m = np.linspace(0.0, sigma_m_max, N)
Z, SIGMA_M = np.meshgrid(z, sigma_m)
sigma_vel = sigma_velocity(Z, SIGMA_M, c)

# Plotting value of mag dispersion
Cyan_map = cyan_map()
fig = plt.figure(figsize=(8, 6))
ax = plt.subplot(111)
image = ax.imshow(sigma_vel, interpolation='bilinear', aspect=0.65, cmap=Cyan_map,
                  extent=[0.0, z_max, 0.0, sigma_m_max], origin='lower')
plt.setp(ax.get_yticklabels()[0], visible=False)
plt.setp(ax.get_yticklabels()[-1], visible=False)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_major_locator(MultipleLocator(0.05))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
axins = inset_axes(ax,
                   width="5%",
                   height="100%",
                   loc='lower right',
                   bbox_to_anchor=(0.08, 0.0, 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = fig.colorbar(image, cax=axins)
cbar.ax.set_ylabel('velocity dispersion (km per s)')
ax.set_xlabel(u'redshift')
ax.set_ylabel(r'$\sigma_m$ (mag)')

# Plotting some values of mag dispersion for different velocity values
ax.plot(z, sigma_magnitude(z, 250, c), color='k', linestyle='-', lw=2, label='250 km per s (Pantheon)')
ax.plot(z, sigma_magnitude(z, 1040, c), color='k', linestyle=':', lw=2.5, label='1040 km per s (Perseus)')
ax.plot(z, sigma_magnitude(z, 1545, c), color='k', linestyle='--', lw=2, label='1545 km per s (ACO 2319)')

# Plotting intrinsic dispersion value for two scatter models on Pantheon
ax.hlines(0.10, 0.0, z_max, color='#edc949', linestyle='solid', lw=2, label='C11 scatter model')
ax.hlines(0.08, 0.0, z_max, color='#edc949', linestyle='--', lw=2, label='G10 scatter model')

ax.axis([0.0, z_max, 0.0, sigma_m_max])
ax.legend(loc='upper right', framealpha=0.0)

plt.show()
