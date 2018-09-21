############################################
# Makes a flux map from observational data
# (which matches the kappa map style)
# and a cube of [kappa, gamma, flux]
# (works only for a smaller flux map as data
# is reduced to its dimensions)
############################################


import os
import sys
import numpy as np
from scipy import interpolate
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits


print("Make data cube...")

if len(sys.argv) > 1 and sys.argv[1] == "default":
    # Tested

    # Create a map with default value at all points
    print("... using default flux value.")

    hdul_kappa = fits.open('./kappa_map.fits')
    hdr_kappa = hdul_kappa[0].header
    data_kappa = hdul_kappa[0].data

    hdul_gamma = fits.open('./gamma_map.fits')
    data_gamma = hdul_gamma[0].data

    one = np.ones(data_kappa.shape)
    default_flux = 1.*10**(-6)*u.erg/(u.cm**2*u.s)

    new_hdr = hdr_kappa
    new_hdr['COMMENT'] = 'Cube of data containing kappa, gamma and flux maps (in this order) for target MACSJ1149.5+2223, made by Emma pipeline for computing transient detection statistics.'
    new_hdr['COMMENT'] = 'This header corresponds to the kappa and gamma maps (computed for an infinite redshift of the source).'
    new_hdr['COMMENT'] = 'Complementary data for the flux: - units: of erg/s/cm**2  - origin: default value mode.'

    hdu_cube = fits.PrimaryHDU([data_kappa, data_gamma, default_flux.value*one])
    hdul_cube = fits.HDUList([hdu_cube])
    hdul_cube.writeto('./data_cube.fits')
else:
    # Tested (possibly add check that corners are inside)

    # Create the map using observational data
    print("... using observational flux data.")

    hdul_kappa = fits.open('./kappa_map.fits')

    # TO TEST
    if not(os.path.isfile('./obs_map.fits')):
        print('Error: No observation map.')

    hdul_obs = fits.open('./obs_map.fits')

    hdr_kappa = hdul_kappa[0].header
    hdr_obs = hdul_obs[0].header  # Units are deg (see CUNIT in header)

    data_kappa = hdul_kappa[0].data
    # print(data_kappa.shape) gives (1000, 1000)
    data_obs = hdul_obs[0].data
    # print(data_obs.shape) gives (3, 6006, 6006)

    hdul_gamma = fits.open('./gamma_map.fits')
    data_gamma = hdul_gamma[0].data

    # Take only red flux into account as it traces cluster galaxies best
    # and get rid of unphysical negative values.
    data_obs_R = (data_obs[0, :, :]+np.absolute(data_obs[0, :, :]))/2
    if np.amin(data_obs_R) >= 0:
        print("Successfully removed unphysical negative values of electron count.")

    wcs_kappa = wcs.WCS(hdr_kappa)
    wcs_obs = wcs.WCS(hdr_obs)

    x_kappa, y_kappa = data_kappa.shape
    color, x_obs, y_obs = data_obs.shape

    print("Matching the corners.")
    wcs_obs_R = wcs.WCS(hdr_obs, naxis=[1, 2])
    # Maybe add a check that corners are within kappa map ...
    corners = wcs_obs_R.calc_footprint()

    catalog_along_x = SkyCoord.from_pixel(np.arange(x_kappa), 0, wcs_kappa)
    catalog_along_y = SkyCoord.from_pixel(0, np.arange(y_kappa), wcs_kappa)

    kx = []
    ky = []

    for i in range(len(corners)):
        obs_ra, obs_dec = corners[i]
        obs_coord = SkyCoord(ra=obs_ra, dec=obs_dec, unit=(u.degree, u.degree))
        kx.append(int(obs_coord.match_to_catalog_sky(catalog_along_x)[0]))
        ky.append(int(obs_coord.match_to_catalog_sky(catalog_along_y)[0]))

    # Check with Christoph:
    print("Interpolation starts.")
    if kx[0] == kx[1] and ky[0] == ky[3]:
        x = np.arange(x_obs)
        y = np.arange(y_obs)
        f = interpolate.interp2d(x, y, data_obs_R, kind='linear')

        space_x = kx[3]-kx[0]+1
        space_y = ky[1]-ky[0]+1
        x_new = np.round(np.linspace(0, x_obs-1, num=space_x)).astype(int)
        y_new = np.round(np.linspace(0, y_obs-1, num=space_y)).astype(int)
        z_new = f(x_new, y_new)  # Check this !!!

        new_kappa = np.zeros((space_x, space_y))
        new_gamma = np.zeros((space_x, space_y))
        new_flux = z_new

        for i in range(space_x):
            for j in range(space_y):
                new_kappa[i, j] = data_kappa[kx[0]+i, ky[0]+j]
                new_gamma[i, j] = data_gamma[kx[0]+i, ky[0]+j]

        hst_factor = hdr_obs['PHOTFLAM']
        band_width = 2511  # in angstrom, from ACS camera documentation.
        new_flux = hst_factor*band_width*new_flux

        new_hdr = hdr_kappa
        new_hdr['NAXIS1'], new_hdr['NAXIS2'] = new_kappa.shape
        new_hdr['CRPIX1'] = new_hdr['CRPIX1'] - kx[0]
        new_hdr['CRPIX2'] = new_hdr['CRPIX2'] - ky[0]
        new_hdr['COMMENT'] = 'Cube of data containing kappa, gamma and flux maps (in this order) for target MACSJ1149.5+2223, made by Emma pipeline for computing transient detection statistics.'
        new_hdr['COMMENT'] = 'This header corresponds to the kappa and gamma maps (computed for an infinite redshift of the source).'
        new_hdr['COMMENT'] = 'Complementary data for the flux: - units: of erg/s/cm**2  - origin: Remy data (i.e. HST image stripped of blue background) in band F814W.'

        hdu_cube = fits.PrimaryHDU([new_kappa, new_gamma, new_flux], header=new_hdr)
        hdul_cube = fits.HDUList([hdu_cube])
        hdul_cube.writeto('./data_cube.fits')
    else:
        print("Error: Observational image should be straight.")

    hdul_kappa.close()
    hdul_gamma.close()
    hdul_obs.close()
