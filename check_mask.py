import numpy as np
import astropy.io.fits as fits
import glob
from gal_radii_pb import correct_rgc
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

sr_masks = ['masks/mask_hi_25.fits', 'masks/mask_hi_18.fits', 'masks/mask_hi_08.fits']
sr_mask_lims= [(18,25),(8,18),(0,8)]

def do_check(masklist = sr_masks, masklims = sr_mask_lims):
    for i, mask in enumerate(masklist):
        lower, upper = masklims[i]

        # read in mask to be checked, set its NaNs to 0
        maskdat, maskhead = fits.getdata(mask, header=True)
        maskdat[np.isnan(maskdat)]=0

        # read mask WCS and construct array of coordinates, compute distance
        maskwcs = WCS(maskhead)
        pix_indices = np.indices(maskdat.shape)
        coords = SkyCoord.from_pixel(pix_indices[1],pix_indices[0],maskwcs,origin=0)
        galacto_dist_kpc = correct_rgc(coords).value # using fact that defaults of this function are for M31!

        # construct a new mask with specified distances and check to see if it's same as original
        newmask = np.zeros(maskdat.shape)
        newmask[np.logical_and(galacto_dist_kpc > lower, galacto_dist_kpc< upper)]=1.0
        check = np.array_equal(maskdat, newmask)
        print '%s %s' % (mask, check)
        if not check:
            pix_diff = 
        
        # write out the new mask (with the old header) to a fits file
        newname = mask[:-5] + '.check.fits'
        fits.writeto(newname, newmask, maskhead)
    return
