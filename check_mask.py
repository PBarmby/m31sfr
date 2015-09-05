import numpy as np
import astropy.io.fits as fits
import astropy.units as u
import glob
from gal_radii_pb import correct_rgc
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle, Distance

# hard code some filenames here to save typing
hi_masks = ['masks/mask_hi_25.fits', 'masks/mask_hi_18.fits', 'masks/mask_hi_08.fits']
co_masks = ['masks/mask_co_25.fits', 'masks/mask_co_18.fits', 'masks/mask_co_08.fits']
hi750_masks = ['masks/mask_750_hi_25.fits', 'masks/mask_750_hi_18.fits', 'masks/mask_750_hi_08.fits']
sr_mask_lims= [(18,25),(8,18),(0,8)]

def do_check(masklist = hi_masks, masklims = sr_mask_lims, save_new_mask=False):
    '''check that masks which are supposed to select only certain galactocentric distances
       do this correctly
       input: 
       masklist: list of masks to check
       masklims: list of tuples with distance limits corresponding to each mask
       save_new_mask: write out new masks?'''
    if len(masklist) != len(masklims):
        print 'first 2 arguments must have same length'
        return

    for i, mask in enumerate(masklist):
        lower, upper = masklims[i]

        # read in mask to be checked, set its NaNs to 0
        maskdat, maskhead = fits.getdata(mask, header=True)
        maskdat[np.isnan(maskdat)]=0

        # read mask WCS and construct array of coordinates
        maskwcs = WCS(maskhead)
        pix_indices = np.indices(maskdat.shape)
        coords = SkyCoord.from_pixel(pix_indices[1],pix_indices[0],maskwcs,origin=0)
        # compute galactocentric distance using same parameters for galaxy as Sahar did
        galacto_dist = correct_rgc(coords, glx_ctr=SkyCoord(10.724382145*u.deg, 41.332205842*u.deg), 
                                           glx_PA=Angle('38d'), 
                                           glx_dist=Distance(780, unit=u.kpc),
                                           glx_incl=Angle('75.883d'))
        galacto_dist_kpc = galacto_dist.to(u.kpc).value

        # construct a new mask with specified distances and check to see if it's same as original
        newmask = np.zeros(maskdat.shape)
        newmask[np.logical_and(galacto_dist_kpc > lower, galacto_dist_kpc< upper)]=1.0
        check = np.array_equal(maskdat, newmask)
        print '%s %s' % (mask, check)
        if not check:
            matches = np.count_nonzero(maskdat==newmask)
            print '%d matches out of %d pix (%f)' % (matches, maskdat.size, float(matches)/maskdat.size)
        
        if save_new_mask: # write out the new mask (with the old header) to a fits file
            newname = mask[:-5] + '.check.fits'
            fits.writeto(newname, newmask, maskhead)
    return
