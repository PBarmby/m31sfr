import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
import glob

sr_masks = ['masks/mask_hi_25.fits', 'masks/mask_hi_18.fits', 'masks/mask_hi_08.fits']
# mean values from eg ./ks/hi/tir/lt8/IndividualFit.Rout -- medians are almost the same
slopes=[0.4879,1.9492,0.09117] 
ints=[-9.9556,-10.0130,-9.07312]

def doplot(sfr='sfr_fir_reg_hi.fits',gas='hi_per_mass_sun_pc2.fits',masklist=sr_masks):
    fig, ax = plt.subplots()
    ax.set_ylabel('SFR(FIR)')
    ax.set_xlabel('Sigma_HI')
    ax.set_xlim(-1.15,1.15)
    ax.set_ylim(-12,-7)
    xgas=np.arange(-1,1,0.1)
    sfrdat = fits.getdata(sfr)
    gasdat = fits.getdata(gas)
    for i,mask in enumerate(masklist):
        maskdat = fits.getdata(mask)
        ax.plot(np.log10(gasdat[~np.isnan(maskdat)]),np.log10(sfrdat[~np.isnan(maskdat)]),marker='.',ms=1,label= mask)
        ax.plot(xgas, xgas*slopes[i]+ints[i], ls='solid',marker='None',color='k')
    ax.legend()
    fig.show()
    return
