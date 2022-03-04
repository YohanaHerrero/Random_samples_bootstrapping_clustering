from __future__ import division
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import matplotlib.pyplot as plt
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import fits
import math
from matplotlib.ticker import ScalarFormatter
import glob
import itertools
from matplotlib.lines import Line2D
from sklearn.utils import resample

#I read the fits file to obtain the data
event_filename = get_pkg_data_filename('60fields.fits')
events = Table.read(event_filename, hdu=1)
hdul = fits.open(event_filename)
data = hdul[1].data  

#extract the colums from the table
dec = data['DEC']
ra = data['RA']
redshift = data['Z']

#some specific selection
dec_sel = dec[:1593]
ra_sel = ra[:1593]
redshift_sel = redshift[:1593]

#redshift selection from 3 to 6
select = (redshift_sel >= 3 ) & (redshift_sel <= 6.)
zf = redshift_sel[select]
decf = dec_sel[select] 
raf = ra_sel[select] 


#I create 100 pseudo-samples that slightly differ from the real data set and save them in txt files
#I perturb the real data set by re-sampling with the bootstraping technique
file_name = "bootstrap_data68fields_3.3z6.{i}.txt"
ra_dec_z = np.vstack([raf, decf, zf]) 
i=0
while i<100:
    ra_dec_z_trans = ra_dec_z.T  
    gal_bootstrapped = resample(ra_dec_z_trans)
    new_ra_dec_z = gal_bootstrapped.T  
    ra_new=new_ra_dec_z[0]
    dec_new=new_ra_dec_z[1]
    z_new=new_ra_dec_z[2]
    header = " RA, DEC, z"
    np.savetxt(file_name.format(i=i), np.array([ra_new, dec_new, z_new]).T, delimiter='\t', header=header, fmt="%s")
    i = i+1

#read the created 100 files and compute the clustering in each of the 100 pseudo-samples 

#clustering statistic = K-estimator (see Adelberger et al. 2005)
ax = plt.figure().add_subplot(111) 
files = glob.glob('*.txt')
kab_full_array = np.array([])
for f in files:
    #read data
    data = np.loadtxt(f)
    RAf=data[:,0]
    DECf=data[:,1]
    Zf=data[:,2]
    
    #calculate comoving separations
    zij = np.array([]) 
    co = np.array([])
    cos = np.array([])  
    for k, zk in enumerate(Zf):
        for l, zl in enumerate(Zf[k+1:]):
            co = np.append(co, cosmo.comoving_distance(zk).value)
            cos = np. append(cos, cosmo.comoving_distance(zl).value)
    zij = np.append(zij, abs(co-cos))*0.7 
    
    #calculate transverse separations
    phi = np.array([])
    d = np.array([])
    r = np.array([])
    for i, deci in enumerate(DECf):
        d = np.append(d, deci-DECf[i+1:])
        for j, rai in enumerate(RAf[i+1:]):
            r = np.append(r, RAf[i]-rai)
    phi = np.append(phi, np.sqrt((r)**2+(d)**2)*math.pi/180)
    rij = np.array([])
    pairs = itertools.combinations(Zf, 2)
    rzav = np.array([])
    for pair in pairs:
        rzav= np.append(rzav, cosmo.comoving_distance((pair[0] + pair[1])/2).value)
    rij = np.append(rij, rzav * phi )*0.7 
    
    #calculate the K-estimator
    kab = np.array([])
    bins=np.array([0.045,0.19,0.375,0.605,0.89,1.25,1.7,2.26,2.95,3.79,4.85,6.15, 7.76,9.8,12.25])
    err = np.array([])
    binp = np.array([])    
    for k, bini in enumerate(bins):
        if k < len(bins)-1:                                                         
            idxtrans = (rij >= bini) & (rij < (bini+bins[k+1]))                                    
            idxlos1 = (zij > 0) & (zij < 7)                                                    
            idxlos2 = (zij > 0) & (zij < 45)                                                
            kab = np.append(kab, sum(idxtrans & idxlos1)/sum(idxtrans & idxlos2))     
            binp = np.append(binp, bini + (bins[k+1]-bini)/2)
    
    plt.plot(binp, kab, color='lightgray', alpha=0.9, linewidth=0.5, zorder=1) 
    kab_full_array = np.append(kab_full_array,kab) 
    
kab_splitted = np.split(kab_full_array, 100) 
k_mean = np.mean(kab_splitted, axis=0) 
err_k_mean = np.std(kab_splitted, axis=0) #standard deviation of the 100 kab = error from bootstrapping approach

ax.scatter(binp, k_mean, s=10, c = 'red', marker='o')
ax.errorbar(binp, k_mean, yerr=err_k_mean, xerr=None, c = 'red', ls='None', barsabove=True, capsize=2, elinewidth=1)

#we calculate the K estimator for the real sample
zij = np.array([]) 
co = np.array([])
cos = np.array([])  
for k, zk in enumerate(zf):
    for l, zl in enumerate(zf[k+1:]):
        co = np.append(co, cosmo.comoving_distance(zk).value)
        cos = np. append(cos, cosmo.comoving_distance(zl).value)
zij = np.append(zij, abs(co-cos))*0.7 

phi = np.array([])
d = np.array([])
r = np.array([])
for i, deci in enumerate(decf):
    d = np.append(d, deci-decf[i+1:])
    for j, rai in enumerate(raf[i+1:]):
        r = np.append(r, raf[i]-rai)
phi = np.append(phi, np.sqrt((r)**2+(d)**2)*math.pi/180) 

rij = np.array([])
pairs = itertools.combinations(zf, 2)
rzav = np.array([])
for pair in pairs:
    rzav= np.append(rzav, cosmo.comoving_distance((pair[0] + pair[1])/2).value)
rij = np.append(rij, rzav * phi )*0.7 

kab = np.array([])
bins=np.array([0.045,0.19,0.375,0.605,0.89,1.25,1.7,2.26,2.95,3.79,4.85,6.15, 7.76,9.8,12.25])
err = np.array([])
binp = np.array([])
for k, bini in enumerate(bins):
    if k < len(bins)-1:                                                         
        idxtrans = (rij >= bini) & (rij < (bini+bins[k+1]))                                    
        idxlos1 = (zij > 0) & (zij < 7)                                                    
        idxlos2 = (zij > 0) & (zij < 45)   
        err = np.append(err, math.sqrt(sum(idxtrans & idxlos1))/sum(idxtrans & idxlos2)) #Poisson error                                                  
        kab = np.append(kab, sum(idxtrans & idxlos1)/sum(idxtrans & idxlos2))     
        binp = np.append(binp, bini + (bins[k+1]-bini)/2)
   
ax.scatter(binp, kab, s=10, c = 'blue', marker='o')
ax.errorbar(binp, kab, yerr=err, xerr=None, c = 'blue', ls='None', barsabove=True, capsize=2, elinewidth=1)

horiz_line = np.array([0.2 for m in range(len(kab))])
ax.plot(binp, horiz_line, 'k-', linewidth = 1)
plt.xlabel(r'$R_{ij}$ [$h^{-1}$Mpc]', fontsize=14)
plt.ylabel(r'$K^{0,7}_{7,35}$', fontsize=14)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_tick_params(direction='in', which='both')
ax.yaxis.set_tick_params(direction='in', which='both')
plt.tick_params(labelsize = 'large')
blue_dot = Line2D([0], [0], marker='o', color='blue', linestyle = 'None', label='real sample (with Poisson)')
red_dot = Line2D([0], [0], marker='o', color='red', linestyle = 'None', label='mean boots (with bootstrapping)')
lightgray_line = Line2D([0], [0], color='lightgray', label='boots')
plt.legend(handles=[blue_dot, red_dot,lightgray_line], loc = 'upper right').get_frame().set_edgecolor('black')
ax.set_xscale('log') 
for axis in [ax.xaxis, ax.yaxis]:
     axis.set_major_formatter(ScalarFormatter())
plt.tight_layout()                                                                      
plt.grid(False) 
#plt.savefig('K-estimator in 100 boostrapping samples', dpi=500)
plt.show()
