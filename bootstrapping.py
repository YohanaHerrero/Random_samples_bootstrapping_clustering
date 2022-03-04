from __future__ import division
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import fits
import math
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
import glob
import itertools
from matplotlib.lines import Line2D
from scipy.sparse import coo_matrix
from sklearn.utils import resample


#CALCULATION OF ERRORS FOLLOWING DURKALEK2015
######################## BOOTSTRAP + RESAMPLING METHOD ###########################


#load and open file the 60 fields catalogue
event_filename1 = get_pkg_data_filename('60fields.fits')
events1 = Table.read(event_filename1, hdu=1)
hdul1 = fits.open(event_filename1)
data1 = hdul1[1].data  


#extract the colums from the table
d1 = data1['DEC']
r1 = data1['RA']
redshift1 = data1['Z']
id1 = data1['UNIQUE_ID']
#ID selection, until ID row=1602, all the ID start with 1, after they start with 2
decsel1 = d1[:1593]
rasel1 = r1[:1593]
redshiftsel1 = redshift1[:1593]
id1 = id1[:1593]

#Right Ascension correction (differentials)
RA_offset1 = (rasel1 - np.mean(rasel1))*np.cos(decsel1*math.pi/180)
RA1 = np.mean(rasel1) + RA_offset1


#redshift selection from 3 to 6
select1 = (redshiftsel1 >= 3.3 ) & (redshiftsel1 <= 6.)
z1 = redshiftsel1[select1]
dec1 = decsel1[select1] 
ra1 = RA1[select1] 


#load and open file the merged catalogue
event_filename2 = get_pkg_data_filename('merged_catalog_e40_v0.9.fits')
events2 = Table.read(event_filename2, hdu=1)
hdul2 = fits.open(event_filename2)
data2 = hdul2[1].data  


#extract the colums from the table
d2 = data2['DEC']
r2 = data2['RA']
redshift2 = data2['REDSHIFT']
id2 = data2['ID']

#ID selection, until ID row=700, all the ID start with 1, after they start with 2
decsel2 = d2[:701]
rasel2 = r2[:701]
redshiftsel2 = redshift2[:701]

#Right Ascension correction (differentials)
RA_offset2 = (rasel2 - np.mean(rasel2))*np.cos(decsel2*math.pi/180)
RA2 = np.mean(rasel2) + RA_offset2


#redshift selection from 3 to 6
select2 = (redshiftsel2 >= 3.3 ) & (redshiftsel2 <= 6.)
z2 = redshiftsel2[select2]
dec2 = decsel2[select2] 
ra2 = RA2[select2] 


#Galaxies with id3
declsel2 = d2[956:]
rassel2 = r2[956:]
redshifttsel2 = redshift2[956:]

#Right Ascension correction (differentials)
RA_offset2 = (rassel2 - np.mean(rassel2))*np.cos(declsel2*math.pi/180)
RAS2 = np.mean(rassel2) + RA_offset2
#redshift selection from 3 to 6 for id3
select2 = (redshifttsel2 >= 3.3 ) & (redshifttsel2 <= 6.)
zs2 = redshifttsel2[select2]
decl2 = declsel2[select2] 
ras2 = RAS2[select2] 

raf = np.hstack((ra1,ra2,ras2))
decf = np.hstack((dec1,dec2,decl2))
zf = np.hstack((z1,z2,zs2))


'''
#plot 68 fields
cm = plt.cm.get_cmap('jet') 
fig = plt.figure().add_subplot(111)
plt.scatter(raf,decf, s=10, c=zf, marker='o', cmap=cm)
plt.gca().invert_xaxis()
plt.colorbar().set_label('z')
plt.clim(3.3, 6)  
fig.xaxis.set_ticks_position('both')
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_tick_params(direction='in', which='both')
fig.yaxis.set_tick_params(direction='in', which='both')
plt.xlabel("RA", fontsize=14)
plt.ylabel("DEC", fontsize=14)
plt.tight_layout()
plt.show()
'''

print('number of real gal', len(raf))
'''
#First: I create 100 pseudo-samples that slightly differ from the real data set (I should perturb the real data set by re-sampling with bootstrap)
#we create 100 pseudo-samples txt files for the real data in which we will compute the K estimator

file_name = "bootstrap_data68fields_3.3z6.{i}.txt"
ra_dec_z = np.vstack([raf, decf, zf]) #for the bootstrap of each gal (ra,dec,z of each gal at a time)
i=0
while i<100:
    #bootstrap of each gal (ra,dec,z of each gal at a time)
    ra_dec_z_trans = ra_dec_z.T #now each row has the ra, dec and z of each gal. each row is a different gal with the 3 values
    gal_bootstrapped = resample(ra_dec_z_trans)
    new_ra_dec_z = gal_bootstrapped.T #now each column is one galaxy, so each row is again ra, dec and z (3 rows in total)
    ra_new=new_ra_dec_z[0]
    dec_new=new_ra_dec_z[1]
    z_new=new_ra_dec_z[2]
    
    #for the other methods of bootstrapping
    #ra_new = resample(raf, replace=True, n_samples=len(raf)) #n_samples is the lenght that I want for the new array, random_state=0 or 1 (an int) gives us the same DD_new every time we run the code
    #not_taken = [x for x in raf if x not in ra_new] #the element/s of ra/dec/z that were not selected for the creation of wp_new this time
    #dec_new = resample(decf, replace=True, n_samples=len(decf))
    #z_new = resample(zf, replace=True, n_samples=len(zf))
    header = " RA, DEC, z"
    np.savetxt(file_name.format(i=i), np.array([ra_new, dec_new, z_new]).T, delimiter='\t', header=header, fmt="%s")
    i = i+1

#we read the created 100 files and we compute the correlation function for each of the 100 pseudo-samples 

#we plot ra-dec of the 100 bootstrap samples
files = glob.glob('*.txt')
cm = plt.cm.get_cmap('jet') 
fig = plt.figure().add_subplot(111)

for f in files:
    data = np.loadtxt(f)
    ra=data[:,0]
    dec=data[:,1]
    z=data[:,2]
    plt.scatter(ra,dec, s=10, c=z, marker='o', cmap=cm)
    
plt.clim(3.3, 6) 
plt.gca().invert_xaxis()
plt.colorbar().set_label('z')
fig.xaxis.set_ticks_position('both')
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_tick_params(direction='in', which='both')
fig.yaxis.set_tick_params(direction='in', which='both')
plt.title("100_bootstrap_samples")
plt.xlabel("RA", fontsize=14)
plt.ylabel("DEC", fontsize=14)
plt.tight_layout()
plt.savefig('RA-DEC_100_bootstrap_samples',dpi=500)
#plt.show()
'''

#we calculate the K estimator for the 100 bootstrap samples
ax = plt.figure().add_subplot(111) #only define it once, otherwise I will not plot all curves together
files = glob.glob('*.txt')
kab_full_array = np.array([])
k_1bin=np.array([])
err_1bin=np.array([])
for f in files:
    data = np.loadtxt(f)
    RAf=data[:,0]
    DECf=data[:,1]
    Zf=data[:,2]
    zij = np.array([]) 
    co = np.array([])
    cos = np.array([])  
    for k, zk in enumerate(Zf):
        for l, zl in enumerate(Zf[k+1:]):
            co = np.append(co, cosmo.comoving_distance(zk).value)
            cos = np. append(cos, cosmo.comoving_distance(zl).value)
    zij = np.append(zij, abs(co-cos))*0.7 #radial distance #the unit will be h^(-1)Mpc
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
    rij = np.append(rij, rzav * phi )*0.7 #transverse distance #the unit will be h^(-1)Mpc
    
    idx=(rij<3.5)
    idxlos1 = (zij > 0) & (zij < 7)
    idxlos2 = (zij > 0) & (zij < 35) 
    k_1bin = np.append(k_1bin, sum(idx & idxlos1)/sum(idx & idxlos2))
    err_1bin = np.append(err_1bin, math.sqrt(sum(idx & idxlos1))/sum(idx & idxlos2))
    
    kab = np.array([])
    bins=np.array([0.045,0.19,0.375,0.605,0.89,1.25,1.7,2.26,2.95,3.79,4.85,6.15, 7.76,9.8,12.25])
    err = np.array([])
    binp = np.array([])
    
    for k, bini in enumerate(bins):
        idxtrans = (rij >= bini) & (rij < (bini+1.125))
        idxlos1 = (zij > 0) & (zij < 7)
        idxlos2 = (zij > 0) & (zij < 35) 
            
        kab = np.append(kab, sum(idxtrans & idxlos1)/sum(idxtrans & idxlos2))
        binp = np.append(binp, bini + 1.125/2)
    
    plt.plot(binp, kab, color='lightgray', alpha=0.9, linewidth=0.5, zorder=1) #I plot all the mocks without error bars #zorder=1 makes the gray lines to stay behind the blue and red curves
    kab_full_array = np.append(kab_full_array,kab) #creo un array con todos los kab juntos, no una matriz #hstack, matrix, concatenate would do the same as np.append
    
kab_splitted = np.split(kab_full_array, 100) #we created a matrix where each row is one of the kab calculated arrays
k_mean = np.mean(kab_splitted, axis=0) #we calculate the average of each column of kab_splitted (average of the first element of each kab, average of the second element of each kab, etc)
err_k_mean = np.std(kab_splitted, axis=0) #standard deviation of the 100 kab
print(kab_splitted)
print('k mean', k_mean)
print('err mean', err_k_mean) 

#Print out the K and err measured in just one bin (rij<3.5)
k_1bin_mean = np.mean(k_1bin)
err_1bin_mean = np.mean(err_1bin)
print('K measured in rij<3.5 Mpc', k_1bin)
print('err of K measured in rij<3.5 Mpc', err_1bin)
print('Mean K measured in rij<3.5 Mpc', k_1bin_mean)
print('Mean err of K measured in rij<3.5 Mpc', err_1bin_mean)
   
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
zij = np.append(zij, abs(co-cos))*0.7 #radial distance #the unit will be h^(-1)Mpc


phi = np.array([])
d = np.array([])
r = np.array([])
for i, deci in enumerate(decf):
    d = np.append(d, deci-decf[i+1:])
    for j, rai in enumerate(raf[i+1:]):
        r = np.append(r, raf[i]-rai)
phi = np.append(phi, np.sqrt((r)**2+(d)**2)*math.pi/180) #approx. of dtheta = arccos{ sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2) }
 

rij = np.array([])
pairs = itertools.combinations(zf, 2)
rzav = np.array([])

for pair in pairs:
    rzav= np.append(rzav, cosmo.comoving_distance((pair[0] + pair[1])/2).value)
rij = np.append(rij, rzav * phi )*0.7 #transverse distance #the unit will be h^(-1)Mpc

kab = np.array([])
bins=np.array([0.045,0.19,0.375,0.605,0.89,1.25,1.7,2.26,2.95,3.79,4.85,6.15, 7.76,9.8,12.25])
err = np.array([])
binp = np.array([])
for k, bini in enumerate(bins):
     idxtrans = (rij >= bini) & (rij < (bini+1.125))
     idxlos1 = (zij > 0) & (zij < 7)
     idxlos2 = (zij > 0) & (zij < 35) 
            
     kab = np.append(kab, sum(idxtrans & idxlos1)/sum(idxtrans & idxlos2))
     err = np.append(err, math.sqrt(sum(idxtrans & idxlos1))/sum(idxtrans & idxlos2))
     binp = np.append(binp, bini + 1.125/2)
   
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
blue_dot = Line2D([0], [0], marker='o', color='blue', linestyle = 'None', label='MUSE Wide')
red_dot = Line2D([0], [0], marker='o', color='red', linestyle = 'None', label='mean mocks')
lightgray_line = Line2D([0], [0], color='lightgray', label='mocks')
plt.legend(handles=[blue_dot, red_dot,lightgray_line], loc = 'upper right').get_frame().set_edgecolor('black')
ax.set_xscale('log') 
for axis in [ax.xaxis, ax.yaxis]:
     axis.set_major_formatter(ScalarFormatter())
plt.tight_layout()                                                                      
plt.grid(False) 
plt.savefig('K-estimator68fields_h_7-35_100_mocks_ra,dec,z_boots_at_a_time.png', dpi=500)
plt.show()

'''
#K7:40 fitting for gamma=1.9
r0val = np.arange(0.2, 2., 0.1)
xi_ij_a = np.array([])
xi_ij_b = np.array([])
err1 = np.array([])
k1 = np.array([])
integ1 = np.array([])
integ2 = np.array([])

idx = (rij > .2) & (rij < 2.)
transin = rij[idx]
for j, r0valj in enumerate(r0val):
    for h, transinh in enumerate(transin):
        integ1 = integrate.quad(lambda Z: (math.sqrt(transinh**2 + Z**2) / r0valj)**(-1.9),7,40)[0]
        integ2 = integrate.quad(lambda Z: (math.sqrt(transinh**2 + Z**2) / r0valj)**(-1.9),0,7)[0] 
        xi_ij_a = np.append(xi_ij_a, (1/(40-7)) * integ1)
        xi_ij_b = np.append(xi_ij_b, (1/(7-0)) * integ2)
                        
    k1 = np.append(k1, ((7-0)*sum(1+xi_ij_b[j*len(transin):]))/((7-0)*sum(1+xi_ij_b[j*len(transin):])+(40-7)*sum(1+xi_ij_a[j*len(transin):])))
print('k1: ', k1)

kmeasured = sum(idx & idxlos1)/sum(idx & idxlos2) #I count the pairs in a single bin rij < 5
err1 = np.append(err1, math.sqrt(sum(idx & idxlos1))/sum(idx & idxlos2))     
difference = np.abs(k1-kmeasured) #Difference between the actual K and the expectation value
differencelow = np.abs(k1-(kmeasured-err1))
differenceup = np.abs(k1-(kmeasured+err1))
min_dif = np.argmin(difference) #find the min element position in the array difference
r0opt = r0val[min_dif]
min_dif_up = np.argmin(differenceup)
r0up = r0val[min_dif_up]
min_dif_low = np.argmin(differencelow)
r0low = r0val[min_dif_low]

#from the average of the 100 pseudo samples and its deviation, we calculate .........not really, covariance matrix firs.......... the difference to our real correlation function curve. This will be the error bar estimation that we will use


#plot a histogram for the r0 values
'''
