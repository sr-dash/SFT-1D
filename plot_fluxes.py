import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
import matplotlib.style
#plt.ion()
## Plotting canvas properties.
params = {'legend.fontsize': 12,
          'axes.labelsize': 12,
          'axes.titlesize': 12,
          'xtick.labelsize' :10,
          'ytick.labelsize': 10,
          'grid.color': 'k',
          'grid.linestyle': ':',
          'grid.linewidth': 0.5,
          'mathtext.fontset' : 'stix',
          'mathtext.rm'      : 'DejaVu serif',
          'font.family'      : 'DejaVu serif',
          'font.serif'       : "Times New Roman", # or "Times"          
         }
matplotlib.rcParams.update(params)

from scipy.io import netcdf_file
import os
import pickle
import glob
import datetime
from datetime import datetime as dt1
from datetime import timedelta
import f90nml

# Read the parameter file to retrive the values
nml = f90nml.read(os.getcwd()+'/initial_Parameters.nml')
eta = nml['user']['eta']
peak_lat = nml['user']['peak_lat']
out_freq = nml['user']['output_freq']

# Set the time range for the plot.
# datetime(2010,6,17)+timedelta(days=5987)
t_sft = np.arange(dt1(2010,6,17), dt1(2026,11,7), timedelta(days=1)).astype(dt1)


# Function to convert time string to fractional year.
def frac_year(time_string):
	SECONDS_IN_DAY = 60.0*60.0*24.0
	SECONDS_IN_YEAR = SECONDS_IN_DAY*365.25
	t1 = datetime.datetime.strptime(time_string, '%Y-%m-%d') #stime.parse_time(time_string)
	t1_diff = t1 - datetime.datetime(t1.year,1,1,0,0,0)
	frac_year = (t1_diff.days + t1_diff.seconds/(60*60*24.0))/365.25
	return t1.year + frac_year

# Create a plot dir to dave the results.
PLOTPATH = os.getcwd()+'/plots'
if not os.path.exists(PLOTPATH):
    os.makedirs(PLOTPATH)

# Read the butterfly diagram file from the output directory.
bfly1 = glob.glob(os.getcwd()+'/output_files/bfly_%3d_*.nc'%int(eta))[0]

bfly_eta = glob.glob(os.getcwd()+'/output_files/bfly_advfluxes_%3d_*.nc'%int(eta))[0]
bfly_adv = glob.glob(os.getcwd()+'/output_files/bfly_resfluxes_%3d_*.nc'%int(eta))[0]

fh2 = netcdf_file(bfly1)
bfly = fh2.variables['bfly'].data.copy()
sth = fh2.variables['lat'].data.copy()
time = fh2.variables['time'].data.copy()
fh2.close()

fh3 = netcdf_file(bfly_eta)
bfly_eta = fh3.variables['bfly'].data.copy()
sth_eta = fh3.variables['lat'].data.copy()
time_eta = fh3.variables['time'].data.copy()
fh3.close()

fh4 = netcdf_file(bfly_adv)
bfly_adv = fh4.variables['bfly'].data.copy()
sth_adv = fh4.variables['lat'].data.copy()
time_adv = fh4.variables['time'].data.copy()
fh4.close()

# Read HMI butterfly diagram data from the file.
picklefile_hmi_bfly = open('./hmi_bfly_new.p','rb')
bflyhmi = pickle.load(picklefile_hmi_bfly)
timearr = pickle.load(picklefile_hmi_bfly)
latbfly = pickle.load(picklefile_hmi_bfly)
picklefile_hmi_bfly.close()

# Plot the butterfly diagram for the SFT simulation output
fig = plt.figure(figsize=[12,12])
ax1 = plt.subplot(411)

vel = glob.glob(os.getcwd()+'/output_files/MC_vel*.dat')[0]
v1 = np.loadtxt(vel)
L1 = 6.96e5
# print('%2.1f'%(np.max(v1[:,1])*L1*1E3))
pm = ax1.pcolormesh(time,np.rad2deg(np.arcsin(sth)),bfly,cmap='bwr',vmax=10,vmin=-10)
ax1.set_xlim([time[0],time[-1]])  
# ax1.set_xlim([-90,90])
ax1.set_xlabel('Years')
ax1.set_ylabel('Latitude (degrees)')
ax1.axvline(x = frac_year('2025-11-08'), c='brown',ls='--',alpha=0.8)
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.15)
fig.colorbar(pm, cax=cax, orientation='vertical',label=r'B$_r$ [G]')
ax1.set_title(r'$\eta$ = %3d km$^2$/s, V0 = %2.1f m/s'%(int(eta),np.max(v1[:,1])*L1*1E3))
ax1.text(0.01,0.9,r'SFT B$_r$ butterfly diagram',transform=ax1.transAxes,
        fontsize=10,color='brown')


ax2 = plt.subplot(412)
pm2 = ax2.pcolormesh(time,np.rad2deg(np.arcsin(sth)),bfly_eta,cmap='bwr',vmax=1e-7,vmin=-1e-7)
ax2.set_xlim([time[0],time[-1]])  
# ax1.set_xlim([-90,90])
ax2.set_xlabel('Years')
ax2.set_ylabel('Latitude (degrees)')
ax2.axvline(x = frac_year('2025-11-08'), c='brown',ls='--',alpha=0.8)
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.15)
fig.colorbar(pm2, cax=cax, orientation='vertical',label=r'F$_{resistive}$ [G/s]')
# ax2.set_title(r'$\eta$ = %3d km$^2$/s, V0 = %2.1f m/s'%(int(eta),np.max(v1[:,1])*L1*1E3))
ax2.text(0.01,0.9,r'Resistive flux butterfly diagram',transform=ax2.transAxes,
        fontsize=10,color='brown')

ax3 = plt.subplot(413)
# print('%2.1f'%(np.max(v1[:,1])*L1*1E3))
pm3 = ax3.pcolormesh(time,np.rad2deg(np.arcsin(sth)),bfly_adv,cmap='bwr',vmax=1e-7,vmin=-1e-7)
ax3.set_xlim([time[0],time[-1]])  
# ax1.set_xlim([-90,90])
ax3.set_xlabel('Years')
ax3.set_ylabel('Latitude (degrees)')
ax3.axvline(x = frac_year('2025-11-08'), c='brown',ls='--',alpha=0.8)
divider = make_axes_locatable(ax3)
cax = divider.append_axes('right', size='5%', pad=0.15)
fig.colorbar(pm3, cax=cax, orientation='vertical',label=r'F$_{advective}$ [G/s]')
# ax3.set_title(r'$\eta$ = %3d km$^2$/s, V0 = %2.1f m/s'%(int(eta),np.max(v1[:,1])*L1*1E3))
ax3.text(0.01,0.9,r'Advective flux butterfly diagram',transform=ax3.transAxes,
        fontsize=10,color='brown')


ax4 = plt.subplot(414)
im4 = ax4.pcolormesh(timearr,np.rad2deg(np.arcsin(latbfly)),bflyhmi,cmap='bwr',vmax=10,vmin=-10)
# ax2.set_xlim([timearr[0],timearr[-1]])
ax4.set_xlim([time[0],time[-1]]) 
ax4.set_xlabel('Years')
ax4.set_ylabel('Latitude (degrees)')
divider = make_axes_locatable(ax4)
cax = divider.append_axes('right', size='5%', pad=0.15)
fig.colorbar(im4, cax=cax, orientation='vertical',label='B$_r$ [G]')
ax4.text(0.01,0.9,r'HMI B$_r$ butterfly diagram',transform=ax4.transAxes,
        fontsize=10,color='brown')

plt.savefig(PLOTPATH+'/bfly_all_bipoles_fluxes%3d_%3d.png'%(int(eta),np.max(v1[:,1])*L1*1E4),
            dpi=300,transparent=False,bbox_inches='tight')
plt.show()