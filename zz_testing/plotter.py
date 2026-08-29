
"""
maybe pick a vertical level, plot the temperatures (as colored scatter points) on the map, make another plot of the climatology T at that same level, and another plot of the uncertainty, and another plot of (prof T - clim T) / uncertainty (sigma, not sigma^2).  for the last plot we should expect to see values between about -3 to +3 (+/- 3 standard deviations from the assumed uncertainty).
"""

import pdb
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import argparse
import cartopy.crs as ccrs

parser = argparse.ArgumentParser()
parser.add_argument("verticallevel", type=int)
args = parser.parse_args()

vertical_level = args.verticallevel

#vertical_level = 0

dot_size = 1
#dot_size = 5

# SUPER hacky, only works when run interactivley in ipython
input_file = [p.name for p in Path(input_directory).iterdir() if p.is_file()][0]

T_min = min(MITprof_ds['prof_T'].min().item(), MITprof_ds['prof_Tclim'].min().item())
T_max = max(MITprof_ds['prof_T'].max().item(), MITprof_ds['prof_Tclim'].max().item())

fig,axes = plt.subplots(2,2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,8))

for ax in axes.flat:
    ax.coastlines(color='black', linewidth=0.15)
    ax.patch.set_facecolor('#D9D9D9')
    ax.gridlines(draw_labels=True)
    ax.set_aspect('equal', anchor='C')
    #ax.set_adjustable("datalim")

prof_valid_indices = MITprof_ds['prof_T'][:,vertical_level].notnull()

a1 = axes[0,0].scatter(MITprof_ds['prof_lon'][prof_valid_indices], MITprof_ds['prof_lat'][prof_valid_indices], s=dot_size, c=MITprof_ds['prof_T'][:,vertical_level][prof_valid_indices], vmin=T_min, vmax=T_max)
#a1 = axes[0,0].scatter(MITprof_ds['prof_lon'], MITprof_ds['prof_lat'], s=dot_size, c=MITprof_ds['prof_T'][:,vertical_level], vmin=T_min, vmax=T_max)
fig.colorbar(a1, ax=axes[0,0])
axes[0,0].set_title("prof_T")

a2 = axes[0,1].scatter(MITprof_ds['prof_lon'][prof_valid_indices], MITprof_ds['prof_lat'][prof_valid_indices], s=dot_size, c=MITprof_ds['prof_Tclim'][:,vertical_level][prof_valid_indices], vmin=T_min, vmax=T_max)
#a2 = axes[0,1].scatter(MITprof_ds['prof_lon'], MITprof_ds['prof_lat'], s=dot_size, c=MITprof_ds['prof_Tclim'][:,vertical_level], vmin=T_min, vmax=T_max)
fig.colorbar(a2, ax=axes[0,1])
axes[0,1].set_title("prof_Tclim")

a3 = axes[1,0].scatter(MITprof_ds['prof_lon'][prof_valid_indices], MITprof_ds['prof_lat'][prof_valid_indices], s=dot_size, c=MITprof_ds['prof_Tuncertainty'][:,vertical_level][prof_valid_indices])
fig.colorbar(a3, ax=axes[1,0])
axes[1,0].set_title("prof_Tuncertainty")


values = (MITprof_ds['prof_T'][:,vertical_level] - MITprof_ds['prof_Tclim'][:,vertical_level])**2/MITprof_ds['prof_Tuncertainty'][:,vertical_level]**2
cbar_min = 0
cbar_max = 16
a4 = axes[1,1].scatter(MITprof_ds['prof_lon'], MITprof_ds['prof_lat'], s=dot_size, c=values, vmin=cbar_min, vmax=cbar_max)
axes[1,1].set_title(f"(prof_T - prof_Tclim)**2/prof_Tuncertainty**2")
cbar = fig.colorbar(a4, ax=axes[1,1], extend='max')

quantiles_fractions = [0.1, 0.15, 0.18, 0.25, 0.5, 0.75, 0.9] 
quantiles = np.quantile(values[values.notnull()].data, quantiles_fractions)
quantiles_strings = [f'{int(100 * qval)}%' for qval in quantiles_fractions]

text_color = "black"
text_va = "center"
cbar_norm = colors.Normalize(vmin=cbar_min, vmax=cbar_max)

for ii in range(len(quantiles)):

    quantile_norm = cbar_norm(quantiles[ii])

    cbar.ax.hlines(y=quantile_norm, xmin=0, xmax=1, color=text_color, zorder=3, linewidth=0.2,
                   transform=cbar.ax.transAxes, clip_on=True)

    cbar.ax.text(x=0.5, y=quantile_norm+0.01, s=quantiles_strings[ii], color=text_color,
                 va=text_va, ha='center', fontsize='xx-small',
                 transform=cbar.ax.transAxes, clip_on=True)


num_survivors = len(values[values.notnull()][values[values.notnull()] < cbar_max])
num_total = len(values[values.notnull()])
valid_fraction = num_survivors/num_total

fig.suptitle(f"input file: {input_file}\nlevel: {vertical_level} (0 - {MITprof_ds['prof_T'].shape[-1]-1})\nsurvivors: {num_survivors}/{num_total} ({valid_fraction*100:.2f}%)\n(less than {int(np.sqrt(cbar_max))} times uncertainty from climatology)", ha='left', x=0.42)
#fig.suptitle(f"input file: {input_file}\nlevel: {vertical_level+1} (1 - {MITprof_ds['prof_T'].shape[-1]})\nsurvivors: {num_survivors}/{num_total} ({valid_fraction*100:.2f}%)\n(less than {int(np.sqrt(cbar_max))} times uncertainty from climatology)", ha='left', x=0.42)


plt.show()
