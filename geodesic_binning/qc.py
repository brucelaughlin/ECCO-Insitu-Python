import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.collections import PatchCollection
import numpy as np

def qc(lons, lats):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()

    plt.plot(lons, lats, 'k.', transform=ccrs.PlateCarree())
    ax.gridlines(draw_labels=True)
    plt.show()


def pp(geodesic_bin_data):

    key_variable = "T"
    key_depth = "00"
    key_anomaly_var = "std"

    xmin,ymin = 1000, 1000
    xmax,ymax = -1000, -1000

    value_min = 1000
    value_max = -1000

    for index in geodesic_bin_data[key_variable][key_depth].keys():
        gbd_polygon = geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        gbd_anomaly_var = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var]
        if gbd_polygon.size > 0:
           
            if value_min > gbd_anomaly_var:
                value_min = gbd_anomaly_var
                
            if value_max < gbd_anomaly_var:
                value_max = gbd_anomaly_var

    norm = colors.Normalize(vmin=value_min, vmax=value_max)
    cmap = plt.colormaps.get_cmap('jet')
    #cmap = cm.get_cmap('viridis')


    # ax = plt.axes(projection=ccrs.PlateCarree())
    # ax.coastlines()

    value_list = []
    patch_list = []

    for index in geodesic_bin_data[key_variable][key_depth].keys():
        gbd_polygon = geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        #color_val = cmap(norm(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var]))
        val = norm(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var])
        value_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var])
        if gbd_polygon.size > 0:
            #ax.add_patch(
                #patches.Polygon(gbd_polygon, closed=True, facecolor=color_val, alpha=0.3, edgecolor='black'))
                #patches.Polygon(gbd_polygon, closed=True, facecolor='green', alpha=0.3, edgecolor='black'))

            patch_list.append(patches.Polygon(gbd_polygon, closed=True))
            value_list.append(val)
            #value_list.append(color_val)
            
            if xmin > np.min(gbd_polygon[:,0]):
                xmin = np.min(gbd_polygon[:,0])

            if xmax < np.max(gbd_polygon[:,0]):
                xmax = np.max(gbd_polygon[:,0])

            if ymin > np.min(gbd_polygon[:,1]):
                ymin = np.min(gbd_polygon[:,1])

            if ymax < np.max(gbd_polygon[:,1]):
                ymax = np.max(gbd_polygon[:,1])

            #break

    col = PatchCollection(patch_list, cmap=cmap, transform=ccrs.PlateCarree())
    col.set_array(value_list) 

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    ax.add_collection(col)
    plt.colorbar(col, ax=ax)  
    ax.gridlines(draw_labels=True)

    #profiles_ds = xr.open_dataset(profile_file)
    #plt.plot(profiles_ds['prof_lon'].values, profiles_ds['prof_lat'].values, 'k.', transform=ccrs.PlateCarree())

    plt.show()


