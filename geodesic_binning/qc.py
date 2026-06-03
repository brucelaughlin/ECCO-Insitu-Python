import pdb
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection
import numpy as np

def qc(lons, lats):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()

    plt.plot(lons, lats, 'k.', transform=ccrs.PlateCarree())
    ax.gridlines(draw_labels=True)
    plt.show()


def pp(geodesic_bin_data):

    # -------------------------
    # Hardcoded test parameters
    # -------------------------
    key_variable = "T"
    key_depth = "00"
    key_anomaly_var_edge = "std"
    key_anomaly_var_face = "mean"
    # -------------------------

    #line_width_max = 100
    #line_width_max = 10
    line_width_max = 5

    dummyMegaNumber = 1e30

    xmin,ymin = dummyMegaNumber, dummyMegaNumber
    xmax,ymax = -dummyMegaNumber, -dummyMegaNumber

    value_min_edge = dummyMegaNumber
    value_max_edge = -dummyMegaNumber
    value_min_face = dummyMegaNumber
    value_max_face = -dummyMegaNumber

    count_min = dummyMegaNumber 
    count_max = 0
    count_list = []

    for index in geodesic_bin_data[key_variable][key_depth].keys():

        if geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"].size == 0:
            continue

        count_list.append(geodesic_bin_data[key_variable][key_depth][index]["count"])

        if geodesic_bin_data[key_variable][key_depth][index]["count"] > 1:
            if geodesic_bin_data[key_variable][key_depth][index]["count"] < count_min:
                count_min = geodesic_bin_data[key_variable][key_depth][index]["count"] 
            if geodesic_bin_data[key_variable][key_depth][index]["count"] > count_max:
                count_max = geodesic_bin_data[key_variable][key_depth][index]["count"] 

        if geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge] < value_min_edge:
            value_min_edge = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge]
        if geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge] > value_max_edge:
            value_max_edge = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge]

        if geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face] < value_min_face:
            value_min_face = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face]
        if geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face] > value_max_face:
            value_max_face = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face]

    '''
    count_width = count_max - count_min
    if count_width == 0:
        count_width = 1
    '''
    count_mean = np.mean(count_list)

    count_log2 = np.log2(count_list) 

    #print(f"{count_max}, {count_min}, {count_width}")


    norm_face = mcolors.CenteredNorm(vcenter=0)
    #norm_face = mcolors.TwoSlopeNorm(vmin=value_min_face, vcenter=0, vmax=value_max_face)
    cmap_face = cm.get_cmap('RdBu_r')
    #cmap_face = cm.get_cmap('RdBu')
    norm_edge = mcolors.Normalize(vmin=value_min_edge, vmax=value_max_edge)
    #norm_edge = mcolors.Normalize(vmin=count_min, vmax=count_max)
    cmap_edge = cm.get_cmap('viridis')
    #cmap_edge = cm.get_cmap('viridis_r')


    value_list = []
    patch_list = []
    count_list = []
    count_relative_list = []
    linewidths = []
    edgecolors = []


    for index in geodesic_bin_data[key_variable][key_depth].keys():
        if geodesic_bin_data[key_variable][key_depth][index]["count"] == 1 and key_anomaly_var_edge == "std":
            continue
        gbd_polygon = geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        value_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face])
        #value_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge])
        count_list.append(geodesic_bin_data[key_variable][key_depth][index]['count'])
        #count_relative_list.append(geodesic_bin_data[key_variable][key_depth][index]['count']/count_width)
        if gbd_polygon.size > 0:
            
            line_width_pre = np.log(geodesic_bin_data[key_variable][key_depth][index]['count'])
            #line_width_pre = np.log2(geodesic_bin_data[key_variable][key_depth][index]['count'])
            if line_width_pre > 5:
                print(line_width_pre)
            #line_width_pre = geodesic_bin_data[key_variable][key_depth][index]['count']/count_mean
            '''
            #line_width_pre = line_width_max * (geodesic_bin_data[key_variable][key_depth][index]['count'] - count_min)/count_mean
            #line_width_pre = line_width_max * (geodesic_bin_data[key_variable][key_depth][index]['count'] - count_min)/count_width
            if line_width_pre > line_width_max/2:
                print(f"{index}: {line_width_pre}")
                linewidths.append(line_width_max/3) # No idea why this fixes the weird circle bug around the single thickest border
            elif line_width_pre == 0:
                linewidths.append(1)
            elif line_width_pre < 1:
                linewidths.append(2)
            else:
                linewidths.append(line_width_pre)
            '''

            linewidths.append(line_width_pre)

            #linewidths.append(1)

            #print(f"{count_width}, {geodesic_bin_data[key_variable][key_depth][index]['count']}, {linewidths[-1]}")
            #linewidths.append(line_width_max * geodesic_bin_data[key_variable][key_depth][index]['count']/count_width)
            edgecolors.append(cmap_edge(norm_edge(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge])))
            patch_list.append(patches.Polygon(gbd_polygon, closed=True))
            #patch_list.append(patches.Polygon(gbd_polygon, closed=True, linewidth=linewidth, edgecolor='black'))
            #patch_list.append(patches.Polygon(gbd_polygon, closed=True, linewidth=linewidth, edgecolor=edgecolor))

            #patch_list.append(patches.Polygon(gbd_polygon, closed=True), linewidth=patch_border_width_max)
            #patch_list.append(patches.Polygon(gbd_polygon, closed=True))
            if xmin > np.min(gbd_polygon[:,0]):
                xmin = np.min(gbd_polygon[:,0])

            if xmax < np.max(gbd_polygon[:,0]):
                xmax = np.max(gbd_polygon[:,0])

            if ymin > np.min(gbd_polygon[:,1]):
                ymin = np.min(gbd_polygon[:,1])

            if ymax < np.max(gbd_polygon[:,1]):
                ymax = np.max(gbd_polygon[:,1])

            #break

    
    count_range= np.max(count_list) - np.min(count_list)

    counts_normalized = count_list/count_range


    quantile_cutoff = 0.9
    #quantile_cutoff = 0.75
    val_max = np.quantile(np.array(value_list), quantile_cutoff)

    quartiles = [0.25, 0.5, 0.75]
    quartile_values = np.quantile(np.array(value_list), quartiles)
    #quartile_labels = [f"{quartile*100}%" for quartile in quartiles]
    #quartile_labels = [f"{quart*100}%" for val, quart in zip(quartiles, quartile_values)]

    #pdb.set_trace()

    col = PatchCollection(patch_list, cmap=cmap_face, norm=norm_face, linewidths=linewidths, edgecolors=edgecolors, antialiased=True, transform=ccrs.PlateCarree())
    #col = PatchCollection(patch_list, cmap=cmap_face, norm=norm_face, linewidths=linewidths, edgecolors=edgecolors, transform=ccrs.PlateCarree())
    #col = PatchCollection(patch_list, cmap=cmap_face, norm=norm_face, linewidths=linewidths, edgecolor="black", transform=ccrs.PlateCarree())
    #col = PatchCollection(patch_list, cmap=cmap_face, norm=norm_face, transform=ccrs.PlateCarree())
    col.set_array(value_list) 

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(color='black', linewidth=0.1)

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    #ax.set_xlim(-20,0)
    #ax.set_ylim(60,81)

    collection_ax = ax.add_collection(col)
    #collection_ax.set_clim(np.min(value_list), np.max(value_list))
    #collection_ax.set_clim(np.min(value_list), val_max)

    cbar_shrink = 0.5
    cbar_pad = 0.12
    cbar_aspect = 10

    cbar = plt.colorbar(col, ax=ax, shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, label='anomaly mean')
    #cbar = plt.colorbar(col, ax=ax, shrink=0.5, pad=0.12, aspect=10, extend="both")
    #cbar = plt.colorbar(col, ax=ax, shrink=0.5, pad=0.12, aspect=10, ticks=quartile_values)

    cbar.ax.set_ylim(np.min(value_list), np.max(value_list))

    print(np.min(value_list))
    print(np.max(value_list))

    plt.colorbar(mpl.cm.ScalarMappable(norm=norm_edge, cmap=cmap_edge),
             ax=ax, orientation='vertical', label='anomaly std', shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, location='left')

    #default_ticks = list(cbar.ax.get_yticks())
    #new_ticks = sorted(default_ticks + list(quartile_values))
    #cbar.set_ticks(new_ticks)

    #collection_ax.set_clim(np.min(value_list), np.max(value_list))

    '''
    for qval in quartile_values:
        cbar.ax.axhline(qval, color='white', zorder=3)
    '''


    #ax.gridlines(draw_labels=True)

    plt.show()


