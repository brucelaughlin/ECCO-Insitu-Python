from matplotlib.lines import Line2D
import textwrap
import math
import pdb
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection
import numpy as np

'''
def qc(lons, lats):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()

    plt.plot(lons, lats, 'k.', transform=ccrs.PlateCarree())
    ax.gridlines(draw_labels=True)
    plt.show()
'''


def pp(geodesic_bin_data, num_bins):

    # -------------------------
    # Hardcoded test parameters
    # -------------------------
    key_variable = "T"
    key_depth = "00"
    key_anomaly_var_edge = "std"
    key_anomaly_var_face = "mean"
    variable_units = "($^\circ$C)"
    # -------------------------

    dummyMegaNumber = 1e30

    xmin,ymin = dummyMegaNumber, dummyMegaNumber
    xmax,ymax = -dummyMegaNumber, -dummyMegaNumber

    value_min_edge = dummyMegaNumber
    value_max_edge = -dummyMegaNumber
    value_min_face = dummyMegaNumber
    value_max_face = -dummyMegaNumber

    count_min = dummyMegaNumber 
    count_max = 0

    for index in geodesic_bin_data[key_variable][key_depth].keys():

        if geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"].size == 0:
            continue

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


    norm_face = mcolors.CenteredNorm(vcenter=0)
    #norm_face = mcolors.TwoSlopeNorm(vmin=value_min_face, vcenter=0, vmax=value_max_face)
    cmap_face = cm.get_cmap('PRGn')
    #cmap_face = cm.get_cmap('RdBu_r')
    #cmap_face = cm.get_cmap('RdBu')
    norm_edge = mcolors.Normalize(vmin=value_min_edge, vmax=value_max_edge)
    #norm_edge = mcolors.Normalize(vmin=count_min, vmax=count_max)
    #cmap_edge = cm.get_cmap('viridis')
    #cmap_edge = cm.get_cmap('viridis_r')
    cmap_edge = cm.get_cmap('cividis_r')
    #cmap_edge = cm.get_cmap('cividis')

    linewidth_floor = 0
    #original_linewidth_max = 0.5
    original_linewidth_max = 1
    #original_linewidth_max = 0.25

    original_scale = 0.01
    #original_scale = 0.05
    #original_scale = 0.1

    anomaly_var_face_list = []
    anomaly_var_edge_list= []
    patch_list = []
    count_list = []
    count_relative_list = []
    linewidths = []
    edgecolors = []


    for index in geodesic_bin_data[key_variable][key_depth].keys():
        if geodesic_bin_data[key_variable][key_depth][index]["count"] == 1 and key_anomaly_var_edge == "std":
            continue
        gbd_polygon = geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        #anomaly_var_face_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face])
        #count_list.append(geodesic_bin_data[key_variable][key_depth][index]['count'])
        if gbd_polygon.size > 0:

            anomaly_var_face_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face])
            anomaly_var_edge_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge])
            count_list.append(geodesic_bin_data[key_variable][key_depth][index]['count'])
            
            linewidth_pre = original_scale * geodesic_bin_data[key_variable][key_depth][index]['count']

            linewidths.append(linewidth_pre)

            edgecolors.append(cmap_edge(norm_edge(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge])))
            patch_list.append(patches.Polygon(gbd_polygon, closed=True))

            if xmin > np.min(gbd_polygon[:,0]):
                xmin = np.min(gbd_polygon[:,0])

            if xmax < np.max(gbd_polygon[:,0]):
                xmax = np.max(gbd_polygon[:,0])

            if ymin > np.min(gbd_polygon[:,1]):
                ymin = np.min(gbd_polygon[:,1])

            if ymax < np.max(gbd_polygon[:,1]):
                ymax = np.max(gbd_polygon[:,1])

            #break


    original_linewidths_raw = np.array(linewidths)
    original_linewidths_plot = np.clip(original_linewidths_raw, linewidth_floor, original_linewidth_max)

    face_array = np.array(anomaly_var_face_list)
    count_array = np.array(count_list)
    edge_array = np.array(anomaly_var_edge_list)

    col = PatchCollection(patch_list, cmap=cmap_face, norm=norm_face, linewidths=original_linewidths_plot, edgecolors=edgecolors, transform=ccrs.PlateCarree(), joinstyle='miter')

    col.set_array(face_array) 
    #col.set_array(anomaly_var_face_list) 

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(color='black', linewidth=0.1)
    ax.patch.set_facecolor('#D9D9D9')
    #ax.patch.set_facecolor('darkgrey')
    #ax.patch.set_facecolor('lightgrey')

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    #ax.set_xlim(-20,0)
    #ax.set_ylim(60,81)

    collection_ax = ax.add_collection(col)

    cbar_shrink = 0.5
    cbar_pad = 0.2
    #cbar_pad = 0.12
    cbar_aspect = 10

    cbar = plt.colorbar(col, ax=ax, shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, label=rf'{key_variable} anomaly mean {variable_units}')

    cbar.ax.set_ylim(np.min(anomaly_var_face_list), np.max(anomaly_var_face_list))

    print(np.min(anomaly_var_face_list))
    print(np.max(anomaly_var_face_list))

    cbar_std = plt.colorbar(mpl.cm.ScalarMappable(norm=norm_edge, cmap=cmap_edge),
             ax=ax, orientation='vertical', label=rf'{key_variable} anomaly std {variable_units}', shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, location='left')


    quartiles = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    #quartiles = [0.25, 0.5, 0.75, 0.9]
    #quartiles = [0.25, 0.5, 0.75]
    quartiles_edgecolors = np.quantile(edge_array, quartiles)

    quartiles_strings = [rf'$\downarrow${int(100 * qval)}%' for qval in quartiles]
    #quartiles_strings = [rf'{int(100 * qval)}%$\downarrow$' for qval in quartiles]

    for q_dex in range(len(quartiles_edgecolors)):
    #for qval in quartiles_edgecolors:
        #cbar_std.ax.axhline(qval color='white', zorder=3)
        cbar_std.ax.axhline(quartiles_edgecolors[q_dex], color='white', zorder=3)
        cbar_std.ax.text(x=0.5, y=quartiles_edgecolors[q_dex], s=quartiles_strings[q_dex], color='black',
             va='center', ha='center', fontsize='xx-small')
             #va='center', ha='left', fontsize='xx-small')

    original_range_x = ax.get_xlim()[1] - ax.get_xlim()[0]  
    original_range_y = ax.get_ylim()[1] - ax.get_ylim()[0]  
    original_area = original_range_x * original_range_y

    def scale_with_zoom(axes):

        current_range_x = ax.get_xlim()[1] - ax.get_xlim()[0]  
        current_range_y = ax.get_ylim()[1] - ax.get_ylim()[0]  
        current_area = current_range_x * current_range_y

        scale_factor = np.sqrt(original_area / current_area)

        new_linewidth_max = scale_factor 

        new_linewidths = np.clip(original_linewidths_raw * scale_factor, linewidth_floor, new_linewidth_max)

        col.set_linewidths(new_linewidths)

        custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, ax)

        if custom_handles != 0:
            ax.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                      handles=custom_handles, title=f"{legend_title}")

        else:
            plt.gca().get_legend().remove()


    ax.callbacks.connect('xlim_changed', scale_with_zoom)
    ax.callbacks.connect('ylim_changed', scale_with_zoom)


    # Only runs for first plot, before zooming
    custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, ax)
    if custom_handles != 0:
        ax.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                  handles=custom_handles, title=f"{legend_title}")

    ax.gridlines(draw_labels=True)

    plt.title(f"variable: {key_variable }\ndepth level: {key_depth }\nnum bins populated: {len(count_array)}/{num_bins}\nnum profiles binned: {np.sum(count_array)}\n\n")

    caption_string = ("As you zoom, bin edge widths will scale with profile counts. Within each bin, face color corresponds to anomaly mean, edge color "
                      "corresponds to anomaly standard deviation (std), and edge width corresponds to the number of profiles within the bin.")

    wrap_width = 80

    caption_string_wrapped = "\n".join(textwrap.wrap(caption_string, width=wrap_width))

    plt.figtext(0.46, 0.15, caption_string_wrapped, ha="center", fontsize="small", style="italic")

    plt.show()



def make_handles_and_titles(patch_list, count_array, ax):

    visible_patch_mask = np.zeros(len(patch_list)).astype(bool)

    for patch_dex in range(len(patch_list)):

        patch_coords = patch_list[patch_dex].get_xy()

        for coord_dex in range(patch_coords.shape[0]):

            if (patch_coords[coord_dex,0] > ax.get_xlim()[0]
                and patch_coords[coord_dex,0] < ax.get_xlim()[1] 
                and patch_coords[coord_dex,1] > ax.get_ylim()[0] 
                and patch_coords[coord_dex,1] < ax.get_ylim()[1]):

                    visible_patch_mask[patch_dex] = True 
                    break

    if np.sum(visible_patch_mask) == 0:
        return 0, 0

    else:
        count_min_zoom = np.min(count_array[visible_patch_mask])
        count_max_zoom = np.max(count_array[visible_patch_mask])

        if count_min_zoom > 1:
           legend_min_string = "profiles"
        else:
           legend_min_string = "profile"

        if count_max_zoom > 1:
           legend_max_string = "profiles"
        else:
           legend_max_string = "profile"

        if len(np.unique(count_array[visible_patch_mask])) == 1: 
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"{count_min_zoom}"),
            ]
            legend_title = "profiles per patch"

        elif len(np.unique(count_array[visible_patch_mask])) == 2: 
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"min {count_min_zoom}"),
                Line2D([0], [0], color='gray', label=f"max {count_max_zoom}")
            ]
            legend_title = "profiles per patch"

        else:
            count_mid_zoom = int(np.median(np.sort(count_array[visible_patch_mask])))
            if count_mid_zoom > 1:
               legend_mid_string = "profiles"
            else:
               legend_mid_string = "profile"

            custom_handles = [
                Line2D([0], [0], color='gray', label=f"min {count_min_zoom}"),
                Line2D([0], [0], color='gray', label=f"med {count_mid_zoom}"),
                Line2D([0], [0], color='gray', label=f"max {count_max_zoom}")
            ]
            legend_title = "profiles per patch"

        return custom_handles, legend_title




