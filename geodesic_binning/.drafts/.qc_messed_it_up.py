from mpl_toolkits.axes_grid1 import make_axes_locatable
from functools import partial
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
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point
from shapely import get_coordinates as ShapelyCoordinates
from sklearn.cluster import KMeans


def pp(geodesic_bin_data, num_bins):

    # -------------------------
    # Hardcoded test parameters
    # -------------------------
    key_variable = "T"
    #key_depth = "05"
    key_depth = "00"
    key_anomaly_var_edge = "std"
    key_anomaly_var_face = "mean"
    key_anomaly_var_raw_values= "values"
    variable_units = "($^\circ$C)"
    scale_threshold = 10
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

    anomaly_var_raw_values_list = []
    anomaly_var_face_list = []
    anomaly_var_edge_list= []
    patch_list = []
    count_list = []
    count_relative_list = []
    linewidths = []
    edgecolors_list = []


    for index in geodesic_bin_data[key_variable][key_depth].keys():
        if geodesic_bin_data[key_variable][key_depth][index]["count"] == 1 and key_anomaly_var_edge == "std":
            continue
        gbd_polygon = geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        if gbd_polygon.size > 0:

            anomaly_var_raw_values_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_raw_values])

            anomaly_var_face_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face])
            anomaly_var_edge_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge])
            count_list.append(geodesic_bin_data[key_variable][key_depth][index]['count'])
            
            linewidth_pre = original_scale * geodesic_bin_data[key_variable][key_depth][index]['count']

            linewidths.append(linewidth_pre)

            edgecolors_list.append(cmap_edge(norm_edge(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge])))
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

    count_array = np.array(count_list)

    collection = PatchCollection(patch_list, cmap=cmap_face, norm=norm_face, linewidths=original_linewidths_plot, edgecolors=edgecolors_list, transform=ccrs.PlateCarree(), joinstyle='miter')

    collection.set_array(np.array(anomaly_var_face_list))

    fig = plt.figure(layout="constrained")

    #ax = plt.axes([0.25, 0.25, 0.5, 0.5], projection=ccrs.PlateCarree())
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.coastlines(color='black', linewidth=0.15)
    ax.patch.set_facecolor('#D9D9D9')

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    #ax.set_adjustable('datalim')
    ax.set_box_aspect(1)
    #ax.set_aspect('equal', adjustable='box')
    #ax.set_aspect('equal')

    fig = plt.gcf()
    fig.canvas.draw()

    collection_ax = ax.add_collection(collection)

    cbar_shrink = 0.5
    cbar_pad = 0.2
    cbar_aspect = 10

    cbar = plt.colorbar(collection, ax=ax, shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, label=rf'{key_variable} anomaly {variable_units}'+'\n\n(no extensions shown)')


    cbar_min, cbar_max = find_colorbar_limits(cbar, anomaly_var_face_list)
    cbar.ax.set_ylim(cbar_min, cbar_max)

    cbar_std = plt.colorbar(mpl.cm.ScalarMappable(norm=norm_edge, cmap=cmap_edge),
             ax=ax, orientation='vertical', label=rf'{key_variable} anomaly std {variable_units}', shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, location='left')

    quartiles = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    quartiles_edgecolors = np.quantile(np.array(anomaly_var_edge_list), quartiles)

    quartiles_strings = [rf'$\downarrow${int(100 * qval)}%' for qval in quartiles]

    for q_dex in range(len(quartiles_edgecolors)):
        cbar_std.ax.axhline(quartiles_edgecolors[q_dex], color='white', zorder=3, linewidth=0.2)
        cbar_std.ax.text(x=0.5, y=quartiles_edgecolors[q_dex], s=quartiles_strings[q_dex], color='black',
             va='center', ha='center', fontsize='xx-small')

    original_range_x = ax.get_xlim()[1] - ax.get_xlim()[0]  
    original_range_y = ax.get_ylim()[1] - ax.get_ylim()[0]  
    original_area = original_range_x * original_range_y


    bound_callback = partial(scale_with_zoom, colorbar=cbar, original_area=original_area, original_linewidths_raw=original_linewidths_raw, linewidth_floor=linewidth_floor, count_array=count_array, patch_list=patch_list, anomaly_var_raw_values_list=anomaly_var_raw_values_list, cmap_face=cmap_face, norm_face=norm_face, edgecolors_list=edgecolors_list, anomaly_var_face_list=anomaly_var_face_list, scale_threshold=scale_threshold)


    ax.callbacks.connect('xlim_changed', bound_callback)
    ax.callbacks.connect('ylim_changed', bound_callback)

    #right = ax.get_xlim()[1]
    #top = ax.get_ylim()[1] 
    data_to_axes_space = ax.transLimits.transform((xmax, ymax))
    norm_xmax, norm_ymax = data_to_axes_space[0], data_to_axes_space[1]

    # Only runs for first plot, before zooming
    custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, ax)
    if custom_handles != 0:
        ax.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                  handles=custom_handles, title=f"{legend_title}")
                  #handles=custom_handles, title=f"{legend_title}", bbox_to_anchor=(norm_xmax, norm_ymax), bbox_transform=ax.transAxes)


    ax.gridlines(draw_labels=True)

    plt.title(f"variable: {key_variable }\ndepth level: {key_depth }\nnum bins populated: {len(count_array)}/{num_bins}\nnum profiles binned: {np.sum(count_array)}\n\n")

    caption_string = ("Within each bin, face color corresponds to anomaly value, and edge color corresponds to anomaly standard deviation (std).\n\n"
                      "At low-moderate zoom levels, bin face color represents bin anomaly mean, and bin edge width scales linearly with bin profile count.\n\n"
                      "At higher zoom levels, bins are sub-divided into equal-area polygons representing individual profiles within the bin "
                      "(locations ignored) , with face colors "
                      "indicating profile anomaly values and edge colors still representing overall bin anomaly standard deviation.  Note that these"
                      "sub-polygons do not indicate profile location within a bin")

    '''
    caption = cbar.ax.text(
            0.5, -0.25, caption_string, ha='center', va='top', wrap=True, style='italic',
            transform=ax.transAxes,
            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none')
                  )
    '''

    '''
    caption = cbar.ax.annotate(
            caption_string,
            xy=(0.5, 0),
            xycoords='axes fraction',
            xytext=(0,-0.2),
            textcoords='offset points',
            ha='center', va='top',
            style='italic',
            wrap=True,
            transform=ax.transAxes,
            bbox=dict(boxstyle='square,pad=0', fc='none', ec='none'),
            )
    '''

    caption = plt.figtext(0.46, 0.08, caption_string, ha="center", fontsize="small", style="italic")
    #plt.figtext(0.46, 0.08, caption_string_wrapped, ha="center", fontsize="small", style="italic", wrap=True)

    wrap_width = 0.7 
    caption._get_wrap_line_width = lambda: fig.bbox.width * wrap_width

    plt.show()


def determine_sub_polygons(axes, patch_list, count_array, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list):

    visible_patch_mask = np.zeros(len(patch_list)).astype(bool)

    for patch_dex in range(len(patch_list)):
        patch_coords = patch_list[patch_dex].get_xy()
        for coord_dex in range(patch_coords.shape[0]):
            if (patch_coords[coord_dex,0] > axes.get_xlim()[0]
                and patch_coords[coord_dex,0] < axes.get_xlim()[1] 
                and patch_coords[coord_dex,1] > axes.get_ylim()[0] 
                and patch_coords[coord_dex,1] < axes.get_ylim()[1]):
                    visible_patch_mask[patch_dex] = True 
                    break
    if np.sum(visible_patch_mask) == 0:
        return 0, 0

    counts_zoom = count_array[visible_patch_mask]

    patches_zoom = [patch_list[ii] for ii in range(len(patch_list)) if visible_patch_mask[ii]]
    edgecolors_zoom = [edgecolors_list[ii] for ii in range(len(edgecolors_list)) if visible_patch_mask[ii]]
    anomalies_zoom = [anomaly_var_raw_values_list[ii] for ii in range(len(anomaly_var_raw_values_list)) if visible_patch_mask[ii]]

    mini_patches_list = []
    mini_patches_anomaly_list = []
    mini_patches_edgecolors_list = []

    for patch_dex in range(len(patches_zoom)):

        patch_vertices = patches_zoom[patch_dex].get_xy()
        num_profiles = counts_zoom[patch_dex]

        orig_poly = ShapelyPolygon(patch_vertices)

        # VIBING OUT
        num_samples = 2000
        minx, miny, maxx, maxy = orig_poly.bounds
        points = []

        while len(points) < num_samples:
            p = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
            if orig_poly.contains(p):
                points.append([p.x, p.y])

        random_points_array = np.array(points)

        kmeans = KMeans(n_clusters=num_profiles, n_init=10, random_state=42)
        labels = kmeans.fit_predict(random_points_array)

        for profile_index in range(num_profiles):
            cluster_points = random_points_array[labels == profile_index]
            polygon_coords = ShapelyCoordinates(ShapelyPolygon(cluster_points).convex_hull)
            mini_patches_list.append(patches.Polygon(polygon_coords, closed=True))

        mini_patches_anomaly_list += anomalies_zoom[patch_dex]
        mini_patches_edgecolors_list += [edgecolors_zoom[patch_dex]] * num_profiles



    subcol = PatchCollection(mini_patches_list, cmap=cmap_face, norm=norm_face, linewidths=1, edgecolors=mini_patches_edgecolors_list, transform=ccrs.PlateCarree(), joinstyle='miter')

    subcol.set_array(np.array(mini_patches_anomaly_list)) 

    return subcol, mini_patches_anomaly_list




def scale_with_zoom(axes, colorbar, original_area, original_linewidths_raw, linewidth_floor, count_array, patch_list, anomaly_var_face_list, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list, scale_threshold):

    current_range_x = axes.get_xlim()[1] - axes.get_xlim()[0]  
    current_range_y = axes.get_ylim()[1] - axes.get_ylim()[0]  
    current_area = current_range_x * current_range_y

    scale_factor = np.sqrt(original_area / current_area)

    data_to_axes_space = axes.transLimits.transform((axes.get_xlim()[1], axes.get_ylim()[1]))
    norm_xmax, norm_ymax = data_to_axes_space[0], data_to_axes_space[1]


    current_collection=axes.collections[-1]
    if type(current_collection) is PatchCollection:
        current_collection.remove() # erase the old collection/plot, start fresh.  maybe unecessary, but just want to get this working for now

    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()

    if scale_factor < scale_threshold:

        #---------------------------------------------------
        collection = PatchCollection(patch_list, cmap=cmap_face, norm=norm_face, edgecolors=edgecolors_list, transform=ccrs.PlateCarree(), joinstyle='miter')
        collection.set_array(np.array(anomaly_var_face_list))
        collection_ax = axes.add_collection(collection)
        #---------------------------------------------------

        new_linewidth_max = scale_factor # Random choice, but seems to do the job
        new_linewidths = np.clip(original_linewidths_raw * scale_factor, linewidth_floor, new_linewidth_max)
        collection.set_linewidths(new_linewidths)

        custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, axes)
        if custom_handles != 0:
            axes.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                      #handles=custom_handles, title=f"{legend_title}", bbox_to_anchor=(norm_xmax, norm_ymax), bbox_transform=axes.transAxes)
                      #handles=custom_handles, title=f"{legend_title}", loc='upper right', bbox_to_anchor=(right, top), bbox_transform=axes.transData)
                      handles=custom_handles, title=f"{legend_title}")

        cbar_min, cbar_max = find_colorbar_limits(colorbar, anomaly_var_face_list)
        colorbar.ax.set_ylim(cbar_min, cbar_max)

    else:
        collection, new_anomaly_var_face_list = determine_sub_polygons(axes, patch_list, count_array, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list)

        if collection != 0:
            custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, axes)
            axes.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                      #handles=custom_handles, title=f"{legend_title}", bbox_to_anchor=(norm_xmax, norm_ymax), bbox_transform=axes.transAxes)
                      #handles=custom_handles, title=f"{legend_title}", loc='upper right', bbox_to_anchor=(right, top), bbox_transform=axes.transData)
                      handles=custom_handles, title=f"{legend_title}")


            collection_ax = axes.add_collection(collection)
            collection.set_linewidths(1)
            cbar_min, cbar_max = find_colorbar_limits(colorbar, new_anomaly_var_face_list)
            colorbar.ax.set_ylim(cbar_min, cbar_max)



def find_colorbar_limits(colorbar, value_list):

    cbar_min = np.min(value_list)
    cbar_max = np.max(value_list)

    extend_up = False
    extend_down = False

    if cbar_max > colorbar.norm.vmax:
        cbar_max = colorbar.norm.vmax
        extend_up = True
    if cbar_min < colorbar.norm.vmin:
        cbar_min = colorbar.norm.vmin
        extend_down = True

    return cbar_min, cbar_max
    #return cbar_min, cbar_max, extend_down, extend_up




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
            legend_title = "profiles per patch:"

        elif len(np.unique(count_array[visible_patch_mask])) == 2: 
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"min {count_min_zoom}"),
                Line2D([0], [0], color='gray', label=f"max {count_max_zoom}")
            ]
            legend_title = "profiles per patch:"

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
            legend_title = "profiles per patch:"

        return custom_handles, legend_title




