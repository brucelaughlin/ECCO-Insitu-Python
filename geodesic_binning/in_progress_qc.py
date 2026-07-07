from pathlib import Path
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
import xarray as xr


def pp(geodesic_bin_data, num_bins, geodesic_file, profile_file, num_subpolygons_max):

    """

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

    '''
    count_min = dummyMegaNumber 
    count_max = 0
    '''

    profiles_lats = []
    profiles_lons = []
    bin_counts = []

    num_zero_area_bins = 0

    for index in geodesic_bin_data[key_variable][key_depth].keys():

        # Had to add this bc of the profiles_lons/lats, which aren't specific to any bins 
        if not isinstance(geodesic_bin_data[key_variable][key_depth][index], dict):
            continue

        if geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"].size == 0:
            num_zero_area_bins += 1

        '''
        if geodesic_bin_data[key_variable][key_depth][index]["count"] > 1:
            if geodesic_bin_data[key_variable][key_depth][index]["count"] < count_min:
                count_min = geodesic_bin_data[key_variable][key_depth][index]["count"] 
            if geodesic_bin_data[key_variable][key_depth][index]["count"] > count_max:
                count_max = geodesic_bin_data[key_variable][key_depth][index]["count"] 
        '''

        if geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge] < value_min_edge:
            value_min_edge = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge]
        if geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge] > value_max_edge:
            value_max_edge = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge]

        if geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face] < value_min_face:
            value_min_face = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face]
        if geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face] > value_max_face:
            value_max_face = geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face]

        bin_counts.append(geodesic_bin_data[key_variable][key_depth][index]['count'])

    profiles_lats = geodesic_bin_data[key_variable][key_depth]['profiles_lats']
    profiles_lons = geodesic_bin_data[key_variable][key_depth]['profiles_lons']

    print(f'debug: {num_zero_area_bins} zero-area bins encountered')

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
        #if geodesic_bin_data[key_variable][key_depth][index]["count"] == 1 and key_anomaly_var_edge == "std":
        #    continue

        # Had to add this bc of the profiles_lons/lats, which aren't specific to any bins 
        if not isinstance(geodesic_bin_data[key_variable][key_depth][index], dict):
            continue

        gbd_polygon = geodesic_bin_data[key_variable][key_depth][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        if gbd_polygon.size > 0:

            anomaly_var_raw_values_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_raw_values])

            anomaly_var_face_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_face])
            anomaly_var_edge_list.append(geodesic_bin_data[key_variable][key_depth][index][key_anomaly_var_edge])
            count_list.append(geodesic_bin_data[key_variable][key_depth][index]['count'])
            
            if geodesic_bin_data[key_variable][key_depth][index]['count'] == 1:
                linewidth_pre = 0
            else:
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


    """

    # Now, all basic information needed to make patch collections should be loaded from a file generated by functions in the utilities module.
    # Note: I am gonna save with npz.  so do all patch and patchcollection making here

    patch_list.append(patches.Polygon(gbd_polygon, closed=True))

    patch_collection = PatchCollection(patch_list, transform=ccrs.PlateCarree(), joinstyle='miter')
    patch_collection.set_array(np.array(anomaly_var_face_list))
    patch_collection.set_linewidths(original_linewidths_plot)
    patch_collection.set_edge_colors(edgecolors=edgecolors_list)
    patch_collection.set_cmap(cmap_face)
    patch_collection.set_norm(norm_face)




    fig_width, fig_height = 14, 6
    fig = plt.figure(figsize=(fig_width, fig_height), facecolor='lightskyblue', layout='constrained')

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(color='black', linewidth=0.15)
    ax.patch.set_facecolor('#D9D9D9')

    xmin, ymin, xmax, ymax = collection.get_extents().extents
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    fig.add_axes(ax)
    ax.add_collection(collection)

    original_range_x = ax.get_xlim()[1] - ax.get_xlim()[0]  
    original_range_y = ax.get_ylim()[1] - ax.get_ylim()[0]  
    original_area = original_range_x * original_range_y

    # sure, why not
    cbar_shrink = 0.5
    cbar_pad = 0.2
    cbar_aspect = 10

    cbar = plt.colorbar(collection, ax=ax, shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, label=rf'{key_variable} anomaly mean {variable_units}'+'\n\n(no extensions shown)')

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

    # Only runs for first plot, before zooming
    custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, ax)
    if custom_handles != 0:
        ax.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                  handles=custom_handles, title=f"{legend_title}")


    # Hardcoding "globe_area" to be that of the globe, since otherwise it may be impossible to get the sub-polygons
    # to plot (ie when only a few nearby bins are populated, the first plot will already be somewhat zoomed in, and
    # "original_area" will be small to begin with and thus we won't get past the ratio threshold when zooming further in 
    # before the program freaks out bc we've zoomed too far in to plot a whole bin polygon. 
    globe_range_x = 360
    globe_range_y = 180 
    globe_area = globe_range_x * globe_range_y

    #pdb.set_trace()

    scale_factor = np.sqrt(globe_area/original_area)

    if scale_factor > scale_threshold: 

        ax.collections[-1].remove()

        if plt.gca().get_legend() is not None:
            plt.gca().get_legend().remove()

        collection, new_anomaly_var_face_list = determine_sub_polygons(ax, patch_list, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list, num_subpolygons_max)

        if collection != 0:
            custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, ax)
            if custom_handles != 0:
                ax.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                          handles=custom_handles, title=f"{legend_title}")


            ax.scatter(profiles_lons, profiles_lats, c='red', s=1, zorder=10)
            ax.add_collection(collection)
            cbar_min, cbar_max = find_colorbar_limits(cbar, new_anomaly_var_face_list)
            cbar.ax.set_ylim(cbar_min, cbar_max)


    # !!!!!!!!!!
    # !!!!!!!!!!
    # !!!!!!!!!!
    # More hacky stuff to only remove an axes collection when the last collection added was a patch collection.
    # Used in my hack to remove the previous patch collection, in case we zoom in and would see it underneath the sub-polygons.
    # Note that this means the patch collection must be the last collection added to the axes (ie after calling ax.scatter(), ax.coastlines(), etc)
    num_artists_original = len(ax.collections)
    # !!!!!!!!!!
    # !!!!!!!!!!
    # !!!!!!!!!!

    bound_callback = partial(scale_with_zoom, num_artists_original=num_artists_original, colorbar=cbar, globe_area=globe_area, original_linewidths_raw=original_linewidths_raw, linewidth_floor=linewidth_floor, count_array=count_array, patch_list=patch_list, anomaly_var_raw_values_list=anomaly_var_raw_values_list, cmap_face=cmap_face, norm_face=norm_face, edgecolors_list=edgecolors_list, anomaly_var_face_list=anomaly_var_face_list, scale_threshold=scale_threshold, num_subpolygons_max=num_subpolygons_max, profiles_lons=profiles_lons, profiles_lats=profiles_lats)


    ax.callbacks.connect('xlim_changed', bound_callback)
    ax.callbacks.connect('ylim_changed', bound_callback)


    ax.gridlines(draw_labels=True)


    suptitle_string = (
            f"\nprofile_file: {geodesic_bin_data["profile_file_stem"]}\n"
            f"geodesic_bin_file: {geodesic_bin_data["geodesic_bin_file_stem"]}\n"
            f"variable: {key_variable }\n"
            f"depth level: {key_depth }\n"
            f"num bins populated: {len(count_array)}/{geodesic_bin_data["num_geodesic_bins"]}\n"
            f"num profiles binned: {np.sum(count_array)}\n\n"
            )

    fig.suptitle(suptitle_string)

    caption_string = (
            "Within each bin, face color corresponds to anomaly value, and edge color corresponds to anomaly standard deviation (std).\n"
            "At low-moderate zoom levels, bin face color represents bin anomaly mean, and bin edge width scales linearly with bin profile count.\n"
            "At higher zoom levels, bins are sub-divided into equal-area polygons representing individual profiles, with face colors "
            "indicating profile anomaly values and edge colors still representing overall bin anomaly standard deviation."
            )

    wrap_width = 100

    caption_string_wrapped_list = [textwrap.fill(paragraph, width=wrap_width) for paragraph in caption_string.split('\n')]
    caption_string_wrapped = '\n'.join(caption_string_wrapped_list)

    # VIBING
    ax.annotate(
        caption_string_wrapped,
        xy=(0.5, -0.15),             # Position relative to your data or axes
        xycoords='axes fraction',    # Placed relative to the axes container
        ha='center',
        va='top',
        annotation_clip=False
    )
    
    ax.set_aspect('equal', anchor='C')

    plt.show()


def determine_sub_polygons(axes, patch_list, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list, num_subpolygons_max):
#def determine_sub_polygons(axes, patch_list, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list, linewidths_list, num_subpolygons_max):

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

    patches_zoom = [patch_list[ii] for ii in range(len(patch_list)) if visible_patch_mask[ii]]
    edgecolors_zoom = [edgecolors_list[ii] for ii in range(len(patch_list)) if visible_patch_mask[ii]]
    anomalies_zoom = [anomaly_var_raw_values_list[ii] for ii in range(len(patch_list)) if visible_patch_mask[ii]]
    #linewidths_zoom = [linewidths_list[ii] for ii in range(len(patch_list)) if visible_patch_mask[ii]]
    profiles_per_patch_zoom = [len(anomaly_list) for anomaly_list in anomalies_zoom]

    mini_patches_list = []
    mini_patches_anomaly_list = []
    mini_patches_edgecolors_list = []
    mini_patches_linewidths_list = []

    for patch_dex in range(len(patches_zoom)):

        num_profiles = profiles_per_patch_zoom[patch_dex]

        '''
        if num_profiles > 10000000:
        #if num_profiles > num_subpolygons_max:

            mini_patches_list.append(patches_zoom[patch_dex])
            mini_patches_anomaly_list.append(np.mean(anomalies_zoom[patch_dex]))
            mini_patches_edgecolors_list.append(edgecolors_zoom[patch_dex])
            mini_patches_linewidths_list.append(linewidths_zoom[patch_dex])

        else:
        '''

        patch_vertices = patches_zoom[patch_dex].get_xy()
        orig_poly = ShapelyPolygon(patch_vertices)

        # STRANGE ISSUE WITH A SINGLE POLYGON NEAR THE EQUATOR - IT HAD ZERO AREA, FROZE THE PROGRAM
        # SETTING A PDB TRACE HERE LETS US LOOK AT "orig_poly.area" and "orig_poly.exterior", WHICH REVEAL TROUBLING ZEROS!!!
        # is this bc of the proximity to the equator??

        # FOR NOW, SWEEPING THIS UNDER THE RUG!!!

        # Note - First I'd run this without the trace, but with the "print(num_profiles)", so I could see the last profile count
        # printed before the program freezes.  Then I'd add the trace, and hit "c" until reaching the problematic polygon.
        # That's where "orig_poly.area" was 0, and original_poly.exterior had zeros for all lat values...

        # OK this definitely has to do with the equator.  Brain tired but should be fixable or maybe fine to skip problematic bins

        # OK, we don't even need to zoom past the threshold to see the bug - bins spanning the equator are wonky, only plotting on one
        # side of it (in the current case, the southern side)

        #print(num_profiles)

        if orig_poly.area == 0:
            #pdb.set_trace()
            print("A polygon with zero area was encountered... did one of the bins you zoomed in on span the equator?  Still haven't sorted this bug out...")
            continue

        # THIS NEEDS TO BE CALCULATED AND STORED IN A FILE THAT'S LOADED, REALLY DUMB TO WASTE TIME DOING THIS EVERY TIME THE AXIS LIMITS CHANGE
        # WHILE ZOOMED IN

        # VIBING OUT
        num_samples = 100000
        #num_samples = 10000
        #num_samples = 5000
        #num_samples = 2000
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
            try:
                polygon_coords = ShapelyCoordinates(ShapelyPolygon(cluster_points).convex_hull)
            except:
                pdb.set_trace()
            mini_patches_list.append(patches.Polygon(polygon_coords, closed=True))

        mini_patches_anomaly_list += anomalies_zoom[patch_dex]
        mini_patches_edgecolors_list += [edgecolors_zoom[patch_dex]] * num_profiles
        if num_profiles == 1:
            mini_patches_linewidths_list.append(0)
        else:
            mini_patches_linewidths_list += [1] * num_profiles

    try:
        subcol = PatchCollection(mini_patches_list, cmap=cmap_face, norm=norm_face, edgecolors=mini_patches_edgecolors_list, transform=ccrs.PlateCarree(), joinstyle='miter')
    except:
        pdb.set_trace()

    subcol.set_array(np.array(mini_patches_anomaly_list)) 

    subcol.set_linewidths(mini_patches_linewidths_list)

    return subcol, mini_patches_anomaly_list




def scale_with_zoom(axes, num_artists_original, colorbar, globe_area, original_linewidths_raw, linewidth_floor, count_array, patch_list, anomaly_var_face_list, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list, scale_threshold, num_subpolygons_max, profiles_lons, profiles_lats):
#def scale_with_zoom(axes, num_artists_original, colorbar, original_area, original_linewidths_raw, linewidth_floor, count_array, patch_list, anomaly_var_face_list, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list, scale_threshold, num_subpolygons_max, profiles_lons, profiles_lats):

    # My hack to remove the previous patch collection, in case we zoom in and would see it underneath the sub-polygons.
    # Note that this means the patch collection must be the last collection added to the axes (ie after calling ax.scatter(), ax.coastlines(), etc)
    #if len(axes.collections) == num_artists_original:
    if len(axes.collections) == num_artists_original + 1:
        axes.collections[-1].remove()
        axes.collections[-1].remove()
    elif len(axes.collections) == num_artists_original:
        axes.collections[-1].remove()

    current_range_x = axes.get_xlim()[1] - axes.get_xlim()[0]  
    current_range_y = axes.get_ylim()[1] - axes.get_ylim()[0]  
    current_area = current_range_x * current_range_y

    scale_factor = np.sqrt(globe_area / current_area)
    #scale_factor = np.sqrt(original_area / current_area)
    #print(f"zoom scale factor: {scale_factor}")

    '''
    current_collection=axes.collections[-1]
    if type(current_collection) is PatchCollection:
        current_collection.remove() # erase the old collection/plot, start fresh.  maybe unecessary, but just want to get this working for now
    '''

    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()

    new_linewidth_max = scale_factor # Random choice, but seems to do the job
    new_linewidths = np.clip(original_linewidths_raw * scale_factor, linewidth_floor, new_linewidth_max)

    if scale_factor < scale_threshold:

        #---------------------------------------------------
        collection = PatchCollection(patch_list, cmap=cmap_face, norm=norm_face, edgecolors=edgecolors_list, transform=ccrs.PlateCarree(), joinstyle='miter')
        collection.set_array(np.array(anomaly_var_face_list))
        axes.add_collection(collection)
        #---------------------------------------------------

        collection.set_linewidths(new_linewidths)

        custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, axes)
        if custom_handles != 0:
            axes.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                      handles=custom_handles, title=f"{legend_title}")

        cbar_min, cbar_max = find_colorbar_limits(colorbar, anomaly_var_face_list)
        colorbar.ax.set_ylim(cbar_min, cbar_max)

    else:
        collection, new_anomaly_var_face_list = determine_sub_polygons(axes, patch_list, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list, num_subpolygons_max)
        #collection, new_anomaly_var_face_list = determine_sub_polygons(axes, patch_list, anomaly_var_raw_values_list, cmap_face, norm_face, edgecolors_list, new_linewidths, num_subpolygons_max)

        if collection != 0:
            custom_handles, legend_title = make_handles_and_titles(patch_list, count_array, axes)
            axes.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                      handles=custom_handles, title=f"{legend_title}")

            axes.add_collection(collection)
            cbar_min, cbar_max = find_colorbar_limits(colorbar, new_anomaly_var_face_list)
            colorbar.ax.set_ylim(cbar_min, cbar_max)

        # Haven't figured out why no profile coords are plotting in the one bottom antarctic bin in 
        # /Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2015_PFL_A__ncei_step_10.nc

        '''
        scatter_lons = []
        scatter_lats = []

        for profile_dex in range(len(profiles_lons)):
            if (profiles_lons[profile_dex] > axes.get_xlim()[0]
                and profiles_lons[profile_dex] < axes.get_xlim()[1] 
                and profiles_lats[profile_dex] > axes.get_ylim()[0] 
                and profiles_lats[profile_dex] < axes.get_ylim()[1]):

                scatter_lons.append(profiles_lons[profile_dex])
                scatter_lats.append(profiles_lats[profile_dex])

        pdb.set_trace()
        axes.scatter(scatter_lons,scatter_lats,c='red',s=1,zorder=10)
        '''

        axes.scatter(profiles_lons, profiles_lats, c='red', s=1, zorder=10)

            

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
