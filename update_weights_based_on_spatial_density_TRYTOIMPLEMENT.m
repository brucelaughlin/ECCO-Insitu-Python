function update_weights_based_on_spatial_density(run_code)

% function update_weights_based_on_spatial_density(run_code)
%
% function update_weights_based_on_spatial_density(run_code)
% This script calcuates scaling factors for profile weights
% such that when n profiles appear within a particular geodesic bin
% within a particular month their weights are scaled by 1/n.
% 
% Consequently, the forcing of the adjoint model from the 
% in-situ profile model-data misfit will be relatively reduced where
% multiple profiles appear within the same geodesic bin and month.
%
% For example, if during one particular month there are 100 profiles in a
% bin, each profile will contribute 1/100 of its former forcing to the
% adjoint model.  


% Filename; update_zero_weight_points_on_prepared_profiles.m

%  ** former filename :
% Date Created: 2018-11-25
% Last Modified: 2020-06-30

% notes:
% 2020-06-30 added run code for 2020-06-30
%            changed out netcdf files are read using new matlab syntax
%            changed write routine to write_profile_structure_to_netcdf_2020
% 2018-12-21 added figure making
% 2018-01-16 fixed bug where it wasn't writing the last year
%            also added SH and NH figures for obs count
%            also cleaned up how to profile scalings are applied (at the
%            very end of the code)  shouldn't change any results.  may 
%            speed things up in some limited cases.
              
close all;

% defaults
make_figs = 1
fig_dir = [];
file_glob = '*nc'

%%
plotting_resolution = .25 % degrees
fprintf('global \n')
grs = prep_for_geodesic_bins_to_latlon_mapping(10242, -1, ...
        plotting_resolution, 0);

fprintf('NH \n')
grs_n = prep_for_geodesic_bins_to_latlon_mapping(10242, -2, ...
    plotting_resolution, 0);

fprintf('SH \n')
grs_s = prep_for_geodesic_bins_to_latlon_mapping(10242, -3, ...
    plotting_resolution, 0);

fprintf('loading geo bins 10242 \n')
[lats_g_10k lons_g_10k] = load_geo_bins(10242);

clear scaling;

switch run_code
     case '20200630_all'
    
        run_code
        scaling_years = 1992:2019
        bin_code = 10242
        
        % scaling factor output directory
        scaling_factor_output_dir = '/home/ifenty/data/observations/insitu/scaling_factors'
        
        scaling_factor_output_name = ...
            [run_code '_geodesic_10242_scalar_weight_factors_all_profs_' ...
            num2str(scaling_years(1)) '_' num2str(scaling_years(end)) '.mat']
        
        new_profile_suffix = '_sp_scal.nc';
        profile_output_dir = '/home/ifenty/data11/ifenty/observations/insitu/20200630/llc_90/post-spatial-scaling/'
        mkdir(profile_output_dir);
        
        prof_dir = '/home/ifenty/data/observations/insitu/compilation/llc90/20200630/pre-spatial-scaling/'
        
        apply_scaling_factors_to_weight_fields = 1
        
        make_figs =1 
        fig_dir = [profile_output_dir '/figures']
        mkdir(fig_dir)
        
        file_glob = '*nc'
        
     case '20190220_all'
    
        run_code
        scaling_years = 1992:2017
        bin_code = 10242
        
        %% scaling factor output directory
        scaling_factor_output_dir = '/home/ifenty/data/observations/insitu/scaling_factors'
        
        scaling_factor_output_name = ...
            [run_code '_geodesic_10242_scalar_weight_factors_all_profs_' ...
            num2str(scaling_years(1)) '_' num2str(scaling_years(end)) '.mat']
        
        new_profile_suffix = '_sp_scal.nc';
        profile_output_dir = '/home/ifenty/data11/ifenty/observations/insitu/20190220/llc_90/post-spatial-scaling/'
        mkdir(profile_output_dir);
        
        prof_dir = '/home/ifenty/data/observations/insitu/compilation/llc90/20190220/pre-spatial-scaling/'
        
        apply_scaling_factors_to_weight_fields = 1
        
        make_figs =1 
        fig_dir = [profile_output_dir '/figures']
        mkdir(fig_dir)
        
        file_glob = '*nc'
        
        
    case '20190131_all'
    
        run_code
        scaling_years = 1992:2017
        bin_code = 10242
        
        %% scaling factor output directory
        scaling_factor_output_dir = '/home/ifenty/data/observations/insitu/scaling_factors'
        
        scaling_factor_output_name = ...
            [run_code '_geodesic_10242_scalar_weight_factors_all_profs_' ...
            num2str(scaling_years(1)) '_' num2str(scaling_years(end)) '.mat']
        
        new_profile_suffix = '_sp_scal.nc';
        profile_output_dir = '/home/ifenty/data11/ifenty/observations/insitu/20190131/llc90/post-spatial-scaling/'
        mkdir(profile_output_dir);
        
        prof_dir = '/home/ifenty/data/observations/insitu/compilation/llc90/20190131/pre-spatial-scaling/'
        
        apply_scaling_factors_to_weight_fields = 1
        
        make_figs =1 
        fig_dir = [profile_output_dir '/figures']
        mkdir(fig_dir)
        
        file_glob = '*nc'
        
   
end

mkdir (profile_output_dir)
mkdir(fig_dir)

num_scaling_years = length(scaling_years);

year_range = [num2str(scaling_years(1)) '-' num2str(scaling_years(end))]

%% initialize arrays

prof_lat_all = []
prof_lon_all = []
prof_YYYYMMDD_all = []
prof_bin_id_a_all = []


% loop through all input directories and aggregate stats
cd(prof_dir)
files = dir([file_glob])
num_files = length(files)

clear lats lons dates prof_a prof_b
files_all = {};
for f = 1:num_files
    
    fprintf('filename %s \n', files(f).name)
    
    files_all = {files_all{:} files(f).name};
    
    prof_lat  = ncread(files(f).name,'prof_lat');
    prof_lon  = ncread(files(f).name,'prof_lon');
    prof_YYYYMMDD  = ncread(files(f).name,'prof_YYYYMMDD');
    prof_bin_id_a  = ncread(files(f).name,'prof_bin_id_a');
    %%
    prof_lat_f{f} = prof_lat;
    prof_lon_f{f} = prof_lon;
    prof_YYYYMMDD_f{f} = prof_YYYYMMDD;
    prof_bin_id_a_f{f} = prof_bin_id_a;
    
    %%
    prof_lat_all = [prof_lat_all(:)' prof_lat(:)'];
    prof_lon_all = [prof_lon_all(:)' prof_lon(:)'];
    prof_YYYYMMDD_all = [prof_YYYYMMDD_all(:)' prof_YYYYMMDD(:)'];
    prof_bin_id_a_all = [prof_bin_id_a_all(:)' prof_bin_id_a(:)'];
    
    clear prof_lat prof_lon prof_YYYYMMDD prof_bin_id_a prof_bin_id_b;
end



%%
% Count in 1x1 degree bins for reference
[lon_mat, lat_mat] = meshgrid(-180:180,-90:90);

min_lon = -180;
min_lat = -90;

tmp_lat = round(prof_lat_all);
tmp_lon = round(prof_lon_all);

tmp_lat_i = tmp_lat - min_lat + 1;
tmp_lon_i = tmp_lon - min_lon +1;

prof_count = zeros(size(lon_mat));

%%     
tmp = floor(length(tmp_lat_i)/100);
for i = 1:length(tmp_lat_i)
    prof_count(tmp_lat_i(i), tmp_lon_i(i))=prof_count(tmp_lat_i(i), tmp_lon_i(i)) + 1;
    if mod(i,tmp)==0
        fprintf('progress: %f pct" \n', 100*i/length(tmp_lat_i))
    end
end
fprintf('finished counting profiles\n')
    
%%
if make_figs
    %%
    figure(10);clf
    pltstmp(gcf, 0, [run_code])

    subplot(211);
    imagesc(prof_count);axis xy;
    c=colormap(jet(1000));c(1,:) = 0;colormap(c)
    %caxis([0 10.^(2*median(log10(prof_count(:))))])
    caxis([0 2*median((prof_count(prof_count > 0)))]);
    colorbar;
    title('count in 1x1 degree boxes')
    
    subplot(212);
    imagesc(log10(prof_count));axis xy;
    c=colormap(jet(1000));c(1,:) = 0;colormap(c)
    %caxis([0 (2*median(log10(prof_count(:))))])
    caxis auto
    colorbar;
    title('log10 count in 1x1 degree boxes')
    
    if length(fig_dir) > 0
        paperx=16;papery=16;prep_figure_for_print;
        fig_name = ['profile_count_' run_code '-' year_range '.png']
        fig_name = [fig_dir '/' fig_name]
        print(gcf,'-dpng','-r300',fig_name);
    end
end

%%    
count_10k = zeros(length(lats_g_10k), length(scaling_years),12);

fprintf('counting profiles in each month/year \n')
for y = 1:length(scaling_years)
    for m = 1:12
            
        start_date = scaling_years(y)*1e4 +     m*1e2;
        end_date   = scaling_years(y)*1e4 + (m+1)*1e2;
            
        fprintf('start & end date %10.0f %10.0f \n', start_date, end_date)
            
        ins = find(prof_YYYYMMDD_all  >= start_date & ...
                   prof_YYYYMMDD_all  < end_date);
            %%
        if length(ins) > 0
            [a,b]=hist(prof_bin_id_a_all(ins),[.5:10241.5]);
            
            count_10k(:,y,m) = a';
        end
    end
end

%%

% an array with the # of profiles in each bin for each month over the 
% number of scaling yeras
count_10k_long = reshape(count_10k,[10242, length(scaling_years)*12]);

% take sum in the bin dimension.  result is # of profiles total in each
% bin over the scaling period.
count_10k_sum = sum(count_10k_long,2);

% ins where there are no profiles
ins_no_profiles = find(count_10k == 0);

% bins with no profiles ever
ins_never_profiles = find(count_10k_sum==0);

% nan the count where we never have profiles
count_10k_sum(ins_never_profiles) = NaN;

['finished count']


%% WEIGHT FACTOR

% determine weight factor, alpha(i,t) = 1/n(i,t)  
% where n = num of profiles in geodesic bin 'i' for each simulation month 't'
%    where 'i' is bin id

scaling_factor = zeros(size(count_10k));

% ins where at least one profile
ins2 = find(count_10k >= 1);
% ins where no profiles
ins3 = find(count_10k == 0);

% the scaling factor is 1/# of profiles in that month for each bin
scaling_factor(ins2) = 1./count_10k(ins2);

% set factor to be 1 where there are no profiles
scaling_factor(ins_no_profiles) = 1;

% 
% % add a little more scaling factor the argentine basin 
% switch run_code
%      case '20200630_all'
%          
%         % distance from center of argentine basine
%         d_from_AP = ((-45 - lats_g_10k).^2 + (-50 - lons_g_10k).^2).^0.5;
% 
%         scaled_d = d_from_AP/70 - 0.25;
%         regional_factor = 2-1 ./ (1 + exp(-15*scaled_d));
%         figure(99);clf;
%         reverse_map_grid_to_latlon(grs.lon, grs.lat, grs.latlon_bins,...
%             regional_factor,grs.latlon_land_ins,2,0,10);;caxis([0 2]);;colorbar
% 
%         
%         orig_scaling_factor = scaling_factor;
%         
%         
%         for y = 1:length(scaling_years)
%             for m = 1:12
% 
%                 scaling_factor(:,y,m) = scaling_factor(:,y,m).*regional_factor;
%             end
%         end
%         fprintf(['finished applying regional scaling factor\n'])      
%     otherwise
%         fprintf(['no regional scaling factor applied\n'])
% end
  


% create scaling factor data structure.  save for later!
scaling.geo_bin = bin_code;
scaling.scaling_factor = scaling_factor;
scaling.files = files_all;
scaling.scaling_years = scaling_years;
scaling.calculation_date = date;

cd(scaling_factor_output_dir)
save(scaling_factor_output_name, 'scaling')
fprintf('finished calculating and saving scaling factor\n')



%%




%% plot count
if make_figs
    c=colormap(jet(1000));
    c(1,:)=0;
    figure(30);clf;
    pltstmp(gcf, 0, [run_code])    
    [lon, lat, field]=reverse_map_grid_to_latlon(grs.lon, grs.lat, ...
        grs.latlon_bins,count_10k_sum, grs.latlon_land_ins, 2,0,10);
    caxis([0 10*median(count_10k_sum(count_10k_sum>0))])
    colormap(c);
    colorbar
    title(['count in 10242 bins ' year_range])
    set(gcf, 'Position',[0 0 2000 1000])

    if length(fig_dir) > 0
        paperx=20;papery=10;prep_figure_for_print;
        fig_name = ['profile_count_10242_bins_' run_code '-' year_range '_GLOBAL.png']
        fig_name = [fig_dir '/' fig_name]
        print(gcf,'-dpng','-r100',fig_name);
    end
    %%
    
    figure(31);clf;
    pltstmp(gcf, 0, [run_code])    
    [lon, lat, field]=reverse_map_grid_to_latlon(grs.lon, grs.lat, ...
        grs.latlon_bins,count_10k_sum, grs.latlon_land_ins, 3.1,0,10);
    caxis([0 10*median(count_10k_sum(count_10k_sum>0))])
    colormap(c);
    colorbar
    title(['count in 10242 bins ' year_range])    
    set(gcf, 'Position',[0 0 1000 1000])

    if length(fig_dir) > 0
        paperx=10;papery=10;prep_figure_for_print;
        fig_name = ['profile_count_10242_bins_' run_code '-' year_range '_SH.png']
        fig_name = [fig_dir '/' fig_name]
        print(gcf,'-dpng','-r300',fig_name);
    end
    %%
    close (32)
    figure(32);clf;
    pltstmp(gcf, 0, [run_code])    
    [lon, lat, field]=reverse_map_grid_to_latlon(grs.lon, grs.lat, ...
        grs.latlon_bins,count_10k_sum, grs.latlon_land_ins, 3,0,10);
    caxis([0 10*median(count_10k_sum(count_10k_sum>0))])
    colormap(c);
    colorbar
    title(['count in 10242 bins ' year_range])
    set(gcf, 'Position',[0 0 1000 1000])

    if length(fig_dir) > 0
        paperx=10;papery=10;prep_figure_for_print;
        fig_name = ['profile_count_10242_bins_' run_code '-' year_range '_NH.png']
        fig_name = [fig_dir '/' fig_name]
        print(gcf,'-dpng','-r100',fig_name);
    end
    
    %%
    
    close(33)
    figure(33);clf;
    [lon, lat, field]=reverse_map_grid_to_latlon(grs.lon, grs.lat, ...
        grs.latlon_bins,log10(count_10k_sum), grs.latlon_land_ins, 2,0,10);
    cmax = ceil(max(log10(count_10k_sum)))
    caxis([1 min(cmax,4)]);
    colormap(c);
    colorbar
    set(gcf, 'Position',[0 0 2000 1000])
    set(gcf,'color','w');
    title(['log10(number of profiles in 10242 bins ' year_range])
    if length(fig_dir) > 0
        paperx=20;papery=10;prep_figure_for_print;
        fig_name = ['profile_count_10242_bins_log_10_' run_code '-' year_range '.png']
        fig_name = [fig_dir '/' fig_name]
        print(gcf,'-dpng','-r100',fig_name);
    end
   
end


%% plot probability of having one or more profiles in a given month

count_10k_long = reshape(count_10k,[10242, length(scaling_years)*12]);
count_10k_p = sum(count_10k_long >= 1,2)./size(count_10k_long,2);

count_10k_p(ins_never_profiles) = NaN;

if make_figs
    
    close(31)
    figure(31);clf;
    [lon, lat, field]=reverse_map_grid_to_latlon(grs.lon, grs.lat, ...
        grs.latlon_bins,count_10k_p, grs.latlon_land_ins, 2,0,10);
    caxis([0 1]);
    colormap(c);
    colorbar
    set(gcf, 'Position',[0 0 2000 1000])

    title(['probability of one or more profiles in simulation month ' year_range])
    if length(fig_dir) > 0
        paperx=20;papery=10;prep_figure_for_print;
        fig_name = ['probability_of_one_or_more_profs_in_month_' run_code '-' year_range '.png']
        fig_name = [fig_dir '/' fig_name]
        print(gcf,'-dpng','-r100',fig_name);
    end
    
    figure(32);clf;
    pltstmp(gcf, 0, [run_code])    
    [lon, lat, field]=reverse_map_grid_to_latlon(grs.lon, grs.lat, ...
        grs.latlon_bins,count_10k_p, grs.latlon_land_ins, 3.1,0,10);
    caxis([0 1]);
    colormap(c);
    colorbar
    set(gcf, 'Position',[0 0 1000 1000])

    title(['probability of one or more profiles in simulation month ' year_range])
    if length(fig_dir) > 0
        paperx=10;papery=10;prep_figure_for_print;
        fig_name = ['probability_of_one_or_more_profs_in_month_' run_code '-' year_range '_SH.png']
        fig_name = [fig_dir '/' fig_name]
        print(gcf,'-dpng','-r100',fig_name);
    end
    
    figure(33);clf;
    pltstmp(gcf, 0, [run_code])    
    [lon, lat, field]=reverse_map_grid_to_latlon(grs.lon, grs.lat, ...
        grs.latlon_bins,count_10k_p, grs.latlon_land_ins, 3,0,10);
    caxis([0 1]);
    colormap(c);
    colormap(c);
    colorbar
    title(['probability of one or more profiles in simulation month ' year_range])
    set(gcf, 'Position',[0 0 1000 1000])

    title(['probability of one or more profiles in simulation month ' year_range])
    if length(fig_dir) > 0
        paperx=10;papery=10;prep_figure_for_print;
        fig_name = ['probability_of_one_or_more_profs_in_month_' run_code '-' year_range '_NH.png']
        fig_name = [fig_dir '/' fig_name]
        print(gcf,'-dpng','-r100',fig_name);
    end
    
end


%% plot mean scaling factors

if make_figs
    %%
    for i = 0:1
        remove_missing_data = i
        tmp=reshape(scaling_factor,[10242, length(scaling_years)*12]);

        if remove_missing_data
            tmp(ins_no_profiles) = NaN;
            suffix = 'remove missing'
        else
            suffix = 'keep missing'
        end

        tmp2=mynanmean(tmp,2);

        tmp2(ins_never_profiles) = 0;

        figure(4+i*100);clf;
        [lon, lat, field]=reverse_map_grid_to_latlon(grs.lon, grs.lat, ...
            grs.latlon_bins,tmp2, grs.latlon_land_ins, 2,0,10);
        caxis([0 1]);
        colormap(c);colorbar
        title(['mean scaling factor ' suffix])        
        if length(fig_dir) > 0
            paperx=16;papery=16;prep_figure_for_print;
            fig_name = ['mean_scaling_factor_2_' run_code '_' suffix '-' year_range '.png']
            fig_name = [fig_dir '/' fig_name]
            print(gcf,'-dpng','-r300',fig_name);
        end

        figure(5+i*100);clf;
        [lon, lat, field]=reverse_map_grid_to_latlon(grs_n.lon, grs_n.lat, ...
            grs_n.latlon_bins,tmp2, grs_n.latlon_land_ins, 3,0,10);
        caxis([0 1]);
        colormap(c);colorbar
        title(['mean scaling factor ' suffix])
        if length(fig_dir) > 0
            paperx=16;papery=16;prep_figure_for_print;
            fig_name = ['mean_scaling_factor_3_' run_code  '_' suffix '-' year_range '.png']
            fig_name = [fig_dir '/' fig_name]
            print(gcf,'-dpng','-r300',fig_name);
        end
        
        figure(6+i*100);clf;
        [lon, lat, field]=reverse_map_grid_to_latlon(grs_s.lon, grs_s.lat, ...
            grs_s.latlon_bins,tmp2, grs_s.latlon_land_ins, 3.1,0,10);
        caxis([0 1]);
        colormap(c);
        colorbar
        title(['mean scaling factor ' suffix])        
        if length(fig_dir) > 0
            paperx=16;papery=16;prep_figure_for_print;
            fig_name = ['mean_scaling_factor_31_' run_code  '_' suffix '-' year_range '.png']
            fig_name = [fig_dir '/' fig_name]
            print(gcf,'-dpng','-r300',fig_name);
        end
    end
end % make figs

%%
if apply_scaling_factors_to_weight_fields
    
    scaling_scaling_years = scaling_years(1):scaling_years(end);
    fprintf('Scaling years : %10.f, %10.f \n', scaling_years(1), scaling_years(end))
    
    
    %% APPLY SCALING FACTORS to weight fields
    for i = 1:length(files_all)
        %%
        cd(prof_dir)
        mitprof = make_MITprof_structure_from_MITprof_netcdf_file(files_all{i});
        %mitprof= MITprof_read(files_all{i});
        %num_profs = length(mitprof.prof_lon);
        
        %%
        prof_year = floor(mitprof.prof_YYYYMMDD./1e4);
        prof_mon  = floor((mitprof.prof_YYYYMMDD - prof_year*1e4)./1e2);
        prof_year_i = prof_year - scaling.scaling_years(1) + 1;
        
        %%
        mitprof_new = mitprof;
        num_profs = length(mitprof.prof_date);
        % set default spatial scaling factor to -9999
        mitprof_new.prof_spatial_scaling_factor = ones(num_profs,1).*(-9999);

        fprintf('\nApplying weight factor to : %s \t %i profiles \n\n', files_all{i}, num_profs)

        process_profiles = 0;
        
        if length(unique(prof_year_i))==1
            if (prof_year(1) >= scaling_years(1) & prof_year(1) <= scaling_years(end))
                fprintf('all the profiles of this file have the same year\n')
                fprintf('and we are applying spatial scale factors on them\n')
                process_profiles = 1;
            else
                fprintf('all the profiles of this file have the same year\n')
                fprintf('but we are not applying spatial scale factors on them\n')
                process_profiles = 0;
            end
        else
            tmp = intersect(prof_year, scaling_scaling_years);
            if length(tmp) > 0
                fprintf('profiles in this file have different scaling_years\n')
                fprintf('and we are applying spatial scale factors to some of them\n')
                process_profiles = 1;
            else
                fprintf('profiles in this file have different scaling_years\n')
                fprintf('and we are applying spatial scale factors to NONE of them\n')
                process_profiles = 0;
            end
        end
        
        
        if process_profiles

            % loop through all profiles in this file
            for p = 1:num_profs

                % pull out its bin
                bin = mitprof.prof_bin_id_a(p);

                % verify that it falls within the scaling_years that we care about
                if (prof_year(p) >= scaling_years(1) & prof_year(p) <= scaling_years(end))
                    sf = scaling.scaling_factor(bin, prof_year_i(p), prof_mon(p));
                else
                    % if not set its scaling factor to ZERO
                    %['we are not looking at this particular year! - setting sf = 0']
                    %files_all{i};
                    %prof_year(p);
                    %prof_year_i(p);
                    %['valid scaling_years']
                    %[scaling_years(1) scaling_years(end)];
                    sf = 0.0;
                end

                % one way or another, apply the scaling factor to profile p.
                mitprof_new.prof_spatial_scaling_factor(p) = sf;
                mitprof_new.prof_Tweight(p,:) = mitprof.prof_Tweight(p,:).*sf;
                mitprof_new.prof_Sweight(p,:) = mitprof.prof_Sweight(p,:).*sf;

                % output intermediate progress
                if mod(p,1000)==0
                    fprintf('progress : %3.0f \n', (100*p/num_profs))
                end
            end

            % stop if missing scaling factor
            if length(find(mitprof_new.prof_spatial_scaling_factor == -9999 >0))
                fprintf('found -9999 spatial scaling factor! \n')
                stop
                
            end


            % check to see if we have at leats one nonzero good weight
            % if so, then we save.  otherwise we skip the file.
            if sum(mitprof_new.prof_Tweight(:) + mitprof_new.prof_Sweight(:)) > 0

                if length(mitprof_new.prof_YYYYMMDD) > 0

                    scaling_factor_output_name = [files_all{i}(1:end-3) new_profile_suffix];


                    fileOut=[profile_output_dir '/' scaling_factor_output_name];
                    fprintf('\n... writing output in %s\n',scaling_factor_output_name);

                    write_profile_structure_to_netcdf_2020(mitprof_new,fileOut);
                end
            else
                fprintf('\n---> not a single nonzero weight remains after spatial scaling\n')
            end
        else
            fprintf('not processing this file : %s \n ', files_all{i})
        end
            
    end
end
%%
% 
% 
%  case '20190128_all'
%     
%         run_code
%         scaling_years = 1992:2017
%         bin_code = 10242
%         
%         %% scaling factor output directory
%         scaling_factor_output_dir = '/home/ifenty/data/observations/insitu/scaling_factors'
%         
%         scaling_factor_output_name = ...
%             [run_code '_geodesic_10242_scalar_weight_factors_all_profs_' ...
%             num2str(scaling_years(1)) '_' num2str(scaling_years(end)) '.mat']
%         
%         new_profile_suffix = '_sp_scal.nc';
%         profile_output_dir = '/home/ifenty/data11/ifenty/observations/insitu/20190128/llc90/post-spatial-scaling'
%         mkdir(profile_output_dir);
%         
%         prof_dir = '/home/ifenty/data/observations/insitu/compilation/llc90/20190128/'
%         
%         apply_scaling_factors_to_weight_fields = 1
%         
%         make_figs =1 
%         fig_dir = [profile_output_dir '/figures']
%         mkdir(fig_dir)
%         
%         file_glob = '*'
%         
%       case '20181221_ctd'
%     
%         run_code
%         scaling_years = 1992:2017
%         bin_code = 10242
%         
%         %% scaling factor output directory
%         scaling_factor_output_dir = '/home/ifenty/data/observations/insitu/scaling_factors'
%         
%         scaling_factor_output_name = '20181221_geodesic_10242_scalar_weight_factors_all_profs_1992_2017.mat'
%         
%         new_profile_suffix = '_sp_scal.nc';
%         profile_output_dir = '/home/ifenty/data11/ifenty/observations/insitu/20181221/llc90/post-spatial-scaling'
%         mkdir(profile_output_dir);
%         
%         %% put all profiles together in one directory
%         prof_dir = '/home/ifenty/data/observations/insitu/compilation/llc90/20181221/pre-spatial-scaling/'
%         
%         apply_scaling_factors_to_weight_fields = 1
%         
%         make_figs =1 
%         fig_dir = [profile_output_dir '/figures']
%         mkdir(fig_dir)
%         
%         file_glob = 'CTD*'
%         
%      case '20181221_all'
%     
%         run_code
%         scaling_years = 1992:2017
%         bin_code = 10242
%         
%         %% scaling factor output directory
%         scaling_factor_output_dir = '/home/ifenty/data/observations/insitu/scaling_factors'
%         
%         scaling_factor_output_name = ...
%             [run_code '_geodesic_10242_scalar_weight_factors_all_profs_' ...
%             num2str(scaling_years(1)) '_' num2str(scaling_years(end)) '.mat']
%         
%         new_profile_suffix = '_sp_scal.nc';
%         profile_output_dir = '/home/ifenty/data11/ifenty/observations/insitu/20181221/llc90/post-spatial-scaling'
%         mkdir(profile_output_dir);
%         
%         prof_dir = '/home/ifenty/data/observations/insitu/compilation/llc90/20181221/pre-spatial-scaling/'
%         
%         apply_scaling_factors_to_weight_fields = 1
%         
%         make_figs =1 
%         fig_dir = [profile_output_dir '/figures']
%         mkdir(fig_dir)
%         
%         file_glob = '*'
%  
%         
%     case '20190116_2015_2016'
%         % because my code had a bug where it wasn't writing the output of
%         % the LAST year.  
%         scaling_years = 2015:2016
%         
%         bin_code = 10242
%         
%         %% scaling factor output directory
%         scaling_factor_output_dir = '/home/ifenty/data/observations/insitu/scaling_factors'
%         
%         scaling_factor_output_name = ...
%             [run_code '_geodesic_10242_scalar_weight_factors_all_profs_' ...
%             num2str(scaling_years(1)) '_' num2str(scaling_years(end)) '.mat']
%         
%         new_profile_suffix = '_sp_scal.nc';
%         profile_output_dir = '/home/ifenty/data11/ifenty/observations/insitu/20190116b/llc90/post-spatial-scaling'
%         mkdir(profile_output_dir);
%         
%         prof_dir = '/home/ifenty/data/observations/insitu/compilation/llc90/20181221/pre-spatial-scaling/tmp/'
%         
%         apply_scaling_factors_to_weight_fields = 1
%         
%         make_figs = 0
%         fig_dir = [profile_output_dir '/figures']
%         
%         file_glob = '*'
%            