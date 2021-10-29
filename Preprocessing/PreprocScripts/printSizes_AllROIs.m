% Count how many voxels are in each ROI (using nifti files)
% print out tables for each hemisphere and subject.

% MMH 4/15/18
%%
clear 
close all

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
root = mypath(1:filesepinds(end-nDirsUp+1));

sample_dir = fullfile(root,'Samples');
figpath = fullfile(root,'figs');

sublist = [2,3,4,5,6,7];

% the first 14 areas here are the retinotopic and motor ROIs that we
% included in our main analyses (note they're slightly out of order from how they appear
% in figures). 
% the next 6  'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA' are
% areas in the multiple-demand network that we ran some initial analyses on
% and didn't include in main analyses.
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
   'S1','M1','Premotor',...
   'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA'};

nROIs = length(ROI_names);

vis_inds=[1:11];
motor_inds = [12:14];
mdloc_inds=[15:20];
all_inds = [1:5,10,11,6:9,12:14];
hemi_names = {'lh','rh'};

%% print sizes of all visual (retinotopic) ROIs

all_sizes = zeros(nROIs,length(hemi_names), length(sublist));

col_names = [];
for ss=1:length(sublist)
    
    substr = sprintf('S%02d',sublist(ss));
    fn2load = fullfile(sample_dir,sprintf('SampleFile_%s.mat',substr));
    fprintf('loading from %s\n',fn2load);
    load(fn2load, 'ROIs');    
    col_names = [col_names, {[hemi_names{1} '_' substr], [hemi_names{2} '_' substr]}];
    
    for vv=1:nROIs
        for hh=1:length(hemi_names)
            
            
            [rowind1,colind1] = find(strcmp(reshape({ROIs.name},2,nROIs),sprintf('%s_%s',hemi_names{hh},ROI_names{vv})));
            if ~isempty(rowind1)
                col_inds = [colind1]; % column is the region
                row_inds = [rowind1];   % row is the hemisphere

                all_sizes(vv,hh,ss) = numel(ROIs(row_inds,col_inds).voxel_inds);
            end
               
        end
        
    end
end

%% print sizes of ROIs in the order reported in paper

allareas_sizes = all_sizes(all_inds,:,:);
my_areas_all = ROI_names(all_inds);
fprintf('\nALL ROI SIZES:\n\n');   
tab = array2table(reshape(allareas_sizes, length(my_areas_all), length(hemi_names)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas_all);
disp(tab)

writetable(tab, fullfile(figpath, 'ROI_size_table.csv'), 'WriteRowNames', true)

%% print sizes of motor (digit localizer) ROIs

motor_sizes = all_sizes(motor_inds,:,:);
my_areas_motor = ROI_names(motor_inds);
fprintf('\nMOTOR ROI SIZES:\n\n');
tab = array2table(reshape(motor_sizes, length(my_areas_motor), length(hemi_names)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas_motor);
disp(tab)


%% print sizes of retinotopic ROIs
          
visual_sizes = all_sizes(vis_inds,:,:);
my_areas_vis = ROI_names(vis_inds);
fprintf('\nRETINOTOPIC ROI SIZES:\n\n');   
tab = array2table(reshape(visual_sizes, length(my_areas_vis), length(hemi_names)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas_vis);
disp(tab)


%% print sizes of motor (digit localizer) ROIs

motor_sizes = all_sizes(motor_inds,:,:);
my_areas_motor = ROI_names(motor_inds);
fprintf('\nMOTOR ROI SIZES:\n\n');
tab = array2table(reshape(motor_sizes, length(my_areas_motor), length(hemi_names)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas_motor);
disp(tab)

%% print sizes of MD localizer ROIs

mdloc_sizes = all_sizes(mdloc_inds,:,:);
my_areas_mdloc = ROI_names(mdloc_inds);
fprintf('\nMULTIPLE-DEMAND ROI SIZES:\n\n');
tab = array2table(reshape(mdloc_sizes, length(my_areas_mdloc), length(hemi_names)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas_mdloc);
disp(tab)
