% Loop through the trees and employ compute_treemetrics.m to extract tree 
% metrics for each tree.
%
% Input:   - rootResults   % Root folder where to store the results [character array]
%          - h_butt        % Stump height
%          - h_interval    % diameter measuring interval
%          - h_layer       % point-cloud layer height for diameter measurements
%          - min_points    % min number of points for circle fitting
%          - p             % stem curve smoothing parameter (cubic spline)
%
% Output:
%         - treeMetrics   % Table of tree metrics 
%
% (c) Tuomas Yrttimaa / Science4Trees @ UEF School of Forest Sciences 2021
% ---------------------------------------------------------------------------
% 

function [treeMetrics,crownMetrics] = pc_treeMetrics(rootResults,h_butt,h_interval,h_layer,min_points,p,computestemmetrics,computecrownmetrics,detectbranches)
    hh = waitbar(0,''); hh.Name = 'Point Cloud Tools: Compute tree metrics';
    waitbar(15/100,hh,'Computing tree metrics, please wait...');
    starttime = datetime;
    filenames = struct2table(dir(fullfile([rootResults,'/pointclouds/stem'], '**', '*.ply')));
    filenames = filenames.name;
    numtrees = size(filenames,1);
    tree_ids = zeros([1,size(filenames,1)]);
    for i = 1:numtrees
        tmp = filenames{i};
        tmp = split(tmp,'_');
        tmp = split(tmp(1),'tree');
        tree_ids(i) = str2double(tmp{2});
    end
   
    treeMetrics = array2table(zeros([max(tree_ids) 11]));
    treeMetrics.Properties.VariableNames = {'tree_id','X','Y','dbh','ba','h','vol','d6','d05','vlog','vpulp'};
    crownMetrics = array2table(zeros([max(tree_ids) 9]));
    crownMetrics.Properties.VariableNames = {'tree_id','X','Y','vol_convhull','vol_voxel','area2d','area3d','perimeter','crownBaseHeight'};


    parfor tree_id = 1:max(tree_ids)
        try
            if computestemmetrics == true
                treemetr = compute_treemetrics(tree_id,rootResults,h_butt,h_interval,h_layer,min_points,p);
                treeMetrics(tree_id,:) = treemetr;
            end
            if computecrownmetrics == true
                crownmetr = compute_crownmetrics(tree_id,rootResults);
                crownMetrics(tree_id,:) = crownmetr;
            end
            if detectbranches == true
                pcStem = pcread(strcat(rootResults,sprintf('/pointclouds/stem/tree%i_stempoints.ply',tree_id)));
                pcCanopy = pcread(strcat(rootResults,sprintf('/pointclouds/canopy/tree%i_canopypoints_tls.ply',tree_id)));
                pc = pointCloud(cat(1,pcStem.Location,pcCanopy.Location));
                [~,~,~,~]=pc_detect_branches(tree_id,pc,rootResults);
            end
        catch
            warning('Cannot compute metrics for tree no. %i!', tree_id);
            continue
        end
        
    end
    treeMetrics(treeMetrics.dbh == 0,:) = [];
    crownMetrics(crownMetrics.vol_convhull == 0,:) = [];

    IDX = rangesearch(table2array(treeMetrics(:,2:3)),table2array(treeMetrics(:,2:3)),0.3);
    dbl = [];
    for i = 1:size(IDX,1)
        id = IDX{i};
        if length(id) > 1
            if ~ismember(id(1),dbl)
                dbl = cat(2,dbl,id(2:end));
            end
        end
    end
    treeMetrics(dbl,:) = [];

    
    savedir = strcat(rootResults,'/stemmetrics');
    if isfolder(savedir) == false; mkdir(savedir);end
    writetable(treeMetrics,strcat(savedir,'/stemmetrics.csv'),'Delimiter',';');

    savedir = strcat(rootResults,'/crownmetrics');
    if isfolder(savedir) == false; mkdir(savedir);end
    writetable(crownMetrics,strcat(savedir,'/crownmetrics.csv'),'Delimiter',';');

    endtime = datetime;
    cprintf('comment',sprintf('   %i Trees extracted from TLS data in %s\n',height(treeMetrics),endtime-starttime));
    close(hh);
end
