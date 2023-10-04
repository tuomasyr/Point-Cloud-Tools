% POINT CLOUD PROCESSING TOOLS
%  1. Import point cloud data to matlab     (pc_import.m)
%  2. Detect tree crown segments            (pc_detect_trees_crowns_v2.m)
%  3. Extract tree-wise point clouds        (pc_tree_extraction.m)
%  4. Detect trees and classify point cloud (pc_classify_v2.m)
%  5. Compute tree metrics                  (pc_treeMetrics.m, compute_treemetrics.m, compute_crownmetrics.m)
%  6. Detect and segment branches           (pc_detect_branches.m, run within prev. step) 
%
%  See the related README-document and function descriptions for more details and requirements
%  for the data structure etc.
%
% (c) Tuomas Yrttimaa, School of Forest Sciences, University of Eastern Finland 2021-2023
% Contact: tuomas.yrttimaa@uef.fi
%------------------------------------------------------------------------

    %  Some input parameters:

    % Root folder for point cloud data and the results.
    rootTLSData = 'enter_your_directory';
    rootResults = 'enter_your_directory';
    % Filenames
    filenameTLS = 'filename.las';

    % create parallel pool
    parpool('local',8)

    %%% 1. IMPORT POINT CLOUD DATA  
    [PCtls, plotExtent] = pc_import(rootTLSData,filenameTLS);
      %  Visualize 
      %   PCtlssample = pcdownsample(PCtls,'gridAverage',0.1);   % Downsampling (10 cm spacing)
      %   PCtlssample = pcdenoise(PCtlssample,'NumNeighbor',20); % Remove noise
      %   figure, pcshow(PCtlssample);


    %%% 2. DETECT TREE CROWN SEGMENTS 
        chmPixelSize = 0.2;          % raster cell size in meters
        minCrownArea = 2;            % min crown area in m^2
        viewFigures = false;         % view figures
        minTreeHeight = 2;

        treePgons = pc_detect_tree_crowns_v2(PCtls,chmPixelSize,minCrownArea,minTreeHeight,viewFigures,rootResults);

    %%% 3. EXTRACT TREE-WISE POINT CLOUDS
        pc_tree_extraction(rootTLSData,filenameTLS);
        
        %clearvars PCtls

    %%% 4. DETECT TREES and CLASSIFY STEM POINTS and CROWN POINTS
      % Saving the classified point clouds into Results folder
        pc_classify_v2(rootTLSData,rootResults);
        rmdir([rootTLSData,'/tmp'],'s');

      % Construct merged and classified point cloud.
        % Class 3 = stem points, class 4 = canopy points
        % TreeID stored as user data in las-files
        lastools(['txt2las -i "',rootResults,'/pointclouds/stem/*.ply" -merged -parse xyzu -set_scale 0.0001 0.0001 0.0001 -set_classification 3 -o stempoints.las']); 
        lastools(['txt2las -i "',rootResults,'/pointclouds/canopy/*.ply" -merged -parse xyzu -set_scale 0.0001 0.0001 0.0001 -set_classification 4 -o canopypoints.las']);

        lastools(['lasmerge -i stempoints.las canopypoints.las  -odir "',rootTLSData,'" -o classified.laz']);
        delete('stempoints.las');delete('canopypoints.las');

        % Visualize. Index by user data to extract pc for a specific tree.
        % lastools(['lasview -i "',rootTLSData,'/classified.laz" -color_by_classification -keep_user_data 56']);
        % lastools(['lasview -i "',rootTLSData,'/classified.laz" -color_by_classification']);

    %%% 5. COMPUTE TREE METRICS
         h_butt = 0.2;       % Stump height
         h_interval = 0.2;   % diameter measuring interval
         h_layer = 0.05;      % point-cloud layer height
         min_points = 20;    % min number of points for circle fitting
         p = 0.5;           % stem curve smoothing parameter (cubic spline)
         computestemmetrics = true; % true/false whether to compute stem metrics
         computecrownmetrics = true; % true/false whether to compute crown metrics
         detectbranches = true; % true/false whether to detect branches
         [stemMetrics,crownMetrics] = pc_treeMetrics(rootResults,h_butt,h_interval,h_layer,min_points,p,computestemmetrics,computecrownmetrics,detectbranches);

         
