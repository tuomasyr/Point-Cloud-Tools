% Detect trees from TLS point clouds.
%  1. Point cloud to CHM raster using lascanopy-tool in LAStools. This
%     requires a licenced LAStools software (see lastools.m & https://rapidlasso.com/lastools).
%
%  2. Individual tree detection in R (see detect_trees_chm.R). This
%     requires the R software and the following R-packages installed:
%     - 'raster'      (https://www.rdocumentation.org/packages/raster)
%     - 'ForestTools' (https://www.rdocumentation.org/packages/ForestTools)
%     - 'sp'          (https://www.rdocumentation.org/packages/sp)
%     - 'rgdal'       (https://www.rdocumentation.org/packages/rgdal)
%
% Input:  - plot_id       % plot ID [double]
%         - chmPixelSize  % size of the CHM pixel in meters [double]
%         - minCrownArea  % minimum crown segment area in m^2 [double]
%         - plotCorners   % plot corner XY-coordinates [4x2 double array]
%         - rootALSData   % Directory for the ALS data [1x1 character array]
%         - filenameALS   % name of the ALS data file [1x1 character array]
%         - Rdir          % Directory for R.exe [1x1 character array]
%         - viewplot      % illustrate tree detection [logical]
%
% Output: - treePgons     % tree crown segment polygons [Mx1 cell array]
%
% (c) Tuomas Yrttimaa / Science4Trees @ UEF School of Forest Sciences 2019
% ---------------------------------------------------------------------------
function [treePgons] = pc_detect_trees_chm(chmPixelSize,minCrownArea,plotExtent,rootTLSData,filenameTLS,Rdir,viewplot) 
    starttime = datetime;
    hh = waitbar(0,''); hh.Name = 'Point Cloud Tools: Tree detection';
    perc = 5; waitbar(perc/100,hh,sprintf('%d%% Generating CHM...',perc));
    warning('off','MATLAB:polyshape:repairedBySimplify');
    
    % Make a directory if the tmp folder not yet exist
    saveDir = './tmp/treeloc'; if isfolder(saveDir) == false, mkdir(saveDir); end
    saveDir = './tmp/crownseg'; if isfolder(saveDir) == false, mkdir(saveDir); end
    saveDir = './tmp/figures'; if isfolder(saveDir) == false, mkdir(saveDir); end
    saveDir = './tmp/raster'; if isfolder(saveDir) == false, mkdir(saveDir); end
    
    % Create TLS-CHM
    lastools(sprintf(['lascanopy -i "',rootTLSData,'/',filenameTLS,'" -step %f -drop_z_below 1.3 -drop_z_above 45 -max -odir ./tmp/raster -o chm.tif -etrs89 -utm 35north'],chmPixelSize));

    % Export parameters and execute tree detection in R
    perc = 25; waitbar(perc/100,hh,sprintf('%d%% Segmenting tree crowns...',perc));
    a = system(['"',Rdir,'" CMD BATCH detect_trees_chm.R']);
    if a == 1, open detect_trees_chm.Rout; error('Failed to execute the R-command. Please check the R-out file');end

    % Import tree locs and segments
    perc = 55; waitbar(perc/100,hh,sprintf('%d%% Processing tree segmentation...',perc));
    treeseg = shaperead('./tmp/crownseg/crownseg_above.shp');

    treePgons = cell(2*size(treeseg,1),1);
    numTrees = 0;
    treeLocations = zeros(2*size(treePgons,1),2);
    for i = 1:size(treeseg,1)
        tmp = treeseg(i);
        R = regions(rmholes(polyshape(tmp.X,tmp.Y,'Simplify',false)));
        for j = 1:size(R,1)
            numTrees = numTrees + 1;
            if area(R(j)) >= minCrownArea
                treePgons{numTrees} = convhull(R(j));
                [treeLocations(numTrees,1),treeLocations(numTrees,2)] = centroid(convhull(R(j)));
            else
                continue
            end
        end
    end
    treePgons = treePgons(~cellfun('isempty',treePgons));
    treeLocations = treeLocations(treeLocations(:,1)>0,:);
    
    plotCorners = [plotExtent.min_x plotExtent.max_y; plotExtent.max_x plotExtent.max_y; plotExtent.max_x plotExtent.min_y;plotExtent.min_x plotExtent.min_y];
    bufferi = polybuffer(polyshape(plotCorners(:,1),plotCorners(:,2)),2.5);
    inID = inpolygon(treeLocations(:,1),treeLocations(:,2),bufferi.Vertices(:,1),bufferi.Vertices(:,2));
    treeLocations = treeLocations(inID,:);
    treePgonstmp = cell([size(treeLocations,1),1]);
    treeLocationstmp = [];
    numTrees = 0;
    emptyspace = polyshape(plotCorners(:,1),plotCorners(:,2));
    for i = 1:size(treePgons,1)
        poly = treePgons{i};
        if overlaps(poly,polyshape(cat(1,plotCorners(:,1),plotCorners(1,1)),cat(1,plotCorners(:,2),plotCorners(1,2))))
            emptyspace = subtract(emptyspace,poly);
            tmptree = treeLocations(inpolygon(treeLocations(:,1),treeLocations(:,2),poly.Vertices(:,1),poly.Vertices(:,2)),:);
            if ~isempty(tmptree)
               numTrees = numTrees + 1;
               treePgonstmp{numTrees} = poly;
               treeLocationstmp = cat(1,treeLocationstmp,tmptree);
            end
        end
    end
    treePgons = treePgonstmp(~cellfun('isempty',treePgonstmp));
    treeLocations = treeLocationstmp;

    R = regions(emptyspace);
    for j = 1:size(R,1)
        if area(R(j)) >= minCrownArea
            numTrees = numTrees + 1;
            poly = R(j);
            %treePgons{numTrees} = convhull(polybuffer(polyshape(poly.Vertices(:,1),poly.Vertices(:,2)),0.1));
            treePgons{numTrees} = polybuffer(polyshape(poly.Vertices(:,1),poly.Vertices(:,2)),0.7);
            [treeLocations(numTrees,1),treeLocations(numTrees,2)] = centroid(convhull(R(j)));
        else
            continue
        end
    end
    
    savedir = './tmp/treepoly'; if isfolder(savedir) == false, mkdir(savedir); end
    for i = 1:size(treePgons,1)
       poly = treePgons{i};
       writematrix(cat(1,poly.Vertices,poly.Vertices(1,:)),[savedir,sprintf('/treepoly_%i.txt',i)],'Delimiter','\t');
    end
    
    perc = 95; waitbar(perc/100,hh,sprintf('%d%% All done!',perc));
    
    % view results
    if viewplot == true
       winopen('./tmp/figures/chm_treeseg_above.png');
       figure, hold on, axis equal;
       plot(polyshape(plotCorners(:,1),plotCorners(:,2)),'FaceAlpha',0.3,'FaceColor',[0.15 0.15 0.15],'EdgeColor',[0.15 0.15 0.15]);
       for k = 1:size(treePgons,1)
          poly = treePgons{k};
          [centx,centy] = centroid(poly);
          plot(poly,'FaceColor',[26 204 62]./255,'EdgeColor',[0 163 32]./255);
          text(centx,centy,sprintf('%i',k),'Color','r');
          plot(centx,centy,'r+'),
       end
        hold off;
    end

    perc = 95; waitbar(perc/100,hh,sprintf('%d%% All done!',perc));
    endtime = datetime; runtime = endtime-starttime;
    cprintf('*comment',   sprintf('   %i trees detected in: %s\n',size(treePgons,1),runtime));
    close(hh);
end