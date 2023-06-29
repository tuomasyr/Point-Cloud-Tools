% Point cloud processing method for identifying and segmenting individual branches of trees. 
%-------------------------------------------------------------------------------------------
% This method has been implemented for southern boreal
% tree species and its performance has been validated with Scots pine 
% (Pinus sylvestris L.) trees in Yrttimaa et al. (2023):
% - - (citation, TBD).
% This function inputs a point cloud representing an individual tree, with ground removed and 
% height-normalized (i.e., Z-coordinates indicating height above the ground). 
% INPUT:
% - plot_id       <- sample plot identifier [integer]
% - tree_id       <- tree identifier [integer]
% - pc            <- point cloud, representing an individual tree, ground removed, height normalized [3-D pointCloud Object]
% - resultsfolder <- directory of a folder to store results [string]
%
% OUTPUT:
% - pcstemsurf    <- points representing stem surface [3-D pointCloud Object]
% - pcbranches    <- points of individual branches, segmentation ID stored as an Intensity value [3-D pointCloud Object]
% - pcrestofcrown <- points not assigned a branch segmentation [3-D pointCloud Object]
% - branchlist    <- a table of branches and their characteristics: id, 3D-location (cartesian & cylinder coordinates), horizontal length, vertical length, actual length, insertion angle,
%
% Required functions:
% - CircleFitByTaubin.m: Nikolai Chernov (2023). Circle Fit (Taubin method) (https://www.mathworks.com/matlabcentral/fileexchange/22678-circle-fit-taubin-method), MATLAB Central File Exchange.
% - euclideanDistanceTwoPointClouds.m: Audrey Cheong (2023). Euclidean distance between two point clouds (https://www.mathworks.com/matlabcentral/fileexchange/59377-euclidean-distance-between-two-point-clouds), MATLAB Central File Exchange.
%
% The function has been implemented in MATLAB R2022b (version 9.3)
% Required toolboxes:
% - Statistics and Machine Learning Toolbox (version 12.4)
% - Computer Vision Toolbox (version 10.3)
% - Lidar Toolbox (version 2.2)
%-----------------------------------------------------------------------------------------
% (c) Tuomas Yrttimaa // University of Eastern Finland // School of Forest Sciences (2023)
%-----------------------------------------------------------------------------------------
function [pcstemsurf,pcbranches,pcrestofcrown,branchlist]=pc_detect_branches(plot_id,tree_id,pc,resultsfolder)
    
    % Downsample and denoise the point cloud
    pc = pcdownsample(pc,'gridAverage',0.001);
    pc = pcdenoise(pc,"NumNeighbors",6);

    % Extract candidate stem points
    pccandstem = pcdownsample(pc,'gridAverage',0.005);
    pccandstem.Normal = pcnormals(pccandstem,30);
    pccandstem = pointCloud(pc.Location(abs(pc.Normal(:,3))< 0.005 & pc.Location(:,3) < pc.ZLimits(2)*0.8,:));
    pccandstem = pcdenoise(pccandstem,"NumNeighbors",5);
    %figure,pcshow(pccandstem)
    hints = 0.1:1:pccandstem.ZLimits(2);
    inliers = [];
    for i = 1:size(hints,2)-1
        pcstempart = pointCloud(pccandstem.Location(pccandstem.Location(:,3)>=hints(i) & pccandstem.Location(:,3)<hints(i+1),:));
        try
            [~,inlrinds] = pcfitcylinder(pcstempart,0.03,[0 0 1],'Confidence',99.99);
            inliers = cat(1,inliers,pcstempart.Location(inlrinds,:));
        catch
            try
                [~,inlrinds] = pcfitcylinder(pcstempart,0.04,[0 0 1],'Confidence',99.99);
                inliers = cat(1,inliers,pcstempart.Location(inlrinds,:));
            catch
                continue
            end
        end
    end
    pccandstem = pointCloud(inliers);
    pccandstem = pcdenoise(pccandstem,"NumNeighbors",5);

    % Extract stem orientation by circle fitting (CircleFitByTaubin.m):  
    % G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar Space Curves Defined By Implicit Equations, With Applications To Edge And Range Image Segmentation", 
    % IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
    hmin = 0.1; hmax = max(pccandstem.ZLimits);
    heightint = hmin:0.2:hmax;
    stemorient = [];
    for i = 1:size(heightint,2)-1
        stempc = pointCloud(pccandstem.Location(pccandstem.Location(:,3)>heightint(i) & pccandstem.Location(:,3)<=heightint(i+1),:));
        if stempc.Count > 20
            par = CircleFitByTaubin(stempc.Location(:,1:2));
            stemorient = cat(1,stemorient,cat(2,par(1:2),heightint(i)+0.1));
        end
    end
    if size(stemorient,1) > 10
        stemorientintX = interp1(stemorient(:,3),stemorient(:,1),0.1:0.05:pc.ZLimits(2)*1,'linear','extrap');
        stemorientintY = interp1(stemorient(:,3),stemorient(:,2),0.1:0.05:pc.ZLimits(2)*1,'linear','extrap');
        stemorientint = transpose(cat(1,stemorientintX,stemorientintY,0.1:0.05:pc.ZLimits(2)*1));
    else
        error('Cannot extract stem orientation axis - not enough observations!')
    end
    
    clear inliers inlrinds pccandstem
    %figure, scatter3(stemorientint(:,1),stemorientint(:,2),stemorientint(:,3)), axis equal

    % STEM SURFACE CLASSIFICATION
    % Cartesian-to-cylinder coordinate transformation
    pc_cylcoord = [];
    pc_cartcoord = [];
    for i = 1:size(stemorientint,1)-1
        pcslice = pc.Location(pc.Location(:,3)>=stemorientint(i,3) & pc.Location(:,3)<stemorientint(i+1,3),:);
        stemcent = mean(stemorientint(stemorientint(:,3)<=stemorientint(i,3)+0.05 & stemorientint(:,3)>=stemorientint(i,3)-0.05,:),1);
        if ~isempty(pcslice) && ~isempty(stemcent)
            y = pcslice(:,2)-stemcent(2);
            x = pcslice(:,1)-stemcent(1);
            theta = atan2(y,x);
            rho = sqrt(x.^2 + y.^2);
            pc_cylcoord = cat(1,pc_cylcoord,cat(2,theta,pcslice(:,3),rho.*4)); % Multiply rho by 4 to make branches more distinct
            pc_cartcoord = cat(1,pc_cartcoord,pcslice);
        end
    end
    pc_cylcoord = pointCloud(pc_cylcoord);
    %figure, pcshowpair(pc_cylcoord_stem,pc_cylcoord_scrown)
    pc_cartcoord = pointCloud(pc_cartcoord);
    [pc_cylcoord,inlierIndices,~] = pcdenoise(pc_cylcoord,'NumNeighbors',6);
    pc_cartcoord = select(pc_cartcoord,inlierIndices);
    
    % Segment stem surface <-- treat it as ground
    stemsurfPtsIdx = segmentGroundSMRF(pc_cylcoord,'ElevationThreshold',0.12,'MaxWindowRadius',5,'SlopeThreshold',0.8,'gridResolution',1);
    %nonstemsurfpc = select(pc_cylcoord,~stemsurfPtsIdx);
    nonstemsurfpc_cart = select(pc_cartcoord,~stemsurfPtsIdx);
    %stemsurfpc = select(pc_cylcoord,stemsurfPtsIdx);
    stemsurfpc_cart = select(pc_cartcoord,stemsurfPtsIdx);
    clear stemsurfPtsIdx pc_cylcoord pc_cartcoord pc pcstem
    %nonstemsurfpc = pointCloud(cat(2,nonstemsurfpc.Location(:,1:2),nonstemsurfpc.Location(:,3)./4));
    %stemsurfpc = pointCloud(cat(2,stemsurfpc.Location(:,1:2),stemsurfpc.Location(:,3)./4));
    %figure, pcshowpair(stemsurfpc,nonstemsurfpc), axis on;
    %figure, pcshowpair(stemsurfpc_cart,nonstemsurfpc_cart),axis on;
    
    %figure, pcshowpair(stemsurfpc_cart,pccrown);
    
    
    % Find stem height threshold
    hmin = 0.1; hmax = max(stemsurfpc_cart.ZLimits);
    heightint = hmin:0.2:hmax;
    stemdiam = [];
    for i = 1:size(heightint,2)-1
        stempc = pointCloud(stemsurfpc_cart.Location(stemsurfpc_cart.Location(:,3)>heightint(i) & stemsurfpc_cart.Location(:,3)<=heightint(i+1),:));
        if stempc.Count > 20
            par = CircleFitByTaubin(stempc.Location(:,1:2));
            stemdiam = cat(1,stemdiam,cat(2,par(3)*2,heightint(i)+0.1));
        end
    end
    
    stemdiam_smooth = stemdiam;
    stemdiam_smooth = rmoutliers(stemdiam_smooth,1);
    
    stemdiam_smooth_clip = stemdiam_smooth(stemdiam_smooth(:,2)>0.2*hmax & stemdiam_smooth(:,1) < mean(stemdiam_smooth(stemdiam_smooth(:,2)<=0.2*hmax,1)),:);
    stemdiam_smooth_clipped = stemdiam_smooth(stemdiam_smooth(:,2)<=0.2*hmax,:);
    
    % Use RANSAC polyfit to determine inliers/outliers
    sampleSize = 6; % number of points to sample per trial
    maxDistance = 0.0005; % max allowable distance for inliers
    
    points = cat(2,stemdiam_smooth_clip(:,2),stemdiam_smooth_clip(:,1));
    fitLineFcn = @(points) polyfit(points(:,1),points(:,2),2); % Fit a second-degree polynomial
    evalLineFcn = @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);
    
    [~, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
    
    modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1);
    
    inlierPts = points(inlierIdx,:);
    maxdiamdist = 2.1; % Max allowable vertical distance between two consecutive diameters
    interv = stemdiam_smooth_clip;
    interv(1,3) = 0;
    for ii = 1:size(interv,1)-1
        interv(ii,3) = interv(ii+1,2) - interv(ii,2);
    end
    maxhh = min(interv(interv(:,3) > maxdiamdist,2));
    
    if ~isempty(maxhh);inlierPts = inlierPts(inlierPts(:,1)<=maxhh,:);end
    
    maxstemh = max(inlierPts(:,1));
    
    % Plot stem curve
    x = [min(inlierPts(:,1)) max(inlierPts(:,1))]; y = modelInliers(1)*x + modelInliers(2);
    figure('visible','off'),plot(stemdiam_smooth(:,2),stemdiam_smooth(:,1).*100,'r.'), hold on, plot(x, y.*100, 'g-'), plot(inlierPts(:,1),inlierPts(:,2).*100,'go'),plot(inlierPts(:,1),inlierPts(:,2).*100,'k.'),plot(stemdiam_smooth_clipped(:,2),stemdiam_smooth_clipped(:,1).*100,'k.'); xlabel('stem height [m]'), ylabel('stem diameter [cm]')
    
    if ~isfolder(strcat(resultsfolder,sprintf('/%i/stemfigs',plot_id))); mkdir(strcat(resultsfolder,sprintf('/%i/stemfigs',plot_id)));end
    saveas(gcf,strcat(resultsfolder,sprintf('/%i/stemfigs/stemfig_%i.png',plot_id,tree_id)))
    
    % write classified point clouds
    pcstemsurf = pointCloud(stemsurfpc_cart.Location(stemsurfpc_cart.Location(:,3)<=maxstemh,:));
    pccrown = pointCloud(cat(1,nonstemsurfpc_cart.Location,stemsurfpc_cart.Location(stemsurfpc_cart.Location(:,3)>maxstemh,:)));
    %figure, pcshow(pcstemsurf)
    %figure, pcshowpair(pcstemsurf,pccrown)
    
    pcstemsurf = pcdownsample(pcstemsurf,"gridAverage",0.001);
    pccrown = pcdownsample(pccrown,"gridAverage",0.001);
    pcstemsurf.Color = uint8([170 115 30] .* ones(size(pcstemsurf.Location)));
    pcwrite(pcstemsurf,strcat(resultsfolder,sprintf('/%i/pointclouds/stem/%i_tree%i_stempoints.ply',plot_id,plot_id,tree_id)));
    pcwrite(pccrown,strcat(resultsfolder,sprintf('/%i/pointclouds/canopy/%i_tree%i_canopypoints.ply',plot_id,plot_id,tree_id)));

    %pcstemsurf = pcread(strcat(resultsfolder,sprintf('/%i/pointclouds/stem/%i_tree%i_stempoints.ply',plot_id,plot_id,tree_id)));
    %pccrown = pcread(strcat(resultsfolder,sprintf('/%i/pointclouds/canopy/%i_tree%i_canopypoints.ply',plot_id,plot_id,tree_id)));
    
    pcstemsurf = pcdownsample(pcstemsurf,'gridAverage',0.002);
    pccrown = pcdownsample(pccrown,'gridAverage',0.002);

    clear stemsurfpc_cart nonstemsurfpc_cart
    
    % Extract preliminary branch origins (knots) at points where stem and crown intersect
    hmin = 1; hmax = max(pcstemsurf.ZLimits);
    heightint = hmin:0.2:hmax;
    knotpc = [];
    %figure, pcshowpair(pcstemsurf,pccrown)
    for i = 1:size(heightint,2)-1
        pcslice = pccrown.Location(pccrown.Location(:,3)>heightint(i) & pccrown.Location(:,3)<=heightint(i+1),:);
        stempc = pcstemsurf.Location(pcstemsurf.Location(:,3)>heightint(i) & pcstemsurf.Location(:,3)<=heightint(i+1),:);
        %figure, pcshowpair(pcslice,stempc)
        if ~isempty(pcslice) && ~isempty(stempc)
            [dist,idx] = euclideanDistanceTwoPointClouds(stempc,pcslice);
            idx = idx(dist<0.06);
            if ~isempty(idx)
                knotpts = pcslice(dist<0.06,:);
                knotpc = cat(1,knotpc,knotpts);
                [dist,idx] = euclideanDistanceTwoPointClouds(knotpts,stempc);
                idx = idx(dist<0.05);
                if ~isempty(idx)
                    knotpts = stempc(dist<0.01,:);
                    knotpc = cat(1,knotpc,knotpts);
                end
            end
        end
    end
    knotpc = pointCloud(knotpc);
    knotpc = pcdenoise(knotpc,"NumNeighbors",10);
    
    %figure, pcshowpair(pcstemsurf,knotpc)
    
    % Cluster individual knots and define their location at stem surface
    [labels,numClusters] = pcsegdist(knotpc,0.03,'ParallelNeighborSearch',true);
    knotclustpc = cat(2,knotpc.Location,double(labels));
    %clrtrplts = colorcube(max(double(labels)));
    %clrtrplts = clrtrplts(randperm(max(double(labels))),:);
    %figure, pcshow(knotpc.Location,labels), colormap(clrtrplts);
    
    knotlocs = [];
    for j = 1:numClusters
        pcCluster = pointCloud(knotclustpc(knotclustpc(:,4) == j,1:3));
        stempc = pcstemsurf.Location(pcstemsurf.Location(:,3)>pcCluster.ZLimits(1)-0.15 & pcstemsurf.Location(:,3)<pcCluster.ZLimits(1)+0.15,:);
        %figure, pcshowpair(pcstemsurf,pcCluster)
        clustdim = cat(2,pcCluster.XLimits(2)-pcCluster.XLimits(1),pcCluster.YLimits(2)-pcCluster.YLimits(1),pcCluster.ZLimits(2)-pcCluster.ZLimits(1));
        if pcCluster.Count > 5 && sum(clustdim) < 0.5
            [dist,~] = euclideanDistanceTwoPointClouds(pcCluster.Location,stempc);
            if ~isempty(dist)
                knotlocs = cat(1,knotlocs,mean(stempc(dist<0.05,:),1));
            end
        else
            continue
        end
    end
    knotlocs2 = knotlocs(knotlocs(:,3)>0,:);
    knotlocs3 = [];
    while ~isempty(knotlocs2)
        knot = knotlocs2(1,:);
        [dist,~] = euclideanDistanceTwoPointClouds(knot,knotlocs2);
        nghbrknts = knotlocs2(dist<0.05,:);
        if ~isempty(nghbrknts)
            knotlocs3 = cat(1,knotlocs3,mean(nghbrknts,1));
            knotlocs2(dist<0.05,:) = [];
        end
    end
    clear knotclustpc
    
    % Branch detection #1
    branches = cell(size(knotlocs3,1),1);
    pcrown = cat(2,pccrown.Location,transpose(1:pccrown.Count));
    for kk = 1:size(knotlocs3,1)
        knotloc = knotlocs3(kk,:);
        branchpc = knotloc;
        crownslice = pcrown(pcrown(:,3)<knotloc(3)+1.5 & pcrown(:,3)>knotloc(3)-1.5,:);
        ids = crownslice(:,4);
        crownslice = pointCloud(crownslice(:,1:3));
        [labels,numClusters] = pcsegdist(crownslice,0.05,'ParallelNeighborSearch',false);
        branchclustpc = cat(2,crownslice.Location,double(labels));
        lbls = [];
        for j = 1:numClusters
            clustpc = branchclustpc(branchclustpc(:,4)==j,1:3);
            [dist,~] = euclideanDistanceTwoPointClouds(knotloc,clustpc);
            dist = min(dist);
            if dist < 0.06
                branchpc = cat(1,branchpc,clustpc);
                lbls = cat(1,lbls,j);
            end
        end
        branchpc = pointCloud(branchpc);
        inds = ismember(labels,lbls);
        remvids = ids(inds);
        pcrown = pcrown(~ismember(pcrown(:,4),remvids),:);
        branchpc.Intensity = ones(branchpc.Count,1).*kk;
        branches{kk} = branchpc;
        %fprintf('%.2f %%\n',kk/size(knotlocs3,1)*100)
    end
    clear inds labels
    %toc
    branchpc = [];
    colors = [];
    for kk = 1:size(knotlocs3,1)
        branch = branches{kk};
        branchpc = cat(1,branchpc,branch.Location);
        colors = cat(1,colors,branch.Intensity);
    end
    branchpc = pointCloud(branchpc);
    branchpc.Intensity = colors;
    colortrplts = colorcube(max(branchpc.Intensity));
    branchpc.Color = colortrplts(branchpc.Intensity,:);
    
    %figure, pcshow(pcstemsurf),hold on; pcshow(branchpc); hold on; scatter3(knotlocs3(:,1),knotlocs3(:,2),knotlocs3(:,3),60,'w','filled');
    warning('off','MATLAB:alphaShape:DupPointsBasicWarnId');
    branchpcsample = pcdownsample(branchpc,"gridAverage",0.01);
    branchAlpha = alphaShape(branchpcsample.Location(:,1),branchpcsample.Location(:,2),branchpcsample.Location(:,3),0.1);
    crownpc = pointCloud(pccrown.Location(logical(abs(~inShape(branchAlpha,pccrown.Location))),:));
    
    [labels,numClusters] = pcsegdist(crownpc,0.04,'ParallelNeighborSearch',true);
    crownclustpc = cat(2,crownpc.Location,double(labels));
    
    knotlocs4 = [];
    addbranches = [];
    count = 0;
    for j = 1:numClusters
        pcCluster = pointCloud(crownclustpc(crownclustpc(:,4) == j,1:3));
        stempc = pcstemsurf.Location(pcstemsurf.Location(:,3)>pcCluster.ZLimits(1)-0.1 & pcstemsurf.Location(:,3)<=pcCluster.ZLimits(2)+0.1,:);
        %figure, pcshowpair(pcstemsurf,pcCluster)
        %figure, pcshowpair(pointCloud(stempc),pcCluster)
        clustdim = cat(2,pcCluster.XLimits(2)-pcCluster.XLimits(1),pcCluster.YLimits(2)-pcCluster.YLimits(1),pcCluster.ZLimits(2)-pcCluster.ZLimits(1));
        hortzdim = sqrt(clustdim(1)^2+clustdim(2)^2);
        if pcCluster.Count > 50 && hortzdim > 0.2
            [dist,~] = euclideanDistanceTwoPointClouds(pcCluster.Location,stempc);
            [minval,idx] = min(dist);
            if ~isempty(idx) && minval < 0.1
                count = count + 1;
                knotlocs4 = cat(1,knotlocs4,stempc(idx,:));
                addbranches = cat(1,addbranches,cat(2,pcCluster.Location,count*ones([pcCluster.Count,1])));
            else
                continue
            end
        else
            continue
        end
        %fprintf('%.2f %%\n',j/numClusters*100);
    end
    clear crownclustpc crownpc labels
    if ~isempty(addbranches)
        addbranchpc = pointCloud(addbranches(:,1:3));
        addbranchpc.Intensity = addbranches(:,4)+max(branchpc.Intensity);
        branch_pc = pointCloud(cat(1,branchpc.Location,addbranchpc.Location));
        branch_pc.Intensity = cat(1,branchpc.Intensity,addbranchpc.Intensity);
        colortrplts = colorcube(max(branch_pc.Intensity));
        colortrplts = colortrplts(randperm(max(branch_pc.Intensity)),:);
        branch_pc.Color = colortrplts(branch_pc.Intensity,:);
    else
        branch_pc = branchpc;
    end
    if ~isempty(knotlocs4)
        knotlocs3 = cat(1,knotlocs3,knotlocs4);
    end
    clear branchpc branchAlpha
    
    %figure, pcshow(pcstemsurf), hold on, pcshow(branch_pc); scatter3(knotlocs3(:,1),knotlocs3(:,2),knotlocs3(:,3),60,'w','filled');
    
    
    % Convert knot locations from cartesian to cylinder coordinates
    knotlocs_cylcoords = [];
    for i = 1:size(knotlocs3,1)
        knotloc = knotlocs3(i,:);
        stemcent = mean(stemorientint(stemorientint(:,3)<=knotloc(3)+0.04 & stemorientint(:,3)>=knotloc(3)-0.04,:),1);
        if ~isempty(knotloc) && ~isempty(stemcent)
            y = knotloc(2)-stemcent(2);
            x = knotloc(1)-stemcent(1);
            theta = atan2(y,x);
            rho = sqrt(x^2 + y^2);
            knotlocs_cylcoords = cat(1,knotlocs_cylcoords,[knotloc theta rho knotloc(3)]);
        end
    end
    knotlocs_cylcoords(isnan(knotlocs_cylcoords(:,5)),:) = [];
    knotlocs_cylcoords = array2table(knotlocs_cylcoords,"VariableNames",{'x','y','z','theta','rho','h'});
    knotlocs_cylcoords.h_rel = knotlocs_cylcoords.h./max(pccrown.ZLimits);
    
    % Branch detection #2
    branchpc_cylcoord = [];
    colors_cylcoords = [];
    intensity_cylcoords = [];
    for i = 1:size(stemorientint,1)-1
        pcslice = branch_pc.Location(branch_pc.Location(:,3)>=stemorientint(i,3) & branch_pc.Location(:,3)<stemorientint(i+1,3),:);
        pccolor = branch_pc.Color(branch_pc.Location(:,3)>=stemorientint(i,3) & branch_pc.Location(:,3)<stemorientint(i+1,3),:);
        pcintensity = branch_pc.Intensity(branch_pc.Location(:,3)>=stemorientint(i,3) & branch_pc.Location(:,3)<stemorientint(i+1,3),:);
        stemcent = mean(stemorientint(stemorientint(:,3)<=stemorientint(i,3)+0.05 & stemorientint(:,3)>=stemorientint(i,3)-0.05,:),1);
        if ~isempty(pcslice) && ~isempty(stemcent)
            y = pcslice(:,2)-stemcent(2);
            x = pcslice(:,1)-stemcent(1);
            theta = atan2(y,x);
            rho = sqrt(x.^2 + y.^2);
            branchpc_cylcoord = cat(1,branchpc_cylcoord,cat(2,theta,pcslice(:,3),rho));
            colors_cylcoords = cat(1,colors_cylcoords,pccolor);
            intensity_cylcoords = cat(1,intensity_cylcoords,pcintensity);
        end
    end
    branchpc_cylcoord = pointCloud(branchpc_cylcoord);
    branchpc_cylcoord.Color = colors_cylcoords;
    branchpc_cylcoord.Intensity = intensity_cylcoords;
    clear branch_pc intensity_cylcoords colors_cylcoords
    %figure, pcshow(branchpc_cylcoord); hold on; scatter3(knotlocs_cylcoords.theta,knotlocs_cylcoords.h,knotlocs_cylcoords.rho,60,'w','filled'); xlabel('theta [rad]'),ylabel('h [m]')
    
    warning('off','MATLAB:polyshape:repairedBySimplify');
    branchclustpc = cat(2,branchpc_cylcoord.Location,double(branchpc_cylcoord.Intensity));
    branchlist_cylcoord = cell([size(knotlocs_cylcoords,1),1]);
    branchlist = [];
    count = 0;
    for i = 1:max(double(branchpc_cylcoord.Intensity))
        clustpc = pointCloud(branchclustpc(branchclustpc(:,4) == i,1:3));
        if ~isempty(clustpc.Location)
            disttheta = euclideanDistanceTwoPointClouds(clustpc.Location(:,1),knotlocs_cylcoords.theta);
            [dist,~] = euclideanDistanceTwoPointClouds(clustpc.Location(:,2:3),[knotlocs_cylcoords.h,knotlocs_cylcoords.rho]);
            %[dist,~] = euclideanDistanceTwoPointClouds(clustpc.Location,[knotlocs_cylcoords.theta,knotlocs_cylcoords.h,knotlocs_cylcoords.rho]);
            kntlcs = knotlocs_cylcoords(dist<0.08 & rad2deg(disttheta) < 8,:);
            if ~isempty(kntlcs) && clustpc.Count > 4
                if size(kntlcs,1) == 1
                    count = count + 1;
                    branchlist_cylcoord{count} = clustpc;
                    branchlist = cat(1,branchlist,[count table2array(kntlcs)]);
                elseif size(kntlcs,1) > 1
                    clustconvh = polyshape(clustpc.Location(convhull(clustpc.Location(:,1:2)),1:2));
                    [Xlim,Ylim] = boundingbox(clustconvh);
                    bbox = polyshape([Xlim(1) Xlim(1) Xlim(2) Xlim(2)],[Ylim(2) Ylim(1) Ylim(1) Ylim(2)]);
                    %figure, plot(bbox),hold on; plot(clustconvh);
    
                    inknot = inpolygon(kntlcs.theta,kntlcs.h,bbox.Vertices(:,1),bbox.Vertices(:,2));
                    kntlcs = kntlcs(inknot,:);
                    if size(kntlcs,1) == 1
                        count = count + 1;
                        branchlist_cylcoord{count} = clustpc;
                        branchlist = cat(1,branchlist,[count table2array(kntlcs)]);
                    elseif size(kntlcs,1) > 1
                        disth = mean(sqrt((diff(kntlcs.h)).^2));
                        distang = rad2deg(mean(sqrt((diff(kntlcs.theta)).^2)));
                        if disth < 0.03 && distang < 8
                            kntlcs = mean(table2array(kntlcs),1);
                            count = count + 1;
                            branchlist_cylcoord{count} = clustpc;
                            branchlist = cat(1,branchlist,[count kntlcs]);
                        else
                            [xq,yq] = meshgrid(Xlim(1):0.1:Xlim(2),Ylim(1):0.01:Ylim(2));
                            xlist = []; ylist = [];
                            for k = 1:size(xq,2)
                                xlist = cat(1,xlist,xq(:,k)); ylist = cat(1,ylist,yq(:,k));
                            end
                            in = inpolygon(xlist,ylist,bbox.Vertices(:,1),bbox.Vertices(:,2));
                            gridpoints = cat(2,xlist,ylist); gridpoints = gridpoints(in,:);
                            if ~isempty(gridpoints)
    
    
                                % Group the gridpoints by finding closest knot for each gridpoint
                                knotsegcell = cell(size(kntlcs,1));
                                for k = 1:size(gridpoints,1)
                                    point = gridpoints(k,:);
                                    kdtree = ExhaustiveSearcher(cat(2,kntlcs.theta,kntlcs.h));
                                    ni = knnsearch(kdtree,point);
                                    temp = knotsegcell{ni}; temp = cat(1,temp,point); knotsegcell{ni} = temp;
                                end
    
                                % Generate split crown segments from the convex hull of the grouped gridpoints
                                branchsegments = cell(size(kntlcs,1));
                                for k = 1:size(knotsegcell,1)
                                    try
                                        K = convhull(knotsegcell{k}(:,1),knotsegcell{k}(:,2), 'simplify',true);
                                        branchsegmt = polyshape(knotsegcell{k}(K,1),knotsegcell{k}(K,2));
                                        %figure, plot(branchsegmt), hold on; plot(kntlc(:,1),kntlc(:,2),'*')
                                        branchsegments{k} = branchsegmt;
                                        kntlc = kntlcs(inpolygon(kntlcs.theta,kntlcs.h,branchsegmt.Vertices(:,1),branchsegmt.Vertices(:,2)),:);
                                        branch = pointCloud(clustpc.Location(inpolygon(clustpc.Location(:,1),clustpc.Location(:,2),branchsegmt.Vertices(:,1),branchsegmt.Vertices(:,2)),:));
                                        if ~isempty(kntlc)
                                            count = count + 1;
                                            branchlist_cylcoord{count} = branch;
                                            branchlist = cat(1,branchlist,[count table2array(kntlc)]);
                                        end
                                    catch
                                        continue
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    branchlist_cylcoord = branchlist_cylcoord(~cellfun('isempty',branchlist_cylcoord));
    branchlist = array2table(branchlist,'VariableNames',{'id','x','y','z','theta','rho','h','h_rel'});
    clear branchpc_cylcoord branchclustpc
    
    knotlist2 = branchlist(branchlist.h > 0,:);
    knotlist3 = [];
    removeids = [];
    while ~isempty(knotlist2)
        knot = knotlist2(1,:);
        [dist,~] = euclideanDistanceTwoPointClouds([knot.x knot.y knot.z],[knotlist2.x knotlist2.y knotlist2.z]);
        nghbrknts = knotlist2(dist<0.05,:);
        if size(nghbrknts,1) == 1
            knotlist3 = cat(1,knotlist3,table2array(nghbrknts));
            knotlist2(dist<0.05,:) = [];
        elseif size(nghbrknts,1) > 1
            tobemerged = [];
            brnchids = nghbrknts.id;
            for jj = 1:size(nghbrknts,1)
                pcc = branchlist_cylcoord{brnchids(jj)};
                tobemerged = cat(1,tobemerged,pcc.Location);
            end
            branchlist_cylcoord{brnchids(1)} = pointCloud(tobemerged);
            knotinfo = mean(table2array(nghbrknts),1);
            knotinfo(1) = brnchids(1);
            removeids = cat(1,removeids,brnchids(2:end));
            knotlist3 = cat(1,knotlist3,knotinfo);
            knotlist2(dist<0.05,:) = [];
        end
    end
    branchlist_cylcoord = branchlist_cylcoord(knotlist3(:,1));
    branchlist = array2table(knotlist3,'VariableNames',{'id','x','y','z','theta','rho','h','h_rel'});
    
    branchpc2 = [];
    colors2 = [];
    branchendlist = [];
    branchmidlist = [];
    removeids = [];
    for kk = 1:size(branchlist_cylcoord,1)
        branch = branchlist_cylcoord{kk};
        if ~isempty(branch.Location)
            kntlc = table2array(branchlist(kk,["theta" "h" "rho"]));
            id = table2array(branchlist(kk,"id"));
            branchend = median(branch.Location(branch.Location(:,3) == branch.ZLimits(2),:),1);
            branchmid = median(branch.Location(branch.Location(:,3) < branch.ZLimits(1)+0.3,:),1);
            branchendlist = cat(1,branchendlist,[id branchend]);
            branchmidlist = cat(1,branchmidlist,[id branchmid]);
            branchlist.horzlength(kk) = branch.ZLimits(2) - kntlc(3);
            branchpc2 = cat(1,branchpc2,branch.Location);
            colors2 = cat(1,colors2,id*ones([branch.Count,1]));
        else
            removeids = cat(1,removeids,kk);
        end
    end
    branchlist(removeids,:) = [];
    clear branchclustpc branchlist_cylcoord
    branchpc2 = pointCloud(branchpc2);
    colortrplts2 = colorcube(max(colors2));
    colortrplts2 = colortrplts2(randperm(max(colors2)),:);
    branchpc2.Color = colortrplts2(colors2,:);
    branchpc2.Intensity = colors2;
    clear colors colors2
    %stemsurfpc = pointCloud(cat(2,stemsurfpc.Location(:,1),stemsurfpc.Location(:,2),stemsurfpc.Location(:,3)./4));
    %stemsurfpc.Color = ones([stemsurfpc.Count,1])*[0.7 0.5 0];
    %figure, pcshow(branchpc2); hold on; pcshow(stemsurfpc); scatter3(knotlist.theta,knotlist.h,knotlist.rho,60,'w','filled');
    %figure, pcshow(branchpc2); hold on; pcshow(stemsurfpc); scatter3(knotlocs_cylcoords.theta,knotlocs_cylcoords.h,-knotlocs_cylcoords.rho,60,'w','filled');
    
    
    %tic
    branchpc3 = [];
    colors3 = [];
    intensity3 = [];
    branchendlist_cyl2cart = [];
    branchmidlist_cyl2cart = [];
    %branchlist_cyl2cart = cell([max(knotlist.id),1]);
    %idd = 0;
    for i = 1:size(stemorientint,1)-1
        pcslice = branchpc2.Location(branchpc2.Location(:,2)>=stemorientint(i,3) & branchpc2.Location(:,2)<stemorientint(i+1,3),:);
        pccolor = branchpc2.Color(branchpc2.Location(:,2)>=stemorientint(i,3) & branchpc2.Location(:,2)<stemorientint(i+1,3),:);
        pcintensity = branchpc2.Intensity(branchpc2.Location(:,2)>=stemorientint(i,3) & branchpc2.Location(:,2)<stemorientint(i+1,3),:);
        branchend = branchendlist(branchendlist(:,3)>=stemorientint(i,3) & branchendlist(:,3)<stemorientint(i+1,3),:);
        branchmid = branchmidlist(branchmidlist(:,3)>=stemorientint(i,3) & branchmidlist(:,3)<stemorientint(i+1,3),:);
        stemcent = mean(stemorientint(stemorientint(:,3)<=stemorientint(i,3)+0.05 & stemorientint(:,3)>=stemorientint(i,3)-0.05,:),1);
        if ~isempty(pcslice) && ~isempty(stemcent)
            x = pcslice(:,3).*sin(-pcslice(:,1)+pi/2) + stemcent(1);
            y = pcslice(:,3).*cos(-pcslice(:,1)+pi/2) + stemcent(2);
            branchpc3 = cat(1,branchpc3,cat(2,x,y,pcslice(:,2)));
            colors3 = cat(1,colors3,pccolor);
            intensity3 = cat(1,intensity3,pcintensity);
            if ~isempty(branchend)
                x = branchend(:,4).*sin(-branchend(:,2)+pi/2) + stemcent(1);
                y = branchend(:,4).*cos(-branchend(:,2)+pi/2) + stemcent(2);
                branchendlist_cyl2cart = cat(1,branchendlist_cyl2cart,cat(2,branchend(:,1),x,y,branchend(:,3)));
            end
            if ~isempty(branchmid)
                x = branchmid(:,4).*sin(-branchmid(:,2)+pi/2) + stemcent(1);
                y = branchmid(:,4).*cos(-branchmid(:,2)+pi/2) + stemcent(2);
                branchmidlist_cyl2cart = cat(1,branchmidlist_cyl2cart,cat(2,branchmid(:,1),x,y,branchmid(:,3)));
            end
        end
        %fprintf('%.2f %%\n',(i/size(stemorientint,1))*100);
    end
    branchpc3 = pointCloud(branchpc3);
    branchpc3.Color = colors3;
    branchpc3.Intensity = intensity3;
    clear branchpc2 colors3 intensity3
    %toc
    
    %figure, pcshow(pcstemsurf),hold on; pcshow(branchpc3); scatter3(knotlist.x,knotlist.y,knotlist.z,60,'w','filled');
    
    branchlist.branchendX = zeros([size(branchlist,1),1]);
    branchlist.branchendY = zeros([size(branchlist,1),1]);
    branchlist.branchendZ = zeros([size(branchlist,1),1]);
    branchlist.branchmidX = zeros([size(branchlist,1),1]);
    branchlist.branchmidY = zeros([size(branchlist,1),1]);
    branchlist.branchmidZ = zeros([size(branchlist,1),1]);
    for kk = 1:size(branchlist,1)
        kntlc = branchlist(kk,:);
        if ismember(kntlc.id,branchendlist_cyl2cart(:,1))
            branchlist.branchendX(kk) = branchendlist_cyl2cart(branchendlist_cyl2cart(:,1)==kntlc.id,2);
            branchlist.branchendY(kk) = branchendlist_cyl2cart(branchendlist_cyl2cart(:,1)==kntlc.id,3);
            branchlist.branchendZ(kk) = branchendlist_cyl2cart(branchendlist_cyl2cart(:,1)==kntlc.id,4);
        end
        if ismember(kntlc.id,branchmidlist_cyl2cart(:,1))
            branchlist.branchmidX(kk) = branchmidlist_cyl2cart(branchmidlist_cyl2cart(:,1)==kntlc.id,2);
            branchlist.branchmidY(kk) = branchmidlist_cyl2cart(branchmidlist_cyl2cart(:,1)==kntlc.id,3);
            branchlist.branchmidZ(kk) = branchmidlist_cyl2cart(branchmidlist_cyl2cart(:,1)==kntlc.id,4);
        end
    end
    branchlist.length = sqrt((branchlist.x-branchlist.branchendX).^2+(branchlist.y-branchlist.branchendY).^2+(branchlist.z-branchlist.branchendZ).^2);
    branchlist.vertlength = abs(branchlist.branchendZ-branchlist.z);
    branchlist.insertangle_end = 90-rad2deg(atan((branchlist.branchendZ-branchlist.z)./(branchlist.horzlength)));
    branchlist.insertangle_mid = 90-rad2deg(atan((branchlist.branchmidZ-branchlist.z)./sqrt((branchlist.x-branchlist.branchmidX).^2 + (branchlist.y-branchlist.branchmidY).^2)));
    
    branchlist(branchlist.horzlength < 0.01,:) = [];
    branchlist(branchlist.horzlength < 0.1 & branchlist.h_rel > 0.4,:) = [];
    branchlist(branchlist.horzlength < 0.05 & branchlist.vertlength > 0.2,:) = [];
    branchlist(branchlist.insertangle_mid < 15 | branchlist.insertangle_mid > 155,:) = [];
    
    branchpointclouds = cell([size(branchlist,1),1]);
    %figure, pcshow(pcstemsurf), hold on, pcshow(branchpc3);scatter3(knotlist.x,knotlist.y,knotlist.z,40,'w','filled');
    branchpc4 = [];
    intensitylist = [];
    for i = 1:size(branchlist,1)
        branch = pointCloud(branchpc3.Location(branchpc3.Intensity==branchlist.id(i),:));
        branch.Intensity = branchpc3.Intensity(branchpc3.Intensity==branchlist.id(i),:);
        intensitylist = cat(1,intensitylist,branch.Intensity);
        branchpointclouds{i} = branch;
        branchlist.id(i) = i;
        branchpc4 = cat(1,branchpc4,cat(2,branch.Location,i*ones([branch.Count,1])));
    end
    clear branchpc3
    
    ids = branchpc4(:,4);
    branchpc4 = pointCloud(branchpc4(:,1:3));
    colortrplts = colorcube(max(ids));
    colortrplts = colortrplts(randperm(max(ids)),:);
    branchpc4.Color = colortrplts(ids,:);
    branchpc4.Intensity = ids;
    clear ids
    
    branchpcsample = pcdownsample(branchpc4,"gridAverage",0.01);
    branchAlpha = alphaShape(branchpcsample.Location(:,1),branchpcsample.Location(:,2),branchpcsample.Location(:,3),0.1);
    crownpc = pointCloud(pccrown.Location(logical(abs(~inShape(branchAlpha,pccrown.Location))),:));
    topcrownpc = pointCloud(crownpc.Location(crownpc.Location(:,3)>=branchpc4.ZLimits(2)-0.4,:));
    crownpc = pointCloud(crownpc.Location(crownpc.Location(:,3)<branchpc4.ZLimits(2)-0.4,:));
    
    crownpc.Color = [1 1 1].*ones([crownpc.Count,1]);
    %figure, pcshow(pcstemsurf), hold on, pcshow(branchpc4.Location,branchpc4.Color.*0.85);pcshow(crownpc);scatter3(knotlist.x,knotlist.y,knotlist.z,40,'yellow');
    %figure, pcshow(pcstemsurf), hold on, pcshowpair(branchpc4,crownpc); plot(branchAlpha, 'FaceAlpha',0.1,'FaceColor','g','EdgeColor','r')
    clear branchAlpha
    
    try
        [labels,~] = pcsegdist(crownpc,0.06,'ParallelNeighborSearch',true);
        lbls = unique(labels);
        c = arrayfun(@(x)length(find(labels == x)), unique(labels), 'Uniform', false);
        c = cell2mat(c);
        lbls = double(lbls(c > 50));
    
        crownpc = cat(2,crownpc.Location,double(labels));
        %clrtrplts = colorcube(max(double(labels)));
        %clrtrplts = clrtrplts(randperm(max(double(labels))),:);
        %figure, pcshow(pcstemsurf), hold on;pcshow(crownpc.Location,labels), colormap(clrtrplts);scatter3(knotlist.x,knotlist.y,knotlist.z,40,'w','filled');
        clear labels
        pcrestofcrown = [];
        for j = 1:size(lbls,1)
            pcCluster = pointCloud(crownpc(crownpc(:,4) == lbls(j),1:3));
            branchespc = branchpc4.Location(branchpc4.Location(:,3) > pcCluster.ZLimits(1)-0.2 & branchpc4.Location(:,3) < pcCluster.ZLimits(2)+0.2,:);
            branchesint = branchpc4.Intensity(branchpc4.Location(:,3) > pcCluster.ZLimits(1)-0.2 & branchpc4.Location(:,3) < pcCluster.ZLimits(2)+0.2,:);
            [dist,~] = euclideanDistanceTwoPointClouds(pcCluster.Location,branchespc);
            [minval,idx] = min(dist);
            if ~isempty(idx) && minval < 0.06
                branchid = branchesint(idx);
                branch = branchpointclouds{branchid};
                if ~isempty(branch.Location)
                    new_branchpc = pointCloud(cat(1,branch.Location,pcCluster.Location));
                    new_branchpc.Intensity = max(branch.Intensity) .* ones([new_branchpc.Count,1]);
                    branchpointclouds{branchid} = new_branchpc;
                    %fprintf('added points to branch # %i\n',branchid);
                end
            else
                pcrestofcrown = cat(1,pcrestofcrown,cat(2,pcCluster.Location));
            end
            %fprintf('%.2f %%\n',j/size(lbls,1)*100);
        end
        clear crownpc dist
        pcrestofcrown = pointCloud(cat(1,pcrestofcrown,topcrownpc.Location));
        pcrestofcrown.Color = [1 1 1].*ones([pcrestofcrown.Count,1]);
    
        pcbranches = [];
        %intensitylist2 = [];
        intensitylist1 = [];
        for i = 1:size(branchlist,1)
            branch = branchpointclouds{i};
            pcbranches = cat(1,pcbranches,branch.Location);
            %intensitylist2 = cat(1,intensitylist2,branch.Intensity);
            intensitylist1 = cat(1,intensitylist1,i*ones([branch.Count,1]));
        end
        pcbranches = pointCloud(pcbranches);
        pcbranches.Intensity = intensitylist1;
        clrtrplts = colorcube(max(intensitylist1));
        clrtrplts = clrtrplts(randperm(max(intensitylist1)),:);
        pcbranches.Color = clrtrplts(intensitylist1,:);
        clear intensitylist1 branchpc4
    catch
        pcbranches = branchpc4;
        pcrestofcrown = topcrownpc;
    end
    %figure, pcshow(pcstemsurf), hold on, pcshow(branchpc5); pcshow(restofcrown); scatter3(knotlist.x,knotlist.y,knotlist.z,50,'w','filled');
    %figure, pcshow(pcstemsurf), hold on, pcshowpair(branchpc5,restofcrown);
    
    pcwrite(pcbranches,strcat(resultsfolder,sprintf('/%i/pointclouds/canopy/%i_tree%i_branches.ply',plot_id,plot_id,tree_id)));
    pcwrite(pcrestofcrown,strcat(resultsfolder,sprintf('/%i/pointclouds/canopy/%i_tree%i_restofcrown.ply',plot_id,plot_id,tree_id)));
    
    savedir = strcat(resultsfolder,sprintf('/%i/branchlists',plot_id));
    if ~isfolder(savedir);mkdir(savedir);end
    writetable(branchlist,strcat(savedir,sprintf('/%i_tree%i_branchlist.csv',plot_id,tree_id)))
    
    try
        figure('visible','off'), set(gcf,'position',[50,50,1500,900]);
        pcstemsurf.Color = uint8([170 115 30] .* ones(size(pcstemsurf.Location)));
        pccrown.Color = uint8([110 190 100] .* ones(size(pccrown.Location)));
    
        subplot(1,4,1);
        scatter3(pcstemsurf.Location(:,1),pcstemsurf.Location(:,2),pcstemsurf.Location(:,3),1,double(pcstemsurf.Color(1,:))./255), view(0,0); title('a)'); axis equal; axis on; zlabel('tree height [m]'); zlim([0 max(pccrown.ZLimits)*1.03]),xlim([min(pccrown.XLimits) max(pccrown.XLimits)])
        %subplot(1,3,1); pcshowpair(pcstemsurf,pccrown,BackgroundColor="white"), view(0,0); axis on; zlabel('tree height [m]')
    
        subplot(1,4,2);
        scatter3(pccrown.Location(:,1),pccrown.Location(:,2),pccrown.Location(:,3),1,double(pccrown.Color(1,:))./255); view(0,0); title('b)');axis equal; axis on; zlabel('tree height [m]'); zlim([0 max(pccrown.ZLimits)*1.03]),xlim([min(pccrown.XLimits) max(pccrown.XLimits)])
        %subplot(1,3,1); pcshowpair(pcstemsurf,pccrown,BackgroundColor="white"), view(0,0); axis on; zlabel('tree height [m]')
        subplot(1,4,3);
        scatter3(pcstemsurf.Location(:,1),pcstemsurf.Location(:,2),pcstemsurf.Location(:,3),1,double(pcstemsurf.Color(1,:))./255), hold on;
        scatter3(pccrown.Location(:,1),pccrown.Location(:,2),pccrown.Location(:,3),1,double(pccrown.Color(1,:))./255);
        scatter3(branchlist.x,branchlist.y,branchlist.z,50,'k','filled'); view(0,0); title('c)');axis equal; axis on; zlabel('tree height [m]'); zlim([0 max(pccrown.ZLimits)*1.03])
    
        subplot(1,4,4);
        scatter3(pcstemsurf.Location(:,1),pcstemsurf.Location(:,2),pcstemsurf.Location(:,3),1,double(pcstemsurf.Color(1,:))./255), hold on;
        scatter3(pcrestofcrown.Location(:,1),pcrestofcrown.Location(:,2),pcrestofcrown.Location(:,3),1,[0.85 0.85 0.85].*ones([pcrestofcrown.Count,1]));
        scatter3(pcbranches.Location(:,1),pcbranches.Location(:,2),pcbranches.Location(:,3),1,double(pcbranches.Color)./255.*0.75);
        view(0,0); axis equal; axis on; title('d)');zlabel('tree height [m]'); zlim([0 max(pccrown.ZLimits)*1.03])%xlim([0.9*min(pccrown.XLimits) 1.1*max(pccrown.XLimits)])
        scatter3(branchlist.x,branchlist.y,branchlist.z,50,'k','filled'); view(0,0); axis equal; axis on; zlabel('tree height [m]'); zlim([0 max(pccrown.ZLimits)*1.03])
        savedir = strcat(resultsfolder,sprintf('/%i/figures/branchclassif',plot_id));
        if ~isfolder(savedir);mkdir(savedir);end
        saveas(gcf,strcat(savedir,sprintf('/%i_tree%i_branchclassif.png',plot_id,tree_id)))
    catch
        warning('Could not visualize branch segmentation for tree # %i (out of memory)',tree_id)
        return
    end
end
