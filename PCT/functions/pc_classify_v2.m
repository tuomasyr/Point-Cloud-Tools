% Point cloud classifier. Classifies TLS point cloud into stem points and
% crown points.
%
% Input: - rootTLSData   % Root folder for the point cloud data [character array]
%        - rootResults   % Root folder where to store the results [character array]
%
% Output:
%        - TLS stem points for each tree (.txt)
%        - TLS & ALS crown points for each tree (.txt)
%
% (c) Tuomas Yrttimaa / Science4Trees @ UEF School of Forest Sciences 2021
% ---------------------------------------------------------------------------

function pc_classify_v2(rootTLSData,rootResults)
% --------------------------------------------------------------------------------------------    
    % starttime = datetime;
    % Create folder for the results (if not already exist)
    if isfolder(rootResults) == false, mkdir(rootResults); end
    currdir = cd; paramdir = [currdir,'/tmp/param'];
    if isfolder(paramdir) == false, mkdir(paramdir);end
    f = fopen([paramdir,'/numtrees.txt'], 'w');fprintf(f, '%d\n', 0);fclose(f);
    list = struct2table(dir(fullfile([rootTLSData,'/tmp'], '**', '*.las')));
    filenames = list.name;
    parfor_progress(size(filenames,1));
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
    
    tic
    parfor i = 1:size(filenames,1)
        try
            PCtls = lasdata([rootTLSData,'/tmp/',filenames{i}]);
            PCtls = pointCloud(table2array(table(PCtls.x,PCtls.y,PCtls.z)));
            
            warning('off','vision:ransac:maxTrialsReached');warning('off','MATLAB:polyshape:repairedBySimplify');warning('off','MATLAB:alphaShape:DupPointsBasicWarnId'); warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
            parfor_progress;
            %--------------------------------------------------------   
            % 1. DETECT TREE STEMS
            %--------------------------------------------------------
            try
                % Delineate point cloud to detect tree stems (1 m < z < 5 m)
                pc = pcdownsample(pointCloud(PCtls.Location(PCtls.Location(:,3)<= 4,:)),'gridAverage',0.005);
                pc = pointCloud(pc.Location(pc.Location(:,3) >= 1,:));
                pc = pcdenoise(pc,'NumNeighbors',20);

                % Compute surface normals and roughly filter out points on
                % non-vertical surfaces
                normals = pcnormals(pc,30);
                pc = pointCloud(pc.Location(abs(normals(:,3)) < 0.010,:));
                pc = pcdenoise(pc,'NumNeighbors',30);
            catch
                continue
            end

            % Segment the point cloud into clusters and keep only clusters that
            % have enough points and vertical dimensions.
            pctrees = cell([10,1]);  % clusters representing different trees are stored here
            minDistance = [0.3 0.5]; % use 2 different parameters for clustering
            count = 0;
            for k = minDistance
                count = count + 1;
                pc = pcdenoise(pc,'NumNeighbors',10);
                [labels,numClusters] = pcsegdist(pc,k);
                clusterPC = cat(2,pc.Location,double(labels));
                % figure, pcshow(pc.Location,labels), colormap(hsv(numClusters));
                inliers = [];
                for j = 1:numClusters
                    pcCluster = pointCloud(clusterPC(clusterPC(:,4) == j,1:3));
                    pcCluster = pcdenoise(pcCluster);

                    if (pcCluster.Count >= 100 && abs(diff([pcCluster.ZLimits(1) pcCluster.ZLimits(2)])) > 1) == true
                        inliers = cat(1,inliers,pcCluster.Location);
                        if count == size(minDistance,2)
                            pctrees{j} = pcCluster;
                        end
                    else
                        continue
                    end
                end
                if ~isempty(inliers)
                    pc = pointCloud(inliers);
                else
                    continue
                end
            end
            pctrees = pctrees(~cellfun('isempty',pctrees));
            % 'pctrees' now includes candidate stem point clusters that each
            % represent individual tree stem. If 'pctrees' is empty, the crown
            % segment contains no eligible point cloud structures for tree
            % detection.

            % Next, the candidate stem point clusters are filtered to confirm
            % tree detection by applying RANSAC-cylinder filtering procedure.
            if ~isempty(pctrees)
                cylmodels = cell([size(pctrees,1),1]);
                cylcents = [];
                    % figure, hold on;
                for j = 1:size(pctrees,1)
                    % The RANSAC-cylinder filtering procedure is repreated for
                    % robustness. Vertical reference vector is used to guide
                    % the procedure to fit vertical cylinders.
                    iter = 100;
                    maxDistance = 0.02;
                    referenceVector = [0 0 1];
                    maxNumTrials = 5000;
                    paramlist = zeros([iter,7]);
                    inlierlist = zeros([iter,1]);
                    pc = pctrees{j}; pc.Normal = pcnormals(pc,20);
                    for k = 1:iter
                        try
                            [model,inlierIndices] = pcfitcylinder(pc,maxDistance,referenceVector,maxNumTrials);
                            % Number of inlier points is stored for later use.
                            inliers = select(pc,inlierIndices);
                            inlierlist(k) = inliers.Count;
                            % Cylinder parameters for each iteration are stored for later use.
                            paramlist(k,:) = model.Parameters;
                        catch
                            continue
                        end
                    end
                    % Radius of the fitted cylinders are compared and clear
                    % outliers are removed.
                    TF = isoutlier(paramlist(:,7));
                    paramlist(TF,:) = []; inlierlist(TF,:) = [];
                    % If there is enough inlier points, average cylinder model
                    % parameters are computed and the cylinder model object is
                    % stored for later use.
                    if mean(inlierlist) > 100
                        cylmodel = cylinderModel(mean(paramlist));
                        cylmodels{j} = cylmodel;
                        cylcents = cat(1,cylcents,cylmodel.Center);
                            % plot(cylmodel); plot3(pctrees{j}.Location(:,1),pctrees{j}.Location(:,2),pctrees{j}.Location(:,3),'r.');
                            % xlim([cylmodel.Center(1)-5 cylmodel.Center(1)+5]);ylim([cylmodel.Center(2)-5 cylmodel.Center(2)+5]);
                    else
                        continue
                    end
                end
                    % hold off;
                cylmodels = cylmodels(~cellfun('isempty',cylmodels));
                if ~isempty(cylmodels)
                    % Now 'cylmodels' contain cylinder models that represent each
                    % detected tree stem inside the CHM-derived crown segment.

                    % If multiple trees have been detected from one crown segment,
                    % the crown segment is split in parts.

                    % First generate a set of grid points to fill the crown segment.
                    convhulpoints = select(PCtls,convhull(PCtls.Location(:,1:2)));
                    crownseg = polyshape(convhulpoints.Location(:,1),convhulpoints.Location(:,2));
                    [Xlim,Ylim] = boundingbox(crownseg);
                    [xq,yq] = meshgrid(Xlim(1):0.1:Xlim(2),Ylim(1):0.1:Ylim(2));
                    xlist = []; ylist = [];
                    for k = 1:size(xq,2)
                        xlist = cat(1,xlist,xq(:,k)); ylist = cat(1,ylist,yq(:,k));
                    end
                    in = inpolygon(xlist,ylist,crownseg.Vertices(:,1),crownseg.Vertices(:,2));
                    gridpoints = cat(2,xlist,ylist); gridpoints = gridpoints(in,:);

                    % Group the gridpoints by finding closest tree for each gridpoint
                    treesegcell = cell(size(cylmodels));
                    for k = 1:size(gridpoints,1)
                        point = gridpoints(k,:);
                        trees = cat(2,cylcents(:,1),cylcents(:,2));
                        kdtree = ExhaustiveSearcher(trees);
                        nn = knnsearch(kdtree,point);
                        temp = treesegcell{nn}; temp = cat(1,temp,point); treesegcell{nn} = temp;
                    end

                    % Generate split crown segments from the convex hull of the grouped gridpoints
                    treecrownsegments = cell(size(cylmodels));
                    for k = 1:size(treesegcell,1)
                        K = convhull(treesegcell{k}(:,1),treesegcell{k}(:,2), 'simplify',true);
                        treecrownsegments{k} = polyshape(treesegcell{k}(K,1),treesegcell{k}(K,2));
                    end
                    

                    %----------------------------------------------------------  
                    % 2. CLASSIFY POINT CLOUD INTO STEM POINTS AND CROWN POINTS
                    %---------------------------------------------------------- 
                    % Once the lower part of the tree stem (up to 5 m) has detected,
                    % the point cloud is classified into stem points and crown
                    % points by using an iterative cylinder filtering approach.
                    for j = 1:size(cylmodels,1)
                        % Stem is characterized as a series of cylinders and the
                        % parameters of the cylinders are stored in 'treemodelparamlist'
                        treemodelparamlist = [];
                        cylmodelpar = cylmodels{j}.Parameters;
                        cylmodelpar(3) = 0; cylmodel = cylinderModel(cylmodelpar);
                        treemodelparamlist = cat(1,treemodelparamlist,cylmodel.Parameters);

                        % Crown segments are used to delineate tree-wise point
                        % cloud with a point-in-polygon approach
                        crownseg = treecrownsegments{j};
                        pc = PCtls; 
                        if ~isempty(pc.Location)
                            pc = pointCloud(pc.Location(inpolygon(pc.Location(:,1),pc.Location(:,2),crownseg.Vertices(:,1),crownseg.Vertices(:,2)),:));
                            if ~isempty(pc.Location)
                                % Points within the surface of the first cylinder (z < 4 m)
                                % are delineated and classified as stem points
                                inPoints = []; htresh = cylmodel.Parameters(6); hmax = pc.ZLimits(2);
                                tmp = pc.Location(pc.Location(:,3) < htresh,:);
                                flagBelowVec = (tmp(:,3) <= cylmodel.Parameters(6));
                                flagAboveVec = (tmp(:,3) >= cylmodel.Parameters(3));
                                radialDistanceSquaredVec = (tmp(:,1)-mean([cylmodel.Parameters(1) cylmodel.Parameters(4)])).^2 + (tmp(:,2)-mean([cylmodel.Parameters(2) cylmodel.Parameters(5)])).^2;
                                flagInsideVec = (radialDistanceSquaredVec <= 1.7*cylmodel.Parameters(7)^2);
                                flagIsIn = (flagBelowVec & flagAboveVec & flagInsideVec);

                                tmp = tmp(flagIsIn,:);
                                inPoints = cat(1,inPoints,tmp);
                                if ~isempty(inPoints)
                                    inPoints = pointCloud(inPoints);
                                end

                                % Next, the remaining point cloud above the 5 m treshold
                                % is split into horizontal slices of 50 cm high.
                                hstep = 4:0.5:ceil(hmax);

                                % From each slice of points, cylindrical stem surfaces are
                                % searched by applying a guided cylinder filtering
                                % procedure where information on the previous cylinder is
                                % used for more robust RANSAC cylinder fit.

                                % 'infotab' contains information on the previous cylinder
                                % that is used as a guide for the next one. It is also
                                % investigated whether the stem diverges or not, which
                                % means that there might be more than one 'guide cylinders'
                                infotab = table('Size',[length(hstep),2],'variableNames',{'divnum','cylcell'},'variableTypes',{'double','cell'});
                                infotab.cylcell = cell([length(hstep),5]);
                                infotab.divnum = ones([length(hstep),1]); % No stem diverge by default
                                maxdiv = 3;
                                tmpcell = cell([1,1]); tmpcell{1} = cylmodel;
                                infotab.cylcell(1,1) = tmpcell;
                                noCylinder = 0;

                                for k = 1:length(hstep)-1
                                    % First get the guide information
                                    numdiverges = infotab.divnum(k); % number of diverging stems at this slice
                                    numdiverges(numdiverges > maxdiv) = maxdiv;
                                    cyls = infotab.cylcell(k,:);     % guide cylinders
                                    cyls = cyls(~cellfun('isempty',cyls));
                                    cylcount = 0;
                                    for q = 1:numdiverges
                                        % Extract the point cloud slice
                                        lowlim = hstep(k); highlim = hstep(k+1);
                                        tmp = pointCloud(pc.Location(pc.Location(:,3) >= lowlim & pc.Location(:,3) < highlim,:));

                                        if ~isempty(tmp.Location)
                                            % Downsample to 5 mm grid
                                            tmp = pcdownsample(tmp,'gridAverage',0.005);
                                            % Use the previous cylinder as a guide. If no cylindrical structures 
                                            % were detected at the previous slice, use the stem origin as a guide.
                                            try
                                                prevcyl = cyls{q};
                                                noCylinder = 0;
                                            catch
                                                noCylinder = noCylinder + 1;
                                                if noCylinder < 2
                                                    try
                                                        prevcyl = clustcylmodel;
                                                    catch
                                                        continue
                                                    end
                                                else
                                                    continue
                                                end
                                            end

                                            % Define region of interest (ROI) based on the XY-location and radius of the 
                                            % previous cylinder. Use more tolerance for the first slice.
                                            if k == 1
                                                roi = [prevcyl.Center(1)-2*prevcyl.Radius prevcyl.Center(1)+2*prevcyl.Radius prevcyl.Center(2)-2*prevcyl.Radius prevcyl.Center(2)+2*prevcyl.Radius lowlim highlim];
                                            else
                                                roi = [prevcyl.Center(1)-1.5*prevcyl.Radius prevcyl.Center(1)+1.5*prevcyl.Radius prevcyl.Center(2)-1.5*prevcyl.Radius prevcyl.Center(2)+1.5*prevcyl.Radius lowlim highlim];
                                            end

                                            % Select points within ROI
                                            try
                                                tmp = select(tmp,findPointsInROI(tmp,roi));
                                            catch
                                                continue
                                            end

                                            % Cluster the point cloud slice to filter out non-vertical surfaces
                                            if ~isempty(tmp.Location)
                                                normals = pcnormals(tmp,40);
                                                tmp = pointCloud(tmp.Location(abs(normals(:,3)) < 0.05,:));

                                                if ~isempty(tmp.Location)
                                                    tmp = pcdenoise(tmp,'NumNeighbors',40);
                                                    [labels,numClusters] = pcsegdist(tmp,0.3);
                                                    clusterPC = cat(2,tmp.Location,double(labels));

                                                    for m = 1:numClusters
                                                        clusterPoints = clusterPC(clusterPC(:,4) == m,:);
                                                        pcCluster = pointCloud(clusterPoints(:,1:3));
                                                        pcCluster.Normal = pcnormals(pcCluster,20);

                                                        % Apply RANSAC-cylinder filtering. Store cylinder parameters for later use.
                                                        if pcCluster.Count >= 40
                                                            iter = 300; maxDistance = 0.02;
                                                            referenceVector = [0 0 1]; maxNumTrials = 5000;
                                                            paramlist = zeros([iter,7]);inlierlist = zeros([iter,1]);

                                                            for n = 1:iter
                                                                try
                                                                    [model,inlierIndices] = pcfitcylinder(pcCluster,maxDistance,referenceVector,maxNumTrials);
                                                                    inliers = select(pcCluster,inlierIndices);
                                                                    inlierlist(n) = inliers.Count;
                                                                    paramlist(n,:) = model.Parameters;
                                                                catch
                                                                    continue
                                                                end
                                                            end

                                                            % Remove clear outliers
                                                            try
                                                                TF = isoutlier(paramlist(:,7));
                                                                paramlist(TF,:) = []; inlierlist(TF,:) = [];
                                                            catch
                                                                continue
                                                            end

                                                            % Compute average cylinder parameters
                                                            param = mean(paramlist);
                                                            % Fix cylinder height
                                                            param(3) = lowlim; param(6) = highlim;
                                                            % Take stem diverge into account by extending the cylinder downwards
                                                            if (k > 1 && numdiverges > infotab.divnum(k-1)) == true
                                                                param(3) = lowlim - 0.5;
                                                                roi(5) = lowlim - 0.5;
                                                            end
                                                            clustcylmodel = cylinderModel(param);

                                                            % Check that the fitted cylinder meets the conditions
                                                            mod_ori = [0 0 abs(clustcylmodel.Orientation(3))];
                                                            theta = acosd(dot(mod_ori,[0 0 1])/(norm(mod_ori)*norm([0 0 1])));
                                                            locdiff = sqrt((clustcylmodel.Center(1)-prevcyl.Center(1))^2+(clustcylmodel.Center(2)-prevcyl.Center(2))^2);

                                                            if (abs(theta) < 25 && clustcylmodel.Radius >= 0.015 && clustcylmodel.Radius <= 0.35 && mean(inlierlist) > 30 && locdiff < 0.3) == true
                                                                % Store cylinder parameters
                                                                treemodelparamlist = cat(1,treemodelparamlist,clustcylmodel.Parameters);

                                                                % Investigate if the stem diverges. Delineate a 15 cm high point cloud slice
                                                                % on top of the 50 cm slice and segment it into clusters with 5 cm apart from each other.
                                                                topslice = pointCloud(pcCluster.Location(pcCluster.Location(:,3) >= pcCluster.ZLimits(2) - 0.15,:));
                                                                if ~isempty(topslice.Location)
                                                                    try
                                                                        topslice = pcdenoise(topslice);
                                                                        [labelstop,numClusterstop] = pcsegdist(topslice,0.05);
                                                                        clusterPCtop = cat(2,topslice.Location,double(labelstop));
                                                                    catch
                                                                        continue
                                                                    end

                                                                    % Ensure that the clusters are cylincrical using the RANSAC-cylinder filtering procedure
                                                                    topcylmodelclust = cell([1,numClusterstop]);
                                                                    for p = 1:numClusterstop
                                                                        clusterPoints = clusterPCtop(clusterPCtop(:,4) == p,:);
                                                                        try
                                                                            pcCluster = pointCloud(clusterPoints(:,1:3));
                                                                            pcCluster.Normal = pcnormals(pcCluster,20);
                                                                        catch
                                                                            continue
                                                                        end
                                                                        if pcCluster.Count >= 20
                                                                            iter = 300; maxDistance = 0.02;
                                                                            referenceVector = [0 0 1]; maxNumTrials = 5000;
                                                                            paramlistclust = zeros([iter,7]);inlierlistclust = zeros([iter,1]);
                                                                            for n = 1:iter
                                                                                try
                                                                                    [model,inlierIndices] = pcfitcylinder(pcCluster,maxDistance,referenceVector,maxNumTrials);
                                                                                    inliers = select(pcCluster,inlierIndices);
                                                                                    inlierlistclust(n) = inliers.Count;
                                                                                    paramlistclust(n,:) = model.Parameters;
                                                                                catch
                                                                                    continue
                                                                                end
                                                                            end
                                                                        end
                                                                        % Remove outliers
                                                                        try
                                                                            TF = isoutlier(paramlistclust(:,7));
                                                                            paramlistclust(TF,:) = []; inlierlistclust(TF,:) = [];
                                                                        catch
                                                                            continue
                                                                        end

                                                                        % Compute average cylinder model parameters.
                                                                        cyl = cylinderModel(mean(paramlistclust));
                                                                        locdiff = sqrt((cyl.Center(1)-prevcyl.Center(1))^2+(cyl.Center(2)-prevcyl.Center(2))^2);
                                                                        mod_ori1 = cyl.Orientation; mod_ori2 = clustcylmodel.Orientation;
                                                                        theta = acosd(dot(mod_ori1,mod_ori2)/(norm(mod_ori1)*norm(mod_ori2)));

                                                                        % Check that the fitted cylinder meets the conditions, and save the cylindermodel
                                                                        if (abs(theta) <= 35 && cyl.Radius >= 0.02 && (cyl.Radius <= 1.3*prevcyl.Radius || cyl.Radius <= 0.10) && locdiff <= 2.5*prevcyl.Radius && mean(inlierlistclust) > 20) == true
                                                                            topcylmodelclust{p} = cyl;
                                                                        end
                                                                    end
                                                                    topcylmodelclust = topcylmodelclust(~cellfun('isempty',topcylmodelclust));
                                                                    divnum = size(topcylmodelclust,2);

                                                                    % Store info to guide the stem points classification in the next slice
                                                                    if divnum > 1
                                                                        cylcount = cylcount + divnum;
                                                                        infotab.cylcell(k+1,cylcount:(cylcount+divnum-1)) = topcylmodelclust;
                                                                    else
                                                                        cylcount = cylcount + 1;
                                                                        temp = cell([1,1]); temp{1} = clustcylmodel;
                                                                        infotab.cylcell(k+1,cylcount) = temp;
                                                                    end
                                                                    temp = infotab.cylcell(k+1,:);
                                                                    infotab.divnum(k+1) = length(temp(~cellfun('isempty',temp)));

                                                                    % Delineate stem points that are located within the cylinder surface
                                                                    try
                                                                        tmp = select(pc,findPointsInROI(pc,roi));
                                                                        flagBelowVec = (tmp.Location(:,3) <= clustcylmodel.Parameters(6));
                                                                        flagAboveVec = (tmp.Location(:,3) >= clustcylmodel.Parameters(3));
                                                                        radialDistanceSquaredVec = (tmp.Location(:,1)-mean([clustcylmodel.Parameters(1) clustcylmodel.Parameters(4)])).^2 + (tmp.Location(:,2)-mean([clustcylmodel.Parameters(2) clustcylmodel.Parameters(5)])).^2;
                                                                        flagInsideVec = (radialDistanceSquaredVec <= 1.2*clustcylmodel.Parameters(7)^2);
                                                                        flagIsIn = (flagBelowVec & flagAboveVec & flagInsideVec);
                                                                        tmp = tmp.Location(flagIsIn,:);
                                                                        inPoints = pointCloud(cat(1,inPoints.Location,tmp));
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
                                end

                                % Extract tree number from a txt file
                                f = fopen([paramdir,'/numtrees.txt'], 'a');fprintf(f, '1\n');fclose(f);
                                f = fopen([paramdir,'/numtrees.txt'], 'r');progress = fscanf(f, '%d');fclose(f);
                                numTrees = (length(progress)-1);

                                % Save the cylinder parameters
                                savedir = [rootResults,'/treecylinders'];
                                if isfolder(savedir) == false, mkdir(savedir);end
                                writematrix(treemodelparamlist,sprintf([savedir,'/cylinderparameters_tree%i.txt'],numTrees),'Delimiter','\t');

                                % Save stem points
                                if ~isempty(inPoints.Location)
                                    savedir = strcat(rootResults,'/pointclouds/stem');
                                    if isfolder(savedir) == false; mkdir(savedir);end
                                    inPoints.Intensity = uint8(numTrees.*ones([inPoints.Count,1]));
                                    pcwrite(inPoints,strcat(savedir,sprintf('/tree%i_stempoints.ply',numTrees)))
                                end 

                                % Cluster the point cloud for crown points extraction
                                try
                                    [labels,numClusters] = pcsegdist(pc,0.7);
                                    %figure, pcshow(pc.Location,labels), colormap(hsv(numClusters));
                                    limits = zeros([numClusters,6]);
                                    guidecyl = infotab.cylcell{1};
                                    for k = 1:numClusters
                                       tmpc = pointCloud(pc.Location(labels==k,:));
                                       horiz_dev = sqrt((guidecyl.Center(1)-mean(tmpc.Location(:,1)))^2+(guidecyl.Center(2)-mean(tmpc.Location(:,2)))^2);
                                       limits(k,:) = [tmpc.ZLimits tmpc.Count abs(tmpc.ZLimits(2)-tmpc.ZLimits(1)) k horiz_dev];
                                    end
                                    limits(limits(:,3) < 1000,:) = [];
                                    limits(limits(:,1) > 2,:) = [];
                                    limits(limits(:,6) > 1,:) = [];
                                    limits = sortrows(limits,4,'descend');
                                    crownpoints = pointCloud(pc.Location(labels==limits(1,5),:)); 
                                catch
                                    crownpoints = pc;
                                end
                                
                                % Extract non-stem points
                                outPoints = pcdenoise(pointCloud(crownpoints.Location(logical(abs(inShape(alphaShape(inPoints.Location),crownpoints.Location)-1)),:)),'NumNeighbors',10);

                                % Construct alphashape to filter out non-stem points from undergrowth vegetation
                                stemloc = treemodelparamlist(1,1:2);
                                t = 0:pi/10:2*pi;
                                [X,Y,Z] = cylinder(0.1+0.4*t.^1.6);
                                xlist = []; ylist = []; zlist = [];
                                for k = 1:size(X,2)
                                    xlist = cat(1,xlist,X(:,k)); ylist = cat(1,ylist,Y(:,k));zlist = cat(1,zlist,Z(:,k));
                                end
                                hmax = max(hmax,pc.ZLimits(2));
                                xlist = xlist + stemloc(1); ylist = ylist + stemloc(2); zlist = zlist * hmax;

                                crownAlpha = alphaShape(xlist,ylist,zlist,15);
                                % hold on; plot(crownAlpha, 'FaceAlpha',0.1,'FaceColor','g','EdgeColor','r')
                                
                                % Save crown points
                                try
                                    canopypoints = pointCloud(outPoints.Location(logical(abs(inShape(crownAlpha,outPoints.Location))),:));
                                catch
                                    warning('No crown points for this tree!');
                                end
                                
                                
                                % View preliminary point classification
                                % figure, pcshowpair(inPoints,canopypointsTLS), view(90,0);
                                % hold on, plot(crownAlpha, 'FaceAlpha',0.1,'FaceColor','g','EdgeColor','r'), pcshow(canopypointsALS),hold off;
                                % set(gca,'color','w');set(gcf,'color','w');
                                % set(gca, 'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'ZColor', [0.15 0.15 0.15])
                                
                                try
                                    inPoints.Color = uint8([170 115 30] .* ones(size(inPoints.Location)));
                                    canopypoints.Color = uint8([110 190 100] .* ones(size(canopypoints.Location)));
                                    savedir = [rootResults,'/figures/pointclassification'];
                                    if isfolder(savedir) == false; mkdir(savedir);end
                                    % Illustrate point cloud classification
                                    myfunctions.pcshow_parfor(inPoints,canopypoints)
                                    filename = sprintf([savedir,sprintf('/pc_treeID_%i.png',numTrees)]);
                                    saveas(gcf,filename);
                                    close all
                                catch
                                    warning('Cannot illustrate point cloud for tree no. %i!',numTrees);
                                end

                                %winopen(filename);

                                

                                % Save crown points
                                if ~isempty(canopypoints.Location)
                                    savedir = strcat(rootResults,'/pointclouds/canopy');
                                    if isfolder(savedir) == false; mkdir(savedir); end
                                    canopypoints.Intensity = uint8(numTrees.*ones([canopypoints.Count,1]));
                                    pcwrite(canopypoints,strcat(savedir,sprintf('/tree%i_canopypoints_tls.ply',numTrees)))
                                end
                            end
                        end
                    end
                end
            end
        catch
            continue
        end
    end
    toc
    parfor_progress(0);
end
