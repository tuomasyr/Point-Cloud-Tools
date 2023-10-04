function treePgons = pc_detect_tree_crowns_v2(ptCloud,gridRes,minCrownArea,minTreeHeight,viewplot,rootResults) 
    
    warning('off','MATLAB:polyshape:repairedBySimplify');
   
    normalizedPoints = ptCloud.Location;
    % Visualize normalized points
    if viewplot == true; figure; pcshow(normalizedPoints); title("Point Cloud with Normalized Elevation"); end
    ptCloud = pointCloud(normalizedPoints);

    % Generate CHM
    [canopyModel,~,~] = pc2dem(pointCloud(normalizedPoints),gridRes,CornerFillMethod="max");
    % Clip invalid and negative CHM values to zero
    canopyModel(isnan(canopyModel) | canopyModel<0) = 0;
    % Perform gaussian smoothing to remove noise effects
    H = fspecial("gaussian",[3 3],1);
    canopyModel = imfilter(canopyModel,H,'replicate','same');
    % Visualize CHM
    if viewplot == true; figure; imagesc(canopyModel); title('Canopy Height Model'); axis off; colormap(parula); end
    
    % Detect tree tops
    [treeTopRowId,treeTopColId] = helperDetectTreeTops(canopyModel,gridRes,minTreeHeight);
    % Visualize treetops
     if viewplot == true; figure; imagesc(canopyModel);hold on;plot(treeTopColId,treeTopRowId,"rx",MarkerSize=3);title("CHM with Detected Tree Tops");axis off; colormap("gray"); end
    
    %treecoords = cat(2,xlimits(1) + gridRes*treeTopColId,ylimits(1) + gridRes*treeTopRowId);
    
    % Segment individual trees
    label2D = helperSegmentTrees(canopyModel,treeTopRowId,treeTopColId,minTreeHeight);
    % Identify row and column id of each point in label2D and transfer labels
    % to each point
    rowId = ceil((ptCloud.Location(:,2) - ptCloud.YLimits(1))/gridRes) + 1;
    colId = ceil((ptCloud.Location(:,1) - ptCloud.XLimits(1))/gridRes) + 1;
    ind = sub2ind(size(label2D),rowId,colId);
    label3D = label2D(ind);
    % Extract valid labels and corresponding points
    validSegIds = label3D ~= 0;
    ptVeg = select(ptCloud,validSegIds);
    veglabel3D = label3D(validSegIds);
    % Assign color to each label
    numColors = max(veglabel3D);
    colorMap = randi([0 255],numColors,3)/255;
    labelColors = label2rgb(veglabel3D,colorMap,OutputFormat="triplets");
    % Visualize tree segments
    if viewplot == true; figure; pcshow(ptVeg.Location,labelColors); title("Individual Tree Segments"); view(3); end
    
    % Extract crown segments
    labels = unique(veglabel3D);
    treePgons = cell(2*size(labels,1),1);
    numTrees = 0;
    %treeLocations = zeros(2*size(treePgons,1),2);
    savedir = strcat(rootResults,'/treepoly'); if isfolder(savedir) == false, mkdir(savedir); end
    for id = transpose(labels)
        treecloud = double(ptVeg.Location(veglabel3D == id,1:2));
        if size(treecloud,1) > 10
            try
                R = regions(rmholes(polyshape(treecloud(convhull(treecloud),:), "Simplify",false)));
                for j = 1:size(R,1)
                    numTrees = numTrees + 1;
                    if area(R(j)) >= minCrownArea
                        poly = R(j);
                        writematrix(round(cat(1,poly.Vertices,poly.Vertices(1,:)),2),strcat(savedir,sprintf('/treepoly_%i.txt',numTrees)),'Delimiter','\t');
                        pause(1);
                        treePgons{numTrees} = convhull(R(j));
                        %[treeLocations(numTrees,1),treeLocations(numTrees,2)] = centroid(convhull(poly));
                    else
                        continue
                    end
                end
            catch
                continue
            end
        else
            continue
        end
    end
    treePgons = treePgons(~cellfun('isempty',treePgons));
    savedir = strcat(rootResults,'/treepoly'); if isfolder(savedir) == false, mkdir(savedir); end
    save(strcat(savedir,'/treePgons.mat'),'treePgons')
    for treeid = 1:max(size(treePgons))
       poly = treePgons{treeid};
       writematrix(round(cat(1,poly.Vertices,poly.Vertices(1,:)),2),strcat(savedir,sprintf('/treepoly_%i.txt',treeid)),'Delimiter','\t');
    end
end 