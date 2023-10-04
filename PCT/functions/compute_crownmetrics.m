function crownMetrics = compute_crownmetrics(tree_id,rootResults)
    warning('off','MATLAB:polyshape:repairedBySimplify');

    pctls = pcread(strcat(rootResults,sprintf('/%i/pointclouds/canopy/%i_tree%i_canopypoints_tls.ply',plot_id,plot_id,tree_id)));
    pctls = pointCloud(pctls.Location(pctls.Location(:,3)>0.5,:));
    
    if pctls.Count == 0
        warning('Cannot compute crown metrics for tree %i',tree_id)
    else
        crownMetrics.tree_id = tree_id;
        
        pc_tls_top = pointCloud(pctls.Location(pctls.Location(:,3) > 0.7 * pctls.ZLimits(2),:));
        pc_tls_bottom = pointCloud(pctls.Location(pctls.Location(:,3) <= 0.7 * pctls.ZLimits(2),:));
        theta = (0:50-1)*(2*pi/50);
        try
            pcstem = pcread(strcat(rootResults,sprintf('/pointclouds/stem/tree%i_stempoints.ply',tree_id)));
            xy_base = mean(pcstem.Location(pcstem.Location(:,3) < pcstem.ZLimits(1)+1.5,1:2));
            crownMetrics.X = xy_base(1);
            crownMetrics.Y = xy_base(2);
        catch
            xy_base = mean(pctls.Location(pctls.Location(:,3) < pctls.ZLimits(1)+1.5,1:2));
            crownMetrics.X = xy_base(1);
            crownMetrics.Y = xy_base(2);
        end
        x = xy_base(1) + 0.2*cos(theta);
        y = xy_base(2) + 0.2*sin(theta);
        pc = pointCloud(cat(1,pc_tls_top.Location,(pc_tls_bottom.Location(not(inpolygon(pc_tls_bottom.Location(:,1),pc_tls_bottom.Location(:,2),x,y)),:))));            
            
        %figure, pcshowpair(pc_tls_out,pc_als_out)
        heights = (pc.ZLimits(1):0.2:pc.ZLimits(2));
        features = zeros([size(heights,1),5]);
        for j = 1:(length(heights)-1)
            bottom = heights(j);
            top = heights(j+1);
            slice = pc.Location(pc.Location(:,3)<= top & pc.Location(:,3)>= bottom,1:2);
            features(j,1) = bottom;
            features(j,2) = size(slice,1);
            if size(slice,1) > 20
                poly = polyshape(slice(convhull(slice),:));
                features(j,3) = poly.area;
                features(j,4) = poly.perimeter;
                features(j,5) = poly.perimeter/poly.area;
            else
                continue
            end
        end

        features(:,6) = smooth(features(:,5));

        if pc.ZLimits(2) < 14
            min_area = 2;
            max_perimeter2area_ratio = 3;
        else
            min_area = 5;
            max_perimeter2area_ratio = 2;
        end
        
        features = features(features(:,3) >= min_area & features(:,6) < max_perimeter2area_ratio,:);
        if ~isempty(features)
            crownBaseHeight = min(features(:,1));
        else
            crownBaseHeight = 1.3;
        end
        pc = pointCloud(pc.Location(pc.Location(:,3) >= crownBaseHeight,:));
        [~,vol] = convhull(pc.Location);
        [~,surf_area] = area3d(pc.Location);
        
        voxel = pcdownsample(pc,'gridAverage',0.1);
        
        [chull2d,area] = convhull(pc.Location(:,1:2));
        
        poly = polyshape(pc.Location(chull2d,1:2));
        crownMetrics.vol_convhull = vol;
        crownMetrics.area3d = surf_area;
        crownMetrics.vol_voxel = voxel.Count/1000;
        crownMetrics.area2d = area;
        crownMetrics.crownBaseHeight = crownBaseHeight;
        crownMetrics.perimeter = poly.perimeter; 
        
        savedir = strcat(rootResults,'/crownmetrics');
        if isfolder(savedir) == false; mkdir(savedir);end
        writetable(crownMetrics,sprintf(strcat(savedir,'/crownmetrics_%i.csv'),tree_id),'Delimiter',';');
        
        try
            savedir = strcat(rootResults,'/crownpolygons');
            if isfolder(savedir) == false; mkdir(savedir);end
            shape.Geometry = 'Polygon';
            shape.X = cat(1,poly.Vertices(:,1),poly.Vertices(1,1));
            shape.Y = cat(1,poly.Vertices(:,2),poly.Vertices(1,2));
            shapewrite(shape,sprintf(strcat(savedir,'/crownpolygon_tree%i.shp'),tree_id));
        catch
            warning('cannot extract crown polygon for tree %i',tree_id);
        end
    end
end
