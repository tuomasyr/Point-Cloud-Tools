
% Compute tree metrics for each tree

function treemetr = compute_treemetrics(tree_id,rootResults,h_butt,h_interval,h_layer,min_points,p)
    currdir = cd; paramdir = [currdir,'\tmp\param'];
    if isfolder(paramdir) == false
        mkdir(paramdir)
    end
    proceed = false; trials = 0;
    warning('off','MATLAB:polyshape:repairedBySimplify');
    warning('off','signal:findpeaks:largeMinPeakHeight');
    warning('off','MATLAB:detrend:PolyNotUnique');
    treecylinder = myfunctions.import_treecylinder(sprintf([rootResults,'/treecylinders/cylinderparameters_tree%i.txt'],tree_id));
    prelimdiam = 2*treecylinder.diam(1)*100;
    treeloc = [treecylinder.x_bottom(1) treecylinder.y_bottom(1)];
    while proceed == false && trials < 6
        treemetr = array2table(zeros(1,11));
        treemetr.Properties.VariableNames = {'treeID','X','Y','dbh','ba','h','vol','d6','d05','vlog','vpulp'};
        treemetr.treeID = tree_id;

        % Import point cloud data
        try
            pcStem = myfunctions.importclassifiedpcfile(rootResults,tree_id,'stem');
            pcCanopy = myfunctions.importclassifiedpcfile(rootResults,tree_id,'canopy');
        catch
            warning('Cannot load TLS points for tree no. %i!', tree_id);
            return
        end

        if prelimdiam < 14
            poly = polybuffer(treeloc,'points',1);
            pcCanopy = pointCloud(pcCanopy.Location(inpolygon(pcCanopy.Location(:,1),pcCanopy.Location(:,2),poly.Vertices(:,1),poly.Vertices(:,2)),:));
        end
        hmax = max(pcCanopy.ZLimits); hmin = max(pcStem.ZLimits);
        hrange = hmin:0.75:hmax;
        crownskeleton = zeros([size(hrange,2),8]);
        for ind = 1:size(hrange,2)
            try
                h = hrange(ind);
                hlayer = pcCanopy.Location(pcCanopy.Location(:,3) > h & pcCanopy.Location(:,3) < h+0.75,:);
                if size(hlayer,1) > 9
                    convexh = polyshape(hlayer(convhull(hlayer(:,1),hlayer(:,2)),1:2));
                    [centx,centy] = convexh.centroid;
                    if ind > 1
                        distprev = sqrt((centx-prev(1))^2+(centy-prev(2))^2);
                        distbottom = sqrt((centx-treeloc(1))^2+(centy-treeloc(2))^2);
                        crownskeleton(ind,:) = [centx centy h size(hlayer,1) convexh.area convexh.perimeter distprev distbottom];
                    else
                        distbottom = sqrt((centx-treeloc(1))^2+(centy-treeloc(2))^2);
                        crownskeleton(ind,:) = [centx centy h size(hlayer,1) convexh.area convexh.perimeter 0 distbottom];
                    end
                    prev = [centx,centy];
                end
            catch
                continue
            end
        end
        crownskeleton = crownskeleton(crownskeleton(:,3) > 0,:);
        crownskeleton(:,7) = smoothdata(detrend(crownskeleton(:,7)),'smoothingfactor',0.10);
        crownskeleton(:,8) = smoothdata(detrend(crownskeleton(:,8)),'smoothingfactor',0.10);
        if prelimdiam < 14
            maxdif = 0.2;
        else
            maxdif = 0.5;
        end
        try
            
            [~,locs] = findpeaks(crownskeleton(:,7),'MinPeakHeight',maxdif);
            
            if ~isempty(locs)
                hEnd = crownskeleton(min(locs),3);
                pcCanopy = pointCloud(pcCanopy.Location(pcCanopy.Location(:,3) < hEnd,:));
                %crownskeleton = crownskeleton(crownskeleton(:,3) <= hend,:);
            end
        catch
            hEnd = max([pcStem.ZLimits(2) pcCanopy.ZLimits(2)]);
        end
        
        
        
        hEnd = max([pcStem.ZLimits(2) pcCanopy.ZLimits(2)]);
        if hEnd <= 4
            treemetr.h = hEnd;
        else
            starth = prctile(pcCanopy.Location(pcCanopy.Location(:,3) > 3,3),50);
            klist = (max([starth pcCanopy.ZLimits(1)]):0.75:hEnd);
            k = 1;
            while k <= (size(klist,2))
                hlayer = pointCloud(pcCanopy.Location(pcCanopy.Location(:,3)>(klist(k)-0.75) & pcCanopy.Location(:,3)<(klist(k)),:));
                convexh = polyshape(hlayer.Location(convhull(hlayer.Location(:,1),hlayer.Location(:,2)),1:2));
                minpointcount = 10*convexh.area;
                if ~isempty(hlayer.Location)
                    if hlayer.Count >= minpointcount
                       treemetr.h = hlayer.ZLimits(2);
                    end
                    if k == size(klist,2)
                        treemetr.h = pcCanopy.ZLimits(2);
                    elseif hlayer.Count < minpointcount && hlayer.ZLimits(2) < min(pcCanopy.ZLimits(2),18)
                        break
                    elseif hlayer.ZLimits(1) > 7 && convexh.area < 0.4
                        break
                    end
                end
                k = k + 1;
            end
        end
        if treemetr.h == 0
            treemetr.h = hEnd;
        end

        % Measure stem diameters
        h_start = h_butt;   % height for first diameter measurement
        h_end = treemetr.h; % tree height
        d_dif_underdbh = 0.2;  % relative diameter difference limit below 1.3 m height
        d_dif_overdbh = 0.12;  % relative diameter difference limit above 1.3 m height

        layer_heights = (h_start:h_interval:h_end);
        circle_list = zeros([length(layer_heights) 4]);
        for j = 1:length(layer_heights)
            layer_h = layer_heights(j);
            layer_low = layer_h - h_layer/2;
            layer_high = layer_h + h_layer/2;
            try
                layer_points = pcStem.Location(pcStem.Location(:,3) < layer_high & pcStem.Location(:,3) > layer_low,:);
                points = [layer_points(:,1) layer_points(:,2)];
            catch
                fprintf('No points at height %f\n',layer_low);
                continue
            end
            % Diameters are measured by fitting a circle to the points using CircleFitByTaubin.m
            if length(points) >= min_points
                par = CircleFitByTaubin(points);
                if par(3) <= 0.5
                    innercircle = [par(1:2) par(3)/2];
                    % figure, scatter(points(:,1),points(:,2),'k.'), hold on, axis equal;
                    % viscircles(cat(1,par(1:2),innercircle(1:2)),cat(1,par(3),innercircle(3))), hold off;
                    radialdistsq = (points(:,1) - innercircle(1)).^2 + (points(:,2) - innercircle(2)).^2;
                    isIn = radialdistsq <= innercircle(3)^2;
                    if sum(isIn) < 3
                        circle_list(j,1:3) = par;
                        circle_list(j,4) = layer_heights(j);
                    end
                end
            end
        end
        circle_list = circle_list(circle_list(:,3) > 0 & circle_list(:,3) < 0.4,:);
        dhlist = [circle_list(:,4) circle_list(:,3)*2];
        dhlist = sortrows(dhlist,1);
        circle_list = sortrows(circle_list,4);
        dhcell = cell([1 1]); dhcell{1} = dhlist;

        % OUTLIER REMOVAL
        % The stem is split into sections (2.5m), and diameters
        % deviating more than three median absolute deviations away
        % from the median diameter of the section are removed.
        splits = round(h_end/2.5,0);
        splitInterval = h_end/splits;
        for s = 1:splits
            splitinterval = splitInterval*s;
            inlierpisteet = [];
            tt = splits/s;
            for t = 1:tt
               pisteet = dhlist(dhlist(:,1)> (t*splitinterval-splitinterval) & dhlist(:,1)<= (t*splitinterval),:);
               TF = isoutlier(pisteet(:,2));
               pisteet(TF, :) = [];
               inlierpisteet = cat(1,inlierpisteet,pisteet);
            end
            dhlist = inlierpisteet;
        end

        % Comparison of diameter to mean of three previous (or
        % three closest at the bottom) diameters. A diameter is
        % considered as an outlier and removed if the relative
        % diameter difference exceeds the set limits.
        try
            underdbh = size(dhlist(dhlist(:,1)<= 1.3,:),1); % Number of diameters < 1.3 m
            remove_ids = [];            % ids of the removed diameters are collected here
            for j = 1:length(dhlist(:,1))
                    if j <= underdbh || underdbh < 3
                        if j == 1    % diameter is compared with 2 next diameters 
                            d = dhlist(j,2);
                            d_prev = mean([dhlist((j+1),2) dhlist((j+2),2)]);
                            d_diff = abs(d - d_prev)/d_prev;
                            if d_diff >= d_dif_underdbh
                                remove_ids = cat(2,remove_ids,j);
                                dhlist(j,2) = 1.02*d_prev;
                            end
                        elseif j == 2   % diameter is compared with neighbouring diameters 
                            d = dhlist(j,2);
                            d_prev = mean([dhlist((j-1),2) dhlist((j+1),2)]);
                            d_diff = abs(d - d_prev)/d_prev;
                            if d_diff >= d_dif_underdbh
                                remove_ids = cat(2,remove_ids,j);
                                dhlist(j,2) = d_prev;
                            end
                        elseif j == 3   % diameter is compared with 2 previous diameters
                            d = dhlist(j,2);
                            d_prev = mean([dhlist((j-1),2) dhlist((j-2),2)]);
                            d_diff = abs(d - d_prev)/d_prev;
                            if d_diff >= d_dif_underdbh
                                remove_ids = cat(2,remove_ids,j);
                                dhlist(j,2) = 0.98*d_prev;
                            end
                        else            % diameter is compared with 3 previous diameters
                            d = dhlist(j,2);
                            d_prev = mean([dhlist((j-1),2) dhlist((j-2),2) dhlist((j-3),2)]);
                            d_diff = abs(d - d_prev)/d_prev;
                            if d_diff >= d_dif_underdbh
                                remove_ids = cat(2,remove_ids,j);
                                dhlist(j,2) = 0.98*d_prev;
                            end
                        end
                    else                % diameter is compared with 3 previous diameters
                        d = dhlist(j,2);
                        d_prev = mean([dhlist((j-1),2) dhlist((j-2),2) dhlist((j-3),2)]);
                        d_diff = abs(d - d_prev);
                        if d_diff >= d_dif_overdbh
                            remove_ids = cat(2,remove_ids,j);
                            dhlist(j,2) = 0.98*d_prev;
                        end
                    end
            end
            if ~isempty(remove_ids)
                dhlist(remove_ids, :) = []; % outliers are removed
            end
        catch
            warning('Error in d-h-filtering for tree no. %i!',tree_id);
            return
        end


        % STEM CURVE
        % a smoothing cubic spline curve is fitted to the measured stem 
        % diameters to level the unevenness of the measurements and to interpolate the
        % missing diameters:

        heights = dhlist(:,1);
        diameters = dhlist(:,2);
        % dbhcheck = mean(dhlist(dhlist(:,1)> 0.8 & dhlist(:,1) < 2,2));
        % if (dbhcheck < 0.10 && abs(heights(end)-h_end) > 14) == true
        %     h_end = heights(end);
        % end
        % height and diameter at the top of the tree is set:
        heights = cat(1,heights,h_end);
        diameters = cat(1,diameters,0);
        sizediff = length(heights) - length(diameters);
        if sizediff < 0
                diameters = cat(1,diameters,zeros([1 abs(sizediff)]));
        elseif sizediff > 0
                diameters = diameters(1:end-sizediff);
        end
        h_smooth = 0:0.1:h_end;
        % Taper curve at 10 cm intervals from ground level to the tree top
        if (min(heights) < 2.5 && h_end > 3.5 && size(diameters,1) > 8)
            try
                d_smooth = csaps(heights,diameters,p,h_smooth);
                f = fopen([paramdir,'/d_smooth.txt'], 'w');fprintf(f, '%f\t', d_smooth);fclose(f);
            catch
                warning('Error in d-h-smoothing for tree no. %i!', tree_id);
                return
            end
        else
            warning('Not enough d-h-observations for tree no. %i!', tree_id);
            return
        end
        try
            f = fopen([paramdir,'/d_smooth.txt'], 'r');d_smooth = transpose(fscanf(f, '%f'));fclose(f);
            sizediff = length(h_smooth) - length(d_smooth);
            if sizediff < 0
                d_smooth = cat(1,d_smooth,zeros([1 abs(sizediff)]));
            elseif sizediff > 0
                d_smooth = d_smooth(1:end-sizediff);
            end
            taper = transpose(cat(1,h_smooth,d_smooth));
            taper(taper(:,2) < 0,:) = [];
            savedir = [rootResults,'/stemcurves'];
            if isfolder(savedir) == false, mkdir(savedir);end
            writetable(array2table(taper,'VariableNames',{'h','d'}),[savedir,sprintf('/tree%i_stemcurve.txt',tree_id)],'Delimiter','\t');
            proceed = true;
        catch
            proceed = false; trials = trials + 1;
            fprintf('trials: %i/6\n',trials)
        end
    end
    if proceed == false
        warning('Cannot construct taper curve for tree no. %i!',tree_id);
        return
    end
    
    % Construct cylinder model
    savedir = [rootResults,'/cylindermodels'];
    if isfolder(savedir) == false, mkdir(savedir);end
    h_max = max(taper(:,1));
    z_smooth = 0:0.1:h_max;
    pcyl = 0.8;
    cyllist = circle_list;
    cyllist(~ismember(cyllist(:,4),heights),:) = [];
    xxx = cyllist(:,1); yyy = cyllist(:,2); zzz = cyllist(:,4);
    x_smooth = csaps(zzz,xxx,pcyl,z_smooth);
    y_smooth = csaps(zzz,yyy,pcyl,z_smooth);
    diam = taper(taper(:,1)<=(h_max+0.01),2);
    sizediff = length(diam) - length(x_smooth);
    if sizediff == 0
        cylindermodel = cat(2,x_smooth',y_smooth',z_smooth',diam,0.1*ones([length(z_smooth),1]));
        writetable(array2table(cylindermodel,'VariableNames',{'x','y','z','d','h'}),[savedir,sprintf('/tree%i_cylinders.txt',tree_id)],'Delimiter','\t');
    elseif sizediff < 0
        diam = cat(1,diam,zeros([abs(sizediff),1]));
        cylindermodel = cat(2,x_smooth',y_smooth',z_smooth',diam,0.1*ones([length(z_smooth),1]));
        writetable(array2table(cylindermodel,'VariableNames',{'x','y','z','d','h'}),[savedir,sprintf('/tree%i_cylinders.txt',tree_id)],'Delimiter','\t');
    elseif sizediff > 0
        diam = diam(1:end-sizediff);
        cylindermodel = cat(2,x_smooth',y_smooth',z_smooth',diam,0.1*ones([length(z_smooth),1]));
        writetable(array2table(cylindermodel,'VariableNames',{'x','y','z','d','h'}),[savedir,sprintf('/tree%i_cylinders.txt',tree_id)],'Delimiter','\t');
    end
    % Get tree location
    treemetr.X = cylindermodel(cylindermodel(:,3) == 1.3,1);
    treemetr.Y = cylindermodel(cylindermodel(:,3) == 1.3,2);

    % Calculate volume and store tree metrics
    vol = 0; vlog = 0; vpulp = 0;
    for j = 2:length(taper(:,1))
        cyl_h = taper(j,1)-taper((j-1),1);
        cyl_r = mean([taper(j,2) taper((j-1),2)])/2;
        cyl_vol = pi*cyl_r^2*cyl_h;
        vol = vol + cyl_vol;
        if cyl_r >= 0.08, vlog = vlog + cyl_vol; end
        if (cyl_r <= 0.08 && cyl_r >= 0.03), vpulp = vpulp + cyl_vol; end
    end
    
    dbh = taper(round(taper(:,1),1) == 1.3,2)*100;
    treemetr.dbh = dbh;
    treemetr.d6 = taper(round(taper(:,1),1) == 6,2)*100;
    treemetr.d05 = taper(round(taper(:,1),1) == round(max(taper(:,1))*0.5,1),2)*100;
    treemetr.ba = pi/4*(dbh/100)^2;
    treemetr.h = max(taper(:,1));
    treemetr.vol = vol;
    treemetr.vlog = vlog;
    treemetr.vpulp = vpulp;
    
    pcStem.Color = uint8([170 115 30] .* ones(size(pcStem.Location)));
    pcCanopy.Color = uint8([110 190 100] .* ones(size(pcCanopy.Location)));
    
    savedir = [rootResults,'/figures/pointclassification'];
    if isfolder(savedir) == false; mkdir(savedir);end
    % Illustrate point cloud classification
    myfunctions.pcshow_parfor(pcStem,pcCanopy)
    filename = sprintf([savedir,sprintf('/pc_treeID_%i.png',tree_id)]);
    saveas(gcf,filename);
    close all
    
    % Save stem curve figure
    savedir = [rootResults,'/figures/stemcurves'];
    if isfolder(savedir) == false; mkdir(savedir); end
    ff = figure('visible','off');
    plot(dhcell{1}(:,1), dhcell{1}(:,2),'r.'), title(sprintf('tree #%i',tree_id)),hold on;
    plot(heights, diameters,'k.'), hold on;
    plot(0:0.1:h_end,d_smooth,'b-'),ylim([0 (max(d_smooth)+0.1)]), hold off
    legend('outliers','filtered','smoothed stem curve');
    filename = sprintf([savedir,sprintf('/stem_curve_tree%i.png',tree_id)]);
    saveas(ff,filename)
end