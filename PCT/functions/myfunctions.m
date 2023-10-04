
% Some functions for point cloud processing

classdef myfunctions
    methods(Static)
        function pcshow_parfor(stem,crown)
           points1 = stem.Location;
           pointcol1 = double(stem.Color(1,:))./255;
           points2 = crown.Location;
           pointcol2 = double(crown.Color(1,:))./255;
           figure, set(gcf,'Position',[50 50 450 900]);
           scatter3(points1(:,1),points1(:,2),points1(:,3),1,pointcol1), hold on;
           scatter3(points2(:,1),points2(:,2),points2(:,3),1,pointcol2); view(90,0);axis equal;
           xlim(crown.XLimits), ylim(crown.YLimits), zlim([0 crown.ZLimits(2)]);
           xlabel('X (m)'), ylabel('Y (m)'), zlabel('Height (m)');
           zticks(0:2:floor(crown.ZLimits(2)));
           set(gca,'color','w');set(gcf,'color','w');set(gca, 'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'ZColor', [0.15 0.15 0.15]);
        end

        function treecylinder = import_treecylinder(filename)
            dataLines = [1, Inf];
            opts = delimitedTextImportOptions("NumVariables", 7);
            opts.DataLines = dataLines;
            opts.Delimiter = "\t";
            opts.VariableNames = ["x_bottom", "y_bottom", "h_bottom", "x_top", "y_top", "h_top", "diam"];
            opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            treecylinder = readtable(filename, opts);
        end

        function pc = importclassifiedpcfile(rootData,tree_id,pcType)
            if strcmp(pcType,'stem') == true
                filename = sprintf([rootData,'/pointclouds/stem/tree%i_stempoints.txt'],tree_id);
            elseif strcmp(pcType,'canopy') == true
                filename = sprintf([rootData,'/pointclouds/canopy/tree%i_canopypoints.txt'],tree_id);
            end
            dataLines = [1, Inf];
            opts = delimitedTextImportOptions("NumVariables", 4);
            opts.DataLines = dataLines;
            opts.Delimiter = "\t";
            opts.VariableNames = ["X", "Y", "Z", "ID"];
            opts.VariableTypes = ["double", "double", "double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            pc = readtable(filename, opts);
            pc = pointCloud(table2array(pc(:,1:3)));

        end
    end
end
