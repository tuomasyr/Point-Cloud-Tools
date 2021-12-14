% Extract TLS and UAV point clouds representing individual trees detected from the
% TLS point cloud. Points inside the tree segment are considered to form a 
% point cloud representing each tree.
%
% Input:  - TLSpoints     % TLS data [pointCloud object]
%         - ALSpoints     % ALS data [pointCloud object]
%         - treeLocations % XY-locations for each detected tree [Mx2 double array]
%         - treePgons     % tree crown segment polygons [Mx1 cell array]
%
% Output: - TLSPCs        % Tree-wise extracted TLS point cloud objects [Mx1 cell array]
%         - ALSPCs        % Tree-wise extracted ALS point cloud objects [Mx1 cell array]
%
% (c) Tuomas Yrttimaa / Science4Trees @ UEF School of Forest Sciences 2019
% ---------------------------------------------------------------------------
function pc_tree_extraction(rootTLSData,filenameTLS)
    hh = waitbar(0,''); hh.Name = 'Point Cloud Tools: Extract Tree Segments';
    waitbar(15/100,hh,'Extracting tree segments, please wait...');
    starttime = datetime;
    savedir = [rootTLSData,'/tmp']; if isfolder(savedir) == false, mkdir(savedir); end
    list = struct2table(dir(fullfile('./tmp/treepoly', '**', '*.txt')));
    filenames = list.name;
    
    parfor_progress(size(filenames,1));
    parfor i = 1:size(filenames,1)
        parfor_progress;
        try
            poly = filenames{i};
            lastools(['lasclip -i "',rootTLSData,'/',filenameTLS,'" -poly ./tmp/treepoly/',poly,' -odir "',savedir,'" -o ',sprintf('treeseg_%i.las',i)]);
        catch
            continue
        end
    end
    parfor_progress(0);
    endtime = datetime; runtime = endtime-starttime;
    cprintf('comment',   sprintf('   Trees delineated in: %s\n',runtime));
    close(hh);
end

