% Import point cloud data into MATLAB.
%
% Input: 
%         - rootTLSData      % directory for TLS data [1x1 character]
%         - filenameTLS      % name of TLS data file [1x1 character]
% 
% Output: - PCtls            % TLS point cloud [pointCloud object]
%         - plotExtent       % min and max XY [1x4 table]
%
% (c) Tuomas Yrttimaa / Science4Trees @ UEF School of Forest Sciences 2019
% ---------------------------------------------------------------------------
function [PCtls,plotExtent] = pc_import(rootTLSData,filenameTLS)
    starttime = datetime;
    hh = waitbar(0,''); hh.Name = 'Point Cloud Tools: Import point cloud data';
  
    waitbar(5/100,hh,'Loading TLS data, please wait...');
   
    % Read the las-files
    lasReader = lasFileReader(strcat(rootTLSData,'/',filenameTLS));
    PCtls = readPointCloud(lasReader);
    %PCtls = pointCloud(PCtls.Location(PCtls.Location(:,3)< 40 & PCtls.Location(:,3)>=0,:));
    clearvars c1;

    % min-max values for xy-coordinates
    min_x = PCtls.XLimits(1); max_x = PCtls.XLimits(2);
    min_y = PCtls.YLimits(1); max_y = PCtls.YLimits(2);
    plotExtent = table(min_x,max_x,min_y,max_y);
   
    close(hh);
    endtime = datetime;
    cprintf('comment',sprintf('   Point cloud data loaded in %s\n',endtime-starttime));
end