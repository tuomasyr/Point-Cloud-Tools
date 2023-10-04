% Function to run LAStools commands on MATLAB.
% 
% LAStools command to be executed is entered as a character vector
%
% The use of this function requires licensed LAStools (see rapidlasso.com)
% -------------------------------------------------------------------------
% example: lastools('lasinfo -i pointcloud.las')
% -------------------------------------------------------------------------

function lastools(las_string)
    ld = 'C:/LAStools/bin/';
    command = [ld,las_string];
    system(command);
end