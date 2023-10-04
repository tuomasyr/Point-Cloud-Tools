# Point-Cloud-Tools

Dr. Tuomas Yrttimaa, School of Forest Sciences, University of Eastern Finland // tuomas.yrttimaa@uef.fi

## Disclaimer
This is a readme document for the users of automatic point cloud processing tools to characterize trees initially developed and described in detail in Yrttimaa et al. ([2019](https://doi.org/10.3390/rs11121423), [2020](https://doi.org/10.1016/j.isprsjprs.2020.08.017)). The tools have been implemented in [MATLAB](https://se.mathworks.com/products/matlab.html) but the workflow requires the user to also have a licensed [LAStools software](https://rapidlasso.com/) installed on the computer. The tools can be used free of charge under the [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/). Appropriately citing the original research articles as well as this document is preferred when used within academia. No guarantees or technical support related to the use and applicability of the tools are provided by default. For scientific or commercial cooperation, please contact the author.

To cite Yrttimaa T. (2021). Automatic Point Cloud Processing Tools to Characterize Trees (Point-Cloud-Tools). Zenodo. [https://doi.org/10.5281/zenodo.5779288](https://doi.org/10.5281/zenodo.5779288).

## Contents
The point cloud processing tools consist of a main script (mainscript.m) and a group of functions (see Table 1) executed subsequently to:
1) import point cloud data (pc_import.m),
2) detect tree crown segments (pc_detect_tree_crowns_v2.m),
3) extract point clouds for individual trees or groups of trees (pc_tree_extraction.m),
4) classify individual tree point clouds into stem and non-stem points (pc_classify_v2),
5) compute tree metrics for each tree (pc_treeMetrics.m, compute_treemetrics.m, compute_crownmetrics.m), and
6) detect and segment branches of individual trees (pc_detect_branches.m, run within previous step)

As an input, the tools require a single height-normalized point cloud representing the tree community of interest. By default, the functions read .las/laz-files that can be obtained from practically any point cloud file format using the open source tools of LAStools software. The height normalization procedure can be carried out using the LAStools software and following the workflow presented in e.g. Ritter et al. (2017). 

As an output, the point cloud processing tools return for each tree, identified by a treeID number: 
1) a classified point cloud with stem points and crown points labelled separately, 
2) a stem taper curve (i.e. diameters along the stem from base to top of the tree)
3) figures to illustrate point cloud classification and stem taper curve estimation
3) cylinder models representing the stem dimensions and volume,
4) stem metrics, or attributes characterizing the stem properties such as the XY-location, diameter (dbh) and basal area (ba) at the breast height, tree height (h), stem volume (vol), diameter at 6 m height d6), diameter at 50% of tree height (d05), logwood volume (vlog, i.e. volume of the stem section with diameter > 15 cm) and pulpwood volume (vpulp i.e. volume of the stem section with diameter < 15 cm but > 8 cm).
5) crown metrics, or attributes characterizing the crown properties such as the XY-location, crown volume based on 3D convex hull (vol_convhull) and voxelization (vol_voxel), projection area (area2d), surface area (area3d), perimeter (perimeter) and crown base height (crownbaseheight) in m, m^2 or m^3.
6) point cloud clusters representing individual branches, information of their characteristics (e.g. location of branch origin on stem surface, branching angle, length).

The point cloud processing tools have been designed to process one multi-scan terrestrial laser scanning point cloud covering a sample plot of 32 m x 32 m in size, including a variety of different-sized trees. Performance of the point cloud-based methods have been investigated in southern boreal forests, and extensive investigations of its feasibility in characterizing trees is summarized in [Yrttimaa (2021)](https://doi.org/10.14214/df.314). When applied outside the southern boreal forest zone, the user may need to adjust some of the computing parameters. 

## Description of the contents

- **mainscript.m**
Define parameter values, directories, and filenames, run the processing functions.

- **pc_import.m**
Import height-normalized terrestrial laser scanning point cloud data as .las/laz-format and extract point cloud extent.

- **pc_detect_tree_crowns_v2.m**
Detect tree crown segments to partition the point cloud into smaller parts in a meaningful way. Outputs crown segment polygons.

- **pc_tree_extraction.m**
Delineate the sample plot point cloud based on the detected crown segments using the point-in-polygon approach. This uses the [lasclip](https://rapidlasso.com/lastools/lasclip/)-tool from the LAStools software and outputs .las-files of individual trees.

- **pc_classify_v2.m**
This function executes the point cloud classification for each individual tree point cloud. It is an iterative procedure beginning from the base of the tree and proceeding towards the tree top, aiming at separating points representing the tree stem (i.e., stem points) from points not representing the tree stem (i.e., non-stem points). Smooth, cylindrical and vertical surfaces with vertical continuity are the point cloud characteristics that are to distinct a tree stem from branches and foliage. For more details, see [Yrttimaa et al. (2020)](https://doi.org/10.1016/j.isprsjprs.2020.08.017). Outputs stem points and crown points as .ply-files for each tree. The LAStools software can be used to merge tree-wise extracted point clouds together.
 
 - **pc_treeMetrics.m**
Loops through all the detected trees and uses the functions **compute_treemetrics.m**, **compute_crownmetrics.m** and **pc_detect_branches.m** to extract stem and crown metrics (described above) as well as to detect and segment individual branches (implemented and validated for Scots pine trees only!) 

- **compute_treemetrics.m**
Tree metrics for each extracted tree is computed using this function. It uses the function [CircleFitByTaubin.m](https://www.mathworks.com/matlabcentral/fileexchange/22678-circle-fit-taubin-method) (Chernov 2021) to fit circles into horizontal point cloud slices to measure stem cross-section diameters. Taper curve is stored for each tree to quantify its characteristics. For more details, see [Yrttimaa et al. (2019)](https://doi.org/10.3390/rs11121423).

- **compute_crownmetrics.m**
Computes metrics characterizing the crown structure of individual trees. THe procedure is described in [Yrttimaa et al. (2022)](https://doi.org/10.1016/j.foreco.2022.120303).

- **pc_detect_branches.m**
Identifies and segments individual branches and returns information of their characteristics.

- **lastools.m**
This function enables the user to run LAStools commands on MATLAB. User’s own LAStools bin directory must be defined in order to use this function.

- **myfunctions.m**
Contains some functions to be used in the data processing workflow.

## References
Chernov N. 2021. Circle Fit (Taubin method) (https://www.mathworks.com/matlabcentral/fileexchange/22678-circle-fit-taubin-method), MATLAB Central File Exchange. Retrieved December 14, 2021.

Ritter T, Schwarz M, Tockner A, Leisch F, Nothdurft A. 2017. Automatic mapping of forest stands based on three-dimensional point clouds derived from terrestrial laser-scanning. Forests, 8(8), 265. https://doi.org/10.3390/f8080265

Yrttimaa T, Saarinen N, Kankare V, Liang X, Hyyppä J, Holopainen M, Vastaranta M. 2019. Investigating the feasibility of multi-scan terrestrial laser scanning to characterize tree communities in southern boreal forests. Remote Sensing, 11(12), 1423. https://doi.org/10.3390/rs11121423

Yrttimaa T, Saarinen N, Kankare V, Hynynen J, Huuskonen S, Holopainen M, Hyyppä J., Vastaranta, M. 2020. Performance of terrestrial laser scanning to characterize managed Scots pine (Pinus sylvestris L.) stands is dependent on forest structural variation. ISPRS Journal of Photogrammetry and Remote Sensing, 168, 277-287. https://doi.org/10.1016/j.isprsjprs.2020.08.017

Yrttimaa T. 2021. Characterizing tree communities in space and time using point clouds. Dissertationes Forestales 314. 52 p. https://doi.org/10.14214/df.314

Yrttimaa T, Luoma V, Saarinen N, Kankare V, Junttila S, Holopainen M, Hyyppä J, Vastaranta M, 2022. Exploring tree growth allometry using two-date terrestrial laser scanning. Forest Ecology and Management, 518, p.120303. https://doi.org/10.1016/j.foreco.2022.120303
