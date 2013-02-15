This file briefly describes format and sources of 2D dataset files for SLAM++.

There are several conventions for the 2D graph SLAM datasets, each with its own
shortcomings, some of them requiring additional processing for incremental
scenarios. In order to run our code with several already existing datasets we
use a unified graph SLAM dataset format.

---------------------------------------------------
- SLAM++ format specifications
---------------------------------------------------
vertices:
    VERTEX_SE2/VERTEX2 <node_id> <X> <Y> <Theta>

edges:
    EDGE_SE2/EDGE2/ODOMETRY <node_id_from> <node_id_to> <X> <Y> <Theta> <XX> <XY> <XT> <YY> <YT> <TT>
	where node_id_from < node_id_to

Variable initialization entries are supproted and marked by tag VERTEX_SE2 or VERTEX2.
Edges are marked by tags EDGE_SE2, EDGE2 or ODOMETRY.
Angles are in radians.
The information matrix associated to each edge is NOT in squared form. The
information is in upper-triangular form (xx xy xy yy yt yy). The matrix is
stored row by row from top to bottom, with left to right column order. Example:

|0 1 2|
|  3 4|
|    5|

---------------------------------------------------
- Datasets
---------------------------------------------------
The following text contains information about sources of dataset files
available to download.


intel [4]
---------------------------------------------------
download: http://sourceforge.net/projects/slam-plus-plus/files/data/intel.txt/download

manhattanOlson3500 [2]
---------------------------------------------------
modification: the information values have been squared
download: http://sourceforge.net/projects/slam-plus-plus/files/data/manhattanOlson3500.txt/download

10kHog-man  [1]
---------------------------------------------------
modification: order of information values has been changed
download: http://sourceforge.net/projects/slam-plus-plus/files/data/10kHog-man.txt/download

10k [1]
---------------------------------------------------
modification: order of information values has been changed
download: http://sourceforge.net/projects/slam-plus-plus/files/data/10k.txt/download

100k [1]
---------------------------------------------------
modification: order of information values has been changed
download: http://sourceforge.net/projects/slam-plus-plus/files/data/100k.txt/download

city10k [3]
---------------------------------------------------
modification: the information values have been squared
download: http://sourceforge.net/projects/slam-plus-plus/files/data/city10k.txt/download

cityTrees10k [3]
---------------------------------------------------
modification: the information values have been squared
download: http://sourceforge.net/projects/slam-plus-plus/files/data/cityTrees10k.txt/download

Victoria park [4]
---------------------------------------------------
modification: the information values have been squared
download: http://sourceforge.net/projects/slam-plus-plus/files/data/victoria-park.txt/download

Killian court [5]
---------------------------------------------------
modification: the information values have been squared and the order changed
download: http://sourceforge.net/projects/slam-plus-plus/files/data/killian-court.txt/download


[1] G. Grisetti, C. Stachniss, S. Grzonka, and W. Burgard, "A tree parameterization for efficiently computing maximum likelihood maps using gradient descent," in Robotics: Science and Systems (RSS), June 2007.
[2] E. Olson, "Robust and efficient robot mapping," Ph.D. dissertation, Massachusetts Institute of Technology, 2008.
[3] M. Kaess, A. Ranganathan, and F. Dellaert, "iSAM: Fast incremental smoothing and mapping with efficient data association," in IEEE Intl. Conf. on Robotics and Automation (ICRA), Rome, Italy, April 2007, pp. 1670-1677.
[4] A. Howard and N. Roy, "The robotics data set repository (Radish)," 2003. [Online]. Available: http://radish.sourceforge.net/
[5] M. Bosse, P. Newman, J. Leonard, and S. Teller, "Simultaneous localization and map building in large-scale cyclic environments using the Atlas framework," Intl. J. of Robotics Research, vol. 23, no. 12, pp. 1113-1139, Dec 2004.
