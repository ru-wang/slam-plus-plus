=== Building ===

This library was tested on Windows (using MSVC 2008 and above),
Mac, Linux and BSD. There is an optimizer issue in g++ 4.2 or earlier
which causes the code to crash at runtime.

You can download thestable sources in a zip file from SourceForge,
or checkout the newest sources from the svn:

$ mkdir slam
$ cd slam
$ svn co svn+ssh://svn.code.sf.net/p/slam-plus-plus/code/

There is a CMakeFile. To be able to change the code and commit back
to the svn, do an out-of-source build, like this:

$ cd build
$ cmake -i ..
$ <press enter 100x, real fast>
$ make
$ ../bin/slam_plus_plus -i ../data/manhattanOlson3500.txt --pose-only

You can also use fast parallel build, like this:

$ make -j 8

And that should take care of the build. Don't use CMake for Visual Studio,
there is a link conflict in eigen, there are some files in different
directories that generate the same .obj files, overwriting each other.
Either download the configured workspace from SorceForge, or if you must,
generate, open the workspace in Visual Studio, open Solution Explorer and
in project eigen3, right-click on complex_double.cpp, go to "C++" and
"Output Files", and change "Object File Name" from "$(IntDir)\" to
"$(IntDir)\lapack\" (the second slash is important). Now click on the
second complex_double.cpp and repeat the procedure, only put "$(IntDir)\blas\".
This needs to be repeated also for complex_single.cpp, double.cpp and
single.cpp.

Note that with CMake, there is no way to do multiple target x86 / x64
builds in Visual Studio.

Also, it's good to tweak optimizations in Visual Studio, you can turn
"Full Optimization" on, enable link-time code generation, and importantly
enable OpenMP support.

=== Data ===

For the default optimization framework, it is mostly safe to use iSAM
(http://people.csail.mit.edu/kaess/isam/) datafiles.
