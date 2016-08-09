=== Building ===

This library was tested on Windows (using MSVC 2008 and above),
Mac, Linux and BSD. There is an optimizer issue in g++ 4.2 or earlier
which causes the code to crash at runtime.

You can download thestable sources in a zip file from SourceForge,
or checkout the newest sources from the svn:

$ mkdir slam
$ cd slam
$ svn checkout svn://svn.code.sf.net/p/slam-plus-plus/code/trunk .

There is a CMakeFile. To be able to change the code and commit back
to the svn, do an out-of-source build, like this:

$ cd build
$ cmake ..

This will configure the project without any configuration. To change
configuration, run cmake -i .. instead of the last line above (or at
a later time, should a change in the configuration be needed). One
interesting option is to specify the default linear solver. Supernodal
CHOLMOD or block Cholesky are the fastest, CSparse is slightly slower
and simplical CHOLMOD is the slowest. Another option is to enable GPU
acceleration support, which currently applies to the Schur complement
solver only (specify -us when running SLAM ++; requires CUDA and CULA
toolkits).

On Mac, you might want to configure your C and C++ compilers to be the
GNU ones (e.g. from the MacPorts project) rather than Clang which does
not support OpenMP nowadays and your code will be somewhat slower. You
can do that by using:

$ cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.7 \
$	-D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.7 ..

Where you might want to change the version to the latest one that you
have installed. But it will build and run correctly even without that.

$ make 

Will finish the building process. You should now be able to run by typing:

$ ../bin/SLAM_plus_plus --help 

You can also use fast parallel build, like this:

$ make -j

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

Data can be downloaded from SourceForge, at:
http://sourceforge.net/projects/slam-plus-plus/files/data/
