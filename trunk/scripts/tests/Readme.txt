This is for internal unit testing and possibly also benchmarking.

The make_unittest.sh is a simple script that records a reference
result. This is only supposed to be run once, using a highly stable
version of SLAM ++. It prints parts of the second script, which
verifies these results using the current version of SLAM ++ and
prints any differences.

unit_tests.sh [-v|--verbose|-q|--quiet] [-rt|--robust-timing]
	[-i|--data-path <path-to-SLAM_plus_plus/data>]

runs the unit tests with the specified options. The ground truth
part of the script is not supposed to be modified under normal
conditions.

In case the data path is not specified, the script tries to read
the SLAMPP_DATA environment variable and if that is not set, falls
back to the default "../data" (which works, assuming this is being
run from build).
