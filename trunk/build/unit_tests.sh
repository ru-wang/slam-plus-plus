#!/bin/sh

#rm -f src/slam_app/CMakeFiles/slam_plus_plus.dir/Solve3DImpl.o src/slam_app/CMakeFiles/slam_plus_plus.dir/Solve3DPoseOnlyImpl.o src/slam_app/CMakeFiles/slam_plus_plus.dir/SolveBAImpl.o src/slam_app/CMakeFiles/slam_plus_plus.dir/SolveBAIntrinsicsImpl.o src/slam_app/CMakeFiles/slam_plus_plus.dir/SolveBAStereoImpl.o src/slam_app/CMakeFiles/slam_plus_plus.dir/SolveSpheronImpl.o
#make slam_plus_plus
# need to rebuild?

data_path="../data"
#data_path="$HOME/my-projects/SLAM_plus_plus/data"

datasets[0]="$data_path/intel.txt"
md5sums[0]="9fc7f05950e0a5b57ed2663a2c0f3db4"
command[0]=" -po "
chi2s[0]=17 # 17.06
iters[0]=3

datasets[1]="$data_path/manhattanOlson3500.txt"
md5sums[1]="f9e150f35a47bbaf8433efd25eb32c4b"
command[1]=" -po "
chi2s[1]=146 # 146.08
iters[1]=5

datasets[2]="$data_path/10K.graph"
md5sums[2]="df8fb5bb0f14ed0e0d51090e0a796cca"
command[2]=" -po "
chi2s[2]=303 # 303.17
iters[2]=5

datasets[3]="$data_path/10KHOG-MAN_g2o.txt"
md5sums[3]="50c4b37aef5347f13ddb059a1d9c46df"
command[3]=" -po "
iters[3]=5
chi2s[3]=171545 # 171545.45

datasets[4]="$data_path/victoria_park_original.graph"
md5sums[4]="b64e9f2eb9f322eea96a1354991e1de0"
command[4]=" -nsp 1 -s "
iters[4]=3424
chi2s[4]=140 # 140.09

datasets[5]="$data_path/sphere2500.txt"
md5sums[5]="4a78a67f9ce065bde01e55590a3d5fed"

command[5]=" -po "
chi2s[5]=728 # 727.72
iters[5]=5

datasets[6]="$data_path/parking-garage.txt"
md5sums[6]="1e8889f61e4f4e09451ebfdd921eb565"
command[6]=" -po "
chi2s[6]=1 # 1.46
iters[6]=5

datasets[7]="$data_path/w100K_sort.txt"
md5sums[7]="c595eafbaeabe72185aa9ddfa42f471d"
command[7]=" -po "
chi2s[7]=8685 # 8685.07
iters[7]=5

datasets[8]="$data_path/venice871.g2o"
md5sums[8]="54f3f9dd06e781f0c69abd6911967498"
command[8]=" -us "
chi2s[8]=234013899 # 234013899.18
iters[8]=5

processed=0
problems=0
index=0
for ds in "${datasets[@]}"
do
	if [ ! -f $ds ]; then
		echo "warning: $ds not found"
		index=$(( $index + 1 ))
		continue
	fi

	processed=$(( $processed + 1 ))

	ds_sum="`md5sum $ds | cut -d ' ' -f 1`"
	stored_sum="${md5sums[$index]}"
	if [ $ds_sum != $stored_sum ]; then
		echo "$ds - bad checksum - $ds_sum != $stored_sum"
	fi
	# make sure the dataset is the correct dataset

	rm -f slampp.log

	`../bin/slam_plus_plus -i $ds ${command[$index]} > slampp.log 2>&1`
	#cat slampp.log
	# run slam++

	num_iters=`cat slampp.log | egrep "solver took [0-9]+ iterations" | cut -d ' ' -f 3 | tail -n 1`
	chi2=`cat slampp.log | grep "denormalized chi2 error:" | cut -d ' ' -f 4 | tail -n 1`
	# get the number of iterations and chi2

	chi2_round=`echo $chi2 | gawk '{ print int($0 + .5); }'`

	#echo "$ds: $num_iters iterations, chi2: $chi2 ($chi2_round)"
	# verbose - dont want

	if [[ -z $num_iters || -z $chi2 ]]; then
		echo "$ds - SLAM++ did not finish (`cat slampp.log | grep error | wc -l` error(s))"
		problems=$(( $problems + 1 ))
	else
		if [ $num_iters -ne ${iters[$index]} ]; then
			echo "$ds - bad number of iterations - $num_iters != ${iters[$index]}"
			problems=$(( $problems + 1 ))
		fi
		if [ $chi2_round -ne ${chi2s[$index]} ]; then
			echo "$ds - bad chi2 - $chi2 ($chi2_round) != ${chi2s[$index]}"
			problems=$(( $problems + 1 ))
		fi
	fi
	# make sure that

	#echo "chi2s[$index]=$chi2_round # $chi2"
	#echo "iters[$index]=$num_iters"
	# print reference values

	rm -f slampp.log

	index=$(( $index + 1 ))
done

if [ $problems -eq 0 ]; then
	echo "no problems occured, all $processed tests finished successfully"
else
	echo "finished; processed $processed datasets, there were $problems problem(s)"
fi
