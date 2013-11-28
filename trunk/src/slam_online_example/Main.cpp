/*
								+-----------------------------------+
								|                                   |
								|   ***  Online SLAM Example  ***   |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2013  |
								|                                   |
								|             Main.cpp              |
								|                                   |
								+-----------------------------------+
*/

/**
 *	@file src/slam_online_example/Main.cpp
 *	@brief contains the main() function of the online example SLAM program
 *	@author -tHE SWINe-
 *	@date 2013-11-27
 */

#include <stdio.h> // printf
#include "slam/LinearSolver_UberBlock.h"
#include "slam/LinearSolver_CholMod.h" // linear solvers (only one is required)
#include "slam/ConfigSolvers.h" // nonlinear graph solvers
#include "slam/SE2_Types.h" // SE(2) types

#if defined(_WIN32) || defined(_WIN64)
#define NOMINMAX
#include <windows.h> // to show the result image, otherwise not required
#endif // _WIN32) || _WIN64

/**
 *	@brief main
 *
 *	@param[in] n_arg_num is number of commandline arguments
 *	@param[in] p_arg_list is the list of commandline arguments
 *
 *	@return Returns 0 on success, -1 on failure.
 */
int main(int UNUSED(n_arg_num), const char **UNUSED(p_arg_list))
{
	typedef MakeTypelist(CVertexPose2D) TVertexTypelist;
	typedef MakeTypelist(CEdgePose2D) TEdgeTypelist;

	typedef CFlatSystem<CVertexPose2D, TVertexTypelist, CEdgePose2D, TEdgeTypelist> CSystemType;

	typedef CLinearSolver_UberBlock<CSystemType::_TyHessianMatrixBlockList> CLinearSolverType;

	CSystemType system;
	CNonlinearSolver_FastL<CSystemType, CLinearSolverType> solver(system, 0, 1); // solve each 1

	Eigen::Matrix3d information;
	information <<
		45,  0,  0,
		 0, 45,  0,
		 0,  0, 45;
	// prepare the information matrix (all edges have the same)

	system.r_Get_Vertex<CVertexPose2D>(0, Eigen::Vector3d(0, 0, M_PI / 4));
	// initialize the first vertex (the robot starts rotated 45 debrees)

	solver.Incremental_Step(system.r_Add_Edge(CEdgePose2D(0, 1,
		Eigen::Vector3d(1.02765, -0.00745597, 0.00483283), information, system))); // replace (1.02765, -0.00745597, 0.00483283) with readings from the ROS stack
	// get the solution every step

	system.r_Vertex_Pool(); // returns container with all the vertices which are optimized
	system.Dump("vertices.txt");  // puts the same vertices in a text file
	system.Plot2D("afterFirstStep.tga", plot_quality::plot_Printing); // plot in print quality
	// get a snapshot of the system

	solver.Delay_Optimization();
	// when adding a bunch of measurements, maybe we don't care about incremental solutions in between

	solver.Incremental_Step(system.r_Add_Edge(CEdgePose2D(1, 2,
		Eigen::Vector3d(-0.0120155, 1.00436, 1.56679), information, system)));
	solver.Incremental_Step(system.r_Add_Edge(CEdgePose2D(2, 3,
		Eigen::Vector3d(0.0440195, 0.988477, -1.56353), information, system)));
	solver.Incremental_Step(system.r_Add_Edge(CEdgePose2D(3, 0,
		Eigen::Vector3d(0.0139844, -1.02344, -0.00780158), information, system)));
	// put edges in the system, the vertices are created and initialized

	solver.Enable_Optimization();
	// re-enable optimization, and optimize right away if it needs to

	system.Plot2D("result.tga", plot_quality::plot_Printing); // plot in print quality

#if defined(_WIN32) || defined(_WIN64)
	ShellExecute(0, "open", "result.tga", 0, 0, SW_SHOW);
#endif // _WIN32) || _WIN64
	// on windows, we can open the file with the results in a window

	return 0;
}

/*
 *	end-of-file
 */
