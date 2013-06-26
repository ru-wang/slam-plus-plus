/*
								+----------------------------------+
								|                                  |
								|  ***  3D solvers instances  ***  |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2013  |
								|                                  |
								|         Solve3DImpl.cpp          |
								|                                  |
								+----------------------------------+
*/

/**
 *	@file src/slam/Solve3DImpl.cpp
 *	@brief contains instantiation of SE(3) solver templates
 *	@author -tHE SWINe-
 *	@date 2013-06-14
 */

#include "slam/Main.h"
#include "slam/ConfigSolvers.h" // only included in files that actually need the solvers (slow to compile)
#include "slam/SE3_Types.h"

/**
 *	@def __SE3_ENABLED
 *	@brief if defined, the solver specializations for SE(3) are compiled
 */
#define __SE3_ENABLED

int n_Run_SE3_Solver(TCommandLineArgs t_args) // throw(std::runtime_error, std::bad_alloc)
{
#ifdef __SE3_ENABLED
	typedef MakeTypelist_Safe((CVertexPose3D)) TVertexTypelist_SE3;
	typedef MakeTypelist_Safe((CEdgePose3D)) TEdgeTypelist_SE3;
	// define types of vertices, edges

	typedef CFlatSystem<CVertexPose3D/*CSEBaseVertex*/, TVertexTypelist_SE3,
		CEdgePose3D/*CSEBaseEdge*/, TEdgeTypelist_SE3> CSystemType;
	// make a system permitting SE(3) vertex and edge types

	typedef CSolverCaller<CSystemType, CSE3OnlyPoseEdgeTraits, CParseLoop> CSpecializedSolverCaller;
	// specify how the nonlinear solver should be called

	return CTypelistForEach<CCompiledSolverList,
		CSpecializedSolverCaller>::Run(CSpecializedSolverCaller(t_args)).n_Result();
#else // __BA_ENABLED
	fprintf(stderr, "error: SE(3) is disabled\n");
	return -1;
#endif // __SE3_ENABLED
}
