/*
								+----------------------------------+
								|                                  |
								|  ***  2D solvers instances  ***  |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2013  |
								|                                  |
								|         Solve2DImpl.cpp          |
								|                                  |
								+----------------------------------+
*/

/**
 *	@file src/slam/Solve2DImpl.cpp
 *	@brief contains instantiation of SE(2) solver templates
 *	@author -tHE SWINe-
 *	@date 2013-06-14
 */

#include "slam/Main.h"
#include "slam/ConfigSolvers.h" // only included in files that actually need the solvers (slow to compile)
#include "slam/SE2_Types.h"

/**
 *	@def __SE2_ENABLED
 *	@brief if defined, the solver specializations for SE(2) are compiled
 */
#define __SE2_ENABLED

int n_Run_SE2_Solver(TCommandLineArgs t_args) // throw(std::runtime_error, std::bad_alloc)
{
#ifdef __SE2_ENABLED
	typedef MakeTypelist2(CVertexPose2D, CVertexLandmark2D) TVertexTypelist;
	typedef MakeTypelist2(CEdgePose2D, CEdgePoseLandmark2D) TEdgeTypelist;
	// define types of vertices, edges

	typedef CFlatSystem<CSEBaseVertex, TVertexTypelist,
		/*CVertexTypeTraits,*/ CSEBaseEdge, TEdgeTypelist> CSystemType;
	// make a system permitting SE(2) vertex and edge types

	typedef CSolverCaller<CSystemType, CSE2EdgeTraits, CParseLoop> CSpecializedSolverCaller;
	// specify how the nonlinear solver should be called

	return CTypelistForEach<CCompiledSolverList,
		CSpecializedSolverCaller>::Run(CSpecializedSolverCaller(t_args)).n_Result();
#else // __SE2_ENABLED
	fprintf(stderr, "error: SE(2) is disabled\n");
	return -1;
#endif // __SE2_ENABLED
}
