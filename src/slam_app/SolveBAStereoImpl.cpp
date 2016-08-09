/*
								+----------------------------------+
								|                                  |
								|  ***  BA solvers instances  ***  |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2014  |
								|                                  |
								|      SolveBAStereoImpl.cpp       |
								|                                  |
								+----------------------------------+
*/

/**
 *	@file src/slam_app/SolveBAStereoImpl.cpp
 *	@brief contains instantiation of stereo bundle adjustment solver templates
 *	@author -tHE SWINe-
 *	@date 2014-05-23
 */

/**
 *	@def __BA_STEREO_ENABLED
 *	@brief if defined, the solver specializations for stereo BA are compiled
 */
#define __BA_STEREO_ENABLED

#include "slam_app/Main.h"
#ifdef __BA_STEREO_ENABLED
#include "slam/ConfigSolvers.h" // only included in files that actually need the solvers (slow to compile)
#include "slam/BA_Types.h"
#endif // __BA_STEREO_ENABLED

int n_Run_BA_Stereo_Solver(TCommandLineArgs t_args) // throw(std::runtime_error, std::bad_alloc)
{
#ifdef __BA_STEREO_ENABLED
	_ASSERTE(t_args.b_use_BAS);

	typedef MakeTypelist_Safe((CVertexSCam, CVertexXYZ)) TVertexTypelist_BA;
	typedef MakeTypelist_Safe((CEdgeP2SC3D)) TEdgeTypelist_BA;
	// define types of vertices, edges

	typedef CFlatSystem<CBaseVertex, TVertexTypelist_BA, CEdgeP2SC3D, TEdgeTypelist_BA> CSystemType;
	// make a system permitting BA vertex and edge types

	if(t_args.n_solver_choice == nlsolver_Lambda)
		t_args.n_solver_choice = nlsolver_LambdaLM;
	// use Levenberg-Marquardt for bundle adjustment, GN not good enough

	typedef CSolverCaller<CSystemType, CBASEdgeTraits,
		CBASVertexTraits, CParseLoop> CSpecializedSolverCaller;
	// specify how the nonlinear solver should be called

	return CTypelistForEach<CCompiledSolverList,
		CSpecializedSolverCaller>::Run(CSpecializedSolverCaller(t_args)).n_Result();
#else // __BA_STEREO_ENABLED
	fprintf(stderr, "error: stereo BA is disabled\n");
	return -1;
#endif // __BA_STEREO_ENABLED
}
