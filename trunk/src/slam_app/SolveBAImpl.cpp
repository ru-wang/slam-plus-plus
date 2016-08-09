/*
								+----------------------------------+
								|                                  |
								|  ***  BA solvers instances  ***  |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2013  |
								|                                  |
								|         SolveBAImpl.cpp          |
								|                                  |
								+----------------------------------+
*/

/**
 *	@file src/slam_app/SolveBAImpl.cpp
 *	@brief contains instantiation of bundle adjustment solver templates
 *	@author -tHE SWINe-
 *	@date 2013-06-14
 */

/**
 *	@def __BA_ENABLED
 *	@brief if defined, the solver specializations for BA are compiled
 */
#define __BA_ENABLED

#include "slam_app/Main.h"
#ifdef __BA_ENABLED
#include "slam/ConfigSolvers.h" // only included in files that actually need the solvers (slow to compile)
#include "slam/BA_Types.h"
#include "slam/SE3_Types.h"
#endif // __BA_ENABLED

int n_Run_BA_Solver(TCommandLineArgs t_args) // throw(std::runtime_error, std::bad_alloc)
{
#ifdef __BA_ENABLED
	_ASSERTE(t_args.b_use_BA);
	typedef MakeTypelist_Safe((CVertexCam, CVertexXYZ)) TVertexTypelist_BA;
	typedef MakeTypelist_Safe((CEdgeP2C3D)) TEdgeTypelist_BA;
	// define types of vertices, edges

	typedef CFlatSystem<CBaseVertex, TVertexTypelist_BA, CEdgeP2C3D, TEdgeTypelist_BA> CSystemType;
	// make a system permitting BA vertex and edge types

	if(t_args.n_solver_choice == nlsolver_Lambda)
		t_args.n_solver_choice = nlsolver_LambdaLM;
	// use Levenberg-Marquardt for bundle adjustment, GN not good enough

	typedef CSolverCaller<CSystemType, CBAEdgeTraits,
		CBAVertexTraits, CParseLoop> CSpecializedSolverCaller;
	// specify how the nonlinear solver should be called

	return CTypelistForEach<CCompiledSolverList,
		CSpecializedSolverCaller>::Run(CSpecializedSolverCaller(t_args)).n_Result();
#else // __BA_ENABLED
	fprintf(stderr, "error: BA is disabled\n");
	return -1;
#endif // __BA_ENABLED
}
