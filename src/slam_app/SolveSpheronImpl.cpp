/*
								+-----------------------------------+
								|                                   |
								|   ***  Spherical BA solver  ***   |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2013  |
								|                                   |
								|       SolveSpheronImpl.cpp        |
								|                                   |
								+-----------------------------------+
*/

/**
 *	@file src/slam_app/SolveSpheronImpl.cpp
 *	@brief contains instantiation of spherical bundle adjustment solver templates
 *	@author isolony
 *	@date 2014-07-02
 */

/**
 *	@def __SPHERON_ENABLED
 *	@brief if defined, the solver specializations for BA are compiled
 */
#define __SPHERON_ENABLED

#include "slam_app/Main.h"
#ifdef __SPHERON_ENABLED
#include "slam/ConfigSolvers.h" // only included in files that actually need the solvers (slow to compile)
#include "slam/BA_Types.h"
//#include "slam/SE3_Types.h"
#endif // __SPHERON_ENABLED

int n_Run_Spheron_Solver(TCommandLineArgs t_args) // throw(std::runtime_error, std::bad_alloc)
{
#ifdef __SPHERON_ENABLED
	typedef MakeTypelist_Safe((CVertexSpheron, CVertexXYZ)) TVertexTypelist_BA;
	typedef MakeTypelist_Safe((CEdgeSpheronXYZ)) TEdgeTypelist_BA;
	// define types of vertices, edges

	typedef CFlatSystem<CBaseVertex, TVertexTypelist_BA, CEdgeSpheronXYZ, TEdgeTypelist_BA> CSystemType;
	// make a system permitting BA vertex and edge types

	if(t_args.n_solver_choice == nlsolver_Lambda)
		t_args.n_solver_choice = nlsolver_LambdaLM;
	// use Levenberg-Marquardt for bundle adjustment, GN not good enough

	typedef CSolverCaller<CSystemType, CSpheronEdgeTraits,
		CSpheronVertexTraits, CParseLoop> CSpecializedSolverCaller;
	// specify how the nonlinear solver should be called

	return CTypelistForEach<CCompiledSolverList,
		CSpecializedSolverCaller>::Run(CSpecializedSolverCaller(t_args)).n_Result();
#else // __SPHERON_ENABLED
	fprintf(stderr, "error: Spheron is disabled\n");
	return -1;
#endif // __SPHERON_ENABLED
}
