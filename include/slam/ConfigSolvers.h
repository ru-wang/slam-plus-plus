/*
								+-----------------------------------+
								|                                   |
								|  ***  Solvers configuration  ***  |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2013  |
								|                                   |
								|          ConfigSolvers.h          |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __SOLVERS_CONFIGURATION_INCLUDED
#define __SOLVERS_CONFIGURATION_INCLUDED

/**
 *	@file include/slam/ConfigSolvers.h
 *	@author -tHE SWINe-
 *	@brief a file, containing solvers configuration (intended for enabling / disabling solvers)
 *	@date 2013
 */

#ifndef __SLAMPP_CONFIGURATION_INCLUDED
#define __SE_TYPES_SUPPORT_A_SOLVERS
#define __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
#define __SE_TYPES_SUPPORT_L_SOLVERS
#endif // !__SLAMPP_CONFIGURATION_INCLUDED
// if the config is not included, just enable everything (default)

/**
 *	@def __SLAM_COUNT_ITERATIONS_AS_VERTICES
 *	@brief if defined, the incremental updates are scheduled per N vertices,
 *		if not defined, they are scheduled per N edges (presumably faster in batch each 1,
 *		practically slower in batch each 10)
 */
#if !defined(__SLAM_COUNT_ITERATIONS_AS_EDGES) && !defined(__SLAM_COUNT_ITERATIONS_AS_VERTICES)
#define __SLAM_COUNT_ITERATIONS_AS_VERTICES
#endif // !__SLAM_COUNT_ITERATIONS_AS_EDGES && !__SLAM_COUNT_ITERATIONS_AS_VERTICES

#include "slam/NonlinearSolver_A.h"
#include "slam/NonlinearSolver_Lambda.h"
#include "slam/NonlinearSolver_Lambda_LM.h"
//#include "slam/NonlinearSolver_L.h"
#include "slam/NonlinearSolver_FastL.h"
// hint - disable or enable solvers at will

#include "slam/SolverTraits.h"
// some helper classes

/**
 *	@def _TySolverA_Name
 *	@brief name of the nonlinear solver A (not using typedef as it is a template)
 */
#ifdef __SE_TYPES_SUPPORT_A_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED
#define _TySolverA_Name CNonlinearSolver_A
#else // __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED
#define _TySolverA_Name CSolverNotIncluded
#endif // __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED
#else // __SE_TYPES_SUPPORT_A_SOLVERS
#define _TySolverA_Name CSolverNotSupported
#endif // __SE_TYPES_SUPPORT_A_SOLVERS

/**
 *	@def _TySolverLambda_Name
 *	@brief name of the Lambda nonlinear solver (not using typedef as it is a template)
 */
#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
#define _TySolverLambda_Name CNonlinearSolver_Lambda
#else // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
#define _TySolverLambda_Name CSolverNotIncluded
#endif // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
#else // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
#define _TySolverLambda_Name CSolverNotSupported
#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS

/**
 *	@def _TySolverLambdaLM_Name
 *	@brief name of the Lambda nonlinear solver with LM (not using typedef as it is a template)
 */
#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
#define _TySolverLambdaLM_Name CNonlinearSolver_Lambda_LM
#else // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
#define _TySolverLambdaLM_Name CSolverNotIncluded
#endif // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
#else // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
#define _TySolverLambdaLM_Name CSolverNotSupported
#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS

/**
 *	@def _TySolverL_Name
 *	@brief name of the nonlinear solver L (not using typedef as it is a template)
 */
#ifdef __SE_TYPES_SUPPORT_L_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_L_INCLUDED
#define _TySolverL_Name CNonlinearSolver_L
#else // __NONLINEAR_BLOCKY_SOLVER_L_INCLUDED
#define _TySolverL_Name CSolverNotIncluded
#endif // __NONLINEAR_BLOCKY_SOLVER_L_INCLUDED
#else // __SE_TYPES_SUPPORT_L_SOLVERS
#define _TySolverL_Name CSolverNotSupported
#endif // __SE_TYPES_SUPPORT_L_SOLVERS

/**
 *	@def _TySolverFastL_Name
 *	@brief name of the FastL nonlinear solver (not using typedef as it is a template)
 */
#ifdef __SE_TYPES_SUPPORT_L_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED
#define _TySolverFastL_Name CNonlinearSolver_FastL
#else // __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED
#define _TySolverFastL_Name CSolverNotIncluded
#endif // __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED
#else // __SE_TYPES_SUPPORT_L_SOLVERS
#define _TySolverFastL_Name CSolverNotSupported
#endif // __SE_TYPES_SUPPORT_L_SOLVERS

/**
 *	@brief nonlinear solver type
 */
enum ENonlinearSolverType {
	nlsolver_A, /**< @brief nonlinear solver A */
	nlsolver_Lambda, /**< @brief nonlinear solver lambda */
	nlsolver_LambdaLM, /**< @brief nonlinear solver lambda with Levenberg-Marquardt */
	nlsolver_L, /**< @brief nonlinear solver L */
	nlsolver_FastL /**< @brief nonlinear progressively reordering solver L */
};

/**
 *	@brief a list of all the available and unavailable solvers
 */
typedef MakeTypelist_Safe((
	CSolverTypeIdPair<nlsolver_A, _TySolverA_Name>,
	CSolverTypeIdPair<nlsolver_Lambda, _TySolverLambda_Name>,
	CSolverTypeIdPair<nlsolver_LambdaLM, _TySolverLambdaLM_Name>,
	CSolverTypeIdPair<nlsolver_L, _TySolverL_Name>,
	CSolverTypeIdPair<nlsolver_FastL, _TySolverFastL_Name>
	)) CCompiledSolverList;
// hint - add new solvers here

#endif // __SOLVERS_CONFIGURATION_INCLUDED
