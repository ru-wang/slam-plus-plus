/*
								+-----------------------------------+
								|                                   |
								|      ***  Configuration  ***      |
								|                                   |
								|   Copyright  � -tHE SWINe- 2012   |
								|                                   |
								|             Config.h              |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __PRECOMPILED_HEADER_INCLUDED
#define __PRECOMPILED_HEADER_INCLUDED

/**
 *	@file include/slam/Config.h
 *	@author -tHE SWINe-
 *	@brief main application configuration file
 *	@note Also a precompiled header (makes the compilation faster in windows).
 *	@date 2012
 */

#include <csparse/cs.hpp>

/*#define EIGEN_NO_DEBUG
#define EIGEN_NO_STATIC_ASSERT*/
// to compile with the experimental allocator

#include "eigen/Eigen/Dense"
#include <stdio.h>
#include <time.h>
#if defined(_WIN32) || defined(_WIN64)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <shellapi.h> // ShellExecute()
#endif // _WIN32 || _WIN64
#include "slam/Debug.h"
#include "slam/Unused.h"
#include "slam/Parser.h"
#include "slam/Timer.h"

/**
 *	@def __PIMP_MY_WINDOW
 *	@brief if defined, looks if there is a second monitor and moves the console window there
 *	@note This only works on windows platforms.
 */
#define __PIMP_MY_WINDOW

/**
 *	@def __COMPILE_LIBRARY_TESTS
 *	@brief if defined, csparse and eigen test functions are compiled (use to debug linking problems)
 */
//#define __COMPILE_LIBRARY_TESTS

/**
 *	@def __USE_NATIVE_CHOLESKY
 *	@brief if defined, the native linear solver is used instead of CSparse
 */
//#define __USE_NATIVE_CHOLESKY

/**
 *	@def __USE_CHOLMOD
 *	@brief if defined, CHOLMOD is used as linear solver instead of CSparse
 *		(only if __USE_NATIVE_CHOLESKY is not defined)
 */
//#define __USE_CHOLMOD

/**
 *	@def __CHOLMOD_SHORT
 *	@brief if defined, CHOLMOD is forced to use 32-bit indices (it is only effective when
 *		in x64 mode, if __USE_CHOLMOD is defined and __USE_NATIVE_CHOLESKY is not)
 */
//#define __CHOLMOD_SHORT

/**
 *	@def __USE_CXSPARSE
 *	@brief if defined, CXSparse is used as linear solver instead of CSparse
 *		(only if __USE_CHOLMOD or __USE_NATIVE_CHOLESKY is not defined)
 */
//#define __USE_CXSPARSE

/**
 *	@def __CXSPARSE_SHORT
 *	@brief if defined, CXSparse is forced to use 32-bit indices (it is only effective when in x64
 *		mode, if __USE_CXSPARSE is defined and __USE_CHOLMOD or __USE_NATIVE_CHOLESKY are not)
 */
//#define __CXSPARSE_SHORT

/**
 *	@def __SLAM_COUNT_ITERATIONS_AS_VERTICES
 *	@brief if defined, the incremental updates are scheduled per N vertices,
 *		if not defined, they are scheduled per N edges (presumably faster in batch each 1,
 *		practically slower in batch each 10)
 */
//#define __SLAM_COUNT_ITERATIONS_AS_VERTICES

/**
 *	@def __LINEAR_SOLVER_OVERRIDE
 *	@brief enables linear solver selection from commandline (0 = CSparse, 1 = CXSparse, 2 = CHOLMOD, 3 = native)
 */
#ifdef __LINEAR_SOLVER_OVERRIDE
//#pragma message("__LINEAR_SOLVER_OVERRIDE = " __LINEAR_SOLVER_OVERRIDE)
#ifdef __USE_CHOLMOD
#undef __USE_CHOLMOD
#endif // __USE_CHOLMOD
#ifdef __USE_CXSPARSE
#undef __USE_CXSPARSE
#endif // __USE_CXSPARSE
#if __LINEAR_SOLVER_OVERRIDE == 0
// use CSparse
#elif __LINEAR_SOLVER_OVERRIDE == 1
#define __USE_CXSPARSE
#elif __LINEAR_SOLVER_OVERRIDE == 2
#define __USE_CHOLMOD
#else // __LINEAR_SOLVER_OVERRIDE == 0

#endif // __LINEAR_SOLVER_OVERRIDE == 0
#endif // __LINEAR_SOLVER_OVERRIDE
// add override option to easily compile from linux

/**
 *	@def __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_A
 *	@brief if defined, SE(2) types implement functions, required by CNonlinearSolver_A
 */
#define __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_A

/**
 *	@def __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_LAMBDA
 *	@brief if defined, SE(2) types implement functions, required
 *		by CNonlinearSolver_L and CNonlinearSolver_Lambda
 */
#define __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_LAMBDA

/**
 *	@def __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_L
 *	@brief if defined, SE(2) types implement functions, required by CNonlinearSolver_L
 */
#define __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_L

#include "slam/Segregated.h"
#include "slam/BlockMatrix.h"
//#include "slam/Tetris.h" // to speed up builds
#include "slam/BlockBench.h" // don't need this right now
#include "slam/BlockUnit.h"
#include "slam/FlatSystem.h"
#include "slam/ParseLoop.h"
#include "slam/SE2_Types.h"
#include "slam/LinearSolver_CSparse.h"
#include "slam/LinearSolver_CXSparse.h"
#include "slam/LinearSolver_CholMod.h"
// rarely change (pch optimization)

#endif // __PRECOMPILED_HEADER_INCLUDED
