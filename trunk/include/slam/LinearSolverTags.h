/*
								+----------------------------------+
								|                                  |
								|   ***  Linear solver tags  ***   |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2012  |
								|                                  |
								|        LinearSolverTags.h        |
								|                                  |
								+----------------------------------+
*/

#pragma once
#ifndef __LINEAR_SOLVER_TAGS_INCLUDED
#define __LINEAR_SOLVER_TAGS_INCLUDED

/**
 *	@file include/slam/LinearSolverTags.h
 *	@brief linear solver tags
 *	@author -tHE SWINe-
 *	@date 2012-09-04
 *
 *	Each solver class needs to have public type named _Tag,
 *	which is one of the types defined below. Based on that,
 *	the nonlinear solvers can then decide which functions
 *	the linear solver implements and what is the best
 *	approach to use.
 *
 */

/**
 *	@brief default tag for all the linear solvers
 */
class CBasicLinearSolverTag {};

/**
 *	@brief solvers supporting blockwise solutions
 *
 *	The usage of such solver is supposed to be:
 *
 *@code
 *	UberBlockMatrix A;
 *	Eigen::VectorXd b;
 *	Solver s;
 *
 *	s.SymbolicDecomposition_Blocky(A);
 *	s.Solve_PosDef_Blocky(A, b);
 *@endcode
 */
class CBlockwiseLinearSolverTag {};

#endif // __LINEAR_SOLVER_TAGS_INCLUDED
