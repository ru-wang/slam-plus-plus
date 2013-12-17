/*
								+-----------------------------------+
								|                                   |
								| *** L factor nonlinear solver *** |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2013  |
								|                                   |
								|      NonlinearSolver_FastL.h      |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED
#define __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED

/**
 *	@file include/slam/NonlinearSolver_FastL.h
 *	@brief nonlinear blocky solver with progressive reordering, working above the L factor matrix
 *	@author -tHE SWINe-
 *	@date 2013-01-28
 */

#include "slam/FlatSystem.h"
#include "slam/OrderingMagic.h"
#include "slam/IncrementalPolicy.h"
#include "slam/Marginals.h"

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
 *	@brief if defined, solution is calcualted at each step, otherwise it is calculated
 *		only at the specified intervals (but is ready to be calculated at each step quickly)
 */
#define __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
 *	@brief if defined, R is always used to calculate updates, even in the subsequent iterations,
 *		so that when the solver finishes, R is up to date and does not need to be refreshed again
 */
#define __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
 *	@brief if defined, permutation folding is verified by comparing parts
 *		of lambda (lambda00 and lambda11)
 *	@note This slows down quite a bit, should not be enabled for production builds.
 */
//#define __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
 *	@brief enables writes of diagnostic data (timing samples, ...)
 */
//#define __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2
 *	@brief enables writes of chi2 errors at each step
 */
//#define __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE
 *	@brief if defined, chi2 is calculated at the last edge, just before introducing a new vertex,
 *		that gives a different chi2 than calculating chi2 just after introducing a new vertex (default)
 */
//#define __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY
 *	@brief enables writes of density of R given different ordering strategies
 */
//#define __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
 *	@brief if defined, enables writes of timing of different R update variants
 */
//#define __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
 *	@brief replaces timer implementation with dummy timer to avoid costy runtime
 *		library calls to get the time (more costy on windows than on linux)
 */
//#define __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING

#ifndef _DEBUG

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
 *	@brief enables computation of dense cholesky for small loops (a speed optimization)
 */
#define __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY

#endif // !_DEBUG

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_PARALLEL_MATMULT_THRESH
 *	@brief matrix size threshold for parallel multiplication
 *		in R update (in blocks, used to calculate \f$R_{10}R_{10}^T\f$ and \f$R_{11}R_{11}^T\f$)
 */
#define __NONLINEAR_SOLVER_FAST_L_PARALLEL_MATMULT_THRESH 200

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
 *	@brief dump RSS 2013 matrix animation data
 */
//#define __NONLINEAR_SOLVER_FAST_L_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

/**
 *	@brief namespace with fast L solver utilities
 */
namespace __fL_util {

/**
 *	@brief a simple native solver predicate
 *	@tparam CLinearSolver is linear solver type
 */
template <class CLinearSolver>
class CIsNativeSolver {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		b_result = false /**< @brief result of comparison */
	};
};

/**
 *	@brief a simple native solver predicate (specialization for native solver)
 *	@tparam CBlockMatrixSizeList is list of block matrix sized
 */
template <class CBlockMatrixSizeList>
class CIsNativeSolver<CLinearSolver_UberBlock<CBlockMatrixSizeList> > {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		b_result = true /**< @brief result of comparison */
	};
};

} // ~__fL_util

/**
 *	@brief nonlinear blocky solver working above the L factor matrix
 *
 *	@tparam CSystem is the system type
 *	@tparam CLinearSolver is a linear solver
 */
template <class CSystem, class CLinearSolver, class CAMatrixBlockSizes = typename CSystem::_TyJacobianMatrixBlockList>
class CNonlinearSolver_FastL {
public:
	typedef CSystem _TySystem; /**< @brief system type */
	typedef CLinearSolver _TyLinearSolver; /**< @brief linear solver type */

	typedef typename CSystem::_TyBaseVertex _TyBaseVertex; /**< @brief the data type for storing vertices */
	typedef typename CSystem::_TyVertexTypelist _TyVertexTypelist; /**< @brief list of vertex types */
	typedef typename CSystem::_TyBaseEdge _TyBaseEdge; /**< @brief the data type for storing measurements */
	typedef typename CSystem::_TyEdgeTypelist _TyEdgeTypelist; /**< @brief list of edge types */

	typedef typename CSystem::_TyVertexMultiPool _TyVertexMultiPool; /**< @brief vertex multipool type */
	typedef typename CSystem::_TyEdgeMultiPool _TyEdgeMultiPool; /**< @brief edge multipool type */

	typedef typename CLinearSolver::_Tag _TySolverTag; /**< @brief linear solver tag */
	typedef CLinearSolverWrapper<_TyLinearSolver, _TySolverTag> _TyLinearSolverWrapper; /**< @brief wrapper for linear solvers (shields solver capability to solve blockwise) */

	typedef typename CUniqueTypelist<CAMatrixBlockSizes>::_TyResult _TyAMatrixBlockSizes; /**< @brief possible block matrices, that can be found in A */
	typedef typename __fbs_ut::CBlockMatrixTypesAfterPreMultiplyWithSelfTranspose<
		_TyAMatrixBlockSizes>::_TyResult _TyLambdaMatrixBlockSizes; /**< @brief possible block matrices, found in lambda and R */

	/**
	 *	@brief some run-time constants, stored as enum
	 */
	enum {
		b_Is_PoseOnly_SLAM = CTypelistLength<_TyAMatrixBlockSizes>::n_result == 1, /**< @brief determines if we're doing pose-only SLAM (10k) */
		b_Have_NativeSolver = __fL_util::CIsNativeSolver<_TyLinearSolver>::b_result /**< @brief determines if the native linear solver is being used */
	};

	/**
	 *	@brief solver interface properties, stored as enum (see also CSolverTraits)
	 */
	enum {
		solver_HasDump = true, /**< @brief timing statistics support flag */
		solver_HasChi2 = true, /**< @brief Chi2 error calculation support flag */
		solver_HasMarginals = true, /**< @brief marginal covariance support flag */
		solver_HasGaussNewton = true, /**< @brief Gauss-Newton support flag */
		solver_HasLevenberg = false, /**< @brief Levenberg-Marquardt support flag */
		solver_HasGradient = false, /**< @brief gradient-based linear solving support flag */
		solver_HasSchur = false, /**< @brief Schur complement support flag */
		solver_HasDelayedOptimization = true, /**< @brief delayed optimization support flag */
		solver_IsPreferredBatch = false, /**< @brief preferred batch solver flag */
		solver_IsPreferredIncremental = true, /**< @brief preferred incremental solver flag */
		solver_ExportsJacobian = false, /**< @brief interface for exporting jacobian system matrix flag */
		solver_ExportsHessian = false, /**< @brief interface for exporting hessian system matrix flag */
		solver_ExportsFactor = false /**< @brief interface for exporting factorized system matrix flag */
	};

protected:
	CSystem &m_r_system; /**< @brief reference to the system */
	CLinearSolver m_linear_solver; /**< @brief linear solver */
	CLinearSolver m_linear_solver2; /**< @brief linear solver for calculating cholesky of R and cholesky of R increment */

	std::vector<size_t> m_chol_etree; /**< @brief reusable e-tree storage */
	std::vector<size_t> m_chol_ereach_stack; /**< @brief reusable workspace for Cholesky */
	std::vector<size_t> m_chol_bitfield; /**< @brief reusable workspace for Cholesky */

	CUberBlockMatrix m_R; /**< @brief the R matrix (built / updated incrementally) */

	bool m_b_had_loop_closure; /**< @brief (probable) loop closure flag */
	bool m_b_first_iteration_use_R; /**< @brief flag for using the R matrix or rather lambda in the first iteration of nonlinear optimization */
	bool m_b_R_up_to_date; /**< @brief dirty flag for the R matrix (required to keep track after lambda updates and linearization point changes) */
	size_t m_n_last_full_R_update_size; /**< @brief the last number of block columns in R when it was fully updated */
	std::vector<size_t> m_R_row_lookup_table; /**< @brief row lookup table for R (used by b_Refresh_R() and Refresh_R11()) */
	//size_t m_n_big_loop_threshold; /**< @brief threshold for what is considered a "big" loop (incrementing R is avoided) */

	CMatrixOrdering m_lambda_ordering; /**< @brief lambda block ordering calculator (CAMD wrapper) */
	const size_t *m_p_lambda_block_ordering; /**< @brief lambda block ordering (only valid if m_b_R_up_to_date is set) */ // todo - convert all those to size_t
	size_t m_n_lambda_block_ordering_size; /**< @brief lambda block ordering size */
	CUberBlockMatrix m_lambda_perm; /**< @brief the reordered reference to the lambda matrix */
	CUberBlockMatrix m_lambda; /**< @brief the lambda matrix (built / updated incrementally) */

	CMatrixOrdering m_lambda11_ordering; /**< @brief lambda11 block ordering calculator (CAMD wrapper) */
	const size_t *m_p_lambda11_block_ordering; /**< @brief lambda block ordering (only valid if m_b_R_up_to_date is set) */ // todo - convert all those to size_t
	size_t m_n_lambda_block11_ordering_size; /**< @brief lambda block ordering size */

	CFirstLastElementOrderingConstraint m_lambda11_constraint; /**< @brief incremental lambda ordering constraint */
	CLastElementOrderingConstraint m_lambda_constraint; /**< @brief global lambda ordering constraint */
	CMatrixOrdering m_lambda_alt_ordering; /**< @brief secondary lambda ordering, calculated from m_lambda_perm */
	CNFirst1LastElementOrderingConstraint m_lambda_alt_constraint; /**< @brief constraint for the secondary lambda ordering */

	Eigen::VectorXd m_v_dx; /**< @brief dx vector */
	Eigen::VectorXd m_v_d; /**< @brief d vector */
	Eigen::VectorXd m_v_perm_temp; /**< @brief temporary storage for the permutation vector, same dimension as d and dx */
	size_t m_n_verts_in_lambda; /**< @brief number of vertices already in lambda */
	size_t m_n_edges_in_lambda; /**< @brief number of edges already in lambda */
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
	size_t m_n_last_optimized_vertex_num;
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
	size_t m_n_step; /**< @brief counter of incremental steps modulo m_n_nonlinear_solve_threshold */
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
	size_t m_n_linear_solve_threshold; /**< @brief step threshold for linear solve */
	size_t m_n_nonlinear_solve_threshold; /**< @brief step threshold for nonlinear solve */
	size_t m_n_nonlinear_solve_max_iteration_num; /**< @brief maximal number of iterations in incremental nonlinear solve */
	double m_f_nonlinear_solve_error_threshold; /**< @brief error threshold in incremental nonlinear solve */ // t_odo - document these in elementwise A and R
	bool m_b_verbose; /**< @brief verbosity flag */

	size_t m_n_real_step; /**< @brief counter of incremental steps (no modulo) */

	bool m_b_system_dirty; /**< @brief system updated without relinearization flag */
	bool m_b_linearization_dirty; /**< @brief system matrices updated but relinearization point was not set flag */

#ifndef __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
	typedef CVoidTimerSampler _TyTimeSampler; /**< @brief timer sampler type */
#else // !__NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
	typedef CTimerSampler _TyTimeSampler; /**< @brief timer sampler type */
#endif // !__NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
	typedef _TyTimeSampler::_TySample _TyTime; /**< @brief time type */

	bool m_b_inhibit_optimization; /**< @brief optimization enable falg */
	size_t m_n_iteration_num; /**< @brief number of linear solver iterations */
	_TyTime m_f_chol_time; /**< @brief time spent in Choleski() section */
	_TyTime m_f_norm_time; /**< @brief time spent in norm calculation section */
	_TyTime m_f_vert_upd_time; /**< @brief time spent in updating the vertices */

	size_t m_n_full_forwardsubst_num; /**< @brief number of d updates performed using full R forward substitution */
	size_t m_n_resumed_forwardsubst_num; /**< @brief number of d updates performed using resumed R forward substitution */
	size_t m_n_resumed_perm_forwardsubst_num; /**< @brief number of d updates performed using resumed R forward substitution, while being permutated in the updated area */
	size_t m_n_R_optim_num; /**< @brief number of system optimizations performed using R backsubstitution */
	size_t m_n_lambda_optim_num; /**< @brief number of system optimizations performed using cholsol(lambda) */
	size_t m_n_Rup_num; /**< @brief number of R increments */
	size_t m_n_omega_update_num; /**< @brief number of R increments calculated using omega */
	size_t m_n_lambda_update_num; /**< @brief number of R increments calculated using lambda */
	size_t m_n_full_R_num; /**< @brief number of R updates */
	_TyTime m_f_lambda_refresh_time; /**< @brief time spent in updating and allocating lambda */
	_TyTime m_f_rhs_time; /**< @brief time spent in updating right-hand side vector */
	_TyTime m_f_ordering_time; /**< @brief time spent calculating ordering of lambda */
	_TyTime m_f_fullR_d; /**< @brief time spent in updating d while doing full R */
	_TyTime m_f_r11_omega_calc_time; /**< @brief time spent calculating omega (R increment) */
	_TyTime m_f_r11_omega_slice_time; /**< @brief time spent in slicing \f$R_{11}\f$ (R increment) */
	_TyTime m_f_r11_omega_ata_time; /**< @brief time spent calculating \f$R_{11}^TR_{11}\f$ (R increment) */
	_TyTime m_f_r11_omega_add_time; /**< @brief time spent adding \f$R_{11}^TR_{11}\f$ + omega (R increment) */
	_TyTime m_f_r11_lambda_slice_time; /**< @brief time spent in slicing lambda11 and \f$R_{01}\f$ (R increment) */
	_TyTime m_f_r11_lambda_ata_time; /**< @brief time spent calculating \f$R_{01}^TR_{01}\f$ (R increment) */
	_TyTime m_f_r11_lambda_add_time; /**< @brief time spent adding \f$R_{01}^TR_{01} + \Lambda_{11}\f$ (R increment) */
	_TyTime m_f_Rupdate_time; /**< @brief time spent calculating cholesky of new \f$R_{11}\f$ (R increment) */
	_TyTime m_f_d_time; /**< @brief time spent updating d (right hand side vector) */
	_TyTime m_f_backsubst_time; /**< @brief time spent in backsubstitution (solving for R / d) */
	_TyTime m_f_fullR_cholesky; /**< @brief time spent in calculating cholesky (R update) */
				
	size_t m_n_resumed_chol_num; /**< @brief number of times the resumed Cholesky was used */
	size_t m_n_blocks_above_num; /**< @brief number of times there were blocks above lambda_11 */
	size_t m_n_limited_search_num; /**< @brief number of times there were blocks above lambda_11 but only a smaller submatrix was sufficient for permutation calculation */
	_TyTime m_f_ordering_fold_time; /**< @brief time spent folding two orderings */
	_TyTime m_f_repermute_time; /**< @brief time spent repermuting lambda matrix with incremented ordering */
	_TyTime m_f_Rslice_time; /**< @brief time spent slicing R for resumed Cholesky */
	_TyTime m_f_etree_time; /**< @brief time spent calculating the elimination tree */
	_TyTime m_f_resumed_chol_time; /**< @brief time spent in resumed Cholesky */
	_TyTime m_f_ordering11_time; /**< @brief time spent in calculating the incremental ordering */
	_TyTime m_f_ordering11_part_time; /**< @brief time spent in calculating the incremental ordering, only the small or inflated lambda_11 cases */
	_TyTime m_f_ordering11_full_time; /**< @brief time spent in calculating the incremental ordering, only the full lambda_perm cases */

	std::vector<size_t> lambda_perm_frontline; /**< @brief cached frontline of the lambda_perm matrix */

	TMarginalsComputationPolicy m_t_marginals_config; /**< @brief marginal covariance computation configuration */
	_TyTime m_f_marginals_time; /**< @brief time spent in calculating marginal covariances */
	CMarginalCovariance m_marginals; /***< @brief marginals cache */

	CTimer m_shared_timer; /**< @brief timer object */

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
	size_t m_n_loop_size_cumsum; /**< @brief cumulative sum of loops processed so far */
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS

public:
	/**
	 *	@brief initializes the nonlinear solver
	 *
	 *	@param[in] r_system is the system to be optimized
	 *		(it is only referenced, not copied - must not be deleted)
	 *	@param[in] n_linear_solve_threshold is the step threshold
	 *		for linear solver to be called (0 = disable)
	 *	@param[in] n_nonlinear_solve_threshold is the step threshold
	 *		for nonlinear solver to be called (0 = disable)
	 *	@param[in] n_nonlinear_solve_max_iteration_num is maximal
	 *		number of iterations in nonlinear solver
	 *	@param[in] f_nonlinear_solve_error_threshold is error threshold
	 *		for the nonlinear solver
	 *	@param[in] b_verbose is verbosity flag
	 *	@param[in] linear_solver is linear solver instance
	 *	@param[in] b_use_schur is Schur complement trick flag (not supported)
	 *
	 *	@deprecated This is deprecated version of the constructor, use constructor
	 *		with TIncrementalSolveSetting instead.
	 */
	CNonlinearSolver_FastL(CSystem &r_system, size_t n_linear_solve_threshold,
		size_t n_nonlinear_solve_threshold, size_t n_nonlinear_solve_max_iteration_num = 5,
		double f_nonlinear_solve_error_threshold = .01, bool b_verbose = false,
		CLinearSolver linear_solver = CLinearSolver(), bool UNUSED(b_use_schur) = false)
		:m_r_system(r_system), m_linear_solver(linear_solver),
		m_linear_solver2(linear_solver), m_b_had_loop_closure(false),
		m_b_first_iteration_use_R(true), m_b_R_up_to_date(true), m_n_last_full_R_update_size(0),
		m_p_lambda_block_ordering(0), m_n_lambda_block_ordering_size(0),
		m_p_lambda11_block_ordering(0), m_n_lambda_block11_ordering_size(0),
		m_n_verts_in_lambda(0), m_n_edges_in_lambda(0),
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_last_optimized_vertex_num(0),
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_step(0),
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_linear_solve_threshold(n_linear_solve_threshold),
		m_n_nonlinear_solve_threshold(n_nonlinear_solve_threshold),
		m_n_nonlinear_solve_max_iteration_num(n_nonlinear_solve_max_iteration_num),
		m_f_nonlinear_solve_error_threshold(f_nonlinear_solve_error_threshold),
		m_b_verbose(b_verbose), m_n_real_step(0), m_b_system_dirty(false),
		m_b_linearization_dirty(false), m_b_inhibit_optimization(false),
		m_n_iteration_num(0), m_f_chol_time(0), m_f_norm_time(0), m_f_vert_upd_time(0),
		m_n_full_forwardsubst_num(0), m_n_resumed_forwardsubst_num(0),
		m_n_resumed_perm_forwardsubst_num(0), m_n_R_optim_num(0), m_n_lambda_optim_num(0),
		m_n_Rup_num(0), m_n_omega_update_num(0), m_n_lambda_update_num(0), m_n_full_R_num(0),
		m_f_lambda_refresh_time(0), m_f_rhs_time(0), m_f_ordering_time(0), m_f_fullR_d(0),
		m_f_r11_omega_calc_time(0), m_f_r11_omega_slice_time(0), m_f_r11_omega_ata_time(0),
		m_f_r11_omega_add_time(0), m_f_r11_lambda_slice_time(0), m_f_r11_lambda_ata_time(0),
		m_f_r11_lambda_add_time(0), m_f_Rupdate_time(0), m_f_d_time(0),
		m_f_backsubst_time(0), m_f_fullR_cholesky(0),
		m_n_resumed_chol_num(0), m_n_blocks_above_num(0), m_n_limited_search_num(0),
		m_f_ordering_fold_time(0), m_f_repermute_time(0), m_f_Rslice_time(0), m_f_etree_time(0),
		m_f_resumed_chol_time(0), m_f_ordering11_time(0), m_f_ordering11_part_time(0),
		m_f_ordering11_full_time(0), m_f_marginals_time(0)
	{
		_ASSERTE(!m_n_nonlinear_solve_threshold || !m_n_linear_solve_threshold); // only one of those

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
		m_n_loop_size_cumsum = 0;
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
	}

	/**
	 *	@brief initializes the nonlinear solver
	 *
	 *	@param[in] r_system is the system to be optimized
	 *		(it is only referenced, not copied - must not be deleted)
	 *	@param[in] t_incremental_config is incremental solving configuration
	 *	@param[in] t_marginals_config is marginal covariance calculation configuration
	 *	@param[in] b_verbose is verbosity flag
	 *	@param[in] linear_solver is linear solver instance
	 *	@param[in] b_use_schur is Schur complement trick flag (not supported)
	 */
	CNonlinearSolver_FastL(CSystem &r_system,
		TIncrementalSolveSetting t_incremental_config = TIncrementalSolveSetting(),
		TMarginalsComputationPolicy t_marginals_config = TMarginalsComputationPolicy(),
		bool b_verbose = false,
		CLinearSolver linear_solver = CLinearSolver(), bool UNUSED(b_use_schur) = false)
		:m_r_system(r_system), m_linear_solver(linear_solver),
		m_linear_solver2(linear_solver), m_b_had_loop_closure(false),
		m_b_first_iteration_use_R(true), m_b_R_up_to_date(true), m_n_last_full_R_update_size(0),
		m_p_lambda_block_ordering(0), m_n_lambda_block_ordering_size(0),
		m_p_lambda11_block_ordering(0), m_n_lambda_block11_ordering_size(0),
		m_n_verts_in_lambda(0), m_n_edges_in_lambda(0),
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_last_optimized_vertex_num(0),
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_step(0),
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_linear_solve_threshold(t_incremental_config.t_linear_freq.n_period),
		m_n_nonlinear_solve_threshold(t_incremental_config.t_nonlinear_freq.n_period),
		m_n_nonlinear_solve_max_iteration_num(t_incremental_config.n_max_nonlinear_iteration_num),
		m_f_nonlinear_solve_error_threshold(t_incremental_config.f_nonlinear_error_thresh),
		m_b_verbose(b_verbose), m_n_real_step(0), m_b_system_dirty(false),
		m_b_linearization_dirty(false), m_b_inhibit_optimization(false),
		m_n_iteration_num(0), m_f_chol_time(0), m_f_norm_time(0), m_f_vert_upd_time(0),
		m_n_full_forwardsubst_num(0), m_n_resumed_forwardsubst_num(0),
		m_n_resumed_perm_forwardsubst_num(0), m_n_R_optim_num(0), m_n_lambda_optim_num(0),
		m_n_Rup_num(0), m_n_omega_update_num(0), m_n_lambda_update_num(0), m_n_full_R_num(0),
		m_f_lambda_refresh_time(0), m_f_rhs_time(0), m_f_ordering_time(0), m_f_fullR_d(0),
		m_f_r11_omega_calc_time(0), m_f_r11_omega_slice_time(0), m_f_r11_omega_ata_time(0),
		m_f_r11_omega_add_time(0), m_f_r11_lambda_slice_time(0), m_f_r11_lambda_ata_time(0),
		m_f_r11_lambda_add_time(0), m_f_Rupdate_time(0), m_f_d_time(0),
		m_f_backsubst_time(0), m_f_fullR_cholesky(0),
		m_n_resumed_chol_num(0), m_n_blocks_above_num(0), m_n_limited_search_num(0),
		m_f_ordering_fold_time(0), m_f_repermute_time(0), m_f_Rslice_time(0), m_f_etree_time(0),
		m_f_resumed_chol_time(0), m_f_ordering11_time(0), m_f_ordering11_part_time(0),
		m_f_ordering11_full_time(0), m_t_marginals_config(t_marginals_config), m_f_marginals_time(0)
	{
		_ASSERTE(!m_n_nonlinear_solve_threshold || !m_n_linear_solve_threshold); // only one of those

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
		m_n_loop_size_cumsum = 0;
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
	}

	/**
	 *	@brief destructor (only required if ? is defined)
	 */
	inline ~CNonlinearSolver_FastL()
	{}

	/**
	 *	@brief displays performance info on stdout
	 *	@param[in] f_total_time is total time taken by everything (can be -1 to ommit)
	 */
	void Dump(double f_total_time = -1) const
	{
		printf("solver took " PRIsize " iterations\n", m_n_iteration_num); // debug, to be able to say we didn't botch it numerically
#ifdef __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
		if(m_t_marginals_config.b_calculate)
			printf("solver spent %f seconds in calculating the marginals\n", m_f_marginals_time);

		double f_serial_time = m_f_backsubst_time + m_f_chol_time + m_f_norm_time + m_f_vert_upd_time;
		if(f_total_time > 0) {
			printf("solver spent %f seconds in parallelizable section (updating R)\n",
				f_total_time - f_serial_time);
		}
		double f_total_resumed_up = m_f_resumed_chol_time + m_f_Rslice_time +
			m_f_etree_time + m_f_ordering_fold_time + m_f_repermute_time;
		double f_total_omega_up = m_f_r11_omega_calc_time + m_f_r11_omega_slice_time +
			m_f_r11_omega_ata_time + m_f_r11_omega_add_time;
		double f_total_lambda_up = m_f_r11_lambda_slice_time +
			m_f_r11_lambda_ata_time + m_f_r11_lambda_add_time;
		double f_l_upd_time = f_total_resumed_up + f_total_omega_up +
			f_total_lambda_up + m_f_Rupdate_time + m_f_ordering11_time +
			m_f_ordering11_part_time + m_f_ordering11_full_time;
		double f_measured_parallel_time = m_f_lambda_refresh_time + m_f_rhs_time + m_f_ordering_time +
			m_f_fullR_d + m_f_fullR_cholesky + f_l_upd_time + m_f_d_time;
		printf("measured parallel time: %f, disparity: %f; out of which:\n", f_measured_parallel_time,
			f_total_time - f_serial_time - f_measured_parallel_time);
		printf("\t   ,\\: %f\n", m_f_lambda_refresh_time);
		printf("\t  rhs: %f\n", m_f_rhs_time);
		printf("\torder: %f\n", m_f_ordering_time);
		printf("\tfullR: %f (ran " PRIsize " times)\n", m_f_fullR_d + m_f_fullR_cholesky, m_n_full_R_num);
		printf("\tout of which:\n");
		printf("\t\t chol: %f\n", m_f_fullR_cholesky);
		printf("\t\t    d: %f\n", m_f_fullR_d);
		printf("\tR update: %f (ran " PRIsize " times)\n", f_l_upd_time, m_n_Rup_num);
		printf("\t\tordfu: %f (blocks above " PRIsize " times)\n", m_f_ordering11_full_time, m_n_blocks_above_num);
		printf("\t\tordli: %f (ran " PRIsize " times)\n", m_f_ordering11_part_time, m_n_limited_search_num);
		printf("\t\tordsm: %f (ran " PRIsize " times)\n", m_f_ordering11_time, m_n_Rup_num - m_n_blocks_above_num - m_n_limited_search_num);
		printf("\t\tresum: %f (ran " PRIsize " times)\n", f_total_resumed_up, m_n_resumed_chol_num);
		printf("\t\t\tofold: %f\n", m_f_ordering_fold_time);
		printf("\t\t\trperm: %f\n", m_f_repermute_time);
		printf("\t\t\tR cut: %f\n", m_f_Rslice_time);
		printf("\t\t\tetree: %f\n", m_f_etree_time);
		printf("\t\t\t chol: %f\n", m_f_resumed_chol_time);
		printf("\t\t  add: %f (ran " PRIsize " times)\n", f_total_omega_up + f_total_lambda_up +
			m_f_Rupdate_time, m_n_Rup_num - m_n_resumed_chol_num);
		printf("\t\t\tomega: %f (ran " PRIsize " times)\n", f_total_omega_up, m_n_omega_update_num);
		printf("\t\t\t\t calc: %f\n", m_f_r11_omega_calc_time);
		printf("\t\t\t\tslice: %f\n", m_f_r11_omega_slice_time);
		printf("\t\t\t\t Rata: %f\n", m_f_r11_omega_ata_time);
		printf("\t\t\t\tR11up: %f\n", m_f_r11_omega_add_time);
		printf("\t\t\t   ,\\: %f (ran " PRIsize " times)\n", f_total_lambda_up, m_n_lambda_update_num);
		printf("\t\t\t\tslice: %f\n", m_f_r11_lambda_slice_time);
		printf("\t\t\t\t Rata: %f\n", m_f_r11_lambda_ata_time);
		printf("\t\t\t\tR11up: %f\n", m_f_r11_lambda_add_time);
		printf("\t\t\t  Rup: %f // cholesky and fill\n", m_f_Rupdate_time);
		printf("\t    d: %f (resumed " PRIsize ", p-resumed " PRIsize ", full "
			PRIsize ")\n", m_f_d_time, m_n_resumed_forwardsubst_num,
			m_n_resumed_perm_forwardsubst_num, m_n_full_forwardsubst_num);
		printf("solver spent %f seconds in serial section\n", f_serial_time);
		printf("out of which:\n");
		printf("\t chol: %f (ran " PRIsize " times)\n", m_f_chol_time, m_n_lambda_optim_num);
		printf("\tbksub: %f (ran " PRIsize " times)\n", m_f_backsubst_time, m_n_R_optim_num);
		printf("\t norm: %f\n", m_f_norm_time);
		printf("\tv-upd: %f\n", m_f_vert_upd_time);
		/*printf("in unrelated news, small cholesky ran " PRIsize " times\n", m_n_dense_cholesky_num);
		printf("\t dense: %f\n", m_f_dense_cholesky_time);
		printf("\tsparse: %f\n", m_f_sparse_cholesky_time);*/ // dont want to do it runtime
#else // __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
		printf("it took: %f\n", f_total_time);
#endif // __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
	}

	/**
	 *	@brief calculates chi-squared error
	 *	@return Returns chi-squared error.
	 *	@note This only works with systems with edges of one degree of freedom
	 *		(won't work for e.g. systems with both poses and landmarks).
	 */
	inline double f_Chi_Squared_Error() /*const*/
	{
		if(m_r_system.r_Edge_Pool().b_Empty())
			return 0;
#ifndef __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
		if(m_lambda.n_BlockColumn_Num() != m_r_system.r_Vertex_Pool().n_Size()) {
			Optimize(0, 0);
			// optimize but don't allow iterations - just updates lambda, d and R
			// in order to be able to generate approximate solutions on request
		}
#endif // !__NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
		if(m_b_linearization_dirty) {
			if(!CalculateOneTimeDx())
				return -1;
			PushValuesInGraphSystem(m_v_dx);
			double f_chi2 = m_r_system.r_Edge_Pool().For_Each(CSum_ChiSquareError()) /
				(m_r_system.r_Edge_Pool().n_Size() - m_r_system.r_Edge_Pool()[0].n_Dimension());
			PushValuesInGraphSystem(-m_v_dx); // !!
			return f_chi2;
		} else {
			return m_r_system.r_Edge_Pool().For_Each(CSum_ChiSquareError()) /
				(m_r_system.r_Edge_Pool().n_Size() - m_r_system.r_Edge_Pool()[0].n_Dimension());
		}
	}

	/**
	 *	@brief calculates denormalized chi-squared error
	 *	@return Returns denormalized chi-squared error.
	 *	@note This doesn't perform the final division by (number of edges - degree of freedoms).
	 */
	inline double f_Chi_Squared_Error_Denorm() /*const*/
	{
		if(m_r_system.r_Edge_Pool().b_Empty())
			return 0;
#ifndef __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
		if(m_lambda.n_BlockColumn_Num() != m_r_system.r_Vertex_Pool().n_Size()) {
			Optimize(0, 0);
			// optimize but don't allow iterations - just updates lambda, d and R
			// in order to be able to generate approximate solutions on request
		}
#endif // !__NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
		if(m_b_linearization_dirty) {
			typedef CNonlinearSolver_FastL<CSystem, CLinearSolver, CAMatrixBlockSizes> TThisType;
			if(!CalculateOneTimeDx())
				return -1;
			PushValuesInGraphSystem(m_v_dx);
			double f_chi2 = m_r_system.r_Edge_Pool().For_Each(CSum_ChiSquareError());
			PushValuesInGraphSystem(-m_v_dx); // !!
			return f_chi2;
		} else
			return m_r_system.r_Edge_Pool().For_Each(CSum_ChiSquareError());
	}

	/**
	 *	@brief delays optimization upon Incremental_Step()
	 */
	inline void Delay_Optimization()
	{
		m_b_inhibit_optimization = true;
	}

	/**
	 *	@brief enables optimization upon Incremental_Step()
	 *
	 *	This is default behavior. In case it was disabled (by Delay_Optimization()),
	 *	and optimization is required, this will also run the optimization.
	 */
	inline void Enable_Optimization()
	{
		if(m_b_inhibit_optimization) {
			m_b_inhibit_optimization = false;

			TryOptimize(m_n_nonlinear_solve_max_iteration_num, m_f_nonlinear_solve_error_threshold);
			// optimize
		}
	}

	/**
	 *	@brief gets marginal covariances
	 *	@return Returns reference to the marginal covariances object.
	 */
	inline CMarginalCovariance &r_MarginalCovariance()
	{
		return m_marginals;
	}

	/**
	 *	@brief gets marginal covariances
	 *	@return Returns const reference to the marginal covariances object.
	 */
	inline const CMarginalCovariance &r_MarginalCovariance() const
	{
		return m_marginals;
	}

	/**
	 *	@brief incremental optimization function
	 *	@param[in] r_last_edge is the last edge that was added to the system
	 *	@note This function throws std::bad_alloc and std::runtime_error (when R is not-pos-def).
	 */
	void Incremental_Step(_TyBaseEdge &r_last_edge) // throw(std::bad_alloc, std::runtime_error)
	{
		size_t n_vertex_num = m_r_system.r_Vertex_Pool().n_Size();
		if(!m_b_had_loop_closure) {
			_ASSERTE(r_last_edge.n_Vertex_Id(0) != r_last_edge.n_Vertex_Id(1));
			size_t n_first_vertex = std::min(r_last_edge.n_Vertex_Id(0), r_last_edge.n_Vertex_Id(1));
			//size_t n_vertex_num = m_r_system.r_Vertex_Pool().n_Size();
			m_b_had_loop_closure = (n_first_vertex < n_vertex_num - 2);
			_ASSERTE(m_b_had_loop_closure || std::max(r_last_edge.n_Vertex_Id(0),
				r_last_edge.n_Vertex_Id(1)) == n_vertex_num - 1);
			/*if(m_b_had_loop_closure) {
				printf("" PRIsize ", " PRIsize " (out of " PRIsize " and " PRIsize ")\n", n_vertex_num, n_first_vertex,
					r_last_edge.n_Vertex_Id(0), r_last_edge.n_Vertex_Id(1));
			}*/ // debug
		}
		// detect loop closures (otherwise the edges are initialized based on measurement and error would be zero)

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
		{
			FILE *p_fw = fopen("timeSteps_R.txt", (m_n_real_step > 0)? "a" : "w");
			fprintf(p_fw, "" PRIsize ";%f;" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";"
				PRIsize ";" PRIsize "\n", m_n_real_step, m_shared_timer.f_Time(),
				m_n_Rup_num, m_n_full_R_num, m_n_R_optim_num, m_n_lambda_optim_num,
				std::max(r_last_edge.n_Vertex_Id(0), r_last_edge.n_Vertex_Id(1)) -
				std::min(r_last_edge.n_Vertex_Id(0), r_last_edge.n_Vertex_Id(1)),
				m_n_loop_size_cumsum);
			fclose(p_fw);
		}
		// dump time per step
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2
		bool b_new_vert = m_lambda.n_BlockColumn_Num() < m_r_system.r_Vertex_Pool().n_Size();

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE
		if(b_new_vert) {
			FILE *p_fw = fopen("chi2perVert.txt", (m_n_real_step > 0)? "a" : "w");
			double f_chi2 = 0;
			do {
				if(m_r_system.r_Edge_Pool().b_Empty())
					break;
#ifndef __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
				if(m_lambda.n_BlockColumn_Num() != m_r_system.r_Vertex_Pool().n_Size()) {
#pragma error "this might not work, before lambda is too small and after it is too big"
					Optimize(0, 0);
					// optimize but don't allow iterations - just updates lambda, d and R
					// in order to be able to generate approximate solutions on request
				}
#endif // !__NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
				if(m_b_linearization_dirty) {
					typedef CNonlinearSolver_FastL<CSystem, CLinearSolver, CAMatrixBlockSizes> TThisType;
					TThisType *p_this = const_cast<TThisType*>(this); // don't want to put 'mutable' arround everything
					if(!p_this->CalculateOneTimeDx(1)) // this is one smaller
						break;
					p_this->m_r_system.r_Vertex_Pool().For_Each_Parallel(0,
						m_r_system.r_Vertex_Pool().n_Size() - 1, CUpdateEstimates(m_v_dx)); // ignore the last vertex
					f_chi2 = m_r_system.r_Edge_Pool().For_Each(0,
						m_r_system.r_Edge_Pool().n_Size() - 1, CSum_ChiSquareError()); // ignore the last edge
					p_this->m_r_system.r_Vertex_Pool().For_Each_Parallel(0,
						m_r_system.r_Vertex_Pool().n_Size() - 1, CUpdateEstimates(-m_v_dx)); // !!
					break;
				} else {
					f_chi2 = m_r_system.r_Edge_Pool().For_Each(0,
						m_r_system.r_Edge_Pool().n_Size() - 1, CSum_ChiSquareError());
				}
			} while(0);
			// calculate chi2, excluding the last edge (the current one)

			fprintf(p_fw, "%f\n", f_chi2);
			fclose(p_fw);
		}
		// dump chi2
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2

		TryOptimize(m_n_nonlinear_solve_max_iteration_num, m_f_nonlinear_solve_error_threshold);
		// optimize

		if(m_t_marginals_config.b_calculate) {
			// todo - handle freq settings
			// todo - handle policies, now it defaults to full everytime

			_TyTimeSampler timer(m_shared_timer);

			CMarginals::Calculate_DenseMarginals_Fast(m_marginals.r_Matrix(), m_R);
			// calculate the thing

			timer.Accum_DiffSample(m_f_marginals_time);
		}
		// now R is up to date, can get marginals

#ifndef __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE
#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2
		if(b_new_vert) {
			FILE *p_fw = fopen("chi2perVert.txt", (m_n_real_step > 0)? "a" : "w");
			fprintf(p_fw, "%f\n", f_Chi_Squared_Error_Denorm());
			fclose(p_fw);
		}
		// dump chi2
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY
		if(b_new_vert)
			Dump_RDensity();
		// dump nnz of R, given different ordering strategies
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY

		++ m_n_real_step;
		// used only above
	}

protected:
	/**
	 *	@brief optimization function with optimization decision
	 *
	 *	@param[in] n_max_iteration_num is the maximal number of iterations
	 *	@param[in] f_min_dx_norm is the residual norm threshold
	 *
	 *	@note This function throws std::bad_alloc and std::runtime_error (when R is not-pos-def).
	 */
	void TryOptimize(size_t n_max_iteration_num, double f_min_dx_norm) // throw(std::bad_alloc, std::runtime_error)
	{
		bool b_optimization_triggered = false;

		size_t n_vertex_num = m_r_system.r_Vertex_Pool().n_Size();
		if(!m_b_inhibit_optimization) {
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
			size_t n_new_vertex_num = n_vertex_num - m_n_last_optimized_vertex_num;
			if(!n_new_vertex_num) // evidently, this backfires; need to have another m_n_last_optimized_vertex_num where it would remember when was the system last extended and check if there are new vertices since *then* // fixed now
				return; // no new vertices; don't go in ... (otherwise 2x slower on molson35, for obvious reasons)
			// the optimization periods are counted in vertices, not in edges (should save work on 10k, 100k)
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
			++ m_n_step;
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES

#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
			if(m_n_nonlinear_solve_threshold && n_new_vertex_num >= m_n_nonlinear_solve_threshold) {
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
			if(m_n_nonlinear_solve_threshold && m_n_step == m_n_nonlinear_solve_threshold) {
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
				m_n_last_optimized_vertex_num = n_vertex_num;
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
				m_n_step = 0;
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
				// do this first, in case Optimize() threw

				b_optimization_triggered = true;
				// nonlinear optimization
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
			} else if(m_n_linear_solve_threshold && n_new_vertex_num >= m_n_linear_solve_threshold) {
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
			} else if(m_n_linear_solve_threshold && m_n_step % m_n_linear_solve_threshold == 0) {
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
				_ASSERTE(!m_n_nonlinear_solve_threshold);
				m_n_last_optimized_vertex_num = n_vertex_num;
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
				if(!m_n_nonlinear_solve_threshold)
					m_n_step = 0;
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
				// do this first, in case Optimize() threw

				n_max_iteration_num = 1;
				f_min_dx_norm = m_f_nonlinear_solve_error_threshold; // right?
				b_optimization_triggered = true;
				// simple optimization
			}
		}

		if(!b_optimization_triggered) {
#ifdef __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
			size_t n_vertex_num = m_r_system.r_Vertex_Pool().n_Size();
			if(m_lambda.n_BlockColumn_Num() == n_vertex_num)
				return;
			// there is enough vertices in lambda, none would be added

			Optimize(0, 0); // big todo - remove this in order to be faster for each 100; move it to function that does approx solutions on request
			// optimize but don't allow iterations - just updates lambda, d and R
			// in order to be able to generate approximate solutions on request
#endif // __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
		} else
			Optimize(n_max_iteration_num, f_min_dx_norm);
	}

	/**
	 *	@brief refreshes system matrices lambda and R
	 *	@return Returns true if optimization should take place, otherwise returns false
	 */
	bool RefreshLambdaR()
	{
		_TyTimeSampler timer(m_shared_timer);

		const size_t n_variables_size = m_r_system.n_VertexElement_Num();
		const size_t n_measurements_size = m_r_system.n_EdgeElement_Num();
		if(n_variables_size > n_measurements_size) {
			fprintf(stderr, "warning: the system is underspecified\n");
			return false;
		}
		if(!n_measurements_size)
			return false; // nothing to solve (but no results need to be generated so it's ok)
		// can't solve in such conditions

		// note that n_order_min can be possibly used in Refresh_Lambda()
		// to minimize number of vertices that require update of hessian blocks
		// note that it needs to be the one with permutation? or does it? // todo - try that after it works

		Extend_LambdaR(m_n_verts_in_lambda, m_n_edges_in_lambda);
		// recalculated all the jacobians inside Extend_LambdaR(), also extend R structurally

		if(!m_b_system_dirty)
			Refresh_Lambda(0/*m_n_verts_in_lambda*/, m_n_edges_in_lambda); // calculate only for new edges // @todo - but how to mark affected vertices?
		else
			Refresh_Lambda(); // calculate for entire system, rebuild R from scratch
		// refresh lambda (does not fully refresh permutated lambda, even though it is a reference matrix)

		m_v_dx.resize(n_variables_size);
		m_v_perm_temp.resize(n_variables_size);
		if(m_b_R_up_to_date && !m_b_system_dirty)
			m_v_d.conservativeResize(n_variables_size); // b_Refresh_R() also refreshes d (rhs), needs dx as temp // !!
		else
			m_v_d.resize(n_variables_size); // in case we're about to rebuild R from scratch, don't care about contents of d
		// resize the helper vectors

		timer.Accum_DiffSample(m_f_lambda_refresh_time);

		{
			if(m_b_R_up_to_date && // can increment R only if up to date
			   !m_b_system_dirty) // avoidance of big incremental updates of R is inside b_Refresh_R() - can only decide if ordering is known
				m_b_first_iteration_use_R = b_Refresh_R(0/*m_n_verts_in_lambda*/, m_n_edges_in_lambda); // calculate only for new edges // @todo - but how to mark affected vertices?
			else
				m_b_first_iteration_use_R = b_Refresh_R(); // calculate for entire system, rebuild R from scratch

			m_b_R_up_to_date = m_b_first_iteration_use_R;
			// in case R is not used, it will fall behind
		}
		// refresh R

		m_b_system_dirty = false;
		m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
		m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
		_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
			m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
			m_lambda.n_Column_Num() == n_variables_size);
		_ASSERTE(!m_b_R_up_to_date || (m_R.n_Row_Num() == m_R.n_Column_Num() &&
			m_R.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
			m_R.n_Column_Num() == n_variables_size)); // lambda is square, blocks on either side = number of vertices
		// need to have lambda and perhaps also R

		return true;
	}

	/**
	 *	@brief updates the m_v_dx vector from the current R and d (or lambda eta)
	 *	@param[in] n_ignore_vertices is number of vertices at the end of the system to be ignored
	 *	@return Returns true on success, false on failure (numerical issues).
	 */
	bool CalculateOneTimeDx(size_t n_ignore_vertices = 0)
	{
		_ASSERTE(m_b_linearization_dirty); // this should only be called in case the linearization point was not updated

		if(m_b_R_up_to_date && m_b_first_iteration_use_R) { // Optimize() clears m_b_R_up_to_date but not m_b_first_iteration_use_R at the same time
			_ASSERTE(m_b_R_up_to_date);
			// we have R and can use it efficiently

			{
				bool b_cholesky_result;
				{
					_ASSERTE(m_p_lambda_block_ordering);
					m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
						m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
					b_cholesky_result = m_R.UpperTriangular_Solve_FBS<_TyLambdaMatrixBlockSizes>(&m_v_perm_temp(0), m_v_perm_temp.rows());
					m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_dx(0), &m_v_perm_temp(0), m_v_dx.rows(),
						m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
					// dx = R'/d // note this never fails (except if R is null)
				}
				// R solves with permutation (note that m_v_d is not modified!)
				// calculate cholesky, reuse block ordering if the linear solver supports it

				if(!b_cholesky_result)
					return false;

#ifdef _DEBUG
				for(size_t i = 0, n_variables_size = m_v_dx.rows(); i < n_variables_size; ++ i) {
					if(_isnan(m_v_dx(i)))
						fprintf(stderr, "warning: p_dx[" PRIsize "] = NaN (file \'%s\', line " PRIsize ")\n", i, __FILE__, __LINE__);
				}
				// detect the NaNs, if any (warn, but don't modify)
#endif // _DEBUG
			}
		} else {
			if(!n_ignore_vertices)
				Collect_RightHandSide_Vector(m_v_dx);
			else {
				m_r_system.r_Vertex_Pool().For_Each_Parallel(0,
					m_r_system.r_Vertex_Pool().n_Size() - n_ignore_vertices,
					CCollect_RightHandSide_Vector(m_v_dx));
			}
			// collects the right-hand side vector

			{
				bool b_cholesky_result;
				{
					Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
					b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda, v_eta); // p_dx = eta = lambda / eta
					// dont reuse block ordering
				}
				// lambda is good without permutation (there is one inside and we save copying eta arround)
				// calculate cholesky, reuse block ordering if the linear solver supports it

				if(!b_cholesky_result)
					return false;

#ifdef _DEBUG
				for(size_t i = 0, n_variables_size = m_v_dx.rows(); i < n_variables_size; ++ i) {
					if(_isnan(m_v_dx(i)))
						fprintf(stderr, "warning: p_dx[" PRIsize "] = NaN (file \'%s\', line " PRIsize ")\n", i, __FILE__, __LINE__);
				}
				// detect the NaNs, if any (warn, but don't modify)
#endif // _DEBUG
			}
		}
		// just solve and check NaNs in debug, nothing more
		// can't update timers as some part of pipeline is not run, don't update counters neither

		return true;
		// the result is in m_v_dx
	}

public:
	/**
	 *	@brief final optimization function
	 *
	 *	@param[in] n_max_iteration_num is the maximal number of iterations
	 *	@param[in] f_min_dx_norm is the residual norm threshold
	 *
	 *	@note This function throws std::bad_alloc and std::runtime_error (when R is not-pos-def).
	 */
	void Optimize(size_t n_max_iteration_num = 5, double f_min_dx_norm = .01) // throw(std::bad_alloc, std::runtime_error)
	{
		if(!RefreshLambdaR())
			return;
		// decide whether to optimize or not

		if(!n_max_iteration_num) {
			m_b_linearization_dirty = true;
			return;
		}
		// in case we're not required to optimize, do nothing
		// (the user can still request solution, R is in good shape)

		if(m_b_had_loop_closure)
			m_b_had_loop_closure = false;
		else {
			m_b_linearization_dirty = true;
			return; // nothing to optimize, dx would be zero
		}
		// handle loop closures a bit differently

#if 0
		static bool b_had_lambda_up = false;
		if(m_b_first_iteration_use_R && b_had_lambda_up) {
			b_had_lambda_up = false;
			Check_RLambdaTracking(); // seems to be working now
		}
#endif // 0
		// make sure lambda and R contain the same system

		bool b_verbose = m_b_verbose;

		for(size_t n_iteration = 0; n_iteration < n_max_iteration_num; ++ n_iteration) {
			++ m_n_iteration_num;
			// debug

			if(b_verbose) {
				if(n_max_iteration_num == 1)
					printf("\n=== incremental optimization step ===\n\n");
				else
					printf("\n=== nonlinear optimization: iter #" PRIsize " ===\n\n", n_iteration);
			}
			b_verbose = m_b_verbose; // restore
			// verbose

			if(m_b_R_up_to_date/*m_b_first_iteration_use_R && !n_iteration*/) { // always when we have R
				++ m_n_R_optim_num;

				_ASSERTE(m_b_R_up_to_date);
				// we have R and can use it efficiently

				_TyTimeSampler timer(m_shared_timer);

				{
					bool b_utsolve_result;
					{
						_ASSERTE(m_p_lambda_block_ordering);
						m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
							m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
						b_utsolve_result = m_R.UpperTriangular_Solve_FBS<_TyLambdaMatrixBlockSizes>(&m_v_perm_temp(0), m_v_perm_temp.rows());
						m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_dx(0), &m_v_perm_temp(0), m_v_dx.rows(),
							m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
						// dx = R'/d // note this never fails (except if R is null)

						if(m_b_verbose) {
							printf("%s", (b_utsolve_result)? "backsubstitution succeeded\n" :
								"backsubstitution failed\n");
						}
					}
					// R solves with permutation (note that m_v_d is not modified!)

					// calculate cholesky, reuse block ordering if the linear solver supports it

					timer.Accum_DiffSample(m_f_backsubst_time);

#ifdef _DEBUG
					for(size_t i = 0, n_variables_size = m_r_system.n_VertexElement_Num();
					   i < n_variables_size; ++ i) {
						if(_isnan(m_v_dx(i)))
							fprintf(stderr, "warning: p_dx[" PRIsize "] = NaN (file \'%s\', line " PRIsize ")\n", i, __FILE__, __LINE__);
					}
					// detect the NaNs, if any (warn, but don't modify)
#endif // _DEBUG

					double f_residual_norm = 0;
					if(b_utsolve_result) {
						f_residual_norm = m_v_dx.norm(); // Eigen likely uses SSE and OpenMP
						if(m_b_verbose)
							printf("residual norm: %.4f\n", f_residual_norm);
					}
					// calculate residual norm

					timer.Accum_DiffSample(m_f_norm_time);

					if(f_residual_norm <= f_min_dx_norm) {
						m_b_linearization_dirty = true;
						break;
					}
					// in case the error is low enough, quit (saves us recalculating the hessians)

					if(b_utsolve_result) {
						/*printf("just optimized using R\n");*/
						PushValuesInGraphSystem(m_v_dx); // note this kills R
						m_b_system_dirty = true;
						m_b_R_up_to_date = false; // !!

						m_b_linearization_dirty = false;

						timer.Accum_DiffSample(m_f_vert_upd_time);
					}
					// update the system (in parallel)

					if(!b_utsolve_result)
						break;
					// in case cholesky failed, quit
				}
			} else {
				_TyTimeSampler timer(m_shared_timer);

				if(n_iteration && m_b_system_dirty) {
					Refresh_Lambda(); // want only lambda, leave R behind
					m_b_system_dirty = false;
					m_b_R_up_to_date = false; // lambda not dirty anymore, but R still is

					timer.Accum_DiffSample(m_f_lambda_refresh_time);

#ifdef __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
					m_b_R_up_to_date = b_Refresh_R(0, 0, n_iteration > 1); // refresh R as well
					// suppress reordering in iterations 2 and above
					// (already reordered in iteration 1, won't get any better)
					// note that RHS vector is updated inside

					_TyTime f_dummy_sample = 0;
					timer.Accum_DiffSample(f_dummy_sample); // b_Refresh_R() contains timing inside

					if(m_b_R_up_to_date) {
						-- n_iteration;
						-- m_n_iteration_num;
						b_verbose = false; // suppress the banner
						continue;
					}
					// try again, this time with R
#endif // __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
				}
				// no need to rebuild lambda, just refresh the values that are being referenced

				++ m_n_lambda_optim_num;
				// we fall back to lambda

				Collect_RightHandSide_Vector(m_v_dx);
				// collects the right-hand side vector

				timer.Accum_DiffSample(m_f_rhs_time);

				{
					bool b_cholesky_result;
					{
						Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
						if((/*m_b_first_iteration_use_R &&*/ n_max_iteration_num > 2) ||
						   (!m_b_first_iteration_use_R && n_max_iteration_num > 1)) {
							do {
								if(n_iteration == ((m_b_first_iteration_use_R)? 1 : 0) &&
								   !_TyLinearSolverWrapper::FinalBlockStructure(m_linear_solver, m_lambda)) {
									b_cholesky_result = false;
									break;
								}
								// prepare symbolic decomposition, structure of lambda won't change in the next steps

								b_cholesky_result = _TyLinearSolverWrapper::Solve(m_linear_solver, m_lambda, v_eta);
								// p_dx = eta = lambda / eta
							} while(0);
						} else
							b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda, v_eta); // p_dx = eta = lambda / eta

						if(m_b_verbose)
							printf("%s", (b_cholesky_result)? "Cholesky succeeded\n" : "Cholesky failed\n");
					}
					// lambda is good without permutation (there is one inside and we save copying eta arround)
					// calculate cholesky, reuse block ordering if the linear solver supports it

					timer.Accum_DiffSample(m_f_chol_time);

#ifdef _DEBUG
					for(size_t i = 0, n_variables_size = m_r_system.n_VertexElement_Num();
					   i < n_variables_size; ++ i) {
						if(_isnan(m_v_dx(i)))
							fprintf(stderr, "warning: p_dx[" PRIsize "] = NaN (file \'%s\', line " PRIsize ")\n", i, __FILE__, __LINE__);
					}
					// detect the NaNs, if any (warn, but don't modify)
#endif // _DEBUG

					double f_residual_norm = 0;
					if(b_cholesky_result) {
						f_residual_norm = m_v_dx.norm(); // Eigen likely uses SSE and OpenMP
						if(m_b_verbose)
							printf("residual norm: %.4f\n", f_residual_norm);
					}
					// calculate residual norm

					timer.Accum_DiffSample(m_f_norm_time);
					// timing breakup

					if(f_residual_norm <= f_min_dx_norm) {
						m_b_linearization_dirty = true;
						break;
					}
					// in case the error is low enough, quit (saves us recalculating the hessians)

					if(b_cholesky_result) {
						/*printf("just optimized using lambda\n");*/

						PushValuesInGraphSystem(m_v_dx);

						timer.Accum_DiffSample(m_f_vert_upd_time);

						m_b_system_dirty = true;
						m_b_R_up_to_date = false;

						m_b_linearization_dirty = false;
					}
					// update the system (in parallel)

					if(!b_cholesky_result)
						break;
					// in case cholesky failed, quit
				}
			}
		}

#ifdef __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
		_ASSERTE(m_b_system_dirty || m_b_R_up_to_date);
		// make sure that R is indeed kept up-to-date, unless the solver
		// was stopped by reaching the maximum number of iterations
#endif // __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
	}

protected:
	/**
	 *	@brief function object that calls lambda hessian block allocation for all edges
	 */
	class CAlloc_LambdaRBlocks { // t_odo - R probably only allocates blocks on vertices; retain old version of this functor with lambda only for edges
	protected:
		CUberBlockMatrix &m_r_lambda; /**< @brief reference to the lambda matrix (out) */
		CUberBlockMatrix &m_r_R; /**< @brief reference to the R matrix (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_lambda is reference to the lambda matrix
		 *	@param[in] r_R is reference to the lambda matrix
		 */
		inline CAlloc_LambdaRBlocks(CUberBlockMatrix &r_lambda, CUberBlockMatrix &r_R)
			:m_r_lambda(r_lambda), m_r_R(r_R)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyVertexOrEdge is vertex or edge type
		 *	@param[in,out] r_vertex_or_edge is vertex or edge to have hessian blocks allocated in lambda
		 *	@note This function throws std::bad_alloc.
		 */
		template <class _TyVertexOrEdge>
		inline void operator ()(_TyVertexOrEdge &r_vertex_or_edge) // throw(std::bad_alloc)
		{
			r_vertex_or_edge.Alloc_HessianBlocks(m_r_lambda);
			r_vertex_or_edge.Alloc_LBlocks(m_r_R); // t_odo - alloc R blocks as well
		}
	};

	/**
	 *	@brief function object that calls R factor block allocation for all vertices
	 */
	class CAlloc_RBlocks {
	protected:
		CUberBlockMatrix &m_r_R; /**< @brief reference to the R matrix (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_R is reference to the lambda matrix
		 */
		inline CAlloc_RBlocks(CUberBlockMatrix &r_R)
			:m_r_R(r_R)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyVertex is vertex type
		 *	@param[in,out] r_vertex is vertex to have hessian blocks allocated in R
		 *	@note This function throws std::bad_alloc.
		 */
		template <class _TyVertex>
		inline void operator ()(_TyVertex &r_vertex) // throw(std::bad_alloc)
		{
			r_vertex.Alloc_LBlocks(m_r_R); // t_odo - alloc R blocks as well
		}
	};

	/**
	 *	@brief function object that calls lambda hessian block allocation for all edges
	 */
	class CAlloc_LambdaBlocks {
	protected:
		CUberBlockMatrix &m_r_lambda; /**< @brief reference to the lambda matrix (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_lambda is reference to the lambda matrix
		 */
		inline CAlloc_LambdaBlocks(CUberBlockMatrix &r_lambda)
			:m_r_lambda(r_lambda)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyVertexOrEdge is vertex or edge type
		 *	@param[in,out] r_vertex_or_edge is vertex or edge to have hessian blocks allocated in lambda
		 *	@note This function throws std::bad_alloc.
		 */
		template <class _TyVertexOrEdge>
		inline void operator ()(_TyVertexOrEdge &r_vertex_or_edge) // throw(std::bad_alloc)
		{
			r_vertex_or_edge.Alloc_HessianBlocks(m_r_lambda);
		}
	};

	/**
	 *	@brief function object that calls omega hessian block allocation and evaluation for all edges
	 */
	class CCalculateOmega {
	protected:
		CUberBlockMatrix &m_r_omega; /**< @brief reference to the omega matrix (out) */
		size_t m_n_min_elem_order; /**< @brief minimal order of vertex to be filled in omega (in elements) */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in] r_omega is reference to the omega matrix
		 *	@param[in] n_min_elem_order is minimal order of vertex to be filled in omega (in elements)
		 */
		inline CCalculateOmega(CUberBlockMatrix &r_omega, size_t n_min_elem_order)
			:m_r_omega(r_omega), m_n_min_elem_order(n_min_elem_order)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is an edge type
		 *	@param[in] r_edge is edge to have its hessian blocks added to omega
		 *	@note This function throws std::bad_alloc.
		 */
		template <class _TyVertexOrEdge>
		inline void operator ()(const _TyVertexOrEdge &r_edge) // throw(std::bad_alloc)
		{
			r_edge.Calculate_Omega(m_r_omega, m_n_min_elem_order);
		}
	};

	/**
	 *	@brief function object that calls b vector calculation for all edges
	 */
	class CCollect_RightHandSide_Vector {
	protected:
		Eigen::VectorXd &m_r_b; /**< @brief reference to the right-hand side vector (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_b is reference to the right-hand side vector
		 */
		inline CCollect_RightHandSide_Vector(Eigen::VectorXd &r_b)
			:m_r_b(r_b)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in,out] r_t_edge is edge to output its part R error vector
		 */
		template <class _TyEdge>
		inline void operator ()(_TyEdge &r_t_edge) // throw(std::bad_alloc)
		{
			r_t_edge.Get_RightHandSide_Vector(m_r_b);
		}
	};

	/**
	 *	@brief calculates the right-hand side vector
	 *	@param[out] r_v_b is the right-hand side vector (allocated by the caller)
	 */
	inline void Collect_RightHandSide_Vector(Eigen::VectorXd &r_v_b)
	{
		m_r_system.r_Vertex_Pool().For_Each_Parallel(CCollect_RightHandSide_Vector(r_v_b)); // can do this in parallel
		// collect b
	}

	/**
	 *	@brief creates the lambda matrix from scratch
	 */
	inline void AddEntriesInSparseSystem() // throw(std::bad_alloc)
	{
		// note: don't worry about very large matrices being built at once,
		// this will most likely only be any good for incremental

		m_lambda.Clear();
		const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
		if(!m_lambda.Append_Block(Eigen::MatrixXd(r_t_uf.transpose() * r_t_uf), 0, 0))
			throw std::bad_alloc();
		// add unary factor (actually UF^T * UF, but it's the same matrix)

		m_R.Clear();
		m_r_system.r_Edge_Pool().For_Each(CAlloc_LambdaBlocks(m_lambda));
		m_r_system.r_Vertex_Pool().For_Each(CAlloc_LambdaRBlocks(m_lambda, m_R)); // can stay, there is no ordering to be applied
		// add all the hessian blocks
	}

	/**
	 *	@brief incrementally updates the lambda matrix structure (must not be empty)
	 *
	 *	@param[in] n_skip_vertices is number of vertices before the first vertex that changes
	 *	@param[in] n_skip_edges is number of edges before the first edge that changes
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	inline void UpdateSparseSystem(size_t n_skip_vertices, size_t n_skip_edges) // throw(std::bad_alloc)
	{
		//throw std::runtime_error("UpdateSparseSystem not implemented for R"); // t_odo

		_ASSERTE(m_lambda.n_Row_Num() > 0 && m_lambda.n_Column_Num() == m_lambda.n_Row_Num()); // make sure lambda is not empty
		//m_r_system.r_Edge_Pool().For_Each(n_skip_edges,
		//	m_r_system.r_Edge_Pool().n_Size(), CAlloc_LambdaRBlocks(m_lambda, m_R)); // no, edges do not add R blocks
		m_r_system.r_Edge_Pool().For_Each(n_skip_edges,
			m_r_system.r_Edge_Pool().n_Size(), CAlloc_LambdaBlocks(m_lambda));
		m_r_system.r_Vertex_Pool().For_Each(n_skip_vertices,
			m_r_system.r_Vertex_Pool().n_Size(), CAlloc_LambdaRBlocks(m_lambda, m_R)); // will not work if ordering is applied (but it mostly isn't, the increments follow identity ordering)
		// add the hessian blocks of the new edges
	}

	/**
	 *	@brief incrementally updates the lambda matrix structure (can be empty)
	 *
	 *	@param[in] n_vertices_already_in_lambda is number of vertices before the first vertex that changes
	 *	@param[in] n_edges_already_in_lambda is number of edges before the first edge that changes
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	inline void Extend_LambdaR(size_t n_vertices_already_in_lambda, size_t n_edges_already_in_lambda) // throw(std::bad_alloc)
	{
		if(!n_vertices_already_in_lambda && !n_edges_already_in_lambda)
			AddEntriesInSparseSystem(); // works for empty
		else
			UpdateSparseSystem(n_vertices_already_in_lambda, n_edges_already_in_lambda); // does not work for empty
		// create block matrix lambda
	}

#if 0
	/**
	 *	@brief checks if R == chol(lambda), prints the norm of the difference to stdout
	 */
	void Check_RLambdaTracking() const
	{
		CUberBlockMatrix RtR_upper;
		m_R.PreMultiplyWithSelfTransposeTo(RtR_upper, true);
		//cs *p_R = m_R.p_Convert_to_Sparse();
		cs *p_lam = m_lambda.p_Convert_to_Sparse();
		//cs *p_Rt = cs_transpose(p_R, 1);
		cs *p_RtR = RtR_upper.p_Convert_to_Sparse();//cs_multiply(p_Rt, p_R);
		cs *p_diff = cs_add(p_RtR, p_lam, 1, -1);
		double f_norm = cs_norm(p_diff);
		//cs_spfree(p_R);
		cs_spfree(p_lam);
		//cs_spfree(p_Rt);
		cs_spfree(p_RtR);
		cs_spfree(p_diff);
		// calculate norm (R*R' - lambda)

		printf("R - lambda tracking: %f\n", f_norm);
	}
#endif // 0

	/**
	 *	@brief calculates the new \f$R_{11}\f$ matrix
	 *
	 *	@param[in] n_order_min is the minimum vertex that changes in R (zero-based index in blocks)
	 *	@param[in] n_order_max is number of column blocks in R (in blocks)
	 *	@param[in] r_R11_new is matrix, containing the new \f$R_{11}\f$, before calculating cholesky of it
	 *
	 *	@note This may modify / damage r_R11_new as it is no longer needed and it is *not* a reference
	 *		to a part of R, in case that enables speed optimizations.
	 *	@note This function throws std::bad_alloc and std::runtime_error (when not-pos-def).
	 */
	void Refresh_R11(size_t n_order_min, size_t n_order_max, CUberBlockMatrix &r_R11_new) // throw(std::bad_alloc, std::runtime_error)
	{
#ifdef __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
		_ASSERTE(r_R11_new.n_Row_Num() == r_R11_new.n_Column_Num());
		if(r_R11_new.n_Column_Num() < 150) {
			if(!r_R11_new.Cholesky_Dense_FBS<_TyLambdaMatrixBlockSizes, 15>()) // 15, not 150 (would yield >1024 template depth)
				throw std::runtime_error("Cholesky_Dense_FBS() failed to increment R");
			// use statically sized matrices up to 30x30, then dynamically allocated matrices up to 150x150

			m_R.From_Matrix(n_order_min, n_order_min, r_R11_new);
			// put r_R11_new to m_R
		} else {
#else // __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
		{
#endif // __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
			if(!b_Have_NativeSolver) {
				CUberBlockMatrix R11_old;
				m_R.SliceTo(R11_old, n_order_min, n_order_max, n_order_min, n_order_max, true); // get R11 as well, need to clear the blocks first
				// todo - make a ClearBlocks() function to a) delete or b) memset(0) blocks in a rectangular area

				R11_old.Scale_FBS_Parallel<_TyLambdaMatrixBlockSizes>(0);
				// clears the data in the update area (it would be better to erase the blocks, but there is the fap) // todo - see indices of the blocks in R and see if these could be efficiently erased right now (while keeping the structure)
				// note that even if we thrash memory taken by (some) R11 blocks,
				// it will be recollected once a full update takes place
			}
			// only have to clear for the old solvers, the native solver does it automatically

			if(!m_linear_solver2.Factorize_PosDef_Blocky(m_R, r_R11_new,
			   m_R_row_lookup_table, n_order_min, n_order_min))
				throw std::runtime_error("Factorize_PosDef_Blocky() failed to increment R");
		}
	}

	/**
	 *	@brief calculates the new R matrix incrementally using lambda or omega
	 *
	 *	@param[in] n_refresh_from_edge is the first edge that changes
	 *	@param[in] n_order_min is the minimum vertex that changes in R (zero-based index in blocks)
	 *
	 *	@note This function throws std::bad_alloc and std::runtime_error (when not-pos-def).
	 */
	inline void Refresh_R_IncR11(size_t n_refresh_from_edge, size_t n_order_min) // throw(std::bad_alloc, std::runtime_error)
	{
		_ASSERTE(m_b_R_up_to_date);
		// make sure R is up to date with lambda and we can actually increment

		_ASSERTE(n_refresh_from_edge > 0);
		// make sure that lambda didn't rebuild the ordering completely

		const size_t n_order_max = m_lambda.n_BlockColumn_Num();
		// makes sure that max is fixed at the end

		size_t n_min_vertex = m_p_lambda_block_ordering[n_order_min];
		bool b_identity_perm = true;
		for(size_t i = n_order_min; i < n_order_max; ++ i) {
			if(m_p_lambda_block_ordering[i] != i) {
				b_identity_perm = false;
				break;
			}
		}
		if(!b_identity_perm) {
			for(size_t i = n_order_min + 1; i < n_order_max; ++ i)
				n_min_vertex = std::min(n_min_vertex, m_p_lambda_block_ordering[i]);
		}
		// see if the ordering is identity ordering

#ifndef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		if(n_order_max - n_order_min >= 88) // disable this for timing bench
			b_identity_perm = false;
#endif // !__NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		// yet another threshold

		_TyTimeSampler timer(m_shared_timer);

		bool b_omega_available = b_identity_perm;

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		b_identity_perm = true; // i want both branches to run
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

		CUberBlockMatrix R11TR11;
		if(b_identity_perm) {
			++ m_n_omega_update_num;
			// count them

			CUberBlockMatrix omega, R11;
			size_t n_elem_order_min = m_lambda_perm.n_BlockColumn_Base(n_order_min);
			m_r_system.r_Edge_Pool().For_Each(n_refresh_from_edge, m_r_system.r_Edge_Pool().n_Size(),
				CCalculateOmega(omega, n_elem_order_min));
			omega.CheckIntegrity();

			timer.Accum_DiffSample(m_f_r11_omega_calc_time);

			m_R.SliceTo(R11, n_order_min, n_order_max, n_order_min, n_order_max, true); // row(0 - min) x col(min - max)
			// calculate the omega matrix (ho, ho, ho) and slice R11

			timer.Accum_DiffSample(m_f_r11_omega_slice_time);

			if(n_order_max - n_order_min >= __NONLINEAR_SOLVER_FAST_L_PARALLEL_MATMULT_THRESH) // big one // t_odo this never runs, the limit for using R is also 100
				R11.PreMultiplyWithSelfTransposeTo_FBS_Parallel<_TyLambdaMatrixBlockSizes>(R11TR11, true); // calculate R11^T * R11 // t_odo - use FBS and maybe also parallel // t_odo - only need lower diagonal, what to do?
			else
				R11.PreMultiplyWithSelfTransposeTo_FBS<_TyLambdaMatrixBlockSizes>(R11TR11, true); // calculate R11^T * R11 // t_odo - use FBS and maybe also parallel // t_odo - only need lower diagonal, what to do?
			// calculate R11TR11

			timer.Accum_DiffSample(m_f_r11_omega_ata_time);

			bool UNUSED(b_result) = omega.AddTo_FBS<_TyLambdaMatrixBlockSizes>(R11TR11); // todo - maybe also parallel
			_ASSERTE(b_result); // if the block order in omega was wrong, this would fail
			// calculate R11TR11_new = R11TR11 + omega
			// note this uses faster addition algorithm

			timer.Accum_DiffSample(m_f_r11_omega_add_time);
		}

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		b_identity_perm = false; // i want both branches to run
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

		CUberBlockMatrix R01TR01;
		if(!b_identity_perm) {
			CUberBlockMatrix lambda11, R01;
			++ m_n_lambda_update_num;
			// count them

			m_R.SliceTo(R01, 0, n_order_min, n_order_min, n_order_max, true); // row(0 - min) x col(min - max)
			m_lambda_perm.SliceTo(lambda11, n_order_min, n_order_max, n_order_min, n_order_max, true);

			timer.Accum_DiffSample(m_f_r11_lambda_slice_time);

			if(n_order_max - n_order_min >= __NONLINEAR_SOLVER_FAST_L_PARALLEL_MATMULT_THRESH) // big one // t_odo this never runs, the limit for using R is also 100
				R01.PreMultiplyWithSelfTransposeTo_FBS_Parallel<_TyLambdaMatrixBlockSizes>(R01TR01, true); // calculate R01^T * R01 // t_odo - use FBS and maybe also parallel // t_odo - only need lower diagonal, what to do?
			else
				R01.PreMultiplyWithSelfTransposeTo_FBS<_TyLambdaMatrixBlockSizes>(R01TR01, true); // calculate R01^T * R01 // t_odo - use FBS and maybe also parallel // t_odo - only need lower diagonal, what to do?

			timer.Accum_DiffSample(m_f_r11_lambda_ata_time);

			lambda11.AddTo_FBS<_TyLambdaMatrixBlockSizes>(R01TR01, -1, 1); // t_odo - use FBS // todo - maybe also parallel
			// calculates R01TR01 = -R01TR01 + lambda11 (note the "-1, 1" is correct, the opposite way it crashes)
			// note that lambda11 is upper diagonal, as well as R01TR01

			timer.Accum_DiffSample(m_f_r11_lambda_add_time);
		}

		Refresh_R11(n_order_min, n_order_max, (b_omega_available)? R11TR11 : R01TR01);

		timer.Accum_DiffSample(m_f_Rupdate_time);

		Refresh_d_IncR11(n_refresh_from_edge, n_order_min); // use the function, do not repeat code, it is ...
		// note that this contains its own timing inside
	}

	/**
	 *	@brief calculates the new right-hand-side vector, does it incrementally were possible
	 *
	 *	@param[in] n_refresh_from_edge is the first edge that changes
	 *	@param[in] n_order_min is the minimum vertex that changes in R (zero-based index in blocks)
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	inline void Refresh_d_IncR11(size_t UNUSED(n_refresh_from_edge), size_t n_order_min) // throw(std::bad_alloc)
	{
		_ASSERTE(m_b_R_up_to_date);
		// make sure R is up to date with lambda and we can actually increment

		_ASSERTE(n_refresh_from_edge > 0);
		// make sure that lambda didn't rebuild the ordering completely

		_TyTimeSampler timer(m_shared_timer);

		const size_t n_order_max = m_lambda.n_BlockColumn_Num();
		// makes sure that max is fixed at the end

		size_t n_min_vertex = m_p_lambda_block_ordering[n_order_min];
		bool b_identity_perm = true;
		for(size_t i = n_order_min; i < n_order_max; ++ i) {
			if(m_p_lambda_block_ordering[i] != i) {
				b_identity_perm = false;
				break;
			}
		}
		if(!b_identity_perm) {
			for(size_t i = n_order_min + 1; i < n_order_max; ++ i)
				n_min_vertex = std::min(n_min_vertex, m_p_lambda_block_ordering[i]);
		}
		const bool b_is_identity_perm = b_identity_perm; // make a copy

		{
			if(b_is_identity_perm) {
				++ m_n_resumed_forwardsubst_num;

				_ASSERTE(m_v_d.rows() == m_lambda.n_Column_Num());
				m_r_system.r_Vertex_Pool().For_Each_Parallel(n_order_min,
					m_r_system.r_Vertex_Pool().n_Size(), CCollect_RightHandSide_Vector(m_v_d));
				// collect part of b to the lower part of d (this is inside of Collect_RightHandSide_Vector())

				timer.Accum_DiffSample(m_f_rhs_time);

				m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				m_R.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(
					&m_v_perm_temp(0), m_v_perm_temp.rows(), n_order_min);
				m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				// "resumed forward substitution"
			} else {
				++ m_n_resumed_perm_forwardsubst_num; // a different category

				_ASSERTE(m_v_d.rows() == m_lambda.n_Column_Num());
				{
					if(!n_min_vertex) {
						-- m_n_resumed_perm_forwardsubst_num;
						++ m_n_full_forwardsubst_num; // this is really a full one

						m_r_system.r_Vertex_Pool().For_Each_Parallel(n_min_vertex,
							m_r_system.r_Vertex_Pool().n_Size(), CCollect_RightHandSide_Vector(m_v_d));
					} else {
						const size_t *p_order_inv = m_lambda_ordering.p_Get_Ordering();
						CCollect_RightHandSide_Vector collector(m_v_d);
#ifdef _OPENMP
						_ASSERTE(n_order_max <= INT_MAX);
						const int n = int(n_order_max);
						#pragma omp parallel for default(shared) if(n - int(n_min_vertex) >= 50)
						for(int i = int(n_min_vertex); i < n; ++ i)
#else // _OPENMP
						for(size_t i = n_min_vertex; i < n_order_max; ++ i)
#endif // _OPENMP
							collector(m_r_system.r_Vertex_Pool()[p_order_inv[i]]);
						// can do this in parallel as well
					}
				}
				// do this instead

				timer.Accum_DiffSample(m_f_rhs_time);

				// so ... updating vertices n_order_min to n_order_max (a contiguous range?)
				// after permutation, they go to vector[m_p_lambda_block_ordering[n_order_min to n_order_max]]
				// so the whole part of the _permuted_ vector from m_p_lambda_block_ordering[n_order_min] to
				// the end must be updated, and *no more*

				m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num()); // dest[p ++] = *src ++
				m_R.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(
					&m_v_perm_temp(0), m_v_perm_temp.rows(), n_min_vertex);
				m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num()); // *dest ++ = src[p ++]
				// "resumed forward substitution"
			}
			// convert eta to d (d = eta/R)
		}
		// update d incrementally as well

		timer.Accum_DiffSample(m_f_d_time);
	}

	/**
	 *	@brief calculates the new R matrix from scratch as Cholesky of lambda
	 *	@note This function throws std::bad_alloc and std::runtime_error (when not-pos-def).
	 */
	inline void Refresh_R_FullR() // throw(std::bad_alloc, std::runtime_error)
	{
		_TyTimeSampler timer(m_shared_timer);

		m_n_last_full_R_update_size = m_lambda.n_BlockColumn_Num();
		// R will have the same size once updated ...

		if(!b_Have_NativeSolver) {
			if(b_Is_PoseOnly_SLAM) { // this is known at compile-time, should optimize the unused branch away
				m_R.Clear();
				m_r_system.r_Vertex_Pool().For_Each(CAlloc_RBlocks(m_R)); // won't work with VP problems, need to set correct ordering to vertices
			} else {
				//m_R.Clear(); // already inside PermuteTo()
				CUberBlockMatrix t_new_R;
				m_r_system.r_Vertex_Pool().For_Each(CAlloc_RBlocks(t_new_R));
				t_new_R.PermuteTo(m_R, m_p_lambda_block_ordering, m_n_lambda_block_ordering_size);
			}
			// only need to do there things if not using the native solver
		}
		// do the right thing and thrash R

		if(!m_linear_solver2.Factorize_PosDef_Blocky(m_R, m_lambda_perm, m_R_row_lookup_table, 0, 0))
			throw std::runtime_error("Factorize_PosDef_Blocky() failed to calculate full R");
		// factorize (uses cached cholesky, saves some time on allocation of workspace memory)

		timer.Accum_DiffSample(m_f_fullR_cholesky);

		{
			Collect_RightHandSide_Vector(m_v_dx);
			// collects the right-hand side vector (eta)

			timer.Accum_DiffSample(m_f_rhs_time);

			//++ m_n_full_forwardsubst_num;
			// do not count it here, we know how many times we did full R, it is the same count

			_ASSERTE(m_p_lambda_block_ordering);
			m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_dx(0), m_v_dx.rows(),
				m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
			m_R.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(&m_v_perm_temp(0), m_v_perm_temp.rows());
			m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
				m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
			// d = eta = eta/R
		}
		// convert eta to d

		timer.Accum_DiffSample(m_f_fullR_d);
	}

	/**
	 *	@brief refreshes the R matrix either from (pert of) lambda or from omega
	 *
	 *	@param[in] n_referesh_from_vertex is zero-based index of the first vertex that changes (unused)
	 *	@param[in] n_refresh_from_edge is zero-based index of the first edge that changes
	 *	@param[in] b_supress_reorder is 
	 *
	 *	@return Returns true if R was refreshed, false if it decided to take lambda fallback instead.
	 *
	 *	@note This function throws std::bad_alloc and std::runtime_error (when not-pos-def).
	 */
	inline bool b_Refresh_R(size_t UNUSED(n_referesh_from_vertex) = 0,
		size_t n_refresh_from_edge = 0, bool b_supress_reorder = false) // throw(std::bad_alloc, std::runtime_error)
	{
		_TyTimeSampler timer(m_shared_timer);

		// note that lambda is now up to date

		bool b_force_reorder = !n_refresh_from_edge && !b_supress_reorder; // if rebuilding whole lambda, it would be shame not to reorder
		// flag for forcing reorder

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		if(false) { // no optimizations for R up variants timing
#else // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		if(!b_supress_reorder) { // if allow optimizations ...
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			if(!b_force_reorder && m_lambda.n_BlockColumn_Num() > m_n_last_full_R_update_size + 10) { // it's always dense at the beginning // hlamfb 10
				size_t n_nnz = m_R.n_Storage_Size();
				float f_area = float(m_R.n_Column_Num()) * m_R.n_Column_Num();
				float f_density = n_nnz / f_area;
				if(f_density > 0.02f) {
					b_force_reorder = true;
					//printf("R became too dense (%.2f %%), forcing reorder\n", f_density * 100); // verbose
				}
			}
			// permit 2% density in R, then rebuild
		}
		// these two should just set b_force_reorder

		if(b_force_reorder) {
			// only calculate a new ordering on full refresh or if forced

			//printf("build new ordering ...\n");

			m_lambda_ordering.p_BlockOrdering(m_lambda,
				m_lambda_constraint.p_Get(m_lambda.n_BlockColumn_Num()),
				m_lambda.n_BlockColumn_Num(), true); // constrained blocky, calculate inverse as well
			m_p_lambda_block_ordering = m_lambda_ordering.p_Get_InverseOrdering(); // todo - make sure that the last vertex is a pose (otherwise we need to modify the constraint to select the pose, not the landmark)
			// get blockwise and elementwise ordering ...

			if(!b_Have_NativeSolver)
				m_R_row_lookup_table.clear(); // unused in native solver
			// can't reuse lookup of R's rows since these change with ordering
		} else if(m_lambda.n_BlockColumn_Num() > m_lambda_perm.n_BlockColumn_Num()) {
			// simply appends ordering with a new value (identity ordering at the end)

			//printf("extend ordering by " PRIsize "\n", m_lambda.n_BlockColumn_Num() - m_lambda_perm.n_BlockColumn_Num());

			m_p_lambda_block_ordering = m_lambda_ordering.p_InvertOrdering(
				m_lambda_ordering.p_ExtendBlockOrdering_with_Identity(m_lambda.n_BlockColumn_Num()),
				m_lambda.n_BlockColumn_Num());
			// get blockwise and elementwise ordering ...
		}
		m_n_lambda_block_ordering_size = m_lambda.n_BlockColumn_Num();
		// refresh/update the ordering (update means append with identity)

		_ASSERTE(CMatrixOrdering::b_IsValidOrdering(m_p_lambda_block_ordering, m_lambda.n_BlockColumn_Num()));
		// make sure that the ordering is good

		m_lambda.Permute_UpperTriangular_To(m_lambda_perm, m_p_lambda_block_ordering,
			m_lambda.n_BlockColumn_Num(), true);
		// make a reordered version of lambda (*always* changes if lambda changes)

		_ASSERTE(m_n_lambda_block_ordering_size == m_r_system.r_Vertex_Pool().n_Size());
		size_t n_order_min;
		if(b_force_reorder)
			n_order_min = 0; // a new ordering? from the ground up then ...
		else if(b_supress_reorder)
			n_order_min = 0; // the whole system changed numerically? from the ground up then ...
		else {
			n_order_min = m_p_lambda_block_ordering[m_n_lambda_block_ordering_size - 1];
			for(size_t i = n_refresh_from_edge, n = m_r_system.r_Edge_Pool().n_Size(); i < n; ++ i) {
				size_t n_order_v0 = m_p_lambda_block_ordering[m_r_system.r_Edge_Pool()[i].n_Vertex_Id(0)];
				size_t n_order_v1 = m_p_lambda_block_ordering[m_r_system.r_Edge_Pool()[i].n_Vertex_Id(1)]; // note that these are ids, but these equal order at the moment
				n_order_min = std::min(n_order_min, std::min(n_order_v0, n_order_v1));
			}
			//printf("loop size: " PRIsize "\n", m_n_lambda_block_ordering_size - 1 - n_order_min); // debug
		}
		// calculate min vertex order that needs to be updated (within the ordering!)

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
		m_n_loop_size_cumsum += m_n_lambda_block_ordering_size - 1 - n_order_min;
		// calculate how much loops did it process
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS

		timer.Accum_DiffSample(m_f_ordering_time);
		// stats

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		static bool b_first_time_dump = true;
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

		if(n_order_min > 0) {
			// if !n_order_min, L01 is empty and R merely equals chol(lambda)

			if(m_n_edges_in_lambda == m_r_system.r_Edge_Pool().n_Size()) {
				_ASSERTE(m_n_verts_in_lambda == m_lambda.n_BlockColumn_Num());
				_ASSERTE(m_n_verts_in_lambda == m_r_system.r_Vertex_Pool().n_Size());
				return true;
			}
			// this is final optimization, no need to refresh, there is no new edge / vertex

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			double p_inc_upd_times_start[] = {
				m_f_Rupdate_time,
				m_f_d_time
			};
			double p_omega_upd_times_start[] = {
				m_f_r11_omega_calc_time,
				m_f_r11_omega_slice_time,
				m_f_r11_omega_ata_time,
				m_f_r11_omega_add_time
			};
			double p_lambda_upd_times_start[] = {
				m_f_r11_lambda_slice_time,
				m_f_r11_lambda_ata_time,
				m_f_r11_lambda_add_time
			};
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

			const size_t n_order_max = m_lambda.n_BlockColumn_Num();

			// now: a) extend the reorder area to contain all the blocks in frontline (so lambda10 is null)
			//      b) freezee the few blocks in the frontline using constrained ordering? will that help?
			//      c) ?

			bool b_limited_search_region = false;
			bool b_blocks_above = false;

			m_lambda_perm.Get_UpperTriangular_BlockFrontline(lambda_perm_frontline);
			// get the frontline

			_ASSERTE(lambda_perm_frontline.size() == n_order_max);
			for(size_t i = n_order_min + 1; i < n_order_max /*- 1*/; ++ i) { // not sure why dont i check the last one - todo
				if(lambda_perm_frontline[i] < n_order_min) {
					b_blocks_above = true;
					break;
				}
			}
			// see if there are blocks above (except in the first or last column which are fixed)

#ifdef __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
			CUberBlockMatrix lambda11; // todo - make it a member?
#endif // __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
			if(!b_blocks_above) { // 4.2 seconds for full ordering always -> 2.3 seconds
				// this is the insufficient ordering only on lambda11
				// (causes a horrible fill-in if there are any blocks above)

#ifndef __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
				CUberBlockMatrix lambda11; // local is enough
#endif // !__NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
				m_lambda_perm.SliceTo(lambda11, n_order_min,
					n_order_max, n_order_min, n_order_max, true);
				_ASSERTE(lambda11.b_Square() && lambda11.b_SymmetricLayout());

				m_p_lambda11_block_ordering = m_lambda11_ordering.p_BlockOrdering(lambda11,
					m_lambda11_constraint.p_Get(lambda11.n_BlockColumn_Num()), lambda11.n_BlockColumn_Num());
				m_n_lambda_block11_ordering_size = lambda11.n_BlockColumn_Num();
			} else {
				size_t n_frontline_min = lambda_perm_frontline[n_order_min];
				for(size_t i = n_order_min + 1; i < n_order_max; ++ i) {
					if(n_frontline_min > lambda_perm_frontline[i])
						n_frontline_min = lambda_perm_frontline[i];
				}
				// find a minimal frontline value over the updated submatrix

				size_t n_context_min = 0;
				if(n_frontline_min > 0) {
					for(size_t i = n_order_min; i > 0;) { // note that the loop terminates itself, the condition is not required (i think)
						-- i;
						if(n_frontline_min >= i) {
							n_context_min = i;
							// in case frontline is at the same level as a submatrix, we can order using a smaller matrix

							break;
						}
						//if(n_frontline_min > lambda_perm_frontline[i])
						//	n_frontline_min = lambda_perm_frontline[i]; // do not do this, will enlarge the search window too much, in order to include elements that can not be reordered anyway (seems like a good idea now)
					}
					// see if we can expand the ordering submatrix a little, to avoid using the whole matrix

#ifdef _DEBUG
					if(n_context_min > 0) {
						//for(size_t i = n_context_min + 1; i < n_order_max - 1; ++ i) // forbids first / last // do not do this
						for(size_t i = n_order_min + 1; i < n_order_max - 1; ++ i) // in case the expanding over previous area is ignored
							_ASSERTE(lambda_perm_frontline[i] >= n_context_min);
						/*printf("can order on " PRIsize " x " PRIsize " instead of " PRIsize " x " PRIsize
							" (update is " PRIsize " x " PRIsize ")\n", n_order_max - n_context_min,
							n_order_max - n_context_min, n_order_max, n_order_max, n_order_max - n_order_min,
							n_order_max - n_order_min);*/
					}
					// debug - see if we can use it and tell us
#endif // _DEBUG
				}

				m_n_lambda_block11_ordering_size = n_order_max - n_order_min;
				_ASSERTE(m_lambda_perm.n_BlockColumn_Num() == m_lambda.n_BlockColumn_Num());

				if(n_context_min > n_order_max / 8 /*std::min(size_t(100), n_order_max / 16)*/) { // t_odo - need to dump n_context_min, probably citytrees have one at 0 or something, that can not be ordered away (make it go away by forcing it somewhere else?l)
					// this is prefix-constrained ordering on part of lambda perm (not full, not update)
					// this works rather well

					b_limited_search_region = true;
					// say we did it

					_ASSERTE(n_order_min > n_context_min);
					size_t n_skip_blocks = n_context_min; // how much to slice
					size_t n_skip_order = n_order_min - n_context_min; // how much to replace by diagonal
					size_t n_order_size = n_order_max - n_context_min; // size of the lambda11 part
					_ASSERTE(n_order_size > n_skip_order);

					const size_t *p_order = m_lambda_alt_ordering.p_BlockOrdering_MiniSkirt(m_lambda_perm,
						n_skip_blocks, n_order_min, m_lambda_alt_constraint.p_Get(n_order_size,
						n_skip_order), n_order_size);
					// get ordering on the submatrix of permuted lambda, without calculating the submatrix,
					// also with less items on the 

#ifdef _DEBUG
					for(size_t i = 0; i < n_skip_order; ++ i)
						_ASSERTE(p_order[i] == i);
					// the prefix should be identity

					std::vector<bool> coverage;
					coverage.resize(m_n_lambda_block11_ordering_size, false);
					// make sure we produce a valid ordering
#endif // _DEBUG

					size_t *p_order11 = (size_t*)p_order + n_skip_order;
					for(size_t i = 0; i < m_n_lambda_block11_ordering_size; ++ i) {
						_ASSERTE(p_order11[i] >= n_skip_order);
						p_order11[i] -= n_skip_order;
						_ASSERTE(p_order11[i] < m_n_lambda_block11_ordering_size);

#ifdef _DEBUG
						_ASSERTE(!coverage[p_order11[i]]); // no repeated ones
						coverage[p_order11[i]] = true;
						// make sure we produce a valid ordering
#endif // _DEBUG
					}
					// just subtract to get ordering on lambda11

#ifdef _DEBUG
					_ASSERTE(std::find(coverage.begin(), coverage.end(), false) == coverage.end());
					// make sure all elements of lambda11 are covered
#endif // _DEBUG

					m_p_lambda11_block_ordering = p_order11;
				} else {
					// this is prefix-constrained ordering on the full lambda (perm)
					// this works rather well

					const size_t *p_order = m_lambda_alt_ordering.p_BlockOrdering_MiniSkirt(m_lambda_perm,
						0, n_order_min, m_lambda_alt_constraint.p_Get(m_lambda.n_BlockColumn_Num(),
						n_order_min), m_lambda.n_BlockColumn_Num());
					// just give diagonal matrix all the way from 0 to n_order_min, then the actual
					// sparsity pattern till n_order_max
					// get ordering on the permuted lambda
					// gets the full ordering on lambda, making sure that the prefix is the same and only the lambda11 suffix changes

#ifdef _DEBUG
					for(size_t i = 0; i < n_order_min; ++ i)
						_ASSERTE(p_order[i] == i);
					// the prefix should be identity

					std::vector<bool> coverage;
					coverage.resize(m_n_lambda_block11_ordering_size, false);
					// make sure we produce a valid ordering
#endif // _DEBUG

					size_t *p_order11 = (size_t*)p_order + n_order_min;
					for(size_t i = 0; i < m_n_lambda_block11_ordering_size; ++ i) {
						_ASSERTE(p_order11[i] >= n_order_min);
						p_order11[i] -= n_order_min;
						_ASSERTE(p_order11[i] < m_n_lambda_block11_ordering_size);

#ifdef _DEBUG
						_ASSERTE(!coverage[p_order11[i]]); // no repeated ones
						coverage[p_order11[i]] = true;
						// make sure we produce a valid ordering
#endif // _DEBUG
					}
					// just subtract to get ordering on lambda11

#ifdef _DEBUG
					_ASSERTE(std::find(coverage.begin(), coverage.end(), false) == coverage.end());
					// make sure all elements of lambda11 are covered
#endif // _DEBUG

					m_p_lambda11_block_ordering = p_order11;
					//m_p_lambda11_block_ordering = m_lambda11_ordering.p_InvertOrdering(p_order11,
					//	m_n_lambda_block11_ordering_size); // inverse
					// invert it for merging with the original ordering
				}

#ifdef __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
				//CUberBlockMatrix lambda11; // above; scope / define hell
				m_lambda_perm.SliceTo(lambda11, n_order_min,
					n_order_max, n_order_min, n_order_max, true);
				_ASSERTE(lambda11.b_Square() && lambda11.b_SymmetricLayout());
				// will need this for verification purposes, otherwise not required by this method
#endif // __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
			}

			bool b_identity_ordering = true;
			for(size_t i = 0; i < m_n_lambda_block11_ordering_size; ++ i) {
				if(m_p_lambda11_block_ordering[i] != i) {
					b_identity_ordering = false;
					break;
				}
			}
			// calculates new block ordering for lambda11 (the updated area of lambda)

#ifdef __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
			CUberBlockMatrix lambda00_p, lambda11_p;
			m_lambda_perm.SliceTo(lambda00_p, 0, n_order_min, 0, n_order_min, false); // make a deep copy
			lambda11.Permute_UpperTriangular_To(lambda11_p, m_p_lambda11_block_ordering,
				m_n_lambda_block11_ordering_size, false); // make a deep copy
			// copy lambda 00 and lambda 11 (don't care about lambda 01, it is hard to permute correctly at this point)
#endif // __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING

			_TyTime f_ordering11_time = 0;
			timer.Accum_DiffSample(f_ordering11_time);

			if(!b_identity_ordering) {
				if(!b_Have_NativeSolver)
					m_R_row_lookup_table.clear(); // unused in native solver
				// !! we defined a new ordering

				const size_t *p_order;
				m_p_lambda_block_ordering = m_lambda_ordering.p_InvertOrdering(p_order =
					m_lambda_ordering.p_ExtendBlockOrdering_with_SubOrdering(n_order_min,
					m_p_lambda11_block_ordering, m_n_lambda_block11_ordering_size), m_lambda.n_BlockColumn_Num());
				_ASSERTE(m_n_lambda_block_ordering_size == n_order_min + m_n_lambda_block11_ordering_size);
				// update the ordering (update means append with lambda11 sub-block ordering)
				// this is quick, no bottleneck in here

				_ASSERTE(CMatrixOrdering::b_IsValidOrdering(m_p_lambda_block_ordering,
					m_lambda.n_BlockColumn_Num()));
				// make sure that the ordering is good

				timer.Accum_DiffSample(m_f_ordering_fold_time);

				m_lambda.Permute_UpperTriangular_To(m_lambda_perm, m_p_lambda_block_ordering,
					m_lambda.n_BlockColumn_Num(), true, n_order_min, true);
				// note that this does leave allocated blocks, but in the next round, lambda will reperm and free those

#ifdef __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
				CUberBlockMatrix lambda00_r, lambda11_r;
				m_lambda_perm.SliceTo(lambda00_r, 0, n_order_min, 0, n_order_min, false);
				m_lambda_perm.SliceTo(lambda11_r, n_order_min, n_order_max, n_order_min, n_order_max, false);
				// make copies of the new permutated lambda; should be identical to what was intended

				lambda00_p.AddTo(lambda00_r, -1);
				lambda11_p.AddTo(lambda11_r, -1);

				double f_diff0 = lambda00_r.f_Norm();
				double f_diff1 = lambda11_r.f_Norm();
				printf(" %g/%g", f_diff0, f_diff1);
#endif // __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING

				timer.Accum_DiffSample(m_f_repermute_time);

				m_R.SliceTo(m_R, n_order_min, n_order_min, true); // 3 seconds -> 0.3 seconds

				timer.Accum_DiffSample(m_f_Rslice_time);

				if(m_chol_bitfield.capacity() < n_order_max) {
					m_chol_bitfield.clear();
					m_chol_bitfield.reserve(std::max(n_order_max, 2 * m_chol_bitfield.capacity()));
				}
				m_chol_bitfield.resize(n_order_max, 0);

				m_lambda_perm.Build_EliminationTree(m_chol_etree, m_chol_ereach_stack); // use ereach stack as workspace
				_ASSERTE(m_chol_ereach_stack.size() == n_order_max);
				// build an elimination tree

				timer.Accum_DiffSample(m_f_etree_time);

				if(!m_R.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(m_lambda_perm, m_chol_etree,
				   m_chol_ereach_stack, m_chol_bitfield, n_order_min)) // todo - do incremental etree as well, might save considerable time
					throw std::runtime_error("error: got not pos def in incR section anyways"); // does not really happen
				// calcualte updated R11 and R10 using resumed Cholesky

				++ m_n_resumed_chol_num;
				timer.Accum_DiffSample(m_f_resumed_chol_time);

				Refresh_d_IncR11(n_refresh_from_edge, n_order_min); // timing inside
				// all that is left to do is to refresh d
			} else {
				Refresh_R_IncR11(n_refresh_from_edge, n_order_min); // timing inside
				// run the "fast" refresh of R
			}
			// choose between progressive reordering and "fast" update to R

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
			char p_s_filename[256];

			CUberBlockMatrix lambda_perm;
			m_lambda.Permute_UpperTriangular_To(lambda_perm, m_p_lambda_block_ordering,
				m_lambda.n_BlockColumn_Num(), true);
			// need to reperm, may only have a part of lambda_perm, effectively selecting everything above n_order_min as a new nnz

			size_t n_verts_in_lambda = m_lambda.n_BlockColumn_Num();
			sprintf(p_s_filename, "rss2013/%05d_6_lambda-perm.tga", n_verts_in_lambda);
			//if(n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			//	lambda_perm.Rasterize_Symmetric(p_s_filename, (n_verts_in_lambda < 750 * 6 / _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 3); // do not really need lambdas right now

			sprintf(p_s_filename, "rss2013/%05d_7_R.tga", n_verts_in_lambda);
			//if(n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			{
				int n_ss; // scalar size
				TBmp *p_img = m_R.p_Rasterize(lambda_perm, false, 0, n_ss = ((n_verts_in_lambda < 750 * 6 / _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 3 : 2)); // highlight fill-in
				if(p_img) {
					CTgaCodec::Save_TGA(p_s_filename, *p_img, false, true);

					sprintf(p_s_filename, "rss2013/%05d_8_R_marked.tga", n_verts_in_lambda);
					//printf("drawing\n");
					int n_line0 = (n_ss - 1) * m_R.n_BlockColumn_Base(n_order_min);
					//int n_line1 = (n_ss - 1) * m_R.n_BlockColumn_Base(n_context_min);
					p_img->DrawLine(n_line0, n_line0, n_line0, p_img->n_height, 0xfff38630U, 8);
					p_img->DrawLine(n_line0, n_line0, p_img->n_width, n_line0, 0xfff38630U, 8); // draw line that separates order min
					for(int y = n_line0, h = p_img->n_height; y < h; ++ y) {
						for(int x = n_line0, w = p_img->n_width; x < w; ++ x) {
							if(p_img->p_buffer[x + w * y] == 0xffffffffU) // white?
								p_img->p_buffer[x + w * y] = 0xffffffbbU; // yellowish
						}
					}
					//p_img->DrawLine(n_line1, n_line1, n_line1, p_img->n_height, 0xff00ff00U, 3);
					//p_img->DrawLine(n_line1, n_line1, p_img->n_width, n_line1, 0xff00ff00U, 3); // draw line that separates context min
					//printf("saving\n");
					CTgaCodec::Save_TGA(p_s_filename, *p_img, false, true);
					//printf("deleting\n");
					p_img->Delete();
				} else
					fprintf(stderr, "error: not enough memory to rasterize the matrix\n");
			}

			sprintf(p_s_filename, "rss2013/%05d_9_stats.txt", n_verts_in_lambda);
			FILE *p_fw;
			if((p_fw = fopen(p_s_filename, "w"))) {
				fprintf(p_fw, PRIsize "\n", m_lambda_perm.n_Block_Num()); // save density of lambda
				fprintf(p_fw, PRIsize "\n", m_R.n_Block_Num()); // save density of R
				fprintf(p_fw, PRIsize "\n", n_order_min); // save density of R
				fclose(p_fw);
			}
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

			if(b_blocks_above) {
				if(b_limited_search_region) {
					++ m_n_limited_search_num;
					m_f_ordering11_part_time += f_ordering11_time;
				} else {
					++ m_n_blocks_above_num;
					m_f_ordering11_full_time += f_ordering11_time;
				}
			} else
				m_f_ordering11_time += f_ordering11_time;

			++ m_n_Rup_num;
			// count incremental R updates

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			double p_inc_upd_times[] = {
				m_f_Rupdate_time,
				m_f_d_time
			};
			double f_inc_upd_sum = 0;
			for(int i = 0; i < 2; ++ i) {
				p_inc_upd_times[i] -= p_inc_upd_times_start[i];
				f_inc_upd_sum += p_inc_upd_times[i];
			}
			double p_omega_upd_times[] = {
				m_f_r11_omega_calc_time,
				m_f_r11_omega_slice_time,
				m_f_r11_omega_ata_time,
				m_f_r11_omega_add_time
			};
			double f_omega_upd_sum = f_inc_upd_sum;
			for(int i = 0; i < 4; ++ i) {
				p_omega_upd_times[i] -= p_omega_upd_times_start[i];
				f_omega_upd_sum += p_omega_upd_times[i];
			}
			double p_lambda_upd_times[] = {
				m_f_r11_lambda_slice_time,
				m_f_r11_lambda_ata_time,
				m_f_r11_lambda_add_time
			};
			double f_lambda_upd_sum = f_inc_upd_sum;
			for(int i = 0; i < 3; ++ i) {
				p_lambda_upd_times[i] -= p_lambda_upd_times_start[i];
				f_lambda_upd_sum += p_lambda_upd_times[i];
			}
			// calculate times

			size_t n_loop_size = m_n_lambda_block_ordering_size - 1 - n_order_min;
			bool b_had_omega_upd = f_omega_upd_sum > f_inc_upd_sum; // omega update only if it can
			_ASSERTE(f_lambda_upd_sum > 0); // lambda update always

			double f_full_R_start = m_shared_timer.f_Time();
			Refresh_R_FullR();
			double f_full_R_time = m_shared_timer.f_Time() - f_full_R_start;
			// measure full R as well

			FILE *p_fw = fopen("Rup_variants_time.txt", (b_first_time_dump)? "w" : "a");
			if(b_first_time_dump) {
				fprintf(p_fw, "verts-in-R;loop-size;full-R-time;lambda-up-time;lambda-slice-time;"
					"lambda-ata-time;lambda-add-time;omega-up-time;omega-calc-time;"
					"omega-slice-time;omega-ata-time;omega-add-time\n");
			}
			fprintf(p_fw, "" PRIsize ";" PRIsize ";%f;%f;%f;%f;%f", m_R.n_BlockColumn_Num(), n_loop_size, f_full_R_time,
				f_lambda_upd_sum, p_lambda_upd_times[0], p_lambda_upd_times[1], p_lambda_upd_times[2]);
			if(b_had_omega_upd) {
				fprintf(p_fw, ";%f;%f;%f;%f;%f\n", f_omega_upd_sum, p_omega_upd_times[0],
					p_omega_upd_times[1], p_omega_upd_times[2], p_omega_upd_times[3]);
			} else {
				fprintf(p_fw, ";;;;;\n");
			}
			fclose(p_fw);
			// print timing to a file
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		} else {
			//printf("doing full R\n"); // debug

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			double f_full_R_start = m_shared_timer.f_Time();
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

			++ m_n_full_R_num;
			// R is not up to date, need to rebuild from scratch

			Refresh_R_FullR();
			// do the "full" R = chol(lambda)

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			double f_full_R_time = m_shared_timer.f_Time() - f_full_R_start;
			// measure full R

			size_t n_loop_size = m_n_lambda_block_ordering_size - 1 - n_order_min;

			FILE *p_fw = fopen("Rup_variants_time.txt", (b_first_time_dump)? "w" : "a");
			if(b_first_time_dump) {
				fprintf(p_fw, "verts-in-R;loop-size;full-R-time;lambda-up-time;lambda-slice-time;"
					"lambda-ata-time;lambda-add-time;omega-up-time;omega-calc-time;"
					"omega-slice-time;omega-ata-time;omega-add-time\n");
			}
			fprintf(p_fw, "" PRIsize ";" PRIsize ";%f;;;;", m_R.n_BlockColumn_Num(), n_loop_size, f_full_R_time); // no lambda upd
			fprintf(p_fw, ";;;;;\n"); // no omega upd
			fclose(p_fw);
			// print timing to a file
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
			char p_s_filename[256];

			size_t n_verts_in_lambda = m_lambda.n_BlockColumn_Num();
			sprintf(p_s_filename, "rss2013/%05d_6_lambda-perm.tga", n_verts_in_lambda);
			//if(n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			//	lambda_perm.Rasterize_Symmetric(p_s_filename, (n_verts_in_lambda < 750 * 6 / _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 3); // do not really need lambdas right now

			sprintf(p_s_filename, "rss2013/%05d_7_R.tga", n_verts_in_lambda);
			//if(n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			{
				int n_ss; // scalar size
				TBmp *p_img = m_R.p_Rasterize(m_lambda_perm, false, 0, n_ss = ((n_verts_in_lambda < 750 * 6 / _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 3 : 2)); // highlight fill-in
				if(p_img) {
					CTgaCodec::Save_TGA(p_s_filename, *p_img, false, true);

					sprintf(p_s_filename, "rss2013/%05d_8_R_marked.tga", n_verts_in_lambda);
					//printf("drawing\n");
					int n_line0 = 0;//(n_ss - 1) * m_R.n_BlockColumn_Base(n_order_min);
					//int n_line1 = (n_ss - 1) * m_R.n_BlockColumn_Base(n_context_min);
					p_img->DrawLine(n_line0, n_line0, n_line0, p_img->n_height, 0xfff38630U, 8);
					p_img->DrawLine(n_line0, n_line0, p_img->n_width, n_line0, 0xfff38630U, 8); // draw line that separates order min
					for(int y = n_line0, h = p_img->n_height; y < h; ++ y) {
						for(int x = n_line0, w = p_img->n_width; x < w; ++ x) {
							if(p_img->p_buffer[x + w * y] == 0xffffffffU) // white?
								p_img->p_buffer[x + w * y] = 0xffffffbbU;//0xffddddddU; // grayish // yellowish
						}
					}
					//p_img->DrawLine(n_line1, n_line1, n_line1, p_img->n_height, 0xff00ff00U, 3);
					//p_img->DrawLine(n_line1, n_line1, p_img->n_width, n_line1, 0xff00ff00U, 3); // draw line that separates context min
					//printf("saving\n");
					CTgaCodec::Save_TGA(p_s_filename, *p_img, false, true);
					//printf("deleting\n");
					p_img->Delete();
				} else
					fprintf(stderr, "error: not enough memory to rasterize the matrix\n");
			}

			sprintf(p_s_filename, "rss2013/%05d_9_stats.txt", n_verts_in_lambda);
			FILE *p_fw;
			if((p_fw = fopen(p_s_filename, "w"))) {
				fprintf(p_fw, PRIsize "\n", m_lambda_perm.n_Block_Num()); // save density of lambda
				fprintf(p_fw, PRIsize "\n", m_R.n_Block_Num()); // save density of R
				fclose(p_fw);
			}
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
		}

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		b_first_time_dump = false;
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

		return true;
	}

	/**
	 *	@brief refreshes the lambda matrix by recalculating edge hessians
	 *
	 *	@param[in] n_referesh_from_vertex is zero-based index of the first vertex that changes
	 *	@param[in] n_refresh_from_edge is zero-based index of the first edge that changes
	 */
	inline void Refresh_Lambda(size_t n_referesh_from_vertex = 0, size_t n_refresh_from_edge = 0)
	{
		if(n_refresh_from_edge) {
			m_r_system.r_Edge_Pool().For_Each_Parallel(n_refresh_from_edge,
				m_r_system.r_Edge_Pool().n_Size(), CCalculate_Lambda()); // this is not ok, the update will contain some contributions that are not supposed to be there
			// in order to fix this, we need to modify the call to vertex update (the vertices need to specifically ignore contributions from some edges)
		} else
			m_r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Lambda()); // this is ok, the update is going to be calculated correctly
		if(n_referesh_from_vertex) {
			m_r_system.r_Vertex_Pool().For_Each_Parallel(n_referesh_from_vertex,
				m_r_system.r_Vertex_Pool().n_Size(), CCalculate_Lambda());
		} else
			m_r_system.r_Vertex_Pool().For_Each_Parallel(CCalculate_Lambda()); // this is currently always from vert 0
		// can do this in parallel

		if(!n_referesh_from_vertex) {
			const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
			m_lambda.t_FindBlock(0, 0).noalias() += r_t_uf.transpose() * r_t_uf;
			// for lambda yes
		}
		// add unary factor (gets overwritten by the first vertex' block)
	}

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY

	/**
	 *	@brief dumps density of R factor, given different ordering strategies
	 *	@note This is only available if __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY is defined.
	 */
	void Dump_RDensity()
	{
		size_t n_nnz_ideal, n_nnz_ideal_elem;
		{
			cs *A = m_lambda.p_Convert_to_Sparse();
			css *S;
			S = cs_schol(1, A); // use AMD
			csn *N = cs_chol(A, S);
			cs *L = N->L;

			n_nnz_ideal_elem = L->p[L->n]; // @todo - count this by bylocks

			CUberBlockMatrix blockL;
			m_lambda.CopyLayoutTo(blockL);
			std::vector<size_t> workspace;
			blockL.From_Sparse(0, 0, L, false, workspace);
			n_nnz_ideal = blockL.n_NonZero_Num();

			cs_sfree(S);
			cs_spfree(N->U);
			cs_free(N->pinv);
			cs_free(N->B);
			cs_spfree(L);
			cs_spfree(A);
			// calculate cholesky with elementwise AMD
		}
		// "ideal" ordering

		size_t n_nnz_blocky, n_nnz_blocky_elem;
		{
			CMatrixOrdering mord;
			mord.p_BlockOrdering(m_lambda, true); // unconstrained, calculate inverse right away
			const size_t *p_order = mord.p_Get_InverseOrdering();
			// unconstrained blocky ordering on lambda

			CUberBlockMatrix lord;
			m_lambda.Permute_UpperTriangular_To(lord, p_order,
				m_lambda.n_BlockColumn_Num(), true);
			// order the matrix

			cs *A = lord.p_Convert_to_Sparse();
			css *S;
			S = cs_schol(0, A); // use AMD
			csn *N = cs_chol(A, S);
			cs *L = N->L;

			n_nnz_blocky_elem = L->p[L->n]; // @todo - count this by bylocks

			CUberBlockMatrix _blockL, blockL;
			m_lambda.CopyLayoutTo(_blockL);
			_blockL.Permute_UpperTriangular_To(blockL, p_order, m_lambda.n_BlockColumn_Num());
			std::vector<size_t> workspace;
			blockL.From_Sparse(0, 0, L, false, workspace);
			n_nnz_blocky = blockL.n_NonZero_Num();

			cs_sfree(S);
			cs_spfree(N->U);
			cs_free(N->pinv);
			cs_free(N->B);
			cs_spfree(L);
			cs_spfree(A);
			// calculate cholesky with natural ordering (= no ordering)
		}
		// blocky ordering

		size_t n_nnz_blocky_constr, n_nnz_blocky_constr_elem;
		{
			CLastElementOrderingConstraint constr;
			const size_t *p_constraint = constr.p_Get(m_lambda.n_BlockColumn_Num());
			CMatrixOrdering mord;
			mord.p_BlockOrdering(m_lambda, p_constraint, m_lambda.n_BlockColumn_Num(), true);
			const size_t *p_order = mord.p_GetInverseOrdering();
			// unconstrained blocky ordering on lambda

			CUberBlockMatrix lord;
			m_lambda.Permute_UpperTriangular_To(lord, p_order,
				m_lambda.n_BlockColumn_Num(), true);
			// order the matrix

			cs *A = lord.p_Convert_to_Sparse();
			css *S;
			S = cs_schol(0, A); // use AMD
			csn *N = cs_chol(A, S);
			cs *L = N->L;

			n_nnz_blocky_constr_elem = L->p[L->n]; // @todo - count this by bylocks

			CUberBlockMatrix _blockL, blockL;
			m_lambda.CopyLayoutTo(_blockL);
			_blockL.Permute_UpperTriangular_To(blockL, p_order, m_lambda.n_BlockColumn_Num());
			std::vector<size_t> workspace;
			blockL.From_Sparse(0, 0, L, false, workspace);
			n_nnz_blocky_constr = blockL.n_NonZero_Num();

			cs_sfree(S);
			cs_spfree(N->U);
			cs_free(N->pinv);
			cs_free(N->B);
			cs_spfree(L);
			cs_spfree(A);
			// calculate cholesky with natural ordering (= no ordering)
		}
		// constrained blocky ordering

		size_t n_nnz_actual = m_R.n_NonZero_Num();
		// actual NNZ

		FILE *p_fw = fopen("RDensityByOrdering.txt", (m_n_real_step > 0)? "a" : "w");
		if(!m_n_real_step) {
			fprintf(p_fw, "block-cols;amd;blocky-amd;blocky-constrained-amd;"
				"amd-blocks;blocky-amd-blocks;blocky-constrained-amd-blocks;actual-R-nnz-blocks\n");
		}
		fprintf(p_fw, "" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize "\n", m_lambda.n_BlockColumn_Num(), n_nnz_ideal_elem,
			n_nnz_blocky_elem, n_nnz_blocky_constr_elem, n_nnz_ideal, n_nnz_blocky,
			n_nnz_blocky_constr, n_nnz_actual);
		fclose(p_fw);
	}

#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY

	/**
	 *	@brief function object that calculates and accumulates chi^2 from all the edges
	 */
	class CSum_ChiSquareError {
	protected:
		double m_f_sum; /**< @brief a running sum of chi-square errors */

	public:
		/**
		 *	@brief default constructor
		 */
		inline CSum_ChiSquareError()
			:m_f_sum(0)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in] r_t_edge is edge to output its part R error vector
		 */
		template <class _TyEdge>
		inline void operator ()(const _TyEdge &r_t_edge)
		{
			m_f_sum += r_t_edge.f_Chi_Squared_Error();
		}

		/**
		 *	@brief gets the current value of the accumulator
		 *	@return Returns the current value of the accumulator.
		 */
		inline operator double() const
		{
			return m_f_sum;
		}
	};

	/**
	 *	@brief function object that updates states of all the vertices
	 */
	class CUpdateEstimates {
	protected:
		const Eigen::VectorXd &m_r_dx; /**< @brief vector of differences */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_dx is reference to the vector of differences
		 */
		inline CUpdateEstimates(const Eigen::VectorXd &r_dx)
			:m_r_dx(r_dx)
		{}

		/**
		 *	@brief updates vertex state
		 *	@tparam _TyVertex is type of vertex
		 *	@param[in,out] r_t_vertex is reference to vertex, being updated
		 */
		template <class _TyVertex>
		inline void operator ()(_TyVertex &r_t_vertex) const
		{
			r_t_vertex.Operator_Plus(m_r_dx);
		}
	};

	/**
	 *	@brief updates states of all the vertices with the error vector
	 */
	inline void PushValuesInGraphSystem(const Eigen::VectorXd &r_v_dx)
	{
		m_r_system.r_Vertex_Pool().For_Each_Parallel(CUpdateEstimates(r_v_dx)); // can do this in parallel
	}

	/**
	 *	@brief function object that calculates hessians in all the edges
	 */
	class CCalculate_Lambda {
	public:
		/**
		 *	@brief function operator
		 *	@tparam _TyVertexOrEdge is edge type
		 *	@param[in] r_t_vertex_or_edge is vertex or edge to update it's hessians
		 */
		template <class _TyVertexOrEdge>
		inline void operator ()(_TyVertexOrEdge &r_t_vertex_or_edge) const
		{
			r_t_vertex_or_edge.Calculate_Hessians();
		}
	};

	CNonlinearSolver_FastL(const CNonlinearSolver_FastL &UNUSED(r_solver)); /**< @brief the object is not copyable */
	CNonlinearSolver_FastL &operator =(const CNonlinearSolver_FastL &UNUSED(r_solver)) { return *this; } /**< @brief the object is not copyable */
};

#endif // __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED
