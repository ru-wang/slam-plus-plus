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

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
 *	@brief if defined, solution is calcualted at each step, otherwise it is calculated
 *		only at the specified intervals (but is ready to be calculated at each step quickly)
 */
#define __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
 *	@brief if defined, L is always used to calculate updates, even in the subsequent iterations,
 *		so that when the solver finishes, L is up to date and does not need to be refreshed again
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
 *	@brief enables writes of density of L given different ordering strategies
 */
//#define __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY

/**
 *	@def __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
 *	@brief if defined, enables writes of timing of different L update variants
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
 *	@def __NONLIVEAR_SOLVER_FAST_L_PARALLEL_MATMULT_THRESH
 *	@brief matrix size threshold for parallel multiplication
 *		in L update (in blocks, used to calculate L10L10T and L11L11T)
 */
#define __NONLIVEAR_SOLVER_FAST_L_PARALLEL_MATMULT_THRESH 200

/**
 *	@def __NONLINEAR_SOLVER_FASTL_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
 *	@brief dump RSS 2013 matrix animation data
 */
//#define __NONLINEAR_SOLVER_FASTL_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

//extern int n_dummy_param; // remove me

#include "slam/OrderingMagic.h"

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

	typedef typename CUniqueTypelist<CAMatrixBlockSizes>::_TyResult _TyAMatrixBlockSizes; /**< @brief possible block matrices, that can be found in A */
	typedef typename __fbs_ut::CBlockSizesAfterPreMultiplyWithSelfTranspose<
		_TyAMatrixBlockSizes>::_TyResult _TyLambdaMatrixBlockSizes; /**< @brief possible block matrices, found in lambda and L */

	/**
	 *	@brief some run-time constants, stored as enum
	 */
	enum {
		b_Is_PoseOnly_SLAM = CTypelistLength<_TyAMatrixBlockSizes>::n_result == 1, /**< @brief determines if we're doing pose-only SLAM (10k) */
		b_Have_NativeSolver = CEqualType<CLinearSolver_UberBlock<_TyLambdaMatrixBlockSizes>,
			 _TyLinearSolver>::b_result /**< @brief determines if the native linear solver is being used */
	};

protected:
	CSystem &m_r_system; /**< @brief reference to the system */
	CLinearSolver m_linear_solver; /**< @brief linear solver */
	CLinearSolver m_linear_solver2; /**< @brief linear solver for calculating cholesky of L and cholesky of L increment */

	std::vector<size_t> m_chol_etree; /**< @brief reusable e-tree storage */
	std::vector<size_t> m_chol_ereach_stack; /**< @brief reusable workspace for Cholesky */
	std::vector<size_t> m_chol_bitfield; /**< @brief reusable workspace for Cholesky */

	CUberBlockMatrix m_L; /**< @brief the L matrix (built / updated incrementally) */

	bool m_b_had_loop_closure; /**< @brief (probable) loop closure flag */
	bool m_b_first_iteration_use_L; /**< @brief flag for using the L matrix or rather lambda in the first iteration of nonlinear optimization */
	bool m_b_L_up_to_date; /**< @brief dirty flag for the L matrix (required to keep track after lambda updates and linearization point changes) */
	size_t m_n_last_full_L_update_size; /**< @brief the last number of block columns in L when it was fully updated */
	std::vector<size_t> m_L_row_lookup_table; /**< @brief row lookup table for L (used by b_Refresh_L() and Refresh_L11()) */
	//size_t m_n_big_loop_threshold; /**< @brief threshold for what is considered a "big" loop (incrementing L is avoided) */

	CMatrixOrdering m_lambda_ordering; /**< @brief lambda block ordering calculator (CAMD wrapper) */
	const size_t *m_p_lambda_block_ordering; /**< @brief lambda block ordering (only valid if m_b_L_up_to_date is set) */ // todo - convert all those to size_t
	size_t m_n_lambda_block_ordering_size; /**< @brief lambda block ordering size */
	CUberBlockMatrix m_lambda_perm; /**< @brief the reordered reference to the lambda matrix */
	CUberBlockMatrix m_lambda; /**< @brief the lambda matrix (built / updated incrementally) */

	CMatrixOrdering m_lambda11_ordering; /**< @brief lambda11 block ordering calculator (CAMD wrapper) */
	const size_t *m_p_lambda11_block_ordering; /**< @brief lambda block ordering (only valid if m_b_L_up_to_date is set) */ // todo - convert all those to size_t
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
	double m_f_nonlinear_solve_error_threshold; /**< @brief error threshold in incremental nonlinear solve */ // t_odo - document these in elementwise A and L
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

	size_t m_n_iteration_num; /**< @brief number of linear solver iterations */
	_TyTime m_f_chol_time; /**< @brief time spent in Choleski() section */
	_TyTime m_f_norm_time; /**< @brief time spent in norm calculation section */
	_TyTime m_f_vert_upd_time; /**< @brief time spent in updating the vertices */

	size_t m_n_full_forwardsubst_num; /**< @brief number of d updates performed using full L forward substitution */
	size_t m_n_resumed_forwardsubst_num; /**< @brief number of d updates performed using resumed L forward substitution */
	size_t m_n_resumed_perm_forwardsubst_num; /**< @brief number of d updates performed using resumed L forward substitution, while being permutated in the updated area */
	size_t m_n_L_optim_num; /**< @brief number of system optimizations performed using L backsubstitution */
	size_t m_n_lambda_optim_num; /**< @brief number of system optimizations performed using cholsol(lambda) */
	size_t m_n_Lup_num; /**< @brief number of L increments */
	size_t m_n_omega_update_num; /**< @brief number of L increments calculated using omega */
	size_t m_n_lambda_update_num; /**< @brief number of L increments calculated using lambda */
	size_t m_n_full_L_num; /**< @brief number of L updates */
	_TyTime m_f_lambda_refresh_time; /**< @brief time spent in updating and allocating lambda */
	_TyTime m_f_rhs_time; /**< @brief time spent in updating right-hand side vector */
	_TyTime m_f_ordering_time; /**< @brief time spent calculating ordering of lambda */
	_TyTime m_f_fullL_d; /**< @brief time spent in updating d while doing full L */
	_TyTime m_f_l11_omega_calc_time; /**< @brief time spent calculating omega (L increment) */
	_TyTime m_f_l11_omega_slice_time; /**< @brief time spent in slicing L11 (L increment) */
	_TyTime m_f_l11_omega_ata_time; /**< @brief time spent calculating L11TL11 (L increment) */
	_TyTime m_f_l11_omega_add_time; /**< @brief time spent adding L11TL11 + omega (L increment) */
	_TyTime m_f_l11_lambda_slice_time; /**< @brief time spent in slicing lambda11 and L01 (L increment) */
	_TyTime m_f_l11_lambda_ata_time; /**< @brief time spent calculating L01TL01 (L increment) */
	_TyTime m_f_l11_lambda_add_time; /**< @brief time spent adding L01TL01 + lambda11 (L increment) */
	_TyTime m_f_lupdate_time; /**< @brief time spent calculating cholesky of new L11 (L increment) */
	_TyTime m_f_d_time; /**< @brief time spent updating d (right hand side vector) */
	_TyTime m_f_backsubst_time; /**< @brief time spent in backsubstitution (solving for L / d) */
	_TyTime m_f_fullL_cholesky; /**< @brief time spent in calculating cholesky (L update) */
				
	size_t m_n_resumed_chol_num; /**< @brief number of times the resumed Cholesky was used */
	size_t m_n_blocks_above_num; /**< @brief number of times there were blocks above lambda_11 */
	size_t m_n_limited_search_num; /**< @brief number of times there were blocks above lambda_11 but only a smaller submatrix was sufficient for permutation calculation */
	_TyTime m_f_ordering_fold_time; /**< @brief time spent folding two orderings */
	_TyTime m_f_repermute_time; /**< @brief time spent repermuting lambda matrix with incremented ordering */
	_TyTime m_f_Lslice_time; /**< @brief time spent slicing L for resumed Cholesky */
	_TyTime m_f_etree_time; /**< @brief time spent calculating the elimination tree */
	_TyTime m_f_resumed_chol_time; /**< @brief time spent in resumed Cholesky */
	_TyTime m_f_ordering11_time; /**< @brief time spent in calculating the incremental ordering */
	_TyTime m_f_ordering11_part_time; /**< @brief time spent in calculating the incremental ordering, only the small or inflated lambda_11 cases */
	_TyTime m_f_ordering11_full_time; /**< @brief time spent in calculating the incremental ordering, only the full lambda_perm cases */

	CTimer m_shared_timer; /**< @brief timer object */

	/**
	 *	@brief wrapper for linear solvers (shields solver capability to solve blockwise)
	 *	@tparam CSolverTag is linear solver tag
	 */
	template <class _CSystem, class _CLinearSolver, class CSolverTag>
	class CLinearSolverWrapper {
	public:
		/**
		 *	@brief estabilishes final block matrix structure before solving iteratively
		 *
		 *	@param[in] r_solver is linear solver
		 *	@param[in] r_lambda is the block matrix
		 *
		 *	@return Always returns true.
		 */
		static inline bool FinalBlockStructure(CLinearSolver &r_solver,
			const CUberBlockMatrix &r_lambda)
		{
			return true;
		}

		/**
		 *	@brief calculates ordering, solves a system
		 *
		 *	@param[in] r_solver is linear solver
		 *	@param[in] r_lambda is the block matrix
		 *	@param[in] r_v_eta is the right side vector
		 *
		 *	@return Returns true on success, false on failure.
		 */
		static inline bool Solve(CLinearSolver &r_solver,
			const CUberBlockMatrix &r_lambda, Eigen::VectorXd &r_v_eta)
		{
			return r_solver.Solve_PosDef(r_lambda, r_v_eta);
		}
	};

	/**
	 *	@brief wrapper for linear solvers (specialization for CBlockwiseLinearSolverTag sovers)
	 */
	template <class _CSystem, class _CLinearSolver>
	class CLinearSolverWrapper<_CSystem, _CLinearSolver, CBlockwiseLinearSolverTag> {
	public:
		/**
		 *	@brief estabilishes final block matrix structure before solving iteratively
		 *
		 *	@param[in] r_solver is linear solver
		 *	@param[in] r_lambda is the block matrix
		 *
		 *	@return Always returns true.
		 */
		static inline bool FinalBlockStructure(CLinearSolver &r_solver,
			const CUberBlockMatrix &UNUSED(r_lambda))
		{
			r_solver.Clear_SymbolicDecomposition();
			// will trigger automatic recalculation and saves one needless converting lambda to cs*

			return true;
		}

		/**
		 *	@brief solves a system, reusing the previously calculated block ordering
		 *
		 *	@param[in] r_solver is linear solver
		 *	@param[in] r_lambda is the block matrix
		 *	@param[in] r_v_eta is the right side vector
		 *
		 *	@return Returns true on success, false on failure.
		 */
		static inline bool Solve(CLinearSolver &r_solver,
			const CUberBlockMatrix &r_lambda, Eigen::VectorXd &r_v_eta)
		{
			return r_solver.Solve_PosDef_Blocky(r_lambda, r_v_eta);
		}
	};

	typedef CLinearSolverWrapper<CSystem, CLinearSolver, _TySolverTag> _TyLinearSolverWrapper; /**< @brief wrapper for linear solvers (shields solver capability to solve blockwise) */

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
	 *	@param[in] linear_solver is linear solver instance for solving for lambda
	 *	@param[in] linear_solver2 is linear solver instance for updating L
	 */
	CNonlinearSolver_FastL(CSystem &r_system, size_t n_linear_solve_threshold,
		size_t n_nonlinear_solve_threshold, size_t n_nonlinear_solve_max_iteration_num,
		double f_nonlinear_solve_error_threshold, bool b_verbose,
		CLinearSolver linear_solver = CLinearSolver(),
		CLinearSolver linear_solver2 = CLinearSolver())
		:m_r_system(r_system), m_linear_solver(linear_solver),
		m_linear_solver2(linear_solver2), m_b_had_loop_closure(false),
		m_b_first_iteration_use_L(true), m_b_L_up_to_date(true), m_n_last_full_L_update_size(0),
		//m_n_big_loop_threshold((n_dummy_param > 0)? n_dummy_param : SIZE_MAX),
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
		m_b_linearization_dirty(false), m_n_iteration_num(0),
		m_f_chol_time(0), m_f_norm_time(0), m_f_vert_upd_time(0),
		m_n_full_forwardsubst_num(0), m_n_resumed_forwardsubst_num(0),
		m_n_resumed_perm_forwardsubst_num(0), m_n_L_optim_num(0), m_n_lambda_optim_num(0),
		m_n_Lup_num(0), m_n_omega_update_num(0), m_n_lambda_update_num(0), m_n_full_L_num(0),
		m_f_lambda_refresh_time(0), m_f_rhs_time(0), m_f_ordering_time(0), m_f_fullL_d(0),
		m_f_l11_omega_calc_time(0), m_f_l11_omega_slice_time(0), m_f_l11_omega_ata_time(0),
		m_f_l11_omega_add_time(0), m_f_l11_lambda_slice_time(0), m_f_l11_lambda_ata_time(0),
		m_f_l11_lambda_add_time(0), m_f_lupdate_time(0), m_f_d_time(0),
		m_f_backsubst_time(0), m_f_fullL_cholesky(0),
		m_n_resumed_chol_num(0), m_n_blocks_above_num(0), m_n_limited_search_num(0),
		m_f_ordering_fold_time(0), m_f_repermute_time(0), m_f_Lslice_time(0), m_f_etree_time(0),
		m_f_resumed_chol_time(0), m_f_ordering11_time(0), m_f_ordering11_part_time(0),
		m_f_ordering11_full_time(0)
	{
		_ASSERTE(!m_n_nonlinear_solve_threshold || !m_n_linear_solve_threshold); // only one of those

		/*printf("hardcoded: ");
		__fbs_ut::CDumpBlockMatrixTypelist<_TyAMatrixBlockSizes>::Print();
		printf("from system: ");
		__fbs_ut::CDumpBlockMatrixTypelist<typename _TySystem::_TyJacobianMatrixBlockList>::Print();*/
		// debug (now these are the same)

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
		m_n_loop_size_cumsum = 0;
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
	}

	/**
	 *	@brief destructor (only required if ? is defined)
	 */
	inline ~CNonlinearSolver_FastL()
	{
	}

	/**
	 *	@brief displays performance info on stdout
	 *	@param[in] f_total_time is total time taken by everything (can be -1 to ommit)
	 */
	void Dump(double f_total_time = -1) const
	{
		printf("solver took " PRIsize " iterations\n", m_n_iteration_num); // debug, to be able to say we didn't botch it numerically
#ifdef __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
		double f_serial_time = m_f_backsubst_time + m_f_chol_time + m_f_norm_time + m_f_vert_upd_time;
		if(f_total_time > 0) {
			printf("solver spent %f seconds in parallelizable section (updating L)\n",
				f_total_time - f_serial_time);
		}
		double f_total_resumed_up = m_f_resumed_chol_time + m_f_Lslice_time +
			m_f_etree_time + m_f_ordering_fold_time + m_f_repermute_time;
		double f_total_omega_up = m_f_l11_omega_calc_time + m_f_l11_omega_slice_time +
			m_f_l11_omega_ata_time + m_f_l11_omega_add_time;
		double f_total_lambda_up = m_f_l11_lambda_slice_time +
			m_f_l11_lambda_ata_time + m_f_l11_lambda_add_time;
		double f_l_upd_time = f_total_resumed_up + f_total_omega_up +
			f_total_lambda_up + m_f_lupdate_time + m_f_ordering11_time +
			m_f_ordering11_part_time + m_f_ordering11_full_time;
		double f_measured_parallel_time = m_f_lambda_refresh_time + m_f_rhs_time + m_f_ordering_time +
			m_f_fullL_d + m_f_fullL_cholesky + f_l_upd_time + m_f_d_time;
		printf("measured parallel time: %f, disparity: %f; out of which:\n", f_measured_parallel_time,
			f_total_time - f_serial_time - f_measured_parallel_time);
		printf("\t   ,\\: %f\n", m_f_lambda_refresh_time);
		printf("\t  rhs: %f\n", m_f_rhs_time);
		printf("\torder: %f\n", m_f_ordering_time);
		printf("\tfullL: %f (ran " PRIsize " times)\n", m_f_fullL_d + m_f_fullL_cholesky, m_n_full_L_num);
		printf("\tout of which:\n");
		printf("\t\t chol: %f\n", m_f_fullL_cholesky);
		printf("\t\t    d: %f\n", m_f_fullL_d);
		printf("\tL update: %f (ran " PRIsize " times)\n", f_l_upd_time, m_n_Lup_num);
		printf("\t\tordfu: %f (blocks above " PRIsize " times)\n", m_f_ordering11_full_time, m_n_blocks_above_num);
		printf("\t\tordli: %f (ran " PRIsize " times)\n", m_f_ordering11_part_time, m_n_limited_search_num);
		printf("\t\tordsm: %f (ran " PRIsize " times)\n", m_f_ordering11_time, m_n_Lup_num - m_n_blocks_above_num - m_n_limited_search_num);
		printf("\t\tresum: %f (ran " PRIsize " times)\n", f_total_resumed_up, m_n_resumed_chol_num);
		printf("\t\t\tofold: %f\n", m_f_ordering_fold_time);
		printf("\t\t\trperm: %f\n", m_f_repermute_time);
		printf("\t\t\tL cut: %f\n", m_f_Lslice_time);
		printf("\t\t\tetree: %f\n", m_f_etree_time);
		printf("\t\t\t chol: %f\n", m_f_resumed_chol_time);
		printf("\t\t  add: %f (ran " PRIsize " times)\n", f_total_omega_up + f_total_lambda_up +
			m_f_lupdate_time, m_n_Lup_num - m_n_resumed_chol_num);
		printf("\t\t\tomega: %f (ran " PRIsize " times)\n", f_total_omega_up, m_n_omega_update_num);
		printf("\t\t\t\t calc: %f\n", m_f_l11_omega_calc_time);
		printf("\t\t\t\tslice: %f\n", m_f_l11_omega_slice_time);
		printf("\t\t\t\t Lata: %f\n", m_f_l11_omega_ata_time);
		printf("\t\t\t\tL11up: %f\n", m_f_l11_omega_add_time);
		printf("\t\t\t   ,\\: %f (ran " PRIsize " times)\n", f_total_lambda_up, m_n_lambda_update_num);
		printf("\t\t\t\tslice: %f\n", m_f_l11_lambda_slice_time);
		printf("\t\t\t\t Lata: %f\n", m_f_l11_lambda_ata_time);
		printf("\t\t\t\tL11up: %f\n", m_f_l11_lambda_add_time);
		printf("\t\t\t  Lup: %f // cholesky and fill\n", m_f_lupdate_time);
		printf("\t    d: %f (resumed " PRIsize ", p-resumed " PRIsize ", full "
			PRIsize ")\n", m_f_d_time, m_n_resumed_forwardsubst_num,
			m_n_resumed_perm_forwardsubst_num, m_n_full_forwardsubst_num);
		printf("solver spent %f seconds in serial section\n", f_serial_time);
		printf("out of which:\n");
		printf("\t chol: %f (ran " PRIsize " times)\n", m_f_chol_time, m_n_lambda_optim_num);
		printf("\tbksub: %f (ran " PRIsize " times)\n", m_f_backsubst_time, m_n_L_optim_num);
		printf("\t norm: %f\n", m_f_norm_time);
		printf("\tv-upd: %f\n", m_f_vert_upd_time);
		/*printf("in unrelated news, small cholesky ran " PRIsize " times\n", m_n_dense_cholesky_num);
		printf("\t dense: %f\n", m_f_dense_cholesky_time);
		printf("\tsparse: %f\n", m_f_sparse_cholesky_time);*/ // dont want to do it runtime
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
			// optimize but don't allow iterations - just updates lambda, d and L
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
			// optimize but don't allow iterations - just updates lambda, d and L
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
	 *	@brief incremental optimization function
	 *	@param[in] r_last_edge is the last edge that was added to the system
	 *	@note This function throws std::bad_alloc and std::runtime_error (when L is not-pos-def).
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
			FILE *p_fw = fopen("timeSteps_L.txt", (m_n_real_step > 0)? "a" : "w");
			fprintf(p_fw, "" PRIsize ";%f;" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";"
				PRIsize ";" PRIsize "\n", m_n_real_step, m_shared_timer.f_Time(),
				m_n_Lup_num, m_n_full_L_num, m_n_L_optim_num, m_n_lambda_optim_num,
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
					// optimize but don't allow iterations - just updates lambda, d and L
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
		if(b_new_vert) {
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
				m_lambda.Permute_UpperTriangluar_To(lord, p_order,
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
				_blockL.Permute_UpperTriangluar_To(blockL, p_order, m_lambda.n_BlockColumn_Num());
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
				m_lambda.Permute_UpperTriangluar_To(lord, p_order,
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
				_blockL.Permute_UpperTriangluar_To(blockL, p_order, m_lambda.n_BlockColumn_Num());
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

			size_t n_nnz_actual = m_L.n_NonZero_Num();
			// actual NNZ

			FILE *p_fw = fopen("LDensityByOrdering.txt", (m_n_real_step > 0)? "a" : "w");
			if(!m_n_real_step) {
				fprintf(p_fw, "block-cols;amd;blocky-amd;blocky-constrained-amd;"
					"amd-blocks;blocky-amd-blocks;blocky-constrained-amd-blocks;actual-L-nnz-blocks\n");
			}
			fprintf(p_fw, "" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize ";" PRIsize "\n", m_lambda.n_BlockColumn_Num(), n_nnz_ideal_elem,
				n_nnz_blocky_elem, n_nnz_blocky_constr_elem, n_nnz_ideal, n_nnz_blocky,
				n_nnz_blocky_constr, n_nnz_actual);
			fclose(p_fw);
		}
		// dump nnz of L, given different ordering strategies
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
	 *	@note This function throws std::bad_alloc and std::runtime_error (when L is not-pos-def).
	 */
	void TryOptimize(size_t n_max_iteration_num, double f_min_dx_norm) // throw(std::bad_alloc, std::runtime_error)
	{
		bool b_optimization_triggered = false;

		size_t n_vertex_num = m_r_system.r_Vertex_Pool().n_Size();
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

		if(!b_optimization_triggered) {
#ifdef __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
			size_t n_vertex_num = m_r_system.r_Vertex_Pool().n_Size();
			if(m_lambda.n_BlockColumn_Num() == n_vertex_num)
				return;
			// there is enough vertices in lambda, none would be added

			Optimize(0, 0); // big todo - remove this in order to be faster for each 100; move it to function that does approx solutions on request
			// optimize but don't allow iterations - just updates lambda, d and L
			// in order to be able to generate approximate solutions on request
#endif // __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
		} else
			Optimize(n_max_iteration_num, f_min_dx_norm);
	}

	/**
	 *	@brief refreshes system matrices lambda and L
	 *	@return Returns true if optimization should take place, otherwise returns false
	 */
	bool RefreshLambdaL()
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

		Extend_LambdaL(m_n_verts_in_lambda, m_n_edges_in_lambda);
		// recalculated all the jacobians inside Extend_LambdaL(), also extend L structurally

		{
			//_TyTimeSampler timer1(m_shared_timer);

			if(!m_b_system_dirty)
				Refresh_Lambda(0/*m_n_verts_in_lambda*/, m_n_edges_in_lambda); // calculate only for new edges // @todo - but how to mark affected vertices?
			else
				Refresh_Lambda(); // calculate for entire system, rebuild L from scratch

			/*double f_refresh_time = 0;
			timer1.Accum_DiffSample(f_refresh_time);
			printf("lambda block width: " PRIsize ", lambda block num: " PRIsize " "
				"(sparsity %.2f %%); update size (edges): " PRIsize " (took %.2f msec)\n",
				m_lambda.n_BlockColumn_Num(), m_lambda.n_Block_Num(),
				float(m_lambda.n_Block_Num() * 9) / (m_lambda.n_Row_Num() *
				m_lambda.n_Column_Num()) * 100, m_lambda.n_BlockColumn_Num() -
				m_n_verts_in_lambda, f_refresh_time * 1000);*/ // debug
		}
		// refresh lambda (does not fully refresh permutated lambda, even though it is a reference matrix)

		m_v_dx.resize(n_variables_size);
		m_v_perm_temp.resize(n_variables_size);
		if(m_b_L_up_to_date && !m_b_system_dirty)
			m_v_d.conservativeResize(n_variables_size); // b_Refresh_L() also refreshes d (rhs), needs dx as temp // !!
		else
			m_v_d.resize(n_variables_size); // in case we're about to rebuild L from scratch, don't care about contents of d
		// resize the helper vectors

		timer.Accum_DiffSample(m_f_lambda_refresh_time);

		{
			if(m_b_L_up_to_date && // can increment L only if up to date
			   !m_b_system_dirty) // avoidance of big incremental updates of L is inside b_Refresh_L() - can only decide if ordering is known
				m_b_first_iteration_use_L = b_Refresh_L(0/*m_n_verts_in_lambda*/, m_n_edges_in_lambda); // calculate only for new edges // @todo - but how to mark affected vertices?
			else
				m_b_first_iteration_use_L = b_Refresh_L(); // calculate for entire system, rebuild L from scratch

			m_b_L_up_to_date = m_b_first_iteration_use_L;
			// in case L is not used, it will fall behind
		}
		// refresh L

		m_b_system_dirty = false;
		m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
		m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
		_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
			m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
			m_lambda.n_Column_Num() == n_variables_size);
		_ASSERTE(!m_b_L_up_to_date || (m_L.n_Row_Num() == m_L.n_Column_Num() &&
			m_L.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
			m_L.n_Column_Num() == n_variables_size)); // lambda is square, blocks on either side = number of vertices
		// need to have lambda and perhaps also L

		return true;
	}

	/**
	 *	@brief updates the m_v_dx vector from the current L and d (or lambda eta)
	 *	@param[in] n_ignore_vertices is number of vertices at the end of the system to be ignored
	 *	@return Returns true on success, false on failure (numerical issues).
	 */
	bool CalculateOneTimeDx(size_t n_ignore_vertices = 0)
	{
		_ASSERTE(m_b_linearization_dirty); // this should only be called in case the linearization point was not updated

		if(m_b_L_up_to_date && m_b_first_iteration_use_L) { // Optimize() clears m_b_L_up_to_date but not m_b_first_iteration_use_L at the same time
			_ASSERTE(m_b_L_up_to_date);
			// we have L and can use it efficiently

			{
				bool b_cholesky_result;
				{
					_ASSERTE(m_p_lambda_block_ordering);
					m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
						m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
					b_cholesky_result = m_L.UpperTriangular_Solve_FBS<_TyLambdaMatrixBlockSizes>(&m_v_perm_temp(0), m_v_perm_temp.rows());
					m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_dx(0), &m_v_perm_temp(0), m_v_dx.rows(),
						m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
					// dx = L'/d // note this never fails (except if L is null)
				}
				// L solves with permutation (note that m_v_d is not modified!)
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
#if 1
				{
					Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
					b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda, v_eta); // p_dx = eta = lambda / eta
					// dont reuse block ordering
				}
				// lambda is good without permutation (there is one inside and we save copying eta arround)
#else // 1
				{
					_ASSERTE(m_p_lambda_block_ordering);
					m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_dx(0), m_v_dx.rows(),
						m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
					// permute the dx vector to eta (stored in m_v_perm_temp)

					b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda_perm, m_v_perm_temp); // p_dx = eta = lambda / eta
					// solve the permutated lambda (dont reuse block ordering)

					m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_dx(0), &m_v_perm_temp(0), m_v_dx.rows(),
						m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
					// permute the result in eta back to bx
				}
#endif // 1
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
	 *	@note This function throws std::bad_alloc and std::runtime_error (when L is not-pos-def).
	 */
	void Optimize(size_t n_max_iteration_num, double f_min_dx_norm) // throw(std::bad_alloc, std::runtime_error)
	{
		if(!RefreshLambdaL())
			return;
		// decide whether to optimize or not

		if(!n_max_iteration_num) {
			m_b_linearization_dirty = true;
			return;
		}
		// in case we're not required to optimize, do nothing
		// (the user can still request solution, L is in good shape)

		if(m_b_had_loop_closure)
			m_b_had_loop_closure = false;
		else {
			m_b_linearization_dirty = true;
			return; // nothing to optimize, dx would be zero
		}
		// handle loop closures a bit differently

#if 0
		static bool b_had_lambda_up = false;
		if(m_b_first_iteration_use_L && b_had_lambda_up) {
			b_had_lambda_up = false;
			Check_LLambdaTracking(); // seems to be working now
		}
#endif // 0
		// make sure lambda and L contain the same system

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

			if(m_b_L_up_to_date/*m_b_first_iteration_use_L && !n_iteration*/) { // always when we have L
				++ m_n_L_optim_num;

				_ASSERTE(m_b_L_up_to_date);
				// we have L and can use it efficiently

				_TyTimeSampler timer(m_shared_timer);

				{
					bool b_utsolve_result;
					{
						_ASSERTE(m_p_lambda_block_ordering);
						m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
							m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
						b_utsolve_result = m_L.UpperTriangular_Solve_FBS<_TyLambdaMatrixBlockSizes>(&m_v_perm_temp(0), m_v_perm_temp.rows());
						m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_dx(0), &m_v_perm_temp(0), m_v_dx.rows(),
							m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
						// dx = L'/d // note this never fails (except if L is null)

						if(m_b_verbose) {
							printf("%s", (b_utsolve_result)? "backsubstitution rulez!\n" :
								"backsubstitution failed!\n");
						}
					}
					// L solves with permutation (note that m_v_d is not modified!)

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
						/*printf("just optimized using L\n");*/
						PushValuesInGraphSystem(m_v_dx); // note this kills L
						m_b_system_dirty = true;
						m_b_L_up_to_date = false; // !!

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
					Refresh_Lambda(); // want only lambda, leave L behind
					m_b_system_dirty = false;
					m_b_L_up_to_date = false; // lambda not dirty anymore, but L still is

					timer.Accum_DiffSample(m_f_lambda_refresh_time);

#ifdef __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
					m_b_L_up_to_date = b_Refresh_L(0, 0, n_iteration > 1); // refresh L as well
					// suppress reordering in iterations 2 and above
					// (already reordered in iteration 1, won't get any better)
					// note that RHS vector is updated inside

					_TyTime f_dummy_sample = 0;
					timer.Accum_DiffSample(f_dummy_sample); // b_Refresh_L() contains timing inside

					if(m_b_L_up_to_date) {
						-- n_iteration;
						-- m_n_iteration_num;
						b_verbose = false; // suppress the banner
						continue;
					}
					// try again, this time with L
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
#if 1
					{
						Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
						if((/*m_b_first_iteration_use_L &&*/ n_max_iteration_num > 2) ||
						   (!m_b_first_iteration_use_L && n_max_iteration_num > 1)) {
							do {
								if(n_iteration == ((m_b_first_iteration_use_L)? 1 : 0) &&
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
							printf("%s", (b_cholesky_result)? "Cholesky rulez!\n" : "optim success: 0\n");
					}
					// lambda is good without permutation (there is one inside and we save copying eta arround)
#else // 1
					{
						_ASSERTE(m_p_lambda_block_ordering);
						m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_dx(0), m_v_dx.rows(),
							m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
						// permute the dx vector to eta (stored in m_v_perm_temp)


						if((/*m_b_first_iteration_use_L &&*/ n_max_iteration_num > 2) ||
						   (!m_b_first_iteration_use_L && n_max_iteration_num > 1)) {
							do {
								if(n_iteration == ((m_b_first_iteration_use_L)? 1 : 0) &&
								   !_TyLinearSolverWrapper::FinalBlockStructure(m_linear_solver, m_lambda_perm)) {
									b_cholesky_result = false;
									break;
								}
								// prepare symbolic decomposition, structure of lambda won't change in the next steps

								b_cholesky_result = _TyLinearSolverWrapper::Solve(m_linear_solver, m_lambda_perm, m_v_perm_temp);
								// p_dx = eta = lambda / eta
							} while(0);
						} else
							b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda_perm, m_v_perm_temp); // p_dx = eta = lambda / eta
						// solve the permutated lambda

						if(m_b_verbose)
							printf("%s", (b_cholesky_result)? "perv Cholesky rulez!\n" : "perm optim success: 0\n");

						m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_dx(0), &m_v_perm_temp(0), m_v_dx.rows(),
							m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
						// permute the result in eta back to bx
					}
					// todo - this branch can be improved by using native cholesky directly on the matrix
					// (and therefore reusing the ordering without calculating a new one)
#endif // 1
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
						m_b_L_up_to_date = false;

#if 0
						b_had_lambda_up = true; // debug
#endif // 0

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
		_ASSERTE(m_b_system_dirty || m_b_L_up_to_date);
		// make sure that L is indeed kept up-to-date, unless the solver
		// was stopped by reaching the maximum number of iterations
#endif // __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
	}

protected:
	/**
	 *	@brief function object that calls lambda hessian block allocation for all edges
	 */
	class CAlloc_LambdaLBlocks { // t_odo - L probably only allocates blocks on vertices; retain old version of this functor with lambda only for edges
	protected:
		CUberBlockMatrix &m_r_lambda; /**< @brief reference to the lambda matrix (out) */
		CUberBlockMatrix &m_r_L; /**< @brief reference to the L matrix (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_lambda is reference to the lambda matrix
		 *	@param[in] r_L is reference to the lambda matrix
		 */
		inline CAlloc_LambdaLBlocks(CUberBlockMatrix &r_lambda, CUberBlockMatrix &r_L)
			:m_r_lambda(r_lambda), m_r_L(r_L)
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
			r_vertex_or_edge.Alloc_LBlocks(m_r_L); // t_odo - alloc L blocks as well
		}
	};

	/**
	 *	@brief function object that calls L factor block allocation for all vertices
	 */
	class CAlloc_LBlocks {
	protected:
		CUberBlockMatrix &m_r_L; /**< @brief reference to the L matrix (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_L is reference to the lambda matrix
		 */
		inline CAlloc_LBlocks(CUberBlockMatrix &r_L)
			:m_r_L(r_L)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyVertex is vertex type
		 *	@param[in,out] r_vertex is vertex to have hessian blocks allocated in L
		 *	@note This function throws std::bad_alloc.
		 */
		template <class _TyVertex>
		inline void operator ()(_TyVertex &r_vertex) // throw(std::bad_alloc)
		{
			r_vertex.Alloc_LBlocks(m_r_L); // t_odo - alloc L blocks as well
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

		m_L.Clear();
		/*Eigen::MatrixXd t_uf_L = (r_t_uf.transpose() * r_t_uf).llt().matrixU();
		if(!m_L.Append_Block(t_uf_L, 0, 0))
			throw std::bad_alloc();*/ // UF is kept in lambda, it is propagated to L additively
		// add unary factor to L (actually cholesky(UF^T * UF), but it's the same matrix) // todo - is this right?

		//m_r_system.r_Edge_Pool().For_Each(CAlloc_LambdaLBlocks(m_lambda, m_L)); // no, edges do not add L blocks
		m_r_system.r_Edge_Pool().For_Each(CAlloc_LambdaBlocks(m_lambda));
		m_r_system.r_Vertex_Pool().For_Each(CAlloc_LambdaLBlocks(m_lambda, m_L)); // can stay, there is no ordering to be applied
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
		//throw std::runtime_error("UpdateSparseSystem not implemented for L"); // t_odo

		_ASSERTE(m_lambda.n_Row_Num() > 0 && m_lambda.n_Column_Num() == m_lambda.n_Row_Num()); // make sure lambda is not empty
		//m_r_system.r_Edge_Pool().For_Each(n_skip_edges,
		//	m_r_system.r_Edge_Pool().n_Size(), CAlloc_LambdaLBlocks(m_lambda, m_L)); // no, edges do not add L blocks
		m_r_system.r_Edge_Pool().For_Each(n_skip_edges,
			m_r_system.r_Edge_Pool().n_Size(), CAlloc_LambdaBlocks(m_lambda));
		m_r_system.r_Vertex_Pool().For_Each(n_skip_vertices,
			m_r_system.r_Vertex_Pool().n_Size(), CAlloc_LambdaLBlocks(m_lambda, m_L)); // will not work if ordering is applied (but it mostly isn't, the increments follow identity ordering)
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
	inline void Extend_LambdaL(size_t n_vertices_already_in_lambda, size_t n_edges_already_in_lambda) // throw(std::bad_alloc)
	{
		if(!n_vertices_already_in_lambda && !n_edges_already_in_lambda)
			AddEntriesInSparseSystem(); // works for empty
		else
			UpdateSparseSystem(n_vertices_already_in_lambda, n_edges_already_in_lambda); // does not work for empty
		// create block matrix lambda
	}

#ifdef _DEBUG
	/**
	 *	@brief checks if L == chol(lambda), prints the norm of the difference to stdout
	 */
	void Check_LLambdaTracking() const
	{
		CUberBlockMatrix LtL_upper;
		m_L.PreMultiplyWithSelfTransposeTo(LtL_upper, true);
		//cs *p_L = m_L.p_Convert_to_Sparse();
		cs *p_lam = m_lambda.p_Convert_to_Sparse();
		//cs *p_Lt = cs_transpose(p_L, 1);
		cs *p_LtL = LtL_upper.p_Convert_to_Sparse();//cs_multiply(p_Lt, p_L);
		cs *p_diff = cs_add(p_LtL, p_lam, 1, -1);
		double f_norm = cs_norm(p_diff);
		//cs_spfree(p_L);
		cs_spfree(p_lam);
		//cs_spfree(p_Lt);
		cs_spfree(p_LtL);
		cs_spfree(p_diff);
		// calculate norm (L*L' - lambda)

		printf("L - lambda tracking: %f\n", f_norm);
	}
#endif // _DEBUG

	/**
	 *	@brief calculates the new L11 matrix
	 *
	 *	@param[in] n_order_min is the minimum vertex that changes in L (zero-based index in blocks)
	 *	@param[in] n_order_max is number of column blocks in L (in blocks)
	 *	@param[in] r_L11_new is matrix, containing the new L11, before calculating cholesky of it
	 *
	 *	@note This may modify / damage r_L11_new as it is no longer needed and it is *not* a reference
	 *		to a part of L, in case that enables speed optimizations.
	 *	@note This function throws std::bad_alloc and std::runtime_error (when not-pos-def).
	 */
	void Refresh_L11(size_t n_order_min, size_t n_order_max, CUberBlockMatrix &r_L11_new) // throw(std::bad_alloc, std::runtime_error)
	{
#ifdef __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
		_ASSERTE(r_L11_new.n_Row_Num() == r_L11_new.n_Column_Num());
		if(r_L11_new.n_Column_Num() < 150) {
			if(!r_L11_new.Cholesky_Dense_FBS<_TyLambdaMatrixBlockSizes, 15>()) // 15, not 150 (would yield >1024 template depth)
				throw std::runtime_error("Cholesky_Dense_FBS() failed to increment L");
			// use statically sized matrices up to 30x30, then dynamically allocated matrices up to 150x150

			m_L.From_Matrix(n_order_min, n_order_min, r_L11_new);
			// put r_L11_new to m_L
		} else {
#else // __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
		{
#endif // __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
			if(!b_Have_NativeSolver) {
				CUberBlockMatrix L11_old;
				m_L.SliceTo(L11_old, n_order_min, n_order_max, n_order_min, n_order_max, true); // get L11 as well, need to clear the blocks first
				// todo - make a ClearBlocks() function to a) delete or b) memset(0) blocks in a rectangular area

				L11_old.Scale_FBS_Parallel<_TyLambdaMatrixBlockSizes>(0);
				// clears the data in the update area (it would be better to erase the blocks, but there is the fap) // todo - see indices of the blocks in L and see if these could be efficiently erased right now (while keeping the structure)
				// note that even if we thrash memory taken by (some) L11 blocks,
				// it will be recollected once a full update takes place
			}
			// only have to clear for the old solvers, the native solver does it automatically

			if(!m_linear_solver2.Factorize_PosDef_Blocky(m_L, r_L11_new,
			   m_L_row_lookup_table, n_order_min, n_order_min))
				throw std::runtime_error("Factorize_PosDef_Blocky() failed to increment L");
		}
	}

	/**
	 *	@brief calculates the new L matrix incrementally using lambda or omega
	 *
	 *	@param[in] n_refresh_from_edge is the first edge that changes
	 *	@param[in] n_order_min is the minimum vertex that changes in L (zero-based index in blocks)
	 *
	 *	@note This function throws std::bad_alloc and std::runtime_error (when not-pos-def).
	 */
	inline void Refresh_L_IncL11(size_t n_refresh_from_edge, size_t n_order_min) // throw(std::bad_alloc, std::runtime_error)
	{
		_ASSERTE(m_b_L_up_to_date);
		// make sure L is up to date with lambda and we can actually increment

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
		const bool b_is_identity_perm = b_identity_perm; // make a copy
		//if(!b_Is_PoseOnly_SLAM)
		//	b_identity_perm = false; // don't allow omega upd for VP // it *is* buggy
		// see if the ordering is identity ordering

#ifndef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		if(n_order_max - n_order_min >= 88) // disable this for timing bench
			b_identity_perm = false;
#endif // !__NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		// yet another threshold

		_TyTimeSampler timer(m_shared_timer);

		bool b_omega_available = b_identity_perm;
		//b_identity_perm = true; // i want both branches to run

		CUberBlockMatrix L11TL11;
		if(b_identity_perm) {
			++ m_n_omega_update_num;
			// count them

			CUberBlockMatrix omega, L11;
			size_t n_elem_order_min = m_lambda_perm.n_BlockColumn_Base(n_order_min);
			m_r_system.r_Edge_Pool().For_Each(n_refresh_from_edge, m_r_system.r_Edge_Pool().n_Size(),
				CCalculateOmega(omega, n_elem_order_min));
			omega.CheckIntegrity();

			timer.Accum_DiffSample(m_f_l11_omega_calc_time);

			m_L.SliceTo(L11, n_order_min, n_order_max, n_order_min, n_order_max, true); // row(0 - min) x col(min - max)
			// calculate the omega matrix (ho, ho, ho) and slice L11

			timer.Accum_DiffSample(m_f_l11_omega_slice_time);

			if(n_order_max - n_order_min >= __NONLIVEAR_SOLVER_FAST_L_PARALLEL_MATMULT_THRESH) // big one // t_odo this never runs, the limit for using L is also 100
				L11.PreMultiplyWithSelfTransposeTo_FBS_Parallel<_TyLambdaMatrixBlockSizes>(L11TL11, true); // calculate L11^T * L11 // t_odo - use FBS and maybe also parallel // t_odo - only need lower diagonal, what to do?
			else
				L11.PreMultiplyWithSelfTransposeTo_FBS<_TyLambdaMatrixBlockSizes>(L11TL11, true); // calculate L11^T * L11 // t_odo - use FBS and maybe also parallel // t_odo - only need lower diagonal, what to do?
			// calculate L11TL11

			timer.Accum_DiffSample(m_f_l11_omega_ata_time);

			bool UNUSED(b_result) = omega.AddTo_FBS<_TyLambdaMatrixBlockSizes>(L11TL11); // todo - maybe also parallel
			_ASSERTE(b_result); // if the block order in omega was wrong, this would fail
			// calculate L11TL11_new = L11TL11 + omega
			// note this uses faster addition algorithm

			timer.Accum_DiffSample(m_f_l11_omega_add_time);
		} /*else {
			f_omega_up_calc_end = f_omega_start;
			f_omega_up_slice_end = f_omega_start;
			f_omega_up_ata_end = f_omega_start;
			f_omega_up_end = f_omega_start;
			// omega update not performed - zero time
		}*/ // .. and do nothing :)

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		b_identity_perm = false; // i want both branches to run
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

		CUberBlockMatrix L01TL01;
		if(!b_identity_perm) {
			CUberBlockMatrix lambda11, L01;
			++ m_n_lambda_update_num;
			// count them

			m_L.SliceTo(L01, 0, n_order_min, n_order_min, n_order_max, true); // row(0 - min) x col(min - max)
			m_lambda_perm.SliceTo(lambda11, n_order_min, n_order_max, n_order_min, n_order_max, true);

			timer.Accum_DiffSample(m_f_l11_lambda_slice_time);

			if(n_order_max - n_order_min >= __NONLIVEAR_SOLVER_FAST_L_PARALLEL_MATMULT_THRESH) // big one // t_odo this never runs, the limit for using L is also 100
				L01.PreMultiplyWithSelfTransposeTo_FBS_Parallel<_TyLambdaMatrixBlockSizes>(L01TL01, true); // calculate L01^T * L01 // t_odo - use FBS and maybe also parallel // t_odo - only need lower diagonal, what to do?
			else
				L01.PreMultiplyWithSelfTransposeTo_FBS<_TyLambdaMatrixBlockSizes>(L01TL01, true); // calculate L01^T * L01 // t_odo - use FBS and maybe also parallel // t_odo - only need lower diagonal, what to do?

			timer.Accum_DiffSample(m_f_l11_lambda_ata_time);

			lambda11.AddTo_FBS<_TyLambdaMatrixBlockSizes>(L01TL01, -1, 1); // t_odo - use FBS // todo - maybe also parallel
			// calculates L01TL01 = -L01TL01 + lambda11 (note the "-1, 1" is correct, the opposite way it crashes)
			// note that lambda11 is upper diagonal, as well as L01TL01

			timer.Accum_DiffSample(m_f_l11_lambda_add_time);
		} /*else {
			f_lambda_up_slice_end = f_omega_up_end;
			f_lambda_up_ata_end = f_omega_up_end;
			f_lambda_up_end = f_omega_up_end;
			// lambda update not taken
		}*/ // do nothing :)

		//double f_l11_omega_time = f_omega_up_end - f_omega_start;
		//double f_l11_lambda_time = f_lambda_up_end - f_omega_up_end; // unused

		Refresh_L11(n_order_min, n_order_max, (b_omega_available)? L11TL11 : L01TL01);

		timer.Accum_DiffSample(m_f_lupdate_time);

#if 1
		Refresh_d_IncL11(n_refresh_from_edge, n_order_min); // use the function, do not repeat code, it is ...
		// note that this contains its own timing inside
#else // 1
		{
			if(b_is_identity_perm) {
				++ m_n_resumed_forwardsubst_num;

				_ASSERTE(m_v_d.rows() == m_lambda.n_Column_Num());
				m_r_system.r_Vertex_Pool().For_Each_Parallel(n_order_min,
					m_r_system.r_Vertex_Pool().n_Size(), CCollect_RightHandSide_Vector(m_v_d));
				// collect part of b to the lower part of d (this is inside of Collect_RightHandSide_Vector())

				m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				m_L.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(
					&m_v_perm_temp(0), m_v_perm_temp.rows(), n_order_min);
				m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				// "resumed forward substitution"
			} else {
				++ m_n_full_forwardsubst_num;

				/*Collect_RightHandSide_Vector(m_v_dx);
				// collects the right-hand side vector (eta)

				_ASSERTE(m_p_lambda_block_ordering);
				m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_dx(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				m_L.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(&m_v_perm_temp(0), m_v_perm_temp.rows());
				m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());*/
				// standard forward substitution 

				n_min_vertex = 0;
				// t_odo - figure out how to resume this one

				_ASSERTE(m_v_d.rows() == m_lambda.n_Column_Num());
				m_r_system.r_Vertex_Pool().For_Each_Parallel(n_min_vertex,
					m_r_system.r_Vertex_Pool().n_Size(), CCollect_RightHandSide_Vector(m_v_d));
				// collect part of b to the lower part of d (this is inside of Collect_RightHandSide_Vector())

				m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				m_L.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(
					&m_v_perm_temp(0), m_v_perm_temp.rows(), n_min_vertex);
				m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				// "resumed forward substitution"
			}
			// convert eta to d (d = eta/L)
		}
		// update d incrementally as well
#endif // 1
	}

	/**
	 *	@brief calculates the new right-hand-side vector, does it incrementally were possible
	 *
	 *	@param[in] n_refresh_from_edge is the first edge that changes
	 *	@param[in] n_order_min is the minimum vertex that changes in L (zero-based index in blocks)
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	inline void Refresh_d_IncL11(size_t UNUSED(n_refresh_from_edge), size_t n_order_min) // throw(std::bad_alloc)
	{
		_ASSERTE(m_b_L_up_to_date);
		// make sure L is up to date with lambda and we can actually increment

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
				m_L.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(
					&m_v_perm_temp(0), m_v_perm_temp.rows(), n_order_min);
				m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				// "resumed forward substitution"
			} else {
				++ m_n_resumed_perm_forwardsubst_num; // a different category

				/*Collect_RightHandSide_Vector(m_v_dx);
				// collects the right-hand side vector (eta)

				_ASSERTE(m_p_lambda_block_ordering);
				m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_dx(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
				m_L.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(&m_v_perm_temp(0), m_v_perm_temp.rows());
				m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());*/
				// standard forward substitution

#if 0
				Eigen::VectorXd v_d = m_v_d;
				{
					Eigen::VectorXd v_eta = m_v_d, v_perm_temp = m_v_d; // no need to initialize, just resize (lazy me)
					Collect_RightHandSide_Vector(v_eta);
					// collects the right-hand side vector (eta)

					_ASSERTE(m_p_lambda_block_ordering);
					m_lambda_perm.InversePermute_LeftHandSide_Vector(&v_perm_temp(0), &v_eta(0), v_eta.rows(),
						m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
					m_L.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(&v_perm_temp(0), v_perm_temp.rows());
					m_lambda_perm.Permute_LeftHandSide_Vector(&v_d(0), &v_perm_temp(0), v_d.rows(),
						m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
					// standard forward substitution
				}
				// "ground truth" solution
				// resumed backsubstitution debugging
#endif // 0

				/*Eigen::VectorXd v_d_copy = m_v_d;
				CCollect_RightHandSide_Vector collector(v_d_copy);
				for(size_t i = n_order_min; i < n_order_max; ++ i)
					collector(m_r_system.r_Vertex_Pool()[m_p_lambda_block_ordering[i]]);*/
				// update a copy of d, only partially

				_ASSERTE(m_v_d.rows() == m_lambda.n_Column_Num());
				/*m_r_system.r_Vertex_Pool().For_Each_Parallel(n_min_vertex,
					m_r_system.r_Vertex_Pool().n_Size(), CCollect_RightHandSide_Vector(m_v_d));*/
				// collect part of b to the lower part of d (this is inside of Collect_RightHandSide_Vector())
				// this also works well

				//n_min_vertex = 0; // disables p-resumed backsubst (testing)

				/*{
					CCollect_RightHandSide_Vector collector(m_v_d);
					for(size_t i = n_order_min; i < n_order_max; ++ i)
						collector(m_r_system.r_Vertex_Pool()[m_p_lambda_block_ordering[i]]);
				}*/ // noes. not enough in some cases
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

				//n_min_vertex = 0;
				// todo - figure out how to resume this one

				/*double f_delta = (v_d_copy - m_v_d).norm();
				if(f_delta > 1e-9) {
					fprintf(stderr, "warning: difference between incrementally "
						"collected RHS and full RHS: %f\n", f_delta);
				}*/

#if 0
				printf("min-vertex: " PRIsize " order-min: " PRIsize " order-max: "
					PRIsize "\n", n_min_vertex, n_order_min, n_order_max);
				// resumed backsubstitution debugging
#endif // 0

				// so ... updating vertices n_order_min to n_order_max (a contiguous range?)
				// after permutation, they go to vector[m_p_lambda_block_ordering[n_order_min to n_order_max]]
				// so the whole part of the _permuted_ vector from m_p_lambda_block_ordering[n_order_min] to
				// the end must be updated, and *no more*

				m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_d(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num()); // dest[p ++] = *src ++
				m_L.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(
					&m_v_perm_temp(0), m_v_perm_temp.rows(), n_min_vertex);
				m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
					m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num()); // *dest ++ = src[p ++]
				// "resumed forward substitution"

#if 0
				double f_delta_d = (v_d - m_v_d).norm();
				if(f_delta_d > 1e-9) {
					fprintf(stderr, "warning: difference between d and ground truth: %f\n", f_delta_d);
					printf("refreshed:");
					for(size_t i = n_order_min; i < n_order_max; ++ i)
						printf(" " PRIsize, m_p_lambda_block_ordering[i]);
					printf("\n");
					printf("should have refreshed:");
					const size_t *p_order_inv = m_lambda_ordering.p_Get_Ordering();
					for(size_t i = n_min_vertex; i < n_order_max; ++ i)
						printf(" " PRIsize, m_p_lambda_block_ordering[p_order_inv[i]]);
					printf("\n");
				}
				// resumed backsubstitution debugging
#endif // 0
				// todo - this can be carried to solver L, if needed
			}
			// convert eta to d (d = eta/L)
		}
		// update d incrementally as well

		timer.Accum_DiffSample(m_f_d_time);
	}

	/**
	 *	@brief calculates the new L matrix from scratch as Cholesky of lambda
	 *	@note This function throws std::bad_alloc and std::runtime_error (when not-pos-def).
	 */
	inline void Refresh_L_FullL() // throw(std::bad_alloc, std::runtime_error)
	{
		_TyTimeSampler timer(m_shared_timer);

		m_n_last_full_L_update_size = m_lambda.n_BlockColumn_Num();
		// L will have the same size once updated ...

		if(!b_Have_NativeSolver) {
			if(b_Is_PoseOnly_SLAM) { // this is known at compile-time, should optimize the unused branch away
				m_L.Clear();
				m_r_system.r_Vertex_Pool().For_Each(CAlloc_LBlocks(m_L)); // won't work with VP problems, need to set correct ordering to vertices
			} else {
				//m_L.Clear(); // already inside PermuteTo()
				CUberBlockMatrix t_new_L;
				m_r_system.r_Vertex_Pool().For_Each(CAlloc_LBlocks(t_new_L));
				t_new_L.PermuteTo(m_L, m_p_lambda_block_ordering, m_n_lambda_block_ordering_size);
			}
			// only need to do there things if not using the native solver
		}
		// do the right thing and thrash L

#if 0
		m_L_row_lookup_table.clear();
		CUberBlockMatrix t_new_L, t_L, t_L2;
		m_r_system.r_Vertex_Pool().For_Each(CAlloc_LBlocks(t_new_L));
		t_new_L.PermuteTo(t_L, m_p_lambda_block_ordering, m_n_lambda_block_ordering_size);
		CLinearSolver_CSparse backup_solver;
		if(!backup_solver.Factorize_PosDef_Blocky(t_L, m_lambda_perm, m_L_row_lookup_table, 0, 0))
			throw std::runtime_error("Factorize_PosDef_Blocky() failed to calculate full L");

		if(!m_linear_solver2.Factorize_PosDef_Blocky(t_L2, m_lambda_perm, m_L_row_lookup_table, 0, 0))
			throw std::runtime_error("Factorize_PosDef_Blocky() failed to calculate full L");
		// factorize (uses cached cholesky, saves some time on allocation of workspace memory)

		if(!t_L.b_EqualLayout(t_L2))
			throw std::runtime_error("blocky L has different layout");
		m_L.Swap(t_L);
#else
		if(!m_linear_solver2.Factorize_PosDef_Blocky(m_L, m_lambda_perm, m_L_row_lookup_table, 0, 0))
			throw std::runtime_error("Factorize_PosDef_Blocky() failed to calculate full L");
		// factorize (uses cached cholesky, saves some time on allocation of workspace memory)
#endif
		// just checking, not it

		timer.Accum_DiffSample(m_f_fullL_cholesky);

		{
			Collect_RightHandSide_Vector(m_v_dx);
			// collects the right-hand side vector (eta)

			timer.Accum_DiffSample(m_f_rhs_time);

			//++ m_n_full_forwardsubst_num;
			// do not count it here, we know how many times we did full L, it is the same count

			_ASSERTE(m_p_lambda_block_ordering);
			m_lambda_perm.InversePermute_LeftHandSide_Vector(&m_v_perm_temp(0), &m_v_dx(0), m_v_dx.rows(),
				m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
			m_L.UpperTriangularTranspose_Solve_FBS<_TyLambdaMatrixBlockSizes>(&m_v_perm_temp(0), m_v_perm_temp.rows());
			m_lambda_perm.Permute_LeftHandSide_Vector(&m_v_d(0), &m_v_perm_temp(0), m_v_dx.rows(),
				m_p_lambda_block_ordering, m_lambda_perm.n_BlockRow_Num());
			// d = eta = eta/L
		}
		// convert eta to d

		timer.Accum_DiffSample(m_f_fullL_d);
	}

	std::vector<size_t> lambda_perm_frontline;

	/**
	 *	@brief refreshes the L matrix either from (pert of) lambda or from omega
	 *
	 *	@param[in] n_referesh_from_vertex is zero-based index of the first vertex that changes (unused)
	 *	@param[in] n_refresh_from_edge is zero-based index of the first edge that changes
	 *	@param[in] b_supress_reorder is 
	 *
	 *	@return Returns true if L was refreshed, false if it decided to take lambda fallback instead.
	 *
	 *	@note This function throws std::bad_alloc and std::runtime_error (when not-pos-def).
	 */
	inline bool b_Refresh_L(size_t UNUSED(n_referesh_from_vertex) = 0,
		size_t n_refresh_from_edge = 0, bool b_supress_reorder = false) // throw(std::bad_alloc, std::runtime_error)
	{
		_TyTimeSampler timer(m_shared_timer);

		// note that lambda is now up to date

		bool b_force_reorder = !n_refresh_from_edge && !b_supress_reorder; // if rebuilding whole lambda, it would be shame not to reorder
		// flag for forcing reorder

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		if(false) { // no optimizations for L up variants timing
#else // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		if(!b_supress_reorder) { // if allow optimizations ...
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			/*if(!b_force_reorder) {
				size_t n_last = m_r_system.r_Vertex_Pool().n_Size() - 1;
				size_t n_order_min_projection = (n_last < m_n_lambda_block_ordering_size)?
					m_p_lambda_block_ordering[n_last] : n_last;
				for(size_t i = n_refresh_from_edge, n = m_r_system.r_Edge_Pool().n_Size(); i < n; ++ i) {
					size_t n_order_v0 = m_r_system.r_Edge_Pool()[i].n_Vertex_Id(0);
					size_t n_order_v1 = m_r_system.r_Edge_Pool()[i].n_Vertex_Id(1); // note that these are ids, but these equal order at the moment
					n_order_v0 = (n_order_v0 < m_n_lambda_block_ordering_size)?
						m_p_lambda_block_ordering[n_order_v0] : n_order_v0;
					n_order_v1 = (n_order_v1 < m_n_lambda_block_ordering_size)?
						m_p_lambda_block_ordering[n_order_v1] : n_order_v1;
					n_order_min_projection = std::min(n_order_min_projection, std::min(n_order_v0, n_order_v1));
				}
				size_t n_order_max = m_lambda.n_BlockColumn_Num();
				// project where will n_order_min end up after ordering is updated

				size_t n_projected_update_size = n_order_max - n_order_min_projection;
				if(n_projected_update_size > m_n_big_loop_threshold) // todo - remove m_n_big_loop_threshold
					b_force_reorder = true;
			}*/
			// calculates how big would the next edge loop be, assuming that the permutation
			// is just extended with identity permutation and forces full L update if it was large

			if(!b_force_reorder && m_lambda.n_BlockColumn_Num() > m_n_last_full_L_update_size + 10) { // it's always dense at the beginning // hlamfb 10
				size_t n_nnz = m_L.n_Storage_Size();
				float f_area = float(m_L.n_Column_Num()) * m_L.n_Column_Num();
				float f_density = n_nnz / f_area;
				if(f_density > 0.02f) {
					b_force_reorder = true;
					//printf("L became too dense (%.2f %%), forcing reorder\n", f_density * 100); // verbose
				}
			}
			// permit 2% density in L, then rebuild
		}
		// these two should just set b_force_reorder
		
		/*if(!b_force_reorder) {
			m_lambda_perm.Get_UpperTriangular_BlockFrontline(lambda_perm_frontline);
			// get the frontline

			// now: a) extend the reorder area to contain all the blocks in frontline (so lambda10 is null)
			//      b) freezee the few blocks in the frontline using constrained ordering? will that help?
			//      c) ?

			const size_t n_order_max = m_lambda.n_BlockColumn_Num();
			
			size_t n_order_min = m_p_lambda_block_ordering[m_n_lambda_block_ordering_size - 1];
			for(size_t i = n_refresh_from_edge, n = m_r_system.r_Edge_Pool().n_Size(); i < n; ++ i) {
				size_t n_order_v0 = m_p_lambda_block_ordering[m_r_system.r_Edge_Pool()[i].n_Vertex_Id(0)];
				size_t n_order_v1 = m_p_lambda_block_ordering[m_r_system.r_Edge_Pool()[i].n_Vertex_Id(1)]; // note that these are ids, but these equal order at the moment
				n_order_min = std::min(n_order_min, std::min(n_order_v0, n_order_v1));
			}

			bool b_limited_search_region = false;
			size_t n_blocks_above = 0;
			_ASSERTE(lambda_perm_frontline.size() == n_order_max);
			for(size_t i = n_order_min + 1; i < n_order_max - 1; ++ i) {
				if(lambda_perm_frontline[i] < n_order_min) {
					++ n_blocks_above;
					/*b_blocks_above = true;
					break;* /
				}
			}
			if(n_blocks_above) {
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
						if(n_frontline_min > lambda_perm_frontline[i])
							n_frontline_min = lambda_perm_frontline[i];
					}
					// see if we can expand the ordering submatrix a little, to avoid using the whole matrix
				}

				if(n_context_min <= n_order_max / 8)
					b_force_reorder = true;
			}
		}*/
		// force reordering instead of doing heavily constrained ordering (slow)
		// not a good idea! it is spending much more time doing the full cholesky
		// maybe there would be a solution in doing some other reorder (larger reorder),
		// and constraining only a part of the matrix (a larger part, still). that would
		// give the advantage of reorder while still not doing full cholesky (but how
		// to calculate the desired size of reorder?)

		if(b_force_reorder) {
			// only calculate a new ordering on full refresh or if forced

			//printf("build new ordering ...\n");

			/*m_p_lambda_block_ordering = m_lambda_ordering.p_InvertOrdering(
				m_lambda_ordering.p_HybridBlockOrdering(m_lambda, m_r_system.r_Edge_Pool().n_Size(),
				m_lambda_constraint.p_Get(m_lambda.n_BlockColumn_Num()), m_lambda.n_BlockColumn_Num()), m_lambda.n_BlockColumn_Num());*/ // uses hybrid aat calculation
			m_lambda_ordering.p_BlockOrdering(m_lambda,
				m_lambda_constraint.p_Get(m_lambda.n_BlockColumn_Num()),
				m_lambda.n_BlockColumn_Num(), true); // constrained blocky, calculate inverse as well
			m_p_lambda_block_ordering = m_lambda_ordering.p_Get_InverseOrdering();
				//m_lambda_ordering.p_BlockOrdering(m_lambda), m_lambda.n_BlockColumn_Num()); // todo - make sure that the last vertex is a pose (otherwise we need to modify the constraint to select the pose, not the landmark)
			// get blockwise and elementwise ordering ...

			if(!b_Have_NativeSolver)
				m_L_row_lookup_table.clear(); // unused in native solver
			// can't reuse lookup of L's rows since these change with ordering
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

		m_lambda.Permute_UpperTriangluar_To(m_lambda_perm, m_p_lambda_block_ordering,
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
			// if !n_order_min, L01 is empty and L merely equals chol(lambda)

			if(m_n_edges_in_lambda == m_r_system.r_Edge_Pool().n_Size()) {
				_ASSERTE(m_n_verts_in_lambda == m_lambda.n_BlockColumn_Num());
				_ASSERTE(m_n_verts_in_lambda == m_r_system.r_Vertex_Pool().n_Size());
				return true;
			}
			// this is final optimization, no need to refresh, there is no new edge / vertex

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			double p_inc_upd_times_start[] = {
				m_f_lupdate_time,
				m_f_d_time
			};
			double p_omega_upd_times_start[] = {
				m_f_l11_omega_calc_time,
				m_f_l11_omega_slice_time,
				m_f_l11_omega_ata_time,
				m_f_l11_omega_add_time
			};
			double p_lambda_upd_times_start[] = {
				m_f_l11_lambda_slice_time,
				m_f_l11_lambda_ata_time,
				m_f_l11_lambda_add_time
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

			//size_t n_blocks_above = 0; // unused, not sure why i thought i needed it
			_ASSERTE(lambda_perm_frontline.size() == n_order_max);
			for(size_t i = n_order_min + 1; i < n_order_max /*- 1*/; ++ i) { // not sure why dont i check the last one - todo
				if(lambda_perm_frontline[i] < n_order_min) {
					//++ n_blocks_above;
					b_blocks_above = true;
					break;
				}
			}
			//if(n_blocks_above)
			//	b_blocks_above = true;
			// see if there are blocks above (except in the first or last column which are fixed)

			/*size_t n_frontline_min = m_lambda_perm.n_Get_UpperTriangular_BlockFrontline_Minimum(
				n_order_min + 1, n_order_max);
			b_blocks_above = n_frontline_min < n_order_min;*/
			// do the same without representing the frontline (works worse, that gotta be the last one ...)

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
				//m_p_lambda11_block_ordering = m_lambda11_ordering.p_InvertOrdering(
				//	m_lambda11_ordering.p_BlockOrdering(lambda11,
				//	m_lambda11_constraint.p_Get(lambda11.n_BlockColumn_Num()),
				//	lambda11.n_BlockColumn_Num()), lambda11.n_BlockColumn_Num()); // inverse (unused at the moment)
				m_n_lambda_block11_ordering_size = lambda11.n_BlockColumn_Num();
				//const size_t *p_lam11_inv = m_lambda11_ordering.p_InvertOrdering(m_p_lambda11_block_ordering,
				//	m_n_lambda_block11_ordering_size); // the order of the update // no need for this one
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

#if 0
				/*fprintf(stderr, "size: " PRIsize ", order-min: " PRIsize ", context-min: " PRIsize "\n",
					n_order_max, n_order_min, n_context_min);*/
				if(n_context_min < 10 || (n_order_max > 6000 && n_context_min <= n_order_max / 8)) { // in case we will do full
					char p_s_filename[256];
					sprintf(p_s_filename, /*"/media/Data/rssdumps/"*/"size%06" _PRIsize "_order-min" PRIsize "_context-min" PRIsize ".tga",
						n_order_max, n_order_min, n_context_min);
					//printf("rast to \'%s\'\n", p_s_filename);
					const int n_ss = 3;
					TBmp *p_img = m_lambda_perm.p_Rasterize(0, n_ss);
					if(p_img) {
						//printf("drawing\n");
						int n_block_size = _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime;
						int n_line0 = (n_ss - 1) * n_block_size * n_order_min;
						int n_line1 = (n_ss - 1) * n_block_size * n_context_min;
						p_img->DrawLine(n_line0, 0, n_line0, p_img->n_height, 0xffff0000U, 3);
						p_img->DrawLine(0, n_line0, p_img->n_width, n_line0, 0xffff0000U, 3); // draw line that separates order min
						p_img->DrawLine(n_line1, 0, n_line1, p_img->n_height, 0xff00ff00U, 3);
						p_img->DrawLine(0, n_line1, p_img->n_width, n_line1, 0xff00ff00U, 3); // draw line that separates context min
						//printf("saving\n");
						CTgaCodec::Save_TGA(p_s_filename, *p_img, false, true);
						//printf("deleting\n");
						p_img->Delete();
					} else
						fprintf(stderr, "error: not enough memory to rasterize the matrix\n");
				}
				// debug - dump the matrices
#endif // 0

				if(n_context_min > n_order_max / 8 /*std::min(size_t(100), n_order_max / 16)*/) { // todo - need to dump n_context_min, probably citytrees have one at 0 or something, that can not be ordered away (make it go away by forcing it somewhere else?l)
					// this is prefix-constrained ordering on part of lambda perm (not full, not update)
					// this works rather well

					b_limited_search_region = true;
					// say we did it

					/*printf("optimized order on " PRIsize " x " PRIsize " instead of " PRIsize " x " PRIsize
						" (update is " PRIsize " x " PRIsize ")\n", n_order_max - n_context_min,
						n_order_max - n_context_min, n_order_max, n_order_max, n_order_max - n_order_min,
						n_order_max - n_order_min);*/
					// debug

#if 1
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
#else // 1
#ifndef __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
					CUberBlockMatrix lambda11; // local is enough
#endif // !__NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
					m_lambda_perm.SliceTo(lambda11, n_context_min,
						n_order_max, n_context_min, n_order_max, true);
					_ASSERTE(lambda11.b_Square() && lambda11.b_SymmetricLayout());
					// slice bigger part of lambda with a sufficient context

					_ASSERTE(n_order_min > n_context_min);
					size_t n_skip_order = n_order_min - n_context_min;
					const size_t *p_order = m_lambda_alt_ordering.p_BlockOrdering(lambda11,
						m_lambda_alt_constraint.p_Get(lambda11.n_BlockColumn_Num(),
						n_skip_order), lambda11.n_BlockColumn_Num());
					// get ordering on the submatrix of permuted lambda
#endif // 1

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
					//m_p_lambda11_block_ordering = m_lambda11_ordering.p_InvertOrdering(p_order11,
					//	m_n_lambda_block11_ordering_size); // inverse
					// invert it for merging with the original ordering
				} else {
					// this is prefix-constrained ordering on the full lambda (perm)
					// this works rather well

					/*const size_t *p_order = m_lambda_alt_ordering.p_HybridBlockOrdering(m_lambda_perm,
						m_r_system.r_Edge_Pool().n_Size(), m_lambda_alt_constraint.p_Get(
						m_lambda.n_BlockColumn_Num(), n_order_min), m_lambda.n_BlockColumn_Num());*/
#if 1
					const size_t *p_order = m_lambda_alt_ordering.p_BlockOrdering_MiniSkirt(m_lambda_perm,
						0, n_order_min, m_lambda_alt_constraint.p_Get(m_lambda.n_BlockColumn_Num(),
						n_order_min), m_lambda.n_BlockColumn_Num());
					// just give diagonal matrix all the way from 0 to n_order_min, then the actual
					// sparsity pattern till n_order_max
#else // 1
					const size_t *p_order = m_lambda_alt_ordering.p_BlockOrdering(m_lambda_perm,
						m_lambda_alt_constraint.p_Get(m_lambda.n_BlockColumn_Num(), n_order_min),
						m_lambda.n_BlockColumn_Num());
#endif // 1
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
			lambda11.Permute_UpperTriangluar_To(lambda11_p, m_p_lambda11_block_ordering,
				m_n_lambda_block11_ordering_size, false); // make a deep copy
			// copy lambda 00 and lambda 11 (don't care about lambda 01, it is hard to permute correctly at this point)
#endif // __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING

			_TyTime f_ordering11_time = 0;
			timer.Accum_DiffSample(f_ordering11_time);

			if(!b_identity_ordering) {
				if(!b_Have_NativeSolver)
					m_L_row_lookup_table.clear(); // unused in native solver
				// !! we defined a new ordering

				const size_t *p_order;
				m_p_lambda_block_ordering = m_lambda_ordering.p_InvertOrdering(p_order =
					m_lambda_ordering.p_ExtendBlockOrdering_with_SubOrdering(n_order_min,
					m_p_lambda11_block_ordering, m_n_lambda_block11_ordering_size), m_lambda.n_BlockColumn_Num());
				_ASSERTE(m_n_lambda_block_ordering_size == n_order_min + m_n_lambda_block11_ordering_size);
				// update the ordering (update means append with lambda11 sub-block ordering)
				// this is quick, no bottleneck in here

				timer.Accum_DiffSample(m_f_ordering_fold_time);

				/*printf("order = {");
				for(size_t i = 0; i < m_lambda.n_BlockColumn_Num(); ++ i)
					printf((i)? ", " PRIsize : PRIsize, p_order[i]);
				printf("}, order-min = " PRIsize "\n", n_order_min);*/ // debug

				m_lambda.Permute_UpperTriangluar_To(m_lambda_perm, m_p_lambda_block_ordering,
					m_lambda.n_BlockColumn_Num(), true, n_order_min); // only writes a small square n_order_min - n_order_max // saves 1.45 seconds -> 0.47 seconds
				// make a new reordered version of lambda (*always* changes if lambda changes)
				// this is not very fast, a significant portion of m_lambda_perm could be reused
				// note that the resumed cholesky will not reference anything but lambda11, no need to fill the prev columns
				// need to be careful about use of lambda_perm elsewhere, though (but elsewhere mostly L is used,
				// noone cares about lambda_perm, except maybe for its dimensions)

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

				//CUberBlockMatrix L_copy;
				//m_L.SliceTo(L_copy, 0, n_order_min, 0, n_order_min, false);
				//L_copy.Swap(m_L);
				// throw away L11 and L10
				// this is *very* slow
				//m_L.SliceTo(m_L, n_order_min, n_order_min); // 3 seconds -> 2 seconds, still quite slow
				m_L.SliceTo(m_L, n_order_min, n_order_min, true); // 3 seconds -> 0.3 seconds

				timer.Accum_DiffSample(m_f_Lslice_time);

				if(m_chol_bitfield.capacity() < n_order_max) {
					m_chol_bitfield.clear();
					m_chol_bitfield.reserve(std::max(n_order_max, 2 * m_chol_bitfield.capacity()));
				}
				m_chol_bitfield.resize(n_order_max, 0);

				m_lambda_perm.Build_EliminationTree(m_chol_etree, m_chol_ereach_stack); // use ereach stack as workspace
				_ASSERTE(m_chol_ereach_stack.size() == n_order_max);
				// build an elimination tree

				timer.Accum_DiffSample(m_f_etree_time);

				if(!m_L.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(m_lambda_perm, m_chol_etree,
				   m_chol_ereach_stack, m_chol_bitfield, n_order_min)) { // todo - do incremental etree as well, might save considerable time
					m_lambda.Rasterize("npd0_lambda.tga");
					m_lambda_perm.Rasterize("npd1_lambda_perm.tga");
					m_L.Rasterize("npd2_L.tga");
					fprintf(stderr, "error: got not pos def in incL section\n"); // big todo - throw
				}
				// calcualte updated L11 and L10 using resumed Cholesky

#if 0
				if(m_L.n_BlockColumn_Num() == 1750) {
					{
						cs *p_lbl = m_lambda_perm.p_BlockStructure_to_Sparse();
						cs *p_Lbl = m_L.p_BlockStructure_to_Sparse();
						CDebug::Dump_SparseMatrix("L_at_1750.tga", p_Lbl, p_lbl, 3);
						printf("lambda11 at 1750 is " PRIsize " x " PRIsize " blocks\n",
							m_n_lambda_block11_ordering_size, m_n_lambda_block11_ordering_size);
						cs_spfree(p_Lbl);
						cs_spfree(p_lbl);
					}
					{
						CUberBlockMatrix L11;
						L11.CholeskyOf(lambda11_p);
						cs *p_lbl = lambda11_p.p_BlockStructure_to_Sparse();
						cs *p_Lbl = L11.p_BlockStructure_to_Sparse();
						CDebug::Dump_SparseMatrix("L11_at_1750.tga", p_Lbl, p_lbl, 3);
						printf("lambda11 at 1750 is " PRIsize " x " PRIsize " blocks\n",
							m_n_lambda_block11_ordering_size, m_n_lambda_block11_ordering_size);
						cs_spfree(p_Lbl);
						cs_spfree(p_lbl);
					}
					/*m_p_lambda_block_ordering = m_lambda_ordering.p_InvertOrdering(
						m_lambda_ordering.p_HybridBlockOrdering(m_lambda,
						m_r_system.r_Edge_Pool().n_Size()), m_lambda.n_BlockColumn_Num());*/
					m_lambda_ordering.p_BlockOrdering(m_lambda, true);
					m_p_lambda_block_ordering = m_lambda_ordering.p_Get_InverseOrdering();
					m_lambda.Permute_UpperTriangluar_To(m_lambda_perm, m_p_lambda_block_ordering,
						m_lambda.n_BlockColumn_Num(), true);
					m_L.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(m_lambda_perm);
					{
						cs *p_lbl = m_lambda_perm.p_BlockStructure_to_Sparse();
						cs *p_Lbl = m_L.p_BlockStructure_to_Sparse();
						CDebug::Dump_SparseMatrix("L_at_1750_fullordered.tga", p_Lbl, p_lbl, 3);
						printf("lambda11 at 1750 is " PRIsize " x " PRIsize " blocks\n",
							m_n_lambda_block11_ordering_size, m_n_lambda_block11_ordering_size);
						cs_spfree(p_Lbl);
						cs_spfree(p_lbl);
					}
					exit(-1);
				}
				// debug (manhattan olson 3500, need verify perm folding)
#endif // 0

				++ m_n_resumed_chol_num;
				timer.Accum_DiffSample(m_f_resumed_chol_time);

				Refresh_d_IncL11(n_refresh_from_edge, n_order_min); // timing inside
				// all that is left to do is to refresh d
			} else {
				Refresh_L_IncL11(n_refresh_from_edge, n_order_min); // timing inside
				// run the "fast" refresh of L
			}
			// choose between progressive reordering and "fast" update to L


#ifdef __NONLINEAR_SOLVER_FASTL_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
			char p_s_filename[256];

			CUberBlockMatrix lambda_perm;
			m_lambda.Permute_UpperTriangluar_To(lambda_perm, m_p_lambda_block_ordering,
				m_lambda.n_BlockColumn_Num(), true);
			// need to reperm, may only have a part of lambda_perm, effectively selecting everything above n_order_min as a new nnz

			size_t n_verts_in_lambda = m_lambda.n_BlockColumn_Num();
			sprintf(p_s_filename, "rss2013/%05d_6_lambda-perm.tga", n_verts_in_lambda);
			//if(n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			//	lambda_perm.Rasterize_Symmetric(p_s_filename, (n_verts_in_lambda < 750 * 6 / _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 3); // do not really need lambdas right now

			sprintf(p_s_filename, "rss2013/%05d_7_L.tga", n_verts_in_lambda);
			//if(n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			{
				int n_ss; // scalar size
				TBmp *p_img = m_L.p_Rasterize(lambda_perm, false, 0, n_ss = ((n_verts_in_lambda < 750 * 6 / _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 3 : 2)); // highlight fill-in
				if(p_img) {
					CTgaCodec::Save_TGA(p_s_filename, *p_img, false, true);

					sprintf(p_s_filename, "rss2013/%05d_8_L_marked.tga", n_verts_in_lambda);
					//printf("drawing\n");
					int n_line0 = (n_ss - 1) * m_L.n_BlockColumn_Base(n_order_min);
					//int n_line1 = (n_ss - 1) * m_L.n_BlockColumn_Base(n_context_min);
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
				fprintf(p_fw, PRIsize "\n", m_L.n_Block_Num()); // save density of L
				fprintf(p_fw, PRIsize "\n", n_order_min); // save density of L
				fclose(p_fw);
			}
#endif // __NONLINEAR_SOLVER_FASTL_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

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

			++ m_n_Lup_num;
			// count incremental L updates

#ifdef __NONLINEAR_SOLVER_FAST_L_PROOF_OF_CONCEPT
			double f_proper_end = m_shared_timer.f_Time();

			L_copy.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(m_lambda_perm, n_order_min); // todo - do incremental etree as well, might save considerable time

			double f_resumed_end = m_shared_timer.f_Time();
			f_resumed_end -= f_proper_end;
			f_proper_end -= f_proper_start; // f_proper_start got lost somewhere

			m_L.AddTo(L_copy, -1);
			double f_difference = L_copy.f_Norm();
			printf("L - L_ref = %g (resumed %f , proper %f , speedup %f )\n", // the otherwise completely unacceptable spaces are used to simplify import into excel spreadsheet
				f_difference, f_resumed_end, f_proper_end, f_proper_end / f_resumed_end);
			// test the concept of resumed cholesky
#endif // __NONLINEAR_SOLVER_FAST_L_PROOF_OF_CONCEPT

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			double p_inc_upd_times[] = {
				m_f_lupdate_time,
				m_f_d_time
			};
			double f_inc_upd_sum = 0;
			for(int i = 0; i < 2; ++ i) {
				p_inc_upd_times[i] -= p_inc_upd_times_start[i];
				f_inc_upd_sum += p_inc_upd_times[i];
			}
			double p_omega_upd_times[] = {
				m_f_l11_omega_calc_time,
				m_f_l11_omega_slice_time,
				m_f_l11_omega_ata_time,
				m_f_l11_omega_add_time
			};
			double f_omega_upd_sum = f_inc_upd_sum;
			for(int i = 0; i < 4; ++ i) {
				p_omega_upd_times[i] -= p_omega_upd_times_start[i];
				f_omega_upd_sum += p_omega_upd_times[i];
			}
			double p_lambda_upd_times[] = {
				m_f_l11_lambda_slice_time,
				m_f_l11_lambda_ata_time,
				m_f_l11_lambda_add_time
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

			double f_full_L_start = m_shared_timer.f_Time();
			Refresh_L_FullL();
			double f_full_L_time = m_shared_timer.f_Time() - f_full_L_start;
			// measure full L as well

			FILE *p_fw = fopen("Lup_variants_time.txt", (b_first_time_dump)? "w" : "a");
			if(b_first_time_dump) {
				fprintf(p_fw, "verts-in-L;loop-size;full-L-time;lambda-up-time;lambda-slice-time;"
					"lambda-ata-time;lambda-add-time;omega-up-time;omega-calc-time;"
					"omega-slice-time;omega-ata-time;omega-add-time\n");
			}
			fprintf(p_fw, "" PRIsize ";" PRIsize ";%f;%f;%f;%f;%f", m_L.n_BlockColumn_Num(), n_loop_size, f_full_L_time,
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
			//printf("doing full L\n"); // debug

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			double f_full_L_start = m_shared_timer.f_Time();
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

			++ m_n_full_L_num;
			// L is not up to date, need to rebuild from scratch

			Refresh_L_FullL();
			// do the "full" L = chol(lambda)

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
			double f_full_L_time = m_shared_timer.f_Time() - f_full_L_start;
			// measure full L

			size_t n_loop_size = m_n_lambda_block_ordering_size - 1 - n_order_min;

			FILE *p_fw = fopen("Lup_variants_time.txt", (b_first_time_dump)? "w" : "a");
			if(b_first_time_dump) {
				fprintf(p_fw, "verts-in-L;loop-size;full-L-time;lambda-up-time;lambda-slice-time;"
					"lambda-ata-time;lambda-add-time;omega-up-time;omega-calc-time;"
					"omega-slice-time;omega-ata-time;omega-add-time\n");
			}
			fprintf(p_fw, "" PRIsize ";" PRIsize ";%f;;;;", m_L.n_BlockColumn_Num(), n_loop_size, f_full_L_time); // no lambda upd
			fprintf(p_fw, ";;;;;\n"); // no omega upd
			fclose(p_fw);
			// print timing to a file
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES


#ifdef __NONLINEAR_SOLVER_FASTL_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
			char p_s_filename[256];

			size_t n_verts_in_lambda = m_lambda.n_BlockColumn_Num();
			sprintf(p_s_filename, "rss2013/%05d_6_lambda-perm.tga", n_verts_in_lambda);
			//if(n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			//	lambda_perm.Rasterize_Symmetric(p_s_filename, (n_verts_in_lambda < 750 * 6 / _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 3); // do not really need lambdas right now

			sprintf(p_s_filename, "rss2013/%05d_7_L.tga", n_verts_in_lambda);
			//if(n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			{
				int n_ss; // scalar size
				TBmp *p_img = m_L.p_Rasterize(m_lambda_perm, false, 0, n_ss = ((n_verts_in_lambda < 750 * 6 / _TyLambdaMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 3 : 2)); // highlight fill-in
				if(p_img) {
					CTgaCodec::Save_TGA(p_s_filename, *p_img, false, true);

					sprintf(p_s_filename, "rss2013/%05d_8_L_marked.tga", n_verts_in_lambda);
					//printf("drawing\n");
					int n_line0 = 0;//(n_ss - 1) * m_L.n_BlockColumn_Base(n_order_min);
					//int n_line1 = (n_ss - 1) * m_L.n_BlockColumn_Base(n_context_min);
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
				fprintf(p_fw, PRIsize "\n", m_L.n_Block_Num()); // save density of L
				fclose(p_fw);
			}
#endif // __NONLINEAR_SOLVER_FASTL_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
		}

#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
		b_first_time_dump = false;
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES

#ifdef _DEBUG
		//m_L.Rasterize("lslam7_LAfterOpt.tga");
#endif // _DEBUG

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
