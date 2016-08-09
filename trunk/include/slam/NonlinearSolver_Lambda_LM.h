/*
								+-----------------------------------+
								|                                   |
								| ***  Lambda nonlinear solver  *** |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2012  |
								|                                   |
								|    NonlinearSolver_Lambda_LM.h    |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
#define __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED

/**
 *	@file include/slam/NonlinearSolver_Lambda_LM.h
 *	@brief nonlinear blocky solver working above the lambda matrix, with Levenberg Marquardt
 *	@author -tHE SWINe-
 *	@date 2012-09-13
 */

#include "slam/FlatSystem.h"
#include "slam/LinearSolver_Schur.h"
#include "slam/OrderingMagic.h"
#include "slam/NonlinearSolver_Lambda.h"

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
 *	@brief enables writes of chi2 errors at each step
 */
//#define __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
 *	@brief if defined, the number of iterations is quadrupled and after
 *		every full iteration there are three landmark settle iterations
 */
//#define __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS

/**
 *	@brief trust region heuristic interface
 *	@tparam CLambdaLM_Solver is a specialization of CNonlinearSolver_Lambda_LM
 *	@note This just enforces that all the functions are declated and in correct form.
 */
template <class CLambdaLM_Solver>
class CTrustRegion_Interface {
public:
	/**
	 *	@brief calculates error metric for LM
	 *	@param[in] r_solver is reference to the solver (solver state should be all valid at this point)
	 *	@return Returns error metric for LM.
	 */
	virtual double f_Error(const CLambdaLM_Solver &r_solver) = 0;

	/**
	 *	@brief calculates initial damping for LM
	 *	@param[in] m_r_system is reference to the system (solver state is not completely up to date)
	 *	@return Returns the initial damping for LM.
	 */
	virtual double f_InitialDamping(const typename CLambdaLM_Solver::_TySystem& m_r_system) = 0;

	/**
	 *	@brief decides on stepping forward / backward after updating the system
	 *
	 *	@param[in,out] r_f_last_error is reference to error before the update
	 *	@param[in] f_error is error after the update
	 *	@param[in,out] r_f_alpha is the current damping factor
	 *	@param[in] r_lambda is reference to the system matrix before the update
	 *	@param[in] r_solver is reference to the solver
	 *	@param[in] r_v_dx is reference to the dx vector
	 *	@param[in] r_v_rhs is reference to the right-hand-side vector in this iteration
	 *
	 *	@return Returns true if the step was a success, false if the step needs to be rolled back.
	 */
	virtual bool Aftermath(double &r_f_last_error, double f_error, double &r_f_alpha, const CUberBlockMatrix &r_lambda,
		const CLambdaLM_Solver &r_solver, const Eigen::VectorXd &r_v_dx, const Eigen::VectorXd &r_v_rhs) = 0;

	/**
	 *	@brief adds the damping factor to the system matrix
	 *
	 *	@param[in,out] r_lambda is reference to the system matrix before the update
	 *	@param[in] f_alpha is the current damping factor
	 *	@param[in] n_first_vertex is zero-based index of the first column of lambda to refresh
	 *	@param[in] n_last_vertex is zero-based index of the one past the last column of lambda to refresh
	 */
	virtual void ApplyDamping(CUberBlockMatrix &r_lambda, double f_alpha, size_t n_first_vertex, size_t n_last_vertex) = 0;
};

/**
 *	@brief baseline LM algorithm
 *	@tparam CLambdaLM_Solver is a specialization of CNonlinearSolver_Lambda_LM
 */
template <class CLambdaLM_Solver>
class CLevenbergMarquardt_Baseline : public CTrustRegion_Interface<CLambdaLM_Solver> {
protected:
	double m_nu; /**< @brief damping step value */ // reused between nonlinear solver iterations

public:
	/**
	 *	@copydoc CTrustRegion_Interface::f_Error
	 */
	double f_Error(const CLambdaLM_Solver &r_solver)
	{
		return r_solver.f_Chi_Squared_Error_Denorm();
	}

	/**
	 *	@copydoc CTrustRegion_Interface::f_InitialDamping
	 */
	double f_InitialDamping(const typename CLambdaLM_Solver::_TySystem& m_r_system)
	{
		double f_alpha = 0.0;
		m_nu = 2.0;
		const double tau = 1e-3;
		/*for(size_t i = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_last = m_lambda.n_BlockColumn_Block_Num(i) - 1;
			CUberBlockMatrix::_TyMatrixXdRef block = m_lambda.t_BlockAt(n_last, i);
			for(size_t j = 0, m = block.rows(); j < m; ++ j)
				f_alpha = std::max(block(j, j), f_alpha);
		}*/
		for(size_t i = 0, n = m_r_system.n_Edge_Num(); i < n; ++ i) {
			/*typename CLambdaLM_Solver::_TyEdgeMultiPool::_TyConstBaseRef e = m_r_system.r_Edge_Pool()[i]; // mouthful
			f_alpha = std::max(e.f_Max_VertexHessianDiagValue(), f_alpha);*/
			f_alpha = std::max(m_r_system.r_Edge_Pool()[i].f_Max_VertexHessianDiagValue(), f_alpha); // simple
		}

//#define __FIND_BOTTLENECK_EDGE
#ifdef __FIND_BOTTLENECK_EDGE
		double f_average = 0;
		for(size_t i = 0, n = m_r_system.n_Edge_Num(); i < n; ++ i)
			f_average += m_r_system.r_Edge_Pool()[i].f_Max_VertexHessianDiagValue();
		f_average /= m_r_system.n_Edge_Num();
		double f_thresh = f_alpha / 1000000;
		for(size_t i = 0, n = m_r_system.n_Edge_Num(); i < n; ++ i) {
			double f_max_diag = m_r_system.r_Edge_Pool()[i].f_Max_VertexHessianDiagValue();
			if(f_max_diag > f_thresh) {
				double f_ratio = f_max_diag / f_average;
				fprintf(stderr, "warning: edge " PRIsize " exceeds the average hessian "
					"diagonal value by %d orders of magnitude\n",
					i, int(log(f_ratio) / log(10.0) + .5));
				fprintf(stderr, "warning: edge " PRIsize ": (", i);
				for(size_t j = 0, m = m_r_system.r_Edge_Pool()[i].n_Vertex_Num(); j < m; ++ j)
					printf((j)? ", " PRIsize : PRIsize, m_r_system.r_Edge_Pool()[i].n_Vertex_Id(j));
				printf(")\n");

				//const_cast<typename CLambdaLM_Solver::_TySystem&>(m_r_system).r_Edge_Pool()[i].Calculate_Hessians_v2(); // just for debuggingh purposes, to see how the big number comes into existence
				const_cast<typename CLambdaLM_Solver::_TySystem&>(m_r_system).r_Edge_Pool().For_Each(i, i + 1,
					typename CLambdaLM_Solver::_TyLambdaOps::CCalculate_Hessians_v2());
			}
		}
#endif // __FIND_BOTTLENECK_EDGE
		// debug - see where the extreme derivatives are produced

		f_alpha *= tau;
		//f_alpha = 1.44; // copy from g2o for 10khogman

		return f_alpha;
	}

	/**
	 *	@copydoc CTrustRegion_Interface::Aftermath
	 */
	bool Aftermath(double &r_f_last_error, double f_error, double &r_f_alpha, const CUberBlockMatrix &r_lambda,
		const CLambdaLM_Solver &r_solver, const Eigen::VectorXd &r_v_dx, const Eigen::VectorXd &r_v_rhs)
	{
		double rho = (r_f_last_error - f_error) / (r_v_dx.transpose()).dot(r_f_alpha * r_v_dx + r_v_rhs);
		if(rho > 0) {
			r_f_alpha *= std::max(1 / 3.0, 1.0 - pow((2 * rho - 1), 3));
			m_nu = 2;

			r_f_last_error = f_error;
			// step taken: update the last error to this error (not all strategies
			// may want to do this - hence it is a part of TR policy)

			return true; // good step
		} else {
			r_f_alpha *= m_nu; // g2o has ni
			m_nu *= 2;

			return false; // fail
		}
	}

	/**
	 *	@copydoc CTrustRegion_Interface::ApplyDamping
	 */
	void ApplyDamping(CUberBlockMatrix &r_lambda, double f_alpha, size_t n_first_vertex, size_t n_last_vertex)
	{
		_ASSERTE(n_first_vertex <= n_last_vertex);
		_ASSERTE(n_last_vertex <= r_lambda.n_BlockColumn_Num());
		for(size_t i = n_first_vertex; i < n_last_vertex; ++ i) {
			size_t n_last = r_lambda.n_BlockColumn_Block_Num(i) - 1;
			CUberBlockMatrix::_TyMatrixXdRef block = r_lambda.t_Block_AtColumn(i, n_last);
			for(size_t j = 0, m = block.rows(); j < m; ++ j) // could use unroll (but wouldn't be able to use SSE anyway)
				block(j, j) += f_alpha;
		}
	}
};

/**
 *	@brief nonlinear blocky solver working above the lambda matrix
 *
 *	@tparam CSystem is optimization system type
 *	@tparam CLinearSolver is linear solver type
 *	@tparam CAMatrixBlockSizes is list of block sizes in the Jacobian matrix
 *	@tparam CLambdaMatrixBlockSizes is list of block sizes in the information (Hessian) matrix
 *	@tparam CTRAlgorithm is trust region algorithm
 */
template <class CSystem, class CLinearSolver, class CAMatrixBlockSizes = typename CSystem::_TyJacobianMatrixBlockList,
	class CLambdaMatrixBlockSizes = typename CSystem::_TyHessianMatrixBlockList,
	template <class> class CTRAlgorithm = CLevenbergMarquardt_Baseline>
class CNonlinearSolver_Lambda_LM {
public:
	typedef CTRAlgorithm<CNonlinearSolver_Lambda_LM<CSystem, CLinearSolver,
		CAMatrixBlockSizes, CLambdaMatrixBlockSizes, CTRAlgorithm> > _TyTRAlgorithm; /**< @brief trust-region algorithm type */

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

	typedef /*typename CUniqueTypelist<*/CAMatrixBlockSizes/*>::_TyResult*/ _TyAMatrixBlockSizes; /**< @brief possible block matrices, that can be found in A */
	typedef /*typename*/ CLambdaMatrixBlockSizes/*fbs_ut::CBlockMatrixTypesAfterPreMultiplyWithSelfTranspose<
		_TyAMatrixBlockSizes>::_TySizeList*/ _TyLambdaMatrixBlockSizes; /**< @brief possible block matrices, found in lambda and R */

	typedef typename CChooseType<lambda_utils::CLambdaOps<_TyLambdaMatrixBlockSizes>,
		lambda_utils::CLambdaOps2<_TyLambdaMatrixBlockSizes>, !base_iface::lambda_ReductionPlan_v2>::_TyResult _TyLambdaOps; /**< @brief implementation of operations for filling the lambda matrix */
	typedef typename _TyLambdaOps::_TyReductionPlan _TyReductionPlan; /**< @brief reduction plan implementation */

	/**
	 *	@brief solver interface properties, stored as enum (see also CSolverTraits)
	 */
	enum {
		solver_HasDump = true, /**< @brief timing statistics support flag */
		solver_HasChi2 = true, /**< @brief Chi2 error calculation support flag */
		solver_HasMarginals = true, /**< @brief marginal covariance support flag */
		solver_HasGaussNewton = false, /**< @brief Gauss-Newton support flag */
		solver_HasLevenberg = true, /**< @brief Levenberg-Marquardt support flag */
		solver_HasGradient = false, /**< @brief gradient-based linear solving support flag */
		solver_HasSchur = true, /**< @brief Schur complement support flag */
		solver_HasDelayedOptimization = false, /**< @brief delayed optimization support flag */
		solver_IsPreferredBatch = true, /**< @brief preferred batch solver flag */
		solver_IsPreferredIncremental = false, /**< @brief preferred incremental solver flag */
		solver_ExportsJacobian = false, /**< @brief interface for exporting jacobian system matrix flag */
		solver_ExportsHessian = false, /**< @brief interface for exporting hessian system matrix flag */
		solver_ExportsFactor = false /**< @brief interface for exporting factorized system matrix flag */
	};

protected:
	CSystem &m_r_system; /**< @brief reference to the system */
	CLinearSolver m_linear_solver; /**< @brief linear solver */
	CLinearSolver_Schur<CLinearSolver, _TyAMatrixBlockSizes, CSystem> m_schur_solver; /**< @brief linear solver with Schur trick */

	_TyTRAlgorithm m_TR_algorithm; /**< @brief Levenberg-Marquardt algorithm */

	CUberBlockMatrix m_lambda; /**< @brief the lambda matrix (built / updated incrementally) */
	_TyReductionPlan m_reduction_plan; /**< @brief lambda incremental reduction plan */
	Eigen::VectorXd m_v_dx; /**< @brief dx vector */
	Eigen::VectorXd m_v_rhs; /**< @brief right hand side vector */
	Eigen::VectorXd m_v_saved_state; /**< @brief saved state of the vertices (for LM step back) */
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
	bool m_b_use_schur; /**< @brief use schur complement flag */

	size_t m_n_real_step; /**< @brief counter of incremental steps (no modulo) */

	bool m_b_system_dirty; /**< @brief system updated without relinearization flag */

	size_t m_n_iteration_num; /**< @brief number of linear solver iterations */
	double m_f_chi2_time; /**< @brief time spent in Choleski() section */
	double m_f_lambda_time; /**< @brief time spent updating lambda */
	double m_f_rhs_time; /**< @brief time spent in right-hand-side calculation */
	double m_f_linsolve_time; /**< @brief time spent in luinear solving (Cholesky / Schur) */
	double m_f_norm_time; /**< @brief time spent in norm calculation section */
	double m_f_damping_time; /**< @brief time spent in calculating initial damping */
	double m_f_dampingupd_time; /**< @brief time spent in updating the damping */
	double m_f_sysupdate_time; /**< @brief time spent in norm calculation section */

	TMarginalsComputationPolicy m_t_marginals_config; /**< @brief marginal covariance computation configuration */
	CMarginalCovariance m_marginals; /**< @brief marginals cache */
	double m_f_extra_chol_time; /**< @brief time spent in calculating extra Cholesky factorization for marginal covariances */
	double m_f_marginals_time; /**< @brief time spent in calculating marginal covariances (batch) */
	double m_f_incmarginals_time; /**< @brief time spent in calculating marginal covariances (update) */
	size_t m_n_incmarginals_num; /**< @brief number of times the marginals update ran instead of batch recalculation */

	CTimer m_timer; /**< @brief timer object */

	bool m_b_had_loop_closure; /**< @brief (probable) loop closure flag */

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
	 *	@param[in] b_use_schur is Schur complement trick flag
	 *
	 *	@deprecated This is deprecated version of the constructor, use constructor
	 *		with TIncrementalSolveSetting instead.
	 */
	CNonlinearSolver_Lambda_LM(CSystem &r_system, size_t n_linear_solve_threshold,
		size_t n_nonlinear_solve_threshold, size_t n_nonlinear_solve_max_iteration_num = 5,
		double f_nonlinear_solve_error_threshold = .01, bool b_verbose = false,
		CLinearSolver linear_solver = CLinearSolver(), bool b_use_schur = true)
		:m_r_system(r_system), m_linear_solver(linear_solver),
		m_schur_solver(linear_solver), m_n_verts_in_lambda(0), m_n_edges_in_lambda(0),
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_last_optimized_vertex_num(0),
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_step(0),
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_linear_solve_threshold(n_linear_solve_threshold),
		m_n_nonlinear_solve_threshold(n_nonlinear_solve_threshold),
		m_n_nonlinear_solve_max_iteration_num(n_nonlinear_solve_max_iteration_num),
		m_f_nonlinear_solve_error_threshold(f_nonlinear_solve_error_threshold),
		m_b_verbose(b_verbose), m_b_use_schur(b_use_schur), m_n_real_step(0), m_b_system_dirty(false),
		m_n_iteration_num(0), m_f_chi2_time(0), m_f_lambda_time(0), m_f_rhs_time(0), m_f_linsolve_time(0),
		m_f_norm_time(0), m_f_damping_time(0), m_f_dampingupd_time(0), m_f_sysupdate_time(0),
		m_b_had_loop_closure(false), m_f_extra_chol_time(0), m_f_marginals_time(0),
		m_f_incmarginals_time(0), m_n_incmarginals_num(0)
	{
		_ASSERTE(!m_n_nonlinear_solve_threshold || !m_n_linear_solve_threshold); // only one of those
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
	 *	@param[in] b_use_schur is Schur complement trick flag
	 */
	CNonlinearSolver_Lambda_LM(CSystem &r_system,
		TIncrementalSolveSetting t_incremental_config = TIncrementalSolveSetting(),
		TMarginalsComputationPolicy t_marginals_config = TMarginalsComputationPolicy(),
		bool b_verbose = false,
		CLinearSolver linear_solver = CLinearSolver(), bool b_use_schur = true)
		:m_r_system(r_system), m_linear_solver(linear_solver),
		m_schur_solver(linear_solver), m_n_verts_in_lambda(0), m_n_edges_in_lambda(0),
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_last_optimized_vertex_num(0),
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_step(0),
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_linear_solve_threshold(t_incremental_config.t_linear_freq.n_period),
		m_n_nonlinear_solve_threshold(t_incremental_config.t_nonlinear_freq.n_period),
		m_n_nonlinear_solve_max_iteration_num(t_incremental_config.n_max_nonlinear_iteration_num),
		m_f_nonlinear_solve_error_threshold(t_incremental_config.f_nonlinear_error_thresh),
		m_b_verbose(b_verbose), m_b_use_schur(b_use_schur), m_n_real_step(0), m_b_system_dirty(false),
		m_n_iteration_num(0), m_f_chi2_time(0), m_f_lambda_time(0), m_f_rhs_time(0),
		m_f_linsolve_time(0), m_f_norm_time(0), m_f_damping_time(0), m_f_dampingupd_time(0),
		m_f_sysupdate_time(0), m_b_had_loop_closure(false), m_t_marginals_config(t_marginals_config),
		m_f_extra_chol_time(0), m_f_marginals_time(0), m_f_incmarginals_time(0), m_n_incmarginals_num(0)
	{
		_ASSERTE(!m_n_nonlinear_solve_threshold || !m_n_linear_solve_threshold); // only one of those

		/*if(t_marginals_config.b_calculate) // not supported at the moment
			throw std::runtime_error("CNonlinearSolver_Lambda_LM does not support marginals calculation");*/

		if(t_marginals_config.b_calculate) {
			if(t_marginals_config.t_increment_freq.n_period != t_incremental_config.t_nonlinear_freq.n_period &&
			   t_marginals_config.t_increment_freq.n_period != t_incremental_config.t_linear_freq.n_period) {
				throw std::runtime_error("in CNonlinearSolver_Lambda, the marginals must"
					" be updated with the same frequency as the system");
			}
			// unfortunately, yes

			/*if(t_marginals_config.n_incremental_policy != (mpart_LastColumn | mpart_Diagonal)) {
				throw std::runtime_error("in CNonlinearSolver_Lambda, the marginals update"
					" policy must be mpart_LastColumn | mpart_Diagonal");
			}
			if(t_marginals_config.n_incremental_policy != t_marginals_config.n_relinearize_policy) {
				throw std::runtime_error("in CNonlinearSolver_Lambda, the marginals "
					" incremental and relinearize update policy must be the same");
			}*/ // these are now implemented
			if(t_marginals_config.n_cache_miss_policy != mpart_Nothing) {
				throw std::runtime_error("in CNonlinearSolver_Lambda, the marginals cache"
					" miss policy is not supported at the moment, sorry for inconvenience");
			}
			// nothing else is implemented so far
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
	 *	@brief displays performance info on stdout
	 *	@param[in] f_total_time is total time taken by everything (can be -1 to ommit)
	 */
	void Dump(double f_total_time = -1) const
	{
		printf("solver took " PRIsize " iterations\n", m_n_iteration_num); // debug, to be able to say we didn't botch it numerically
		double f_parallel_time = m_f_lambda_time + m_f_rhs_time;
		double f_all_time = f_parallel_time + m_f_chi2_time + m_f_linsolve_time -
			m_f_norm_time + m_f_damping_time + m_f_dampingupd_time + m_f_sysupdate_time +
			m_f_extra_chol_time + m_f_marginals_time + m_f_incmarginals_time;
		printf("solver spent %f seconds in parallelizable section (updating lambda; disparity %g seconds)\n",
			f_parallel_time, (f_total_time > 0)? f_total_time - f_all_time : 0);
		printf("out of which:\n");
		printf("\tlambda: %f\n", m_f_lambda_time);
		printf("\t   rhs: %f\n", m_f_rhs_time);
		if(m_t_marginals_config.b_calculate) {
			printf("solver spent %f seconds in marginals\n"
				"\t chol: %f\n"
				"\tmargs: %f\n"
				"\t incm: %f (ran " PRIsize " times)\n",
				m_f_extra_chol_time + m_f_marginals_time + m_f_incmarginals_time,
				m_f_extra_chol_time, m_f_marginals_time,
				m_f_incmarginals_time, m_n_incmarginals_num);
		}
		printf("solver spent %f seconds in serial section\n", f_all_time - f_parallel_time);
		printf("out of which:\n");
		printf("\t  chi2: %f\n", m_f_chi2_time);
		printf("\t  damp: %f\n", m_f_damping_time);
		printf("\tlinsol: %f\n", m_f_linsolve_time);
		printf("\t  norm: %f\n", m_f_norm_time);
		printf("\tsysupd: %f\n", m_f_sysupdate_time);
		printf("\tdamupd: %f\n", m_f_dampingupd_time);
	}

	/**
	 *	@brief writes system matrix for art purposes
	 *
	 *	@param[in] p_s_filename is output file name (.tga)
	 *	@param[in] n_scalar_size is size of one scalar, in pixels
	 *
	 *	@return Returns true on success, false on failure.
	 */
	bool Dump_SystemMatrix(const char *p_s_filename, int n_scalar_size = 5)
	{
		try {
			_TyLambdaOps::Extend_Lambda(m_r_system, m_reduction_plan, m_lambda,
				m_n_verts_in_lambda, m_n_edges_in_lambda);
			if(!m_b_system_dirty)
				_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda, 0, m_n_edges_in_lambda);
			else
				_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda);
			m_b_system_dirty = false;
			m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
			m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size();
			_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
				m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
				m_lambda.n_Column_Num() == m_r_system.n_VertexElement_Num()); // lambda is square, blocks on either side = number of vertices
			// need to have lambda
		} catch(std::bad_alloc&) {
			return false;
		}

		return m_lambda.Rasterize(p_s_filename, n_scalar_size);
	}

	/**
	 *	@brief writes system matrix in matrix market for benchmarking purposes
	 *
	 *	@param[in] p_s_filename is output file name (.mtx)
	 *
	 *	@return Returns true on success, false on failure.
	 */
	bool Save_SystemMatrix_MM(const char *p_s_filename) const
	{
		char p_s_layout_file[256];
		strcpy(p_s_layout_file, p_s_filename);
		if(strrchr(p_s_layout_file, '.'))
			*(char*)strrchr(p_s_layout_file, '.') = 0;
		strcat(p_s_layout_file, ".bla");
		// only really required for landmark datasets

		return m_lambda.Save_MatrixMarket(p_s_filename, p_s_layout_file, "lambda matrix for SLAM problem");
	}

	/**
	 *	@brief calculates chi-squared error
	 *	@return Returns chi-squared error.
	 *	@note This only works with systems with edges of one degree of freedom
	 *		(won't work for e.g. systems with both poses and landmarks).
	 */
	double f_Chi_Squared_Error() const
	{
		return _TyLambdaOps::f_Chi_Squared_Error(m_r_system);
	}

	/**
	 *	@brief calculates denormalized chi-squared error
	 *	@return Returns denormalized chi-squared error.
	 *	@note This doesn't perform the final division by (number of edges - degree of freedoms).
	 */
	double f_Chi_Squared_Error_Denorm() const
	{
		return _TyLambdaOps::f_Chi_Squared_Error_Denorm(m_r_system);
	}

	/**
	 *	@brief incremental optimization function
	 *	@param[in] r_last_edge is the last edge that was added to the system
	 *	@note This function throws std::bad_alloc.
	 */
	void Incremental_Step(_TyBaseEdge &UNUSED(r_last_edge)) // throw(std::bad_alloc)
	{
		size_t n_vertex_num = m_r_system.r_Vertex_Pool().n_Size();
		if(!m_b_had_loop_closure) {
			_ASSERTE(!m_r_system.r_Edge_Pool().b_Empty());
			typename _TyEdgeMultiPool::_TyConstBaseRef r_edge =
				m_r_system.r_Edge_Pool()[m_r_system.r_Edge_Pool().n_Size() - 1];
			// get a reference to the last edge interface

			if(r_edge.n_Vertex_Num() > 1) { // unary factors do not cause classical loop closures
				_ASSERTE(r_edge.n_Vertex_Id(0) != r_edge.n_Vertex_Id(1));
				size_t n_first_vertex = std::min(r_edge.n_Vertex_Id(0), r_edge.n_Vertex_Id(1));
				m_b_had_loop_closure = (n_first_vertex < n_vertex_num - 2);
				_ASSERTE(m_b_had_loop_closure || std::max(r_edge.n_Vertex_Id(0),
					r_edge.n_Vertex_Id(1)) == n_vertex_num - 1);
			} else {
				size_t n_first_vertex = r_edge.n_Vertex_Id(0);
				m_b_had_loop_closure = (n_first_vertex < n_vertex_num - 1);
			}
			// todo - encapsulate this code, write code to support hyperedges as well, use it
		}
		// detect loop closures (otherwise the edges are initialized based on measurement and error would be zero)

		bool b_new_vert = false, b_ran_opt = false;

		++ m_n_real_step;

#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
		size_t n_new_vertex_num = n_vertex_num - m_n_last_optimized_vertex_num;
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

			if(m_b_had_loop_closure) {
				b_ran_opt = true;
				m_b_had_loop_closure = false;
				Optimize(m_n_nonlinear_solve_max_iteration_num, m_f_nonlinear_solve_error_threshold);
			}
			// nonlinear optimization

			b_new_vert = true;
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

			if(m_b_had_loop_closure) {
				b_ran_opt = true;
				m_b_had_loop_closure = false;
				Optimize(1, 0); // only if there was a loop (ignores possibly high residual after single step optimization)
			}
			// simple optimization

			b_new_vert = true;
		}

		if(b_new_vert && !b_ran_opt && m_t_marginals_config.b_calculate)
			Optimize(0, 0);
		// run optimization in order to calculate marginals after each vertex

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
		if(b_new_vert) {
			FILE *p_fw = fopen("chi2perVert.txt", (m_n_real_step > 0)? "a" : "w");
			fprintf(p_fw, "%f\n", f_Chi_Squared_Error_Denorm());
			fclose(p_fw);
		}
		// dump chi2
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
	}

	/**
	 *	@brief final optimization function
	 *
	 *	@param[in] n_max_iteration_num is the maximal number of iterations
	 *	@param[in] f_min_dx_norm is the residual norm threshold
	 */
	void Optimize(size_t n_max_iteration_num = 5, double f_min_dx_norm = .01) // throw(std::bad_alloc)
	{
		CTimerSampler timer(m_timer);

		const size_t n_variables_size = m_r_system.n_VertexElement_Num();
		const size_t n_measurements_size = m_r_system.n_EdgeElement_Num();
		if(n_variables_size > n_measurements_size) {
			if(n_measurements_size)
				fprintf(stderr, "warning: the system is underspecified\n");
			else
				fprintf(stderr, "warning: the system contains no edges at all: nothing to optimize\n");
			//return;
		}
		if(!n_measurements_size)
			return; // nothing to solve (but no results need to be generated so it's ok)
		// can't solve in such conditions

		_TyLambdaOps::Extend_Lambda(m_r_system, m_reduction_plan, m_lambda,
			m_n_verts_in_lambda, m_n_edges_in_lambda); // recalculated all the jacobians inside Extend_Lambda()
		m_v_dx.resize(n_variables_size, 1);
		m_v_saved_state.resize(n_variables_size, 1);
		// allocate more memory

		if(!m_b_system_dirty)
			_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda, 0, m_n_edges_in_lambda); // calculate only for new edges // @todo - but how to mark affected vertices? // simple test if edge id is greater than m_n_edges_in_lambda, the vertex needs to be recalculated
		else
			_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda); // calculate for entire system
		//m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
		//m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // not yet, the damping was not applied
		//m_b_system_dirty = false; // cannot do that yet, the damping was not applied
		// lambda is (partially) updated, but (partially) without damping - don't give it to m_TR_algorithm as it would quite increase complexity

		timer.Accum_DiffSample(m_f_lambda_time);

		if(m_lambda.n_BlockColumn_Num() < m_r_system.r_Vertex_Pool().n_Size()) {
			fprintf(stderr, "warning: waiting for more edges\n");
			m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
			m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // probably necessary here
			m_b_system_dirty = true; // to correctly adjust damping to the whole matrix
			return;
		}
		// waiting for more edges

#ifdef __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
		printf("warning: landmark settle iterations enabled\n");
#endif // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS

		double f_alpha = m_TR_algorithm.f_InitialDamping(m_r_system); // lambda is not ready at this point yet, can only use the edges
		if(m_b_verbose) {
			printf("alpha: %f\n", f_alpha);

			/*std::vector<float> pix_err(m_r_system.r_Edge_Pool().n_Size());
			for(size_t i = 0, n = m_r_system.r_Edge_Pool().n_Size(); i < n; ++ i) {
				//const CEdgeP2C3D &r_edge = m_r_system.r_Edge_Pool().r_At<CEdgeP2C3D>(i);
				pix_err[i] = float(m_r_system.r_Edge_Pool().r_At<_TyBaseEdge>(i).f_Reprojection_Error());
			}
			size_t n_med = pix_err.size() / 2;
			std::nth_element(pix_err.begin(), pix_err.begin() + n_med, pix_err.end());
			printf("median reprojection error: %.2f px\n", pix_err[n_med]);*/
			// debug - print median reprojection error
		}

		double f_errorx = m_TR_algorithm.f_Error(*this);
		//printf("init chi: %f\n", f_errorx);

		timer.Accum_DiffSample(m_f_damping_time);

		if(!m_b_system_dirty)
			Apply_Damping(0, m_n_edges_in_lambda, f_alpha); // calculate only for new edges // @todo - but how to mark affected vertices?
		else
			Apply_Damping(0, 0, f_alpha); // for the entire system

		m_b_system_dirty = false;
		m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
		m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
		_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
			m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
			m_lambda.n_Column_Num() == n_variables_size); // lambda is square, blocks on either side = number of vertices
		// need to have lambda

		timer.Accum_DiffSample(m_f_lambda_time);

		double f_last_error = m_TR_algorithm.f_Error(*this);
		// this may not be called until lambda is finished

		timer.Accum_DiffSample(m_f_chi2_time);

		if(m_b_verbose) {
			size_t n_sys_size = m_r_system.n_Allocation_Size();
			size_t n_rp_size = m_reduction_plan.n_Allocation_Size();
			size_t n_lam_size = m_lambda.n_Allocation_Size();
			printf("memory_use(sys: %.2f MB, redplan: %.2f MB, ,\\: %.2f MB)\n",
				n_sys_size / 1048576.0, n_rp_size / 1048576.0, n_lam_size / 1048576.0);
		}
		// print memory use statistics

		int fail = 10;
#ifdef __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
		for(size_t n_iteration = 0; n_iteration < n_max_iteration_num * 4; ++ n_iteration) {
#else // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
		for(size_t n_iteration = 0; n_iteration < n_max_iteration_num; ++ n_iteration) {
#endif // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
			++ m_n_iteration_num;
			// debug

			if(m_b_verbose) {
				if(n_max_iteration_num == 1)
					printf("\n=== incremental optimization step ===\n\n");
				else {
#ifdef __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
					if(!(n_iteration % 4))
						printf("\n=== nonlinear optimization: iter #" PRIsize " ===\n\n", n_iteration / 4);
					else {
						printf("\n=== nonlinear optimization: iter #" PRIsize
							", settle iter #" PRIsize " ===\n\n", n_iteration / 4, n_iteration % 4);
					}
#else // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
					printf("\n=== nonlinear optimization: iter #" PRIsize " ===\n\n", n_iteration);
#endif // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
				}
			}
			// verbose

			if(n_iteration && m_b_system_dirty) {
				_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda, 0, 0);
				Apply_Damping(0, 0, f_alpha);
				m_b_system_dirty = false;

				timer.Accum_DiffSample(m_f_lambda_time);
			}
			// no need to rebuild lambda, just refresh the values that are being referenced

			_TyLambdaOps::Collect_RightHandSide_Vector(m_r_system, m_reduction_plan, m_v_dx);
			// collects the right-hand side vector

			m_v_rhs = m_v_dx; // copy intended
			// save the original rhs for step estimation

			timer.Accum_DiffSample(m_f_rhs_time);

			bool b_cholesky_result = LinearSolve(n_iteration, n_max_iteration_num);
			// calculate cholesky, reuse block ordering if the linear solver supports it

			timer.Accum_DiffSample(m_f_linsolve_time);

			if(!b_cholesky_result)
				break;
			// in case cholesky failed, quit

			double f_residual_norm = 0;
			if(b_cholesky_result) {
				f_residual_norm = m_v_dx.norm(); // Eigen likely uses SSE and OpenMP
				if(m_b_verbose)
					printf("residual norm: %.4f\n", f_residual_norm);
			}

			// calculate residual norm

			timer.Accum_DiffSample(m_f_norm_time);
			// timing breakup

			if(f_residual_norm <= f_min_dx_norm/* && n_iteration % 4 == 0*/)
				break;
			// in case the error is low enough, quit (saves us recalculating the hessians)

#ifdef __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
			if(n_iteration % 4 == 0) // save state at the beginning!
				base_iface::CSolverOps_Base::Save_State(m_r_system, m_v_saved_state);
#else // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
			base_iface::CSolverOps_Base::Save_State(m_r_system, m_v_saved_state);
#endif // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
			// god save the vertices

			base_iface::CSolverOps_Base::PushValuesInGraphSystem(m_r_system, m_v_dx);
			// update the system (in parallel)

			timer.Accum_DiffSample(m_f_sysupdate_time);

#ifdef __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
			if(n_iteration % 4 < 3) { // evaluate chi2 at the end of the relaxation (this is maybe wrong, maybe a bad step is caught unnecessarily late)
				m_b_system_dirty = true;
				continue;
			}
#endif // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS

			double f_error = m_TR_algorithm.f_Error(*this);
			if(m_b_verbose) {
				printf("chi2: %f\n", f_error);

				/*std::vector<float> pix_err(m_r_system.r_Edge_Pool().n_Size());
				for(size_t i = 0, n = m_r_system.r_Edge_Pool().n_Size(); i < n; ++ i) {
					//const CEdgeP2C3D &r_edge = m_r_system.r_Edge_Pool().r_At<CEdgeP2C3D>(i);
					pix_err[i] = float(m_r_system.r_Edge_Pool().r_At<_TyBaseEdge>(i).f_Reprojection_Error());
				}
				size_t n_med = pix_err.size() / 2;
				std::nth_element(pix_err.begin(), pix_err.begin() + n_med, pix_err.end());
				printf("median reprojection error: %.2f px\n", pix_err[n_med]);
				// debug - print median reprojection error*/
			}

			timer.Accum_DiffSample(m_f_chi2_time);

			bool b_good_step = m_TR_algorithm.Aftermath(f_last_error, f_error, f_alpha, m_lambda, *this, m_v_dx, m_v_rhs);
			if(!b_good_step) {
				fprintf(stderr, "warning: chi2 rising\n");
				// verbose

				base_iface::CSolverOps_Base::Load_State(m_r_system, m_v_saved_state);
				// restore saved vertives

				if(fail > 0) {
					-- fail;
					n_max_iteration_num ++;
				}
				// increase the number of iterations, up to a certain limit
			} else
				m_marginals.DisableUpdate(); // linearization point just changed, all the marginals will change - need full recalc

			timer.Accum_DiffSample(m_f_dampingupd_time);
			m_b_system_dirty = true; // even though we saved state, we need to change alpha
		}
		// optimize the system

		if(m_t_marginals_config.b_calculate) {
			bool b_batch = !m_n_linear_solve_threshold && !m_n_nonlinear_solve_threshold;
			// are we running batch?

			if(b_batch) {
				m_linear_solver.Free_Memory();
				m_schur_solver.Free_Memory();
			}
			// unable to reuse these, free memory

			if(m_b_verbose && b_batch)
				printf("\n=== calculating marginals ===\n\n");
			// todo - handle freq settings
			// todo - handle policies

			if(m_b_verbose && b_batch)
				printf("refreshing lambda with null damping\n");

			m_b_system_dirty = true;
			if(f_alpha > 0 || m_b_system_dirty) {
				f_alpha = 0; // !! otherwise the marginals are something else
				_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda, 0, 0);
				Apply_Damping(0, 0, f_alpha);
				//m_b_system_dirty = false; // this will break something

				timer.Accum_DiffSample(m_f_lambda_time);
			}
			_ASSERTE(m_n_verts_in_lambda == m_r_system.r_Vertex_Pool().n_Size());
			_ASSERTE(m_n_edges_in_lambda == m_r_system.r_Edge_Pool().n_Size());
			_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
				m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
				m_lambda.n_Column_Num() == n_variables_size);
			// need to update or will end up with forever bad marginals!

			CUberBlockMatrix R;
			//R.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(m_lambda); // this makes ugly dense factor, dont do that

			//m_linear_solver.Factorize_PosDef_Blocky(R, m_lambda, std::vector<size_t>()); // dense as well, no ordering inside

			if(m_b_verbose && b_batch)
				printf("calculating fill-reducing ordering\n");

			CMatrixOrdering mord;
			if((m_marginals.b_CanUpdate() && (m_t_marginals_config.n_incremental_policy &
			   mpart_LastColumn) == mpart_LastColumn) || // can tell for sure if incremental is going to be used
			   (m_t_marginals_config.n_relinearize_policy & mpart_LastColumn) == mpart_LastColumn) { // never know if we fallback to batch, though
				CLastElementOrderingConstraint leoc;
				mord.p_BlockOrdering(m_lambda, leoc.p_Get(m_lambda.n_BlockColumn_Num()),
					m_lambda.n_BlockColumn_Num(), true); // constrain the last column to be the last column (a quick fix) // todo - handle this properly, will be unable to constrain like this in fast R (well, actually ...) // todo - see what is this doing to the speed
			} else
				mord.p_BlockOrdering(m_lambda, true); // unconstrained; the last column may be anywhere (could reuse R from the linear solver here - relevant also in batch (e.g. on venice))
			const size_t *p_order = mord.p_Get_InverseOrdering();
			{
				CUberBlockMatrix lambda_perm; // not needed afterwards

				if(m_b_verbose && b_batch)
					printf("forming lambda perm\n");

				m_lambda.Permute_UpperTriangular_To(lambda_perm, p_order, mord.n_Ordering_Size(), true);

				if(m_b_verbose && b_batch)
					printf("calculating Cholesky\n");

				if(!R.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(lambda_perm))
					throw std::runtime_error("fatal error: R.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(lambda_perm) failed");
				// note that now the marginals are calculated with ordering: need to count with that, otherwise those are useless!

				if(m_b_verbose && b_batch)
					printf("memory_use(R: %.2f MB)\n", R.n_Allocation_Size() / 1048576.0);
				// verbose
			}
			if(m_b_verbose && b_batch)
				printf("calculating marginals (%s)\n", (m_marginals.b_CanUpdate())? "incrementally" : "batch");

			// todo - reuse what the linear solver calculated, if we have it (not if schur, )
			// todo - think of what happens if using schur ... have to accelerate the dense margs differently
			//		probably the whole mindset of having R is wrong, it would be much better to leave
			//		it up to the linear solver to solve for the columns

			/*printf("debug: matrix size: " PRIsize " / " PRIsize " (" PRIsize " nnz)\n",
				lambda_perm.n_BlockColumn_Num(), lambda_perm.n_Column_Num(), lambda_perm.n_NonZero_Num());
			float f_avg_block_size = float(lambda_perm.n_Column_Num()) / lambda_perm.n_BlockColumn_Num();
			printf("debug: diagonal nnz: " PRIsize "\n", size_t(lambda_perm.n_BlockColumn_Num() *
				(f_avg_block_size * f_avg_block_size)));
			printf("debug: factor size: " PRIsize " / " PRIsize " (" PRIsize " nnz)\n",
				R.n_BlockColumn_Num(), R.n_Column_Num(), R.n_NonZero_Num());*/
			// see how much we compute, compared to g2o

			timer.Accum_DiffSample(m_f_extra_chol_time);

			size_t n_add_edge_num = m_r_system.r_Edge_Pool().n_Size() - m_marginals.n_Edge_Num();
			bool b_incremental = m_marginals.b_CanUpdate() && CMarginals::b_PreferIncremental(m_r_system,
				m_marginals.r_SparseMatrix(), m_lambda, R, mord, m_marginals.n_Edge_Num(),
				m_t_marginals_config.n_incremental_policy);
//#ifdef __SE_TYPES_SUPPORT_L_SOLVERS
			if(b_incremental) { // otherwise just update what we have
				CUberBlockMatrix &r_m = const_cast<CUberBlockMatrix&>(m_marginals.r_SparseMatrix()); // watch out, need to call Swap_SparseMatrix() afterwards
				if(!CMarginals::Update_BlockDiagonalMarginals_FBS<false>(m_r_system, r_m, m_lambda,
				   R, mord, m_marginals.n_Edge_Num(), m_t_marginals_config.n_incremental_policy)) {
#ifdef _DEBUG
					fprintf(stderr, "warning: Update_BlockDiagonalMarginals_FBS() had a numerical issue:"
						" restarting with Calculate_DenseMarginals_Recurrent_FBS() instead\n");
#endif // _DEBUG
					b_incremental = false;
					// failed, will enter the batch branch below, that will not have a numerical issue
				} else {
					m_marginals.Swap_SparseMatrix(r_m); // now the marginals know that the matrix changed

					timer.Accum_DiffSample(m_f_incmarginals_time);
					++ m_n_incmarginals_num;
				}
			}
/*#else // __SE_TYPES_SUPPORT_L_SOLVERS
#pragma message("warning: the fast incremental marginals not available: __SE_TYPES_SUPPORT_L_SOLVERS not defined")
			b_incremental = false;
#endif // __SE_TYPES_SUPPORT_L_SOLVERS*/
			if(!b_incremental) { // if need batch marginals
				CUberBlockMatrix margs_ordered;
				CMarginals::Calculate_DenseMarginals_Recurrent_FBS<_TyLambdaMatrixBlockSizes>(margs_ordered, R,
					mord, m_t_marginals_config.n_relinearize_policy, false);
				// calculate the thing

				{
					CUberBlockMatrix empty;
					R.Swap(empty);
				}
				// delete R, don't need it and it eats a lot of memory

				if(m_b_verbose && b_batch)
					printf("reordering the marginals\n");

				CUberBlockMatrix &r_m = const_cast<CUberBlockMatrix&>(m_marginals.r_SparseMatrix()); // watch out, need to call Swap_SparseMatrix() afterwards
				margs_ordered.Permute_UpperTriangular_To(r_m, mord.p_Get_Ordering(),
					mord.n_Ordering_Size(), false); // no share! the original will be deleted
				m_marginals.Swap_SparseMatrix(r_m); // now the marginals know that the matrix changed
				// take care of having the correct permutation there

				m_marginals.EnableUpdate();
				// now the marginals are current and can be updated until the linearization point is changed again

				timer.Accum_DiffSample(m_f_marginals_time);
			}

			m_marginals.Set_Edge_Num(m_r_system.r_Edge_Pool().n_Size());
			// now all those edges are in the marginals

			FILE *p_fw;
			if((p_fw = fopen("marginals.txt", "w"))) {
				for(size_t i = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
					size_t n_order = m_lambda.n_BlockColumn_Base(i);
					size_t n_dimension = m_lambda.n_BlockColumn_Column_Num(i);
					// get col

					CUberBlockMatrix::_TyConstMatrixXdRef block =
						m_marginals.r_SparseMatrix().t_FindBlock(n_order, n_order);
					// get block

					//fprintf(p_fw, "block_%d_%d = ", int(i), int(i));
					//CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, block);
					// prints the matrix

					_ASSERTE(block.rows() == block.cols() && block.cols() == n_dimension);
					for(size_t i = 0; i < n_dimension; ++ i)
						fprintf(p_fw, (i)? " %lf" : "%lf", block(i, i));
					fprintf(p_fw, "\n");
					// print just the diagonal, one line per every vertex
				}
				fclose(p_fw);
			}
			// dump diagonal blocks of the marginals to a file
		}
		// now R is up to date, can get marginals
	}

protected:
	/**
	 *	@brief solves m_lambda \ m_v_dx, result left in m_v_dx
	 *
	 *	@param[in] n_iteration is nonlinear solver iteration (zero-based index)
	 *	@param[in] n_max_iteration_num is maximum nonlinear iteration count
	 *
	 *	@return Returns true if the factorization succeeded, otherwise returns false.
	 */
	bool LinearSolve(size_t n_iteration, size_t n_max_iteration_num)
	{
		bool b_cholesky_result;
		{
			Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta

			if(!m_b_use_schur) { // use cholesky
				if(n_max_iteration_num > 1) {
					do {
						if(!n_iteration &&
						   !_TyLinearSolverWrapper::FinalBlockStructure(m_linear_solver, m_lambda)) {
							b_cholesky_result = false;

							break;
						}
						// prepare symbolic factorization, structure of lambda won't change in the next steps
						b_cholesky_result = _TyLinearSolverWrapper::Solve(m_linear_solver, m_lambda, v_eta);
						// p_dx = eta = lambda / eta
					} while(0);
				} else
					b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda, v_eta); // p_dx = eta = lambda / eta
			} else { // use Schur complement

				bool b_marginalize_out_poses = false;
				// set to true to optimize only landmarks

#ifdef __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS
				b_marginalize_out_poses = true;
				if(n_iteration % 4 == 0)
					b_marginalize_out_poses = false;
#endif // __NONLINEAR_SOLVER_LAMBDA_LM_LANDMARK_SETTLE_ITERATIONS

				if(!n_iteration) {
					bool b_force_guided_ordering = b_marginalize_out_poses; // or set this to true to make sure that it will be the poses and not landmarks or a mixture thereof that will be marginalized out in Schur_Solve_MarginalPoses()
					m_schur_solver.SymbolicDecomposition_Blocky(m_lambda, b_force_guided_ordering);
				}
				// calculate the ordering once, it does not change

				if(!b_marginalize_out_poses)
					b_cholesky_result = m_schur_solver.Solve_PosDef_Blocky(m_lambda, v_eta); // as usual
				else
					b_cholesky_result = m_schur_solver.Solve_PosDef_Blocky_MarginalPoses(m_lambda, v_eta); // then m_v_dx contains change only for the landmarks and the poses remain constant
				// Schur
			}

			if(m_b_verbose) {
				printf("%s %s", (m_b_use_schur)? "Schur" : "Cholesky",
					(b_cholesky_result)? "succeeded\n" : "failed\n");
			}

#ifdef _DEBUG
			for(size_t i = 0, n_variables_size = m_v_dx.cols(); i < n_variables_size; ++ i) {
				if(_isnan(m_v_dx(i))) {
					fprintf(stderr, "warning: p_dx[" PRIsize "] = NaN (file \'%s\', line "
						PRIsize ")\n", i, __FILE__, size_t(__LINE__));
				}
			}
			// detect the NaNs, if any (warn, but don't modify)
#endif // _DEBUG
		}
		// calculate cholesky, reuse block ordering if the linear solver supports it

		return b_cholesky_result;
	}

	/**
	 *	@brief calls the TR algorithm to apply damping to the lambda matrix
	 *
	 *	@param[in] n_refresh_from_vertex is zero-based index of the first vertex to refresh (unused)
	 *	@param[in] n_refresh_from_edge is zero-based index of the first edge to refresh
	 *	@param[in] f_alpha is value of the damping factor
	 *
	 *	@note The damping is applied incrementally.
	 */
	void Apply_Damping(size_t n_refresh_from_vertex, size_t n_refresh_from_edge, double f_alpha)
	{
#ifdef __LAMBDA_USE_V2_REDUCTION_PLAN
		if(n_refresh_from_edge) {
			for(size_t e = n_refresh_from_edge, n = m_r_system.n_Edge_Num(); e < n; ++ e) { // todo - parallel
				typename CSystem::_TyConstEdgeRef r_edge = m_r_system.r_Edge_Pool()[e];
				size_t n_v0 = r_edge.n_Vertex_Id(0);
				size_t n_v1 = r_edge.n_Vertex_Id(1);
				/*if(n_v1 == n_v0 + 1) // can't really expect that to happen in BA
					m_TR_algorithm.ApplyDamping(m_lambda, f_alpha, n_v0, n_v1 + 1);
				else*/ {
					m_TR_algorithm.ApplyDamping(m_lambda, f_alpha, n_v0, n_v0 + 1); // only damp the new or updated vertices
					m_TR_algorithm.ApplyDamping(m_lambda, f_alpha, n_v1, n_v1 + 1); // only damp the new or updated vertices
				}
			}
			// refresh only vertices belonging to the new edges
		} else {
			m_TR_algorithm.ApplyDamping(m_lambda, f_alpha, n_refresh_from_vertex, m_lambda.n_BlockColumn_Num());
			// refresh the full diagonal
		}
#else // __LAMBDA_USE_V2_REDUCTION_PLAN
		m_TR_algorithm.ApplyDamping(m_lambda, f_alpha, n_refresh_from_vertex, m_lambda.n_BlockColumn_Num());
		// always all of the vertices (all of them updated)
#endif // __LAMBDA_USE_V2_REDUCTION_PLAN
	}

	CNonlinearSolver_Lambda_LM(const CNonlinearSolver_Lambda_LM &UNUSED(r_solver)); /**< @brief the object is not copyable */
	CNonlinearSolver_Lambda_LM &operator =(const CNonlinearSolver_Lambda_LM &UNUSED(r_solver)) { return *this; } /**< @brief the object is not copyable */
};

#endif // !__NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
