/*
								+-----------------------------------+
								|                                   |
								| ***  Lambda nonlinear solver  *** |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2012  |
								|                                   |
								|     NonlinearSolver_Lambda.h      |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
#define __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED

/**
 *	@file include/slam/NonlinearSolver_Lambda.h
 *	@brief nonlinear blocky solver working above the lambda matrix
 *	@author -tHE SWINe-
 *	@date 2012-09-13
 */

#include "slam/FlatSystem.h"
#include <iostream> // SOSO

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
 *	@brief enables writes of chi2 errors at each step
 */
//#define __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2

/**
 *	@brief nonlinear blocky solver working above the lambda matrix
 *
 *	@tparam CSystem is the system type
 *	@tparam CLinearSolver is a linear solver
 */
template <class CSystem, class CLinearSolver, class CAMatrixBlockSizes = typename CSystem::_TyJacobianMatrixBlockList>
class CNonlinearSolver_Lambda_LM {
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

protected:
	CSystem &m_r_system; /**< @brief reference to the system */
	CLinearSolver m_linear_solver; /**< @brief linear solver */
	CLinearSolver_Schur<CLinearSolver, _TyAMatrixBlockSizes> m_schur_solver; /**< @brief linear solver with Schur trick */

	CUberBlockMatrix m_lambda; /**< @brief the lambda matrix (built / updated incrementally) */
	Eigen::VectorXd m_v_dx; /**< @brief dx vector */
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
	double m_f_serial_time; /**< @brief time spent in serial section */
	double m_f_ata_time; /**< @brief time spent in A^T * A section */
	double m_f_premul_time; /**< @brief time spent in b * A section */
	double m_f_chol_time; /**< @brief time spent in Choleski() section */
	double m_f_norm_time; /**< @brief time spent in norm calculation section */

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
	 */
	CNonlinearSolver_Lambda_LM(CSystem &r_system, size_t n_linear_solve_threshold,
		size_t n_nonlinear_solve_threshold, size_t n_nonlinear_solve_max_iteration_num,
		double f_nonlinear_solve_error_threshold, bool b_verbose,
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
		m_b_verbose(b_verbose), m_n_real_step(0), m_b_system_dirty(false),
		m_n_iteration_num(0), m_f_serial_time(0), m_f_ata_time(0), m_f_premul_time(0),
		m_f_chol_time(0), m_f_norm_time(0), m_b_had_loop_closure(false), m_b_use_schur(b_use_schur)
	{
		_ASSERTE(!m_n_nonlinear_solve_threshold || !m_n_linear_solve_threshold); // only one of those
	}

	/**
	 *	@brief displays performance info on stdout
	 *	@param[in] f_total_time is total time taken by everything (can be -1 to ommit)
	 */
	void Dump(double f_total_time = -1) const
	{
		printf("solver took " PRIsize " iterations\n", m_n_iteration_num); // debug, to be able to say we didn't botch it numerically
		if(f_total_time > 0)
			printf("solver spent %f seconds in parallelizable section (updating lambda)\n", f_total_time - m_f_serial_time);
		printf("solver spent %f seconds in serial section\n", m_f_serial_time);
		printf("out of which:\n");
		printf("\t  ata: %f\n", m_f_ata_time);
		printf("\tgaxpy: %f\n", m_f_premul_time);
		printf("\t chol: %f\n", m_f_chol_time);
		printf("\t norm: %f\n", m_f_norm_time);
		printf("\ttotal: %f\n", m_f_ata_time + m_f_premul_time + m_f_chol_time + m_f_norm_time);
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
			Extend_Lambda(m_n_verts_in_lambda, m_n_edges_in_lambda); // recalculated all the jacobians inside Extend_Lambda()
			if(!m_b_system_dirty)
				Refresh_Lambda(0/*m_n_verts_in_lambda*/, m_n_edges_in_lambda); // calculate only for new edges // @todo - but how to mark affected vertices? // simple test if edge id is greater than m_n_edges_in_lambda, the vertex needs to be recalculated
			else
				Refresh_Lambda(); // calculate for entire system
			m_b_system_dirty = false;
			m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
			m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
			//const size_t n_variables_size = m_r_system.n_VertexElement_Num();
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
	bool Save_SystemMatrix_MM(const char *p_s_filename)
	{
		size_t n_rows = m_lambda.n_Row_Num();
		size_t n_columns = m_lambda.n_Column_Num();
		size_t n_nonzeros = m_lambda.n_NonZero_Num();
		if(CTypelistLength<_TyAMatrixBlockSizes>::n_result > 1) { // only really required for landmark datasets
			char p_s_layout_file[256];
			strcpy(p_s_layout_file, p_s_filename);
			if(strrchr(p_s_layout_file, '.'))
				*(char*)strrchr(p_s_layout_file, '.') = 0;
			strcat(p_s_layout_file, "_block-layout.txt");
			FILE *p_fw;
			if(!(p_fw = fopen(p_s_layout_file, "w")))
				return false;
			fprintf(p_fw, PRIsize " x " PRIsize " (" PRIsize ")\n", n_rows, n_columns, n_nonzeros);
			fprintf(p_fw, PRIsize " x " PRIsize " (" PRIsize ")\n", m_lambda.n_BlockRow_Num(),
				m_lambda.n_BlockColumn_Num(), m_lambda.n_Block_Num());
			for(size_t i = 0, n = m_lambda.n_BlockRow_Num(); i < n; ++ i)
				fprintf(p_fw, PRIsize " ", m_lambda.n_BlockRow_Base(i));
			fprintf(p_fw, PRIsize "\n", n_rows);
			for(size_t i = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i)
				fprintf(p_fw, PRIsize " ", m_lambda.n_BlockColumn_Base(i));
			fprintf(p_fw, PRIsize "\n", n_columns);
			fclose(p_fw);
		}
		// save the block layout for the benchmarks

		FILE *p_fw;
		if(!(p_fw = fopen(p_s_filename, "w")))
			return false;
		fprintf(p_fw, "%%%%MatrixMarket matrix coordinate real symmetric\n");
		fprintf(p_fw, "%%-------------------------------------------------------------------------------\n");
		fprintf(p_fw, "%% SLAM_plus_plus matrix dump\n");
		fprintf(p_fw, "%% kind: lambda matrix for SLAM problem\n");
		fprintf(p_fw, "%%-------------------------------------------------------------------------------\n");
		fprintf(p_fw, "" PRIsize " " PRIsize " " PRIsize "\n", n_rows, n_columns, n_nonzeros);
		for(size_t n_col = 0, n_column_blocks = m_lambda.n_BlockColumn_Num();
		   n_col < n_column_blocks; ++ n_col) {
			size_t n_col_base = m_lambda.n_BlockColumn_Base(n_col);
			size_t n_col_width = m_lambda.n_BlockColumn_Column_Num(n_col);
			size_t n_block_num = m_lambda.n_BlockColumn_Block_Num(n_col);
			for(size_t j = 0; j < n_block_num; ++ j) {
				size_t n_row = m_lambda.n_Block_Row(n_col, j);
				size_t n_row_base = m_lambda.n_BlockRow_Base(n_row);
				size_t n_row_height = m_lambda.n_BlockRow_Row_Num(n_row);

				CUberBlockMatrix::_TyMatrixXdRef t_block = m_lambda.t_BlockAt(j, n_col);
				_ASSERTE(t_block.rows() == n_row_height && t_block.cols() == n_col_width);
				// get a block

				for(size_t k = 0; k < n_row_height; ++ k) {
					for(size_t l = 0; l < n_col_width; ++ l) {
						fprintf(p_fw, "" PRIsize " " PRIsize " %f\n", n_row_base + 1 + k,
							n_col_base + 1 + l, t_block(k, l));
					}
				}
				// print a block (includes the nulls, but so does n_NonZero_Num())
			}
			// for all blocks in a column
		}
		// for all columns

		if(ferror(p_fw)) {
			fclose(p_fw);
			return false;
		}
		fclose(p_fw);

		return true;
	}

	/**
	 *	@brief calculates chi-squared error
	 *	@return Returns chi-squared error.
	 *	@note This only works with systems with edges of one degree of freedom
	 *		(won't work for e.g. systems with both poses and landmarks).
	 */
	inline double f_Chi_Squared_Error() const
	{
		if(m_r_system.r_Edge_Pool().b_Empty())
			return 0;
		return m_r_system.r_Edge_Pool().For_Each(CSum_ChiSquareError()) /
			(m_r_system.r_Edge_Pool().n_Size() - m_r_system.r_Edge_Pool()[0].n_Dimension());
	}

	/**
	 *	@brief calculates denormalized chi-squared error
	 *	@return Returns denormalized chi-squared error.
	 *	@note This doesn't perform the final division by (number of edges - degree of freedoms).
	 */
	inline double f_Chi_Squared_Error_Denorm() const
	{
		if(m_r_system.r_Edge_Pool().b_Empty())
			return 0;
		return m_r_system.r_Edge_Pool().For_Each(CSum_ChiSquareError());
	}

	/**
	 *	@brief incremental optimization function
	 *	@param[in] r_last_edge is the last edge that was added to the system
	 */
	void Incremental_Step(_TyBaseEdge &r_last_edge) // throw(std::bad_alloc)
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

		/*FILE *p_fw = fopen("timeSteps_lambda.txt", (m_n_real_step > 0)? "a" : "w");
		fprintf(p_fw, "" PRIsize ";%f\n", m_n_real_step, m_timer.f_Time());
		fclose(p_fw);*/
		// dump time per step

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
		bool b_new_vert = false;
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2

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
				m_b_had_loop_closure = false;
				Optimize(m_n_nonlinear_solve_max_iteration_num, m_f_nonlinear_solve_error_threshold);
			}
			// nonlinear optimization

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
			b_new_vert = true;
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
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
				m_b_had_loop_closure = false;
				Optimize(1, 0); // only if there was a loop (ignores possibly high residual after single step optimization)
			}
			// simple optimization

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
			b_new_vert = true;
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
		}

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
	void Optimize(size_t n_max_iteration_num, double f_min_dx_norm) // throw(std::bad_alloc)
	{
		const size_t n_variables_size = m_r_system.n_VertexElement_Num();
		const size_t n_measurements_size = m_r_system.n_EdgeElement_Num();
		if(n_variables_size > n_measurements_size) {
			fprintf(stderr, "warning: the system is underspecified\n");
			return;
		}
		if(!n_measurements_size)
			return; // nothing to solve (but no results need to be generated so it's ok)
		// can't solve in such conditions
		Extend_Lambda(m_n_verts_in_lambda, m_n_edges_in_lambda); // recalculated all the jacobians inside Extend_Lambda()

		if(!m_b_system_dirty) {
			m_r_system.r_Edge_Pool().For_Each_Parallel(m_n_edges_in_lambda,
				m_r_system.r_Edge_Pool().n_Size(), CCalculate_Hessians());
		} else
			m_r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Hessians());

		double alpha = 0.0;
		const double tau = 1e-3;
		/*for(size_t i = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_last = m_lambda.n_BlockColumn_Block_Num(i) - 1;
			CUberBlockMatrix::_TyMatrixXdRef block = m_lambda.t_BlockAt(n_last, i);
			for(size_t j = 0, m = block.rows(); j < m; ++ j)
				alpha = std::max(block(j, j), alpha);
		}*/
		for(size_t i = 0, n = m_r_system.n_Edge_Num(); i < n; ++ i) {
			const _TyBaseEdge &e = m_r_system.r_Edge_Pool()[i];
			alpha = std::max(e.f_Max_VertexHessianDiagValue(), alpha);
		}
		alpha *= tau;
		//alpha = 1.44; // copy from g2o for 10khogman
		std::cout << "alfa: " << alpha << std::endl;

		if(!m_b_system_dirty)
			Refresh_Lambda(0/*m_n_verts_in_lambda*/, m_n_edges_in_lambda, alpha); // calculate only for new edges // @todo - but how to mark affected vertices?
		else
			Refresh_Lambda(0, 0, alpha); // calculate for entire system
		m_b_system_dirty = false;
		m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
		m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
		_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
			m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
			m_lambda.n_Column_Num() == n_variables_size); // lambda is square, blocks on either side = number of vertices
		// need to have lambda

		m_v_dx.resize(n_variables_size, 1);

#if 0
		{
			// test Schur
			size_t n_cameras;
			std::vector<size_t> order(m_lambda.n_BlockColumn_Num());
			n_cameras = m_schur_solver.n_Schur_Ordering(m_lambda, &order[0], order.size());
			size_t n_points = m_lambda.n_BlockColumn_Num() - n_cameras; // the rest

			Collect_RightHandSide_Vector(m_v_dx); // !!
			Eigen::VectorXd schur_sol(n_variables_size);
			bool b_cholesky_result = m_schur_solver.Schur_Solve(m_lambda,
				m_v_dx, schur_sol, &order[0], order.size(), n_cameras);

			//Collect_RightHandSide_Vector(m_v_dx);
			//bool b_cholesky_result;
			{
				Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
				b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda, v_eta); // p_dx = eta = lambda / eta
				if(m_b_verbose)
					printf("%s", (b_cholesky_result)? "Cholesky rulez!\n" : "optim success: 0\n");
			}
			// calculate cholesky

			double f_error = (m_v_dx - schur_sol).norm();
			printf("error is %g\n", f_error);
			if(n_cameras * 6 + n_points * 3 == schur_sol.rows()) { // 6 and 3 is only good for BA
				f_error = (m_v_dx.head(n_cameras * 6) - schur_sol.head(n_cameras * 6)).norm();
				printf("error in head is %g\n", f_error);
				f_error = (m_v_dx.tail(n_points * 3) - schur_sol.tail(n_points * 3)).norm();
				printf("error in tail is %g\n", f_error);
			}
			// see how precise we are
		}
#endif // 0
		// test schur complement by comparing with the result of Cholesky

		double last_error = f_Chi_Squared_Error_Denorm();
		//std::cout << "initial chi2: " << last_error << std::endl;

		//initialize vertex saver
		Eigen::VectorXd v_saved_state;
		v_saved_state.resize(n_variables_size, 1);

		if(m_b_use_schur)
			m_schur_solver.SymbolicDecomposition_Blocky(m_lambda);
		// calculate the ordering once, it does not change

		//n_max_iteration_num = 20;
		int fail = 10;
		for(size_t n_iteration = 0; n_iteration < n_max_iteration_num; ++ n_iteration) {
			//std::cout << "---------------" << std::endl;
			++ m_n_iteration_num;
			// debug

			if(m_b_verbose) {
				if(n_max_iteration_num == 1)
					printf("\n=== incremental optimization step ===\n\n");
				else
					printf("\n=== nonlinear optimization: iter #" PRIsize " ===\n\n", n_iteration);
			}
			// verbose

			if(n_iteration && m_b_system_dirty) {
				Refresh_Lambda(0, 0, alpha);
				m_b_system_dirty = false;
			}
			// no need to rebuild lambda, just refresh the values that are being referenced

			Collect_RightHandSide_Vector(m_v_dx);
			// collects the right-hand side vector

#ifdef _DEBUG
			/*_ASSERTE(m_lambda.n_Row_Num() == n_variables_size); // should be the same
			// collect errors (in parallel)

			{
				CUberBlockMatrix A;
				{
					const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
					if(!A.Append_Block(r_t_uf, 0, 0))
						throw std::bad_alloc();
					// add unary factor

					m_r_system.r_Edge_Pool().For_Each(CAlloc_JacobianBlocks(A));
					// add all the hessian blocks

					m_r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Jacobians());
					// calculate the values as well
				}
				// calculate A

				Eigen::VectorXd v_error(n_measurements_size);
				Collect_R_Errors(v_error);
				Eigen::VectorXd v_eta(n_variables_size);
				v_eta.setZero(); // eta = 0
				Eigen::VectorXd &v_b = v_error; // b = Rz * p_errors
				_ASSERTE(v_eta.rows() == n_variables_size);
				_ASSERTE(v_b.rows() == n_measurements_size);
				A.PostMultiply_Add_FBS<_TyAMatrixBlockSizes>(&v_eta(0), n_variables_size, &v_b(0),
					n_measurements_size); // works (fast parallel post-multiply)
				// calculates eta

				_ASSERTE(v_eta.rows() == m_v_dx.rows());
				for(size_t i = 0, n = v_eta.rows(); i < n; ++ i) {
					if(fabs(m_v_dx(i) - v_eta(i)) > 1e-5)
						fprintf(stderr, "error: RHS vectors do not match\n");
				}
				// compare
			}*/
			// calculate eta = A^T * b
#endif // _DEBUG

			double f_serial_start = m_timer.f_Time();

			{
				{}
				// no AtA :(

				double f_ata_end = f_serial_start;//m_timer.f_Time();

				{}
				// no mul :(

				double f_mul_end = f_ata_end;//m_timer.f_Time();

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
								// prepare symbolic decomposition, structure of lambda won't change in the next steps
								b_cholesky_result = _TyLinearSolverWrapper::Solve(m_linear_solver, m_lambda, v_eta);
								// p_dx = eta = lambda / eta
							} while(0);
						} else
							b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda, v_eta); // p_dx = eta = lambda / eta
					} else { // use Schur complement
						b_cholesky_result = m_schur_solver.Solve_PosDef_Blocky(m_lambda, v_eta);
						// Schur
					}

					if(m_b_verbose) {
						printf("%s", (b_cholesky_result)? ((m_b_use_schur)? "Schur rulez!\n" :
							"Cholesky rulez!\n") : "optim success: 0\n");
					}
				}
				// calculate cholesky, reuse block ordering if the linear solver supports it

				double f_chol_end = m_timer.f_Time();

#ifdef _DEBUG
				for(size_t i = 0; i < n_variables_size; ++ i) {
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

				double f_norm_end = m_timer.f_Time();
				double f_serial_time = f_norm_end - f_serial_start;
				m_f_ata_time += f_ata_end - f_serial_start;
				m_f_premul_time += f_mul_end - f_ata_end;
				m_f_chol_time += f_chol_end - f_mul_end;
				m_f_norm_time += f_norm_end - f_chol_end;
				m_f_serial_time += f_serial_time;
				// timing breakup

				if(f_residual_norm <= f_min_dx_norm)
					break;
				// in case the error is low enough, quit (saves us recalculating the hessians)

				//god save the vertices
				m_r_system.r_Vertex_Pool().For_Each_Parallel(CSaveState(v_saved_state));

				if(b_cholesky_result) {
					PushValuesInGraphSystem(m_v_dx);

					m_b_system_dirty = true;
				}
				// update the system (in parallel)

				if(!b_cholesky_result)
					break;
				// in case cholesky failed, quit

				double f_error = f_Chi_Squared_Error_Denorm();
				std::cout << "chi2: " << f_error << std::endl;
				if(f_error <= last_error) {
					alpha = alpha / 10; // g2o calculates scale
					last_error = f_error;
				} else {
					fprintf(stderr, "warning: chi2 rising\n");
					alpha = alpha * 10; // g2o has ni
					if(fail > 0) {
						-- fail;
						n_max_iteration_num ++;
						// restore saved vertives
						m_r_system.r_Vertex_Pool().For_Each_Parallel(CLoadState(v_saved_state));
					}
				}
				m_b_system_dirty = true;
			}
		}
	}

protected:
	/**
	 *	@brief function object that calls lambda hessian block allocation for all edges
	 */
	class CAlloc_HessianBlocks {
	protected:
		CUberBlockMatrix &m_r_lambda; /**< @brief reference to the lambda matrix (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_lambda is reference to the lambda matrix
		 */
		inline CAlloc_HessianBlocks(CUberBlockMatrix &r_lambda)
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
	 */
	inline void Collect_RightHandSide_Vector(Eigen::VectorXd &r_v_b)
	{
		m_r_system.r_Vertex_Pool().For_Each_Parallel(CCollect_RightHandSide_Vector(r_v_b)); // can do this in parallel
		// collect b
	}

	/**
	 *	@brief function object that saves state of all the vertices
	 */
	class CSaveState {
	protected:
		Eigen::VectorXd &m_r_state; /**< @brief reference to the state vector (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_state is reference to the state vector
		 */
		inline CSaveState(Eigen::VectorXd &r_state)
			:m_r_state(r_state)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyVertex is vertex type
		 *	@param[in,out] r_t_vertex is vertex to output its part the state vector
		 */
		template <class _TyVertex>
		inline void operator ()(_TyVertex &r_t_vertex)
		{
			r_t_vertex.SaveState(m_r_state);
		}
	};

	/**
	 *	@brief function object that saves state of all the vertices
	 */
	class CLoadState {
	protected:
		Eigen::VectorXd &m_r_state; /**< @brief reference to the state vector (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_state is reference to the state vector
		 */
		inline CLoadState(Eigen::VectorXd &r_state)
			:m_r_state(r_state)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyVertex is vertex type
		 *	@param[in,out] r_t_vertex is vertex to output its part the state vector
		 */
		template <class _TyVertex>
		inline void operator ()(_TyVertex &r_t_vertex)
		{
			r_t_vertex.LoadState(m_r_state);
		}
	};

#if 0
	/**
	 *	@brief function object that copies ordering from the edges
	 */
	class CGetCumsums {
	protected:
		std::vector<size_t>::iterator m_p_cumsum_it; /**< @brief cumsum iterator (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_row_cumsum_list is cumsum iterator
		 */
		inline CGetCumsums(std::vector<size_t> &r_row_cumsum_list)
			:m_p_cumsum_it(r_row_cumsum_list.begin())
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in,out] r_t_edge is edge to output its cumsum
		 */
		template <class _TyEdge>
		inline void operator ()(const _TyEdge &r_t_edge)
		{
			size_t n_order = r_t_edge.n_Order();
			_ASSERTE(n_order > 0); // don't want 0, but that is assigned to the unary factor
			*m_p_cumsum_it = n_order;
			++ m_p_cumsum_it;
		}

		/**
		 *	@brief conversion back to cumsum iterator
		 *	@return Returns the value of cumsum iterator.
		 */
		inline operator std::vector<size_t>::iterator() const
		{
			return m_p_cumsum_it;
		}
	};
#endif // 0 // unused, need to operate on edges, not on vertices

	/**
	 *	@brief creates the lambda matrix from scratch
	 */
	inline void AddEntriesInSparseSystem() // throw(std::bad_alloc)
	{
#if 0
		if(m_r_system.r_Edge_Pool().n_Size() > 1000) { // wins 2.42237 - 2.48938 = .06701 seconds on 10k.graph, likely more on larger graphs
			//printf("building large matrix from scratch ...\n"); // debug
			std::vector<size_t> row_cumsum_list(m_r_system.r_Edge_Pool().n_Size());
			/*std::vector<size_t>::iterator p_end_it =*/
				m_r_system.r_Edge_Pool().For_Each(CGetCumsums(row_cumsum_list));
			//_ASSERTE(p_end_it == row_cumsum_list.end());
			// collect cumsums

			CUberBlockMatrix tmp(row_cumsum_list.begin(),
				row_cumsum_list.end(), m_r_system.r_Vertex_Pool().n_Size());
			m_lambda.Swap(tmp);
			// use this one instead

			// todo - see if there are some row_reindex on 100k, fix it by collecting
			// cumsums and building matrix with that (proven to be faster before)
		} else
#endif // 0 // todo - need to write function that gets cumsums from vertices (it's not difficult)
		{
			//printf("building small matrix from scratch ...\n"); // debug
			m_lambda.Clear();
			// ...
		}

		const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
		if(!m_lambda.Append_Block(Eigen::MatrixXd(r_t_uf.transpose() * r_t_uf), 0, 0))
			throw std::bad_alloc();
		// add unary factor (actually UF^T * UF, but it's the same matrix)

		m_r_system.r_Edge_Pool().For_Each(CAlloc_HessianBlocks(m_lambda));
		m_r_system.r_Vertex_Pool().For_Each(CAlloc_HessianBlocks(m_lambda));
		// add all the hessian blocks

		//printf("building lambda from scratch finished\n"); // debug

		/*double alpha = 0.01; // ela?
		for(size_t i = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_last = m_lambda.n_BlockColumn_Block_Num(i) - 1;
			CUberBlockMatrix::_TyMatrixXdRef block = m_lambda.t_BlockAt(n_last, i);
			for(size_t j = 0, m = block.rows(); j < m; ++ j)
				block(j, j) += alpha;
		}*/
	}

	/**
	 *	@brief incrementally updates the lambda matrix structure (must not be empty)
	 */
	inline void UpdateSparseSystem(size_t n_skip_vertices, size_t n_skip_edges) // throw(std::bad_alloc)
	{
		_ASSERTE(m_lambda.n_Row_Num() > 0 && m_lambda.n_Column_Num() == m_lambda.n_Row_Num()); // make sure lambda is not empty
		m_r_system.r_Edge_Pool().For_Each(n_skip_edges,
			m_r_system.r_Edge_Pool().n_Size(), CAlloc_HessianBlocks(m_lambda));
		m_r_system.r_Vertex_Pool().For_Each(n_skip_vertices,
			m_r_system.r_Vertex_Pool().n_Size(), CAlloc_HessianBlocks(m_lambda));
		// add the hessian blocks of the new edges
	}

	/**
	 *	@brief incrementally updates the lambda matrix structure (can be empty)
	 */
	inline void Extend_Lambda(size_t n_vertices_already_in_lambda, size_t n_edges_already_in_lambda) // throw(std::bad_alloc)
	{
		if(!n_vertices_already_in_lambda && !n_edges_already_in_lambda)
			AddEntriesInSparseSystem(); // works for empty
		else
			UpdateSparseSystem(n_vertices_already_in_lambda, n_edges_already_in_lambda); // does not work for empty
		// create block matrix lambda

#ifdef _DEBUG
		/*{
			CUberBlockMatrix A;
			const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
			if(!A.Append_Block(r_t_uf, 0, 0))
				throw std::bad_alloc();
			// add unary factor

			m_r_system.r_Edge_Pool().For_Each(CAlloc_JacobianBlocks(A));
			// add all the hessian blocks

			CUberBlockMatrix lambda_ref;
			A.PreMultiplyWithSelfTransposeTo(lambda_ref, true); // only upper diag!
			// calculate lambda = AtA

			if(!m_lambda.b_EqualStructure(lambda_ref)) {
				lambda_ref.Rasterize("lambda1_reference_structure.tga");
				m_lambda.Rasterize("lambda0_structure.tga");
			}

			_ASSERTE(m_lambda.b_EqualStructure(lambda_ref));
			// make sure the matrix has the same structure
		}*/
#endif // _DEBUG
	}

	/**
	 *	@brief refreshes the lambda matrix by recalculating edge hessians
	 */
	inline void Refresh_Lambda(size_t n_referesh_from_vertex = 0,
		size_t n_referesh_from_edge = 0, double alpha = 0.01)
	{
		if(n_referesh_from_edge) {
			m_r_system.r_Edge_Pool().For_Each_Parallel(n_referesh_from_edge,
				m_r_system.r_Edge_Pool().n_Size(), CCalculate_Hessians());
		} else {
			m_r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Hessians());
		}
		if(n_referesh_from_vertex) {
			m_r_system.r_Vertex_Pool().For_Each_Parallel(n_referesh_from_vertex,
				m_r_system.r_Vertex_Pool().n_Size(), CCalculate_Hessians());
		} else {
			m_r_system.r_Vertex_Pool().For_Each_Parallel(CCalculate_Hessians());
		}
		// can do this in parallel

		if(!n_referesh_from_vertex) {
			const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
			m_lambda.t_FindBlock(0, 0).noalias() += r_t_uf.transpose() * r_t_uf;
		}
		// add unary factor (gets overwritten by the first vertex' block)

		// ela?
		for(size_t i = n_referesh_from_vertex, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_last = m_lambda.n_BlockColumn_Block_Num(i) - 1;
			CUberBlockMatrix::_TyMatrixXdRef block = m_lambda.t_BlockAt(n_last, i);
			for(size_t j = 0, m = block.rows(); j < m; ++ j)
				block(j, j) += alpha;
		}

#ifdef _DEBUG
		/*{
			CUberBlockMatrix A;
			const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
			if(!A.Append_Block(r_t_uf, 0, 0))
				throw std::bad_alloc();
			// add unary factor

			m_r_system.r_Edge_Pool().For_Each(CAlloc_JacobianBlocks(A));
			// add all the hessian blocks

			m_r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Jacobians());
			// calculate the values as well

			CUberBlockMatrix lambda_ref;
			A.PreMultiplyWithSelfTransposeTo(lambda_ref, true); // only upper diag!
			// calculate lambda = AtA

			if(!m_lambda.b_Equal(lambda_ref, 1e-3)) {
				m_lambda.Rasterize("lambda2_values.tga");
				lambda_ref.Rasterize("lambda3_reference_values.tga");
				CUberBlockMatrix &diff = lambda_ref;
				m_lambda.AddTo(diff, -1);
				diff.Rasterize("lambda4_diff_values.tga");
				fprintf(stderr, "error: lambda and it's reference have different value\n");
				exit(-1);
			}

			_ASSERTE(m_lambda.b_EqualStructure(lambda_ref));
			// make sure the matrix has the same structure
		}*/
#endif // _DEBUG
	}

#ifdef _DEBUG

	/**
	 *	@brief function object that calls error vector calculation for all edges
	 */
	class CCollect_R_Errors {
	protected:
		Eigen::VectorXd &m_r_R_errors; /**< @brief reference to the R error vector (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_R_errors is reference to the R error vector
		 */
		inline CCollect_R_Errors(Eigen::VectorXd &r_R_errors)
			:m_r_R_errors(r_R_errors)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in,out] r_t_edge is edge to output its part R error vector
		 */
		template <class _TyEdge>
		inline void operator ()(_TyEdge &r_t_edge) // throw(std::bad_alloc)
		{
			r_t_edge.Get_R_Error(m_r_R_errors);
		}
	};

	/**
	 *	@brief calculates the error vector
	 */
	inline void Collect_R_Errors(Eigen::VectorXd &r_v_R_error)
	{
		const Eigen::VectorXd &r_v_err = m_r_system.r_v_Unary_Error();
		r_v_R_error.segment(0, r_v_err.rows()) = r_v_err;
		// add the first error

		m_r_system.r_Edge_Pool().For_Each_Parallel(CCollect_R_Errors(r_v_R_error)); // can do this in parallel
		// collect errors
	}

	/**
	 *	@brief function object that calls hessian block allocation
	 *		for all edges (this is here for debugging purposes only)
	 */
	class CAlloc_JacobianBlocks {
	protected:
		CUberBlockMatrix &m_r_A; /**< @brief reference to the A matrix (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_A is reference to the A matrix
		 */
		inline CAlloc_JacobianBlocks(CUberBlockMatrix &r_A)
			:m_r_A(r_A)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in,out] r_t_edge is edge to have hessian blocks allocated in A
		 *	@note This function throws std::bad_alloc.
		 */
		template <class _TyEdge>
		inline void operator ()(_TyEdge &r_t_edge) // throw(std::bad_alloc)
		{
			r_t_edge.Alloc_JacobianBlocks(m_r_A);
		}
	};

	/**
	 *	@brief function object that calculates hessians
	 *		in all the edges (this is here for debugging purposes only)
	 */
	class CCalculate_Jacobians {
	public:
		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in] r_t_edge is edge to update it's hessians
		 */
		template <class _TyEdge>
		inline void operator ()(_TyEdge &r_t_edge) const
		{
			r_t_edge.Calculate_Jacobians();
		}
	};

#endif // _DEBUG

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
	class CCalculate_Hessians {
	public:
		/**
		 *	@brief function operator
		 *	@tparam _TyVertexOrEdge is vertex or edge type
		 *	@param[in] r_t_vertex_or_edge is vertex or edge to update it's hessians
		 */
		template <class _TyVertexOrEdge>
		inline void operator ()(_TyVertexOrEdge &r_t_vertex_or_edge) const
		{
			r_t_vertex_or_edge.Calculate_Hessians();
		}
	};

	CNonlinearSolver_Lambda_LM(const CNonlinearSolver_Lambda_LM &UNUSED(r_solver)); /**< @brief the object is not copyable */
	CNonlinearSolver_Lambda_LM &operator =(const CNonlinearSolver_Lambda_LM &UNUSED(r_solver)) { return *this; } /**< @brief the object is not copyable */
};

#endif // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
