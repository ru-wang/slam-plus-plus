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
#ifndef __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
#define __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED

/**
 *	@file include/slam/NonlinearSolver_Lambda.h
 *	@brief nonlinear blocky solver working above the lambda matrix
 *	@author -tHE SWINe-
 *	@date 2012-09-13
 */

#include "slam/FlatSystem.h"

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
 *	@brief enables writes of chi2 errors at each step
 */
//#define __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_DUMP_ICRA2013_ANIMATION_DATA
 *	@brief dump ICRA 2013 slam race data
 */
//#define __NONLINEAR_SOLVER_LAMBDA_DUMP_ICRA2013_ANIMATION_DATA

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
 *	@brief dump RSS 2013 matrix animation data
 */
//#define __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

//extern int n_dummy_param; // remove me

/**
 *	@brief nonlinear blocky solver working above the lambda matrix
 *
 *	@tparam CSystem is the system type
 *	@tparam CLinearSolver is a linear solver
 */
template <class CSystem, class CLinearSolver, class CAMatrixBlockSizes = typename CSystem::_TyJacobianMatrixBlockList>
class CNonlinearSolver_Lambda {
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
	double m_f_lambda_time; /**< @brief time spent updating lambda */
	double m_f_rhs_time; /**< @brief time spent gathering the right-hand-side vector */
	double m_f_chol_time; /**< @brief time spent in Choleski() section */
	double m_f_norm_time; /**< @brief time spent in norm calculation section */
	double m_f_vert_upd_time; /**< @brief time spent updating the system */

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
	CNonlinearSolver_Lambda(CSystem &r_system, size_t n_linear_solve_threshold = 0,
		size_t n_nonlinear_solve_threshold = 0, size_t n_nonlinear_solve_max_iteration_num = 5,
		double f_nonlinear_solve_error_threshold = .01, bool b_verbose = false,
		CLinearSolver linear_solver = CLinearSolver(), bool b_use_schur = true)
		:m_r_system(r_system), m_linear_solver(linear_solver), m_schur_solver(linear_solver),
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
		m_b_verbose(b_verbose), m_b_use_schur(b_use_schur), m_n_real_step(0),
		m_b_system_dirty(false), m_n_iteration_num(0), m_f_lambda_time(0), m_f_rhs_time(0),
		m_f_chol_time(0), m_f_norm_time(0), m_f_vert_upd_time(0), m_b_had_loop_closure(false)
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
		printf("solver spent %f seconds in parallelizable section (disparity %f seconds)\n",
			m_f_lambda_time + m_f_rhs_time, (f_total_time > 0)? f_total_time -
			(m_f_lambda_time + m_f_rhs_time + m_f_chol_time + m_f_norm_time + m_f_vert_upd_time) : 0);
		printf("out of which:\n");
		printf("\t   ,\\: %f\n", m_f_lambda_time);
		printf("\t  rhs: %f\n", m_f_rhs_time);
		printf("solver spent %f seconds in serial section\n",
			m_f_chol_time + m_f_norm_time + m_f_vert_upd_time);
		printf("out of which:\n");
		printf("\t chol: %f\n", m_f_chol_time);
		printf("\t norm: %f\n", m_f_norm_time);
		printf("\tv-upd: %f\n", m_f_vert_upd_time);
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

		return m_lambda.Rasterize_Symmetric(p_s_filename, n_scalar_size);
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

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
		n_new_vertex_num = m_n_nonlinear_solve_threshold;
		m_b_had_loop_closure = true;
		// make it trigger
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

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
				Optimize((m_b_had_loop_closure)?
					m_n_nonlinear_solve_max_iteration_num : 0,
					m_f_nonlinear_solve_error_threshold);
				m_b_had_loop_closure = false;
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
			// dump chi2

            //printf("time is %lf\n", m_timer.f_Time()); // print incremental times (ICRA 2013)
#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_ICRA2013_ANIMATION_DATA
			char p_s_filename[256];
			sprintf(p_s_filename, "icra2013/sys%05" _PRIsize ".txt", m_r_system.r_Vertex_Pool().n_Size());
			m_r_system.Dump(p_s_filename);
			//sprintf(p_s_filename, "icra2013/sys%05" _PRIsize ".tga", m_r_system.r_Vertex_Pool().n_Size());
			//m_r_system.Plot2D(p_s_filename); // no plot, takes a ton of memory
			// dump vertices for the ICRA animation
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_ICRA2013_ANIMATION_DATA
		}
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
	}

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
	CMatrixOrdering mord100;
	CUberBlockMatrix lambda_prev;
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

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
			fprintf(stderr, "warning: the system is underspecified\n");
			return;
		}
		if(!n_measurements_size)
			return; // nothing to solve (but no results need to be generated so it's ok)
		// can't solve in such conditions
		Extend_Lambda(m_n_verts_in_lambda, m_n_edges_in_lambda); // recalculated all the jacobians inside Extend_Lambda()
		if(!m_b_system_dirty) {
			Refresh_Lambda(0/*m_n_verts_in_lambda*/, m_n_edges_in_lambda); // calculate only for new edges // @todo - but how to mark affected vertices?
		} else {
			Refresh_Lambda(); // calculate for entire system
		}
		m_b_system_dirty = false;
		m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
		m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
		_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
			m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
			m_lambda.n_Column_Num() == n_variables_size); // lambda is square, blocks on either side = number of vertices
		// need to have lambda

		//Dump_SystemMatrix("lambda_system.tga");

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
		char p_s_filename[256];

		sprintf(p_s_filename, "rss2013/%05d_0_lambda.tga", m_n_edges_in_lambda);
		m_lambda.Rasterize_Symmetric(p_s_filename); // this lambda
		sprintf(p_s_filename, "rss2013/%05d_1_lambda2.tga", m_n_edges_in_lambda);
		m_lambda.Rasterize_Symmetric(lambda_prev, true, p_s_filename); // with changes marked, not symmetric

		CUberBlockMatrix Lprev;
		Lprev.CholeskyOf(lambda_prev);
		// do it now, will damage lambda_prev

		_TyBaseEdge &r_last_edge = m_r_system.r_Edge_Pool()[m_r_system.r_Edge_Pool().n_Size() - 1];
		size_t v0 = r_last_edge.n_Vertex_Id(0);
		size_t v1 = r_last_edge.n_Vertex_Id(1);
		if(lambda_prev.n_BlockColumn_Num()) {
			size_t v = std::min(v0, v1);
			lambda_prev.t_GetBlock_Log(v, v).setZero(); // make sure this is clear
		}

		sprintf(p_s_filename, "rss2013/%05d_2_lambda3.tga", m_n_edges_in_lambda);
		m_lambda.CopyTo(lambda_prev);
		lambda_prev.Scale(0); // zero out the whole thing
		m_lambda.Rasterize_Symmetric(lambda_prev, true, p_s_filename); // all changed, not symmetric

		CUberBlockMatrix L;
		CUberBlockMatrix Lem;
		L.CopyLayoutTo(Lem);
		L.CholeskyOf(m_lambda);
		sprintf(p_s_filename, "rss2013/%05d_3_Lnoord.tga", m_n_edges_in_lambda);
		L.Rasterize(p_s_filename); // no fill-in highlighting
		sprintf(p_s_filename, "rss2013/%05d_4_Lnoord_fill-in.tga", m_n_edges_in_lambda);
		L.Rasterize(m_lambda, false, p_s_filename); // highlight fill-in
		sprintf(p_s_filename, "rss2013/%05d_a_Lnoord_red.tga", m_n_edges_in_lambda);
		L.Rasterize(Lem, false, p_s_filename); // highlight fill-in
		sprintf(p_s_filename, "rss2013/%05d_b_Lnoord_inc.tga", m_n_edges_in_lambda);
		L.Rasterize(Lprev, true, p_s_filename); // highlight fill-in

		m_lambda.CopyTo(lambda_prev);
		// copy the lambda matrix

		std::vector<size_t> cumsum, ccumsum;
		for(size_t i = 0, n_csum = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			n_csum += m_lambda.n_BlockColumn_Column_Num(i);
			cumsum.push_back(n_csum);
		}
		size_t p_cumsum[1] = {1};
		ccumsum.insert(ccumsum.begin(), p_cumsum, p_cumsum + 1);
		CUberBlockMatrix rhs(cumsum.begin(), cumsum.end(), ccumsum.begin(), ccumsum.end());
		m_v_dx.resize(n_variables_size, 1);
		Collect_RightHandSide_Vector(m_v_dx);
		for(size_t i = 0, n_csum = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_w = m_lambda.n_BlockColumn_Column_Num(i);
			if(!rhs.Append_Block(m_v_dx.segment(n_csum, n_w), n_csum, 0))
				throw std::runtime_error("segment size fail");
			n_csum += n_w;
		}
		sprintf(p_s_filename, "rss2013/%05d_5_rhs.tga", m_n_edges_in_lambda);
		rhs.Rasterize(p_s_filename);

		CUberBlockMatrix rhs_nonzero, rhs_hl;
		rhs.CopyLayoutTo(rhs_nonzero);
		rhs.CopyLayoutTo(rhs_hl);
		m_v_dx.setConstant(1);
		for(size_t i = 0, n_csum = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_w = m_lambda.n_BlockColumn_Column_Num(i);
			if(!(i == v0 || i == v1)) {
				if(!rhs_hl.Append_Block(m_v_dx.segment(n_csum, n_w), n_csum, 0))
					throw std::runtime_error("segment size fail");
			}
			n_csum += n_w;
		}
		for(size_t i = 0, n_csum = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_w = m_lambda.n_BlockColumn_Column_Num(i);
			if(!rhs_nonzero.Append_Block(m_v_dx.segment(n_csum, n_w), n_csum, 0))
				throw std::runtime_error("segment size fail");
			n_csum += n_w;
		}
		sprintf(p_s_filename, "rss2013/%05d_9_rhs_add.tga", m_n_edges_in_lambda);
		rhs_nonzero.Rasterize(rhs_hl, true, p_s_filename);
		rhs.Scale(0);
		sprintf(p_s_filename, "rss2013/%05d_6_rhs_red.tga", m_n_edges_in_lambda);
		rhs_nonzero.Rasterize(rhs, true, p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_7_rhs_nnz.tga", m_n_edges_in_lambda);
		rhs_nonzero.Rasterize(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_8_rhs_z.tga", m_n_edges_in_lambda);
		rhs.Rasterize(p_s_filename);

#if 0
		/*sprintf(p_s_filename, "rss2013/%05d_0_lambda.tga", m_n_verts_in_lambda);
		m_lambda.Rasterize_Symmetric(p_s_filename);*/ // probably don't need lambdas

		CUberBlockMatrix lambda_perm;

		CMatrixOrdering mord;
		mord.p_BlockOrdering(m_lambda, true);
		const size_t *p_perm = mord.p_Get_InverseOrdering();

		if(mord100.n_Ordering_Size() < m_lambda.n_BlockColumn_Num()) {
			if(!(m_lambda.n_BlockColumn_Num() % 100))
				mord100.p_BlockOrdering(m_lambda, true); // do full every 100
			else {
				mord100.p_InvertOrdering(mord100.p_ExtendBlockOrdering_with_Identity(
					m_lambda.n_BlockColumn_Num()), m_lambda.n_BlockColumn_Num());
			}
		}
		// maintain an "every 100" ordering as well

		m_lambda.Permute_UpperTriangular_To(lambda_perm, p_perm, mord.n_Ordering_Size(), true);
		sprintf(p_s_filename, "rss2013/%05d_1_lambda-perm.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
		//	lambda_perm.Rasterize_Symmetric(p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 3); // do not really need lambdas right now

		L.CholeskyOf(lambda_perm);
		size_t n_L_good_blocks = L.n_Block_Num();
		sprintf(p_s_filename, "rss2013/%05d_2_L.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			L.Rasterize(lambda_perm, false, p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 3 : 2); // highlight fill-in

		L.CholeskyOf(m_lambda);
		sprintf(p_s_filename, "rss2013/%05d_9_Lnoord.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			L.Rasterize(lambda_perm, false, p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 2); // highlight fill-in

		const size_t *p_perm100 = mord100.p_Get_InverseOrdering();
		m_lambda.Permute_UpperTriangular_To(lambda_perm, p_perm100, mord100.n_Ordering_Size(), true);
		sprintf(p_s_filename, "rss2013/%05d_4_lambda-perm-100.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
		//	lambda_perm.Rasterize_Symmetric(p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 3); // do not really need lambdas right now

		L.CholeskyOf(lambda_perm);
		sprintf(p_s_filename, "rss2013/%05d_5_L100.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			L.Rasterize(lambda_perm, false, p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 3 : 2); // highlight fill-in

		sprintf(p_s_filename, "rss2013/%05d_3_stats.txt", m_n_verts_in_lambda);
		FILE *p_fw;
		if((p_fw = fopen(p_s_filename, "w"))) {
			fprintf(p_fw, PRIsize "\n", lambda_perm.n_Block_Num()); // save density of lambda
			fprintf(p_fw, PRIsize "\n", n_L_good_blocks); // save density of L
			fprintf(p_fw, PRIsize "\n", L.n_Block_Num()); // save density of L at every100
			fclose(p_fw);
		}
#endif // 0
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

#ifdef __ICRA_GET_PRESENTATION_ANIMATION_DATA
		char p_s_filename[256];
		sprintf(p_s_filename, "lambda_%03d.tga", m_n_real_step);
		m_lambda.Rasterize(p_s_filename);
		// dump lambda

		sprintf(p_s_filename, "lambda_%03d.txt", m_n_real_step);
		cs *p_lam = m_lambda.p_BlockStructure_to_Sparse();
		FILE *p_fw = fopen(p_s_filename, "w");
		fprintf(p_fw, "%d x %d\n", p_lam->m, p_lam->n);

		{
			cs *m = p_lam;
			int cols = m->n;
			int rows = m->m;
			int size_i = m->p[cols];
			int col = 0;
			for(int i = 0; i < size_i; ++ i) { // 8x if, 8x add
				int row = m->i[i]; // 7x read
				while(m->p[col + 1] <= i) // 7x if / 10x if (in case of while), the same amount of reads
					col ++; // 4x add
				fprintf(p_fw, "%d, %d\n", row, col);
			}
			// 14 - 17x read, 12x add, 15 - 18x if
		}
		fclose(p_fw);
		sprintf(p_s_filename, "lambda_%03d_ref.txt", m_n_real_step);
		p_fw = fopen(p_s_filename, "w");
		fprintf(p_fw, "%d x %d\n", p_lam->m, p_lam->n);
		{
			cs *m = p_lam;
			int cols = m->n;
			int rows = m->m;
			for(int j = 0; j < cols; ++ j) { // 5x if, 5x add
				int col = j;
				int first_i = m->p[j]; // 4x read
				int last_i = m->p[j + 1]; // 4x read
				for(int i = first_i; i < last_i; ++ i) { // 7+4=11x if, 11x add
					int row = m->i[i]; // 7x read
					fprintf(p_fw, "%d, %d\n", row, col);
				}
			}
			// 15x read, 16x add, 16x if
		}
		fclose(p_fw);
		cs_spfree(p_lam);
#endif // __ICRA_GET_PRESENTATION_ANIMATION_DATA

		m_v_dx.resize(n_variables_size, 1);

		if(m_b_use_schur)
			m_schur_solver.SymbolicDecomposition_Blocky(m_lambda);
		// calculate the ordering once, it does not change

		for(size_t n_iteration = 0; n_iteration < n_max_iteration_num; ++ n_iteration) {

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
				Refresh_Lambda();
				m_b_system_dirty = false;
			}
			// no need to rebuild lambda, just refresh the values that are being referenced

			timer.Accum_DiffSample(m_f_lambda_time);

			Collect_RightHandSide_Vector(m_v_dx);
			// collects the right-hand side vector

			timer.Accum_DiffSample(m_f_rhs_time);

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

			{
				bool b_cholesky_result;
				{
					Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
					if(!m_b_use_schur) {
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
						printf("%s %s", (m_b_use_schur)? "Schur" : "Cholesky",
							(b_cholesky_result)? "succeeded\n" : "failed\n");
					}
				}
				// calculate cholesky, reuse block ordering if the linear solver supports it

				timer.Accum_DiffSample(m_f_chol_time);

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

				timer.Accum_DiffSample(m_f_norm_time);

				if(f_residual_norm <= f_min_dx_norm)
					break;
				// in case the error is low enough, quit (saves us recalculating the hessians)

				if(b_cholesky_result) {
					PushValuesInGraphSystem(m_v_dx);
					m_b_system_dirty = true;

					timer.Accum_DiffSample(m_f_vert_upd_time);
				}
				// update the system (in parallel)

				if(!b_cholesky_result)
					break;
				// in case cholesky failed, quit
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
	inline void Refresh_Lambda(size_t n_referesh_from_vertex = 0, size_t n_refresh_from_edge = 0)
	{
		if(n_refresh_from_edge) {
			m_r_system.r_Edge_Pool().For_Each_Parallel(n_refresh_from_edge,
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

	CNonlinearSolver_Lambda(const CNonlinearSolver_Lambda &UNUSED(r_solver)); /**< @brief the object is not copyable */
	CNonlinearSolver_Lambda &operator =(const CNonlinearSolver_Lambda &UNUSED(r_solver)) { return *this; } /**< @brief the object is not copyable */
};

#endif // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
