/*
								+-----------------------------------+
								|                                   |
								|   ***  A non-linear solver  ***   |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2012  |
								|                                   |
								|        NonlinearSolver_A.h        |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED
#define __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED

/**
 *	@file include/slam/NonlinearSolver_A.h
 *	@brief nonlinear blocky solver working above the A matrix
 *	@author -tHE SWINe-
 *	@date 2012-09-03
 */

#include "slam/FlatSystem.h"

/**
 *	@brief nonlinear blocky solver working above the A matrix
 *
 *	@tparam CSystem is the system type
 *	@tparam CLinearSolver is a linear solver
 */
template <class CSystem, class CLinearSolver, class CAMatrixBlockSizes = typename CSystem::_TyJacobianMatrixBlockList>
class CNonlinearSolver_A {
public:
	typedef CSystem _TySystem; /**< @brief sy	stem type */
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

	CUberBlockMatrix m_A; /**< @brief the A matrix (built incrementally) */
	Eigen::VectorXd m_v_error; /**< @brief error vector */
	Eigen::VectorXd m_v_dx; /**< @brief dx vector */
	size_t m_n_edges_in_A; /**< @brief number of edges already in A */
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
	CNonlinearSolver_A(CSystem &r_system, size_t n_linear_solve_threshold,
		size_t n_nonlinear_solve_threshold, size_t n_nonlinear_solve_max_iteration_num,
		double f_nonlinear_solve_error_threshold, bool b_verbose,
		CLinearSolver linear_solver = CLinearSolver(), bool b_use_schur = true)
		:m_r_system(r_system), m_linear_solver(linear_solver),
		m_schur_solver(linear_solver), m_n_edges_in_A(0),
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
		m_b_system_dirty(false), m_n_iteration_num(0), m_f_serial_time(0), m_f_ata_time(0),
		m_f_premul_time(0), m_f_chol_time(0), m_f_norm_time(0), m_b_had_loop_closure(false)
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
			printf("solver spent %f seconds in parallelizable section (updating A)\n", f_total_time - m_f_serial_time);
		printf("solver spent %f seconds in serial section\n", m_f_serial_time);
		printf("out of which:\n");
		printf("\t  ata: %f\n", m_f_ata_time);
		printf("\tgaxpy: %f\n", m_f_premul_time);
		printf("\t chol: %f\n", m_f_chol_time);
		printf("\t norm: %f\n", m_f_norm_time);
		printf("\ttotal: %f\n", m_f_ata_time + m_f_premul_time + m_f_chol_time + m_f_norm_time);
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
			Extend_A(m_n_edges_in_A); // recalculated all the jacobians inside Extend_A()
			if(!m_b_system_dirty)
				Refresh_A(m_n_edges_in_A); // calculate only for new edges
			else
				Refresh_A(); // calculate for entire system
			m_b_system_dirty = false;
			m_n_edges_in_A = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
			_ASSERTE(m_A.n_Row_Num() >= m_A.n_Column_Num()); // no, A is *not* square, but size of measurements >= vertices
			// need to have A
		} catch(std::bad_alloc&) {
			return false;
		}

		return m_A.Rasterize(p_s_filename, n_scalar_size);
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
		size_t n_rows = m_A.n_Row_Num();
		size_t n_columns = m_A.n_Column_Num();
		size_t n_nonzeros = m_A.n_NonZero_Num();
		FILE *p_fw;
		if(!(p_fw = fopen(p_s_filename, "w")))
			return false;
		fprintf(p_fw, "%%%%MatrixMarket matrix coordinate real general\n");
		fprintf(p_fw, "%%-------------------------------------------------------------------------------\n");
		fprintf(p_fw, "%% SLAM_plus_plus matrix dump\n");
		fprintf(p_fw, "%% kind: A matrix for SLAM problem\n");
		fprintf(p_fw, "%%-------------------------------------------------------------------------------\n");
		fprintf(p_fw, "" PRIsize " " PRIsize " " PRIsize "\n", n_rows, n_columns, n_nonzeros);
		for(size_t n_col = 0, n_column_blocks = m_A.n_BlockColumn_Num();
		   n_col < n_column_blocks; ++ n_col) {
			size_t n_col_base = m_A.n_BlockColumn_Base(n_col);
			size_t n_col_width = m_A.n_BlockColumn_Column_Num(n_col);
			size_t n_block_num = m_A.n_BlockColumn_Block_Num(n_col);
			for(size_t j = 0; j < n_block_num; ++ j) {
				size_t n_row = m_A.n_Block_Row(n_col, j);
				size_t n_row_base = m_A.n_BlockRow_Base(n_row);
				size_t n_row_height = m_A.n_BlockRow_Row_Num(n_row);

				CUberBlockMatrix::_TyMatrixXdRef t_block = m_A.t_BlockAt(j, n_col);
				_ASSERTE(t_block.rows() == n_row_height && t_block.cols() == n_col_width);
				// get a block

				for(size_t k = 0; k < n_row_height; ++ k) {
					for(size_t l = 0; l < n_col_width; ++ l) {
						fprintf(p_fw, "" PRIsize " " PRIsize " %lf\n", n_row_base + 1 + k,
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
			m_b_had_loop_closure = !(n_first_vertex == n_vertex_num - 2);
			_ASSERTE(m_b_had_loop_closure || std::max(r_last_edge.n_Vertex_Id(0),
				r_last_edge.n_Vertex_Id(1)) == n_vertex_num - 1);
		}
		// detect loop closures (otherwise the edges are initialized based on measurement and error would be zero)

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
		}
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

		Extend_A(m_n_edges_in_A); // recalculated all the jacobians inside Extend_A()
		if(!m_b_system_dirty)
			Refresh_A(m_n_edges_in_A); // calculate only for new edges
		else
			Refresh_A(); // calculate for entire system
		m_b_system_dirty = false;
		m_n_edges_in_A = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
		_ASSERTE(m_A.n_Row_Num() >= m_A.n_Column_Num()); // no, A is *not* square, but size of measurements >= vertices
		// need to have A

		//Dump_SystemMatrix("A_system.tga");

		m_v_error.resize(n_measurements_size, 1);
		m_v_dx.resize(n_variables_size, 1);

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
				Refresh_A();
				m_b_system_dirty = false;
			}
			// no need to rebuild A, just refresh the values that are being referenced

			_ASSERTE(m_A.n_Row_Num() == n_measurements_size); // should be the same
			Collect_R_Errors(m_v_error);
			// collect errors (in parallel)

			double f_serial_start = m_timer.f_Time();

			//MakeTypelist_Safe((Eigen::Matrix3d, Eigen::Matrix2d, Eigen::Matrix<double, 2, 3>)) pose_landmark_2d_mats; // yes, instance

			{
				CUberBlockMatrix lambda;
#ifdef _OPENMP
				/*CUberBlockMatrix At, lambdat;
				m_A.TransposeTo(At);
				At.PreMultiplyWithSelfTransposeTo_FBS_Parallel(lambdat, _TyAMatrixBlockSizes(), true);
				lambdat.TransposeTo(lambda);*/
				// does not work, the matrix doesn't have correct size. need to modify the multiplication itself (would have to have PostMultiplyWith() ...) // todo - implement that in free time

				//if(m_A.n_Column_Block_Num() >= 50)
					m_A.PreMultiplyWithSelfTransposeTo_FBS_Parallel<_TyAMatrixBlockSizes>(lambda, true);
				//else
				//	m_A.PreMultiplyWithSelfTransposeTo_FBS(lambda, _TyAMatrixBlockSizes(), true);
#else // _OPENMP
				m_A.PreMultiplyWithSelfTransposeTo_FBS<_TyAMatrixBlockSizes>(lambda, true);
#endif // _OPENMP
				// calculate lambda without calculating At (and only upper diagonal thereof)

				double f_ata_end = m_timer.f_Time();

				{
					Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
					v_eta.setZero(); // eta = 0
					Eigen::VectorXd &v_b = m_v_error; // b = Rz * p_errors
					_ASSERTE(v_eta.rows() == n_variables_size);
					_ASSERTE(v_b.rows() == n_measurements_size);

					//m_A.PostMultiply_Add(&v_eta(0), n_variables_size, &v_b(0), n_measurements_size); // eta^T += b^T * A
					// works (default)

					/*CUberBlockMatrix At;
					m_A.TransposeTo(At);
					At.PreMultiply_Add_FBS(&v_eta(0), n_variables_size, &v_b(0),
						n_measurements_size, transpose_pose_landmark_2d_mats);*/ // eta = A^T * b
					// works (slow pre-multiply with transpose) // t_odo - test it // works

					m_A.PostMultiply_Add_FBS_Parallel<_TyAMatrixBlockSizes>(&v_eta(0), n_variables_size, &v_b(0),
						n_measurements_size); // works (fast parallel post-multiply)
					// use the function with fixed-block sizes
				}
				// calculate eta = A^T * b

				double f_mul_end = m_timer.f_Time();

				bool b_cholesky_result;
				{
					Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
					if(!m_b_use_schur) {
						if(n_max_iteration_num > 1) {
							do {
								if(!n_iteration &&
								   !_TyLinearSolverWrapper::FinalBlockStructure(m_linear_solver, lambda)) {
									b_cholesky_result = false;
									break;
								}
								// prepare symbolic decomposition, structure of lambda won't change in the next steps

								b_cholesky_result = _TyLinearSolverWrapper::Solve(m_linear_solver, lambda, v_eta);
								// p_dx = eta = lambda / eta
							} while(0);
						} else
							b_cholesky_result = m_linear_solver.Solve_PosDef(lambda, v_eta); // p_dx = eta = lambda / eta
					} else { // use Schur complement
						if(!n_iteration)
							m_schur_solver.SymbolicDecomposition_Blocky(lambda);
						b_cholesky_result = m_schur_solver.Solve_PosDef_Blocky(lambda, v_eta);
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

				if(b_cholesky_result) {
					PushValuesInGraphSystem(m_v_dx);
					m_b_system_dirty = true;

					// t_odo - create more efficient system; use FAP to keep vertices and edges (doesn't change addresses, can quickly iterate)
					// t_odo - keep in mind that the solver will have some matrix dimensions fixed, do that the way g2o guys do // did it better
					// t_odo - design the system so that vertices and edges are c++ objects? might present problems with allocators, but still ... might be nice to be able to template stuff away
					// t_odo - keep in mind that we might want to only put stuff in A (make refs to blocks), not to carry out actual calculation (makes management of m_b_system_dirty complex)
					// t_odo - *parallelize* PushValuesInGraphSystem() and p_Collect_R_Errors() and CalcJacobiansAndErrors(), if we can parallelize those, we're missing only about 10 seconds to be the best

					// t_odo - mark vertices affected by p_dx, recalculate error and jacobians in those only. // bad idea: all are affected (would have to threshold)
					// t_odo - find out an efficient way of mapping indices of different granularity to vertices // don't know what that means anymore
					// t_odo - find out how to calculate the set of edges that reference the vertices (could use per-vertex sets and merge those. don't want to use trees, think about sorted arrays) // kept in vertex, no need to sort it
					// t_odo - update lambda and at incrementaly (or use at = transpose(chol(lambda))?) // could also keep and update b // this is offtopic now
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
	 *	@brief function object that calls hessian block allocation for all edges
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
		inline void operator ()(_TyEdge &r_t_edge)
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

	/**
	 *	@brief creates the A matrix from scratch
	 */
	inline void AddEntriesInSparseSystem() // throw(std::bad_alloc)
	{
		if(m_r_system.r_Edge_Pool().n_Size() > 1000) { // wins 2.42237 - 2.48938 = .06701 seconds on 10k.graph, likely more on larger graphs
			//printf("building large matrix from scratch ...\n"); // debug
			std::vector<size_t> row_cumsum_list(m_r_system.r_Edge_Pool().n_Size());
			/*std::vector<size_t>::iterator p_end_it =*/
				m_r_system.r_Edge_Pool().For_Each(CGetCumsums(row_cumsum_list));
			//_ASSERTE(p_end_it == row_cumsum_list.end());
			// collect cumsums

			CUberBlockMatrix tmp(row_cumsum_list.begin(),
				row_cumsum_list.end(), m_r_system.r_Vertex_Pool().n_Size());
			m_A.Swap(tmp);
			// use this one instead

			// todo - see if there are some row_reindex on 100k, fix it by collecting
			// cumsums and building matrix with that (proven to be faster before)
		} else {
			//printf("building small matrix from scratch ...\n"); // debug
			m_A.Clear();
			// ...
		}

		const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
		if(!m_A.Append_Block(r_t_uf, 0, 0))
			throw std::bad_alloc();
		// add unary factor

		m_r_system.r_Edge_Pool().For_Each(CAlloc_JacobianBlocks(m_A));
		// add all the hessian blocks

		//printf("building A from scratch finished\n"); // debug
	}

	/**
	 *	@brief incrementally updates the A matrix structure (must not be empty)
	 *	@param[in] n_skip_edges is number of edges, already in the system
	 *	@note This function throws std::bad_alloc.
	 */
	inline void UpdateSparseSystem(size_t n_skip_edges) // throw(std::bad_alloc)
	{
		_ASSERTE(m_A.n_Row_Num() > 0 && m_A.n_Column_Num() > 0); // make sure A is not empty
		m_r_system.r_Edge_Pool().For_Each(n_skip_edges, m_r_system.r_Edge_Pool().n_Size(), CAlloc_JacobianBlocks(m_A));
		// add the hessian blocks of the new edges
	}

	/**
	 *	@brief incrementally updates the A matrix structure (can be empty)
	 *	@param[in] n_edges_already_in_A is number of edges, already in the system
	 *	@note This function throws std::bad_alloc.
	 */
	inline void Extend_A(size_t n_edges_already_in_A) // throw(std::bad_alloc)
	{
		if(!n_edges_already_in_A)
			AddEntriesInSparseSystem();
		else
			UpdateSparseSystem(n_edges_already_in_A);
		// create block matrix A
	}

	/**
	 *	@brief refreshes the A matrix by recalculating edge hessians
	 *	@param[in] n_refresh_from_edge is zero-based index of the first edge to be refreshed
	 */
	inline void Refresh_A(size_t n_refresh_from_edge = 0)
	{
		if(n_refresh_from_edge) {
			m_r_system.r_Edge_Pool().For_Each_Parallel(n_refresh_from_edge,
				m_r_system.r_Edge_Pool().n_Size(), CCalculate_Jacobians());
		} else
			m_r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Jacobians());
		// can do this in parallel
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

	CNonlinearSolver_A(const CNonlinearSolver_A &UNUSED(r_solver)); /**< @brief the object is not copyable */
	CNonlinearSolver_A &operator =(const CNonlinearSolver_A &UNUSED(r_solver)) { return *this; } /**< @brief the object is not copyable */
};

#endif // __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED
