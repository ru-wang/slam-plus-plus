/*
								+---------------------------------+
								|                                 |
								| ***  Marginals calculation  *** |
								|                                 |
								| Copyright  (c) -tHE SWINe- 2013 |
								|                                 |
								|           Marginals.h           |
								|                                 |
								+---------------------------------+
*/

#pragma once
#ifndef __MARGINALS_INCLUDED
#define __MARGINALS_INCLUDED

/**
 *	@file include/slam/Marginals.h
 *	@date 2013
 *	@author -tHE SWINe-
 *	@brief calculation of marginal covariances in sparse systems
 */

#include "slam/BlockMatrix.h"
#include "slam/Timer.h"
#include "eigen/Eigen/Core"

/**
 *	@brief prototype marginal covariance methods implementation
 */
class CMarginals {
protected:
public:
	static void Calculate_DenseMarginals_Ref(Eigen::MatrixXd &r_marginals, const CUberBlockMatrix &r_R)
	{
		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.resize(n, n);
		Eigen::MatrixXd &R_inv = r_marginals; // R_inv = S
		for(size_t i = 0, n_block_col = -1, n_col_remains = 1; i < n; ++ i) {
			double *p_column = &R_inv.col(i)(0);
			memset(p_column, 0, n * sizeof(double));
			p_column[i] = 1;
			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = r_R.n_BlockColumn_Column_Num(++ n_block_col);
			r_R.UpperTriangular_Solve(p_column, n, n_block_col); // backsub, only the nonzero part of the column (started at (block) column which contains column i, with no loss of generality)
		}
		r_marginals = R_inv * R_inv.transpose(); // C = SS^T, might need some memory
		// calculate the covariance (assume that this is correct)
	}

	static void Calculate_DenseMarginals_Slow(Eigen::MatrixXd &r_marginals, const CUberBlockMatrix &r_R)
	{
		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.resize(n, n);
		r_marginals.setZero(); // !!

		Eigen::MatrixXd R_inv_column(n, 1); // R_inv = S
		double *p_column = &R_inv_column.col(0)(0);
		// get dense column data from the Eigen matrix (actually only need one)

		for(size_t i = 0, n_block_col = -1, n_col_remains = 1; i < n; ++ i) {
			memset(p_column, 0, n * sizeof(double));
			p_column[i] = 1;
			// make a column vector with a single 1 in it

			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = r_R.n_BlockColumn_Column_Num(++ n_block_col);
			// if done before, it avoids both referencing block col 0 in an empty matrix
			// and determining no. of columns in one past the last block column at the end

			size_t UNUSED(n_block_column_size);
			size_t n_block_column = r_R.n_Find_BlockColumn(i, n_block_column_size);
			_ASSERTE(n_block_col == n_block_column); // should be the same
			// get which block column contains column i (optimize this away, probably need to use it when resuming)

			//_ASSERTE(n_block_column_size <= n_diag_band_width); // make this into a run-time check in the production code
			// make sure it is not longer than the diagonal (otherwise we will not have enough backlog to calculate all the off-diagonal elements)

			r_R.UpperTriangular_Solve(p_column, n, n_block_col); // backsub, only the nonzero part of the column (started at (block) column which contains column i, with no loss of generality)
			// this seems to be O(i) divisions + O(nnz) MADs in the given (block) column range
			// that sums up to O(n^2/2) divisions + O(nnz log(nnz))?? MADs ... some quadratute anyways

			std::vector<double> backsub_test(n, 0);
			backsub_test[i] = 1; // !!
			r_R.UpperTriangular_Solve(&backsub_test[0], n); // full backsub
			_ASSERTE(!memcmp(p_column, &backsub_test[0], n * sizeof(double)));
			// make sure that the result is correct

			_ASSERTE((Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_column + i + 1, n - i - 1).norm() == 0)); // double pars required because of the comma in Map params
			// everything below i is zero (true)

			for(size_t k = 0; k <= i; ++ k) {
				for(size_t j = 0; j <= i; ++ j)
					r_marginals(j, k) += p_column[j] * p_column[k]; // it is symmetric, indexing arbitrary
			}
			// accumulate the entries of the covariace matrix. this is O(n^3/2) MADs for the full matrix
			// note that to calculate even only the diagonal, we need full columns
		}
	}

	static void Calculate_DenseMarginals_Fast(Eigen::MatrixXd &r_marginals,
		const CUberBlockMatrix &r_R)
	{
		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.resize(n, n);
		size_t n_block_col = -1, n_col_remains = 1;
		_ASSERTE(n <= INT_MAX);
		int _n = int(n);
		//#pragma omp parallel for // todo: need to redo the indexing logic
		for(int i = 0; i < _n; ++ i) {
			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = r_R.n_BlockColumn_Column_Num(++ n_block_col);
			// if done before, it avoids both referencing block col 0 in an empty matrix
			// and determining no. of columns in one past the last block column at the end

			double *p_column = &r_marginals.col(i)(0);
			memset(p_column, 0, n * sizeof(double));
			p_column[i] = 1;
			// make a column vector with a single 1 in it

			r_R.UpperTriangularTranspose_Solve(p_column, n, n_block_col);
			r_R.UpperTriangular_Solve(p_column, n/*, n_block_col*/);
			// solve for the whole column thing, generates one column at a time
		}
	}

	static void Calculate_DenseMarginals_LastNColumns_Fast(Eigen::MatrixXd &r_marginals,
		const CUberBlockMatrix &r_R, size_t n_column_num)
	{
		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.resize(n, n);
		for(size_t i = 0, n_block_col = -1, n_col_remains = 1; i < n; ++ i) {
			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = r_R.n_BlockColumn_Column_Num(++ n_block_col);
			// if done before, it avoids both referencing block col 0 in an empty matrix
			// and determining no. of columns in one past the last block column at the end

			if(n > n_column_num - 1 && i < n - n_column_num - 1)
				continue;
			// skip the prefix columns (todo - just call r_R.n_Find_BlockColumn())

			double *p_column = &r_marginals.col(i)(0);
			memset(p_column, 0, n * sizeof(double));
			p_column[i] = 1;
			// make a column vector with a single 1 in it

			r_R.UpperTriangularTranspose_Solve(p_column, n, n_block_col);
			r_R.UpperTriangular_Solve(p_column, n/*, n_block_col*/);
			// solve for the whole column thing, generates one column at a time
		}
	}

	static void Calculate_DenseMarginals_LastNColumns_Slow(Eigen::MatrixXd &r_marginals,
		const CUberBlockMatrix &r_R, size_t n_column_num)
	{
		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.resize(n, n);
		r_marginals.setZero(); // !!

		Eigen::MatrixXd R_inv_column(n, 1); // R_inv = S
		double *p_column = &R_inv_column.col(0)(0);
		// get dense column data from the Eigen matrix (actually only need one)

		for(size_t i = 0, n_block_col = -1, n_col_remains = 1; i < n; ++ i) {
			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = r_R.n_BlockColumn_Column_Num(++ n_block_col);
			// if done before, it avoids both referencing block col 0 in an empty matrix
			// and determining no. of columns in one past the last block column at the end

			if(n > n_column_num - 1 && i < n - n_column_num - 1)
				continue;
			// skip the prefix columns (todo - just call

			memset(p_column, 0, n * sizeof(double));
			p_column[i] = 1;
			// make a column vector with a single 1 in it

			size_t UNUSED(n_block_column_size);
			_ASSERTE(n_block_col == r_R.n_Find_BlockColumn(i, n_block_column_size)); // should be the same
			// get which block column contains column i (optimize this away, probably need to use it when resuming)

			//_ASSERTE(n_block_column_size <= n_diag_band_width); // make this into a run-time check in the production code
			// make sure it is not longer than the diagonal (otherwise we will not have enough backlog to calculate all the off-diagonal elements)

			r_R.UpperTriangular_Solve(p_column, n, n_block_col); // backsub, only the nonzero part of the column (started at (block) column which contains column i, with no loss of generality)
			// this seems to be O(i) divisions + O(nnz) MADs in the given (block) column range
			// that sums up to O(n^2/2) divisions + O(nnz log(nnz))?? MADs ... some quadratute anyways

#ifdef _DEBUG
			std::vector<double> backsub_test(n, 0);
			backsub_test[i] = 1; // !!
			r_R.UpperTriangular_Solve(&backsub_test[0], n); // full backsub
			_ASSERTE(!memcmp(p_column, &backsub_test[0], n * sizeof(double)));
			// make sure that the result is correct
#endif // _DEBUG

			_ASSERTE((Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_column + i + 1, n - i - 1).norm() == 0)); // double pars required because of the comma in Map params
			// everything below i is zero (true)

			for(size_t k = n - n_column_num; k <= i; ++ k) {
				for(size_t j = 0; j <= i; ++ j)
					r_marginals(j, k) += p_column[j] * p_column[k]; // it is symmetric, indexing arbitrary
			}
			// accumulate the entries of the covariace matrix. this is O(n^3/2) MADs for the full matrix
			// note that to calculate even only the diagonal, we need full columns
		}
	}

	static void Calculate_DenseMarginals_Recurrent(Eigen::MatrixXd &r_marginals, const CUberBlockMatrix &r_R)
	{
		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.resize(n, n);

		CUberBlockMatrix R_tr;
		R_tr.TransposeOf(r_R);
		// need transpose of R

		for(size_t i = 0; i < n; ++ i) {
			if(!i) { // the last row / col is easy, as it has no refs to the next ones
				Eigen::MatrixXd::ColXpr last_col = r_marginals.col(n - 1);
				size_t lb = r_R.n_BlockColumn_Num() - 1;
				CUberBlockMatrix::_TyMatrixXdRef last_block =
					r_R.t_Block_AtColumn(lb, r_R.n_BlockColumn_Block_Num(lb) - 1);
				last_col.setZero(); // !!
				last_col(n - 1) = 1 / last_block(last_block.rows() - 1, last_block.cols() - 1);
				r_R.UpperTriangular_Solve(&last_col(0), n); // all columns needed
				// calculates the whole last column of C

				r_marginals.row(n - 1).head(n - 1) = r_marginals.col(n - 1).head(n - 1).transpose();
				// copy that also to the last row, to form the full matrix and not just upper-triangular
			} else { // columns with references to the subsequent columns
				i = n - 1 - i;
				// fill the matrix from the back

				size_t n_block_column_size;
				size_t n_block_column = r_R.n_Find_BlockColumn(i, n_block_column_size);
				size_t n_block_column_base = r_R.n_BlockColumn_Base(n_block_column);
				// gets the corresponding block col (can use decrementing strategy like in C_direct)

				CUberBlockMatrix::_TyMatrixXdRef cur_diag_block =
					r_R.t_Block_AtColumn(n_block_column, r_R.n_BlockColumn_Block_Num(n_block_column) - 1);
				double f_R_ii = cur_diag_block(i - n_block_column_base, i - n_block_column_base);
				double f_R_ii_inv = 1 / f_R_ii;
				// get the diagonal element

				double f_diag_sum = 0;
				for(size_t j = i + 1; j < n; ++ j) {
					size_t n_block_row_size;
					size_t n_block_row = r_R.n_Find_BlockColumn(j, n_block_row_size); // R has symmetric layout
					size_t n_block_row_base = r_R.n_BlockColumn_Base(n_block_row);
					// look up the row in R (= column in R_tr)

					CUberBlockMatrix::_TyMatrixXdRef block_i_j = R_tr.t_GetBlock_Log(n_block_row, n_block_column);
					if(!block_i_j.cols())
						continue; // no such block (zero, sparse) // todo - rewrite to form of a loop over *existing* blocks instead of looking up all blocks in dense manner
					double R_i_j = block_i_j(j - n_block_row_base, i - n_block_column_base);
					f_diag_sum += r_marginals(i, j) * R_i_j; // this is bad, accesses the matrix by rows (need transpose)
				}
				r_marginals(i, i) = f_R_ii_inv * (f_R_ii_inv - f_diag_sum);
				// calculate the diagonal element

				for(size_t j = i; j > 0;) {
					-- j; // j is i in the book
					// i is k in the book

					size_t n_block_j_size;
					size_t n_block_j = r_R.n_Find_BlockColumn(j, n_block_j_size);
					size_t n_block_j_base = r_R.n_BlockColumn_Base(n_block_j);

					double f_sum_Rjk_Cik = 0/*, f_sum_part0, f_sum_first_elem*/;
					for(size_t k = j + 1; k < n; ++ k) {
						//if(k == i + 1)
						//	f_sum_part0 = f_sum_Rjk_Cik; // note that the second half of the sum might be actually recurrent and easy to recover from the previous diagonal elements
						// less code this way

						size_t n_block_k_size;
						size_t n_block_k = r_R.n_Find_BlockColumn(k, n_block_k_size);
						size_t n_block_k_base = r_R.n_BlockColumn_Base(n_block_k);

						CUberBlockMatrix::_TyMatrixXdRef block_j_k = R_tr.t_GetBlock_Log(n_block_k, n_block_j);
						if(!block_j_k.cols())
							continue; // no such block (zero, sparse) // todo - rewrite to form of a loop over *existing* blocks 

						f_sum_Rjk_Cik += r_marginals(i, k) * block_j_k(k - n_block_k_base, j - n_block_j_base);
						// product

						//if(k == i)
						//	f_sum_first_elem = C_dep(i, k) * block_j_k(k - n_block_k_base, j - n_block_j_base); // save that as well
					}
					// note that this is a single loop, which skips iterating over element i

					/*{ // debugging of more recurrent formula
						size_t m = i + 1; // wha?
						size_t n_block_m_size;
						size_t n_block_m = R.n_Find_BlockColumn(m, n_block_m_size);
						size_t n_block_m_base = R.n_BlockColumn_Base(n_block_m);
						CUberBlockMatrix::_TyMatrixXdRef diag_block_m =
							R.t_BlockAt(R.n_BlockColumn_Block_Num(n_block_m) - 1, n_block_m);
						double f_R_mm = cur_diag_block(m - n_block_m_base, m - n_block_m_base);
						double f_sum_part1 = -(C_dep(m, m) * f_R_mm - 1 / f_R_mm);
						double f_sum_parted = f_sum_part0 + f_sum_part1;
						double f_sum_err = fabs(f_sum_parted - f_sum_Rjk_Cik);
					}*/

					CUberBlockMatrix::_TyMatrixXdRef block_j_j =
						r_R.t_Block_AtColumn(n_block_j, r_R.n_BlockColumn_Block_Num(n_block_j) - 1);
					double R_j_j = block_j_j(j - n_block_j_base, j - n_block_j_base);
					// get another diagonal element of R

					r_marginals(i, j) = r_marginals(j, i) = f_sum_Rjk_Cik / -R_j_j;
				}
				// calculate C_dep_{n - 1 - i} by recurrent formula

				//printf("%d, error: %g%6s\r", i, (C_dep.col(i) - C.col(i)).norm(), "");

				i = n - 1 - i;
				// go back to i
			}
		}
	}

	static void Calculate_DenseMarginals_LastNCols_Recurrent(Eigen::MatrixXd &r_marginals,
		const CUberBlockMatrix &r_R, size_t n_column_num)
	{
		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.resize(n, n);

		CUberBlockMatrix R_tr;
		R_tr.TransposeOf(r_R);
		// need transpose of R

		for(size_t i = 0; i < std::min(n_column_num, n); ++ i) { // counts in reverse
			if(!i) { // the last row / col is easy, as it has no refs to the next ones
				Eigen::MatrixXd::ColXpr last_col = r_marginals.col(n - 1);
				size_t lb = r_R.n_BlockColumn_Num() - 1;
				CUberBlockMatrix::_TyMatrixXdRef last_block =
					r_R.t_Block_AtColumn(lb, r_R.n_BlockColumn_Block_Num(lb) - 1);
				last_col.setZero(); // !!
				last_col(n - 1) = 1 / last_block(last_block.rows() - 1, last_block.cols() - 1);
				r_R.UpperTriangular_Solve(&last_col(0), n); // all columns needed
				// calculates the whole last column of C

				r_marginals.row(n - 1).head(n - 1) = r_marginals.col(n - 1).head(n - 1).transpose();
				// copy that also to the last row, to form the full matrix and not just upper-triangular
			} else { // columns with references to the subsequent columns
				i = n - 1 - i;
				// fill the matrix from the back

				size_t n_block_column_size;
				size_t n_block_column = r_R.n_Find_BlockColumn(i, n_block_column_size);
				size_t n_block_column_base = r_R.n_BlockColumn_Base(n_block_column);
				// gets the corresponding block col (can use decrementing strategy like in C_direct)

				CUberBlockMatrix::_TyMatrixXdRef cur_diag_block =
					r_R.t_Block_AtColumn(n_block_column, r_R.n_BlockColumn_Block_Num(n_block_column) - 1);
				double f_R_ii = cur_diag_block(i - n_block_column_base, i - n_block_column_base);
				double f_R_ii_inv = 1 / f_R_ii;
				// get the diagonal element

				double f_diag_sum = 0;
				for(size_t j = i + 1; j < n; ++ j) {
					size_t n_block_row_size;
					size_t n_block_row = r_R.n_Find_BlockColumn(j, n_block_row_size); // R has symmetric layout
					size_t n_block_row_base = r_R.n_BlockColumn_Base(n_block_row);
					// look up the row in R (= column in R_tr)

					CUberBlockMatrix::_TyMatrixXdRef block_i_j = R_tr.t_GetBlock_Log(n_block_row, n_block_column);
					if(!block_i_j.cols())
						continue; // no such block (zero, sparse) // todo - rewrite to form of a loop over *existing* blocks instead of looking up all blocks in dense manner
					double R_i_j = block_i_j(j - n_block_row_base, i - n_block_column_base);
					f_diag_sum += r_marginals(i, j) * R_i_j; // this is bad, accesses the matrix by rows (need transpose)
				}
				r_marginals(i, i) = f_R_ii_inv * (f_R_ii_inv - f_diag_sum);
				// calculate the diagonal element

				for(size_t j = i; j > 0;) {
					-- j; // j is i in the book
					// i is k in the book

					size_t n_block_j_size;
					size_t n_block_j = r_R.n_Find_BlockColumn(j, n_block_j_size);
					size_t n_block_j_base = r_R.n_BlockColumn_Base(n_block_j);

					double f_sum_Rjk_Cik = 0/*, f_sum_part0, f_sum_first_elem*/;
					for(size_t k = j + 1; k < n; ++ k) {
						//if(k == i + 1)
						//	f_sum_part0 = f_sum_Rjk_Cik; // note that the second half of the sum might be actually recurrent and easy to recover from the previous diagonal elements
						// less code this way

						size_t n_block_k_size;
						size_t n_block_k = r_R.n_Find_BlockColumn(k, n_block_k_size);
						size_t n_block_k_base = r_R.n_BlockColumn_Base(n_block_k);

						CUberBlockMatrix::_TyMatrixXdRef block_j_k = R_tr.t_GetBlock_Log(n_block_k, n_block_j);
						if(!block_j_k.cols())
							continue; // no such block (zero, sparse) // todo - rewrite to form of a loop over *existing* blocks 

						f_sum_Rjk_Cik += r_marginals(i, k) * block_j_k(k - n_block_k_base, j - n_block_j_base);
						// product

						//if(k == i)
						//	f_sum_first_elem = C_dep(i, k) * block_j_k(k - n_block_k_base, j - n_block_j_base); // save that as well
					}
					// note that this is a single loop, which skips iterating over element i

					/*{ // debugging of more recurrent formula
						size_t m = i + 1; // wha?
						size_t n_block_m_size;
						size_t n_block_m = R.n_Find_BlockColumn(m, n_block_m_size);
						size_t n_block_m_base = R.n_BlockColumn_Base(n_block_m);
						CUberBlockMatrix::_TyMatrixXdRef diag_block_m =
							R.t_BlockAt(R.n_BlockColumn_Block_Num(n_block_m) - 1, n_block_m);
						double f_R_mm = cur_diag_block(m - n_block_m_base, m - n_block_m_base);
						double f_sum_part1 = -(C_dep(m, m) * f_R_mm - 1 / f_R_mm);
						double f_sum_parted = f_sum_part0 + f_sum_part1;
						double f_sum_err = fabs(f_sum_parted - f_sum_Rjk_Cik);
					}*/

					CUberBlockMatrix::_TyMatrixXdRef block_j_j =
						r_R.t_Block_AtColumn(n_block_j, r_R.n_BlockColumn_Block_Num(n_block_j) - 1);
					double R_j_j = block_j_j(j - n_block_j_base, j - n_block_j_base);
					// get another diagonal element of R

					r_marginals(i, j) = r_marginals(j, i) = f_sum_Rjk_Cik / -R_j_j;
				}
				// calculate C_dep_{n - 1 - i} by recurrent formula

				//printf("%d, error: %g%6s\r", i, (C_dep.col(i) - C.col(i)).norm(), "");

				i = n - 1 - i;
				// go back to i
			}
		}
	}

	static void Calculate_DiagonalMarginals(Eigen::MatrixXd &r_marginals, const CUberBlockMatrix &r_R,
		size_t n_diag_band_width = 3)
	{
		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.resize(n, n);

		Eigen::MatrixXd prev_column_buffer(n, n_diag_band_width); // todo - cache align, will likely access the same elements in different vectors (needs to be implemented using a strided map)
		// previous columns

		for(size_t i = 0, n_block_col = -1, n_col_remains = 1; i < n; ++ i) {
			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = r_R.n_BlockColumn_Column_Num(++ n_block_col);
			// if done before, it avoids both referencing block col 0 in an empty matrix
			// and determining no. of columns in one past the last block column at the end

			Eigen::MatrixXd::ColXpr cur_row = prev_column_buffer.col(i % n_diag_band_width);
			cur_row.setZero();
			cur_row(i) = 1; // !!
			r_R.UpperTriangularTranspose_Solve(&cur_row(0), n, n_block_col); // forward substitution; can also skip to the current column
			// calculate a row of the R_inv? maybe? seems outright wrong, but actually gives a rather correct answer.

			_ASSERTE(cur_row.head(i).norm() == 0);
			// everything above i is zero (true)

			size_t n_block_org = r_R.n_BlockColumn_Base(n_block_col);
			r_marginals(i, i) = cur_row.tail(n - i).squaredNorm();
			//for(size_t j = std::max(n_diag_band_width, i) - n_diag_band_width; j < i; ++ j) { // could do two loops; one for i < n_diag_band_width and the one for the rest (low prio, almost doesn't matter)
			for(size_t j = i - n_block_org/*i % n_diag_band_width*/; j < i; ++ j) { // thin block diagonal, works even for mixed-size vertices
				Eigen::MatrixXd::ColXpr prev_row = prev_column_buffer.col(j % n_diag_band_width);
				r_marginals(j, i) = r_marginals(i, j) = prev_row.tail(n - j).dot(cur_row.tail(n - j));
			}
			// calculate banded diagonal; this works well, the banded diagonal is identical to the one in C or C_full
			// todo - store in a smaller matrix (will only contain block cross-covs of the vertices; that involves further ~50% reduction in computation)
		}
	}

	static void Calculate_DiagonalMarginals_Parallel(Eigen::MatrixXd &r_marginals,
		const CUberBlockMatrix &r_R, size_t n_diag_band_width = 3)
	{
		const size_t n = r_R.n_Column_Num(); // in elements
		r_marginals.resize(n, n);

#ifdef _OPENMP
		#pragma omp parallel
		{
			int n_tid = omp_get_thread_num();
			int n_thread_num = omp_get_num_threads();
			size_t n_start = n_tid * (n / n_thread_num);
			size_t n_end = (n_tid + 1 < n_thread_num)? n_start + n / n_thread_num : n;
			// split to bands to be processed in parallel

			Calculate_DiagonalMarginals(r_marginals, r_R, n_start, n_end, n_diag_band_width);
			// process in parallel
		}
#else // _OPENMP
		Calculate_DiagonalMarginals(r_marginals, r_R, n_diag_band_width);
#endif // _OPENMP
	}

	static void Calculate_DiagonalMarginals(Eigen::MatrixXd &r_marginals, const CUberBlockMatrix &r_R,
		size_t n_start_column, size_t n_end_column, size_t n_diag_band_width = 3)
	{
		_ASSERTE(n_start_column <= n_end_column); // should be a valid range

		const size_t n = r_R.n_Column_Num(); // in elements

		r_marginals.conservativeResize(n, n); // might work in a group

		Eigen::MatrixXd prev_column_buffer(n, n_diag_band_width); // todo - cache align, will likely access the same elements in different vectors (needs to be implemented using a strided map)
		// previous columns

		for(size_t i = 0, n_block_col = -1, n_col_remains = 1; i < std::min(n, n_end_column); ++ i) {
			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = r_R.n_BlockColumn_Column_Num(++ n_block_col);
			// if done before, it avoids both referencing block col 0 in an empty matrix
			// and determining no. of columns in one past the last block column at the end

			if(i + n_diag_band_width < n_start_column)
				continue;
			// todo - do a proper binary search for the column

			Eigen::MatrixXd::ColXpr cur_row = prev_column_buffer.col(i % n_diag_band_width);
			cur_row.setZero();
			cur_row(i) = 1; // !!
			r_R.UpperTriangularTranspose_Solve(&cur_row(0), n, n_block_col); // forward substitution; can also skip to the current column
			// calculate a row of the R_inv? maybe? seems outright wrong, but actually gives a rather correct answer.

			_ASSERTE(cur_row.head(i).norm() == 0);
			// everything above i is zero (true)

			size_t n_block_org = r_R.n_BlockColumn_Base(n_block_col);
			r_marginals(i, i) = cur_row.tail(n - i).squaredNorm();
			//for(size_t j = std::max(n_diag_band_width, i) - n_diag_band_width; j < i; ++ j) { // could do two loops; one for i < n_diag_band_width and the one for the rest (low prio, almost doesn't matter)
			for(size_t j = i - n_block_org/*i % n_diag_band_width*/; j < i; ++ j) { // thin block diagonal, works even for mixed-size vertices
				Eigen::MatrixXd::ColXpr prev_row = prev_column_buffer.col(j % n_diag_band_width);
				r_marginals(j, i) = r_marginals(i, j) = prev_row.tail(n - j).dot(cur_row.tail(n - j));
			}
			// calculate banded diagonal; this works well, the banded diagonal is identical to the one in C or C_full
			// todo - store in a smaller matrix (will only contain block cross-covs of the vertices; that involves further ~50% reduction in computation)
		}
	}

	static void Marginals_Test(const CUberBlockMatrix &r_R, size_t n_diag_band_width = 3)
	{
		if(r_R.n_Column_Num() < n_diag_band_width) {
			fprintf(stderr, "error: matrix too small for banded tests (%d x %d): skipping\n",
				int(r_R.n_Row_Num()), int(r_R.n_Column_Num()));
			return;
		}
		// might not be handled correctly

		Eigen::MatrixXd C_ref, C_slow, C_fast, C_rec, C_diag,
			C_diag_para, C_rband_slow, C_rband_fast, C_rband_rec;

		CTimer t;
		double f_time_ref = 0, f_time_slow = 0, f_time_rec = 0,
			f_time_diag = 0, f_time_diag_para = 0, f_time_rbslow = 0,
			f_time_rbfast = 0, f_time_rbrec = 0, f_time_fast = 0;
		size_t n_ref_pass_num = 0, n_slow_pass_num = 0, n_rec_pass_num = 0,
			n_diag_pass_num = 0, n_diag_para_pass_num = 0, n_rbslow_pass_num = 0,
			n_rbfast_pass_num = 0, n_rbrec_pass_num = 0, n_fast_pass_num = 0;

		/*for(;;) {
			double f_start = t.f_Time();
			Calculate_DenseMarginals_Ref(C_ref, r_R);
			++ n_ref_pass_num;
			f_time_ref += t.f_Time() - f_start;
			if((f_time_ref >= 1 && n_ref_pass_num >= 10) || f_time_ref > 4)
				break;
		}
		f_time_ref /= n_ref_pass_num;

		for(;;) {
			double f_start = t.f_Time();
			Calculate_DenseMarginals_Slow(C_slow, r_R);
			++ n_slow_pass_num;
			f_time_slow += t.f_Time() - f_start;
			if((f_time_slow >= 1 && n_slow_pass_num >= 10) || f_time_slow > 4)
				break;
		}
		f_time_slow /= n_slow_pass_num;*/

		for(;;) {
			double f_start = t.f_Time();
			Calculate_DenseMarginals_Fast(C_fast, r_R);
			++ n_fast_pass_num;
			f_time_fast += t.f_Time() - f_start;
			if((f_time_fast >= 1 && n_fast_pass_num >= 10) || f_time_fast > 4)
				break;
		}
		f_time_fast /= n_fast_pass_num;

		/*for(;;) {
			double f_start = t.f_Time();
			Calculate_DenseMarginals_Recurrent(C_rec, r_R);
			++ n_rec_pass_num;
			f_time_rec += t.f_Time() - f_start;
			if((f_time_rec >= 1 && n_rec_pass_num >= 10) || f_time_rec > 4)
				break;
		}
		f_time_rec /= n_rec_pass_num;

		for(;;) {
			double f_start = t.f_Time();
			Calculate_DiagonalMarginals(C_diag, r_R, n_diag_band_width);
			++ n_diag_pass_num;
			f_time_diag += t.f_Time() - f_start;
			if((f_time_diag >= 1 && n_diag_pass_num >= 10) || f_time_diag > 4)
				break;
		}
		f_time_diag /= n_diag_pass_num;

		for(;;) {
			double f_start = t.f_Time();
			Calculate_DiagonalMarginals_Parallel(C_diag_para, r_R, n_diag_band_width);
			++ n_diag_para_pass_num;
			f_time_diag_para += t.f_Time() - f_start;
			if((f_time_diag_para >= 1 && n_diag_para_pass_num >= 10) || f_time_diag_para > 4)
				break;
		}
		f_time_diag_para /= n_diag_para_pass_num;

		for(;;) {
			double f_start = t.f_Time();
			Calculate_DenseMarginals_LastNColumns_Slow(C_rband_slow, r_R, n_diag_band_width);
			++ n_rbslow_pass_num;
			f_time_rbslow += t.f_Time() - f_start;
			if((f_time_rbslow >= 1 && n_rbslow_pass_num >= 10) || f_time_rbslow > 4)
				break;
		}
		f_time_rbslow /= n_rbslow_pass_num;*/

		for(;;) {
			double f_start = t.f_Time();
			Calculate_DenseMarginals_LastNColumns_Fast(C_rband_fast, r_R, n_diag_band_width);
			++ n_rbfast_pass_num;
			f_time_rbfast += t.f_Time() - f_start;
			if((f_time_rbfast >= 1 && n_rbfast_pass_num >= 10) || f_time_rbfast > 4)
				break;
		}
		f_time_rbfast /= n_rbfast_pass_num;

		/*for(;;) {
			double f_start = t.f_Time();
			Calculate_DenseMarginals_LastNCols_Recurrent(C_rband_rec, r_R, n_diag_band_width);
			++ n_rbrec_pass_num;
			f_time_rbrec += t.f_Time() - f_start;
			if((f_time_rbrec >= 1 && n_rbrec_pass_num >= 10) || f_time_rbrec > 4)
				break;
		}
		f_time_rbrec /= n_rbrec_pass_num;*/

		printf("%6s took %.3f msec\n", "ref", f_time_ref * 1000);
		printf("%6s took %.3f msec\n", "slow", f_time_slow * 1000);
		printf("%6s took %.3f msec\n", "fast", f_time_fast * 1000);
		printf("%6s took %.3f msec\n", "rec", f_time_rec * 1000);
		printf("%6s took %.3f msec\n", "diag", f_time_diag * 1000);
		printf("%6s took %.3f msec\n", "diag-p", f_time_diag_para * 1000);
		printf("%6s took %.3f msec\n", "rbslow", f_time_rbslow * 1000);
		printf("%6s took %.3f msec\n", "rbfast", f_time_rbfast * 1000);
		printf("%6s took %.3f msec\n", "rbrec", f_time_rbrec * 1000);
		// print times

		/*printf("norm of C_slow - C_ref is %g\n", (C_slow - C_ref).norm());
		printf("norm of C_fast - C_ref is %g\n", (C_fast - C_ref).norm());
		printf("norm of C_rec - C_ref is %g\n", (C_rec - C_ref).norm());
		printf("norm of the diagonal of C_diag - C_ref is %g\n", (C_diag.diagonal() - C_ref.diagonal()).norm());
		printf("norm of the diagonal of C_diag_para - C_ref is %g\n", (C_diag_para.diagonal() - C_ref.diagonal()).norm());
		printf("norm of the last N columns of C_rband_slow - C_ref is %g\n",
			(C_rband_slow.rightCols(n_diag_band_width) - C_ref.rightCols(n_diag_band_width)).norm());
		printf("norm of the last N columns of C_rband_fast - C_ref is %g\n",
			(C_rband_fast.rightCols(n_diag_band_width) - C_ref.rightCols(n_diag_band_width)).norm());
		printf("norm of the last N columns of C_rband_rec - C_ref is %g\n",
			(C_rband_rec.rightCols(n_diag_band_width) - C_ref.rightCols(n_diag_band_width)).norm());*/
		// print precision
	}

#if 0
	void Marginals_DevelCode()
	{
		const CUberBlockMatrix &lambda = solver.r_Lambda();
		// get system matrix

		cs *p_lam = lambda.p_Convert_to_Sparse();
		FILE *p_fw = fopen("lambda.txt", "w");
		CDebug::Print_SparseMatrix_in_MatlabFormat(p_fw, p_lam, "lambda = ", "' % this is supposed to be upper-tri\n");
		fclose(p_fw);
		cs_spfree(p_lam);
		// dump the optimized lambda

		CUberBlockMatrix R;
		R.CholeskyOf(lambda);
		// take Cholesky

		const size_t n = R.n_Column_Num(); // in elements

		Eigen::MatrixXd R_inv(n, n); // R_inv = S
		for(size_t i = 0, n_block_col = -1, n_col_remains = 1; i < n; ++ i) {
			double *p_column = &R_inv.col(i)(0);
			memset(p_column, 0, n * sizeof(double));
			p_column[i] = 1;
			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = R.n_BlockColumn_Column_Num(++ n_block_col);
			R.UpperTriangular_Solve(p_column, n, n_block_col); // backsub, only the nonzero part of the column (started at (block) column which contains column i, with no loss of generality)
		}
		Eigen::MatrixXd C = R_inv.lazyProduct(R_inv.transpose()); // C = SS^T
		// calculate the covariance (assume that this is correct)

		Eigen::MatrixXd C_full(n, n), C_direct(n, n), C_dep(n, n);
		C_full.setZero();
		C_direct.setZero();
		C_dep.setZero();
		// only the diagonal of C

		size_t n_diag_band_width = 3; // = max column width occuring in the matrix (max vertex dim)
		// size of the diagonal (also the size of the buffer for the past rows of inv(R))

		Eigen::MatrixXd prev_column_buffer(n, n_diag_band_width); // todo - cache align, will likely access the same elements in different vectors
		// previous columns

		CUberBlockMatrix R_tr;
		R.TransposeTo(R_tr);
		// need transpose of R

		//Eigen::MatrixXd R_inv(n, n); // R_inv = S
		for(size_t i = 0, n_block_col = -1, n_col_remains = 1; i < n; ++ i) {
			double *p_column = &R_inv.col(i)(0);
			// get dense column data from the Eigen matrix (should be packed)

			memset(p_column, 0, n * sizeof(double));
			p_column[i] = 1;
			// make a column vector with a single 1 in it

			if(!(-- n_col_remains)) // triggers in the first iteration, loads up column width
				n_col_remains = R.n_BlockColumn_Column_Num(++ n_block_col);
			// if done before, it avoids both referencing block col 0 in an empty matrix
			// and determining no. of columns in one past the last block column at the end

			size_t UNUSED(n_block_column_size);
			size_t n_block_column = R.n_Find_BlockColumn(i, n_block_column_size);
			_ASSERTE(n_block_col == n_block_column); // should be the same
			// get which block column contains column i (optimize this away, probably need to use it when resuming)

			_ASSERTE(n_block_column_size <= n_diag_band_width); // make this into a run-time check in the production code
			// make sure it is not longer than the diagonal (otherwise we will not have enough backlog to calculate all the off-diagonal elements)

			R.UpperTriangular_Solve(p_column, n, n_block_col); // backsub, only the nonzero part of the column (started at (block) column which contains column i, with no loss of generality)
			// this seems to be O(i) divisions + O(nnz) MADs in the given (block) column range
			// that sums up to O(n^2/2) divisions + O(nnz log(nnz))?? MADs ... some quadratute anyways

			std::vector<double> backsub_test(n, 0);
			backsub_test[i] = 1; // !!
			R.UpperTriangular_Solve(&backsub_test[0], n); // full backsub
			_ASSERTE(!memcmp(p_column, &backsub_test[0], n * sizeof(double)));
			// make sure that the result is correct

			_ASSERTE((Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_column + i + 1, n - i - 1).norm() == 0)); // double pars required because of the comma in Map params
			// everything below i is zero (true)

			for(size_t k = 0; k <= i; ++ k) {
				for(size_t j = 0; j <= i; ++ j)
					C_full(j, k) += p_column[j] * p_column[k]; // it is symmetric, indexing arbitrary
			}
			// accumulate the entries of the covariace matrix. this is O(n^3/2) MADs for the full matrix
			// note that to calculate even only the diagonal, we need full columns

			Eigen::MatrixXd::ColXpr cur_row = prev_column_buffer.col(i % n_diag_band_width);
			cur_row.setZero();
			cur_row(i) = 1; // !!
			R.UpperTriangularTranspose_Solve(&cur_row(0), n, n_block_col); // forward substitution; can also skip to the current column
			// calculate a row of the R_inv? maybe? seems outright wrong, but actually gives a rather correct answer.

			_ASSERTE(cur_row.head(i).norm() == 0);
			// everything above i is zero (true)

			size_t n_block_org = R.n_BlockColumn_Base(n_block_col);
			C_direct(i, i) = cur_row.tail(n - i).squaredNorm();
			//for(size_t j = std::max(n_diag_band_width, i) - n_diag_band_width; j < i; ++ j) { // could do two loops; one for i < n_diag_band_width and the one for the rest (low prio, almost doesn't matter)
			for(size_t j = i - n_block_org/*i % n_diag_band_width*/; j < i; ++ j) { // thin block diagonal, works even for mixed-size vertices
				Eigen::MatrixXd::ColXpr prev_row = prev_column_buffer.col(j % n_diag_band_width);
				C_direct(j, i) = C_direct(i, j) = prev_row.tail(n - j).dot(cur_row.tail(n - j));
			}
			// calculate banded diagonal; this works well, the banded diagonal is identical to the one in C or C_full
			// todo - store in a smaller matrix (will only contain block cross-covs of the vertices; that involves further ~50% reduction in computation)

			//C_direct(i, i) = 0;
			//for(size_t j = i; j < n; ++ j)
			//	C_direct(i, i) += p_column[j] * p_column[j]; // O(n^2/2) for the diagonal of the full matrix
			// the diagonal is just the stuff squared. to calculate off-diagonals, we need to cache previous columns as well
			// also the columns are independent and there is no communication involved; this is easily parallelised, even on GPU

			if(!i) { // the last row / col is easy, as it has no refs to the next ones
				Eigen::MatrixXd::ColXpr last_col = C_dep.col(n - 1);
				size_t lb = R.n_BlockColumn_Num() - 1;
				CUberBlockMatrix::_TyMatrixXdRef last_block =
					R.t_BlockAt(R.n_BlockColumn_Block_Num(lb) - 1, lb);
				last_col(n - 1) = 1 / last_block(last_block.rows() - 1, last_block.cols() - 1);
				R.UpperTriangular_Solve(&last_col(0), n); // all columns needed
				// calculates the whole last column of C

				C_dep.row(n - 1).head(n - 1) = C_dep.col(n - 1).head(n - 1).transpose();
				// copy that also to the last row, to form the full matrix and not just upper-triangular
			} else { // columns with references to the subsequent columns
				i = n - 1 - i;
				// fill the matrix from the back

				size_t n_block_column_size;
				size_t n_block_column = R.n_Find_BlockColumn(i, n_block_column_size);
				size_t n_block_column_base = R.n_BlockColumn_Base(n_block_column);
				// gets the corresponding block col (can use decrementing strategy like in C_direct)

				CUberBlockMatrix::_TyMatrixXdRef cur_diag_block =
					R.t_BlockAt(R.n_BlockColumn_Block_Num(n_block_column) - 1, n_block_column);
				double f_R_ii = cur_diag_block(i - n_block_column_base, i - n_block_column_base);
				double f_R_ii_inv = 1 / f_R_ii;
				// get the diagonal element

				double f_diag_sum = 0;
				for(size_t j = i + 1; j < n; ++ j) {
					size_t n_block_row_size;
					size_t n_block_row = R.n_Find_BlockColumn(j, n_block_row_size); // R has symmetric layout
					size_t n_block_row_base = R.n_BlockColumn_Base(n_block_row);
					// look up the row in R (= column in R_tr)

					CUberBlockMatrix::_TyMatrixXdRef block_i_j = R_tr.t_GetBlock_Log(n_block_row, n_block_column);
					if(!block_i_j.cols())
						continue; // no such block (zero, sparse) // todo - rewrite to form of a loop over *existing* blocks instead of looking up all blocks in dense manner
					double R_i_j = block_i_j(j - n_block_row_base, i - n_block_column_base);
					f_diag_sum += C_dep(i, j) * R_i_j; // this is bad, accesses the matrix by rows (need transpose)
				}
				C_dep(i, i) = f_R_ii_inv * (f_R_ii_inv - f_diag_sum);
				for(size_t j = i; j > 0;) {
					-- j; // j is i in the book
					// i is k in the book

					size_t n_block_j_size;
					size_t n_block_j = R.n_Find_BlockColumn(j, n_block_j_size);
					size_t n_block_j_base = R.n_BlockColumn_Base(n_block_j);

					double f_sum_Rjk_Cik = 0/*, f_sum_part0, f_sum_first_elem*/;
					for(size_t k = j + 1; k < n; ++ k) {
						//if(k == i + 1)
						//	f_sum_part0 = f_sum_Rjk_Cik; // note that the second half of the sum might be actually recurrent and easy to recover from the previous diagonal elements
						// less code this way

						size_t n_block_k_size;
						size_t n_block_k = R.n_Find_BlockColumn(k, n_block_k_size);
						size_t n_block_k_base = R.n_BlockColumn_Base(n_block_k);

						CUberBlockMatrix::_TyMatrixXdRef block_j_k = R_tr.t_GetBlock_Log(n_block_k, n_block_j);
						if(!block_j_k.cols())
							continue; // no such block (zero, sparse) // todo - rewrite to form of a loop over *existing* blocks 

						f_sum_Rjk_Cik += C_dep(i, k) * block_j_k(k - n_block_k_base, j - n_block_j_base);
						// product

						//if(k == i)
						//	f_sum_first_elem = C_dep(i, k) * block_j_k(k - n_block_k_base, j - n_block_j_base); // save that as well
					}
					// note that this is a single loop, which skips iterating over element i

					/*{ // debugging of more recurrent formula
						size_t m = i + 1; // wha?
						size_t n_block_m_size;
						size_t n_block_m = R.n_Find_BlockColumn(m, n_block_m_size);
						size_t n_block_m_base = R.n_BlockColumn_Base(n_block_m);
						CUberBlockMatrix::_TyMatrixXdRef diag_block_m =
							R.t_BlockAt(R.n_BlockColumn_Block_Num(n_block_m) - 1, n_block_m);
						double f_R_mm = cur_diag_block(m - n_block_m_base, m - n_block_m_base);
						double f_sum_part1 = -(C_dep(m, m) * f_R_mm - 1 / f_R_mm);
						double f_sum_parted = f_sum_part0 + f_sum_part1;
						double f_sum_err = fabs(f_sum_parted - f_sum_Rjk_Cik);
					}*/

					CUberBlockMatrix::_TyMatrixXdRef block_j_j =
						R.t_BlockAt(R.n_BlockColumn_Block_Num(n_block_j) - 1, n_block_j);
					double R_j_j = block_j_j(j - n_block_j_base, j - n_block_j_base);
					// get another diagonal element of R

					C_dep(i, j) = C_dep(j, i) = f_sum_Rjk_Cik / -R_j_j;
				}
				// calculate C_dep_{n - 1 - i} by recurrent formula

				printf("%d, error: %g%6s\r", i, (C_dep.col(i) - C.col(i)).norm(), "");

				i = n - 1 - i;
				// go back to i
			}
		}
		// calculate inverse of R using backsubstitution

		printf("norm of C - C_full is %g\n", (C_full - C).norm());
		printf("norm of the diagonal of C - C_direct is %g\n", (C_direct.diagonal() - C.diagonal()).norm());
		printf("norm of C - C_dep is %g\n", (C_dep - C).norm());
		printf("norm of the last column of C - C_dep is %g\n", (C_dep.col(n - 1) - C.col(n - 1)).norm());
		printf("norm of the diagonal of C - C_dep is %g\n", (C_dep.diagonal() - C.diagonal()).norm());

		p_fw = fopen("R_inv.txt", "w");
		CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, R_inv, "R_inv = ", "\n");
		fclose(p_fw);
		p_fw = fopen("C.txt", "w");
		CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, C, "C = ", "\n");
		fclose(p_fw);
		p_fw = fopen("C_banded.txt", "w");
		CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, C_direct, "C_banded = ", "\n");
		fclose(p_fw);
		// save for comparison with matlab formulas
	}
#endif // 0
};

/**
 *	@brief marginal covariances cache object
 */
class CMarginalCovariance {
protected:
	Eigen::MatrixXd m_matrix; /**< @brief marginal covariance matrix */

public:
	/**
	 *	@brief gets marginal covariance matrix
	 *
	 *	@return Returns reference to the marginal covariance matrix.
	 *
	 *	@note Depending on the policy used, this matrix might not contain all of the
	 *		marginal covariance values, or they might not be up to date.
	 *		See e.g. Request_Block() for more details.
	 *	@note The space complexity required to store the full matrix can be significant,
	 *		more efficient methods for passing only a part of the matrix will be provided.
	 */
	inline Eigen::MatrixXd &r_Matrix()
	{
		return m_matrix;
	}

	/**
	 *	@brief gets marginal covariance matrix
	 *
	 *	@return Returns const reference to the marginal covariance matrix.
	 *
	 *	@note Depending on the policy used, this matrix might not contain all of the
	 *		marginal covariance values, or they might not be up to date.
	 *		See e.g. Request_Block() for more details.
	 *	@note The space complexity required to store the full matrix can be significant,
	 *		more efficient methods for passing only a part of the matrix will be provided.
	 */
	inline const Eigen::MatrixXd &r_Matrix() const
	{
		return m_matrix;
	}

	/**
	 *	@brief determines whether a specific block in the marginals matrix is up-to-date
	 *
	 *	@param[in] n_block_row is zero-based block row index (index of the first vertex)
	 *	@param[in] n_block_column is zero-based block column index (index of the second vertex)
	 *
	 *	@return Returns true if the queried area is up-to-date, otherwise returns false.
	 *
	 *	@note This function just determines if the queried area is up-to-date but does
	 *		not update it. Thus, it is very fast and can be used for e.g. covariance polling.
	 */
	inline bool b_Block_UpToDate(size_t UNUSED(n_block_row), size_t UNUSED(n_block_column)) const
	{
		return true;
	}

	/**
	 *	@brief determines whether a specific block column in the marginals matrix is up-to-date
	 *	@param[in] n_block_row is zero-based block row index (index of the first vertex)
	 *	@return Returns true if the queried area is up-to-date, otherwise returns false.
	 *	@note This function just determines if the queried area is up-to-date but does
	 *		not update it. Thus, it is very fast and can be used for e.g. covariance polling.
	 */
	inline bool b_Column_UpToDate(size_t UNUSED(n_block_column)) const
	{
		return true;
	}

	/**
	 *	@brief determines whether the block diagonal in the marginals matrix is up-to-date
	 *	@return Returns true if the queried area is up-to-date, otherwise returns false.
	 *	@note This function just determines if the queried area is up-to-date but does
	 *		not update it. Thus, it is very fast and can be used for e.g. covariance polling.
	 */
	inline bool b_Diagonal_UpToDate() const
	{
		return true;
	}

	/**
	 *	@brief determines whether all of the marginals matrix is up-to-date
	 *	@return Returns true if the queried area is up-to-date, otherwise returns false.
	 *	@note This function just determines if the queried area is up-to-date but does
	 *		not update it. Thus, it is very fast and can be used for e.g. covariance polling.
	 */
	inline bool b_FullMatrix_UpToDate() const
	{
		return true;
	}

	/**
	 *	@brief makes sure that a block, corresponding to covariance of two vertices
	 *		is calculated (will perform the calculation if it is not)
	 *
	 *	@param[in] n_block_row is zero-based block row index (index of the first vertex)
	 *	@param[in] n_block_column is zero-based block column index (index of the second vertex)
	 *
	 *	@note Depending on the cache miss policy, this will possibly also calculate
	 *		other parts of the matrix (e.g. the whole block column or the full matrix).
	 *	@note This function throws std::bad_alloc.
	 */
	void Request_Block(size_t UNUSED(n_block_row), size_t UNUSED(n_block_column)) // throw(std::bad_alloc)
	{
		// at the moment, all of the matrix is available, requests have no effect.
		// this will change with more efficient algorithms in the future ...
	}

	/**
	 *	@brief makes sure that a block column, corresponding to covariance of one vertex
	 *		with all the vertices (including itself) is calculated (will perform the
	 *		calculation if it is not)
	 *
	 *	@param[in] n_block_column is zero-based block column index
	 *
	 *	@note Depending on the cache miss policy, this will possibly also calculate
	 *		other parts of the matrix (e.g. the full matrix).
	 *	@note This function throws std::bad_alloc.
	 */
	void Request_BlockColumn(size_t UNUSED(n_block_column)) // throw(std::bad_alloc)
	{
		// at the moment, all of the matrix is available, requests have no effect.
		// this will change with more efficient algorithms in the future ...
	}

	/**
	 *	@brief makes sure that block diagonal portion of the covariance matrix is calculated
	 *		(will perform the calculation if it is not)
	 *
	 *	The block diagonal means blocks at the diagonal (corresponds to auto-covariances of each
	 *	vertex with itself).
	 *
	 *	@note Depending on the cache miss policy, this will possibly also calculate
	 *		other parts of the matrix (e.g. the full matrix).
	 *	@note This function throws std::bad_alloc.
	 */
	void Request_BlockDiagonal() // throw(std::bad_alloc)
	{
		// at the moment, all of the matrix is available, requests have no effect.
		// this will change with more efficient algorithms in the future ...
	}

	/**
	 *	@brief makes sure that everything in the marginal covariance matrix is up-to-date
	 *	@note This function throws std::bad_alloc.
	 */
	void Request_Full() // throw(std::bad_alloc)
	{
		// at the moment, all of the matrix is available, requests have no effect.
		// this will change with more efficient algorithms in the future ...
	}
};

#endif // __MARGINALS_INCLUDED
