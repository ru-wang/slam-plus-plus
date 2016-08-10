/*
								+---------------------------------+
								|                                 |
								|  ***   Über Block Matrix   ***  |
								|                                 |
								| Copyright  (c) -tHE SWINe- 2013 |
								|                                 |
								|       BlockMatrixFBS.inl        |
								|                                 |
								+---------------------------------+
*/

#pragma once
#ifndef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_FUNCTIONS_IMPLEMENTATION_INCLUDED
#define __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_FUNCTIONS_IMPLEMENTATION_INCLUDED

/**
 *	@file include/slam/BlockMatrixFBS.inl
 *	@date 2013
 *	@author -tHE SWINe-
 *	@brief the überblockmatrix fixed block size functions
 *	@note This file is not to be included; it is automatically included from BlockMatrix.h
 */

#include "slam/BlockMatrix.h"

template <class CBlockMatrixTypelist>
void CUberBlockMatrix::PostMultiply_Add_FBS(double *p_dest_vector, size_t UNUSED(n_dest_size),
	const double *p_src_vector, size_t UNUSED(n_src_size)) const
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	PostMultiply_Add(p_dest_vector, n_dest_size, p_src_vector, n_src_size);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

	_ASSERTE(p_dest_vector != p_src_vector);
	_ASSERTE(n_dest_size == m_n_col_num);
	_ASSERTE(n_src_size == m_n_row_num);

	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;
		blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PostMAD<CBlockMatrixTypelist>::PostMAD_OuterLoop(
			r_t_col, p_dest_vector, p_src_vector, m_block_rows_list);
		// execute templated middle loop selector
	}
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
void CUberBlockMatrix::PostMultiply_Add_FBS_Parallel(double *p_dest_vector, size_t UNUSED(n_dest_size),
	const double *p_src_vector, size_t UNUSED(n_src_size)) const
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	PostMultiply_Add_Parallel(p_dest_vector, n_dest_size, p_src_vector, n_src_size);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

	_ASSERTE(p_dest_vector != p_src_vector);
	_ASSERTE(n_dest_size == m_n_col_num);
	_ASSERTE(n_src_size == m_n_row_num);

#ifdef _OPENMP
	// note that this loop can be omp parallelized (just have to use signed integer indexing)
	_ASSERTE(m_block_cols_list.size() <= INT_MAX);
	int n = int(m_block_cols_list.size());
#pragma omp parallel for default(shared)
	for(int i = 0; i < n; ++ i) {
		const TColumn &r_t_col = m_block_cols_list[i];
#else // _OPENMP
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;
#endif // _OPENMP

		blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PostMAD<CBlockMatrixTypelist>::PostMAD_OuterLoop(
			r_t_col, p_dest_vector, p_src_vector, m_block_rows_list);
		// execute templated middle loop selector
	}
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
void CUberBlockMatrix::PreMultiply_Add_FBS(double *p_dest_vector, size_t UNUSED(n_dest_size),
	const double *p_src_vector, size_t UNUSED(n_src_size)) const
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	PreMultiply_Add(p_dest_vector, n_dest_size, p_src_vector, n_src_size);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

	_ASSERTE(p_dest_vector != p_src_vector);
	_ASSERTE(n_dest_size == m_n_row_num);
	_ASSERTE(n_src_size == m_n_col_num);

	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;

		blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PreMAD<CBlockMatrixTypelist>::PreMAD_OuterLoop(
			r_t_col, p_dest_vector, p_src_vector, m_block_rows_list);
		// execute templated middle loop selector
	}
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
inline void CUberBlockMatrix::Scale_FBS(double f_scalar)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	ElementwiseUnaryOp_ZeroInvariant(CScaleBy(f_scalar));
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	ElementwiseUnaryOp_ZeroInvariant_FBS<CBlockMatrixTypelist>(CScaleBy(f_scalar));
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
inline void CUberBlockMatrix::Scale_FBS_Parallel(double f_scalar)
{
	ElementwiseUnaryOp_ZeroInvariant_FBS_Parallel<CBlockMatrixTypelist>(CScaleBy(f_scalar));
}

template <class CBlockMatrixTypelist, class COp>
void CUberBlockMatrix::ElementwiseUnaryOp_ZeroInvariant_FBS(COp op)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	ElementwiseUnaryOp_ZeroInvariant(op);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn &r_t_col = m_block_cols_list[i];
		const size_t n_block_width = r_t_col.n_width;
		for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
			const TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
			size_t n_block_height = m_block_rows_list[r_t_block.first].n_height;
			// get block position, size and data

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseUnaryOp<CBlockMatrixTypelist>::ElementwiseUnary_Loop(
				r_t_block.second, n_block_width * n_block_height, op);
			// perform an operation on the block, with fixed block size
		}
	}
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist, class COp>
void CUberBlockMatrix::ElementwiseUnaryOp_ZeroInvariant_FBS_Parallel(COp op)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	ElementwiseUnaryOp_ZeroInvariant_Parallel(op);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

#ifdef _OPENMP
	_ASSERTE(m_block_cols_list.size() <= INT_MAX);
	int n = int(m_block_cols_list.size());
	#pragma omp parallel for default(shared)
	for(int i = 0; i < n; ++ i) {
		const TColumn &r_t_col = m_block_cols_list[i];
#else // _OPENMP
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn &r_t_col = m_block_cols_list[i];
#endif // OPENMP
		const size_t n_block_width = r_t_col.n_width;
		for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
			const TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
			size_t n_block_height = m_block_rows_list[r_t_block.first].n_height;
			// get block position, size and data

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseUnaryOp<CBlockMatrixTypelist>::ElementwiseUnary_Loop(
				r_t_block.second, n_block_width * n_block_height, op);
			// perform an operation on the block, with fixed block size
		}
	}
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class _CBlockMatrixTypelist>
inline bool CUberBlockMatrix::AddTo_FBS(CUberBlockMatrix &r_dest) const // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp_ZeroInvariant(r_dest, f_Add);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp_ZeroInvariant_FBS<_CBlockMatrixTypelist>(r_dest, f_Add);
	// t_odo - optimize away / generalize as ElementwiseBinaryOp<COp>
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
inline bool CUberBlockMatrix::AddTo_FBS(CUberBlockMatrix &r_dest, double f_factor) const // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp_RightSideZeroInvariant(r_dest, CAddWeighted(f_factor));
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp_RightSideZeroInvariant_FBS<CBlockMatrixTypelist>(r_dest,
		CAddWeighted(f_factor));
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
inline bool CUberBlockMatrix::AddTo_FBS(CUberBlockMatrix &r_dest, double f_factor_dest, double f_factor_this) const // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp(r_dest, CWeightedAddWeighted(f_factor_dest, f_factor_this));
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp_FBS<CBlockMatrixTypelist>(r_dest,
		CWeightedAddWeighted(f_factor_dest, f_factor_this));
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist, class CBinaryOp>
bool CUberBlockMatrix::ElementwiseBinaryOp_ZeroInvariant_FBS(CUberBlockMatrix &r_dest, CBinaryOp op) const // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp_ZeroInvariant(r_dest, op);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);
	r_dest.CheckIntegrity(true);

	// t_odo - decide whether to optimize for addition of matrices with exactly the same block layout (no need to allocate anything, but will anyone use that?) // no
	// t_odo - optimize for matrices with no empty blocks / columns (such as in SLAM), that probably need to be another class // can't see what the optimization was anymore
	// note - this one is a bit hard, as the layout might be completely different in case there are dummy rows / columns,
	//		it is therefore hard to check the layout, and b_EqualStructure() doesn't help. maybe need something like b_CompatibleStructure()

	if(r_dest.m_n_row_num != m_n_row_num || r_dest.m_n_col_num != m_n_col_num)
		return false;
	// the dimensions must be the same

	const std::vector<TRow> &r_row_list_first = m_block_rows_list;
	const std::vector<TColumn> &r_column_list_first = m_block_cols_list;
	std::vector<TRow> &r_row_list_second = r_dest.m_block_rows_list;
	std::vector<TColumn> &r_column_list_second = r_dest.m_block_cols_list;

	std::vector<size_t> row_mapping;
	std::vector<size_t> column_mapping; // t_odo - pack this code to a function
	{
		if(!Build_AdditionLayouts(r_row_list_first, r_column_list_first,
		   r_row_list_second, r_column_list_second, row_mapping, column_mapping))
			return false;
		//r_dest.CheckIntegrity(); // make sure that this operation didn't damage matrix integrity
	}
	// reorganize the destination matrix so that the layout is compatible with this (if required)
	// the (amortized) complexity is linear, about O(row blocks + col blocks)

	for(size_t i = 0, n = r_column_list_first.size(); i < n; ++ i) {
		const TColumn &r_t_col = r_column_list_first[i];
		if(r_t_col.block_list.empty())
			continue;
		// skip empty columns ...

		size_t n_dest = column_mapping[i];
		if(n_dest == size_t(-1))
			return false; // this column was split, the blocks in it are now homeless
		TColumn &r_t_dest_col = r_column_list_second[n_dest];
		// find dest column

		_TyBlockIter p_first_it = r_t_dest_col.block_list.begin();
		// where to begin searching for blocks

		size_t j = 0, m = r_t_col.block_list.size();
		if(!r_t_dest_col.block_list.empty()) {
			for(; j < m; ++ j) {
				const TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
				size_t n_old_row, n_new_row;
				if((n_new_row = row_mapping[n_old_row = r_t_block.first]) == size_t(-1))
					return false; // this row was split, the blocks in it are now homeless
				// get block and its new row

				double *p_value_src = r_t_block.second;
				size_t n_block_size = r_t_col.n_width * m_block_rows_list[n_old_row].n_height;
				// get block data

				_TyBlockIter p_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
					std::lower_bound(p_first_it, r_t_dest_col.block_list.end(), n_new_row, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(p_first_it, r_t_dest_col.block_list.end(), n_new_row, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				// find where to put the block in the column

				bool b_insert_at_end;
				if((b_insert_at_end = (p_block_it == r_t_dest_col.block_list.end())) ||
				   (*p_block_it).first != n_new_row) {
					// a new block

					double *p_value = r_dest.p_Get_DenseStorage(n_block_size);
					blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_UninitializedCopy(
						p_value, p_value_src, n_block_size, op);
					// create a new block, initialize with values

					p_first_it = r_t_dest_col.block_list.insert(p_block_it,
						TColumn::TBlockEntry(n_new_row, p_value));
					// add it to the list, remember iterator

					++ p_first_it;
					// next time search from here (look out for iterator invalidation)

					if(b_insert_at_end) {
						_ASSERTE(p_first_it == r_t_dest_col.block_list.end()); // added to the end
						++ j; // don't forget to count this block as processed
						break; // blocks are sorted, will be always adding to the end from now on
					}
					// in case inserted at end, can continue using the optimized loop below
				} else {
					double *p_value_dest = (*p_block_it).second;
					// get existing block data

					blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop(p_value_dest,
						p_value_src, n_block_size, op);
					// add values to an existing block

					p_first_it = p_block_it + 1;
					// next time, search from here
				}
			}
		}
		for(; j < m; ++ j) {
			const TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
			size_t n_old_row, n_new_row;
			if((n_new_row = row_mapping[n_old_row = r_t_block.first]) == size_t(-1))
				return false; // this row was split, the blocks in it are now homeless

			double *p_value_src = r_t_block.second;
			size_t n_block_size = r_t_col.n_width * m_block_rows_list[n_old_row].n_height;
			// get block data

			_ASSERTE(r_t_dest_col.block_list.empty() || r_t_dest_col.block_list.back().first < n_new_row);
			// make sure the new block comes at the end of the row

			double *p_value = r_dest.p_Get_DenseStorage(n_block_size);
			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_UninitializedCopy(
				p_value, p_value_src, n_block_size, op);
			// create a new block, initialize with values

			r_t_dest_col.block_list.push_back(TColumn::TBlockEntry(n_new_row, p_value));
			// add it to the list
		}
		// merge blocks
	}

	//r_dest.CheckIntegrity(true);
	// make sure that this operation didn't damage matrix integrity

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist, class CBinaryOp>
bool CUberBlockMatrix::ElementwiseBinaryOp_RightSideZeroInvariant_FBS(CUberBlockMatrix &r_dest, CBinaryOp op) const // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp_RightSideZeroInvariant(r_dest, op);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);
	r_dest.CheckIntegrity(true);

	if(r_dest.m_n_row_num != m_n_row_num || r_dest.m_n_col_num != m_n_col_num)
		return false;
	// the dimensions must be the same

	const std::vector<TRow> &r_row_list_first = m_block_rows_list;
	const std::vector<TColumn> &r_column_list_first = m_block_cols_list;
	std::vector<TRow> &r_row_list_second = r_dest.m_block_rows_list;
	std::vector<TColumn> &r_column_list_second = r_dest.m_block_cols_list;

	std::vector<size_t> row_mapping;
	std::vector<size_t> column_mapping; // t_odo - pack this code to a function
	{
		if(!Build_AdditionLayouts(r_row_list_first, r_column_list_first,
		   r_row_list_second, r_column_list_second, row_mapping, column_mapping))
			return false;
		//r_dest.CheckIntegrity(); // make sure that this operation didn't damage matrix integrity
	}
	// reorganize the destination matrix so that the layout is compatible with this (if required)
	// the (amortized) complexity is linear, about O(row blocks + col blocks)

	for(size_t i = 0, n = r_column_list_first.size(); i < n; ++ i) {
		const TColumn &r_t_col = r_column_list_first[i];
		if(r_t_col.block_list.empty())
			continue;
		// skip empty columns ...

		size_t n_dest = column_mapping[i];
		if(n_dest == size_t(-1))
			return false; // this column was split, the blocks in it are now homeless
		TColumn &r_t_dest_col = r_column_list_second[n_dest];
		// find dest column

		_TyBlockIter p_first_it = r_t_dest_col.block_list.begin();
		// where to begin searching for blocks

		size_t j = 0, m = r_t_col.block_list.size();
		if(!r_t_dest_col.block_list.empty()) {
			for(; j < m; ++ j) {
				const TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
				size_t n_old_row, n_new_row;
				if((n_new_row = row_mapping[n_old_row = r_t_block.first]) == size_t(-1))
					return false; // this row was split, the blocks in it are now homeless
				// get block and its new row

				double *p_value_src = r_t_block.second;
				size_t n_block_size = r_t_col.n_width * m_block_rows_list[n_old_row].n_height;
				// get block data

				_TyBlockIter p_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
					std::lower_bound(p_first_it, r_t_dest_col.block_list.end(), n_new_row, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(p_first_it, r_t_dest_col.block_list.end(), n_new_row, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				// find where to put the block in the column

				bool b_insert_at_end;
				if((b_insert_at_end = (p_block_it == r_t_dest_col.block_list.end())) ||
				   (*p_block_it).first != n_new_row) {
					// a new block

					double *p_value = r_dest.p_Get_DenseStorage(n_block_size);
					blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_UninitLeftSide(p_value,
						p_value_src, n_block_size, op);
					// create a new block, initialize with values

					p_first_it = r_t_dest_col.block_list.insert(p_block_it, TColumn::TBlockEntry(n_new_row, p_value));
					// add it to the list, remember iterator

					++ p_first_it;
					// next time search from here (look out for iterator invalidation)

					if(b_insert_at_end) {
						_ASSERTE(p_first_it == r_t_dest_col.block_list.end()); // added to the end
						++ j; // don't forget to count this block as processed
						break; // blocks are sorted, will be always adding to the end from now on
					}
					// in case inserted at end, can continue using the optimized loop below
				} else {
					double *p_value_dest = (*p_block_it).second;
					// get existing block data

					blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop(p_value_dest,
						p_value_src, n_block_size, op);
					// add values to an existing block

					p_first_it = p_block_it + 1;
					// next time, search from here
				}
			}
		}
		for(; j < m; ++ j) {
			const TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
			size_t n_old_row, n_new_row;
			if((n_new_row = row_mapping[n_old_row = r_t_block.first]) == size_t(-1))
				return false; // this row was split, the blocks in it are now homeless

			double *p_value_src = r_t_block.second;
			size_t n_block_size = r_t_col.n_width * m_block_rows_list[n_old_row].n_height;
			// get block data

			_ASSERTE(r_t_dest_col.block_list.empty() || r_t_dest_col.block_list.back().first < n_new_row);
			// make sure the new block comes at the end of the row

			double *p_value = r_dest.p_Get_DenseStorage(n_block_size);
			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_UninitLeftSide(
				p_value, p_value_src, n_block_size, op);
			// create a new block, initialize with values

			r_t_dest_col.block_list.push_back(TColumn::TBlockEntry(n_new_row, p_value));
			// add it to the list
		}
		// merge blocks
	}

	//r_dest.CheckIntegrity(true);
	// make sure that this operation didn't damage matrix integrity

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist, class CBinaryOp>
bool CUberBlockMatrix::ElementwiseBinaryOp_FBS(CUberBlockMatrix &r_dest, CBinaryOp op) const // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return ElementwiseBinaryOp(r_dest, op);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);
	r_dest.CheckIntegrity(true);

	if(r_dest.m_n_row_num != m_n_row_num || r_dest.m_n_col_num != m_n_col_num)
		return false;
	// the dimensions must be the same

	const std::vector<TRow> &r_row_list_first = m_block_rows_list;
	const std::vector<TColumn> &r_column_list_first = m_block_cols_list;
	std::vector<TRow> &r_row_list_second = r_dest.m_block_rows_list;
	std::vector<TColumn> &r_column_list_second = r_dest.m_block_cols_list;

	std::vector<size_t> row_mapping;
	std::vector<size_t> column_mapping; // t_odo - pack this code to a function
	{
		if(!Build_AdditionLayouts(r_row_list_first, r_column_list_first,
		   r_row_list_second, r_column_list_second, row_mapping, column_mapping))
			return false;
		//r_dest.CheckIntegrity(); // make sure that this operation didn't damage matrix integrity
	}
	// reorganize the destination matrix so that the layout is compatible with this (if required)
	// the (amortized) complexity is linear, about O(row blocks + col blocks)

	size_t n_last_dest_col = 0;
	for(size_t i = 0, n = r_column_list_first.size(); i < n; ++ i) {
		const TColumn &r_t_col = r_column_list_first[i];
		if(r_t_col.block_list.empty())
			continue;
		// skip empty columns ...

		size_t n_dest = column_mapping[i];
		if(n_dest == size_t(-1))
			return false; // this column was split, the blocks in it are now homeless
		TColumn &r_t_dest_col = r_column_list_second[n_dest];
		// find dest column

		for(; n_last_dest_col != n_dest; ++ n_last_dest_col) {
			TColumn &r_t_dest_col = r_column_list_second[n_last_dest_col];
			for(size_t j = 0, m = r_t_dest_col.block_list.size(); j < m; ++ j) {
				const TColumn::TBlockEntry &r_t_block = r_t_dest_col.block_list[j];

				double *p_value_dest = r_t_block.second;
				size_t n_block_size_dest = r_t_dest_col.n_width * r_row_list_second[r_t_block.first].n_height;
				// get block data

				blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_NullRightSide(
					p_value_dest, n_block_size_dest, op);
			}
		}
		// don't forget to modify all of the previous columns we skipped in dest

		_ASSERTE(n_last_dest_col == n_dest);
		++ n_last_dest_col;
		// and don't modify this column either

		_TyBlockIter p_first_it = r_t_dest_col.block_list.begin();
		// where to begin searching for blocks

		size_t j = 0, m = r_t_col.block_list.size();
		if(!r_t_dest_col.block_list.empty()) {
			for(; j < m; ++ j) {
				const TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
				size_t n_old_row, n_new_row;
				if((n_new_row = row_mapping[n_old_row = r_t_block.first]) == size_t(-1))
					return false; // this row was split, the blocks in it are now homeless
				// get block and its new row

				double *p_value_src = r_t_block.second;
				size_t n_block_size = r_t_col.n_width * m_block_rows_list[n_old_row].n_height;
				// get block data

				_TyBlockIter p_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
					std::lower_bound(p_first_it, r_t_dest_col.block_list.end(), n_new_row, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(p_first_it, r_t_dest_col.block_list.end(), n_new_row, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				// find where to put the block in the column

				for(; p_first_it != p_block_it; ++ p_first_it) {
					size_t n_dest_row = (*p_first_it).first;
					size_t n_block_size_dest = r_t_col.n_width * r_row_list_second[n_dest_row].n_height;
					double *p_value_dest = (*p_first_it).second;

					blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_NullRightSide(
						p_value_dest, n_block_size_dest, op);
				}
				// don't forget to modify all the blocks in dest between last one and the one being added to

				bool b_insert_at_end;
				if((b_insert_at_end = (p_block_it == r_t_dest_col.block_list.end())) ||
				   (*p_block_it).first != n_new_row) {
					// a new block

					double *p_value = r_dest.p_Get_DenseStorage(n_block_size);
					double *p_value_dest = p_value;
					blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_UninitLeftSide(
						p_value_dest, p_value_src, n_block_size, op);
					// create a new block, initialize with values

					p_first_it = r_t_dest_col.block_list.insert(p_block_it, TColumn::TBlockEntry(n_new_row, p_value));
					// add it to the list, remember iterator

					++ p_first_it;
					// next time search from here (look out for iterator invalidation)

					if(b_insert_at_end) {
						_ASSERTE(p_first_it == r_t_dest_col.block_list.end()); // added to the end
						++ j; // don't forget to count this block as processed
						break; // blocks are sorted, will be always adding to the end from now on
					}
					// in case inserted at end, can continue using the optimized loop below
				} else {
					double *p_value_dest = (*p_block_it).second;
					// get existing block data

					blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop(p_value_dest,
						p_value_src, n_block_size, op);
					// add values to an existing block

					p_first_it = p_block_it + 1;
					// next time, search from here
				}
			}
		}
		{
			_TyBlockIter p_end_it = r_t_dest_col.block_list.end();
			for(; p_first_it != p_end_it; ++ p_first_it) {
				size_t n_dest_row = (*p_first_it).first;
				size_t n_block_size_dest = r_t_col.n_width * r_row_list_second[n_dest_row].n_height;
				double *p_value_dest = (*p_first_it).second;
				blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_NullRightSide(
					p_value_dest, n_block_size_dest, op);
			}
			// don't forget to modify all the blocks in dest between last added and the end of the list
		}
		for(; j < m; ++ j) {
			const TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
			size_t n_old_row, n_new_row;
			if((n_new_row = row_mapping[n_old_row = r_t_block.first]) == size_t(-1))
				return false; // this row was split, the blocks in it are now homeless

			double *p_value_src = r_t_block.second;
			size_t n_block_size = r_t_col.n_width * m_block_rows_list[n_old_row].n_height;
			// get block data

			_ASSERTE(r_t_dest_col.block_list.empty() || r_t_dest_col.block_list.back().first < n_new_row);
			// make sure the new block comes at the end of the row

			double *p_value = r_dest.p_Get_DenseStorage(n_block_size);
			double *p_value_dest = p_value;
			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_UninitLeftSide(
				p_value_dest, p_value_src, n_block_size, op);
			// create a new block, initialize with values

			r_t_dest_col.block_list.push_back(TColumn::TBlockEntry(n_new_row, p_value));
			// add it to the list
		}
		// merge blocks
	}

	for(size_t n = r_column_list_second.size(); n_last_dest_col != n; ++ n_last_dest_col) {
		TColumn &r_t_dest_col = r_column_list_second[n_last_dest_col];
		for(size_t j = 0, m = r_t_dest_col.block_list.size(); j < m; ++ j) {
			const TColumn::TBlockEntry &r_t_block = r_t_dest_col.block_list[j];

			double *p_value_dest = r_t_block.second;
			size_t n_block_size_dest = r_t_dest_col.n_width * r_row_list_second[r_t_block.first].n_height;
			// get block data

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_ElementwiseBinaryOp<CBlockMatrixTypelist>::ElementwiseBinary_Loop_NullRightSide(
				p_value_dest, n_block_size_dest, op);
		}
	}
	// don't forget to modify all of the columns after the last modified one

	//r_dest.CheckIntegrity(true);
	// make sure that this operation didn't damage matrix integrity

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelistA, class CBlockMatrixTypelistB>
inline bool CUberBlockMatrix::ProductOf_FBS(const CUberBlockMatrix &r_A, const CUberBlockMatrix &r_B) // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return r_A.MultiplyToWith(*this, r_B);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return r_A.MultiplyToWith_FBS<CBlockMatrixTypelistA, CBlockMatrixTypelistB>(*this, r_B);
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelistA, class CBlockMatrixTypelistB>
inline bool CUberBlockMatrix::ProductOf_FBS(const CUberBlockMatrix &r_A, const CUberBlockMatrix &r_B, bool b_upper_diag_only) // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return r_A.MultiplyToWith(*this, r_B, b_upper_diag_only);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return r_A.MultiplyToWith_FBS<CBlockMatrixTypelistA, CBlockMatrixTypelistB>(*this, r_B, b_upper_diag_only);
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelistThis, class CBlockMatrixTypelistOther>
bool CUberBlockMatrix::MultiplyToWith_FBS(CUberBlockMatrix &r_dest, const CUberBlockMatrix &r_other) const // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return MultiplyToWith(r_dest, r_other);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);
	r_other.CheckIntegrity(true);

	typedef CBlockMatrixTypelistThis tlA;
	typedef CBlockMatrixTypelistOther tlB;

	if(m_n_col_num != r_other.m_n_row_num)
		return false;
	// check common dimension

	const CUberBlockMatrix &A = *this;
	const CUberBlockMatrix &B = r_other;
	// name the matrices; the operation is r_dest = A * B

	const std::vector<TRow> &r_row_list_A = A.m_block_rows_list;
	const std::vector<TRow> &r_row_list_B = B.m_block_rows_list;
	std::vector<TRow> &r_row_list_dest = r_dest.m_block_rows_list;
	const std::vector<TColumn> &r_col_list_A = A.m_block_cols_list;
	const std::vector<TColumn> &r_col_list_B = B.m_block_cols_list;
	std::vector<TColumn> &r_col_list_dest = r_dest.m_block_cols_list;
	// name the cumsums for easier access

	size_t n_dest_col_num = B.m_n_col_num; // t_odo - these and some others do not follow the naming conventions; take care of that
	size_t n_dest_row_num = A.m_n_row_num;
	// these are the dimensions of the new matrix

	r_dest.Clear();
	r_row_list_dest = r_row_list_A; // copy row layout from this
	r_col_list_dest.resize(r_col_list_B.size());
	//std::for_each(r_col_list_dest.begin(), r_col_list_dest.end(), CColumnCumsumCopy(r_col_list_B.begin())); // copy column layout but not the blocks
	std::transform(r_col_list_B.begin(), r_col_list_B.end(), r_col_list_dest.begin(), t_ColumnCumsumCopy);
	r_dest.m_n_col_num = n_dest_col_num;
	r_dest.m_n_row_num = n_dest_row_num;
	// create layout for the destination matrix (linear time)

	//r_dest.CheckIntegrity();
	// makes sure the dest layout is ok

	std::vector<size_t> reindex_rows_B_to_cols_A;
	{
		std::vector<size_t> reindex_cols_A;
		size_t n_common_size = n_MergeLayout(r_col_list_A, r_row_list_B, reindex_cols_A, reindex_rows_B_to_cols_A);
		// merge matrix layouts (linear time in max(A column blocks, B row blocks) or something like that)

		std::vector<size_t> common(n_common_size, size_t(-2)); // helper array (note the use of -2 !!)
		for(size_t i = 0, n = reindex_cols_A.size(); i < n; ++ i) {
			size_t n_common_A;
			if((n_common_A = reindex_cols_A[i]) == size_t(-1))
				continue;
			_ASSERTE(n_common_A < common.size());
			common[n_common_A] = i;
		}
		// create inverse mapping for columns of A (linear time in number of column blocks of A)

		if(reindex_cols_A.capacity() < reindex_rows_B_to_cols_A.size())
			reindex_cols_A.clear(); // to prevent resize() on the next line from copying the data (that are going to be overwritten anyway)
		reindex_cols_A.resize(reindex_rows_B_to_cols_A.size()); // reuse this array as output (may not actually resize in most cases)
		for(size_t i = 0, n = reindex_rows_B_to_cols_A.size(); i < n; ++ i) {
			size_t n_common_B;
			if((n_common_B = reindex_rows_B_to_cols_A[i]) == size_t(-1)) {
				reindex_cols_A[i] = -1; // !!
				continue;
			}
			_ASSERTE(n_common_B < common.size());
			reindex_cols_A[i] = common[n_common_B];
		}
		reindex_cols_A.swap(reindex_rows_B_to_cols_A); // swap with the array we want output in
		// map inverse mapping of A to B (linear time in number of row blocks of B)
	}
	// merge the common dimension layout (linear time)
	// -1 means the row did not map/exist in B, -2 meants it did not map/exist in A

	// this version have output to transpose matrix, then transpose that to dest, block lookup is not needed

#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
	std::vector<size_t> cols_load_list(r_col_list_dest.size(), 0);
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
	std::vector<std::vector<TColumn::TBlockEntry> > transpose_cols_list(r_row_list_dest.size());
	// list for storing transpose columns (note that the sizes and cumsums are not initialized and are invalid)

	_TyDenseAllocator alloc(r_dest.m_data_pool);
	// get allocator in dest

	size_t n_column_id_B = 0;
	for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(), p_col_B_end_it = r_col_list_B.end();
	   p_col_B_it != p_col_B_end_it; ++ p_col_B_it, ++ n_column_id_B) {
		const TColumn &r_t_column_B = *p_col_B_it;
		// for each column in B (index is n_column_id_B)

		//if(r_t_column_B.block_list.empty()) // better enter decission tree; empty columns are rare in realworld applications
		//	continue;
		// otherwise width mismatch occurs

		if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_MatrixMultiply<tlA, tlB>::MatrixMultiply_OuterLoop(
		   r_t_column_B, n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
		   transpose_cols_list, reindex_rows_B_to_cols_A, alloc))
			return false;
		// the first decision tree for r_t_column_B.n_width (decides all block widths)
#if 0
		for(_TyBlockConstIter p_block_B_it = r_t_column_B.block_list.begin(),
		   p_block_B_end_it = r_t_column_B.block_list.end(); p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
			const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
			size_t n_row_id_B = r_t_block_B.first; // get row of the block
			// for each block in the current column in B

			// the second decision tree for r_row_list_B[n_row_id_B].n_height (decides block B heights, given constant block B width)
			// note that this is the common dimension

			size_t n_bmB_cols;
			_TyMatrixXdRef blockB(r_t_block_B.second, r_row_list_B[n_row_id_B].n_height,
				n_bmB_cols = r_t_column_B.n_width);
			// create map to block B data

			size_t n_column_id_A = reindex_rows_B_to_cols_A[n_row_id_B];
			_ASSERTE(size_t(-1) > size_t(-2)); // just to make sure the next line is correct
			if(n_column_id_A >= size_t(-2)) {
				if(n_column_id_A == -1)
					return false; // didn't map from B to common and we know it was not empty (we have a block here)
				continue; // do not fail, it might also mean that there are no blocks in that column in A and hence the result of multiplication is zero
			}
			_ASSERTE(n_column_id_A < r_col_list_A.size());
			const TColumn &r_t_column_A = r_col_list_A[n_column_id_A];
			_ASSERTE(r_t_column_A.n_width == r_row_list_B[n_row_id_B].n_height);
			// lookup which column in A corresponds with current block row in B

			for(_TyBlockConstIter p_block_A_it = r_t_column_A.block_list.begin(),
			   p_block_A_end_it = r_t_column_A.block_list.end(); p_block_A_it != p_block_A_end_it; ++ p_block_A_it) {
				const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
				size_t n_row_id_A = r_t_block_A.first; // get row of the block
				// for each block in the current column in A

				// the third decision tree for r_row_list_A[n_row_id_A].n_height (decides block A heights, given block A width == block B height)
				// note that r_t_column_A.n_width must equal r_row_list_B[n_row_id_B].n_height (the second decision tree)

				size_t n_bmA_rows;
				_TyMatrixXdRef blockA(r_t_block_A.second,
					n_bmA_rows = r_row_list_A[n_row_id_A].n_height, r_t_column_A.n_width);
				// create map to block A data

				// multiplication of blockA * blockB yields block at (n_row_id_A, n_column_id_B)

				_ASSERTE(n_row_id_A < r_row_list_dest.size());
				_ASSERTE(n_column_id_B < r_col_list_dest.size());
				//size_t n_prod_rows = r_row_list_dest[n_row_id_A].n_height;
				//size_t n_prod_cols = r_col_list_dest[n_column_id_B].n_width; // unused
				_ASSERTE(r_row_list_dest[n_row_id_A].n_height == blockA.rows());
				_ASSERTE(r_col_list_dest[n_column_id_B].n_width == blockB.cols());
				if(blockA.cols() != blockB.rows())
					return false; // make sure the blocks are multiplicable (not able to verify by just merging the layout)
				// basic checks about matrix dimensions

				_ASSERTE(n_row_id_A < transpose_cols_list.size());
				std::vector<TColumn::TBlockEntry> &r_transpose_column = transpose_cols_list[n_row_id_A];
				if(r_transpose_column.empty() || r_transpose_column.back().first < n_column_id_B) {
					double *p_new_block_data = alloc.p_Get_DenseStorage(n_bmA_rows * n_bmB_cols);
					// get storage

					_TyMatrixXdRef block_dest(p_new_block_data, n_bmA_rows, n_bmB_cols);
					block_dest = blockA * blockB;
					// initialize a new block

					r_transpose_column.push_back(TColumn::TBlockEntry(n_column_id_B, p_new_block_data));
					// add it to the list

#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
					_ASSERTE(n_column_id_B < cols_load_list.size());
					++ cols_load_list[n_column_id_B];
					// we have a list of numbers of entries per row, that can be done in linear time and it might help speeding up the final matrix transpose later
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
				} else {
					_ASSERTE(r_transpose_column.back().first == n_column_id_B); // n_column_id_B is monotonically increasing, no matter what, don't have to do log(N) lookup here ;)
					double *p_block_data = r_transpose_column.back().second;
					_TyMatrixXdRef block_dest(p_block_data, n_bmA_rows, n_bmB_cols);
					block_dest.noalias() += blockA * blockB;
					// add to the existing block (run the dot sum)
				}
				// perform the dense multiplication using reference matrices and eigen
			}
		}
#endif // 0
	}
	// performs sparse matrix multiplication (linear time in number of blocks * constant time to insert the blocks)

#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
	_ASSERTE(cols_load_list.size() == r_col_list_dest.size());
	for(size_t i = 0, n = cols_load_list.size(); i < n; ++ i)
		r_col_list_dest[i].block_list.reserve(cols_load_list[i]);
	// allocate block lists
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS

	for(size_t i = 0, n = transpose_cols_list.size(); i < n; ++ i) {
		const std::vector<TColumn::TBlockEntry> &r_col = transpose_cols_list[i];
		for(size_t j = 0, m = r_col.size(); j < m; ++ j) {
			const TColumn::TBlockEntry &r_block = r_col[j];
			_ASSERTE(r_block.first < r_col_list_dest.size());
#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
			_ASSERTE(r_col_list_dest[r_block.first].block_list.capacity() > r_col_list_dest[r_block.first].block_list.size()); // since it was preallocated
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
			r_col_list_dest[r_block.first].block_list.push_back(TColumn::TBlockEntry(i, r_block.second));
		}
	}
	// performs the final transpose (linear time in number of blocks)

	//r_dest.CheckIntegrity(true);
	// makes sure the dest matrix is ok

	// note it might be beneficial (and possibly quite easy) to implement Strassen's algorithm here.
	// would be a nice work for a paper, too.

	// no need to implement more cunning algorithms, though (Wiki):
	// The Coppersmith–Winograd algorithm is frequently used as a building block in other algorithms to
	// prove theoretical time bounds. However, unlike the Strassen algorithm, it is not used in practice
	// because it only provides an advantage for matrices so large that they cannot be processed by modern
	// hardware.

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelistThis, class CBlockMatrixTypelistOther>
bool CUberBlockMatrix::MultiplyToWith_FBS(CUberBlockMatrix &r_dest,
	const CUberBlockMatrix &r_other, bool b_upper_diag_only) const // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return MultiplyToWith(r_dest, r_other);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);
	r_other.CheckIntegrity(true);

	typedef CBlockMatrixTypelistThis tlA;
	typedef CBlockMatrixTypelistOther tlB;

	if(m_n_col_num != r_other.m_n_row_num)
		return false;
	// check common dimension

	const CUberBlockMatrix &A = *this;
	const CUberBlockMatrix &B = r_other;
	// name the matrices; the operation is r_dest = A * B

	const std::vector<TRow> &r_row_list_A = A.m_block_rows_list;
	const std::vector<TRow> &r_row_list_B = B.m_block_rows_list;
	std::vector<TRow> &r_row_list_dest = r_dest.m_block_rows_list;
	const std::vector<TColumn> &r_col_list_A = A.m_block_cols_list;
	const std::vector<TColumn> &r_col_list_B = B.m_block_cols_list;
	std::vector<TColumn> &r_col_list_dest = r_dest.m_block_cols_list;
	// name the cumsums for easier access

	size_t n_dest_col_num = B.m_n_col_num; // t_odo - these and some others do not follow the naming conventions; take care of that
	size_t n_dest_row_num = A.m_n_row_num;
	// these are the dimensions of the new matrix

	r_dest.Clear();
	r_row_list_dest = r_row_list_A; // copy row layout from this
	r_col_list_dest.resize(r_col_list_B.size());
	//std::for_each(r_col_list_dest.begin(), r_col_list_dest.end(), CColumnCumsumCopy(r_col_list_B.begin())); // copy column layout but not the blocks
	std::transform(r_col_list_B.begin(), r_col_list_B.end(), r_col_list_dest.begin(), t_ColumnCumsumCopy);
	r_dest.m_n_col_num = n_dest_col_num;
	r_dest.m_n_row_num = n_dest_row_num;
	// create layout for the destination matrix (linear time)

	//r_dest.CheckIntegrity();
	// makes sure the dest layout is ok

	std::vector<size_t> reindex_rows_B_to_cols_A;
	{
		std::vector<size_t> reindex_cols_A;
		size_t n_common_size = n_MergeLayout(r_col_list_A, r_row_list_B, reindex_cols_A, reindex_rows_B_to_cols_A);
		// merge matrix layouts (linear time in max(A column blocks, B row blocks) or something like that)

		std::vector<size_t> common(n_common_size, size_t(-2)); // helper array (note the use of -2 !!)
		for(size_t i = 0, n = reindex_cols_A.size(); i < n; ++ i) {
			size_t n_common_A;
			if((n_common_A = reindex_cols_A[i]) == size_t(-1))
				continue;
			_ASSERTE(n_common_A < common.size());
			common[n_common_A] = i;
		}
		// create inverse mapping for columns of A (linear time in number of column blocks of A)

		if(reindex_cols_A.capacity() < reindex_rows_B_to_cols_A.size())
			reindex_cols_A.clear(); // to prevent resize() on the next line from copying the data (that are going to be overwritten anyway)
		reindex_cols_A.resize(reindex_rows_B_to_cols_A.size()); // reuse this array as output (may not actually resize in most cases)
		for(size_t i = 0, n = reindex_rows_B_to_cols_A.size(); i < n; ++ i) {
			size_t n_common_B;
			if((n_common_B = reindex_rows_B_to_cols_A[i]) == size_t(-1)) {
				reindex_cols_A[i] = -1; // !!
				continue;
			}
			_ASSERTE(n_common_B < common.size());
			reindex_cols_A[i] = common[n_common_B];
		}
		reindex_cols_A.swap(reindex_rows_B_to_cols_A); // swap with the array we want output in
		// map inverse mapping of A to B (linear time in number of row blocks of B)
	}
	// merge the common dimension layout (linear time)
	// -1 means the row did not map/exist in B, -2 meants it did not map/exist in A

	// this version have output to transpose matrix, then transpose that to dest, block lookup is not needed

#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
	std::vector<size_t> cols_load_list(r_col_list_dest.size(), 0);
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
	std::vector<std::vector<TColumn::TBlockEntry> > transpose_cols_list(r_row_list_dest.size());
	// list for storing transpose columns (note that the sizes and cumsums are not initialized and are invalid)

	_TyDenseAllocator alloc(r_dest.m_data_pool);
	// get allocator in dest

	size_t n_column_id_B = 0;
	if(b_upper_diag_only) {
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(), p_col_B_end_it = r_col_list_B.end();
		   p_col_B_it != p_col_B_end_it; ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			//if(r_t_column_B.block_list.empty()) // better enter decission tree; empty columns are rare in realworld applications
			//	continue;
			// otherwise width mismatch occurs

			if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_MatrixMultiply<tlA, tlB>::MatrixMultiply_OuterLoop_UpperTriag(
			   r_t_column_B, n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
			   transpose_cols_list, reindex_rows_B_to_cols_A, alloc))
				return false;
			// the first decision tree for r_t_column_B.n_width (decides all block widths)
		}
		// performs sparse matrix multiplication (linear time in number of blocks * constant time to insert the blocks)
	} else {
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(), p_col_B_end_it = r_col_list_B.end();
		   p_col_B_it != p_col_B_end_it; ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			//if(r_t_column_B.block_list.empty()) // better enter decission tree; empty columns are rare in realworld applications
			//	continue;
			// otherwise width mismatch occurs

			if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_MatrixMultiply<tlA, tlB>::MatrixMultiply_OuterLoop(
			   r_t_column_B, n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
			   transpose_cols_list, reindex_rows_B_to_cols_A, alloc))
				return false;
		}
		// performs sparse matrix multiplication (linear time in number of blocks * constant time to insert the blocks)
	}

#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
	_ASSERTE(cols_load_list.size() == r_col_list_dest.size());
	for(size_t i = 0, n = cols_load_list.size(); i < n; ++ i)
		r_col_list_dest[i].block_list.reserve(cols_load_list[i]);
	// allocate block lists
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS

	for(size_t i = 0, n = transpose_cols_list.size(); i < n; ++ i) {
		const std::vector<TColumn::TBlockEntry> &r_col = transpose_cols_list[i];
		for(size_t j = 0, m = r_col.size(); j < m; ++ j) {
			const TColumn::TBlockEntry &r_block = r_col[j];
			_ASSERTE(r_block.first < r_col_list_dest.size());
#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
			_ASSERTE(r_col_list_dest[r_block.first].block_list.capacity() > r_col_list_dest[r_block.first].block_list.size()); // since it was preallocated
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
			r_col_list_dest[r_block.first].block_list.push_back(TColumn::TBlockEntry(i, r_block.second));
		}
	}
	// performs the final transpose (linear time in number of blocks)

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
void CUberBlockMatrix::PreMultiplyWithSelfTransposeTo_FBS(CUberBlockMatrix &r_dest,
	bool b_upper_diagonal_only) // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	PreMultiplyWithSelfTransposeTo(r_dest, b_upper_diagonal_only);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

	// if the original multiplication operation was r_dest = A * B, B = this and A = *this^T

	const std::vector<TColumn> &r_row_list_A = m_block_cols_list;
	std::vector<TColumn> col_list_A; // columns of *this^T
	const std::vector<TRow> &r_row_list_B = m_block_rows_list;
	const std::vector<TColumn> &r_col_list_B = m_block_cols_list;
	std::vector<TRow> &r_row_list_dest = r_dest.m_block_rows_list;
	std::vector<TColumn> &r_col_list_dest = r_dest.m_block_cols_list; // t_odo - these do not conform to naming conventions; rename them
	// name the cumsums for easier access

	{
		col_list_A.resize(r_row_list_B.size());
		//std::for_each(col_list_A.begin(), col_list_A.end(), CRowCumsumToColumnCumsum(r_row_list_B.begin()));
		std::transform(r_row_list_B.begin(), r_row_list_B.end(), col_list_A.begin(), t_RowCumsumToColumnCumsum);
		for(size_t i = 0, n = r_col_list_B.size(); i < n; ++ i) {
			const TColumn &r_t_col = r_col_list_B[i];
			for(_TyBlockConstIter p_block_it = r_t_col.block_list.begin(),
			   p_end_it = r_t_col.block_list.end(); p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				col_list_A[r_t_block.first].block_list.push_back(TColumn::TBlockEntry(i, r_t_block.second));
			}
		}
	}
	// calculate fast transpose of this (no rows alloc, no matrix structure,
	// no block transpose, no extra storage, O(n) in number of blocks)

	size_t n_dest_col_num = m_n_col_num; //B.m_n_col_num; // t_odo - these and some others do not follow the naming conventions; take care of that
	size_t n_dest_row_num = m_n_col_num; //A.m_n_row_num; // A = B^T
	// these are the dimensions of the new matrix

	{
		r_dest.Clear();
		r_row_list_dest.resize(r_col_list_B.size()); // copy row layout from this
		r_col_list_dest.resize(r_col_list_B.size()); // column layout is the same, except for the blocks
		//std::for_each(r_row_list_dest.begin(), r_row_list_dest.end(), CColumnCumsumToRowCumsum(r_col_list_B.begin()));
		std::transform(r_col_list_B.begin(), r_col_list_B.end(), r_row_list_dest.begin(), t_ColumnCumsumToRowCumsum);
		//std::for_each(r_col_list_dest.begin(), r_col_list_dest.end(), CColumnCumsumCopy(r_col_list_B.begin())); // copy column layout but not the blocks
		std::transform(r_col_list_B.begin(), r_col_list_B.end(), r_col_list_dest.begin(), t_ColumnCumsumCopy);
		r_dest.m_n_col_num = n_dest_col_num;
		r_dest.m_n_row_num = n_dest_row_num;
		// create layout for the destination matrix (linear time)

		//r_dest.CheckIntegrity();
		// makes sure the dest layout is ok
	}
	// create dest structure

	// no need to reindex / merge, all the indices match and exist

	_TyDenseAllocator alloc(r_dest.m_data_pool);
	// make the allocator object for the templates to use

	{
		for(_TyColumnConstIter p_col_A_it = col_list_A.begin(), p_col_A_end_it = col_list_A.end();
		   p_col_A_it != p_col_A_end_it; ++ p_col_A_it) {
			const TColumn &r_t_column_A = *p_col_A_it;
			// for each column in A

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PreATA<CBlockMatrixTypelist>::PreATA_UpperTriangle_OuterLoop(
				r_t_column_A, r_row_list_A, r_col_list_dest, alloc);
			// wrap the outer loop in a column size decission tree
		}
		// perform the multiplication of the blocks above diagonal, need to use log(n) lookup

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok

		// t_odo test and convert the first loop as well.

		// note that this can run in parallel easily without conflicts, except for the block allocator
		size_t n_column_id_B = 0;
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(), p_col_B_end_it = r_col_list_B.end();
		   p_col_B_it != p_col_B_end_it; ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			if(r_t_column_B.block_list.empty())
				continue;
			// only process non-empty columns

			double *p_block_data =
				blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PreATA<CBlockMatrixTypelist>::p_PreATA_Diagonal_OuterLoop(
				r_t_column_B, alloc, r_row_list_B);
			// produce the column dot product, each block with itself

			_ASSERTE(r_col_list_dest[n_column_id_B].block_list.empty() ||
				r_col_list_dest[n_column_id_B].block_list.back().first < n_column_id_B);
			// makes sure that inserting block at diagonal doesn't violate block ordering

			r_col_list_dest[n_column_id_B].block_list.push_back(TColumn::TBlockEntry(n_column_id_B, p_block_data));
			// add the block at the end of the list
		}
		// generate the blocks at diagonal (const time lookup, fast multiplication)

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok
	}

	if(!b_upper_diagonal_only) {
		for(size_t i = 0, n = r_col_list_dest.size(); i < n; ++ i) {
			const TColumn &r_t_col = r_col_list_dest[i];
			const size_t n_col = i;
			for(_TyBlockConstIter p_block_it = r_t_col.block_list.begin(),
			   p_end_it = r_t_col.block_list.end(); p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				const size_t n_row = r_t_block.first;
				_ASSERTE(n_row <= n_col); // make sure the block is above or at the diagonal
				if(n_row < n_col) {
					// it is guaranteed (by the condition above) that this will not alter block i or any of the
					// blocks to be processed next

					_ASSERTE(r_col_list_dest[n_row].block_list.empty() || r_col_list_dest[n_row].block_list.back().first < n_col);
					// make sure that this won't violate block ordering (since i = n_col is monotonically increasing)

					const size_t n_block_width = r_row_list_dest[n_row].n_height, n_block_height = r_t_col.n_width; // these are dims of the dest block
					double *p_block_data = r_dest.p_Get_DenseStorage(n_block_width * n_block_height);
					r_col_list_dest[n_row].block_list.push_back(TColumn::TBlockEntry(n_col, p_block_data));
					// this block is above the diagonal, at (n_row, n_col), create transpose block at (n_col, n_row)

					const double *p_src_data = r_t_block.second;
					/*for(size_t x = 0; x < n_block_width; ++ x)
						for(size_t y = 0; y < n_block_height; ++ y)
							p_block_data[y + n_block_height * x] = p_src_data[x + n_block_width * y];*/ // f_ixme - does the src / dest ordering matter? is it interchangable? // t_odo - test transposing a matrix // this seems to be correct
					_TyConstMatrixXdRef src(p_src_data, n_block_width, n_block_height);
					_TyMatrixXdRef dest(p_block_data, n_block_height, n_block_width);
					dest.noalias() = src.transpose();
					// copy/transpose the block
				} else {
					_ASSERTE(p_block_it + 1 == p_end_it); // diagonal, gotta be the last block
					break; // still can save the comparison
				}
			}
		}
		// transpose the blocks from the triangle above diagonal to the triangle below diagonal
		// this takes O(n) time in number of blocks above the diagonal, or O(n/2) in number of all blocks
	}
	// no need to optimize this, it is not used in slam, but optimize it later, might transpose blocks using SSE (or at least ommit the two innermost loops)

	//r_dest.CheckIntegrity();
	// makes sure the dest matrix is ok
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
void CUberBlockMatrix::PreMultiplyWithSelfTransposeTo_FBS_Parallel(CUberBlockMatrix &r_dest,
	bool b_upper_diagonal_only) // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	PreMultiplyWithSelfTransposeTo/*_Parallel*/(r_dest, b_upper_diagonal_only); // non-fbs parallel version not written (todo?)
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

#ifdef _OPENMP
	// if the original multiplication operation was r_dest = A * B, B = this and A = *this^T

	const std::vector<TColumn> &r_row_list_A = m_block_cols_list;
	std::vector<TColumn> col_list_A; // columns of *this^T
	const std::vector<TRow> &r_row_list_B = m_block_rows_list;
	const std::vector<TColumn> &r_col_list_B = m_block_cols_list;
	std::vector<TRow> &r_row_list_dest = r_dest.m_block_rows_list;
	std::vector<TColumn> &r_col_list_dest = r_dest.m_block_cols_list; // t_odo - these do not conform to naming conventions; rename them
	// name the cumsums for easier access

	{
		col_list_A.resize(r_row_list_B.size());
		//std::for_each(col_list_A.begin(), col_list_A.end(), CRowCumsumToColumnCumsum(r_row_list_B.begin()));
		std::transform(r_row_list_B.begin(), r_row_list_B.end(), col_list_A.begin(), t_RowCumsumToColumnCumsum);
		for(size_t i = 0, n = r_col_list_B.size(); i < n; ++ i) {
			const TColumn &r_t_col = r_col_list_B[i];
			for(_TyBlockConstIter p_block_it = r_t_col.block_list.begin(),
			   p_end_it = r_t_col.block_list.end(); p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				col_list_A[r_t_block.first].block_list.push_back(TColumn::TBlockEntry(i, r_t_block.second));
			}
		}
	}
	// calculate fast transpose of this (no rows alloc, no matrix structure,

	size_t n_dest_col_num = m_n_col_num; //B.m_n_col_num; // t_odo - these and some others do not follow the naming conventions; take care of that
	size_t n_dest_row_num = m_n_col_num; //A.m_n_row_num; // A = B^T
	// these are the dimensions of the new matrix

	{
		r_dest.Clear();
		r_row_list_dest.resize(r_col_list_B.size()); // copy row layout from this
		r_col_list_dest.resize(r_col_list_B.size()); // column layout is the same, except for the blocks
		//std::for_each(r_row_list_dest.begin(), r_row_list_dest.end(), CColumnCumsumToRowCumsum(r_col_list_B.begin()));
		std::transform(r_col_list_B.begin(), r_col_list_B.end(), r_row_list_dest.begin(), t_ColumnCumsumToRowCumsum);
		//std::for_each(r_col_list_dest.begin(), r_col_list_dest.end(), CColumnCumsumCopy(r_col_list_B.begin())); // copy column layout but not the blocks
		std::transform(r_col_list_B.begin(), r_col_list_B.end(), r_col_list_dest.begin(), t_ColumnCumsumCopy);
		r_dest.m_n_col_num = n_dest_col_num;
		r_dest.m_n_row_num = n_dest_row_num;
		// create layout for the destination matrix (linear time)

		//r_dest.CheckIntegrity();
		// makes sure the dest layout is ok
	}
	// create dest structure

	// no need to reindex / merge, all the indices match and exist

	_TyDenseAllocator alloc(r_dest.m_data_pool);
	// make the allocator object for the templates to use

	/*double f_muinit = m_timer.f_Time();
	printf("_muinit\n");*/
	{
#if 1
		/*double f_loop1 = m_timer.f_Time();
		double f_anal = f_loop1; // ...
		printf("there is 100 %% conflicting columns (assume)\n");
		printf("_loop1\n");*/

		for(_TyColumnConstIter p_col_A_it = col_list_A.begin(), p_col_A_end_it = col_list_A.end();
		   p_col_A_it != p_col_A_end_it; ++ p_col_A_it) {
			const TColumn &r_t_column_A = *p_col_A_it;
			// for each column in A

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PreATA<CBlockMatrixTypelist>::PreATA_UpperTriangle_OuterLoop(
				r_t_column_A, r_row_list_A, r_col_list_dest, alloc);
			// wrap the outer loop in a column size decission tree
		}
		// perform the multiplication of the blocks above diagonal, need to use log(n) lookup
#else
		/*std::vector<CMutexWrap> mutex_list, no_synchronize;
		mutex_list.resize(r_col_list_dest.size() + 1); // one more for allocator lock
		no_synchronize.resize(1);*/

		omp_lock_t t_mutex;
		omp_init_lock(&t_mutex);
		// initialize mutexes

		/*double f_anal = m_timer.f_Time();
		printf("_anal\n");*/

		std::vector<int> refcount_list;
		refcount_list.resize(r_col_list_dest.size(), 0); // initialize to all nulls
		for(_TyColumnConstIter p_col_A_it = col_list_A.begin(), p_col_A_end_it = col_list_A.end();
		   p_col_A_it != p_col_A_end_it; ++ p_col_A_it) {
			const TColumn &r_t_col = *p_col_A_it;
			// for each column in A

			if(r_t_col.block_list.size() < 2)
				continue;
			for(_TyBlockConstIter p_block_it = r_t_col.block_list.begin() + 1,
			   p_end_it = r_t_col.block_list.end(); p_block_it != p_end_it; ++ p_block_it)
				++ refcount_list[(*p_block_it).first];
		}
		// count write-references for columns

		size_t n_conflict_num = 0;
		std::vector<char> conflicting_column_list;
		conflicting_column_list.resize(col_list_A.size(), false);
		{
		_ASSERTE(col_list_A.size() <= INT_MAX);
		int n = int(col_list_A.size());
		#pragma omp parallel for default(shared) reduction(+:n_conflict_num)
		for(int n_column_id_A = 0; n_column_id_A < n; ++ n_column_id_A) {
			const TColumn &r_t_col = col_list_A[n_column_id_A];
			if(r_t_col.block_list.size() < 2)
				continue;
			for(_TyBlockConstIter p_block_it = r_t_col.block_list.begin() + 1,
			   p_end_it = r_t_col.block_list.end(); p_block_it != p_end_it; ++ p_block_it) {
				if(refcount_list[(*p_block_it).first] > 1) {
					conflicting_column_list[n_column_id_A] = true;
					++ n_conflict_num;
					break;
				}
			}
		}}
		// see about columns that write to the same dest-columns (there shouldn't be many)

		printf("there is %lf %% conflicting columns\n", double(n_conflict_num) / col_list_A.size() * 100);
		/*double f_loop1 = m_timer.f_Time();
		printf("_loop1\n");*/

		{
		_ASSERTE(col_list_A.size() <= INT_MAX);
		int n = int(col_list_A.size());
		#pragma omp parallel for default(shared)
		for(int n_column_id_A = 0; n_column_id_A < n; ++ n_column_id_A) {
			const TColumn &r_t_column_A = col_list_A[n_column_id_A];
			// for each column in A

			if(r_t_column_A.block_list.size() < 2 || conflicting_column_list[n_column_id_A])
				continue;
			// only columns with > 1 blocks can produce results

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PreATA<CBlockMatrixTypelist>::PreATA_UpperTriangle_Parallel_OuterLoop(
				r_t_column_A, r_row_list_A, r_col_list_dest, t_mutex, alloc);
			// wrap the outer loop in a column size decission tree
		}}
		// perform the multiplication of the blocks above diagonal, need to use log(n) lookup

		size_t n_column_id_A = 0;
		for(_TyColumnConstIter p_col_A_it = col_list_A.begin(), p_col_A_end_it = col_list_A.end();
		   p_col_A_it != p_col_A_end_it; ++ p_col_A_it, ++ n_column_id_A) {
			const TColumn &r_t_column_A = *p_col_A_it;
			// for each column in A

			if(r_t_column_A.block_list.size() < 2 || !conflicting_column_list[n_column_id_A])
				continue;

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PreATA<CBlockMatrixTypelist>::PreATA_UpperTriangle_OuterLoop(
				r_t_column_A, r_row_list_A, r_col_list_dest, alloc);
			// wrap the outer loop in a column size decission tree
		}
		// perform the multiplication of the blocks above diagonal, need to use log(n) lookup
#endif

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok

		// t_odo test and convert the first loop as well.

		/*double f_loop2 = m_timer.f_Time();
		printf("_loop2\n");*/
		// note that this can run in parallel easily without conflicts, except for the block allocator

		{
		_ASSERTE(r_col_list_B.size() <= INT_MAX);
		int n = int(r_col_list_B.size());

#if 1
		std::vector<double*> block_alloc_list(n);
		for(int n_column_id_B = 0; n_column_id_B < n; ++ n_column_id_B) {
			const TColumn &r_t_column_B = r_col_list_B[n_column_id_B];
			// for each column in B (index is n_column_id_B)

			if(r_t_column_B.block_list.empty())
				continue;
			// only process non-empty columns

			block_alloc_list[n_column_id_B] = r_dest.p_Get_DenseStorage(
				r_t_column_B.n_width * r_t_column_B.n_width);
		}
		// prealloc blocks sequentially
#endif // 1

		#pragma omp parallel for default(shared)
		for(int n_column_id_B = 0; n_column_id_B < n; ++ n_column_id_B) {
			const TColumn &r_t_column_B = r_col_list_B[n_column_id_B];
			// for each column in B (index is n_column_id_B)

			if(r_t_column_B.block_list.empty())
				continue;
			// only process non-empty columns

			double *p_block_data;
			{
#if 1
				p_block_data = block_alloc_list[n_column_id_B];
#else // 1
				omp_set_lock(&t_mutex);
				p_block_data = r_dest.p_Get_DenseStorage(r_t_column_B.n_width * r_t_column_B.n_width);
				omp_unset_lock(&t_mutex); // note this can be moved in front of the loop
#endif // 1
			}
			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_PreATA<CBlockMatrixTypelist>::PreATA_Parallel_Diagonal_OuterLoop(
				r_t_column_B, p_block_data, r_row_list_B);
			// produce the column dot product, each block with itself

			_ASSERTE(r_col_list_dest[n_column_id_B].block_list.empty() ||
				r_col_list_dest[n_column_id_B].block_list.back().first < size_t(n_column_id_B));
			// makes sure that inserting block at diagonal doesn't violate block ordering

			r_col_list_dest[n_column_id_B].block_list.push_back(TColumn::TBlockEntry(n_column_id_B, p_block_data));
			// add the block at the end of the list
		}}
		// generate the blocks at diagonal (const time lookup, fast multiplication)
		// note this is about 30% of work of this function

		/*double f_mudest = m_timer.f_Time();
		printf("_mudestroy\n");
		omp_destroy_lock(&t_mutex);

		double f_end = m_timer.f_Time();
		printf("times:\nmuinit: %lf\nanal: %lf\nloop1: %lf\nloop2: %lf\nmudestroy: %lf\n",
			f_anal - f_muinit, f_loop1 - f_anal, f_loop2 - f_loop1,
			f_mudest - f_loop2, f_end - f_mudest);*/

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok
	}

	if(!b_upper_diagonal_only) {
		for(size_t i = 0, n = r_col_list_dest.size(); i < n; ++ i) {
			const TColumn &r_t_col = r_col_list_dest[i];
			const size_t n_col = i;
			for(_TyBlockConstIter p_block_it = r_t_col.block_list.begin(),
			   p_end_it = r_t_col.block_list.end(); p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				const size_t n_row = r_t_block.first;
				_ASSERTE(n_row <= n_col); // make sure the block is above or at the diagonal
				if(n_row < n_col) {
					// it is guaranteed (by the condition above) that this will not alter block i or any of the
					// blocks to be processed next

					_ASSERTE(r_col_list_dest[n_row].block_list.empty() || r_col_list_dest[n_row].block_list.back().first < n_col);
					// make sure that this won't violate block ordering (since i = n_col is monotonically increasing)

					const size_t n_block_width = r_row_list_dest[n_row].n_height, n_block_height = r_t_col.n_width; // these are dims of the dest block
					double *p_block_data = r_dest.p_Get_DenseStorage(n_block_width * n_block_height);
					r_col_list_dest[n_row].block_list.push_back(TColumn::TBlockEntry(n_col, p_block_data));
					// this block is above the diagonal, at (n_row, n_col), create transpose block at (n_col, n_row)

					const double *p_src_data = r_t_block.second;
					/*for(size_t x = 0; x < n_block_width; ++ x)
						for(size_t y = 0; y < n_block_height; ++ y)
							p_block_data[y + n_block_height * x] = p_src_data[x + n_block_width * y];*/ // f_ixme - does the src / dest ordering matter? is it interchangable? // t_odo - test transposing a matrix // this seems to be correct
					_TyConstMatrixXdRef src(p_src_data, n_block_width, n_block_height);
					_TyMatrixXdRef dest(p_block_data, n_block_height, n_block_width);
					dest.noalias() = src.transpose();
					// copy/transpose the block
				} else {
					_ASSERTE(p_block_it + 1 == p_end_it); // diagonal, gotta be the last block
					break; // still can save the comparison
				}
			}
		}
		// transpose the blocks from the triangle above diagonal to the triangle below diagonal
		// this takes O(n) time in number of blocks above the diagonal, or O(n/2) in number of all blocks
	}
	// no need to optimize this, it is not used in slam, but optimize it later, might transpose blocks using SSE (or at least ommit the two innermost loops)

	//r_dest.CheckIntegrity();
	// makes sure the dest matrix is ok
#else // _OPENMP
	PreMultiplyWithSelfTransposeTo_FBS<CBlockMatrixTypelist>(r_dest, b_upper_diagonal_only);
	// just call serial version (uses less memory than the parallel version)
#endif // _OPENMP
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
void CUberBlockMatrix::InverseOf_BlockDiag_FBS_Parallel(const CUberBlockMatrix &r_A) // throw(std::bad_alloc)
{
	_ASSERTE(r_A.b_SymmetricLayout());
	// inverse only defined for square matrices, this also requires symmetric layout

#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	if(&r_a == this) {
		CUberBlockMatrix tmp; // does not work inplace
		tmp.InverseOf_Symmteric(r_A); // use the regular version
		tmp.Swap(*this);
	} else
		InverseOf_Symmteric(r_A); // use the regular version
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	const size_t n = r_A.m_block_cols_list.size();
	// number of block columns (and rows) in both src and dest matrix

	if(&r_A != this) {
		//Clear(); // CopyLayoutTo() does that
		r_A.CopyLayoutTo(*this);
		// assume the inverse will have the same layout

		bool b_run_in_parallel = (n <= INT_MAX && n > 50);
		// decide whether to run in parallel

		for(size_t i = 0; i < n; ++ i) {
			const TColumn &r_t_col = r_A.m_block_cols_list[i];
			if(r_t_col.block_list.empty())
				continue; // structural rank deficient (allowed)
			_ASSERTE(r_t_col.block_list.size() == 1);
			// contains just a single block (independent from the rest of the matrix)

			TColumn &r_t_dest_col = m_block_cols_list[i];
			_ASSERTE(r_t_dest_col.block_list.empty()); // should be initially empty
			r_t_dest_col.block_list.reserve(1);
			TColumn::TBlockEntry t_block(i,
				p_Get_DenseStorage(r_t_col.n_width * r_t_col.n_width));
			r_t_dest_col.block_list.push_back(t_block);
			// alloc a new (destination) block in this matrix

			if(!b_run_in_parallel) {
				const TColumn::TBlockEntry &r_src_block = r_t_col.block_list.front();
				_ASSERTE(r_src_block.first == i); // block at the diagonal
				const double *p_data = r_src_block.second;
				// make a map of the source block

				blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_BlockwiseUnaryOp<CBlockMatrixTypelist>::template
					BlockwiseUnary_Static_Op_Square<CInvertBlock>(t_block.second, p_data, r_t_col.n_width);
				// calculate inverse of a single block
				// note that in case that there is a rectangular independent (off-diagonal) block,
				// the inverse will fail (inverse is only defined for square matrices)-
			}
		}
		// perform the allocation in series

		if(b_run_in_parallel) {
			_ASSERTE(n <= INT_MAX);
			int _n = int(n);
			#pragma omp parallel for /*if(_n > 50)*/
			for(int i = 0; i < _n; ++ i) {
				const TColumn &r_t_col = r_A.m_block_cols_list[i];
				if(r_t_col.block_list.empty())
					continue; // structural rank deficient (allowed)
				_ASSERTE(r_t_col.block_list.size() == 1);
				// contains just a single block (independent from the rest of the matrix)

				const TColumn::TBlockEntry &r_src_block = r_t_col.block_list.front();
				_ASSERTE(r_src_block.first == i); // block at the diagonal
				const double *p_data = r_src_block.second;
				// make a map of the source block

				TColumn &r_t_dest_col = m_block_cols_list[i];
				_ASSERTE(r_t_dest_col.block_list.size() == 1); // the one we allocated above
				double *p_dest = r_t_dest_col.block_list.front().second;
				// get the destination block in this matrix

				blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_BlockwiseUnaryOp<CBlockMatrixTypelist>::template
					BlockwiseUnary_Static_Op_Square<CInvertBlock>(p_dest, p_data, r_t_col.n_width);
				// calculate inverse of a single block
				// note that in case that there is a rectangular independent (off-diagonal) block,
				// the inverse will fail (inverse is only defined for square matrices)-
			}
		}
		// perform the inverses in parallel unless we chose to do it in series
	} else {
		// working inplace

		_ASSERTE(n <= INT_MAX);
		int _n = int(n);
		#pragma omp parallel for if(_n > 50)
		for(int i = 0; i < _n; ++ i) {
			const TColumn &r_t_col = r_A.m_block_cols_list[i];
			if(r_t_col.block_list.empty())
				continue; // structural rank deficient (allowed)
			_ASSERTE(r_t_col.block_list.size() == 1);
			// contains just a single block (independent from the rest of the matrix)

			const TColumn::TBlockEntry &r_src_block = r_t_col.block_list.front();
			_ASSERTE(r_src_block.first == i); // block at the diagonal
			const double *p_data = r_src_block.second;
			// make a map of the source block

			/*TColumn &r_t_dest_col = m_block_cols_list[i];
			_ASSERTE(r_t_dest_col.block_list.empty()); // should be initially empty
			r_t_dest_col.block_list.reserve(1);
			TColumn::TBlockEntry t_block(i,
				p_Get_DenseStorage(r_t_col.n_width * r_t_col.n_width));
			r_t_dest_col.block_list.push_back(t_block);*/ // working inplace
			// alloc a new (destination) block in this matrix

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_BlockwiseUnaryOp<CBlockMatrixTypelist>::template
				BlockwiseUnary_Static_Op_Square<CInvertBlock>((double*)p_data/*t_block.second*/, // working inplace
				p_data, r_t_col.n_width);
			// calculate inverse of a single block
			// note that in case that there is a rectangular independent (off-diagonal) block,
			// the inverse will fail (inverse is only defined for square matrices)-
		}
	}
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
void CUberBlockMatrix::InverseOf_Symmteric_FBS(const CUberBlockMatrix &r_A, bool b_upper_triangular_source) // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	InverseOf_Symmteric(r_A); // use the regular version
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	// note that we could catch some simple cases where the inverse is not correct,
	// but handling all the cases would be difficult, therefore this function gives
	// absolutely no guarantees about the inverse, and will not fail even when it is
	// evident that the given matrix does not have an inverse. this is intentional.

	r_A.CheckIntegrity(true);

	_ASSERTE(r_A.b_SymmetricLayout()); // actually doesn't need to be symmetric to be invertible
	// makes sure that A is symmetric

	//Clear(); // CopyLayoutTo() does that
	r_A.CopyLayoutTo(*this);
	// assume the inverse will have the same layout

	const size_t n = m_block_cols_list.size();
	// number of block columns (and rows) in both src and dest matrix

	// here, a notion of supernodes is used, where a supernode
	// means a contiguous dense block in the resulting inverse matrix

	std::vector<std::pair<size_t, size_t> > forward_super, backward_super;
	for(size_t i = 0; i < n; ++ i) {
		size_t n_start = i;
		size_t n_end = i;
		// start and end of the supernode

		const TColumn &r_t_col = r_A.m_block_cols_list[i];
		if(!r_t_col.block_list.empty())
			n_end = std::max(n_end, r_t_col.block_list.back().first);
		for(; i < n_end;) {
			++ i;
			const TColumn &r_t_col = r_A.m_block_cols_list[i];
			if(!r_t_col.block_list.empty())
				n_end = std::max(n_end, r_t_col.block_list.back().first);
			// the maximum referenced column
		}
		_ASSERTE(i == n_end);
		// find the maximum referenced column, in all columns between n_start and n_end

		forward_super.push_back(std::make_pair(n_start, n_end));
		// have a new forward supernode
	}
	for(size_t i = n; i > 0;) {
		-- i;
		// here

		size_t n_start = i;
		size_t n_end = i;
		// start and end of the supernode

		const TColumn &r_t_col = r_A.m_block_cols_list[i];
		if(!r_t_col.block_list.empty())
			n_start = std::min(n_start, r_t_col.block_list.front().first);
		for(; i > n_start;) {
			-- i;
			const TColumn &r_t_col = r_A.m_block_cols_list[i];
			if(!r_t_col.block_list.empty())
				n_start = std::min(n_start, r_t_col.block_list.front().first);
			// the minimum referenced column
		}
		_ASSERTE(i == n_start);
		// find the minimum referenced column, in all columns between n_start and n_end

		backward_super.push_back(std::make_pair(n_start, n_end));
		// have a new backward supernode
	}
	// these are always exactly O(n) each, and use up to O(2 * sizeof(size_t) * n) = O(16n) storage each
	// the result is two (sorted) lists of ranges which possibly overlap

	Eigen::MatrixXd dense; // reuse storage; need to expand slices of A to a potentially large dense matrix
	for(size_t i = 0, j = backward_super.size(), n_max_i = forward_super.size(); i < n_max_i;) {
		size_t n_begin = forward_super[i].first;
		size_t n_end = forward_super[i].second;
		// lets begin with a forward-determined supernode

		++ i;
		// skip it, we will process it

		for(;;) {
			//   bb------be
			//               fb------fe
			//	fb > be

			//               bb------be
			//   fb------fe
			//  bb > fe

			bool b_extend = false;

			for(;;) {
				if(j > 0 && !(backward_super[j - 1].first > n_end ||
				   n_begin > backward_super[j - 1].second)) {
					-- j;
					n_begin = std::min(n_begin, backward_super[j].first);
					//_ASSERTE(n_end >= backward_super[j].second); // note that backward probably always extends only the begin // not true
					n_end = std::max(n_end, backward_super[j].second);
					b_extend = true;
					// extend the range
				} else
					break; // range not extended
			}
			// detect overlapping ranges (it is easier to negate disjunct ranges condition)
			// note that the ranges are inclusive

			// maybe it is guaranteed that i < n_max_i, would have to think about it
			for(;;) {
				if(i < n_max_i && !(forward_super[i].first > n_end ||
				   n_begin > forward_super[i].second)) {
					n_begin = std::min(n_begin, forward_super[i].first);
					//_ASSERTE(n_begin <= forward_super[i].first); // note that forward probably always extends only the end // not true
					n_end = std::max(n_end, forward_super[i].second);
					++ i;
					b_extend = true;
					// extend the range
				} else
					break; // range not extended
			}
			// detect overlapping ranges

			if(!b_extend)
				break;
		}
		// merge the sorted ranges in ping-pong fashion (first one extends, then the other extends)

		_ASSERTE((i == n_max_i) == !j);
		// make sure that both of the counter run out at the same time

		// now we know that the supernode is n_start to n_end, we need to grab the blocks from A,
		// put them in a dense matrix, invert it and put it back to this

		if(n_begin == n_end) {
			// a simple case - just a single block gets inverted

			const TColumn &r_t_col = r_A.m_block_cols_list[n_begin];
			if(r_t_col.block_list.empty())
				continue; // structural rank deficient
			_ASSERTE(r_t_col.block_list.size() == 1);
			// contains just a single block (independent from the rest of the matrix)

			const TColumn::TBlockEntry &r_src_block = r_t_col.block_list.front();
			size_t n_row_id = r_src_block.first;
			const double *p_data = r_src_block.second;
			size_t n_row_height = r_A.m_block_cols_list[n_row_id].n_width;
			// make a map of the source block

			_ASSERTE(!b_upper_triangular_source || n_row_id == n_begin);
			// if b_upper_triangular_source is set, make sure that the block is at the diagonal,
			// otherwise this wouldnt work - the inverted block would be bigger by the transposed block

			TColumn &r_t_dest_col = m_block_cols_list[n_begin];
			_ASSERTE(r_t_dest_col.block_list.empty()); // should be initially empty
			r_t_dest_col.block_list.reserve(1);
			TColumn::TBlockEntry t_block(n_row_id,
				p_Get_DenseStorage(n_row_height * r_t_col.n_width));
			r_t_dest_col.block_list.push_back(t_block);
			// alloc a new (destination) block in this matrix

			_ASSERTE(n_row_height == r_t_col.n_width);
			if(n_row_height == r_t_col.n_width) {
				blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_BlockwiseUnaryOp<CBlockMatrixTypelist>::template
					BlockwiseUnary_Static_Op_Square<CInvertBlock>(t_block.second,
					p_data, n_row_height);
			}
			// calculate inverse of a single block
			// note that in case that there is a rectangular independent (off-diagonal) block,
			// the inverse will fail (inverse is only defined for square matrices)
		} else {
			size_t n_first_col = m_block_cols_list[n_begin].n_cumulative_width_sum;
			size_t n_last_col = m_block_cols_list[n_end].n_cumulative_width_sum +
				m_block_cols_list[n_end].n_width; // not inclusive
			// determine elementwise column range

			{
				dense.resize(n_last_col - n_first_col, n_last_col - n_first_col); // no need for conservative
				dense.setZero();
				// allocate and clear a dense matrix

				for(size_t k = n_begin; k <= n_end; ++ k) {
					const TColumn &r_t_col = r_A.m_block_cols_list[k];
					size_t n_dest_column_org = r_t_col.n_cumulative_width_sum - n_first_col; // in elements
					for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
						const TColumn::TBlockEntry &r_block = r_t_col.block_list[j];
						size_t n_row_id = r_block.first;
						const double *p_data = r_block.second;
						// get a block

						_ASSERTE(n_row_id >= n_begin && n_row_id <= n_end);
						// make sure we are inside (otherwise the supernode is not valid and should have been bigger)

						size_t n_row_org = r_A.m_block_cols_list[n_row_id].n_cumulative_width_sum;
						size_t n_row_height = r_A.m_block_cols_list[n_row_id].n_width;
						// is symmetric

						_TyConstMatrixXdRef fill_block(p_data, n_row_height, r_t_col.n_width);
						// make a map

						dense.block(n_row_org - n_first_col, n_dest_column_org,
							n_row_height, r_t_col.n_width) = fill_block;
						// copy the data inside the dense matrix
					}
				}
				// go through all the columns in A and put them in the dense matrix

				if(b_upper_triangular_source)
					dense.triangularView<Eigen::StrictlyLower>() = dense.triangularView<Eigen::StrictlyUpper>().transpose();
				// in case mirroring is needed, do it

				dense = dense.inverse();
				// invert the matrix, making it most likely rather dense
				// todo - now that we are doing this, it would be better to use the template for block combination sizes
				// and calculate the inverse using a fixed size matrix?
			}
			// todo - it would be better to slice the matrix and make a blockwise sparse LU factorization,
			// and use that to calculate the inverse (now we are calculating a rather expensive inverse
			// of a dense (although containing many zeroes) matrix, using dense LU factorization inside Eigen)

			for(size_t k = n_begin; k <= n_end; ++ k) {
				TColumn &r_t_col = m_block_cols_list[k];
				_ASSERTE(r_t_col.block_list.empty()); // should be initially empty
				r_t_col.block_list.resize(n_end + 1 - n_begin);
				// get column, allocate the blocks

				size_t n_dest_column_org = r_t_col.n_cumulative_width_sum - n_first_col; // in elements

				std::vector<TColumn::TBlockEntry>::iterator p_dest_it = r_t_col.block_list.begin();
				for(size_t j = n_begin; j <= n_end; ++ j, ++ p_dest_it) {
					size_t n_row_org = r_A.m_block_cols_list[j].n_cumulative_width_sum;
					size_t n_row_height = m_block_cols_list[j].n_width;
					// is symmetric

					TColumn::TBlockEntry &r_block = *p_dest_it;
					r_block.first = j;
					try {
						r_block.second = p_Get_DenseStorage(n_row_height * r_t_col.n_width);
					} catch(std::bad_alloc &r_exc) {
						r_t_col.block_list.erase(p_dest_it, r_t_col.block_list.end());
						// erase unitialized ones

						throw r_exc; // rethrow
					}
					// fill a block

					_TyMatrixXdRef fill_block(r_block.second, n_row_height, r_t_col.n_width);
					// make a map

					fill_block = dense.block(n_row_org - n_first_col, n_dest_column_org,
						n_row_height, r_t_col.n_width);
					// copy the data back from the dense matrix
				}
			}
			// go through all the columns in the inverse and fill them with inverse data
		}
	}
	// up to O(2n) loops
	// note that this loop could probably be merged with one of the above loops,
	// and only half of the storage would be required (in exchange for slightly messier code)
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
bool CUberBlockMatrix::UpperTriangularTranspose_Solve_FBS(double *p_x, size_t UNUSED(n_vector_size)) const
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return UpperTriangularTranspose_Solve(p_x, n_vector_size);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS

	CheckIntegrity(true);

	_ASSERTE(b_Square());
	// triangular is a special case of square

	_ASSERTE(n_vector_size == m_n_col_num);
	// make sure that the vector's got correct size

	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) { // forward substitution
		const TColumn &r_t_col = *p_col_it;
		if(r_t_col.block_list.empty())
			return false; // a 0 on the diagonal - no solution (is it? csparse would divide by zero and produce a NaN)

		if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_TriangularSolve<CBlockMatrixTypelist>::TriangularSolve_Forward_OuterLoop(r_t_col,
		   p_x, m_block_rows_list))
			return false;
	}

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
bool CUberBlockMatrix::UpperTriangularTranspose_Solve_FBS(double *p_x, size_t UNUSED(n_vector_size),
	const size_t *p_dependent_column, size_t n_dependent_column_num) const
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return UpperTriangularTranspose_Solve(p_x, n_vector_size, p_dependent_column, n_dependent_column_num);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS

	CheckIntegrity(true);

	_ASSERTE(b_Square());
	// triangular is a special case of square

	_ASSERTE(n_vector_size == m_n_col_num);
	// make sure that the vector's got correct size

	_ASSERTE(n_dependent_column_num <= m_block_cols_list.size());

	for(size_t i = 0; i < n_dependent_column_num; ++ i) { // forward substitution
		_ASSERTE(!i || p_dependent_column[i] > p_dependent_column[i - 1]); // must be ordered
		const TColumn &r_t_col = m_block_cols_list[p_dependent_column[i]];

		// the rest of the loop is identical to the full forward-substitution

		if(r_t_col.block_list.empty())
			return false; // a 0 on the diagonal - no solution (is it? csparse would divide by zero and produce a NaN)

		if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_TriangularSolve<CBlockMatrixTypelist>::TriangularSolve_Forward_OuterLoop(r_t_col,
		   p_x, m_block_rows_list))
			return false;
	}

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
bool CUberBlockMatrix::UpperTriangularTranspose_Solve_FBS(double *p_x, size_t UNUSED(n_vector_size), size_t n_skip_columns) const
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return UpperTriangularTranspose_Solve(p_x, n_vector_size, n_skip_columns);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

	//Check_Block_Alignment();
	_ASSERTE(b_Square());
	// triangular is a special case of square

	_ASSERTE(n_vector_size == m_n_col_num);
	// make sure that the vector's got correct size

	_ASSERTE(n_skip_columns <= m_block_cols_list.size());
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin() + n_skip_columns,
	   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) { // forward substitution
		const TColumn &r_t_col = *p_col_it;
		if(r_t_col.block_list.empty())
			return false; // a 0 on the diagonal - no solution (is it? csparse would divide by zero and produce a NaN)

		if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_TriangularSolve<CBlockMatrixTypelist>::TriangularSolve_Forward_OuterLoop(r_t_col,
		   p_x, m_block_rows_list))
			return false;
	}

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
bool CUberBlockMatrix::UpperTriangular_Solve_FBS(double *p_x,
	size_t UNUSED(n_vector_size), size_t n_last_column, size_t n_first_column) const
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return UpperTriangular_Solve(p_x, n_vector_size);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

	//Check_Block_Alignment();
	_ASSERTE(b_Square());
	// triangular is a special case of square

	//_ASSERTE(b_SymmetricLayout()); // note that symmetric layout is not really required, it just must hold that the last block of every column is the diagonal block
	// make sure that the layout is symmetric (this is an optimization, not a prerequisite given by math)

	_ASSERTE(n_vector_size == m_n_col_num);
	// make sure that the vector's got correct size

	_ASSERTE(n_last_column <= n_first_column);
	// it goes backwards so last is before first

	_ASSERTE(n_first_column < m_block_cols_list.size());
	// should be a valid block index

	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin() + (n_first_column + 1),
	   p_end_it = m_block_cols_list.begin() + n_last_column; p_col_it != p_end_it;) { // back substitution
		-- p_col_it;
		// decrement at the beginning of the loop! (but after the comparison)

		const TColumn &r_t_col = *p_col_it;
		if(r_t_col.block_list.empty())
			return false; // a 0 on the diagonal - no solution (is it? csparse would divide by zero and produce a NaN)
		//_ASSERTE(r_t_col.block_list.back().first == n_col); // makes sure that the matrix really is upper diagonal

		if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_TriangularSolve<CBlockMatrixTypelist>::TriangularSolve_Back_OuterLoop(r_t_col,
		   p_x, m_block_rows_list))
			return false;
	}

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CMatrixBlockSizeList, const int n_max_matrix_size>
bool CUberBlockMatrix::Cholesky_Dense_FBS() // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return Cholesky_Dense<Eigen::Dynamic>(); // use dynamic-sized matrices (slow)
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	CheckIntegrity(true);

	_ASSERTE(b_SymmetricLayout());
	// makes sure it has symmetric layout

	if(m_n_col_num <= n_max_matrix_size) {
		bool b_result;
		_ASSERTE(m_n_col_num <= INT_MAX);
		fbs_ut::CWrap3<>::In_SquareMatrixSize_DecisionTree<CMatrixBlockSizeList, n_max_matrix_size>(int(m_n_col_num),
			blockmatrix_detail::CUberBlockMatrix_FBS::CCallDenseCholesky(*this, b_result));
		return b_result;
	} else
		return Cholesky_Dense<Eigen::Dynamic>(); // use dynamic-sized matrices (slow)
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
bool CUberBlockMatrix::CholeskyOf_FBS(const CUberBlockMatrix &r_lambda)
{
	//r_lambda.Check_Block_Alignment();
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return CholeskyOf(r_lambda);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	const size_t n = r_lambda.m_block_cols_list.size();
	// get number of block columns

	std::vector<size_t> etree, ereach_stack, bitfield;
	bitfield.resize(n, 0);

	r_lambda.Build_EliminationTree(etree, ereach_stack); // use ereach stack as workspace
	_ASSERTE(ereach_stack.size() == n);
	// build an elimination tree

	return CholeskyOf_FBS<CBlockMatrixTypelist>(r_lambda, etree, ereach_stack, bitfield);
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
bool CUberBlockMatrix::CholeskyOf_FBS(const CUberBlockMatrix &r_lambda, size_t n_start_on_column)
{
	//r_lambda.Check_Block_Alignment();
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return CholeskyOf(r_lambda, n_start_on_column);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	const size_t n = r_lambda.m_block_cols_list.size();
	// get number of block columns

	std::vector<size_t> etree, ereach_stack, bitfield;
	bitfield.resize(n, 0);

	r_lambda.Build_EliminationTree(etree, ereach_stack); // use ereach stack as workspace
	_ASSERTE(ereach_stack.size() == n);
	// build an elimination tree

	return CholeskyOf_FBS<CBlockMatrixTypelist>(r_lambda,
		etree, ereach_stack, bitfield, n_start_on_column);
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
bool CUberBlockMatrix::CholeskyOf_FBS(const CUberBlockMatrix &r_lambda, const std::vector<size_t> &r_elim_tree,
	std::vector<size_t> &r_workspace, std::vector<size_t> &r_zero_workspace) // throw(std::bad_alloc)
{
	//r_lambda.Check_Block_Alignment();
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return CholeskyOf(r_lambda, r_elim_tree, r_workspace, r_zero_workspace);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	Clear();

	_ASSERTE(r_lambda.b_SymmetricLayout());
	// makes sure that lambda is symmetric

	_ASSERTE(r_lambda.m_n_row_num == r_lambda.m_n_col_num);
	m_n_row_num = m_n_col_num = r_lambda.m_n_col_num;
	m_block_rows_list = r_lambda.m_block_rows_list;
	m_block_cols_list.resize(r_lambda.m_block_cols_list.size());
	// copy the layout

	const size_t n = m_block_cols_list.size();
	// get number of block columns

	std::vector<size_t> &ereach_stack = r_workspace, &bitfield = r_zero_workspace; // alias
	_ASSERTE(ereach_stack.size() >= n);
	_ASSERTE(bitfield.size() >= n);
	// alloc all the whatnots ...

	_TyDenseAllocator alloc(m_data_pool);

	for(size_t j = 0; j < n; ++ j) { // for every column (no sparsity here, L should be full-rank)
#if 1
		const TColumn &r_col_A_j = r_lambda.m_block_cols_list[j];
		if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_Cholesky<CBlockMatrixTypelist>::b_Column_Loop(r_col_A_j.n_width, j, n,
		   m_block_cols_list, r_col_A_j, r_lambda, r_elim_tree, ereach_stack, bitfield, alloc)) {
			//printf("error: not pos def\n"); // not pos def
			Clear(); // otherwise leaving uninit columns behind, CheckIntegrity() will yell
			return false;
		}
#else // 1
		TColumn &r_col_L_j = m_block_cols_list[j];
		const TColumn &r_col_A_j = r_lambda.m_block_cols_list[j];
		const size_t n_col_j_width = r_col_L_j.n_width = r_col_A_j.n_width;
		r_col_L_j.n_cumulative_width_sum = r_col_A_j.n_cumulative_width_sum;
		// get columns

		_ASSERTE(!r_col_A_j.block_list.empty()); // otherwise rank deficient (make it a runtime check?)
		_TyBlockConstIter p_A_block_it =
			r_col_A_j.block_list.begin(), p_A_block_end_it = r_col_A_j.block_list.end();
		// get iterator to blocks of the original matrix

#if 0
		size_t n_ereach_size = r_lambda.n_Build_EReach(j, r_elim_tree, ereach_stack, bitfield);//cs_ereach(p_block_struct, j, p_etree, &s[0], &w[0]);
#else
		size_t n_ereach_first = r_lambda.n_Build_EReach(j, r_elim_tree, ereach_stack, bitfield);//cs_ereach(p_block_struct, j, p_etree, &s[0], &w[0]);
#endif
		// use ereach to compute nonzero pattern

		r_col_L_j.block_list.reserve(n - n_ereach_first/*n_ereach_size*/ + 1); // + 1 amounts for the diagonal
		// reserve space for blocks to avoid reallocation later

		for(size_t u = n_ereach_first, n_highest_k = 0; u < n; ++ u) { // seems to work rather nicely (strange because not every A_up[k, j] is not accessed then - it is likely null) // todo - verify this
			const size_t k = ereach_stack[u]; // use ereach to predict which columns will have nonzero products
			// k generally rises, but it doesn't have to (can be non-monotonic)

			_ASSERTE(k != n_highest_k || u == n_ereach_first); // there should be no column repeated (except for zero, in the first iteration)
			bool b_ereach_mono;
			if((b_ereach_mono = (k >= n_highest_k)))
				n_highest_k = k; // don't remember previous k unless it increased
			// see if the ereach is (locally) increasing

			const TColumn &r_col_L_k = m_block_cols_list[k];
			const size_t n_col_k_width = r_col_L_k.n_width;
			// get column k

			_TyBlockConstIter p_jk_block_it;
			double *p_k_block_data = p_Get_DenseStorage(n_col_k_width * n_col_j_width);
			//_TyMatrixXdRef L_block_kj(p_k_block_data, n_col_k_width, n_col_j_width); // only accessed in FBS section
			if(b_ereach_mono) { // most of the time
				_ASSERTE(r_col_L_j.block_list.empty() || r_col_L_j.block_list.back().first < k); // this is generally not true - ereach doesn't have to be monotonic, might need to insert the block in a middle
				// makes sure that the column doesn't contain the k-th block yet

				r_col_L_j.block_list.push_back(TColumn::TBlockEntry(k, p_k_block_data));
				p_jk_block_it = r_col_L_j.block_list.end() - 1; // it is here
			} else {
				//printf("ereach not mono\n"); // debug
				_ASSERTE(!r_col_L_j.block_list.empty() && r_col_L_j.block_list.back().first > k); // r_col_L_j.block_list.back().first = n_highest_k and n_highest_k > k
				// make sure we're not going to search for the correct position in vain

				_TyBlockIter p_dest_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
					std::lower_bound(r_col_L_j.block_list.begin(), r_col_L_j.block_list.end(), k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(r_col_L_j.block_list.begin(), r_col_L_j.block_list.end(), k, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				// have to search for the insertion position to keep the column sorted

				p_jk_block_it = r_col_L_j.block_list.insert(p_dest_block_it,
					TColumn::TBlockEntry(k, p_k_block_data));
				// insert and remember where it is
			}
			// add a new off-diagonal block to L

			{
				_ASSERTE(p_A_block_it != p_A_block_end_it);
				size_t n_block_row_id;
				if((n_block_row_id = (*p_A_block_it).first) < k) {
					p_A_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						std::lower_bound(p_A_block_it, p_A_block_end_it, k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(p_A_block_it, p_A_block_end_it, k, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				} else if(n_block_row_id > k) {
					p_A_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						std::lower_bound(r_col_A_j.block_list.begin(), p_A_block_end_it, k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(r_col_A_j.block_list.begin(), p_A_block_end_it, k, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					// a rare case where ereach is not monotonically increasing
				}
				_ASSERTE(p_A_block_it != p_A_block_end_it);
			}
			// look up the block in the source matrix

			_ASSERTE(!r_col_L_j.block_list.empty() && (!b_ereach_mono || r_col_L_j.block_list.back().first == k)); // this is generally not true; it might not be the last
			_ASSERTE((*p_jk_block_it).first == k); // this should be always true
			// column j now contains the k-th block

			_TyBlockConstIter
				p_j_block_it = r_col_L_j.block_list.begin(),
				p_j_block_end_it = p_jk_block_it/*r_col_L_j.block_list.end() - 1*/, // this is not end() - 1, might have to look for the block before k
				p_k_block_it = r_col_L_k.block_list.begin(),
				p_k_block_end_it = r_col_L_k.block_list.end();
			// have to loop through both lists and merge the blocks to find the ones, referencing the same rows

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_Cholesky<CBlockMatrixTypelist>::OffDiagonal_Loop(p_k_block_data,
				n_col_j_width, n_col_k_width, p_j_block_it, p_j_block_end_it,
				p_k_block_it, p_k_block_end_it, *p_A_block_it, k, m_block_cols_list);
			// execute cmod and solve by diagonals using FBS

			if((*p_A_block_it).first == k)
				++ p_A_block_it;
			// skip to the next one for the next iteration / for the diagonal
		}
		// complex data dependencies in here, not sure if i like it

		{
			_ASSERTE(p_A_block_it != p_A_block_end_it && (*p_A_block_it).first <= j);
			// it is pointing before or at the diagonal already
			if((*p_A_block_it).first != j) {
				p_A_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
					std::lower_bound(p_A_block_it, p_A_block_end_it, j, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(p_A_block_it, p_A_block_end_it, j, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			}
			_ASSERTE(p_A_block_it != p_A_block_end_it && (*p_A_block_it).first == j); // should always be nonzero
			// find the diagonal block (note that if A is upper diagonal, it is the
			// last block and we can optimize for that)

			_ASSERTE(r_col_L_j.block_list.empty() || r_col_L_j.block_list.back().first < j);
			// makes sure that the column doesn't contain the diagonal block yet

			double *p_diag_data = p_Get_DenseStorage(n_col_j_width * n_col_j_width);
			r_col_L_j.block_list.push_back(TColumn::TBlockEntry(j, p_diag_data));
			// allocates a new block in this column

			if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_Cholesky<CBlockMatrixTypelist>::b_Diagonal_Loop(p_diag_data,
			   n_col_j_width, r_col_L_j.block_list.begin(), r_col_L_j.block_list.end() - 1,
			   (*p_A_block_it).second, m_block_cols_list)) {
				//printf("error: not pos def\n"); // not pos def
				Clear(); // otherwise leaving uninit columns behind, CheckIntegrity() will yell
				return false;
			}
			// execute cdiv using FBS
		}
		// cdiv; reads an entire column and produces diagonal

		_ASSERTE(r_col_L_j.block_list.size() == n - n_ereach_first/*n_ereach_size*/ + 1);
		// make sure we preallocated it correclty
#endif // 1
	}
	// note that Pigglet could probably write it better

	//CheckIntegrity();
	// makes sure the factor matrix (this) is ok

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

template <class CBlockMatrixTypelist>
bool CUberBlockMatrix::CholeskyOf_FBS(const CUberBlockMatrix &r_lambda, const std::vector<size_t> &r_elim_tree,
	std::vector<size_t> &r_workspace, std::vector<size_t> &r_zero_workspace,
	size_t n_start_on_column) // throw(std::bad_alloc)
{
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	return CholeskyOf(r_lambda, r_elim_tree, r_workspace, r_zero_workspace, n_start_on_column);
#else // __UBER_BLOCK_MATRIX_SUPRESS_FBS
	//Check_Block_Alignment(); // resumed, check itself as well ...
	//r_lambda.Check_Block_Alignment();
	//Clear();

	_ASSERTE(r_lambda.b_SymmetricLayout());
	// makes sure that lambda is symmetric

	_ASSERTE(r_lambda.m_n_row_num == r_lambda.m_n_col_num);
	m_n_row_num = m_n_col_num = r_lambda.m_n_col_num;
	m_block_rows_list = r_lambda.m_block_rows_list;
	m_block_cols_list.resize(r_lambda.m_block_cols_list.size()); // todo - make sure that prefix of the layout till n_start_on_column is the same
	// copy the layout

	const size_t n = m_block_cols_list.size();
	// get number of block columns

	std::vector<size_t> &ereach_stack = r_workspace, &bitfield = r_zero_workspace; // alias
	_ASSERTE(ereach_stack.size() >= n);
	_ASSERTE(bitfield.size() >= n);
	// alloc all the whatnots ...

	_ASSERTE(n_start_on_column <= n);

#ifdef _DEBUG
	for(size_t j = 0; j < n_start_on_column; ++ j) {
		TColumn &r_col_L_j = m_block_cols_list[j];
		const TColumn &r_col_A_j = r_lambda.m_block_cols_list[j];
		_ASSERTE(r_col_L_j.n_width == r_col_A_j.n_width);
		_ASSERTE(r_col_L_j.n_cumulative_width_sum == r_col_A_j.n_cumulative_width_sum);
		// just copy column layouts
	}
#endif // _DEBUG
	// makes sure that prefix of the layout till n_start_on_column is the same

	_TyDenseAllocator alloc(m_data_pool);

	for(size_t j = n_start_on_column; j < n; ++ j) { // for every column (no sparsity here, L should be full-rank)
#if 1
		const TColumn &r_col_A_j = r_lambda.m_block_cols_list[j];
		if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_Cholesky<CBlockMatrixTypelist>::b_Column_Loop(r_col_A_j.n_width, j, n,
		   m_block_cols_list, r_col_A_j, r_lambda, r_elim_tree, ereach_stack, bitfield, alloc)) {
			//printf("error: not pos def\n"); // not pos def
			Clear(); // otherwise leaving uninit columns behind, CheckIntegrity() will yell
			return false;
		}
#else // 1
		TColumn &r_col_L_j = m_block_cols_list[j];
		const TColumn &r_col_A_j = r_lambda.m_block_cols_list[j];
		const size_t n_col_j_width = r_col_L_j.n_width = r_col_A_j.n_width;
		r_col_L_j.n_cumulative_width_sum = r_col_A_j.n_cumulative_width_sum;
		// get columns

		_ASSERTE(!r_col_A_j.block_list.empty()); // otherwise rank deficient (make it a runtime check?)
		_TyBlockConstIter p_A_block_it =
			r_col_A_j.block_list.begin(), p_A_block_end_it = r_col_A_j.block_list.end();
		// get iterator to blocks of the original matrix

#if 0
		size_t n_ereach_size = r_lambda.n_Build_EReach(j, r_elim_tree, ereach_stack, bitfield);//cs_ereach(p_block_struct, j, p_etree, &s[0], &w[0]);
#else
		size_t n_ereach_first = r_lambda.n_Build_EReach(j, r_elim_tree, ereach_stack, bitfield);//cs_ereach(p_block_struct, j, p_etree, &s[0], &w[0]);
#endif
		// use ereach to compute nonzero pattern

		r_col_L_j.block_list.reserve(n - n_ereach_first/*n_ereach_size*/ + 1); // + 1 amounts for the diagonal
		// reserve space for blocks to avoid reallocation later

		for(size_t u = n_ereach_first, n_highest_k = 0; u < n; ++ u) { // seems to work rather nicely (strange because not every A_up[k, j] is not accessed then - it is likely null) // todo - verify this
			const size_t k = ereach_stack[u]; // use ereach to predict which columns will have nonzero products
			// k generally rises, but it doesn't have to (can be non-monotonic)

			_ASSERTE(k != n_highest_k || u == n_ereach_first); // there should be no column repeated (except for zero, in the first iteration)
			bool b_ereach_mono;
			if((b_ereach_mono = (k >= n_highest_k)))
				n_highest_k = k; // don't remember previous k unless it increased
			// see if the ereach is (locally) increasing

			const TColumn &r_col_L_k = m_block_cols_list[k];
			const size_t n_col_k_width = r_col_L_k.n_width;
			// get column k

			_TyBlockConstIter p_jk_block_it;
			double *p_k_block_data = p_Get_DenseStorage(n_col_k_width * n_col_j_width);
			//_TyMatrixXdRef L_block_kj(p_k_block_data, n_col_k_width, n_col_j_width); // only accessed in FBS section
			if(b_ereach_mono) { // most of the time
				_ASSERTE(r_col_L_j.block_list.empty() || r_col_L_j.block_list.back().first < k); // this is generally not true - ereach doesn't have to be monotonic, might need to insert the block in a middle
				// makes sure that the column doesn't contain the k-th block yet

				r_col_L_j.block_list.push_back(TColumn::TBlockEntry(k, p_k_block_data));
				p_jk_block_it = r_col_L_j.block_list.end() - 1; // it is here
			} else {
				//printf("ereach not mono\n"); // debug
				_ASSERTE(!r_col_L_j.block_list.empty() && r_col_L_j.block_list.back().first > k); // r_col_L_j.block_list.back().first = n_highest_k and n_highest_k > k
				// make sure we're not going to search for the correct position in vain

				_TyBlockIter p_dest_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
					std::lower_bound(r_col_L_j.block_list.begin(), r_col_L_j.block_list.end(), k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(r_col_L_j.block_list.begin(), r_col_L_j.block_list.end(), k, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				// have to search for the insertion position to keep the column sorted

				p_jk_block_it = r_col_L_j.block_list.insert(p_dest_block_it,
					TColumn::TBlockEntry(k, p_k_block_data));
				// insert and remember where it is
			}
			// add a new off-diagonal block to L

			{
				_ASSERTE(p_A_block_it != p_A_block_end_it);
				size_t n_block_row_id;
				if((n_block_row_id = (*p_A_block_it).first) < k) {
					p_A_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						std::lower_bound(p_A_block_it, p_A_block_end_it, k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(p_A_block_it, p_A_block_end_it, k, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				} else if(n_block_row_id > k) {
					p_A_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						std::lower_bound(r_col_A_j.block_list.begin(), p_A_block_end_it, k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(r_col_A_j.block_list.begin(), p_A_block_end_it, k, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					// a rare case where ereach is not monotonically increasing
				}
				_ASSERTE(p_A_block_it != p_A_block_end_it);
			}
			// look up the block in the source matrix

			_ASSERTE(!r_col_L_j.block_list.empty() && (!b_ereach_mono || r_col_L_j.block_list.back().first == k)); // this is generally not true; it might not be the last
			_ASSERTE((*p_jk_block_it).first == k); // this should be always true
			// column j now contains the k-th block

			_TyBlockConstIter
				p_j_block_it = r_col_L_j.block_list.begin(),
				p_j_block_end_it = p_jk_block_it/*r_col_L_j.block_list.end() - 1*/, // this is not end() - 1, might have to look for the block before k
				p_k_block_it = r_col_L_k.block_list.begin(),
				p_k_block_end_it = r_col_L_k.block_list.end();
			// have to loop through both lists and merge the blocks to find the ones, referencing the same rows

			blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_Cholesky<CBlockMatrixTypelist>::OffDiagonal_Loop(p_k_block_data,
				n_col_j_width, n_col_k_width, p_j_block_it, p_j_block_end_it,
				p_k_block_it, p_k_block_end_it, *p_A_block_it, k, m_block_cols_list);
			// execute cmod and solve by diagonals using FBS

			if((*p_A_block_it).first == k)
				++ p_A_block_it;
			// skip to the next one for the next iteration / for the diagonal
		}
		// complex data dependencies in here, not sure if i like it

		{
			_ASSERTE(p_A_block_it != p_A_block_end_it && (*p_A_block_it).first <= j);
			// it is pointing before or at the diagonal already
			if((*p_A_block_it).first != j) {
				p_A_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
					std::lower_bound(p_A_block_it, p_A_block_end_it, j, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(p_A_block_it, p_A_block_end_it, j, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			}
			_ASSERTE(p_A_block_it != p_A_block_end_it && (*p_A_block_it).first == j); // should always be nonzero
			// find the diagonal block (note that if A is upper diagonal, it is the
			// last block and we can optimize for that)

			_ASSERTE(r_col_L_j.block_list.empty() || r_col_L_j.block_list.back().first < j);
			// makes sure that the column doesn't contain the diagonal block yet

			double *p_diag_data = p_Get_DenseStorage(n_col_j_width * n_col_j_width);
			r_col_L_j.block_list.push_back(TColumn::TBlockEntry(j, p_diag_data));
			// allocates a new block in this column

			if(!blockmatrix_detail::CUberBlockMatrix_FBS::CFBS_Cholesky<CBlockMatrixTypelist>::b_Diagonal_Loop(p_diag_data,
			   n_col_j_width, r_col_L_j.block_list.begin(), r_col_L_j.block_list.end() - 1,
			   (*p_A_block_it).second, m_block_cols_list)) {
				//printf("error: not pos def\n"); // not pos def
				return false;
			}
			// execute cdiv using FBS
		}
		// cdiv; reads an entire column and produces diagonal

		_ASSERTE(r_col_L_j.block_list.size() == n - n_ereach_first/*n_ereach_size*/ + 1);
		// make sure we preallocated it correclty
#endif // 1
	}
	// note that Pigglet could probably write it better

	//CheckIntegrity();
	// makes sure the factor matrix (this) is ok

	return true;
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
}

#endif // !__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_FUNCTIONS_IMPLEMENTATION_INCLUDED
