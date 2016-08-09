/*
								+---------------------------------+
								|                                 |
								|  ***   Über Block Matrix   ***  |
								|                                 |
								| Copyright  (c) -tHE SWINe- 2013 |
								|                                 |
								|         BlockMatrix.inl         |
								|                                 |
								+---------------------------------+
*/

#pragma once
#ifndef __UBER_BLOCK_MATRIX_INLINES_INCLUDED
#define __UBER_BLOCK_MATRIX_INLINES_INCLUDED

/**
 *	@file include/slam/BlockMatrix.inl
 *	@date 2012
 *	@author -tHE SWINe-
 *	@brief the überblockmatrix inline and template function definitions
 */

#include "slam/BlockMatrix.h"

inline size_t CUberBlockMatrix::n_Find_BlockColumn(size_t n_column, size_t &r_n_block_column_num) const
{
	size_t n_column_index;
	if(n_column >= m_n_col_num)
		return -1;
	else {
		_TyColumnConstIter p_col_it;
		_ASSERTE(!m_block_cols_list.empty());
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
		p_col_it = std::upper_bound(m_block_cols_list.begin(),
			m_block_cols_list.end(), n_column, CFindLowerColumn()); // t_odo
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		p_col_it = std::upper_bound(m_block_cols_list.begin(),
			m_block_cols_list.end(), n_column, _FindLowerColumn); // t_odo
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		//_ASSERTE(p_col_it != m_block_cols_list.end()); // t_odo - handle this as proper runtime error // done (the line below)
		if(p_col_it == m_block_cols_list.begin())
			return -1;
		-- p_col_it;
		const TColumn &r_t_col = *p_col_it;

		_ASSERTE(n_column >= r_t_col.n_cumulative_width_sum &&
			n_column < r_t_col.n_cumulative_width_sum + r_t_col.n_width);
		// since n_column is nonnegative, it is inside of the matrix (the condition above),
		// and block columns span the entire matrix area, it can't be not found

		r_n_block_column_num = r_t_col.n_width;
		n_column_index = p_col_it - m_block_cols_list.begin();
		// just find the index, the column in question does exist, and has the right size (ideal case)

		_ASSERTE(n_column_index < m_block_cols_list.size());
		_ASSERTE(m_block_cols_list[n_column_index].n_width == r_n_block_column_num);
		_ASSERTE(n_column >= m_block_cols_list[n_column_index].n_cumulative_width_sum &&
			n_column < m_block_cols_list[n_column_index].n_cumulative_width_sum +
			m_block_cols_list[n_column_index].n_width);
		// make sure that the column reference is really resolved

		return n_column_index;
	}
	// resolve column reference (essentially the same code as for the row reference)
}

inline CUberBlockMatrix::_TyConstMatrixXdRef CUberBlockMatrix::t_Block_AtColumn(size_t n_column_index, size_t n_block_index) const
{
	_ASSERTE(n_column_index < m_block_cols_list.size()); // make sure n_column_index points to a valid column
	_ASSERTE(n_block_index < m_block_cols_list[n_column_index].block_list.size()); // make sure n_block_index selects a block in this column
	const TColumn::TBlockEntry &r_block = m_block_cols_list[n_column_index].block_list[n_block_index];
	size_t n_block_row_num = m_block_rows_list[r_block.first].n_height,
		n_block_column_num = m_block_cols_list[n_column_index].n_width;
	return _TyConstMatrixXdRef(r_block.second, n_block_row_num, n_block_column_num);
}

inline CUberBlockMatrix::_TyMatrixXdRef CUberBlockMatrix::t_Block_AtColumn(size_t n_column_index, size_t n_block_index)
{
	_ASSERTE(n_column_index < m_block_cols_list.size()); // make sure n_column_index points to a valid column
	_ASSERTE(n_block_index < m_block_cols_list[n_column_index].block_list.size()); // make sure n_block_index selects a block in this column
	TColumn::TBlockEntry &r_block = m_block_cols_list[n_column_index].block_list[n_block_index];
	size_t n_block_row_num = m_block_rows_list[r_block.first].n_height,
		n_block_column_num = m_block_cols_list[n_column_index].n_width;
	return _TyMatrixXdRef(r_block.second, n_block_row_num, n_block_column_num);
}

template <const int n_block_row_num, const int n_block_column_num>
inline typename CUberBlockMatrix::CMakeMatrixRef<n_block_row_num, n_block_column_num>::_TyConst
	CUberBlockMatrix::t_Block_AtColumn(size_t n_column_index, size_t n_block_index) const
{
	_ASSERTE(n_column_index < m_block_cols_list.size()); // make sure n_column_index points to a valid column
	_ASSERTE(n_block_index < m_block_cols_list[n_column_index].block_list.size()); // make sure n_block_index selects a block in this column
	const TColumn::TBlockEntry &r_block = m_block_cols_list[n_column_index].block_list[n_block_index];
	_ASSERTE(n_block_row_num == m_block_rows_list[r_block.first].n_height &&
		n_block_column_num == m_block_cols_list[n_column_index].n_width);
	return typename CMakeMatrixRef<n_block_row_num, n_block_column_num>::_TyConst(r_block.second);
}

template <const int n_block_row_num, const int n_block_column_num>
inline typename CUberBlockMatrix::CMakeMatrixRef<n_block_row_num, n_block_column_num>::_Ty
	CUberBlockMatrix::t_Block_AtColumn(size_t n_column_index, size_t n_block_index)
{
	_ASSERTE(n_column_index < m_block_cols_list.size()); // make sure n_column_index points to a valid column
	_ASSERTE(n_block_index < m_block_cols_list[n_column_index].block_list.size()); // make sure n_block_index selects a block in this column
	const TColumn::TBlockEntry &r_block = m_block_cols_list[n_column_index].block_list[n_block_index];
	_ASSERTE(n_block_row_num == m_block_rows_list[r_block.first].n_height &&
		n_block_column_num == m_block_cols_list[n_column_index].n_width);
	return typename CMakeMatrixRef<n_block_row_num, n_block_column_num>::_Ty(r_block.second);
}

inline size_t CUberBlockMatrix::n_Row_Num() const
{
	_ASSERTE((!m_n_row_num && m_block_rows_list.empty()) ||
		(m_n_row_num && m_block_rows_list.back().n_height +
		m_block_rows_list.back().n_cumulative_height_sum == m_n_row_num));
	return m_n_row_num;
}

inline size_t CUberBlockMatrix::n_Column_Num() const
{
	_ASSERTE((!m_n_col_num && m_block_cols_list.empty()) ||
		(m_n_col_num && m_block_cols_list.back().n_width +
		m_block_cols_list.back().n_cumulative_width_sum == m_n_col_num));
	return m_n_col_num;
}

inline size_t CUberBlockMatrix::n_Block_Num() const
{
	size_t n_num = 0;
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i)
		n_num += m_block_cols_list[i].block_list.size();
	return n_num;
}

inline bool CUberBlockMatrix::b_Empty() const
{
	bool b_empty = !m_n_row_num && !m_n_col_num;
	_ASSERTE(!b_empty || m_block_rows_list.empty()); // if empty, the layouts should also be empty
	_ASSERTE(!b_empty || m_block_cols_list.empty()); // if empty, the layouts should also be empty
	_ASSERTE(!b_empty || m_data_pool.empty()); // if empty, the pool should also be empty
	_ASSERTE(!b_empty || !m_n_ref_elem_num); // if empty, m_n_ref_elem_num should also be null
	return b_empty;
}

inline CUberBlockMatrix::_TyMatrixXdRef CUberBlockMatrix::t_FindBlock(size_t n_row, size_t n_column,
	size_t n_block_row_num, size_t n_block_column_num, bool b_alloc_if_not_found /*= true*/,
	bool b_mind_uninitialized /*= true*/) // throw(std::bad_alloc)
{
	double *p_data = p_FindBlock(n_row, n_column, n_block_row_num,
		n_block_column_num, b_alloc_if_not_found, b_mind_uninitialized);
	if(p_data)
		return _TyMatrixXdRef(p_data, n_block_row_num, n_block_column_num);
	else
		return _TyMatrixXdRef(0, 0, 0); // fixme - does this throw an exception or what?
}

inline CUberBlockMatrix::_TyMatrixXdRef CUberBlockMatrix::t_FindBlock(size_t n_row, size_t n_column)
{
	size_t n_block_row_num, n_block_column_num;
	double *p_data = p_FindBlock_ResolveSize(n_row, n_column, n_block_row_num, n_block_column_num);
	if(p_data)
		return _TyMatrixXdRef(p_data, n_block_row_num, n_block_column_num);
	else
		return _TyMatrixXdRef(0, 0, 0); // fixme - does this throw an exception or what?
}

inline CUberBlockMatrix::_TyConstMatrixXdRef CUberBlockMatrix::t_FindBlock(size_t n_row, size_t n_column) const
{
	size_t n_block_row_num, n_block_column_num;
	const double *p_data = p_FindBlock_ResolveSize(n_row, n_column, n_block_row_num, n_block_column_num);
	if(p_data)
		return _TyConstMatrixXdRef(p_data, n_block_row_num, n_block_column_num);
	else
		return _TyConstMatrixXdRef(0, 0, 0); // fixme - does this throw an exception or what?
}

inline CUberBlockMatrix::_TyConstMatrixXdRef CUberBlockMatrix::t_GetBlock_Log(size_t n_row_index, size_t n_column_index) const
{
	_ASSERTE(n_row_index < m_block_rows_list.size());
	_ASSERTE(n_column_index < m_block_cols_list.size());
	size_t n_block_row_num = m_block_rows_list[n_row_index].n_height;
	size_t n_block_column_num = m_block_cols_list[n_column_index].n_width;
	const double *p_data = p_GetBlockData(n_row_index, n_column_index);
	if(p_data)
		return _TyConstMatrixXdRef(p_data, n_block_row_num, n_block_column_num);
	else
		return _TyConstMatrixXdRef(0, 0, 0); // fixme - does this throw an exception or what?
}

inline CUberBlockMatrix::_TyMatrixXdRef CUberBlockMatrix::t_GetBlock_Log(size_t n_row_index, size_t n_column_index)
{
	_ASSERTE(n_row_index < m_block_rows_list.size());
	_ASSERTE(n_column_index < m_block_cols_list.size());
	size_t n_block_row_num = m_block_rows_list[n_row_index].n_height;
	size_t n_block_column_num = m_block_cols_list[n_column_index].n_width;
	double *p_data = p_GetBlockData(n_row_index, n_column_index);
	if(p_data)
		return _TyMatrixXdRef(p_data, n_block_row_num, n_block_column_num);
	else
		return _TyMatrixXdRef(0, 0, 0); // fixme - does this throw an exception or what?
}

inline const double *CUberBlockMatrix::p_GetBlock_Log(size_t n_row_index, size_t n_column_index,
	size_t n_block_row_num, size_t n_block_column_num) const
{
	_ASSERTE(n_row_index < m_block_rows_list.size());
	_ASSERTE(n_column_index < m_block_cols_list.size());
	if(n_block_row_num != m_block_rows_list[n_row_index].n_height ||
	   n_block_column_num != m_block_cols_list[n_column_index].n_width)
		return 0; // size mismatch
	return p_GetBlockData(n_row_index, n_column_index);
}

#ifdef __UBER_BLOCK_MATRIX_IO

inline bool CUberBlockMatrix::Load_MatrixMarket(const char *p_s_filename,
	const char *p_s_layout_filename) // throw(std::bad_alloc)
{
	_ASSERTE(p_s_filename);
	_ASSERTE(p_s_layout_filename); // must always be specified in this version

	if(!Load_BlockLayout(p_s_layout_filename))
		return false;
	// load layout first

	return Load_MatrixMarket_Layout(p_s_filename);
}

#endif // __UBER_BLOCK_MATRIX_IO

inline CUberBlockMatrix::TColumn CUberBlockMatrix::t_ColumnCumsumCopy(const TColumn &r_t_src)
{
	TColumn t_dest;
	t_dest.n_width = r_t_src.n_width; // copy dimensions
	t_dest.n_cumulative_width_sum = r_t_src.n_cumulative_width_sum; // copy cumsums
	_ASSERTE(t_dest.block_list.empty()); // no blocks
	return t_dest;
}

inline CUberBlockMatrix::TRow CUberBlockMatrix::t_ColumnCumsumToRowCumsum(const TColumn &r_t_src)
{
	TRow t_dest;
	t_dest.n_height = r_t_src.n_width; // copy dimensions
	t_dest.n_cumulative_height_sum = r_t_src.n_cumulative_width_sum; // copy cumsums
	// no blocks in dest
	return t_dest;
}

inline CUberBlockMatrix::TColumn CUberBlockMatrix::t_RowCumsumToColumnCumsum(const TRow &r_t_src)
{
	TColumn t_dest;
	t_dest.n_width = r_t_src.n_height; // copy dimensions
	t_dest.n_cumulative_width_sum = r_t_src.n_cumulative_height_sum; // copy cumsums
	_ASSERTE(t_dest.block_list.empty()); // no blocks in src
	return t_dest;
}

template <class _TyDest, class _TyFirst, class _TySecond>
void CUberBlockMatrix::MergeLayout(std::vector<_TyDest> &r_merged_layout, const std::vector<_TyFirst> &r_layout_0,
	const std::vector<_TySecond> &r_layout_1, std::vector<size_t> &r_reindexing_table_0,
	std::vector<size_t> &r_reindexing_table_1) // throw(std::bad_alloc)
{
	r_reindexing_table_0.resize(r_layout_0.size());
	r_reindexing_table_1.resize(r_layout_1.size());
	r_merged_layout.clear(); // note this is always called from inside function where the layout is created, clearing is useless then
	r_merged_layout.reserve(std::max(r_layout_0.size(), r_layout_1.size()));
	// alloc in / out lists

	_ASSERTE(r_layout_0.empty() == r_layout_1.empty()); // both empty or both not
	_ASSERTE(r_layout_0.empty() || r_layout_1.empty() ||
		r_layout_0.back().n_GetCumulative() + r_layout_0.back().n_GetAbsolute() ==
		r_layout_1.back().n_GetCumulative() + r_layout_1.back().n_GetAbsolute()); // same last dimension
	// make sure matrix size is the same

	size_t n_last_cum = 0, i = 0, j = 0, m = r_layout_0.size(), n = r_layout_1.size();
	size_t n_last_first = 0;
	size_t n_last_second = 0;
	while(i < m && j < n) {
		size_t n_cum0 = r_layout_0[i].n_GetCumulative() + r_layout_0[i].n_GetAbsolute();
		size_t n_cum1 = r_layout_1[j].n_GetCumulative() + r_layout_1[j].n_GetAbsolute(); // want the end of the row/col
		size_t n_cum_next = std::min(n_cum0, n_cum1);
		// get cumsums, decide which to advance first

		size_t n_index = r_merged_layout.size();
		if(n_cum0 == n_cum_next) {
			r_reindexing_table_0[i] = (n_last_first == n_last_cum)? n_index : size_t(-1);
			n_last_first = n_cum_next;
			++ i;
		}
		if(n_cum1 == n_cum_next) {
			r_reindexing_table_1[j] = (n_last_second == n_last_cum)? n_index : size_t(-1);
			n_last_second = n_cum_next;
			++ j;
		}
		// build reindexing tables

		_TyDest t_row_or_column;
		t_row_or_column.SetAbsolute(n_cum_next - n_last_cum);
		t_row_or_column.SetCumulative(n_last_cum); // !!
		n_last_cum = n_cum_next;
		// make a new row (or a column, depending on _TyDest)

		r_merged_layout.push_back(t_row_or_column);
		// add it to the list
	}
	_ASSERTE(i == m && j == n); // in case the matrices have the same size, they both end at the same time
}

template <class _TyFirst, class _TySecond>
size_t CUberBlockMatrix::n_MergeLayout(const std::vector<_TyFirst> &r_layout_0,
	const std::vector<_TySecond> &r_layout_1, std::vector<size_t> &r_reindexing_table_0,
	std::vector<size_t> &r_reindexing_table_1) // throw(std::bad_alloc)
{
	r_reindexing_table_0.resize(r_layout_0.size());
	r_reindexing_table_1.resize(r_layout_1.size());
	// alloc out lists

	_ASSERTE(r_layout_0.empty() == r_layout_1.empty()); // both empty or both not
	_ASSERTE(r_layout_0.empty() || r_layout_1.empty() ||
		r_layout_0.back().n_GetCumulative() + r_layout_0.back().n_GetAbsolute() ==
		r_layout_1.back().n_GetCumulative() + r_layout_1.back().n_GetAbsolute()); // same last dimension
	// make sure matrix size is the same

	size_t n_layout_size = 0;
	size_t n_last_cum = 0, i = 0, j = 0, m = r_layout_0.size(), n = r_layout_1.size();
	size_t n_last_first = 0;
	size_t n_last_second = 0;
	while(i < m && j < n) {
		size_t n_cum0 = r_layout_0[i].n_GetCumulative() + r_layout_0[i].n_GetAbsolute();
		size_t n_cum1 = r_layout_1[j].n_GetCumulative() + r_layout_1[j].n_GetAbsolute(); // want the end of the row/col
		size_t n_cum_next = std::min(n_cum0, n_cum1);
		// get cumsums, decide which to advance first

		size_t n_index = n_layout_size;
		if(n_cum0 == n_cum_next) {
			r_reindexing_table_0[i] = (n_last_first == n_last_cum)? n_index : size_t(-1);
			n_last_first = n_cum_next;
			++ i;
		}
		if(n_cum1 == n_cum_next) {
			r_reindexing_table_1[j] = (n_last_second == n_last_cum)? n_index : size_t(-1);
			n_last_second = n_cum_next;
			++ j;
		}
		n_last_cum = n_cum_next;
		// build reindexing tables

		++ n_layout_size;
		// calculate the size of the resulting layout
	}
	_ASSERTE(i == m && j == n); // in case the matrices have the same size, they both end at the same time

	return n_layout_size;
}

inline size_t CUberBlockMatrix::n_RowGet(size_t n_row, size_t &r_n_block_row_num) const
{
	size_t n_row_index;
	if(n_row >= m_n_row_num)
		return -1;
	else {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
		_TyRowConstIter p_row_it = std::upper_bound(m_block_rows_list.begin(),
			m_block_rows_list.end(), n_row, CFindLowerRow()); // t_odo
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		_TyRowConstIter p_row_it = std::upper_bound(m_block_rows_list.begin(),
			m_block_rows_list.end(), n_row, _FindLowerRow); // t_odo
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		//_ASSERTE(p_row_it != m_block_rows_list.end());
		if(p_row_it == m_block_rows_list.begin())
			return -1;
		-- p_row_it;
		const TRow &r_t_row = *p_row_it;

		if(r_t_row.n_cumulative_height_sum == n_row) {
			r_n_block_row_num = r_t_row.n_height;
			n_row_index = p_row_it - m_block_rows_list.begin();
			// just find the index, the row in question does exist, and has the right size (ideal case)

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_row_reref_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		} else {
			return -1;
			// the current row would need to be fragmented
		}
	}
	// resolve row reference

	_ASSERTE(n_row_index < m_block_rows_list.size());
	_ASSERTE(m_block_rows_list[n_row_index].n_height == r_n_block_row_num);
	_ASSERTE(m_block_rows_list[n_row_index].n_cumulative_height_sum == n_row);
	// make sure that the row reference is really resolved

	//_ASSERTE(b_last_row == (n_row_index == m_block_rows_list.size() - 1));
	// makes sure the b_last_row flag is filled correctly

	return n_row_index;
}

inline size_t CUberBlockMatrix::n_ColumnGet(size_t n_column, size_t &r_n_block_column_num) const
{
	size_t n_column_index;
	if(n_column >= m_n_col_num)
		return -1;
	else {
		_TyColumnConstIter p_col_it;
		_ASSERTE(!m_block_cols_list.empty());
		/*if(m_block_cols_list.back().n_cumulative_width_sum == n_column) // saves almost 10% of the time if filling by columns
			p_col_it = m_block_cols_list.end() - 1;
		else*/ {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			p_col_it = std::upper_bound(m_block_cols_list.begin(),
				m_block_cols_list.end(), n_column, CFindLowerColumn()); // t_odo
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			p_col_it = std::upper_bound(m_block_cols_list.begin(),
				m_block_cols_list.end(), n_column, _FindLowerColumn); // t_odo
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			//_ASSERTE(p_col_it != m_block_cols_list.end()); // t_odo - handle this as proper runtime error // done (the line below)
			if(p_col_it == m_block_cols_list.begin())
				return -1;
			-- p_col_it;
		}
		const TColumn &r_t_col = *p_col_it;

		if(r_t_col.n_cumulative_width_sum == n_column) {
			r_n_block_column_num = r_t_col.n_width;
			n_column_index = p_col_it - m_block_cols_list.begin();
			// just find the index, the column in question does exist, and has the right size (ideal case)

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_col_reref_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		} else {
			return -1;
			// handle the subdivision cases
		}
	}
	// resolve column reference (essentially the same code as for the row reference)

	_ASSERTE(n_column_index < m_block_cols_list.size());
	_ASSERTE(m_block_cols_list[n_column_index].n_width == r_n_block_column_num);
	_ASSERTE(m_block_cols_list[n_column_index].n_cumulative_width_sum == n_column);
	// make sure that the column reference is really resolved

	return n_column_index;
}

inline double *CUberBlockMatrix::p_AllocBlockData(size_t n_row_index, size_t n_column_index,
	size_t n_block_row_num, size_t n_block_column_num, bool &r_b_was_a_new_block) // throw(std::bad_alloc)
{
	_ASSERTE(n_column_index < m_block_cols_list.size());
	TColumn &r_t_col = m_block_cols_list[n_column_index];
	return p_AllocBlockData(n_row_index, r_t_col, n_block_row_num,
		n_block_column_num, r_b_was_a_new_block);
}

template <class CMatrixType>
bool CUberBlockMatrix::Append_Block(const CMatrixType &r_t_block, size_t n_row, size_t n_column) // throw(std::bad_alloc)
{
	//CheckIntegrity(); // not here
	// make sure the matrix is ok

	size_t n_row_index, n_column_index;
	if((n_row_index = n_RowAlloc(n_row, r_t_block.rows())) == size_t(-1) ||
	   (n_column_index = n_ColumnAlloc(n_column, r_t_block.cols())) == size_t(-1))
		return false;
	// find row / column

	double *p_data;
	bool b_uninitialized; // ignored
	if(!(p_data = p_AllocBlockData(n_row_index, n_column_index,
	   r_t_block.rows(), r_t_block.cols(), b_uninitialized)))
		return false;
	// allocate a new block / reuse existing one

	Eigen::Map<Eigen::MatrixXd, map_Alignment> dest(p_data, r_t_block.rows(), r_t_block.cols());
	dest = r_t_block; // can't always use memcpy (matrix expressions, matrices with stride, ...)
	// copy the dense data

	return true;
}

template <int n_compile_time_row_num, int n_compile_time_col_num, int n_options>
bool CUberBlockMatrix::Append_Block(const Eigen::Matrix<double, n_compile_time_row_num, n_compile_time_col_num,
	n_options, n_compile_time_row_num, n_compile_time_col_num> &r_t_block, size_t n_row, size_t n_column) // throw(std::bad_alloc)
{
	//CheckIntegrity(); // not here
	// make sure the matrix is ok

	size_t n_row_index, n_column_index;
	if((n_row_index = n_RowAlloc(n_row, r_t_block.rows())) == size_t(-1) ||
	   (n_column_index = n_ColumnAlloc(n_column, r_t_block.cols())) == size_t(-1))
		return false;
	// find row / column

	double *p_data;
	bool b_uninitialized; // ignored
	if(!(p_data = p_AllocBlockData(n_row_index, n_column_index,
	   r_t_block.rows(), r_t_block.cols(), b_uninitialized)))
		return false;
	// allocate a new block / reuse existing one

	memcpy(p_data, &r_t_block(0, 0), r_t_block.rows() * r_t_block.cols() * sizeof(double));
	// copy the dense data

	return true;
}

template <class CMatrixType>
bool CUberBlockMatrix::Append_Block_Log(const CMatrixType &r_t_block, size_t n_row_index, size_t n_column_index) // throw(std::bad_alloc)
{
	_ASSERTE(n_row_index <= m_block_rows_list.size());
	_ASSERTE(n_column_index <= m_block_cols_list.size());

	//CheckIntegrity(); // not here
	// make sure the matrix is ok

	if(n_row_index == m_block_rows_list.size() &&
	   (n_row_index = n_RowAlloc(m_n_row_num, r_t_block.rows())) == size_t(-1))
		return false;
	else if(m_block_rows_list[n_row_index].n_height != r_t_block.rows())
		return false;
	if(n_column_index == m_block_cols_list.size() &&
	   (n_column_index = n_ColumnAlloc(m_n_col_num, r_t_block.cols())) == size_t(-1))
		return false;
	else if(m_block_cols_list[n_column_index].n_width != r_t_block.cols())
		return false;
	// in case either index points one past the last row / column, creates a new one
	// note that n_ColumnAlloc() / n_RowAlloc() check dimensions,
	// if not used, the dimensions must be checked here

	double *p_data;
	bool b_uninitialized; // ignored
	if(!(p_data = p_AllocBlockData(n_row_index, n_column_index,
	   r_t_block.rows(), r_t_block.cols(), b_uninitialized)))
		return false;
	// allocate a new block / reuse existing one

	Eigen::Map<Eigen::MatrixXd, map_Alignment> dest(p_data, r_t_block.rows(), r_t_block.cols());
	dest = r_t_block; // can't always use memcpy (matrix expressions, matrices with stride, ...)
	// copy the dense data

	return false;
}

template <int n_compile_time_row_num, int n_compile_time_col_num, int n_options>
bool CUberBlockMatrix::Append_Block_Log(const Eigen::Matrix<double, n_compile_time_row_num, n_compile_time_col_num,
	n_options, n_compile_time_row_num, n_compile_time_col_num> &r_t_block, size_t n_row_index, size_t n_column_index) // throw(std::bad_alloc)
{
	_ASSERTE(n_row_index <= m_block_rows_list.size());
	_ASSERTE(n_column_index <= m_block_cols_list.size());

	//CheckIntegrity(); // not here
	// make sure the matrix is ok

	if(n_row_index == m_block_rows_list.size() &&
	   (n_row_index = n_RowAlloc(m_n_row_num, r_t_block.rows())) == size_t(-1))
		return false;
	else if(m_block_rows_list[n_row_index].n_height != r_t_block.rows())
		return false;
	if(n_column_index == m_block_cols_list.size() &&
	   (n_column_index = n_ColumnAlloc(m_n_col_num, r_t_block.cols())) == size_t(-1))
		return false;
	else if(m_block_cols_list[n_column_index].n_width != r_t_block.cols())
		return false;
	// in case either index points one past the last row / column, creates a new one
	// note that n_ColumnAlloc() / n_RowAlloc() check dimensions,
	// if not used, the dimensions must be checked here

	double *p_data;
	bool b_uninitialized; // ignored
	if(!(p_data = p_AllocBlockData(n_row_index, n_column_index,
	   r_t_block.rows(), r_t_block.cols(), b_uninitialized)))
		return false;
	// allocate a new block / reuse existing one

	memcpy(p_data, &r_t_block(0, 0), r_t_block.rows() * r_t_block.cols() * sizeof(double));
	// copy the dense data

	return true; // success
}

template <class _TyFixedSizeMatrixOrMap>
void CUberBlockMatrix::Convert_to_Dense(_TyFixedSizeMatrixOrMap &r_dest) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	//r_dest.resize(m_n_row_num, m_n_col_num); // by the caller, also this is specialization for fixed-size matrices
	_ASSERTE(r_dest.rows() == m_n_row_num && r_dest.cols() == m_n_col_num);
	r_dest.setZero();
	// alloc dest matrix

	for(size_t j = 0, n_col_num = m_block_cols_list.size(); j < n_col_num; ++ j) {
		const TColumn &t_col = m_block_cols_list[j];

		size_t n_x = t_col.n_cumulative_width_sum;
		size_t n_width = t_col.n_width;
		// dimensions of the block

		for(size_t i = 0, n_block_num = t_col.block_list.size(); i < n_block_num; ++ i) {
			size_t n_row = t_col.block_list[i].first;
			const double *p_data = t_col.block_list[i].second;

			size_t n_y = m_block_rows_list[n_row].n_cumulative_height_sum;
			size_t n_height = m_block_rows_list[n_row].n_height;
			// dimensions of the block

			for(size_t x = 0; x < n_width; ++ x) {
				for(size_t y = 0; y < n_height; ++ y) {
					double f_elem = p_data[x * n_height + y]; // t_odo - this is probably the other way around // nope, it's this way
					r_dest(n_y + y, n_x + x) = f_elem;
				}
			}
			// iterate through data, sparse fill ...
		}
	}
}

template <class COp>
void CUberBlockMatrix::ElementwiseUnaryOp_ZeroInvariant(COp op)
{
	CheckIntegrity(true);

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn &r_t_col = m_block_cols_list[i];
		const size_t n_block_width = r_t_col.n_width;
		for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
			size_t n_row_idx = r_t_col.block_list[j].first;
			double *p_data = r_t_col.block_list[j].second;
			size_t n_block_height = m_block_rows_list[n_row_idx].n_height;
			// get block position, size and data

			for(const double *p_end = p_data + n_block_width * n_block_height; p_data != p_end; ++ p_data)
				*p_data = op(*p_data);
			// modify the block
		}
	}
}

template <class COp>
void CUberBlockMatrix::ElementwiseUnaryOp_ZeroInvariant_Parallel(COp op)
{
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
			size_t n_row_idx = r_t_col.block_list[j].first;
			double *p_data = r_t_col.block_list[j].second;
			size_t n_block_height = m_block_rows_list[n_row_idx].n_height;
			// get block position, size and data

			for(const double *p_end = p_data + n_block_width * n_block_height; p_data != p_end; ++ p_data)
				*p_data = op(*p_data);
			// modify the block
		}
	}
}

template <class CBinaryOp>
bool CUberBlockMatrix::ElementwiseBinaryOp_ZeroInvariant(CUberBlockMatrix &r_dest, CBinaryOp op) const // throw(std::bad_alloc)
{
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
					memcpy(p_value, p_value_src, n_block_size * sizeof(double));
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

					for(double *p_value_end = p_value_dest + n_block_size;
					   p_value_dest != p_value_end; ++ p_value_src, ++ p_value_dest)
						*p_value_dest = op(*p_value_dest, *p_value_src); // t_odo - replace by op()
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
			memcpy(p_value, p_value_src, n_block_size * sizeof(double));
			// create a new block, initialize with values

			r_t_dest_col.block_list.push_back(TColumn::TBlockEntry(n_new_row, p_value));
			// add it to the list
		}
		// merge blocks
	}

	//r_dest.CheckIntegrity(true);
	// make sure that this operation didn't damage matrix integrity

	return true;
}

template <class CBinaryOp>
bool CUberBlockMatrix::ElementwiseBinaryOp_RightSideZeroInvariant(CUberBlockMatrix &r_dest, CBinaryOp op) const // throw(std::bad_alloc)
{
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
					double *p_value_dest = p_value;
					for(double *p_value_end = p_value_dest + n_block_size;
					   p_value_dest != p_value_end; ++ p_value_src, ++ p_value_dest)
						*p_value_dest = op(0, *p_value_src); // t_odo - replace by op() // note op needs to be zero-invariant (blocks of the second matrix that are missed by blocks of the first matrix are also untouched)
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

					for(double *p_value_end = p_value_dest + n_block_size;
					   p_value_dest != p_value_end; ++ p_value_src, ++ p_value_dest)
						*p_value_dest = op(*p_value_dest, *p_value_src); // t_odo - replace by op()
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
			double *p_value_dest = p_value;
			for(double *p_value_end = p_value_dest + n_block_size;
			   p_value_dest != p_value_end; ++ p_value_src, ++ p_value_dest)
				*p_value_dest = op(0, *p_value_src);  // t_odo - replace by op() // note that op needs to be zero-invariant (blocks of the second matrix that are missed by blocks of the first matrix are also untouched)
			// create a new block, initialize with values

			r_t_dest_col.block_list.push_back(TColumn::TBlockEntry(n_new_row, p_value));
			// add it to the list
		}
		// merge blocks
	}

	//r_dest.CheckIntegrity(true);
	// make sure that this operation didn't damage matrix integrity

	return true;
}

template <class CBinaryOp>
bool CUberBlockMatrix::ElementwiseBinaryOp(CUberBlockMatrix &r_dest, CBinaryOp op) const // throw(std::bad_alloc)
{
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

				for(double *p_value_end = p_value_dest + n_block_size_dest;
				   p_value_dest != p_value_end; ++ p_value_dest)
					*p_value_dest = op(*p_value_dest, 0);
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
					for(double *p_value_end = p_value_dest + n_block_size_dest;
					   p_value_dest != p_value_end; ++ p_value_dest)
						*p_value_dest = op(*p_value_dest, 0);
				}
				// don't forget to modify all the blocks in dest between last one and the one being added to

				bool b_insert_at_end;
				if((b_insert_at_end = (p_block_it == r_t_dest_col.block_list.end())) ||
				   (*p_block_it).first != n_new_row) {
					// a new block

					double *p_value = r_dest.p_Get_DenseStorage(n_block_size);
					double *p_value_dest = p_value;
					for(double *p_value_end = p_value_dest + n_block_size;
					   p_value_dest != p_value_end; ++ p_value_src, ++ p_value_dest)
						*p_value_dest = op(0, *p_value_src); // t_odo - replace by op()
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

					for(double *p_value_end = p_value_dest + n_block_size;
					   p_value_dest != p_value_end; ++ p_value_src, ++ p_value_dest)
						*p_value_dest = op(*p_value_dest, *p_value_src); // t_odo - replace by op()
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
				for(double *p_value_end = p_value_dest + n_block_size_dest;
				   p_value_dest != p_value_end; ++ p_value_dest)
					*p_value_dest = op(*p_value_dest, 0);
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
			for(double *p_value_end = p_value_dest + n_block_size;
			   p_value_dest != p_value_end; ++ p_value_src, ++ p_value_dest)
				*p_value_dest = op(0, *p_value_src);  // t_odo - replace by op()
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

			for(double *p_value_end = p_value_dest + n_block_size_dest;
			   p_value_dest != p_value_end; ++ p_value_dest)
				*p_value_dest = op(*p_value_dest, 0);
		}
	}
	// don't forget to modify all of the columns after the last modified one

	//r_dest.CheckIntegrity(true);
	// make sure that this operation didn't damage matrix integrity

	return true;
}

template <const int MatrixRowsAtCompileTime>
bool CUberBlockMatrix::Cholesky_Dense() // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	//Check_Block_Alignment();

	_ASSERTE(b_SymmetricLayout());
	// makes sure it has symmetric layout

	_ASSERTE(MatrixRowsAtCompileTime == Eigen::Dynamic || n_Row_Num() == MatrixRowsAtCompileTime);
	_ASSERTE(MatrixRowsAtCompileTime == Eigen::Dynamic || n_Column_Num() == MatrixRowsAtCompileTime);
	typedef Eigen::Matrix<double, MatrixRowsAtCompileTime, MatrixRowsAtCompileTime> TMatrixType;

	TMatrixType t_factor(n_Row_Num(), n_Column_Num()); // see what happens here with fixed-size matrices
	t_factor.setZero();
	// alloc dense matrix and clear it

	std::vector<size_t> frontline;
	frontline.resize(n_BlockColumn_Num());
	// alloc vector for frontline blocks (we can observe that cholesky leaves nonzero
	// values in columns above diagonal and under the top-most block of the original matrix)

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn &r_t_col = m_block_cols_list[i];
		if(r_t_col.block_list.empty()) {
			frontline[i] = i; // no block - frontline touches diagonal here
			continue;
		} else
			frontline[i] = r_t_col.block_list.front().first; // use the smallest row id as fronline
		// get frontline

		size_t n_x = r_t_col.n_cumulative_width_sum;
		size_t n_width = r_t_col.n_width;
		// dimensions of the block

		for(size_t i = 0, n_block_num = r_t_col.block_list.size(); i < n_block_num; ++ i) {
			size_t n_row = r_t_col.block_list[i].first;
			const double *p_data = r_t_col.block_list[i].second;

			size_t n_y = m_block_rows_list[n_row].n_cumulative_height_sum;
			size_t n_height = m_block_rows_list[n_row].n_height;
			// dimensions of the block

			for(size_t x = 0; x < n_width; ++ x) {
				for(size_t y = 0; y < n_height; ++ y) {
					double f_elem = p_data[x * n_height + y]; // t_odo - this is probably the other way around // nope, it's this way
					t_factor(n_y + y, n_x + x) = f_elem;
				}
			}
			// iterate through data, sparse fill ...
		}
	}
	// collect the frontline, convert the matrix to dense

	Eigen::LLT<TMatrixType, Eigen::Upper> cholesky(t_factor);
	if(cholesky.info() != Eigen::Success)
		return false; // probably Eigen::NumericalIssue - not pos def
	t_factor = cholesky.matrixU();
	// calculate upper cholesky by eigen

	m_data_pool.clear();
	// clears data pool, invalidates all the existing blocks (faster than
	// sorting the blocks - we already know the correct order)

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		TColumn &r_t_col = m_block_cols_list[i];

		size_t n_x = r_t_col.n_cumulative_width_sum;
		size_t n_width = r_t_col.n_width;
		// dimensions of the block

		size_t n_min_row = frontline[i];
		size_t n_max_row = i + 1; // one past the last one
		// calculates min / max rows that will be nonzero

		if(r_t_col.block_list.capacity() < n_max_row - n_min_row)
			r_t_col.block_list.clear(); // avoid copying data in case it will be reallocated
		r_t_col.block_list.resize(n_max_row - n_min_row);
		// resize the list of blocks

		_TyBlockIter p_block_it = r_t_col.block_list.begin();
		for(size_t n_row = n_min_row; n_row < n_max_row; ++ n_row, ++ p_block_it) {
			TColumn::TBlockEntry &r_t_block = *p_block_it;
			// get block iterator

			size_t n_y = m_block_rows_list[n_row].n_cumulative_height_sum;
			size_t n_height = m_block_rows_list[n_row].n_height;
			// dimensions of the block

			r_t_block.first = n_row;
			double *p_data = r_t_block.second = p_Get_DenseStorage(n_width * n_height);
			// alloc the block

			for(size_t x = 0; x < n_width; ++ x) {
				for(size_t y = 0; y < n_height; ++ y)
					p_data[x * n_height + y] = t_factor(n_y + y, n_x + x); // note this could be faster through FBS and Eigen::Map
			}
			// get cholesky data back
		}
		// fill all the blocks with the data

#ifdef _DEBUG
		for(size_t n_row = 0; n_row < n_min_row; ++ n_row) {
			size_t n_y = m_block_rows_list[n_row].n_cumulative_height_sum;
			size_t n_height = m_block_rows_list[n_row].n_height;
			// dimensions of the block

			for(size_t x = 0; x < n_width; ++ x) {
				for(size_t y = 0; y < n_height; ++ y)
					_ASSERTE(fabs(t_factor(n_y + y, n_x + x)) < 1e-5f); // or directly == 0?
			}
		}
		for(size_t n_row = n_max_row, n_end = m_block_rows_list.size(); n_row < n_end; ++ n_row) {
			size_t n_y = m_block_rows_list[n_row].n_cumulative_height_sum;
			size_t n_height = m_block_rows_list[n_row].n_height;
			// dimensions of the block

			for(size_t x = 0; x < n_width; ++ x) {
				for(size_t y = 0; y < n_height; ++ y)
					_ASSERTE(fabs(t_factor(n_y + y, n_x + x)) < 1e-5f); // or directly == 0?
			}
		}
		// makes sure that all the other entries are really zero
#endif // _DEBUG
	}
	// gets data back from the dense matrix to this matrix

	//CheckIntegrity(true);

	return true;
}

inline bool CUberBlockMatrix::CholeskyOf(const CUberBlockMatrix &r_lambda) // throw(std::bad_alloc)
{
	//r_lambda.CheckIntegrity(true); // inside Build_EliminationTree()
	//r_lambda.Check_Block_Alignment();

	const size_t n = r_lambda.m_block_cols_list.size();
	// get number of block columns

	std::vector<size_t> etree, ereach_stack, bitfield;
	bitfield.resize(n, 0);

	r_lambda.Build_EliminationTree(etree, ereach_stack); // use ereach stack as workspace
	_ASSERTE(ereach_stack.size() == n);
	// build an elimination tree

	return CholeskyOf(r_lambda, etree, ereach_stack, bitfield);
}

inline bool CUberBlockMatrix::CholeskyOf(const CUberBlockMatrix &r_lambda, size_t n_start_on_column) // throw(std::bad_alloc)
{
	//r_lambda.CheckIntegrity(true); // inside Build_EliminationTree()
	//r_lambda.Check_Block_Alignment();

	const size_t n = r_lambda.m_block_cols_list.size();
	// get number of block columns

	std::vector<size_t> etree, ereach_stack, bitfield;
	bitfield.resize(n, 0);

	r_lambda.Build_EliminationTree(etree, ereach_stack); // use ereach stack as workspace
	_ASSERTE(ereach_stack.size() == n);
	// build an elimination tree

	return CholeskyOf(r_lambda, etree, ereach_stack, bitfield, n_start_on_column);
}

#endif // !__UBER_BLOCK_MATRIX_INLINES_INCLUDED
