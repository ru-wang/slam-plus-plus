/*
								+---------------------------------+
								|                                 |
								|  ***   �ber Block Matrix   ***  |
								|                                 |
								| Copyright  (c) -tHE SWINe- 2012 |
								|                                 |
								|         BlockMatrix.cpp         |
								|                                 |
								+---------------------------------+
*/

/**
 *	@file src/slam/BlockMatrix.cpp
 *	@date 2012
 *	@author -tHE SWINe-
 *	@brief the �berblockmatrix class implementation
 */

#include "slam/BlockMatrix.h"

/*
 *								=== CUberBlockMatrix_Base ===
 */

void CUberBlockMatrix_Base::Dump_PerfCounters()
{
#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
	printf("row-reindex-num: %d\n", m_n_row_reindex_num); // nonzero - rows may be split (this is a problem - row reindex references all the blocks)
	printf("rows-list-shift-num: %d\n", m_n_rows_list_shift_num); // equal to m_n_row_reindex_num
	printf("cols-list-shift-num: %d\n", m_n_cols_list_shift_num); // 0 - columns are only appended
	printf("rows-list-realloc-num: %d\n", m_n_rows_list_realloc_num); // 0 - we know the size of the matrix in advance
	printf("cols-list-realloc-num: %d\n", m_n_cols_list_realloc_num); // 0 - we know the size of the matrix in advance
	printf("dummy-row-num: %d\n", m_n_dummy_row_num);
	printf("dummy-col-num: %d\n", m_n_dummy_col_num);
	printf("row-reref-num: %d\n", m_n_row_reref_num); // high numbers here
	printf("col-reref-num: %d\n", m_n_col_reref_num); // high numbers here
	printf("row-append-num: %d\n", m_n_row_append_num);
	printf("col-append-num: %d\n", m_n_col_append_num);
	printf("col-block-search-num: %d\n", m_n_col_block_search_num);
	printf("col-block-reref-num: %d\n", m_n_col_block_reref_num);
#else // __UBER_BLOCK_MATRIX_PERFCOUNTERS
	printf("performance counters not compiled\n");
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
}

void CUberBlockMatrix_Base::Clear()
{
	m_n_row_num = 0;
	m_n_col_num = 0;
	m_block_rows_list.clear();
	m_block_cols_list.clear();

	m_n_ref_elem_num = 0;
	m_data_pool.clear();

	Reset_Perfcounters();
}

void CUberBlockMatrix_Base::Swap(CUberBlockMatrix_Base &r_other)
{
	r_other.m_block_cols_list.swap(m_block_cols_list);
	r_other.m_block_rows_list.swap(m_block_rows_list);
	std::swap(r_other.m_n_col_num, m_n_col_num);
	std::swap(r_other.m_n_row_num, m_n_row_num);
	m_data_pool.swap(r_other.m_data_pool);
	std::swap(m_n_ref_elem_num, r_other.m_n_ref_elem_num);
#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
	std::swap(m_n_row_reindex_num, r_other.m_n_row_reindex_num);
	std::swap(m_n_rows_list_shift_num, r_other.m_n_rows_list_shift_num);
	std::swap(m_n_cols_list_shift_num, r_other.m_n_cols_list_shift_num);
	std::swap(m_n_rows_list_realloc_num, r_other.m_n_rows_list_realloc_num);
	std::swap(m_n_cols_list_realloc_num, r_other.m_n_cols_list_realloc_num);
	std::swap(m_n_dummy_row_num, r_other.m_n_dummy_row_num);
	std::swap(m_n_dummy_col_num, r_other.m_n_dummy_col_num);
	std::swap(m_n_row_reref_num, r_other.m_n_row_reref_num);
	std::swap(m_n_col_reref_num, r_other.m_n_col_reref_num);
	std::swap(m_n_row_append_num, r_other.m_n_row_append_num);
	std::swap(m_n_col_append_num, r_other.m_n_col_append_num);
	std::swap(m_n_col_block_search_num, r_other.m_n_col_block_search_num);
	std::swap(m_n_col_block_reref_num, r_other.m_n_col_block_reref_num);
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
}

void CUberBlockMatrix_Base::CopyTo(CUberBlockMatrix_Base &r_dest, bool b_share_data) const // throw(std::bad_alloc)
{
	_ASSERTE(&r_dest != this); // can't copy to itself (although it would be a no-op)

	CheckIntegrity(true);

#if 0
	SliceTo(r_dest, 0, n_BlockRow_Num(), 0, n_BlockColumn_Num(), b_share_data);
	// can be implemented as that
#else // 0
	_ASSERTE(&r_dest != this); // can't slice to itself

	r_dest.Clear();
	// erase dest ...

	r_dest.m_n_row_num = m_n_row_num;
	r_dest.m_n_col_num = m_n_col_num;
	// set matrix size

	r_dest.m_block_rows_list = m_block_rows_list;
	// set block rows (just copy)

	r_dest.m_block_cols_list.resize(m_block_cols_list.size());
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) { // todo - make iterator loop instead
		TColumn &r_t_dest = r_dest.m_block_cols_list[i];
		const TColumn &r_t_src = m_block_cols_list[i];
		// get src / dest columns

		r_t_dest.n_width = r_t_src.n_width;
		r_t_dest.n_cumulative_width_sum = r_t_src.n_cumulative_width_sum;
		// copy column dimensions

		_TyBlockConstIter p_block_it = r_t_src.block_list.begin();
		_TyBlockConstIter p_end_it = r_t_src.block_list.end();
		r_t_dest.block_list.resize(p_end_it - p_block_it);
		_TyBlockIter p_dest_it = r_t_dest.block_list.begin();
		const size_t n_column_width = r_t_src.n_width;
		if(b_share_data) {
			size_t n_height_sum = 0;
			for(; p_block_it != p_end_it; ++ p_block_it, ++ p_dest_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				// get src block

				*p_dest_it = TColumn::TBlockEntry(r_t_block.first, r_t_block.second);
				// just recalculate row index

				n_height_sum += m_block_rows_list[r_t_block.first].n_height;
			}
			r_dest.m_n_ref_elem_num += r_t_src.n_width * n_height_sum; // !!
		} else {
			for(; p_block_it != p_end_it; ++ p_block_it, ++ p_dest_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				// get src block

				const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
				double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);
				memcpy(p_data, r_t_block.second, n_row_height * n_column_width * sizeof(double));
				// alloc buffer in dest matrix, copy contents

				*p_dest_it = TColumn::TBlockEntry(r_t_block.first, p_data);
				// recalculate row index and use the new buffer
			}
		}
		_ASSERTE(p_dest_it == r_t_dest.block_list.end());
	}
	// set block columns, copy block data (only blocks in the valid column range)

	//r_dest.CheckIntegrity();
	// make sure the dest matrix is ok
#endif // 0
}

/*
 *								=== ~CUberBlockMatrix_Base ===
 */

/*
 *								=== CUberBlockMatrix ===
 */

bool CUberBlockMatrix::b_SymmetricLayout() const
{
	CheckIntegrity();

	if(m_n_row_num != m_n_col_num || m_block_rows_list.size() != m_block_cols_list.size())
		return false;
	_TyRowConstIter p_row_it = m_block_rows_list.begin();
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it, ++ p_row_it) {
		if((*p_col_it).n_cumulative_width_sum != (*p_row_it).n_cumulative_height_sum)
			return false;
		_ASSERTE((*p_col_it).n_width == (*p_row_it).n_height);
		// if cumulative is the same, absolute must be the same as well
	}
	return true;
}

bool CUberBlockMatrix::b_EqualStructure(const CUberBlockMatrix &r_other) const
{
	CheckIntegrity();
	r_other.CheckIntegrity();

	/*if(&r_other == this)
		return true;*/ // don't handle this explicitly, life's not easy
	// simple case

	if(r_other.m_n_col_num != m_n_col_num || r_other.m_n_row_num != m_n_row_num ||
	   r_other.m_block_cols_list.size() != m_block_cols_list.size() ||
	   !(r_other.m_block_rows_list == m_block_rows_list))
		return false;
	// compare dimensions, number of columns and layout of rows

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn r_t_col = m_block_cols_list[i],
			r_t_other_col = r_other.m_block_cols_list[i];
		if(r_t_col.n_width != r_t_other_col.n_width)
			return false;
		_ASSERTE(r_t_col.n_cumulative_width_sum ==
			r_t_other_col.n_cumulative_width_sum); // if widths are equal, cumsums must be equal as well
		// compare block size

		if(r_t_col.block_list.size() != r_t_other_col.block_list.size())
			return false;
		for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
			if(r_t_col.block_list[j].first != r_t_other_col.block_list[j].first)
				return false;
			// make sure there is the same number of blocks
			// and that they reference the same columns
		}
		// compare blocks
	}

	// todo - copy any bug fixes also to b_Equal

	return true;
}

bool CUberBlockMatrix::b_EqualLayout(const CUberBlockMatrix &r_other) const
{
	CheckIntegrity();
	r_other.CheckIntegrity();

	/*if(&r_other == this)
		return true;*/ // don't handle this explicitly, life's not easy
	// simple case

	if(r_other.m_n_col_num != m_n_col_num || r_other.m_n_row_num != m_n_row_num ||
	   r_other.m_block_cols_list.size() != m_block_cols_list.size() ||
	   !(r_other.m_block_rows_list == m_block_rows_list))
		return false;
	// compare dimensions, number of columns and layout of rows

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn r_t_col = m_block_cols_list[i],
			r_t_other_col = r_other.m_block_cols_list[i];
		if(r_t_col.n_width != r_t_other_col.n_width)
			return false;
		// compare block size

		_ASSERTE(r_t_col.n_cumulative_width_sum ==
			r_t_other_col.n_cumulative_width_sum);
		// if widths are equal, cumsums must be equal as well
	}

	// todo - copy any bug fixes also to b_Equal

	return true;
}

bool CUberBlockMatrix::b_Equal(const CUberBlockMatrix &r_other) const
{
	CheckIntegrity(true);
	r_other.CheckIntegrity(true);

	/*if(&r_other == this)
		return true;*/ // don't handle this explicitly, life's not easy
	// simple case

	if(r_other.m_n_col_num != m_n_col_num || r_other.m_n_row_num != m_n_row_num ||
	   r_other.m_block_cols_list.size() != m_block_cols_list.size() ||
	   r_other.m_block_rows_list != m_block_rows_list)
		return false;
	// compare dimensions, number of columns and layout of rows

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn r_t_col = m_block_cols_list[i],
			r_t_other_col = r_other.m_block_cols_list[i];
		if(r_t_col.n_width != r_t_other_col.n_width)
			return false;
		// compare block size

		_ASSERTE(r_t_col.n_cumulative_width_sum ==
			r_t_other_col.n_cumulative_width_sum);
		// if widths are equal, cumsums must be equal as well

		if(r_t_col.block_list.size() != r_t_other_col.block_list.size())
			return false;
		for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
			if(r_t_col.block_list[j].first != r_t_other_col.block_list[j].first)
				return false;
			// make sure there is the same number of blocks
			// and that they reference the same columns

			if(memcmp(r_t_col.block_list[j].second, r_t_other_col.block_list[j].second,
			   r_t_col.n_width * m_block_rows_list[r_t_col.block_list[j].first].n_height *
			   sizeof(double)))
				return false;
			// compare data as well
		}
		// compare blocks
	}

	// todo - copy any bug fixes also to b_EqualStructure

	return true;
}

bool CUberBlockMatrix::b_Equal(const CUberBlockMatrix &r_other, double f_tolerance) const
{
	CheckIntegrity(true);
	r_other.CheckIntegrity(true);

	/*if(&r_other == this)
		return true;*/ // don't handle this explicitly, life's not easy
	// simple case

	if(r_other.m_n_col_num != m_n_col_num || r_other.m_n_row_num != m_n_row_num ||
	   r_other.m_block_cols_list.size() != m_block_cols_list.size() ||
	   r_other.m_block_rows_list != m_block_rows_list)
		return false;
	// compare dimensions, number of columns and layout of rows

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn r_t_col = m_block_cols_list[i],
			r_t_other_col = r_other.m_block_cols_list[i];
		if(r_t_col.n_width != r_t_other_col.n_width)
			return false;
		// compare block size

		_ASSERTE(r_t_col.n_cumulative_width_sum ==
			r_t_other_col.n_cumulative_width_sum);
		// if widths are equal, cumsums must be equal as well

		if(r_t_col.block_list.size() != r_t_other_col.block_list.size())
			return false;
		for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
			if(r_t_col.block_list[j].first != r_t_other_col.block_list[j].first)
				return false;
			// make sure there is the same number of blocks
			// and that they reference the same columns

			const double *p_a = r_t_col.block_list[j].second;
			const double *p_end = p_a + r_t_col.n_width *
				m_block_rows_list[r_t_col.block_list[j].first].n_height;
			const double *p_b = r_t_other_col.block_list[j].second;
			for(; p_a != p_end; ++ p_a, ++ p_b) {
				if(fabs(*p_a - *p_b) > f_tolerance)
					return false;
			}
			// compare data as well, but use tolerance
		}
		// compare blocks
	}

	// todo - copy any bug fixes also to b_EqualStructure

	return true;
}

size_t CUberBlockMatrix::n_NonZero_Num() const
{
	CheckIntegrity(true);

	size_t n_num = 0;
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		size_t n_heights = 0;
		const std::vector<TColumn::TBlockEntry> &r_block_list =
			m_block_cols_list[i].block_list;
		for(size_t j = 0, m = r_block_list.size(); j < m; ++ j)
			n_heights += m_block_rows_list[r_block_list[j].first].n_height;
		// sum up heights of all the blocks in this column

		n_num += n_heights * m_block_cols_list[i].n_width;
		// calculate number of nonzero elements in all the blocks
	}

	return n_num;
}

CUberBlockMatrix::_TyMatrixXdRef CUberBlockMatrix::t_GetBlock_Log(size_t n_row_index,
	size_t n_column_index, size_t n_block_row_num, size_t n_block_column_num,
	bool b_alloc_if_not_found /*= true*/, bool b_mind_uninitialized /*= true*/) // throw(std::bad_alloc)
{
	if(b_alloc_if_not_found) {
		size_t m, n;
		if(n_row_index > (m = m_block_rows_list.size()) ||
		   n_column_index > (n = m_block_cols_list.size()))
			return _TyMatrixXdRef(0, 0, 0);
		// must be inside or just adjacent

		if(n_row_index == m) {
			TRow t_row;
			t_row.n_height = n_block_row_num;
			t_row.n_cumulative_height_sum = m_n_row_num;
			m_block_rows_list.push_back(t_row);
			m_n_row_num += n_block_row_num;
		}
		if(n_column_index == n) {
			TColumn t_col;
			t_col.n_width = n_block_column_num;
			t_col.n_cumulative_width_sum = m_n_col_num;
			m_block_cols_list.push_back(t_col);
			m_n_col_num += n_block_column_num;
		}
		// extend the matrix, if required (p_AllocBlockData requires row and column to exist)

		bool b_uninitialized;
		double *p_data = p_AllocBlockData(n_row_index, n_column_index,
			n_block_row_num, n_block_column_num, b_uninitialized);
		_ASSERTE(p_data); // would throw std::bad_alloc
		if(b_uninitialized && b_mind_uninitialized)
			memset(p_data, 0, n_block_row_num * n_block_column_num * sizeof(double));
		return _TyMatrixXdRef(p_data, n_block_row_num, n_block_column_num);
	} else {
		if(n_row_index >= m_block_rows_list.size() ||
		   n_column_index >= m_block_cols_list.size())
			return _TyMatrixXdRef(0, 0, 0); // fixme - does this throw an exception or what?
		size_t n_found_block_row_num = m_block_rows_list[n_row_index].n_height;
		size_t n_found_block_column_num = m_block_cols_list[n_column_index].n_width;
		double *p_data = (double*)p_GetBlockData(n_row_index, n_column_index);
		if(p_data && n_found_block_row_num == n_block_row_num &&
		   n_found_block_column_num == n_block_column_num)
			return _TyMatrixXdRef(p_data, n_block_row_num, n_block_column_num);
		else
			return _TyMatrixXdRef(0, 0, 0); // fixme - does this throw an exception or what?
	}
}

double *CUberBlockMatrix::p_GetBlock_Log(size_t n_row_index, size_t n_column_index,
	size_t n_block_row_num, size_t n_block_column_num, bool b_alloc_if_not_found /*= true*/,
	bool b_mind_uninitialized /*= true*/) // throw(std::bad_alloc)
{
	if(b_alloc_if_not_found) {
		size_t m, n;
		if(n_row_index > (m = m_block_rows_list.size()) ||
		   n_column_index > (n = m_block_cols_list.size()))
			return 0;
		// must be inside or just adjacent

		if(n_row_index == m) {
			TRow t_row;
			t_row.n_height = n_block_row_num;
			t_row.n_cumulative_height_sum = m_n_row_num;
			m_block_rows_list.push_back(t_row);
			m_n_row_num += n_block_row_num;
		}
		if(n_column_index == n) {
			TColumn t_col;
			t_col.n_width = n_block_column_num;
			t_col.n_cumulative_width_sum = m_n_col_num;
			m_block_cols_list.push_back(t_col);
			m_n_col_num += n_block_column_num;
		}
		// extend the matrix, if required (p_AllocBlockData requires row and column to exist)

		bool b_uninitialized;
		double *p_data = p_AllocBlockData(n_row_index, n_column_index,
			n_block_row_num, n_block_column_num, b_uninitialized);
		if(b_uninitialized && b_mind_uninitialized) {
			_ASSERTE(p_data); // would throw std::bad_alloc
			memset(p_data, 0, n_block_row_num * n_block_column_num * sizeof(double));
		}
		return p_data;
	} else {
		if(n_row_index >= m_block_rows_list.size() ||
		   n_column_index >= m_block_cols_list.size() ||
		   m_block_rows_list[n_row_index].n_height != n_block_row_num ||
		   m_block_cols_list[n_column_index].n_width != n_block_column_num)
			return 0;
		return (double*)p_GetBlockData(n_row_index, n_column_index);
	}
}

double *CUberBlockMatrix::p_GetBlock_Log_Alloc(size_t n_row_index, size_t n_column_index,
	size_t n_block_row_num, size_t n_block_column_num, bool &r_b_is_uninitialized)  // throw(std::bad_alloc)
{
	size_t m, n;
	if(n_row_index > (m = m_block_rows_list.size()) ||
	   n_column_index > (n = m_block_cols_list.size()))
		return 0;
	// must be inside or just adjacent

	if(n_row_index == m) {
		TRow t_row;
		t_row.n_height = n_block_row_num;
		t_row.n_cumulative_height_sum = m_n_row_num;
		m_block_rows_list.push_back(t_row);
		m_n_row_num += n_block_row_num;
	}
	if(n_column_index == n) {
		TColumn t_col;
		t_col.n_width = n_block_column_num;
		t_col.n_cumulative_width_sum = m_n_col_num;
		m_block_cols_list.push_back(t_col);
		m_n_col_num += n_block_column_num;
	}
	// extend the matrix, if required (p_AllocBlockData requires row and column to exist)

	return p_AllocBlockData(n_row_index, n_column_index,
		n_block_row_num, n_block_column_num, r_b_is_uninitialized);
}

double *CUberBlockMatrix::p_FindBlock(size_t n_row, size_t n_column, size_t n_block_row_num,
	size_t n_block_column_num, bool b_alloc_if_not_found /*= true*/,
	bool b_mind_uninitialized /*= true*/) // throw(std::bad_alloc)
{
	//CheckIntegrity(); // not here
	// make sure the matrix is ok

	size_t n_row_index, n_column_index;
	if(b_alloc_if_not_found) {
		if((n_row_index = n_RowAlloc(n_row, n_block_row_num)) == size_t(-1) ||
		   (n_column_index = n_ColumnAlloc(n_column, n_block_column_num)) == size_t(-1))
			return 0;
	} else {
		size_t n_row_rows, n_col_cols;
		if((n_row_index = n_RowGet(n_row, n_row_rows)) == size_t(-1) ||
		   n_row_rows != n_block_row_num ||
		   (n_column_index = n_ColumnGet(n_column, n_col_cols)) == size_t(-1) ||
		   n_col_cols != n_block_column_num)
			return 0;
		// make sure the dimensions are the same
	}
	// alloc row / column

	if(b_alloc_if_not_found) {
		bool b_uninitialized;
		double *p_data = p_AllocBlockData(n_row_index, n_column_index,
			n_block_row_num, n_block_column_num, b_uninitialized);
		if(b_uninitialized && b_mind_uninitialized) {
			_ASSERTE(p_data); // would throw std::bad_alloc
			memset(p_data, 0, n_block_row_num * n_block_column_num * sizeof(double));
		}
		return p_data;
	} else
		return p_GetBlockData(n_row_index, n_column_index);
	// alloc data
}

double *CUberBlockMatrix::p_FindBlock_Alloc(size_t n_row, size_t n_column, size_t n_block_row_num,
	size_t n_block_column_num, bool &r_b_is_uninitialized)
{
	size_t n_row_index, n_column_index;
	if((n_row_index = n_RowAlloc(n_row, n_block_row_num)) == size_t(-1) ||
	   (n_column_index = n_ColumnAlloc(n_column, n_block_column_num)) == size_t(-1))
		return 0;
	return p_AllocBlockData(n_row_index, n_column_index,
		n_block_row_num, n_block_column_num, r_b_is_uninitialized);
}

const double *CUberBlockMatrix::p_FindBlock(size_t n_row, size_t n_column,
	size_t n_block_row_num, size_t n_block_column_num) const
{
	//CheckIntegrity(); not here
	// make sure the matrix is ok

	size_t n_row_index, n_column_index;
	{
		size_t n_row_rows, n_col_cols;
		if((n_row_index = n_RowGet(n_row, n_row_rows)) == size_t(-1) ||
		   n_row_rows != n_block_row_num ||
		   (n_column_index = n_ColumnGet(n_column, n_col_cols)) == size_t(-1) ||
		   n_col_cols != n_block_column_num)
			return 0;
		// make sure the dimensions are the same
	}
	// alloc row / column

	return p_GetBlockData(n_row_index, n_column_index);
	// alloc data
}

const double *CUberBlockMatrix::p_FindBlock_ResolveSize(size_t n_row, size_t n_column,
	size_t &r_n_block_row_num, size_t &r_n_block_column_num) const
{
	//CheckIntegrity(); // not here
	// make sure the matrix is ok

	size_t n_row_index, n_column_index;
	if((n_row_index = n_RowGet(n_row, r_n_block_row_num)) == size_t(-1) ||
	   (n_column_index = n_ColumnGet(n_column, r_n_block_column_num)) == size_t(-1))
		return 0;
	// alloc row / column

	return p_GetBlockData(n_row_index, n_column_index);
	// alloc data
}

size_t CUberBlockMatrix::n_RowAlloc(size_t n_row, size_t n_block_row_num) // throw(std::bad_alloc)
{
	size_t n_row_index;
	if(n_row >= m_n_row_num) {
		if(n_row > m_n_row_num) {
			TRow t_dummy_row;
			t_dummy_row.n_cumulative_height_sum = m_n_row_num;
			t_dummy_row.n_height = n_row - m_n_row_num;
			m_block_rows_list.push_back(t_dummy_row);
#ifdef _DEBUG
			m_n_row_num = n_row; // will be overwritten later ...
#endif // _DEBUG
#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_dummy_row_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		}
		// in case the row is appended further away from the edge,
		// a dummy (empty) row needs to be inserted before

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
		++ m_n_row_append_num;
		if(m_block_rows_list.size() + 1 > m_block_rows_list.capacity())
			++ m_n_rows_list_realloc_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

		TRow t_new_row;
		t_new_row.n_height = n_block_row_num;
		t_new_row.n_cumulative_height_sum = n_row;
		_ASSERTE(m_block_rows_list.empty() || m_block_rows_list.back().n_height +
			m_block_rows_list.back().n_cumulative_height_sum == m_n_row_num); // make sure that the last row either
		m_n_row_num = t_new_row.n_height + t_new_row.n_cumulative_height_sum;
		n_row_index = m_block_rows_list.size(); // remember which
		m_block_rows_list.push_back(t_new_row);
		// make a new row at the end of the matrix to store the block
	} else {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
		_TyRowIter p_row_it = std::upper_bound(m_block_rows_list.begin(),
			m_block_rows_list.end(), n_row, CFindLowerRow()); // t_odo
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		_TyRowIter p_row_it = std::upper_bound(m_block_rows_list.begin(),
			m_block_rows_list.end(), n_row, _FindLowerRow); // t_odo
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		//_ASSERTE(p_row_it != m_block_rows_list.end());
		if(p_row_it == m_block_rows_list.begin())
			return -1;
		-- p_row_it;
		TRow &r_t_row = *p_row_it;

		if(r_t_row.n_cumulative_height_sum == n_row &&
		   r_t_row.n_height == n_block_row_num) {
			n_row_index = p_row_it - m_block_rows_list.begin();
			// just find the index, the row in question does exist, and has the right size (ideal case)

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_row_reref_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		} else {
			_ASSERTE(r_t_row.n_cumulative_height_sum <= n_row); // make sure it's lower bound
			_ASSERTE(r_t_row.n_height > unsigned(n_block_row_num)); // make sure it can possibly fit the matrix

			if(r_t_row.n_height + r_t_row.n_cumulative_height_sum < n_row + n_block_row_num)
				return -1;
			// make sure that the row can contain the block that we'd like to insert

			if(b_RowReferenced(p_row_it - m_block_rows_list.begin()))
				return -1; // the row already contains block of different dimension; can't insert it
			// make sure that the row is not referenced!

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_row_reindex_num;
			++ m_n_rows_list_shift_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

			if(r_t_row.n_cumulative_height_sum == n_row) {
				// the block is inserted at the beginning of a long row

				TRow t_row_after; // not referenced
				t_row_after.n_height = r_t_row.n_height - n_block_row_num;
				t_row_after.n_cumulative_height_sum = r_t_row.n_cumulative_height_sum + n_block_row_num;
				// create a new row to be inserted after the current row

				r_t_row.n_height = n_block_row_num;
				// shorten the current row

				n_row_index = p_row_it - m_block_rows_list.begin();
				// remember row index

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
				if(m_block_rows_list.size() + 1 > m_block_rows_list.capacity())
					++ m_n_rows_list_realloc_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

				m_block_rows_list.insert(p_row_it + 1, 1, t_row_after);
				// insert the new row

				Shift_RowReferences(n_row_index); // fixme: +1 or not? should be.
				// move the indices in the referenced blocks (slow)
			} else if(n_row + n_block_row_num ==
			   r_t_row.n_height + r_t_row.n_cumulative_height_sum) {
				// the block is inserted at the end of a long row

				r_t_row.n_height -= n_block_row_num;
				// shorten the current row

				TRow t_row_after; // referenced
				t_row_after.n_height = n_block_row_num;
				t_row_after.n_cumulative_height_sum = n_row;
				_ASSERTE(t_row_after.n_cumulative_height_sum ==
					r_t_row.n_cumulative_height_sum + r_t_row.n_height);
				// create a new row to be inserted after the current row

				n_row_index = (p_row_it - m_block_rows_list.begin()) + 1;
				// remember row index

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
				if(m_block_rows_list.size() + 1 > m_block_rows_list.capacity())
					++ m_n_rows_list_realloc_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

				m_block_rows_list.insert(p_row_it + 1, 1, t_row_after);
				// insert the new row

				Shift_RowReferences(n_row_index - 1); // fixme: -1 or not? should be.
				// move the indices in the referenced blocks (slow)
			} else {
				// the block is inserted in the middle

				TRow p_row[2];

				TRow &t_row_before = p_row[0];
				t_row_before.n_height = n_row - r_t_row.n_cumulative_height_sum;
				t_row_before.n_cumulative_height_sum = r_t_row.n_cumulative_height_sum;
				// create a new row to be inserted before the current row

				TRow &t_row_after = p_row[1];
				t_row_after.n_cumulative_height_sum = n_row + n_block_row_num;
				t_row_after.n_height = r_t_row.n_height -
					(t_row_after.n_cumulative_height_sum - r_t_row.n_cumulative_height_sum);
				// create a new row to be inserted after the current row

				r_t_row.n_cumulative_height_sum = n_row;
				r_t_row.n_height = n_block_row_num;
				// alter the current row to contain the element in question

				n_row_index = (p_row_it - m_block_rows_list.begin()) + 1;
				// remember row index

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
				if(m_block_rows_list.size() + 2 > m_block_rows_list.capacity())
					++ m_n_rows_list_realloc_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

				std::swap(t_row_before, r_t_row);
				m_block_rows_list.insert(p_row_it + 1, p_row, p_row + 2);
				// insert the new rows

				Shift_RowReferences2(n_row_index - 1, n_row_index + 1);
				// move the indices in the referenced blocks (slow)
			}
			// resolve where the block should be inserted, fragment the current
			// row and renumber row indices in block references
		}
	}
	// resolve row reference

	_ASSERTE(n_row_index < m_block_rows_list.size());
	_ASSERTE(m_block_rows_list[n_row_index].n_height == n_block_row_num);
	_ASSERTE(m_block_rows_list[n_row_index].n_cumulative_height_sum == n_row);
	// make sure that the row reference is really resolved

	//_ASSERTE(b_last_row == (n_row_index == m_block_rows_list.size() - 1));
	// makes sure the b_last_row flag is filled correctly

	//CheckIntegrity(); // leads to exponential cost when inserting N elements in a matrix
	// make sure that adding the row didn't damage matrix integrity

	return n_row_index;
}

size_t CUberBlockMatrix::n_ColumnAlloc(size_t n_column, size_t n_block_column_num) // throw(std::bad_alloc)
{
	size_t n_column_index;
	if(n_column >= m_n_col_num) {
		if(n_column > m_n_col_num) {
			TColumn t_dummy_col;
			t_dummy_col.n_cumulative_width_sum = m_n_col_num;
			t_dummy_col.n_width = n_column - m_n_col_num;
			m_block_cols_list.push_back(t_dummy_col); // not explicitly handled, this does not typically occur in SLAM (works ok, but could be made faster)
#ifdef _DEBUG
			m_n_col_num = n_column; // will be overwritten later ...
#endif // _DEBUG
#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_dummy_col_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		}
		// in case the row is appended further away from the edge,
		// a dummy (empty) column needs to be inserted before

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
		++ m_n_col_append_num;
		if(m_block_cols_list.size() + 1 > m_block_cols_list.capacity())
			++ m_n_cols_list_realloc_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

		TColumn t_new_col; // construct an empty std::vector
		t_new_col.n_width = n_block_column_num;
		t_new_col.n_cumulative_width_sum = n_column;
		//size_t t_new_col_n_width = n_block_column_num;
		//size_t t_new_col_n_cumulative_width_sum = n_column;
		_ASSERTE(m_block_cols_list.empty() || m_block_cols_list.back().n_width +
			m_block_cols_list.back().n_cumulative_width_sum == m_n_col_num); // make sure that the last column either
		m_n_col_num = t_new_col.n_width + t_new_col.n_cumulative_width_sum;
		n_column_index = m_block_cols_list.size(); // remember which

		m_block_cols_list.push_back(t_new_col);
		/*if(m_block_cols_list.capacity() > n_column_index) {
			//m_block_cols_list.resize(n_column_index + 1); // will not realloc & copy many vectors
			//TColumn &r_t_new_col = m_block_cols_list.back();
			TColumn r_t_new_col;
			r_t_new_col.n_width = t_new_col_n_width;
			r_t_new_col.n_cumulative_width_sum = t_new_col_n_cumulative_width_sum;
			m_block_cols_list.push_back(r_t_new_col);
		} else {
			std::vector<TColumn> t_tmp;
			t_tmp.resize(std::max(2 * n_column_index, n_column_index + 1)); // alloc, make sure to double the capacity every time
			//t_tmp.resize(n_column_index + 1); // reserve() resize() does not seem to perform well
			std::vector<TColumn>::iterator p_dst_it = t_tmp.begin();
			for(std::vector<TColumn>::iterator p_src_it = m_block_cols_list.begin(),
			   p_end_it = m_block_cols_list.end(); p_src_it != p_end_it; ++ p_src_it, ++ p_dst_it)
				(*p_src_it).UninitializedSwap(*p_dst_it); // swap the inner vectors instead of copying
			//t_tmp.back() = t_new_col; // copy the last one
			TColumn &r_t_new_col = *p_dst_it;//m_block_cols_list.back();
			r_t_new_col.n_width = t_new_col_n_width;
			r_t_new_col.n_cumulative_width_sum = t_new_col_n_cumulative_width_sum;
			t_tmp.erase(++ p_dst_it, t_tmp.end()); // the rest is for later (not better than the reserve() resize() approach)
			m_block_cols_list.swap(t_tmp); // replace the old list by the new list
		}*/ // this is optimized by uninit_copy on g++, this is just slower
		// make a new column at the end of the matrix to store the block
		// note that this is optimized to avoid copying the inner std::vectors inside TColumn

		// t_odo - insert the new block in this column
	} else {
		_TyColumnIter p_col_it;
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
		TColumn &r_t_col = *p_col_it;

		if(r_t_col.n_cumulative_width_sum == n_column &&
		   r_t_col.n_width == n_block_column_num) {
			n_column_index = p_col_it - m_block_cols_list.begin();
			// just find the index, the column in question does exist, and has the right size (ideal case)

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_col_reref_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		} else {
			_ASSERTE(r_t_col.n_cumulative_width_sum <= n_column); // make sure it's lower bound
			_ASSERTE(r_t_col.n_width > unsigned(n_block_column_num)); // make sure it can possibly fit the matrix (also caught by the runtime condition below)

			if(r_t_col.n_width + r_t_col.n_cumulative_width_sum < n_column + n_block_column_num)
				return -1;
			// make sure that the column can contain the block that we'd like to insert

			if(!r_t_col.block_list.empty())
				return -1; // the column already contains block of different dimension; can't insert it
			// make sure that the column is empty

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_cols_list_shift_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

			if(r_t_col.n_cumulative_width_sum == n_column) {
				// the block is inserted at the beginning of a wide column

				TColumn t_column_after; // not referenced
				t_column_after.n_width = r_t_col.n_width - n_block_column_num;
				t_column_after.n_cumulative_width_sum = r_t_col.n_cumulative_width_sum + n_block_column_num;
				// create a new column to be inserted after the current column

				r_t_col.n_width = n_block_column_num;
				// shorten the current column

				n_column_index = p_col_it - m_block_cols_list.begin();
				// remember column index

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
				if(m_block_cols_list.size() + 1 > m_block_cols_list.capacity())
					++ m_n_cols_list_realloc_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

				m_block_cols_list.insert(p_col_it + 1, 1, t_column_after); // not explicitly handled, this does not typically occur in SLAM (works ok, but could be made faster)
				// insert the new column
			} else if(n_column + n_block_column_num ==
			   r_t_col.n_width + r_t_col.n_cumulative_width_sum) {
				// the block is inserted at the end of a wide column

				r_t_col.n_width -= n_block_column_num;
				// shorten the current column

				TColumn r_t_col_after; // referenced
				r_t_col_after.n_width = n_block_column_num;
				r_t_col_after.n_cumulative_width_sum = n_column;
				_ASSERTE(r_t_col_after.n_cumulative_width_sum ==
					r_t_col.n_cumulative_width_sum + r_t_col.n_width);
				// create a new column to be inserted after the current column

				n_column_index = (p_col_it - m_block_cols_list.begin()) + 1;
				// remember column index

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
				if(m_block_cols_list.size() + 1 > m_block_cols_list.capacity())
					++ m_n_cols_list_realloc_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

				m_block_cols_list.insert(p_col_it + 1, 1, r_t_col_after); // not explicitly handled, this does not typically occur in SLAM (works ok, but could be made faster)
				// insert the new column
			} else {
				// the block is inserted in the middle of a wide column

				TColumn p_col[2];

				TColumn &r_t_col_before = p_col[0];
				r_t_col_before.n_width = n_column - r_t_col.n_cumulative_width_sum;
				r_t_col_before.n_cumulative_width_sum = r_t_col.n_cumulative_width_sum;
				// create a new row to be inserted before the current row

				TColumn &r_t_col_after = p_col[1];
				r_t_col_after.n_cumulative_width_sum = n_column + n_block_column_num;
				r_t_col_after.n_width = r_t_col.n_width -
					(r_t_col_after.n_cumulative_width_sum - r_t_col.n_cumulative_width_sum);
				// create a new row to be inserted after the current row

				r_t_col.n_cumulative_width_sum = n_column;
				r_t_col.n_width = n_block_column_num;
				// alter the current row to contain the element in question

				n_column_index = (p_col_it - m_block_cols_list.begin()) + 1;
				// remember row index

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
				if(m_block_cols_list.size() + 2 > m_block_cols_list.capacity())
					++ m_n_cols_list_realloc_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS

				std::swap(r_t_col_before, r_t_col);
				m_block_cols_list.insert(p_col_it + 1, p_col, p_col + 2); // not explicitly handled, this does not typically occur in SLAM (works ok, but could be made faster)
				// insert the new cols
			}
			// handle the subdivision cases
		}
	}
	// resolve column reference (essentially the same code as for the row reference)

	_ASSERTE(n_column_index < m_block_cols_list.size());
	_ASSERTE(m_block_cols_list[n_column_index].n_width == n_block_column_num);
	_ASSERTE(m_block_cols_list[n_column_index].n_cumulative_width_sum == n_column);
	// make sure that the column reference is really resolved

	//CheckIntegrity(); // leads to exponential cost when inserting N elements in a matrix
	// make sure that adding the column didn't damage matrix integrity

	return n_column_index;
}

bool CUberBlockMatrix::Build_AdditionLayouts(const std::vector<TRow> &r_row_list_first,
	const std::vector<TColumn> &r_column_list_first, std::vector<TRow> &r_row_list_second,
	std::vector<TColumn> &r_column_list_second, std::vector<size_t> &r_reindex_rows_first,
	std::vector<size_t> &r_reindex_columns_first) // throw(std::bad_alloc)
{
	std::vector<TRow> merged_row_list;
	std::vector<size_t> reindex_rows_second;
	MergeLayout(merged_row_list, r_row_list_first, r_row_list_second, r_reindex_rows_first, reindex_rows_second);
	// amortized O(row blocks)

	std::vector<TColumn> merged_column_list;
	std::vector<size_t> reindex_cols_second;
	MergeLayout(merged_column_list, r_column_list_first, r_column_list_second, r_reindex_columns_first, reindex_cols_second);
	// amortized O(col blocks)

	if(r_row_list_second.size() == merged_row_list.size()) { // same size = same layout // fixme - or is it? (asserted below)
		_ASSERTE(reindex_rows_second.size() == merged_row_list.size());
		for(size_t i = 0, n = merged_row_list.size(); i < n; ++ i) {
			_ASSERTE(reindex_rows_second[i] == i); // 1:1 mapping
			_ASSERTE(r_row_list_second[i].n_cumulative_height_sum == merged_row_list[i].n_cumulative_height_sum);
			_ASSERTE(r_row_list_second[i].n_height == merged_row_list[i].n_height);
		}
		// make sure the layout indeed did not change
	} else {
		for(size_t i = 0, n = r_column_list_second.size(); i < n; ++ i) {
			TColumn &r_t_col = r_column_list_second[i];
			for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
				TColumn::TBlockEntry &r_t_block = r_t_col.block_list[j];
				if((r_t_block.first = reindex_rows_second[r_t_block.first]) == size_t(-1))
					return false; // non-empty row was split (different block layout; note that dest is left damaged)
			}
			// for all the blocks in the current column ...
		}
		// reindex block rows in all the columns

		merged_row_list.swap(r_row_list_second);
		// use the merged row layout in the dest matrix
	}
	// resolve rows (zero or O(row blocks))

	if(r_column_list_second.size() == merged_column_list.size()) { // same size = same layout // fixme - or is it? (asserted below)
		_ASSERTE(reindex_cols_second.size() == merged_column_list.size());
		for(size_t i = 0, n = merged_column_list.size(); i < n; ++ i) {
			_ASSERTE(reindex_cols_second[i] == i); // 1:1 mapping
			_ASSERTE(r_column_list_second[i].n_cumulative_width_sum == merged_column_list[i].n_cumulative_width_sum);
			_ASSERTE(r_column_list_second[i].n_width == merged_column_list[i].n_width);
		}
		// make sure the layout indeed did not change
	} else { // the layout somehow changed; scatter old columns in new columns and swap back
		for(size_t i = 0, n = r_column_list_second.size(); i < n; ++ i) {
			TColumn &r_t_col = r_column_list_second[i];
			if(r_t_col.block_list.empty())
				continue;
			// skip empty columns ...

			size_t n_dest = reindex_cols_second[i];
			if(n_dest == size_t(-1))
				return false; // this column was split, the blocks in it are now homeless (different block layout; note that dest is left damaged)
			merged_column_list[n_dest].block_list.swap(r_t_col.block_list);
			// put the column in an reindexed location
		}
		// scatter the columns

		merged_column_list.swap(r_column_list_second);
		// use the merged column layout in the dest matrix
	}
	// resolve columns (zero or O(col blocks))

	return true;
}

void CUberBlockMatrix::PreMultiply_Add(double *p_dest_vector, size_t UNUSED(n_dest_size),
	const double *p_src_vector, size_t UNUSED(n_src_size)) const // @todo - this is largerly untested; test this
{
	CheckIntegrity(true);

	_ASSERTE(p_dest_vector != p_src_vector);
	_ASSERTE(n_dest_size == m_n_row_num);
	_ASSERTE(n_src_size == m_n_col_num);

	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;
		const size_t n_column = r_t_col.n_cumulative_width_sum; // source vector offset
		//Eigen::R_eferencingMatrixXd src(1, r_t_col.n_width, (double*)p_src_vector + n_column); // source vector (column vector)
		_TyVectorXdRef src((double*)p_src_vector + n_column, r_t_col.n_width);
		// forces unaligned memory - no SSE2 (p_src_vector can't be aligned due to block size)

		// note that this loop can be omp parallelized as there will be no two blocks
		// at the same row. it will be also most likely very short (would have to have
		// some kind of test of number of blocks and decide parallel / sequential).
		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
		   p_block_it != p_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			// for each column, for each block ...

			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const size_t n_row = r_t_row.n_cumulative_height_sum; // dest vector offset

			//Eigen::R_eferencingMatrixXd dest(1, r_t_row.n_height, p_dest_vector + n_row); // dest vector (column vector)

			_TyVectorXdRef dest(p_dest_vector + n_row, r_t_row.n_height);
			// forces unaligned memory - no SSE2 (p_dest_vector can't be aligned due to block size)

			_TyMatrixXdRef block((double*)r_t_block.second, r_t_row.n_height, r_t_col.n_width);
			// block

#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			dest += block.lazyProduct(src); // axpy
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			dest.noalias() += block * src; // axpy
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			// perform the multiplication, one block at a time
		}
	}
}

void CUberBlockMatrix::PostMultiply_Add(double *p_dest_vector, size_t UNUSED(n_dest_size),
	const double *p_src_vector, size_t UNUSED(n_src_size)) const // @todo - this is largerly untested; test this
{
	CheckIntegrity(true);

	_ASSERTE(p_dest_vector != p_src_vector);
	_ASSERTE(n_dest_size == m_n_col_num);
	_ASSERTE(n_src_size == m_n_row_num);

	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;
		const size_t n_column = r_t_col.n_cumulative_width_sum; // dest vector offset
		//Eigen::R_eferencingMatrixXd dest(r_t_col.n_width, 1, p_dest_vector + n_column); // dest vector (row vector)
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
		_ASSERTE(r_t_col.n_width == 3);
		_TyRowVector3dRef dest(p_dest_vector + n_column);
		// forces unaligned memory - no SSE2 (p_dest_vector can't be aligned due to block size)
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
		_TyRowVectorXdRef dest(p_dest_vector + n_column, r_t_col.n_width);
		// forces unaligned memory - no SSE2 (p_dest_vector can't be aligned due to block size)
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
		   p_block_it != p_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			// for each column, for each block ...

			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const size_t n_row = r_t_row.n_cumulative_height_sum; // source vector offset

			//Eigen::R_eferencingMatrixXd src(r_t_row.n_height, 1, (double*)p_src_vector + n_row);
			// source vector (row vector)

/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
			_ASSERTE(r_t_col.n_width == 3 && r_t_row.n_height == 3);
			_TyRowVector3dRef src((double*)p_src_vector + n_row);
			// forces unaligned memory - no SSE2 (p_src_vector can't be aligned due to block size)
			_TyMatrix3dRef block((double*)r_t_block.second); // block
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
			_TyRowVectorXdRef src((double*)p_src_vector + n_row, r_t_row.n_height);
			// forces unaligned memory - no SSE2 (p_src_vector can't be aligned due to block size)
			_TyMatrixXdRef block((double*)r_t_block.second,
				r_t_row.n_height, r_t_col.n_width); // block
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS

#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			dest += src.lazyProduct(block); // axpy
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			dest.noalias() += src * block; // axpy
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			// perform the multiplication, one block at a time
		}
	}
}

void CUberBlockMatrix::PostMultiply_Add_Parallel(double *p_dest_vector, size_t UNUSED(n_dest_size),
	const double *p_src_vector, size_t UNUSED(n_src_size)) const // @todo - this is largerly untested; test this
{
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
		const size_t n_column = r_t_col.n_cumulative_width_sum; // dest vector offset
		//Eigen::R_eferencingMatrixXd dest(r_t_col.n_width, 1, p_dest_vector + n_column); // dest vector (row vector)
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
		_ASSERTE(r_t_col.n_width == 3);
		_TyRowVector3dRef dest(p_dest_vector + n_column); // forces unaligned memory - no SSE2 (p_dest_vector can't be aligned due to block size)
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
		_TyRowVectorXdRef dest(p_dest_vector + n_column, r_t_col.n_width); // forces unaligned memory - no SSE2 (p_dest_vector can't be aligned due to block size)
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
		   p_block_it != p_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			// for each column, for each block ...

			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const size_t n_row = r_t_row.n_cumulative_height_sum; // source vector offset

			//Eigen::R_eferencingMatrixXd src(r_t_row.n_height, 1, (double*)p_src_vector + n_row); // source vector (row vector)
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
			_ASSERTE(r_t_col.n_width == 3 && r_t_row.n_height == 3);
			_TyRowVector3dRef src((double*)p_src_vector + n_row); // forces unaligned memory - no SSE2 (p_src_vector can't be aligned due to block size)
			_TyMatrix3dRef block((double*)r_t_block.second); // block
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
			_TyRowVectorXdRef src((double*)p_src_vector + n_row, r_t_row.n_height); // forces unaligned memory - no SSE2 (p_src_vector can't be aligned due to block size)
			_TyMatrixXdRef block((double*)r_t_block.second, r_t_row.n_height, r_t_col.n_width); // block
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS

#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			dest += src.lazyProduct(block); // axpy
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			dest.noalias() += src * block; // axpy
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			// perform the multiplication, one block at a time
		}
	}
}

#ifdef __UBER_BLOCK_MATRIX_HAVE_CSPARSE

cs *CUberBlockMatrix::p_Convert_to_Sparse_UpperTriangular(cs *p_alloc /*= 0*/) const // t_odo - will need a *32 version of this attrocity
{
	if(!p_alloc) {
		if(!(p_alloc = cs_spalloc(m_n_row_num, m_n_col_num,
		   m_n_ref_elem_num + m_data_pool.size(), 1, 0))) // this will not realloc (will alloc a bit more, in fact, but should not be too much more)
			return 0;
	} else {
		if((p_alloc->nzmax < 0 || unsigned(p_alloc->nzmax) < m_n_ref_elem_num +
		   m_data_pool.size()) && !cs_sprealloc(p_alloc, std::max(size_t(p_alloc->nzmax * 2),
		   m_n_ref_elem_num + m_data_pool.size()))) // this will not realloc (will alloc a bit more, in fact, but should not be too much more)
			return 0;
		if((p_alloc->n < 0 || size_t(p_alloc->n) < m_n_col_num) &&
		   !(p_alloc->p = (csi*)realloc(p_alloc->p, (m_n_col_num + 1) * sizeof(csi)))) // cs_sprealloc doesn't realloc column array
			return 0;
		p_alloc->m = m_n_row_num;
		p_alloc->n = m_n_col_num;
	}
	_ASSERTE(p_alloc->nzmax > 0);
	_ASSERTE(unsigned(p_alloc->nzmax) >= m_n_ref_elem_num + m_data_pool.size());
	// prepare dest matrix

	csi *p_column = p_alloc->p, *p_index = p_alloc->i;
	double *p_value = p_alloc->x;
	// access values

	size_t n_nonzero_element_num = 0, n_cur_column = 0;
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;
		const size_t n_width = r_t_col.n_width;
		_ASSERTE(n_width > 0);
		// get column

		if(r_t_col.block_list.empty()) {
			for(size_t x = 0; x < n_width; ++ x, ++ n_cur_column, ++ p_column)
				*p_column = n_nonzero_element_num;
			// fill the column structures (empty columns)
		} else {
			const TColumn::TBlockEntry &r_t_last_block = r_t_col.block_list.back();
			const TRow &r_t_last_row = m_block_rows_list[r_t_last_block.first];
			// get last block in the column

			if(r_t_last_row.n_cumulative_height_sum + r_t_last_row.n_height < n_cur_column) {
				// the column is completely off-diagonal, can use the fast code from p_Convert_to_Sparse()

				size_t n_column_row_num = n_nonzero_element_num;
				for(_TyBlockConstIter p_block_it =
				   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
				   p_block_it != p_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					const TRow &r_t_row = m_block_rows_list[r_t_block.first];
					const double *p_data = r_t_block.second;
					// for each column, for each block ...

					n_nonzero_element_num += r_t_row.n_height;
					_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
					// will write r_t_row.n_height values

					for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
					   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
						*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
						*p_index = n_row;
					}
					// write values to the arrays
				}
				n_column_row_num = n_nonzero_element_num - n_column_row_num; // no need to keep two counters in the loop
				// first loop collects number of rows in the column

				_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
				_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

				for(size_t i = 0; i < n_width; ++ i, ++ p_column)
					*p_column = n_nonzero_element_num - n_column_row_num + i * n_column_row_num;
				n_cur_column += n_width; // !!
				// fill the column structures (works on empty columns as well)

				if(n_width > 1) {
					p_value -= n_column_row_num;
					p_index -= n_column_row_num;
					// make those zero-based indexable

					_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num - n_column_row_num);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num - n_column_row_num);

					n_nonzero_element_num += (n_width - 1) * n_column_row_num;
					_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
					// will write n_column_row_num values, n_width - 1 times (once per each of the rest of element-columns)

					for(_TyBlockConstIter p_block_it =
					   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
					   p_block_it != p_end_it; ++ p_block_it) {
						const TColumn::TBlockEntry &r_t_block = *p_block_it;
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							for(size_t x = 1; x < n_width; ++ x) {
								p_value[x * n_column_row_num] = p_data[x * r_t_row.n_height + i];
								p_index[x * n_column_row_num] = n_row;
							}
						}
						// write values to the arrays
					}
					// process element-columns 1 to n_width - 1

					_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num -
						n_column_row_num * (n_width - 1));
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num -
						n_column_row_num * (n_width - 1));

					p_value += n_column_row_num * (n_width - 1);
					p_index += n_column_row_num * (n_width - 1);
					// shift after the last value used

					_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);
				}
				// second loop processes all the element-columns at once
			} else if(false && (r_t_last_row.n_cumulative_height_sum <= n_cur_column &&
			   r_t_last_row.n_cumulative_height_sum + r_t_last_row.n_height >=
			   n_cur_column + n_width)) { // fie! this code is actually slower than the more generic version in the last branch (well, it is simpler)
				// the column ends with one block that covers the entire diagonal (will happen in symmetric matrices)
				// number of rows in one column = number of rows at first column + index of row (elementwise)

				size_t n_column0_row_num = n_nonzero_element_num;
				for(_TyBlockConstIter p_block_it =
				   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end() - 1;
				   p_block_it != p_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					const TRow &r_t_row = m_block_rows_list[r_t_block.first];
					const double *p_data = r_t_block.second;
					// for each column, for each block ...

					_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column);
					//if(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column) {
						n_nonzero_element_num += r_t_row.n_height;
						_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
						// will write r_t_row.n_height values

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
							*p_index = n_row;
						}
						// write values to the arrays
					/*} else { // note this applies only and exactly on the last block - could separate the loops to avoid condition
						n_nonzero_element_num += n_cur_column - r_t_row.n_cumulative_height_sum + 1;
						_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
						// will write r_t_row.n_height values

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   n_row <= n_cur_column; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
							*p_index = n_row;
						}
						// write values to the arrays
					}*/
				}
				{
					const TColumn::TBlockEntry &r_t_block = r_t_col.block_list.back();
					const TRow &r_t_row = m_block_rows_list[r_t_block.first];
					const double *p_data = r_t_block.second;
					// for each column, for each block ...

					_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height > n_cur_column);
					{ // note this applies only and exactly on the last block - could separate the loops to avoid condition
						n_nonzero_element_num += n_cur_column -
							r_t_row.n_cumulative_height_sum + 1;
						_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
						// will write r_t_row.n_height values

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   n_row <= n_cur_column; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
							*p_index = n_row;
						}
						// write values to the arrays
					}
				}
				n_column0_row_num = n_nonzero_element_num - n_column0_row_num; // no need to keep two counters in the loop
				// first loop collects number of rows in the column

				_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
				_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

				n_nonzero_element_num -= n_column0_row_num;
				for(size_t i = 0; i < n_width; ++ i, ++ p_column) {
					*p_column = n_nonzero_element_num;
					n_nonzero_element_num += n_column0_row_num + i;
				}
				//n_nonzero_element_num += n_column0_row_num;
				_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
				// fill the column structures (works on empty columns as well)

				if(n_width > 1) {
					p_value -= n_column0_row_num;
					p_index -= n_column0_row_num;
					// make those zero-based indexable

					//_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num - n_column0_row_num);
					//_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num - n_column0_row_num);
					// already updated n_nonzero_element_num, can't check

					for(_TyBlockConstIter p_block_it =
					   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end() - 1;
					   p_block_it != p_end_it; ++ p_block_it) {
						const TColumn::TBlockEntry &r_t_block = *p_block_it;
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						_ASSERTE(r_t_row.n_cumulative_height_sum +
							r_t_row.n_height <= n_cur_column);
						/*if(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column) {*/
							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
								for(size_t x = 1; x < n_width; ++ x) {
									p_value[x * n_column0_row_num + x - 1] =
										p_data[x * r_t_row.n_height + i];
									p_index[x * n_column0_row_num + x - 1] = n_row;
								}
							}
							// write values to the arrays
						/*} else {
							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   n_row < n_cur_column + n_width; ++ i, ++ n_row, ++ p_value, ++ p_index) {
								size_t n_first = (n_row > n_cur_column)? n_row - n_cur_column : 1;
								for(size_t x = n_first; x < n_width; ++ x) {
									p_value[x * n_column0_row_num + x - 1] = p_data[x * r_t_row.n_height + i];
									p_index[x * n_column0_row_num + x - 1] = n_row;
								}
							}
							// write values to the arrays
						}*/
					}
					{
						const TColumn::TBlockEntry &r_t_block = r_t_col.block_list.back();
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						_ASSERTE(r_t_row.n_cumulative_height_sum +
							r_t_row.n_height > n_cur_column);
						{
							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   n_row < n_cur_column + n_width; ++ i, ++ n_row,
							   ++ p_value, ++ p_index) {
								size_t n_first = (n_row > n_cur_column)? n_row - n_cur_column : 1;
								for(size_t x = n_first; x < n_width; ++ x) {
									p_value[x * n_column0_row_num + x - 1] =
										p_data[x * r_t_row.n_height + i];
									p_index[x * n_column0_row_num + x - 1] = n_row;
								}
							}
							// write values to the arrays
						}
					}
					// process element-columns 1 to n_width - 1

					size_t n_other_area = 0;
					for(size_t i = 0; i < n_width - 1; ++ i)
						n_other_area += n_column0_row_num + i; // todo this can be optimized away
					//n_other_area -= n_width * r_t_col.block_list.size(); // area without the first rows of blocks

					_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num - n_other_area);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num - n_other_area);

					p_value += n_other_area;
					p_index += n_other_area;
					// shift after the last value used

					_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);
				}
				// second loop processes all the element-columns at once

				n_cur_column += n_width; // !!
			} else {
				// the geometry is complex, will handle it a slow way (should never happen in symmetric layouts)
				// note the code below was tested, though

				for(size_t x = 0; x < n_width; ++ x, ++ n_cur_column, ++ p_column) {
					*p_column = n_nonzero_element_num;
					// fill the column structure

					//size_t n_column_row_num = n_nonzero_element_num;
					for(_TyBlockConstIter p_block_it =
					   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end() /*- 1*/;
					   p_block_it != p_end_it; ++ p_block_it) {
						const TColumn::TBlockEntry &r_t_block = *p_block_it;
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						//_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column); // not true, there might be >1 column below diagonal, tha's why it's complex
						if(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column) {
							n_nonzero_element_num += r_t_row.n_height;
							_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
							// will write r_t_row.n_height values

							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
								*p_value = p_data[x * r_t_row.n_height + i];
								*p_index = n_row;
							}
							// write values to the arrays
						} else if(n_cur_column + 1 > r_t_row.n_cumulative_height_sum) { // if there's anything to write at all ...
							n_nonzero_element_num += n_cur_column -
								r_t_row.n_cumulative_height_sum + 1;
							_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
							// will write as many values as it fits above the diagonal (and on it)

							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   n_row <= n_cur_column; ++ i, ++ n_row, ++ p_value, ++ p_index) {
								*p_value = p_data[x * r_t_row.n_height + i];
								*p_index = n_row;
							}
							// write values to the arrays
						}
					}
					/*{
						const TColumn::TBlockEntry &r_t_block = r_t_col.block_list.back();
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height > n_cur_column); // the last block crosses diagonal
						n_nonzero_element_num += n_cur_column - r_t_row.n_cumulative_height_sum + 1;
						_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
						// will write as many values as it fits above the diagonal (and on it)

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   n_row <= n_cur_column; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							*p_value = p_data[x * r_t_row.n_height + i];
							*p_index = n_row;
						}
						// write values to the arrays
					}*/
					//n_column_row_num = n_nonzero_element_num - n_column_row_num; // no need to keep two counters in the loop
					// first loop collects number of rows in the column
					// note that this loop is rather nasty. it would be nicer to see if the last block is either a) completely off diagonal, b) covering the whole diagonal (block diagonal = matrix diagonal) c) other. a simpler code could be used for a) and b) (the innermost loop could contain x, blocks could be iterated just once, all column pointers could be written in a short loop) and the detection is rather simple // todo

					_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);
				}
				// per each elementwise-column
			}
		}
	}
	// writes out the matrix in compressed column form directly without triplet form and cs_compress()

	_ASSERTE(p_column == p_alloc->p + p_alloc->n); // make sure it's pointing to the last column
	_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num);
	*p_column = n_nonzero_element_num;
	_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
	_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

	return p_alloc;
}

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)

/* allocate a sparse matrix (triplet form or compressed-column form) */
static cs *cs_spalloc32(csi m, csi n, csi nzmax, csi values, csi triplet)
{
	cs *A = (cs*)cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
	if (!A) return (NULL) ;                 /* out of memory */
	A->m = m ;                              /* define dimensions and nzmax */
	A->n = n ;
	A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
	A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
	A->p = (csi*)cs_malloc (triplet ? nzmax : n+1, sizeof (int)) ;
	A->i = (csi*)cs_malloc (nzmax, sizeof (int)) ;
	A->x = (double*)(values ? cs_malloc (nzmax, sizeof (double)) : NULL) ;
	return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
static csi cs_sprealloc32(cs *A, csi nzmax)
{
	csi ok, oki, okj = 1, okx = 1 ;
	if (!A) return (0) ;
	if (nzmax <= 0) nzmax = (CS_CSC (A)) ? (A->p [A->n]) : A->nz ;
	A->i = (csi*)cs_realloc (A->i, nzmax, sizeof (int), &oki) ;
	if (CS_TRIPLET (A)) A->p = (csi*)cs_realloc (A->p, nzmax, sizeof (int), &okj) ;
	if (A->x) A->x = (double*)cs_realloc (A->x, nzmax, sizeof (double), &okx) ;
	ok = (oki && okj && okx) ;
	if (ok) A->nzmax = nzmax ;
	return (ok) ;
}

#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)

cs *CUberBlockMatrix::p_Convert_to_Sparse_UpperTriangular32(cs *p_alloc /*= 0*/) const // t_odo - will need a *32 version of this attrocity
{
	if(!p_alloc) {
		if(!(p_alloc = cs_spalloc32(m_n_row_num, m_n_col_num,
		   m_n_ref_elem_num + m_data_pool.size(), 1, 0))) // this will not realloc (will alloc a bit more, in fact, but should not be too much more)
			return 0;
	} else {
		if((p_alloc->nzmax < 0 || unsigned(p_alloc->nzmax) < m_n_ref_elem_num +
		   m_data_pool.size()) && !cs_sprealloc32(p_alloc, std::max(size_t(p_alloc->nzmax * 2),
		   m_n_ref_elem_num + m_data_pool.size()))) // this will not realloc (will alloc a bit more, in fact, but should not be too much more)
			return 0;
		if((p_alloc->n < 0 || size_t(p_alloc->n) < m_n_col_num) &&
		   !(p_alloc->p = (csi*)realloc(p_alloc->p, (m_n_col_num + 1) * sizeof(int)))) // cs_sprealloc doesn't realloc column array
			return 0;
		p_alloc->m = m_n_row_num;
		p_alloc->n = m_n_col_num;
	}
	_ASSERTE(p_alloc->nzmax > 0);
	_ASSERTE(unsigned(p_alloc->nzmax) >= m_n_ref_elem_num + m_data_pool.size());
	// prepare dest matrix

	int *p_column = (int*)p_alloc->p, *p_index = (int*)p_alloc->i;
	double *p_value = p_alloc->x;
	// access values

	size_t n_nonzero_element_num = 0, n_cur_column = 0;
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;
		const size_t n_width = r_t_col.n_width;
		_ASSERTE(n_width > 0);
		// get column

		if(r_t_col.block_list.empty()) {
			_ASSERTE(n_nonzero_element_num <= INT_MAX);
			for(size_t x = 0; x < n_width; ++ x, ++ n_cur_column, ++ p_column)
				*p_column = int(n_nonzero_element_num);
			// fill the column structures (empty columns)
		} else {
			const TColumn::TBlockEntry &r_t_last_block = r_t_col.block_list.back();
			const TRow &r_t_last_row = m_block_rows_list[r_t_last_block.first];
			// get last block in the column

			if(r_t_last_row.n_cumulative_height_sum + r_t_last_row.n_height < n_cur_column) {
				// the column is completely off-diagonal, can use the fast code from p_Convert_to_Sparse()

				size_t n_column_row_num = n_nonzero_element_num;
				for(_TyBlockConstIter p_block_it =
				   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
				   p_block_it != p_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					const TRow &r_t_row = m_block_rows_list[r_t_block.first];
					const double *p_data = r_t_block.second;
					// for each column, for each block ...

					n_nonzero_element_num += r_t_row.n_height;
					_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
					// will write r_t_row.n_height values

					for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
					   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
						*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
						_ASSERTE(n_row <= INT_MAX);
						*p_index = int(n_row);
					}
					// write values to the arrays
				}
				n_column_row_num = n_nonzero_element_num - n_column_row_num; // no need to keep two counters in the loop
				// first loop collects number of rows in the column

				_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
				_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

				for(size_t i = 0; i < n_width; ++ i, ++ p_column) {
					_ASSERTE(n_nonzero_element_num - n_column_row_num +
						i * n_column_row_num <= INT_MAX);
					*p_column = int(n_nonzero_element_num -
						n_column_row_num + i * n_column_row_num);
				}
				n_cur_column += n_width; // !!
				// fill the column structures (works on empty columns as well)

				if(n_width > 1) {
					p_value -= n_column_row_num;
					p_index -= n_column_row_num;
					// make those zero-based indexable

					_ASSERTE(p_index == ((int*)p_alloc->i) +
						n_nonzero_element_num - n_column_row_num);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num - n_column_row_num);

					n_nonzero_element_num += (n_width - 1) * n_column_row_num;
					_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
					// will write n_column_row_num values, n_width - 1 times (once per each of the rest of element-columns)

					for(_TyBlockConstIter p_block_it =
					   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
					   p_block_it != p_end_it; ++ p_block_it) {
						const TColumn::TBlockEntry &r_t_block = *p_block_it;
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							for(size_t x = 1; x < n_width; ++ x) {
								p_value[x * n_column_row_num] = p_data[x * r_t_row.n_height + i];
								_ASSERTE(n_row <= INT_MAX);
								p_index[x * n_column_row_num] = int(n_row);
							}
						}
						// write values to the arrays
					}
					// process element-columns 1 to n_width - 1

					_ASSERTE(p_index == ((int*)p_alloc->i) +
						n_nonzero_element_num - n_column_row_num * (n_width - 1));
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num -
						n_column_row_num * (n_width - 1));

					p_value += n_column_row_num * (n_width - 1);
					p_index += n_column_row_num * (n_width - 1);
					// shift after the last value used

					_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);
				}
				// second loop processes all the element-columns at once
			} else if(false && (r_t_last_row.n_cumulative_height_sum <= n_cur_column &&
			   r_t_last_row.n_cumulative_height_sum + r_t_last_row.n_height >=
			   n_cur_column + n_width)) { // fie! this code is actually slower than the more generic version in the last branch (well, it is simpler)
				// the column ends with one block that covers the entire diagonal (will happen in symmetric matrices)
				// number of rows in one column = number of rows at first column + index of row (elementwise)

				size_t n_column0_row_num = n_nonzero_element_num;
				for(_TyBlockConstIter p_block_it =
				   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end() - 1;
				   p_block_it != p_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					const TRow &r_t_row = m_block_rows_list[r_t_block.first];
					const double *p_data = r_t_block.second;
					// for each column, for each block ...

					_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column);
					//if(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column) {
						n_nonzero_element_num += r_t_row.n_height;
						_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
						// will write r_t_row.n_height values

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
							_ASSERTE(n_row <= INT_MAX);
							*p_index = int(n_row);
						}
						// write values to the arrays
					/*} else { // note this applies only and exactly on the last block - could separate the loops to avoid condition
						n_nonzero_element_num += n_cur_column -
							r_t_row.n_cumulative_height_sum + 1;
						_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
						// will write r_t_row.n_height values

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   n_row <= n_cur_column; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
							*p_index = n_row;
						}
						// write values to the arrays
					}*/
				}
				{
					const TColumn::TBlockEntry &r_t_block = r_t_col.block_list.back();
					const TRow &r_t_row = m_block_rows_list[r_t_block.first];
					const double *p_data = r_t_block.second;
					// for each column, for each block ...

					_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height > n_cur_column);
					{ // note this applies only and exactly on the last block - could separate the loops to avoid condition
						n_nonzero_element_num += n_cur_column - r_t_row.n_cumulative_height_sum + 1;
						_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
						// will write r_t_row.n_height values

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   n_row <= n_cur_column; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
							_ASSERTE(n_row <= INT_MAX);
							*p_index = int(n_row);
						}
						// write values to the arrays
					}
				}
				n_column0_row_num = n_nonzero_element_num - n_column0_row_num; // no need to keep two counters in the loop
				// first loop collects number of rows in the column

				_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
				_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

				n_nonzero_element_num -= n_column0_row_num;
				for(size_t i = 0; i < n_width; ++ i, ++ p_column) {
					_ASSERTE(n_nonzero_element_num <= INT_MAX);
					*p_column = int(n_nonzero_element_num);
					n_nonzero_element_num += n_column0_row_num + i;
				}
				//n_nonzero_element_num += n_column0_row_num;
				_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
				// fill the column structures (works on empty columns as well)

				if(n_width > 1) {
					p_value -= n_column0_row_num;
					p_index -= n_column0_row_num;
					// make those zero-based indexable

					//_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num - n_column0_row_num);
					//_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num - n_column0_row_num);
					// already updated n_nonzero_element_num, can't check

					for(_TyBlockConstIter p_block_it =
					   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end() - 1;
					   p_block_it != p_end_it; ++ p_block_it) {
						const TColumn::TBlockEntry &r_t_block = *p_block_it;
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						_ASSERTE(r_t_row.n_cumulative_height_sum +
							r_t_row.n_height <= n_cur_column);
						/*if(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column) {*/
							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
								for(size_t x = 1; x < n_width; ++ x) {
									p_value[x * n_column0_row_num + x - 1] =
										p_data[x * r_t_row.n_height + i];
									_ASSERTE(n_row <= INT_MAX);
									p_index[x * n_column0_row_num + x - 1] = int(n_row);
								}
							}
							// write values to the arrays
						/*} else {
							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   n_row < n_cur_column + n_width; ++ i, ++ n_row,
							   ++ p_value, ++ p_index) {
								size_t n_first = (n_row > n_cur_column)? n_row - n_cur_column : 1;
								for(size_t x = n_first; x < n_width; ++ x) {
									p_value[x * n_column0_row_num + x - 1] =
										p_data[x * r_t_row.n_height + i];
									p_index[x * n_column0_row_num + x - 1] = n_row;
								}
							}
							// write values to the arrays
						}*/
					}
					{
						const TColumn::TBlockEntry &r_t_block = r_t_col.block_list.back();
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						_ASSERTE(r_t_row.n_cumulative_height_sum +
							r_t_row.n_height > n_cur_column);
						{
							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   n_row < n_cur_column + n_width; ++ i, ++ n_row,
							   ++ p_value, ++ p_index) {
								size_t n_first = (n_row > n_cur_column)? n_row - n_cur_column : 1;
								for(size_t x = n_first; x < n_width; ++ x) {
									p_value[x * n_column0_row_num + x - 1] =
										p_data[x * r_t_row.n_height + i];
									_ASSERTE(n_row <= INT_MAX);
									p_index[x * n_column0_row_num + x - 1] = int(n_row);
								}
							}
							// write values to the arrays
						}
					}
					// process element-columns 1 to n_width - 1

					size_t n_other_area = 0;
					for(size_t i = 0; i < n_width - 1; ++ i)
						n_other_area += n_column0_row_num + i; // todo this can be optimized away
					//n_other_area -= n_width * r_t_col.block_list.size(); // area without the first rows of blocks

					_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num - n_other_area);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num - n_other_area);

					p_value += n_other_area;
					p_index += n_other_area;
					// shift after the last value used

					_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);
				}
				// second loop processes all the element-columns at once

				n_cur_column += n_width; // !!
			} else {
				// the geometry is complex, will handle it a slow way (should never happen in symmetric layouts)
				// note the code below was tested, though

				for(size_t x = 0; x < n_width; ++ x, ++ n_cur_column, ++ p_column) {
					_ASSERTE(n_nonzero_element_num <= INT_MAX);
					*p_column = int(n_nonzero_element_num);
					// fill the column structure

					//size_t n_column_row_num = n_nonzero_element_num;
					for(_TyBlockConstIter p_block_it =
					   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end() /*- 1*/;
					   p_block_it != p_end_it; ++ p_block_it) {
						const TColumn::TBlockEntry &r_t_block = *p_block_it;
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						//_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column); // not true, there might be >1 column below diagonal, tha's why it's complex
						if(r_t_row.n_cumulative_height_sum + r_t_row.n_height <= n_cur_column) {
							n_nonzero_element_num += r_t_row.n_height;
							_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
							// will write r_t_row.n_height values

							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
								*p_value = p_data[x * r_t_row.n_height + i];
								_ASSERTE(n_row <= INT_MAX);
								*p_index = int(n_row);
							}
							// write values to the arrays
						} else if(n_cur_column + 1 > r_t_row.n_cumulative_height_sum) { // if there's anything to write at all ...
							n_nonzero_element_num += n_cur_column -
								r_t_row.n_cumulative_height_sum + 1;
							_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
							// will write as many values as it fits above the diagonal (and on it)

							for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
							   n_row <= n_cur_column; ++ i, ++ n_row, ++ p_value, ++ p_index) {
								*p_value = p_data[x * r_t_row.n_height + i];
								_ASSERTE(n_row <= INT_MAX);
								*p_index = int(n_row);
							}
							// write values to the arrays
						}
					}
					/*{
						const TColumn::TBlockEntry &r_t_block = r_t_col.block_list.back();
						const TRow &r_t_row = m_block_rows_list[r_t_block.first];
						const double *p_data = r_t_block.second;
						// for each column, for each block ...

						_ASSERTE(r_t_row.n_cumulative_height_sum +
							r_t_row.n_height > n_cur_column); // the last block crosses diagonal
						n_nonzero_element_num += n_cur_column -
							r_t_row.n_cumulative_height_sum + 1;
						_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
						// will write as many values as it fits above the diagonal (and on it)

						for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
						   n_row <= n_cur_column; ++ i, ++ n_row, ++ p_value, ++ p_index) {
							*p_value = p_data[x * r_t_row.n_height + i];
							*p_index = n_row;
						}
						// write values to the arrays
					}*/
					//n_column_row_num = n_nonzero_element_num - n_column_row_num; // no need to keep two counters in the loop
					// first loop collects number of rows in the column
					// note that this loop is rather nasty. it would be nicer to see if the last block is either a) completely off diagonal, b) covering the whole diagonal (block diagonal = matrix diagonal) c) other. a simpler code could be used for a) and b) (the innermost loop could contain x, blocks could be iterated just once, all column pointers could be written in a short loop) and the detection is rather simple // todo

					_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
					_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);
				}
				// per each elementwise-column
			}
		}
	}
	// writes out the matrix in compressed column form directly without triplet form and cs_compress()

	_ASSERTE(p_column == ((int*)p_alloc->p) + p_alloc->n); // make sure it's pointing to the last column
	_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num);
	_ASSERTE(n_nonzero_element_num <= INT_MAX);
	*p_column = int(n_nonzero_element_num);
	_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
	_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

	return p_alloc;
}

#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64

#ifdef _OPENMP

template <class _TIntType>
inline void CUberBlockMatrix::From_Sparse_ParallelInnerLoop(const cs *p_sparse,
	/*omp_lock_t &UNUSED(t_mutex),*/ TColumn &r_t_col, const size_t n_base_row_base,
	const size_t n_base_col_base, const std::vector<size_t> &r_workspace) // throw(std::bad_alloc)
{
	size_t n_base_col = r_t_col.n_cumulative_width_sum;
	size_t n_width = r_t_col.n_width;
	// get column

	const _TIntType *const p_cs_col_ptr = (const _TIntType*)p_sparse->p;
	const _TIntType *const p_cs_row_index = (const _TIntType*)p_sparse->i;
	const double *const p_cs_x = p_sparse->x;
	// make sure that these are also constant pointers

	size_t p_col_ptr_column[16];
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
	size_t p_col_ptr[16], p_col_end[16]; // only up to 16x16 blocks // todo - resolve this (wrap with one more loop)
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
	const _TIntType *p_col_ptr[16], *p_col_end[16];
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
	size_t n_col_ptr_num = 0;
	const _TIntType *p_cs_col = p_cs_col_ptr + n_base_col - n_base_col_base; // mess
	size_t n_min_row = m_n_row_num;
	for(size_t i = 0; i < n_width; ++ i) {
		if(p_cs_col[i] != p_cs_col[i + 1]) {
			_ASSERTE(p_cs_row_index[p_cs_col[i]] >= 0 && p_cs_row_index[p_cs_col[i]] <= SIZE_MAX);
			if(size_t(p_cs_row_index[p_cs_col[i]]) < n_min_row)
				n_min_row = p_cs_row_index[p_cs_col[i]] ;
			_ASSERTE(n_col_ptr_num < 16); // todo
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			p_col_ptr[n_col_ptr_num] = p_cs_col[i];
			p_col_end[n_col_ptr_num] = p_cs_col[i + 1];
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			p_col_ptr[n_col_ptr_num] = p_cs_row_index + p_cs_col[i];
			p_col_end[n_col_ptr_num] = p_cs_row_index + p_cs_col[i + 1];
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
#ifdef _DEBUG
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			for(size_t j = p_cs_col[i] + 1, m = p_cs_col[i + 1]; j < m; ++ j)
				_ASSERTE(p_cs_row_index[j - 1] < p_cs_row_index[j]);
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			for(const _TIntType *j = p_cs_row_index + p_cs_col[i] + 1,
			   *m = p_cs_row_index + p_cs_col[i + 1]; j < m; ++ j)
				_ASSERTE(*(j - 1) < *j);
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			// make sure the row indices are sorted - otherwise this function does not work
#endif // _DEBUG
			p_col_ptr_column[n_col_ptr_num] = i;
			++ n_col_ptr_num;
		}
	}
	n_min_row += n_base_row_base;
	// see on which row the sparse data begin in this block column
	// complexity O(N) in block column width

	if(!n_col_ptr_num)
		return;
	// handle empty columns

	size_t n_min_row_id = r_workspace[n_min_row];
	// get row id

	_TyBlockIter p_block_it = r_t_col.block_list.begin();
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
	p_block_it = std::lower_bound(p_block_it, r_t_col.block_list.end(),
		n_min_row_id, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	p_block_it = std::lower_bound(p_block_it, r_t_col.block_list.end(),
		n_min_row_id, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	// get block iterator - note the use of lower_bound() !!
	// complexity O(log N) in average column fill

	for(;;) {
		const TRow &r_t_min_row = m_block_rows_list[n_min_row_id];
		size_t n_height = r_t_min_row.n_height;
		size_t n_base_row = r_t_min_row.n_cumulative_height_sum;
		size_t n_last_row = n_base_row + n_height;
		// get row

		if(p_block_it != r_t_col.block_list.end()) {
			if((*p_block_it).first < n_min_row_id) {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
				p_block_it = std::lower_bound(p_block_it, r_t_col.block_list.end(),
					n_min_row_id, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				p_block_it = std::lower_bound(p_block_it, r_t_col.block_list.end(),
					n_min_row_id, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				// get block iterator - note the use of lower_bound() !!
				// complexity O(log N) in average column fill
			}
			_ASSERTE((*p_block_it).first >= n_min_row_id);
			if((*p_block_it).first > n_min_row_id) {
				// the block with n_min_row_id does not exist in this column yet

				//omp_set_lock(&t_mutex);
				/*double *p_block;
				#pragma omp critical
				{
					p_block = p_Get_DenseStorage(n_width * n_height);
				}*/
				//omp_unset_lock(&t_mutex);
				//memset(p_block, 0, n_width * n_height * sizeof(double));
				p_block_it = r_t_col.block_list.insert(p_block_it,
					TColumn::TBlockEntry(n_min_row_id, 0));
				// put the new block before the existing block
			}
		} else {
			_ASSERTE(r_t_col.block_list.empty() || r_t_col.block_list.back().first < n_min_row_id);
			// the last block in this column has lower row id

			//omp_set_lock(&t_mutex);
			/*double *p_block;
			#pragma omp critical
			{
				p_block = p_Get_DenseStorage(n_width * n_height);
			}*/
			//omp_unset_lock(&t_mutex);
			//memset(p_block, 0, n_width * n_height * sizeof(double));
			r_t_col.block_list.push_back(TColumn::TBlockEntry(n_min_row_id, 0));
			p_block_it = r_t_col.block_list.end() - 1;
			// put the new block at the end of the list
		}
		if(!(*p_block_it).second) { // null means we inserted a new block and need to alloc (there must be only a single synchronization point)
			#pragma omp critical
			{
				(*p_block_it).second = p_Get_DenseStorage(n_width * n_height);
			}
			memset((*p_block_it).second, 0, n_width * n_height * sizeof(double));
		}
		double *p_block = (*p_block_it).second;
		_ASSERTE((*p_block_it).first == n_min_row_id);
		// make sure the block is up to date

		size_t n_min_next_row = m_n_row_num;
		for(size_t i = 0; i < n_col_ptr_num; ++ i) {
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			size_t n_col_ptr = p_col_ptr[i]; // read the column pointer
			size_t n_col_end = p_col_end[i]; // antialiass
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			const _TIntType *n_col_ptr = p_col_ptr[i]; // read the column pointer
			const _TIntType *n_col_end = p_col_end[i]; // antialiass
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			_ASSERTE(n_col_ptr < n_col_end);
			double *p_block_pointer = p_block + (n_height * p_col_ptr_column[i]);
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			_ASSERTE(p_cs_row_index[n_col_ptr] + n_base_row_base >= n_base_row); // n_base_row is chosen on minima of pointers; all must be above or equal
			if(p_cs_row_index[n_col_ptr] + n_base_row_base < n_last_row) for(;;) {
				size_t n_row_in_block = p_cs_row_index[n_col_ptr] + n_base_row_base - n_base_row;
				p_block_pointer[n_row_in_block] = p_cs_x[n_col_ptr];
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			_ASSERTE(*n_col_ptr + n_base_row_base >= n_base_row); // n_base_row is chosen on minima of pointers; all must be above or equal
			size_t n_min_col_row;
			if((n_min_col_row = *n_col_ptr + n_base_row_base) < n_last_row) {
				for(;;) {
					size_t n_row_in_block = *n_col_ptr + n_base_row_base - n_base_row;
					p_block_pointer[n_row_in_block] = p_cs_x[n_col_ptr - p_cs_row_index];
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
					if(++ n_col_ptr == n_col_end) {
						// this column does not contain any more elements

						if(-- n_col_ptr_num) {
							p_col_ptr_column[i] = p_col_ptr_column[n_col_ptr_num];
							p_col_ptr[i] = p_col_ptr[n_col_ptr_num];
							p_col_end[i] = p_col_end[n_col_ptr_num];
							-- i; // !!
						}
						// swap for column from the end of the list (order doesn't really matter)

						break;
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
					} else if((n_row_in_block = p_cs_row_index[n_col_ptr] +
						n_base_row_base) >= n_last_row) {
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
					} else if((n_row_in_block = *n_col_ptr + n_base_row_base) >= n_last_row) {
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
						// this will be written to the next block(s)

						p_col_ptr[i] = n_col_ptr;
						// store state for the next iteration

						if(n_min_next_row > n_row_in_block)
							n_min_next_row = n_row_in_block;
						// remember next min row to write to

						break;
					}
				}
			} else {
				// this column does not have any elements in this row

				if(n_min_next_row > n_min_col_row)
					n_min_next_row = n_min_col_row;
				// remember next min row to write to
			}
		}
		// go through column pointers and write inside the current block
		// complexity O(N) in number of elements to write

		if(!n_col_ptr_num)
			return;
		// all rows in these columns completely written

		++ p_block_it; // note this is not final - the data might skip >1 block (handled at the beginning of this loop)
		_ASSERTE(n_min_next_row < m_n_row_num);
		n_min_row_id = r_workspace[n_min_next_row]; // t_odo - this is suboptimal, we should know what is the next minimal row to write and use that
		// shift to the next row / block
	}
}

#endif // _OPENMP

template <class _TIntType>
inline void CUberBlockMatrix::From_Sparse_InnerLoop(const cs *p_sparse, TColumn &r_t_col,
	const size_t n_base_row_base, const size_t n_base_col_base,
	const std::vector<size_t> &r_workspace) // throw(std::bad_alloc)
{
	size_t n_base_col = r_t_col.n_cumulative_width_sum;
	size_t n_width = r_t_col.n_width;
	// get column

	const _TIntType *const p_cs_col_ptr = (const _TIntType*)p_sparse->p;
	const _TIntType *const p_cs_row_index = (const _TIntType*)p_sparse->i;
	const double *const p_cs_x = p_sparse->x;
	// make sure that these are also constant pointers

	size_t p_col_ptr_column[16];
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
	size_t p_col_ptr[16], p_col_end[16]; // only up to 16x16 blocks // todo - resolve this (wrap with one more loop)
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
	const _TIntType *p_col_ptr[16], *p_col_end[16];
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
	size_t n_col_ptr_num = 0;
	const _TIntType *p_cs_col = p_cs_col_ptr + n_base_col - n_base_col_base; // mess
	size_t n_min_row = m_n_row_num;
	for(size_t i = 0; i < n_width; ++ i) {
		if(p_cs_col[i] != p_cs_col[i + 1]) {
			_ASSERTE(p_cs_row_index[p_cs_col[i]] >= 0 && p_cs_row_index[p_cs_col[i]] <= SIZE_MAX);
			if(size_t(p_cs_row_index[p_cs_col[i]]) < n_min_row)
				n_min_row = p_cs_row_index[p_cs_col[i]] ;
			_ASSERTE(n_col_ptr_num < 16); // todo
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			p_col_ptr[n_col_ptr_num] = p_cs_col[i];
			p_col_end[n_col_ptr_num] = p_cs_col[i + 1];
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			p_col_ptr[n_col_ptr_num] = p_cs_row_index + p_cs_col[i];
			p_col_end[n_col_ptr_num] = p_cs_row_index + p_cs_col[i + 1];
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
#ifdef _DEBUG
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			for(size_t j = p_cs_col[i] + 1, m = p_cs_col[i + 1]; j < m; ++ j)
				_ASSERTE(p_cs_row_index[j - 1] < p_cs_row_index[j]);
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			for(const _TIntType *j = p_cs_row_index + p_cs_col[i] + 1,
			   *m = p_cs_row_index + p_cs_col[i + 1]; j < m; ++ j)
				_ASSERTE(*(j - 1) < *j);
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			// make sure the row indices are sorted - otherwise this function does not work
#endif // _DEBUG
			p_col_ptr_column[n_col_ptr_num] = i;
			++ n_col_ptr_num;
		}
	}
	n_min_row += n_base_row_base;
	// see on which row the sparse data begin in this block column
	// complexity O(N) in block column width

	if(!n_col_ptr_num)
		return;
	// handle empty columns

	size_t n_min_row_id = r_workspace[n_min_row];
	// get row id

	_TyBlockIter p_block_it = r_t_col.block_list.begin();
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
	p_block_it = std::lower_bound(p_block_it, r_t_col.block_list.end(),
		n_min_row_id, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	p_block_it = std::lower_bound(p_block_it, r_t_col.block_list.end(),
		n_min_row_id, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	// get block iterator - note the use of lower_bound() !!
	// complexity O(log N) in average column fill

	for(;;) {
		const TRow &r_t_min_row = m_block_rows_list[n_min_row_id];
		size_t n_height = r_t_min_row.n_height;
		size_t n_base_row = r_t_min_row.n_cumulative_height_sum;
		size_t n_last_row = n_base_row + n_height;
		// get row

		if(p_block_it != r_t_col.block_list.end()) {
			if((*p_block_it).first < n_min_row_id) {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
				p_block_it = std::lower_bound(p_block_it, r_t_col.block_list.end(),
					n_min_row_id, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				p_block_it = std::lower_bound(p_block_it, r_t_col.block_list.end(),
					n_min_row_id, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				// get block iterator - note the use of lower_bound() !!
				// complexity O(log N) in average column fill
			}
			_ASSERTE((*p_block_it).first >= n_min_row_id);
			if((*p_block_it).first > n_min_row_id) {
				// the block with n_min_row_id does not exist in this column yet

				double *p_block = p_Get_DenseStorage(n_width * n_height);
				memset(p_block, 0, n_width * n_height * sizeof(double));
				p_block_it = r_t_col.block_list.insert(p_block_it,
					TColumn::TBlockEntry(n_min_row_id, p_block));
				// put the new block before the existing block
			}
		} else {
			_ASSERTE(r_t_col.block_list.empty() || r_t_col.block_list.back().first < n_min_row_id);
			// the last block in this column has lower row id

			double *p_block = p_Get_DenseStorage(n_width * n_height);
			memset(p_block, 0, n_width * n_height * sizeof(double));
			r_t_col.block_list.push_back(TColumn::TBlockEntry(n_min_row_id, p_block));
			p_block_it = r_t_col.block_list.end() - 1;
			// put the new block at the end of the list
		}
		double *p_block = (*p_block_it).second;
		_ASSERTE((*p_block_it).first == n_min_row_id);
		// make sure the block is up to date

		size_t n_min_next_row = m_n_row_num;
		for(size_t i = 0; i < n_col_ptr_num; ++ i) {
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			size_t n_col_ptr = p_col_ptr[i]; // read the column pointer
			size_t n_col_end = p_col_end[i]; // antialiass
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			const _TIntType *n_col_ptr = p_col_ptr[i]; // read the column pointer
			const _TIntType *n_col_end = p_col_end[i]; // antialiass
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			_ASSERTE(n_col_ptr < n_col_end);
			double *p_block_pointer = p_block + (n_height * p_col_ptr_column[i]);
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			_ASSERTE(p_cs_row_index[n_col_ptr] + n_base_row_base >= n_base_row); // n_base_row is chosen on minima of pointers; all must be above or equal
			if(p_cs_row_index[n_col_ptr] + n_base_row_base < n_last_row) for(;;) {
				size_t n_row_in_block = p_cs_row_index[n_col_ptr] + n_base_row_base - n_base_row;
				p_block_pointer[n_row_in_block] = p_cs_x[n_col_ptr];
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
			_ASSERTE(*n_col_ptr + n_base_row_base >= n_base_row); // n_base_row is chosen on minima of pointers; all must be above or equal
			size_t n_min_col_row;
			if((n_min_col_row = *n_col_ptr + n_base_row_base) < n_last_row) {
				for(;;) {
					size_t n_row_in_block = *n_col_ptr + n_base_row_base - n_base_row;
					p_block_pointer[n_row_in_block] = p_cs_x[n_col_ptr - p_cs_row_index];
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
					if(++ n_col_ptr == n_col_end) {
						// this column does not contain any more elements

						if(-- n_col_ptr_num) {
							p_col_ptr_column[i] = p_col_ptr_column[n_col_ptr_num];
							p_col_ptr[i] = p_col_ptr[n_col_ptr_num];
							p_col_end[i] = p_col_end[n_col_ptr_num];
							-- i; // !!
						}
						// swap for column from the end of the list (order doesn't really matter)

						break;
#ifdef __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
					} else if((n_row_in_block = p_cs_row_index[n_col_ptr] +
						n_base_row_base) >= n_last_row) {
#else // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
					} else if((n_row_in_block = *n_col_ptr + n_base_row_base) >= n_last_row) {
#endif // __UBER_BLOCK_MATRIX_FROM_SPARSE_ITERATE_INTEGERS
						// this will be written to the next block(s)

						p_col_ptr[i] = n_col_ptr;
						// store state for the next iteration

						if(n_min_next_row > n_row_in_block)
							n_min_next_row = n_row_in_block;
						// remember next min row to write to

						break;
					}
				}
			} else {
				// this column does not have any elements in this row

				if(n_min_next_row > n_min_col_row)
					n_min_next_row = n_min_col_row;
				// remember next min row to write to
			}
		}
		// go through column pointers and write inside the current block
		// complexity O(N) in number of elements to write

		if(!n_col_ptr_num)
			return;
		// all rows in these columns completely written

		++ p_block_it; // note this is not final - the data might skip >1 block (handled at the beginning of this loop)
		_ASSERTE(n_min_next_row < m_n_row_num);
		n_min_row_id = r_workspace[n_min_next_row]; // t_odo - this is suboptimal, we should know what is the next minimal row to write and use that
		// shift to the next row / block
	}
}

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)

bool CUberBlockMatrix::From_Sparse32(size_t n_base_row_id, size_t n_base_column_id,
	const cs *p_sparse, bool b_transpose /*= false*/, std::vector<size_t> &r_workspace) // throw(std::bad_alloc)
{
	if(b_transpose)
		return false;
	// for the moment we only handle the not-transpose case

	size_t n_base_row_base = m_block_rows_list[n_base_row_id].n_cumulative_height_sum;
	size_t n_base_col_base = m_block_cols_list[n_base_column_id].n_cumulative_width_sum;
	// get offset of csparse data

	if(p_sparse->m + n_base_row_base > m_n_row_num || p_sparse->n + n_base_col_base > m_n_col_num)
		return false;
	// make sure the sparse data fit inside the existing layout

	if(r_workspace.size() > m_n_row_num)
		r_workspace.clear(); // couldn't be possibly filled from here; clear
	if(r_workspace.size() < m_n_row_num) {
		size_t n_start = r_workspace.size();
		r_workspace.resize(m_n_row_num);
		if(!n_start) {
			size_t n_row_index = 0;
			for(_TyRowConstIter p_row_it = m_block_rows_list.begin(), p_end_it =
			   m_block_rows_list.end(); p_row_it != p_end_it; ++ p_row_it, ++ n_row_index) {
				for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
				   n_start < n; ++ n_start)
					r_workspace[n_start] = n_row_index;
			}
			// fill workspace with row index lookup
		} else {
			_TyRowIter p_row_it;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_start, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_start, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(p_row_it != m_block_rows_list.begin())
				-- p_row_it;
			_ASSERTE((*p_row_it).n_cumulative_height_sum <= n_start);
			_ASSERTE((*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height > n_start);
			// find where to start

			size_t n_row_index = p_row_it - m_block_rows_list.begin();
			for(_TyRowConstIter p_end_it = m_block_rows_list.end();
			   p_row_it != p_end_it; ++ p_row_it, ++ n_row_index) {
				for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
				   n_start < n; ++ n_start)
					r_workspace[n_start] = n_row_index;
			}
			// fill workspace with row index lookup

#ifdef _DEBUG
			{
				size_t n_row_index = 0;
				for(_TyRowConstIter p_row_it = m_block_rows_list.begin(),
				   p_end_it = m_block_rows_list.end(); p_row_it != p_end_it;
				   ++ p_row_it, ++ n_row_index) {
					for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
					   n_start < n; ++ n_start)
						_ASSERTE(r_workspace[n_start] == n_row_index);
				}
			}
#endif // _DEBUG
		}
	}
	// creates a lookup table; complexity O(N) in number of element-rows
	// todo - put this inside a function

	for(_TyColumnIter p_col_it = m_block_cols_list.begin() + n_base_column_id,
	   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) {
		TColumn &r_t_col = *p_col_it;
		From_Sparse_InnerLoop<int>(p_sparse, r_t_col, n_base_row_base,
			n_base_col_base, r_workspace);
		// a simple one-at-a-time approach
	}

	//CheckIntegrity(true);
	// ...

	return true;
}

bool CUberBlockMatrix::From_Sparse32_Parallel(size_t n_base_row_id, size_t n_base_column_id,
	const cs *p_sparse, bool b_transpose /*= false*/, std::vector<size_t> &r_workspace) // throw(std::bad_alloc)
{
	if(b_transpose)
		return false;
	// for the moment we only handle the not-transpose case

	size_t n_base_row_base = m_block_rows_list[n_base_row_id].n_cumulative_height_sum;
	size_t n_base_col_base = m_block_cols_list[n_base_column_id].n_cumulative_width_sum;
	// get offset of csparse data

	if(p_sparse->m + n_base_row_base > m_n_row_num || p_sparse->n + n_base_col_base > m_n_col_num)
		return false;
	// make sure the sparse data fit inside the existing layout

	if(r_workspace.size() > m_n_row_num)
		r_workspace.clear(); // couldn't be possibly filled from here; clear
	if(r_workspace.size() < m_n_row_num) {
		size_t n_start = r_workspace.size();
		r_workspace.resize(m_n_row_num);
		if(!n_start) {
			size_t n_row_index = 0;
			for(_TyRowConstIter p_row_it = m_block_rows_list.begin(), p_end_it =
			   m_block_rows_list.end(); p_row_it != p_end_it; ++ p_row_it, ++ n_row_index) {
				for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
				   n_start < n; ++ n_start)
					r_workspace[n_start] = n_row_index;
			}
			// fill workspace with row index lookup
		} else {
			_TyRowIter p_row_it;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_start, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_start, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(p_row_it != m_block_rows_list.begin())
				-- p_row_it;
			_ASSERTE((*p_row_it).n_cumulative_height_sum <= n_start);
			_ASSERTE((*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height > n_start);
			// find where to start

			size_t n_row_index = p_row_it - m_block_rows_list.begin();
			for(_TyRowConstIter p_end_it = m_block_rows_list.end();
			   p_row_it != p_end_it; ++ p_row_it, ++ n_row_index) {
				for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
				   n_start < n; ++ n_start)
					r_workspace[n_start] = n_row_index;
			}
			// fill workspace with row index lookup

#ifdef _DEBUG
			{
				size_t n_row_index = 0;
				for(_TyRowConstIter p_row_it = m_block_rows_list.begin(),
				   p_end_it = m_block_rows_list.end(); p_row_it != p_end_it;
				   ++ p_row_it, ++ n_row_index) {
					for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
					   n_start < n; ++ n_start)
						_ASSERTE(r_workspace[n_start] == n_row_index);
				}
			}
#endif // _DEBUG
		}
	}
	// creates a lookup table; complexity O(N) in number of element-rows
	// todo - put this inside a function

#ifdef _OPENMP
	if(m_block_cols_list.size() >= 512) {
		//omp_lock_t t_mutex;
		//omp_init_lock(&t_mutex);
		// initialize the mutex

		_ASSERTE(n_base_column_id <= INT_MAX);
		_ASSERTE(m_block_cols_list.size() <= INT_MAX);
		int n_column_num = int(m_block_cols_list.size());
		#pragma omp parallel for default(shared)
		for(int n_column = int(n_base_column_id); n_column < n_column_num; ++ n_column) {
			TColumn &r_t_col = m_block_cols_list[n_column];
			From_Sparse_ParallelInnerLoop<int>(p_sparse, /*t_mutex,*/ r_t_col, n_base_row_base,
				n_base_col_base, r_workspace);
			// need to be a function in order to make per-thread stack with the static arrays
		}
		// fill the matrix in parallel

		//omp_destroy_lock(&t_mutex);
		// destroy the mutex
	} else {
		for(_TyColumnIter p_col_it = m_block_cols_list.begin() + n_base_column_id,
		   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) {
			TColumn &r_t_col = *p_col_it;
			From_Sparse_InnerLoop<int>(p_sparse, r_t_col, n_base_row_base,
				n_base_col_base, r_workspace);
			// a simple one-at-a-time approach
		}
	}
#else // _OPENMP
	for(_TyColumnIter p_col_it = m_block_cols_list.begin() + n_base_column_id,
	   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) {
		TColumn &r_t_col = *p_col_it;
		From_Sparse_InnerLoop<int>(p_sparse, r_t_col, n_base_row_base,
			n_base_col_base, r_workspace);
		// a simple one-at-a-time approach
	}
#endif // _OPENMP

	//CheckIntegrity(true);
	// ...

	return true;
}

#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64

bool CUberBlockMatrix::From_Sparse_Parallel(size_t n_base_row_id, size_t n_base_column_id,
	const cs *p_sparse, bool b_transpose /*= false*/, std::vector<size_t> &r_workspace) // throw(std::bad_alloc)
{
	if(b_transpose)
		return false;
	// for the moment we only handle the not-transpose case

	size_t n_base_row_base = m_block_rows_list[n_base_row_id].n_cumulative_height_sum;
	size_t n_base_col_base = m_block_cols_list[n_base_column_id].n_cumulative_width_sum;
	// get offset of csparse data

	if(p_sparse->m + n_base_row_base > m_n_row_num || p_sparse->n + n_base_col_base > m_n_col_num)
		return false;
	// make sure the sparse data fit inside the existing layout

	if(r_workspace.size() > m_n_row_num)
		r_workspace.clear(); // couldn't be possibly filled from here; clear
	if(r_workspace.size() < m_n_row_num) {
		size_t n_start = r_workspace.size();
		r_workspace.resize(m_n_row_num);
		if(!n_start) {
			size_t n_row_index = 0;
			for(_TyRowConstIter p_row_it = m_block_rows_list.begin(), p_end_it =
			   m_block_rows_list.end(); p_row_it != p_end_it; ++ p_row_it, ++ n_row_index) {
				for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
				   n_start < n; ++ n_start)
					r_workspace[n_start] = n_row_index;
			}
			// fill workspace with row index lookup
		} else {
			_TyRowIter p_row_it;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_start, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_start, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(p_row_it != m_block_rows_list.begin())
				-- p_row_it;
			_ASSERTE((*p_row_it).n_cumulative_height_sum <= n_start);
			_ASSERTE((*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height > n_start);
			// find where to start

			size_t n_row_index = p_row_it - m_block_rows_list.begin();
			for(_TyRowConstIter p_end_it = m_block_rows_list.end();
			   p_row_it != p_end_it; ++ p_row_it, ++ n_row_index) {
				for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
				   n_start < n; ++ n_start)
					r_workspace[n_start] = n_row_index;
			}
			// fill workspace with row index lookup

#ifdef _DEBUG
			{
				size_t n_row_index = 0;
				for(_TyRowConstIter p_row_it = m_block_rows_list.begin(),
				   p_end_it = m_block_rows_list.end(); p_row_it != p_end_it;
				   ++ p_row_it, ++ n_row_index) {
					for(size_t i = (*p_row_it).n_cumulative_height_sum,
					   n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height; i < n; ++ i)
						_ASSERTE(r_workspace[i] == n_row_index);
				}
			}
#endif // _DEBUG
		}
	}
	// creates a lookup table; complexity O(N) in number of element-rows
	// todo - put this inside a function

#ifdef _OPENMP
	if(m_block_cols_list.size() >= 512) {
		/*omp_lock_t t_mutex;
		omp_init_lock(&t_mutex);*/
		// initialize the mutex

		_ASSERTE(n_base_column_id <= INT_MAX);
		_ASSERTE(m_block_cols_list.size() <= INT_MAX);
		int n_column_num = int(m_block_cols_list.size());
		#pragma omp parallel for default(shared)
		for(int n_column = int(n_base_column_id); n_column < n_column_num; ++ n_column) {
			TColumn &r_t_col = m_block_cols_list[n_column];
			From_Sparse_ParallelInnerLoop<csi>(p_sparse, /*t_mutex,*/ r_t_col, n_base_row_base,
				n_base_col_base, r_workspace);
			// need to be a function in order to make per-thread stack with the static arrays
		}
		// fill the matrix in parallel

		//omp_destroy_lock(&t_mutex);
		// destroy the mutex
	} else {
		for(_TyColumnIter p_col_it = m_block_cols_list.begin() + n_base_column_id,
		   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) {
			TColumn &r_t_col = *p_col_it;
			From_Sparse_InnerLoop<csi>(p_sparse, r_t_col, n_base_row_base,
				n_base_col_base, r_workspace);
			// a simple one-at-a-time approach
		}
	}
#else // _OPENMP
	for(_TyColumnIter p_col_it = m_block_cols_list.begin() + n_base_column_id,
	   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) {
		TColumn &r_t_col = *p_col_it;
		From_Sparse_InnerLoop<csi>(p_sparse, r_t_col, n_base_row_base,
			n_base_col_base, r_workspace);
		// a simple one-at-a-time approach
	}
#endif // _OPENMP

	//CheckIntegrity(true);
	// ...

	return true;
}

bool CUberBlockMatrix::From_Sparse(size_t n_base_row_id, size_t n_base_column_id,
	const cs *p_sparse, bool b_transpose /*= false*/, std::vector<size_t> &r_workspace) // throw(std::bad_alloc)
{
	if(b_transpose)
		return false;
	// for the moment we only handle the not-transpose case

	_ASSERTE(n_base_row_id < m_block_rows_list.size());
	_ASSERTE(n_base_column_id < m_block_cols_list.size());
	// improve user-readability of assertion failures which would otherwise occur in the next two lines

	size_t n_base_row_base = m_block_rows_list[n_base_row_id].n_cumulative_height_sum;
	size_t n_base_col_base = m_block_cols_list[n_base_column_id].n_cumulative_width_sum;
	// get offset of csparse data

	if(p_sparse->m + n_base_row_base > m_n_row_num || p_sparse->n + n_base_col_base > m_n_col_num)
		return false;
	// make sure the sparse data fit inside the existing layout

	if(r_workspace.size() > m_n_row_num)
		r_workspace.clear(); // couldn't be possibly filled from here; clear
	if(r_workspace.size() < m_n_row_num) {
		size_t n_start = r_workspace.size();
		r_workspace.resize(m_n_row_num);
		if(!n_start) {
			size_t n_row_index = 0;
			for(_TyRowConstIter p_row_it = m_block_rows_list.begin(), p_end_it =
			   m_block_rows_list.end(); p_row_it != p_end_it; ++ p_row_it, ++ n_row_index) {
				for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
				   n_start < n; ++ n_start)
					r_workspace[n_start] = n_row_index;
			}
			// fill workspace with row index lookup
		} else {
			_TyRowIter p_row_it;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_start, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_start, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(p_row_it != m_block_rows_list.begin())
				-- p_row_it;
			_ASSERTE((*p_row_it).n_cumulative_height_sum <= n_start);
			_ASSERTE((*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height > n_start);
			// find where to start

			size_t n_row_index = p_row_it - m_block_rows_list.begin();
			for(_TyRowConstIter p_end_it = m_block_rows_list.end();
			   p_row_it != p_end_it; ++ p_row_it, ++ n_row_index) {
				for(size_t n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height;
				   n_start < n; ++ n_start)
					r_workspace[n_start] = n_row_index;
			}
			// fill workspace with row index lookup

#ifdef _DEBUG
			{
				size_t n_row_index = 0;
				for(_TyRowConstIter p_row_it = m_block_rows_list.begin(),
				   p_end_it = m_block_rows_list.end(); p_row_it != p_end_it;
				   ++ p_row_it, ++ n_row_index) {
					for(size_t i = (*p_row_it).n_cumulative_height_sum,
					   n = (*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height; i < n; ++ i)
						_ASSERTE(r_workspace[i] == n_row_index);
				}
			}
#endif // _DEBUG
		}
	}
	// creates a lookup table; complexity O(N) in number of element-rows
	// todo - put this inside a function

	for(_TyColumnIter p_col_it = m_block_cols_list.begin() + n_base_column_id,
	   p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) {
		TColumn &r_t_col = *p_col_it;
		From_Sparse_InnerLoop<csi>(p_sparse, r_t_col, n_base_row_base,
			n_base_col_base, r_workspace);
		// a simple one-at-a-time approach
	}

#if 0 // the code below has high runtime complexity
	if(b_transpose) {
		_TyRowIter p_row_it;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
		p_row_it = std::upper_bound(m_block_rows_list.begin(),
			m_block_rows_list.end(), n_base_row, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		p_row_it = std::upper_bound(m_block_rows_list.begin(),
			m_block_rows_list.end(), n_base_row, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		if(p_row_it != m_block_rows_list.begin())
			-- p_row_it;

		for(csi i = 0, n = p_sparse->n; i < n; ++ i) { // todo this probably could be done in parallel by looping over row blocks and each thread working with just a single column (need to apply the transpose trick)
			size_t n_row = n_base_row + i;
			// column that is about to be filled

			if(p_cs_col_ptr[i] == p_cs_col_ptr[i + 1])
				continue;
			// skip empty columns

			_ASSERTE((*p_row_it).n_cumulative_height_sum <= n_row); // current column should always be to the left
			if((*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height <= n_row) { // we need to move to the next column
				if(p_row_it != m_block_rows_list.end()) {
					++ p_row_it; // try the immediate neighbour
					_ASSERTE((*p_row_it).n_cumulative_height_sum <= n_row); // current column should always be to the left
					if((*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height <= n_row) { // not enough, need to seek
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						p_row_it = std::upper_bound(p_row_it,
							m_block_rows_list.end(), n_base_column, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						p_row_it = std::upper_bound(p_row_it,
							m_block_rows_list.end(), n_base_column, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						_ASSERTE(p_row_it != m_block_rows_list.begin());
						-- p_row_it; // for sure doesn't equal begin()
					}
				} else
					throw std::runtime_error("From_Sparse() inserted row out of range"); // block structure must exist
			}
			TRow &r_t_row = *p_row_it;
			_ASSERTE(r_t_row.n_cumulative_height_sum <= n_row); // current column should always be to the left
			_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height > n_row); // make sure the chosen column is the right one
			// resolve row

			size_t n_first_column = n_base_column + p_cs_row_index[p_cs_col_ptr[i]];

			_TyColumnIter p_col_it;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			p_col_it = std::upper_bound(m_block_cols_list.begin(),
				m_block_cols_list.end(), n_first_column, CFindLowerColumn());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			p_col_it = std::upper_bound(m_block_cols_list.begin(),
				m_block_cols_list.end(), n_first_column, _FindLowerColumn);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(p_col_it != m_block_cols_list.begin())
				-- p_col_it;

			for(csi j = p_cs_col_ptr[i], n_row_end = p_cs_col_ptr[i + 1]; j < n_row_end; ++ j) {
				size_t n_column = n_base_row + p_cs_row_index[j];
				double f_value = p_cs_x[j];
				// get triplet form record

				_ASSERTE((*p_col_it).n_cumulative_width_sum <= n_column); // current column should always be to the left
				if((*p_col_it).n_cumulative_width_sum + (*p_col_it).n_width <= n_column) { // we need to move to the next column
					if(p_col_it != m_block_cols_list.end()) {
						++ p_col_it; // try the immediate neighbour
						_ASSERTE((*p_col_it).n_cumulative_width_sum <= n_column); // current column should always be to the left
						if((*p_col_it).n_cumulative_width_sum + (*p_col_it).n_width <= n_column) { // not enough, need to seek
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
							p_col_it = std::upper_bound(p_col_it,
								m_block_cols_list.end(), n_column, CFindLowerColumn());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
							p_col_it = std::upper_bound(p_col_it,
								m_block_cols_list.end(), n_column, _FindLowerColumn);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
							_ASSERTE(p_col_it != m_block_cols_list.begin());
							-- p_col_it; // for sure doesn't equal begin()
						}
					} else
						throw std::runtime_error("From_Sparse() inserted column out of range"); // block structure must exist
				}
				TColumn &r_t_col = *p_col_it;
				size_t n_column_id = p_col_it - m_block_cols_list.begin();
				_ASSERTE(r_t_col.n_cumulative_width_sum <= n_column); // current column should always be to the left
				_ASSERTE(r_t_col.n_cumulative_width_sum + r_t_col.n_width > n_column); // make sure the chosen column is the right one
				// resolve column

				bool b_new_block;
				double *p_block_data = p_AllocBlockData(p_row_it - m_block_rows_list.begin(),
					n_column_id, r_t_row.n_height, r_t_col.n_width, b_new_block);
				if(b_new_block)
					memset(p_block_data, 0, r_t_row.n_height * r_t_col.n_width * sizeof(double));
				// get initialized block // t_odo - cache the last block iterator, seek from that // can't

				p_block_data[(n_column - r_t_col.n_cumulative_width_sum) * r_t_row.n_height +
					n_row - r_t_row.n_cumulative_height_sum] = f_value;
				// write data in the block
			}
		}
	} else {
		_TyColumnIter p_col_it;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
		p_col_it = std::upper_bound(m_block_cols_list.begin(),
			m_block_cols_list.end(), n_base_column, CFindLowerColumn());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		p_col_it = std::upper_bound(m_block_cols_list.begin(),
			m_block_cols_list.end(), n_base_column, _FindLowerColumn);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		if(p_col_it != m_block_cols_list.begin())
			-- p_col_it;

		for(csi i = 0, n = p_sparse->n; i < n; ++ i) {
			size_t n_column = n_base_column + i;
			// column that is about to be filled

			if(p_cs_col_ptr[i] == p_cs_col_ptr[i + 1])
				continue;
			// skip empty columns

			_ASSERTE((*p_col_it).n_cumulative_width_sum <= n_column); // current column should always be to the left
			if((*p_col_it).n_cumulative_width_sum + (*p_col_it).n_width <= n_column) { // we need to move to the next column
				if(p_col_it != m_block_cols_list.end()) {
					++ p_col_it; // try the immediate neighbour
					_ASSERTE((*p_col_it).n_cumulative_width_sum <= n_column); // current column should always be to the left
					if((*p_col_it).n_cumulative_width_sum + (*p_col_it).n_width <= n_column) { // not enough, need to seek
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						p_col_it = std::upper_bound(p_col_it,
							m_block_cols_list.end(), n_base_column, CFindLowerColumn());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						p_col_it = std::upper_bound(p_col_it,
							m_block_cols_list.end(), n_base_column, _FindLowerColumn);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						_ASSERTE(p_col_it != m_block_cols_list.begin());
						-- p_col_it; // for sure doesn't equal begin()
					}
				} else
					throw std::runtime_error("From_Sparse() inserted column out of range"); // block structure must exist
			}
			TColumn &r_t_col = *p_col_it;
			size_t n_column_id = p_col_it - m_block_cols_list.begin();
			_ASSERTE(r_t_col.n_cumulative_width_sum <= n_column); // current column should always be to the left
			_ASSERTE(r_t_col.n_cumulative_width_sum + r_t_col.n_width > n_column); // make sure the chosen column is the right one
			// resolve column

			size_t n_first_row = n_base_row + p_cs_row_index[p_cs_col_ptr[i]];

			_TyRowIter p_row_it;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_first_row, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			p_row_it = std::upper_bound(m_block_rows_list.begin(),
				m_block_rows_list.end(), n_first_row, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(p_row_it != m_block_rows_list.begin())
				-- p_row_it;

			for(csi j = p_cs_col_ptr[i], n_row_end = p_cs_col_ptr[i + 1]; j < n_row_end; ++ j) {
				size_t n_row = n_base_row + p_cs_row_index[j];
				double f_value = p_cs_x[j];
				// get triplet form record

				_ASSERTE((*p_row_it).n_cumulative_height_sum <= n_row); // current column should always be to the left
				if((*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height <= n_row) { // we need to move to the next column
					if(p_row_it != m_block_rows_list.end()) {
						++ p_row_it; // try the immediate neighbour
						_ASSERTE((*p_row_it).n_cumulative_height_sum <= n_row); // current column should always be to the left
						if((*p_row_it).n_cumulative_height_sum + (*p_row_it).n_height <= n_row) { // not enough, need to seek
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
							p_row_it = std::upper_bound(p_row_it,
								m_block_rows_list.end(), n_base_column, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
							p_row_it = std::upper_bound(p_row_it,
								m_block_rows_list.end(), n_base_column, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
							_ASSERTE(p_row_it != m_block_rows_list.begin());
							-- p_row_it; // for sure doesn't equal begin()
						}
					} else
						throw std::runtime_error("From_Sparse() inserted row out of range"); // block structure must exist
				}
				TRow &r_t_row = *p_row_it;
				_ASSERTE(r_t_row.n_cumulative_height_sum <= n_row); // current column should always be to the left
				_ASSERTE(r_t_row.n_cumulative_height_sum + r_t_row.n_height > n_row); // make sure the chosen column is the right one
				// resolve row

				bool b_new_block;
				double *p_block_data = p_AllocBlockData(p_row_it - m_block_rows_list.begin(),
					n_column_id, r_t_row.n_height, r_t_col.n_width, b_new_block);
				if(b_new_block)
					memset(p_block_data, 0, r_t_row.n_height * r_t_col.n_width * sizeof(double));
				// get initialized block // todo - cache the last block iterator, seek from that

				p_block_data[(n_column - r_t_col.n_cumulative_width_sum) * r_t_row.n_height +
					n_row - r_t_row.n_cumulative_height_sum] = f_value;
				// write data in the block
			}
		}
	}
#endif // 0

	//CheckIntegrity(true);

	return true;
}

cs *CUberBlockMatrix::p_Convert_to_Sparse(cs *p_alloc /*= 0*/) const
{
	CheckIntegrity(true);

	if(!p_alloc) {
		if(!(p_alloc = cs_spalloc(m_n_row_num, m_n_col_num,
		   m_n_ref_elem_num + m_data_pool.size(), 1, 0))) // this will not realloc (will alloc a bit more, in fact, but should not be too much more)
			return 0;
	} else {
		if((p_alloc->nzmax < 0 || unsigned(p_alloc->nzmax) < m_n_ref_elem_num +
		   m_data_pool.size()) && !cs_sprealloc(p_alloc, std::max(size_t(p_alloc->nzmax * 2),
		   m_n_ref_elem_num + m_data_pool.size()))) // this will not realloc (will alloc a bit more, in fact, but should not be too much more)
			return 0;
		if((p_alloc->n < 0 || size_t(p_alloc->n) < m_n_col_num) &&
		   !(p_alloc->p = (csi*)realloc(p_alloc->p, (m_n_col_num + 1) * sizeof(csi)))) // cs_sprealloc doesn't realloc column array
			return 0;
		p_alloc->m = m_n_row_num;
		p_alloc->n = m_n_col_num;
	}
	_ASSERTE(p_alloc->nzmax > 0);
	_ASSERTE(unsigned(p_alloc->nzmax) >= m_n_ref_elem_num + m_data_pool.size());
	// prepare dest matrix

	csi *p_column = p_alloc->p, *p_index = p_alloc->i;
	double *p_value = p_alloc->x;
	// access values

	size_t n_nonzero_element_num = 0;
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;
		const size_t n_width = r_t_col.n_width;
		_ASSERTE(n_width > 0);
		// get column

		size_t n_column_row_num = n_nonzero_element_num;
		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
		   p_block_it != p_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const double *p_data = r_t_block.second;
			// for each column, for each block ...

			n_nonzero_element_num += r_t_row.n_height;
			_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
			// will write r_t_row.n_height values

			for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
			   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
				*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
				*p_index = n_row;
			}
			// write values to the arrays
		}
		n_column_row_num = n_nonzero_element_num - n_column_row_num; // no need to keep two counters in the loop
		// first loop collects number of rows in the column

		_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
		_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

		for(size_t i = 0; i < n_width; ++ i, ++ p_column)
			*p_column = n_nonzero_element_num - n_column_row_num + i * n_column_row_num;
		// fill the column structures (works on empty columns as well)

		if(n_width > 1) {
			p_value -= n_column_row_num;
			p_index -= n_column_row_num;
			// make those zero-based indexable

			_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num - n_column_row_num);
			_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num - n_column_row_num);

			n_nonzero_element_num += (n_width - 1) * n_column_row_num;
			_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
			// will write n_column_row_num values, n_width - 1 times (once per each of the rest of element-columns)

			for(_TyBlockConstIter p_block_it =
			   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
			   p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				const TRow &r_t_row = m_block_rows_list[r_t_block.first];
				const double *p_data = r_t_block.second;
				// for each column, for each block ...

				for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
				   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
					for(size_t x = 1; x < n_width; ++ x) {
						p_value[x * n_column_row_num] = p_data[x * r_t_row.n_height + i];
						p_index[x * n_column_row_num] = n_row;
					}
				}
				// write values to the arrays
			}
			// process element-columns 1 to n_width - 1

			_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num -
				n_column_row_num * (n_width - 1));
			_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num -
				n_column_row_num * (n_width - 1));

			p_value += n_column_row_num * (n_width - 1);
			p_index += n_column_row_num * (n_width - 1);
			// shift after the last value used

			_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
			_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);
		}
		// second loop processes all the element-columns at once
	}
	// writes out the matrix in compressed column form directly without triplet form and cs_compress()

	_ASSERTE(p_column == p_alloc->p + p_alloc->n); // make sure it's pointing to the last column
	_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num);
	*p_column = n_nonzero_element_num;
	_ASSERTE(p_index == p_alloc->i + n_nonzero_element_num);
	_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

	return p_alloc;
}

cs *CUberBlockMatrix::p_BlockStructure_to_Sparse(cs *p_alloc /*= 0*/) const
{
	CheckIntegrity();

	if(!p_alloc) {
		if(!(p_alloc = cs_spalloc(m_block_rows_list.size(),
		   m_block_cols_list.size(), n_Block_Num(), 0, 0))) // this will not realloc
			return 0;
	} else {
		if(p_alloc->nzmax < 1 && !cs_sprealloc(p_alloc,
		   std::max(size_t(p_alloc->nzmax * 2), n_Block_Num()))) // this will not realloc
			return 0;
		if((p_alloc->n < 0 || size_t(p_alloc->n) < m_block_cols_list.size()) &&
		   !(p_alloc->p = (csi*)realloc(p_alloc->p, (m_block_cols_list.size() + 1) * sizeof(csi)))) // cs_sprealloc doesn't realloc column array
			return 0;
		p_alloc->m = m_block_rows_list.size();
		p_alloc->n = m_block_cols_list.size();
	}
	_ASSERTE(p_alloc->nzmax > 0);
	// prepare dest matrix (note it might realloc if allocated outside)

	csi *p_column = p_alloc->p, *p_index = p_alloc->i;
	// access values

	size_t n_nonzero_block_num = 0, n_column_index = 0;
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it;
	   ++ p_col_it, ++ n_column_index) {
		const TColumn &r_t_col = *p_col_it;
		// get column

		size_t n_block_num = r_t_col.block_list.size();
		if(n_nonzero_block_num + n_block_num > size_t(p_alloc->nzmax)) {
			_ASSERTE(n_nonzero_block_num == p_index - p_alloc->i);
			_ASSERTE(n_column_index == p_column - p_alloc->p);
			// make sure the pointers are on track

			if(!cs_sprealloc(p_alloc, std::max(size_t(p_alloc->nzmax * 2),
			   n_nonzero_block_num + n_block_num)))
				return 0;
			_ASSERTE(n_nonzero_block_num + n_block_num <= size_t(p_alloc->nzmax)); // will fit now
			// realloc

			p_index = p_alloc->i + n_nonzero_block_num;
			//p_column = p_alloc->p + n_column_index; // does not even change
			_ASSERTE(p_column == p_alloc->p + n_column_index);
			// reconstruct pointers after realloc, implement cholsol
		}
		// make sure all the blocks will fit

		*p_column = n_nonzero_block_num;
		++ p_column;
		// write column pointer

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
		   p_block_it != p_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			*p_index = r_t_block.first; // row
			++ p_index;
		}
		// go through all the blocks, fill the row numbers

		n_nonzero_block_num += n_block_num;
	}

	double *p_value = p_alloc->x;
	if(p_value) {
		for(const double *p_val_end = p_value + n_nonzero_block_num; p_value != p_val_end; ++ p_value)
			*p_value = 1;
	}
	// someone may want values, fill them all at once to avoid branching inside the above loop(s)

	_ASSERTE(p_column == p_alloc->p + p_alloc->n); // make sure it points to the end
	*p_column = n_nonzero_block_num;
	// finalize the matrix layout

	return p_alloc;
}

cs *CUberBlockMatrix::p_BlockStructure_to_Sparse_Apron(size_t n_min_block,
	size_t n_diag_block_num, cs *p_alloc /*= 0*/) const // todo - write a forced 32-bit version
{
	CheckIntegrity();

	_ASSERTE(m_block_rows_list.size() >= n_min_block);
	_ASSERTE(m_block_cols_list.size() >= n_min_block); // the cut must not be larger than the matrix
	_ASSERTE(n_min_block <= n_diag_block_num); // n_diag_block_num must be after n_min_block
	_ASSERTE(m_block_rows_list.size() >= n_diag_block_num);
	_ASSERTE(m_block_cols_list.size() >= n_diag_block_num); // as well
	if(!p_alloc) {
		if(!(p_alloc = cs_spalloc(m_block_rows_list.size() - n_min_block,
		   m_block_cols_list.size() - n_min_block, n_Block_Num(), 0, 0))) // this will not realloc
			return 0;
	} else {
		if(p_alloc->nzmax < 1 && !cs_sprealloc(p_alloc,
		   std::max(size_t(p_alloc->nzmax * 2), n_Block_Num()))) // this will not realloc
			return 0;
		if((p_alloc->n < 0 || size_t(p_alloc->n) < m_block_cols_list.size()) &&
		   !(p_alloc->p = (csi*)realloc(p_alloc->p, (m_block_cols_list.size() + 1) * sizeof(csi)))) // cs_sprealloc doesn't realloc column array
			return 0;
		p_alloc->m = m_block_rows_list.size() - n_min_block;
		p_alloc->n = m_block_cols_list.size() - n_min_block;
	}
	_ASSERTE(p_alloc->nzmax > 0);
	// prepare dest matrix (note it might realloc if allocated outside)

	csi *p_column = p_alloc->p, *p_index = p_alloc->i;
	// access values

	size_t n_nonzero_block_num = 0, n_column_index = n_min_block;
	for(; n_column_index < n_diag_block_num; ++ n_column_index) {
		const size_t n_block_num = 1;
		if(n_nonzero_block_num + n_block_num > size_t(p_alloc->nzmax)) {
			_ASSERTE(n_nonzero_block_num == p_index - p_alloc->i);
			_ASSERTE(n_column_index - n_min_block == p_column - p_alloc->p);
			// make sure the pointers are on track

			if(!cs_sprealloc(p_alloc, std::max(size_t(p_alloc->nzmax * 2),
			   n_nonzero_block_num + n_block_num)))
				return 0;
			_ASSERTE(n_nonzero_block_num + n_block_num <= size_t(p_alloc->nzmax)); // will fit now
			// realloc

			p_index = p_alloc->i + n_nonzero_block_num;
			//p_column = p_alloc->p + (n_column_index - n_min_block); // does not even change
			_ASSERTE(p_column == p_alloc->p + (n_column_index - n_min_block));
			// reconstruct pointers after realloc, implement cholsol
		}
		// make sure all the blocks will fit

		*p_column = n_nonzero_block_num;
		++ p_column;
		// write column pointer

		{
			*p_index = n_column_index - n_min_block; // row; no underflow here
			++ p_index;
		}
		// go through all the blocks, fill only the diagonal entries

		n_nonzero_block_num += n_block_num;
	}
	// fill-out the mini-skirt

	if(n_min_block > 0) {
		for(_TyColumnConstIter p_col_it = m_block_cols_list.begin() + n_diag_block_num,
		   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it;
		   ++ p_col_it, ++ n_column_index) {
			const TColumn &r_t_col = *p_col_it;
			// get column

			size_t n_block_num = r_t_col.block_list.size();
			if(n_nonzero_block_num + n_block_num > size_t(p_alloc->nzmax)) {
				_ASSERTE(n_nonzero_block_num == p_index - p_alloc->i);
				_ASSERTE(n_column_index - n_min_block == p_column - p_alloc->p);
				// make sure the pointers are on track

				if(!cs_sprealloc(p_alloc, std::max(size_t(p_alloc->nzmax * 2),
				   n_nonzero_block_num + n_block_num)))
					return 0;
				_ASSERTE(n_nonzero_block_num + n_block_num <= size_t(p_alloc->nzmax)); // will fit now
				// realloc

				p_index = p_alloc->i + n_nonzero_block_num;
				//p_column = p_alloc->p + (n_column_index - n_min_block); // does not even change
				_ASSERTE(p_column == p_alloc->p + (n_column_index - n_min_block));
				// reconstruct pointers after realloc, implement cholsol
			}
			// make sure all the blocks will fit

			*p_column = n_nonzero_block_num;
			++ p_column;
			// write column pointer

			for(_TyBlockConstIter p_block_it =
			   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
			   p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				_ASSERTE(r_t_block.first >= n_min_block); // could underflow if we did not skip carefully
				*p_index = r_t_block.first - n_min_block; // row
				++ p_index;
			}
			// go through all the blocks, fill the row numbers

			n_nonzero_block_num += n_block_num;
		}
	} else {
		for(_TyColumnConstIter p_col_it = m_block_cols_list.begin() + n_diag_block_num,
		   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it;
		   ++ p_col_it, ++ n_column_index) {
			const TColumn &r_t_col = *p_col_it;
			// get column

			size_t n_block_num = r_t_col.block_list.size();
			if(n_nonzero_block_num + n_block_num > size_t(p_alloc->nzmax)) {
				_ASSERTE(n_nonzero_block_num == p_index - p_alloc->i);
				_ASSERTE(n_column_index - n_min_block == p_column - p_alloc->p);
				// make sure the pointers are on track

				if(!cs_sprealloc(p_alloc, std::max(size_t(p_alloc->nzmax * 2),
				   n_nonzero_block_num + n_block_num)))
					return 0;
				_ASSERTE(n_nonzero_block_num + n_block_num <= size_t(p_alloc->nzmax)); // will fit now
				// realloc

				p_index = p_alloc->i + n_nonzero_block_num;
				//p_column = p_alloc->p + (n_column_index - n_min_block); // does not even change
				_ASSERTE(p_column == p_alloc->p + (n_column_index - n_min_block));
				// reconstruct pointers after realloc, implement cholsol
			}
			// make sure all the blocks will fit

			*p_column = n_nonzero_block_num;
			++ p_column;
			// write column pointer

			for(_TyBlockConstIter p_block_it =
			   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
			   p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				*p_index = r_t_block.first; // row
				++ p_index;
			}
			// go through all the blocks, fill the row numbers

			n_nonzero_block_num += n_block_num;
		}
	}
	// fill-out the actual blocks

	double *p_value = p_alloc->x;
	if(p_value) {
		for(const double *p_val_end = p_value + n_nonzero_block_num; p_value != p_val_end; ++ p_value)
			*p_value = 1;
	}
	// someone may want values, fill them all at once to avoid branching inside the above loop(s)

	_ASSERTE(p_column == p_alloc->p + p_alloc->n); // make sure it points to the end
	*p_column = n_nonzero_block_num;
	// finalize the matrix layout

	return p_alloc;
}

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)

cs *CUberBlockMatrix::p_BlockStructure_to_Sparse32(cs *p_alloc /*= 0*/) const
{
	CheckIntegrity();

	if(!p_alloc) {
		if(!(p_alloc = cs_spalloc32(m_block_rows_list.size(),
		   m_block_cols_list.size(), n_Block_Num(), 0, 0))) // this will not realloc
			return 0;
	} else {
		if(p_alloc->nzmax < 1) {
			if(!cs_sprealloc32(p_alloc, std::max(size_t(p_alloc->nzmax * 2), n_Block_Num()))) // this will not realloc
				return 0;
		}
		if((p_alloc->n < 0 || size_t(p_alloc->n) < m_block_cols_list.size()) &&
		   !(p_alloc->p = (csi*)realloc(p_alloc->p, (m_block_cols_list.size() + 1) * sizeof(int)))) // cs_sprealloc doesn't realloc column array
			return 0;
		p_alloc->m = m_block_rows_list.size();
		p_alloc->n = m_block_cols_list.size();
	}
	_ASSERTE(p_alloc->nzmax > 0);
	// prepare dest matrix (note it might realloc if allocated outside)

	int *p_column = (int*)p_alloc->p, *p_index = (int*)p_alloc->i;
	// access values

	size_t n_nonzero_block_num = 0, n_column_index = 0;
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it;
	   ++ p_col_it, ++ n_column_index) {
		const TColumn &r_t_col = *p_col_it;
		// get column

		size_t n_block_num = r_t_col.block_list.size();
		if(n_nonzero_block_num + n_block_num > size_t(p_alloc->nzmax)) {
			_ASSERTE(n_nonzero_block_num == p_index - ((int*)p_alloc->i));
			_ASSERTE(n_column_index == p_column - ((int*)p_alloc->p));
			// make sure the pointers are on track

			if(!cs_sprealloc32(p_alloc, (std::max(size_t(p_alloc->nzmax * 2),
			   n_nonzero_block_num + n_block_num))))
				return 0;
			_ASSERTE(n_nonzero_block_num + n_block_num <= size_t(p_alloc->nzmax)); // will fit now
			// realloc

			p_index = ((int*)p_alloc->i) + n_nonzero_block_num;
			//p_column = ((int*)p_alloc->p) + n_column_index;
			_ASSERTE(p_column == ((int*)p_alloc->p) + n_column_index);
			// reconstruct pointers after realloc, implement cholsol
		}
		// make sure all the blocks will fit

		*p_column = int(n_nonzero_block_num);
		++ p_column;
		// write column pointer

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
		   p_block_it != p_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			*p_index = int(r_t_block.first); // row
			++ p_index;
		}
		// go through all the blocks, fill the row numbers

		n_nonzero_block_num += n_block_num;
	}

	double *p_value = p_alloc->x;
	if(p_value) {
		for(const double *p_val_end = p_value + n_nonzero_block_num; p_value != p_val_end; ++ p_value)
			*p_value = 1;
	}
	// someone may want values, fill them all at once to avoid branching inside the above loop(s)

	_ASSERTE(p_column == ((int*)p_alloc->p) + p_alloc->n); // make sure it points to the end
	*p_column = int(n_nonzero_block_num);
	// finalize the matrix layout

	return p_alloc;
}

cs *CUberBlockMatrix::p_Convert_to_Sparse32(cs *p_alloc /*= 0*/) const
{
	CheckIntegrity(true);

	if(!p_alloc) {
		if(!(p_alloc = cs_spalloc32(m_n_row_num, m_n_col_num,
		   m_n_ref_elem_num + m_data_pool.size(), 1, 0))) // this will not realloc (will alloc a bit more, in fact, but should not be too much more)
			return 0;
	} else {
		if((p_alloc->nzmax < 0 || unsigned(p_alloc->nzmax) < m_n_ref_elem_num +
		   m_data_pool.size()) && !cs_sprealloc32(p_alloc,
		   (std::max(size_t(p_alloc->nzmax * 2), m_n_ref_elem_num + m_data_pool.size())))) // this will not realloc (will alloc a bit more, in fact, but should not be too much more)
			return 0;
		if((p_alloc->n < 0 || size_t(p_alloc->n) < m_n_col_num) &&
		   !(p_alloc->p = (csi*)realloc(p_alloc->p, (m_n_col_num + 1) * sizeof(int)))) // cs_sprealloc doesn't realloc column array
			return 0;
		p_alloc->m = m_n_row_num;
		p_alloc->n = m_n_col_num;
	}
	_ASSERTE(p_alloc->nzmax > 0);
	_ASSERTE(unsigned(p_alloc->nzmax) >= m_n_ref_elem_num + m_data_pool.size());
	// prepare dest matrix

	int *p_column = (int*)p_alloc->p, *p_index = (int*)p_alloc->i;
	double *p_value = p_alloc->x;
	// access values

	size_t n_nonzero_element_num = 0;
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;
		const size_t n_width = r_t_col.n_width;
		_ASSERTE(n_width > 0);
		// get column

		size_t n_column_row_num = n_nonzero_element_num;
		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
		   p_block_it != p_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const double *p_data = r_t_block.second;
			// for each column, for each block ...

			n_nonzero_element_num += r_t_row.n_height;
			_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
			// will write r_t_row.n_height values

			for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
			   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
				*p_value = p_data[i]; // x * r_t_row.n_height + i, x = 0
				*p_index = int(n_row);
			}
			// write values to the arrays
		}
		n_column_row_num = n_nonzero_element_num - n_column_row_num; // no need to keep two counters in the loop
		// first loop collects number of rows in the column

		_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
		_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

		for(size_t i = 0; i < n_width; ++ i, ++ p_column)
			*p_column = int(n_nonzero_element_num - n_column_row_num + i * n_column_row_num);
		// fill the column structures (works on empty columns as well)

		if(n_width > 1) {
			p_value -= n_column_row_num;
			p_index -= n_column_row_num;
			// make those zero-based indexable

			_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num - n_column_row_num);
			_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num - n_column_row_num);

			n_nonzero_element_num += (n_width - 1) * n_column_row_num;
			_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num); // make sure we don't go out of bounds
			// will write n_column_row_num values, n_width - 1 times (once per each of the rest of element-columns)

			for(_TyBlockConstIter p_block_it =
			   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
			   p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				const TRow &r_t_row = m_block_rows_list[r_t_block.first];
				const double *p_data = r_t_block.second;
				// for each column, for each block ...

				for(size_t i = 0, n_row = r_t_row.n_cumulative_height_sum;
				   i < r_t_row.n_height; ++ i, ++ n_row, ++ p_value, ++ p_index) {
					for(size_t x = 1; x < n_width; ++ x) {
						p_value[x * n_column_row_num] = p_data[x * r_t_row.n_height + i];
						p_index[x * n_column_row_num] = int(n_row);
					}
				}
				// write values to the arrays
			}
			// process element-columns 1 to n_width - 1

			_ASSERTE(p_index == ((int*)p_alloc->i) +
				n_nonzero_element_num - n_column_row_num * (n_width - 1));
			_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num -
				n_column_row_num * (n_width - 1));

			p_value += n_column_row_num * (n_width - 1);
			p_index += n_column_row_num * (n_width - 1);
			// shift after the last value used

			_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
			_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);
		}
		// second loop processes all the element-columns at once
	}
	// writes out the matrix in compressed column form directly without triplet form and cs_compress()

	_ASSERTE(p_column == ((int*)p_alloc->p) + p_alloc->n); // make sure it's pointing to the last column
	_ASSERTE(unsigned(p_alloc->nzmax) >= n_nonzero_element_num);
	*p_column = int(n_nonzero_element_num);
	_ASSERTE(p_index == ((int*)p_alloc->i) + n_nonzero_element_num);
	_ASSERTE(p_value == p_alloc->x + n_nonzero_element_num);

	return p_alloc;
}

#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64

cs *CUberBlockMatrix::p_Convert_to_Sparse_Debug() const
{
	CheckIntegrity(true);

	cs *p_trip;
	if(!(p_trip = cs_spalloc(m_n_row_num, m_n_col_num, 1024, 1, 1)))
		return 0;

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
					if(/*f_elem != 0 &&*/ !cs_entry(p_trip, n_y + y, n_x + x, f_elem)) {
						cs_spfree(p_trip);
						return 0;
					}
				}
			}
			// iterate through data, sparse fill ...
		}
	}

	cs *p_csc;
	if(!(p_csc = cs_compress(p_trip))) {
		cs_spfree(p_trip);
		return 0;
	}

	cs_spfree(p_trip);

	return p_csc;
}

cs *CUberBlockMatrix::p_Convert_to_Sparse_UpperDiagonal_Debug() const
{
	CheckIntegrity(true);

	cs *p_trip;
	if(!(p_trip = cs_spalloc(m_n_row_num, m_n_col_num, 1024, 1, 1)))
		return 0;

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
					if(n_y + y <= n_x + x && /*f_elem != 0 &&*/
					   !cs_entry(p_trip, n_y + y, n_x + x, f_elem)) { // row <= column
						cs_spfree(p_trip);
						return 0;
					}
				}
			}
			// iterate through data, sparse fill ...
		}
	}

	cs *p_csc;
	if(!(p_csc = cs_compress(p_trip))) {
		cs_spfree(p_trip);
		return 0;
	}

	cs_spfree(p_trip);

	return p_csc;
}

#endif // __UBER_BLOCK_MATRIX_HAVE_CSPARSE

void CUberBlockMatrix::TransposeTo(CUberBlockMatrix &r_dest) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	_ASSERTE(&r_dest != this); // can't transpose to itself

	r_dest.Clear();
	// erase dest ...

	r_dest.m_n_row_num = m_n_col_num;
	r_dest.m_n_col_num = m_n_row_num;
	// set matrix size

	r_dest.m_block_cols_list.resize(m_block_rows_list.size());
	/*std::for_each(r_dest.m_block_cols_list.begin(), r_dest.m_block_cols_list.end(),
		CRowCumsumToColumnCumsum(m_block_rows_list.begin()));*/ // the loop below, but using std::for_each
	std::transform(m_block_rows_list.begin(), m_block_rows_list.end(),
		r_dest.m_block_cols_list.begin(), t_RowCumsumToColumnCumsum);
	/*for(size_t i = 0, n = m_block_rows_list.size(); i < n; ++ i) {
		r_dest.m_block_cols_list[i].n_cumulative_width_sum = m_block_rows_list[i].n_cumulative_height_sum;
		r_dest.m_block_cols_list[i].n_width = m_block_rows_list[i].n_height;
		// t_odo - prealloc lists? // probably not - it would take cumsum O(blocks) = O(rows * avgColFil) and it would save O(rows * log(avgColFil)) // todo - try it anyway and see what happens
	}*/
	// copy the columns structure

	r_dest.m_block_rows_list.resize(m_block_cols_list.size());
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const TColumn &r_t_col = m_block_cols_list[i];
		r_dest.m_block_rows_list[i].n_height = r_t_col.n_width;
		r_dest.m_block_rows_list[i].n_cumulative_height_sum = r_t_col.n_cumulative_width_sum;
		// copy row structure

		size_t n_row_idx = i;
		size_t n_block_height = r_t_col.n_width;
		for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
			size_t n_col_idx = r_t_col.block_list[j].first;
			const double *p_src_data = r_t_col.block_list[j].second;
			size_t n_block_width = m_block_rows_list[n_col_idx].n_height;
			// get block position, size and data

			double *p_data = r_dest.p_Get_DenseStorage(n_block_width * n_block_height);
			// get dense storage for the block

			/*for(size_t x = 0; x < n_block_width; ++ x)
				for(size_t y = 0; y < n_block_height; ++ y)
					p_data[y + n_block_height * x] = p_src_data[x + n_block_width * y];*/ // f_ixme - does the src / dest ordering matter? is it interchangable? // t_odo - test transposing a matrix // this seems to be correct // now it's correct
			// this works, but it doesn't if the indexing is reversed; not very trustworthy code

			_TyMatrixXdRef src((double*)p_src_data, n_block_width, n_block_height),
				dest(p_data, n_block_height, n_block_width);
			dest = src.transpose();
			// copy/transpose the block using Eigen // todo - run tests, see if it's slower

			r_dest.m_block_cols_list[n_col_idx].block_list.push_back(
				TColumn::TBlockEntry(n_row_idx, p_data));
			// add the block to the destinat list of blocks
		}
	}
	// copy the rows structure, scatter the block to the column structures

	// t_odo - transpose the matrix to r_dest

	//r_dest.CheckIntegrity();
	// make sure that the destination matrix was assembled correctly
}

void CUberBlockMatrix::PreMultiplyWithSelfTransposeTo(CUberBlockMatrix &r_dest,
	bool b_upper_diagonal_only /*= false*/) const // throw(std::bad_alloc)
{
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
		std::transform(r_row_list_B.begin(), r_row_list_B.end(),
			col_list_A.begin(), t_RowCumsumToColumnCumsum);
		for(size_t i = 0, n = r_col_list_B.size(); i < n; ++ i) {
			const TColumn &r_t_col = r_col_list_B[i];
			for(_TyBlockConstIter p_block_it =
			   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
			   p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				col_list_A[r_t_block.first].block_list.push_back(TColumn::TBlockEntry(i,
					r_t_block.second));
			}
		}
	}
	// calculate fast transpose of this (no rows alloc, no matrix structure,
	// no block transpose, no extra storage, O(n) in number of blocks) // @todo - make a variant of this function that takes existing transpose, should caller have it at his disposal (low priority)

	size_t n_dest_col_num = m_n_col_num; //B.m_n_col_num; // t_odo - these and some others do not follow the naming conventions; take care of that
	size_t n_dest_row_num = m_n_col_num; //A.m_n_row_num; // A = B^T
	// these are the dimensions of the new matrix

	{
		r_dest.Clear();
		r_row_list_dest.resize(r_col_list_B.size()); // copy row layout from this
		r_col_list_dest.resize(r_col_list_B.size()); // column layout is the same, except for the blocks
		//std::for_each(r_row_list_dest.begin(), r_row_list_dest.end(), CColumnCumsumToRowCumsum(r_col_list_B.begin()));
		std::transform(r_col_list_B.begin(), r_col_list_B.end(),
			r_row_list_dest.begin(), t_ColumnCumsumToRowCumsum);
		//std::for_each(r_col_list_dest.begin(), r_col_list_dest.end(), CColumnCumsumCopy(r_col_list_B.begin())); // copy column layout but not the blocks
		std::transform(r_col_list_B.begin(), r_col_list_B.end(),
			r_col_list_dest.begin(), t_ColumnCumsumCopy);
		r_dest.m_n_col_num = n_dest_col_num;
		r_dest.m_n_row_num = n_dest_row_num;
		// create layout for the destination matrix (linear time)

		//r_dest.CheckIntegrity();
		// makes sure the dest layout is ok
	}
	// create dest structure

	// no need to reindex / merge, all the indices match and exist

#ifndef __UBER_BLOCK_MATRIX_HYBRID_AT_A
#ifndef __UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
	{
		// this version have straight output to dest, block lookup is done in log(N) time

		size_t n_column_id_B = 0;
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(),
		   p_col_B_end_it = r_col_list_B.end(); p_col_B_it != p_col_B_end_it;
		   ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			for(_TyBlockConstIter p_block_B_it =
			   r_t_column_B.block_list.begin(), p_block_B_end_it = r_t_column_B.block_list.end();
			   p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
				const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
				size_t n_row_id_B = r_t_block_B.first; // get row of the block
				// for each block in the current column in B

				size_t n_bmB_cols;
				_TyMatrixXdRef blockB(r_t_block_B.second, r_row_list_B[n_row_id_B].n_height,
					n_bmB_cols = r_t_column_B.n_width);
				// create map to block B data

				size_t n_column_id_A = n_row_id_B;
				_ASSERTE(n_column_id_A < col_list_A.size());
				const TColumn &r_t_column_A = col_list_A[n_column_id_A];
				_ASSERTE(r_t_column_A.n_width == r_row_list_B[n_row_id_B].n_height);
				// lookup which column in A corresponds with current block row in B

				for(_TyBlockConstIter p_block_A_it =
				   r_t_column_A.block_list.begin(), p_block_A_end_it =
				   r_t_column_A.block_list.end(); p_block_A_it != p_block_A_end_it;
				   ++ p_block_A_it) {
					const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
					size_t n_row_id_A = r_t_block_A.first; // get row of the block
					// for each block in the current column in A

					// multiplication of blockA * blockB yields block at (n_row_id_A, n_column_id_B)

					if(n_row_id_A > n_column_id_B)
						break; // these only get higher, might as well break instead of continue
					// this block would've ended up at lower diagonal, don't want that

					size_t n_bmA_rows; // @todo - is this calculated correctly? might fail with rectangular blocks
					_TyMatrixXdRef blockA_not_transposed(r_t_block_A.second,
						r_t_column_A.n_width, n_bmA_rows = r_row_list_A[n_row_id_A].n_width);
					// create map to block A data

					_ASSERTE(n_row_id_A < r_row_list_dest.size());
					_ASSERTE(n_column_id_B < r_col_list_dest.size());
					//size_t n_prod_rows = r_row_list_dest[n_row_id_A].n_height;
					//size_t n_prod_cols = r_col_list_dest[n_column_id_B].n_width; // unused
					_ASSERTE(r_row_list_dest[n_row_id_A].n_height == blockA_not_transposed.cols());
					_ASSERTE(r_col_list_dest[n_column_id_B].n_width == blockB.cols());
					_ASSERTE(blockA_not_transposed.rows() == blockB.rows()); // make sure the blocks are multiplicable (not able to verify by just merging the layout, but should be since we're multiplying A^T * A)
					// basic checks about matrix dimensions

					bool b_uninitialized;
					double *p_new_block_data = r_dest.p_AllocBlockData(n_row_id_A,
						n_column_id_B, n_bmA_rows, n_bmB_cols, b_uninitialized); // note this clears block data to 0 if not there yet // t_odo - create version with boolean flag saying whether it cleared or not
					// get storage (log(N) time :(, no good way arround it)

					_TyMatrixXdRef block_dest(p_new_block_data, n_bmA_rows, n_bmB_cols);
					if(b_uninitialized)
						block_dest.noalias() = blockA_not_transposed.transpose() * blockB; // initialize a new block
					else
						block_dest.noalias() += blockA_not_transposed.transpose() * blockB; // add to the existing block (run the dot sum)
					// perform the dense multiplication using reference matrices and eigen
				}
			}
		}
		// performs sparse matrix multiplication (linear time in number of blocks * log time to lookup and insert the blocks)

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok
	}
	// performs the multiplication, produces only upper / diagonal portion of the matrix
#else // !__UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
	{
		// this version have output to transpose matrix, then transpose that to dest, block lookup is not needed

#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
		std::vector<size_t> cols_load_list(r_col_list_dest.size(), 0);
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
		std::vector<std::vector<TColumn::TBlockEntry> > transpose_cols_list(r_row_list_dest.size());
		// list for storing transpose columns (note that the sizes and cumsums are not initialized and are invalid)

		size_t n_column_id_B = 0;
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(),
		   p_col_B_end_it = r_col_list_B.end(); p_col_B_it != p_col_B_end_it;
		   ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			for(_TyBlockConstIter p_block_B_it =
			   r_t_column_B.block_list.begin(), p_block_B_end_it = r_t_column_B.block_list.end();
			   p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
				const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
				size_t n_row_id_B = r_t_block_B.first; // get row of the block
				// for each block in the current column in B

				size_t n_bmB_cols;
				_TyMatrixXdRef blockB(r_t_block_B.second,
					r_row_list_B[n_row_id_B].n_height, n_bmB_cols = r_t_column_B.n_width);
				// create map to block B data

				size_t n_column_id_A = n_row_id_B;
				_ASSERTE(n_column_id_A < col_list_A.size());
				const TColumn &r_t_column_A = col_list_A[n_column_id_A];
				_ASSERTE(r_t_column_A.n_width == r_row_list_B[n_row_id_B].n_height);
				// lookup which column in A corresponds with current block row in B

				for(_TyBlockConstIter p_block_A_it =
				   r_t_column_A.block_list.begin(), p_block_A_end_it =
				   r_t_column_A.block_list.end(); p_block_A_it != p_block_A_end_it;
				   ++ p_block_A_it) {
					const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
					size_t n_row_id_A = r_t_block_A.first; // get row of the block
					// for each block in the current column in A

					// multiplication of blockA * blockB yields block at (n_row_id_A, n_column_id_B)

					if(n_row_id_A > n_column_id_B)
						break; // these only get higher, might as well break instead of continue
					// this block would've ended up at lower diagonal, don't want that

					size_t n_bmA_rows; // @todo - is this calculated correctly? might fail with rectangular blocks
					_TyMatrixXdRef blockA_not_transposed(r_t_block_A.second,
						r_t_column_A.n_width, n_bmA_rows = r_row_list_A[n_row_id_A].n_width);
					// create map to block A data

					_ASSERTE(n_row_id_A < r_row_list_dest.size());
					_ASSERTE(n_column_id_B < r_col_list_dest.size());
					//size_t n_prod_rows = r_row_list_dest[n_row_id_A].n_height;
					//size_t n_prod_cols = r_col_list_dest[n_column_id_B].n_width; // unused
					_ASSERTE(r_row_list_dest[n_row_id_A].n_height == blockA_not_transposed.cols());
					_ASSERTE(r_col_list_dest[n_column_id_B].n_width == blockB.cols());
					_ASSERTE(blockA_not_transposed.rows() == blockB.rows()); // make sure the blocks are multiplicable (not able to verify by just merging the layout, but should be since we're multiplying A^T * A)
					// basic checks about matrix dimensions

					_ASSERTE(n_row_id_A < transpose_cols_list.size());
					std::vector<TColumn::TBlockEntry> &r_transpose_column =
						transpose_cols_list[n_row_id_A];
					if(r_transpose_column.empty() ||
					   r_transpose_column.back().first < n_column_id_B) {
						double *p_new_block_data =
							r_dest.p_Get_DenseStorage(n_bmA_rows * n_bmB_cols);
						// get storage

						_TyMatrixXdRef block_dest(p_new_block_data, n_bmA_rows, n_bmB_cols);
						block_dest.noalias() = blockA_not_transposed.transpose() * blockB;
						// initialize a new block

						r_transpose_column.push_back(TColumn::TBlockEntry(n_column_id_B,
							p_new_block_data));
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
						block_dest.noalias() += blockA_not_transposed.transpose() * blockB;
						// add to the existing block (run the dot sum)
					}
					// perform the dense multiplication using reference matrices and eigen
				}
			}
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
				_ASSERTE(r_col_list_dest[r_block.first].block_list.capacity() >
					r_col_list_dest[r_block.first].block_list.size()); // since it was preallocated
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
				r_col_list_dest[r_block.first].block_list.push_back(
					TColumn::TBlockEntry(i, r_block.second));
			}
		}
		// performs the final transpose (linear time in number of blocks)

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok
	}
#endif // !__UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
#else // !__UBER_BLOCK_MATRIX_HYBRID_AT_A
	{
		for(_TyColumnConstIter p_col_A_it =
		   col_list_A.begin(), p_col_A_end_it = col_list_A.end();
		   p_col_A_it != p_col_A_end_it; ++ p_col_A_it) {
			const TColumn &r_t_column_A = *p_col_A_it;
			const size_t n_bm_cols = r_t_column_A.n_width;
			// for each column in A

			for(_TyBlockConstIter p_block_A_it =
			   r_t_column_A.block_list.begin(), p_block_A_end_it = r_t_column_A.block_list.end();
			   p_block_A_it != p_block_A_end_it;) { // increments at the beginning of the inner loop
				const TColumn::TBlockEntry &r_t_block_A0 = *p_block_A_it;
				size_t n_row_id_A0 = r_t_block_A0.first; // get row of the block
				// for each block in the current column in A

/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
				_ASSERTE(n_bm_cols == 3 && r_row_list_A[n_row_id_A0].n_width == 3);
				const size_t n_bmA0_rows = 3; // this is of transposed block
				_TyMatrix3dRef blockA0_not_transposed(r_t_block_A0.second);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
				size_t n_bmA0_rows; // this is of transposed block
				_TyMatrixXdRef blockA0_not_transposed(r_t_block_A0.second,
					n_bm_cols, n_bmA0_rows = r_row_list_A[n_row_id_A0].n_width); // t_odo - this is mixing rows and cols, transpose and straight, definitely buggy // nope, not buggy at all
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
				// create map to block A0 data
				// this block is oriented as in original (B) matrix, the data is not transposed either

				for(_TyBlockConstIter p_block_A1_it =
				   ++ p_block_A_it; p_block_A1_it != p_block_A_end_it; ++ p_block_A1_it) {
					const TColumn::TBlockEntry &r_t_block_A1 = *p_block_A1_it;
					size_t n_row_id_A1 = r_t_block_A1.first; // get row of the block
					// for each next block below the current block in the current column in A

					// multiplication of blockA0 * blockA1 yields block at (n_row_id_A0, n_row_id_A1),
					// guaranteed to be above diagonal

/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
					_ASSERTE(n_bm_cols == 3 && r_row_list_A[n_row_id_A1].n_width == 3);
					const size_t n_bmA1_cols = 3; // this is of straight block
					_TyMatrix3dRef blockA1(r_t_block_A1.second);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
					size_t n_bmA1_cols; // this is of straight block
					_TyMatrixXdRef blockA1(r_t_block_A1.second,
						n_bm_cols, n_bmA1_cols = r_row_list_A[n_row_id_A1].n_width); // t_odo - this is mixing rows and cols, transpose and straight, definitely buggy // no, not buggy at all
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
					// create map to block A1 data
					// this block is oriented as in original (B) matrix, the data is not transposed either

					_ASSERTE(n_row_id_A0 < r_row_list_dest.size());
					_ASSERTE(n_row_id_A1 < r_col_list_dest.size());
					//size_t n_prod_rows = r_row_list_dest[n_row_id_A0].n_height;
					//size_t n_prod_cols = r_col_list_dest[n_row_id_A1].n_width; // unused
					_ASSERTE(r_row_list_dest[n_row_id_A0].n_height ==
						blockA0_not_transposed.cols());
					_ASSERTE(r_col_list_dest[n_row_id_A1].n_width == blockA1.cols());
					_ASSERTE(blockA0_not_transposed.rows() == blockA1.rows()); // make sure the blocks are multiplicable (not able to verify by just merging the layout, but should be since we're multiplying A^T * A)
					// basic checks about matrix dimensions

					/*bool b_uninitialized;
					double *p_new_block_data = r_dest.p_AllocBlockData(n_row_id_A0,
						n_row_id_A1, n_bmA0_rows, n_bmA1_cols, b_uninitialized); // note this clears block data to 0 if not there yet // t_odo - create version with boolean flag saying whether it cleared or not
					// get storage (log(N) time :(, no good way arround it)

//#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
//					_ASSERTE(n_bmA0_rows == 3 && n_bmA1_cols == 3);
//					_TyMatrix3dRef block_dest(p_new_block_data);
//#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
					_TyMatrixXdRef block_dest(p_new_block_data, n_bmA0_rows, n_bmA1_cols);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
					if(b_uninitialized)
						block_dest.noalias() = blockA0_not_transposed.transpose() * blockA1; // initialize a new block
					else
						block_dest.noalias() += blockA0_not_transposed.transpose() * blockA1; // add to the existing block (run the dot sum)
					// perform the dense multiplication using reference matrices and eigen*/

					{
						TColumn &r_t_dest_col = r_col_list_dest[n_row_id_A1];
						_TyBlockIter p_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
							std::lower_bound(r_t_dest_col.block_list.begin(),
							r_t_dest_col.block_list.end(), n_row_id_A0, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
							std::lower_bound(r_t_dest_col.block_list.begin(),
							r_t_dest_col.block_list.end(), n_row_id_A0, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						// find where to put the block in the column

						if((p_block_it == r_t_dest_col.block_list.end()) ||
						   (*p_block_it).first != n_row_id_A0) {
							// a new block

							double *p_value = r_dest.p_Get_DenseStorage(n_bmA0_rows * n_bmA1_cols);
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
							_ASSERTE(n_bmA0_rows == 3 && n_bmA1_cols == 3);
							_TyMatrix3dRef block_dest(p_value);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
							_TyMatrixXdRef block_dest(p_value, n_bmA0_rows, n_bmA1_cols);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS

#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
							block_dest = blockA0_not_transposed.transpose().lazyProduct(blockA1);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
							block_dest.noalias() = blockA0_not_transposed.transpose() * blockA1;
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
							// create a new block, initialize with values

							r_t_dest_col.block_list.insert(p_block_it,
								TColumn::TBlockEntry(n_row_id_A0, p_value));
							// add it to the list
						} else {
							double *p_value_dest = (*p_block_it).second;
							// get existing block data

/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
							_ASSERTE(n_bmA0_rows == 3 && n_bmA1_cols == 3);
							_TyMatrix3dRef block_dest(p_value_dest);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
							_TyMatrixXdRef block_dest(p_value_dest, n_bmA0_rows, n_bmA1_cols);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS

#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
							block_dest += blockA0_not_transposed.transpose().lazyProduct(blockA1);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
							block_dest.noalias() += blockA0_not_transposed.transpose() * blockA1;
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
							// add values to an existing block
						}
					}
					// column is directly indexable, just need to sort blocks inside it
				}
			}
		}
		// perform the multiplication of the blocks above diagonal, need to use log(n) lookup

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok

		size_t n_column_id_B = 0;
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(),
		   p_col_B_end_it = r_col_list_B.end(); p_col_B_it != p_col_B_end_it;
		   ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			if(r_t_column_B.block_list.empty())
				continue;
			// only process non-empty columns

			const size_t n_bm_cols = r_t_column_B.n_width; // B^T * B has shape of B.cols by B.cols
			double *p_block_data = r_dest.p_Get_DenseStorage(n_bm_cols * n_bm_cols);
			memset(p_block_data, 0, n_bm_cols * n_bm_cols * sizeof(double)); // or handle the first block differently
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
			_ASSERTE(n_bm_cols == 3 && n_bm_cols == 3);
			_TyMatrix3dRef block_dest(p_block_data);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
			_TyMatrixXdRef block_dest(p_block_data, n_bm_cols, n_bm_cols);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
			// allocate new block

			for(_TyBlockConstIter p_block_B_it =
			   r_t_column_B.block_list.begin(), p_block_B_end_it = r_t_column_B.block_list.end();
			   p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
				const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
				size_t n_row_id_B = r_t_block_B.first; // get row of the block
				// for each block in the current column in B

/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
				_ASSERTE(r_row_list_B[n_row_id_B].n_height == 3 && n_bm_cols == 3);
				_TyMatrix3dRef blockB(r_t_block_B.second);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
				_TyMatrixXdRef blockB(r_t_block_B.second,
					r_row_list_B[n_row_id_B].n_height, n_bm_cols);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
				// create map to block B data

#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
				block_dest += blockB.transpose().lazyProduct(blockB);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
				block_dest.noalias() += blockB.transpose() * blockB;
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
				// multiply block by its transpose (hope eigen can really optimize this)
			}
			// produce the column dot product, each block with itself

			_ASSERTE(r_col_list_dest[n_column_id_B].block_list.empty() ||
				r_col_list_dest[n_column_id_B].block_list.back().first < n_column_id_B);
			// makes sure that inserting block at diagonal doesn't violate block ordering

			r_col_list_dest[n_column_id_B].block_list.push_back(
				TColumn::TBlockEntry(n_column_id_B, p_block_data));
			// add the block at the end of the list
		}
		// generate the blocks at diagonal (linear lookup, fast multiplication)

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok
	}
#endif // !__UBER_BLOCK_MATRIX_HYBRID_AT_A

	if(!b_upper_diagonal_only) {
		for(size_t i = 0, n = r_col_list_dest.size(); i < n; ++ i) {
			const TColumn &r_t_col = r_col_list_dest[i];
			const size_t n_col = i;
			for(_TyBlockConstIter p_block_it =
			   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
			   p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				const size_t n_row = r_t_block.first;
				_ASSERTE(n_row <= n_col); // make sure the block is above or at the diagonal
				if(n_row < n_col) {
					// it is guaranteed (by the condition above) that this will not alter block i or any of the
					// blocks to be processed next

					_ASSERTE(r_col_list_dest[n_row].block_list.empty() ||
						r_col_list_dest[n_row].block_list.back().first < n_col);
					// make sure that this won't violate block ordering (since i = n_col is monotonically increasing)

					const size_t n_block_width = r_row_list_dest[n_row].n_height,
						n_block_height = r_t_col.n_width; // these are dims of the dest block
					double *p_block_data =
						r_dest.p_Get_DenseStorage(n_block_width * n_block_height);
					r_col_list_dest[n_row].block_list.push_back(
						TColumn::TBlockEntry(n_col, p_block_data));
					// this block is above the diagonal, at (n_row, n_col), create transpose block at (n_col, n_row)

					const double *p_src_data = r_t_block.second;
					/*for(size_t x = 0; x < n_block_width; ++ x)
						for(size_t y = 0; y < n_block_height; ++ y)
							p_block_data[y + n_block_height * x] = p_src_data[x + n_block_width * y];*/ // fixme - does the src / dest ordering matter? is it interchangable? // t_odo - test transposing a matrix // this seems to be correct
					_TyMatrixXdRef src((double*)p_src_data, n_block_width, n_block_height),
						dest(p_block_data, n_block_height, n_block_width);
					dest = src.transpose();
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

	//r_dest.CheckIntegrity();
	// makes sure the dest matrix is ok
}

bool CUberBlockMatrix::MultiplyToWith(CUberBlockMatrix &r_dest,
	const CUberBlockMatrix &r_other) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);
	r_other.CheckIntegrity(true);

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

	{
		r_dest.Clear();
		r_row_list_dest = r_row_list_A; // copy row layout from this
		r_col_list_dest.resize(r_col_list_B.size());
		//std::for_each(r_col_list_dest.begin(), r_col_list_dest.end(), CColumnCumsumCopy(r_col_list_B.begin())); // copy column layout but not the blocks
		std::transform(r_col_list_B.begin(), r_col_list_B.end(),
			r_col_list_dest.begin(), t_ColumnCumsumCopy);
		r_dest.m_n_col_num = n_dest_col_num;
		r_dest.m_n_row_num = n_dest_row_num;
		// create layout for the destination matrix (linear time)

		//r_dest.CheckIntegrity();
		// makes sure the dest layout is ok
	}

	std::vector<size_t> reindex_rows_B_to_cols_A;
	{
		{
			std::vector<size_t> reindex_cols_A;
			size_t n_common_size = n_MergeLayout(r_col_list_A,
				r_row_list_B, reindex_cols_A, reindex_rows_B_to_cols_A);
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
	}
	// merge the common dimension layout (linear time)
	// -1 means the row did not map/exist in B, -2 meants it did not map/exist in A

	{
#ifndef __UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
		// this version have straight output to dest, block lookup is done in log(N) time

		size_t n_column_id_B = 0;
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(),
		   p_col_B_end_it = r_col_list_B.end(); p_col_B_it != p_col_B_end_it;
		   ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			for(_TyBlockConstIter
			   p_block_B_it = r_t_column_B.block_list.begin(),
			   p_block_B_end_it = r_t_column_B.block_list.end();
			   p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
				const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
				size_t n_row_id_B = r_t_block_B.first; // get row of the block
				// for each block in the current column in B

				size_t n_bmB_cols;
				//Eigen::R_eferencingMatrixXd blockB(r_row_list_B[n_row_id_B].n_height,
				//	n_bmB_cols = r_t_column_B.n_width, r_t_block_B.second);
				_TyMatrixXdRef blockB(r_t_block_B.second,
					r_row_list_B[n_row_id_B].n_height, n_bmB_cols = r_t_column_B.n_width);
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

				for(_TyBlockConstIter p_block_A_it =
				   r_t_column_A.block_list.begin(), p_block_A_end_it =
				   r_t_column_A.block_list.end(); p_block_A_it != p_block_A_end_it;
				   ++ p_block_A_it) {
					const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
					size_t n_row_id_A = r_t_block_A.first; // get row of the block
					// for each block in the current column in A

					size_t n_bmA_rows;
					//Eigen::R_eferencingMatrixXd blockA(n_bmA_rows = r_row_list_A[n_row_id_A].n_height,
					//	r_t_column_A.n_width, r_t_block_A.second);
					_TyMatrixXdRef blockA(r_t_block_A.second, n_bmA_rows =
						r_row_list_A[n_row_id_A].n_height, r_t_column_A.n_width);
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

					{
						bool b_uninitialized;
						double *p_new_block_data = r_dest.p_AllocBlockData(n_row_id_A,
							n_column_id_B, n_bmA_rows, n_bmB_cols, b_uninitialized); // note this clears block data to 0 if not there yet // t_odo - create version with boolean flag saying whether it cleared or not
						// get storage (log(N) time :(, no good way arround it)

						//Eigen::R_eferencingMatrixXd block_dest(n_bmA_rows, n_bmB_cols, p_new_block_data);
						_TyMatrixXdRef block_dest(p_new_block_data,
							n_bmA_rows, n_bmB_cols);
						if(b_uninitialized)
							block_dest.noalias() = blockA * blockB; // initialize a new block
						else
							block_dest.noalias() += blockA * blockB; // add to the existing block (run the dot sum)
						// perform the dense multiplication using reference matrices and eigen
					}
				}
			}
		}
		// performs sparse matrix multiplication (linear time in number of blocks * log time to lookup and insert the blocks)

		// note that with proper loop reorganization, it should be possible to create
		// a linear time version (the resulting blocks would come out in sorted order) // todo - investigate the possibility
		// evidently you have to traverse column in B and a row in A at the same time. traversing row in A is log(N) again,
		// but if A is transposed, then it's linear. the multiplication then reduces to merging two sorted lists, right?

		// or transpose the result, that will allow for direct modification of this loop and will be very simple
#else // !__UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
		// this version have output to transpose matrix, then transpose that to dest, block lookup is not needed

#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
		std::vector<size_t> cols_load_list(r_col_list_dest.size(), 0);
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
		std::vector<std::vector<TColumn::TBlockEntry> >
			transpose_cols_list(r_row_list_dest.size());
		// list for storing transpose columns (note that the sizes and cumsums are not initialized and are invalid)

		size_t n_column_id_B = 0;
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(),
		   p_col_B_end_it = r_col_list_B.end(); p_col_B_it != p_col_B_end_it;
		   ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			for(_TyBlockConstIter p_block_B_it =
			   r_t_column_B.block_list.begin(), p_block_B_end_it = r_t_column_B.block_list.end();
			   p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
				const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
				size_t n_row_id_B = r_t_block_B.first; // get row of the block
				// for each block in the current column in B

				size_t n_bmB_cols;
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
				_ASSERTE(r_row_list_B[n_row_id_B].n_height == 3);
				_ASSERTE(r_t_column_B.n_width == 3);
				n_bmB_cols = 3;
				_TyMatrix3dRef blockB(r_t_block_B.second);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
				_TyMatrixXdRef blockB(r_t_block_B.second, r_row_list_B[n_row_id_B].n_height,
					n_bmB_cols = r_t_column_B.n_width);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
				// create map to block B data

				size_t n_column_id_A = reindex_rows_B_to_cols_A[n_row_id_B];
				_ASSERTE(size_t(-1) > size_t(-2)); // just to make sure the next line is correct
				if(n_column_id_A >= size_t(-2)) {
					if(n_column_id_A == size_t(-1))
						return false; // didn't map from B to common and we know it was not empty (we have a block here)
					continue; // do not fail, it might also mean that there are no blocks in that column in A and hence the result of multiplication is zero
				}
				_ASSERTE(n_column_id_A < r_col_list_A.size());
				const TColumn &r_t_column_A = r_col_list_A[n_column_id_A];
				_ASSERTE(r_t_column_A.n_width == r_row_list_B[n_row_id_B].n_height);
				// lookup which column in A corresponds with current block row in B

				for(_TyBlockConstIter p_block_A_it =
				   r_t_column_A.block_list.begin(), p_block_A_end_it =
				   r_t_column_A.block_list.end(); p_block_A_it != p_block_A_end_it;
				   ++ p_block_A_it) {
					const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
					size_t n_row_id_A = r_t_block_A.first; // get row of the block
					// for each block in the current column in A

					size_t n_bmA_rows;
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
					n_bmA_rows = 3;
					_TyMatrix3dRef blockA(r_t_block_A.second);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
					_TyMatrixXdRef blockA(r_t_block_A.second,
						n_bmA_rows = r_row_list_A[n_row_id_A].n_height, r_t_column_A.n_width);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
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
					std::vector<TColumn::TBlockEntry> &r_transpose_column =
						transpose_cols_list[n_row_id_A];
					if(r_transpose_column.empty() ||
					   r_transpose_column.back().first < n_column_id_B) {
						double *p_new_block_data =
							r_dest.p_Get_DenseStorage(n_bmA_rows * n_bmB_cols);
						// get storage

						_TyMatrixXdRef block_dest(p_new_block_data, n_bmA_rows, n_bmB_cols);
						block_dest = blockA * blockB;
						// initialize a new block

						r_transpose_column.push_back(
							TColumn::TBlockEntry(n_column_id_B, p_new_block_data));
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
#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						block_dest += blockA.lazyProduct(blockB);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						block_dest.noalias() += blockA * blockB;
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						// add to the existing block (run the dot sum)
					}
					// perform the dense multiplication using reference matrices and eigen
				}
			}
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
				_ASSERTE(r_col_list_dest[r_block.first].block_list.capacity() >
					r_col_list_dest[r_block.first].block_list.size()); // since it was preallocated
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
				r_col_list_dest[r_block.first].block_list.push_back(
					TColumn::TBlockEntry(i, r_block.second));
			}
		}
		// performs the final transpose (linear time in number of blocks)
#endif // !__UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok
	}

	// note it might be beneficial (and possibly quite easy) to implement Strassen's algorithm here.
	// would be a nice work for a paper, too.

	// no need to implement more cunning algorithms, though (Wiki):
	// The Coppersmith�Winograd algorithm is frequently used as a building block in other algorithms to
	// prove theoretical time bounds. However, unlike the Strassen algorithm, it is not used in practice
	// because it only provides an advantage for matrices so large that they cannot be processed by modern
	// hardware.

	return true;
}

bool CUberBlockMatrix::MultiplyToWith(CUberBlockMatrix &r_dest,
	const CUberBlockMatrix &r_other, bool b_upper_diag_only) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);
	r_other.CheckIntegrity(true);

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

	{
		r_dest.Clear();
		r_row_list_dest = r_row_list_A; // copy row layout from this
		r_col_list_dest.resize(r_col_list_B.size());
		//std::for_each(r_col_list_dest.begin(), r_col_list_dest.end(), CColumnCumsumCopy(r_col_list_B.begin())); // copy column layout but not the blocks
		std::transform(r_col_list_B.begin(), r_col_list_B.end(),
			r_col_list_dest.begin(), t_ColumnCumsumCopy);
		r_dest.m_n_col_num = n_dest_col_num;
		r_dest.m_n_row_num = n_dest_row_num;
		// create layout for the destination matrix (linear time)

		//r_dest.CheckIntegrity();
		// makes sure the dest layout is ok
	}

	std::vector<size_t> reindex_rows_B_to_cols_A;
	{
		{
			std::vector<size_t> reindex_cols_A;
			size_t n_common_size = n_MergeLayout(r_col_list_A,
				r_row_list_B, reindex_cols_A, reindex_rows_B_to_cols_A);
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
	}
	// merge the common dimension layout (linear time)
	// -1 means the row did not map/exist in B, -2 meants it did not map/exist in A

	{
#ifndef __UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
		// this version have straight output to dest, block lookup is done in log(N) time

		size_t n_column_id_B = 0;
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(),
		   p_col_B_end_it = r_col_list_B.end(); p_col_B_it != p_col_B_end_it;
		   ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			for(_TyBlockConstIter
			   p_block_B_it = r_t_column_B.block_list.begin(),
			   p_block_B_end_it = r_t_column_B.block_list.end();
			   p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
				const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
				size_t n_row_id_B = r_t_block_B.first; // get row of the block
				// for each block in the current column in B

				size_t n_bmB_cols;
				//Eigen::R_eferencingMatrixXd blockB(r_row_list_B[n_row_id_B].n_height,
				//	n_bmB_cols = r_t_column_B.n_width, r_t_block_B.second);
				_TyMatrixXdRef blockB(r_t_block_B.second,
					r_row_list_B[n_row_id_B].n_height, n_bmB_cols = r_t_column_B.n_width);
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

				for(_TyBlockConstIter p_block_A_it =
				   r_t_column_A.block_list.begin(), p_block_A_end_it =
				   r_t_column_A.block_list.end(); p_block_A_it != p_block_A_end_it;
				   ++ p_block_A_it) {
					const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
					size_t n_row_id_A = r_t_block_A.first; // get row of the block
					// for each block in the current column in A

					if(b_upper_diag_only && n_row_id_A > n_column_id_B)
						break; // these only get higher, might as well break instead of continue
					// this block would've ended up at lower diagonal, don't want that

					size_t n_bmA_rows;
					//Eigen::R_eferencingMatrixXd blockA(n_bmA_rows = r_row_list_A[n_row_id_A].n_height,
					//	r_t_column_A.n_width, r_t_block_A.second);
					_TyMatrixXdRef blockA(r_t_block_A.second, n_bmA_rows =
						r_row_list_A[n_row_id_A].n_height, r_t_column_A.n_width);
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

					{
						bool b_uninitialized;
						double *p_new_block_data = r_dest.p_AllocBlockData(n_row_id_A,
							n_column_id_B, n_bmA_rows, n_bmB_cols, b_uninitialized); // note this clears block data to 0 if not there yet // t_odo - create version with boolean flag saying whether it cleared or not
						// get storage (log(N) time :(, no good way arround it)

						//Eigen::R_eferencingMatrixXd block_dest(n_bmA_rows, n_bmB_cols, p_new_block_data);
						_TyMatrixXdRef block_dest(p_new_block_data,
							n_bmA_rows, n_bmB_cols);
						if(b_uninitialized)
							block_dest.noalias() = blockA * blockB; // initialize a new block
						else
							block_dest.noalias() += blockA * blockB; // add to the existing block (run the dot sum)
						// perform the dense multiplication using reference matrices and eigen
					}
				}
			}
		}
		// performs sparse matrix multiplication (linear time in number of blocks * log time to lookup and insert the blocks)

		// note that with proper loop reorganization, it should be possible to create
		// a linear time version (the resulting blocks would come out in sorted order) // todo - investigate the possibility
		// evidently you have to traverse column in B and a row in A at the same time. traversing row in A is log(N) again,
		// but if A is transposed, then it's linear. the multiplication then reduces to merging two sorted lists, right?

		// or transpose the result, that will allow for direct modification of this loop and will be very simple
#else // !__UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
		// this version have output to transpose matrix, then transpose that to dest, block lookup is not needed

#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
		std::vector<size_t> cols_load_list(r_col_list_dest.size(), 0);
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
		std::vector<std::vector<TColumn::TBlockEntry> >
			transpose_cols_list(r_row_list_dest.size());
		// list for storing transpose columns (note that the sizes and cumsums are not initialized and are invalid)

		size_t n_column_id_B = 0;
		for(_TyColumnConstIter p_col_B_it = r_col_list_B.begin(),
		   p_col_B_end_it = r_col_list_B.end(); p_col_B_it != p_col_B_end_it;
		   ++ p_col_B_it, ++ n_column_id_B) {
			const TColumn &r_t_column_B = *p_col_B_it;
			// for each column in B (index is n_column_id_B)

			for(_TyBlockConstIter p_block_B_it =
			   r_t_column_B.block_list.begin(), p_block_B_end_it = r_t_column_B.block_list.end();
			   p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
				const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
				size_t n_row_id_B = r_t_block_B.first; // get row of the block
				// for each block in the current column in B

				size_t n_bmB_cols;
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
				_ASSERTE(r_row_list_B[n_row_id_B].n_height == 3);
				_ASSERTE(r_t_column_B.n_width == 3);
				n_bmB_cols = 3;
				_TyMatrix3dRef blockB(r_t_block_B.second);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
				_TyMatrixXdRef blockB(r_t_block_B.second, r_row_list_B[n_row_id_B].n_height,
					n_bmB_cols = r_t_column_B.n_width);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
				// create map to block B data

				size_t n_column_id_A = reindex_rows_B_to_cols_A[n_row_id_B];
				_ASSERTE(size_t(-1) > size_t(-2)); // just to make sure the next line is correct
				if(n_column_id_A >= size_t(-2)) {
					if(n_column_id_A == size_t(-1))
						return false; // didn't map from B to common and we know it was not empty (we have a block here)
					continue; // do not fail, it might also mean that there are no blocks in that column in A and hence the result of multiplication is zero
				}
				_ASSERTE(n_column_id_A < r_col_list_A.size());
				const TColumn &r_t_column_A = r_col_list_A[n_column_id_A];
				_ASSERTE(r_t_column_A.n_width == r_row_list_B[n_row_id_B].n_height);
				// lookup which column in A corresponds with current block row in B

				for(_TyBlockConstIter p_block_A_it =
				   r_t_column_A.block_list.begin(), p_block_A_end_it =
				   r_t_column_A.block_list.end(); p_block_A_it != p_block_A_end_it;
				   ++ p_block_A_it) {
					const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
					size_t n_row_id_A = r_t_block_A.first; // get row of the block
					// for each block in the current column in A

					if(b_upper_diag_only && n_row_id_A > n_column_id_B)
						break; // these only get higher, might as well break instead of continue
					// this block would've ended up at lower diagonal, don't want that

					size_t n_bmA_rows;
/*#ifdef __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
					n_bmA_rows = 3;
					_TyMatrix3dRef blockA(r_t_block_A.second);
#else // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS*/
					_TyMatrixXdRef blockA(r_t_block_A.second,
						n_bmA_rows = r_row_list_A[n_row_id_A].n_height, r_t_column_A.n_width);
//#endif // __UBER_BLOCK_MATRIX_FORCE_3x3_BLOCKS
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
					std::vector<TColumn::TBlockEntry> &r_transpose_column =
						transpose_cols_list[n_row_id_A];
					if(r_transpose_column.empty() ||
					   r_transpose_column.back().first < n_column_id_B) {
						double *p_new_block_data =
							r_dest.p_Get_DenseStorage(n_bmA_rows * n_bmB_cols);
						// get storage

						_TyMatrixXdRef block_dest(p_new_block_data, n_bmA_rows, n_bmB_cols);
						block_dest = blockA * blockB;
						// initialize a new block

						r_transpose_column.push_back(
							TColumn::TBlockEntry(n_column_id_B, p_new_block_data));
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
#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						block_dest += blockA.lazyProduct(blockB);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						block_dest.noalias() += blockA * blockB;
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						// add to the existing block (run the dot sum)
					}
					// perform the dense multiplication using reference matrices and eigen
				}
			}
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
				_ASSERTE(r_col_list_dest[r_block.first].block_list.capacity() >
					r_col_list_dest[r_block.first].block_list.size()); // since it was preallocated
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
				r_col_list_dest[r_block.first].block_list.push_back(
					TColumn::TBlockEntry(i, r_block.second));
			}
		}
		// performs the final transpose (linear time in number of blocks)
#endif // !__UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR

		//r_dest.CheckIntegrity();
		// makes sure the dest matrix is ok
	}

	// note it might be beneficial (and possibly quite easy) to implement Strassen's algorithm here.
	// would be a nice work for a paper, too.

	// no need to implement more cunning algorithms, though (Wiki):
	// The Coppersmith�Winograd algorithm is frequently used as a building block in other algorithms to
	// prove theoretical time bounds. However, unlike the Strassen algorithm, it is not used in practice
	// because it only provides an advantage for matrices so large that they cannot be processed by modern
	// hardware.

	return true;
}

TBmp *CUberBlockMatrix::p_Rasterize(TBmp *p_storage /*= 0*/, int n_scalar_size /*= 5*/) const
{
	CheckIntegrity(true);

	const uint32_t n_background = 0xffffffffU,
		n_row_col_border = 0xff808080U,
		n_nonzero = 0xff8080ffU,
		//n_nonzero_new = 0xffff0000U,
		n_zero_new = 0xff00ff00U;
	// colors of the individual areas

	const uint32_t n_nonzero_border = 0x80000000U | ((n_nonzero >> 1) & 0x7f7f7f7fU);
	//const uint32_t n_nonzero_new_border = 0x80000000U | ((n_nonzero_new >> 1) & 0x7f7f7f7fU),
		//n_zero_new_border = 0x80000000U | ((n_zero_new >> 1) & 0x7f7f7f7fU);
	// colors of borders (just darker)

	size_t m = m_n_row_num;
	size_t n = m_n_col_num; // in fact, it's the other way arround, but lets not confuse things

	if(m == SIZE_MAX || n == SIZE_MAX || m > INT_MAX || n > INT_MAX ||
	   m * (n_scalar_size - 1) + 1 > INT_MAX || n * (n_scalar_size - 1) + 1 > INT_MAX ||
	   uint64_t(m * (n_scalar_size - 1) + 1) * (m * (n_scalar_size - 1) + 1) > INT_MAX)
		return 0;
	/*if(m == SIZE_MAX)
		throw std::runtime_error("m == SIZE_MAX, would lead to infinite loop");
	if(n == SIZE_MAX)
		throw std::runtime_error("n == SIZE_MAX, would lead to infinite loop");*/

	TBmp *p_bitmap;
	if(p_storage && p_storage->n_width == int(n * (n_scalar_size - 1) + 1) &&
	   p_storage->n_height == int(m * (n_scalar_size - 1) + 1))
		p_bitmap = p_storage; // use the existing storage
	else {
		if(!(p_bitmap = TBmp::p_Alloc(int(n * (n_scalar_size - 1) + 1),
		   int(m * (n_scalar_size - 1) + 1))))
			return 0;
		// need to alloc a new bitmap

		if(p_storage) {
			p_storage->Free();
			*p_storage = *p_bitmap;
			delete p_bitmap;
			p_bitmap = p_storage;
		}
		// reuse storage for bitmap structure (but use buffer of the new bitmap)
	}
	// alloc storage

	p_bitmap->Clear(n_background);

#if 1
	for(size_t j = 0, n_col_num = m_block_cols_list.size(); j < n_col_num; ++ j) {
		const TColumn &t_col = m_block_cols_list[j];
		for(size_t i = 0; i <= m; ++ i) {
			p_bitmap->PutPixel(int(t_col.n_cumulative_width_sum * (n_scalar_size - 1)),
				int(i * (n_scalar_size - 1)), n_row_col_border);
		}
	}
	for(size_t i = 0; i <= m; ++ i) {
		p_bitmap->PutPixel(int(n * (n_scalar_size - 1)),
			int(i * (n_scalar_size - 1)), n_row_col_border);
	}
	for(size_t j = 0, n_row_num = m_block_rows_list.size(); j < n_row_num; ++ j) {
		const TRow &t_row = m_block_rows_list[j];
		for(size_t i = 0; i <= n; ++ i) {
			p_bitmap->PutPixel(int(i * (n_scalar_size - 1)), int(t_row.n_cumulative_height_sum *
				(n_scalar_size - 1)), n_row_col_border);
		}
	}
	for(size_t i = 0; i <= n; ++ i) {
		p_bitmap->PutPixel(int(i * (n_scalar_size - 1)), int(m *
			(n_scalar_size - 1)), n_row_col_border);
	}
	// rasterize divisions between the rows / columns
#endif // 0

	for(size_t j = 0, n_col_num = m_block_cols_list.size(); j < n_col_num; ++ j) {
		const TColumn &t_col = m_block_cols_list[j];

		size_t n_x = t_col.n_cumulative_width_sum * (n_scalar_size - 1);
		size_t n_width = t_col.n_width;
		// dimensions of the block

		for(size_t i = 0, n_block_num = t_col.block_list.size(); i < n_block_num; ++ i) {
			size_t n_row = t_col.block_list[i].first;
			const double *p_data = t_col.block_list[i].second;

			size_t n_y = m_block_rows_list[n_row].n_cumulative_height_sum * (n_scalar_size - 1);
			size_t n_height = m_block_rows_list[n_row].n_height;
			// dimensions of the block

			{
				for(size_t dy = 1; dy < n_height * (n_scalar_size - 1); ++ dy) {
					for(size_t dx = 1; dx < n_width * (n_scalar_size - 1); ++ dx) {
						bool b_is_border = false;
						if(n_scalar_size > 1) {
							b_is_border = ((dx - 1) % size_t(n_scalar_size - 1) ==
								(size_t(n_scalar_size - 1) - 1)) ||
								((dy - 1) % size_t(n_scalar_size - 1) ==
								(size_t(n_scalar_size - 1) - 1));
							// makes the matrices look more symmetric
						}
						size_t n_col = (dx - 1) / size_t(n_scalar_size - 1);
						size_t n_row = (dy - 1) / size_t(n_scalar_size - 1);
						_ASSERTE(n_col < n_width);
						_ASSERTE(n_row < n_height);
						p_bitmap->PutPixel(int(n_x + dx), int(n_y + dy), ((b_is_border ||
							fabs(p_data[n_col * n_height + n_row]) > 1e-10)?
							n_nonzero : n_zero_new));
					}
				}
				// fills the area (colors for zero/nonzero entries)

				for(size_t dy = 1; dy < n_height; ++ dy) {
					for(size_t dx = 1; dx < n_width; ++ dx) {
						p_bitmap->PutPixel(int(n_x + dx * (n_scalar_size - 1)),
							int(n_y + dy * (n_scalar_size - 1)), n_nonzero_border);
					}
				}
				// draws the separation between the block elements (as dots in the corners)

				p_bitmap->DrawLine_SP(float(n_x), float(n_y),
					float(n_x + n_width * (n_scalar_size - 1)), float(n_y), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_x), float(n_y + n_height * (n_scalar_size - 1)),
					float(n_x + n_width * (n_scalar_size - 1)),
					float(n_y + n_height * (n_scalar_size - 1)), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_x), float(n_y),
					float(n_x), float(n_y + n_height * (n_scalar_size - 1)), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_x + n_width * (n_scalar_size - 1)), float(n_y),
					float(n_x + n_width * (n_scalar_size - 1)),
					float(n_y + n_height * (n_scalar_size - 1)), n_nonzero_border);
				// draws the borders
			}
			// draw a square into the bitmap (containing the whole block) // t_odo - separate scalars somehow
		}
	}

	return p_bitmap;
}

TBmp *CUberBlockMatrix::p_Rasterize_Symmetric(TBmp *p_storage /*= 0*/, int n_scalar_size /*= 5*/) const
{
	CheckIntegrity(true);

	const uint32_t n_background = 0xffffffffU,
		n_row_col_border = 0xff808080U,
		n_nonzero = 0xff8080ffU,
		//n_nonzero_new = 0xffff0000U,
		n_zero_new = 0xff00ff00U;
	// colors of the individual areas

	const uint32_t n_nonzero_border = 0x80000000U | ((n_nonzero >> 1) & 0x7f7f7f7fU);
	//const uint32_t n_nonzero_new_border = 0x80000000U | ((n_nonzero_new >> 1) & 0x7f7f7f7fU),
		//n_zero_new_border = 0x80000000U | ((n_zero_new >> 1) & 0x7f7f7f7fU);
	// colors of borders (just darker)

	size_t m = m_n_row_num;
	size_t n = m_n_col_num; // in fact, it's the other way arround, but lets not confuse things

	if(m != n || !b_SymmetricLayout())
		return 0;
	// must be square and symmetric

	if(m == SIZE_MAX || n == SIZE_MAX || m > INT_MAX || n > INT_MAX ||
	   m * (n_scalar_size - 1) + 1 > INT_MAX || n * (n_scalar_size - 1) + 1 > INT_MAX ||
	   uint64_t(m * (n_scalar_size - 1) + 1) * (m * (n_scalar_size - 1) + 1) > INT_MAX)
		return 0;
	/*if(m == SIZE_MAX)
		throw std::runtime_error("m == SIZE_MAX, would lead to infinite loop");
	if(n == SIZE_MAX)
		throw std::runtime_error("n == SIZE_MAX, would lead to infinite loop");*/

	TBmp *p_bitmap;
	if(p_storage && p_storage->n_width == int(n * (n_scalar_size - 1) + 1) &&
	   p_storage->n_height == int(m * (n_scalar_size - 1) + 1))
		p_bitmap = p_storage; // use the existing storage
	else {
		if(!(p_bitmap = TBmp::p_Alloc(int(n * (n_scalar_size - 1) + 1),
		   int(m * (n_scalar_size - 1) + 1))))
			return 0;
		// need to alloc a new bitmap

		if(p_storage) {
			p_storage->Free();
			*p_storage = *p_bitmap;
			delete p_bitmap;
			p_bitmap = p_storage;
		}
		// reuse storage for bitmap structure (but use buffer of the new bitmap)
	}
	// alloc storage

	p_bitmap->Clear(n_background);

#if 1
	for(size_t j = 0, n_col_num = m_block_cols_list.size(); j < n_col_num; ++ j) {
		const TColumn &t_col = m_block_cols_list[j];
		for(size_t i = 0; i <= m; ++ i) {
			p_bitmap->PutPixel(int(t_col.n_cumulative_width_sum * (n_scalar_size - 1)),
				int(i * (n_scalar_size - 1)), n_row_col_border);
		}
	}
	for(size_t i = 0; i <= m; ++ i) {
		p_bitmap->PutPixel(int(n * (n_scalar_size - 1)),
			int(i * (n_scalar_size - 1)), n_row_col_border);
	}
	for(size_t j = 0, n_row_num = m_block_rows_list.size(); j < n_row_num; ++ j) {
		const TRow &t_row = m_block_rows_list[j];
		for(size_t i = 0; i <= n; ++ i) {
			p_bitmap->PutPixel(int(i * (n_scalar_size - 1)), int(t_row.n_cumulative_height_sum *
				(n_scalar_size - 1)), n_row_col_border);
		}
	}
	for(size_t i = 0; i <= n; ++ i) {
		p_bitmap->PutPixel(int(i * (n_scalar_size - 1)), int(m *
			(n_scalar_size - 1)), n_row_col_border);
	}
	// rasterize divisions between the rows / columns
#endif // 0

	for(size_t j = 0, n_col_num = m_block_cols_list.size(); j < n_col_num; ++ j) {
		const TColumn &t_col = m_block_cols_list[j];

		size_t n_x = t_col.n_cumulative_width_sum * (n_scalar_size - 1);
		size_t n_width = t_col.n_width;
		// dimensions of the block

		for(size_t i = 0, n_block_num = t_col.block_list.size(); i < n_block_num; ++ i) {
			size_t n_row = t_col.block_list[i].first;
			const double *p_data = t_col.block_list[i].second;

			size_t n_y = m_block_rows_list[n_row].n_cumulative_height_sum * (n_scalar_size - 1);
			size_t n_height = m_block_rows_list[n_row].n_height;
			// dimensions of the block

			{
				for(size_t dy = 1; dy < n_height * (n_scalar_size - 1); ++ dy) {
					for(size_t dx = 1; dx < n_width * (n_scalar_size - 1); ++ dx) {
						bool b_is_border = false;
						if(n_scalar_size > 1) {
							b_is_border = ((dx - 1) % size_t(n_scalar_size - 1) ==
								(size_t(n_scalar_size - 1) - 1)) ||
								((dy - 1) % size_t(n_scalar_size - 1) ==
								(size_t(n_scalar_size - 1) - 1));
							// makes the matrices look more symmetric
						}
						size_t n_col = (dx - 1) / size_t(n_scalar_size - 1);
						size_t n_row = (dy - 1) / size_t(n_scalar_size - 1);
						_ASSERTE(n_col < n_width);
						_ASSERTE(n_row < n_height);
						uint32_t n_color = ((b_is_border || fabs(p_data[n_col * n_height +
							n_row]) > 1e-10)? n_nonzero : n_zero_new);
						p_bitmap->PutPixel(int(n_x + dx), int(n_y + dy), n_color);
						p_bitmap->PutPixel(int(n_y + dy), int(n_x + dx), n_color); // transpose
					}
				}
				// fills the area (colors for zero/nonzero entries)

				for(size_t dy = 1; dy < n_height; ++ dy) {
					for(size_t dx = 1; dx < n_width; ++ dx) {
						p_bitmap->PutPixel(int(n_x + dx * (n_scalar_size - 1)),
							int(n_y + dy * (n_scalar_size - 1)), n_nonzero_border);
						p_bitmap->PutPixel(int(n_y + dy * (n_scalar_size - 1)),
							int(n_x + dx * (n_scalar_size - 1)), n_nonzero_border); // transpose
					}
				}
				// draws the separation between the block elements (as dots in the corners)

				p_bitmap->DrawLine_SP(float(n_x), float(n_y),
					float(n_x + n_width * (n_scalar_size - 1)), float(n_y), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_x), float(n_y + n_height * (n_scalar_size - 1)),
					float(n_x + n_width * (n_scalar_size - 1)),
					float(n_y + n_height * (n_scalar_size - 1)), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_x), float(n_y),
					float(n_x), float(n_y + n_height * (n_scalar_size - 1)), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_x + n_width * (n_scalar_size - 1)), float(n_y),
					float(n_x + n_width * (n_scalar_size - 1)),
					float(n_y + n_height * (n_scalar_size - 1)), n_nonzero_border);
				// draws the borders

				p_bitmap->DrawLine_SP(float(n_y), float(n_x),
					float(n_y), float(n_x + n_width * (n_scalar_size - 1)), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_y + n_height * (n_scalar_size - 1)), float(n_x),
					float(n_y + n_height * (n_scalar_size - 1)),
					float(n_x + n_width * (n_scalar_size - 1)), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_y), float(n_x),
					float(n_y + n_height * (n_scalar_size - 1)), float(n_x), n_nonzero_border);
				p_bitmap->DrawLine_SP(float(n_y), float(n_x + n_width * (n_scalar_size - 1)),
					float(n_y + n_height * (n_scalar_size - 1)),
					float(n_x + n_width * (n_scalar_size - 1)), n_nonzero_border);
				// and the transpose borders
			}
			// draw a square into the bitmap (containing the whole block) // t_odo - separate scalars somehow
		}
	}

	return p_bitmap;
}

TBmp *CUberBlockMatrix::p_Rasterize(const CUberBlockMatrix &r_prev_state,
	bool b_handle_changed_as_new, TBmp *p_storage /*= 0*/, int n_scalar_size /*= 5*/) const
{
	CheckIntegrity(true);

	_ASSERTE(r_prev_state.m_n_row_num <= m_n_row_num && r_prev_state.m_n_col_num <= m_n_col_num);
	_ASSERTE(r_prev_state.m_block_rows_list.size() <= m_block_rows_list.size() &&
		r_prev_state.m_block_cols_list.size() <= m_block_cols_list.size());
#ifdef _DEBUG
	for(size_t i = 0, n = r_prev_state.m_block_rows_list.size(); i < n; ++ i) {
		_ASSERTE(r_prev_state.m_block_rows_list[i].n_height ==
			m_block_rows_list[i].n_height);
		// cumulative heights will be the same if the heights are the same
	}
	for(size_t i = 0, n = r_prev_state.m_block_cols_list.size(); i < n; ++ i) {
		_ASSERTE(r_prev_state.m_block_cols_list[i].n_width ==
			m_block_cols_list[i].n_width);
		// cumulative widths will be the same if the widths are the same
	}
#endif // _DEBUG
	// make sure that r_prev_state is a prefix of this (at least the block layout)

	const uint32_t n_background = 0xffffffffU,
		n_row_col_border = 0xff808080U,
		n_nonzero = 0xff8080ffU,
		n_nonzero_new = 0xffff0000U,
		n_zero_new = 0xff00ff00U;
	// colors of the individual areas

	const uint32_t n_nonzero_border = 0x80000000U | ((n_nonzero >> 1) & 0x7f7f7f7fU);
	const uint32_t n_nonzero_new_border = 0x80000000U | ((n_nonzero_new >> 1) & 0x7f7f7f7fU);
		//n_zero_new_border = 0x80000000U | ((n_zero_new >> 1) & 0x7f7f7f7fU);
	// colors of borders (just darker)

	size_t m = m_n_row_num;
	size_t n = m_n_col_num; // in fact, it's the other way arround, but lets not confuse things

	if(m == SIZE_MAX || n == SIZE_MAX || m > INT_MAX || n > INT_MAX ||
	   m * (n_scalar_size - 1) + 1 > INT_MAX || n * (n_scalar_size - 1) + 1 > INT_MAX ||
	   uint64_t(m * (n_scalar_size - 1) + 1) * (m * (n_scalar_size - 1) + 1) > INT_MAX)
		return 0;
	/*if(m == SIZE_MAX)
		throw std::runtime_error("m == SIZE_MAX, would lead to infinite loop");
	if(n == SIZE_MAX)
		throw std::runtime_error("n == SIZE_MAX, would lead to infinite loop");*/

	TBmp *p_bitmap;
	if(p_storage && p_storage->n_width == int(n * (n_scalar_size - 1) + 1) &&
	   p_storage->n_height == int(m * (n_scalar_size - 1) + 1))
		p_bitmap = p_storage; // use the existing storage
	else {
		if(!(p_bitmap = TBmp::p_Alloc(int(n * (n_scalar_size - 1) + 1),
		   int(m * (n_scalar_size - 1) + 1))))
			return 0;
		// need to alloc a new bitmap

		if(p_storage) {
			p_storage->Free();
			*p_storage = *p_bitmap;
			delete p_bitmap;
			p_bitmap = p_storage;
		}
		// reuse storage for bitmap structure (but use buffer of the new bitmap)
	}
	// alloc storage

	p_bitmap->Clear(n_background);

#if 1
	for(size_t j = 0, n_col_num = m_block_cols_list.size(); j < n_col_num; ++ j) {
		const TColumn &t_col = m_block_cols_list[j];
		for(size_t i = 0; i <= m; ++ i) {
			p_bitmap->PutPixel(int(t_col.n_cumulative_width_sum * (n_scalar_size - 1)),
				int(i * (n_scalar_size - 1)), n_row_col_border);
		}
	}
	for(size_t i = 0; i <= m; ++ i) {
		p_bitmap->PutPixel(int(n * (n_scalar_size - 1)),
			int(i * (n_scalar_size - 1)), n_row_col_border);
	}
	for(size_t j = 0, n_row_num = m_block_rows_list.size(); j < n_row_num; ++ j) {
		const TRow &t_row = m_block_rows_list[j];
		for(size_t i = 0; i <= n; ++ i) {
			p_bitmap->PutPixel(int(i * (n_scalar_size - 1)), int(t_row.n_cumulative_height_sum *
				(n_scalar_size - 1)), n_row_col_border);
		}
	}
	for(size_t i = 0; i <= n; ++ i) {
		p_bitmap->PutPixel(int(i * (n_scalar_size - 1)), int(m *
			(n_scalar_size - 1)), n_row_col_border);
	}
	// rasterize divisions between the rows / columns
#endif // 0

	size_t n_prev_col_num = r_prev_state.m_block_cols_list.size();
	for(size_t j = 0, n_col_num = m_block_cols_list.size(); j < n_col_num; ++ j) {
		const TColumn &t_col = m_block_cols_list[j];

		size_t n_x = t_col.n_cumulative_width_sum * (n_scalar_size - 1);
		size_t n_width = t_col.n_width;
		// dimensions of the block

		size_t n_prev_block_lbound = 0;
		size_t UNUSED(n_prev_block_num) = (j < n_prev_col_num)?
			r_prev_state.m_block_cols_list[j].block_list.size() : 0;
		// prepare a lower bound and number of blocks in the same column in the "previous" matrix

		for(size_t i = 0, n_block_num = t_col.block_list.size(); i < n_block_num; ++ i) {
			size_t n_row = t_col.block_list[i].first;
			const double *p_data = t_col.block_list[i].second;

			const double *p_prev_data = 0;
			bool b_new_block = true;
			if(j < n_prev_col_num) {
				const TColumn &t_prev_col = r_prev_state.m_block_cols_list[j];
				_ASSERTE(n_prev_block_num == t_prev_col.block_list.size()); // make sure it is up to date
				_ASSERTE((n_prev_block_lbound < n_prev_block_num && // n_prev_block_lbound is a valid index
					(!n_prev_block_lbound || (n_prev_block_lbound && // and either it is 0 and we know nothing
					t_prev_col.block_list[n_prev_block_lbound - 1].first < n_row))) || // or it is greater than 0 and it must be a lower bound
					(n_prev_block_lbound == n_prev_block_num && // or we are past the previous matrix
					(t_prev_col.block_list.empty() || // and it is either empty ... or the last block is indeed above the current one
					t_prev_col.block_list.back().first <= n_row))); // make sure that the lower bound is correct
				// get a column in the prev matrix

				_TyBlockConstIter p_prev_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
					std::lower_bound(t_prev_col.block_list.begin(),
					t_prev_col.block_list.end(), n_row, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(t_prev_col.block_list.begin(),
					t_prev_col.block_list.end(), n_row, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				// look for a block with the same row as the current block

				n_prev_block_lbound = p_prev_block_it - t_prev_col.block_list.begin();
				// remember the new lower bound

				if(p_prev_block_it != t_prev_col.block_list.end() &&
				   (*p_prev_block_it).first == n_row) {
					b_new_block = false;
					p_prev_data = (*p_prev_block_it).second;

					++ n_prev_block_lbound;
					// this one was already used
				}
				// see if we found
			}
			// find out if there is the same block in the prev matrix

			size_t n_y = m_block_rows_list[n_row].n_cumulative_height_sum * (n_scalar_size - 1);
			size_t n_height = m_block_rows_list[n_row].n_height;
			// dimensions of the block

			{
				if(b_handle_changed_as_new) {
					bool b_changed = (p_prev_data)? memcmp(p_prev_data,
						p_data, n_width * n_height * sizeof(double)) != 0 : false;
					if(b_changed)
						b_new_block = true;
				}
				// draw changed blocks as new ones

				uint32_t n_nnz_color = (b_new_block)? n_nonzero_new : n_nonzero;

				for(size_t dy = 1; dy < n_height * (n_scalar_size - 1); ++ dy) {
					for(size_t dx = 1; dx < n_width * (n_scalar_size - 1); ++ dx) {
						bool b_is_border = false;
						if(n_scalar_size > 1) {
							b_is_border = ((dx - 1) % size_t(n_scalar_size - 1) ==
								(size_t(n_scalar_size - 1) - 1)) ||
								((dy - 1) % size_t(n_scalar_size - 1) ==
								(size_t(n_scalar_size - 1) - 1));
							// makes the matrices look more symmetric
						}
						size_t n_col = (dx - 1) / size_t(n_scalar_size - 1);
						size_t n_row = (dy - 1) / size_t(n_scalar_size - 1);
						_ASSERTE(n_col < n_width);
						_ASSERTE(n_row < n_height);
						double f_data = fabs(p_data[n_col * n_height + n_row]);
						double f_prev_data = (p_prev_data)?
							fabs(p_prev_data[n_col * n_height + n_row]) : 0;
						uint32_t n_color = n_nnz_color;
						if(!b_new_block && !b_is_border && f_prev_data <= 1e-10 && f_data > 1e-10)
							n_color = n_nonzero_new; // highlight new nonzeros inside preexisting blocks
						p_bitmap->PutPixel(int(n_x + dx), int(n_y + dy),
							(b_is_border || f_data > 1e-10)? n_color : n_zero_new);
					}
				}
				// fills the area (colors for zero/nonzero entries)

				uint32_t n_border_color = (b_new_block)? n_nonzero_new_border : n_nonzero_border;

				for(size_t dy = 1; dy < n_height; ++ dy) {
					for(size_t dx = 1; dx < n_width; ++ dx) {
						p_bitmap->PutPixel(int(n_x + dx * (n_scalar_size - 1)),
							int(n_y + dy * (n_scalar_size - 1)), n_border_color);
					}
				}
				// draws the separation between the block elements (as dots in the corners)

				p_bitmap->DrawLine_SP(float(n_x), float(n_y),
					float(n_x + n_width * (n_scalar_size - 1)), float(n_y), n_border_color);
				p_bitmap->DrawLine_SP(float(n_x), float(n_y + n_height * (n_scalar_size - 1)),
					float(n_x + n_width * (n_scalar_size - 1)),
					float(n_y + n_height * (n_scalar_size - 1)), n_border_color);
				p_bitmap->DrawLine_SP(float(n_x), float(n_y),
					float(n_x), float(n_y + n_height * (n_scalar_size - 1)), n_border_color);
				p_bitmap->DrawLine_SP(float(n_x + n_width * (n_scalar_size - 1)), float(n_y),
					float(n_x + n_width * (n_scalar_size - 1)),
					float(n_y + n_height * (n_scalar_size - 1)), n_border_color);
				// draws the borders
			}
			// draw a square into the bitmap (containing the whole block)
		}
	}

	return p_bitmap;
}

bool CUberBlockMatrix::Rasterize(const char *p_s_filename, int n_scalar_size /*= 5*/) const
{
	TBmp *p_bitmap;
	if(!(p_bitmap = p_Rasterize(0, n_scalar_size)))
		return false;

	bool b_result = CTgaCodec::Save_TGA(p_s_filename, *p_bitmap, false, true);

	p_bitmap->Delete();

	return b_result;
}

bool CUberBlockMatrix::Rasterize_Symmetric(const char *p_s_filename, int n_scalar_size /*= 5*/) const
{
	TBmp *p_bitmap;
	if(!(p_bitmap = p_Rasterize_Symmetric(0, n_scalar_size)))
		return false;

	bool b_result = CTgaCodec::Save_TGA(p_s_filename, *p_bitmap, false, true);

	p_bitmap->Delete();

	return b_result;
}

bool CUberBlockMatrix::Rasterize_Symmetric(const CUberBlockMatrix &r_prev_state,
	bool b_handle_changed_as_new, const char *p_s_filename, int n_scalar_size /*= 5*/) const
{
	TBmp *p_bitmap;
	if(!(p_bitmap = p_Rasterize(r_prev_state, b_handle_changed_as_new, 0, n_scalar_size)))
		return false;

	if(p_bitmap->n_width != p_bitmap->n_height) { // !!
		p_bitmap->Delete();
		return false;
	}

	for(int y = 0, w = p_bitmap->n_width, h = p_bitmap->n_height; y < h; ++ y)
		for(int x = 0; x < y; ++ x)
			p_bitmap->p_buffer[x + y * w] = p_bitmap->p_buffer[y + x * w];
	// transpose the upper half to the lower half (in the image only)

	bool b_result = CTgaCodec::Save_TGA(p_s_filename, *p_bitmap, false, true);

	p_bitmap->Delete();

	return b_result;
}

bool CUberBlockMatrix::Rasterize(const CUberBlockMatrix &r_prev_state,
	bool b_handle_changed_as_new, const char *p_s_filename, int n_scalar_size /*= 5*/) const
{
	TBmp *p_bitmap;
	if(!(p_bitmap = p_Rasterize(r_prev_state, b_handle_changed_as_new, 0, n_scalar_size)))
		return false;

	bool b_result = CTgaCodec::Save_TGA(p_s_filename, *p_bitmap, false, true);

	p_bitmap->Delete();

	return b_result;
}

bool CUberBlockMatrix::b_RowReferenced(size_t n_row_index) const
{
	_ASSERTE(n_row_index < m_block_rows_list.size()); // make sure it's a valid index

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		const std::vector<TColumn::TBlockEntry> &r_block_list =
			m_block_cols_list[i].block_list;
		if(std::find_if(r_block_list.begin(), r_block_list.end(),
		   CReferencesRow(n_row_index)) != r_block_list.end())
			return true; // it's referenced
	}
	// go through all columns
	// note this could be faster if the column with the same index as the
	// row was checked first (elements tend to lie on the diagonal)

	return false; // not referenced
}

void CUberBlockMatrix::Shift_RowReferences(size_t n_row_index)
{
	_ASSERTE(!b_RowReferenced(n_row_index));
	/*_ASSERTE((n_row_index && !b_RowReferenced(n_row_index - 1)) ||
		(n_row_index != m_block_rows_list.size() && !b_RowReferenced(n_row_index + 1)));*/ // doesn't work like that, the row indices are not shifted with expanding row vector, can't check for this.
	// it's supposed to be index of the unused row that was just fragmented

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) { // todo - for each
		std::vector<TColumn::TBlockEntry> &r_block_list =
			m_block_cols_list[i].block_list;
		std::for_each(r_block_list.begin(), r_block_list.end(),
			CIncrementRowIndex<1>(n_row_index));
	}
	// go through all columns
}

void CUberBlockMatrix::Shift_RowReferences2(size_t n_row_index0, size_t UNUSED(n_row_index1))
{
	_ASSERTE(n_row_index1 == n_row_index0 + 2);
	_ASSERTE(!b_RowReferenced(n_row_index0));
	/*_ASSERTE(!b_RowReferenced(n_row_index0 + 1));
	_ASSERTE(!b_RowReferenced(n_row_index1));*/ // doesn't work like that, the row indices are not shifted with expanding row vector, can't check for this.
	// there is only one call to this function, when a long and unused row
	// needs to be fragmented. make sure that is the case.

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) { // todo - for each
		std::vector<TColumn::TBlockEntry> &r_block_list =
			m_block_cols_list[i].block_list;
		std::for_each(r_block_list.begin(), r_block_list.end(),
			CIncrementRowIndex<2>(n_row_index0));
	}
	// go through all columns
}

double *CUberBlockMatrix::p_AllocBlockData(size_t n_row_index, TColumn &r_t_col,
	size_t n_block_row_num, size_t n_block_column_num, bool &r_b_was_a_new_block) // throw(std::bad_alloc)
{
	_ASSERTE(n_row_index < m_block_rows_list.size());
	_ASSERTE(m_block_rows_list[n_row_index].n_height == n_block_row_num);
	_ASSERTE(r_t_col.n_width == n_block_column_num);
	double *p_data;
	if(r_t_col.block_list.empty() || r_t_col.block_list.back().first < n_row_index) { // no blocks / last block is before (should be true for 'nice' matrix filling methods - where blocks are added to columns in already sorted fashion)
		p_data = p_Get_DenseStorage(n_block_row_num * n_block_column_num);
		r_b_was_a_new_block = true;
		r_t_col.block_list.push_back(TColumn::TBlockEntry(n_row_index, p_data)); // just add to the end
	} else {
		_TyBlockIter p_block_dest_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::lower_bound(r_t_col.block_list.begin(), r_t_col.block_list.end(),
			n_row_index, CCompareBlockRow()); // t_odo - use lower bound instead of upper. this will never find blocks it's looking for // t_odo - make sure this works
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::lower_bound(r_t_col.block_list.begin(), r_t_col.block_list.end(),
			n_row_index, CompareBlockRow); // t_odo - use lower bound instead of upper. this will never find blocks it's looking for // t_odo - make sure this works
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		if(p_block_dest_it == r_t_col.block_list.end() ||
		   (*p_block_dest_it).first != n_row_index) {
			p_data = p_Get_DenseStorage(n_block_row_num * n_block_column_num);
			r_b_was_a_new_block = true;
			_ASSERTE(p_block_dest_it == r_t_col.block_list.end() ||
				(*p_block_dest_it).first > n_row_index); // make sure it will not unorder
			r_t_col.block_list.insert(p_block_dest_it, TColumn::TBlockEntry(n_row_index, p_data)); // a new block
		} else {
			p_data = (*p_block_dest_it).second; // CCompareBlockRow is ok, probably points to the right block // t_odo - fixme
			// the block at the given position already exists, no need to add another one

			r_b_was_a_new_block = false; // !! don't forget to set it

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_col_block_reref_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		}
		// t_odo - is there any preference between upper and lower bound? // yes, if we're looking for an element that possibly exists, always use lower_bound (guaranteed to return it, if it exists)

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
		++ m_n_col_block_search_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
	}
	// allocate a new block

	//CheckIntegrity(); // leads to exponential cost when inserting N elements in a matrix
	// make sure that allocating the block didn't damage matrix integrity

	return p_data;
}

const double *CUberBlockMatrix::p_GetBlockData(size_t n_row_index, size_t n_column_index) const
{
	const TColumn &r_t_col = m_block_cols_list[n_column_index];
	const double *p_data;
	if(r_t_col.block_list.empty() || r_t_col.block_list.back().first < n_row_index) // no blocks / last block is before (should be true for 'nice' matrix filling methods - where blocks are added to columns in already sorted fashion)
		return 0;
	else {
		_TyBlockConstIter p_block_dest_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::lower_bound(r_t_col.block_list.begin(), r_t_col.block_list.end(),
			n_row_index, CCompareBlockRow()); // t_odo - use lower bound instead of upper. this will never find blocks it's looking for // t_odo - make sure this works
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::lower_bound(r_t_col.block_list.begin(), r_t_col.block_list.end(),
			n_row_index, CompareBlockRow); // t_odo - use lower bound instead of upper. this will never find blocks it's looking for // t_odo - make sure this works
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		if(p_block_dest_it == r_t_col.block_list.end() || (*p_block_dest_it).first != n_row_index)
			return 0;
		else {
			p_data = (*p_block_dest_it).second; // CCompareBlockRow is ok, probably points to the right block // t_odo - fixme
			// the block at the given position already exists, no need to add another one

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_col_block_reref_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		}
		// t_odo - is there any preference between upper and lower bound? // yes, if we're looking for an element that possibly exists, always use lower_bound (guaranteed to return it, if it exists)

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
		++ m_n_col_block_search_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
	}
	// allocate a new block

	//CheckIntegrity(); // leads to exponential cost when inserting N elements in a matrix
	// make sure that allocating the block didn't damage matrix integrity

	return p_data;
}

double *CUberBlockMatrix::p_GetBlockData(size_t n_row_index, size_t n_column_index)
{
	const TColumn &r_t_col = m_block_cols_list[n_column_index];
	double *p_data;
	if(r_t_col.block_list.empty() || r_t_col.block_list.back().first < n_row_index) // no blocks / last block is before (should be true for 'nice' matrix filling methods - where blocks are added to columns in already sorted fashion)
		return 0;
	else {
		_TyBlockConstIter p_block_dest_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::lower_bound(r_t_col.block_list.begin(), r_t_col.block_list.end(),
			n_row_index, CCompareBlockRow()); // t_odo - use lower bound instead of upper. this will never find blocks it's looking for // t_odo - make sure this works
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::lower_bound(r_t_col.block_list.begin(), r_t_col.block_list.end(),
			n_row_index, CompareBlockRow); // t_odo - use lower bound instead of upper. this will never find blocks it's looking for // t_odo - make sure this works
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		if(p_block_dest_it == r_t_col.block_list.end() || (*p_block_dest_it).first != n_row_index)
			return 0;
		else {
			p_data = (*p_block_dest_it).second; // CCompareBlockRow is ok, probably points to the right block // t_odo - fixme
			// the block at the given position already exists, no need to add another one

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
			++ m_n_col_block_reref_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
		}
		// t_odo - is there any preference between upper and lower bound? // yes, if we're looking for an element that possibly exists, always use lower_bound (guaranteed to return it, if it exists)

#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
		++ m_n_col_block_search_num;
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
	}
	// allocate a new block

	//CheckIntegrity(); // leads to exponential cost when inserting N elements in a matrix
	// make sure that allocating the block didn't damage matrix integrity

	return p_data;
}

void CUberBlockMatrix::PasteTo(CUberBlockMatrix &r_dest,
	size_t n_dest_row_index, size_t n_dest_column_index) // throw(std::bad_alloc)
{
	_ASSERTE(&r_dest != this); // can't append to itself

	CheckIntegrity(true);

	_ASSERTE(r_dest.m_block_rows_list.size() > n_dest_row_index);
	_ASSERTE(r_dest.m_block_cols_list.size() > n_dest_column_index);
	_ASSERTE(r_dest.m_block_rows_list.size() >= n_dest_row_index + m_block_rows_list.size());
	_ASSERTE(r_dest.m_block_cols_list.size() >= n_dest_column_index + m_block_cols_list.size());
	// make sure there is space for the slice

#ifdef _DEBUG
	size_t n_base_dest_row = r_dest.m_block_rows_list[n_dest_row_index].n_cumulative_height_sum;
	size_t n_base_dest_col = r_dest.m_block_cols_list[n_dest_column_index].n_cumulative_width_sum;
#endif // _DEBUG
	// todo - check row layout as well

	/*r_dest.Rasterize("paste0_original.tga");
	Rasterize("paste1_insert.tga");*/
	// debug

	const size_t m = m_block_rows_list.size();
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) { // todo - make iterator loop instead
		TColumn &r_t_dest = r_dest.m_block_cols_list[i + n_dest_column_index];
		const TColumn &r_t_src = m_block_cols_list[i];
		const size_t n_width = r_t_src.n_width;
		// get src / dest columns

		_ASSERTE(r_t_src.block_list.empty() || (r_t_dest.n_width == r_t_src.n_width &&
			r_t_dest.n_cumulative_width_sum == r_t_src.n_cumulative_width_sum + n_base_dest_col));
		// make sure that the columns match

		_TyBlockIter p_dst_range_begin =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::lower_bound(r_t_dest.block_list.begin(), r_t_dest.block_list.end(),
			n_dest_row_index, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::lower_bound(r_t_dest.block_list.begin(), r_t_dest.block_list.end(),
			n_dest_row_index, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		_TyBlockConstIter p_dst_range_end = // note it might be faster to just check for the end in iterating the range, especially if the number of blocks is very low
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::upper_bound(r_t_dest.block_list.begin(), r_t_dest.block_list.end(),
			n_dest_row_index + m, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::upper_bound(r_t_dest.block_list.begin(), r_t_dest.block_list.end(),
			n_dest_row_index + m, _CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		// find where the first block comes to (note upper_bound() used to find the end - one past the last)

		_TyBlockConstIter p_src_block_it = r_t_src.block_list.begin(),
			p_src_block_end_it = r_t_src.block_list.end();
		// get a range of blocks to insert

		if(p_dst_range_begin != p_dst_range_end) {
			for(; p_src_block_it != p_src_block_end_it; ++ p_src_block_it) {
				size_t n_src_row = (*p_src_block_it).first;
				size_t n_dest_row = n_src_row + n_dest_row_index;

				_ASSERTE(p_dst_range_begin < p_dst_range_end);
				while((*p_dst_range_begin).first < n_dest_row) {
					size_t n_row_height = r_dest.m_block_rows_list[
						(*p_dst_range_begin).first].n_height;
					memset((*p_dst_range_begin).second, 0,
						n_width * n_row_height * sizeof(double));
					if(++ p_dst_range_begin == p_dst_range_end)
						break;
				}
				if(p_dst_range_begin == p_dst_range_end)
					break;
				// zero all the blocks prior to the next inserted block // todo - do that? add a flag!

				if((*p_dst_range_begin).first == n_dest_row) {
					size_t n_row_height = m_block_rows_list[n_src_row].n_height;
					_ASSERTE(n_row_height == r_dest.m_block_rows_list[n_dest_row].n_height);
					_ASSERTE(m_block_rows_list[n_src_row].n_cumulative_height_sum + n_base_dest_row ==
						r_dest.m_block_rows_list[n_dest_row].n_cumulative_height_sum);
					memcpy((*p_dst_range_begin).second, (*p_src_block_it).second,
						n_width * n_row_height * sizeof(double));
					// there is a block that is rewritten - just copy data

					if(++ p_dst_range_begin == p_dst_range_end) {
						++ p_src_block_it; // !!
						// tricky - still, worth not moving this condition up as
						// the branch below always knows it is not at the end
						break;
					}
					// can happen
				} else {
					_ASSERTE((*p_dst_range_begin).first > n_dest_row);
					size_t n_row_height = m_block_rows_list[n_src_row].n_height;
					_ASSERTE(n_row_height == r_dest.m_block_rows_list[n_dest_row].n_height);
					_ASSERTE(m_block_rows_list[n_src_row].n_cumulative_height_sum + n_base_dest_row ==
						r_dest.m_block_rows_list[n_dest_row].n_cumulative_height_sum);
					double *p_new_block = r_dest.p_Get_DenseStorage(n_width * n_row_height);
					memcpy(p_new_block, (*p_src_block_it).second,
						n_width * n_row_height * sizeof(double));
					size_t n_remaining = p_dst_range_end - p_dst_range_begin;
					p_dst_range_begin = r_t_dest.block_list.insert(p_dst_range_begin,
						TColumn::TBlockEntry(n_dest_row, p_new_block)) + 1;
					p_dst_range_end = p_dst_range_begin + n_remaining; // mind pointer invalidation
					_ASSERTE(p_dst_range_begin < r_t_dest.block_list.end());
					_ASSERTE(p_dst_range_end <= r_t_dest.block_list.end()); // should still end up with valid iterators
					// a new block; alloc data, copy, insert to the list
				}
			}
		}
		// merge blocks (solves collisions and sorts the blocks)

		if(p_dst_range_begin == p_dst_range_end) {
			if(p_dst_range_begin == r_t_dest.block_list.end()) {
				// no more dest blocks

				r_t_dest.block_list.reserve(r_t_dest.block_list.size() +
					(p_src_block_end_it - p_src_block_it));
				// avoid reallocs (note p_dst_range_begin is invalidated)

				for(; p_src_block_it != p_src_block_end_it; ++ p_src_block_it) {
					size_t n_src_row = (*p_src_block_it).first;
					size_t n_dest_row = n_src_row + n_dest_row_index;

					size_t n_row_height = m_block_rows_list[n_src_row].n_height;
					_ASSERTE(n_row_height == r_dest.m_block_rows_list[n_dest_row].n_height);
					_ASSERTE(m_block_rows_list[n_src_row].n_cumulative_height_sum + n_base_dest_row ==
						r_dest.m_block_rows_list[n_dest_row].n_cumulative_height_sum);
					double *p_new_block = r_dest.p_Get_DenseStorage(n_width * n_row_height);
					memcpy(p_new_block, (*p_src_block_it).second,
						n_width * n_row_height * sizeof(double));

					r_t_dest.block_list.push_back(
						TColumn::TBlockEntry(n_dest_row, p_new_block)); // can use push_back() here
					// a new block; alloc data, copy, insert to the list
				}
			} else {
				// dest blocks below inserted range

				size_t n_begin_off = p_dst_range_begin - r_t_dest.block_list.begin();
				r_t_dest.block_list.reserve(r_t_dest.block_list.size() +
					(p_src_block_end_it - p_src_block_it));
				p_dst_range_begin = r_t_dest.block_list.begin() + n_begin_off; // mind iterator invalidation
				// avoid reallocs

				for(; p_src_block_it != p_src_block_end_it; ++ p_src_block_it) {
					size_t n_src_row = (*p_src_block_it).first;
					size_t n_dest_row = n_src_row + n_dest_row_index;

					size_t n_row_height = m_block_rows_list[n_src_row].n_height;
					_ASSERTE(n_row_height == r_dest.m_block_rows_list[n_dest_row].n_height);
					_ASSERTE(m_block_rows_list[n_src_row].n_cumulative_height_sum + n_base_dest_row ==
						r_dest.m_block_rows_list[n_dest_row].n_cumulative_height_sum);
					double *p_new_block = r_dest.p_Get_DenseStorage(n_width * n_row_height);
					memcpy(p_new_block, (*p_src_block_it).second,
						n_width * n_row_height * sizeof(double));

					p_dst_range_begin = r_t_dest.block_list.insert(p_dst_range_begin,
						TColumn::TBlockEntry(n_dest_row, p_new_block)) + 1;
					// a new block; alloc data, copy, insert to the list
				}
			}
		}
		// all the other dest blocks are below, just insert
	}
	// todo - avoid memset and memcpy and use matrix refs instead

	/*r_dest.Rasterize("paste2_modified.tga");
	Rasterize("paste3_insert-after.tga");*/
	// debug

	//r_dest.CheckIntegrity();
	// make sure the dest matrix is ok
}

void CUberBlockMatrix::CopyLayoutTo(CUberBlockMatrix &r_dest) const
{
	_ASSERTE(&r_dest != this); // can't copy layout to itself

	CheckIntegrity();

#if 0
	SliceTo(r_dest, 0, n_BlockRow_Num(), 0, n_BlockColumn_Num(), true);
	// copy all, reference the data

	for(size_t i = 0, n = n_BlockColumn_Num(); i < n; ++ i)
		r_dest.m_block_cols_list[i].block_list.clear();
	r_dest.m_data_pool.clear();
	r_dest.m_n_ref_elem_num = 0;
	// remove blocks, free data
#else // 0
	r_dest.Clear();
	// erase dest ...

	r_dest.m_n_row_num = m_n_row_num;
	r_dest.m_n_col_num = m_n_col_num;
	// set matrix size

	r_dest.m_block_rows_list = m_block_rows_list;
	// set block rows (just copy)

	r_dest.m_block_cols_list.resize(m_block_cols_list.size());
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) { // todo - make iterator loop instead
		TColumn &r_t_dest = r_dest.m_block_cols_list[i];
		const TColumn &r_t_src = m_block_cols_list[i];
		// get src / dest columns

		r_t_dest.n_width = r_t_src.n_width;
		r_t_dest.n_cumulative_width_sum = r_t_src.n_cumulative_width_sum;
		// copy column dimensions
	}
	// set block columns, copy block data (only blocks in the valid column range)

	//r_dest.CheckIntegrity();
	// make sure the dest matrix is ok
#endif // 0
}

void CUberBlockMatrix::SetZero()
{
	m_data_pool.clear();
	m_n_ref_elem_num = 0;
	// free any allocated blocks, discount referenced data

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i)
		m_block_cols_list[i].block_list.clear();
	// clear all the block lists

	CheckIntegrity();
}

void CUberBlockMatrix::SetIdentity() // throw(std::bad_alloc)
{
	_ASSERTE(b_SymmetricLayout());

	m_data_pool.clear();
	m_n_ref_elem_num = 0;
	// free any allocated blocks, discount referenced data

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		TColumn &r_t_column = m_block_cols_list[i];
		r_t_column.block_list.clear();
		r_t_column.block_list.reserve(1);
		// delete any previous contents, make sure there will be space for one block

		TColumn::TBlockEntry t_block;
		t_block.first = i; // a diagonal block
		t_block.second = p_Get_DenseStorage(r_t_column.n_width * r_t_column.n_width);
		r_t_column.block_list.push_back(t_block);
		// allocate a single diagonal block

		_TyMatrixXdRef block(t_block.second, r_t_column.n_width, r_t_column.n_width);
		block.setIdentity();
		// set identity
	}

	CheckIntegrity(true);
}

void CUberBlockMatrix::SliceTo(CUberBlockMatrix &r_dest, size_t n_block_row_num,
	size_t n_block_column_num, bool b_share_data /*= false*/) // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	_ASSERTE(!n_block_row_num == !n_block_column_num); // make sure they are either both empty or both full
	_ASSERTE(n_block_row_num <= m_block_rows_list.size());
	_ASSERTE(n_block_column_num <= m_block_cols_list.size());

	r_dest.m_n_col_num = (n_block_column_num)? m_block_cols_list[n_block_column_num - 1].n_cumulative_width_sum +
		m_block_cols_list[n_block_column_num - 1].n_width : 0;
	r_dest.m_n_row_num = (n_block_row_num)? m_block_rows_list[n_block_row_num - 1].n_cumulative_height_sum +
		m_block_rows_list[n_block_row_num - 1].n_height : 0;
	// don't forget to set the size

	if(/*&r_dest == this ||*/ !b_share_data) {
		_TyPool data_pool;
		_TyDenseAllocator alloc(data_pool);
		// a temp pool

		if(&r_dest != this) {
			r_dest.m_block_cols_list.clear();
			r_dest.m_block_cols_list.insert(r_dest.m_block_cols_list.begin(),
				m_block_cols_list.begin(), m_block_cols_list.begin() + n_block_column_num);
			r_dest.m_block_rows_list.clear();
			r_dest.m_block_rows_list.insert(r_dest.m_block_rows_list.begin(),
				m_block_rows_list.begin(), m_block_rows_list.begin() + n_block_row_num);
			// clear and copy
		} else {
			m_block_rows_list.resize(n_block_row_num);
			m_block_cols_list.resize(n_block_column_num);
			// just throw awawy leftovers
		}
		// initialize block cols / rows

		for(_TyColumnIter p_col_it = r_dest.m_block_cols_list.begin(),
		   p_end_it = r_dest.m_block_cols_list.begin() + n_block_column_num;
		   p_col_it != p_end_it; ++ p_col_it) {
			TColumn &r_t_col = *p_col_it;
			size_t n_width = r_t_col.n_width;
			for(_TyBlockIter p_block_it = r_t_col.block_list.begin(), p_block_end_it =
			   r_t_col.block_list.end(); p_block_it != p_block_end_it; ++ p_block_it) {
				size_t n_row_id = (*p_block_it).first;
				if(n_row_id < n_block_row_num) {
					size_t n_height = m_block_rows_list[n_row_id].n_height;
					double *p_new = alloc.p_Get_DenseStorage(n_width * n_height);
					memcpy(p_new, (*p_block_it).second, n_width * n_height * sizeof(double));
					(*p_block_it).second = p_new;
					// relocate the block to the new data pool
				} else {
					r_t_col.block_list.erase(p_block_it, p_block_end_it);
					// this is the end of this column; delete any leftover blocks
					break;
				}
			}
		}
		// go through all the columns in the matrix, vacate the ones below thresh, delete the rest

		r_dest.m_data_pool.swap(data_pool);
		// use the new pool, effectively delete (free) the old data

		r_dest.m_n_ref_elem_num = 0;
		// now all the blocks are here
	} else {
		// this is largerly untested // todo test this

		if(&r_dest != this) {
			r_dest.m_block_cols_list.clear();
			r_dest.m_block_cols_list.insert(r_dest.m_block_cols_list.begin(),
				m_block_cols_list.begin(), m_block_cols_list.begin() + n_block_column_num);
			r_dest.m_block_rows_list.clear();
			r_dest.m_block_rows_list.insert(r_dest.m_block_rows_list.begin(),
				m_block_rows_list.begin(), m_block_rows_list.begin() + n_block_row_num);
			// clear and copy
		} else {
			// do this later
		}
		// initialize block cols / rows

		if(&r_dest != this) {
			r_dest.m_n_ref_elem_num = 0;
			for(_TyColumnIter p_col_it = r_dest.m_block_cols_list.begin(),
			   p_end_it = r_dest.m_block_cols_list.begin() + n_block_column_num;
			   p_col_it != p_end_it; ++ p_col_it) {
				TColumn &r_t_col = *p_col_it;
				size_t n_width = r_t_col.n_width, n_height_sum = 0;
				for(_TyBlockIter p_block_it = r_t_col.block_list.begin(), p_block_end_it =
				   r_t_col.block_list.end(); p_block_it != p_block_end_it; ++ p_block_it) {
					size_t n_row_id = (*p_block_it).first;
					if(n_row_id < n_block_row_num) {
						size_t n_height = m_block_rows_list[n_row_id].n_height;
						n_height_sum += n_height;
						// just count the blocks
					} else {
						r_t_col.block_list.erase(p_block_it, p_block_end_it);
						// this is the end of this column; delete any leftover blocks
						break;
					}
				}
				r_dest.m_n_ref_elem_num += n_width * n_height_sum;
			}
			// go through all the columns in the matrix, vacate the ones below thresh, delete the rest
		} else {
			size_t n_min_index = 0;
			for(_TyColumnIter p_col_it = r_dest.m_block_cols_list.begin(),
			   p_end_it = r_dest.m_block_cols_list.begin() + n_block_column_num;
			   p_col_it != p_end_it; ++ p_col_it) {
				TColumn &r_t_col = *p_col_it;
				size_t n_width = r_t_col.n_width;
				for(_TyBlockIter p_block_it = r_t_col.block_list.begin(), p_block_end_it =
				   r_t_col.block_list.end(); p_block_it != p_block_end_it; ++ p_block_it) {
					size_t n_row_id = (*p_block_it).first;
					if(n_row_id < n_block_row_num) {
						size_t n_height = m_block_rows_list[n_row_id].n_height;
						size_t n_index = m_data_pool.index_of((*p_block_it).second);
						if(n_index == size_t(-1)) // !! need to handle those, otherwise we'll get really strange results
							continue; // a block, referenced from a different matrix
						n_index += n_width * n_height;
						if(n_min_index < n_index)
							n_min_index = n_index;
						// relocate the block to the new data pool
					} else {
						r_t_col.block_list.erase(p_block_it, p_block_end_it);
						// this is the end of this column; delete any leftover blocks
						break;
					}
				}
			}
			// go through all the blocks and record indices that are in use

#ifdef _DEBUG
			size_t n_dead_block_num = 0;
			size_t n_dead_blocks_size = 0;
			for(_TyColumnIter p_col_it = r_dest.m_block_cols_list.begin(),
			   p_end_it = r_dest.m_block_cols_list.begin() + n_block_column_num;
			   p_col_it != p_end_it; ++ p_col_it) {
				TColumn &r_t_col = *p_col_it;
				size_t n_width = r_t_col.n_width;
				for(_TyBlockIter p_block_it = r_t_col.block_list.begin(), p_block_end_it =
				   r_t_col.block_list.end(); p_block_it != p_block_end_it; ++ p_block_it) {
					size_t n_row_id = (*p_block_it).first;
					if(n_row_id >= n_block_row_num) {
						size_t n_height = m_block_rows_list[n_row_id].n_height;
						size_t n_index = m_data_pool.index_of((*p_block_it).second);
						n_index += n_width * n_height;
						if(n_min_index > n_index) {
							++ n_dead_block_num;
							n_dead_blocks_size += n_width * n_height;
						}
					}
				}
			}
			for(_TyColumnIter p_col_it = r_dest.m_block_cols_list.begin() +
			   n_block_column_num, p_end_it = r_dest.m_block_cols_list.end();
			   p_col_it != p_end_it; ++ p_col_it) {
				TColumn &r_t_col = *p_col_it;
				size_t n_width = r_t_col.n_width;
				for(_TyBlockIter p_block_it = r_t_col.block_list.begin(), p_block_end_it =
				   r_t_col.block_list.end(); p_block_it != p_block_end_it; ++ p_block_it) {
					size_t n_row_id = (*p_block_it).first;
					size_t n_height = m_block_rows_list[n_row_id].n_height;
					size_t n_index = m_data_pool.index_of((*p_block_it).second);
					n_index += n_width * n_height;
					if(n_min_index > n_index) {
						++ n_dead_block_num;
						n_dead_blocks_size += n_width * n_height;
					}
				}
			}
			// go through all the deleted blocks and see if any of those remain allocated

			/*if(n_dead_block_num) {
				fprintf(stderr, "warning: CUberBlockMatrix::SliceTo() leaves " PRIsize
					" allocated blocks (" PRIsize " elems)\n", n_dead_block_num, n_dead_blocks_size); // debug
			}*/ // silence!
			// warn the user about this operation
#endif // _DEBUG

			m_block_rows_list.resize(n_block_row_num);
			m_block_cols_list.resize(n_block_column_num);
			// just throw awawy leftovers

			_ASSERTE(m_data_pool.size() >= n_min_index); // note that n_min_index points *after* the last block, not at its origin
#if defined(_DEBUG) && defined(__UBER_BLOCK_MATRIX_BUFFER_BOUNDS_CHECK)
			if(n_min_index) { // could end up empty
				n_min_index += CBlockBoundsDebugInitializer<_TyPool>::blockBoundMarker_Element_Num;
				// the pool is not empty, there is a block in there
				// and we want to keep one block bound marker after it.
				// another right thing to do would be to fill the rest of the pool with uninit_page markers
			}
			// account for the bounds check
#endif // _DEBUG && __UBER_BLOCK_MATRIX_BUFFER_BOUNDS_CHECK

			if(pool_MemoryAlignment)
				n_min_index = n_Align_Up(n_min_index, pool_MemoryAlignment / sizeof(double));
			// account for memory alignment

			_ASSERTE(m_data_pool.size() >= n_min_index);
#ifdef _DEBUG
			/*fprintf(stderr, "debug: CUberBlockMatrix::SliceTo() deletes " PRIsize
				" elems)\n", m_data_pool.size() - n_min_index);*/ // debug
#endif // _DEBUG
			m_data_pool.resize(n_min_index);
			// delete stuff from the end of the pool
		}
	}

	//r_dest.CheckIntegrity(true);
}

void CUberBlockMatrix::SliceTo(CUberBlockMatrix &r_dest, size_t n_first_row_index,
	size_t n_last_row_index, size_t n_first_column_index, size_t n_last_column_index,
	bool b_share_data) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	_ASSERTE(&r_dest != this); // can't slice to itself

	_ASSERTE(n_first_column_index < n_last_column_index);
	_ASSERTE(n_first_row_index < n_last_row_index); // note this requires a non-empty matrix
	_ASSERTE(n_last_row_index <= m_block_rows_list.size());
	_ASSERTE(n_last_column_index <= m_block_cols_list.size());
	// make sure the slice indices are sane

	r_dest.Clear();
	// erase dest ...

	r_dest.m_n_row_num = n_BlockRow_Base(n_last_row_index - 1) +
		n_BlockRow_Row_Num(n_last_row_index - 1) - n_BlockRow_Base(n_first_row_index);
	r_dest.m_n_col_num = n_BlockColumn_Base(n_last_column_index - 1) +
		n_BlockColumn_Column_Num(n_last_column_index - 1) - n_BlockColumn_Base(n_first_column_index);
	// set matrix size // note last indices point one past the last, so there must be "- 1"

	r_dest.m_block_rows_list.resize(n_last_row_index - n_first_row_index);
	std::transform(m_block_rows_list.begin() + n_first_row_index,
		m_block_rows_list.begin() + n_last_row_index,
		r_dest.m_block_rows_list.begin(), CRecalcRowCumsum());
	// set block rows (just copy, recalc cumulative sums; this is doing another cumsum, but we might just as well subtract cumsum of the first block, might be slightly faster)

	r_dest.m_block_cols_list.resize(n_last_column_index - n_first_column_index);
	size_t n_base_col_cumsum = m_block_cols_list[n_first_column_index].n_cumulative_width_sum;
	for(size_t i = n_first_column_index; i < n_last_column_index; ++ i) { // todo - make iterator loop instead
		TColumn &r_t_dest = r_dest.m_block_cols_list[i - n_first_column_index];
		const TColumn &r_t_src = m_block_cols_list[i];
		// get src / dest columns

		r_t_dest.n_width = r_t_src.n_width;
		r_t_dest.n_cumulative_width_sum = r_t_src.n_cumulative_width_sum - n_base_col_cumsum;
		// copy/offset column dimensions

		_TyBlockConstIter p_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::lower_bound(r_t_src.block_list.begin(), r_t_src.block_list.end(),
			n_first_row_index, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::lower_bound(r_t_src.block_list.begin(), r_t_src.block_list.end(),
			n_first_row_index, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		_TyBlockConstIter p_end_it = // note it might be faster to just check for the end in iterating the range, especially if the number of blocks is very low
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::lower_bound(p_block_it, r_t_src.block_list.end(),
			n_last_row_index, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::lower_bound(p_block_it, r_t_src.block_list.end(),
			n_last_row_index, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		r_t_dest.block_list.resize(p_end_it - p_block_it);
		_TyBlockIter p_dest_it = r_t_dest.block_list.begin();
		const size_t n_column_width = r_t_src.n_width;
		if(b_share_data) {
			size_t n_height_sum = 0;
			for(; p_block_it != p_end_it; ++ p_block_it, ++ p_dest_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				// get src block

				*p_dest_it = TColumn::TBlockEntry(r_t_block.first -
					n_first_row_index, r_t_block.second);
				// just recalculate row index

				n_height_sum += m_block_rows_list[r_t_block.first].n_height;
			}
			r_dest.m_n_ref_elem_num += r_t_src.n_width * n_height_sum; // !!
		} else {
			for(; p_block_it != p_end_it; ++ p_block_it, ++ p_dest_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				// get src block

				const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
				double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);
				memcpy(p_data, r_t_block.second, n_row_height * n_column_width * sizeof(double));
				// alloc buffer in dest matrix, copy contents

				*p_dest_it = TColumn::TBlockEntry(r_t_block.first -
					n_first_row_index, p_data);
				// recalculate row index and use the new buffer
			}
		}
		_ASSERTE(p_dest_it == r_t_dest.block_list.end());
	}
	// set block columns, copy block data (only blocks in the valid column range)

	//r_dest.CheckIntegrity(true);
	// make sure the dest matrix is ok
}

void CUberBlockMatrix::PermuteTo(CUberBlockMatrix &r_dest, const size_t *p_block_ordering,
	size_t UNUSED(n_ordering_size), bool b_reorder_rows,
	bool b_reorder_columns, bool b_share_data) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	/*b_reorder_rows = true;
	b_reorder_columns = true;*/
	// todo - remove this override

	_ASSERTE(&r_dest != this);
	// can't permute to itself

	_ASSERTE(!b_reorder_rows || n_ordering_size >= m_block_rows_list.size());
	_ASSERTE(!b_reorder_columns || n_ordering_size >= m_block_cols_list.size());
	// make sure that permutation vector is long enough (might actually be longer)

	_ASSERTE(!b_reorder_rows || !b_reorder_columns || b_SymmetricLayout()); // either permutation on one side only, or the matrix is symmetric
	// might fail, b_SymmetricLayout() doesn't cope well with empty blocks ... but then again, permutation does neither

	r_dest.Clear();
	// erase dest ...

	r_dest.m_n_row_num = m_n_row_num;
	r_dest.m_n_col_num = m_n_col_num;
	// set matrix size

	if(b_share_data)
		r_dest.m_n_ref_elem_num = m_data_pool.size() + m_n_ref_elem_num;
	// copy matrix size estimate!

	if(b_reorder_rows) {
		r_dest.m_block_rows_list.resize(m_block_rows_list.size());
		for(size_t i = 0, n = m_block_rows_list.size(); i < n; ++ i) {
			size_t n_dest_index = p_block_ordering[i];
			size_t n_src_index = i;
			r_dest.m_block_rows_list[n_dest_index].n_height =
				m_block_rows_list[n_src_index].n_height;
		}
		for(size_t i = 0, n = m_block_rows_list.size(), n_running_cumsum = 0; i < n; ++ i) {
			size_t n_height = r_dest.m_block_rows_list[i].n_height;
			r_dest.m_block_rows_list[i].n_cumulative_height_sum = n_running_cumsum;
			n_running_cumsum += n_height;
		}
		// copy rows / cols with ordering, refresh cumsums
	} else {
		r_dest.m_block_rows_list.insert(r_dest.m_block_rows_list.begin(),
			m_block_rows_list.begin(), m_block_rows_list.end());
		// just copy
	}
	// set block rows

	r_dest.m_block_cols_list.resize(m_block_cols_list.size());
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		size_t n_src_col, n_dest_col;
		if(b_reorder_columns) {
			n_dest_col = p_block_ordering[i];
			n_src_col = i;
		} else
			n_src_col = n_dest_col = i;
		// calculate ordered indices

		TColumn &r_t_dest = r_dest.m_block_cols_list[n_dest_col];
		const TColumn &r_t_src = m_block_cols_list[n_src_col];
		// get src / dest columns

		const size_t n_column_width = r_t_dest.n_width = r_t_src.n_width;
		// copy/offset column dimensions

		r_t_dest.block_list.resize(r_t_src.block_list.size());
		_TyBlockIter p_dest_it = r_t_dest.block_list.begin();
		_TyBlockConstIter p_block_it = r_t_src.block_list.begin();
		_TyBlockConstIter p_end_it = r_t_src.block_list.end();
		if(b_share_data) {
			if(b_reorder_rows) {
				for(; p_block_it != p_end_it; ++ p_block_it, ++ p_dest_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					// ignore blocks below diagonal

					size_t n_dest_row = p_block_ordering[r_t_block.first];
					// note the permutation is in opposite order

					*p_dest_it = TColumn::TBlockEntry(n_dest_row, r_t_block.second);
				}
			} else {
				for(; p_block_it != p_end_it; ++ p_block_it, ++ p_dest_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					*p_dest_it = TColumn::TBlockEntry(r_t_block.first, r_t_block.second);
				}
			}
		} else {
			if(b_reorder_rows) {
				for(; p_block_it != p_end_it; ++ p_block_it, ++ p_dest_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;

					size_t n_dest_row = p_block_ordering[r_t_block.first];
					// note the permutation is in opposite order

					const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
					_ASSERTE(r_dest.m_block_rows_list[n_dest_row].n_height == n_row_height);
					double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);

					memcpy(p_data, r_t_block.second,
						n_row_height * n_column_width * sizeof(double));
					// alloc buffer in dest matrix, copy contents

					*p_dest_it = TColumn::TBlockEntry(n_dest_row, p_data);
					// recalculate row index and use the new buffer
				}
			} else {
				const size_t n_column_width = r_t_src.n_width;
				for(; p_block_it != p_end_it; ++ p_block_it, ++ p_dest_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					// get src block

					const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
					double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);
					memcpy(p_data, r_t_block.second,
						n_row_height * n_column_width * sizeof(double));
					// alloc buffer in dest matrix, copy contents

					*p_dest_it = TColumn::TBlockEntry(r_t_block.first, p_data);
					// recalculate row index and use the new buffer
				}
			}
		}
		_ASSERTE(p_dest_it == r_t_dest.block_list.end()); // make sure it used all the space
	}
	// set block columns, copy block data (only blocks in the valid column range)

	for(size_t i = 0, n = m_block_cols_list.size(), n_running_width_sum = 0; i < n; ++ i) {
		TColumn &r_t_dest = r_dest.m_block_cols_list[i];
		r_t_dest.n_cumulative_width_sum = n_running_width_sum;
		n_running_width_sum += r_t_dest.n_width;

		if(b_reorder_rows && r_t_dest.block_list.size() > 1)
			std::sort(r_t_dest.block_list.begin(), r_t_dest.block_list.end(), CCompareBlockRow()); // todo - any ideas on how to keep it sorted even after permutation?
		// in case we reordered rows, we need to sort the block list (it is usually short, no worries)
	}
	// fix column cumsums

	//r_dest.CheckIntegrity(true);
	// make sure the dest matrix is ok
}

void CUberBlockMatrix::Permute_UpperTriangular_To(CUberBlockMatrix &r_dest,
	const size_t *p_block_ordering, size_t UNUSED(n_ordering_size), bool b_share_data) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	/*b_reorder_rows = true;
	b_reorder_columns = true;*/
	// todo - remove this override

	_ASSERTE(&r_dest != this);
	// can't permute to itself

	_ASSERTE(n_ordering_size >= m_block_rows_list.size());
	_ASSERTE(n_ordering_size >= m_block_cols_list.size()); // todo - permit >=, or be strict with ==?
	// make sure that permutation vector is long enough (might actually be longer)

	_ASSERTE(b_SymmetricLayout()); // either permutation on one side only, or the matrix is symmetric
	// might fail, b_SymmetricLayout() doesn't cope well with empty blocks ... but then again, permutation does neither

	r_dest.Clear();
	// erase dest ...

	r_dest.m_n_row_num = m_n_row_num;
	r_dest.m_n_col_num = m_n_col_num;
	// set matrix size

	if(b_share_data)
		r_dest.m_n_ref_elem_num = m_data_pool.size() + m_n_ref_elem_num;
	// copy matrix size estimate!

	{
		r_dest.m_block_rows_list.resize(m_block_rows_list.size());
		for(size_t i = 0, n = m_block_rows_list.size(); i < n; ++ i) {
			size_t n_dest_index = p_block_ordering[i];
			size_t n_src_index = i;
			r_dest.m_block_rows_list[n_dest_index].n_height =
				m_block_rows_list[n_src_index].n_height;
		}
		for(size_t i = 0, n = m_block_rows_list.size(), n_running_cumsum = 0; i < n; ++ i) {
			size_t n_height = r_dest.m_block_rows_list[i].n_height;
			r_dest.m_block_rows_list[i].n_cumulative_height_sum = n_running_cumsum;
			n_running_cumsum += n_height;
		}
		// copy rows / cols with ordering, refresh cumsums
	}
	// set block rows

	r_dest.m_block_cols_list.resize(m_block_cols_list.size());
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		size_t n_dest_col = p_block_ordering[i];
		size_t n_src_col = i;
		// calculate ordered indices

		TColumn &r_t_dest = r_dest.m_block_cols_list[n_dest_col];
		const TColumn &r_t_src = m_block_cols_list[n_src_col];
		// get src / dest columns

		const size_t n_column_width = r_t_dest.n_width = r_t_src.n_width;
		// copy/offset column dimensions

		r_t_dest.block_list.resize(r_t_src.block_list.size());
		_TyBlockIter p_dest_it = r_t_dest.block_list.begin();
		_TyBlockConstIter p_block_it = r_t_src.block_list.begin();
		_TyBlockConstIter p_end_it = r_t_src.block_list.end();
		if(b_share_data) {
			for(; p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				if(r_t_block.first > n_src_col)
					continue;
				// ignore blocks below diagonal

				size_t n_dest_row = p_block_ordering[r_t_block.first];
				// note the permutation is in opposite order

				if(n_dest_row <= n_dest_col)
					*p_dest_it = TColumn::TBlockEntry(n_dest_row, r_t_block.second);
				else {
					const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
					_ASSERTE(r_dest.m_block_rows_list[n_dest_row].n_height == n_row_height);
					double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);
					_TyMatrixXdRef src((double*)r_t_block.second, n_row_height, n_column_width), // src is rows, cols
						dest(p_data, n_column_width, n_row_height); // dest is cols, rows
					_ASSERTE(r_dest.m_block_cols_list[n_dest_row].n_width == dest.cols() &&
						r_dest.m_block_rows_list[n_dest_col].n_height == dest.rows()); // make sure that the transpose has correct dimensions
					dest.noalias() = src.transpose();
					TColumn &r_t_dest2 = r_dest.m_block_cols_list[n_dest_row];
					r_t_dest2.block_list.push_back(TColumn::TBlockEntry(n_dest_col, p_data));
					continue; // did not put it there !!
					// in case the block would end below diagonal, we need to make
					// a local transpose of the block
				}

				++ p_dest_it;
			}
		} else {
			for(; p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				if(r_t_block.first > n_src_col)
					continue;
				// ignore blocks below diagonal

				size_t n_dest_row = p_block_ordering[r_t_block.first];
				// note the permutation is in opposite order

				const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
				_ASSERTE(r_dest.m_block_rows_list[n_dest_row].n_height == n_row_height);
				double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);

				if(n_dest_row <= n_dest_col) {
					memcpy(p_data, r_t_block.second,
						n_row_height * n_column_width * sizeof(double));
					// alloc buffer in dest matrix, copy contents

					*p_dest_it = TColumn::TBlockEntry(n_dest_row, p_data);
					// recalculate row index and use the new buffer
				} else {
					_TyMatrixXdRef src((double*)r_t_block.second, n_row_height, n_column_width), // src is rows, cols
						dest(p_data, n_column_width, n_row_height); // dest is cols, rows
					_ASSERTE(r_dest.m_block_cols_list[n_dest_row].n_width == dest.cols() &&
						r_dest.m_block_rows_list[n_dest_col].n_height == dest.rows()); // make sure that the transpose has correct dimensions
					dest.noalias() = src.transpose();
					TColumn &r_t_dest2 = r_dest.m_block_cols_list[n_dest_row];
					r_t_dest2.block_list.push_back(TColumn::TBlockEntry(n_dest_col, p_data));
					continue; // did not put it there !!
					// in case the block would end below diagonal, we need to make
					// a local transpose of the block
				}

				++ p_dest_it;
			}
		}
		r_t_dest.block_list.resize(p_dest_it - r_t_dest.block_list.begin()); // might not use all the space
	}
	// set block columns, copy block data (only blocks in the valid column range)

	for(size_t i = 0, n = m_block_cols_list.size(), n_running_width_sum = 0; i < n; ++ i) {
		TColumn &r_t_dest = r_dest.m_block_cols_list[i];
		r_t_dest.n_cumulative_width_sum = n_running_width_sum;
		n_running_width_sum += r_t_dest.n_width;

		if(r_t_dest.block_list.size() > 1)
			std::sort(r_t_dest.block_list.begin(), r_t_dest.block_list.end(), CCompareBlockRow()); // todo - any ideas on how to keep it sorted even after permutation?
		// in case we reordered rows, we need to sort the block list (it is usually short, no worries)
	}
	// fix column cumsums

	if(b_share_data)
		r_dest.m_n_ref_elem_num -= r_dest.m_data_pool.size();
	// in case data is shared and reflection was needed

	//r_dest.CheckIntegrity(true);
	// make sure the dest matrix is ok
}

void CUberBlockMatrix::Permute_UpperTriangular_To(CUberBlockMatrix &r_dest,
	const size_t *p_block_ordering, size_t UNUSED(n_ordering_size),
	bool b_share_data, size_t n_min_block_row_column,
	bool b_keep_upper_left_part) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	/*b_reorder_rows = true;
	b_reorder_columns = true;*/
	// todo - remove this override

	_ASSERTE(&r_dest != this);
	// can't permute to itself

	_ASSERTE(n_ordering_size >= m_block_rows_list.size());
	_ASSERTE(n_ordering_size >= m_block_cols_list.size()); // todo - permit >=, or be strict with ==?
	// make sure that permutation vector is long enough (might actually be longer)

	_ASSERTE(b_SymmetricLayout()); // either permutation on one side only, or the matrix is symmetric
	// might fail, b_SymmetricLayout() doesn't cope well with empty blocks ... but then again, permutation does neither

	if(b_keep_upper_left_part && n_min_block_row_column > 0) {
		r_dest.SliceTo(r_dest, n_min_block_row_column,
			n_min_block_row_column, true);
		// can only keep the small upper-lef square due to the fact that the permutation is symmetric
	} else
		r_dest.Clear();
	// erase dest ...

	r_dest.m_n_row_num = m_n_row_num;
	r_dest.m_n_col_num = m_n_col_num;
	// set matrix size

	if(b_share_data)
		r_dest.m_n_ref_elem_num = m_data_pool.size() + m_n_ref_elem_num;
	// copy matrix size estimate!

	{
		r_dest.m_block_rows_list.resize(m_block_rows_list.size());
		if(b_keep_upper_left_part) {
			for(size_t i = 0, n = m_block_rows_list.size(); i < n; ++ i) {
				size_t n_dest_index = p_block_ordering[i];
				size_t n_src_index = i;
				if(n_dest_index < n_min_block_row_column) {
					_ASSERTE(r_dest.m_block_rows_list[n_dest_index].n_height ==
						m_block_rows_list[n_src_index].n_height);
					// make sure that the prefix does not change
				} else {
					r_dest.m_block_rows_list[n_dest_index].n_height =
						m_block_rows_list[n_src_index].n_height;
				}
			}
		} else {
			for(size_t i = 0, n = m_block_rows_list.size(); i < n; ++ i) {
				size_t n_dest_index = p_block_ordering[i];
				size_t n_src_index = i;
				r_dest.m_block_rows_list[n_dest_index].n_height =
					m_block_rows_list[n_src_index].n_height;
			}
		}
		for(size_t i = 0, n = m_block_rows_list.size(), n_running_cumsum = 0; i < n; ++ i) {
			size_t n_height = r_dest.m_block_rows_list[i].n_height;
			r_dest.m_block_rows_list[i].n_cumulative_height_sum = n_running_cumsum;
			n_running_cumsum += n_height;
		}
		// copy rows / cols with ordering, refresh cumsums
	}
	// set block rows

	size_t n_min_block_row = (b_keep_upper_left_part)? 0 : n_min_block_row_column;
	// in case we wanted to keep the upper part, we need to reperm anyway as it is going to be shuffled

	r_dest.m_block_cols_list.resize(m_block_cols_list.size());
	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) {
		size_t n_dest_col = p_block_ordering[i];
		size_t n_src_col = i;
		// calculate ordered indices

		TColumn &r_t_dest = r_dest.m_block_cols_list[n_dest_col];
		const TColumn &r_t_src = m_block_cols_list[n_src_col];
		// get src / dest columns

		_ASSERTE(!b_keep_upper_left_part || n_dest_col >= n_min_block_row_column ||
			r_t_dest.n_width == r_t_src.n_width); // the part that was kept must be the same
		const size_t n_column_width = r_t_dest.n_width = r_t_src.n_width;
		// copy/offset column dimensions

		if(n_dest_col < n_min_block_row_column) {
			if(!b_keep_upper_left_part)
				continue;
			// usually we don't want to do anything

			_TyBlockConstIter p_end_it = r_t_src.block_list.end();
			_TyBlockConstIter p_block_it = r_t_src.block_list.begin();
			if(b_share_data) {
				for(; p_block_it != p_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					if(r_t_block.first > n_src_col)
						continue;
					// ignore blocks below diagonal

					size_t n_dest_row = p_block_ordering[r_t_block.first];
					// note the permutation is in opposite order

					_ASSERTE(n_dest_row >= n_min_block_row); // since b_keep_upper_left_part, n_dest_row = 0
					/*if(n_dest_row < n_min_block_row)
						continue;*/
					// ignore blocks below threshold

					if(n_dest_row >= n_min_block_row_column) {
						_ASSERTE(n_dest_row > n_dest_col);
						const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
						_ASSERTE(r_dest.m_block_rows_list[n_dest_row].n_height == n_row_height);
						double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);
						_TyMatrixXdRef src((double*)r_t_block.second, n_row_height, n_column_width), // src is rows, cols
							dest(p_data, n_column_width, n_row_height); // dest is cols, rows
						_ASSERTE(r_dest.m_block_cols_list[n_dest_row].n_width == dest.cols() &&
							r_dest.m_block_rows_list[n_dest_col].n_height == dest.rows()); // make sure that the transpose has correct dimensions
						dest.noalias() = src.transpose();
						TColumn &r_t_dest2 = r_dest.m_block_cols_list[n_dest_row];
						r_t_dest2.block_list.push_back(TColumn::TBlockEntry(n_dest_col, p_data));
						// in case the block would end below diagonal, we need to make
						// a local transpose of the block
					}
				}
			} else {
				for(; p_block_it != p_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					if(r_t_block.first > n_src_col)
						continue;
					// ignore blocks below diagonal

					size_t n_dest_row = p_block_ordering[r_t_block.first];
					// note the permutation is in opposite order

					_ASSERTE(n_dest_row >= n_min_block_row); // since b_keep_upper_left_part, n_dest_row = 0
					/*if(n_dest_row < n_min_block_row)
						continue;*/
					// ignore blocks below threshold

					if(n_dest_row >= n_min_block_row_column) {
						_ASSERTE(n_dest_row > n_dest_col);
						const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
						_ASSERTE(r_dest.m_block_rows_list[n_dest_row].n_height == n_row_height);
						double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);
						_ASSERTE(n_dest_row > n_dest_col);
						_TyMatrixXdRef src((double*)r_t_block.second, n_row_height, n_column_width), // src is rows, cols
							dest(p_data, n_column_width, n_row_height); // dest is cols, rows
						_ASSERTE(r_dest.m_block_cols_list[n_dest_row].n_width == dest.cols() &&
							r_dest.m_block_rows_list[n_dest_col].n_height == dest.rows()); // make sure that the transpose has correct dimensions
						dest.noalias() = src.transpose();
						TColumn &r_t_dest2 = r_dest.m_block_cols_list[n_dest_row];
						r_t_dest2.block_list.push_back(TColumn::TBlockEntry(n_dest_col, p_data));
						// in case the block would end below diagonal, we need to make
						// a local transpose of the block
					}
				}
			}
			// we want to fill some entries that would otherwise end up discarded
			// this is tricky since below-diagonal entries will get transposed to other places
			// so if min row is 0, min col must also be 0 (but only for the transpose blocks)

			continue;
			// and we're done
		}
		// skip some columns

		_TyBlockConstIter p_end_it = r_t_src.block_list.end();
		_TyBlockConstIter p_block_it = r_t_src.block_list.begin();
/*#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::lower_bound(r_t_src.block_list.begin(),
			p_end_it, n_min_block_row_column, CUberBlockMatrix::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::lower_bound(r_t_src.block_list.begin(),
			p_end_it, n_min_block_row_column, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400*/ // can't do this, there's a permutation in there
		size_t n_off = (b_keep_upper_left_part)? r_t_dest.block_list.size() : 0;
		r_t_dest.block_list.resize((p_end_it - p_block_it) + n_off);
		_TyBlockIter p_dest_it = (b_keep_upper_left_part)?
			r_t_dest.block_list.begin() + n_off : r_t_dest.block_list.begin();
		// see how much blocks there are in the selected range

		if(b_share_data) {
			for(; p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				if(r_t_block.first > n_src_col)
					continue;
				// ignore blocks below diagonal

				size_t n_dest_row = p_block_ordering[r_t_block.first];
				// note the permutation is in opposite order

				//if(!b_keep_upper_left_part && n_dest_row < n_min_block_row_column)
				if(n_dest_row < n_min_block_row)
					continue;
				// ignore blocks below threshold

				if(n_dest_row <= n_dest_col)
					*p_dest_it = TColumn::TBlockEntry(n_dest_row, r_t_block.second);
				else {
					const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
					_ASSERTE(r_dest.m_block_rows_list[n_dest_row].n_height == n_row_height);
					double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);
					_TyMatrixXdRef src((double*)r_t_block.second, n_row_height, n_column_width), // src is rows, cols
						dest(p_data, n_column_width, n_row_height); // dest is cols, rows
					_ASSERTE(r_dest.m_block_cols_list[n_dest_row].n_width == dest.cols() &&
						r_dest.m_block_rows_list[n_dest_col].n_height == dest.rows()); // make sure that the transpose has correct dimensions
					dest.noalias() = src.transpose();
					TColumn &r_t_dest2 = r_dest.m_block_cols_list[n_dest_row];
					r_t_dest2.block_list.push_back(TColumn::TBlockEntry(n_dest_col, p_data));
					continue; // did not put it there !!
					// in case the block would end below diagonal, we need to make
					// a local transpose of the block
				}

				++ p_dest_it;
			}
		} else {
			for(; p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				if(r_t_block.first > n_src_col)
					continue;
				// ignore blocks below diagonal

				size_t n_dest_row = p_block_ordering[r_t_block.first];
				// note the permutation is in opposite order

				//if(!b_keep_upper_left_part && n_dest_row < n_min_block_row_column)
				if(n_dest_row < n_min_block_row)
					continue;
				// ignore blocks below threshold

				const size_t n_row_height = m_block_rows_list[r_t_block.first].n_height;
				_ASSERTE(r_dest.m_block_rows_list[n_dest_row].n_height == n_row_height);
				double *p_data = r_dest.p_Get_DenseStorage(n_row_height * n_column_width);

				if(n_dest_row <= n_dest_col) {
					memcpy(p_data, r_t_block.second,
						n_row_height * n_column_width * sizeof(double));
					// alloc buffer in dest matrix, copy contents

					*p_dest_it = TColumn::TBlockEntry(n_dest_row, p_data);
					// recalculate row index and use the new buffer
				} else {
					_TyMatrixXdRef src((double*)r_t_block.second, n_row_height, n_column_width), // src is rows, cols
						dest(p_data, n_column_width, n_row_height); // dest is cols, rows
					_ASSERTE(r_dest.m_block_cols_list[n_dest_row].n_width == dest.cols() &&
						r_dest.m_block_rows_list[n_dest_col].n_height == dest.rows()); // make sure that the transpose has correct dimensions
					dest.noalias() = src.transpose();
					TColumn &r_t_dest2 = r_dest.m_block_cols_list[n_dest_row];
					r_t_dest2.block_list.push_back(TColumn::TBlockEntry(n_dest_col, p_data));
					continue; // did not put it there !!
					// in case the block would end below diagonal, we need to make
					// a local transpose of the block
				}

				++ p_dest_it;
			}
		}
		r_t_dest.block_list.resize(p_dest_it - r_t_dest.block_list.begin()); // might not use all the space
	}
	// set block columns, copy block data (only blocks in the valid column range)

	for(size_t i = 0, n = m_block_cols_list.size(), n_running_width_sum = 0; i < n; ++ i) {
		TColumn &r_t_dest = r_dest.m_block_cols_list[i];
		r_t_dest.n_cumulative_width_sum = n_running_width_sum;
		n_running_width_sum += r_t_dest.n_width;

		if(r_t_dest.block_list.size() > 1)
			std::sort(r_t_dest.block_list.begin(), r_t_dest.block_list.end(), CCompareBlockRow()); // todo - any ideas on how to keep it sorted even after permutation?
		// in case we reordered rows, we need to sort the block list (it is usually short, no worries)
	}
	// fix column cumsums

	if(b_share_data)
		r_dest.m_n_ref_elem_num -= r_dest.m_data_pool.size();
	// in case data is shared and reflection was needed

	r_dest.CheckIntegrity(true);
	// make sure the dest matrix is ok
}

bool CUberBlockMatrix::UpperTriangularTranspose_Solve(double *p_x, size_t UNUSED(n_vector_size)) const // todo add support for block permutation (2nd function), should be faster than permutating x there and back (measure, measure)
{
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

		const TColumn::TBlockEntry &r_t_diag_block = r_t_col.block_list.back();
		const TRow &r_t_diag_row = m_block_rows_list[r_t_diag_block.first];
		if(r_t_col.n_width != r_t_diag_row.n_height ||
		   r_t_col.n_cumulative_width_sum != r_t_diag_row.n_cumulative_height_sum)
			return false;
		const size_t n_width = r_t_col.n_width, n_x = r_t_col.n_cumulative_width_sum;
		// make sure that the last block is on the diagonal, and is square (that does happen in L)

		_TyVectorXdRef diag_x(p_x + n_x, n_width); // note this is always unaligned
		// part of the x vector, corresponding to this column

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_block_end_it = r_t_col.block_list.end() - 1;
		   p_block_it != p_block_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const size_t n_height = r_t_row.n_height, n_y = r_t_row.n_cumulative_height_sum;

			_ASSERTE(n_y + n_height <= n_x); // make sure that the diagonal part of the vector is not affected here

			_TyRowVectorXdRef x_part(p_x + n_y, n_height); // also column vector, also unaligned
			_TyMatrixXdRef block(r_t_block.second, n_height, n_width); // block

#if 0
			for(size_t i = 0; i < n_width; ++ i) {
				for(size_t j = 0; j < n_height; ++ j)
					diag_x(i) -= x_part(j) * block(j, i);
			}
#else // 0
#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			diag_x -= x_part.lazyProduct(block);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			diag_x.noalias() -= x_part * block; // do the same as the loop above, using matrix multiplication
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
#endif // 0
			// calculate x -= column(j) * x(j) in parallel for all columns in this block
		}
		// for every other block in this column (note no reverse order here)
		// note that this loop mostly executes once (optimize for that? how?)

		_TyMatrixXdRef diag_block(r_t_diag_block.second, n_width, n_width);
		for(size_t i = 0; i < n_width; ++ i) {
			double f_old_x = diag_x(i);
			for(size_t j = 0; j < i; ++ j) // elements strictly above diagonal
				f_old_x -= diag_block(j, i) * diag_x(j);
			// forward substitute the diagonal block

#ifdef _DEBUG
			for(size_t j = i + 1; j < n_width; ++ j) // elements strictly under diagonal
				_ASSERTE(diag_block(j, i) == 0);
			// make sure there are only nulls under the diagonal
#endif // _DEBUG

			double f_diag = diag_block(i, i);
			if(f_diag == 0) // what about negative zero? does it cover it? yes.
				return false; // note that if diag_x(i) is zero as well, it probably means "infinite number of solutions", and not "no solution". could still produce one of them.
			diag_x(i) = f_old_x / f_diag;
		}
		// resolve values of the diagonal x
	}

	return true;
}

bool CUberBlockMatrix::UpperTriangularTranspose_Solve(double *p_x,
	size_t UNUSED(n_vector_size), size_t n_skip_columns) const
{
	CheckIntegrity(true);

	_ASSERTE(b_Square());
	// triangular is a special case of square

	_ASSERTE(n_vector_size == m_n_col_num);
	// make sure that the vector's got correct size

	_ASSERTE(n_skip_columns <= m_block_cols_list.size());
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin() +
	   n_skip_columns, p_end_it = m_block_cols_list.end(); p_col_it != p_end_it; ++ p_col_it) { // forward substitution
		const TColumn &r_t_col = *p_col_it;
		if(r_t_col.block_list.empty())
			return false; // a 0 on the diagonal - no solution (is it? csparse would divide by zero and produce a NaN)

		const TColumn::TBlockEntry &r_t_diag_block = r_t_col.block_list.back();
		const TRow &r_t_diag_row = m_block_rows_list[r_t_diag_block.first];
		if(r_t_col.n_width != r_t_diag_row.n_height ||
		   r_t_col.n_cumulative_width_sum != r_t_diag_row.n_cumulative_height_sum)
			return false;
		const size_t n_width = r_t_col.n_width, n_x = r_t_col.n_cumulative_width_sum;
		// make sure that the last block is on the diagonal, and is square (that does happen in L)

		_TyVectorXdRef diag_x(p_x + n_x, n_width); // note this is always unaligned
		// part of the x vector, corresponding to this column

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_block_end_it = r_t_col.block_list.end() - 1;
		   p_block_it != p_block_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const size_t n_height = r_t_row.n_height, n_y = r_t_row.n_cumulative_height_sum;

			_ASSERTE(n_y + n_height <= n_x); // make sure that the diagonal part of the vector is not affected here

			_TyRowVectorXdRef x_part(p_x + n_y, n_height); // also column vector, also unaligned
			_TyMatrixXdRef block(r_t_block.second, n_height, n_width); // block

#if 0
			for(size_t i = 0; i < n_width; ++ i) {
				for(size_t j = 0; j < n_height; ++ j)
					diag_x(i) -= x_part(j) * block(j, i);
			}
#else // 0
#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			diag_x -= x_part.lazyProduct(block);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			diag_x.noalias() -= x_part * block; // do the same as the loop above, using matrix multiplication
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
#endif // 0
			// calculate x -= column(j) * x(j) in parallel for all columns in this block
		}
		// for every other block in this column (note no reverse order here)
		// note that this loop mostly executes once (optimize for that? how?)

		_TyMatrixXdRef diag_block(r_t_diag_block.second, n_width, n_width);
		for(size_t i = 0; i < n_width; ++ i) {
			double f_old_x = diag_x(i);
			for(size_t j = 0; j < i; ++ j) // elements strictly above diagonal
				f_old_x -= diag_block(j, i) * diag_x(j);
			// forward substitute the diagonal block

#ifdef _DEBUG
			for(size_t j = i + 1; j < n_width; ++ j) // elements strictly under diagonal
				_ASSERTE(diag_block(j, i) == 0);
			// make sure there are only nulls under the diagonal
#endif // _DEBUG

			double f_diag = diag_block(i, i);
			if(f_diag == 0) // what about negative zero? does it cover it? yes.
				return false; // note that if diag_x(i) is zero as well, it probably means "infinite number of solutions", and not "no solution". could still produce one of them.
			diag_x(i) = f_old_x / f_diag;
		}
		// resolve values of the diagonal x
	}

	return true;
}

#if 0
bool CUberBlockMatrix::UpperTriangular_Solve(double *p_x, size_t UNUSED(n_vector_size))
{
	CheckIntegrity(true);

	_ASSERTE(b_Square());
	// triangular is a special case of square

	//_ASSERTE(b_SymmetricLayout()); // note that symmetric layout is not really required, it just must hold that the last block of every column is the diagonal block
	// make sure that the layout is symmetric (this is an optimization, not a prerequisite given by math)

	_ASSERTE(n_vector_size == m_n_col_num);
	// make sure that the vector's got correct size

	for(_TyColumnConstIter p_col_it = m_block_cols_list.end(),
	   p_end_it = m_block_cols_list.begin(); p_col_it != p_end_it;) { // back substitution
		-- p_col_it;
		// decrement at the beginning of the loop! (but after the comparison)

		const TColumn &r_t_col = *p_col_it;
		if(r_t_col.block_list.empty())
			return false; // a 0 on the diagonal - no solution (is it? csparse would divide by zero and produce a NaN)
		//_ASSERTE(r_t_col.block_list.back().first == n_col); // makes sure that the matrix really is upper diagonal

		const TColumn::TBlockEntry &r_t_diag_block = r_t_col.block_list.back();
		const TRow &r_t_diag_row = m_block_rows_list[r_t_diag_block.first];
		if(r_t_col.n_width != r_t_diag_row.n_height ||
		   r_t_col.n_cumulative_width_sum != r_t_diag_row.n_cumulative_height_sum)
			return false;
		const size_t n_width = r_t_col.n_width, n_x = r_t_col.n_cumulative_width_sum;
		// make sure that the last block is on the diagonal, and is square (that does happen in L)

		_TyVectorXdRef diag_x(p_x + n_x, n_width); // note this is always unaligned
		_TyMatrixXdRef diag_block(r_t_diag_block.second, n_width, n_width);
		for(size_t i = n_width; i > 0;) {
			-- i; // !!

			double f_diag = diag_block(i, i);
			if(f_diag == 0) // what about negative zero? does it cover it? yes.
				return false; // note that if diag_x(i) is zero as well, it probably means "infinite number of solutions", and not "no solution". could still produce one of them.
			double f_new_x = (diag_x(i) /= f_diag); // f_new_x = diag_x(i)
			for(size_t j = 0; j < i; ++ j) // elements strictly above diagonal
				diag_x(j) -= diag_block(j, i) * f_new_x;
			// backsubstitute the diagonal block

#ifdef _DEBUG
			for(size_t j = i + 1; j < n_width; ++ j) // elements strictly under diagonal
				_ASSERTE(diag_block(j, i) == 0);
			// make sure there are only nulls under the diagonal
#endif // _DEBUG
		}
		// resolve values of the diagonal x

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_block_end_it = r_t_col.block_list.end() - 1;
		   p_block_it != p_block_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const size_t n_height = r_t_row.n_height, n_y = r_t_row.n_cumulative_height_sum;

			_ASSERTE(n_y + n_height <= n_x); // make sure that the diagonal part of the vector is not affected here

			_TyVectorXdRef x_part(p_x + n_y, n_height); // also column vector, also unaligned
			_TyMatrixXdRef block(r_t_block.second, n_height, n_width); // block

#if 0
			for(size_t i = n_width; i > 0;) { // even reversing loop direction gives numerical discrepancies (big dynamic range in the matrix / vector)
				-- i; // !!
				for(size_t j = 0; j < n_height; ++ j)
					x_part(j) -= block(j, i) * diag_x(i);
			}
#else // 0
#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			x_part -= block.lazyProduct(diag_x);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			x_part.noalias() -= block * diag_x; // do the same as the loop above, using matrix multiplication
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
#endif // 0
			// calculate x -= column(j) * x(j) in parallel for all columns in this block
		}
		// for every other block in this column (note no reverse order here)
		// note that this loop mostly executes once (optimize for that? how?)
	}
}
#endif // 0

bool CUberBlockMatrix::UpperTriangular_Solve(double *p_x,
	size_t UNUSED(n_vector_size), size_t n_first_column) const
{
	CheckIntegrity(true);

	_ASSERTE(b_Square());
	// triangular is a special case of square

	//_ASSERTE(b_SymmetricLayout()); // note that symmetric layout is not really required, it just must hold that the last block of every column is the diagonal block
	// make sure that the layout is symmetric (this is an optimization, not a prerequisite given by math)

	_ASSERTE(n_vector_size == m_n_col_num);
	// make sure that the vector's got correct size

	_ASSERTE(n_first_column < m_block_cols_list.size());
	// should be a valid block index

	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin() + (n_first_column + 1),
	   p_end_it = m_block_cols_list.begin(); p_col_it != p_end_it;) { // back substitution
		-- p_col_it;
		// decrement at the beginning of the loop! (but after the comparison)

		const TColumn &r_t_col = *p_col_it;
		if(r_t_col.block_list.empty())
			return false; // a 0 on the diagonal - no solution (is it? csparse would divide by zero and produce a NaN)
		//_ASSERTE(r_t_col.block_list.back().first == n_col); // makes sure that the matrix really is upper diagonal

		const TColumn::TBlockEntry &r_t_diag_block = r_t_col.block_list.back();
		const TRow &r_t_diag_row = m_block_rows_list[r_t_diag_block.first];
		if(r_t_col.n_width != r_t_diag_row.n_height ||
		   r_t_col.n_cumulative_width_sum != r_t_diag_row.n_cumulative_height_sum)
			return false;
		const size_t n_width = r_t_col.n_width, n_x = r_t_col.n_cumulative_width_sum;
		// make sure that the last block is on the diagonal, and is square (that does happen in L)

		_TyVectorXdRef diag_x(p_x + n_x, n_width); // note this is always unaligned
		_TyMatrixXdRef diag_block(r_t_diag_block.second, n_width, n_width);
		for(size_t i = n_width; i > 0;) {
			-- i; // !!

			double f_diag = diag_block(i, i);
			if(f_diag == 0) // what about negative zero? does it cover it? yes.
				return false; // note that if diag_x(i) is zero as well, it probably means "infinite number of solutions", and not "no solution". could still produce one of them.
			double f_new_x = (diag_x(i) /= f_diag); // f_new_x = diag_x(i)
			for(size_t j = 0; j < i; ++ j) // elements strictly above diagonal
				diag_x(j) -= diag_block(j, i) * f_new_x;
			// backsubstitute the diagonal block

#ifdef _DEBUG
			for(size_t j = i + 1; j < n_width; ++ j) // elements strictly under diagonal
				_ASSERTE(diag_block(j, i) == 0);
			// make sure there are only nulls under the diagonal
#endif // _DEBUG
		}
		// resolve values of the diagonal x

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_block_end_it = r_t_col.block_list.end() - 1;
		   p_block_it != p_block_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			const TRow &r_t_row = m_block_rows_list[r_t_block.first];
			const size_t n_height = r_t_row.n_height, n_y = r_t_row.n_cumulative_height_sum;

			_ASSERTE(n_y + n_height <= n_x); // make sure that the diagonal part of the vector is not affected here

			_TyVectorXdRef x_part(p_x + n_y, n_height); // also column vector, also unaligned
			_TyMatrixXdRef block(r_t_block.second, n_height, n_width); // block

#if 0
			for(size_t i = n_width; i > 0;) { // even reversing loop direction gives numerical discrepancies (big dynamic range in the matrix / vector)
				-- i; // !!
				for(size_t j = 0; j < n_height; ++ j)
					x_part(j) -= block(j, i) * diag_x(i);
			}
#else // 0
#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			x_part -= block.lazyProduct(diag_x);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
			x_part.noalias() -= block * diag_x; // do the same as the loop above, using matrix multiplication
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
#endif // 0
			// calculate x -= column(j) * x(j) in parallel for all columns in this block
		}
		// for every other block in this column (note no reverse order here)
		// note that this loop mostly executes once (optimize for that? how?)
	}

	return true;
}

bool CUberBlockMatrix::From_Matrix(size_t n_base_row_id, size_t n_base_column_id,
	const CUberBlockMatrix &r_matrix, bool b_allow_layout_extension) // throw(std::bad_alloc)
{
	r_matrix.CheckIntegrity(true);

	_ASSERTE(&r_matrix != this);
	// can't copy from self

	if(r_matrix.m_block_rows_list.size() + n_base_row_id > m_block_rows_list.size() ||
	   r_matrix.m_block_cols_list.size() + n_base_column_id > m_block_cols_list.size()) {
		if(!b_allow_layout_extension)
			return false; // too many rows / columns
		else {
			throw std::runtime_error("CUberBlockMatrix::From_Matrix()"
				" layout extension not implemented");
			// todo - see if the layout needs to be extended and do it
		}
	}
	// make sure the sparse data fit inside the existing layout

	// todo explicitly test row layout, it cannot be reliably and efficiently tested inside the loop
	// also - should we support subdivision of empty rows? it would make sense in the
	// b_allow_layout_extension context and it would probably be quite simple to implement

	for(size_t i = 0, n = r_matrix.m_block_cols_list.size(); i < n; ++ i) {
		const TColumn &r_t_src_col = r_matrix.m_block_cols_list[i];
		TColumn &r_t_dest_col = m_block_cols_list[i + n_base_column_id];
		if(r_t_src_col.n_width != r_t_dest_col.n_width)
			return false;
		// get source / dest column and make sure the size matches

		if(r_t_src_col.block_list.empty())
			continue;
		// skip empties

		const size_t n_width = r_t_src_col.n_width; // or dest, these are equal

		size_t n_insert_row = r_t_src_col.block_list.front().first + n_base_row_id;
		// see at which row to insert the first block

		_TyBlockIter p_insert_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			std::lower_bound(r_t_dest_col.block_list.begin(),
			r_t_dest_col.block_list.end(), n_insert_row,
			CUberBlockMatrix::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			std::lower_bound(r_t_dest_col.block_list.begin(),
			r_t_dest_col.block_list.end(), n_insert_row, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		// find where to put the block in the column

		_TyBlockConstIter p_block_it =
			r_t_src_col.block_list.begin(), p_end_it = r_t_src_col.block_list.end();
		if(p_insert_it != r_t_dest_col.block_list.end()) {
			size_t n_next_dest_row = (*p_insert_it).first;
			for(;p_block_it != p_end_it; ++ p_block_it) {
				const TColumn::TBlockEntry &r_t_block = *p_block_it;
				// get a block

				//size_t n_src_row = r_t_block.first;
				size_t n_dest_row = r_t_block.first + n_base_row_id;
				_ASSERTE(r_matrix.m_block_rows_list[/*n_src_row*/r_t_block.first].n_height ==
					m_block_rows_list[n_dest_row].n_height);
				// get src / dest row

				size_t n_height = m_block_rows_list[n_dest_row].n_height;

				_ASSERTE(p_insert_it != r_t_dest_col.block_list.end());
				if(n_next_dest_row < n_dest_row) {
					p_insert_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						std::lower_bound(p_insert_it,
						r_t_dest_col.block_list.end(), n_dest_row,
						CUberBlockMatrix::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(p_insert_it,
						r_t_dest_col.block_list.end(), n_dest_row, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					// find where to put the block in the column

					if(p_insert_it == r_t_dest_col.block_list.end())
						break; // this block not added (add to the end in the next loop)
					n_next_dest_row = (*p_insert_it).first;
				}
				_ASSERTE(n_next_dest_row >= n_dest_row); // it must always be after
				if(n_next_dest_row > n_dest_row) {
					size_t n_insert_it = p_insert_it - r_t_dest_col.block_list.begin();
					double *p_data = p_Get_DenseStorage(n_width * n_height);
					memcpy(p_data, r_t_block.second, n_width * n_height * sizeof(double));
					r_t_dest_col.block_list.insert(p_insert_it,
						TColumn::TBlockEntry(n_dest_row, p_data));
					p_insert_it = r_t_dest_col.block_list.begin() + ++ n_insert_it;
					// insert a new block before the row, mind pointer invalidation
				} else {
					_ASSERTE(n_next_dest_row == n_dest_row);
					double *p_data = (*p_insert_it).second;
					memcpy(p_data, r_t_block.second, n_width * n_height * sizeof(double));
					if(++ p_insert_it == r_t_dest_col.block_list.end()) { // is ++ enough to be >= ? no.
						++ p_block_it; // we've just added this block
						break;
					}
					n_next_dest_row = (*p_insert_it).first;
					// copy data to an existing block
				}
			}
		}
		// insert between existing blocks

		for(;p_block_it != p_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			// get a block

			//size_t n_src_row = r_t_block.first;
			size_t n_dest_row = r_t_block.first + n_base_row_id;
			_ASSERTE(r_matrix.m_block_rows_list[/*n_src_row*/r_t_block.first].n_height ==
				m_block_rows_list[n_dest_row].n_height);
			// get src / dest row

			size_t n_height = m_block_rows_list[n_dest_row].n_height;

			double *p_data = p_Get_DenseStorage(n_width * n_height);
			memcpy(p_data, r_t_block.second, n_width * n_height * sizeof(double));
			r_t_dest_col.block_list.push_back(TColumn::TBlockEntry(n_dest_row, p_data));
			// copy data, add to the end
		}
		// insert the new blocks at the end
	}

	//CheckIntegrity(true);

	return true;
}

void CUberBlockMatrix::Convert_to_Dense(Eigen::MatrixXd &r_dest) const // throw(std::bad_alloc)
{
	CheckIntegrity(true);

	r_dest.resize(m_n_row_num, m_n_col_num);
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

double CUberBlockMatrix::f_Norm() const
{
	CheckIntegrity(true);

	double f_norm = 0;
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it; ++ p_col_it) {
		const TColumn &r_t_col = *p_col_it;

		for(_TyBlockConstIter p_block_it =
		   r_t_col.block_list.begin(), p_block_end_it = r_t_col.block_list.end();
		   p_block_it != p_block_end_it; ++ p_block_it) {
			const TColumn::TBlockEntry &r_t_block = *p_block_it;
			size_t n_row_id = r_t_block.first; // get row of the block

			_TyMatrixXdRef block(r_t_block.second,
				m_block_rows_list[n_row_id].n_height, r_t_col.n_width);
			// create map to block data

			f_norm += block.squaredNorm();
		}
	}
	return sqrt(f_norm);
}

void CUberBlockMatrix::Permute_RightHandSide_Vector(double *p_dest, const double *p_src,
	size_t UNUSED(n_vector_length), const size_t *p_permutation,
	size_t UNUSED(n_permutation_length)) const
{
	CheckIntegrity();

	_ASSERTE(p_permutation); // permutation can't be null (csparse allows that, cs_pvec then acts like memcpy)
	_ASSERTE(p_dest != p_src); // can't permute inplace
	_ASSERTE(n_vector_length == m_n_row_num);
	_ASSERTE(n_permutation_length == m_block_rows_list.size());

	for(size_t i = 0, n = m_block_rows_list.size(); i < n; ++ i) { // note this could run in parallel
		size_t p = p_permutation[i];
		_ASSERTE(p < n); // makes sure the permutation contains valid block indices
		size_t n_block_base = m_block_rows_list[p].n_cumulative_height_sum;
		size_t n_block_size = m_block_rows_list[p].n_height;
		_ASSERTE(n_block_size > 0); // can't contain zero-size blocks anyway
		do {
			*p_dest = p_src[n_block_base];
			++ p_dest;
			++ n_block_base;
		} while(-- n_block_size);
	}
}

void CUberBlockMatrix::Permute_LeftHandSide_Vector(double *p_dest, const double *p_src,
	size_t UNUSED(n_vector_length), const size_t *p_permutation,
	size_t UNUSED(n_permutation_length)) const
{
	CheckIntegrity();

	_ASSERTE(p_permutation); // permutation can't be null (csparse allows that, cs_pvec then acts like memcpy)
	_ASSERTE(p_dest != p_src); // can't permute inplace
	_ASSERTE(n_vector_length == m_n_col_num);
	_ASSERTE(n_permutation_length == m_block_cols_list.size());

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) { // note this could run in parallel
		size_t p = p_permutation[i];
		_ASSERTE(p < n); // makes sure the permutation contains valid block indices
		size_t n_block_base = m_block_cols_list[p].n_cumulative_width_sum;
		size_t n_block_size = m_block_cols_list[p].n_width;
		_ASSERTE(n_block_size > 0); // can't contain zero-size blocks anyway
		do {
			*p_dest = p_src[n_block_base];
			++ p_dest;
			++ n_block_base;
		} while(-- n_block_size);
	}
}

void CUberBlockMatrix::InversePermute_RightHandSide_Vector(double *p_dest, const double *p_src,
	size_t UNUSED(n_vector_length), const size_t *p_permutation,
	size_t UNUSED(n_permutation_length)) const
{
	CheckIntegrity();

	_ASSERTE(p_permutation); // permutation can't be null (csparse allows that, cs_pvec then acts like memcpy)
	_ASSERTE(p_dest != p_src); // can't permute inplace
	_ASSERTE(n_vector_length == m_n_row_num);
	_ASSERTE(n_permutation_length == m_block_rows_list.size());

	for(size_t i = 0, n = m_block_rows_list.size(); i < n; ++ i) { // note this could run in parallel
		size_t p = p_permutation[i];
		_ASSERTE(p < n); // makes sure the permutation contains valid block indices
		size_t n_block_base = m_block_rows_list[p].n_cumulative_height_sum;
		size_t n_block_size = m_block_rows_list[p].n_height;
		_ASSERTE(n_block_size > 0); // can't contain zero-size blocks anyway
		do {
			p_dest[n_block_base] = *p_src;
			++ p_src;
			++ n_block_base;
		} while(-- n_block_size);
	}
}

void CUberBlockMatrix::InversePermute_LeftHandSide_Vector(double *p_dest, const double *p_src,
	size_t UNUSED(n_vector_length), const size_t *p_permutation,
	size_t UNUSED(n_permutation_length)) const
{
	CheckIntegrity();

	_ASSERTE(p_permutation); // permutation can't be null (csparse allows that, cs_pvec then acts like memcpy)
	_ASSERTE(p_dest != p_src); // can't permute inplace
	_ASSERTE(n_vector_length == m_n_col_num);
	_ASSERTE(n_permutation_length == m_block_cols_list.size());

	for(size_t i = 0, n = m_block_cols_list.size(); i < n; ++ i) { // note this could run in parallel
		size_t p = p_permutation[i];
		_ASSERTE(p < n); // makes sure the permutation contains valid block indices
		size_t n_block_base = m_block_cols_list[p].n_cumulative_width_sum;
		size_t n_block_size = m_block_cols_list[p].n_width;
		_ASSERTE(n_block_size > 0); // can't contain zero-size blocks anyway
		do {
			p_dest[n_block_base] = *p_src;
			++ p_src;
			++ n_block_base;
		} while(-- n_block_size);
	}
}

void CUberBlockMatrix::Build_EliminationTree(std::vector<size_t> &r_elim_tree,
	std::vector<size_t> &r_workspace) const // throw(std::bad_alloc)
{
	//CheckIntegrity(); // inside b_SymmetricLayout()
	_ASSERTE(b_SymmetricLayout());
	// makes sure that this is a symmetric matrix

	const size_t n = m_block_cols_list.size();
	// get number of block columns

	r_elim_tree.resize(n);
	std::vector<size_t> &r_highest_dependence = r_workspace;
	r_workspace.resize(n);
	// allocate result and workspace

	size_t n_column_id = 0;
	for(_TyColumnConstIter p_col_it = m_block_cols_list.begin(),
	   p_col_end_it = m_block_cols_list.end(); p_col_it != p_col_end_it;
	   ++ p_col_it, ++ n_column_id) {
		const TColumn &r_t_column = *p_col_it;
		// for each column (index is n_column_id)

		r_elim_tree[n_column_id] = size_t(-1);
		r_highest_dependence[n_column_id] = size_t(-1);
		// no information for this column yet

		for(_TyBlockConstIter
		   p_block_it = r_t_column.block_list.begin(),
		   p_block_end_it = r_t_column.block_list.end();
		   p_block_it != p_block_end_it; ++ p_block_it) {
			size_t n_ref_column_id = (*p_block_it).first;
			// get row of the block (a column which depends on this column)

			if(n_ref_column_id < n_column_id) { // uses the upper-diagonal part only
				do {
					size_t n_next_column_id = r_highest_dependence[n_ref_column_id];
					r_highest_dependence[n_ref_column_id] = n_column_id;
					if(n_next_column_id == size_t(-1)) { // can't go any higher, have the dependence
						r_elim_tree[n_ref_column_id] = n_column_id;
						break;
					}
					n_ref_column_id = n_next_column_id;
				} while(n_ref_column_id < n_column_id);
			} else
				break; // no sense iterating over the blocks under the diagonal
		}
	}
	// takes O(nnz(triu(this)) * average etree depth)
}

size_t CUberBlockMatrix::n_Build_EReach(size_t n_column_id,
	const std::vector<size_t> &r_elim_tree, std::vector<size_t> &r_ereach_stack,
	std::vector<size_t> &r_workspace) const
{
	//CheckIntegrity(); // inside b_SymmetricLayout()
	_ASSERTE(b_SymmetricLayout()); // could build an etree in the first place?
	_ASSERTE(r_elim_tree.size() == m_block_cols_list.size());

	const size_t n = m_block_cols_list.size();
	// get number of block columns

	_ASSERTE(n_column_id < n);
	// make sure the column id is valid

	_ASSERTE(r_ereach_stack.size() >= n);
	_ASSERTE(r_workspace.size() >= n);
	// need space for n elemets

#ifdef _DEBUG
	for(size_t i = 0; i < n; ++ i)
		_ASSERTE(!r_workspace[i]);
#endif //_DEBUG
	// make sure that the workspace contains all nulls (or at least n of them)

#if 0
	// this is good and even works for some matrices, but it fails if one
	// column depends on more columns (the paths, although themselves in good
	// order, are concatenated in reverse order)

	size_t n_ereach_size = 0;
	r_workspace[n_column_id] = true;
	const TColumn &r_t_column = m_block_cols_list[n_column_id];
	for(_TyBlockConstIter
	   p_block_it = r_t_column.block_list.begin(),
	   p_block_end_it = r_t_column.block_list.end();
	   p_block_it != p_block_end_it; ++ p_block_it) {
		size_t n_ref_column_id = (*p_block_it).first;
		if(n_ref_column_id > n_column_id)
			break;
		// go through the upper-triangular part

		for(; !r_workspace[n_ref_column_id]; n_ref_column_id =
		   r_elim_tree[n_ref_column_id], ++ n_ereach_size) {
			_ASSERTE(n_ereach_size < n); // if there are no loops, the path will never be longer
			r_ereach_stack[n_ereach_size] = n_ref_column_id;
			r_workspace[n_ref_column_id] = true;
		}
		// walk the elimination tree and put the path onto the stack
	}
	// O(n_ereach_size)

	for(size_t i = 0; i < n_ereach_size; ++ i)
		r_workspace[r_ereach_stack[i]] = false;
	r_workspace[n_column_id] = false;
	// clear the workspace for the next call to this function

	return n_ereach_size; // could return an iterator of r_ereach_stack instead
#else // 0
	size_t n_ereach_first = n;
	r_workspace[n_column_id] = true;
	const TColumn &r_t_column = m_block_cols_list[n_column_id];
	for(_TyBlockConstIter
	   p_block_it = r_t_column.block_list.begin(),
	   p_block_end_it = r_t_column.block_list.end();
	   p_block_it != p_block_end_it; ++ p_block_it) {
		size_t n_ref_column_id = (*p_block_it).first;
		if(n_ref_column_id > n_column_id)
			break;
		// go through the upper-triangular part

		size_t n_path_length = 0;
		for(; !r_workspace[n_ref_column_id]; n_ref_column_id =
		   r_elim_tree[n_ref_column_id], ++ n_path_length) {
			_ASSERTE(n_path_length < n); // if there are no loops, the path will never be longer
			r_ereach_stack[n_path_length] = n_ref_column_id;
			r_workspace[n_ref_column_id] = true;
		}
		// walk the elimination tree and put the path onto the stack

		while(n_path_length)
			r_ereach_stack[-- n_ereach_first] = r_ereach_stack[-- n_path_length];
		// copy to the end of the stack
	}
	// O(n_ereach_size)

	for(size_t i = n_ereach_first; i < n; ++ i)
		r_workspace[r_ereach_stack[i]] = false;
	r_workspace[n_column_id] = false;
	// clear the workspace for the next call to this function

	return n_ereach_first; // could return an iterator of r_ereach_stack instead
#endif // 0
}

bool CUberBlockMatrix::CholeskyOf(const CUberBlockMatrix &r_lambda,
	const std::vector<size_t> &r_elim_tree, std::vector<size_t> &r_workspace,
	std::vector<size_t> &r_zero_workspace) // throw(std::bad_alloc)
{
	Clear();

	// todo - make versions with workspace allocated outside, with precalculated etree (to reuse)
	// and with inplace processing

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

	for(size_t j = 0; j < n; ++ j) { // for every column (no sparsity here, L should be full-rank)
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
		size_t n_ereach_size = r_lambda.n_Build_EReach(j, r_elim_tree, ereach_stack, bitfield);
#else
		size_t n_ereach_first = r_lambda.n_Build_EReach(j, r_elim_tree, ereach_stack, bitfield);
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
			_TyMatrixXdRef L_block_kj(p_k_block_data, n_col_k_width, n_col_j_width);
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
					std::lower_bound(r_col_L_j.block_list.begin(),
					r_col_L_j.block_list.end(), k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(r_col_L_j.block_list.begin(),
					r_col_L_j.block_list.end(), k, CompareBlockRow);
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
						std::lower_bound(r_col_A_j.block_list.begin(),
						p_A_block_end_it, k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(r_col_A_j.block_list.begin(),
						p_A_block_end_it, k, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					// a rare case where ereach is not monotonically increasing
				}
				_ASSERTE(p_A_block_it != p_A_block_end_it);

				if((*p_A_block_it).first == k) {
					L_block_kj = _TyMatrixXdRef((*p_A_block_it).second,
						n_col_k_width, n_col_j_width);
					// here or after the following loop, no matter; L(k, j) not accessed : only L(k - 1, j) is accessed

					++ p_A_block_it;
					// skip to the next one for the next iteration / for the diagonal
				} else
					L_block_kj.setZero(); // sometimes A(k, j) is null
				// todo - move this to FBS as well, both copy and setZero could be faster
			}
			// look up the block in the source matrix

			_ASSERTE(!r_col_L_j.block_list.empty() && (!b_ereach_mono ||
				r_col_L_j.block_list.back().first == k)); // this is generally not true; it might not be the last
			_ASSERTE((*p_jk_block_it).first == k); // this should be always true
			// column j now contains the k-th block

			{
				_TyBlockConstIter
					p_j_block_it = r_col_L_j.block_list.begin(),
					p_j_block_end_it = p_jk_block_it/*r_col_L_j.block_list.end() - 1*/, // this is not end() - 1, might have to look for the block before k
					p_k_block_it = r_col_L_k.block_list.begin();
#ifdef _DEBUG
				_TyBlockConstIter
					p_k_block_end_it = r_col_L_k.block_list.end(); // only used in assertions
#endif // _DEBUG
				// have to loop through both lists and merge the blocks to find the ones, referencing the same rows

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					p_k_block_it != p_k_block_end_it);
				// if the first is a non-empty range, the second should be neither

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					r_col_L_k.block_list.back().first == k);
				// the last block of k-th column is on diagonal (kth row) and it is always present
				// this block will serve as the sentinell for this loop

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					(*(p_j_block_end_it - 1)).first < k); // note that the p_j_block_end_it points at the diagonal block
				// the last block of j-th column (omitting the diagonal) is *at most* at kth row

				for(; p_j_block_it != p_j_block_end_it; ++ p_j_block_it) { // for each block in column j
					size_t n_row_i = (*p_j_block_it).first;
					_ASSERTE(n_row_i < k); // make sure that the sentinell is indeed functional
					_ASSERTE(p_k_block_it != p_k_block_end_it); // should not be pointing at the end in the first place (if we got inside this loop)
					while((*p_k_block_it).first < n_row_i) {
						++ p_k_block_it;
						_ASSERTE(p_k_block_it != p_k_block_end_it); // should never reach the end (sentinell)
					}
					if((*p_k_block_it).first == n_row_i) {
						const size_t n_row_i_height = m_block_cols_list[n_row_i].n_width;
						// an optimistic case, we found blocks at the same row

						_TyMatrixXdRef L_block_ik((*p_k_block_it).second,
							n_row_i_height, n_col_k_width);
						_TyMatrixXdRef L_block_ij((*p_j_block_it).second,
							n_row_i_height, n_col_j_width);

#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						L_block_kj -= L_block_ik.transpose().lazyProduct(L_block_ij);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						L_block_kj.noalias() -= L_block_ik.transpose() * L_block_ij;
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						 // takes blocks from two different columns - need to know which columns have nonzero blocks
					} else {
						// next block in column k is on the next row,
						// we have to skip to the next block in jth row
					}
				}
				// this loop can be probably written in many different ways
				// @todo - investigate which gives the best run time
				// @todo - try to reverse it and loop for p_k_block_it, looking for p_j_block_it (could it save one more iterator?)
			}
			// cmod; causes fill-in in the current column

			_ASSERTE(!m_block_cols_list[k].block_list.empty() &&
				m_block_cols_list[k].block_list.back().first == k); // makes sure that k-th column contains a diagonal block
			_TyMatrixXdRef d(m_block_cols_list[k].block_list.back().second,
				n_col_k_width, n_col_k_width);
			d.triangularView<Eigen::Upper>().transpose().solveInPlace(L_block_kj); // modifies L_block_kj
			// d.marked<Eigen::UpperTriangular>().transpose().solveTriangularInPlace(L_block_kj); // the above line is deprecated, this line should do the trick in the new version(L_block_kj) of Eigen
			// late divide by the diagonal entries (probably unavoidable)
		}
		// complex data dependencies in here, not sure if i like it

		{
			_ASSERTE(r_col_L_j.block_list.empty() || r_col_L_j.block_list.back().first < j);
			// makes sure that the column doesn't contain the diagonal block yet

			double *p_diag_data = p_Get_DenseStorage(n_col_j_width * n_col_j_width);
			r_col_L_j.block_list.push_back(TColumn::TBlockEntry(j, p_diag_data));
			_TyMatrixXdRef L_block_jj(p_diag_data, n_col_j_width, n_col_j_width);
			// allocates a new block in this column

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

			L_block_jj = _TyMatrixXdRef((*p_A_block_it).second, n_col_j_width, n_col_j_width);
			// copy contents of the diagonal block

			for(_TyBlockConstIter p_L_block_j_it =
			   r_col_L_j.block_list.begin(), p_j_block_end_it = r_col_L_j.block_list.end() - 1;
			   p_L_block_j_it != p_j_block_end_it; ++ p_L_block_j_it) {
				_TyMatrixXdRef L_block_ij((*p_L_block_j_it).second,
					m_block_cols_list[(*p_L_block_j_it).first].n_width, n_col_j_width);
				L_block_jj.triangularView<Eigen::Upper>() -=
					L_block_ij.transpose().lazyProduct(L_block_ij);
				// this is an extremely interesting idea, saving 1/2 of diagonal FLOPs (possibly not with SSE)
			}
			// use sparse loop instead
			// calculates square of all the blocks in the current column (has shape
			// column-width * column-width, and the diagonal block has the same shape)

			Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> chol(L_block_jj); // Eigen::LLT only accesses a half of the matrix (upper tri in this case), no need to clear the lower half
			if(chol.info() != Eigen::Success) {
				//printf("error: not pos def\n"); // not pos def
				Clear(); // otherwise leaving uninit columns behind, CheckIntegrity() will yell
				return false;
			}
			L_block_jj = chol.matrixU();
			// calculates cholesky of a square block
		}
		// cdiv; reads an entire column and produces diagonal

		_ASSERTE(r_col_L_j.block_list.size() == n - n_ereach_first/*n_ereach_size*/ + 1);
		// make sure we preallocated it correclty
	}
	// note that Pigglet could probably write it better

	//CheckIntegrity(true);
	// makes sure the factor matrix (this) is ok

	return true;
}

bool CUberBlockMatrix::CholeskyOf(const CUberBlockMatrix &r_lambda,
	const std::vector<size_t> &r_elim_tree, std::vector<size_t> &r_workspace,
	std::vector<size_t> &r_zero_workspace, size_t n_start_on_column) // throw(std::bad_alloc)
{
	//Clear();

	// todo - make versions with workspace allocated outside, with precalculated etree (to reuse)
	// and with inplace processing

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

	for(size_t j = n_start_on_column; j < n; ++ j) { // for every column (no sparsity here, L should be full-rank)
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
		size_t n_ereach_size = r_lambda.n_Build_EReach(j, r_elim_tree, ereach_stack, bitfield);
#else
		size_t n_ereach_first = r_lambda.n_Build_EReach(j, r_elim_tree, ereach_stack, bitfield);
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
			_TyMatrixXdRef L_block_kj(p_k_block_data, n_col_k_width, n_col_j_width);
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
					std::lower_bound(r_col_L_j.block_list.begin(),
					r_col_L_j.block_list.end(), k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					std::lower_bound(r_col_L_j.block_list.begin(),
					r_col_L_j.block_list.end(), k, CompareBlockRow);
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
						std::lower_bound(r_col_A_j.block_list.begin(),
						p_A_block_end_it, k, CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(r_col_A_j.block_list.begin(),
						p_A_block_end_it, k, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					// a rare case where ereach is not monotonically increasing
				}
				_ASSERTE(p_A_block_it != p_A_block_end_it);

				if((*p_A_block_it).first == k) {
					L_block_kj = _TyMatrixXdRef((*p_A_block_it).second,
						n_col_k_width, n_col_j_width);
					// here or after the following loop, no matter; L(k, j) not accessed : only L(k - 1, j) is accessed

					++ p_A_block_it;
					// skip to the next one for the next iteration / for the diagonal
				} else
					L_block_kj.setZero(); // sometimes A(k, j) is null
			}
			// look up the block in the source matrix

			_ASSERTE(!r_col_L_j.block_list.empty() && (!b_ereach_mono ||
				r_col_L_j.block_list.back().first == k)); // this is generally not true; it might not be the last
			_ASSERTE((*p_jk_block_it).first == k); // this should be always true
			// column j now contains the k-th block

			{
				_TyBlockConstIter
					p_j_block_it = r_col_L_j.block_list.begin(),
					p_j_block_end_it = p_jk_block_it/*r_col_L_j.block_list.end() - 1*/, // this is not end() - 1, might have to look for the block before k
					p_k_block_it = r_col_L_k.block_list.begin();
#ifdef _DEBUG
				_TyBlockConstIter
					p_k_block_end_it = r_col_L_k.block_list.end(); // only used in assertions
#endif // _DEBUG
				// have to loop through both lists and merge the blocks to find the ones, referencing the same rows

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					p_k_block_it != p_k_block_end_it);
				// if the first is a non-empty range, the second should be neither

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					r_col_L_k.block_list.back().first == k);
				// the last block of k-th column is on diagonal (kth row) and it is always present
				// this block will serve as the sentinell for this loop

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					(*(p_j_block_end_it - 1)).first < k); // note that the p_j_block_end_it points at the diagonal block
				// the last block of j-th column (omitting the diagonal) is *at most* at kth row

				for(; p_j_block_it != p_j_block_end_it; ++ p_j_block_it) { // for each block in column j
					size_t n_row_i = (*p_j_block_it).first;
					_ASSERTE(n_row_i < k); // make sure that the sentinell is indeed functional
					_ASSERTE(p_k_block_it != p_k_block_end_it); // should not be pointing at the end in the first place (if we got inside this loop)
					while((*p_k_block_it).first < n_row_i) {
						++ p_k_block_it;
						_ASSERTE(p_k_block_it != p_k_block_end_it); // should never reach the end (sentinell)
					}
					if((*p_k_block_it).first == n_row_i) {
						const size_t n_row_i_height = m_block_cols_list[n_row_i].n_width;
						// an optimistic case, we found blocks at the same row

						_TyMatrixXdRef L_block_ik((*p_k_block_it).second,
							n_row_i_height, n_col_k_width);
						_TyMatrixXdRef L_block_ij((*p_j_block_it).second,
							n_row_i_height, n_col_j_width);

#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						L_block_kj -= L_block_ik.transpose().lazyProduct(L_block_ij);
#else // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						L_block_kj.noalias() -= L_block_ik.transpose() * L_block_ij;
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
						 // takes blocks from two different columns - need to know which columns have nonzero blocks
					} else {
						// next block in column k is on the next row,
						// we have to skip to the next block in jth row
					}
				}
				// this loop can be probably written in many different ways
				// @todo - investigate which gives the best run time
				// @todo - try to reverse it and loop for p_k_block_it, looking for p_j_block_it (could it save one more iterator?)
			}
			// cmod; causes fill-in in the current column

			_ASSERTE(!m_block_cols_list[k].block_list.empty() &&
				m_block_cols_list[k].block_list.back().first == k); // makes sure that k-th column contains a diagonal block
			_TyMatrixXdRef d(m_block_cols_list[k].block_list.back().second,
				n_col_k_width, n_col_k_width);
			d.triangularView<Eigen::Upper>().transpose().solveInPlace(L_block_kj); // modifies L_block_kj
			// d.marked<Eigen::UpperTriangular>().transpose().solveTriangularInPlace(L_block_kj); // the above line is deprecated, this line should do the trick in the new version(L_block_kj) of Eigen
			// late divide by the diagonal entries (probably unavoidable)
		}
		// complex data dependencies in here, not sure if i like it

		{
			_ASSERTE(r_col_L_j.block_list.empty() || r_col_L_j.block_list.back().first < j);
			// makes sure that the column doesn't contain the diagonal block yet

			double *p_diag_data = p_Get_DenseStorage(n_col_j_width * n_col_j_width);
			r_col_L_j.block_list.push_back(TColumn::TBlockEntry(j, p_diag_data));
			_TyMatrixXdRef L_block_jj(p_diag_data, n_col_j_width, n_col_j_width);
			// allocates a new block in this column

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
			// last block and we can optimize for that) // RSStodo

			L_block_jj = _TyMatrixXdRef((*p_A_block_it).second, n_col_j_width, n_col_j_width);
			// copy contents of the diagonal block

			for(_TyBlockConstIter p_L_block_j_it =
			   r_col_L_j.block_list.begin(), p_j_block_end_it = r_col_L_j.block_list.end() - 1;
			   p_L_block_j_it != p_j_block_end_it; ++ p_L_block_j_it) {
				_TyMatrixXdRef L_block_ij((*p_L_block_j_it).second,
					m_block_cols_list[(*p_L_block_j_it).first].n_width, n_col_j_width);
				L_block_jj.triangularView<Eigen::Upper>() -=
					L_block_ij.transpose().lazyProduct(L_block_ij);
				// this is an extremely interesting idea, saving 1/2 of diagonal FLOPs (possibly not with SSE)
			}
			// use sparse loop instead // todo - can accum this during the cmod loop, the results will be in cache // RSStodo
			// calculates square of all the blocks in the current column (has shape
			// column-width * column-width, and the diagonal block has the same shape)

			Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> chol(L_block_jj); // Eigen::LLT only accesses a half of the matrix (upper tri in this case), no need to clear the lower half
			if(chol.info() != Eigen::Success) {
				//printf("error: not pos def\n"); // not pos def
				Clear(); // otherwise leaving uninit columns behind, CheckIntegrity() will yell
				return false;
			}
			L_block_jj = chol.matrixU();
			// calculates cholesky of a square block
		}
		// cdiv; reads an entire column and produces diagonal

		_ASSERTE(r_col_L_j.block_list.size() == n - n_ereach_first/*n_ereach_size*/ + 1);
		// make sure we preallocated it correclty
	}
	// note that Pigglet could probably write it better

	//CheckIntegrity(true);
	// makes sure the factor matrix (this) is ok

	return true;
}

void CUberBlockMatrix::Get_UpperTriangular_BlockFrontline(std::vector<size_t> &r_frontline) const // throw(std::bad_alloc)
{
	CheckIntegrity();

	size_t n = m_block_cols_list.size();
	if(r_frontline.capacity() < n) {
		r_frontline.clear();
		r_frontline.reserve(std::max(n, 2 * r_frontline.capacity()));
	}
	r_frontline.resize(n);

	for(size_t i = 0; i < n; ++ i) {
		const TColumn &r_t_col = m_block_cols_list[i];
		if(r_t_col.block_list.empty())
			r_frontline[i] = i; // no block - assume the frontline touches diagonal here
		else
			r_frontline[i] = r_t_col.block_list.front().first; // use the smallest row id as fronline
		// get frontline
	}
}

size_t CUberBlockMatrix::n_Get_UpperTriangular_BlockFrontline_Minimum(size_t n_order_lo, size_t n_order_hi) const
{
	CheckIntegrity();

	_ASSERTE(n_order_hi >= n_order_lo);
	_ASSERTE(n_order_hi <= m_block_cols_list.size());

	size_t n_min = n_order_hi;
	for(++ n_order_lo; n_order_lo < n_order_hi; ++ n_order_lo) {
		const TColumn &r_t_col = m_block_cols_list[n_order_lo];
		size_t n_frontline;
		if(r_t_col.block_list.empty())
			n_frontline = n_order_lo; // no block - assume the frontline touches diagonal here
		else
			n_frontline = r_t_col.block_list.front().first; // use the smallest row id as fronline
		// get frontline

		if(n_min > n_frontline)
			n_min = n_frontline;
	}

	return n_min;
}

void CUberBlockMatrix::InverseOf_Symmteric(const CUberBlockMatrix &r_A) // throw(std::bad_alloc)
{
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

	// inverse of a positive definite matrix = L^-TL^-1, which probably does not help too much

	/*std::vector<size_t> supernode_memberships(n, size_t(-1));
	std::vector<std::pair<size_t, size_t> > supernodes;
	// each supernode has the first column and last column (inclusive)

	for(size_t i = 0; i < n; ++ i) {
		size_t n_start = i;
		size_t n_end = i;
		// start and end of the supernode

		size_t n_col_id = i;
		// to make expressions better readable

		size_t n_node_member = supernode_memberships[n_col_id];
		// which node member are we?

		const TColumn &r_t_col = r_A.m_block_cols_list[i];
		for(size_t j = 0, m = r_t_col.block_list.size(); j < m; ++ j) {
			size_t n_row_id = r_t_col.block_list[j].first;
			// we have a nnz block at (n_row_id, n_col_id)

			if(n_row_id == n_col_id)
				continue;
			// diagonal blocks do not add fill-in

			if(n_row_id < n_col_id) { // above-diagonal block
				size_t n_supernode = supernode_memberships[n_row_id];
				supernodes[n_supernode].second = std::max(supernodes[n_supernode].second, n_col_id); // add itself
				n_start = std::min(n_start, supernodes[n_supernode].first);
				n_end = std::max(n_end, supernodes[n_supernode].second);
				n_node_member = std::min(n_node_member, n_supernode); // lets become a member
				// add this column to an existing supernode
			} else {
				_ASSERTE(n_row_id > n_col_id); // below-diagonal block
				// have to introduce a new supernode, or extend an existing one

				if(n_node_member == size_t(-1)) {
					n_node_member = supernodes.size();
					supernodes.push_back(std::make_pair(n_start, n_end));
				}
				// in case we are not a member of any supernode yet, make a new one

				if(supernode_memberships[n_row_id] != size_t(-1)) {
					// have two supernodes which overlap
				}

				std::pair<size_t, size_t> &r_my_supernode = supernodes[n_node_member];
				r_my_supernode.second = std::max(r_my_supernode.second, n_row_id); // add ref
				n_end = std::max(r_my_supernode.second, n_row_id); // also
				supernode_memberships[n_row_id] = n_node_member; // add ref to the same supernode
			}
		}

		supernode_memberships[n_col_id] = n_node_member;
	}*/
	// creates supernode spans, which potentially overlap
	// overly complicated, need too much memory

	/*for(size_t i = 0; i < n; ++ i) {
		size_t n_start = i;
		size_t n_end = i;
		// start and end of the supernode

		_ASSERTE(r_t_col.block_list.empty() || r_t_col.block_list.front().first >= n_start);
		// make sure that this is the first block of the supernode (not looking back)

		const TColumn &r_t_col = r_A.m_block_cols_list[i];
		if(!r_t_col.block_list.empty())
			n_end = std::max(n_end, r_t_col.block_list.back().first);
		for(;;) {
			for(; i < n_end;) {
				++ i;
				const TColumn &r_t_col = r_A.m_block_cols_list[i];
				if(!r_t_col.block_list.empty()) {
					_ASSERTE(r_t_col.block_list.front().first >= n_start); // make sure it is not looking to the previous supernode
					n_end = std::max(n_end, r_t_col.block_list.back().first);
				}
				// the maximum referenced column
			}
			_ASSERTE(i == n_end);
			// find the maximum referenced column, an all columns between i and n_end

			bool b_expand = false;
			size_t n_new_end = n_end;
			for(int j = i + 1; j < n; ++ j) {
				const TColumn &r_t_col = r_A.m_block_cols_list[i];
				if(!r_t_col.block_list.empty()) {
					n_new_end = std::max(n_new_end, r_t_col.block_list.back().first); // track the new end in case we need to extend
					if(r_t_col.block_list.front().first <= n_end) { // some block in the right references our supernode, need to expand it
						b_expand = true;
						break;
					}
				}
			}
			// now need to loop over the remaining columns and see if there is a ref before begin

			if(!b_expand)
				break;
			// noone references us, the supernode is between n_start and n_end (inclusive)

			n_end = n_new_end;
			// supernode just extended, need to recheck for maximum
			// reference of all columns between i and n_end
		}
		_ASSERTE(i == n_end);
		// does O(n - i) iterations (always, no matter how many times the outer loop repeats)

		// now we know that the supernode is n_start to n_end, we need to grab the blocks from A,
		// put them in a dense matrix, invert it and put it back to this
	}*/
	// does up to O(n) iterations (depends on the size of the supernodes, as it skips all the way to the end)
	// the worst-case complexity of finding the supernodes is therefore O(1/2 * n^2), the best case is O(n).
	// it goes faster on very connected matrices, slower on block diagonal ones. it could probably be done faster ...

	// if we did forward iterations, it would be O(n) in a nested loop (outer loop loops over starting nodes, inner loop finds the end)
	// then we need backward iterations, it would be another O(n) (same as above, only backwards)
	// then we would need to merge the (sorted) intervals, which is probably another worst-case O(n), best case O(1)
	// thus, it could be done in O(2n + 1) to O(3n), but we would need O(2n) storage

	// for w100k, it is 5000000000 / 300000 = 16666.667x faster, that is quite a lot for (worst-case) 800 kB

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
			_TyMatrixXdRef src_block((double*)p_data, n_row_height, r_t_col.n_width);
			// make a map of the source block

			TColumn &r_t_dest_col = m_block_cols_list[n_begin];
			_ASSERTE(r_t_dest_col.block_list.empty()); // should be initially empty
			r_t_dest_col.block_list.reserve(1);
			TColumn::TBlockEntry t_block(n_row_id,
				p_Get_DenseStorage(n_row_height * r_t_col.n_width));
			r_t_dest_col.block_list.push_back(t_block);
			// alloc a new (destination) block in this matrix

			_TyMatrixXdRef dest_block(t_block.second, n_row_height, r_t_col.n_width);
			dest_block = src_block.inverse();
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

						_TyMatrixXdRef fill_block((double*)p_data, n_row_height, r_t_col.n_width);
						// make a map

						dense.block(n_row_org - n_first_col, n_dest_column_org,
							n_row_height, r_t_col.n_width) = fill_block;
						// copy the data inside the dense matrix
					}
				}
				// go through all the columns in A and put them in the dense matrix

				dense = dense.inverse();
				// invert the matrix, making it most likely rather dense
			}
			// todo - it would be better to slice the matrix and make a blockwise sparse LU decomposition,
			// and use that to calculate the inverse (now we are calculating a rather expensive inverse
			// of a dense (although containing many zeroes) matrix, using dense LU decomposition inside Eigen)

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

	// again, orderings to the rescue - ordering the matrix to have the connected
	// blocks contiguous will make the inverse more sparse
}

#ifdef __UBER_BLOCK_MATRIX_IO

static void TrimSpace(std::string &r_s_string)
{
	size_t b = 0, e = r_s_string.length();
	while(e > 0 && isspace((unsigned char)(r_s_string[e - 1])))
		-- e;
	while(b < e && isspace((unsigned char)(r_s_string[b])))
		++ b;
	r_s_string.erase(e);
	r_s_string.erase(0, b);
}

static bool ReadLine(std::string &r_s_line, FILE *p_fr)
{
	r_s_line.erase();
	for(int c = fgetc(p_fr); c != '\n' && c != EOF; c = fgetc(p_fr)) {
		if(ferror(p_fr))
			return false;
		try {
			r_s_line += c; // some string implementations don't have push_back()
		} catch(std::bad_alloc&) {
			return false;
		}
	}
	// read line

	return !ferror(p_fr);
}

// the two above functions are already in CParserBase, but as the block matrix code is otherwise
// independent from the SLAM code, this sin against good style is overlooked here

// note that the matrix IO functions could also easily be implemented as an outside object, which
// would make all this much easier, and CParserBase could depend on that // todo

bool CUberBlockMatrix::Load_BlockLayout(const char *p_s_layout_filename) // throw(std::bad_alloc)
{
	_ASSERTE(p_s_layout_filename);

	size_t n_layout_w, n_layout_h, n_layout_nnz;
	size_t n_layout_bw, n_layout_bh, n_layout_bnnz;
	std::vector<size_t> cols_cumsums;
	std::vector<size_t> rows_cumsums;

	FILE *p_fr = 0;
	do {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
		if(fopen_s(&p_fr, p_s_layout_filename, "r"))
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		if(!(p_fr = fopen(p_s_layout_filename, "r")))
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			return false;
		// open the layout file

#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
		if(fscanf_s(p_fr, PRIsize " x " PRIsize " (" PRIsize ")\n"
		   PRIsize " x " PRIsize " (" PRIsize ")\n", &n_layout_h,
		   &n_layout_w, &n_layout_nnz, &n_layout_bh,
		   &n_layout_bw, &n_layout_bnnz) != 6)
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		if(fscanf(p_fr, PRIsize " x " PRIsize " (" PRIsize ")\n"
		   PRIsize " x " PRIsize " (" PRIsize ")\n", &n_layout_h,
		   &n_layout_w, &n_layout_nnz, &n_layout_bh,
		   &n_layout_bw, &n_layout_bnnz) != 6)
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			break;
		// read layout geometry

		try {
			rows_cumsums.resize(n_layout_bh + 1);
			cols_cumsums.resize(n_layout_bw + 1);
		} catch(std::bad_alloc &r_exc) {
			fclose(p_fr); // clsoe the file
			throw r_exc; // rethrow
		}
		// the first one is zero

		bool b_fail = false;
		for(size_t i = 0, n = n_layout_bh + 1; i < n; ++ i) {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			if(fscanf_s(p_fr, PRIsize, &rows_cumsums[i]) != 1) {
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(fscanf(p_fr, PRIsize, &rows_cumsums[i]) != 1) {
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				b_fail = true;
				break;
			}
		}
		if(b_fail)
			break;
		for(size_t i = 0, n = n_layout_bw + 1; i < n; ++ i) {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			if(fscanf_s(p_fr, PRIsize, &cols_cumsums[i]) != 1) {
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(fscanf(p_fr, PRIsize, &cols_cumsums[i]) != 1) {
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				b_fail = true;
				break;
			}
		}
		if(b_fail)
			break;
		// read row and column cumsums

		fclose(p_fr);
		p_fr = 0;
		// close the file

		if(rows_cumsums.front() != 0 || rows_cumsums.back() != n_layout_h ||
		   cols_cumsums.front() != 0 || cols_cumsums.back() != n_layout_w)
			break;
		// unable to use stored layout: corrupt cumsums

		CUberBlockMatrix bm_lay(rows_cumsums.begin() + 1, rows_cumsums.end(),
			cols_cumsums.begin() + 1, cols_cumsums.end());
		Swap(bm_lay);
		// build an empty matrix with a given layout and swap
		// note that the file is closed now, don't mind bad_alloc here

		return true;
		// success
	} while(0);
	if(p_fr)
		fclose(p_fr);

	return false;
}

bool CUberBlockMatrix::Load_MatrixMarket(const char *p_s_filename, size_t n_block_size,
	bool b_allow_underfilled_last_block) // throw(std::bad_alloc)
{
	_ASSERTE(p_s_filename);

	FILE *p_fr = 0;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
	if(fopen_s(&p_fr, p_s_filename, "r"))
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	if(!(p_fr = fopen(p_s_filename, "r")))
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		return false;
	// open the matrix market file

	bool b_symmetric = false;
	bool b_had_specifier = false;
	bool b_had_header = false;
	size_t n_rows, n_cols, n_nnz = -1;
	std::string s_line;
	while(!feof(p_fr)) {
		if(!ReadLine(s_line, p_fr)) {
			fclose(p_fr);
			return false;
		}
		// read a single line

		size_t n_pos;
		if(!b_had_specifier && (n_pos = s_line.find("%%MatrixMarket")) != std::string::npos) {
			s_line.erase(0, n_pos + strlen("%%MatrixMarket"));
			// get rid of header

			if(s_line.find("matrix") == std::string::npos ||
			   s_line.find("coordinate") == std::string::npos ||
			   s_line.find("real") == std::string::npos) {
				fclose(p_fr);
				return false;
			}
			// must be matrix coordinate real

			if(s_line.find("general") != std::string::npos)
				b_symmetric = false;
			else if(s_line.find("symmetric") != std::string::npos)
				b_symmetric = true;
			else {
				b_symmetric = false;
				/*fclose(p_fr);
				return false;*/ // or assume general
			}
			// either general or symmetric

			b_had_specifier = true;
			continue;
		}
		if((n_pos = s_line.find('%')) != std::string::npos)
			s_line.erase(n_pos);
		TrimSpace(s_line);
		if(s_line.empty())
			continue;
		// trim comments, skip empty lines

		if(!b_had_header) {
			if(!b_had_specifier) {
				fclose(p_fr);
				return false;
			}
			// specifier must come before header

#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			if(sscanf_s(s_line.c_str(), PRIsize " " PRIsize " " PRIsize,
			   &n_rows, &n_cols, &n_nnz) != 3) {
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(sscanf(s_line.c_str(), PRIsize " " PRIsize " " PRIsize,
			   &n_rows, &n_cols, &n_nnz) != 3) {
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				fclose(p_fr);
				return false;
			}
			// read header

			if(n_nnz > n_rows * n_cols) {
				fclose(p_fr);
				return false;
			}
			// sanity check (may also fail on big matrices)

			b_had_header = true;
			break;
		}
	}
	if(ferror(p_fr) || !b_had_header) {
		fclose(p_fr);
		return false;
	}
	fclose(p_fr);
	// read matrix market header in order to know the size of the matrix

	if(!b_allow_underfilled_last_block && (n_rows % n_block_size || n_cols % n_block_size))
		return false;
	size_t n_block_row_num = (n_rows + n_block_size - 1) / n_block_size;
	size_t n_block_col_num = (n_cols + n_block_size - 1) / n_block_size;
	// make sure that the given matrix can fit the given layout

	m_block_rows_list.reserve(n_block_row_num);
	m_block_cols_list.reserve(n_block_col_num);
	// make sure it fits

	Clear();
	// destroy the matrix

	m_block_rows_list.resize(n_block_row_num);
	m_block_cols_list.resize(n_block_col_num);
	// now they should both pass

	for(size_t i = 0, n_cumsum = 0; i < n_block_row_num; ++ i, n_cumsum += n_block_size) {
		m_block_rows_list[i].n_cumulative_height_sum = n_cumsum;
		m_block_rows_list[i].n_height = n_block_size;
	}
	for(size_t i = 0, n_cumsum = 0; i < n_block_col_num; ++ i, n_cumsum += n_block_size) {
		m_block_cols_list[i].n_cumulative_width_sum = n_cumsum;
		m_block_cols_list[i].n_width = n_block_size;
	}
	// build the layout

	m_n_row_num = n_block_row_num * n_block_size;
	m_n_col_num = n_block_col_num * n_block_size;
	// set size

	CheckIntegrity();
	// make sure the matrix is ok

	return Load_MatrixMarket_Layout(p_s_filename);
}

bool CUberBlockMatrix::Load_MatrixMarket_Layout(const char *p_s_filename) // throw(std::bad_alloc)
{
	SetZero();
	// delete any errant blocks

	FILE *p_fr = 0;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
	if(fopen_s(&p_fr, p_s_filename, "r"))
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	if(!(p_fr = fopen(p_s_filename, "r")))
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		return false;
	// open the matrix market file

	bool b_symmetric = false;
	bool b_had_specifier = false;
	bool b_had_header = false;
	size_t n_rows, n_cols, n_nnz = -1, n_read_nnz = 0;
	std::string s_line;
	while(!feof(p_fr)) {
		if(!ReadLine(s_line, p_fr)) {
			fclose(p_fr);
			return false;
		}
		// read a single line

		size_t n_pos;
		if(!b_had_specifier && (n_pos = s_line.find("%%MatrixMarket")) != std::string::npos) {
			s_line.erase(0, n_pos + strlen("%%MatrixMarket"));
			// get rid of header

			if(s_line.find("matrix") == std::string::npos ||
			   s_line.find("coordinate") == std::string::npos ||
			   s_line.find("real") == std::string::npos) {
				fclose(p_fr);
				return false;
			}
			// must be matrix coordinate real

			if(s_line.find("general") != std::string::npos)
				b_symmetric = false;
			else if(s_line.find("symmetric") != std::string::npos)
				b_symmetric = true;
			else {
				b_symmetric = false;
				/*fclose(p_fr);
				return false;*/ // or assume general
			}
			// either general or symmetric

			b_had_specifier = true;
			continue;
		}
		if((n_pos = s_line.find('%')) != std::string::npos)
			s_line.erase(n_pos);
		TrimSpace(s_line);
		if(s_line.empty())
			continue;
		// trim comments, skip empty lines

		if(!b_had_header) {
			if(!b_had_specifier) {
				fclose(p_fr);
				return false;
			}
			// specifier must come before header

#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			if(sscanf_s(s_line.c_str(), PRIsize " " PRIsize " " PRIsize,
			   &n_rows, &n_cols, &n_nnz) != 3) {
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(sscanf(s_line.c_str(), PRIsize " " PRIsize " " PRIsize,
			   &n_rows, &n_cols, &n_nnz) != 3) {
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				fclose(p_fr);
				return false;
			}
			// read header

			if(m_n_row_num != n_rows || m_n_col_num != n_cols) {
				fclose(p_fr);
				return false;
			}
			// make sure that the layout matches the matrix

			if(n_nnz > n_rows * n_cols) {
				fclose(p_fr);
				return false;
			}
			// sanity check (may also fail on big matrices)

			if(n_nnz && (m_block_rows_list.empty() || m_block_cols_list.empty())) {
				fclose(p_fr);
				return false;
			}
			// there are nonzeros, but the matrix layout is empty

			b_had_header = true;
		} else {
			double f_value;
			size_t n_rowi, n_coli;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
			int n = sscanf_s(s_line.c_str(), PRIsize " " PRIsize " %lf", &n_rowi, &n_coli, &f_value);
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			int n = sscanf(s_line.c_str(), PRIsize " " PRIsize " %lf", &n_rowi, &n_coli, &f_value);
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			if(n < 2) {
				fclose(p_fr);
				return false;
			} else if(n == 2)
				f_value = 1; // a binary matrix
			// read a triplet

			-- n_rowi;
			-- n_coli;
			// indices are 1-based, i.e. a(1,1) is the first element

			if(n_rowi >= n_rows || n_coli >= n_cols || n_read_nnz >= n_nnz) {
				fclose(p_fr);
				return false;
			}
			// make sure it figures

			for(int n_pass = 0;; ++ n_pass) {
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
				_TyRowConstIter p_row_it = std::upper_bound(m_block_rows_list.begin(),
					m_block_rows_list.end(), n_rowi, CFindLowerRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				_TyRowConstIter p_row_it = std::upper_bound(m_block_rows_list.begin(),
					m_block_rows_list.end(), n_rowi, _FindLowerRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				_ASSERTE(p_row_it != m_block_rows_list.begin()); // should start with zero, n_rowi is nonnegative
				-- p_row_it;
				const TRow &r_t_row = *p_row_it;
				_ASSERTE(r_t_row.n_cumulative_height_sum <= n_rowi &&
					r_t_row.n_cumulative_height_sum + r_t_row.n_height > n_rowi); // make sure it's it
				//
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
				_TyColumnIter p_col_it = std::upper_bound(m_block_cols_list.begin(),
					m_block_cols_list.end(), n_coli, CFindLowerColumn()); // t_odo
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				_TyColumnIter p_col_it = std::upper_bound(m_block_cols_list.begin(),
					m_block_cols_list.end(), n_coli, _FindLowerColumn); // t_odo
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
				_ASSERTE(p_col_it != m_block_cols_list.begin()); // should start with zero, n_rowi is nonnegative
				-- p_col_it;
				TColumn &r_t_col = *p_col_it;
				_ASSERTE(r_t_col.n_cumulative_width_sum <= n_coli &&
					r_t_col.n_cumulative_width_sum + r_t_col.n_width > n_coli); // make sure it's it
				// find a column and a row (not precise it, but one containing the element in question),
				// note that there is no function for that so far (elementwise access was not needed)

				size_t n_row_id = p_row_it - m_block_rows_list.begin();
				// get row id

				size_t n_block_row_num = r_t_row.n_height;
				size_t n_block_col_num = r_t_col.n_width;
				// get size of the block

				size_t n_block_row_off = n_rowi - r_t_row.n_cumulative_height_sum;
				size_t n_block_col_off = n_coli - r_t_col.n_cumulative_width_sum;
				// where the nonzero belongs in the dense block

				try {
					bool b_uninitialized;
					double *p_data = p_AllocBlockData(n_row_id, r_t_col,
						n_block_row_num, n_block_col_num, b_uninitialized);
					_TyMatrixXdRef t_block(p_data, n_block_row_num, n_block_col_num);
					if(b_uninitialized)
						t_block.setZero(); // clean first, if a new block
					t_block(n_block_row_off, n_block_col_off) = f_value;
				} catch(std::bad_alloc &r_exc) {
					fclose(p_fr);
					throw r_exc; // rethrow
				}
				// find the block and assign the element value

				if(b_symmetric && !n_pass && n_rowi != n_coli)
					std::swap(n_rowi, n_coli);
				else
					break;
				// handle symmetric matrices
			}
			++ n_read_nnz;
			// append the nonzero
		}
	}
	// read the file line by line

	if(ferror(p_fr) || !b_had_header || n_read_nnz != n_nnz) {
		fclose(p_fr);
		return false;
	}
	// make sure no i/o errors occured, and that the entire matrix was read

	fclose(p_fr);
	// close the file

	CheckIntegrity(true);
	// make sure the matrix is ok

	return true;
}

bool CUberBlockMatrix::Save_BlockLayout(const char *p_s_layout_filename) const
{
	CheckIntegrity();

	_ASSERTE(p_s_layout_filename);

	size_t n_rows = n_Row_Num();
	size_t n_columns = n_Column_Num();
	size_t n_nonzeros = n_NonZero_Num();

	FILE *p_fw;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
	if(fopen_s(&p_fw, p_s_layout_filename, "w"))
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	if(!(p_fw = fopen(p_s_layout_filename, "w")))
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		return false;
	fprintf(p_fw, PRIsize " x " PRIsize " (" PRIsize ")\n", n_rows, n_columns, n_nonzeros);
	fprintf(p_fw, PRIsize " x " PRIsize " (" PRIsize ")\n", n_BlockRow_Num(),
		n_BlockColumn_Num(), n_Block_Num());
	for(size_t i = 0, n = n_BlockRow_Num(); i < n; ++ i)
		fprintf(p_fw, PRIsize " ", n_BlockRow_Base(i));
	fprintf(p_fw, PRIsize "\n", n_rows);
	for(size_t i = 0, n = n_BlockColumn_Num(); i < n; ++ i)
		fprintf(p_fw, PRIsize " ", n_BlockColumn_Base(i));
	fprintf(p_fw, PRIsize "\n", n_columns);
	bool b_result = !ferror(p_fw);
	fclose(p_fw);
	// save the block layout

	return b_result;
}

bool CUberBlockMatrix::Save_MatrixMarket(const char *p_s_filename,
	const char *p_s_layout_filename /*= 0*/,
	const char *p_s_kind /*= "general block matrix"*/) const
{
	CheckIntegrity(true);

	_ASSERTE(p_s_filename);

	if(p_s_layout_filename && !Save_BlockLayout(p_s_layout_filename))
		return false;
	// save the block layout

	FILE *p_fw;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
	if(fopen_s(&p_fw, p_s_filename, "w"))
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	if(!(p_fw = fopen(p_s_filename, "w")))
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		return false;
	fprintf(p_fw, "%%%%MatrixMarket matrix coordinate real general\n"); // t_odo - maybe here should also be some keywords
	fprintf(p_fw, "%%-------------------------------------------------------------------------------\n");
	fprintf(p_fw, "%% �berBlockMatrix matrix dump\n");
	if(p_s_kind)
		fprintf(p_fw, "%% kind: %s\n", p_s_kind);
	fprintf(p_fw, "%%-------------------------------------------------------------------------------\n");
	size_t n_rows = n_Row_Num();
	size_t n_columns = n_Column_Num();
	size_t n_nonzeros = n_NonZero_Num();
	fprintf(p_fw, "" PRIsize " " PRIsize " " PRIsize "\n", n_rows, n_columns, n_nonzeros);
	for(size_t n_col = 0, n_column_blocks = n_BlockColumn_Num();
	   n_col < n_column_blocks; ++ n_col) {
		size_t n_col_base = n_BlockColumn_Base(n_col);
		size_t n_col_width = n_BlockColumn_Column_Num(n_col);
		size_t n_block_num = n_BlockColumn_Block_Num(n_col);
		for(size_t j = 0; j < n_block_num; ++ j) {
			size_t n_row = n_Block_Row(n_col, j);
			size_t n_row_base = n_BlockRow_Base(n_row);
			size_t n_row_height = n_BlockRow_Row_Num(n_row);

			CUberBlockMatrix::_TyMatrixXdRef t_block = t_Block_AtColumn(n_col, j);
			_ASSERTE(t_block.rows() == n_row_height && t_block.cols() == n_col_width);
			// get a block

			for(size_t k = 0; k < n_row_height; ++ k) {
				for(size_t l = 0; l < n_col_width; ++ l) {
					fprintf(p_fw, (fabs(t_block(k, l)) > 1)?
						PRIsize " " PRIsize " %.15f\n" :
						PRIsize " " PRIsize " %.15g\n", n_row_base + 1 + k,
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

#endif // __UBER_BLOCK_MATRIX_IO

/*
 *								=== ~CUberBlockMatrix ===
 */
