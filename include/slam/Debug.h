/*
								+-----------------------------------+
								|                                   |
								| ***  Debugging functionality  *** |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2012  |
								|                                   |
								|              Debug.h              |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __DEBUGGING_FUNCTIONALITY_INCLUDED
#define __DEBUGGING_FUNCTIONALITY_INCLUDED

/**
 *	@file include/slam/Debug.h
 *	@brief basic debugging functionality
 *	@author -tHE SWINe-
 *	@date 2012-04-26
 *
 *	@date 2012-08-22
 *
 *	Fixed CSparseMatrixShapedIterator and CSparseMatrixIterator, there was a slight misconception
 *	of CSparse format, the i array is mostly in sorted order, but may not be. This caused some of
 *	matrix elements to be skipped from iteration.
 *
 */

#include <csparse/cs.hpp>
#include <stdio.h>

#ifndef _ASSERTE
#include <assert.h>
/**
 *	@brief basic debug assertion macro
 *	@param[in] x is the condition that is supposed to hold
 */
#define _ASSERTE(x) assert(x)
#endif // !_ASSERTE

#include "slam/Bitmap.h"
#include "slam/Tga.h"

// @todo - at this point the functions below are only compiled in Main.cpp, but should really be moved to Debug.cpp to avoid multiple compilation and link conflicts

/**
 *	@brief implements some basic debugging functionality (printing stuff to stdout, basically)
 */
class CDebug {
public:
	/**
	 *	@brief shaped sparse matrix iterator
	 */
	class CSparseMatrixShapedIterator {
	protected:
		const cs *m_p_matrix; /**< @brief the matrix in question */
		size_t m_n_rows; /**< @brief number of rows of shape to iterate through */
		size_t m_n_columns; /**< @brief number of columns of shape to iterate through */
		size_t m_n_column; /**< @brief current column */
		size_t m_n_row; /**< @brief current row */
		size_t m_n_row_pointer; /**< @brief current row pointer offset (index in m_p_matrix->p) */
		size_t m_n_row_pointer_column_end; /**< @brief maximum value m_n_row_pointer can take at the current column */
		size_t m_n_row_pointer_end; /**< @brief maximum value m_n_row_pointer can take */

	public:
		/**
		 *	@brief constructor
		 *	@param[in] p_matrix is a sparse matrix to iterate over
		 */
		CSparseMatrixShapedIterator(const cs *p_matrix)
			:m_p_matrix(p_matrix), m_n_rows(p_matrix->m), m_n_columns(p_matrix->n),
			m_n_column(0), m_n_row(0), m_n_row_pointer(p_matrix->p[0]),
			m_n_row_pointer_column_end(p_matrix->p[1]),
			m_n_row_pointer_end(p_matrix->p[p_matrix->n])
		{
			if(m_n_row_pointer <= m_n_row_pointer_column_end) {
				size_t n_data_row = m_p_matrix->i[m_n_row_pointer];
				if(m_n_row_pointer == m_n_row_pointer_column_end || n_data_row != m_n_row) {
					size_t p = (m_n_column < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column] : m_n_row_pointer_end,
						e = (m_n_column + 1 < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column + 1] : m_n_row_pointer_end;
					for(size_t i = p; i < e; ++ i) {
						if(size_t(m_p_matrix->i[i]) == m_n_row) {
							m_n_row_pointer = i;
							break;
						}
					}
					// e.g. cs_multiply() produces matrices where i is not sorted :(
				}
			}
		}

		/**
		 *	@brief constructor
		 *	@param[in] p_matrix is a sparse matrix to iterate over
		 *	@param[in] n_rows is override for the number of matrix rows
		 *	@param[in] n_cols is override for the number of matrix columns
		 */
		CSparseMatrixShapedIterator(const cs *p_matrix, size_t n_rows, size_t n_cols)
			:m_p_matrix(p_matrix), m_n_rows(n_rows), m_n_columns(n_cols),
			m_n_column(0), m_n_row(0), m_n_row_pointer(p_matrix->p[0]),
			m_n_row_pointer_column_end(p_matrix->p[1]),
			m_n_row_pointer_end(p_matrix->p[p_matrix->n])
		{
			if(m_n_row_pointer <= m_n_row_pointer_column_end) {
				size_t n_data_row = m_p_matrix->i[m_n_row_pointer];
				if(m_n_row_pointer == m_n_row_pointer_column_end || n_data_row != m_n_row) {
					size_t p = (m_n_column < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column] : m_n_row_pointer_end,
						e = (m_n_column + 1 < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column + 1] : m_n_row_pointer_end;
					for(size_t i = p; i < e; ++ i) {
						if(size_t(m_p_matrix->i[i]) == m_n_row) {
							m_n_row_pointer = i;
							break;
						}
					}
					// e.g. cs_multiply() produces matrices where i is not sorted :(
				}
			}
		}

		/**
		 *	@brief gets current row
		 *	@return Returns current row this iterator points to.
		 */
		size_t n_Row() const
		{
			return m_n_row;
		}

		/**
		 *	@brief gets current column
		 *	@return Returns current column this iterator points to.
		 */
		size_t n_Column() const
		{
			return m_n_column;
		}

		/**
		 *	@brief dereferences the iterator
		 *	@return Returns value under the iterator or 0 in case there
		 *		is no value associated with the current row / column (the matrix is sparse).
		 */
		double operator *() const
		{
			_ASSERTE(m_n_column < m_n_columns && m_n_row < m_n_rows); // make sure the iterator is dereferencable
			if(m_n_row_pointer < m_n_row_pointer_column_end) {
				_ASSERTE(unsigned(m_p_matrix->i[m_n_row_pointer]) >= m_n_row);
				if(size_t(m_p_matrix->i[m_n_row_pointer]) == m_n_row)
					return (m_p_matrix->x)? m_p_matrix->x[m_n_row_pointer] : 1;
			}
			return 0; // no data (outside the matrix or being in "sparse area")
		}

		/**
		 *	@brief dereferences the iterator
		 *	@param[out] r_b_value is set in case there is a value associated
		 *		with the current row / column, otherwise it is cleared
		 *	@return Returns value under the iterator or 0 in case there
		 *		is no value associated with the current row / column (the matrix is sparse).
		 */
		double f_Get(bool &r_b_value) const
		{
			_ASSERTE(m_n_column < m_n_columns && m_n_row < m_n_rows); // make sure the iterator is dereferencable
			if(m_n_row_pointer < m_n_row_pointer_column_end) {
				_ASSERTE(unsigned(m_p_matrix->i[m_n_row_pointer]) >= m_n_row);
				if(size_t(m_p_matrix->i[m_n_row_pointer]) == m_n_row) {
					r_b_value = true;
					return (m_p_matrix->x)? m_p_matrix->x[m_n_row_pointer] : 1;
				}
			}
			r_b_value = false;
			return 0; // no data (outside the matrix or being in "sparse area")
		}

		/**
		 *	@brief increments the iterator (to the next position, regardless of matrix sparsity)
		 */
		void operator ++()
		{
			_ASSERTE(m_n_column < m_n_columns);
			_ASSERTE(m_n_row < m_n_rows); // make sure we're not iterating too far

			if(m_n_column >= unsigned(m_p_matrix->n)) {
				if(++ m_n_row == m_n_rows) {
					m_n_row = 0;
					++ m_n_column;
				}
				// we are below the matrix data, just iterating zeroes
			} else if(m_n_row >= unsigned(m_p_matrix->m)) {
				if(++ m_n_row == m_n_rows) {
					m_n_row = 0;
					++ m_n_column;
					_ASSERTE(m_n_column <= unsigned(m_p_matrix->n));
					/*if(m_n_column <= m_p_matrix->n) */ { // always true // m_p_matrix->n is a valid index in p (sets m_n_row_pointer_offset equal to m_n_row_pointer_end)
						m_n_row_pointer = m_p_matrix->p[m_n_column];
						m_n_row_pointer_column_end = (m_n_column < unsigned(m_p_matrix->n))?
							m_p_matrix->p[m_n_column + 1] : m_n_row_pointer_end;
					}
				}
				// we are right of the matrix data, just iterating zeroes (the same code)
			} else {
				if(++ m_n_row == m_n_rows) {
					m_n_row = 0;
					++ m_n_column;
					_ASSERTE(m_n_column <= unsigned(m_p_matrix->n));
					/*if(m_n_column <= m_p_matrix->n) */ { // always true // m_p_matrix->n is a valid index in p (sets m_n_row_pointer_offset equal to m_n_row_pointer_end)
						m_n_row_pointer = m_p_matrix->p[m_n_column];
						m_n_row_pointer_column_end = (m_n_column < unsigned(m_p_matrix->n))?
							m_p_matrix->p[m_n_column + 1] : m_n_row_pointer_end;
					}
				}
				// shift to the next row / column

				size_t p = (m_n_column < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column] : m_n_row_pointer_end,
					e = (m_n_column + 1 < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column + 1] : m_n_row_pointer_end;
				_ASSERTE(e >= p);
				_ASSERTE(m_n_row_pointer_column_end == e);
				_ASSERTE(m_n_row_pointer >= p && m_n_row_pointer <= e); // note that m_n_row_pointer == e is a valid state in case there is not enough data on the row
				// we are inside the matrix data; just need to find row in the CSC array

				if(m_n_row_pointer <= m_n_row_pointer_column_end) { // have to check if we hit the last nonzero row in the column
					size_t n_data_row = m_p_matrix->i[m_n_row_pointer];
					while(n_data_row < m_n_row && m_n_row_pointer < m_n_row_pointer_column_end)
						n_data_row = m_p_matrix->i[++ m_n_row_pointer];
					if(m_n_row_pointer == m_n_row_pointer_column_end || n_data_row != m_n_row) {
						for(size_t i = p; i < e; ++ i) {
							if(size_t(m_p_matrix->i[i]) == m_n_row) {
								m_n_row_pointer = i;
								break;
							}
						}
						// e.g. cs_multiply() produces matrices where i is not sorted :(
					}
				}
				// makes sure that m_n_row_pointer points to an element at the current or greater
			}

			_ASSERTE((m_n_column < m_n_columns && m_n_row < m_n_rows) ||
				(m_n_column == m_n_columns && !m_n_row)); // make sure we are still inside, or just outside
			if(m_n_column < unsigned(m_p_matrix->n)) {
#ifdef _DEBUG
				size_t p = (m_n_column < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column] : m_n_row_pointer_end,
					e = (m_n_column + 1 < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column + 1] : m_n_row_pointer_end;
#endif // _DEBUG
				_ASSERTE(m_n_row_pointer_column_end == e);
				if(m_n_row < unsigned(m_p_matrix->m)) {
					_ASSERTE(m_n_row_pointer == e ||
						m_n_row <= unsigned(m_p_matrix->i[m_n_row_pointer]));
					// either we are at the end of data, or m_n_row_pointer points at or after current row
				} else {
					_ASSERTE(m_n_row_pointer == e);
					// we are at the end of the row for sure
				}
			} else {
				_ASSERTE(m_n_row_pointer == m_n_row_pointer_end);
				_ASSERTE(m_n_row_pointer_column_end == m_n_row_pointer_end);
				// we are at the end
			}
			// iterator integrity check
		}
	};

	/**
	 *	@brief simple sparse matrix iterator
	 *	@todo Test this one (been using shaped so far).
	 */
	class CSparseMatrixIterator {
	protected:
		const cs *m_p_matrix; /**< @brief the matrix in question */
		size_t m_n_column; /**< @brief current column */
		size_t m_n_row; /**< @brief current row */
		size_t m_n_row_pointer; /**< @brief current row pointer offset (index in m_p_matrix->p) */
		size_t m_n_row_pointer_column_end; /**< @brief maximum value m_n_row_pointer can take at the current column */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] p_matrix is a sparse matrix to iterate over
		 */
		CSparseMatrixIterator(const cs *p_matrix)
			:m_p_matrix(p_matrix), m_n_column(0), m_n_row(0), m_n_row_pointer(p_matrix->p[0]),
			m_n_row_pointer_column_end(p_matrix->p[1])
		{
			if(m_n_row_pointer <= m_n_row_pointer_column_end) {
				size_t n_data_row = m_p_matrix->i[m_n_row_pointer];
				if(m_n_row_pointer == m_n_row_pointer_column_end || n_data_row != m_n_row) {
					size_t p = m_p_matrix->p[m_n_column], e = (m_n_column + 1 < size_t(m_p_matrix->n))? m_p_matrix->p[m_n_column + 1] : p;
					for(size_t i = p; i < e; ++ i) {
						if(size_t(m_p_matrix->i[i]) == m_n_row) {
							m_n_row_pointer = i;
							break;
						}
					}
					// e.g. cs_multiply() produces matrices where i is not sorted :(
				}
			}
		}

		/**
		 *	@copydoc CSparseMatrixShapedIterator::n_Row()
		 */
		size_t n_Row() const
		{
			return m_n_row;
		}

		/**
		 *	@copydoc CSparseMatrixShapedIterator::n_Column()
		 */
		size_t n_Column() const
		{
			return m_n_column;
		}

		/**
		 *	@copydoc CSparseMatrixShapedIterator::operator *()
		 */
		double operator *() const
		{
			_ASSERTE(m_n_column < unsigned(m_p_matrix->n) && m_n_row < unsigned(m_p_matrix->m)); // make sure the iterator is dereferencable
			if(m_n_row_pointer < m_n_row_pointer_column_end) {
				_ASSERTE(unsigned(m_p_matrix->i[m_n_row_pointer]) >= m_n_row);
				if(size_t(m_p_matrix->i[m_n_row_pointer]) == m_n_row)
					return (m_p_matrix->x)? m_p_matrix->x[m_n_row_pointer] : 1;
			}
			return 0; // no data (being in "sparse area")
		}

		/**
		 *	@copydoc CSparseMatrixShapedIterator::operator ++()
		 */
		void operator ++()
		{
			_ASSERTE(m_n_column < unsigned(m_p_matrix->n));
			_ASSERTE(m_n_row < unsigned(m_p_matrix->m)); // make sure we're not iterating too far

			{
				if(++ m_n_row == size_t(m_p_matrix->m)) {
					m_n_row = 0;
					++ m_n_column;
					_ASSERTE(m_n_column <= unsigned(m_p_matrix->n));
					/*if(m_n_column <= m_p_matrix->n) */ { // always true // m_p_matrix->n is a valid index in p (sets m_n_row_pointer_offset equal to the ebd)
						m_n_row_pointer = m_p_matrix->p[m_n_column];
						m_n_row_pointer_column_end = (m_n_column < unsigned(m_p_matrix->n))?
							m_p_matrix->p[m_n_column + 1] : m_n_row_pointer; // m_n_row_pointer == p_matrix->p[p_matrix->n] by then
						_ASSERTE(m_n_column < unsigned(m_p_matrix->n) || m_n_row_pointer == m_p_matrix->p[m_p_matrix->n]); // just to make sure
					}
				}
				// shift to the next row / column

				_ASSERTE(m_n_column <= unsigned(m_p_matrix->n));
				size_t p = m_p_matrix->p[m_n_column], e = (m_n_column + 1 < unsigned(m_p_matrix->n))? m_p_matrix->p[m_n_column + 1] : p;
				_ASSERTE(e >= p);
				_ASSERTE(m_n_row_pointer_column_end == e);
				_ASSERTE(m_n_row_pointer >= p && m_n_row_pointer <= e); // note that m_n_row_pointer == e is a valid state in case there is not enough data on the row
				// we are inside the matrix data; just need to find row in the CSC array

				if(m_n_row_pointer <= m_n_row_pointer_column_end) { // have to check if we hit the last nonzero row in the column
					size_t n_data_row = m_p_matrix->i[m_n_row_pointer];
					while(n_data_row < m_n_row && m_n_row_pointer < m_n_row_pointer_column_end)
						n_data_row = m_p_matrix->i[++ m_n_row_pointer];
					if(m_n_row_pointer == m_n_row_pointer_column_end || n_data_row != m_n_row) {
						for(size_t i = p; i < e; ++ i) {
							if(size_t(m_p_matrix->i[i]) == m_n_row) {
								m_n_row_pointer = i;
								break;
							}
						}
						// e.g. cs_multiply() produces matrices where i is not sorted :(
					}
				}
				// makes sure that m_n_row_pointer points to an element at the current or greater
			}

			_ASSERTE((m_n_column < unsigned(m_p_matrix->n) && m_n_row < unsigned(m_p_matrix->m)) ||
				(m_n_column == m_p_matrix->n && !m_n_row)); // make sure we are still inside, or just outside
			if(m_n_column < unsigned(m_p_matrix->n)) {
				_ASSERTE(m_n_column <= unsigned(m_p_matrix->n));
#ifdef _DEBUG
				size_t p = m_p_matrix->p[m_n_column], e = (m_n_column + 1 < size_t(m_p_matrix->n))? m_p_matrix->p[m_n_column + 1] : p;
#endif // _DEBUG
				_ASSERTE(m_n_row_pointer_column_end == e);
				if(m_n_row < unsigned(m_p_matrix->m)) {
					_ASSERTE(m_n_row_pointer == e ||
						m_n_row <= unsigned(m_p_matrix->i[m_n_row_pointer]));
					// either we are at the end of data, or m_n_row_pointer points at or after current row
				} else {
					_ASSERTE(m_n_row_pointer == e);
					// we are at the end of the row for sure
				}
			} else {
				_ASSERTE(m_n_row_pointer == m_n_row_pointer_column_end);
				_ASSERTE(m_n_row_pointer_column_end == m_p_matrix->p[m_p_matrix->n]);
				// we are at the end
			}
			// iterator integrity check
		}
	};

	/**
	 *	@brief rasterizes a sparse matrix and saves as a .tga images
	 *
	 *	@param[in] p_s_filename is output file name
	 *	@param[in] A is input matrix
	 *	@param[in] A_prev is previous state of the input matrix, for change tracking (can be 0)
	 *	@param[in] n_scalar_size is size of scalar, in pixels
	 *
	 *	@return Returns true on success, false on failure.
	 */
	static bool Dump_SparseMatrix(const char *p_s_filename, const cs *A,
		const cs *A_prev = 0, int n_scalar_size = 5)
	{
		_ASSERTE(CS_CSC(A));
		_ASSERTE(!A_prev || CS_CSC(A_prev));
		// the matrices need to be in compressed column form

		const uint32_t n_background = 0xffffffffU,
			n_nonzero = 0xff8080ffU,
			n_nonzero_new = 0xffff0000U,
			n_zero_new = 0xff00ff00U;
		// colors of the individual areas

		const uint32_t n_nonzero_border = 0x80000000U | ((n_nonzero >> 1) & 0x7f7f7f7fU),
			n_nonzero_new_border = 0x80000000U | ((n_nonzero_new >> 1) & 0x7f7f7f7fU),
			n_zero_new_border = 0x80000000U | ((n_zero_new >> 1) & 0x7f7f7f7fU);
		// colors of borders (just darker)

		size_t m = (A_prev)? std::max(A->m, A_prev->m) : A->m;
		size_t n = (A_prev)? std::max(A->n, A_prev->n) : A->n; // in fact, it's the other way around, but lets not confuse things

		if(m == SIZE_MAX || n == SIZE_MAX || m > INT_MAX || n > INT_MAX ||
		   m * (n_scalar_size - 1) + 1 > INT_MAX || n * (n_scalar_size - 1) + 1 > INT_MAX ||
		   uint64_t(m * (n_scalar_size - 1) + 1) * (m * (n_scalar_size - 1) + 1) > INT_MAX)
			return false;

		TBmp *p_bitmap;

		if(!(p_bitmap = TBmp::p_Alloc(int(n * (n_scalar_size - 1) + 1), int(m * (n_scalar_size - 1) + 1))))
			return false;
		p_bitmap->Clear(n_background);

		if(A_prev) {
			CSparseMatrixShapedIterator p_a_it(A, m, n);
			CSparseMatrixShapedIterator p_a_prev_it(A_prev, m, n);
			for(size_t n_col = 0; n_col < n; ++ n_col) {
				for(size_t n_row = 0; n_row < m; ++ n_row, ++ p_a_it, ++ p_a_prev_it) {
					_ASSERTE(n_row == p_a_it.n_Row() && n_col == p_a_it.n_Column());
					_ASSERTE(n_row == p_a_prev_it.n_Row() && n_col == p_a_prev_it.n_Column());
					// make sure it iterates to the right position

					bool b_is;
					double f_value = p_a_it.f_Get(b_is);// = *p_a_it;
					double f_prev_value = *p_a_prev_it;
					// read the value

					double f_value_thr = (fabs(f_value) < 1e-10)? 0 : f_value;
					double f_prev_value_thr = (fabs(f_prev_value) < 1e-10)? 0 : f_prev_value;
					// threshold

					if(b_is || f_value != 0 || f_prev_value != 0) {
						uint32_t n_fill_color = (f_value_thr != 0 && f_prev_value_thr != 0)? n_nonzero :
											    (f_value_thr != 0 && f_prev_value_thr == 0)? n_nonzero_new :
											    (f_value != 0 && f_prev_value != 0)? n_nonzero : n_zero_new;
						uint32_t n_line_color = (f_value_thr != 0 && f_prev_value_thr != 0)? n_nonzero_border :
											    (f_value_thr != 0 && f_prev_value_thr == 0)? n_nonzero_new_border :
												(f_value != 0 && f_prev_value != 0)? n_nonzero_border : n_zero_new_border;

						size_t n_x = n_col * (n_scalar_size - 1);
						size_t n_y = n_row * (n_scalar_size - 1);
						for(int dy = 1; dy < n_scalar_size - 1; ++ dy)
							for(int dx = 1; dx < n_scalar_size - 1; ++ dx)
								p_bitmap->PutPixel(int(n_x + dx), int(n_y + dy), n_fill_color);
						p_bitmap->DrawLine_SP(float(n_x), float(n_y), float(n_x + n_scalar_size - 1), float(n_y), n_line_color);
						p_bitmap->DrawLine_SP(float(n_x), float(n_y + n_scalar_size - 1), float(n_x + n_scalar_size - 1),
							float(n_y + n_scalar_size - 1), n_line_color);
						p_bitmap->DrawLine_SP(float(n_x), float(n_y), float(n_x), float(n_y + n_scalar_size - 1), n_line_color);
						p_bitmap->DrawLine_SP(float(n_x + n_scalar_size - 1), float(n_y), float(n_x + n_scalar_size - 1),
							float(n_y + n_scalar_size - 1), n_line_color);
					}
					// draw a square into the bitmap
				}
			}
		} else {
			CSparseMatrixShapedIterator p_a_it(A, m, n);
			for(size_t n_col = 0; n_col < n; ++ n_col) {
				for(size_t n_row = 0; n_row < m; ++ n_row, ++ p_a_it) {
					_ASSERTE(n_row == p_a_it.n_Row() && n_col == p_a_it.n_Column());
					// make sure it iterates to the right position

					double f_value = *p_a_it;
					// read the value

					if(f_value != 0) {
						size_t n_x = n_col * (n_scalar_size - 1);
						size_t n_y = n_row * (n_scalar_size - 1);
						for(int dy = 1; dy < n_scalar_size - 1; ++ dy)
							for(int dx = 1; dx < n_scalar_size - 1; ++ dx)
								p_bitmap->PutPixel(int(n_x + dx), int(n_y + dy), n_nonzero);
						p_bitmap->DrawLine_SP(float(n_x), float(n_y), float(n_x + n_scalar_size - 1), float(n_y), n_nonzero_border);
						p_bitmap->DrawLine_SP(float(n_x), float(n_y + n_scalar_size - 1), float(n_x + n_scalar_size - 1),
							float(n_y + n_scalar_size - 1), n_nonzero_border);
						p_bitmap->DrawLine_SP(float(n_x), float(n_y), float(n_x), float(n_y + n_scalar_size - 1), n_nonzero_border);
						p_bitmap->DrawLine_SP(float(n_x + n_scalar_size - 1), float(n_y), float(n_x + n_scalar_size - 1),
							float(n_y + n_scalar_size - 1), n_nonzero_border);
					}
					// draw a square into the bitmap
				}
			}
		}
		// iterate the matrix, draw small tiles

		bool b_result = CTgaCodec::Save_TGA(p_s_filename, *p_bitmap, false);

		p_bitmap->Delete();

		/*for(int j = 0; j < n; ++ j) {
			//printf("    col %d\n", j);
			int y = 0;
			if(j < A->n) {
				for(int p = A->p[j]; p < A->p[j + 1]; ++ p, ++ y) {
					int n_cur_y = A->i[p];
					for(; y < n_cur_y; ++ y)
						;//printf(" %f, ", .0f); // zero entries
					double f_value = (A->x)? A->x[p] : 1; // nonzero entry
				}
				for(; y < m; ++ y)
					;//printf(" %f%s", .0f, (x + 1 == m)? "\n" : ", "); // zero entries
			}
		}*/
		// t_odo - make this an iterator anyway

		// t_odo - rasterize the matrix (A_prev is used to mark difference areas)
		// todo - implement libPNG image writing

		return b_result;
	}

	/**
	 *	@brief prints sparse matrix as a dense matrix
	 *
	 *	@param[in] A is the matrix to be printed
	 *	@param[in] p_s_label is the name of the matrix (can be null)
	 */
	static void Print_SparseMatrix(const cs *A, const char *p_s_label = 0)
	{
		if(p_s_label)
			printf("matrix \'%s\' = ", p_s_label);
		csi m, n, nzmax, /*nz,*/ *Ap, *Ai;
		double *Ax;
		if(!A) {
			printf("(null)\n");
			return;
		}
		cs *Ac = 0;
		if(CS_TRIPLET(A) && !(Ac = cs_compress(A))) {
			printf("[not enough memory]\n");
			return;
		}
		cs *At;
		if(!(At = cs_transpose((Ac)? Ac : A, 1))) {
			if(Ac)
				cs_spfree(Ac);
			printf("[not enough memory]\n");
			return;
		}
		m = At->m; n = At->n; Ap = At->p; Ai = At->i; Ax = At->x;
		nzmax = At->nzmax; //nz = At->nz; // unused

		_ASSERTE(CS_CSC(At));

		{
			_ASSERTE(n <= INT_MAX);
			_ASSERTE(m <= INT_MAX);
			_ASSERTE(nzmax <= INT_MAX);
			_ASSERTE(Ap[n] <= INT_MAX);
			printf("%d-by-%d, nzmax: %d nnz: %d, 1-norm: %g\n", int(n), int(m), int(nzmax),
					int(Ap[n]), cs_norm(A));
			for(csi j = 0; j < n; ++ j) {
				printf("    row %d\n", int(j));
				int x = 0;
				for(csi p = Ap[j]; p < Ap[j + 1]; ++ p, ++ x) {
					csi n_cur_x = Ai[p];
					for(; x < n_cur_x; ++ x)
						printf(" %f, ", .0f);
					printf("%s%f%s", (Ax && Ax[p] < 0)? "" : " ", (Ax)? Ax[p] : 1, (n_cur_x + 1 == m)? "\n" : ", ");
				}
				for(; x < m; ++ x)
					printf(" %f%s", .0f, (x + 1 == m)? "\n" : ", ");
			}
		}

		cs_spfree(At);
		if(Ac)
			cs_spfree(Ac);
	}

	/**
	 *	@brief prints sparse matrix as a dense matrix, in matlab format
	 *
	 *	@param[in] A is the matrix to be printed
	 *	@param[in] p_s_label is the name of the matrix (can be null)
	 */
	static void Print_SparseMatrix_in_MatlabFormat(const cs *A, const char *p_s_label = 0)
	{
		if(p_s_label)
			printf("matrix \'%s\' = ", p_s_label);
		csi m, n, /*nzmax, nz,*/ *Ap, *Ai;
		double *Ax;
		if(!A) {
			printf("(null)\n");
			return;
		}
		cs *Ac = 0;
		if(CS_TRIPLET(A) && !(Ac = cs_compress(A))) {
			printf("[not enough memory]\n");
			return;
		}
		cs *At;
		if(!(At = cs_transpose((Ac)? Ac : A, 1))) {
			if(Ac)
				cs_spfree(Ac);
			printf("[not enough memory]\n");
			return;
		}
		m = At->m; n = At->n; Ap = At->p; Ai = At->i; Ax = At->x;
		//nzmax = At->nzmax; nz = At->nz;

		_ASSERTE(CS_CSC(At));

		{
			printf("[");
			for(csi j = 0; j < n; ++ j) {
				csi x = 0;
				for(csi p = Ap[j]; p < Ap[j + 1]; ++ p, ++ x) {
					csi n_cur_x = Ai[p];
					for(; x < n_cur_x; ++ x)
						printf("%f ", .0f);
					printf("%f%s", (Ax)? Ax[p] : 1, (n_cur_x + 1 == m)? ((j + 1 == n)? "" : "; ") : " ");
				}
				for(; x < m; ++ x)
					printf("%f%s", .0f, (x + 1 == m)? ((j + 1 == n)? "" : "; ") : " ");
			}
			printf("]\n");
		}

		cs_spfree(At);
		if(Ac)
			cs_spfree(Ac);
	}

	/**
	 *	@brief prints sparse matrix as a dense matrix, in matlab format
	 *
	 *	@param[in] p_fw is the output stream
	 *	@param[in] A is the matrix to be printed
	 *	@param[in] p_s_prefix is prefix code before the matrix body (can be null)
	 *	@param[in] p_s_suffix is suffix code after the matrix body (can be null)
	 *
	 *	@return Returns true on success, false on failure.
	 */
	static bool Print_SparseMatrix_in_MatlabFormat(FILE *p_fw,
		const cs *A, const char *p_s_prefix = 0, const char *p_s_suffix = ";\n")
	{
		if(p_s_prefix)
			fprintf(p_fw, "%s", p_s_prefix);

		{
			fprintf(p_fw, "[");
			CSparseMatrixShapedIterator it(A); // @todo - replace this by a simple iterator as soon as it's tested
			for(size_t i = 0, m = A->m; i < m; ++ i) { // iterate rows in the outer loop
				for(size_t j = 0, n = A->n; j < n; ++ j, ++ it) { // iterate columns in the inner loop
					double f_value = *it;
					fprintf(p_fw, " %f", f_value);
				}
				if(i + 1 != m)
					fprintf(p_fw, "; ");
			}
			fprintf(p_fw, "]");
		}

		// "to define a matrix, you can treat it like a column of row vectors"
		/*
		 *	>> A = [ 1 2 3; 3 4 5; 6 7 8]
		 *
		 *	A =
		 *
		 *		 1     2     3
		 *		 3     4     5
		 *		 6     7     8
		 */

		if(p_s_suffix)
			fprintf(p_fw, "%s", p_s_suffix);

		return true;
	}

	/**
	 *	@brief prints dense matrix dimensions
	 *
	 *	@tparam MatrixType is specialization of Eigen::Matrix
	 *
	 *	@param[in] p_fw is the output stream
	 *	@param[in] r_A is the matrix to be printed
	 *	@param[in] p_s_prefix is prefix code before the matrix body (can be null)
	 *	@param[in] p_s_suffix is suffix code after the matrix body (can be null)
	 */
	template <class MatrixType>
	static void Print_DenseMatrix_Dimensions(FILE *p_fw, const MatrixType &r_A,
		const char *p_s_prefix = 0, const char *p_s_suffix = "\n")
	{
		if(p_s_prefix)
			fprintf(p_fw, "%s", p_s_prefix);
		fprintf(p_fw, "dense, " PRIsize " x " PRIsize, r_A.rows(), r_A.cols());
		if(p_s_suffix)
			fprintf(p_fw, "%s", p_s_suffix);
	}

	/**
	 *	@brief prints dense matrix in matlab format
	 *
	 *	@tparam MatrixType is specialization of Eigen::Matrix
	 *
	 *	@param[in] p_fw is the output stream
	 *	@param[in] r_A is the matrix to be printed
	 *	@param[in] p_s_prefix is prefix code before the matrix body (can be null)
	 *	@param[in] p_s_suffix is suffix code after the matrix body (can be null)
	 *	@param[in] p_s_fmt is format string (" %f" by default; the space is important)
	 *
	 *	@return Returns true on success, false on failure.
	 */
	template <class MatrixType>
	static bool Print_DenseMatrix_in_MatlabFormat(FILE *p_fw,
		const MatrixType &r_A, const char *p_s_prefix = 0,
		const char *p_s_suffix = ";\n", const char *p_s_fmt = " %f")
	{
		if(p_s_prefix)
			fprintf(p_fw, "%s", p_s_prefix);

		{
			fprintf(p_fw, "[");
			for(size_t i = 0, m = r_A.rows(); i < m; ++ i) { // iterate rows in the outer loop
				for(size_t j = 0, n = r_A.cols(); j < n; ++ j) { // iterate columns in the inner loop
					double f_value = r_A(i, j);
					fprintf(p_fw, p_s_fmt, f_value);
				}
				if(i + 1 != m)
					fprintf(p_fw, "; ");
			}
			fprintf(p_fw, "]");
		}

		// "to define a matrix, you can treat it like a column of row vectors"
		/*
		 *	>> A = [ 1 2 3; 3 4 5; 6 7 8]
		 *
		 *	A =
		 *
		 *		 1     2     3
		 *		 3     4     5
		 *		 6     7     8
		 */

		if(p_s_suffix)
			fprintf(p_fw, "%s", p_s_suffix);

		return true;
	}

	/**
	 *	@brief prints sparse matrix as a true sparse matrix,
	 *		in matlab format (requires suitesparse)
	 *
	 *	@param[in] p_fw is the output stream
	 *	@param[in] A is the matrix to be printed
	 *	@param[in] p_s_label is a name of the matrix variable
	 *	@param[in] p_s_prefix is prefix code before the matrix body (can be null)
	 *	@param[in] p_s_suffix is suffix code on the next
	 *		line after the matrix body (can be null)
	 *
	 *	@return Returns true on success, false on failure.
	 */
	static bool Print_SparseMatrix_in_MatlabFormat2(FILE *p_fw, const cs *A,
		const char *p_s_label, const char *p_s_prefix = 0, const char *p_s_suffix = 0)
	{
		if(p_s_prefix)
			fprintf(p_fw, "%s", p_s_prefix);

		fprintf(p_fw, "%s = sparse(" PRIdiff ", " PRIdiff ");\n", p_s_label, A->m, A->n);
		for(csi i = 0; i < A->n; ++ i) {
			for(csi p = A->p[i], e = A->p[i + 1]; p < e; ++ p) {
				csi j = A->i[p];
				double x = (A->x)? A->x[p] : 1;
				if(x != 0) { // hard comparison here, try to avoid zero lines introduced by conversion from block matrices
					fprintf(p_fw, (fabs(x) > 1)? "%s(" PRIdiff ", " PRIdiff ") = %f;\n" :
						"%s(" PRIdiff ", " PRIdiff ") = %g;\n", p_s_label, j + 1, i + 1, x);
				}
			}
		}

		// to make sparse unit matrix, one writes:
		/*
		 *	>> A = sparse(2, 2);
		 *	>> A(1, 1) = 1
		 *	>> A(2, 2) = 1
		 *
		 *	A =
		 *		 (1, 1)     1
		 *		 (2, 2)     1
		 */

		if(p_s_suffix)
			fprintf(p_fw, "%s", p_s_suffix);

		return true;
	}

	/**
	 *	@brief prints a dense vector
	 *
	 *	@param[in] b is the vector to be printed
	 *	@param[in] n_vector_length is number of elements in b
	 *	@param[in] p_s_label is the name of the matrix (can be null)
	 */
	static void Print_DenseVector(const double *b, size_t n_vector_length, const char *p_s_label = 0)
	{
		//_ASSERTE(n_vector_length > 0); // doesn't work for empties

		double f_min_abs = (n_vector_length)? fabs(b[0]) : 0;
		for(size_t i = 1; i < n_vector_length; ++ i)
			f_min_abs = (fabs(b[i]) > 0 && (!f_min_abs || f_min_abs > fabs(b[i])))? fabs(b[i]) : f_min_abs;
		int n_log10 = (f_min_abs > 0)? int(ceil(log(f_min_abs) / log(10.0))) : 0;
		// calculate log10 to display the smallest nonzero value as 0 to 1 number

		if(n_log10 < -1)
			n_log10 += 2;
		else if(n_log10 < 0)
			++ n_log10;
		// it's ok to have two first places 0 (in matlab it is, apparently)

		double f_scale = pow(10.0, -n_log10);
		// calculatze scale

		if(p_s_label)
			printf("%s = ", p_s_label);
		if(n_vector_length) {
			printf("vec(" PRIsize ") = 1e%+d * [%.4f", n_vector_length, n_log10, f_scale * b[0]);
			for(size_t i = 1; i < n_vector_length; ++ i)
				printf(", %.4f", f_scale * b[i]);
			printf("]\n");
		} else
			printf("[ ]\n");
		// print with scale
	}

	/**
	 *	@brief prints a dense vector in matlab format
	 *
	 *	@param[in] p_fw is a file where to save the vector
	 *	@param[in] b is the vector to be printed
	 *	@param[in] n_vector_length is number of elements in b
	 *	@param[in] p_s_prefix is prefix code before the vector body (can be null)
	 *	@param[in] p_s_suffix is suffix code after the vector body (can be null)
	 */
	static void Print_DenseVector_in_MatlabFormat(FILE *p_fw, const double *b,
		size_t n_vector_length, const char *p_s_prefix = 0, const char *p_s_suffix = ";\n")
	{
		//_ASSERTE(n_vector_length > 0); // doesn't work for empties

		double f_min_abs = (n_vector_length)? fabs(b[0]) : 0;
		for(size_t i = 1; i < n_vector_length; ++ i)
			f_min_abs = (fabs(b[i]) > 0 && (!f_min_abs || f_min_abs > fabs(b[i])))? fabs(b[i]) : f_min_abs;
		int n_log10 = (f_min_abs > 0)? int(ceil(log(f_min_abs) / log(10.0))) : 0;
		// calculate log10 to display the smallest nonzero value as 0 to 1 number

		if(n_log10 < -1)
			n_log10 += 2;
		else if(n_log10 < 0)
			++ n_log10;
		// it's ok to have two first places 0 (in matlab it is, apparently)

		double f_scale = pow(10.0, -n_log10);
		// calculatze scale

		if(p_s_prefix)
			fprintf(p_fw, "%s", p_s_prefix);
		if(n_vector_length) {
			if(fabs(f_scale * b[0]) > 1)
				fprintf(p_fw, "1e%+d * [%f", n_log10, f_scale * b[0]);
			else
				fprintf(p_fw, "1e%+d * [%g", n_log10, f_scale * b[0]);
			for(size_t i = 1; i < n_vector_length; ++ i)
				fprintf(p_fw, (fabs(f_scale * b[i]) > 1)? " %f" : " %g", f_scale * b[i]);
			fprintf(p_fw, "]");
		} else
			fprintf(p_fw, "[ ]");
		if(p_s_suffix)
			fprintf(p_fw, "%s", p_s_suffix);
		// print with scale
	}
};

#endif // !__DEBUGGING_FUNCTIONALITY_INCLUDED
