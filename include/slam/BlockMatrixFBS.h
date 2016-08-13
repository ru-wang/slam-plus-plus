/*
								+---------------------------------+
								|                                 |
								|  ***   Über Block Matrix   ***  |
								|                                 |
								| Copyright  (c) -tHE SWINe- 2012 |
								|                                 |
								|        BlockMatrixFBS.h         |
								|                                 |
								+---------------------------------+
*/

#pragma once
#ifndef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_OPS_IMPLEMENTATION_INCLUDED
#define __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_OPS_IMPLEMENTATION_INCLUDED

/**
 *	@file include/slam/BlockMatrixFBS.h
 *	@date 2012
 *	@author -tHE SWINe-
 *	@brief the überblockmatrix fixed block size implementation
 *	@note This file is not to be included; it is automatically included from BlockMatrix.h
 */

#include "slam/BlockMatrixBase.h"

namespace blockmatrix_detail {

/**
 *	@brief a dummy class, containing implementation of the FBS functionality (to be dissolved later)
 */
class CUberBlockMatrix_FBS : public CUberBlockMatrix_Base { // todo un-base, just add typedefs for types that are required below (those are public anyway)
public:
	/**
	 *	@brief fixed block size operations template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 */
	template <class CBlockMatrixTypelistA, class CBlockMatrixTypelistB>
	class CFBS_MatrixMultiply : public fbs_ut::deprecate::CFixedBlockSize_BinaryBase<
		CBlockMatrixTypelistA, CBlockMatrixTypelistB> {
	public:
		typedef typename fbs_ut::deprecate::CFixedBlockSize_BinaryBase<CBlockMatrixTypelistA, CBlockMatrixTypelistB>::_TyBlockMatrixTypelistA _TyBlockMatrixTypelistA; /**< @brief list of block sizes, represented as Eigen::Matrix (of the left-hand-side matrix) */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_BinaryBase<CBlockMatrixTypelistA, CBlockMatrixTypelistB>::_TyBlockMatrixTypelistB _TyBlockMatrixTypelistB; /**< @brief list of block sizes, represented as Eigen::Matrix (of the right-hand-side matrix) */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_BinaryBase<CBlockMatrixTypelistA, CBlockMatrixTypelistB>::CRowHeightsListB CRowHeightsListB; /**< @brief list of unique block row heights (of the right-hand-side matrix) */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_BinaryBase<CBlockMatrixTypelistA, CBlockMatrixTypelistB>::CColumnWidthsListB CColumnWidthsListB; /**< @brief list of unique block column widths (of the right-hand-side matrix) */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_BinaryBase<CBlockMatrixTypelistA, CBlockMatrixTypelistB>::CRowHeightsListA CRowHeightsListA; /**< @brief list of unique block row heights (of the left-hand-side matrix) */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_BinaryBase<CBlockMatrixTypelistA, CBlockMatrixTypelistB>::CColumnWidthsListA CColumnWidthsListA; /**< @brief list of unique block column widths (of the left-hand-side matrix) */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_BinaryBase<CBlockMatrixTypelistA, CBlockMatrixTypelistB>::CDimsList_UniqA CDimsList_UniqA; /**< @brief list of block sizes as CCTSize2D (of the left-hand-side matrix; duplicate records removed) */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_BinaryBase<CBlockMatrixTypelistA, CBlockMatrixTypelistB>::CDimsList_UniqB CDimsList_UniqB; /**< @brief list of block sizes as CCTSize2D (of the right-hand-side matrix; duplicate records removed) */

		/**
		 *	@brief (outer) loop of the matrix-matrix multiplication kernel template
		 *
		 *	@tparam _TyGppContextA is template instance context for g++ (compatibility workarround)
		 *	@tparam _TyGppContextB is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthsListB is list of posible column widths of the right-hand matrix
		 */
		template <class _TyGppContextA, class _TyGppContextB, class _CColumnWidthsListB>
		class CMatrixMultiply_OuterLoop {
		protected:
			typedef typename _CColumnWidthsListB::_TyHead CCurrentColumnWidthB; /**< @brief right-hand matrix column width for this decision tree recursion */

		public:
			/**
			 *	@brief wraps the outer loop of matrix-matrix multiplication in block width decision tree
			 *
			 *	@param[in] r_t_column_B is the current block column of the right-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column of the right-hand matrix
			 *	@param[in] r_row_list_B is list of the right-hand matrix block rows
			 *	@param[in] r_col_list_A is list of the left-hand matrix block columns
			 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
			 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
			 *	@param[in] reindex_rows_B_to_cols_A is the mapping functions from right-hand
			 *		matrix block rows to left-hand matrix block columns
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop(const TColumn &r_t_column_B, size_t n_column_id_B,
				const std::vector<TRow> &r_row_list_B,
				const std::vector<TColumn> &r_col_list_A, const std::vector<TRow> &r_row_list_A,
				std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
				const std::vector<size_t> &reindex_rows_B_to_cols_A, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				if(r_t_column_B.n_width == CCurrentColumnWidthB::n_size) {
					return CMatrixMultiply_OuterLoop<_TyGppContextA, _TyGppContextB,
						CTypelist<CCurrentColumnWidthB, CTypelistEnd> >::Loop(r_t_column_B,
						n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
						transpose_cols_list, reindex_rows_B_to_cols_A, alloc);
				} else {
					return CMatrixMultiply_OuterLoop<_TyGppContextA, _TyGppContextB,
						typename _CColumnWidthsListB::_TyTail>::Loop(r_t_column_B,
						n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
						transpose_cols_list, reindex_rows_B_to_cols_A, alloc);
				}
			}

			/**
			 *	@brief wraps the outer loop of matrix-matrix multiplication in block width decision tree
			 *
			 *	@param[in] r_t_column_B is the current block column of the right-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column of the right-hand matrix
			 *	@param[in] r_row_list_B is list of the right-hand matrix block rows
			 *	@param[in] r_col_list_A is list of the left-hand matrix block columns
			 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
			 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
			 *	@param[in] reindex_rows_B_to_cols_A is the mapping functions from right-hand
			 *		matrix block rows to left-hand matrix block columns
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop_UpperTriag(const TColumn &r_t_column_B, size_t n_column_id_B,
				const std::vector<TRow> &r_row_list_B,
				const std::vector<TColumn> &r_col_list_A, const std::vector<TRow> &r_row_list_A,
				std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
				const std::vector<size_t> &reindex_rows_B_to_cols_A, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				if(r_t_column_B.n_width == CCurrentColumnWidthB::n_size) {
					return CMatrixMultiply_OuterLoop<_TyGppContextA, _TyGppContextB,
						CTypelist<CCurrentColumnWidthB, CTypelistEnd> >::Loop_UpperTriag(r_t_column_B,
						n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
						transpose_cols_list, reindex_rows_B_to_cols_A, alloc);
				} else {
					return CMatrixMultiply_OuterLoop<_TyGppContextA, _TyGppContextB,
						typename _CColumnWidthsListB::_TyTail>::Loop_UpperTriag(r_t_column_B,
						n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
						transpose_cols_list, reindex_rows_B_to_cols_A, alloc);
				}
			}

			// todo - debug loop
		};

		/**
		 *	@brief (outer) loop of the matrix-matrix multiplication kernel template
		 *
		 *	@tparam _TyGppContextA is template instance context for g++ (compatibility workarround)
		 *	@tparam _TyGppContextB is template instance context for g++ (compatibility workarround)
		 *	@tparam CLastColumnWidthB is the chosen right-hand matrix column width
		 */
		template <class _TyGppContextA, class _TyGppContextB, class CLastColumnWidthB>
		class CMatrixMultiply_OuterLoop<_TyGppContextA, _TyGppContextB,
			CTypelist<CLastColumnWidthB, CTypelistEnd> > {
		public:
			typedef typename CFilterTypelist2<CRowHeightsListB, CLastColumnWidthB,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_UniqB>::_TyResult _CSelectedRowHeightsListB; /**< @brief row heights list, filtered by the selected column width */

			/**
			 *	@brief performs the outer loop of matrix-matrix multiplication
			 *
			 *	@param[in] r_t_column_B is the current block column of the right-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column of the right-hand matrix
			 *	@param[in] r_row_list_B is list of the right-hand matrix block rows
			 *	@param[in] r_col_list_A is list of the left-hand matrix block columns
			 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
			 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
			 *	@param[in] reindex_rows_B_to_cols_A is the mapping functions from right-hand
			 *		matrix block rows to left-hand matrix block columns
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop(const TColumn &r_t_column_B, size_t n_column_id_B,
				const std::vector<TRow> &r_row_list_B,
				const std::vector<TColumn> &r_col_list_A, const std::vector<TRow> &r_row_list_A,
				std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
				const std::vector<size_t> &reindex_rows_B_to_cols_A, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				_ASSERTE(r_t_column_B.block_list.empty() || r_t_column_B.n_width == CLastColumnWidthB::n_size);

				for(_TyBlockConstIter p_block_B_it = r_t_column_B.block_list.begin(),
				   p_block_B_end_it = r_t_column_B.block_list.end(); p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
					const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
					size_t n_row_id_B = r_t_block_B.first; // get row of the block
					// for each block in the current column in B

					// the second decision tree for r_row_list_B[n_row_id_B].n_height (decides block B heights, given constant block B width)
					// note that this is the common dimension

					const double *p_block_B = r_t_block_B.second;
					size_t n_row_height_B = r_row_list_B[n_row_id_B].n_height;
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

					//if(r_t_column_A.block_list.empty()) // better enter decision tree; empty columns are rare in realworld applications
					//	continue;
					// otherwise block width mismatch occurs

					if(!CMatrixMultiply_MiddleLoop<_TyGppContextA, _TyGppContextB,
					   _CSelectedRowHeightsListB, CLastColumnWidthB>::Loop(r_t_column_A,
						r_row_list_A, transpose_cols_list, n_row_height_B, n_column_id_B,
						p_block_B, alloc))
						return false;
					// call the middle loop
				}

				return true;
			}

			/**
			 *	@brief performs the outer loop of matrix-matrix multiplication
			 *
			 *	@param[in] r_t_column_B is the current block column of the right-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column of the right-hand matrix
			 *	@param[in] r_row_list_B is list of the right-hand matrix block rows
			 *	@param[in] r_col_list_A is list of the left-hand matrix block columns
			 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
			 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
			 *	@param[in] reindex_rows_B_to_cols_A is the mapping functions from right-hand
			 *		matrix block rows to left-hand matrix block columns
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop_UpperTriag(const TColumn &r_t_column_B, size_t n_column_id_B,
				const std::vector<TRow> &r_row_list_B,
				const std::vector<TColumn> &r_col_list_A, const std::vector<TRow> &r_row_list_A,
				std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
				const std::vector<size_t> &reindex_rows_B_to_cols_A, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				_ASSERTE(r_t_column_B.block_list.empty() || r_t_column_B.n_width == CLastColumnWidthB::n_size);

				for(_TyBlockConstIter p_block_B_it = r_t_column_B.block_list.begin(),
				   p_block_B_end_it = r_t_column_B.block_list.end(); p_block_B_it != p_block_B_end_it; ++ p_block_B_it) {
					const TColumn::TBlockEntry &r_t_block_B = *p_block_B_it;
					size_t n_row_id_B = r_t_block_B.first; // get row of the block
					// for each block in the current column in B

					// the second decision tree for r_row_list_B[n_row_id_B].n_height (decides block B heights, given constant block B width)
					// note that this is the common dimension

					const double *p_block_B = r_t_block_B.second;
					size_t n_row_height_B = r_row_list_B[n_row_id_B].n_height;
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

					//if(r_t_column_A.block_list.empty()) // better enter decision tree; empty columns are rare in realworld applications
					//	continue;
					// otherwise block width mismatch occurs

					if(!CMatrixMultiply_MiddleLoop<_TyGppContextA, _TyGppContextB,
					   _CSelectedRowHeightsListB, CLastColumnWidthB>::Loop_UpperTriag(r_t_column_A,
						r_row_list_A, transpose_cols_list, n_row_height_B, n_column_id_B,
						p_block_B, alloc))
						return false;
					// call the middle loop
				}

				return true;
			}

			// todo - debug loop
		};

		/**
		 *	@brief (middle) loop of the matrix-matrix multiplication kernel template
		 *
		 *	@tparam _TyGppContextA is template instance context for g++ (compatibility workarround)
		 *	@tparam _TyGppContextB is template instance context for g++ (compatibility workarround)
		 *	@tparam _CRowHeightsListB is list of posible row heights of the right-hand matrix
		 *	@tparam CColumnWidthB is the chosen right-hand matrix column width
		 */
		template <class _TyGppContextA, class _TyGppContextB, class _CRowHeightsListB, class CColumnWidthB>
		class CMatrixMultiply_MiddleLoop {
		protected:
			typedef typename _CRowHeightsListB::_TyHead CCurrentRowHeightB; /**< @brief right-hand matrix row height for this decision tree recursion */

		public:
			/**
			 *	@brief wraps the outer loop of matrix-matrix multiplication in a block size decision tree
			 *
			 *	@param[in] r_t_column_A is the current block column of the right-hand matrix
			 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
			 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
			 *	@param[in] n_row_height_B is height of the current block in the right-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column in the right-hand matrix
			 *	@param[in] p_block_B is pointer to dense data of the current block in the right-hand matrix
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop(
				const TColumn &r_t_column_A, const std::vector<TRow> &r_row_list_A,
				std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
				size_t n_row_height_B, size_t n_column_id_B, const double *p_block_B,
				_TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				if(n_row_height_B == CCurrentRowHeightB::n_size) {
					return CMatrixMultiply_MiddleLoop<_TyGppContextA, _TyGppContextB,
						CTypelist<CCurrentRowHeightB, CTypelistEnd>, CColumnWidthB>::Loop(
						r_t_column_A, r_row_list_A, transpose_cols_list, n_row_height_B,
						n_column_id_B, p_block_B, alloc);
				} else {
					return CMatrixMultiply_MiddleLoop<_TyGppContextA, _TyGppContextB,
						typename _CRowHeightsListB::_TyTail, CColumnWidthB>::Loop(
						r_t_column_A, r_row_list_A, transpose_cols_list, n_row_height_B,
						n_column_id_B, p_block_B, alloc);
				}
			}

			/**
			 *	@brief wraps the outer loop of matrix-matrix multiplication in a block size decision tree
			 *
			 *	@param[in] r_t_column_A is the current block column of the right-hand matrix
			 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
			 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
			 *	@param[in] n_row_height_B is height of the current block in the right-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column in the right-hand matrix
			 *	@param[in] p_block_B is pointer to dense data of the current block in the right-hand matrix
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop_UpperTriag(
				const TColumn &r_t_column_A, const std::vector<TRow> &r_row_list_A,
				std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
				size_t n_row_height_B, size_t n_column_id_B, const double *p_block_B,
				_TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				if(n_row_height_B == CCurrentRowHeightB::n_size) {
					return CMatrixMultiply_MiddleLoop<_TyGppContextA, _TyGppContextB,
						CTypelist<CCurrentRowHeightB, CTypelistEnd>, CColumnWidthB>::Loop_UpperTriag(
						r_t_column_A, r_row_list_A, transpose_cols_list, n_row_height_B,
						n_column_id_B, p_block_B, alloc);
				} else {
					return CMatrixMultiply_MiddleLoop<_TyGppContextA, _TyGppContextB,
						typename _CRowHeightsListB::_TyTail, CColumnWidthB>::Loop_UpperTriag(
						r_t_column_A, r_row_list_A, transpose_cols_list, n_row_height_B,
						n_column_id_B, p_block_B, alloc);
				}
			}

			// todo - debug loop
		};

		/**
		 *	@brief (middle) loop of the matrix-matrix multiplication kernel template
		 *
		 *	@tparam _TyGppContextA is template instance context for g++ (compatibility workarround)
		 *	@tparam _TyGppContextB is template instance context for g++ (compatibility workarround)
		 *	@tparam CLastRowHeightB is the chosen right-hand matrix row height
		 *	@tparam CColumnWidthB is the chosen right-hand matrix column width
		 */
		template <class _TyGppContextA, class _TyGppContextB, class CLastRowHeightB, class CColumnWidthB>
		class CMatrixMultiply_MiddleLoop<_TyGppContextA, _TyGppContextB,
			CTypelist<CLastRowHeightB, CTypelistEnd>, CColumnWidthB> {
		protected:
			typedef CLastRowHeightB CLastColumnWidthA; /**< @brief column width of the chosen block in the left-hand matrix (supposed to equal row height of the chosen block in the right-hand matrix if the matrix product is defined) */
			typedef typename CFilterTypelist2<CRowHeightsListA, CLastColumnWidthA,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_UniqA/*2015-03-15 bug*/>::_TyResult _CSelectedRowHeightsListA; /**< @brief left-hand matrix row heights list, filtered by the selected right-hand matrix column width */
			// todo - filter the list by possible B widths? will that work or will that force unmultipliable matrices to crash in release because the corresponding dimension is not found? would have to rewrite a bit, but should be faster in the end

		public:
			/**
			 *	@brief performs the outer loop of matrix-matrix multiplication
			 *
			 *	@param[in] r_t_column_A is the current block column of the right-hand matrix
			 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
			 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
			 *	@param[in] n_row_height_B is height of the current block in the right-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column in the right-hand matrix
			 *	@param[in] p_block_B is pointer to dense data of the current block in the right-hand matrix
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop(
				const TColumn &r_t_column_A, const std::vector<TRow> &r_row_list_A,
				std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
				size_t UNUSED(n_row_height_B), size_t n_column_id_B, const double *p_block_B,
				_TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				_ASSERTE(r_t_column_A.block_list.empty() || n_row_height_B == CLastRowHeightB::n_size);
				// not sure if the r_t_column_A.block_list.empty() is required here

				size_t n_column_width_A = r_t_column_A.n_width;
				for(_TyBlockConstIter p_block_A_it =
				   r_t_column_A.block_list.begin(), p_block_A_end_it = r_t_column_A.block_list.end();
				   p_block_A_it != p_block_A_end_it; ++ p_block_A_it) {
					const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
					size_t n_row_id_A = r_t_block_A.first; // get row of the block
					// for each block in the current column in A

					// the third decision tree for r_row_list_A[n_row_id_A].n_height (decides block A heights, given block A width == block B height)
					// note that r_t_column_A.n_width must equal r_row_list_B[n_row_id_B].n_height (the second decision tree)

					size_t n_row_height_A = r_row_list_A[n_row_id_A].n_height;
					const double *p_block_A = r_t_block_A.second;
					// create map to block A data

					//_ASSERTE(n_row_id_A < r_row_list_dest.size());
					//_ASSERTE(r_row_list_dest[n_row_id_A].n_height == n_row_height_A); // todo - uncomment this if the data is available
					_ASSERTE(n_row_id_A < transpose_cols_list.size());
					std::vector<TColumn::TBlockEntry> &r_transpose_column = transpose_cols_list[n_row_id_A];
					if(!CMatrixMultiply_InnerLoop<_TyGppContextA, _TyGppContextB,
					   _CSelectedRowHeightsListA, CLastRowHeightB, CColumnWidthB>::Loop(
					   r_transpose_column, p_block_A, n_row_height_A, n_column_width_A,
					   n_column_id_B, p_block_B, alloc))
						return false;
					// call the inner loop
				}

				return true;
			}

			/**
			 *	@brief performs the outer loop of matrix-matrix multiplication
			 *
			 *	@param[in] r_t_column_A is the current block column of the right-hand matrix
			 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
			 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
			 *	@param[in] n_row_height_B is height of the current block in the right-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column in the right-hand matrix
			 *	@param[in] p_block_B is pointer to dense data of the current block in the right-hand matrix
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop_UpperTriag(
				const TColumn &r_t_column_A, const std::vector<TRow> &r_row_list_A,
				std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
				size_t UNUSED(n_row_height_B), size_t n_column_id_B, const double *p_block_B,
				_TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				_ASSERTE(r_t_column_A.block_list.empty() || n_row_height_B == CLastRowHeightB::n_size);
				// not sure if the r_t_column_A.block_list.empty() is required here

				size_t n_column_width_A = r_t_column_A.n_width;
				for(_TyBlockConstIter p_block_A_it =
				   r_t_column_A.block_list.begin(), p_block_A_end_it = r_t_column_A.block_list.end();
				   p_block_A_it != p_block_A_end_it; ++ p_block_A_it) {
					const TColumn::TBlockEntry &r_t_block_A = *p_block_A_it;
					size_t n_row_id_A = r_t_block_A.first; // get row of the block
					// for each block in the current column in A

					if(n_row_id_A > n_column_id_B) // note that this could be replaced by std::upper_bound on blocks // todo
						break; // these only get higher, might as well break instead of continue
					// this block would've ended up at lower diagonal, don't want that

					// the third decision tree for r_row_list_A[n_row_id_A].n_height (decides block A heights, given block A width == block B height)
					// note that r_t_column_A.n_width must equal r_row_list_B[n_row_id_B].n_height (the second decision tree)

					size_t n_row_height_A = r_row_list_A[n_row_id_A].n_height;
					const double *p_block_A = r_t_block_A.second;
					// create map to block A data

					//_ASSERTE(n_row_id_A < r_row_list_dest.size());
					//_ASSERTE(r_row_list_dest[n_row_id_A].n_height == n_row_height_A); // todo - uncomment this if the data is available
					_ASSERTE(n_row_id_A < transpose_cols_list.size());
					std::vector<TColumn::TBlockEntry> &r_transpose_column = transpose_cols_list[n_row_id_A];
					if(!CMatrixMultiply_InnerLoop<_TyGppContextA, _TyGppContextB,
					   _CSelectedRowHeightsListA, CLastRowHeightB, CColumnWidthB>::Loop(
					   r_transpose_column, p_block_A, n_row_height_A, n_column_width_A,
					   n_column_id_B, p_block_B, alloc))
						return false;
					// call the inner loop
				}

				return true;
			}
			// todo - debug loop
		};

		/**
		 *	@brief (inner) loop of the matrix-matrix multiplication kernel template
		 *
		 *	@tparam _TyGppContextA is template instance context for g++ (compatibility workarround)
		 *	@tparam _TyGppContextB is template instance context for g++ (compatibility workarround)
		 *	@tparam _CRowHeightsListA is left-hand matrix row heights list, filtered by the selected right-hand matrix column width
		 *	@tparam CRowHeightB is the chosen right-hand matrix row height
		 *	@tparam CColumnWidthB is the chosen right-hand matrix column width
		 */
		template <class _TyGppContextA, class _TyGppContextB,
			class _CRowHeightsListA, class CRowHeightB, class CColumnWidthB>
		class CMatrixMultiply_InnerLoop {
		protected:
			typedef typename _CRowHeightsListA::_TyHead CCurrentRowHeightA; /**< @brief left-hand matrix row height for this decision tree recursion */

		public:
			/**
			 *	@brief wraps the inner loop of matrix-matrix multiplication in block size decision tree
			 *
			 *	@param[out] r_transpose_column is the transpose product matrix block column to store the new block
			 *	@param[in] p_block_A is pointer to dense data of the current block in the left-hand matrix
			 *	@param[in] n_row_height_A is height of the current block in the left-hand matrix
			 *	@param[in] n_column_width_A is width of the current block in the left-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column in the right-hand matrix
			 *	@param[in] p_block_B is pointer to dense data of the current block in the right-hand matrix
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop(std::vector<TColumn::TBlockEntry> &r_transpose_column,
				const double *p_block_A, size_t n_row_height_A, size_t n_column_width_A,
				size_t n_column_id_B, const double *p_block_B, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				if(n_row_height_A == CCurrentRowHeightA::n_size) {
					return CMatrixMultiply_InnerLoop<_TyGppContextA, _TyGppContextB,
						CTypelist<CCurrentRowHeightA, CTypelistEnd>, CRowHeightB,
						CColumnWidthB>::Loop(r_transpose_column, p_block_A, n_row_height_A,
						n_column_width_A, n_column_id_B, p_block_B, alloc);
				} else {
					return CMatrixMultiply_InnerLoop<_TyGppContextA, _TyGppContextB,
						typename _CRowHeightsListA::_TyTail, CRowHeightB, CColumnWidthB>::Loop(
						r_transpose_column, p_block_A, n_row_height_A, n_column_width_A,
						n_column_id_B, p_block_B, alloc);
				}
			}

			// todo - debug loop
		};

		/**
		 *	@brief (inner) loop of the matrix-matrix multiplication kernel template
		 *
		 *	@tparam _TyGppContextA is template instance context for g++ (compatibility workarround)
		 *	@tparam _TyGppContextB is template instance context for g++ (compatibility workarround)
		 *	@tparam CLastRowHeightA is the chosen left-hand matrix row height
		 *	@tparam CRowHeightB is the chosen right-hand matrix row height
		 *	@tparam CColumnWidthB is the chosen right-hand matrix column width
		 */
		template <class _TyGppContextA, class _TyGppContextB,
			class CLastRowHeightA, class CRowHeightB, class CColumnWidthB>
		class CMatrixMultiply_InnerLoop<_TyGppContextA, _TyGppContextB,
			CTypelist<CLastRowHeightA, CTypelistEnd>, CRowHeightB, CColumnWidthB> {
		protected:
			typedef CLastRowHeightA CRowHeightA; /**< @brief the chosen left-hand matrix row height (just to provide consistent names) */
			typedef CRowHeightB CColumnWidthA; /**< @brief left-hand matrix row height (supposed to equal the right-hand matrix row height if the matrix product is defined) */

			typedef typename CMakeMatrixRef<CRowHeightA::n_size, CColumnWidthA::n_size>::_Ty _TyBlockAShape; /**< @brief matrix block Eigen::Map type */
			typedef typename CMakeMatrixRef<CRowHeightB::n_size, CColumnWidthB::n_size>::_Ty _TyBlockBShape; /**< @brief matrix block Eigen::Map type */

			typedef typename CMakeMatrixRef<CRowHeightA::n_size, CColumnWidthB::n_size>::_Ty _TyDestBlockShape; /**< @brief matrix block Eigen::Map type */

		public:
			/**
			 *	@brief performs the inner loop of matrix-matrix multiplication
			 *
			 *	@param[out] r_transpose_column is the transpose product matrix block column to store the new block
			 *	@param[in] p_block_A is pointer to dense data of the current block in the left-hand matrix
			 *	@param[in] n_row_height_A is height of the current block in the left-hand matrix
			 *	@param[in] n_column_width_A is width of the current block in the left-hand matrix
			 *	@param[in] n_column_id_B is id of the current block column in the right-hand matrix
			 *	@param[in] p_block_B is pointer to dense data of the current block in the right-hand matrix
			 *	@param[in] alloc is the dense block allocator for the product matrix
			 *
			 *	@return Returns true on success, false on failure (incompatible block layout).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool Loop(std::vector<TColumn::TBlockEntry> &r_transpose_column,
				const double *p_block_A, size_t UNUSED(n_row_height_A), size_t n_column_width_A,
				size_t n_column_id_B, const double *p_block_B, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				_ASSERTE(n_row_height_A == CRowHeightA::n_size);

				// the B block is CRowHeightB x CColumnWidthB
				// the A block is CRowHeightA x CColumnWidthA

				// multiplication of blockA * blockB yields block at (n_row_id_A, n_column_id_B)

				/*_ASSERTE(n_row_id_A < r_row_list_dest.size());
				_ASSERTE(n_column_id_B < r_col_list_dest.size());
				//size_t n_prod_rows = r_row_list_dest[n_row_id_A].n_height;
				//size_t n_prod_cols = r_col_list_dest[n_column_id_B].n_width; // unused
				_ASSERTE(r_row_list_dest[n_row_id_A].n_height == blockA.rows());
				_ASSERTE(r_col_list_dest[n_column_id_B].n_width == blockB.cols());*/ // todo - not sure if we have all of them here, uncomment if possible
				if(n_column_width_A != CColumnWidthA::n_size)
					return false; // make sure the blocks are multiplicable (not able to verify by just merging the layout)
				// basic checks about matrix dimensions

				if(r_transpose_column.empty() || r_transpose_column.back().first < n_column_id_B) {
					double *p_new_block_data = alloc.p_Get_DenseStorage(CRowHeightA::n_size * CColumnWidthB::n_size);
					_TyDestBlockShape block_dest(p_new_block_data);
					// get storage

					_TyBlockAShape blockA((double*)p_block_A);
					_TyBlockBShape blockB((double*)p_block_B);
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
					block_dest = blockA.lazyProduct(blockB); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
					block_dest.noalias() = blockA * blockB;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
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
					_TyDestBlockShape block_dest(p_block_data);
					_TyBlockAShape blockA((double*)p_block_A);
					_TyBlockBShape blockB((double*)p_block_B);
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
					block_dest += blockA.lazyProduct(blockB); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
					block_dest.noalias() += blockA * blockB;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
					// add to the existing block (run the dot sum)
				}
				// perform the dense multiplication using reference matrices and eigen
				// do stuff / calculate matmult

				return true;
			}

			// todo - debug loop
		};

		/**
		 *	@brief performs the outer loop of matrix-matrix multiplication
		 *
		 *	@param[in] r_t_column_B is the current block column of the right-hand matrix
		 *	@param[in] n_column_id_B is id of the current block column of the right-hand matrix
		 *	@param[in] r_row_list_B is list of the right-hand matrix block rows
		 *	@param[in] r_col_list_A is list of the left-hand matrix block columns
		 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
		 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
		 *	@param[in] reindex_rows_B_to_cols_A is the mapping functions from right-hand
		 *		matrix block rows to left-hand matrix block columns
		 *	@param[in] alloc is the dense block allocator for the product matrix
		 *
		 *	@return Returns true on success, false on failure (incompatible block layout).
		 *
		 *	@note This function throws std::bad_alloc.
		 */
		static __forceinline bool MatrixMultiply_OuterLoop(const TColumn &r_t_column_B,
			size_t n_column_id_B, const std::vector<TRow> &r_row_list_B,
			const std::vector<TColumn> &r_col_list_A, const std::vector<TRow> &r_row_list_A,
			std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
			const std::vector<size_t> &reindex_rows_B_to_cols_A, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
		{
			return CMatrixMultiply_OuterLoop<_TyBlockMatrixTypelistA,
				_TyBlockMatrixTypelistB, CColumnWidthsListB>::Loop(r_t_column_B,
				n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
				transpose_cols_list, reindex_rows_B_to_cols_A, alloc);
		}

		/**
		 *	@brief performs the outer loop of matrix-matrix multiplication
		 *
		 *	@param[in] r_t_column_B is the current block column of the right-hand matrix
		 *	@param[in] n_column_id_B is id of the current block column of the right-hand matrix
		 *	@param[in] r_row_list_B is list of the right-hand matrix block rows
		 *	@param[in] r_col_list_A is list of the left-hand matrix block columns
		 *	@param[in] r_row_list_A is list of the left-hand matrix block rows
		 *	@param[out] transpose_cols_list is list of the transpose product matrix block columns
		 *	@param[in] reindex_rows_B_to_cols_A is the mapping functions from right-hand
		 *		matrix block rows to left-hand matrix block columns
		 *	@param[in] alloc is the dense block allocator for the product matrix
		 *
		 *	@return Returns true on success, false on failure (incompatible block layout).
		 *
		 *	@note This function throws std::bad_alloc.
		 *	@note This version produces only the upper-triangular part of the multiplication.
		 */
		static __forceinline bool MatrixMultiply_OuterLoop_UpperTriag(const TColumn &r_t_column_B,
			size_t n_column_id_B, const std::vector<TRow> &r_row_list_B,
			const std::vector<TColumn> &r_col_list_A, const std::vector<TRow> &r_row_list_A,
			std::vector<std::vector<TColumn::TBlockEntry> > &transpose_cols_list,
			const std::vector<size_t> &reindex_rows_B_to_cols_A, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
		{
			return CMatrixMultiply_OuterLoop<_TyBlockMatrixTypelistA,
				_TyBlockMatrixTypelistB, CColumnWidthsListB>::Loop_UpperTriag(r_t_column_B,
				n_column_id_B, r_row_list_B, r_col_list_A, r_row_list_A,
				transpose_cols_list, reindex_rows_B_to_cols_A, alloc);
		}
	};

	/**
	 *	@brief fixed block size operations template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 */
	template <class CBlockMatrixTypelist>
	class CFBS_ElementwiseUnaryOp : public fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist> {
	public:
		//typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CRowHeightsList CRowHeightsList; /**< @brief list of unique block row heights */
		//typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CColumnWidthsList CColumnWidthsList; /**< @brief list of unique block column widths */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CDimsList_Uniq CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */

		typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_Uniq,
			fbs_ut::CTransformDimensionToAreaSize>::_TyResult>::_TyResult CBlockAreasList; /**< @brief list of unique block areas (for elementwise ops) */
		// create typelists to generate the decision tree for block sizes

		// @todo - create blockwise operations (a variant to elementwise, using eigen map and aligned stuff)
		// t_odo - push forceinline and pragmas

		/**
		 *	@brief (inner) loop of the elementwise unary kernel template
		 *
		 *	Loops over block elements and applies unary operator to their values.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CBlockAreasList is list of posible numbers of elements of the blocks
		 *	@tparam CUnaryOperator is unary functor, or a pointer to function type
		 */
		template <class _TyGppContext, class _CBlockAreasList, class CUnaryOperator>
		class CElementwiseUnary_Loop {
		public:
			typedef typename _CBlockAreasList::_TyHead CCurrentBlockArea; /**< @brief block area for this decision tree recursion */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] p_block is current matrix block
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 *	@param[out] r_op is reference to an unary operator
			 *		object instance (or a function)
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(double *p_block,
				size_t n_block_area, CUnaryOperator &r_op)
			{
				if(n_block_area == CCurrentBlockArea::n_size) {
					CElementwiseUnary_Loop<_TyGppContext, CTypelist<CCurrentBlockArea,
						CTypelistEnd>, CUnaryOperator>::Loop(p_block, n_block_area, r_op);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CElementwiseUnary_Loop<_TyGppContext, typename _CBlockAreasList::_TyTail,
						CUnaryOperator>::Loop(p_block, n_block_area, r_op);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(n_block_area == %d) {\n", CCurrentBlockArea::n_size);
				CElementwiseUnary_Loop<_TyGppContext, CTypelist<CCurrentBlockArea,
					CTypelistEnd>, CUnaryOperator>::Debug();
				printf("} else {\n");
				CElementwiseUnary_Loop<_TyGppContext, typename _CBlockAreasList::_TyTail, CUnaryOperator>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief (inner) loop of the elementwise unary kernel template
		 *
		 *	Loops over block elements and applies unary operator to their values.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CLastBlockArea is number of elements of the block being processed
		 *	@tparam CUnaryOperator is unary functor, or a pointer to function type
		 */
		template <class _TyGppContext, class CLastBlockArea, class CUnaryOperator>
		class CElementwiseUnary_Loop<_TyGppContext, CTypelist<CLastBlockArea, CTypelistEnd>, CUnaryOperator> {
		public:
			/**
			 *	@brief loops over block elements and applies unary operator to their values
			 *
			 *	@param[in] p_block is current matrix block
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 *	@param[out] r_op is reference to an unary operator
			 *		object instance (or a function)
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(double *p_block,
				size_t UNUSED(n_block_area), CUnaryOperator &r_op)
			{
				_ASSERTE(n_block_area == CLastBlockArea::n_size);
				double *p_data = p_block; // @todo - at this point, the compiler does not know that the pointer is aligned, and will not use SSE (except in x64 maybe). see what assembly this generates and try using __declspace(align(16)) to make it (__attribute__((aligned(16)) in linux). // not sure that specifying alignment for a pointer will enable SSE. see if eigen does something like this in Map<T,A>
				for(size_t i = 0; i < CLastBlockArea::n_size; ++ i, ++ p_data)
					*p_data = r_op(*p_data);
				// process all the data (should be able to unroll the loop)
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(n_block_area == %d);\n", CLastBlockArea::n_size);
				printf("for(size_t i = 0; i < %d; ++ i)\n", CLastBlockArea::n_size);
				printf("p_block[i] = op(p_block[i]);\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief loops over block elements and applies unary operator to their values
		 *
		 *	@tparam CUnnaryOperator is unary operator object type (or a pointer to a function)
		 *
		 *	@param[in] p_block is pointer ot the current matrix block data
		 *	@param[in] n_block_area is number of elements in the current matrix block
		 *	@param[out] r_op is reference to an unary operator
		 *		object instance (or a function)
		 */
		template <class CUnaryOperator>
		static __forceinline void ElementwiseUnary_Loop(double *p_block,
			size_t n_block_area, CUnaryOperator &r_op)
		{
			CElementwiseUnary_Loop<CBlockMatrixTypelist, CBlockAreasList,
				CUnaryOperator>::Loop(p_block, n_block_area, r_op);
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unindented) pseudocode of the elementwise-op function to stdout
		 *	@tparam CUnaryOperator is unary operator object type (or a pointer to a function)
		 *	@note This function is only available if the
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING macro is defined.
		 */
		template <class CUnaryOperator>
		static void Debug_ElementwiseUnaryOp()
		{
			printf("for(each p_col_it in m_block_cols_list) {\n");
			printf("for(each p_block_it in (*p_col_it).block_list) {\n");
			CElementwiseUnary_Loop<CBlockMatrixTypelist, CBlockAreasList, CUnaryOperator>::Debug();
			printf("}\n");
			printf("}\n\n");
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief fixed block size operations template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 */
	template <class CBlockMatrixTypelist>
	class CFBS_ElementwiseBinaryOp : public fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist> {
	public:
		//typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CRowHeightsList CRowHeightsList; /**< @brief list of unique block row heights */
		//typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CColumnWidthsList CColumnWidthsList; /**< @brief list of unique block column widths */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CDimsList_Uniq CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */

		typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_Uniq,
			fbs_ut::CTransformDimensionToAreaSize>::_TyResult>::_TyResult CBlockAreasList; /**< @brief list of unique block areas (for elementwise ops) */
		// create typelists to generate the decision tree for block sizes

		/**
		 *	@brief (inner) loop of the elementwise binary kernel template
		 *
		 *	Loops over block elements and applies binary operator to their values.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CBlockAreasList is list of posible numbers of elements of the blocks
		 *	@tparam CBinaryOperator is binary functor, or a pointer to function type
		 */
		template <class _TyGppContext, class _CBlockAreasList, class CBinaryOperator>
		class CElementwiseBinary_Loop {
		public:
			typedef typename _CBlockAreasList::_TyHead CCurrentBlockArea; /**< @brief block area for this decision tree recursion */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] p_block_left is current matrix block (left side operand, destination)
			 *	@param[in] p_block_right is current matrix block (right side operand)
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 *	@param[out] r_op is reference to an unary operator
			 *		object instance (or a function)
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(double *p_block_left, const double *p_block_right,
				size_t n_block_area, CBinaryOperator &r_op)
			{
				if(n_block_area == CCurrentBlockArea::n_size) {
					CElementwiseBinary_Loop<_TyGppContext, CTypelist<CCurrentBlockArea,
						CTypelistEnd>, CBinaryOperator>::Loop(p_block_left,
						p_block_right, n_block_area, r_op);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CElementwiseBinary_Loop<_TyGppContext, typename _CBlockAreasList::_TyTail,
						CBinaryOperator>::Loop(p_block_left, p_block_right, n_block_area, r_op);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] p_block_left is current matrix block (left side operand, destination)
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 *	@param[out] r_op is reference to an unary operator
			 *		object instance (or a function)
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop_NullRightSide(double *p_block_left,
				size_t n_block_area, CBinaryOperator &r_op)
			{
				if(n_block_area == CCurrentBlockArea::n_size) {
					CElementwiseBinary_Loop<_TyGppContext, CTypelist<CCurrentBlockArea,
						CTypelistEnd>, CBinaryOperator>::Loop_NullRightSide(p_block_left,
						n_block_area, r_op);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CElementwiseBinary_Loop<_TyGppContext, typename _CBlockAreasList::_TyTail,
						CBinaryOperator>::Loop_NullRightSide(p_block_left, n_block_area, r_op);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] p_block_left is current matrix block (left side operand, destination)
			 *	@param[in] p_block_right is current matrix block (right side operand)
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 *	@param[out] r_op is reference to an unary operator
			 *		object instance (or a function)
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop_UninitLeftSide(double *p_block_left, const double *p_block_right,
				size_t n_block_area, CBinaryOperator &r_op)
			{
				if(n_block_area == CCurrentBlockArea::n_size) {
					CElementwiseBinary_Loop<_TyGppContext, CTypelist<CCurrentBlockArea,
						CTypelistEnd>, CBinaryOperator>::Loop_UninitLeftSide(p_block_left,
						p_block_right, n_block_area, r_op);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CElementwiseBinary_Loop<_TyGppContext, typename _CBlockAreasList::_TyTail,
						CBinaryOperator>::Loop_UninitLeftSide(p_block_left, p_block_right,
						n_block_area, r_op);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] p_block_left is current matrix block (left side operand, destination)
			 *	@param[in] p_block_right is current matrix block (right side operand)
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop_UninitializedCopy(double *p_block_left,
				const double *p_block_right, size_t n_block_area)
			{
				if(n_block_area == CCurrentBlockArea::n_size) {
					CElementwiseBinary_Loop<_TyGppContext, CTypelist<CCurrentBlockArea,
						CTypelistEnd>, CBinaryOperator>::Loop_UninitializedCopy(p_block_left,
						p_block_right, n_block_area);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CElementwiseBinary_Loop<_TyGppContext, typename _CBlockAreasList::_TyTail,
						CBinaryOperator>::Loop_UninitializedCopy(p_block_left,
						p_block_right, n_block_area);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(n_block_area == %d) {\n", CCurrentBlockArea::n_size);
				CElementwiseBinary_Loop<_TyGppContext, CTypelist<CCurrentBlockArea,
					CTypelistEnd>, CBinaryOperator>::Debug();
				printf("} else {\n");
				CElementwiseBinary_Loop<_TyGppContext, typename _CBlockAreasList::_TyTail,
					CBinaryOperator>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief (inner) loop of the elementwise binary kernel template
		 *
		 *	Loops over block elements and applies binary operator to their values.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CLastBlockArea is number of elements of the block being processed
		 *	@tparam CBinaryOperator is binary functor, or a pointer to function type
		 */
		template <class _TyGppContext, class CLastBlockArea, class CBinaryOperator>
		class CElementwiseBinary_Loop<_TyGppContext, CTypelist<CLastBlockArea, CTypelistEnd>, CBinaryOperator> {
		public:
			/**
			 *	@brief loops over block elements and applies binary operator to their values
			 *
			 *	@param[in] p_block_left is current matrix block (left side operand, destination)
			 *	@param[in] p_block_right is current matrix block (right side operand)
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 *	@param[out] r_op is reference to an unary operator
			 *		object instance (or a function)
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(double *p_block_left, const double *p_block_right,
				size_t UNUSED(n_block_area), CBinaryOperator &r_op)
			{
				_ASSERTE(n_block_area == CLastBlockArea::n_size);
				double *p_dest = p_block_left; // @todo - at this point, the compiler does not know that the pointer is aligned, and will not use SSE (except in x64 maybe). see what assembly this generates and try using __declspace(align(16)) to make it (__attribute__((aligned(16)) in linux). // not sure that specifying alignment for a pointer will enable SSE. see if eigen does something like this in Map<T,A>
				const double *p_right = p_block_right;
				for(size_t i = 0; i < CLastBlockArea::n_size; ++ i, ++ p_dest, ++ p_right)
					*p_dest = r_op(*p_dest, *p_right);
				// process all the data (should be able to unroll the loop)
			}

			/**
			 *	@brief loops over block elements and applies binary operator to their values
			 *
			 *	@param[in] p_block_left is current matrix block (left side operand, destination)
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 *	@param[out] r_op is reference to an unary operator
			 *		object instance (or a function)
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop_NullRightSide(double *p_block_left,
				size_t UNUSED(n_block_area), CBinaryOperator &r_op)
			{
				_ASSERTE(n_block_area == CLastBlockArea::n_size);
				double *p_dest = p_block_left; // @todo - at this point, the compiler does not know that the pointer is aligned, and will not use SSE (except in x64 maybe). see what assembly this generates and try using __declspace(align(16)) to make it (__attribute__((aligned(16)) in linux). // not sure that specifying alignment for a pointer will enable SSE. see if eigen does something like this in Map<T,A>
				for(size_t i = 0; i < CLastBlockArea::n_size; ++ i, ++ p_dest)
					*p_dest = r_op(*p_dest, 0); // right hand side is null
				// process all the data (should be able to unroll the loop)
			}

			/**
			 *	@brief loops over block elements and applies binary operator to their values
			 *
			 *	@param[in] p_block_left is current matrix block (left side operand, destination)
			 *	@param[in] p_block_right is current matrix block (right side operand)
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 *	@param[out] r_op is reference to an unary operator
			 *		object instance (or a function)
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop_UninitLeftSide(double *p_block_left, const double *p_block_right,
				size_t UNUSED(n_block_area), CBinaryOperator &r_op)
			{
				_ASSERTE(n_block_area == CLastBlockArea::n_size);
				double *p_dest = p_block_left; // @todo - at this point, the compiler does not know that the pointer is aligned, and will not use SSE (except in x64 maybe). see what assembly this generates and try using __declspace(align(16)) to make it (__attribute__((aligned(16)) in linux). // not sure that specifying alignment for a pointer will enable SSE. see if eigen does something like this in Map<T,A>
				const double *p_right = p_block_right;
				for(size_t i = 0; i < CLastBlockArea::n_size; ++ i, ++ p_dest, ++ p_right)
					*p_dest = r_op(0, *p_right); // left hand side is not initialized; write only (read zero)
				// process all the data (should be able to unroll the loop)
			}

			/**
			 *	@brief loops over block elements and applies binary operator to their values
			 *
			 *	@param[in] p_block_left is current matrix block (left side operand, destination)
			 *	@param[in] p_block_right is current matrix block (right side operand)
			 *	@param[in] n_block_area is number of elements in the current matrix block
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop_UninitializedCopy(double *p_block_left,
				const double *p_block_right, size_t UNUSED(n_block_area))
			{
				_ASSERTE(n_block_area == CLastBlockArea::n_size);
				double *p_dest = p_block_left; // @todo - at this point, the compiler does not know that the pointer is aligned, and will not use SSE (except in x64 maybe). see what assembly this generates and try using __declspace(align(16)) to make it (__attribute__((aligned(16)) in linux). // not sure that specifying alignment for a pointer will enable SSE. see if eigen does something like this in Map<T,A>
				const double *p_right = p_block_right;
				for(size_t i = 0; i < CLastBlockArea::n_size; ++ i, ++ p_dest, ++ p_right)
					*p_dest = *p_right;
				// process all the data (should be able to unroll the loop)
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(n_block_area == %d);\n", CLastBlockArea::n_size);
				printf("for(size_t i = 0; i < %d; ++ i)\n", CLastBlockArea::n_size);
				printf("p_block_left[i] = op(p_block_left[i], p_block_right[i]);\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief loops over block elements and applies binary operator to their values
		 *
		 *	@tparam CBinaryOperator is binary operator object type (or a pointer to a function)
		 *
		 *	@param[in] p_block_left is pointer to the current left-side matrix block data (also destination)
		 *	@param[in] p_block_right is pointer to the current right-side matrix block data
		 *	@param[in] n_block_area is number of elements in the current matrix block
		 *	@param[out] r_op is reference to an binary operator
		 *		object instance (or a function)
		 */
		template <class CBinaryOperator>
		static __forceinline void ElementwiseBinary_Loop(double *p_block_left,
			const double *p_block_right, size_t n_block_area, CBinaryOperator &r_op)
		{
			CElementwiseBinary_Loop<CBlockMatrixTypelist, CBlockAreasList,
				CBinaryOperator>::Loop(p_block_left, p_block_right, n_block_area, r_op);
		}

		/**
		 *	@brief loops over block elements and applies binary operator to their values
		 *
		 *	@tparam CBinaryOperator is binary operator object type (or a pointer to a function)
		 *
		 *	@param[in] p_block_left is pointer to the current left-side matrix block data (also destination)
		 *	@param[in] n_block_area is number of elements in the current matrix block
		 *	@param[out] r_op is reference to an binary operator
		 *		object instance (or a function)
		 */
		template <class CBinaryOperator>
		static __forceinline void ElementwiseBinary_Loop_NullRightSide(double *p_block_left,
			size_t n_block_area, CBinaryOperator &r_op)
		{
			CElementwiseBinary_Loop<CBlockMatrixTypelist, CBlockAreasList,
				CBinaryOperator>::Loop_NullRightSide(p_block_left, n_block_area, r_op);
		}

		/**
		 *	@brief loops over block elements and applies binary operator to their values
		 *
		 *	@tparam CBinaryOperator is binary operator object type (or a pointer to a function)
		 *
		 *	@param[in] p_block_left is pointer to the current left-side matrix block data (also destination)
		 *	@param[in] p_block_right is pointer to the current right-side matrix block data
		 *	@param[in] n_block_area is number of elements in the current matrix block
		 *	@param[out] r_op is reference to an binary operator
		 *		object instance (or a function)
		 */
		template <class CBinaryOperator>
		static __forceinline void ElementwiseBinary_Loop_UninitLeftSide(double *p_block_left,
			const double *p_block_right, size_t n_block_area, CBinaryOperator &r_op)
		{
			CElementwiseBinary_Loop<CBlockMatrixTypelist, CBlockAreasList,
				CBinaryOperator>::Loop_UninitLeftSide(p_block_left, p_block_right, n_block_area, r_op);
		}

		/**
		 *	@brief copies values from the right side block to the left side (no operation involved)
		 *
		 *	@tparam CBinaryOperator is binary operator object type (or a pointer to a function)
		 *
		 *	@param[in] p_block_left is pointer to the current left-side matrix block data (also destination)
		 *	@param[in] p_block_right is pointer to the current right-side matrix block data
		 *	@param[in] n_block_area is number of elements in the current matrix block
		 *	@param[out] r_op is reference to an binary operator
		 *		object instance (or a function)
		 */
		template <class CBinaryOperator>
		static __forceinline void ElementwiseBinary_Loop_UninitializedCopy(double *p_block_left,
			const double *p_block_right, size_t n_block_area, CBinaryOperator &UNUSED(r_op))
		{
			CElementwiseBinary_Loop<CBlockMatrixTypelist, CBlockAreasList,
				CBinaryOperator>::Loop_UninitializedCopy(p_block_left, p_block_right, n_block_area);
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unindented) pseudocode of the elementwise-op function to stdout
		 *	@tparam CBinaryOperator is binary operator object type (or a pointer to a function)
		 *	@note This function is only available if the
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING macro is defined.
		 */
		template <class CBinaryOperator>
		static void Debug_ElementwiseBinaryOp()
		{
			printf("for(each p_col_it in m_block_cols_list) {\n");
			printf("for(each p_block_it in (*p_col_it).block_list) {\n");
			printf("TColumn::TBlockEntry &r_t_block_left = *p_block_it;\n");
			printf("const TColumn::TBlockEntry &r_t_block_right =\n"
				"r_t_other_matrix.t_FindBlock(*p_col_it, *p_block_it);\n"); // very pseudocode
			CElementwiseBinary_Loop<CBlockMatrixTypelist, CBlockAreasList, CBinaryOperator>::Debug();
			printf("}\n");
			printf("}\n\n");
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief fixed block size operations template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 */
	template <class CBlockMatrixTypelist>
	class CFBS_PostMAD : public fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist> {
	public:
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CRowHeightsList CRowHeightsList; /**< @brief list of unique block row heights */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CColumnWidthsList CColumnWidthsList; /**< @brief list of unique block column widths */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CDimsList_Uniq CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */

		/**
		 *	@brief inner loop of the post multiply add kernel
		 *
		 *	Performs a single block-vector multiplication and accumulation.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CRowHeightsList is list of possible row heights
		 *	@tparam CColumnWidth is current column width
		 */
		template <class _TyGppContext, class _CRowHeightsList, class CColumnWidth>
		class CPostMAD_InnerLoop {
		public:
			typedef typename _CRowHeightsList::_TyHead CCurrentRowHeight; /**< @brief row height for this decision tree recursion */
			typedef typename CMakeRowVectorRef<CColumnWidth::n_size>::_Ty _TyDestVectorShape; /**< @brief destination vector Eigen::Map type (note it is unaligned) */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_block is current matrix block
			 *	@param[in] r_t_row is row for the current matrix block
			 *	@param[out] r_dest is destination vector map (only the part
			 *		that aligns with the current matrix block)
			 *	@param[in] p_src_vector is pointer to the (whole) source vector
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_block, const TRow &r_t_row,
				_TyDestVectorShape &r_dest, const double *p_src_vector)
			{
				if(r_t_row.n_height == CCurrentRowHeight::n_size) {
					CPostMAD_InnerLoop<_TyGppContext, CTypelist<CCurrentRowHeight,
						CTypelistEnd>, CColumnWidth>::Loop(r_t_block,
						r_t_row, r_dest, p_src_vector);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPostMAD_InnerLoop<_TyGppContext, typename _CRowHeightsList::_TyTail, CColumnWidth>::Loop(
						r_t_block, r_t_row, r_dest, p_src_vector);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(r_t_row.n_height == %d) {\n", CCurrentRowHeight::n_size);
				CPostMAD_InnerLoop<_TyGppContext, CTypelist<CCurrentRowHeight,
					CTypelistEnd>, CColumnWidth>::Debug();
				printf("} else {\n");
				CPostMAD_InnerLoop<_TyGppContext, typename _CRowHeightsList::_TyTail, CColumnWidth>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief inner loop of the post multiply add kernel
		 *		(specialization for the end-of-the-list (branch-less))
		 *
		 *	Performs a single block-vector multiplication and accumulation.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CLastRowHeight is height of the row at the end of the list (current)
		 *	@tparam CColumnWidth is current column width
		 */
		template <class _TyGppContext, class CLastRowHeight, class CColumnWidth>
		class CPostMAD_InnerLoop<_TyGppContext, CTypelist<CLastRowHeight, CTypelistEnd>, CColumnWidth> {
		public:
			typedef typename CMakeRowVectorRef<CColumnWidth::n_size>::_Ty _TyDestVectorShape; /**< @brief destination vector Eigen::Map type (note it is unaligned) */
			typedef typename CMakeRowVectorRef<CLastRowHeight::n_size>::_Ty _TySrcVectorShape; /**< @brief source vector Eigen::Map type (note it is unaligned) */
			typedef typename CMakeMatrixRef<CLastRowHeight::n_size, CColumnWidth::n_size>::_Ty _TyBlockShape; /**< @brief matrix block Eigen::Map type */

			/**
			 *	@brief performs a single block-vector multiplication and accumulation
			 *
			 *	@param[in] r_t_block is current matrix block
			 *	@param[in] r_t_row is row for the current matrix block
			 *	@param[out] r_dest is destination vector map (only the part
			 *		that aligns with the current matrix block)
			 *	@param[in] p_src_vector is pointer to the (whole) source vector
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_block, const TRow &r_t_row,
				_TyDestVectorShape &r_dest, const double *p_src_vector)
			{
				_ASSERTE(r_t_row.n_height == CLastRowHeight::n_size); // this is the last option, make sure it is this one (omits the branch)
				const size_t n_row = r_t_row.n_cumulative_height_sum; // src vector offset
				// for-each row

				_TySrcVectorShape src((double*)p_src_vector + n_row);
				_TyBlockShape block((double*)r_t_block.second);

#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				r_dest += src.lazyProduct(block); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				r_dest.noalias() += src * block; // xapy
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				// perform the multiplication, one block at a time
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(r_t_row.n_heigh == %d);\n", CLastRowHeight::n_size);
				printf("typedef Eigen::Map<Eigen::Matrix<double, %d, %d> > _TySrcVectorShape;\n",
					1, CLastRowHeight::n_size);
				printf("typedef Eigen::Map<Eigen::Matrix<double, %d, %d> > _TyBlockShape;\n",
					CColumnWidth::n_size, CLastRowHeight::n_size);
				printf("r_dest += src * block; // xapy\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief outer loop of the post multiply add kernel
		 *
		 *	Iterates through all the blocks of the column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthsList is list of possible column widths
		 */
		template <class _TyGppContext, class _CColumnWidthsList>
		class CPostMAD_OuterLoop {
		public:
			typedef typename _CColumnWidthsList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_col is current column
			 *	@param[out] p_dest_vector is pointer to the (whole) destination vector
			 *	@param[in] p_src_vector is pointer to the (whole) source vector
			 *	@param[in] r_row_list is matrix row layout
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn &r_t_col, double *p_dest_vector,
				const double *p_src_vector, const std::vector<TRow> &r_row_list)
			{
				if(r_t_col.n_width == CCurrentColumnWidth::n_size) {
					CPostMAD_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth, CTypelistEnd> >::Loop(
						r_t_col, p_dest_vector, p_src_vector, r_row_list);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPostMAD_OuterLoop<_TyGppContext, typename _CColumnWidthsList::_TyTail>::Loop(
						r_t_col, p_dest_vector, p_src_vector, r_row_list);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(r_t_col.n_width == %d) {\n", CCurrentColumnWidth::n_size);
				CPostMAD_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth, CTypelistEnd> >::Debug();
				printf("} else {\n");
				CPostMAD_OuterLoop<_TyGppContext, typename _CColumnWidthsList::_TyTail>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief outer loop of the post multiply add kernel
		 *		(specialization for the end-of-the-list (branch-less))
		 *
		 *	Iterates through all the blocks of the column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CLastColumnWidth is width of the column at the end of the list (current)
		 */
		template <class _TyGppContext, class CLastColumnWidth>
		class CPostMAD_OuterLoop<_TyGppContext, CTypelist<CLastColumnWidth, CTypelistEnd> > {
		public:
			typedef typename CMakeRowVectorRef<CLastColumnWidth::n_size>::_Ty _TyDestVectorShape; /**< @brief destination vector Eigen::Map type (note it is unaligned) */
			typedef typename CFilterTypelist2<CRowHeightsList, CLastColumnWidth,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_Uniq>::_TyResult _CSelectedRowHeightsList; /**< @brief row heights list, filtered by the selected column width */

			/**
			 *	@brief iterates through all the blocks in the current
			 *		column and calls the inner loop
			 *
			 *	@param[in] r_t_col is current column
			 *	@param[out] p_dest_vector is pointer to the (whole) destination vector
			 *	@param[in] p_src_vector is pointer to the (whole) source vector
			 *	@param[in] r_row_list is matrix row layout
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn &r_t_col, double *p_dest_vector,
				const double *p_src_vector, const std::vector<TRow> &r_row_list)
			{
				_ASSERTE(r_t_col.block_list.empty() || r_t_col.n_width == CLastColumnWidth::n_size); // this is the last option, make sure it is this one (omits the branch)
				const size_t n_column = r_t_col.n_cumulative_width_sum; // dest vector offset
				// for-each column

				_TyDestVectorShape dest(p_dest_vector + n_column); // forces unaligned memory - no SSE2 (p_dest_vector can't be aligned due to block size)
				// gets dest vector as eigen blah

				for(_TyBlockConstIter p_block_it =
				   r_t_col.block_list.begin(), p_end_it = r_t_col.block_list.end();
				   p_block_it != p_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					// for each column, for each block ...

					const TRow &r_t_row = r_row_list[r_t_block.first];
					// get block row ...

					CPostMAD_InnerLoop<_TyGppContext, _CSelectedRowHeightsList,
						CLastColumnWidth>::Loop(r_t_block,
						r_t_row, /*r_t_col,*/ dest, p_src_vector);
					// wrap the inner loop in a block size decision tree
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(r_t_col.n_width == %d);\n", CLastColumnWidth::n_size);
				printf("typedef Eigen::Map<Eigen::Matrix<double, %d, %d> > _TyDestVectorShape;\n",
					1, CLastColumnWidth::n_size);
				printf("for(each p_block_it in r_t_col.block_list) {\n");
				CPostMAD_InnerLoop<_TyGppContext, _CSelectedRowHeightsList, CLastColumnWidth>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief iterates through all the blocks in the current
		 *		column and calls the inner loop
		 *
		 *	@param[in] r_t_col is current column
		 *	@param[out] p_dest_vector is pointer to the (whole) destination vector
		 *	@param[in] p_src_vector is pointer to the (whole) source vector
		 *	@param[in] r_row_list is matrix row layout
		 */
		static __forceinline void PostMAD_OuterLoop(const TColumn &r_t_col,
			double *p_dest_vector, const double *p_src_vector,
			const std::vector<TRow> &r_row_list)
		{
			CPostMAD_OuterLoop<CBlockMatrixTypelist, CColumnWidthsList>::Loop(r_t_col,
				p_dest_vector, p_src_vector, r_row_list);
			// wrap the outer loop body in a column size decision tree
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unindented) pseudocode of the post-multiply-add function to stdout
		 *	@note This function is only available if the
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING macro is defined.
		 */
		static void Debug_PostMAD()
		{
			printf("for(each p_col_it in m_block_cols_list) {\n");
			CPostMAD_OuterLoop<CBlockMatrixTypelist, CColumnWidthsList>::Debug();
			printf("}\n\n");
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief fixed block size operations template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 */
	template <class CBlockMatrixTypelist>
	class CFBS_PreMAD {
	public:
		/**
		 *	@brief inner loop of the post multiply add kernel
		 *	Performs a single block-vector multiplication and accumulation.
		 *	@tparam CColumnWidth is current column width
		 */
		template <class CColumnWidth>
		class CPreMAD_InnerLoop {
		protected:
			const double *m_p_block_data;
			double *m_p_dest_vector;
			const double *m_p_src_vector;

		public:
			CPreMAD_InnerLoop(const double *p_block_data, double *p_dest_vector, const double *p_src_vector)
				:m_p_block_data(p_block_data), m_p_dest_vector(p_dest_vector), m_p_src_vector(p_src_vector)
			{}

			template <class CRowHeight>
			__forceinline void operator ()()
			{
				typedef typename CMakeVectorRef<CColumnWidth::n_size>::_TyConst _TySrcVectorShape; // source vector Eigen::Map type (note it is unaligned)
				typedef typename CMakeVectorRef<CRowHeight::n_size>::_Ty _TyDestVectorShape; // destination vector Eigen::Map type (note it is unaligned)
				typedef typename CMakeMatrixRef<CRowHeight::n_size, CColumnWidth::n_size>::_TyConst _TyBlockShape; // matrix block Eigen::Map type

				_TySrcVectorShape src(m_p_src_vector);
				_TyDestVectorShape dest(m_p_dest_vector);
				_TyBlockShape block(m_p_block_data);

#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				dest += block.lazyProduct(src); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				dest.noalias() += block * src; // axpy
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				// perform the multiplication, one block at a time
			}
		};

		/**
		 *	@brief outer loop of the post multiply add kernel
		 *	Iterates through all the blocks of the column.
		 */
		class CPreMAD_OuterLoop {
		protected:
			const TColumn &m_r_t_col;
			double *m_p_dest_vector;
			const double *m_p_src_vector;
			const std::vector<TRow> &m_r_row_list;

		public:
			CPreMAD_OuterLoop(const TColumn &r_t_col, double *p_dest_vector,
				const double *p_src_vector, const std::vector<TRow> &r_row_list)
				:m_r_t_col(r_t_col), m_p_dest_vector(p_dest_vector),
				m_p_src_vector(p_src_vector), m_r_row_list(r_row_list)
			{}

			template <class CColWidth>
			__forceinline void operator ()()
			{
				_ASSERTE(m_r_t_col.block_list.empty() || m_r_t_col.n_width == CColWidth::n_size); // this is the last option, make sure it is this one (omits the branch)
				const size_t n_column = m_r_t_col.n_cumulative_width_sum; // src vector offset

				const double *p_src = m_p_src_vector + n_column;
				for(_TyBlockConstIter p_block_it =
				   m_r_t_col.block_list.begin(), p_end_it = m_r_t_col.block_list.end();
				   p_block_it != p_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					// for each column, for each block ...

					const TRow &r_t_row = m_r_row_list[r_t_block.first];
					// get block row ...

					const size_t n_row = r_t_row.n_cumulative_height_sum; // dest vector offset
				
					fbs_ut::CWrap3<>::In_RowHeight_DecisionTree_Given_ColumnWidth<CBlockMatrixTypelist,
						CColWidth::n_size>(int(r_t_row.n_height), CPreMAD_InnerLoop<CColWidth>(r_t_block.second,
						m_p_dest_vector + n_row, p_src));
					// wrap the inner loop in a row height decision tree
				}
			}
		};

		/**
		 *	@brief iterates through all the blocks in the current
		 *		column and calls the inner loop
		 *
		 *	@param[in] r_t_col is current column
		 *	@param[out] p_dest_vector is pointer to the (whole) destination vector
		 *	@param[in] p_src_vector is pointer to the (whole) source vector
		 *	@param[in] r_row_list is matrix row layout
		 */
		static __forceinline void PreMAD_OuterLoop(const TColumn &r_t_col,
			double *p_dest_vector, const double *p_src_vector,
			const std::vector<TRow> &r_row_list)
		{
			fbs_ut::CWrap3<fbs_ut::CCallFnOp>::In_ColumnWidth_DecisionTree<CBlockMatrixTypelist>(int(r_t_col.n_width),
				CPreMAD_OuterLoop(r_t_col, p_dest_vector, p_src_vector, r_row_list), !r_t_col.block_list.empty()); // only mind mismatches if not empty
			// use a decision tree
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unindented) pseudocode of the pre-multiply-add function to stdout
		 *	@note This function is only available if the
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING macro is defined.
		 */
		static void Debug_PreMAD()
		{}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief fixed block size operations template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 */
	template <class CBlockMatrixTypelist>
	class CFBS_PreATA : public fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist> {
	public:
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CRowHeightsList CRowHeightsList; /**< @brief list of unique block row heights */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CColumnWidthsList CColumnWidthsList; /**< @brief list of unique block column widths */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CDimsList_Uniq CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */

		/**
		 *	@brief inner loop of the part part of the AtA kernel
		 *		that calculates the upper triangle
		 *
		 *	Performs a single block multiplication.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CSelectedTransposeSecondRowHeightList is list of possible
		 *		row heights in the transposed matrix
		 *	@tparam CTransposeFirstRowHeight is row height of the first row
		 *		(the one selected in the middle loop) in the transpoed matrix
		 *	@tparam CTransposeColumnWidth is column width in the transpoed matrix
		 *		(selected in the outer loop)
		 */
		template <class _TyGppContext, class _CSelectedTransposeSecondRowHeightList,
			class CTransposeFirstRowHeight, class CTransposeColumnWidth>
		class CPreATA_UpperTri_InnerLoop {
		public:
			typedef typename _CSelectedTransposeSecondRowHeightList::_TyHead
				CCurrentTransposeSecondRowHeight; /**< @brief transpose row height for this decision tree recursion */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_first_block is the first block in the current column
			 *	@param[in] r_t_second_block is the second block in the current column
			 *		(guaranteed to be on a row below the first block row)
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_first_block,
				const TColumn::TBlockEntry &r_t_second_block,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, _TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				if(r_transpose_row_list[r_t_second_block.first].n_width ==
				   CCurrentTransposeSecondRowHeight::n_size) {
					CPreATA_UpperTri_InnerLoop<_TyGppContext, CTypelist<CCurrentTransposeSecondRowHeight,
						CTypelistEnd>, CTransposeFirstRowHeight,
						CTransposeColumnWidth>::Loop(r_t_first_block,
						r_t_second_block, r_transpose_row_list, r_col_list_dest, r_allocator);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPreATA_UpperTri_InnerLoop<_TyGppContext,
						typename _CSelectedTransposeSecondRowHeightList::_TyTail,
						CTransposeFirstRowHeight, CTransposeColumnWidth>::Loop(r_t_first_block,
						r_t_second_block, r_transpose_row_list, r_col_list_dest, r_allocator);
					// not CCurrentTransposeSecondRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#ifdef _OPENMP

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_first_block is the first block in the current column
			 *	@param[in] r_t_second_block is the second block in the current column
			 *		(guaranteed to be on a row below the first block row)
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_first_block,
				const TColumn::TBlockEntry &r_t_second_block,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, omp_lock_t &r_t_mutex,
				_TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				if(r_transpose_row_list[r_t_second_block.first].n_width ==
				   CCurrentTransposeSecondRowHeight::n_size) {
					CPreATA_UpperTri_InnerLoop<_TyGppContext, CTypelist<CCurrentTransposeSecondRowHeight,
						CTypelistEnd>, CTransposeFirstRowHeight,
						CTransposeColumnWidth>::Loop(r_t_first_block,
						r_t_second_block, r_transpose_row_list, r_col_list_dest,
						r_t_mutex, r_allocator);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPreATA_UpperTri_InnerLoop<_TyGppContext,
						typename _CSelectedTransposeSecondRowHeightList::_TyTail,
						CTransposeFirstRowHeight, CTransposeColumnWidth>::Loop(r_t_first_block,
						r_t_second_block, r_transpose_row_list, r_col_list_dest,
						r_t_mutex, r_allocator);
					// not CCurrentTransposeSecondRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#endif // _OPENMP

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(r_transpose_row_list[r_t_second_block.first].n_width == %d) {\n",
					CCurrentTransposeSecondRowHeight::n_size);
				CPreATA_UpperTri_InnerLoop<_TyGppContext, CTypelist<CCurrentTransposeSecondRowHeight,
					CTypelistEnd>, CTransposeFirstRowHeight, CTransposeColumnWidth>::Debug();
				printf("} else {\n");
				CPreATA_UpperTri_InnerLoop<_TyGppContext,
					typename _CSelectedTransposeSecondRowHeightList::_TyTail,
					CTransposeFirstRowHeight, CTransposeColumnWidth>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief inner loop of the part part of the AtA kernel that calculates
		 *		the upper triangle (specialization for the end-of-the-list (branch-less))
		 *
		 *	Performs a single block multiplication.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CTransposeSecondRowHeight is height of the row (in the transposed matrix)
		 *		at the end of the list (current)
		 *	@tparam CTransposeFirstRowHeight is row height of the first row
		 *		(the one selected in the middle loop) in the transpoed matrix
		 *	@tparam CTransposeColumnWidth is column width in the transpoed matrix
		 *		(selected in the outer loop)
		 */
		template <class _TyGppContext, class CTransposeSecondRowHeight,
			class CTransposeFirstRowHeight, class CTransposeColumnWidth>
		class CPreATA_UpperTri_InnerLoop<_TyGppContext, CTypelist<CTransposeSecondRowHeight, CTypelistEnd>,
			CTransposeFirstRowHeight, CTransposeColumnWidth> {
		public:
			typedef typename CMakeMatrixRef<CTransposeColumnWidth::n_size,
				CTransposeFirstRowHeight::n_size>::_Ty _TyFirstBlock_NotTransposed; /**< @brief first block matrix Eigen::Map type (before transpose!) */
			typedef typename CMakeMatrixRef<CTransposeColumnWidth::n_size,
				CTransposeSecondRowHeight::n_size>::_Ty _TySecondBlock; /**< @brief second block matrix Eigen::Map type */
			typedef typename CMakeMatrixRef<CTransposeFirstRowHeight::n_size,
				CTransposeSecondRowHeight::n_size>::_TyMatrix _TyDestShape; /**< @brief destination block matrix type */
			typedef typename CMakeMatrixRef<CTransposeFirstRowHeight::n_size,
				CTransposeSecondRowHeight::n_size>::_Ty _TyDestBlock; /**< @brief destination block matrix Eigen::Map type */
			// typedef all block reference matrices (note rows x cols, it's all transposed since
			// we're looking at transpose this as an original, but haven't actually transposed it)

			/**
			 *	@brief performs a single block multiplication
			 *
			 *	@param[in] r_t_first_block is the first block in the current column
			 *	@param[in] r_t_second_block is the second block in the current column
			 *		(guaranteed to be on a row below the first block row)
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_first_block,
				const TColumn::TBlockEntry &r_t_second_block,
				const std::vector<TColumn> &UNUSED(r_transpose_row_list),
				std::vector<TColumn> &r_col_list_dest, _TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				const TColumn::TBlockEntry &r_t_block_A0 = r_t_first_block;
				size_t n_row_id_A0 = r_t_block_A0.first; // get row of the block
				// for each block in the current column in A

				const TColumn::TBlockEntry &r_t_block_A1 = r_t_second_block;
				size_t n_row_id_A1 = r_t_block_A1.first; // get row of the block
				// for each next block below the current block in the current column in A

				_ASSERTE(r_transpose_row_list[n_row_id_A0].n_width ==
					CTransposeFirstRowHeight::n_size); // was checked before
				_ASSERTE(r_transpose_row_list[n_row_id_A1].n_width ==
					CTransposeSecondRowHeight::n_size);
				// makes sure that the sizes do match

				//const size_t n_bm_cols = CTransposeColumnWidth::n_size;
				// column size alias

				const size_t n_bmA0_rows = CTransposeFirstRowHeight::n_size; // this is of transposed block
				_TyFirstBlock_NotTransposed blockA0_not_transposed(r_t_block_A0.second); // t_odo - this is mixing rows and cols, transpose and straight, definitely buggy // i'd say it works
				// create map to block A0 data
				// this block is oriented as in original (B) matrix, the data is not transposed either

				const size_t n_bmA1_cols = CTransposeSecondRowHeight::n_size; // this is of straight block
				_TySecondBlock blockA1(r_t_block_A1.second); // t_odo - this is mixing rows and cols, transpose and straight, definitely buggy // i'd say it works
				// create map to block A1 data
				// this block is oriented as in original (B) matrix, the data is not transposed either

#ifdef _DEBUG
				const std::vector<TColumn> &r_row_list_dest = r_transpose_row_list; // these are the same, just represented as TColumn, not TRow
				_ASSERTE(n_row_id_A0 < r_row_list_dest.size());
				_ASSERTE(n_row_id_A1 < r_col_list_dest.size());
				//size_t n_prod_rows = r_row_list_dest[n_row_id_A0].n_width;
				//size_t n_prod_cols = r_col_list_dest[n_row_id_A1].n_width; // unused
				_ASSERTE(r_row_list_dest[n_row_id_A0].n_width == blockA0_not_transposed.cols());
				_ASSERTE(r_col_list_dest[n_row_id_A1].n_width == blockA1.cols());
				_ASSERTE(blockA0_not_transposed.rows() == blockA1.rows()); // make sure the blocks are multiplicable (not able to verify by just merging the layout, but should be since we're multiplying A^T * A)
#endif // _DEBUG
				// basic checks about matrix dimensions

				{
					TColumn &r_t_dest_col = r_col_list_dest[n_row_id_A1];
					_TyBlockIter p_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						std::lower_bound(r_t_dest_col.block_list.begin(),
						r_t_dest_col.block_list.end(), n_row_id_A0,
						CUberBlockMatrix_Base::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(r_t_dest_col.block_list.begin(),
						r_t_dest_col.block_list.end(), n_row_id_A0, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					// find where to put the block in the column

					if((p_block_it == r_t_dest_col.block_list.end()) ||
					   (*p_block_it).first != n_row_id_A0) {
						// a new block

						double *p_value =
							r_allocator.p_Get_DenseStorage(n_bmA0_rows * n_bmA1_cols);
						_TyDestBlock block_dest(p_value);
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						block_dest = blockA0_not_transposed.transpose().lazyProduct(blockA1); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						block_dest.noalias() = blockA0_not_transposed.transpose() * blockA1;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						// create a new block, initialize with values

						r_t_dest_col.block_list.insert(p_block_it,
							TColumn::TBlockEntry(n_row_id_A0, p_value));
						// add it to the list
					} else {
						double *p_value_dest = (*p_block_it).second;
						// get existing block data

						_TyDestBlock block_dest(p_value_dest);
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						block_dest += blockA0_not_transposed.transpose().lazyProduct(blockA1); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						block_dest.noalias() += blockA0_not_transposed.transpose() * blockA1;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						// add values to an existing block
					}
				}
				// column is directly indexable, just need to sort blocks inside it
			}

#ifdef _OPENMP

			/**
			 *	@brief performs a single block multiplication
			 *
			 *	@param[in] r_t_first_block is the first block in the current column
			 *	@param[in] r_t_second_block is the second block in the current column
			 *		(guaranteed to be on a row below the first block row)
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_first_block,
				const TColumn::TBlockEntry &r_t_second_block,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, omp_lock_t &r_t_mutex,
				_TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				const TColumn::TBlockEntry &r_t_block_A0 = r_t_first_block;
				size_t n_row_id_A0 = r_t_block_A0.first; // get row of the block
				// for each block in the current column in A

				const TColumn::TBlockEntry &r_t_block_A1 = r_t_second_block;
				size_t n_row_id_A1 = r_t_block_A1.first; // get row of the block
				// for each next block below the current block in the current column in A

				_ASSERTE(r_transpose_row_list[n_row_id_A0].n_width ==
					CTransposeFirstRowHeight::n_size); // was checked before
				_ASSERTE(r_transpose_row_list[n_row_id_A1].n_width ==
					CTransposeSecondRowHeight::n_size);
				// makes sure that the sizes do match

				//const size_t n_bm_cols = CTransposeColumnWidth::n_size;
				// column size alias

				const size_t n_bmA0_rows = CTransposeFirstRowHeight::n_size; // this is of transposed block
				_TyFirstBlock_NotTransposed blockA0_not_transposed(r_t_block_A0.second); // t_odo - this is mixing rows and cols, transpose and straight, definitely buggy // i'd say it works
				// create map to block A0 data
				// this block is oriented as in original (B) matrix, the data is not transposed either

				const size_t n_bmA1_cols = CTransposeSecondRowHeight::n_size; // this is of straight block
				_TySecondBlock blockA1(r_t_block_A1.second); // t_odo - this is mixing rows and cols, transpose and straight, definitely buggy // i'd say it works
				// create map to block A1 data
				// this block is oriented as in original (B) matrix, the data is not transposed either

				const std::vector<TColumn> &r_row_list_dest = r_transpose_row_list; // these are the same, just represented as TColumn, not TRow
				_ASSERTE(n_row_id_A0 < r_row_list_dest.size());
				_ASSERTE(n_row_id_A1 < r_col_list_dest.size());
				//size_t n_prod_rows = r_row_list_dest[n_row_id_A0].n_width;
				//size_t n_prod_cols = r_col_list_dest[n_row_id_A1].n_width; // unused
				_ASSERTE(r_row_list_dest[n_row_id_A0].n_width == blockA0_not_transposed.cols());
				_ASSERTE(r_col_list_dest[n_row_id_A1].n_width == blockA1.cols());
				_ASSERTE(blockA0_not_transposed.rows() == blockA1.rows()); // make sure the blocks are multiplicable (not able to verify by just merging the layout, but should be since we're multiplying A^T * A)
				// basic checks about matrix dimensions

				{
					TColumn &r_t_dest_col = r_col_list_dest[n_row_id_A1];

					_TyBlockIter p_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
						std::lower_bound(r_t_dest_col.block_list.begin(),
						r_t_dest_col.block_list.end(), n_row_id_A0,
						CUberBlockMatrix_Base::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						std::lower_bound(r_t_dest_col.block_list.begin(),
						r_t_dest_col.block_list.end(), n_row_id_A0, CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					// find where to put the block in the column

					if((p_block_it == r_t_dest_col.block_list.end()) ||
					   (*p_block_it).first != n_row_id_A0) {
						// a new block

						double *p_value;
						omp_set_lock(&r_t_mutex); // note this could be implemented using OpenMP critical section
						p_value = r_allocator.p_Get_DenseStorage(n_bmA0_rows * n_bmA1_cols);
						omp_unset_lock(&r_t_mutex);
						// allocate, update block with the value (unsafe?)

						r_t_dest_col.block_list.insert(p_block_it,
							TColumn::TBlockEntry(n_row_id_A0, p_value));
						// add it to the list

						_TyDestBlock block_dest(p_value);
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						block_dest = blockA0_not_transposed.transpose().lazyProduct(blockA1); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						block_dest.noalias() = blockA0_not_transposed.transpose() * blockA1;//t_result;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						// create a new block, initialize with values

					} else {
						double *p_value_dest = (*p_block_it).second;
						// get existing block data

						_TyDestBlock block_dest(p_value_dest);
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						block_dest += blockA0_not_transposed.transpose().lazyProduct(blockA1); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						block_dest.noalias() += blockA0_not_transposed.transpose() * blockA1;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
						// add values to an existing block
					}
				}
				// column is directly indexable, just need to sort blocks inside it
			}

#endif // _OPENMP

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(r_transpose_row_list[r_t_first_block.first].n_width"
					" == %d); // was checked before\n", CTransposeFirstRowHeight::n_size);
				printf("_ASSERTE(r_transpose_row_list[r_t_second_block.first].n_width"
					" == %d);\n", CTransposeSecondRowHeight::n_size);
				printf("const size_t n_bm_cols = %d;\n", CTransposeColumnWidth::n_size);
				printf("block_dest += blockA0_not_transposed.transpose() * blockA1;\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief middle loop of the part part of the AtA kernel that calculates the upper triangle
		 *
		 *	Iterates through all the second blocks in the current column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CSelectedTransposeFirstRowHeightList is list of possible row heights in the transposed matrix
		 *		(the recursive class specializations peel this one off, can't be passed to inner loop)
		 *	@tparam _CSelectedTransposeSecondRowHeightList is list of possible row heights in the transposed matrix
		 *	@tparam CTransposeColumnWidth is column width in the transpoed matrix (selected in the outer loop)
		 */
		template <class _TyGppContext, class _CSelectedTransposeFirstRowHeightList,
		class _CSelectedTransposeSecondRowHeightList, class CTransposeColumnWidth>
		class CPreATA_UpperTri_MiddleLoop {
		public:
			typedef typename _CSelectedTransposeFirstRowHeightList::_TyHead
				CCurrentTransposeFirstRowHeight; /**< @brief transpose row height for this decision tree recursion */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_first_block is the first block in the current column
			 *	@param[in] p_second_block_it is the second block iterator, pointing
			 *		to the first of second blocks
			 *	@param[in] p_block_end_it is the current column block list end iterator
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_first_block,
				_TyBlockConstIter p_second_block_it,
				_TyBlockConstIter p_block_end_it,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, _TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				if(r_transpose_row_list[r_t_first_block.first].n_width ==
				   CCurrentTransposeFirstRowHeight::n_size) {
					CPreATA_UpperTri_MiddleLoop<_TyGppContext, CTypelist<CCurrentTransposeFirstRowHeight,
						CTypelistEnd>, _CSelectedTransposeSecondRowHeightList,
						CTransposeColumnWidth>::Loop(r_t_first_block, p_second_block_it,
						p_block_end_it, r_transpose_row_list, r_col_list_dest, r_allocator);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPreATA_UpperTri_MiddleLoop<_TyGppContext,
						typename _CSelectedTransposeFirstRowHeightList::_TyTail,
						_CSelectedTransposeSecondRowHeightList,
						CTransposeColumnWidth>::Loop(r_t_first_block, p_second_block_it,
						p_block_end_it, r_transpose_row_list, r_col_list_dest, r_allocator);
					// not CCurrentTransposeFirstRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#ifdef _OPENMP

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_first_block is the first block in the current column
			 *	@param[in] p_second_block_it is the second block iterator, pointing
			 *		to the first of second blocks
			 *	@param[in] p_block_end_it is the current column block list end iterator
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_first_block,
				_TyBlockConstIter p_second_block_it,
				_TyBlockConstIter p_block_end_it,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, omp_lock_t &r_t_mutex,
				_TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				if(r_transpose_row_list[r_t_first_block.first].n_width ==
				   CCurrentTransposeFirstRowHeight::n_size) {
					CPreATA_UpperTri_MiddleLoop<_TyGppContext, CTypelist<CCurrentTransposeFirstRowHeight,
						CTypelistEnd>, _CSelectedTransposeSecondRowHeightList,
						CTransposeColumnWidth>::Loop(r_t_first_block, p_second_block_it,
						p_block_end_it, r_transpose_row_list, r_col_list_dest, r_t_mutex, r_allocator);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPreATA_UpperTri_MiddleLoop<_TyGppContext,
						typename _CSelectedTransposeFirstRowHeightList::_TyTail,
						_CSelectedTransposeSecondRowHeightList,
						CTransposeColumnWidth>::Loop(r_t_first_block, p_second_block_it,
						p_block_end_it, r_transpose_row_list, r_col_list_dest, r_t_mutex, r_allocator);
					// not CCurrentTransposeFirstRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#endif // _OPENMP

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(r_transpose_row_list[r_t_first_block.first].n_width == %d) {\n",
					CCurrentTransposeFirstRowHeight::n_size);
				CPreATA_UpperTri_MiddleLoop<_TyGppContext, CTypelist<CCurrentTransposeFirstRowHeight,
					CTypelistEnd>, _CSelectedTransposeSecondRowHeightList,
					CTransposeColumnWidth>::Debug();
				printf("} else {\n");
				CPreATA_UpperTri_MiddleLoop<_TyGppContext,
					typename _CSelectedTransposeFirstRowHeightList::_TyTail,
					_CSelectedTransposeSecondRowHeightList, CTransposeColumnWidth>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief middle loop of the part part of the AtA kernel that calculates
		 *		the upper triangle (specialization for the end-of-the-list (branch-less))
		 *
		 *	Iterates through all the second blocks in the current column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CSelectedTransposeFirstRowHeightList is the row height
		 *		(in the transposed matrix) at the end of the list (current);
		 *		(the recursive class specializations peel this one off,
		 *		can't be passed to inner loop)
		 *	@tparam _CSelectedTransposeSecondRowHeightList is list of possible row
		 *		heights in the transposed matrix
		 *	@tparam CTransposeColumnWidth is column width in the transpoed matrix
		 *		(selected in the outer loop)
		 */
		template <class _TyGppContext, class CTransposeFirstRowHeight,
			class _CSelectedTransposeSecondRowHeightList, class CTransposeColumnWidth>
		class CPreATA_UpperTri_MiddleLoop<_TyGppContext, CTypelist<CTransposeFirstRowHeight,
			CTypelistEnd>, _CSelectedTransposeSecondRowHeightList, CTransposeColumnWidth> {
		public:
			/**
			 *	@brief iterates through all the second blocks in the current column
			 *
			 *	@param[in] r_t_first_block is the first block in the current column
			 *	@param[in] p_second_block_it is the second block iterator, pointing
			 *		to the first of second blocks
			 *	@param[in] p_block_end_it is the current column block list end iterator
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_first_block,
				_TyBlockConstIter p_second_block_it,
				_TyBlockConstIter p_block_end_it,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, _TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				_ASSERTE(r_transpose_row_list[r_t_first_block.first].n_width ==
					CTransposeFirstRowHeight::n_size);
				// make sure the block has correct height

				for(; p_second_block_it != p_block_end_it; ++ p_second_block_it) {
					const TColumn::TBlockEntry &r_t_second_block = *p_second_block_it;

					// multiplication of blockA0 * blockA1 yields block at (n_row_id_A0, n_row_id_A1),
					// guaranteed to be above diagonal

					CPreATA_UpperTri_InnerLoop<_TyGppContext, _CSelectedTransposeSecondRowHeightList,
						CTransposeFirstRowHeight, CTransposeColumnWidth>::Loop(r_t_first_block,
						r_t_second_block, r_transpose_row_list, r_col_list_dest, r_allocator);
					// wrap the inner loop in a (second) block size decision tree
				}
				// for each next block below the current block in the current column in A
			}

#ifdef _OPENMP

			/**
			 *	@brief iterates through all the second blocks in the current column
			 *
			 *	@param[in] r_t_first_block is the first block in the current column
			 *	@param[in] p_second_block_it is the second block iterator, pointing
			 *		to the first of second blocks
			 *	@param[in] p_block_end_it is the current column block list end iterator
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn::TBlockEntry &r_t_first_block,
				_TyBlockConstIter p_second_block_it,
				_TyBlockConstIter p_block_end_it,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, omp_lock_t &r_t_mutex,
				_TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				_ASSERTE(r_transpose_row_list[r_t_first_block.first].n_width ==
					CTransposeFirstRowHeight::n_size);
				// make sure the block has correct height

				for(; p_second_block_it != p_block_end_it; ++ p_second_block_it) {
					const TColumn::TBlockEntry &r_t_second_block = *p_second_block_it;

					// multiplication of blockA0 * blockA1 yields block at (n_row_id_A0, n_row_id_A1),
					// guaranteed to be above diagonal

					CPreATA_UpperTri_InnerLoop<_TyGppContext, _CSelectedTransposeSecondRowHeightList,
						CTransposeFirstRowHeight, CTransposeColumnWidth>::Loop(r_t_first_block,
						r_t_second_block, r_transpose_row_list, r_col_list_dest,
						r_t_mutex, r_allocator);
					// wrap the inner loop in a (second) block size decision tree
				}
				// for each next block below the current block in the current column in A
			}

#endif // _OPENMP

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(r_transpose_row_list[r_t_first_block.first].n_width == %d);\n",
					CTransposeFirstRowHeight::n_size);
				printf("for(each p_second_block_it = p_first_block_it + 1 in "
					"r_t_column.block_list) {\nconst TColumn::TBlockEntry &r_t_second_block"
					" = *p_second_block_it;\n");
				printf("// multiplication of blockA0 * blockA1 yields block at (n_row_id_A0,"
					" n_row_id_A1),\n");
				CPreATA_UpperTri_InnerLoop<_TyGppContext, _CSelectedTransposeSecondRowHeightList,
					CTransposeFirstRowHeight, CTransposeColumnWidth>::Debug();
				printf("}\n");
				printf("// for each next block below the current"
					" block in the current column in A\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief outer loop of the part part of the AtA kernel that calculates
		 *		the upper triangle
		 *
		 *	Iterates through all the first blocks in the current column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CTransposeColumnWidthList is list of possible column widths
		 *		in the transpoed matrix
		 */
		template <class _TyGppContext, class _CTransposeColumnWidthList>
		class CPreATA_UpperTri_OuterLoop {
		public:
			typedef typename _CTransposeColumnWidthList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn &r_t_column,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, _TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				if(r_t_column.n_width == CCurrentColumnWidth::n_size) {
					CPreATA_UpperTri_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::Loop(r_t_column, r_transpose_row_list, r_col_list_dest, r_allocator);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPreATA_UpperTri_OuterLoop<_TyGppContext,
						typename _CTransposeColumnWidthList::_TyTail>::Loop(r_t_column,
						r_transpose_row_list, r_col_list_dest, r_allocator);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

#ifdef _OPENMP

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn &r_t_column,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, omp_lock_t &r_t_mutex,
				_TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				if(r_t_column.n_width == CCurrentColumnWidth::n_size) {
					CPreATA_UpperTri_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::Loop(r_t_column, r_transpose_row_list,
						r_col_list_dest, r_t_mutex, r_allocator);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPreATA_UpperTri_OuterLoop<_TyGppContext,
						typename _CTransposeColumnWidthList::_TyTail>::Loop(r_t_column,
						r_transpose_row_list, r_col_list_dest, r_t_mutex, r_allocator);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

#endif // _OPENMP

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(r_t_column.n_width == %d) {\n", CCurrentColumnWidth::n_size);
				CPreATA_UpperTri_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth, CTypelistEnd> >::Debug();
				printf("} else {\n");
				CPreATA_UpperTri_OuterLoop<_TyGppContext, typename _CTransposeColumnWidthList::_TyTail>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief outer loop of the part part of the AtA kernel that calculates
		 *		the upper triangle (specialization for the end-of-the-list (branch-less))
		 *
		 *	Iterates through all the first blocks in the current column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CTransposeColumnWidth is width of the column
		 *		(in the transpoed matrix) at the end of the list (current)
		 */
		template <class _TyGppContext, class CTransposeColumnWidth>
		class CPreATA_UpperTri_OuterLoop<_TyGppContext, CTypelist<CTransposeColumnWidth, CTypelistEnd> > {
		public:
			typedef typename CFilterTypelist2<CColumnWidthsList, CTransposeColumnWidth,
				fbs_ut::CHaveColumnWidthForRowHeight, CDimsList_Uniq>::_TyResult _CSelectedTransposeRowHeightsList; /**< @brief transpose row heights list, filtered by the selected transpose column width */

			/**
			 *	@brief iterates through all the first blocks in the current column
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn &r_t_column,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, _TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				_ASSERTE(r_t_column.block_list.empty() || r_t_column.n_width == CTransposeColumnWidth::n_size);
				// make sure the column has correct width

				for(_TyBlockConstIter p_first_block_it =
				   r_t_column.block_list.begin(), p_block_end_it = r_t_column.block_list.end();
				   p_first_block_it != p_block_end_it;) { // increments at the beginning of the inner loop
					const TColumn::TBlockEntry &r_t_first_block = *p_first_block_it;

					CPreATA_UpperTri_MiddleLoop<_TyGppContext, _CSelectedTransposeRowHeightsList,
						_CSelectedTransposeRowHeightsList, CTransposeColumnWidth>::Loop(
						r_t_first_block, ++ p_first_block_it, p_block_end_it,
						r_transpose_row_list, r_col_list_dest, r_allocator);
					// wrap the middle loop in a (first) block size decision tree
				}
				// for each block in the current column in A
			}

#ifdef _OPENMP

			/**
			 *	@brief iterates through all the first blocks in the current column
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in] r_transpose_row_list is list of transpose matrix rows
			 *	@param[out] r_col_list_dest is list of destination matrix columns
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn &r_t_column,
				const std::vector<TColumn> &r_transpose_row_list,
				std::vector<TColumn> &r_col_list_dest, omp_lock_t &r_t_mutex,
				_TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
			{
				_ASSERTE(r_t_column.block_list.empty() || r_t_column.n_width == CTransposeColumnWidth::n_size);
				// make sure the column has correct width

				for(_TyBlockConstIter p_first_block_it =
				   r_t_column.block_list.begin(), p_block_end_it = r_t_column.block_list.end();
				   p_first_block_it != p_block_end_it;) { // increments at the beginning of the inner loop
					const TColumn::TBlockEntry &r_t_first_block = *p_first_block_it;

					CPreATA_UpperTri_MiddleLoop<_TyGppContext, _CSelectedTransposeRowHeightsList,
						_CSelectedTransposeRowHeightsList, CTransposeColumnWidth>::Loop(
						r_t_first_block, ++ p_first_block_it, p_block_end_it,
						r_transpose_row_list, r_col_list_dest, r_t_mutex, r_allocator);
					// wrap the middle loop in a (first) block size decision tree
				}
				// for each block in the current column in A
			}

#endif // _OPENMP

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(r_t_column.n_width == %d);\n", CTransposeColumnWidth::n_size);
				printf("for(each p_first_block_it in r_t_column.block_list) {\n"
					"const TColumn::TBlockEntry &r_t_first_block = *p_first_block_it;\n");
				CPreATA_UpperTri_MiddleLoop<_TyGppContext, _CSelectedTransposeRowHeightsList,
					_CSelectedTransposeRowHeightsList, CTransposeColumnWidth>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief inner loop of the part part of the AtA kernel that calculates the diagonal
		 *
		 *	Performs a single block multiplication.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CSelectedRowHeightsList is list of possible row heights
		 *	@tparam CColumnWidth is column width (selected in the outer loop)
		 */
		template <class _TyGppContext, class _CSelectedRowHeightsList, class CColumnWidth>
		class CPreATA_Diagonal_InnerLoop {
		public:
			typedef typename _CSelectedRowHeightsList::_TyHead CCurrentRowHeight; /**< @brief row height for this decision tree recursion */
			typedef typename CMakeMatrixRef<CColumnWidth::n_size, CColumnWidth::n_size>::_Ty _TyDestBlock; /**< @brief destination matrix block Eigen::Map type */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[out] r_block_dest is the destination block (dot product accumulator)
			 *	@param[in] r_t_block is source block to be multiplied with its transpose
			 *	@param[in] r_row_list is list of matrix rows
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDestBlock &r_block_dest,
				const TColumn::TBlockEntry &r_t_block, const std::vector<TRow> &r_row_list) // t_odo - follow naming conventions
			{
				if(r_row_list[r_t_block.first].n_height == CCurrentRowHeight::n_size) {
					CPreATA_Diagonal_InnerLoop<_TyGppContext, CTypelist<CCurrentRowHeight, CTypelistEnd>,
						CColumnWidth>::Loop(r_block_dest, r_t_block, r_row_list);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CPreATA_Diagonal_InnerLoop<_TyGppContext, typename _CSelectedRowHeightsList::_TyTail,
						CColumnWidth>::Loop(r_block_dest, r_t_block, r_row_list);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(r_row_list[r_t_block.first].n_height == %d) {\n",
					CCurrentRowHeight::n_size);
				CPreATA_Diagonal_InnerLoop<_TyGppContext, CTypelist<CCurrentRowHeight,
					CTypelistEnd>, CColumnWidth>::Debug();
				printf("} else {\n");
				CPreATA_Diagonal_InnerLoop<_TyGppContext, typename _CSelectedRowHeightsList::_TyTail,
					CColumnWidth>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief inner loop of the part part of the AtA kernel that calculates the diagonal
		 *		(specialization for the end-of-the-list (branch-less))
		 *
		 *	Performs a single block multiplication.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CRowHeight is height of the row at the end of the list (current)
		 *	@tparam CColumnWidth is column width (selected in the outer loop)
		 */
		template <class _TyGppContext, class CRowHeight, class CColumnWidth>
		class CPreATA_Diagonal_InnerLoop<_TyGppContext, CTypelist<CRowHeight, CTypelistEnd>, CColumnWidth> {
		public:
			typedef typename CMakeMatrixRef<CRowHeight::n_size, CColumnWidth::n_size>::_Ty _TyBlock; /**< @brief source matrix block Eigen::Map type */
			typedef typename CMakeMatrixRef<CColumnWidth::n_size, CColumnWidth::n_size>::_Ty _TyDestBlock; /**< @brief destination matrix block Eigen::Map type */

			/**
			 *	@brief performs a single block multiplication
			 *
			 *	@param[out] r_block_dest is the destination block (dot product accumulator)
			 *	@param[in] r_t_block is source block to be multiplied with its transpose
			 *	@param[in] r_row_list is list of matrix rows
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDestBlock &r_block_dest,
				const TColumn::TBlockEntry &r_t_block,
				const std::vector<TRow> &UNUSED(r_row_list))
			{
				const TColumn::TBlockEntry &r_t_block_B = r_t_block;
				size_t UNUSED(n_row_id_B) = r_t_block_B.first; // get row of the block
				// for each block in the current column in B

				_ASSERTE(r_row_list[n_row_id_B].n_height == CRowHeight::n_size);
				// make sure that the row has correct size

				_TyBlock blockB(r_t_block_B.second);
				// create map to block B data

#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				r_block_dest += blockB.transpose().lazyProduct(blockB); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				r_block_dest.noalias() += blockB.transpose() * blockB;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				// multiply block by its transpose (hope eigen can really optimize this)
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(r_row_list[r_t_block.first].n_height == %d);\n",
					CRowHeight::n_size);
				printf("block_dest += block.transpose() * block;\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief outer loop of the part part of the AtA kernel that calculates the diagonal
		 *
		 *	Iterates through all the blocks in the column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthList is list of possible column widths
		 */
		template <class _TyGppContext, class _CColumnWidthList>
		class CPreATA_Diagonal_OuterLoop {
		public:
			typedef typename _CColumnWidthList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *	@param[in] r_row_list is list of matrix rows
			 *
			 *	@return Returns pointer to allocated and initialized block data.
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline double *p_Loop(const TColumn &r_t_column,
				_TyDenseAllocator &r_allocator, const std::vector<TRow> &r_row_list) // throw(std::bad_alloc)
			{
				if(r_t_column.n_width == CCurrentColumnWidth::n_size) {
					return CPreATA_Diagonal_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::p_Loop(r_t_column, r_allocator, r_row_list);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					return CPreATA_Diagonal_OuterLoop<_TyGppContext,
						typename _CColumnWidthList::_TyTail>::p_Loop(r_t_column,
						r_allocator, r_row_list);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in] p_block_data is pointer to preallocated block data (in the destination matrix)
			 *	@param[in] r_row_list is list of matrix rows
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn &r_t_column,
				double *p_block_data, const std::vector<TRow> &r_row_list)
			{
				if(r_t_column.n_width == CCurrentColumnWidth::n_size) {
					return CPreATA_Diagonal_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::Loop(r_t_column, p_block_data, r_row_list);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					return CPreATA_Diagonal_OuterLoop<_TyGppContext,
						typename _CColumnWidthList::_TyTail>::Loop(r_t_column,
						p_block_data, r_row_list);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(r_t_column.n_width == %d) {\n", CCurrentColumnWidth::n_size);
				CPreATA_Diagonal_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth, CTypelistEnd> >::Debug();
				printf("} else {\n");
				CPreATA_Diagonal_OuterLoop<_TyGppContext, typename _CColumnWidthList::_TyTail>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief outer loop of the part part of the AtA kernel that calculates the diagonal
		 *		(specialization for the end-of-the-list (branch-less))
		 *
		 *	Iterates through all the blocks in the column.
		 *
		 *	@tparam CColumnWidth is width of the column at the end of the list (current)
		 */
		template <class _TyGppContext, class CColumnWidth>
		class CPreATA_Diagonal_OuterLoop<_TyGppContext, CTypelist<CColumnWidth, CTypelistEnd> > {
		public:
			typedef typename CMakeMatrixRef<CColumnWidth::n_size, CColumnWidth::n_size>::_Ty _TyDestBlock; /**< @brief destination matrix block Eigen::Map type */
			typedef typename CFilterTypelist2<CRowHeightsList, CColumnWidth,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_Uniq>::_TyResult _CSelectedRowHeightsList; /**< @brief row heights list, filtered by the selected column width */

			/**
			 *	@brief iterates through all the blocks in the column
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in] r_allocator is block allocator, associated with the destination matrix
			 *	@param[in] r_row_list is list of matrix rows
			 *
			 *	@return Returns pointer to allocated and initialized block data.
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline double *p_Loop(const TColumn &r_t_column,
				_TyDenseAllocator &r_allocator, const std::vector<TRow> &r_row_list) // throw(std::bad_alloc)
			{
				_ASSERTE(r_t_column.block_list.empty() || r_t_column.n_width == CColumnWidth::n_size);
				// makes sure the column has the right size

				const size_t n_bm_cols = CColumnWidth::n_size; // B^T * B has shape of B.cols by B.cols
				double *p_block_data = r_allocator.p_Get_DenseStorage(n_bm_cols * n_bm_cols);
				//memset(p_block_data, 0, n_bm_cols * n_bm_cols * sizeof(double)); // or handle the first block differently
				_TyDestBlock block_dest(p_block_data);
				block_dest.setZero(); // might use SSE
				// allocate new block

				for(_TyBlockConstIter p_block_it = r_t_column.block_list.begin(),
				   p_block_end_it = r_t_column.block_list.end(); p_block_it != p_block_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;

					CPreATA_Diagonal_InnerLoop<_TyGppContext, _CSelectedRowHeightsList, CColumnWidth>::Loop(block_dest,
						r_t_block, r_row_list);
					// wrap the inner loop in a block size decision tree

					// t_odo - remove B from varname, call innerloop
					// t_odo - write debug, debug
					// t_odo - test, run VP, run speed comparisons with CSparse
				}
				// calculate dot product of all the blocks in this column

				return p_block_data;
			}

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in] p_block_data is pointer to preallocated block data (in the destination matrix)
			 *	@param[in] r_row_list is list of matrix rows
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(const TColumn &r_t_column,
				double *p_block_data, const std::vector<TRow> &r_row_list)
			{
				_ASSERTE(r_t_column.block_list.empty() || r_t_column.n_width == CColumnWidth::n_size);
				// makes sure the column has the right size

				//const size_t n_bm_cols = CColumnWidth::n_size; // B^T * B has shape of B.cols by B.cols
				//double *p_block_data = r_allocator.p_Get_DenseStorage(n_bm_cols * n_bm_cols); // allocated outside
				//memset(p_block_data, 0, n_bm_cols * n_bm_cols * sizeof(double)); // or handle the first block differently
				_TyDestBlock block_dest(p_block_data);
				block_dest.setZero(); // might use SSE
				// allocate new block

				for(_TyBlockConstIter p_block_it = r_t_column.block_list.begin(),
				   p_block_end_it = r_t_column.block_list.end(); p_block_it != p_block_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;

					CPreATA_Diagonal_InnerLoop<_TyGppContext, _CSelectedRowHeightsList, CColumnWidth>::Loop(block_dest,
						r_t_block, r_row_list);
					// wrap the inner loop in a block size decision tree

					// t_odo - remove B from varname, call innerloop
					// t_odo - write debug, debug
					// t_odo - test, run VP, run speed comparisons with CSparse
				}
				// calculate dot product of all the blocks in this column
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(r_t_column.n_width == %d);\n", CColumnWidth::n_size);
				printf("_TyDestBlock block_dest(r_allocator.p_Get_DenseStorage(%d * %d));\n",
					CColumnWidth::n_size, CColumnWidth::n_size);
				printf("for(each p_block_it in r_t_column.block_list) {\n"
					"const TColumn::TBlockEntry &r_t_block = *p_block_it;\n");
				CPreATA_Diagonal_InnerLoop<_TyGppContext, _CSelectedRowHeightsList, CColumnWidth>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief iterates through all the first blocks in the current column
		 *
		 *	@param[in] r_t_column is the current column
		 *	@param[in] r_transpose_row_list is list of transpose matrix rows
		 *	@param[out] r_col_list_dest is list of destination matrix columns
		 *	@param[in] r_allocator is block allocator, associated with the destination matrix
		 *
		 *	@note This function throws std::bad_alloc.
		 */
		static __forceinline void PreATA_UpperTriangle_OuterLoop(const TColumn &r_t_column,
			const std::vector<TColumn> &r_transpose_row_list,
			std::vector<TColumn> &r_col_list_dest, _TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
		{
			CPreATA_UpperTri_OuterLoop<CBlockMatrixTypelist, CRowHeightsList>::Loop(r_t_column,
				r_transpose_row_list, r_col_list_dest, r_allocator);
			// wrap the outer loop body in a transpose column size decision tree
		}

#ifdef _OPENMP

		/**
		 *	@brief iterates through all the first blocks in the current column
		 *
		 *	@param[in] r_t_column is the current column
		 *	@param[in] r_transpose_row_list is list of transpose matrix rows
		 *	@param[out] r_col_list_dest is list of destination matrix columns
		 *	@param[in] r_allocator is block allocator, associated with the destination matrix
		 *
		 *	@note This function throws std::bad_alloc.
		 */
		static __forceinline void PreATA_UpperTriangle_Parallel_OuterLoop(const TColumn &r_t_column,
			const std::vector<TColumn> &r_transpose_row_list,
			std::vector<TColumn> &r_col_list_dest, omp_lock_t &r_t_mutex,
			_TyDenseAllocator &r_allocator) // throw(std::bad_alloc)
		{
			CPreATA_UpperTri_OuterLoop<CBlockMatrixTypelist, CRowHeightsList>::Loop(r_t_column,
				r_transpose_row_list, r_col_list_dest, r_t_mutex, r_allocator);
			// wrap the outer loop body in a transpose column size decision tree
		}

#endif // _OPENMP

		/**
		 *	@brief iterates through all the blocks in the column
		 *
		 *	@param[in] r_t_column is the current column
		 *	@param[in] r_allocator is block allocator, associated with the destination matrix
		 *	@param[in] r_row_list is list of matrix rows
		 *
		 *	@return Returns pointer to a block, calculated in the loop
		 *		(allocated inside, but not added to the destination matrix).
		 *
		 *	@note This function throws std::bad_alloc.
		 */
		static __forceinline double *p_PreATA_Diagonal_OuterLoop(const TColumn &r_t_column,
			_TyDenseAllocator &r_allocator, const std::vector<TRow> &r_row_list) // throw(std::bad_alloc)
		{
			return CPreATA_Diagonal_OuterLoop<CBlockMatrixTypelist, CColumnWidthsList>::p_Loop(r_t_column,
				r_allocator, r_row_list);
			// wrap the outer loop body in a column size decision tree
		}

		/**
		 *	@brief iterates through all the blocks in the column
		 *
		 *	@param[in] r_t_column is the current column
		 *	@param[in] p_block_data is pointer to preallocated block data (in the destination matrix)
		 *	@param[in] r_row_list is list of matrix rows
		 */
		static __forceinline void PreATA_Parallel_Diagonal_OuterLoop(const TColumn &r_t_column,
			double *p_block_data, const std::vector<TRow> &r_row_list)
		{
			CPreATA_Diagonal_OuterLoop<CBlockMatrixTypelist,
				CColumnWidthsList>::Loop(r_t_column, p_block_data, r_row_list);
			// wrap the outer loop body in a column size decision tree
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unindented) pseudocode of the pre-multiply with
		 *		self transpose function to stdout
		 *	@note This function is only available if the
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING macro is defined.
		 */
		static void Debug_PreMultiplyWithSelfTransposeTo()
		{
			printf("for(each p_col_A_it in col_list_A) {\n"
				"const TColumn &r_t_column = *p_col_A_it; // for each column in A\n");
			CPreATA_UpperTri_OuterLoop<CBlockMatrixTypelist, CRowHeightsList>::Debug();
			printf("}\n"
				"// perform the multiplication of the blocks above diagonal,"
				" need to use log(n) lookup\n\n");

			printf("size_t n_column_id_B = 0;\n"
				"for(each p_col_B_it in cols_B, ++ n_column_id_B) {\n"
				"const TColumn &r_t_column_B = *p_col_B_it;\n"
				"if(r_t_column_B.block_list.empty())\n"
				"continue;\n"
				"// only process non-empty columns\n");
			CPreATA_Diagonal_OuterLoop<CBlockMatrixTypelist, CColumnWidthsList>::Debug();
			printf("}\n"
				"// perform the multiplication of the blocks on the diagonal,"
				" (const time lookup, fast multiplication)\n\n");
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief fixed block size Cholesky template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 *	@todo write documentations for members, sort members, write debug functions
	 */
	template <class CBlockMatrixTypelist>
	class CFBS_Cholesky : public fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist> {
	public:
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CColumnWidthsList CColumnWidthsList; /**< @brief list of unique block column widths */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CDimsList_Uniq CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */

		/**
		 *	@brief off-diagonal loop of the Cholesky factorization kernel
		 *
		 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, k)
		 *	@param[in] n_col_j_width is width of the j-th column, in elements
		 *	@param[in] n_col_k_width is width of the k-th column, in elements
		 *		(equal to the number of rows of the dest block)
		 *	@param[in] p_j_block_it is iterator to blocks in j-th column
		 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
		 *	@param[in] p_k_block_it is iterator to blocks in k-th column
		 *	@param[in] p_k_block_end_it is iterator pointing one past the last block in k-th column
		 *	@param[in] t_block_A is input block at A(j, k), same size as the dest block
		 *	@param[in] k is index of the previous column
		 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
		 */
		static __forceinline void OffDiagonal_Loop(double *p_dest_block, size_t n_col_j_width,
			size_t n_col_k_width, _TyBlockConstIter p_j_block_it, _TyBlockConstIter p_j_block_end_it,
			_TyBlockConstIter p_k_block_it, _TyBlockConstIter p_k_block_end_it,
			TColumn::TBlockEntry t_block_A, size_t k, const std::vector<TColumn> &r_block_cols_list)
		{
			CCholesky_OffDiagonalLoop<CBlockMatrixTypelist, CColumnWidthsList>::Loop(p_dest_block,
				n_col_j_width, n_col_k_width, p_j_block_it, p_j_block_end_it, p_k_block_it,
				p_k_block_end_it, t_block_A, k, r_block_cols_list);
		}

		/**
		 *	@brief off-diagonal loop of the Cholesky factorization kernel
		 *
		 *	Performs calculation of the off-diagonal elements in one column of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthList is list of possible column widths
		 */
		template <class _TyGppContext, class _CColumnWidthList>
		class CCholesky_OffDiagonalLoop {
		public:
			typedef typename _CColumnWidthList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Wraps the inside of the loop in a decision tree for j-th column width.
			 *
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, k)
			 *	@param[in] n_col_j_width is width of the j-th column, in elements
			 *	@param[in] n_col_k_width is width of the k-th column, in elements
			 *		(equal to the number of rows of the dest block)
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_k_block_it is iterator to blocks in k-th column
			 *	@param[in] p_k_block_end_it is iterator pointing one past the last block in k-th column
			 *	@param[in] t_block_A is input block at A(j, k), same size as the dest block
			 *	@param[in] k is index of the previous column
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(double *p_dest_block, size_t n_col_j_width,
				size_t n_col_k_width, _TyBlockConstIter p_j_block_it,
				_TyBlockConstIter p_j_block_end_it, _TyBlockConstIter p_k_block_it,
				_TyBlockConstIter p_k_block_end_it, TColumn::TBlockEntry t_block_A, size_t k,
				const std::vector<TColumn> &r_block_cols_list)
			{
				if(n_col_j_width == CCurrentColumnWidth::n_size) {
					CCholesky_OffDiagonalLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::Loop(p_dest_block, n_col_j_width, n_col_k_width,
						p_j_block_it, p_j_block_end_it, p_k_block_it, p_k_block_end_it,
						t_block_A, k, r_block_cols_list);
				} else {
					CCholesky_OffDiagonalLoop<_TyGppContext, typename
						_CColumnWidthList::_TyTail>::Loop(p_dest_block, n_col_j_width,
						n_col_k_width, p_j_block_it, p_j_block_end_it, p_k_block_it,
						p_k_block_end_it, t_block_A, k, r_block_cols_list);
				}
			}
		};

		/**
		 *	@brief off-diagonal loop of the Cholesky factorization kernel
		 *
		 *	Performs calculation of the off-diagonal elements in one column of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CColumnWidth_j is selected column width for this decision tree branch
		 */
		template <class _TyGppContext, class CColumnWidth_j>
		class CCholesky_OffDiagonalLoop<_TyGppContext, CTypelist<CColumnWidth_j, CTypelistEnd> > {
		public:
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_TyMatrix _TyDiagMatrix; /**< @brief a dense matrix with the same shape as the diagonal block */

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Just calls the next level of the decision tree.
			 *
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, k)
			 *	@param[in] n_col_j_width is width of the j-th column, in elements
			 *	@param[in] n_col_k_width is width of the k-th column, in elements
			 *		(equal to the number of rows of the dest block)
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_k_block_it is iterator to blocks in k-th column
			 *	@param[in] p_k_block_end_it is iterator pointing one past the last block in k-th column
			 *	@param[in] t_block_A is input block at A(j, k), same size as the dest block
			 *	@param[in] k is index of the previous column
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(double *p_dest_block, size_t UNUSED(n_col_j_width),
				size_t n_col_k_width, _TyBlockConstIter p_j_block_it,
				_TyBlockConstIter p_j_block_end_it, _TyBlockConstIter p_k_block_it,
				_TyBlockConstIter p_k_block_end_it, TColumn::TBlockEntry t_block_A, size_t k,
				const std::vector<TColumn> &r_block_cols_list)
			{
				// this is now empty, but one of the loops could be inside here
				// (did i not put it here because of some data dependences?)

				_ASSERTE(n_col_j_width == CColumnWidth_j::n_size);
				CCholesky_OffDiagonalLoop1<_TyGppContext, CColumnWidthsList,
					CColumnWidth_j>::Loop(p_dest_block, n_col_k_width,
					p_j_block_it, p_j_block_end_it, p_k_block_it, p_k_block_end_it,
					t_block_A, k, r_block_cols_list);
				// decide on the second column width
			}

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Just calls the next level of the decision tree.
			 *
			 *	@param[out] r_column_dot is the running dor product of the column blocks
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, k)
			 *	@param[in] n_col_k_width is width of the k-th column, in elements
			 *		(equal to the number of rows of the dest block)
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_k_block_it is iterator to blocks in k-th column
			 *	@param[in] p_k_block_end_it is iterator pointing one past the last block in k-th column
			 *	@param[in] t_block_A is input block at A(j, k), same size as the dest block
			 *	@param[in] k is index of the previous column
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDiagMatrix &r_column_dot, double *p_dest_block,
				size_t n_col_k_width, _TyBlockConstIter p_j_block_it,
				_TyBlockConstIter p_j_block_end_it, _TyBlockConstIter p_k_block_it,
				_TyBlockConstIter p_k_block_end_it, TColumn::TBlockEntry t_block_A, size_t k,
				const std::vector<TColumn> &r_block_cols_list)
			{
				// this is now empty, but one of the loops could be inside here
				// (did i not put it here because of some data dependences?)

				CCholesky_OffDiagonalLoop1<_TyGppContext, CColumnWidthsList,
					CColumnWidth_j>::Loop(r_column_dot, p_dest_block, n_col_k_width,
					p_j_block_it, p_j_block_end_it, p_k_block_it, p_k_block_end_it,
					t_block_A, k, r_block_cols_list);
				// decide on the second column width
			}
		};

		/**
		 *	@brief off-diagonal loop of the Cholesky factorization kernel
		 *
		 *	Performs calculation of the off-diagonal elements in one column of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthList is list of possible column widths
		 *	@tparam CColumnWidth_j is width of the j-th column
		 */
		template <class _TyGppContext, class _CColumnWidthList, class CColumnWidth_j>
		class CCholesky_OffDiagonalLoop1 {
		public:
			typedef typename _CColumnWidthList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_TyMatrix _TyDiagMatrix; /**< @brief a dense matrix with the same shape as the diagonal block */

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Wraps the inside of the loop in a decision tree for k-th column width.
			 *
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, k)
			 *	@param[in] n_col_k_width is width of the k-th column, in elements
			 *		(equal to the number of rows of the dest block)
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_k_block_it is iterator to blocks in k-th column
			 *	@param[in] p_k_block_end_it is iterator pointing one past the last block in k-th column
			 *	@param[in] t_block_A is input block at A(j, k), same size as the dest block
			 *	@param[in] k is index of the previous column
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(double *p_dest_block,
				size_t n_col_k_width, _TyBlockConstIter UNUSED(p_j_block_it),
				_TyBlockConstIter p_j_block_end_it, _TyBlockConstIter p_k_block_it,
				_TyBlockConstIter p_k_block_end_it, TColumn::TBlockEntry t_block_A, size_t k,
				const std::vector<TColumn> &r_block_cols_list)
			{
				if(n_col_k_width == CCurrentColumnWidth::n_size) {
					CCholesky_OffDiagonalLoop1<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd>, CColumnWidth_j>::Loop(p_dest_block,
						n_col_k_width, p_j_block_it, p_j_block_end_it, p_k_block_it,
						p_k_block_end_it, t_block_A, k, r_block_cols_list);
				} else {
					CCholesky_OffDiagonalLoop1<_TyGppContext, typename _CColumnWidthList::_TyTail,
						CColumnWidth_j>::Loop(p_dest_block, n_col_k_width,
						p_j_block_it, p_j_block_end_it, p_k_block_it, p_k_block_end_it,
						t_block_A, k, r_block_cols_list);
				}
			}

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Wraps the inside of the loop in a decision tree for k-th column width.
			 *
			 *	@param[out] r_column_dot is the running dor product of the column blocks
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, k)
			 *	@param[in] n_col_k_width is width of the k-th column, in elements
			 *		(equal to the number of rows of the dest block)
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_k_block_it is iterator to blocks in k-th column
			 *	@param[in] p_k_block_end_it is iterator pointing one past the last block in k-th column
			 *	@param[in] t_block_A is input block at A(j, k), same size as the dest block
			 *	@param[in] k is index of the previous column
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDiagMatrix &r_column_dot, double *p_dest_block,
				size_t n_col_k_width, _TyBlockConstIter UNUSED(p_j_block_it),
				_TyBlockConstIter p_j_block_end_it, _TyBlockConstIter p_k_block_it,
				_TyBlockConstIter p_k_block_end_it, TColumn::TBlockEntry t_block_A, size_t k,
				const std::vector<TColumn> &r_block_cols_list)
			{
				if(n_col_k_width == CCurrentColumnWidth::n_size) {
					CCholesky_OffDiagonalLoop1<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd>, CColumnWidth_j>::Loop(r_column_dot, p_dest_block,
						n_col_k_width, p_j_block_it, p_j_block_end_it, p_k_block_it,
						p_k_block_end_it, t_block_A, k, r_block_cols_list);
				} else {
					CCholesky_OffDiagonalLoop1<_TyGppContext, typename _CColumnWidthList::_TyTail,
						CColumnWidth_j>::Loop(r_column_dot, p_dest_block, n_col_k_width,
						p_j_block_it, p_j_block_end_it, p_k_block_it, p_k_block_end_it,
						t_block_A, k, r_block_cols_list);
				}
			}
		};

		/**
		 *	@brief off-diagonal loop of the Cholesky factorization kernel
		 *
		 *	Performs calculation of the off-diagonal elements in one column of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CColumnWidth_k is width of the k-th column
		 *	@tparam CColumnWidth_j is width of the j-th column
		 */
		template <class _TyGppContext, class CColumnWidth_k, class CColumnWidth_j>
		class CCholesky_OffDiagonalLoop1<_TyGppContext, CTypelist<CColumnWidth_k, CTypelistEnd>, CColumnWidth_j> {
		public:
			typedef typename CMakeMatrixRef<CColumnWidth_k::n_size, CColumnWidth_j::n_size>::_Ty _TyDestBlock; /**< @brief destination block shape */
			typedef typename CMakeMatrixRef<CColumnWidth_k::n_size, CColumnWidth_k::n_size>::_Ty _TyDiagBlock_k; /**< @brief shape of the diagonal block at the k-th column */
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_TyMatrix _TyDiagMatrix; /**< @brief a dense matrix with the same shape as the diagonal block */

			/**
			 *	@brief occurence of an item in a list predicate
			 *
			 *	@tparam _TyItem is an item
			 *	@tparam _TyList is a typelist, tested for the occurence of the given item
			 */
			template <class _TyItem, class _TyList>
			class CContainsItem {
			public:
				/**
				 *	@brief the result, stored as enum
				 */
				enum {
					b_result = CFindTypelistItem<_TyList, _TyItem>::b_result /**< @brief the result; set if the list contains the item, otherwise cleared */
					// unfortunately this has reversed order of arguments (list, needle) instead of (needle, list)
				};
			};

			typedef typename CFilterTypelist2<CColumnWidthsList, CColumnWidth_k,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_Uniq>::_TyResult
				_CSelectedRowHeightsList_on_k; /**< @brief list of row heights in column k */
			typedef typename CFilterTypelist2<CColumnWidthsList, CColumnWidth_j,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_Uniq>::_TyResult
				_CSelectedRowHeightsList_on_j; /**< @brief list of row heights in column j */
			typedef typename CFilterTypelist<_CSelectedRowHeightsList_on_k, _CSelectedRowHeightsList_on_j,
				CContainsItem>::_TyResult _CSelectedRowHeightsList; /**< @brief list of row heights that can occur in both columns at the same time (optimization for matrices with >1 block sizes) */

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Performs the actual loop.
			 *
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, k)
			 *	@param[in] n_col_k_width is width of the k-th column, in elements
			 *		(equal to the number of rows of the dest block)
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_k_block_it is iterator to blocks in k-th column
			 *	@param[in] p_k_block_end_it is iterator pointing one past the last block in k-th column
			 *	@param[in] t_block_A is input block at A(j, k), same size as the dest block
			 *	@param[in] k is index of the previous column
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(double *p_dest_block,
				size_t UNUSED(n_col_k_width), _TyBlockConstIter p_j_block_it,
				_TyBlockConstIter p_j_block_end_it, _TyBlockConstIter p_k_block_it,
				_TyBlockConstIter p_k_block_end_it, TColumn::TBlockEntry t_block_A, size_t k,
				const std::vector<TColumn> &r_block_cols_list)
			{
				_ASSERTE(n_col_k_width == CColumnWidth_k::n_size);
				// make sure this is the correct block size

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					p_k_block_it != p_k_block_end_it);
				// if the first is a non-empty range, the second should be neither

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					/*r_col_L_k.block_list.back()*/(*(p_k_block_end_it - 1)).first == k);
				// the last block of k-th column is on diagonal (kth row) and it is always present
				// this block will serve as the sentinell for this loop

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					(*(p_j_block_end_it - 1)).first < k); // note that the p_j_block_end_it points at the diagonal block
				// the last block of j-th column (omitting the diagonal) is *at most* at kth row

				_TyDestBlock L_block_kj(p_dest_block);
				// the destination block

				if(t_block_A.first == k)
					L_block_kj = _TyDestBlock(t_block_A.second);
				else
					L_block_kj.setZero();
				// copy data from the source matrix

				{
					for(; p_j_block_it != p_j_block_end_it; ++ p_j_block_it) { // for each block in column j
						size_t n_row_i = (*p_j_block_it).first;
						_ASSERTE(n_row_i < k); // make sure that the sentinell is indeed functional
						_ASSERTE(p_k_block_it != p_k_block_end_it); // should not be pointing at the end in the first place (if we got inside this loop)
						while((*p_k_block_it).first < n_row_i) {
							++ p_k_block_it;
							_ASSERTE(p_k_block_it != p_k_block_end_it); // should never reach the end (sentinell)
						}
						if((*p_k_block_it).first == n_row_i) {
							const size_t n_row_i_height = r_block_cols_list[n_row_i].n_width;
							// an optimistic case, we found blocks at the same row

							CCholesky_OffDiagonalProduct<_TyGppContext, _CSelectedRowHeightsList/*CColumnWidthsList*/,
								CColumnWidth_k, CColumnWidth_j>::Loop(L_block_kj,
								(*p_k_block_it).second, (*p_j_block_it).second, n_row_i_height);
							// takes blocks from two different columns - need to know which columns have nonzero blocks
						} else {
							// next block in column k is on the next row,
							// we have to skip to the next block in jth row
						}
					}
					// this loop can be probably written in many different ways
					// todo - see if we could loop p_k_block_it instead and loose one parameter (p_j_block_end_it)
				}
				// cmod; causes fill-in in the current column

				-- p_k_block_end_it;
				_ASSERTE((*p_k_block_end_it).first == k); // makes sure that k-th column contains a diagonal block
				_TyDiagBlock_k d((*p_k_block_end_it).second);
				d.template triangularView<Eigen::Upper>().transpose().solveInPlace(L_block_kj); // modifies L_block_kj
				// d.marked<Eigen::UpperTriangular>().transpose().solveTriangularInPlace(L_block_kj); // the above line is deprecated, this line should do the trick in the new version(L_block_kj) of Eigen
				// late divide by the diagonal entries (probably unavoidable)
			}

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Performs the actual loop.
			 *
			 *	@param[out] r_column_dot is the running dor product of the column blocks
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, k)
			 *	@param[in] n_col_k_width is width of the k-th column, in elements
			 *		(equal to the number of rows of the dest block)
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_k_block_it is iterator to blocks in k-th column
			 *	@param[in] p_k_block_end_it is iterator pointing one past the last block in k-th column
			 *	@param[in] t_block_A is input block at A(j, k), same size as the dest block
			 *	@param[in] k is index of the previous column
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDiagMatrix &r_column_dot, double *p_dest_block,
				size_t UNUSED(n_col_k_width), _TyBlockConstIter p_j_block_it,
				_TyBlockConstIter p_j_block_end_it, _TyBlockConstIter p_k_block_it,
				_TyBlockConstIter p_k_block_end_it, TColumn::TBlockEntry t_block_A, size_t k,
				const std::vector<TColumn> &r_block_cols_list)
			{
				_ASSERTE(n_col_k_width == CColumnWidth_k::n_size);
				// make sure this is the correct block size

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					p_k_block_it != p_k_block_end_it);
				// if the first is a non-empty range, the second should be neither

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					/*r_col_L_k.block_list.back()*/(*(p_k_block_end_it - 1)).first == k);
				// the last block of k-th column is on diagonal (kth row) and it is always present
				// this block will serve as the sentinell for this loop

				_ASSERTE(p_j_block_it == p_j_block_end_it ||
					(*(p_j_block_end_it - 1)).first < k); // note that the p_j_block_end_it points at the diagonal block
				// the last block of j-th column (omitting the diagonal) is *at most* at kth row

				_TyDestBlock L_block_kj(p_dest_block);
				// the destination block

				if(t_block_A.first == k)
					L_block_kj = _TyDestBlock(t_block_A.second);
				else
					L_block_kj.setZero();
				// copy data from the source matrix

				{
					for(; p_j_block_it != p_j_block_end_it; ++ p_j_block_it) { // for each block in column j
						size_t n_row_i = (*p_j_block_it).first;
						_ASSERTE(n_row_i < k); // make sure that the sentinell is indeed functional
						_ASSERTE(p_k_block_it != p_k_block_end_it); // should not be pointing at the end in the first place (if we got inside this loop)
						while((*p_k_block_it).first < n_row_i) {
							++ p_k_block_it;
							_ASSERTE(p_k_block_it != p_k_block_end_it); // should never reach the end (sentinell)
						}
						if((*p_k_block_it).first == n_row_i) {
							const size_t n_row_i_height = r_block_cols_list[n_row_i].n_width;
							// an optimistic case, we found blocks at the same row

							CCholesky_OffDiagonalProduct<_TyGppContext, _CSelectedRowHeightsList/*CColumnWidthsList*/,
								CColumnWidth_k, CColumnWidth_j>::Loop(L_block_kj,
								(*p_k_block_it).second, (*p_j_block_it).second, n_row_i_height);
							// takes blocks from two different columns - need to know which columns have nonzero blocks
						} else {
							// next block in column k is on the next row,
							// we have to skip to the next block in jth row
						}
					}
					// this loop can be probably written in many different ways
					// todo - see if we could loop p_k_block_it instead and loose one parameter (p_j_block_end_it)
				}
				// cmod; causes fill-in in the current column

				-- p_k_block_end_it;
				_ASSERTE((*p_k_block_end_it).first == k); // makes sure that k-th column contains a diagonal block
				_TyDiagBlock_k d((*p_k_block_end_it).second);
				d.template triangularView<Eigen::Upper>().transpose().solveInPlace(L_block_kj); // modifies L_block_kj
				// d.marked<Eigen::UpperTriangular>().transpose().solveTriangularInPlace(L_block_kj); // the above line is deprecated, this line should do the trick in the new version(L_block_kj) of Eigen
				// late divide by the diagonal entries (probably unavoidable)

#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				//r_column_dot += L_block_kj.transpose().lazyProduct(L_block_kj);
				r_column_dot.template triangularView<Eigen::Upper>() += L_block_kj.transpose().lazyProduct(L_block_kj); // save FLOPs? unlikely with SSE; it *is* actually worth it!
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				//r_column_dot.noalias() += L_block_kj.transpose() * L_block_kj;
				r_column_dot.template triangularView<Eigen::Upper>().noalias() += L_block_kj.transpose() * L_block_kj; // save FLOPs? unlikely with SSE; it *is* actually worth it!
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
			}
		};

		/**
		 *	@brief off-diagonal product, one step of the Cholesky factorization kernel
		 *
		 *	Performs calculation of a product of two off-diagonal blocks in the same
		 *	row at two different columns of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CRowHeightList is a list of possible occupied row heights (given the column widths)
		 *	@tparam CColumnWidth_k is width of the k-th column
		 *	@tparam CColumnWidth_j is width of the j-th column
		 */
		template <class _TyGppContext, class _CRowHeightList, class CColumnWidth_k, class CColumnWidth_j>
		class CCholesky_OffDiagonalProduct {
		public:
			typedef typename _CRowHeightList::_TyHead CCurrentRowHeight; /**< @brief row height for this decision tree recursion */
			typedef typename CMakeMatrixRef<CColumnWidth_k::n_size, CColumnWidth_j::n_size>::_Ty _TyDestBlock; /**< @brief destination block shape */

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Wraps the product in a decision tree for k-th column width.
			 *
			 *	@param[in,out] L_block_kj is the output block at L(j, k)
			 *	@param[in] p_block_ik is the output block at L(i, k)
			 *	@param[in] p_block_ij is the output block at L(i, j)
			 *	@param[in] n_row_i_height is height of the i-th row, in elements
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDestBlock L_block_kj, const double *p_block_ik,
				const double *p_block_ij, size_t n_row_i_height)
			{
				if(n_row_i_height == CCurrentRowHeight::n_size) {
					CCholesky_OffDiagonalProduct<_TyGppContext, CTypelist<CCurrentRowHeight,
						CTypelistEnd>, CColumnWidth_k, CColumnWidth_j>::Loop(L_block_kj,
						p_block_ik, p_block_ij, n_row_i_height);
				} else {
					CCholesky_OffDiagonalProduct<_TyGppContext, typename
						_CRowHeightList::_TyTail, CColumnWidth_k, CColumnWidth_j>::Loop(
						L_block_kj, p_block_ik, p_block_ij, n_row_i_height);
				}
			}
		};

		/**
		 *	@brief off-diagonal product, one step of the Cholesky factorization kernel
		 *
		 *	Performs calculation of a product of two off-diagonal blocks in the same
		 *	row at two different columns of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CRowHeight_i is height of the i-th row
		 *	@tparam CColumnWidth_k is width of the k-th column
		 *	@tparam CColumnWidth_j is width of the j-th column
		 */
		template <class _TyGppContext, class CRowHeight_i, class CColumnWidth_k, class CColumnWidth_j>
		class CCholesky_OffDiagonalProduct<_TyGppContext, CTypelist<CRowHeight_i,
			CTypelistEnd>, CColumnWidth_k, CColumnWidth_j> {
		public:
			typedef typename CMakeMatrixRef<CColumnWidth_k::n_size, CColumnWidth_j::n_size>::_Ty _TyDestBlock; /**< @brief destination block shape */
			typedef typename CMakeMatrixRef<CRowHeight_i::n_size, CColumnWidth_j::n_size>::_Ty _TyBlock_ij; /**< @brief right-hand-side source block shape */
			typedef typename CMakeMatrixRef<CRowHeight_i::n_size, CColumnWidth_k::n_size>::_Ty _TyBlock_ik; /**< @brief left-hand-side source block shape */

			/**
			 *	@brief off-diagonal loop of the Cholesky factorization kernel
			 *
			 *	Performs the product.
			 *
			 *	@param[in,out] L_block_kj is the output block at L(j, k)
			 *	@param[in] p_block_ik is the output block at L(i, k)
			 *	@param[in] p_block_ij is the output block at L(i, j)
			 *	@param[in] n_row_i_height is height of the i-th row, in elements
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDestBlock L_block_kj, const double *p_block_ik,
				const double *p_block_ij, size_t UNUSED(n_row_i_height))
			{
				_ASSERTE(n_row_i_height == CRowHeight_i::n_size);

				_TyBlock_ik L_block_ik((double*)p_block_ik);
				_TyBlock_ij L_block_ij((double*)p_block_ij);

#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				L_block_kj -= L_block_ik.transpose().lazyProduct(L_block_ij);
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				L_block_kj.noalias() -= L_block_ik.transpose() * L_block_ij;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
			}
		};

		/**
		 *	@brief outer loop of the Cholesky factorization
		 *
		 *	@param[in] n_col_j_width is width of the j-th block column
		 *	@param[in] j is zero-based index the j-th block column
		 *	@param[in] n is number of block columns in the matrix
		 *	@param[in] r_block_cols_list is list of block columns in the factorized matrix
		 *	@param[in] r_col_A_j is reference to the j-th column in the original matrix
		 *	@param[in] r_lambda is reference to original matrix
		 *	@param[in] r_elim_tree is reference to the elimination tree
		 *	@param[in] ereach_stack is temporary workspace
		 *	@param[in] bitfield is temporary workspace
		 *	@param[in] alloc is reference to the block allocator for the factorized matrix
		 *
		 *	@return Returns true on success, false on failure (not positive definite).
		 *
		 *	@note This function throws std::bad_alloc.
		 */
		static __forceinline bool b_Column_Loop(const size_t n_col_j_width, const size_t j, const size_t n,
			std::vector<TColumn> &r_block_cols_list, const TColumn &r_col_A_j,
			const CUberBlockMatrix &r_lambda, const std::vector<size_t> &r_elim_tree,
			std::vector<size_t> &ereach_stack, std::vector<size_t> &bitfield, _TyDenseAllocator &alloc) // throw(std::bad_alloc) // todo - doc, throws, ...
		{
			return CCholesky_ColumnLoop<CBlockMatrixTypelist,
				CColumnWidthsList>::b_Loop(n_col_j_width, j, n, r_block_cols_list, r_col_A_j,
				r_lambda, r_elim_tree, ereach_stack, bitfield, alloc);
		}

		/**
		 *	@brief column loop of cholesky factorization
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthList is list of possible column widths
		 */
		template <class _TyGppContext, class _CColumnWidthList>
		class CCholesky_ColumnLoop {
		public:
			typedef typename _CColumnWidthList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief outer loop of the Cholesky factorization
			 *
			 *	Wraps the body of the loop in block column size decision tree.
			 *
			 *	@param[in] n_col_j_width is width of the j-th block column
			 *	@param[in] j is zero-based index the j-th block column
			 *	@param[in] n is number of block columns in the matrix
			 *	@param[in] r_block_cols_list is list of block columns in the factorized matrix
			 *	@param[in] r_col_A_j is reference to the j-th column in the original matrix
			 *	@param[in] r_lambda is reference to original matrix
			 *	@param[in] r_elim_tree is reference to the elimination tree
			 *	@param[in] ereach_stack is temporary workspace
			 *	@param[in] bitfield is temporary workspace
			 *	@param[in] alloc is reference to the block allocator for the factorized matrix
			 *
			 *	@return Returns true on success, false on failure (not positive definite).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool b_Loop(const size_t n_col_j_width, const size_t j, const size_t n,
				std::vector<TColumn> &r_block_cols_list, const TColumn &r_col_A_j,
				const CUberBlockMatrix &r_lambda, const std::vector<size_t> &r_elim_tree,
				std::vector<size_t> &ereach_stack, std::vector<size_t> &bitfield, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				if(n_col_j_width == CCurrentColumnWidth::n_size) {
					return CCholesky_ColumnLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::b_Loop(n_col_j_width, j, n, r_block_cols_list, r_col_A_j,
						r_lambda, r_elim_tree, ereach_stack, bitfield, alloc);
				} else {
					return CCholesky_ColumnLoop<_TyGppContext, typename
						_CColumnWidthList::_TyTail>::b_Loop(n_col_j_width, j, n, r_block_cols_list,
						r_col_A_j, r_lambda, r_elim_tree, ereach_stack, bitfield, alloc);
				}
			}
		};

		/**
		 *	@brief column loop of cholesky factorization (specialization for the chosen column width)
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CColumnWidth_j is width of the current column
		 */
		template <class _TyGppContext, class CColumnWidth_j>
		class CCholesky_ColumnLoop<_TyGppContext, CTypelist<CColumnWidth_j, CTypelistEnd> > {
		public:
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_Ty _TyDiagBlock; /**< @brief reference to the diagonal block */
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_TyMatrix _TyDiagMatrix; /**< @brief a dense matrix with the same shape as the diagonal block */

			/**
			 *	@brief outer loop of the Cholesky factorization
			 *
			 *	Contains a specialized instance the loop body for the given column width.
			 *
			 *	@param[in] n_col_j_width is width of the j-th block column (unused)
			 *	@param[in] j is zero-based index the j-th block column
			 *	@param[in] n is number of block columns in the matrix
			 *	@param[in] r_block_cols_list is list of block columns in the factorized matrix
			 *	@param[in] r_col_A_j is reference to the j-th column in the original matrix
			 *	@param[in] r_lambda is reference to original matrix
			 *	@param[in] r_elim_tree is reference to the elimination tree
			 *	@param[in] ereach_stack is temporary workspace
			 *	@param[in] bitfield is temporary workspace
			 *	@param[in] alloc is reference to the block allocator for the factorized matrix
			 *
			 *	@return Returns true on success, false on failure (not positive definite).
			 *
			 *	@note This function throws std::bad_alloc.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool b_Loop(const size_t UNUSED(n_col_j_width), const size_t j, const size_t n,
				std::vector<TColumn> &r_block_cols_list, const TColumn &r_col_A_j,
				const CUberBlockMatrix &r_lambda, const std::vector<size_t> &r_elim_tree,
				std::vector<size_t> &ereach_stack, std::vector<size_t> &bitfield, _TyDenseAllocator &alloc) // throw(std::bad_alloc)
			{
				_ASSERTE(n_col_j_width == CColumnWidth_j::n_size);

				TColumn &r_col_L_j = r_block_cols_list[j];
				r_col_L_j.n_width = CColumnWidth_j::n_size;
				r_col_L_j.n_cumulative_width_sum = r_col_A_j.n_cumulative_width_sum;
				// get & copy columns

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

				_TyDiagMatrix column_dot;
				column_dot.setZero();
				// accumulate column dot product on the stack

				for(size_t u = n_ereach_first, n_highest_k = 0; u < n; ++ u) { // seems to work rather nicely (strange because not every A_up[k, j] is not accessed then - it is likely null) // todo - verify this
					const size_t k = ereach_stack[u]; // use ereach to predict which columns will have nonzero products
					// k generally rises, but it doesn't have to (can be non-monotonic)

					_ASSERTE(k != n_highest_k || u == n_ereach_first); // there should be no column repeated (except for zero, in the first iteration)
					bool b_ereach_mono;
					if((b_ereach_mono = (k >= n_highest_k)))
						n_highest_k = k; // don't remember previous k unless it increased
					// see if the ereach is (locally) increasing

					const TColumn &r_col_L_k = r_block_cols_list[k];
					const size_t n_col_k_width = r_col_L_k.n_width;
					// get column k

					_TyBlockConstIter p_jk_block_it;
					double *p_k_block_data = alloc.p_Get_DenseStorage(n_col_k_width * CColumnWidth_j::n_size);
					//_TyMatrixXdRef L_block_kj(p_k_block_data, n_col_k_width, CColumnWidth_j::n_size); // only accessed in FBS section
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
							std::lower_bound(r_col_L_j.block_list.begin(), r_col_L_j.block_list.end(), k, CUberBlockMatrix_Base::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
							std::lower_bound(r_col_L_j.block_list.begin(), r_col_L_j.block_list.end(), k, CUberBlockMatrix_Base::CompareBlockRow);
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
								std::lower_bound(p_A_block_it, p_A_block_end_it, k, CUberBlockMatrix_Base::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
								std::lower_bound(p_A_block_it, p_A_block_end_it, k, CUberBlockMatrix_Base::CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
						} else if(n_block_row_id > k) {
							p_A_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
								std::lower_bound(r_col_A_j.block_list.begin(), p_A_block_end_it, k, CUberBlockMatrix_Base::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
								std::lower_bound(r_col_A_j.block_list.begin(), p_A_block_end_it, k, CUberBlockMatrix_Base::CompareBlockRow);
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

					/*CFBS_Cholesky<CBlockMatrixTypelist>::OffDiagonal_Loop(p_k_block_data,
						CColumnWidth_j::n_size, n_col_k_width, p_j_block_it, p_j_block_end_it,
						p_k_block_it, p_k_block_end_it, *p_A_block_it, k, r_block_cols_list);*/ // this is how it's called from BlockMatrix.h
					CCholesky_OffDiagonalLoop<CBlockMatrixTypelist, CTypelist<CColumnWidth_j,
						CTypelistEnd> >::Loop(column_dot, p_k_block_data, n_col_k_width,
						p_j_block_it, p_j_block_end_it, p_k_block_it, p_k_block_end_it,
						*p_A_block_it, k, r_block_cols_list); // save one decision tree level, calculate column_dot
					// execute cmod and solve by diagonals using FBS

					if((*p_A_block_it).first == k)
						++ p_A_block_it;
					// skip to the next one for the next iteration / for the diagonal
				}
				// complex data dependencies in here, not sure if i like it

				bool b_result;
				{
					_ASSERTE(p_A_block_it != p_A_block_end_it && (*p_A_block_it).first <= j);
					// it is pointing before or at the diagonal already
					if((*p_A_block_it).first != j) {
						p_A_block_it =
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
							std::lower_bound(p_A_block_it, p_A_block_end_it, j, CUberBlockMatrix_Base::CCompareBlockRow());
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
							std::lower_bound(p_A_block_it, p_A_block_end_it, j, CUberBlockMatrix_Base::CompareBlockRow);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
					}
					_ASSERTE(p_A_block_it != p_A_block_end_it && (*p_A_block_it).first == j); // should always be nonzero
					// find the diagonal block (note that if A is upper diagonal, it is the
					// last block and we can optimize for that)

					_ASSERTE(r_col_L_j.block_list.empty() || r_col_L_j.block_list.back().first < j);
					// makes sure that the column doesn't contain the diagonal block yet

					double *p_diag_data = alloc.p_Get_DenseStorage(CColumnWidth_j::n_size * CColumnWidth_j::n_size);
					r_col_L_j.block_list.push_back(TColumn::TBlockEntry(j, p_diag_data));
					// allocates a new block in this column

					/*b_result = CFBS_Cholesky<CBlockMatrixTypelist>::b_Diagonal_Loop(p_diag_data,
						CColumnWidth_j::n_size, r_col_L_j.block_list.begin(), r_col_L_j.block_list.end() - 1,
						(*p_A_block_it).second, r_block_cols_list);*/ // this is how it's called from BlockMatrix.h
					/*b_result = CCholesky_DiagonalLoop<CBlockMatrixTypelist, CTypelist<CColumnWidth_j,
						CTypelistEnd> >::b_Loop(p_diag_data, CColumnWidth_j::n_size,
						r_col_L_j.block_list.begin(), r_col_L_j.block_list.end() - 1,
						(*p_A_block_it).second, r_block_cols_list);*/ // this calulates column_dot inside, requires another decision tree
					{
						_TyDiagBlock L_block_jj(p_diag_data);
						L_block_jj = _TyDiagBlock((double*)(*p_A_block_it).second);
						// get the dest block and initialize it with values from A

						L_block_jj -= column_dot;
						// calculated while calculating cmods

						Eigen::LLT<_TyDiagMatrix, Eigen::Upper> chol(L_block_jj); // Eigen::LLT only accesses a half of the matrix (upper tri in this case), no need to clear the lower half
						b_result = (chol.info() == Eigen::Success);
						L_block_jj = chol.matrixU(); // t_odo - make sure it will not throw if chol.info() != Eigen::Success // does not, in 3.1.3
						// calculates cholesky of a square block
					}
					// execute cdiv using FBS
				}
				// cdiv; reads an entire column and produces diagonal

				_ASSERTE(r_col_L_j.block_list.size() == n - n_ereach_first/*n_ereach_size*/ + 1);
				// make sure we preallocated it correclty

				return b_result;
			}
		};

		/**
		 *	@brief diagonal loop of the Cholesky factorization kernel
		 *
		 *	Performs calculation of the diagonal element in one column of Cholesky factor.
		 *
		 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, j), which is square
		 *	@param[in] n_col_j_width is width of the j-th column, in elements
		 *	@param[in] p_j_block_it is iterator to blocks in j-th column
		 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
		 *	@param[in] p_block_A is input block at A(j, j), same size as the dest block
		 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
		 *
		 *	@return Returns true on success, false on failure (the matrix is not positive definite).
		 */
		static __forceinline bool b_Diagonal_Loop(double *p_dest_block, size_t n_col_j_width,
			_TyBlockConstIter p_j_block_it, _TyBlockConstIter p_j_block_end_it,
			const double *p_block_A, const std::vector<TColumn> &r_block_cols_list)
		{
			return CCholesky_DiagonalLoop<CBlockMatrixTypelist,
				CColumnWidthsList>::b_Loop(p_dest_block, n_col_j_width,
				p_j_block_it, p_j_block_end_it, p_block_A, r_block_cols_list);
		}

		/**
		 *	@brief diagonal loop of the Cholesky factorization kernel
		 *
		 *	Performs calculation of the diagonal element in one column of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthList is list of possible column widths
		 */
		template <class _TyGppContext, class _CColumnWidthList>
		class CCholesky_DiagonalLoop {
		public:
			typedef typename _CColumnWidthList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief diagonal loop of the Cholesky factorization kernel
			 *
			 *	Wraps the diagonal loop in column width decision tree.
			 *
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, j), which is square
			 *	@param[in] n_col_j_width is width of the j-th column, in elements
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_block_A is input block at A(j, j), same size as the dest block
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 *
			 *	@return Returns true on success, false on failure (the matrix is not positive definite).
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool b_Loop(double *p_dest_block, size_t n_col_j_width,
				_TyBlockConstIter p_j_block_it, _TyBlockConstIter p_j_block_end_it,
				const double *p_block_A, const std::vector<TColumn> &r_block_cols_list)
			{
				if(n_col_j_width == CCurrentColumnWidth::n_size) {
					return CCholesky_DiagonalLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::b_Loop(p_dest_block, n_col_j_width, p_j_block_it,
						p_j_block_end_it, p_block_A, r_block_cols_list);
				} else {
					return CCholesky_DiagonalLoop<_TyGppContext, typename
						_CColumnWidthList::_TyTail>::b_Loop(p_dest_block, n_col_j_width,
						p_j_block_it, p_j_block_end_it, p_block_A, r_block_cols_list);
				}
			}
		};

		/**
		 *	@brief diagonal loop of the Cholesky factorization kernel
		 *
		 *	Performs calculation of the diagonal element in one column of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CColumnWidth_j is width of the j-th column
		 */
		template <class _TyGppContext, class CColumnWidth_j>
		class CCholesky_DiagonalLoop<_TyGppContext, CTypelist<CColumnWidth_j, CTypelistEnd> > {
		public:
			typedef typename CFilterTypelist2<CColumnWidthsList, CColumnWidth_j,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_Uniq>::_TyResult _CSelectedRowHeightsList_on_j; /**< @brief list of row heights in column j */
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_Ty _TyDiagBlock; /**< @brief reference to the diagonal block */
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_TyMatrix _TyDiagMatrix; /**< @brief a dense matrix with the same shape as the diagonal block */

			/**
			 *	@brief diagonal loop of the Cholesky factorization kernel
			 *
			 *	Performs calculation of the diagonal element in one column of Cholesky factor.
			 *
			 *	@param[out] p_dest_block is the (uninitialized) output block at L(j, j), which is square
			 *	@param[in] n_col_j_width is width of the j-th column, in elements
			 *	@param[in] p_j_block_it is iterator to blocks in j-th column
			 *	@param[in] p_j_block_end_it is iterator pointing one past the last block in j-th column
			 *	@param[in] p_block_A is input block at A(j, j), same size as the dest block
			 *	@param[in] r_block_cols_list is block column layout of the (symmetric) matrix
			 *
			 *	@return Returns true on success, false on failure (the matrix is not positive definite).
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool b_Loop(double *p_dest_block, size_t UNUSED(n_col_j_width),
				_TyBlockConstIter p_j_block_it, _TyBlockConstIter p_j_block_end_it,
				const double *p_block_A, const std::vector<TColumn> &r_block_cols_list)
			{
				_ASSERTE(n_col_j_width == CColumnWidth_j::n_size);

				_TyDiagBlock L_block_jj(p_dest_block);
				L_block_jj = _TyDiagBlock((double*)p_block_A);
				// get the dest block and initialize it with values from A

				for(; p_j_block_it != p_j_block_end_it; ++ p_j_block_it) {
					TColumn::TBlockEntry t_L_ij = *p_j_block_it;
					size_t n_row_i_height = r_block_cols_list[t_L_ij.first].n_width;
					CCholesky_DiagonalProduct<_TyGppContext, _CSelectedRowHeightsList_on_j/*CColumnWidthsList*/,
						CColumnWidth_j>::Loop(L_block_jj, t_L_ij.second, n_row_i_height);
				}
				// use sparse loop instead
				// calculates square of all the blocks in the current column (has shape
				// column-width * column-width, and the diagonal block has the same shape)

				Eigen::LLT<_TyDiagMatrix, Eigen::Upper> chol(L_block_jj); // Eigen::LLT only accesses a half of the matrix (upper tri in this case), no need to clear the lower half
				if(chol.info() != Eigen::Success)
					return false;
				L_block_jj = chol.matrixU();
				// calculates cholesky of a square block

				return true;
			}
		};

		/**
		 *	@brief diagonal product from the Cholesky factorization kernel
		 *
		 *	Performs calculation of the diagonal element in one column of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CRowHeightList is a list of possible row heights
		 *	@tparam CColumnWidth_j is width of the j-th column
		 */
		template <class _TyGppContext, class _CRowHeightList, class CColumnWidth_j>
		class CCholesky_DiagonalProduct {
		public:
			typedef typename _CRowHeightList::_TyHead CCurrentRowHeight; /**< @brief row height for this decision tree recursion */
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_Ty _TyDestBlock; /**< @brief destination block shape */

			/**
			 *	@brief diagonal product from the Cholesky factorization kernel
			 *
			 *	Wraps the product in a row height decision tree.
			 *
			 *	@param[in,out] L_block_jj is the output block at L(j, j)
			 *	@param[in] p_block_ij is the output block at L(i, j)
			 *	@param[in] n_row_i_height is height of the i-th row, in elements
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDestBlock L_block_jj,
				const double *p_block_ij, size_t n_row_i_height)
			{
				if(n_row_i_height == CCurrentRowHeight::n_size) {
					CCholesky_DiagonalProduct<_TyGppContext, CTypelist<CCurrentRowHeight,
						CTypelistEnd>, CColumnWidth_j>::Loop(L_block_jj,
						p_block_ij, n_row_i_height);
				} else {
					CCholesky_DiagonalProduct<_TyGppContext, typename
						_CRowHeightList::_TyTail, CColumnWidth_j>::Loop(
						L_block_jj, p_block_ij, n_row_i_height);
				}
			}
		};

		/**
		 *	@brief diagonal product from the Cholesky factorization kernel
		 *
		 *	Performs calculation of the diagonal element in one column of Cholesky factor.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CRowHeight_i is height of the i-th row
		 *	@tparam CColumnWidth_j is width of the j-th column
		 */
		template <class _TyGppContext, class CRowHeight_i, class CColumnWidth_j>
		class CCholesky_DiagonalProduct<_TyGppContext, CTypelist<CRowHeight_i,
			CTypelistEnd>, CColumnWidth_j> {
		public:
			typedef typename CMakeMatrixRef<CColumnWidth_j::n_size, CColumnWidth_j::n_size>::_Ty _TyDestBlock; /**< @brief destination block shape */
			typedef typename CMakeMatrixRef<CRowHeight_i::n_size, CColumnWidth_j::n_size>::_Ty _TyBlock_ij; /**< @brief source block shape */

			/**
			 *	@brief diagonal product from the Cholesky factorization kernel
			 *
			 *	Performs the product.
			 *
			 *	@param[in,out] L_block_jj is the output block at L(j, j)
			 *	@param[in] p_block_ij is the output block at L(i, j)
			 *	@param[in] n_row_i_height is height of the i-th row, in elements
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Loop(_TyDestBlock L_block_jj,
				const double *p_block_ij, size_t UNUSED(n_row_i_height))
			{
				_ASSERTE(n_row_i_height == CRowHeight_i::n_size);

				_TyBlock_ij L_block_ij((double*)p_block_ij);

#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				L_block_jj -= L_block_ij.transpose().lazyProduct(L_block_ij);
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				L_block_jj.noalias() -= L_block_ij.transpose() * L_block_ij;
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
			}
		};
	};

	/**
	 *	@brief fixed block size operations template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 */
	template <class CBlockMatrixTypelist>
	class CFBS_TriangularSolve : public fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist> {
	public:
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CRowHeightsList CRowHeightsList; /**< @brief list of unique block row heights */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CColumnWidthsList CColumnWidthsList; /**< @brief list of unique block column widths */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CDimsList_Uniq CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */

		/**
		 *	@brief inner loop of the forward / back substitution kernels,
		 *		calculating solution of a system with triangular matrix
		 *
		 *	Performs forward/back substitution on non-diagonal blocks in the current column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CRowHeightList is list of possible row heights
		 *	@tparam CColumnWidth is selected column width
		 */
		template <class _TyGppContext, class _CRowHeightList, class CColumnWidth>
		class CTriangularSolve_InnerLoop {
		public:
			typedef typename CMakeVectorRef<CColumnWidth::n_size>::_Ty _TyVector_ColPart; /**< @brief part of the vector, corresponding to this column */
			//typedef typename CMakeRowVectorRef<CColumnWidth::n_size>::_Ty _TyRowVector_ColPart; /**< @brief part of the vector, corresponding to this column */
			typedef typename _CRowHeightList::_TyHead CCurrentRowHeight; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief wraps forward-substitution loop body in row height decision tree
			 *
			 *	@param[in,out] r_diag_x is the part of the left-hand side vector, correspoinding to the current column
			 *	@param[in] p_x_row is the part of the left-hand side vector, correspoinding to the current row
			 *	@param[in] p_block_data is pointer to the current block data
			 *	@param[in] n_row_height is the current row height, in elements
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void ForwardLoop(_TyVector_ColPart &r_diag_x, const double *p_x_row,
				const double *p_block_data, size_t n_row_height)
			{
				if(n_row_height == CCurrentRowHeight::n_size) {
					CTriangularSolve_InnerLoop<_TyGppContext, CTypelist<CCurrentRowHeight,
						CTypelistEnd>, CColumnWidth>::ForwardLoop(r_diag_x, p_x_row, p_block_data, n_row_height);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CTriangularSolve_InnerLoop<_TyGppContext, typename _CRowHeightList::_TyTail,
						CColumnWidth>::ForwardLoop(r_diag_x, p_x_row, p_block_data, n_row_height);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

			/**
			 *	@brief wraps back-substitution loop body in row height decision tree
			 *
			 *	@param[in] r_diag_x is the part of the left-hand side vector, correspoinding to the current column
			 *	@param[in,out] p_x_row is the part of the left-hand side vector, correspoinding to the current row
			 *	@param[in] p_block_data is pointer to the current block data
			 *	@param[in] n_row_height is the current row height, in elements
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void BackLoop(const _TyVector_ColPart &r_diag_x, double *p_x_row,
				const double *p_block_data, size_t n_row_height)
			{
				if(n_row_height == CCurrentRowHeight::n_size) {
					CTriangularSolve_InnerLoop<_TyGppContext, CTypelist<CCurrentRowHeight,
						CTypelistEnd>, CColumnWidth>::BackLoop(r_diag_x, p_x_row, p_block_data, n_row_height);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CTriangularSolve_InnerLoop<_TyGppContext, typename _CRowHeightList::_TyTail,
						CColumnWidth>::BackLoop(r_diag_x, p_x_row, p_block_data, n_row_height);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(n_row_height == %d) {\n", CCurrentRowHeight::n_size);
				CTriangularSolve_InnerLoop<_TyGppContext, CTypelist<CCurrentRowHeight,
					CTypelistEnd>, CColumnWidth>::Debug();
				printf("} else {\n");
				CTriangularSolve_InnerLoop<_TyGppContext,
					typename _CRowHeightList::_TyTail, CColumnWidth>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief inner loop of the forward / back substitution kernels,
		 *		calculating solution of a system with triangular matrix
		 *
		 *	Performs forward/back substitution on non-diagonal blocks in the current column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CRowHeight is selected row height
		 *	@tparam CColumnWidth is selected column width
		 */
		template <class _TyGppContext, class CRowHeight, class CColumnWidth>
		class CTriangularSolve_InnerLoop<_TyGppContext, CTypelist<CRowHeight, CTypelistEnd>, CColumnWidth> {
		public:
			typedef typename CMakeRowVectorRef<CRowHeight::n_size>::_Ty _TyRowVector_RowPart; /**< @brief part of the vector, corresponding to this column as a row vector */
			typedef typename CMakeVectorRef<CRowHeight::n_size>::_Ty _TyVector_RowPart; /**< @brief part of the vector, corresponding to this column as a column vector */
			typedef typename CMakeMatrixRef<CRowHeight::n_size, CColumnWidth::n_size>::_Ty _TyBlock; /**< @brief off-diagonal matrix block Eigen::Map type */
			//typedef typename CMakeRowVectorRef<CColumnWidth::n_size>::_Ty _TyRowVector_ColPart; /**< @brief part of the vector, corresponding to this column */
			typedef typename CMakeVectorRef<CColumnWidth::n_size>::_Ty _TyVector_ColPart; /**< @brief part of the vector, corresponding to this column */

			/**
			 *	@brief performs forward-substitution on a single matrix block
			 *
			 *	@param[in,out] r_diag_x is the part of the left-hand side vector, correspoinding to the current column
			 *	@param[in] p_x_row is the part of the left-hand side vector, correspoinding to the current row
			 *	@param[in] p_block_data is pointer to the current block data
			 *	@param[in] n_row_height is the current row height, in elements
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void ForwardLoop(_TyVector_ColPart &r_diag_x, const double *p_x_row,
				const double *p_block_data, size_t UNUSED(n_row_height))
			{
				_ASSERTE(n_row_height == CRowHeight::n_size);
				// make sure that's really it

				_TyRowVector_RowPart x_part((double*)p_x_row); // also column vector, also unaligned
				_TyBlock block((double*)p_block_data); // block

#if 0 // t_odo - clear this
				for(size_t i = 0; i < CColumnWidth::n_size; ++ i) {
					for(size_t j = 0; j < CRowHeight::n_size; ++ j)
						r_diag_x(i) -= x_part(j) * block(j, i);
				}
#else // 0
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				r_diag_x -= x_part.lazyProduct(block); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				r_diag_x.noalias() -= x_part * block; // do the same as the loop above, using matrix multiplication
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
#endif // 0
				// calculate x -= column(j) * x(j) in parallel for all columns in this block
			}

			/**
			 *	@brief performs back-substitution on a single matrix block
			 *
			 *	@param[in] r_diag_x is the part of the left-hand side vector, correspoinding to the current column
			 *	@param[in,out] p_x_row is the part of the left-hand side vector, correspoinding to the current row
			 *	@param[in] p_block_data is pointer to the current block data
			 *	@param[in] n_row_height is the current row height, in elements
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void BackLoop(const _TyVector_ColPart &r_diag_x, double *p_x_row,
				const double *p_block_data, size_t UNUSED(n_row_height))
			{
				_ASSERTE(n_row_height == CRowHeight::n_size);
				// make sure that's really it

				_TyRowVector_RowPart x_part(p_x_row); // also column vector, also unaligned
				_TyBlock block((double*)p_block_data); // block

#if 0 // t_odo - clear this
				for(size_t i = CColumnWidth::n_size; i > 0;) { // even reversing loop direction gives numerical discrepancies (big dynamic range in the matrix / vector)
					-- i; // !!
					for(size_t j = 0; j < CRowHeight::n_size; ++ j)
						x_part(j) -= block(j, i) * r_diag_x(i);
				}
#else // 0
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				x_part -= block.lazyProduct(r_diag_x); // fbsla
#else // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
				x_part.noalias() -= block * r_diag_x; // do the same as the loop above, using matrix multiplication
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
#endif // 0
				// calculate x -= column(j) * x(j) in parallel for all columns in this block
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(n_row_height == %d);\n", CRowHeight::n_size);
				printf("r_diag_x.noalias() -= x_part.transpose() * block; // if forward substution\n");
				printf("x_part.noalias() -= block * r_diag_x; // if back substution\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief outer loop of the forward / back substitution kernels,
		 *		calculating solution of a system with triangular matrix
		 *
		 *	Iterates through all the blocks in the column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthList is list of possible column widths
		 */
		template <class _TyGppContext, class _CColumnWidthList>
		class CTriangularSolve_OuterLoop {
		public:
			typedef typename _CColumnWidthList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief wraps forward-substitution loop body in block size decision tree
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in,out] p_x is the left-hand side vector (becomes right-hand side as the system is solved)
			 *	@param[in] r_row_list is list of matrix rows
			 *
			 *	@return Returns true on success, false in case the system doesn't have exactly one solution.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool ForwardLoop(const TColumn &r_t_column,
				double *p_x, const std::vector<TRow> &r_row_list)
			{
				if(r_t_column.n_width == CCurrentColumnWidth::n_size) {
					return CTriangularSolve_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::ForwardLoop(r_t_column, p_x, r_row_list);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					return CTriangularSolve_OuterLoop<_TyGppContext, typename
						_CColumnWidthList::_TyTail>::ForwardLoop(r_t_column, p_x, r_row_list);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

			/**
			 *	@brief wraps back-substitution loop body in block size decision tree
			 *
			 *	@param[in] r_t_column is the current column
			 *	@param[in,out] p_x is the left-hand side vector (becomes right-hand side as the system is solved)
			 *	@param[in] r_row_list is list of matrix rows
			 *
			 *	@return Returns true on success, false in case the system doesn't have exactly one solution.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool BackLoop(const TColumn &r_t_column,
				double *p_x, const std::vector<TRow> &r_row_list)
			{
				if(r_t_column.n_width == CCurrentColumnWidth::n_size) {
					return CTriangularSolve_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd> >::BackLoop(r_t_column, p_x, r_row_list);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					return CTriangularSolve_OuterLoop<_TyGppContext, typename
						_CColumnWidthList::_TyTail>::BackLoop(r_t_column, p_x, r_row_list);
					// not CCurrentColumnWidth::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("if(r_t_column.n_width == %d) {\n", CCurrentColumnWidth::n_size);
				CTriangularSolve_OuterLoop<_TyGppContext, CTypelist<CCurrentColumnWidth, CTypelistEnd> >::Debug();
				printf("} else {\n");
				CTriangularSolve_OuterLoop<_TyGppContext, typename _CColumnWidthList::_TyTail>::Debug();
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief outer loop of the forward / back substitution kernels,
		 *		calculating solution of a system with triangular matrix
		 *
		 *	Iterates through all the blocks in the column.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthList is list of possible column widths
		 */
		template <class _TyGppContext, class CColumnWidth>
		class CTriangularSolve_OuterLoop<_TyGppContext, CTypelist<CColumnWidth, CTypelistEnd> > {
		public:
			typedef typename CMakeRowVectorRef<CColumnWidth::n_size>::_Ty _TyRowVector_ColPart; /**< @brief part of the vector, corresponding to this column */
			typedef typename CMakeVectorRef<CColumnWidth::n_size>::_Ty _TyVector_ColPart; /**< @brief part of the vector, corresponding to this column */
			typedef typename CMakeMatrixRef<CColumnWidth::n_size, CColumnWidth::n_size>::_Ty _TyDiagBlock; /**< @brief diagonal matrix block Eigen::Map type */
			typedef typename CFilterTypelist2<CRowHeightsList, CColumnWidth,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_Uniq>::_TyResult _CSelectedRowHeightsList; /**< @brief row heights list, filtered by the selected column width */

			/**
			 *	@brief performs forward-substitution on a single column
			 *
			 *	@param[in] r_t_col is the current column
			 *	@param[in,out] p_x is the left-hand side vector (becomes right-hand side as the system is solved)
			 *	@param[in] r_row_list is list of matrix rows
			 *
			 *	@return Returns true on success, false in case the system doesn't have exactly one solution.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool ForwardLoop(const TColumn &r_t_col,
				double *p_x, const std::vector<TRow> &r_row_list)
			{
				_ASSERTE(r_t_col.n_width == CColumnWidth::n_size);
				// make sure this is the size

				const TColumn::TBlockEntry &r_t_diag_block = r_t_col.block_list.back();
				const TRow &r_t_diag_row = r_row_list[r_t_diag_block.first];
				if(CColumnWidth::n_size != r_t_diag_row.n_height ||
				   r_t_col.n_cumulative_width_sum != r_t_diag_row.n_cumulative_height_sum)
					return false;
				const size_t n_x = r_t_col.n_cumulative_width_sum;
				// make sure that the last block is on the diagonal, and is square (that does happen in L)

				_TyVector_ColPart diag_x(p_x + n_x); // note this is always unaligned
				// part of the x vector, corresponding to this column

				for(_TyBlockConstIter p_block_it =
				   r_t_col.block_list.begin(), p_block_end_it = r_t_col.block_list.end() - 1;
				   p_block_it != p_block_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					const TRow &r_t_row = r_row_list[r_t_block.first];
					const size_t n_height = r_t_row.n_height, n_y = r_t_row.n_cumulative_height_sum;
					_ASSERTE(n_y + n_height <= n_x); // make sure that the diagonal part of the vector is not affected here

					CTriangularSolve_InnerLoop<_TyGppContext, _CSelectedRowHeightsList,
						CColumnWidth>::ForwardLoop(diag_x, p_x + n_y, r_t_block.second, n_height);
					// call the inner loop
				}
				// for every other block in this column (note no reverse order here)
				// note that this loop mostly executes once (optimize for that? how?)

				_TyDiagBlock diag_block(r_t_diag_block.second);
				for(size_t i = 0; i < CColumnWidth::n_size; ++ i) {
					double f_old_x = diag_x(i);
					for(size_t j = 0; j < i; ++ j) // elements strictly above diagonal
						f_old_x -= diag_block(j, i) * diag_x(j);
					// forward substitute the diagonal block

#ifdef _DEBUG
					for(size_t j = i + 1; j < CColumnWidth::n_size; ++ j) // elements strictly under diagonal
						_ASSERTE(diag_block(j, i) == 0);
					// make sure there are only nulls under the diagonal
#endif // _DEBUG

					double f_diag = diag_block(i, i);
#ifdef __UBER_BLOCK_MATRIX_TRIANGULAR_SOLVE_NEVER_FAILS
					_ASSERTE(f_diag != 0);
#else // __UBER_BLOCK_MATRIX_TRIANGULAR_SOLVE_NEVER_FAILS
					if(f_diag == 0) // what about negative zero? does it cover it? yes.
						return false; // note that if diag_x(i) is zero as well, it probably means "infinite number of solutions", and not "no solution". could still produce one of them.
#endif // __UBER_BLOCK_MATRIX_TRIANGULAR_SOLVE_NEVER_FAILS
					diag_x(i) = f_old_x / f_diag;
				}
				// resolve values of the diagonal x

				return true;
			}

			/**
			 *	@brief performs back-substitution on a single column
			 *
			 *	@param[in] r_t_col is the current column
			 *	@param[in,out] p_x is the left-hand side vector (becomes right-hand side as the system is solved)
			 *	@param[in] r_row_list is list of matrix rows
			 *
			 *	@return Returns true on success, false in case the system doesn't have exactly one solution.
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline bool BackLoop(const TColumn &r_t_col,
				double *p_x, const std::vector<TRow> &r_row_list)
			{
				const TColumn::TBlockEntry &r_t_diag_block = r_t_col.block_list.back();
				const TRow &r_t_diag_row = r_row_list[r_t_diag_block.first];
				if(CColumnWidth::n_size != r_t_diag_row.n_height ||
				   r_t_col.n_cumulative_width_sum != r_t_diag_row.n_cumulative_height_sum)
					return false;
				const size_t n_x = r_t_col.n_cumulative_width_sum;
				// make sure that the last block is on the diagonal, and is square (that does happen in L)

				_TyVector_ColPart diag_x(p_x + n_x); // note this is always unaligned
				_TyDiagBlock diag_block(r_t_diag_block.second);
				for(size_t i = CColumnWidth::n_size; i > 0;) {
					-- i; // !!

					double f_diag = diag_block(i, i);
#ifdef __UBER_BLOCK_MATRIX_TRIANGULAR_SOLVE_NEVER_FAILS
					_ASSERTE(f_diag != 0);
#else // __UBER_BLOCK_MATRIX_TRIANGULAR_SOLVE_NEVER_FAILS
					if(f_diag == 0) // what about negative zero? does it cover it? yes.
						return false; // note that if diag_x(i) is zero as well, it probably means "infinite number of solutions", and not "no solution". could still produce one of them.
#endif // __UBER_BLOCK_MATRIX_TRIANGULAR_SOLVE_NEVER_FAILS
					double f_new_x = (diag_x(i) /= f_diag); // f_new_x = diag_x(i)
					for(size_t j = 0; j < i; ++ j) // elements strictly above diagonal
						diag_x(j) -= diag_block(j, i) * f_new_x;
					// backsubstitute the diagonal block

#ifdef _DEBUG
					for(size_t j = i + 1; j < CColumnWidth::n_size; ++ j) // elements strictly under diagonal
						_ASSERTE(diag_block(j, i) == 0);
					// make sure there are only nulls under the diagonal
#endif // _DEBUG
				}
				// resolve values of the diagonal x

				for(_TyBlockConstIter p_block_it =
				   r_t_col.block_list.begin(), p_block_end_it = r_t_col.block_list.end() - 1;
				   p_block_it != p_block_end_it; ++ p_block_it) {
					const TColumn::TBlockEntry &r_t_block = *p_block_it;
					const TRow &r_t_row = r_row_list[r_t_block.first];
					const size_t n_height = r_t_row.n_height, n_y = r_t_row.n_cumulative_height_sum;
					_ASSERTE(n_y + r_t_row.n_height <= n_x); // make sure that the diagonal part of the vector is not affected here

					CTriangularSolve_InnerLoop<_TyGppContext, _CSelectedRowHeightsList,
						CColumnWidth>::BackLoop(diag_x, p_x + n_y, r_t_block.second, n_height);
					// call the inner loop
				}
				// for every other block in this column (note no reverse order here)
				// note that this loop mostly executes once (optimize for that? how?)

				return true;
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

			/**
			 *	@brief prints (unindented) pseudocode of the loop to stdout
			 *	@note This is only available if
			 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
			 */
			static void Debug()
			{
				printf("_ASSERTE(r_t_column.n_width == %d);\n", CColumnWidth::n_size);
				printf("_TyDiagBlock diag_block(r_t_column.block_list.back().second);\n");
				printf("// get diag_x and backsubstitute if doing that\n");
				printf("for(each p_block_it in r_t_column.block_list) {\n");
				printf("if(is_diagonal(*p_block_it))\ncontinue;\n"); // not diagonal block, that is handle separately
				CTriangularSolve_InnerLoop<_TyGppContext, _CSelectedRowHeightsList, CColumnWidth>::Debug();
				printf("}\n");
				printf("// forward substitute diag_x if doing that\n");
				printf("}\n");
			}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		// t_odo - distribute and document the below ifdef, add ifdef around the whole FBS stuff
		// t_odo - take the function above and put it in the matrix body
		// t_odo - remove typedefs to CUberBlockMatrix_Base types at the beginning (especially the allocator type)

		/**
		 *	@brief performs forward-substitution on a single column
		 *
		 *	@param[in] r_t_column is the current column
		 *	@param[in,out] p_x is the left-hand side vector (becomes right-hand side as the system is solved)
		 *	@param[in] r_row_list is list of matrix rows
		 *
		 *	@return Returns true on success, false in case the system doesn't have exactly one solution.
		 */
		static __forceinline bool TriangularSolve_Forward_OuterLoop(const TColumn &r_t_column,
				double *p_x, const std::vector<TRow> &r_row_list)
		{
			return CTriangularSolve_OuterLoop<CBlockMatrixTypelist,
				CColumnWidthsList>::ForwardLoop(r_t_column, p_x, r_row_list);
		}

		/**
		 *	@brief performs back-substitution on a single column
		 *
		 *	@param[in] r_t_column is the current column
		 *	@param[in,out] p_x is the left-hand side vector (becomes right-hand side as the system is solved)
		 *	@param[in] r_row_list is list of matrix rows
		 *
		 *	@return Returns true on success, false in case the system doesn't have exactly one solution.
		 */
		static __forceinline bool TriangularSolve_Back_OuterLoop(const TColumn &r_t_column,
				double *p_x, const std::vector<TRow> &r_row_list)
		{
			return CTriangularSolve_OuterLoop<CBlockMatrixTypelist,
				CColumnWidthsList>::BackLoop(r_t_column, p_x, r_row_list);
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unindented) pseudocode of the triangular solver function to stdout
		 *	@note This function is only available if the
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING macro is defined.
		 */
		static void Debug_TriangularSolve()
		{
			printf("for(each p_col_it in m_block_cols_list) {\n");
			CTriangularSolve_OuterLoop<CBlockMatrixTypelist, CColumnWidthsList>::Debug();
			printf("}\n\n");
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief fixed block size operations template class
	 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
	 *		matrices with known compile-time sizes
	 */
	template <class CBlockMatrixTypelist>
	class CFBS_BlockwiseUnaryOp : public fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist> {
	public:
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CRowHeightsList CRowHeightsList; /**< @brief list of unique block row heights */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CColumnWidthsList CColumnWidthsList; /**< @brief list of unique block column widths */
		typedef typename fbs_ut::deprecate::CFixedBlockSize_UnaryBase<CBlockMatrixTypelist>::CDimsList_Uniq CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */

		// todo - make loops with >1 blocks (like if i want to process all blocks of the matrix in a loop)
		// todo - make loops with a single operand (src = dest), now it is unary operation, but where src != dest

		/**
		 *	@brief (inner) loop of the blockwise unary kernel template
		 *
		 *	Processes a single block by unary operation.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CColumnWidthList is list of posible widths of the blocks
		 *	@tparam CUnaryOperator is unary functor, or a pointer to function type
		 */
		template <class _TyGppContext, class _CColumnWidthList, class CUnaryOperator>
		class CBlockwiseUnary_Loop {
		public:
			typedef typename _CColumnWidthList::_TyHead CCurrentColumnWidth; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[out] p_block_dest is current matrix block
			 *	@param[in] p_block_src is current matrix block
			 *	@param[in] n_block_row_num is number of elements in the current matrix block
			 *	@param[in] n_block_column_num is number of elements in the current matrix block
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Static_Op(double *p_block_dest, const double *p_block_src,
				size_t n_block_row_num, size_t n_block_column_num)
			{
				if(n_block_column_num == CCurrentColumnWidth::n_size) {
					CBlockwiseUnary_Loop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd>, CUnaryOperator>::Static_Op(p_block_dest, p_block_src,
						n_block_row_num, n_block_column_num);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CBlockwiseUnary_Loop<_TyGppContext, typename _CColumnWidthList::_TyTail,
						CUnaryOperator>::Static_Op(p_block_dest, p_block_src, n_block_row_num,
						n_block_column_num);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

			/**
			 *	@brief wraps loop body in block size decision tree
			 *
			 *	@param[out] p_block_dest is current matrix block
			 *	@param[in] p_block_src is current matrix block
			 *	@param[in] n_block_row_column_num is number of rows and columns in the current matrix block
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Static_Op_Square(double *p_block_dest, const double *p_block_src,
				size_t n_block_row_column_num)
			{
				if(n_block_row_column_num == CCurrentColumnWidth::n_size) {
					CBlockwiseUnary_Loop<_TyGppContext, CTypelist<CCurrentColumnWidth,
						CTypelistEnd>, CUnaryOperator>::Static_Op_Square(p_block_dest,
						p_block_src, n_block_row_column_num);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CBlockwiseUnary_Loop<_TyGppContext, typename _CColumnWidthList::_TyTail,
						CUnaryOperator>::Static_Op_Square(p_block_dest, p_block_src,
						n_block_row_column_num);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief (inner) loop of the blockwise unary kernel template
		 *
		 *	Processes a single block by unary operation.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CColumnWidth is the selected width of the block
		 *	@tparam CUnaryOperator is unary functor, or a pointer to function type
		 */
		template <class _TyGppContext, class CColumnWidth, class CUnaryOperator>
		class CBlockwiseUnary_Loop<_TyGppContext, CTypelist<CColumnWidth, CTypelistEnd>, CUnaryOperator> {
		public:
			typedef typename CMakeMatrixRef<CColumnWidth::n_size, CColumnWidth::n_size>::_Ty _TySquareBlock; /**< @brief  matrix diagonal block Eigen::Map type */
			typedef typename CFilterTypelist2<CRowHeightsList, CColumnWidth,
				fbs_ut::CHaveRowHeightForColumnWidth, CDimsList_Uniq>::_TyResult _CSelectedRowHeightsList; /**< @brief row heights list, filtered by the selected column width */

			/**
			 *	@brief wraps operation body in block size decision tree
			 *
			 *	@param[out] p_block_dest is current matrix block
			 *	@param[in] p_block_src is current matrix block
			 *	@param[in] n_block_row_num is number of elements in the current matrix block
			 *	@param[in] n_block_column_num is number of elements in the current matrix block
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Static_Op(double *p_block_dest, const double *p_block_src,
				size_t n_block_row_num, size_t UNUSED(n_block_column_num))
			{
				_ASSERTE(n_block_column_num == CColumnWidth::n_size);
				CBlockwiseUnary_Loop2<_TyGppContext, _CSelectedRowHeightsList,
					CColumnWidth, CUnaryOperator>::Static_Op(p_block_dest, p_block_src,
					n_block_row_num, n_block_column_num);
				// execute inner loop (use the specialization for end-of-the-list that actually
				// contains the implementaiton - avoids having the same code in two places)
			}

			/**
			 *	@brief performs unary operation on a single block
			 *
			 *	@param[out] p_block_dest is current matrix block
			 *	@param[in] p_block_src is current matrix block
			 *	@param[in] n_block_row_column_num is number of rows and columns in the current matrix block
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Static_Op_Square(double *p_block_dest, const double *p_block_src,
				size_t UNUSED(n_block_row_column_num))
			{
				_ASSERTE(n_block_row_column_num == CColumnWidth::n_size);

				_TySquareBlock dest(p_block_dest);
				_TySquareBlock src((double*)p_block_src);
				CUnaryOperator::template Do<_TySquareBlock>(dest, src);
				// do the operation on the two blocks (dest, src)
			}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
		};

		/**
		 *	@brief (inner-most) loop of the blockwise unary kernel template
		 *
		 *	Processes a single block by unary operation.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam _CRowHeightList is list of posible heights of the blocks
		 *	@tparam CColumnWidth is the selected width of the block
		 *	@tparam CUnaryOperator is unary functor, or a pointer to function type
		 */
		template <class _TyGppContext, class _CRowHeightList, class CColumnWidth, class CUnaryOperator>
		class CBlockwiseUnary_Loop2 {
		public:
			typedef typename _CRowHeightList::_TyHead CCurrentRowHeight; /**< @brief column width for this decision tree recursion */

			/**
			 *	@brief wraps operation body in block size decision tree
			 *
			 *	@param[out] p_block_dest is current matrix block
			 *	@param[in] p_block_src is current matrix block
			 *	@param[in] n_block_row_num is number of elements in the current matrix block
			 *	@param[in] n_block_column_num is number of elements in the current matrix block
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Static_Op(double *p_block_dest, const double *p_block_src,
				size_t n_block_row_num, size_t UNUSED(n_block_column_num))
			{
				if(n_block_row_num == CCurrentRowHeight::n_size) {
					CBlockwiseUnary_Loop2<_TyGppContext, CTypelist<CCurrentRowHeight,
						CTypelistEnd>, CColumnWidth, CUnaryOperator>::Static_Op(p_block_dest,
						p_block_src, n_block_row_num, n_block_column_num);
					// execute inner loop (use the specialization for end-of-the-list that actually
					// contains the implementaiton - avoids having the same code in two places)
				} else {
					CBlockwiseUnary_Loop2<_TyGppContext, typename _CRowHeightList::_TyTail,
						CColumnWidth, CUnaryOperator>::Static_Op(p_block_dest, p_block_src,
						n_block_row_num, n_block_column_num);
					// not CCurrentRowHeight::n_size, gotta be one of the rest of the list
				}
			}
		};

		/**
		 *	@brief (inner-most) loop of the blockwise unary kernel template
		 *
		 *	Processes a single block by unary operation.
		 *
		 *	@tparam _TyGppContext is template instance context for g++ (compatibility workarround)
		 *	@tparam CRowHeight is the selected height of the block
		 *	@tparam CColumnWidth is the selected width of the block
		 *	@tparam CUnaryOperator is unary functor, or a pointer to function type
		 */
		template <class _TyGppContext, class CRowHeight, class CColumnWidth, class CUnaryOperator>
		class CBlockwiseUnary_Loop2<_TyGppContext,
			CTypelist<CRowHeight, CTypelistEnd>, CColumnWidth, CUnaryOperator> {
		public:
			typedef typename CMakeMatrixRef<CRowHeight::n_size, CColumnWidth::n_size>::_Ty _TyBlock; /**< @brief  matrix block Eigen::Map type */

			/**
			 *	@brief performs unary operation on a single block
			 *
			 *	@param[out] p_block_dest is current matrix block
			 *	@param[in] p_block_src is current matrix block
			 *	@param[in] n_block_row_num is number of elements in the current matrix block
			 *	@param[in] n_block_column_num is number of elements in the current matrix block
			 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
			#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
			static __forceinline void Static_Op(double *p_block_dest, const double *p_block_src,
				size_t UNUSED(n_block_row_num), size_t UNUSED(n_block_column_num))
			{
				_ASSERTE(n_block_row_num == CRowHeight::n_size);

				_TyBlock dest(p_block_dest);
				_TyBlock src((double*)p_block_src);
				CUnaryOperator::template Do<_TyBlock>(dest, src);
				// do the operation on the two blocks (dest, src)
			}
		};

		/**
		 *	@brief loops over block elements and applies unary operator to their values
		 *
		 *	@tparam CUnnaryOperator is unary operator object type, must implement
		 *		static function Do, parametrizable by Eigen::Map specialization
		 *
		 *	@param[out] p_block_dest is pointer ot the current matrix block data
		 *	@param[in] p_block_src is pointer ot the current matrix block data
		 *	@param[in] n_block_row_num is number of elements in the current matrix block
		 *	@param[in] n_block_column_num is number of elements in the current matrix block
		 */
		template <class CUnaryOperator>
		static __forceinline void BlockwiseUnary_Static_Op(double *p_block_dest,
			const double *p_block_src, size_t n_block_row_num,
			size_t n_block_column_num)
		{
			CBlockwiseUnary_Loop<CBlockMatrixTypelist, CColumnWidthsList,
				CUnaryOperator>::Static_Op(p_block_dest, p_block_src,
				n_block_row_num, n_block_column_num);
		}

		/**
		 *	@brief loops over (square) block elements and applies unary operator to their values
		 *
		 *	@tparam CUnnaryOperator is unary operator object type, must implement
		 *		static function Do, parametrizable by Eigen::Map specialization
		 *
		 *	@param[out] p_block_dest is pointer ot the current matrix block data
		 *	@param[in] p_block_src is pointer ot the current matrix block data
		 *	@param[in] n_block_row_column_num is number of elements in the current matrix block
		 */
		template <class CUnaryOperator>
		static __forceinline void BlockwiseUnary_Static_Op_Square(double *p_block_dest,
			const double *p_block_src, size_t n_block_row_column_num)
		{
			CBlockwiseUnary_Loop<CBlockMatrixTypelist, CColumnWidthsList,
				CUnaryOperator>::Static_Op_Square(p_block_dest, p_block_src,
				n_block_row_column_num);
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		// todo - write that as well

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief a simple std::pair lookalike, only the types are references
	 *
	 *	@tparam _TyFirst is the first type (can be a reference type)
	 *	@tparam _TySecond is the second type (can be a reference type)
	 */
	template <class _TyFirst, class _TySecond>
	class ref_pair {
	public:
		_TyFirst first; /**< @brief the first item in the pair */
		_TySecond second; /**< @brief the second item in the pair */

		/**
		 *	@brief default constructor
		 *
		 *	@param[in] _first is the first item in the pair
		 *	@param[in] _second is the second item in the pair
		 */
		inline ref_pair(_TyFirst _first, _TySecond _second)
			:first(_first), second(_second)
		{}

		/**
		 *	@brief copy-constructor
		 *	@param[in] r_other is the other pair to initialize from
		 */
		inline ref_pair(ref_pair &r_other)
			:first(r_other.first), second(r_other.second)
		{}
	};

	// todo - use decision trees all the way through this file, record the compile times

	/**
	 *	@brief simple static functor
	 */
	class CCallDenseCholesky {
	protected:
		CUberBlockMatrix &m_r_matrix;
		bool &m_r_b_result;

	public:
		CCallDenseCholesky(CUberBlockMatrix &r_matrix, bool &r_b_result)
			:m_r_matrix(r_matrix), m_r_b_result(r_b_result)
		{}

		/**
		 *	@brief function operator-like function
		 */
		template <class CMatrixSize>
		inline void operator ()()
		{
			m_r_b_result = m_r_matrix.Cholesky_Dense<CMatrixSize::n_size>();
			// call debse cholesky and store result in context
		}
	};
};

} // ~blockmatrix_detail

#endif // !__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_OPS_IMPLEMENTATION_INCLUDED
