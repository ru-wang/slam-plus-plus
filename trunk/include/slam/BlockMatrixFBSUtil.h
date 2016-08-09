/*
								+---------------------------------+
								|                                 |
								|  ***   Über Block Matrix   ***  |
								|                                 |
								| Copyright  (c) -tHE SWINe- 2012 |
								|                                 |
								|      BlockMatrixFBSUtil.h       |
								|                                 |
								+---------------------------------+
*/

#pragma once
#ifndef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_OPS_UTILITIES_INCLUDED
#define __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_OPS_UTILITIES_INCLUDED

/**
 *	@file include/slam/BlockMatrixFBSUtil.h
 *	@date 2012
 *	@author -tHE SWINe-
 *	@brief the überblockmatrix fixed block size utility classes
 *	@note This file is not to be included; it is automatically included from BlockMatrix.h
 */

#include "slam/BlockMatrixBase.h"

/**
 *	@brief fixed block size utility classes and functions
 */
namespace fbs_ut {

/**
 *	@brief compile-time constant booelan
 *	@tparam _n_size is binary flag
 */
template <bool _b_flag>
class CCTFlag {
public:
	/**
	 *	@brief copy of template parameters as enums
	 */
	enum {
		b_flag = _b_flag /**< @brief the flag value */
	};
};

/**
 *	@brief compile-time constant scalar
 *	@tparam _n_size is scalar size
 */
template <int _n_size>
class CCTSize {
public:
	/**
	 *	@brief copy of template parameters as enums
	 */
	enum {
		n_size = _n_size /**< @brief size */
	};
};

/**
 *	@brief compile-time constant 2D vector
 *
 *	@tparam _n_row_num is number of rows (y component)
 *	@tparam _n_column_num is number of columns (x component)
 */
template <int _n_row_num, int _n_column_num>
class CCTSize2D {
public:
	/**
	 *	@brief copy of template parameters as enums
	 */
	enum {
		n_row_num = _n_row_num, /**< @brief number of rows */
		n_column_num = _n_column_num /**< @brief number of columns */
	};
};

/**
 *	@brief binary addition operating on scalars
 *
 *	@tparam _T1 is the first operand (specialization of CCTSize template)
 *	@tparam _T2 is the second operand (specialization of CCTSize template)
 */
template <class _T1, class _T2>
class CBinaryScalarAdd {
public:
	typedef CCTSize<_T1::n_size + _T2::n_size> _TyResult; /**< @brief result of the addition */
};

/**
 *	@brief binary maximum of two compile-time scalars
 *
 *	@tparam _T1 is the first operand (specialization of CCTSize template)
 *	@tparam _T2 is the second operand (specialization of CCTSize template)
 */
template <class _T1, class _T2>
class CBinaryScalarMax {
protected:
	/**
	 *	@brief intermediates stored as an enum
	 */
	enum {
		s0 = _T1::n_size, /**< @brief value of the first scalar */
		s1 = _T2::n_size, /**< @brief value of the second scalar */
		n_result = (s0 > s1)? s0 : s1 /**< @brief the maximum (in here to avoid template problems when using the greater-than operator) */
	};

public:
	typedef CCTSize<n_result> _TyResult; /**< @brief result of the addition */
};

/**
 *	@brief binary less than or equal comparison operating on scalars
 *
 *	@tparam _T1 is the first operand (specialization of CCTSize template)
 *	@tparam _T2 is the second operand (specialization of CCTSize template)
 */
template <class _T1, class _T2>
class CCompareScalar_LEqual {
public:
	/**
	 *	@brief result stored as an enum
	 */
	enum {
		b_result = size_t(_T1::n_size) <= size_t(_T2::n_size) /**< @brief result of the comparison */
	};
	// converting to size_t disables g++ warning: comparison between 'enum CCTSize<X>::<anonymous>' and 'enum CCTSize<Y>::<anonymous>'
};

/**
 *	@brief binary less than comparison operating on scalars
 *
 *	@tparam _T1 is the first operand (specialization of CCTSize template)
 *	@tparam _T2 is the second operand (specialization of CCTSize template)
 */
template <class _T1, class _T2>
class CCompareScalar_Less {
public:
	/**
	 *	@brief result stored as an enum
	 */
	enum {
		b_result = size_t(_T1::n_size) < size_t(_T2::n_size) /**< @brief result of the comparison */
	};
	// converting to size_t disables g++ warning: comparison between 'enum CCTSize<X>::<anonymous>' and 'enum CCTSize<Y>::<anonymous>'
};

/**
 *	@brief binary less than comparison operating on size 2D (number of columns
 *		being more significant, number of rows les significant)
 *
 *	@tparam _T1 is the first operand (specialization of CCTSize2D template)
 *	@tparam _T2 is the second operand (specialization of CCTSize2D template)
 */
template <class _T1, class _T2>
class CCompareSize2D_Less {
public:
	/**
	 *	@brief result stored as an enum
	 */
	enum {
		b_result = size_t(_T1::n_column_num) < size_t(_T2::n_column_num) ||
			(size_t(_T1::n_column_num) == size_t(_T2::n_column_num) &&
			size_t(_T1::n_row_num) < size_t(_T2::n_row_num)) /**< @brief result of the comparison */
	};
	// converting to size_t disables g++ warning: comparison between 'enum CCTSize2D<X>::<anonymous>' and 'enum CCTSize2D<Y>::<anonymous>'
};

/**
 *	@brief gets dimension for a given vertex
 *	@tparam CVertex is vertex type
 */
template <class CVertex>
class CGetVertexDimension {
public:
	typedef fbs_ut::CCTSize<CVertex::n_dimension> _TyResult; /**< @brief vertex dimension */
};

/**
 *	@brief extracts number of columns from compile-time constant 2D vector
 *	@tparam CDimensionType is a compile-time constant 2D vector
 */
template <class CDimensionType>
class CTransformDimensionColumnsToSize {
public:
	typedef CCTSize<CDimensionType::n_column_num> _TyResult; /**< @brief number of columns, extracted from CDimensionType */
};

/*template <int _n_row_num, int _n_column_num>
class CTransformDimensionColumnsToSize<CCTSize2D<_n_row_num, _n_column_num> > {
public:
	typedef CCTSize<_n_column_num> _TyResult;
};*/

/**
 *	@brief extracts number of rows from compile-time constant 2D vector
 *	@tparam CDimensionType is a compile-time constant 2D vector
 */
template <class CDimensionType>
class CTransformDimensionRowsToSize {
public:
	typedef CCTSize<CDimensionType::n_row_num> _TyResult; /**< @brief number of rows, extracted from CDimensionType */
};

/*template <int _n_row_num, int _n_column_num>
class CTransformDimensionRowsToSize<CCTSize2D<_n_row_num, _n_column_num> > {
public:
	typedef CCTSize<_n_row_num> _TyResult;
};*/

/**
 *	@brief calculates area of a rectangle from compile-time constant 2D vector
 *	@tparam CDimensionType is a compile-time constant 2D vector
 */
template <class CDimensionType>
class CTransformDimensionToAreaSize {
public:
	typedef CCTSize<CDimensionType::n_column_num * CDimensionType::n_row_num> _TyResult; /**< @brief number of columns times number of rows, extracted from CDimensionType */
};

/**
 *	@brief converts eigen matrix type to simple CCTSize2D
 *	@tparam CMatrixType is specialization of Eigen::Matrix
 */
template <class CMatrixType>
class CEigenToDimension {
public:
#if !defined(_MSC_VER) || defined(__MWERKS__) // the below version does not work with Intellisense, the one in else branch does
	typedef CCTSize2D<CMatrixType::RowsAtCompileTime,
		CMatrixType::ColsAtCompileTime> _TyResult; /**< @brief size of the matrix, represented as CCTSize2D */
#else // !_MSC_VER || __MWERKS__
	typedef Eigen::internal::traits<CMatrixType> _TyMatrixTraits; /**< @brief matrix traits */
	typedef CCTSize2D<_TyMatrixTraits::RowsAtCompileTime,
		_TyMatrixTraits::ColsAtCompileTime> _TyResult; /**< @brief size of the matrix, represented as CCTSize2D */
#endif // !_MSC_VER || __MWERKS__
};

/**
 *	@brief converts eigen matrix type to simple CCTSize2D
 *		(specialization for input already being CCTSize2D)
 *
 *	@tparam _n_row_num is number of rows (y component)
 *	@tparam _n_column_num is number of columns (x component)
 */
template <int _n_row_num, int _n_column_num>
class CEigenToDimension<CCTSize2D<_n_row_num, _n_column_num> > {
public:
	typedef CCTSize2D<_n_row_num, _n_column_num> _TyResult; /**< @brief size of the matrix, represented as CCTSize2D */
};

#if 0 // discarded as unnecessary complication, Eigen matrix types have RowsAtCompileTime and ColsAtCompileTime enum

/**
 *	@brief specialization for matrices of double (the only supported type)
 *
 *	@tparam n_compile_time_row_num number of rows
 *	@tparam n_compile_time_col_num number of columns
 *	@tparam n_options a combination of either Eigen::RowMajor or Eigen::ColMajor,
 *		and of either Eigen::AutoAlign or Eigen::DontAlign
 *	@tparam n_max_row_num maximum number of rows
 *	@tparam n_max_col_num maximum number of columns
 */
template <int n_compile_time_row_num, int n_compile_time_col_num,
	int n_options, int n_max_row_num, int n_max_col_num>
class CEigenToDimension<Eigen::Matrix<double, n_compile_time_row_num,
	n_compile_time_col_num, n_options, n_max_row_num, n_max_col_num> > {
public:
	typedef CCTSize2D<n_compile_time_row_num, n_compile_time_col_num> _TyResult; /**< @brief size of the matrix, represented as CCTSize2D */
};

#endif // 0

/**
 *	@brief transposes shape of Eigen::Matrix
 *	@tparam CMatrixType is specialization of Eigen::Matrix
 */
template <class CMatrixType>
class CEigenToTransposeEigen;

/**
 *	@brief specialization for matrices of double (the only supported type)
 *
 *	@tparam n_compile_time_row_num number of rows
 *	@tparam n_compile_time_col_num number of columns
 *	@tparam n_options a combination of either Eigen::RowMajor or Eigen::ColMajor,
 *		and of either Eigen::AutoAlign or Eigen::DontAlign
 *	@tparam n_max_row_num maximum number of rows
 *	@tparam n_max_col_num maximum number of columns
 */
template <int n_compile_time_row_num, int n_compile_time_col_num,
	int n_options, int n_max_row_num, int n_max_col_num>
class CEigenToTransposeEigen<Eigen::Matrix<double, n_compile_time_row_num,
	n_compile_time_col_num, n_options, n_max_row_num, n_max_col_num> > {
public:
	typedef Eigen::Matrix<double, n_compile_time_col_num,
		n_compile_time_row_num, n_options, n_max_col_num, n_max_row_num> _TyResult; /**< @brief transposed shape of the original matrix */
};

/**
 *	@brief compile-time assertion (_FBS functions do not work
 *		with dynamic-size block matrices)
 *	@tparam b_expression is the expression being asserted
 */
template <const bool b_expression>
class CMATRIX_BLOCK_DIMENSIONS_MUST_BE_KNOWN_AT_COMPILE_TIME;

/**
 *	@brief compile-time assertion (specialization for assertion passing)
 *	@tparam b_expression is the expression being asserted
 */
template <>
class CMATRIX_BLOCK_DIMENSIONS_MUST_BE_KNOWN_AT_COMPILE_TIME<true> {};

/**
 *	@brief predicate for filtering row height list by selected column width
 *
 *	@tparam CRowHeight is row height (value from the list that is being filtered)
 *	@tparam CColumnWidth is column width (reference value)
 */
template <class CDimsList_Uniq, class CRowHeight, class CColumnWidth>
class CHaveRowHeightForColumnWidth {
public:
	typedef CCTSize2D<CRowHeight::n_size, CColumnWidth::n_size> _TyNeedle; /**< @brief hypothetical matrix size */

	/**
	 *	@brief result stored as enum
	 */
	enum {
		b_result = CFindTypelistItem<CDimsList_Uniq, _TyNeedle>::b_result /**< @brief predicate result (looks for the matrix size in the list) */
	};
};

/**
 *	@brief predicate for filtering column width list by selected row height
 *
 *	@tparam CColumnWidth is column width (value from the list that is being filtered)
 *	@tparam CRowHeight is row height (reference value)
 */
template <class CDimsList_Uniq, class CColumnWidth, class CRowHeight>
class CHaveColumnWidthForRowHeight {
public:
	typedef CCTSize2D<CRowHeight::n_size, CColumnWidth::n_size> _TyNeedle; /**< @brief hypothetical matrix size */

	/**
	 *	@brief result stored as enum
	 */
	enum {
		b_result = CFindTypelistItem<CDimsList_Uniq, _TyNeedle>::b_result /**< @brief predicate result (looks for the matrix size in the list) */
	};
};

/**
 *	@brief conversion of a pair of CCTSize back to an Eigen matrix type with
 *		known compile-time sizes (the result can be found in the _TyResult type)
 *
 *	@tparam _T1 is a specialization of CCTSize, containing number of matrix rows
 *	@tparam _T2 is a specialization of CCTSize, containing number of matrix columns
 */
template <class _T1, class _T2>
class CMakeMatrixSizeType {
public:
	typedef Eigen::Matrix<double, _T1::n_size, _T2::n_size> _TyResult; /**< @brief the resulting Eigen matrix type */
};

/**
 *	@brief converts compile-time constant 2D vector to Eigen matrix type
 *	@tparam CDimensionType is a compile-time constant 2D vector
 */
template <class CDimensionType>
class CDimensionToEigen {
public:
	typedef Eigen::Matrix<double, CDimensionType::n_row_num, CDimensionType::n_column_num> _TyResult; /**< @brief Eigen matrix type, corresponding to the input block size */
};

/**
 *	@brief conversion of a pair of CCTSize to CCTSize2D matrix type with
 *		known compile-time sizes (the result can be found in the _TyResult type)
 *
 *	@tparam _T1 is a specialization of CCTSize, containing number of matrix rows
 *	@tparam _T2 is a specialization of CCTSize, containing number of matrix columns
 */
template <class _T1, class _T2>
class CMakeCTSize2DType {
public:
	typedef CCTSize2D<_T1::n_size, _T2::n_size> _TyResult; /**< @brief the resulting 2D size type */
};

/**
 *	@brief calculates size of product size in form L^T * R
 *
 *	@tparam CSizeLeft is size of the left matrix
 *	@tparam CSizeRight is size of the right matrix
 */
template <class CSizeLeft, class CSizeRight>
class CMakeATAProductSize {
protected:
	/**
	 *	@brief intermediates, stored as enum
	 */
	enum {
		n_left_col_num = CSizeLeft::n_row_num, /**< @brief number of columns of the tranposed left block */
		n_left_row_num = CSizeLeft::n_column_num, /**< @brief number of rows of the tranposed left block */
		// note that this is transposed (the product is A^T * A)

		n_right_row_num = CSizeRight::n_row_num, /**< @brief number of rows of the right block */
		n_right_col_num = CSizeRight::n_column_num, /**< @brief number of columns of the right block */

		b_can_multiply = (int(n_left_col_num) == int(n_right_row_num)) /**< @brief common dimension check */
		// int cast to avoid "warning: comparison between 'enum fbs_ut::CCTSize2D<X, Y>::<anonymous>' and 'enum fbs_ut::CCTSize2D<Z, W>::<anonymous>' [-Wenum-compare]"
	};

public:
	typedef CCTSize2D<(b_can_multiply)? n_left_row_num : -1,
		(b_can_multiply)? n_right_col_num : -1> _TyResult; /**< @brief the resulting block size, or CCTSize<-1, -1> if not multiplicable */
};

/**
 *	@brief nonnegative block size predicate
 *	@tparam CSize is size of the block
 */
template <class CSize>
class CIsNotNegsize {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		b_result = (CSize::n_row_num > 0 && CSize::n_column_num > 0) /**< @brief result; true if the size is a valid block size */
	};
};

/**
 *	@brief square block predicate
 *	@tparam CSize is size of the block
 */
template <class CSize>
class CIsSquare {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		b_result = (int(CSize::n_row_num) == int(CSize::n_column_num)) /**< @brief result; true if the size is square */
	};
};

/**
 *	@brief runtime 2D size comparator
 *	@tparam C2DSize is size of the block
 */
template <class C2DSize>
class CRuntimeCompareSize2D {
public:
	/**
	 *	@brief equality predicate
	 *	@param[in] t_size is pair with the number of rows and columns (in this order)
	 *	@return Returns true if the size specified by C2DSize is equal to t_size.
	 */
	static inline bool b_Equal(std::pair<size_t, size_t> t_size)
	{
		return size_t(C2DSize::n_row_num) == t_size.first && size_t(C2DSize::n_column_num) == t_size.second;
	}

	/**
	 *	@brief less-than predicate
	 *	@param[in] t_size is pair with the number of rows and columns (in this order)
	 *	@return Returns true if the size specified by C2DSize is less than t_size.
	 *	@note The comparison is the same as in CCompareSize2D_Less.
	 */
	static inline bool b_LessThan(std::pair<size_t, size_t> t_size)
	{
		return size_t(C2DSize::n_row_num) < t_size.first || (size_t(C2DSize::n_row_num) == t_size.first &&
			size_t(C2DSize::n_column_num) < t_size.second);
	}

	/**
	 *	@brief greater-than predicate
	 *	@param[in] t_size is pair with the number of rows and columns (in this order)
	 *	@return Returns true if the size specified by C2DSize is greater than t_size.
	 *	@note The comparison has the same ordering on first / second as in CCompareSize2D_Less.
	 */
	static inline bool b_GreaterThan(std::pair<size_t, size_t> t_size)
	{
		return size_t(C2DSize::n_row_num) > t_size.first || (size_t(C2DSize::n_row_num) == t_size.first &&
			size_t(C2DSize::n_column_num) > t_size.second);
	}
};

/**
 *	@brief runtime size comparator
 *	@tparam CScalar is instance of CCTSize
 */
template <class CScalar>
class CRuntimeCompareScalar {
public:
	/**
	 *	@brief equality predicate
	 *	@param[in] n_size is value of the scalar
	 *	@return Returns true if the value specified by CScalar is equal to n_size.
	 */
	static inline bool b_Equal(size_t n_size)
	{
		return size_t(CScalar::n_size) == n_size;
	}

	/**
	 *	@brief less-than predicate
	 *	@param[in] n_size is value of the scalar
	 *	@return Returns true if the value specified by CScalar is less than n_size.
	 *	@note The comparison is the same as in CCompareSize2D_Less.
	 */
	static inline bool b_LessThan(size_t n_size)
	{
		return size_t(CScalar::n_size) < n_size;
	}

	/**
	 *	@brief greater-than predicate
	 *	@param[in] n_size is value of the scalar
	 *	@return Returns true if the value specified by CScalar is greater than n_size.
	 *	@note The comparison has the same ordering on first / second as in CCompareSize2D_Less.
	 */
	static inline bool b_GreaterThan(size_t n_size)
	{
		return size_t(CScalar::n_size) > n_size;
	}
};

/**
 *	@brief calculates a typelist of block sizes after the PreMultiplyWithSelfTranspose()
 *		operation (the result can be found in the _TyResult type)
 *	@tparam CBlockSizeTypelist is typelist, containing
 *		list of block sizes as CCTSize2D
 */
template <class CBlockSizeTypelist>
class CBlockSizesAfterPreMultiplyWithSelfTranspose {
protected:
	typedef CBlockSizeTypelist CDimsList; /**< @brief list of block sizes as CCTSize2D */
	typedef typename CUniqueTypelist<CDimsList>::_TyResult CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */
	//typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_Uniq,
	//	CTransformDimensionRowsToSize>::_TyResult>::_TyResult CRowHeightsList; /**< @brief list of unique block row heights */
	typedef typename CUniqueTypelist<typename CCarthesianProductTypelist<CDimsList_Uniq,
		CDimsList_Uniq, CMakeATAProductSize>::_TyResult>::_TyResult CProdList; /**< @brief list of resulting block sizes and CCTSize<-1, -1> */

public:
	typedef typename CFilterTypelist1<CProdList, CIsNotNegsize>::_TyResult _TyResult; /**< @brief the resulting typelist */
	//typedef typename CUniqueTypelist<typename CCarthesianProductTypelist<CRowHeightsList,
	//	CRowHeightsList, CMakeCTSize2DType>::_TyResult>::_TyResult _TyResult; /**< @brief the resulting typelist */
};

/**
 *	@brief calculates a typelist of block sizes after the PreMultiplyWithSelfTranspose()
 *		operation (the result can be found in the _TyResult type)
 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
 *		matrices with known compile-time sizes
 */
template <class CBlockMatrixTypelist>
class CBlockMatrixTypesAfterPreMultiplyWithSelfTranspose {
public:
	typedef typename CBlockSizesAfterPreMultiplyWithSelfTranspose<CBlockMatrixTypelist>::_TyResult _TySizeList; /**< @brief the resulting size list (as CCTSize2D) */
	typedef typename CTransformTypelist<_TySizeList, CDimensionToEigen>::_TyResult _TyResult; /**< @brief the resulting typelist */
};

/**
 *	@brief prints typelist containing Eigen matrices to stdout (debugging utility)
 *
 *	Use as follows:
 *	@code
 *	typedef MakeTypelist_Safe((Eigen::Matrix2d, Eigen::Matrix3d)) matrices;
 *	fbs_ut::CDumpBlockMatrixTypelist<matrices>::Print();
 *	@endcode
 *
 *	@tparam CBlockMatrixTypelist is typelist containing Eigen matrices
 */
template <class CBlockMatrixTypelist>
class CDumpBlockMatrixTypelist {
public:
	/**
	 *	@brief prints the given typelist to stdout
	 *	@note This assumes that scalar type is double, but in order to make
	 *		the implementation simple, it is only assumed, the scalar type is ignored.
	 */
	static void Print()
	{
		typedef typename CBlockMatrixTypelist::_TyHead THead;
		printf("Eigen::Matrix<double, %d, %d> ", THead::RowsAtCompileTime, THead::ColsAtCompileTime);
		CDumpBlockMatrixTypelist<typename CBlockMatrixTypelist::_TyTail>::Print();
	}
};

/**
 *	@brief prints typelist containing Eigen matrices to stdout
 *		(debugging utility; specialization for the end of the list)
 */
template <>
class CDumpBlockMatrixTypelist<CTypelistEnd> {
public:
	/**
	 *	@brief prints the given typelist to stdout
	 */
	static void Print()
	{
		printf("\n");
	}
};

/**
 *	@brief prints typelist containing CCTSize sizes to stdout (debugging utility)
 *
 *	Use as follows:
 *	@code
 *	typedef MakeTypelist_Safe((fbs_ut::CCTSize<1>, fbs_ut::CCTSize<2>)) sizes;
 *	fbs_ut::CDumpCTSizeTypelist<sizes>::Print();
 *	@endcode
 *
 *	@tparam CCTSizeTypelist is typelist containing CCTSize sizes
 */
template <class CCTSizeTypelist>
class CDumpCTSizeTypelist {
public:
	/**
	 *	@brief prints the given typelist to stdout
	 *	@note This assumes that scalar type is double, but in order to make
	 *		the implementation simple, it is only assumed, the scalar type is ignored.
	 */
	static void Print()
	{
		typedef typename CCTSizeTypelist::_TyHead THead;
		printf("CCTsize<%d> ", THead::n_size);
		CDumpCTSizeTypelist<typename CCTSizeTypelist::_TyTail>::Print();
	}
};

/**
 *	@brief prints typelist containing CCTSize sizes to stdout
 *		(debugging utility; specialization for the end of the list)
 */
template <>
class CDumpCTSizeTypelist<CTypelistEnd> {
public:
	/**
	 *	@brief prints the given typelist to stdout
	 */
	static void Print()
	{
		printf("\n");
	}
};

/**
 *	@brief makes a decission tree from a sorted typelist
 *
 *	@tparam CList is a sorted list of scalars (type CCTSize)
 *	@tparam CContext is a type of context, required by the operation
 *	@tparam CFunctor is a functor template that specializes for a given size,
 *		having a global static function Do(CContext) which is called in the
 *		leafs of the decision tree
 */
template <class CList, class CContext, template <const int> class CFunctor>
class CMakeDecisionTree {
protected:
	/**
	 *	@brief makes a branch of a decission tree from a sorted typelist
	 *
	 *	@tparam _CList is a sorted list of scalars (type CCTSize)
	 *	@tparam _CContext is a type of context, required by the operation
	 *	@tparam _CFunctor is a functor template that specializes for a given size,
	 *		having a global static function Do(CContext) which is called in the
	 *		leafs of the decision tree
	 *	@tparam n_begin is a (zero based) index of the starting element of CList
	 *		this branch decides over
	 *	@tparam n_length is number of elements of CList this branch decides over
	 */
	template <class _CList, class _CContext, template <const int> class _CFunctor,
		const int n_begin, const int n_length>
	class CDecisionTree {
	public:
		/**
		 *	@brief decision indices stored as enums
		 */
		enum {
			n_half = (n_length + 1) / 2, /**< @brief half of the length (round up, comparison is less than) */
			n_pivot = n_begin + n_half, /**< @brief zero-based index of the pivot element */
			n_rest = n_length - n_half /**< @brief the remaining half of the length for the right child */
		};

		typedef typename CTypelistItemAt<_CList, n_pivot>::_TyResult _TyPivot; /**< @brief pivot type */

		/**
		 *	@brief the executive function
		 *
		 *	@param[in] n_size is the parameter to be decided over
		 *	@param[in,out] context is the operation context (task specific)
		 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
		static __forceinline void Do(const size_t n_size, _CContext context)
		{
			if(n_size < _TyPivot::n_size)
				CDecisionTree<_CList, _CContext, _CFunctor, n_begin, n_half>::Do(n_size, context);
			else
				CDecisionTree<_CList, _CContext, _CFunctor, n_pivot, n_rest>::Do(n_size, context);
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unidented) decission tree pseudocode to stdout
		 *	@note This is only available if
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
		 */
		static void Debug()
		{
			//printf("if(n_size < i[%d]) {\n", n_pivot);
			printf("if(n_size < %d) {\n", _TyPivot::n_size);
			CDecisionTree<_CList, _CContext, _CFunctor, n_begin, n_half>::Debug();
			printf("} else {\n");
			CDecisionTree<_CList, _CContext, _CFunctor, n_pivot, n_rest>::Debug();
			printf("}\n");
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief makes a branch of a decission tree from a sorted typelist
	 *		(specialization for leafs)
	 *
	 *	@tparam CList is a sorted list of scalars (type CCTSize)
	 *	@tparam CContext is a type of context, required by the operation
	 *	@tparam CFunctor is a functor template that specializes for a given size,
	 *		having a global static function Do(CContext) which is called in the
	 *		leafs of the decision tree
	 *	@tparam n_begin is a (zero based) index of the starting element of CList
	 *		this branch decides over
	 */
	template <class _CList, class _CContext, template <const int> class _CFunctor, const int n_begin>
	class CDecisionTree<_CList, _CContext, _CFunctor, n_begin, 1> {
	public:
		typedef typename CTypelistItemAt<_CList, n_begin>::_TyResult _TyPivot; /**< @brief pivot type */

		/**
		 *	@brief the executive function
		 *
		 *	@param[in] n_size is the parameter to be decided over
		 *	@param[in,out] context is the operation context (task specific)
		 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
		static __forceinline void Do(const size_t UNUSED(n_size), _CContext context)
		{
			_ASSERTE(n_size == _TyPivot::n_size);
			_CFunctor<_TyPivot::n_size>::Do(context);
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unidented) decission tree pseudocode to stdout
		 *	@note This is only available if
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
		 */
		static void Debug()
		{
			//printf("_ASSERTE(n_size == i[%d]);\n", n_begin);
			printf("_ASSERTE(n_size == %d);\n", _TyPivot::n_size);
			printf("CFunctor<%d>::Do(context);\n", _TyPivot::n_size);
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

public:
	/**
	 *	@brief the executive function
	 *
	 *	@param[in] n_size is the parameter to be decided over
	 *	@param[in,out] context is the operation context (task specific)
	 */
	static __forceinline void Do(const size_t n_size, CContext context)
	{
		CDecisionTree<CList, CContext, CFunctor, 0,
			CTypelistLength<CList>::n_result>::Do(n_size, context);
	}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

	/**
	 *	@brief prints (unidented) decission tree pseudocode to stdout
	 *	@note This is only available if
	 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
	 */
	static void Debug()
	{
		CDecisionTree<CList, CContext, CFunctor, 0, CTypelistLength<CList>::n_result>::Debug();
	}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
};

/**
 *	@brief makes a decission tree from a sorted typelist
 *
 *	@tparam CList is a sorted list of scalars (type CCTSize)
 *	@tparam CContext is a type of context, required by the operation
 *	@tparam CContext2 is a type of context, passed as a second template parameter for the operation
 *	@tparam CFunctor is a functor template that specializes for a given size and
 *		for context which is equal to CContext2, having a global static function
 *		Do(CContext) which is called in the leafs of the decision tree
 */
template <class CList, class CContext, class CContext2,
	template <const int, class> class CFunctor>
class CMakeDecisionTree2 {
protected:
	/**
	 *	@brief makes a branch of a decission tree from a sorted typelist
	 *
	 *	@tparam _CList is a sorted list of scalars (type CCTSize)
	 *	@tparam _CContext is a type of context, required by the operation
	 *	@tparam _CContext2 is a type of context, passed as a second
	 *		template parameter for the operation
	 *	@tparam _CFunctor is a functor template that specializes for a given size,
	 *		having a global static function Do(CContext) which is called in the
	 *		leafs of the decision tree
	 *	@tparam n_begin is a (zero based) index of the starting element of CList
	 *		this branch decides over
	 *	@tparam n_length is number of elements of CList this branch decides over
	 */
	template <class _CList, class _CContext, class _CContext2,
		template <const int, class> class _CFunctor,
		const int n_begin, const int n_length>
	class CDecisionTree {
	public:
		/**
		 *	@brief decision indices stored as enums
		 */
		enum {
			n_half = (n_length + 1) / 2, /**< @brief half of the length (round up, comparison is less than) */
			n_pivot = n_begin + n_half, /**< @brief zero-based index of the pivot element */
			n_rest = n_length - n_half /**< @brief the remaining half of the length for the right child */
		};

		typedef typename CTypelistItemAt<_CList, n_pivot>::_TyResult _TyPivot; /**< @brief pivot type */

		/**
		 *	@brief the executive function
		 *
		 *	@param[in] n_size is the parameter to be decided over
		 *	@param[in,out] context is the operation context (task specific)
		 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
		static __forceinline void Do(const size_t n_size, _CContext context)
		{
			if(n_size < _TyPivot::n_size) {
				CDecisionTree<_CList, _CContext, _CContext2,
					_CFunctor, n_begin, n_half>::Do(n_size, context);
			} else {
				CDecisionTree<_CList, _CContext, _CContext2,
					_CFunctor, n_pivot, n_rest>::Do(n_size, context);
			}
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unidented) decission tree pseudocode to stdout
		 *	@note This is only available if
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
		 */
		static void Debug()
		{
			//printf("if(n_size < i[%d]) {\n", n_pivot);
			printf("if(n_size < %d) {\n", _TyPivot::n_size);
			CDecisionTree<_CList, _CContext, _CContext2, _CFunctor, n_begin, n_half>::Debug();
			printf("} else {\n");
			CDecisionTree<_CList, _CContext, _CContext2, _CFunctor, n_pivot, n_rest>::Debug();
			printf("}\n");
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

	/**
	 *	@brief makes a branch of a decission tree from a sorted typelist
	 *		(specialization for leafs)
	 *
	 *	@tparam CList is a sorted list of scalars (type CCTSize)
	 *	@tparam CContext is a type of context, required by the operation
	 *	@tparam _CContext2 is a type of context, passed as a second
	 *		template parameter for the operation
	 *	@tparam CFunctor is a functor template that specializes for a given size,
	 *		having a global static function Do(CContext) which is called in the
	 *		leafs of the decision tree
	 *	@tparam n_begin is a (zero based) index of the starting element of CList
	 *		this branch decides over
	 */
	template <class _CList, class _CContext, class _CContext2,
		template <const int, class> class _CFunctor, const int n_begin>
	class CDecisionTree<_CList, _CContext, _CContext2, _CFunctor, n_begin, 1> {
	public:
		typedef typename CTypelistItemAt<_CList, n_begin>::_TyResult _TyPivot; /**< @brief pivot type */

		/**
		 *	@brief the executive function
		 *
		 *	@param[in] n_size is the parameter to be decided over
		 *	@param[in,out] context is the operation context (task specific)
		 */
#if defined(_MSC_VER) && !defined(__MWERKS__)
#pragma inline_recursion(on)
#endif // _MSC_VER && !__MWERKS__
		static __forceinline void Do(const size_t UNUSED(n_size), _CContext context)
		{
			_ASSERTE(n_size == _TyPivot::n_size);
			_CFunctor<_TyPivot::n_size, _CContext2>::Do(context);
		}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

		/**
		 *	@brief prints (unidented) decission tree pseudocode to stdout
		 *	@note This is only available if
		 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
		 */
		static void Debug()
		{
			//printf("_ASSERTE(n_size == i[%d]);\n", n_begin);
			printf("_ASSERTE(n_size == %d);\n", _TyPivot::n_size);
			printf("CFunctor<%d, context2>::Do(context);\n", _TyPivot::n_size);
		}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
	};

public:
	/**
	 *	@brief the executive function
	 *
	 *	@param[in] n_size is the parameter to be decided over
	 *	@param[in,out] context is the operation context (task specific)
	 */
	static __forceinline void Do(const size_t n_size, CContext context)
	{
		CDecisionTree<CList, CContext, CContext2, CFunctor, 0,
			CTypelistLength<CList>::n_result>::Do(n_size, context);
	}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

	/**
	 *	@brief prints (unidented) decission tree pseudocode to stdout
	 *	@note This is only available if
	 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
	 */
	static void Debug()
	{
		CDecisionTree<CList, CContext, CContext2, CFunctor, 0,
			CTypelistLength<CList>::n_result>::Debug();
	}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
};

/**
 *	@brief a simple square size predicate
 *	@tparam _TySize is a specialization of CCTSize2D
 */
template <class _TySize>
class CIsSize2DSquare {
public:
	/**
	 *	@brief result stored as enum
	 */
	enum {
		b_result = (_TySize::n_row_num == _TySize::n_column_num)? 1 : 0 /**< @brief nonzero if _TySize is a square, otherwise zero */
	};
};

/**
 *	@brief makes a decission tree for possible block matrix sizes,
 *		given a list of block sizes and a maximal size
 *
 *	@tparam CBlockMatrixTypelist is a list of possible block sizes as Eigen::Matrix
 *	@tparam n_max_matrix_size is maximum matrix size, in elements
 */
template <class CBlockMatrixTypelist, const int n_max_matrix_size>
class CMakeSquareMatrixSizeDecisionTree {
protected:
	typedef CBlockMatrixTypelist _TyBlockMatrixTypelist; /**< @brief list of block sizes, represented as Eigen::Matrix */
	typedef typename CTransformTypelist<_TyBlockMatrixTypelist,
		CEigenToDimension>::_TyResult CDimsList; /**< @brief list of block sizes as CCTSize2D */
	typedef typename CFilterTypelist1<CDimsList, CIsSize2DSquare>::_TyResult CSquareDimsList; /**< @brief list of square block sizes as CCTSize2D */
	typedef typename CTransformTypelist<CSquareDimsList, CTransformDimensionColumnsToSize>::_TyResult CWidthList; /**< @brief list of square block widths as CCTSize */
	typedef typename CBinaryCombinationTypelist<CWidthList, CBinaryScalarAdd,
		CCompareScalar_LEqual, CCTSize<n_max_matrix_size> >::_TyResult _TyPossibleWidthList; /**< @brief list of all possible matrix sizes */
	typedef typename CSortTypelist<_TyPossibleWidthList, CCompareScalar_Less>::_TyResult _TySortedList; /**< @brief sorted list of all possible matrix sizes */

	/**
	 *	@brief an empty functor, used for debugging
	 *	@tparam n_size is matrix size
	 */
	template <const int n_size>
	class CDoNothing {
	public:
		/**
		 *	@brief an empty function
		 *	@param[in] context is a context (unused)
		 */
		static void Do(int UNUSED(context))
		{}
	};

public:
	/**
	 *	@brief the executive function
	 *
	 *	@param[in] n_size is the parameter to be decided over
	 *	@param[in,out] context is the operation context (task specific)
	 */
	template <template <const int> class CFunctor, class CContext>
	static __forceinline void Do(const size_t n_size, CContext context)
	{
		CMakeDecisionTree<_TySortedList, CContext, CFunctor>::Do(n_size, context);
	}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

	/**
	 *	@brief prints (unidented) decission tree pseudocode to stdout
	 *	@note This is only available if
	 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
	 */
	static void Debug()
	{
		CMakeDecisionTree<_TySortedList, int, CDoNothing>::Debug();
	}

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING
};

/**
 *	@brief fixed block size operations template class
 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
 *		matrices with known compile-time sizes
 */
template <class CBlockMatrixTypelistA, class CBlockMatrixTypelistB>
class CFixedBlockSize_BinaryBase {
public:
	typedef CBlockMatrixTypelistA _TyBlockMatrixTypelistA; /**< @brief list of block sizes, represented as Eigen::Matrix */
	typedef typename CTransformTypelist<_TyBlockMatrixTypelistA,
		CEigenToDimension>::_TyResult CDimsListA; /**< @brief list of block sizes as CCTSize2D */
	typedef typename CUniqueTypelist<CDimsListA>::_TyResult CDimsList_UniqA; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_UniqA,
		CTransformDimensionColumnsToSize>::_TyResult>::_TyResult CColumnWidthsListA; /**< @brief list of unique block column widths */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_UniqA,
		CTransformDimensionRowsToSize>::_TyResult>::_TyResult CRowHeightsListA; /**< @brief list of unique block row heights */

	typedef CBlockMatrixTypelistB _TyBlockMatrixTypelistB; /**< @brief list of block sizes, represented as Eigen::Matrix */
	typedef typename CTransformTypelist<_TyBlockMatrixTypelistB,
		CEigenToDimension>::_TyResult CDimsListB; /**< @brief list of block sizes as CCTSize2D */
	typedef typename CUniqueTypelist<CDimsListB>::_TyResult CDimsList_UniqB; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_UniqB,
		CTransformDimensionColumnsToSize>::_TyResult>::_TyResult CColumnWidthsListB; /**< @brief list of unique block column widths */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_UniqB,
		CTransformDimensionRowsToSize>::_TyResult>::_TyResult CRowHeightsListB; /**< @brief list of unique block row heights */
	// create typelists to generate the decision tree for block sizes

	typedef typename ubm_static_check::CStaticAssert<
		!CFindTypelistItem<CColumnWidthsListA, CCTSize<Eigen::Dynamic> >::b_result &&
		!CFindTypelistItem<CRowHeightsListA, CCTSize<Eigen::Dynamic> >::b_result &&
		!CFindTypelistItem<CColumnWidthsListB, CCTSize<Eigen::Dynamic> >::b_result &&
		!CFindTypelistItem<CRowHeightsListB, CCTSize<Eigen::Dynamic> >::b_result>::MATRIX_BLOCK_DIMENSIONS_MUST_BE_KNOWN_AT_COMPILE_TIME CAssert0; /**< @brief make sure that all the blocks have the sizes fixed and known at compile-time */
};

/**
 *	@brief fixed block size operations template class
 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
 *		matrices with known compile-time sizes
 */
template <class CBlockMatrixTypelist>
class CFixedBlockSize_UnaryBase {
public:
	typedef CBlockMatrixTypelist _TyBlockMatrixTypelist; /**< @brief list of block sizes, represented as Eigen::Matrix */
	typedef typename CTransformTypelist<_TyBlockMatrixTypelist,
		CEigenToDimension>::_TyResult CDimsList; /**< @brief list of block sizes as CCTSize2D */
	typedef typename CUniqueTypelist<CDimsList>::_TyResult CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_Uniq,
		CTransformDimensionColumnsToSize>::_TyResult>::_TyResult CColumnWidthsList; /**< @brief list of unique block column widths */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_Uniq,
		CTransformDimensionRowsToSize>::_TyResult>::_TyResult CRowHeightsList; /**< @brief list of unique block row heights */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_Uniq,
		CTransformDimensionToAreaSize>::_TyResult>::_TyResult CBlockAreasList; /**< @brief list of unique block areas (for elementwise ops) */
	// create typelists to generate the decision tree for block sizes

	typedef typename ubm_static_check::CStaticAssert<
		!CFindTypelistItem<CColumnWidthsList, CCTSize<Eigen::Dynamic> >::b_result &&
		!CFindTypelistItem<CRowHeightsList, CCTSize<Eigen::Dynamic> >::b_result>::MATRIX_BLOCK_DIMENSIONS_MUST_BE_KNOWN_AT_COMPILE_TIME CAssert0; /**< @brief make sure that all the blocks have the sizes fixed and known at compile-time */
};

/**
 *	@brief simple utility class for sorting block dimensions
 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
 *		matrices with known compile-time sizes
 */
template <class CBlockMatrixTypelist>
class CSortBlockDims {
public:
	typedef typename CSortTypelist<typename CUniqueTypelist<typename
		CTransformTypelist<CBlockMatrixTypelist, CEigenToDimension>::_TyResult>::_TyResult,
		CCompareSize2D_Less>::_TyResult _TyResult; /**< @brief unique sorted list of block matrix dimensions */
};

/**
 *	@brief simple utility for writing clean FBS loops
 *
 *	This would be used as follows:
 *	@code
 *	class CMy_FBS_Loop {
 *	public:
 *		struct TContext { // call context
 *			int params; // parameters to the loop, return type
 *
 *			inline TContext(int _params) // provide constructor
 *				:params(_params)
 *			{}
 *		};
 *
 *		template <const int n_fixed_size>
 *		struct TLoopBody { // loop body for a given fixed size
 *			static inline void Do(TContext t_ctx) // loop implementation
 *			{
 *				// implement the loop here
 *			}
 *		};
 *
 *	public:
 *		template <class CBlockMatrixTypelist> // list of possible block sizes (will fail on other sizes!)
 *		static inline void Loop(size_t n_size, int params) // loop interface; add parameters here
 *		{
 *			fbs_ut::CWrap<TLoopBody>::template
 *				In_ColumnWidth_DecisionTree<CBlockMatrixTypelist>(n_size, TContext(params));
 *			// wrap this loop in column width decision tree
 *		}
 *	};
 *	@endcode
 *
 *	@tparam CLoopType is a loop type; must be parametrizable by block size and
 *		must implement a <tt>Do(context)</tt> function
 *
 *	@note In case the loop body requires static context as well
 *		(such as block size list), use CWrap2 instead.
 */
template <template <const int n_size> class CLoopType>
class CWrap {
public:
	/**
	 *	@brief wraps a function in scalar size decision tree
	 *
	 *	@tparam CScalarSizeTypelist is list of all possible sizes (as CCTSize)
	 *	@tparam _TyContext is type of context that the wrapped function requires
	 *
	 *	@param[in] n_scalar_size is size at invokation (must be one of sizes in CScalarSizeTypelist)
	 *	@param[in] c is value of the context to be passed to the wrapped function
	 *
	 *	@note This executes operations at run time.
	 */
	template <class CScalarSizeTypelist, class _TyContext>
	static __forceinline void In_ScalarSize_DecisionTree(int n_scalar_size, _TyContext c)
	{
		typedef typename CSortTypelist<typename CUniqueTypelist<CScalarSizeTypelist>::_TyResult,
			fbs_ut::CCompareScalar_Less>::_TyResult _TySortedSizes;
		// sorted list of unique sizes

		CMakeDecisionTree<_TySortedSizes, _TyContext, CLoopType>::Do(n_scalar_size, c);
		// enter the decision tree
	}

	/**
	 *	@brief wraps a function in block width decision tree
	 *
	 *	@tparam CBlockMatrixTypelist is list of all possible matrix block sizes (as Eigen::Matrix with compile-time sizes)
	 *	@tparam _TyContext is type of context that the wrapped function requires
	 *
	 *	@param[in] n_column_width is size at invokation (must be one of column widths in CBlockMatrixTypelist)
	 *	@param[in] c is value of the context to be passed to the wrapped function
	 *
	 *	@note This executes operations at run time.
	 */
	template <class CBlockMatrixTypelist, class _TyContext>
	static __forceinline void In_ColumnWidth_DecisionTree(int n_column_width, _TyContext c)
	{
		typedef typename CTransformTypelist<typename CSortBlockDims<CBlockMatrixTypelist>::_TyResult,
			CTransformDimensionColumnsToSize>::_TyResult _TySortedColumnWidths;
		// sorted list of unique block column widths

		CMakeDecisionTree<_TySortedColumnWidths, _TyContext, CLoopType>::Do(n_column_width, c);
		// enter the decision tree
	}

	/**
	 *	@brief wraps a function in block height decision tree, given column width
	 *
	 *	@tparam CBlockMatrixTypelist is list of all possible matrix block sizes (as Eigen::Matrix with compile-time sizes)
	 *	@tparam n_column_width is column width, in elements (must equal one of column widths in CBlockMatrixTypelist)
	 *	@tparam _TyContext is type of context that the wrapped function requires
	 *
	 *	@param[in] n_row_height is size at invokation (must be one of row heights in CBlockMatrixTypelist)
	 *	@param[in] c is value of the context to be passed to the wrapped function
	 *
	 *	@note This executes operations at run time.
	 */
	template <class CBlockMatrixTypelist, const int n_column_width, class _TyContext>
	static __forceinline void In_RowHeight_DecisionTree_Given_ColumnWidth(int n_row_height, _TyContext c)
	{
		typedef typename CSortBlockDims<CBlockMatrixTypelist>::_TyResult _TyDimsList;
		typedef typename CTransformTypelist<_TyDimsList,
			CTransformDimensionRowsToSize>::_TyResult _TySortedRowHeights;
		// sorted list of unique block row heights

		typedef typename CFilterTypelist2<_TySortedRowHeights, fbs_ut::CCTSize<n_column_width>,
			fbs_ut::CHaveRowHeightForColumnWidth, _TyDimsList>::_TyResult _TySelectedRowHeights;
		// t_odo - apply column width filter here

		CMakeDecisionTree<_TySelectedRowHeights, _TyContext, CLoopType>::Do(n_row_height, c);
		// enter the decision tree
	}

	// todo - RowHeight
	// todo - conditional (ColumnWidth_Given_RowHeight and RowHeight_Given_ColumnWidth)
};

/**
 *	@brief simple utility for writing clean FBS loops
 *
 *	This would be used as follows:
 *	@code
 *	class CMy_FBS_Loop {
 *	public:
 *		struct TContext { // call context
 *			int params; // parameters to the loop, return type
 *
 *			inline TContext(int _params) // provide constructor
 *				:params(_params)
 *			{}
 *		};
 *
 *		template <const int n_fixed_size, class CBlockSizes> // CBlockSizes filled from the CContext2 argument
 *		struct TLoopBody { // loop body for a given fixed size
 *			static inline void Do(TContext t_ctx) // loop implementation
 *			{
 *				// implement the loop here
 *			}
 *		};
 *
 *	public:
 *		template <class CBlockMatrixTypelist> // list of possible block sizes (will fail on other sizes!)
 *		static inline void Loop(size_t n_size, int params) // loop interface; add parameters here
 *		{
 *			fbs_ut::CWrap2<TLoopBody, CBlockMatrixTypelist>::template
 *				In_ColumnWidth_DecisionTree<CBlockMatrixTypelist>(n_size, TContext(params));
 *			// wrap this loop in column width decision tree,
 *			// and also pass CBlockMatrixTypelist to TLoopBody
 *		}
 *	};
 *	@endcode
 *
 *	@tparam CLoopType is a loop type; must be parametrizable by block size and
 *		static context and must implement a <tt>Do(context)</tt> function
 *	@tparam CContext2 is a type of context, passed as a second template parameter for the operation
 *
 *	@note In case the loop body does not require the static context, use CWrap instead.
 */
template <template <const int, class> class CLoopType, class CContext2>
class CWrap2 {
public:
	/**
	 *	@brief wraps a function in scalar size decision tree
	 *
	 *	@tparam CScalarSizeTypelist is list of all possible sizes (as CCTSize)
	 *	@tparam _TyContext is type of context that the wrapped function requires
	 *
	 *	@param[in] n_scalar_size is size at invokation (must be one of sizes in CScalarSizeTypelist)
	 *	@param[in] c is value of the context to be passed to the wrapped function
	 *
	 *	@note This executes operations at run time.
	 */
	template <class CScalarSizeTypelist, class _TyContext>
	static __forceinline void In_ScalarSize_DecisionTree(int n_scalar_size, _TyContext c)
	{
		typedef typename CSortTypelist<typename CUniqueTypelist<CScalarSizeTypelist>::_TyResult,
			fbs_ut::CCompareScalar_Less>::_TyResult _TySortedSizes;
		// sorted list of unique sizes

		CMakeDecisionTree2<_TySortedSizes, _TyContext, CContext2, CLoopType>::Do(n_scalar_size, c);
		// enter the decision tree
	}

	/**
	 *	@brief wraps a function in block width decision tree
	 *
	 *	@tparam CBlockMatrixTypelist is list of all possible matrix block sizes (as Eigen::Matrix with compile-time sizes)
	 *	@tparam _TyContext is type of context that the wrapped function requires
	 *
	 *	@param[in] n_column_width is size at invokation (must be one of column widths in CBlockMatrixTypelist)
	 *	@param[in] c is value of the context to be passed to the wrapped function
	 *
	 *	@note This executes operations at run time.
	 */
	template <class CBlockMatrixTypelist, class _TyContext>
	static __forceinline void In_ColumnWidth_DecisionTree(int n_column_width, _TyContext c)
	{
		typedef typename CTransformTypelist<typename CSortBlockDims<CBlockMatrixTypelist>::_TyResult,
			CTransformDimensionColumnsToSize>::_TyResult _TySortedColumnWidths;
		// sorted list of unique block column widths

		CMakeDecisionTree2<_TySortedColumnWidths, _TyContext, CContext2, CLoopType>::Do(n_column_width, c);
		// enter the decision tree
	}

	/**
	 *	@brief wraps a function in block height decision tree, given column width
	 *
	 *	@tparam CBlockMatrixTypelist is list of all possible matrix block sizes (as Eigen::Matrix with compile-time sizes)
	 *	@tparam n_column_width is column width, in elements (must equal one of column widths in CBlockMatrixTypelist)
	 *	@tparam _TyContext is type of context that the wrapped function requires
	 *
	 *	@param[in] n_row_height is size at invokation (must be one of row heights in CBlockMatrixTypelist)
	 *	@param[in] c is value of the context to be passed to the wrapped function
	 *
	 *	@note This executes operations at run time.
	 */
	template <class CBlockMatrixTypelist, const int n_column_width, class _TyContext>
	static __forceinline void In_RowHeight_DecisionTree_Given_ColumnWidth(int n_row_height, _TyContext c)
	{
		typedef typename CSortBlockDims<CBlockMatrixTypelist>::_TyResult _TyDimsList;
		typedef typename CTransformTypelist<_TyDimsList,
			CTransformDimensionRowsToSize>::_TyResult _TySortedRowHeights;
		// sorted list of unique block row heights

		typedef typename CFilterTypelist2<_TySortedRowHeights, fbs_ut::CCTSize<n_column_width>,
			fbs_ut::CHaveRowHeightForColumnWidth, _TyDimsList>::_TyResult _TySelectedRowHeights;
		// t_odo - apply column width filter here

		CMakeDecisionTree2<_TySelectedRowHeights, _TyContext, CContext2, CLoopType>::Do(n_row_height, c);
		// enter the decision tree
	}

	// todo - RowHeight
	// todo - conditional (ColumnWidth_Given_RowHeight and RowHeight_Given_ColumnWidth)
};

/**
 *	@brief helper object for copying typelist data to an array
 *	@tparam T is data type of the output array elements
 */
template <class T = int>
class CCopyCTSizesToArray {
protected:
	T *m_p_dest; /**< @brief pointer to the next output element */

public:
	/**
	 *	@brief default constructor; initializes the output pointer
	 *	@param[in] p_dest is pointer to the destination array (must be allocated to sufficient number of elements)
	 */
	CCopyCTSizesToArray(T *p_dest)
		:m_p_dest(p_dest)
	{}

	/**
	 *	@brief function operator; writes a single number in the output array
	 *	@tparam CDimension is a specialization of CCTSize (specifies the value to be written)
	 */
	template <class CDimension>
	void operator ()()
	{
		*m_p_dest = CDimension::n_size;
		++ m_p_dest;
	}

	/**
	 *	@brief cast to pointer
	 *	@return Returns pointer to the next output element.
	 */
	operator const T *() const
	{
		return m_p_dest;
	}
};

/**
 *	@brief copies typelist of CCTSize to an array
 *
 *	@tparam CSizeList is typelist of CCTSize specializations
 *	@tparam T is data type of the output array elements
 *	@tparam n_size is size of the output array (must match the length of the list)
 *
 *	@param[out] p_dest is array of integers, filled with the data from the input typelist
 *
 *	@note This executes operations at run time.
 */
template <class CSizeList, class T, const size_t n_size>
void Copy_CTSizes_to_Array(T (&p_dest)[n_size])
{
	const T *p_end = CTypelistForEach<CSizeList, CCopyCTSizesToArray<T> >::Run(CCopyCTSizesToArray<T>(p_dest));
	_ASSERTE(p_end == p_dest + n_size); // make sure the size was correct
}

/**
 *	@brief copies typelist of CCTSize to an array
 *
 *	@tparam CSizeList is typelist of CCTSize specializations
 *	@tparam T is data type of the output array elements
 *
 *	@param[out] p_dest is array of integers, filled with the data from the input typelist
 *	@param[in] n_size is size of the output array (must match the length of the list; used only in debug)
 *
 *	@note This executes operations at run time.
 */
template <class CSizeList, class T>
void Copy_CTSizes_to_Array(T *p_dest, const size_t UNUSED(n_size))
{
	const T *p_end = CTypelistForEach<CSizeList, CCopyCTSizesToArray<T> >::Run(CCopyCTSizesToArray<T>(p_dest));
	_ASSERTE(p_end == p_dest + n_size); // make sure the size was correct
}

} // ~fbs_ut

#include "slam/Self.h"
// it is recommended to use DECLARE_SELF and _TySelf in conjunction with fbs_ut::CWrap

#endif // !__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_OPS_UTILITIES_INCLUDED
