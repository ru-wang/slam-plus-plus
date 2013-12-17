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
namespace __fbs_ut {

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
		// converting to size_t disables g++ warning: comparison between 'enum CCTSize<X>::<anonymous>' and 'enum CCTSize<Y>::<anonymous>'
	};
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
		// converting to size_t disables g++ warning: comparison between 'enum CCTSize<X>::<anonymous>' and 'enum CCTSize<Y>::<anonymous>'
	};
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
	typedef Eigen::Matrix<double, CDimensionType::n_row_num, CDimensionType::n_column_num> _TyResult;
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
	typedef CCTSize2D<_T1::n_size, _T2::n_size> _TyResult; /**< @brief the resulting Eigen matrix type */
};

/**
 *	@brief calculates a typelist of block sizes after the PreMultiplyWithSelfTranspose()
 *		operation (the result can be found in the _TyResult type)
 *	@tparam CBlockMatrixTypelist is typelist, containing Eigen
 *		matrices with known compile-time sizes
 */
template <class CBlockMatrixTypelist>
class CBlockMatrixTypesAfterPreMultiplyWithSelfTranspose {
protected:
	typedef CBlockMatrixTypelist _TyBlockMatrixTypelist; /**< @brief list of block sizes, represented as Eigen::Matrix */
	typedef typename CTransformTypelist<_TyBlockMatrixTypelist,
		CEigenToDimension>::_TyResult CDimsList; /**< @brief list of block sizes as CCTSize2D */
	typedef typename CUniqueTypelist<CDimsList>::_TyResult CDimsList_Uniq; /**< @brief list of block sizes as CCTSize2D (duplicate records removed) */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_Uniq,
		CTransformDimensionRowsToSize>::_TyResult>::_TyResult CRowHeightsList; /**< @brief list of unique block row heights */

public:
	typedef typename CUniqueTypelist<typename CCarthesianProductTypelist<CRowHeightsList,
		CRowHeightsList, CMakeMatrixSizeType>::_TyResult>::_TyResult _TyResult; /**< @brief the resulting typelist */
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
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList_Uniq,
		CTransformDimensionRowsToSize>::_TyResult>::_TyResult CRowHeightsList; /**< @brief list of unique block row heights */

public:
	typedef typename CUniqueTypelist<typename CCarthesianProductTypelist<CRowHeightsList,
		CRowHeightsList, CMakeCTSize2DType>::_TyResult>::_TyResult _TyResult; /**< @brief the resulting typelist */
};

/**
 *	@brief prints typelist containing Eigen matrices to stdout (debugging utility)
 *
 *	Use as follows:
 *	@code
 *	typedef MakeTypelist_Safe((Eigen::Matrix2d, Eigen::Matrix3d)) matrices;
 *	__fbs_ut::CDumpBlockMatrixTypelist<matrices>::Print();
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
 *	typedef MakeTypelist_Safe((__fbs_ut::CCTSize<1>, __fbs_ut::CCTSize<2>)) sizes;
 *	__fbs_ut::CDumpCTSizeTypelist<sizes>::Print();
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
class CMakeDecissionTree {
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
	class CDecissionTree {
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
				CDecissionTree<_CList, _CContext, _CFunctor, n_begin, n_half>::Do(n_size, context);
			else
				CDecissionTree<_CList, _CContext, _CFunctor, n_pivot, n_rest>::Do(n_size, context);
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
			CDecissionTree<_CList, _CContext, _CFunctor, n_begin, n_half>::Debug();
			printf("} else {\n");
			CDecissionTree<_CList, _CContext, _CFunctor, n_pivot, n_rest>::Debug();
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
	class CDecissionTree<_CList, _CContext, _CFunctor, n_begin, 1> {
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
		CDecissionTree<CList, CContext, CFunctor, 0,
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
		CDecissionTree<CList, CContext, CFunctor, 0, CTypelistLength<CList>::n_result>::Debug();
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
		CMakeDecissionTree<_TySortedList, CContext, CFunctor>::Do(n_size, context);
	}

#ifdef __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING

	/**
	 *	@brief prints (unidented) decission tree pseudocode to stdout
	 *	@note This is only available if
	 *		__UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_DEBUGGING is defined.
	 */
	static void Debug()
	{
		CMakeDecissionTree<_TySortedList, int, CDoNothing>::Debug();
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

	typedef __ubm_static_check::CStaticAssert<sizeof(CMATRIX_BLOCK_DIMENSIONS_MUST_BE_KNOWN_AT_COMPILE_TIME<
		!CFindTypelistItem<CColumnWidthsListA, CCTSize<Eigen::Dynamic> >::b_result &&
		!CFindTypelistItem<CRowHeightsListA, CCTSize<Eigen::Dynamic> >::b_result &&
		!CFindTypelistItem<CColumnWidthsListB, CCTSize<Eigen::Dynamic> >::b_result &&
		!CFindTypelistItem<CRowHeightsListB, CCTSize<Eigen::Dynamic> >::b_result>)> _TyAssert0; /**< @brief make sure that all the blocks have the sizes fixed and known at compile-time */
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

	typedef __ubm_static_check::CStaticAssert<sizeof(CMATRIX_BLOCK_DIMENSIONS_MUST_BE_KNOWN_AT_COMPILE_TIME<
		!CFindTypelistItem<CColumnWidthsList, CCTSize<Eigen::Dynamic> >::b_result &&
		!CFindTypelistItem<CRowHeightsList, CCTSize<Eigen::Dynamic> >::b_result>)> _TyAssert0; /**< @brief make sure that all the blocks have the sizes fixed and known at compile-time */
};

} // ~__fbs_ut

#endif // __UBER_BLOCK_MATRIX_FIXED_BLOCK_SIZE_OPS_UTILITIES_INCLUDED
