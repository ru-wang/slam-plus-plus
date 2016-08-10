/*
								+----------------------------------+
								|                                  |
								|      ***  BA Marginals  ***      |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2016  |
								|                                  |
								|          BAMarginals.h           |
								|                                  |
								+----------------------------------+
*/

/**
 *	@file include/ba_marginals_example/BAMarginals.h
 *	@brief helper classes for covariance recovery in Schur-complemented systems
 *	@author -tHE SWINe-
 *	@date 2016-02-15
 */

#ifndef __BA_MARGINALS_INCLUDED
#define __BA_MARGINALS_INCLUDED

#include "slam/Marginals.h"

//template <class CDestMatrix> // everything depends on Derived0, it is ok to be a subclass like this

/**
 *	@brief multiply add operation with a sparse matrix and a block vector
 */
class CBlockVectorMAD_Impl {
public:
	struct TInnerContext {
		double *p_dest;
		const double *p_block;

		inline TInnerContext(double *_p_dest, const double *_p_block)
			:p_dest(_p_dest), p_block(_p_block)
		{}
	};

	template <const int n_row_height, class CColumnWidth>
	class CInnerLoop { // not dependent on CBlockSizeList; don't make it an inner class of COuterLoop
	public:
		enum {
			n_column_width = CColumnWidth::n_size
		};

		typedef typename CUberBlockMatrix::CMakeMatrixRef<n_column_width, n_column_width>::_Ty _TyDestMap;
		typedef typename CUberBlockMatrix::CMakeMatrixRef<n_row_height, n_column_width>::_TyConst _TyBlockMap;

	public:
		static inline void Do(TInnerContext t_ctx)
		{
			_TyBlockMap fbs_block(t_ctx.p_block);
			_TyDestMap(t_ctx.p_dest) += fbs_block.transpose().lazyProduct(fbs_block); // mad
		}
	};

	struct TOuterContext {
		double *p_dest;
		const CUberBlockMatrix &r_block_vector;

		inline TOuterContext(double *_p_dest, const CUberBlockMatrix &_r_block_vector)
			:p_dest(_p_dest), r_block_vector(_r_block_vector)
		{}
	};

	template <const int n_column_width, class CBlockSizeList>
	class COuterLoop {
	public:
		static inline void Do(TOuterContext t_ctx)
		{
			_ASSERTE(n_column_width == t_ctx.r_block_vector.n_BlockColumn_Column_Num(0));
			for(size_t i = 0, n = t_ctx.r_block_vector.n_BlockColumn_Block_Num(0); i < n; ++ i) {
				CUberBlockMatrix::_TyMatrixXdRef t_block =
					const_cast<CUberBlockMatrix&>(t_ctx.r_block_vector).t_Block_AtColumn(0, i);
				// hack - need to cast, const_ref does not expose its pointer to data

				size_t n_row_height = t_block.rows();
				_ASSERTE(n_column_width == t_block.cols());

				fbs_ut::CWrap2<CInnerLoop, fbs_ut::CCTSize<n_column_width> >::template
					In_RowHeight_DecisionTree_Given_ColumnWidth<CBlockSizeList,
					n_column_width>(int(n_row_height), TInnerContext(t_ctx.p_dest, t_block.data()/*&t_block(0, 0)*/));
				// mad
			}
			// for each nnz block
		}
	};

/*public:
	template <class CBlockSizeList>
	static void BlockVector_PreMultiplyWithSelfTranspose_Add_FBS(double *p_dest,
		const CUberBlockMatrix &r_block_vector) // g++ is unable to reference TOuterLoop from the outside for some reason
	{
		fbs_ut::CWrap2<COuterLoop, CBlockSizeList>::template In_ColumnWidth_DecisionTree<CBlockSizeList>(
			r_block_vector.n_Column_Num(), TOuterContext(p_dest, r_block_vector));
		// decide over vector width
	}*/
};

template <class CBlockSizeList, class CDestMatrix>
inline void BlockVector_PreMultiplyWithSelfTranspose_Add_FBS(CDestMatrix &r_dest,
	const CUberBlockMatrix &r_block_vector)
{
	_ASSERTE(r_block_vector.n_BlockColumn_Num() == 1); // it is a block column-vector
	_ASSERTE(r_block_vector.n_Column_Num() == r_dest.rows()); // the number of columns in the block vector matches the size of the destination
	_ASSERTE(r_dest.rows() == r_dest.cols()); // the output is square
	_ASSERTE(!(CDestMatrix::Flags & Eigen::RowMajor)); // must be column-major otherwise the conversion to pointer strips data
	//_ASSERTE(CDestMatrix::MaxColsAtCompileTime == Eigen::Dynamic ||
	//	CDestMatrix::MaxColsAtCompileTime == r_dest.cols() ||
	//	r_dest.cols() == 1); // the stride must be tight or there is a single col and it does not matter
	_ASSERTE(r_dest.cols() <= 1 || &r_dest(0, 1) == &r_dest(0, 0) + r_dest.rows()); // the stride must be tight or there is a single col and it does not matter

	// do not zero r_dest

	//CBlockVectorMAD_Impl::template BlockVector_PreMultiplyWithSelfTranspose_Add_FBS<CBlockSizeList>(r_dest.data(), r_block_vector);

	fbs_ut::CWrap2<CBlockVectorMAD_Impl::COuterLoop, CBlockSizeList>::template
		In_ColumnWidth_DecisionTree<CBlockSizeList>(int(r_block_vector.n_Column_Num()),
		CBlockVectorMAD_Impl::TOuterContext(r_dest.data(), r_block_vector));
	// decide over vector width
}

void Calculate_UpperTriangularTransposeSolve_Bases(CUberBlockMatrix &S_bases,
	const CUberBlockMatrix &S, /*const*/ cs *p_St, size_t n_column,
	Eigen::MatrixXd &r_workspace, std::vector<size_t> &r_workspace1)
{
	_ASSERTE(S.b_EqualLayout(S_bases)); // will yield a matrix with the same sparsity structure
	//S.CopyLayoutTo(S_bases);

	_ASSERTE(n_column < S.n_BlockColumn_Num()); // make sure the column is inside the matrix

	_ASSERTE(p_St->m == S.n_BlockRow_Num() && p_St->n == S.n_BlockColumn_Num()); // should be the same matrix

	_ASSERTE(sizeof(csi) == sizeof(size_t));
	r_workspace1.resize(2 * p_St->n);
	// alloc workspace

	{
		//cs *p_St = cs_transpose(p_S, 0);
		// need p_St

		/*csi p_col[2] = {0, 1};
		csi n_row = 0;
		cs B;
		B.m = p_St->m;
		B.n = 1;
		B.p = p_col;
		B.i = &n_row;
		B.x = 0;
		B.nzmax = 1;
		B.nz = -1;*/
		// prepare a single entry CSC matrix

		//Eigen::Matrix<double, Eigen::Dynamic, 6> S_dense_basis(S.n_Row_Num(), 6); // todo - implement matrix solving and try to make this row major
		// unit basis matrix

		//for(size_t n_column = 0, n = S.n_BlockColumn_Num(); n_column < n; ++ n_column) { // t_odo - do this in parallel (probably explicit matrix reduction rather than locking)
		{ //size_t n_column = n_column;
			size_t w = S.n_BlockColumn_Column_Num(n_column); // t_odo - FBS it
			//_ASSERTE(w == 6);
			//n_row = n_column;

			//size_t n_first_dep_col = cs_reach(p_St, &B, 0, (csi*)&r_workspace1.front(), 0); // modifies p_St but then puts it back
			size_t n_first_dep_col = cs_dfs(n_column, p_St, p_St->n, (csi*)&r_workspace1.front(),
				(csi*)&r_workspace1.front() + p_St->n, 0); // no need for B // todo - reimplement this directly on block matrices, try to avoid needing the transpose
			size_t n_dep_col_num = p_St->n - n_first_dep_col;
			const size_t *p_dep_col = &r_workspace1[n_first_dep_col];
			for(size_t j = 0; j < n_dep_col_num; ++ j)
				CS_MARK(p_St->p, p_dep_col[j]); // restore col. pointers after calling cs_dfs
			// get the list of columns of S that affect block U_Dinv_{n_column, *}

			r_workspace.resize(S.n_Row_Num(), w);
			// alloc workspace

			Eigen::MatrixXd &S_dense_basis = r_workspace; // just rename
			//Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Aligned>
			//	S_dense_basis(r_workspace.data(), S.n_Row_Num(), w);
			// map the workspace as (column-wise) fixed-size matrix

			S_dense_basis.setZero();
			S_dense_basis.middleRows(S_bases.n_BlockColumn_Base(n_column), w).setIdentity();
			// create a vector of zeros

			for(size_t c = 0; c < w; ++ c) {
				//S.UpperTriangularTranspose_Solve(&U_Dinv_i_permd.col(c)(0), U_Dinv_i_permd.rows(), p_dep_col, n_dep_col_num);
				S.UpperTriangularTranspose_Solve/*_FBS<SC_BlockSizes>*/(&S_dense_basis.col(c)(0),
					S_dense_basis.rows(), p_dep_col, n_dep_col_num);
			}
			// sparse sparse UTTSolve

			for(size_t j = 0; j < n_dep_col_num; ++ j) {
				size_t n_row = p_dep_col[j];
				size_t y = S_bases.n_BlockColumn_Base(n_row);
				size_t h = S_bases.n_BlockColumn_Column_Num(n_row); // those are rows but S_bases is symmetric
#ifdef _DEBUG
				if(j > 0) {
					const size_t r_prev = p_dep_col[j - 1];
					size_t y_prev = S_bases.n_BlockColumn_Base(r_prev);
					size_t h_prev = S_bases.n_BlockColumn_Column_Num(r_prev);
					size_t e_prev = y_prev + h_prev;
					_ASSERTE(S_dense_basis.middleRows(e_prev, y - e_prev).squaredNorm() == 0); // make sure there are zeros between consecutive (nonadjacent) blocks
				} else if(y > 0)
					_ASSERTE(S_dense_basis.topRows(y).squaredNorm() == 0); // make sure there are zeros above the first block
				if(j + 1 == n_dep_col_num)
					_ASSERTE(S_dense_basis.bottomRows(S_dense_basis.rows() - (y + h)).squaredNorm() == 0); // make sure there are zeros till the end
				// make sure that there are only zeros in between the elements
#endif // _DEBUG
				//_ASSERTE(h == 6);
				//_ASSERTE(S_bases.n_BlockColumn_Column_Num(n_row) == 6); // t_odo - FBS it
				S_bases.t_GetBlock_Log(n_row, n_column, h, w) =
					S_dense_basis.middleRows(S_bases.n_BlockColumn_Base(n_row), h); // this is transposed (transpose the block as well?): each row is a single basis; this only works if the structure of S is symmetric
			}
			// sparse fill the bases matrix
		}

		//S_bases.Rasterize("S_bases.tga", 3); // ...
	}
}

/**
 *	@brief FBS implementation for the upper triangular transpose solve of the sparse bases matrix
 */
class CUTTSolve_Bases_Impl {
public:
	struct TInnerContext {
		CUberBlockMatrix &r_dest;
		const size_t n_column;
		const Eigen::MatrixXd &r_src;
		const size_t n_row;

		inline TInnerContext(CUberBlockMatrix &_r_dest, size_t _n_column,
			const Eigen::MatrixXd &_r_src, size_t _n_row)
			:r_dest(_r_dest), n_column(_n_column), r_src(_r_src),
			n_row(_n_row)
		{}
	};

	template <const int n_row_height, class CColumnWidth>
	class CInnerLoop {
	public:
		enum {
			n_column_width = CColumnWidth::n_size
		};

	public:
		static inline void Do(TInnerContext t_ctx)
		{
			_ASSERTE(t_ctx.r_src.rows() == t_ctx.r_dest.n_Row_Num());
			_ASSERTE(n_column_width == t_ctx.r_dest.n_BlockColumn_Column_Num(t_ctx.n_column));

			Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, n_column_width>, Eigen::Aligned>
				S_dense_basis(t_ctx.r_src.data(), t_ctx.r_src.rows(), n_column_width);
			// map the source as (column-wise) fixed-size matrix

			Eigen::Map<Eigen::Matrix<double, n_row_height,
				n_column_width> > dest_block(t_ctx.r_dest.p_GetBlock_Log(t_ctx.n_row,
				t_ctx.n_column, n_row_height, n_column_width, true, false));
			dest_block = S_dense_basis.template middleRows<n_row_height>(t_ctx.r_dest.n_BlockColumn_Base(t_ctx.n_row));
		}
	};

	struct TOuterContext {
		CUberBlockMatrix &r_S_bases;
		const size_t n_column;
		const CUberBlockMatrix &r_S;
		Eigen::MatrixXd &r_workspace;
		const size_t *p_dep_col;
		const size_t n_dep_num;

		inline TOuterContext(CUberBlockMatrix &_r_S_bases, size_t _n_column,
			const CUberBlockMatrix &_r_S, Eigen::MatrixXd &_r_workspace,
			const size_t *_p_dep_col, size_t _n_dep_num)
			:r_S_bases(_r_S_bases), n_column(_n_column), r_S(_r_S),
			r_workspace(_r_workspace), p_dep_col(_p_dep_col), n_dep_num(_n_dep_num)
		{}
	};

	template <const int n_column_width, class CBlockSizeList>
	class COuterLoop {
	public:
		static inline void Do(TOuterContext t_ctx)
		{
			t_ctx.r_workspace.resize(t_ctx.r_S.n_Row_Num(), n_column_width);
			// alloc workspace

			Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, n_column_width>, Eigen::Aligned>
				S_dense_basis(t_ctx.r_workspace.data(), t_ctx.r_S.n_Row_Num(), n_column_width);
			// map the workspace as (column-wise) fixed-size matrix

			S_dense_basis.setZero();
			S_dense_basis.template middleRows<n_column_width>(t_ctx.r_S.n_BlockColumn_Base(t_ctx.n_column)).setIdentity();
			// create a vector of zeros

			for(size_t c = 0; c < n_column_width; ++ c) { // todo - make a version of UpperTriangularTranspose_Solve_FBS for vectors
				t_ctx.r_S.UpperTriangularTranspose_Solve_FBS<CBlockSizeList>(&S_dense_basis.col(c)(0),
					S_dense_basis.rows(), t_ctx.p_dep_col, t_ctx.n_dep_num);
			}
			// sparse sparse UTTSolve

			for(size_t j = 0; j < t_ctx.n_dep_num; ++ j) {
				size_t n_row = t_ctx.p_dep_col[j];
				size_t h = t_ctx.r_S.n_BlockColumn_Column_Num(n_row); // those are rows but S_bases is symmetric
#ifdef _DEBUG
				size_t y = t_ctx.r_S.n_BlockColumn_Base(n_row);
				if(j > 0) {
					const size_t r_prev = t_ctx.p_dep_col[j - 1];
					size_t y_prev = t_ctx.r_S.n_BlockColumn_Base(r_prev);
					size_t h_prev = t_ctx.r_S.n_BlockColumn_Column_Num(r_prev);
					size_t e_prev = y_prev + h_prev;
					_ASSERTE(S_dense_basis.middleRows(e_prev, y - e_prev).squaredNorm() == 0); // make sure there are zeros between consecutive (nonadjacent) blocks
				} else if(y > 0)
					_ASSERTE(S_dense_basis.topRows(y).squaredNorm() == 0); // make sure there are zeros above the first block
				if(j + 1 == t_ctx.n_dep_num)
					_ASSERTE(S_dense_basis.bottomRows(S_dense_basis.rows() - (y + h)).squaredNorm() == 0); // make sure there are zeros till the end
				// make sure that there are only zeros in between the elements
#endif // _DEBUG
				//_ASSERTE(h == 6);
				//_ASSERTE(S_bases.n_BlockColumn_Column_Num(n_row) == 6); // t_odo - FBS it
				//S_bases.t_GetBlock_Log(n_row, n_column, h, w) =
				//	S_dense_basis.middleRows<6>(S_bases.n_BlockColumn_Base(n_row)); // this is transposed (transpose the block as well?): each row is a single basis; this only works if the structure of S is symmetric

				fbs_ut::CWrap2<CInnerLoop, fbs_ut::CCTSize<n_column_width> >::template
					In_RowHeight_DecisionTree_Given_ColumnWidth<CBlockSizeList,
					n_column_width>(int(h), TInnerContext(t_ctx.r_S_bases,
					t_ctx.n_column, t_ctx.r_workspace, n_row));
				// use SSE to copy stuff around
			}
			// sparse fill the bases matrix
		}
	};
};

template <class CBlockSizeList>
void Calculate_UpperTriangularTransposeSolve_Bases_FBS(CUberBlockMatrix &S_bases,
	const CUberBlockMatrix &S, /*const*/ cs *p_St, size_t n_column,
	Eigen::MatrixXd &r_workspace, std::vector<size_t> &r_workspace1)
{
	_ASSERTE(S.b_EqualLayout(S_bases)); // will yield a matrix with the same sparsity structure
	//S.CopyLayoutTo(S_bases);

	_ASSERTE(n_column < S.n_BlockColumn_Num()); // make sure the column is inside the matrix

	_ASSERTE(p_St->m == S.n_BlockRow_Num() && p_St->n == S.n_BlockColumn_Num()); // should be the same matrix

	_ASSERTE(sizeof(csi) == sizeof(size_t));
	r_workspace1.resize(2 * p_St->n);
	// alloc workspace

	size_t w = S.n_BlockColumn_Column_Num(n_column); // t_odo - FBS it
	//_ASSERTE(w == 6);
	//n_row = n_column;

	//size_t n_first_dep_col = cs_reach(p_St, &B, 0, (csi*)&r_workspace1.front(), 0); // modifies p_St but then puts it back
	size_t n_first_dep_col = cs_dfs(n_column, p_St, p_St->n, (csi*)&r_workspace1.front(),
		(csi*)&r_workspace1.front() + p_St->n, 0); // no need for B // todo - reimplement this directly on block matrices, try to avoid needing the transpose
	size_t n_dep_col_num = p_St->n - n_first_dep_col;
	const size_t *p_dep_col = &r_workspace1[n_first_dep_col];
	for(size_t j = 0; j < n_dep_col_num; ++ j)
		CS_MARK(p_St->p, p_dep_col[j]); // restore col. pointers after calling cs_dfs
	// get the list of columns of S that affect block U_Dinv_{n_column, *}
	// note that this is FBS-independent

	fbs_ut::CWrap2<CUTTSolve_Bases_Impl::COuterLoop, CBlockSizeList>::template
		In_ColumnWidth_DecisionTree<CBlockSizeList>(int(w),
		CUTTSolve_Bases_Impl::TOuterContext(S_bases,
		n_column, S, r_workspace, p_dep_col, n_dep_col_num));
}

#endif // !__BA_MARGINALS_INCLUDED
