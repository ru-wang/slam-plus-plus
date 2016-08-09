/*
								+-----------------------------------+
								|                                   |
								| ***  Lambda nonlinear solver  *** |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2012  |
								|                                   |
								|     NonlinearSolver_Lambda.h      |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
#define __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED

/**
 *	@file include/slam/NonlinearSolver_Lambda.h
 *	@brief nonlinear blocky solver working above the lambda matrix
 *	@author -tHE SWINe-
 *	@date 2012-09-13
 */

#include "slam/FlatSystem.h"
#include "slam/LinearSolver_Schur.h"
#include "slam/OrderingMagic.h"
#include "slam/IncrementalPolicy.h"
#include "slam/Marginals.h"

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
 *	@brief enables writes of chi2 errors at each step
 */
//#define __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_DUMP_ICRA2013_ANIMATION_DATA
 *	@brief dump ICRA 2013 slam race data
 */
//#define __NONLINEAR_SOLVER_LAMBDA_DUMP_ICRA2013_ANIMATION_DATA

/**
 *	@def __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
 *	@brief dump RSS 2013 matrix animation data
 */
//#define __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

/**
 *	@brief utilities for lambda solvers
 */
namespace lambda_utils {

/**
 *	@brief static assertion helper
 *	@brief b_expression is expression being asserted
 */
template <bool b_expression>
class CReductionPlanAssert {
public:
	typedef void BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST; /**< @brief static assertion tag */
};

/**
 *	@brief static assertion helper (specialization for assertion failed)
 */
template <>
class CReductionPlanAssert<false> {};

/**
 *	@brief calculates std::set memory allocation size estimate
 *	@tparam P is set payload type
 *	@param[in] r_set is the set to estimate memory size of
 *	@return Returns the approximate size of the given set, in bytes.
 */
template <class P>
static size_t n_Set_Allocation_Size(const std::set<P> &r_set)
{
	return (sizeof(P) + sizeof(int) + 2 * sizeof(void*)) * // size of a node (payload, red/black, children)
		((r_set.size() * 2) / 3) + sizeof(std::set<P>); // number of leafs + other nodes + size of the struct
}

/**
 *	@brief calculates std::map memory allocation size estimate
 *
 *	@tparam K is map key type
 *	@tparam P is map payload type
 *
 *	@param[in] r_map is the map to estimate memory size of
 *
 *	@return Returns the approximate size of the given map, in bytes.
 */
template <class K, class P>
static size_t n_Map_Allocation_Size(const std::map<K, P> &r_map)
{
	return (sizeof(K) + sizeof(P) + sizeof(int) + 2 * sizeof(void*)) * // size of a node (key, payload, red/black, children)
		((r_map.size() * 2) / 3) + sizeof(std::map<K, P>); // number of leafs + other nodes + size of the struct
}

// todo - consider not storing the reductions in a pool, instead store them in the map and keep a vector of pointers for parallel processing?

/**
 *	@brief right hand side vector reduction plan
 *
 *	This takes care of summing up right hand side (residual) vector contributions from different edges.
 *	In the v1 reduction plan this was done by the vertex class, each vertex had to keep a list of referencing
 *	edges and each edge contained a vector for the r.h.s. contribution. This does not improve efficiency but
 *	seems like a slightly cleaner solution.
 *
 *	@tparam CDimsList is list of block sizes in lambda matrix
 */
template <class CDimsList> // todo - need clear!
class CVectorReductionPlan { // t_odo - could we possibly move this to nl solver lambda.h? // seems like we did
public:
	typedef typename CTransformTypelist<CDimsList, fbs_ut::CEigenToDimension>::_TyResult _TyBlockSizeList; /**< @brief list of block sizes in lambda */
	typedef typename CUniqueTypelist<typename CFilterTypelist1<_TyBlockSizeList,
		fbs_ut::CIsSquare>::_TyResult>::_TyResult _TyVertDimsList; /**< @brief list of square block sizes (corresponding to vertices) */
	typedef typename CTransformTypelist<_TyVertDimsList,
		fbs_ut::CTransformDimensionColumnsToSize>::_TyResult _TyDimensionList; /**< @brief list of vertex dimensions (as fbs_ut::CCTSize, not 2D) */

	/**
	 *	@brief reduction plan parameters, stored as enum
	 */
	enum {
		n_reductor_num = CTypelistLength<_TyDimensionList>::n_result, /**< @brief number of different reduction sizes */
		n_pool_page_size = 4096, /**< @brief reduction pool page size */
		n_pool_memory_align = 0 /**< @brief reduction pool memory alignment */ // see no benefit in alignment right now, this stores the TReduction elements, those do not need to be aligned
	};

protected:
	/**
	 *	@brief reduction description
	 *	@note This does not store reduction dimension; that is given by the index of the pool that this is found in.
	 */
	struct TReduction {
		size_t n_offset; /**< @brief offset in the right hand side vector */
		std::vector<const double*> src_list; /**< @brief list of sources */
	};

	typedef forward_allocated_pool<TReduction,
		n_pool_page_size, n_pool_memory_align> _TyPool; /**< @brief pool for storing reduction info */ // provides random access and const pointers
	typedef std::map<size_t, TReduction*> _TyReductionMap; /**< @brief map of reductions, ordered by dest block address */ // note that the dest block address is duplicated here, maybe could use std::set with an appropriate comparison

	_TyPool m_p_reduction_pool[n_reductor_num]; /**< @brief list of reduction pools, one per each vertex dimension */
	CUberBlockMatrix m_p_matrix[n_reductor_num]; /**< @brief list of data stores for the per-edge contributions, one per each vertex dimension */ // storage only (could we do the same with just a pool and an allocator from CUberBlockMatrix, that would be better, p_Get_DenseStorage() could be protected again)
	_TyReductionMap m_p_reduction[n_reductor_num]; /**< @brief list of reductions, one per each vertex dimension */ // need to have them sorted by size for loop unrolling

	/**
	 *	@brief reduction function object
	 *	@note This is used to calculate reduction of all the vertices' r.h.s. vectors at once.
	 */
	class CReduce {
	protected:
		const CVectorReductionPlan<CDimsList> *m_p_this; /**< @brief pointer to the parent reduction plan */
		Eigen::VectorXd &m_r_dest; /**< @brief destination r.h.s. vector */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in] p_this is pointer to the parent reduction plan
		 *	@param[out] r_dest is destination r.h.s. vector (filled once the function operator is invoked)
		 */
		inline CReduce(const CVectorReductionPlan<CDimsList> *p_this, Eigen::VectorXd &r_dest)
			:m_p_this(p_this), m_r_dest(r_dest)
		{}

		/**
		 *	@brief function operator; performs all the reductions of one vertex dimension
		 *	@tparam C1DSize is vertex size (a specialization of fbs_ut::CCTSize)
		 */
		template <class C1DSize>
		inline void operator ()()
		{
			typedef CFindTypelistItem<_TyDimensionList, C1DSize> CSearch; // look for the size
			typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
			// the search happens at compile-time, each call to t_GetTempBlock is already linked to use a specific index

			typedef typename Eigen::VectorBlock<Eigen::VectorXd, C1DSize::n_size> CDestMap;
			typedef typename CUberBlockMatrix::CMakeMatrixRef<C1DSize::n_size, 1>::_TyConst CSrcMap;

			const _TyPool &reduction_pool = m_p_this->m_p_reduction_pool[CSearch::n_index];
			size_t _n = reduction_pool.size();
			_ASSERTE(_n <= INT_MAX);
			int n = int(_n);
			#pragma omp parallel for if(n > 50) // todo - dynamic schedule and below as well
			for(int i = 0; i < n; ++ i) {
				const TReduction &r_red = reduction_pool[i];
				CDestMap dest_map = m_r_dest.segment<C1DSize::n_size>(r_red.n_offset);
				CSrcMap src0_map(r_red.src_list.front());
				dest_map = src0_map; // can this be used? can the first block still use the original block inside the matrix, without having to use a temporary? probably not. todo
				for(size_t j = 1, m = r_red.src_list.size(); j < m; ++ j)
					dest_map += CSrcMap(r_red.src_list[j]);
				// reduce
			}
		}
	};

	/**
	 *	@brief reduction function object
	 *	@note This is used to calculate reduction of a subset of vertices.
	 */
	class CReduceRange {
	protected:
		const CVectorReductionPlan<CDimsList> *m_p_this; /**< @brief pointer to the parent reduction plan */
		Eigen::VectorXd &m_r_dest; /**< @brief destination r.h.s. vector */
		size_t m_n_begin; /**< @brief zero-based index of the first element of the r.h.s. vector to calculate (not vertex id) */
		size_t m_n_end; /**< @brief zero-based index of one past the last element of the r.h.s. vector to calculate (not vertex id) */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in] p_this is pointer to the parent reduction plan
		 *	@param[out] r_dest is destination r.h.s. vector (filled once the function operator is invoked)
		 *	@param[in] n_begin is zero-based index of the first element of the r.h.s. vector to calculate (not vertex id)
		 *	@param[in] n_end is zero-based index of one past the last element of the r.h.s. vector to calculate (not vertex id)
		 */
		inline CReduceRange(const CVectorReductionPlan<CDimsList> *p_this,
			Eigen::VectorXd &r_dest, size_t n_begin, size_t n_end)
			:m_p_this(p_this), m_r_dest(r_dest), m_n_begin(n_begin), m_n_end(n_end)
		{}

		/**
		 *	@brief function operator; performs the selected reductions of one vertex dimension
		 *	@tparam C1DSize is vertex size (a specialization of fbs_ut::CCTSize)
		 */
		template <class C1DSize>
		inline void operator ()()
		{
			typedef CFindTypelistItem<_TyDimensionList, C1DSize> CSearch; // look for the size
			typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
			// the search happens at compile-time, each call to t_GetTempBlock is already linked to use a specific index

			typedef typename Eigen::VectorBlock<Eigen::VectorXd, C1DSize::n_size> CDestMap;
			typedef typename CUberBlockMatrix::CMakeMatrixRef<C1DSize::n_size, 1>::_TyConst CSrcMap;

			const _TyReductionMap &reduction_map = m_p_this->m_p_reduction[CSearch::n_index];
			typename _TyReductionMap::const_iterator p_begin_it = reduction_map.lower_bound(m_n_begin);
			typename _TyReductionMap::const_iterator p_end_it = reduction_map.upper_bound(m_n_end);
			// find range of reductions of the selected size to perform

			for(; p_begin_it != p_end_it; ++ p_begin_it) {
				const TReduction &r_red = *(*p_begin_it).second;
				CDestMap dest_map = m_r_dest.segment<C1DSize::n_size>(r_red.n_offset);
				CSrcMap src0_map(r_red.src_list.front());
				dest_map = src0_map; // can this be used? can the first block still use the original block inside the matrix, without having to use a temporary? probably not. // not sure what was wrong with it, this works.
				for(size_t j = 1, m = r_red.src_list.size(); j < m; ++ j)
					dest_map += CSrcMap(r_red.src_list[j]);
				// reduce
			}
			// can't do this easily in parallel
		}
	};

public:
	/**
	 *	@brief destructor; performs consistency checks in debug
	 */
	~CVectorReductionPlan()
	{
#ifdef _DEBUG
		int p_dims_list[n_reductor_num];
		fbs_ut::Copy_CTSizes_to_Array<_TyDimensionList>(p_dims_list);
		// convert the typelist to an array so that we can index it at runtime

		std::vector<std::pair<size_t, int> > allocated_segments;
		for(int i = 0; i < n_reductor_num; ++ i) {
			const int n_dim = p_dims_list[i];
			//_ASSERTE(!i || n_dim > p_dims_list[i - 1]); // it is not sorted, only unique, would have to check differently
			const _TyReductionMap &reduction_map = m_p_reduction[i];
			for(typename _TyReductionMap::const_iterator p_begin_it = reduction_map.begin(),
			   p_end_it = reduction_map.end(); p_begin_it != p_end_it; ++ p_begin_it) {
				const TReduction &r_red = *(*p_begin_it).second;
				std::pair<size_t, int> segment(r_red.n_offset, n_dim);
				allocated_segments.push_back(segment);
			}
		}
		// collect allocated segments as (offset, dimension) pairs

		std::sort(allocated_segments.begin(), allocated_segments.end());
		// sort the segments (they are sorted for each dimension, merging them explicitly
		// would be more efficient but also error prone; this is debug, correctness is paramount)

		_ASSERTE(std::unique(allocated_segments.begin(), allocated_segments.end()) == allocated_segments.end());
		// make sure there are no duplicates

		_ASSERTE(allocated_segments.empty() || !allocated_segments.front().first);
		// make sure that the first segment starts at zero

		for(size_t i = 1, n = allocated_segments.size(); i < n; ++ i)
			_ASSERTE(allocated_segments[i].first == allocated_segments[i - 1].first + allocated_segments[i - 1].second);
		// make sure that the next segment starts where the previous one ends
#endif // _DEBUG
	}

	/**
	 *	@brief gets the maximum reduced dimension
	 *	@return Returns the maximum reduced dimension, in elements.
	 *	@note This is also the expected size of the r.h.s. vector.
	 */
	size_t n_Max_Dimension() const
	{
		int p_dims_list[n_reductor_num];
		fbs_ut::Copy_CTSizes_to_Array<_TyDimensionList>(p_dims_list);
		// convert the typelist to an array so that we can index it at runtime

		size_t n_max = 0;
		for(int i = 0; i < n_reductor_num; ++ i) {
			if(!m_p_reduction[i].empty()) {
				typename _TyReductionMap::const_iterator p_back_it = -- m_p_reduction[i].end(); // no .back() on this thing
				const TReduction &r_red = *(*p_back_it).second;
				size_t n_last = r_red.n_offset + p_dims_list[i];
				if(n_max < n_last)
					n_max = n_last;
			}
		}
		// find the one past the last element that will be written

		return n_max;
	}

	/**
	 *	@brief calculates the size of this object in memory
	 *	@return Returns the size of this object (and of all associated
	 *		arrays or buffers) in memory, in bytes.
	 */
	size_t n_Allocation_Size() const
	{
		size_t n_size = sizeof(CVectorReductionPlan<CDimsList>);
		size_t n_data_size = 0, n_maps_size = 0, n_pools_size = 0,
			n_vectors_size = 0, n_vectors_slack = 0;
		for(int i = 0; i < n_reductor_num; ++ i) {
			n_data_size += m_p_matrix[i].n_Allocation_Size() - sizeof(CUberBlockMatrix);
			n_maps_size += n_Map_Allocation_Size(m_p_reduction[i]) - sizeof(m_p_reduction[i]);
			n_pools_size += m_p_reduction_pool[i].capacity() * sizeof(TReduction) +
				m_p_reduction_pool[i].page_num() * sizeof(TReduction*);
			for(size_t j = 0, m = m_p_reduction_pool[i].size(); j < m; ++ j) {
				const std::vector<const double*> src_list = m_p_reduction_pool[i][j].src_list;
				n_vectors_size += src_list.capacity() * sizeof(const double*);
				n_vectors_slack += (src_list.capacity() - src_list.size()) * sizeof(const double*);
			}
		}
		return n_size + n_data_size + n_maps_size + n_pools_size + n_vectors_size;
	}

	/**
	 *	@brief gets a temporary vector assigned for reduction
	 *
	 *	@tparam n_dimension is size of the requested vector (must match one of vertex dimensions)
	 *
	 *	@param[in] n_vector_offset is offset in the r.h.s. vector, in elements
	 *
	 *	@return Returns pointer to the assigned memory (guaranteed not to change,
	 *		deleted at the end of the lifetime of this object).
	 *
	 *	@note Note that the reduction conflicts are unchecked here (could have multiple block
	 *		sizes reducing in the same destination or could have overlapping blocks).
	 *	@note This function throws std::bad_alloc.
	 */
	template <const int n_dimension>
	double *p_Get_ReductionBlock(size_t n_vector_offset) // throw(std::bad_alloc)
	{
		typedef fbs_ut::CCTSize<n_dimension> C1DBlockSize;
		typedef CFindTypelistItem<_TyDimensionList, C1DBlockSize> CSearch; // look for the size
		typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
		// the search happens at compile-time, each call to p_GetTempBlock is already linked to use a specific index

		_TyReductionMap &r_reduction_list = m_p_reduction[CSearch::n_index];
		_TyPool &r_pool = m_p_reduction_pool[CSearch::n_index];
		typename _TyReductionMap::iterator p_red_it = r_reduction_list.find(n_vector_offset);
		if(p_red_it == r_reduction_list.end()) {
			r_pool.resize(r_pool.size() + 1);
			TReduction *p_red = &(*(r_pool.end() - 1));
			p_red->n_offset = n_vector_offset; // remember
			_ASSERTE(p_red->src_list.empty()); // should be initialized and empty
			p_red_it = r_reduction_list.insert(std::pair<const size_t,
				TReduction*>(n_vector_offset, p_red)).first;
		}
		TReduction &r_red = *(*p_red_it).second;
		r_red.src_list.reserve(r_red.src_list.size() + 1);
		CUberBlockMatrix &r_matrix = m_p_matrix[CSearch::n_index];
		double *p_block = r_matrix.p_Get_DenseStorage(n_dimension);
		r_red.src_list.push_back(p_block); // should not throw now

		return p_block; // store result for reduction here
	}

	/**
	 *	@brief reduce all the blocks (runs in parallel, where available)
	 *	@param[out] r_dest is the destination vector (must be allocated by the caller)
	 *	@note This does not check the integrity of the reduction; if initialized
	 *		incorrectly, some parts of the vector can be left uninitialized.
	 */
	void ReduceAll(Eigen::VectorXd &r_dest) const
	{
		_ASSERTE(r_dest.rows() == n_Max_Dimension()); // check the size of the vector
		CTypelistForEach<_TyDimensionList, CReduce>::Run(CReduce(this, r_dest));
	}

	/**
	 *	@brief reduce all the blocks (runs in parallel, where available)
	 *
	 *	@param[out] r_dest is the destination vector (must be allocated by the caller)
	 *	@param[in] n_begin is zero-based index of the first element of the r.h.s. vector to calculate (not vertex id)
	 *	@param[in] n_end is zero-based index of one past the last element of the r.h.s. vector to calculate (not vertex id)
	 *
	 *	@note This does not check if the begin / end boundaries match vertex boundaries. Upper bound function
	 *		is used to find the nearest conservative boundaries (slightly more elements are updated if the
	 *		begin / end is misaligned).
	 *	@note This does not check the integrity of the reduction; if initialized
	 *		incorrectly, some parts of the vector can be left uninitialized.
	 */
	void ReduceRange(Eigen::VectorXd &r_dest, size_t n_begin, size_t n_end) const // ReduceSingle() would be easier to implement and could run in parallel
	{
		_ASSERTE(r_dest.rows() == n_Max_Dimension()); // check the size of the vector
		CTypelistForEach<_TyDimensionList, CReduceRange>::Run(CReduceRange(this, r_dest, n_begin, n_end));
	}

	/**
	 *	@brief reduce a single block
	 *
	 *	@tparam n_dimension is size of the requested vector (must match one of vertex dimensions)
	 *
	 *	@param[out] r_dest is the destination vector (must be allocated by the caller)
	 *	@param[in] n_order is order of the vertex to reduce (offset in the r.h.s. vector, in elements)
	 *
	 *	@note This assumes that a vertex with the given order is in the reduction plan and that the
	 *		dimension is correct (only checked in debug).
	 */
	template <const int n_dimension>
	void Reduce_Single(Eigen::VectorXd &r_dest, size_t n_order) const
	{
		_ASSERTE(r_dest.rows() == n_Max_Dimension()); // check the size of the vector

		typedef fbs_ut::CCTSize<n_dimension> C1DBlockSize;
		typedef CFindTypelistItem<_TyDimensionList, C1DBlockSize> CSearch; // look for the size
		typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
		// the search happens at compile-time, each call to p_GetTempBlock is already linked to use a specific index

		const _TyReductionMap &r_reduction_list = m_p_reduction[CSearch::n_index];
		typename _TyReductionMap::const_iterator p_it = r_reduction_list.find(n_order);
		_ASSERTE(p_it != r_reduction_list.end()); // it should not

		typedef typename Eigen::VectorBlock<Eigen::VectorXd, n_dimension> CDestMap;
		typedef typename CUberBlockMatrix::CMakeMatrixRef<n_dimension, 1>::_TyConst CSrcMap;

		const TReduction &r_red = *(*p_it).second;
		CDestMap dest_map = r_dest.segment<n_dimension>(r_red.n_offset);
		CSrcMap src0_map(r_red.src_list.front());
		dest_map = src0_map; // can this be used? can the first block still use the original block inside the matrix, without having to use a temporary? probably not. todo
		for(size_t j = 1, m = r_red.src_list.size(); j < m; ++ j)
			dest_map += CSrcMap(r_red.src_list[j]);
		// reduce
	}
};

/**
 *	@brief matrix reduction key type selection and traits
 *
 *	For reduced block identification, one can use either pointer to the original
 *	block in lambda or its coordinates. Using coordinates has minor advantages for
 *	the edge system, as most of the edges do not have the pointer to the original
 *	block, instead they store a pointer to a reduced block and getting the original
 *	then requiers <tt>O(log n)</tt> lookup.
 *
 *	@tparam b_use_block_coord_as_reduction_key is key type selector
 */
template <bool b_use_block_coord_as_reduction_key>
class CMatrixReductionKey_Traits {
public:
	typedef const double *TKey; /**< @brief key type */

	/**
	 *	@brief utility function; distills key from block description
	 *
	 *	@param[in] n_row is zero-based block row (unused)
	 *	@param[in] n_col is zero-based block column (unused)
	 *	@param[in] p_address is block address
	 *
	 *	@return Returns key of the given block.
	 */
	static TKey t_MakeKey(size_t UNUSED(n_row), size_t UNUSED(n_col), const double *p_address)
	{
		return p_address;
	}
};

/**
 *	@brief matrix reduction key type selection and traits (specialization for coordinate-based keys)
 */
template <>
class CMatrixReductionKey_Traits<true> {
public:
	typedef std::pair<size_t, size_t> TKey; /**< @brief key type */

	/**
	 *	@brief utility function; distills key from block description
	 *
	 *	@param[in] n_row is zero-based block row
	 *	@param[in] n_col is zero-based block column
	 *	@param[in] p_address is block address (unused)
	 *
	 *	@return Returns key of the given block.
	 */
	static TKey t_MakeKey(size_t n_row, size_t n_col, const double *UNUSED(p_address))
	{
		return std::make_pair(n_row, n_col);
	}
};

/**
 *	@brief parallel reduction plan for efficiently calculating and updating the hessian matrix
 *
 *	@tparam CDimsList is list of hessian matrix block dimensions, as fbs_ut::CCTSize2D
 *
 *	@note This uses block dimensions to differentiate between blocks, assumes that there will
 *		be no conflicts between blocks of different dimensions (does not check for that).
 *	@todo Redesign the pointers to be objects that wrap the pointer and remove the illusion of
 *		being constant (the shared pointers may change upon block conflict without the knowledge
 *		of the object owning it).
 */
template <class CDimsList>
class CMatrixReductionPlan { // todo - need clear!
public:
	typedef typename CTransformTypelist<CDimsList, fbs_ut::CEigenToDimension>::_TyResult _TyDimensionList; /**< @brief list of hessian matrix block dimensions, as fbs_ut::CCTSize2D */

	/**
	 *	@brief parameters, stored as enum
	 */
	enum {
		n_reductor_num = CTypelistLength<_TyDimensionList>::n_result, /**< @brief number of different reduction sizes */
		n_pool_page_size = 4096, /**< @brief reduction pool page size */
		n_pool_memory_align = 0, /**< @brief reduction pool memory alignment */ // see no benefit in alignment right now, this stores the TReduction elements, those do not need to be aligned
		b_use_block_coord_keys = 1 /**< @brief if set, use <tt>(row, col)</tt> coordinates instead of pointers to identify the blocks */
	};

	typedef typename CMatrixReductionKey_Traits<b_use_block_coord_keys != 0>::TKey TKey; /**< @brief block key type */

	/**
	 *	@brief reduction description
	 *	@note This does not store block dimensions; that is given by the index of the pool that this is found in.
	 */
	struct TReduction {
		double *p_dest; /**< @brief destination block in the lambda matrix */
		std::vector<const double*> src_list; /**< @brief list of reduced blocks */
	};

	typedef forward_allocated_pool<TReduction,
		n_pool_page_size, n_pool_memory_align> _TyPool; /**< @brief pool for storing reduction info */ // provides random access and const pointers
	typedef std::map<TKey, TReduction*> _TyReductionMap; /**< @brief map of reductions, ordered by dest block address / coords */ // sorted by dest block address / coords
	typedef std::map<const double*, double**> _TyOwnerLookup; /**< @brief reverse block lookup */ // sorted by dest block address

protected:
	_TyOwnerLookup m_p_owner_lookup[n_reductor_num]; /**< @brief list of reverse block lookups, one per each vertex dimension @note This only contains records for the un-conflicted owners; it eliminates copying the Jacobian contribution from per-edge memory to the matrix (size 1 reduction). @note Only the off-diagonal block owners are stored here. The diagonal ones where conflicts are anticipated do not have size 1 reduction elimination. @note In SLAM, there can be some conflicts in the off-diagonal blocks (multiple edges between the same vertices, e.g. from different sensors, re-observation, etc.). In BA, there typically aren't. */
	_TyPool m_p_reduction_pool[n_reductor_num]; /**< @brief list of reduction pools, one per each vertex dimension */ // actually need that for parallel processing
	_TyReductionMap m_p_reduction_list[n_reductor_num]; /**< @brief list of reductions, one per each vertex dimension */
	CUberBlockMatrix m_p_matrix[n_reductor_num]; /**< @brief list of data stores for the per-edge contributions, one per each vertex dimension */ // storage only (could we do the same with just a pool and an allocator from CUberBlockMatrix, that would be better, p_Get_DenseStorage() could be protected again)

	/**
	 *	@brief reduction function object
	 */
	class CReduce {
	protected:
		const CMatrixReductionPlan<CDimsList> *m_p_this; /**< @brief pointer to the parent reduction plan */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] p_this is pointer to the parent reduction plan
		 */
		inline CReduce(const CMatrixReductionPlan<CDimsList> *p_this)
			:m_p_this(p_this)
		{}

		/**
		 *	@brief function operator; performs all the reductions of a single Jacobian dimension
		 *	@tparam C2DSize is Jacobian size (a specialization of fbs_ut::CCTSize2D)
		 */
		template <class C2DSize>
		inline void operator ()()
		{
			typedef CFindTypelistItem<_TyDimensionList, C2DSize> CSearch; // look for the size
			typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
			// the search happens at compile-time, each call to t_GetTempBlock is already linked to use a specific index

			typedef typename CUberBlockMatrix::CMakeMatrixRef<C2DSize::n_row_num,
				C2DSize::n_column_num>::_Ty CDestMap;
			typedef typename CUberBlockMatrix::CMakeMatrixRef<C2DSize::n_row_num,
				C2DSize::n_column_num>::_TyConst CSrcMap;

			const _TyPool &reduction_pool = m_p_this->m_p_reduction_pool[CSearch::n_index];
			size_t _n = reduction_pool.size();
			_ASSERTE(_n <= INT_MAX);
			int n = int(_n);
			#pragma omp parallel for if(n > 50) // todo - export the parallel thresh as an arg, maybe consider dynamic schedule
			for(int i = 0; i < n; ++ i) {
				const TReduction &r_red = reduction_pool[i];
				CDestMap dest_map((double*)r_red.p_dest);
				CSrcMap src0_map(r_red.src_list.front());
				dest_map = src0_map; // can this be used? can the first block still use the original block inside the matrix, without having to use a temporary? probably not. todo
				for(size_t j = 1, m = r_red.src_list.size(); j < m; ++ j)
					dest_map += CSrcMap(r_red.src_list[j]);
				// reduce
			}
		}
	};

public:
	/**
	 *	@brief calculates the size of this object in memory
	 *	@return Returns the size of this object (and of all associated
	 *		arrays or buffers) in memory, in bytes.
	 */
	size_t n_Allocation_Size() const
	{
		size_t n_size = sizeof(CMatrixReductionPlan<CDimsList>);
		size_t n_data_size = 0, n_maps_size = 0, n_pools_size = 0,
			n_vectors_size = 0, n_vectors_slack = 0;
		for(int i = 0; i < n_reductor_num; ++ i) {
			n_data_size += m_p_matrix[i].n_Allocation_Size() - sizeof(CUberBlockMatrix);
			n_maps_size += n_Map_Allocation_Size(m_p_owner_lookup[i]) - sizeof(m_p_owner_lookup[i]);
			n_maps_size += n_Map_Allocation_Size(m_p_reduction_list[i]) - sizeof(m_p_reduction_list[i]);
			n_pools_size += m_p_reduction_pool[i].capacity() * sizeof(TReduction) +
				m_p_reduction_pool[i].page_num() * sizeof(TReduction*);
			for(size_t j = 0, m = m_p_reduction_pool[i].size(); j < m; ++ j) {
				const std::vector<const double*> src_list = m_p_reduction_pool[i][j].src_list;
				n_vectors_size += src_list.capacity() * sizeof(const double*);;
				n_vectors_slack += (src_list.capacity() - src_list.size()) * sizeof(const double*);
			}
		}
		return n_size + n_data_size + n_maps_size + n_pools_size + n_vectors_size;
	}

	/**
	 *	@brief gets a single reduction block in case we know it is most likely going to be conflicted (like the blocks at the diagonal)
	 *
	 *	@tparam C2DSize is block size as fbs_ut::CCTSize2D
	 *
	 *	@param[in] n_row is zero-based block row
	 *	@param[in] n_col is zero-based block column
	 *	@param[in] p_reduction_dest is pointer to the first element of the block inside the matrix
	 *
	 *	@return Returns pointer to a temporary block from which the values will be reduced to the specified destination.
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	template <class C2DSize>
	inline double *p_Diagonal_GetTempBlock(size_t n_row, size_t n_col, double *p_reduction_dest) // throw(std::bad_alloc)
	{
		//_ASSERTE(n_row == n_col); // it does not have to be diagonal, it is just a type of block where collision is anticipated in most of the blocks

		TKey t_key = t_MakeKey(n_row, n_col, p_reduction_dest);
		typedef CFindTypelistItem<_TyDimensionList, C2DSize> CSearch; // look for the size
		typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
		// the search happens at compile-time, each call to p_GetTempBlock is already linked to use a specific index

		_TyReductionMap &r_reduction_list = m_p_reduction_list[CSearch::n_index];
		CUberBlockMatrix &r_storage = m_p_matrix[CSearch::n_index];
		_TyPool &r_pool = m_p_reduction_pool[CSearch::n_index];
		return p_GetSingle(r_reduction_list, r_storage, r_pool, C2DSize::n_row_num *
			C2DSize::n_column_num, t_key, p_reduction_dest); // just push those three on the stack and go
	}

	/**
	 *	@brief gets a single reduction block in case we know it is most likely going to be unique (like the off-diagonal blocks)
	 *
	 *	@tparam C2DSize is block size as fbs_ut::CCTSize2D
	 *
	 *	@param[in] p_reduction_dest is pointer to the first element of the block inside the matrix
	 *	@param[in] p_owner_storage is pointer to the pointer to address of the block inside the
	 *		owner object (is liable to be changed upon conflict)
	 *
	 *	@return Returns pointer to the original block (no conflict occured yet).
	 *
	 *	@note This should be used in case the block did not exist in the matrix (<tt>b_uninitialized</tt> is set).
	 *	@note This assumes that there is only a single pointer to the block stored, which can be replaced
	 *		by the pointer to a pointer, passed as the second argument. In case this does not apply, this
	 *		reductor cannot be used.
	 *	@note This function throws std::bad_alloc.
	 */
	template <class C2DSize>
	inline double *p_OffDiagonal_GetTempBlock(double *p_reduction_dest, double **p_owner_storage) // throw(std::bad_alloc)
	{
#ifdef _DEBUG
		typedef CFindTypelistItem<_TyDimensionList, C2DSize> CSearch; // look for the size
		typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
		// the search happens at compile-time, each call to p_GetTempBlock is already linked to use a specific index

		_TyReductionMap &r_reduction_list = m_p_reduction_list[CSearch::n_index];
		_ASSERTE(b_use_block_coord_keys || r_reduction_list.find(t_MakeKey(0,
			0, p_reduction_dest)) == r_reduction_list.end()); // if b_use_block_coord_keys is set, we can't verify anything
		// make sure the block is not there
#endif // _DEBUG

		Set_BlockOwner<C2DSize>(p_reduction_dest, p_owner_storage);
		return p_reduction_dest;
	}

	/**
	 *	@brief gets a single reduction block in case we know it is most likely going to be unique (like the off-diagonal blocks)
	 *
	 *	@tparam C2DSize is block size as fbs_ut::CCTSize2D
	 *
	 *	@param[in] n_row is zero-based block row
	 *	@param[in] n_col is zero-based block column
	 *	@param[in] p_reduction_dest is pointer to the first element of the block inside the matrix
	 *
	 *	@return Returns pointer to a temporary block from which the values will be reduced to the specified destination.
	 *
	 *	@note This should be used in case the block did already exist in the matrix (<tt>b_uninitialized</tt> is not set).
	 *	@note This function throws std::bad_alloc.
	 */
	template <class C2DSize>
	inline double *p_OffDiagonal_GetTempBlock(size_t n_row, size_t n_col, double *p_reduction_dest)
	{
		TKey t_key = t_MakeKey(n_row, n_col, p_reduction_dest);
		std::pair<double*, double*> t_storage = t_GetTempBlock<C2DSize>(t_key, p_reduction_dest);
		if(t_storage.second) {
			double **p_owner_variable = p_Get_BlockOwner<C2DSize>(p_reduction_dest, true);
			_ASSERTE(p_owner_variable != 0);
			// we are replacing the owner, there should be one registered

			*p_owner_variable = t_storage.first;
			// the original owner should use the first temp block

			memcpy(t_storage.first, p_reduction_dest, C2DSize::n_row_num *
				C2DSize::n_column_num * sizeof(double));
			// this block might already have significant value, need to store it in the temp block

			return t_storage.second;
			// the second reducer should use the second temp block
		} else {
			_ASSERTE(!p_Get_BlockOwner<C2DSize>(p_reduction_dest, false)); // there is no owner anymore, we already replaced it
			return t_storage.first;
			// nothing special, just another block
		}
	}

	/**
	 *	@brief reduce all the blocks (runs in parallel, where available)
	 */
	void ReduceAll() const
	{
		/*eigen::MatrixXd::identity();
		for(size_t i = 0; i < n_reductor_num; ++ i) { // todo - do typelist_foreach first, to have fixed block size ;) (will probably lose for loop scheduling though)
			_TyReductionMap &r_reduction_list = m_p_reduction_list[i];
			_TyReductionMap::iterator p_it = r_reduction_list.begin();
			_TyReductionMap::iterator p_end_it = r_reduction_list.end();
			for(; p_it != p_end_it; ++ p_it) { // duh; how to parallelize that? can't.
			}
		}*/
		// can do everything in parallel, need to see which strategy is the fastest

		CTypelistForEach<_TyDimensionList, CReduce>::Run(CReduce(this));
	}

	/**
	 *	@brief reduce a single block
	 *
	 *	Reduce a single block address, this is callable by an edge, probably does
	 *	good enough job in incremental updates where a single edge / a couple edges
	 *	is added. It will recalculate a couple of vertex' sums, but it should still
	 *	be less than all of them, as it is now.
	 *
	 *	@tparam C2DSize is block size as fbs_ut::CCTSize2D
	 *	@param[in] t_key is block key (address / row col coordinate)
	 *	@note Note that this can not run in parallel over the edges. to do it in parallel,
	 *		one would need to collect the *unique* block addresses and then run in parallel.
	 */
	template <class C2DSize>
	void ReduceSingle(TKey t_key) const
	{
		typedef CFindTypelistItem<_TyDimensionList, C2DSize> CSearch; // look for the size
		typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
		// the search happens at compile-time, each call to p_GetTempBlock is already linked to use a specific index

		const _TyReductionMap &r_reduction_list = m_p_reduction_list[CSearch::n_index];
		// get the list

		typename _TyReductionMap::const_iterator p_red_it = r_reduction_list.find(t_key);
		if(p_red_it == r_reduction_list.end())
			return; // the block is not there (may happen with off-diagonal blocks, which do not have a conflict yet and do not need to be reduced)

		const TReduction &r_red = *(*p_red_it).second;
		_ASSERTE(!r_red.src_list.empty()); // should never be (but there can be 1, in case of a lonely vertex where collision was expected but did not occur (yet))
		// get the particular reduction

		typename CUberBlockMatrix::CMakeMatrixRef<C2DSize::n_row_num,
			C2DSize::n_column_num>::_Ty dest_map((double*)r_red.p_dest);
		typename CUberBlockMatrix::CMakeMatrixRef<C2DSize::n_row_num,
			C2DSize::n_column_num>::_TyConst src0_map(r_red.src_list.front());
		dest_map = src0_map;
		for(size_t i = 1, n = r_red.src_list.size(); i < n; ++ i) {
			typename CUberBlockMatrix::CMakeMatrixRef<C2DSize::n_row_num,
				C2DSize::n_column_num>::_TyConst src_map(r_red.src_list[i]);
			dest_map += src_map;
		}
		// reduce
	}

	/**
	 *	@brief utility function for making key
	 *
	 *	@param[in] n_row is zero-based block row
	 *	@param[in] n_col is zero-based block column
	 *	@param[in] p_address is block address
	 *
	 *	@return Returns key of the given block.
	 *
	 *	@note Depending on which type of key is used, some of the arguments are ignored.
	 */
	static inline TKey t_MakeKey(size_t n_row, size_t n_col, const double *p_address)
	{
		return CMatrixReductionKey_Traits<b_use_block_coord_keys != 0>::t_MakeKey(n_row, n_col, p_address);
	}

protected:
	/**
	 *	@brief gets one or two blocks, based on whether there would be a conflict
	 *
	 *	Get a temp block, and in case p_reduction_dest was not present in the list yet,
	 *	alloc also the second block for the original owner of p_reduction_dest (which
	 *	wrote directly to the matrix, bypassing the reduction).
	 *
	 *	@tparam C2DSize is block size as fbs_ut::CCTSize2D
	 *
	 *	@param[in] t_key is block key (address / row col coordinate)
	 *	@param[in] p_reduction_dest is pointer to the first element of the block inside the matrix
	 *
	 *	@return Returns a pair of reduced block storages, first is always a valid pointer and
	 *		second may be set to null if the block was already shared before.
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	template <class C2DSize>
	inline std::pair<double*, double*> t_GetTempBlock(TKey t_key, double *p_reduction_dest) // throw(std::bad_alloc)
	{
		typedef CFindTypelistItem<_TyDimensionList, C2DSize> CSearch; // look for the size
		typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
		// the search happens at compile-time, each call to t_GetTempBlock is already linked to use a specific index

		_TyReductionMap &r_reduction_list = m_p_reduction_list[CSearch::n_index];
		CUberBlockMatrix &r_storage = m_p_matrix[CSearch::n_index];
		_TyPool &r_pool = m_p_reduction_pool[CSearch::n_index];
		return t_GetPair(r_reduction_list, r_storage, r_pool, C2DSize::n_row_num *
			C2DSize::n_column_num, t_key, p_reduction_dest); // just push those three on the stack and go
		//_ASSERTE(t_pair.second && (m_p_owner_lookup[CSearch::n_index].find(p_reduction_dest) !=
		//	m_p_owner_lookup[CSearch::n_index].end())); // if there is second, there was a reduction conflict and the block should be owned (but not vice versa)
		//return t_pair;
	}

	/**
	 *	@brief sets pointer to the block owner pointer storage
	 *
	 *	This is used to bypass size 1 reductions in the blocks where the conflicts are unlikely.
	 *	This adds a record of the original pointer in the matrix, and a pointer to the pointer
	 *	to this block in the edge class, so that when a conflict occurs, the pointer in the original
	 *	owner edge can be modified. This mandates that pointers to the blocks be stored only in
	 *	a single instance (as multiple pointer instances can't be modified like this).
	 *
	 *	@tparam C2DSize is block size as fbs_ut::CCTSize2D
	 *
	 *	@param[in] p_block_addr is pointer to the first element of the block inside the matrix
	 *	@param[in] p_owner is pointer to the pointer to address of the block inside the
	 *		owner object (is liable to be changed upon conflict)
	 *
	 *	@note This should be used in case the block did not exist in the matrix (<tt>b_uninitialized</tt> is set).
	 *	@note This assumes that there is only a single pointer to the block stored, which can be replaced
	 *		by the pointer to a pointer, passed as the second argument. In case this does not apply, this
	 *		reductor cannot be used.
	 *	@note This function throws std::bad_alloc.
	 */
	template <class C2DSize>
	void Set_BlockOwner(const double *p_block_addr, double **p_owner) // throw(std::bad_alloc)
	{
		typedef CFindTypelistItem<_TyDimensionList, C2DSize> CSearch; // look for the size
		typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
		// the search happens at compile-time, each call to t_GetTempBlock is already linked to use a specific index

		_TyOwnerLookup &r_lookup_map = m_p_owner_lookup[CSearch::n_index];
		_ASSERTE(r_lookup_map.find(p_block_addr) == r_lookup_map.end()); // make sure that there is no owner so far
		r_lookup_map[p_block_addr] = p_owner;
		// stores the original block woner
	}

	/**
	 *	@brief gets a pointer to block owner pointer
	 *
	 *	@tparam C2DSize is block size as fbs_ut::CCTSize2D
	 *
	 *	@param[in] p_block_addr is pointer to the first element of the block inside the matrix
	 *	@param[in] b_clear_owner is clear owner flag (if set, and the order exists, it is )
	 *
	 *	@return Returns pointer to the pointer inside the owner of the specified block, or null if the
	 *		block does not have an owner (or has multiple owners and the first owner was cleared).
	 */
	template <class C2DSize>
	double **p_Get_BlockOwner(const double *p_block_addr, bool b_clear_owner = true)
	{
		typedef CFindTypelistItem<_TyDimensionList, C2DSize> CSearch; // look for the size
		typedef typename CReductionPlanAssert<CSearch::b_result>::BLOCK_SIZE_NOT_PRESENT_IN_THE_LIST _TyAssert0; // make sure it is there
		// the search happens at compile-time, each call to t_GetTempBlock is already linked to use a specific index

		_TyOwnerLookup &r_lookup_map = m_p_owner_lookup[CSearch::n_index];
		_TyOwnerLookup::iterator p_it = r_lookup_map.find(p_block_addr);
		if(p_it == r_lookup_map.end())
			return 0; // no such owner
		double **p_owner = (*p_it).second;
		if(b_clear_owner)
			r_lookup_map.erase(p_it); // remove, now all edges will use the reductor anyway
		return p_owner;
		// find the owner
	}

	/**
	 *	@brief gets a single reduced block
	 *
	 *	@param[in,out] r_reduction_list is reduction map
	 *	@param[in,out] r_storage is matrix where the reduced blocks are allocated
	 *	@param[in,out] r_pool is a pool for reduction description objects
	 *	@param[in] n_size is size of the requested block, in elements
	 *	@param[in] t_key is block key (address / row col coordinate)
	 *	@param[in] p_reduction_dest is pointer to the first element of the block inside the matrix
	 *
	 *	@return Returns a pointer to reduced block storage.
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	double *p_GetSingle(_TyReductionMap &r_reduction_list, CUberBlockMatrix &r_storage,
		_TyPool &r_pool, int n_size, TKey t_key, double *p_reduction_dest) // throw(std::bad_alloc)
	{
		typename _TyReductionMap::iterator p_red_it = r_reduction_list.find(t_key);
		if(p_red_it == r_reduction_list.end()) {
			r_pool.resize(r_pool.size() + 1);
			TReduction *p_red = &(*(r_pool.end() - 1));
			p_red->p_dest = p_reduction_dest; // remember
			_ASSERTE(p_red->src_list.empty()); // should be initialized and empty
			p_red_it = r_reduction_list.insert(std::pair<const TKey,
				TReduction*>(t_key, p_red)).first;
		}
		std::vector<const double*> &r_list = (*p_red_it).second->src_list;
		double *p_ptr = r_storage.p_Get_DenseStorage(n_size);
		r_list.push_back(p_ptr);
		return p_ptr;
	}

	/**
	 *	@brief gets one or two blocks, based on whether there would be a conflict
	 *
	 *	@param[in,out] r_reduction_list is reduction map
	 *	@param[in,out] r_storage is matrix where the reduced blocks are allocated
	 *	@param[in,out] r_pool is a pool for reduction description objects
	 *	@param[in] n_size is size of the requested block, in elements
	 *	@param[in] t_key is block key (address / row col coordinate)
	 *	@param[in] p_reduction_dest is pointer to the first element of the block inside the matrix
	 *
	 *	@return Returns a pair of reduced block storages, first is always a valid pointer and
	 *		second may be set to null if the block was already shared before.
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	std::pair<double*, double*> t_GetPair(_TyReductionMap &r_reduction_list,
		CUberBlockMatrix &r_storage, _TyPool &r_pool, int n_size, TKey t_key, double *p_reduction_dest) // throw(std::bad_alloc)
	{
		typename _TyReductionMap::iterator p_red_it = r_reduction_list.find(t_key);
		if(p_red_it == r_reduction_list.end()) {
			r_pool.resize(r_pool.size() + 1);
			TReduction *p_red = &(*(r_pool.end() - 1));
			p_red->p_dest = p_reduction_dest; // remember
			_ASSERTE(p_red->src_list.empty()); // should be initialized and empty
			p_red_it = r_reduction_list.insert(std::pair<const TKey,
				TReduction*>(t_key, p_red)).first;
		}
		std::vector<const double*> &r_list = (*p_red_it).second->src_list;
		if(!r_list.empty()) {
			double *p_ptr = r_storage.p_Get_DenseStorage(n_size);
			r_list.push_back(p_ptr);
			return std::make_pair(p_ptr, (double*)0);
			// just another block
		} else {
			double *p_ptr0 = r_storage.p_Get_DenseStorage(n_size);
			double *p_ptr1 = r_storage.p_Get_DenseStorage(n_size); // must call twice, want each of them aligned
			r_list.push_back(p_ptr0);
			r_list.push_back(p_ptr1);
			return std::make_pair(p_ptr0, p_ptr1);
			// the first time around: alloc two blocks
		}
	}
};

/**
 *	@brief wrapper reduction plans for lambda and the right-hand-side vector
 *
 *	@tparam CDimsList is list of hessian matrix block dimensions, as fbs_ut::CCTSize2D
 */
template <class CDimsList>
class CLambdaReductionPlan {
public:
	typedef CVectorReductionPlan<CDimsList> CRHSReductor; /**< @brief right hand side vector reduction plan type */
	typedef CMatrixReductionPlan<CDimsList> CLambdaReductor; /**< @brief lambda reduction plan type */

protected:
	CVectorReductionPlan<CDimsList> m_vec_plan; /**< @brief right hand side vector reduction plan */
	CMatrixReductionPlan<CDimsList> m_mat_plan; /**< @brief lambda reduction plan */

public:
	/**
	 *	@brief gets right hand side vector reduction plan
	 *	@return Returns a reference to the right hand side vector reduction plan.
	 */
	inline CRHSReductor &r_RHS_ReductionPlan()
	{
		return m_vec_plan;
	}

	/**
	 *	@brief gets lambda reduction plan
	 *	@return Returns a reference to the lambda reduction plan.
	 */
	inline CLambdaReductor &r_Lambda_ReductionPlan()
	{
		return m_mat_plan;
	}

	/**
	 *	@brief gets right hand side vector reduction plan
	 *	@return Returns a const reference to the right hand side vector reduction plan.
	 */
	inline const CRHSReductor &r_RHS_ReductionPlan() const
	{
		return m_vec_plan;
	}

	/**
	 *	@brief gets lambda reduction plan
	 *	@return Returns a const reference to the lambda reduction plan.
	 */
	inline const CLambdaReductor &r_Lambda_ReductionPlan() const
	{
		return m_mat_plan;
	}

	/**
	 *	@brief calculates the size of this object in memory
	 *	@return Returns the size of this object (and of all associated
	 *		arrays or buffers) in memory, in bytes.
	 */
	size_t n_Allocation_Size() const
	{
		return m_vec_plan.n_Allocation_Size() + m_mat_plan.n_Allocation_Size();
	}
};

/**
 *	@brief v1 lambda solver utility function
 *	@tparam CDimsList is list of lambda matrix block sizes, as fbs_ut::CCTSize2D
 */
template <class CDimsList>
class CLambdaOps : public base_iface::CSolverOps_Base {
public:
	typedef CDimsList _TyLambdaMatrixBlockSizes; /**< @brief list of lambda matrix block sizes, as fbs_ut::CCTSize2D */
	struct _TyReductionPlan {}; /**< @brief reduction plan type (an empty class) */ // the v1 lambda did not need reduction plan, it was stored in the edges / vertices

	/**
	 *	@brief function object that calls lambda hessian block allocation for all edges
	 */
	class CAlloc_HessianBlocks {
	protected:
		CUberBlockMatrix &m_r_lambda; /**< @brief reference to the lambda matrix (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_lambda is reference to the lambda matrix
		 */
		inline CAlloc_HessianBlocks(CUberBlockMatrix &r_lambda)
			:m_r_lambda(r_lambda)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyVertexOrEdge is vertex or edge type
		 *	@param[in,out] r_vertex_or_edge is vertex or edge to have hessian blocks allocated in lambda
		 *	@note This function throws std::bad_alloc.
		 */
		template <class _TyVertexOrEdge>
		inline void operator ()(_TyVertexOrEdge &r_vertex_or_edge) // throw(std::bad_alloc)
		{
			r_vertex_or_edge.Alloc_HessianBlocks(m_r_lambda);
		}
	};

	/**
	 *	@brief function object that calculates hessians in all the edges
	 */
	class CCalculate_Hessians {
	public:
		/**
		 *	@brief function operator
		 *	@tparam _TyVertexOrEdge is vertex or edge type
		 *	@param[in] r_t_vertex_or_edge is vertex or edge to update it's hessians
		 */
		template <class _TyVertexOrEdge>
		inline void operator ()(_TyVertexOrEdge &r_t_vertex_or_edge) const
		{
			r_t_vertex_or_edge.Calculate_Hessians();
		}
	};

	/**
	 *	@brief function object that calls b vector calculation for all edges
	 */
	class CCollect_RightHandSide_Vector {
	protected:
		Eigen::VectorXd &m_r_b; /**< @brief reference to the right-hand side vector (out) */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_b is reference to the right-hand side vector
		 */
		inline CCollect_RightHandSide_Vector(Eigen::VectorXd &r_b)
			:m_r_b(r_b)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyVertex is vertex type
		 *	@param[in,out] r_t_vertex is a vertex to output its part R error vector
		 */
		template <class _TyVertex>
		inline void operator ()(const _TyVertex &r_t_vertex) // throw(std::bad_alloc)
		{
			r_t_vertex.Get_RightHandSide_Vector(m_r_b);
		}
	};

public:
	/**
	 *	@brief incrementally updates the lambda matrix structure (can be empty)
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in,out] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan (unused in v1)
	 *	@param[in,out] r_lambda is reference to the lambda matrix
	 *	@param[in] n_vertices_already_in_lambda is number of vertices which are already in the matrix
	 *	@param[in] n_edges_already_in_lambda is number of edges which are already in the matrix
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	template <class CSystem>
	static inline void Extend_Lambda(CSystem &r_system, _TyReductionPlan &UNUSED(r_reduction_plan),
		CUberBlockMatrix &r_lambda, size_t n_vertices_already_in_lambda, size_t n_edges_already_in_lambda) // throw(std::bad_alloc)
	{
		if(!n_vertices_already_in_lambda && !n_edges_already_in_lambda)
			AddEntriesInSparseSystem(r_system, r_lambda); // works for empty
		else
			UpdateSparseSystem(r_system, r_lambda, n_vertices_already_in_lambda, n_edges_already_in_lambda); // does not work for empty
		// create block matrix lambda

#ifdef _DEBUG
		/*{
			CUberBlockMatrix A;
			const Eigen::MatrixXd &r_t_uf = r_system.r_t_Unary_Factor();
			if(!A.Append_Block(r_t_uf, 0, 0))
				throw std::bad_alloc();
			// add unary factor

			r_system.r_Edge_Pool().For_Each(CAlloc_JacobianBlocks(A));
			// add all the hessian blocks

			CUberBlockMatrix lambda_ref;
			A.PreMultiplyWithSelfTransposeTo(lambda_ref, true); // only upper diag!
			// calculate lambda = AtA

			if(!r_lambda.b_EqualStructure(lambda_ref)) {
				lambda_ref.Rasterize("lambda1_reference_structure.tga");
				r_lambda.Rasterize("lambda0_structure.tga");
			}

			_ASSERTE(r_lambda.b_EqualStructure(lambda_ref));
			// make sure the matrix has the same structure
		}*/
#endif // _DEBUG
	}

	/**
	 *	@brief refreshes the lambda matrix by recalculating edge hessians
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in,out] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan (unused in v1)
	 *	@param[in,out] r_lambda is reference to the lambda matrix
	 *	@param[in] n_referesh_from_vertex is zero-based index of the first vertex to refresh
	 *	@param[in] n_refresh_from_edge is zero-based index of the first edge to refresh
	 */
	template <class CSystem>
	static inline void Refresh_Lambda(CSystem &r_system, _TyReductionPlan &UNUSED(r_reduction_plan),
		CUberBlockMatrix &r_lambda, size_t n_referesh_from_vertex = 0, size_t n_refresh_from_edge = 0)
	{
		if(n_refresh_from_edge) {
			r_system.r_Edge_Pool().For_Each_Parallel(n_refresh_from_edge,
				r_system.r_Edge_Pool().n_Size(), CCalculate_Hessians());
		} else {
			r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Hessians());
		}
		if(n_referesh_from_vertex) {
			r_system.r_Vertex_Pool().For_Each_Parallel(n_referesh_from_vertex,
				r_system.r_Vertex_Pool().n_Size(), CCalculate_Hessians());
		} else {
			r_system.r_Vertex_Pool().For_Each_Parallel(CCalculate_Hessians());
		}
		// can do this in parallel

		if(!CSystem::null_UnaryFactor && !n_referesh_from_vertex) {
			const Eigen::MatrixXd &r_t_uf = r_system.r_t_Unary_Factor();
			r_lambda.t_FindBlock(0, 0).noalias() += r_t_uf.transpose() * r_t_uf;
		}
		// add unary factor (gets overwritten by the first vertex' block)
#ifdef _DEBUG
		/*{
			CUberBlockMatrix A;
			const Eigen::MatrixXd &r_t_uf = r_system.r_t_Unary_Factor();
			if(!A.Append_Block(r_t_uf, 0, 0))
				throw std::bad_alloc();
			// add unary factor

			r_system.r_Edge_Pool().For_Each(CAlloc_JacobianBlocks(A));
			// add all the hessian blocks

			r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Jacobians());
			// calculate the values as well

			CUberBlockMatrix lambda_ref;
			A.PreMultiplyWithSelfTransposeTo(lambda_ref, true); // only upper diag!
			// calculate lambda = AtA

			if(!r_lambda.b_Equal(lambda_ref, 1e-3)) {
				r_lambda.Rasterize("lambda2_values.tga");
				lambda_ref.Rasterize("lambda3_reference_values.tga");
				CUberBlockMatrix &diff = lambda_ref;
				r_lambda.AddTo(diff, -1);
				diff.Rasterize("lambda4_diff_values.tga");
				fprintf(stderr, "error: lambda and it's reference have different value\n");
				exit(-1);
			}

			_ASSERTE(r_lambda.b_EqualStructure(lambda_ref));
			// make sure the matrix has the same structure
		}*/
#endif // _DEBUG
	}

	/**
	 *	@brief calculates the right-hand side vector
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan (unused in v1)
	 *	@param[in,out] r_v_b is reference to the r.h.s. vector (needs to be
	 *		allocated by the caller to the appropriate dimension)
	 */
	template <class CSystem>
	static inline void Collect_RightHandSide_Vector(const CSystem &r_system,
		const _TyReductionPlan &UNUSED(r_reduction_plan), Eigen::VectorXd &r_v_b)
	{
		r_system.r_Vertex_Pool().For_Each_Parallel(CCollect_RightHandSide_Vector(r_v_b)); // can do this in parallel
		// collect b
	}

	/**
	 *	@brief calculates a segment of the right-hand side vector, corresponding to a range of vertices
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan (unused in v1)
	 *	@param[in,out] r_v_b is reference to the r.h.s. vector (needs to be
	 *		allocated by the caller to the appropriate dimension)
	 *	@param[in] n_begin is zero-based index of the first vertex to calculate the r.h.s. for
	 *	@param[in] n_end is zero-based index of one past the last vertex to calculate the r.h.s. for
	 */
	template <class CSystem>
	static inline void Collect_RightHandSide_Vector(const CSystem &r_system,
		const _TyReductionPlan &UNUSED(r_reduction_plan), Eigen::VectorXd &r_v_b, size_t n_begin, size_t n_end)
	{
		r_system.r_Vertex_Pool().For_Each_Parallel(n_begin, n_end,
			CCollect_RightHandSide_Vector(r_v_b)); // can do this in parallel
		// collect b
	}

	/**
	 *	@brief calculates a segment of the right-hand side vector, corresponding to a single vertex
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan (unused in v1)
	 *	@param[in,out] r_v_b is reference to the r.h.s. vector (needs to be
	 *		allocated by the caller to the appropriate dimension)
	 *	@param[in] n_vertex is zero-based index of the vertex to calculate the r.h.s. for
	 */
	template <class CSystem>
	static inline void Collect_RightHandSide_Vector(const CSystem &r_system,
		const _TyReductionPlan &r_reduction_plan, Eigen::VectorXd &r_v_b, size_t n_vertex)
	{
		r_system.r_Vertex_Pool().For_Each(n_vertex, n_vertex + 1,
			CCollect_RightHandSide_Vector(r_v_b)); // this may not be the most efficient way
		// collect b
	}

protected:
	/**
	 *	@brief creates the lambda matrix from scratch
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in,out] r_system is optimized system
	 *	@param[in,out] r_lambda is reference to the lambda matrix
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	template <class CSystem>
	static inline void AddEntriesInSparseSystem(CSystem &r_system, CUberBlockMatrix &r_lambda) // throw(std::bad_alloc)
	{
#if 0
		if(r_system.r_Edge_Pool().n_Size() > 1000) { // wins 2.42237 - 2.48938 = .06701 seconds on 10k.graph, likely more on larger graphs
			//printf("building large matrix from scratch ...\n"); // debug
			std::vector<size_t> row_cumsum_list(r_system.r_Edge_Pool().n_Size());
			/*std::vector<size_t>::iterator p_end_it =*/
				r_system.r_Edge_Pool().For_Each(CGetCumsums(row_cumsum_list));
			//_ASSERTE(p_end_it == row_cumsum_list.end());
			// collect cumsums

			CUberBlockMatrix tmp(row_cumsum_list.begin(),
				row_cumsum_list.end(), r_system.r_Vertex_Pool().n_Size());
			r_lambda.Swap(tmp);
			// use this one instead

			// todo - see if there are some row_reindex on 100k, fix it by collecting
			// cumsums and building matrix with that (proven to be faster before)
		} else
#endif // 0 // todo - need to write function that gets cumsums from vertices (it's not difficult)
		{
			//printf("building small matrix from scratch ...\n"); // debug
			r_lambda.Clear();
			// ...
		}

		if(!CSystem::null_UnaryFactor) {
			const Eigen::MatrixXd &r_t_uf = r_system.r_t_Unary_Factor();
			if(!r_lambda.Append_Block(Eigen::MatrixXd(r_t_uf.transpose() * r_t_uf), 0, 0))
				throw std::bad_alloc();
		}
		// add unary factor (actually UF^T * UF, but it's the same matrix)

		// note that the unary error cannot be easily added without introducing a dummy
		// edge that would add itself as a reference to vertex 0
		if(r_system.r_v_Unary_Error().squaredNorm() > 0)
			throw std::runtime_error("unary error is not supported by the v1 reduction plan");
		// this is slightly obsolete so we will not support it for now

		r_system.r_Edge_Pool().For_Each(CAlloc_HessianBlocks(r_lambda));
		r_system.r_Vertex_Pool().For_Each(CAlloc_HessianBlocks(r_lambda));
		// add all the hessian blocks

		//printf("building lambda from scratch finished\n"); // debug
	}

	/**
	 *	@brief incrementally updates the lambda matrix structure (must not be empty)
	 */
	template <class CSystem>
	static inline void UpdateSparseSystem(CSystem &r_system,
		CUberBlockMatrix &r_lambda, size_t n_skip_vertices, size_t n_skip_edges) // throw(std::bad_alloc)
	{
		_ASSERTE(r_lambda.n_Row_Num() > 0 && r_lambda.n_Column_Num() == r_lambda.n_Row_Num()); // make sure lambda is not empty
		r_system.r_Edge_Pool().For_Each(n_skip_edges,
			r_system.r_Edge_Pool().n_Size(), CAlloc_HessianBlocks(r_lambda));
		r_system.r_Vertex_Pool().For_Each(n_skip_vertices,
			r_system.r_Vertex_Pool().n_Size(), CAlloc_HessianBlocks(r_lambda));
		// add the hessian blocks of the new edges
	}
};

/**
 *	@brief v2 lambda solver utility function
 *	@tparam CDimsList is list of lambda matrix block sizes, as fbs_ut::CCTSize2D
 */
template <class CDimsList>
class CLambdaOps2 : public base_iface::CSolverOps_Base {
public:
	typedef CDimsList _TyLambdaMatrixBlockSizes; /**< @brief list of lambda matrix block sizes, as fbs_ut::CCTSize2D */
	typedef CLambdaReductionPlan<CDimsList> _TyReductionPlan; /**< @brief reduction plan type */

	/**
	 *	@brief function object that calls lambda hessian block allocation for all edges
	 */
	class CAlloc_HessianBlocks_v2 {
	protected:
		CUberBlockMatrix &m_r_lambda; /**< @brief reference to the lambda matrix (out) */
		_TyReductionPlan &m_r_redplan; /**< @brief reference to the reduction plan */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in,out] r_lambda is reference to the lambda matrix (modified once the function operator is invoked)
		 *	@param[in,out] r_redplan is reference to the reduction plan (modified once the function operator is invoked)
		 */
		inline CAlloc_HessianBlocks_v2(CUberBlockMatrix &r_lambda, _TyReductionPlan &r_redplan)
			:m_r_lambda(r_lambda), m_r_redplan(r_redplan)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in,out] r_t_edge is edge to have hessian blocks allocated in lambda
		 *	@note This function throws std::bad_alloc.
		 */
		template <class _TyEdge>
		inline void operator ()(_TyEdge &r_t_edge) // throw(std::bad_alloc)
		{
			r_t_edge.Alloc_HessianBlocks_v2(m_r_lambda, m_r_redplan);
		}
	};

	/**
	 *	@brief function object that calculates hessians in all the edges
	 */
	class CCalculate_Hessians_v2 {
	public:
		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in] r_t_edge is edge to update it's hessians
		 */
		template <class _TyEdge>
		inline void operator ()(_TyEdge &r_t_edge) const
		{
			r_t_edge.Calculate_Hessians_v2();
		}
	};

	/**
	 *	@brief function object that calculates hessians in the selected edges
	 */
	class CUpdate_Hessians_v2 {
	protected:
		const CUberBlockMatrix &m_r_lambda; /**< @brief reference to the lambda matrix */
		_TyReductionPlan &m_r_redplan; /**< @brief reference to the reduction plan */
		bool m_b_recalc; /**< @brief Jacobian recalculation flag */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in] r_lambda is reference to the lambda matrix
		 *	@param[in,out] r_redplan is reference to the reduction plan (executed once the function operator is invoked)
		 *	@param[in] b_recalc is Jacobian recalculation flag (if set, the edge Jacobians are recalculated prior to the reduction)
		 */
		inline CUpdate_Hessians_v2(const CUberBlockMatrix &r_lambda, _TyReductionPlan &r_redplan, bool b_recalc)
			:m_r_lambda(r_lambda), m_r_redplan(r_redplan), m_b_recalc(b_recalc)
		{}

		/**
		 *	@brief function operator
		 *	@tparam _TyEdge is edge type
		 *	@param[in] r_t_edge is edge to update it's hessians
		 */
		template <class _TyEdge>
		inline void operator ()(_TyEdge &r_t_edge) const
		{
			if(m_b_recalc) // may choose to do that in parallel earlier
				r_t_edge.Calculate_Hessians_v2(); // calculate
			r_t_edge.Reduce_Hessians_v2(m_r_lambda, m_r_redplan); // reduce (may overwrite an earlier reduce
			// of a shared jacobian, but that one was incomplete because r_t_edge.Calculate_Hessians_v2()
			// was not called yet on this edge)
		}
	};

	/**
	 *	@brief unary factor helper
	 *
	 *	The size of UF depends on the dimension of the first vertex (unknown),
	 *	it needs to be found at runtime in the block size typelist. This is a
	 *	callback for CTypelistItemBFind.
	 */
	class CAddUnaryFactor {
	protected:
		_TyReductionPlan &m_r_rp; /**< @brief reference to lambda and RHS reductor */
		const Eigen::MatrixXd &m_r_t_uf; /**< @brief const reference to the unary factor */
		const Eigen::VectorXd &m_r_t_uerr; /**< @brief const reference to the unary error */
		double *m_p_uf_block; /**< @brief pointer to the hessian block of the first vertex in lambda */
		size_t m_n_vertes_id; /**< @brief id of the anchor vertex the unary factor will be applied to */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in] r_rp is reference to lambda / RHS reductors
		 *	@param[in] r_t_uf is const reference to the unary factor
		 *	@param[in] r_t_uerr is const reference to the unary error
		 *	@param[in] p_uf_block is pointer to the hessian block of the first vertex in lambda
		 *	@param[in] n_vertes_id is id of the anchor vertex the unary factor will be applied to
		 */
		inline CAddUnaryFactor(_TyReductionPlan &r_rp, const Eigen::MatrixXd &r_t_uf,
			const Eigen::VectorXd &r_t_uerr, double *p_uf_block, size_t n_vertes_id)
			:m_r_rp(r_rp), m_r_t_uf(r_t_uf), m_r_t_uerr(r_t_uerr), m_p_uf_block(p_uf_block),
			m_n_vertes_id(n_vertes_id)
		{
			_ASSERTE(p_uf_block);
		}

		/**
		 *	@brief callback operator; adds the unary factor ro the reduction plan
		 *	@tparam C2DSize is size of the unary factor
		 */
		template <class C2DSize>
		void operator ()()
		{
			double *p_temp = m_r_rp.r_Lambda_ReductionPlan().template
				p_Diagonal_GetTempBlock<C2DSize>(m_n_vertes_id, m_n_vertes_id, m_p_uf_block);
			// get a temp reduction block for the unary factor

			typename CUberBlockMatrix::CMakeMatrixRef<C2DSize::n_row_num,
				C2DSize::n_column_num>::_Ty dest(p_temp);
			dest.noalias() = m_r_t_uf.transpose() * m_r_t_uf;
			// copy UF

			if(m_r_t_uerr.squaredNorm() > 0) {
				double *p_temp_vec = m_r_rp.r_RHS_ReductionPlan().template
					p_Get_ReductionBlock<C2DSize::n_row_num>(m_n_vertes_id);
				// get a temp reduction block for the unary error

				typename CUberBlockMatrix::CMakeMatrixRef<C2DSize::n_row_num, 1>::_Ty vec(p_temp_vec);
				vec = m_r_t_uerr;
				// copy the error
			}
			// if norm of the unary error is nonzero (only special applications, it usualy is zero)
		}

		/**
		 *	@brief adds unary factor to the reduction plan
		 *
		 *	@param[in] r_rp is reference to lambda reductor
		 *	@param[in] r_t_uf is const reference to the unary factor
		 *	@param[in] r_t_uerr is const reference to the unary error
		 *	@param[in] r_lambda is reference to lambda
		 *	@param[in] n_vertes_id is id of the anchor vertex the unary factor will be applied to
		 *
		 *	@note This can only be used when the structure of lambda is fully allocated.
		 */
		static void Add_UnaryFactor(_TyReductionPlan &r_rp, const Eigen::MatrixXd &r_t_uf,
			const Eigen::VectorXd &r_t_uerr, CUberBlockMatrix &r_lambda, size_t n_vertes_id)
		{
			double *p_UF_block = r_lambda.p_GetBlock_Log(n_vertes_id, n_vertes_id, r_t_uf.rows(), r_t_uf.cols(), true, false);
			CAddUnaryFactor add_uf(r_rp, r_t_uf, r_t_uerr, p_UF_block, n_vertes_id);
			CTypelistItemBFind<typename CSortTypelist<_TyLambdaMatrixBlockSizes,
				fbs_ut::CCompareSize2D_Less>::_TyResult, fbs_ut::CRuntimeCompareSize2D,
				std::pair<size_t, size_t>, CAddUnaryFactor>::FindExisting(std::make_pair(r_t_uf.rows(),
				r_t_uf.cols()), add_uf);
		}

		/**
		 *	@brief adds unary factor to the reduction plan
		 *
		 *	@param[in] r_rp is reference to lambda reductor
		 *	@param[in] r_t_uf is const reference to the unary factor
		 *	@param[in] r_t_uerr is const reference to the unary error
		 *	@param[in] r_lambda is reference to lambda
		 *	@param[in] n_vertes_id is id of the anchor vertex the unary factor will be applied to
		 *	@param[in] n_vertes_order is order of the vertex with id n_vertes_id
		 */
		static void Add_UnaryFactor(_TyReductionPlan &r_rp, const Eigen::MatrixXd &r_t_uf,
			const Eigen::VectorXd &r_t_uerr, CUberBlockMatrix &r_lambda, size_t n_vertes_id, size_t n_vertes_order)
		{
			double *p_UF_block = r_lambda.p_FindBlock(n_vertes_order, n_vertes_order,
				r_t_uf.rows(), r_t_uf.cols(), true, false);
			CAddUnaryFactor add_uf(r_rp, r_t_uf, r_t_uerr, p_UF_block, n_vertes_id);
			CTypelistItemBFind<typename CSortTypelist<_TyLambdaMatrixBlockSizes,
				fbs_ut::CCompareSize2D_Less>::_TyResult, fbs_ut::CRuntimeCompareSize2D,
				std::pair<size_t, size_t>, CAddUnaryFactor>::FindExisting(std::make_pair(r_t_uf.rows(),
				r_t_uf.cols()), add_uf);
		}
	};

public:
	/**
	 *	@brief incrementally updates the lambda matrix structure (can be empty)
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in,out] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan
	 *	@param[in,out] r_lambda is reference to the lambda matrix
	 *	@param[in] n_vertices_already_in_lambda is number of vertices which are already in the matrix
	 *	@param[in] n_edges_already_in_lambda is number of edges which are already in the matrix
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	template <class CSystem>
	static inline void Extend_Lambda(CSystem &r_system,
		_TyReductionPlan &r_reduction_plan, CUberBlockMatrix &r_lambda,
		size_t n_vertices_already_in_lambda, size_t n_edges_already_in_lambda) // throw(std::bad_alloc)
	{
		if(!n_vertices_already_in_lambda && !n_edges_already_in_lambda)
			AddEntriesInSparseSystem(r_system, r_reduction_plan, r_lambda); // works for empty
		else {
			UpdateSparseSystem(r_system, r_reduction_plan, r_lambda,
				n_vertices_already_in_lambda, n_edges_already_in_lambda); // does not work for empty
		}
		// create block matrix lambda
	}

	/**
	 *	@brief refreshes the lambda matrix by recalculating edge hessians
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in,out] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan
	 *	@param[in,out] r_lambda is reference to the lambda matrix
	 *	@param[in] n_refresh_from_vertex is zero-based index of the first vertex to refresh (unused)
	 *	@param[in] n_refresh_from_edge is zero-based index of the first edge to refresh
	 */
	template <class CSystem>
	static inline void Refresh_Lambda(CSystem &r_system, _TyReductionPlan &r_reduction_plan,
		CUberBlockMatrix &r_lambda, size_t UNUSED(n_refresh_from_vertex) = 0, size_t n_refresh_from_edge = 0)
	{
		size_t n_edge_num = r_system.r_Edge_Pool().n_Size();
		size_t n_new_edge_num = n_edge_num - n_refresh_from_edge;
		const size_t n_parallel_thresh = 50;
		if(n_refresh_from_edge) {
			if(n_new_edge_num > n_parallel_thresh)
				r_system.r_Edge_Pool().For_Each_Parallel(n_refresh_from_edge, n_edge_num, CCalculate_Hessians_v2(), 0); // always run in parallel
		} else
			r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Hessians_v2());
		// can do this in parallel

		if(n_refresh_from_edge) {
			if(n_new_edge_num > n_parallel_thresh) {
				r_system.r_Edge_Pool().For_Each/*_Parallel*/(n_refresh_from_edge,
					r_system.r_Edge_Pool().n_Size(), CUpdate_Hessians_v2(r_lambda, r_reduction_plan, false)/*, 0*/); // always run in parallel, the hessians are already calculated
				// reduce only, can do this in parallel as all the results of the reduction will be the
				// same for all the threads and write conflicts should not matter (maybe on ARM or other
				// strange archs this might be a problem)
			} else {
				r_system.r_Edge_Pool().For_Each(n_refresh_from_edge,
					r_system.r_Edge_Pool().n_Size(), CUpdate_Hessians_v2(r_lambda, r_reduction_plan, true)); // not in parallel, recalculate the hessians as well
				// calculate and reduce, *not* in parallel
			}
		} else
			r_reduction_plan.r_Lambda_ReductionPlan().ReduceAll(); // simple, parallel
		// run the reduction plan
	}

	/**
	 *	@brief calculates the right-hand side vector
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in] r_system is optimized system (unused)
	 *	@param[in] r_reduction_plan is reduction plan
	 *	@param[in,out] r_v_b is reference to the r.h.s. vector (needs to be
	 *		allocated by the caller to the appropriate dimension)
	 */
	template <class CSystem>
	static inline void Collect_RightHandSide_Vector(const CSystem &UNUSED(r_system),
		const _TyReductionPlan &r_reduction_plan, Eigen::VectorXd &r_v_b)
	{
		r_reduction_plan.r_RHS_ReductionPlan().ReduceAll(r_v_b);
		// collect b
	}

	/**
	 *	@brief calculates a segment of the right-hand side vector, corresponding to a range of vertices
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan
	 *	@param[in,out] r_v_b is reference to the r.h.s. vector (needs to be
	 *		allocated by the caller to the appropriate dimension)
	 *	@param[in] n_begin is zero-based index of the first vertex to calculate the r.h.s. for
	 *	@param[in] n_end is zero-based index of one past the last vertex to calculate the r.h.s. for
	 */
	template <class CSystem>
	static inline void Collect_RightHandSide_Vector(const CSystem &r_system,
		const _TyReductionPlan &r_reduction_plan, Eigen::VectorXd &r_v_b, size_t n_begin, size_t n_end)
	{
		if(n_end - n_begin > 50) {
			// todo - parallel implementation using ReduceSingle
		}

		n_begin = r_system.r_Vertex_Pool()[n_begin].n_Order();
		typename CSystem::_TyConstVertexRef last = r_system.r_Vertex_Pool()[n_end - 1];
		n_end = last.n_Order() + last.n_Dimension();
		// need to convert from vertex indices to element indices

		r_reduction_plan.r_RHS_ReductionPlan().ReduceRange(r_v_b, n_begin, n_end);
		// collect b
	}

protected:
	/**
	 *	@brief reduction context
	 */
	struct TSingleRHSReduceCtx {
		const _TyReductionPlan &r_reduction_plan; /**< @brief reference to the reduction plan */
		Eigen::VectorXd &r_v_b; /**< @brief reference to the r.h.s. vector */
		size_t n_order; /**< @brief order of the reduced vertex */

		/**
		 *	@brief default constructor
		 *
		 *	@param[in] _r_reduction_plan is reference to the reduction plan
		 *	@param[in] _r_v_b is reference to the r.h.s. vector
		 *	@param[in] _n_order is order of the reduced vertex
		 */
		inline TSingleRHSReduceCtx(const _TyReductionPlan &_r_reduction_plan, Eigen::VectorXd &_r_v_b, size_t _n_order)
			:r_reduction_plan(_r_reduction_plan), r_v_b(_r_v_b), n_order(_n_order)
		{}
	};

	/**
	 *	@brief signle vertex r.h.s. reduction functor
	 *	@tparam n_dimension is dimension of the given vertex
	 */
	template <const int n_dimension>
	class CSingleRHSReducer {
	public:
		/**
		 *	@brief performs a single vertex r.h.s. reduction
		 *	@param[in] t_ctx is reduction context (reduction plan, the r.h.s. vector and the vertex order)
		 */
		static inline void Do(TSingleRHSReduceCtx t_ctx)
		{
			t_ctx.r_reduction_plan.r_RHS_ReductionPlan().template Reduce_Single<n_dimension>(t_ctx.r_v_b, t_ctx.n_order);
		}
	};

public:
	/**
	 *	@brief calculates a segment of the right-hand side vector, corresponding to a single vertex
	 *
	 *	@tparam CSystem is optimized system type
	 *
	 *	@param[in] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan
	 *	@param[in,out] r_v_b is reference to the r.h.s. vector (needs to be
	 *		allocated by the caller to the appropriate dimension)
	 *	@param[in] n_vertex is zero-based index of the vertex to calculate the r.h.s. for
	 */
	template <class CSystem>
	static inline void Collect_RightHandSide_Vector(const CSystem &r_system,
		const _TyReductionPlan &r_reduction_plan, Eigen::VectorXd &r_v_b, size_t n_vertex)
	{
		typename CSystem::_TyConstVertexRef vertex = r_system.r_Vertex_Pool()[n_vertex];
		size_t n_order = vertex.n_Order();
		int n_dimension = int(vertex.n_Dimension());

		fbs_ut::CWrap<CSingleRHSReducer>::template In_ScalarSize_DecisionTree<typename
			_TyReductionPlan::CRHSReductor::_TyDimensionList,
			TSingleRHSReduceCtx>(n_dimension, TSingleRHSReduceCtx(r_reduction_plan, r_v_b, n_order));
	}

protected:
	/**
	 *	@brief creates the lambda matrix from scratch
	 *
	 *	@param[in,out] r_system is optimized system
	 *	@param[in] r_reduction_plan is reduction plan
	 *	@param[in,out] r_lambda is reference to the lambda matrix
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	template <class CSystem>
	static inline void AddEntriesInSparseSystem(CSystem &r_system,
		_TyReductionPlan &r_reduction_plan, CUberBlockMatrix &r_lambda) // throw(std::bad_alloc)
	{
		r_lambda.Clear();
		// ...

		r_system.r_Edge_Pool().For_Each(CAlloc_HessianBlocks_v2(r_lambda, r_reduction_plan));
		// add all the hessian blocks

		if(!CSystem::null_UnaryFactor) {
			size_t n_first_vertex_id = r_system.r_Edge_Pool()[0].n_Vertex_Id(0);
			size_t n_first_vertex_order = r_system.r_Vertex_Pool()[n_first_vertex_id].n_Order();
			// get id of the first vertex (usually zero)

			const Eigen::MatrixXd &r_t_uf = r_system.r_t_Unary_Factor();
			/*if(!r_lambda.Append_Block_Log(Eigen::MatrixXd(r_t_uf.transpose() * r_t_uf),
			   n_first_vertex_id, n_first_vertex_id))
				throw std::bad_alloc();*/ // no need to do this anymore
			// add unary factor (actually UF^T * UF, but it's the same matrix)

			CAddUnaryFactor::Add_UnaryFactor(r_reduction_plan, r_t_uf,
				r_system.r_v_Unary_Error(), r_lambda, n_first_vertex_id, n_first_vertex_order);
		}
		// add unary factor to the reductor so that we don't have to care about ir anymore
		// do this after the edges, as we are using Append_Block_Log() and hence we need
		// the row / column with n_first_vertex_id to exist (and if n_first_vertex_id > 0,
		// then it might not)
	}

	/**
	 *	@brief incrementally updates the lambda matrix structure (must not be empty)
	 */
	template <class CSystem>
	static inline void UpdateSparseSystem(CSystem &r_system,
		_TyReductionPlan &r_reduction_plan, CUberBlockMatrix &r_lambda,
		size_t n_skip_vertices, size_t n_skip_edges) // throw(std::bad_alloc)
	{
		_ASSERTE(r_lambda.n_Row_Num() > 0 && r_lambda.n_Column_Num() == r_lambda.n_Row_Num()); // make sure lambda is not empty
		r_system.r_Edge_Pool().For_Each(n_skip_edges,
			r_system.r_Edge_Pool().n_Size(), CAlloc_HessianBlocks_v2(r_lambda, r_reduction_plan));
		// add the hessian blocks of the new edges
	}
};

} // ~lambda_utils

/**
 *	@brief nonlinear blocky solver working above the lambda matrix
 *
 *	@tparam CSystem is optimization system type
 *	@tparam CLinearSolver is linear solver type
 *	@tparam CAMatrixBlockSizes is list of block sizes in the Jacobian matrix
 *	@tparam CLambdaMatrixBlockSizes is list of block sizes in the information (Hessian) matrix
 */
template <class CSystem, class CLinearSolver, class CAMatrixBlockSizes = typename CSystem::_TyJacobianMatrixBlockList,
	class CLambdaMatrixBlockSizes = typename CSystem::_TyHessianMatrixBlockList>
class CNonlinearSolver_Lambda {
public:
	typedef CSystem _TySystem; /**< @brief system type */
	typedef CLinearSolver _TyLinearSolver; /**< @brief linear solver type */

	typedef typename CSystem::_TyBaseVertex _TyBaseVertex; /**< @brief the data type for storing vertices */
	typedef typename CSystem::_TyVertexTypelist _TyVertexTypelist; /**< @brief list of vertex types */
	typedef typename CSystem::_TyBaseEdge _TyBaseEdge; /**< @brief the data type for storing measurements */
	typedef typename CSystem::_TyEdgeTypelist _TyEdgeTypelist; /**< @brief list of edge types */

	typedef typename CSystem::_TyVertexMultiPool _TyVertexMultiPool; /**< @brief vertex multipool type */
	typedef typename CSystem::_TyEdgeMultiPool _TyEdgeMultiPool; /**< @brief edge multipool type */

	typedef typename CLinearSolver::_Tag _TySolverTag; /**< @brief linear solver tag */
	typedef CLinearSolverWrapper<_TyLinearSolver, _TySolverTag> _TyLinearSolverWrapper; /**< @brief wrapper for linear solvers (shields solver capability to solve blockwise) */

	typedef /*typename CUniqueTypelist<*/CAMatrixBlockSizes/*>::_TyResult*/ _TyAMatrixBlockSizes; /**< @brief possible block matrices, that can be found in A */
	typedef /*typename*/ CLambdaMatrixBlockSizes /*fbs_ut::CBlockMatrixTypesAfterPreMultiplyWithSelfTranspose<
		_TyAMatrixBlockSizes>::_TySizeList*/ _TyLambdaMatrixBlockSizes; /**< @brief possible block matrices, found in lambda and R */

	typedef typename CChooseType<lambda_utils::CLambdaOps<_TyLambdaMatrixBlockSizes>,
		lambda_utils::CLambdaOps2<_TyLambdaMatrixBlockSizes>, !base_iface::lambda_ReductionPlan_v2>::_TyResult _TyLambdaOps; /**< @brief implementation of operations for filling the lambda matrix */
	typedef typename _TyLambdaOps::_TyReductionPlan _TyReductionPlan; /**< @brief reduction plan implementation */

	/**
	 *	@brief solver interface properties, stored as enum (see also CSolverTraits)
	 */
	enum {
		solver_HasDump = true, /**< @brief timing statistics support flag */
		solver_HasChi2 = true, /**< @brief Chi2 error calculation support flag */
		solver_HasMarginals = true, /**< @brief marginal covariance support flag */
		solver_HasGaussNewton = true, /**< @brief Gauss-Newton support flag */
		solver_HasLevenberg = false, /**< @brief Levenberg-Marquardt support flag */
		solver_HasGradient = false, /**< @brief gradient-based linear solving support flag */
		solver_HasSchur = true, /**< @brief Schur complement support flag */
		solver_HasDelayedOptimization = false, /**< @brief delayed optimization support flag */
		solver_IsPreferredBatch = true, /**< @brief preferred batch solver flag */
		solver_IsPreferredIncremental = false, /**< @brief preferred incremental solver flag */
		solver_ExportsJacobian = false, /**< @brief interface for exporting jacobian system matrix flag */
		solver_ExportsHessian = true, /**< @brief interface for exporting hessian system matrix flag */
		solver_ExportsFactor = false /**< @brief interface for exporting factorized system matrix flag */
	};

protected:
	CSystem &m_r_system; /**< @brief reference to the system */
	CLinearSolver m_linear_solver; /**< @brief linear solver */
	CLinearSolver_Schur<CLinearSolver, _TyAMatrixBlockSizes, CSystem> m_schur_solver; /**< @brief linear solver with Schur trick */

	CUberBlockMatrix m_lambda; /**< @brief the lambda matrix (built / updated incrementally) */
	_TyReductionPlan m_reduction_plan; /**< @brief lambda incremental reduction plan */
	Eigen::VectorXd m_v_dx; /**< @brief dx vector */
	size_t m_n_verts_in_lambda; /**< @brief number of vertices already in lambda */
	size_t m_n_edges_in_lambda; /**< @brief number of edges already in lambda */
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
	double m_f_lambda_time; /**< @brief time spent updating lambda */
	double m_f_rhs_time; /**< @brief time spent gathering the right-hand-side vector */
	double m_f_chol_time; /**< @brief time spent in Choleski() section */
	double m_f_norm_time; /**< @brief time spent in norm calculation section */
	double m_f_vert_upd_time; /**< @brief time spent updating the system */

	TMarginalsComputationPolicy m_t_marginals_config; /**< @brief marginal covariance computation configuration */
	CMarginalCovariance m_marginals; /**< @brief marginals cache */
	double m_f_extra_chol_time; /**< @brief time spent in calculating extra Cholesky factorization for marginal covariances */
	double m_f_marginals_time; /**< @brief time spent in calculating marginal covariances (batch) */
	double m_f_incmarginals_time; /**< @brief time spent in calculating marginal covariances (update) */
	size_t m_n_incmarginals_num; /**< @brief number of times the marginals update ran instead of batch recalculation */

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
	 *
	 *	@deprecated This is deprecated version of the constructor, use constructor
	 *		with TIncrementalSolveSetting instead.
	 */
	CNonlinearSolver_Lambda(CSystem &r_system, size_t n_linear_solve_threshold,
		size_t n_nonlinear_solve_threshold, size_t n_nonlinear_solve_max_iteration_num = 5,
		double f_nonlinear_solve_error_threshold = .01, bool b_verbose = false,
		CLinearSolver linear_solver = CLinearSolver(), bool b_use_schur = true)
		:m_r_system(r_system), m_linear_solver(linear_solver), m_schur_solver(linear_solver),
		m_n_verts_in_lambda(0), m_n_edges_in_lambda(0),
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
		m_b_system_dirty(false), m_n_iteration_num(0), m_f_lambda_time(0), m_f_rhs_time(0),
		m_f_chol_time(0), m_f_norm_time(0), m_f_vert_upd_time(0), m_b_had_loop_closure(false),
		m_f_extra_chol_time(0), m_f_marginals_time(0), m_f_incmarginals_time(0), m_n_incmarginals_num(0)
	{
		_ASSERTE(!m_n_nonlinear_solve_threshold || !m_n_linear_solve_threshold); // only one of those
	}

	/**
	 *	@brief initializes the nonlinear solver
	 *
	 *	@param[in] r_system is the system to be optimized
	 *		(it is only referenced, not copied - must not be deleted)
	 *	@param[in] t_incremental_config is incremental solving configuration
	 *	@param[in] t_marginals_config is marginal covariance calculation configuration
	 *	@param[in] b_verbose is verbosity flag
	 *	@param[in] linear_solver is linear solver instance
	 *	@param[in] b_use_schur is Schur complement trick flag
	 */
	CNonlinearSolver_Lambda(CSystem &r_system,
		TIncrementalSolveSetting t_incremental_config = TIncrementalSolveSetting(),
		TMarginalsComputationPolicy t_marginals_config = TMarginalsComputationPolicy(),
		bool b_verbose = false,
		CLinearSolver linear_solver = CLinearSolver(), bool b_use_schur = false)
		:m_r_system(r_system), m_linear_solver(linear_solver), m_schur_solver(linear_solver),
		m_n_verts_in_lambda(0), m_n_edges_in_lambda(0),
#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_last_optimized_vertex_num(0),
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_step(0),
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		m_n_linear_solve_threshold(t_incremental_config.t_linear_freq.n_period),
		m_n_nonlinear_solve_threshold(t_incremental_config.t_nonlinear_freq.n_period),
		m_n_nonlinear_solve_max_iteration_num(t_incremental_config.n_max_nonlinear_iteration_num),
		m_f_nonlinear_solve_error_threshold(t_incremental_config.f_nonlinear_error_thresh),
		m_b_verbose(b_verbose), m_b_use_schur(b_use_schur), m_n_real_step(0),
		m_b_system_dirty(false), m_n_iteration_num(0), m_f_lambda_time(0), m_f_rhs_time(0),
		m_f_chol_time(0), m_f_norm_time(0), m_f_vert_upd_time(0), m_b_had_loop_closure(false),
		m_t_marginals_config(t_marginals_config), m_f_extra_chol_time(0), m_f_marginals_time(0),
		m_f_incmarginals_time(0), m_n_incmarginals_num(0)
	{
		_ASSERTE(!m_n_nonlinear_solve_threshold || !m_n_linear_solve_threshold); // only one of those
		//_ASSERTE(!t_marginals_config.b_calculate); // not supported at the moment // already partially supported

		if(t_marginals_config.b_calculate) {
			if(t_marginals_config.t_increment_freq.n_period != t_incremental_config.t_nonlinear_freq.n_period &&
			   t_marginals_config.t_increment_freq.n_period != t_incremental_config.t_linear_freq.n_period) {
				throw std::runtime_error("in CNonlinearSolver_Lambda, the marginals must"
					" be updated with the same frequency as the system");
			}
			// unfortunately, yes

			/*if(t_marginals_config.n_incremental_policy != (mpart_LastColumn | mpart_Diagonal)) {
				throw std::runtime_error("in CNonlinearSolver_Lambda, the marginals update"
					" policy must be mpart_LastColumn | mpart_Diagonal");
			}
			if(t_marginals_config.n_incremental_policy != t_marginals_config.n_relinearize_policy) {
				throw std::runtime_error("in CNonlinearSolver_Lambda, the marginals "
					" incremental and relinearize update policy must be the same");
			}*/ // these are now implemented
			if(t_marginals_config.n_cache_miss_policy != mpart_Nothing) {
				throw std::runtime_error("in CNonlinearSolver_Lambda, the marginals cache"
					" miss policy is not supported at the moment, sorry for inconvenience");
			}
			// nothing else is implemented so far
		}
	}

	/**
	 *	@brief gets the current system matrix (not necessarily up-to-date)
	 *	@return Returns const reference to the current system matrix.
	 */
	inline const CUberBlockMatrix &r_Lambda() const
	{
		return m_lambda;
	}

	/**
	 *	@brief gets marginal covariances
	 *	@return Returns reference to the marginal covariances object.
	 */
	inline CMarginalCovariance &r_MarginalCovariance()
	{
		return m_marginals;
	}

	/**
	 *	@brief gets marginal covariances
	 *	@return Returns const reference to the marginal covariances object.
	 */
	inline const CMarginalCovariance &r_MarginalCovariance() const
	{
		return m_marginals;
	}

	/**
	 *	@brief displays performance info on stdout
	 *	@param[in] f_total_time is total time taken by everything (can be -1 to ommit)
	 */
	void Dump(double f_total_time = -1) const
	{
		printf("solver took " PRIsize " iterations\n", m_n_iteration_num); // debug, to be able to say we didn't botch it numerically
		printf("solver spent %f seconds in parallelizable section (disparity %f seconds)\n",
			m_f_lambda_time + m_f_rhs_time, (f_total_time > 0)? f_total_time -
			(m_f_lambda_time + m_f_rhs_time + m_f_chol_time + m_f_norm_time +
			m_f_vert_upd_time + m_f_extra_chol_time + m_f_marginals_time + m_f_incmarginals_time) : 0);

		printf("out of which:\n");
		printf("\t   ,\\: %f\n", m_f_lambda_time);
		printf("\t  rhs: %f\n", m_f_rhs_time);
		if(m_t_marginals_config.b_calculate) {
			printf("solver spent %f seconds in marginals\n"
				"\t chol: %f\n"
				"\tmargs: %f\n"
				"\t incm: %f (ran " PRIsize " times)\n",
				m_f_extra_chol_time + m_f_marginals_time + m_f_incmarginals_time,
				m_f_extra_chol_time, m_f_marginals_time,
				m_f_incmarginals_time, m_n_incmarginals_num);
		}
		printf("solver spent %f seconds in serial section\n",
			m_f_chol_time + m_f_norm_time + m_f_vert_upd_time);
		printf("out of which:\n");
		printf("\t chol: %f\n", m_f_chol_time);
		printf("\t norm: %f\n", m_f_norm_time);
		printf("\tv-upd: %f\n", m_f_vert_upd_time);
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
			_TyLambdaOps::Extend_Lambda(m_r_system, m_reduction_plan, m_lambda,
				m_n_verts_in_lambda, m_n_edges_in_lambda); // recalculated all the jacobians inside Extend_Lambda()
			if(!m_b_system_dirty) {
				_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda,
					m_n_verts_in_lambda, m_n_edges_in_lambda); // calculate only for new edges // @todo - but how to mark affected vertices? // simple test if edge id is greater than m_n_edges_in_lambda, the vertex needs to be recalculated
			} else
				_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda); // calculate for entire system
			m_b_system_dirty = false;
			m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
			m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
			//const size_t n_variables_size = m_r_system.n_VertexElement_Num();
			_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
				m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
				m_lambda.n_Column_Num() == m_r_system.n_VertexElement_Num()); // lambda is square, blocks on either side = number of vertices
		} catch(std::bad_alloc&) {
			return false;
		}
		// need to have lambda

		return m_lambda.Rasterize_Symmetric(p_s_filename, n_scalar_size);
	}

	/**
	 *	@brief writes system matrix in matrix market for benchmarking purposes
	 *	@param[in] p_s_filename is output file name (.mtx)
	 *	@return Returns true on success, false on failure.
	 */
	bool Save_SystemMatrix_MM(const char *p_s_filename) const
	{
		char p_s_layout_file[256];
		strcpy(p_s_layout_file, p_s_filename);
		if(strrchr(p_s_layout_file, '.'))
			*(char*)strrchr(p_s_layout_file, '.') = 0;
		strcat(p_s_layout_file, ".bla");
		// only really required for landmark datasets

		return m_lambda.Save_MatrixMarket(p_s_filename, p_s_layout_file, "lambda matrix for SLAM problem");
	}

	/**
	 *	@brief calculates chi-squared error
	 *	@return Returns chi-squared error.
	 *	@note This only works with systems with edges of one degree of freedom
	 *		(won't work for e.g. systems with both poses and landmarks).
	 */
	double f_Chi_Squared_Error() const
	{
		return _TyLambdaOps::f_Chi_Squared_Error(m_r_system);
	}

	/**
	 *	@brief calculates denormalized chi-squared error
	 *	@return Returns denormalized chi-squared error.
	 *	@note This doesn't perform the final division by (number of edges - degree of freedoms).
	 */
	double f_Chi_Squared_Error_Denorm() const
	{
		return _TyLambdaOps::f_Chi_Squared_Error_Denorm(m_r_system);
	}

	/**
	 *	@brief incremental optimization function
	 *	@param[in] r_last_edge is the last edge that was added to the system
	 *	@note This function throws std::bad_alloc.
	 */
	void Incremental_Step(_TyBaseEdge &UNUSED(r_last_edge)) // throw(std::bad_alloc)
	{
		size_t n_vertex_num = m_r_system.r_Vertex_Pool().n_Size();
		if(!m_b_had_loop_closure) {
			_ASSERTE(!m_r_system.r_Edge_Pool().b_Empty());
			typename _TyEdgeMultiPool::_TyConstBaseRef r_edge =
				m_r_system.r_Edge_Pool()[m_r_system.r_Edge_Pool().n_Size() - 1];
			// get a reference to the last edge interface

			if(r_edge.n_Vertex_Num() > 1) { // unary factors do not cause classical loop closures
				_ASSERTE(r_edge.n_Vertex_Id(0) != r_edge.n_Vertex_Id(1));
				size_t n_first_vertex = std::min(r_edge.n_Vertex_Id(0), r_edge.n_Vertex_Id(1));
				m_b_had_loop_closure = (n_first_vertex < n_vertex_num - 2);
				_ASSERTE(m_b_had_loop_closure || std::max(r_edge.n_Vertex_Id(0),
					r_edge.n_Vertex_Id(1)) == n_vertex_num - 1);
			} else {
				size_t n_first_vertex = r_edge.n_Vertex_Id(0);
				m_b_had_loop_closure = (n_first_vertex < n_vertex_num - 1);
			}
			// todo - encapsulate this code, write code to support hyperedges as well, use it
		}
		// detect loop closures (otherwise the edges are initialized based on measurement and error would be zero)

		++ m_n_real_step;

#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
		size_t n_new_vertex_num = n_vertex_num - m_n_last_optimized_vertex_num;
		// the optimization periods are counted in vertices, not in edges (should save work on 10k, 100k)
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		++ m_n_step;
		size_t n_new_edge_num = m_n_step;
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES

		IncrementalPreOptDemo();

		bool b_new_vert = false, b_ran_opt = false;

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
				b_ran_opt = true;
				Optimize((m_b_had_loop_closure)?
					m_n_nonlinear_solve_max_iteration_num : 0,
					m_f_nonlinear_solve_error_threshold);
				m_b_had_loop_closure = false;
			}
			// nonlinear optimization

			b_new_vert = true;
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
				b_ran_opt = true;
				m_b_had_loop_closure = false;
				Optimize(1, 0); // only if there was a loop (ignores possibly high residual after single step optimization)
			}
			// simple optimization

			b_new_vert = true;
		}

		if(b_new_vert && !b_ran_opt && m_t_marginals_config.b_calculate)
			Optimize(0, 0);
		// run optimization in order to calculate marginals after each vertex

#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
		IncrementalPostOptDemo(n_new_vertex_num);
#else // __SLAM_COUNT_ITERATIONS_AS_VERTICES
		IncrementalPostOptDemo(n_new_edge_num);
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
	}

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
	CMatrixOrdering mord100;
	CUberBlockMatrix lambda_prev;
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

	/**
	 *	@brief final optimization function
	 *
	 *	@param[in] n_max_iteration_num is the maximal number of iterations
	 *	@param[in] f_min_dx_norm is the residual norm threshold
	 */
	void Optimize(size_t n_max_iteration_num = 5, double f_min_dx_norm = .01) // throw(std::bad_alloc)
	{
		CTimerSampler timer(m_timer);

		const size_t n_variables_size = m_r_system.n_VertexElement_Num();
		const size_t n_measurements_size = m_r_system.n_EdgeElement_Num();
		if(n_variables_size > n_measurements_size) {
			if(n_measurements_size)
				fprintf(stderr, "warning: the system is underspecified\n");
			else
				fprintf(stderr, "warning: the system contains no edges at all: nothing to optimize\n");
			//return;
		}
		if(!n_measurements_size)
			return; // nothing to solve (but no results need to be generated so it's ok)
		// can't solve in such conditions
		_TyLambdaOps::Extend_Lambda(m_r_system, m_reduction_plan,
			m_lambda, m_n_verts_in_lambda, m_n_edges_in_lambda); // recalculated all the jacobians inside Extend_Lambda()
		if(!m_b_system_dirty) {
			_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda,
				m_n_verts_in_lambda, m_n_edges_in_lambda); // calculate only for new edges // t_odo - but how to mark affected vertices?
		} else
			_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda); // calculate for entire system
		m_b_system_dirty = false;
		m_n_verts_in_lambda = m_r_system.r_Vertex_Pool().n_Size();
		m_n_edges_in_lambda = m_r_system.r_Edge_Pool().n_Size(); // right? // yes.
		// need to have lambda

		if(m_lambda.n_BlockColumn_Num() < m_r_system.r_Vertex_Pool().n_Size()) {
			timer.Accum_DiffSample(m_f_lambda_time);
			fprintf(stderr, "warning: waiting for more edges\n");
			return;
		}
		// waiting for more edges

		_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
			m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
			m_lambda.n_Column_Num() == n_variables_size); // lambda is square, blocks on either side = number of vertices
		// invariants

		Dump_DemoData();

		m_v_dx.resize(n_variables_size, 1);

		if(m_b_use_schur)
			m_schur_solver.SymbolicDecomposition_Blocky(m_lambda);
		// calculate the ordering once, it does not change

		if(m_b_verbose) {
			size_t n_sys_size = m_r_system.n_Allocation_Size();
			size_t n_rp_size = m_reduction_plan.n_Allocation_Size();
			size_t n_lam_size = m_lambda.n_Allocation_Size();
			printf("memory_use(sys: %.2f MB, redplan: %.2f MB, Lambda: %.2f MB)\n",
				n_sys_size / 1048576.0, n_rp_size / 1048576.0, n_lam_size / 1048576.0);
		}
		// print memory use statistics

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
				_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda);
				m_b_system_dirty = false;
			}
			// no need to rebuild lambda, just refresh the values that are being referenced

			timer.Accum_DiffSample(m_f_lambda_time);

			_TyLambdaOps::Collect_RightHandSide_Vector(m_r_system, m_reduction_plan, m_v_dx);
			// collects the right-hand side vector

			timer.Accum_DiffSample(m_f_rhs_time);

#ifdef _DEBUG
			/*_ASSERTE(m_lambda.n_Row_Num() == n_variables_size); // should be the same
			// collect errors (in parallel)

			{
				CUberBlockMatrix A;
				{
					const Eigen::MatrixXd &r_t_uf = m_r_system.r_t_Unary_Factor();
					if(!A.Append_Block(r_t_uf, 0, 0))
						throw std::bad_alloc();
					// add unary factor

					m_r_system.r_Edge_Pool().For_Each(CAlloc_JacobianBlocks(A));
					// add all the hessian blocks

					m_r_system.r_Edge_Pool().For_Each_Parallel(CCalculate_Jacobians());
					// calculate the values as well
				}
				// calculate A

				Eigen::VectorXd v_error(n_measurements_size);
				Collect_R_Errors(v_error);
				Eigen::VectorXd v_eta(n_variables_size);
				v_eta.setZero(); // eta = 0
				Eigen::VectorXd &v_b = v_error; // b = Rz * p_errors
				_ASSERTE(v_eta.rows() == n_variables_size);
				_ASSERTE(v_b.rows() == n_measurements_size);
				A.PostMultiply_Add_FBS<_TyAMatrixBlockSizes>(&v_eta(0), n_variables_size, &v_b(0),
					n_measurements_size); // works (fast parallel post-multiply)
				// calculates eta

				_ASSERTE(v_eta.rows() == m_v_dx.rows());
				for(size_t i = 0, n = v_eta.rows(); i < n; ++ i) {
					if(fabs(m_v_dx(i) - v_eta(i)) > 1e-5)
						fprintf(stderr, "error: RHS vectors do not match\n");
				}
				// compare
			}*/
			// calculate eta = A^T * b
#endif // _DEBUG

			{
				bool b_cholesky_result;
				{
					Eigen::VectorXd &v_eta = m_v_dx; // dx is calculated inplace from eta
					if(!m_b_use_schur) {
						if(n_max_iteration_num > 1) {
							do {
								if(!n_iteration &&
								   !_TyLinearSolverWrapper::FinalBlockStructure(m_linear_solver, m_lambda)) {
									b_cholesky_result = false;
									break;
								}
								// prepare symbolic factorization, structure of lambda won't change in the next steps
								b_cholesky_result = _TyLinearSolverWrapper::Solve(m_linear_solver, m_lambda, v_eta);
								// p_dx = eta = lambda / eta
							} while(0);
						} else
							b_cholesky_result = m_linear_solver.Solve_PosDef(m_lambda, v_eta); // p_dx = eta = lambda / eta
					} else { // use Schur complement
						b_cholesky_result = m_schur_solver.Solve_PosDef_Blocky(m_lambda, v_eta);
						// Schur
					}

					if(m_b_verbose) {
						printf("%s %s", (m_b_use_schur)? "Schur" : "Cholesky",
							(b_cholesky_result)? "succeeded\n" : "failed\n");
					}
				}
				// calculate cholesky, reuse block ordering if the linear solver supports it

				timer.Accum_DiffSample(m_f_chol_time);

				double f_residual_norm = 0;
				if(b_cholesky_result) {
					f_residual_norm = m_v_dx.norm(); // Eigen likely uses SSE and OpenMP
					if(m_b_verbose)
						printf("residual norm: %.4f\n", f_residual_norm);
				}
				// calculate residual norm

				timer.Accum_DiffSample(m_f_norm_time);

				if(f_residual_norm <= f_min_dx_norm)
					break;
				// in case the error is low enough, quit (saves us recalculating the hessians)

				if(b_cholesky_result) {
					_TyLambdaOps::PushValuesInGraphSystem(m_r_system, m_v_dx);
					m_b_system_dirty = true;

					//printf("debug: just updated the linearization point\n");
					m_marginals.DisableUpdate(); // linearization point just changed, all the marginals will change - need full recalc

					timer.Accum_DiffSample(m_f_vert_upd_time);
				}
				// update the system (in parallel)

				if(!b_cholesky_result)
					break;
				// in case cholesky failed, quit
			}
		}
		// optimize the system

		if(m_t_marginals_config.b_calculate) {
			// todo - handle freq settings
			// todo - handle policies

			if(m_b_system_dirty) {
				_TyLambdaOps::Refresh_Lambda(m_r_system, m_reduction_plan, m_lambda);
				m_b_system_dirty = false; // does this break something? seems not.

				timer.Accum_DiffSample(m_f_lambda_time);
			}
			_ASSERTE(m_n_verts_in_lambda == m_r_system.r_Vertex_Pool().n_Size());
			_ASSERTE(m_n_edges_in_lambda == m_r_system.r_Edge_Pool().n_Size());
			_ASSERTE(m_lambda.n_Row_Num() == m_lambda.n_Column_Num() &&
				m_lambda.n_BlockColumn_Num() == m_r_system.r_Vertex_Pool().n_Size() &&
				m_lambda.n_Column_Num() == n_variables_size);
			// need to update or will end up with forever bad marginals!

			CUberBlockMatrix R;
			//R.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(m_lambda); // this makes ugly dense factor, dont do that

			//m_linear_solver.Factorize_PosDef_Blocky(R, m_lambda, std::vector<size_t>()); // dense as well, no ordering inside

			CMatrixOrdering mord;
			if((m_marginals.b_CanUpdate() && (m_t_marginals_config.n_incremental_policy &
			   mpart_LastColumn) == mpart_LastColumn) || // can tell for sure if incremental is going to be used
			   (m_t_marginals_config.n_relinearize_policy & mpart_LastColumn) == mpart_LastColumn) { // never know if we fallback to batch, though
				CLastElementOrderingConstraint leoc;
				mord.p_BlockOrdering(m_lambda, leoc.p_Get(m_lambda.n_BlockColumn_Num()),
					m_lambda.n_BlockColumn_Num(), true); // constrain the last column to be the last column (a quick fix) // todo - handle this properly, will be unable to constrain like this in fast R (well, actually ...) // todo - see what is this doing to the speed
			} else
				mord.p_BlockOrdering(m_lambda, true); // unconstrained; the last column may be anywhere (could reuse R from the linear solver here - relevant also in batch (e.g. on venice))
			const size_t *p_order = mord.p_Get_InverseOrdering();
			CUberBlockMatrix lambda_perm;
			m_lambda.Permute_UpperTriangular_To(lambda_perm, p_order, mord.n_Ordering_Size(), true);
			if(!R.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(lambda_perm))
				throw std::runtime_error("fatal error: R.CholeskyOf_FBS<_TyLambdaMatrixBlockSizes>(lambda_perm) failed");
			// note that now the marginals are calculated with ordering: need to count with that, otherwise those are useless!

			// todo - reuse what the linear solver calculated, if we have it (not if schur, )
			// todo - think of what happens if using schur ... have to accelerate the dense margs differently
			//		probably the whole mindset of having R is wrong, it would be much better to leave
			//		it up to the linear solver to solve for the columns

			/*printf("debug: matrix size: " PRIsize " / " PRIsize " (" PRIsize " nnz)\n",
				lambda_perm.n_BlockColumn_Num(), lambda_perm.n_Column_Num(), lambda_perm.n_NonZero_Num());
			float f_avg_block_size = float(lambda_perm.n_Column_Num()) / lambda_perm.n_BlockColumn_Num();
			printf("debug: diagonal nnz: " PRIsize "\n", size_t(lambda_perm.n_BlockColumn_Num() *
				(f_avg_block_size * f_avg_block_size)));
			printf("debug: factor size: " PRIsize " / " PRIsize " (" PRIsize " nnz)\n",
				R.n_BlockColumn_Num(), R.n_Column_Num(), R.n_NonZero_Num());*/
			// see how much we compute, compared to g2o

			timer.Accum_DiffSample(m_f_extra_chol_time);

			size_t n_add_edge_num = m_r_system.r_Edge_Pool().n_Size() - m_marginals.n_Edge_Num();
			bool b_incremental = m_marginals.b_CanUpdate() && CMarginals::b_PreferIncremental(m_r_system,
				m_marginals.r_SparseMatrix(), m_lambda, R, mord, m_marginals.n_Edge_Num(),
				m_t_marginals_config.n_incremental_policy);
//#ifdef __SE_TYPES_SUPPORT_L_SOLVERS
			if(b_incremental) { // otherwise just update what we have
				CUberBlockMatrix &r_m = const_cast<CUberBlockMatrix&>(m_marginals.r_SparseMatrix()); // watch out, need to call Swap_SparseMatrix() afterwards
				if(!CMarginals::Update_BlockDiagonalMarginals_FBS<false>(m_r_system, r_m, m_lambda,
				   R, mord, m_marginals.n_Edge_Num(), m_t_marginals_config.n_incremental_policy)) {
#ifdef _DEBUG
					fprintf(stderr, "warning: Update_BlockDiagonalMarginals_FBS() had a numerical issue:"
						" restarting with Calculate_DenseMarginals_Recurrent_FBS() instead\n");
#endif // _DEBUG
					b_incremental = false;
					// failed, will enter the batch branch below, that will not have a numerical issue
				} else {
					m_marginals.Swap_SparseMatrix(r_m); // now the marginals know that the matrix changed

					timer.Accum_DiffSample(m_f_incmarginals_time);
					++ m_n_incmarginals_num;
				}
			}
/*#else // __SE_TYPES_SUPPORT_L_SOLVERS
#pragma message("warning: the fast incremental marginals not available: __SE_TYPES_SUPPORT_L_SOLVERS not defined")
			b_incremental = false;
#endif // __SE_TYPES_SUPPORT_L_SOLVERS*/
			if(!b_incremental) { // if need batch marginals
				CUberBlockMatrix margs_ordered;
				CMarginals::Calculate_DenseMarginals_Recurrent_FBS<_TyLambdaMatrixBlockSizes>(margs_ordered, R,
					mord, m_t_marginals_config.n_relinearize_policy, false); // todo - use this in other solvers as well
				// calculate the thing

				CUberBlockMatrix &r_m = const_cast<CUberBlockMatrix&>(m_marginals.r_SparseMatrix()); // watch out, need to call Swap_SparseMatrix() afterwards
				margs_ordered.Permute_UpperTriangular_To(r_m, mord.p_Get_Ordering(),
					mord.n_Ordering_Size(), false); // no share! the original will be deleted
				m_marginals.Swap_SparseMatrix(r_m); // now the marginals know that the matrix changed
				// take care of having the correct permutation there

				m_marginals.EnableUpdate();
				// now the marginals are current and can be updated until the linearization point is changed again

				timer.Accum_DiffSample(m_f_marginals_time);
			}

			m_marginals.Set_Edge_Num(m_r_system.r_Edge_Pool().n_Size());
			// now all those edges are in the marginals

#if 0
			try {
				CUberBlockMatrix margs_ref, margs_untg;
				CMarginals::Calculate_DenseMarginals_Recurrent_FBS<_TyLambdaMatrixBlockSizes>(margs_ref, R);
				margs_ref.Permute_UpperTriangular_To(margs_untg, mord.p_Get_Ordering(),
					mord.n_Ordering_Size(), true); // ref is ok, it will be a short lived matrix
				Eigen::MatrixXd margs_buffer;
				margs_untg.Convert_to_Dense(margs_buffer);
				Eigen::VectorXd diag_rec = margs_buffer.diagonal();
				// recursive does not hurt

				CUberBlockMatrix R_unord;
				R_unord.CholeskyOf(m_lambda); // sloow
				CMarginals::Calculate_DenseMarginals_Ref(margs_buffer, R_unord); // sloow
				//CMarginals::Calculate_DenseMarginals_Fast/*_Parallel_FBS<_TyLambdaMatrixBlockSizes>*/(margs_buffer, R,
				//	mord.p_Get_Ordering(), mord.n_Ordering_Size()); // why does this not work!?
				double f_full_norm = margs_buffer.norm();
				double f_full_max = margs_buffer.maxCoeff();
				double f_diag_norm = margs_buffer.diagonal().norm();
				double f_diag_max = margs_buffer.diagonal().maxCoeff();

				Eigen::VectorXd diag_ref = margs_buffer.diagonal();
				//if(!margs_untg.b_EqualLayout(m_marginals.r_SparseMatrix()))
				//	printf("error: the marginals have different block layout\n");
				double f_all_err = CMarginals::f_IncompleteDifference(margs_buffer, m_marginals.r_SparseMatrix());
				m_marginals.r_SparseMatrix().Convert_to_Dense(margs_buffer);
				double f_error = (diag_ref - margs_buffer.diagonal()).norm();
				printf("debug: vertex " PRIsize ": added " PRIsize
					" edges: marginals diagonal tracking: %g (weighted: %g, fulltrack: %g, mnorm: %g, mmax: %g, dnorm: %g, dmax: %g, recerr: %g)\n",
					m_lambda.n_BlockColumn_Num(), n_add_edge_num, f_error, f_error / f_diag_norm, f_all_err,
					f_full_norm, f_full_max, f_diag_norm, f_diag_max, (diag_rec - diag_ref).norm());
				if(f_error > 100 /*|| diag_ref.maxCoeff() < 250*/) { // the 250 thing is only good for debugging intel.txt
					printf("\tthis calls for special attention ... this is called as (%d, %g)\n",
						int(n_max_iteration_num), f_min_dx_norm);
					Eigen::VectorXd diag_ref2 = margs_buffer.diagonal(); // now there is m_marginals.r_SparseMatrix() in there
					CUberBlockMatrix R_unord;
					R_unord.CholeskyOf(m_lambda); // sloow
					CMarginals::Calculate_DenseMarginals_Ref(margs_buffer, R_unord); // sloow
					double f_error = (diag_ref2 - margs_buffer.diagonal()).norm();
					printf("debug again: vertex " PRIsize ": added " PRIsize
						" edges: marginals diagonal tracking: %g (max elem: %g, difference of refs: %g)\n",
						m_lambda.n_BlockColumn_Num(), n_add_edge_num, f_error,
						margs_buffer.diagonal().maxCoeff(), (margs_buffer.diagonal() - diag_ref).norm());

					m_marginals.DisableUpdate(); // see how long this lasts
				}
			} catch(std::bad_alloc&) {
				fprintf(stderr, "warning: unable to verify marginals (bad_alloc)\n");
			}
			// the new detailed comparison code

			// the old comparison code is below
#if 0
			try {
				CUberBlockMatrix margs_ref, margs_untg;
				CMarginals::Calculate_DenseMarginals_Recurrent_FBS<_TyLambdaMatrixBlockSizes>(margs_ref, R);
				margs_ref.Permute_UpperTriangular_To(margs_untg, mord.p_Get_Ordering(),
					mord.n_Ordering_Size(), true); // ref is ok, it will be a short lived matrix
				Eigen::MatrixXd margs_buffer;
				margs_untg.Convert_to_Dense(margs_buffer);
				Eigen::VectorXd diag_ref = margs_buffer.diagonal();
				if(!margs_untg.b_EqualLayout(m_marginals.r_SparseMatrix()))
					printf("error: the marginals have different block layout\n");
				m_marginals.r_SparseMatrix().Convert_to_Dense(margs_buffer);
				double f_error = (diag_ref - margs_buffer.diagonal()).norm();
				printf("debug: vertex " PRIsize ": added " PRIsize
					" edges: marginals diagonal tracking: %g (%g weighted, max elem: %g)\n",
					m_lambda.n_BlockColumn_Num(), n_add_edge_num, f_error, f_error / diag_ref.norm(), diag_ref.maxCoeff());
				if(f_error > 100 /*|| diag_ref.maxCoeff() < 250*/) { // the 250 thing is only good for debugging intel.txt
					printf("\tthis calls for special attention ... this is called as (%d, %g)\n",
						int(n_max_iteration_num), f_min_dx_norm);
					Eigen::VectorXd diag_ref2 = margs_buffer.diagonal(); // now there is m_marginals.r_SparseMatrix() in there
					CUberBlockMatrix R_unord;
					R_unord.CholeskyOf(m_lambda); // sloow
					CMarginals::Calculate_DenseMarginals_Ref(margs_buffer, R_unord); // sloow
					double f_error = (diag_ref2 - margs_buffer.diagonal()).norm();
					printf("debug again: vertex " PRIsize ": added " PRIsize
						" edges: marginals diagonal tracking: %g (max elem: %g, difference of refs: %g)\n",
						m_lambda.n_BlockColumn_Num(), n_add_edge_num, f_error,
						margs_buffer.diagonal().maxCoeff(), (margs_buffer.diagonal() - diag_ref).norm());

					m_marginals.DisableUpdate(); // see how long this lasts
				}
			} catch(std::bad_alloc&) {
				fprintf(stderr, "warning: unable to verify marginals (bad_alloc)\n");
			}
			// calculate marginals again and subtract the diagonal to see how much off it gets
			// (i'm lazy, the whole diagonal blocks should be compared instead, that's why the dense bad_alloc business)
#endif // 0
#endif // 0

			/*FILE *p_fw;
			if((p_fw = fopen("marginals.txt", "w"))) {
				for(size_t i = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
					size_t n_order = m_lambda.n_BlockColumn_Base(i);
					size_t n_dimension = m_lambda.n_BlockColumn_Column_Num(i);
					// get col

					Eigen::MatrixXd block = m_marginals.r_Matrix().block(n_order, n_order, n_dimension, n_dimension);
					// get block

					fprintf(p_fw, "block_%d_%d = ", int(i), int(i));
					CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, block);
					// prints the matrix
				}
				fclose(p_fw);
			}*/
			// dump diagonal blocks of the marginals to a file
		}
		// now R is up to date, can get marginals
	}

protected:
	/**
	 *	@brief called from Incremental_Step() before optimizing, used to generate data for various demos
	 */
	void IncrementalPreOptDemo() const
	{
#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
		n_new_vertex_num = m_n_nonlinear_solve_threshold;
		m_b_had_loop_closure = true;
		// make it trigger
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

		/*{
			char p_s_name[256];
			sprintf(p_s_name, "edgeinfo_%03" _PRIsize ".m", m_r_system.r_Edge_Pool().n_Size());
			FILE *p_fw;
			if((p_fw = fopen(p_s_name, "w"))) {
				int n_edge_num = int(m_r_system.r_Edge_Pool().n_Size());
				fprintf(p_fw, "num_edges = %d;\n", n_edge_num);
				for(int i = 0; i < n_edge_num; ++ i) {
					const CEdgePose2D &e = *(const CEdgePose2D*)&m_r_system.r_Edge_Pool()[i];
					Eigen::Matrix3d t_jacobian0, t_jacobian1;
					Eigen::Vector3d v_expectation, v_error;
					e.Calculate_Jacobians_Expectation_Error(t_jacobian0, t_jacobian1, v_expectation, v_error);
					fprintf(p_fw, "edge_%d = [%d %d]; %% (zero-based) indices of vertices\n", i, e.n_Vertex_Id(0), e.n_Vertex_Id(1));
					fprintf(p_fw, "sigma_inv_%d = ", i);
					CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, e.t_Sigma_Inv());
					fprintf(p_fw, "jacobian_%d_0 = ", i);
					CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, t_jacobian0);
					fprintf(p_fw, "jacobian_%d_1 = ", i);
					CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, t_jacobian1);
				}

				fclose(p_fw);
			}
		}*/
		// save also image of the system as matlab
	}

	/**
	 *	@brief called from Optimize(), used to generate data for various demos
	 */
	void Dump_DemoData()
	{
		//Save_SystemMatrix_MM("lambda_system.mtx");
		//Dump_SystemMatrix("lambda_system.tga");

		//char p_s_filename[256];
		//sprintf(p_s_filename, "rss2013/%05d_0_lambda.tga", m_n_edges_in_lambda);
		//m_lambda.Rasterize_Symmetric(p_s_filename); // this lambda

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
		char p_s_filename[256];

		sprintf(p_s_filename, "rss2013/%05d_0_lambda.tga", m_n_edges_in_lambda);
		m_lambda.Rasterize_Symmetric(p_s_filename); // this lambda
		sprintf(p_s_filename, "rss2013/%05d_1_lambda2.tga", m_n_edges_in_lambda);
		m_lambda.Rasterize_Symmetric(lambda_prev, true, p_s_filename); // with changes marked, not symmetric

		CUberBlockMatrix Lprev;
		Lprev.CholeskyOf(lambda_prev);
		// do it now, will damage lambda_prev

		typename _TyEdgeMultiPool::_TyConstBaseRef r_last_edge =
			m_r_system.r_Edge_Pool()[m_r_system.r_Edge_Pool().n_Size() - 1];
		size_t v0 = r_last_edge.n_Vertex_Id(0);
		size_t v1 = r_last_edge.n_Vertex_Id(1);
		if(lambda_prev.n_BlockColumn_Num()) {
			size_t v = std::min(v0, v1);
			lambda_prev.t_GetBlock_Log(v, v).setZero(); // make sure this is clear
		}

		sprintf(p_s_filename, "rss2013/%05d_2_lambda3.tga", m_n_edges_in_lambda);
		m_lambda.CopyTo(lambda_prev);
		lambda_prev.Scale(0); // zero out the whole thing
		m_lambda.Rasterize_Symmetric(lambda_prev, true, p_s_filename); // all changed, not symmetric

		CUberBlockMatrix L;
		CUberBlockMatrix Lem;
		L.CopyLayoutTo(Lem);
		L.CholeskyOf(m_lambda);
		sprintf(p_s_filename, "rss2013/%05d_3_Lnoord.tga", m_n_edges_in_lambda);
		L.Rasterize(p_s_filename); // no fill-in highlighting
		sprintf(p_s_filename, "rss2013/%05d_4_Lnoord_fill-in.tga", m_n_edges_in_lambda);
		L.Rasterize(m_lambda, false, p_s_filename); // highlight fill-in
		sprintf(p_s_filename, "rss2013/%05d_a_Lnoord_red.tga", m_n_edges_in_lambda);
		L.Rasterize(Lem, false, p_s_filename); // highlight fill-in
		sprintf(p_s_filename, "rss2013/%05d_b_Lnoord_inc.tga", m_n_edges_in_lambda);
		L.Rasterize(Lprev, true, p_s_filename); // highlight fill-in

		m_lambda.CopyTo(lambda_prev);
		// copy the lambda matrix

		std::vector<size_t> cumsum, ccumsum;
		for(size_t i = 0, n_csum = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			n_csum += m_lambda.n_BlockColumn_Column_Num(i);
			cumsum.push_back(n_csum);
		}
		size_t p_cumsum[1] = {1};
		ccumsum.insert(ccumsum.begin(), p_cumsum, p_cumsum + 1);
		CUberBlockMatrix rhs(cumsum.begin(), cumsum.end(), ccumsum.begin(), ccumsum.end());
		m_v_dx.resize(n_variables_size, 1);
		Collect_RightHandSide_Vector(m_v_dx);
		for(size_t i = 0, n_csum = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_w = m_lambda.n_BlockColumn_Column_Num(i);
			if(!rhs.Append_Block(m_v_dx.segment(n_csum, n_w), n_csum, 0))
				throw std::runtime_error("segment size fail");
			n_csum += n_w;
		}
		sprintf(p_s_filename, "rss2013/%05d_5_rhs.tga", m_n_edges_in_lambda);
		rhs.Rasterize(p_s_filename);

		CUberBlockMatrix rhs_nonzero, rhs_hl;
		rhs.CopyLayoutTo(rhs_nonzero);
		rhs.CopyLayoutTo(rhs_hl);
		m_v_dx.setConstant(1);
		for(size_t i = 0, n_csum = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_w = m_lambda.n_BlockColumn_Column_Num(i);
			if(!(i == v0 || i == v1)) {
				if(!rhs_hl.Append_Block(m_v_dx.segment(n_csum, n_w), n_csum, 0))
					throw std::runtime_error("segment size fail");
			}
			n_csum += n_w;
		}
		for(size_t i = 0, n_csum = 0, n = m_lambda.n_BlockColumn_Num(); i < n; ++ i) {
			size_t n_w = m_lambda.n_BlockColumn_Column_Num(i);
			if(!rhs_nonzero.Append_Block(m_v_dx.segment(n_csum, n_w), n_csum, 0))
				throw std::runtime_error("segment size fail");
			n_csum += n_w;
		}
		sprintf(p_s_filename, "rss2013/%05d_9_rhs_add.tga", m_n_edges_in_lambda);
		rhs_nonzero.Rasterize(rhs_hl, true, p_s_filename);
		rhs.Scale(0);
		sprintf(p_s_filename, "rss2013/%05d_6_rhs_red.tga", m_n_edges_in_lambda);
		rhs_nonzero.Rasterize(rhs, true, p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_7_rhs_nnz.tga", m_n_edges_in_lambda);
		rhs_nonzero.Rasterize(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_8_rhs_z.tga", m_n_edges_in_lambda);
		rhs.Rasterize(p_s_filename);

#if 0
		/*sprintf(p_s_filename, "rss2013/%05d_0_lambda.tga", m_n_verts_in_lambda);
		m_lambda.Rasterize_Symmetric(p_s_filename);*/ // probably don't need lambdas

		CUberBlockMatrix lambda_perm;

		CMatrixOrdering mord;
		mord.p_BlockOrdering(m_lambda, true);
		const size_t *p_perm = mord.p_Get_InverseOrdering();

		if(mord100.n_Ordering_Size() < m_lambda.n_BlockColumn_Num()) {
			if(!(m_lambda.n_BlockColumn_Num() % 100))
				mord100.p_BlockOrdering(m_lambda, true); // do full every 100
			else {
				mord100.p_InvertOrdering(mord100.p_ExtendBlockOrdering_with_Identity(
					m_lambda.n_BlockColumn_Num()), m_lambda.n_BlockColumn_Num());
			}
		}
		// maintain an "every 100" ordering as well

		m_lambda.Permute_UpperTriangular_To(lambda_perm, p_perm, mord.n_Ordering_Size(), true);
		sprintf(p_s_filename, "rss2013/%05d_1_lambda-perm.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
		//	lambda_perm.Rasterize_Symmetric(p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 3); // do not really need lambdas right now

		L.CholeskyOf(lambda_perm);
		size_t n_L_good_blocks = L.n_Block_Num();
		sprintf(p_s_filename, "rss2013/%05d_2_L.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			L.Rasterize(lambda_perm, false, p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 3 : 2); // highlight fill-in

		L.CholeskyOf(m_lambda);
		sprintf(p_s_filename, "rss2013/%05d_9_Lnoord.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			L.Rasterize(lambda_perm, false, p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 2); // highlight fill-in

		const size_t *p_perm100 = mord100.p_Get_InverseOrdering();
		m_lambda.Permute_UpperTriangular_To(lambda_perm, p_perm100, mord100.n_Ordering_Size(), true);
		sprintf(p_s_filename, "rss2013/%05d_4_lambda-perm-100.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
		//	lambda_perm.Rasterize_Symmetric(p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 5 : 3); // do not really need lambdas right now

		L.CholeskyOf(lambda_perm);
		sprintf(p_s_filename, "rss2013/%05d_5_L100.tga", m_n_verts_in_lambda);
		//if(m_n_verts_in_lambda > size_t(n_dummy_param)) // continue from before
			L.Rasterize(lambda_perm, false, p_s_filename, (m_n_verts_in_lambda < 750 * 6 / _TyAMatrixBlockSizes::_TyHead::ColsAtCompileTime)? 3 : 2); // highlight fill-in

		sprintf(p_s_filename, "rss2013/%05d_3_stats.txt", m_n_verts_in_lambda);
		FILE *p_fw;
		if((p_fw = fopen(p_s_filename, "w"))) {
			fprintf(p_fw, PRIsize "\n", lambda_perm.n_Block_Num()); // save density of lambda
			fprintf(p_fw, PRIsize "\n", n_L_good_blocks); // save density of L
			fprintf(p_fw, PRIsize "\n", L.n_Block_Num()); // save density of L at every100
			fclose(p_fw);
		}
#endif // 0
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

#ifdef __ICRA_GET_PRESENTATION_ANIMATION_DATA
		char p_s_filename[256];
		sprintf(p_s_filename, "lambda_%03d.tga", m_n_real_step);
		m_lambda.Rasterize(p_s_filename);
		// dump lambda

		sprintf(p_s_filename, "lambda_%03d.txt", m_n_real_step);
		cs *p_lam = m_lambda.p_BlockStructure_to_Sparse();
		FILE *p_fw = fopen(p_s_filename, "w");
		fprintf(p_fw, "%d x %d\n", p_lam->m, p_lam->n);

		{
			cs *m = p_lam;
			int cols = m->n;
			int rows = m->m;
			int size_i = m->p[cols];
			int col = 0;
			for(int i = 0; i < size_i; ++ i) { // 8x if, 8x add
				int row = m->i[i]; // 7x read
				while(m->p[col + 1] <= i) // 7x if / 10x if (in case of while), the same amount of reads
					col ++; // 4x add
				fprintf(p_fw, "%d, %d\n", row, col);
			}
			// 14 - 17x read, 12x add, 15 - 18x if
		}
		fclose(p_fw);
		sprintf(p_s_filename, "lambda_%03d_ref.txt", m_n_real_step);
		p_fw = fopen(p_s_filename, "w");
		fprintf(p_fw, "%d x %d\n", p_lam->m, p_lam->n);
		{
			cs *m = p_lam;
			int cols = m->n;
			int rows = m->m;
			for(int j = 0; j < cols; ++ j) { // 5x if, 5x add
				int col = j;
				int first_i = m->p[j]; // 4x read
				int last_i = m->p[j + 1]; // 4x read
				for(int i = first_i; i < last_i; ++ i) { // 7+4=11x if, 11x add
					int row = m->i[i]; // 7x read
					fprintf(p_fw, "%d, %d\n", row, col);
				}
			}
			// 15x read, 16x add, 16x if
		}
		fclose(p_fw);
		cs_spfree(p_lam);
#endif // __ICRA_GET_PRESENTATION_ANIMATION_DATA

		/*m_lambda.Rasterize_Symmetric("lambda.tga");
		{
			FILE *p_fw = fopen("system.dot", "w");
			fprintf(p_fw, "digraph system_in_A {\n");
			printf("\trankdir = LR;\n");
			fprintf(p_fw, "\tnode [shape = doublecircle];");
			for(size_t i = 0, n = m_r_system.r_Vertex_Pool().n_Size(); i < n; ++ i) {
				const _TyBaseVertex &r_vertex = m_r_system.r_Vertex_Pool()[i];
				if(r_vertex.n_Dimension() == 2) {
					char p_s_vname[64];
					sprintf(p_s_vname, "L" PRIsize, i);
					fprintf(p_fw, " %s", p_s_vname);
				}
			}
			printf("\n");
			fprintf(p_fw, "\tnode [shape = circle];\n");
			fprintf(p_fw, "\t%s -> %s;\n", "UF", "V0");
			for(size_t i = 0, n = m_r_system.r_Edge_Pool().n_Size(); i < n; ++ i) {
				const _TyBaseEdge &r_edge = m_r_system.r_Edge_Pool()[i];
				int n_v0 = r_edge.n_Vertex_Id(0);
				int n_v1 = r_edge.n_Vertex_Id(1);
				char p_s_vname0[64], p_s_vname1[64];
				sprintf(p_s_vname0, "%c" PRIsize, (m_r_system.r_Vertex_Pool()[n_v0].n_Dimension() == 2)? 'L' : 'V', n_v0);
				sprintf(p_s_vname1, "%c" PRIsize, (m_r_system.r_Vertex_Pool()[n_v1].n_Dimension() == 2)? 'L' : 'V', n_v1);
				fprintf(p_fw, "\t%s -> %s;\n", p_s_vname0, p_s_vname1);
			}
			fprintf(p_fw, "}\n");
			fclose(p_fw);
		}*/
		// save a matrix and make a .dot file of the graph
	}

	/**
	 *	@brief called from Incremental_Step() after optimizing, used to generate data for various demos
	 */
	void IncrementalPostOptDemo(size_t n_new_vertex_num) const
	{
#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
		if(n_new_vertex_num) {
			FILE *p_fw = fopen("chi2perVert.txt", (m_n_real_step > 0)? "a" : "w");
			fprintf(p_fw, "%f\n", f_Chi_Squared_Error_Denorm());
			fclose(p_fw);
			// dump chi2

            //printf("time is %lf\n", m_timer.f_Time()); // print incremental times (ICRA 2013)
#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_ICRA2013_ANIMATION_DATA
			char p_s_filename[256];
			sprintf(p_s_filename, "icra2013/sys%05" _PRIsize ".txt", m_r_system.r_Vertex_Pool().n_Size());
			m_r_system.Dump(p_s_filename);
			//sprintf(p_s_filename, "icra2013/sys%05" _PRIsize ".tga", m_r_system.r_Vertex_Pool().n_Size());
			//m_r_system.Plot2D(p_s_filename); // no plot, takes a ton of memory
			// dump vertices for the ICRA animation
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_ICRA2013_ANIMATION_DATA
		}
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2

		/*{
			char p_s_name[256], p_s_name1[256];
			sprintf(p_s_name, "lambda_%03" _PRIsize ".tga", m_r_system.r_Edge_Pool().n_Size());
			m_lambda.Rasterize(p_s_name);
			sprintf(p_s_name, "lambda_%03" _PRIsize ".mtx", m_r_system.r_Edge_Pool().n_Size());
			sprintf(p_s_name1, "lambda_%03" _PRIsize ".bla", m_r_system.r_Edge_Pool().n_Size());
			m_lambda.Save_MatrixMarket(p_s_name, p_s_name1);
		}*/
		// save matrix image and matrix market
	}

	CNonlinearSolver_Lambda(const CNonlinearSolver_Lambda &UNUSED(r_solver)); /**< @brief the object is not copyable */
	CNonlinearSolver_Lambda &operator =(const CNonlinearSolver_Lambda &UNUSED(r_solver)) { return *this; } /**< @brief the object is not copyable */
};

#endif // !__NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
