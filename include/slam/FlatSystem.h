/*
								+----------------------------------+
								|                                  |
								| ***  Pool-based flat system  *** |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2012  |
								|                                  |
								|           FlatSystem.h           |
								|                                  |
								+----------------------------------+
*/

#pragma once
#ifndef __FLAT_SYSTEM_TEMPLATE_INCLUDED
#define __FLAT_SYSTEM_TEMPLATE_INCLUDED

/**
 *	@file include/slam/FlatSystem.h
 *	@brief pool-based flat system template
 *	@author -tHE SWINe-
 *	@date 2012-09-03
 *
 *	This is the new take on implementation of reusable system, this time with templates.
 *	This is very simple to use, one merely needs to list all the types handled by the
 *	system in advance, specialize the system template, and then almost all of the calls
 *	to edge and vertex member functions can be statically linked (not virtual). It is
 *	also very simple to modify SE(2) solver to SE(3) or anything else, really.
 *
 *	t_odo - see if we can implement template token list and have:
 *		* vertex_pose_token, vertex_landmark_token,
 *		* edge_pose_pose_token, edge_landmark_landmark_token
 *	and to create lists of pools from that
 *	and to create vertex / edge addition functions featuring those tokens
 *	then there is just a vector of edge pointers and that shields everything out of that
 *		* yes, this was indeed implemented
 *
 *	t_odo - figure out what else the edge needs to implement and add it, implement ASLAM solver
 *		and test it with the new structure
 *	@todo - describe the "virtual call evasion" technique used here
 *	t_odo - figure out what operations can be possibly offloaded to solver to simplify the edge
 *			(and vertex) classes // created base edge and vertex implementation templates
 *	t_odo - handle the unary factor bussiness nicely
 *	t_odo - abstract the linear solvers, test cholmod performance
 *
 *	@date 2012-09-12
 *
 *	Fixed linux build issues, namely missing template argument in CMultiPool specialization
 *	for a single type, added namespace for types contained in CFlatSystem, and some other minor bugs.
 *
 *	@date 2013-05-23
 *
 *	Added minimal parallel threshold in For_Each_Parallel(), maybe that was not originally intended
 *	(the caller was supposed to decide whether they want parallel), but it ended up used in a way
 *	that requires this.
 *
 *	@date 2013-11-01
 *
 *	Added support for fixed vertices (they must have b_IsFixed() function, and the vertices must
 *	be determined on initialization, changing to fixed / not fixed later on is not supported at the
 *	moment).
 *
 */

#include "eigen/Eigen/Dense"
#include "slam/TypeList.h"
#include "slam/Segregated.h"
#include "slam/BlockMatrix.h"

/**
 *	@def __BASE_TYPES_USE_ID_ADDRESSING
 *	@brief if defined, the faster matrix addressing using vertex and edge id's is used (blockwise)
 *		otherwise addressing based on vertex and edge orders (elementwise) is used
 *	@note This must be defined here as the CFlatSystem depends on it.
 */
//#define __BASE_TYPES_USE_ID_ADDRESSING

/*#if !defined(_WIN32) && !defined(_WIN64)
#define __USE_ALIGNED_MULTIPOOL
#endif // !_WIN32 && !_WIN64*/ // aligned SE(2) types do not work, don't align, it's slower
// decide whether to align pool data // todo - document it and add manual overrides

/**
 *	@brief very basic implementation of unary factor initialization that uses unit UF
 */
class CBasicUnaryFactorFactory {
public:
	/**
	 *	@brief initializes unary factor, based on the first edge introduced to the system
	 *
	 *	@param[out] r_t_unary_factor is the unary factor matrix
	 *	@param[out] r_v_unary_error is the error vector associated with the first vertex
	 *	@param[in] r_edge is the first edge in the system (value not used)
	 */
	template <class CEdge>
	inline void operator ()(Eigen::MatrixXd &r_t_unary_factor,
		Eigen::VectorXd &r_v_unary_error, const CEdge &UNUSED(r_edge))
	{
		size_t n_edge_dimension = CEdge::n_vertex0_dimension;
		r_t_unary_factor.resize(n_edge_dimension, n_edge_dimension); // unary factor is a unit matrix
		r_t_unary_factor.setIdentity();
		r_v_unary_error = Eigen::VectorXd(n_edge_dimension);
		r_v_unary_error.setZero(); // no error on the first vertex
	}
};

/**
 *	@brief implementation of unary factor initialization that takes cholesky
 *		of the information matrix of the first edge as unary factor
 */
class CProportionalUnaryFactorFactory {
public:
	/**
	 *	@brief initializes unary factor, based on the first edge introduced to the system
	 *
	 *	@param[out] r_t_unary_factor is the unary factor matrix
	 *	@param[out] r_v_unary_error is the error vector associated with the first vertex
	 *	@param[in] r_edge is the first edge in the system
	 */
	template <class CEdge>
	inline void operator ()(Eigen::MatrixXd &r_t_unary_factor,
		Eigen::VectorXd &r_v_unary_error, const CEdge &r_edge)
	{
		size_t n_edge_dimension = r_edge.n_Dimension();
		r_t_unary_factor = r_edge.t_Sigma_Inv().llt().matrixU(); // unary factor is cholesky of the information matrix
		r_v_unary_error = Eigen::VectorXd(n_edge_dimension);
		r_v_unary_error.setZero(); // no error on the first vertex
	}
};

/**
 *	@brief multiple data type (heterogenous), yet statically typed pools
 *	@todo - put this to multipool.h or something
 */
namespace __multipool {

/**
 *	@brief list of pools per every data type in the list
 *	@note While this works as a simple multipool, it isn't entirely practical without accessors.
 */
template <class CListType, const int n_pool_page_size, const int n_pool_memory_align>
class CPoolList {
public:
	typedef typename CListType::_TyHead TPayloadType; /**< @brief payload data type */
	typedef forward_allocated_pool<TPayloadType, n_pool_page_size, n_pool_memory_align> _TyFirstPool; /**< @brief pool of the first data type being stored */

protected:
	_TyFirstPool m_pool; /**< @brief contains vertex pool of a given type */
	CPoolList<typename CListType::_TyTail, n_pool_page_size, n_pool_memory_align> m_recurse; /**< @brief contains vertex pools of all the types following in the list */

public:
	/**
	 *	@brief adds a new element to the pool
	 *
	 *	@param[out] r_p_storage is pointer to the place where the element is stored
	 *	@param[in] t_payload is initialization value for the new element
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	template <class CGenericPalyoadType>
	inline void AddElement(CGenericPalyoadType *&r_p_storage, CGenericPalyoadType t_payload)
	{
		m_recurse.AddElement(r_p_storage, t_payload); // not this type of payload, perhaps recurse can contain it
		// note that if the compiler generates errors here, the type geing inserted is most likely not on the list
		// @todo - make static assertion here that would print a human-readable message as well
	}

	/**
	 *	@brief adds a new element to the pool (version specialized for the payload type of this pool)
	 *
	 *	@param[out] r_p_storage is pointer to the place where the element is stored
	 *	@param[in] t_payload is initialization value for the new element
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	inline void AddElement(TPayloadType *&r_p_storage, TPayloadType t_payload)
	{
		m_pool.push_back(t_payload); // our type of payload
		r_p_storage = &*(m_pool.end() - 1); // return pointer to the payload
	}

	/**
	 *	@brief gets access to the data pool, containing the first type in the list
	 *	@return Returns reference to the first data pool.
	 */
	inline _TyFirstPool &r_Pool()
	{
		return m_pool;
	}

	/**
	 *	@brief gets access to the data pool, containing the first type in the list
	 *	@return Returns const reference to the first data pool.
	 */
	inline const _TyFirstPool &r_Pool() const
	{
		return m_pool;
	}
};

/**
 *	@brief list of pools per every data type in the list (the list terminator specialization)
 */
template <const int n_pool_page_size, const int n_pool_memory_align>
class CPoolList<CTypelistEnd, n_pool_page_size, n_pool_memory_align> {};

/**
 *	@brief heterogenous (but static) pool container
 *
 *	@tparam CBaseType is base type of items being stored in the multipool
 *	@tparam CTypelist is list of types being stored in the multipool
 *	@tparam n_pool_page_size is pool page size (per each type being stored)
 */
template <class CBaseType, class CTypelist, const size_t n_pool_page_size = 1024>
class CMultiPool {
public:
	/**
	 *	@brief configuration, stored as enum
	 */
	enum {
		pool_PageSize = n_pool_page_size, /**< @brief size of pool page, in elements */
#ifdef __USE_ALIGNED_MULTIPOOL
		pool_MemoryAlign = 16 /**< @brief memory alignment of data stored in the pool */
#else // __USE_ALIGNED_MULTIPOOL
		pool_MemoryAlign = 0 /**< @brief memory alignment of data stored in the pool */ // or 16 for aligned types and SSE // t_odo
#endif // __USE_ALIGNED_MULTIPOOL
	};

	typedef CBaseType _TyBaseType; /**< @brief base type of items being stored in the multipool */
	typedef CTypelist _TyTypelist; /**< @brief list of types being stored in the multipool */

protected:

	std::vector<_TyBaseType*> m_uniform_list;  /**< @brief list of pointers with uniform data types (use virtual functions) */
	CPoolList<_TyTypelist, pool_PageSize, pool_MemoryAlign> m_pool_list;  /**< @brief list of all the pools for all the data types */

	/**
	 *	@brief function object that transparently inserts pointer dereferencing
	 *		between std::for_each and client function object
	 *	@tparam COp is client function object
	 */
	template <class COp>
	class CDereference {
	protected:
		COp m_op; /**< @brief client function object */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] op is instance of client function object
		 */
		inline CDereference(COp op)
			:m_op(op)
		{}

		/**
		 *	@brief function operator; dereferences the element
		 *		and passes it to the client function object
		 *	@param[in] p_elem is pointer to the element of the multipool
		 */
		inline void operator ()(_TyBaseType *p_elem)
		{
			m_op(*p_elem);
		}

		/**
		 *	@brief conversion to the client function object
		 *	@return Returns instance of client function object,
		 *		passed to the constructor.
		 */
		inline operator COp() const
		{
			return m_op;
		}
	};

	/*template <class COp>
	class CConstDereference {
	protected:
		COp m_op;

	public:
		inline CConstDereference(COp op)
			:m_op(op)
		{}

		inline void operator ()(const _TyBaseType *p_elem)
		{
			m_op(*p_elem);
		}

		inline operator COp() const
		{
			return m_op;
		}
	};*/

public:
	/**
	 *	@brief adds an element to the multipool (at the end)
	 *	@tparam CGenericPalyoadType is type of the element being inserted
	 *		(must be on the list, otherwise the call will result in compile errors)
	 *	@param[in] t_payload is the element being inserted
	 *	@return Returns reference to the inserted element
	 *		(the address is guaranteed not to change).
	 *	@note This function throws std::bad_alloc.
	 */
	template <class CGenericPalyoadType>
	CGenericPalyoadType &r_Add_Element(CGenericPalyoadType t_payload) // throw(std::bad_alloc)
	{
		CGenericPalyoadType *p_storage;
		m_pool_list.AddElement(p_storage, t_payload);
		m_uniform_list.push_back(static_cast<CBaseType*>(p_storage));
		return *p_storage;
	}

	/**
	 *	@brief determines whether the multipool is empty
	 *	@return Returns true in case the multipool is empty, otherwise returns false.
	 */
	inline bool b_Empty() const
	{
		return m_uniform_list.empty();
	}

	/**
	 *	@brief gets number of elements
	 *	@return Returns number of elements stored in the multipool.
	 */
	inline size_t n_Size() const
	{
		return m_uniform_list.size();
	}

	/**
	 *	@brief gets an element
	 *	@param[in] n_index is zero-based element index
	 *	@return Returns a reference to the selected element.
	 */
	inline _TyBaseType &operator [](size_t n_index)
	{
		return *m_uniform_list[n_index];
	}

	/**
	 *	@brief gets an element
	 *	@param[in] n_index is zero-based element index
	 *	@return Returns a const reference to the selected element.
	 */
	inline const _TyBaseType &operator [](size_t n_index) const
	{
		return *m_uniform_list[n_index];
	}

	/**
	 *	@brief iterates over all the elements in this multipool and performs operation on each
	 *	@tparam COp is client function object
	 *	@param[in] op is instance of the client function object
	 *	@return Returns instance of the client function object,
	 *		after performing all the operations (can e.g. perform reduction).
	 */
	template <class COp>
	COp For_Each(COp op)
	{
		return std::for_each(m_uniform_list.begin(),
			m_uniform_list.end(), CDereference<COp>(op));
	}

	/**
	 *	@brief iterates over all the elements in this multipool and performs operation on each
	 *
	 *	@tparam COp is client function object
	 *
	 *	@param[in] op is instance of the client function object
	 *	@param[in] n_parallel_thresh is threshold for parallelized processing
	 *
	 *	@note This performs the iteration in parallel, the function object must be reentrant.
	 *	@note This does not return the function object (avoids synchronization
	 *		and explicit reduction).
	 */
	template <class COp>
	void For_Each_Parallel(COp op, const int n_parallel_thresh = 50)
	{
#ifdef _OPENMP
		_ASSERTE(m_uniform_list.size() <= INT_MAX);
		const int n = int(m_uniform_list.size());
		#pragma omp parallel for default(shared) if(n >= n_parallel_thresh)
		for(int i = 0; i < n; ++ i)
			op(*m_uniform_list[i]);
#else // _OPENMP
		std::for_each(m_uniform_list.begin(), m_uniform_list.end(), CDereference<COp>(op));
#endif // _OPENMP
	}

	/**
	 *	@brief iterates over a selected range of elements
	 *		in this multipool and performs operation on each
	 *
	 *	@tparam COp is client function object
	 *
	 *	@param[in] n_first is zero-based index of the first element to be processed
	 *	@param[in] n_last is zero-based index of one after the last element to be processed
	 *	@param[in] op is instance of the client function object
	 *
	 *	@return Returns instance of the client function object,
	 *		after performing all the operations (can e.g. perform reduction).
	 */
	template <class COp>
	COp For_Each(size_t n_first, size_t n_last, COp op)
	{
		_ASSERTE(n_first <= n_last);
		_ASSERTE(n_last <= m_uniform_list.size());
		return std::for_each(m_uniform_list.begin() + n_first,
			m_uniform_list.begin() + n_last, CDereference<COp>(op));
	}

	/**
	 *	@brief iterates over a selected range of elements
	 *		in this multipool and performs operation on each
	 *
	 *	@tparam COp is client function object
	 *	@param[in] n_parallel_thresh is threshold for parallelized processing
	 *
	 *	@param[in] n_first is zero-based index of the first element to be processed
	 *	@param[in] n_last is zero-based index of one after the last element to be processed
	 *	@param[in] op is instance of the client function object
	 *
	 *	@note This performs the iteration in parallel, the function object must be reentrant.
	 *	@note This does not return the function object (avoids synchronization
	 *		and explicit reduction).
	 */
	template <class COp>
	void For_Each_Parallel(size_t n_first, size_t n_last, COp op, const int n_parallel_thresh = 50)
	{
		_ASSERTE(n_first <= n_last);
		_ASSERTE(n_last <= m_uniform_list.size());
#ifdef _OPENMP
		_ASSERTE(n_last <= INT_MAX);
		const int n = int(n_last);
		#pragma omp parallel for default(shared) if(n - int(n_first) >= n_parallel_thresh)
		for(int i = int(n_first); i < n; ++ i)
			op(*m_uniform_list[i]);
#else // _OPENMP
		std::for_each(m_uniform_list.begin() + n_first, m_uniform_list.begin() + n_last, CDereference<COp>(op));
#endif // _OPENMP
	}

	/*template <class COp>
	void For_Each(COp op) const
	{
		std::for_each(m_uniform_list.begin(), m_uniform_list.end(), CConstDereference<COp>(op));
	}

	template <class COp>
	void For_Each(size_t n_first, size_t n_last, COp op) const
	{
		_ASSERTE(n_first <= n_last);
		_ASSERTE(n_last <= m_uniform_list.size());
		std::for_each(m_uniform_list.begin() + n_first, m_uniform_list.begin() + n_last, CConstDereference<COp>(op));
	}*/
};

/**
 *	@brief heterogenous (but static) pool container (specialization
 *		for only a single type being stored)
 *
 *	@tparam CBaseType is base type of items being stored in the multipool
 *	@tparam CDerivedType is the only type of items being stored in the multipool
 *		(derived from CBaseType, might be equal to CBaseType)
 *	@tparam n_pool_page_size is pool page size (per each type being stored)
 *
 *	@todo Try storing pointers in std::vector for faster parallel for_each(),
 *		and try to write fap::for_each() and see what is faster.
 */
template <class CBaseType, class CDerivedType, const size_t n_pool_page_size>
class CMultiPool<CBaseType, CTypelist<CDerivedType, CTypelistEnd>, n_pool_page_size> {
public:
	/**
	 *	@brief configuration, stored as enum
	 */
	enum {
		pool_PageSize = n_pool_page_size, /**< @brief size of pool page, in elements */
#ifdef __USE_ALIGNED_MULTIPOOL
		pool_MemoryAlign = 16 /**< @brief memory alignment of data stored in the pool */
#else // __USE_ALIGNED_MULTIPOOL
		pool_MemoryAlign = 0 /**< @brief memory alignment of data stored in the pool */ // or 16 for aligned types and SSE // t_odo
#endif // __USE_ALIGNED_MULTIPOOL
	};

	typedef CBaseType _TyBaseType; /**< @brief base type of items being stored in the multipool */
	typedef CDerivedType _TyDerivedType; /**< @brief the only type of items being stored in the multipool */

protected:
	forward_allocated_pool<_TyDerivedType, pool_PageSize, pool_MemoryAlign> m_pool; /**< @brief pool for the single type being handled */
	// (static function binding, no space wasted for base type pointers)

public:
	/**
	 *	@brief adds an element to the multipool (at the end)
	 *	@param[in] t_payload is the element being inserted
	 *	@return Returns reference to the inserted element
	 *		(the address is guaranteed not to change).
	 *	@note This function throws std::bad_alloc.
	 */
	_TyDerivedType &r_Add_Element(_TyDerivedType t_payload) // throw(std::bad_alloc)
	{
		m_pool.push_back(t_payload);
		return *(m_pool.end() - 1);
	}

	/**
	 *	@copydoc CMultiPool::b_Empty()
	 */
	inline bool b_Empty() const
	{
		return m_pool.empty();
	}

	/**
	 *	@copydoc CMultiPool::n_Size()
	 */
	inline size_t n_Size() const
	{
		return m_pool.size();
	}

	/**
	 *	@copydoc CMultiPool::operator [](size_t)
	 */
	inline _TyBaseType &operator [](size_t n_index)
	{
		return static_cast<_TyBaseType&>(m_pool[n_index]);
	}

	/**
	 *	@copydoc CMultiPool::operator [](size_t) const
	 */
	inline const _TyBaseType &operator [](size_t n_index) const
	{
		return static_cast<const _TyBaseType&>(m_pool[n_index]);
	}

	/**
	 *	@copydoc CMultiPool::For_Each(COp)
	 */
	template <class COp>
	COp For_Each(COp op)
	{
		return std::for_each(m_pool.begin(), m_pool.end(), op);
	}

	/**
	 *	@copydoc CMultiPool::For_Each_Parallel(COp,const int)
	 */
	template <class COp>
	void For_Each_Parallel(COp op, const int n_parallel_thresh = 50)
	{
#ifdef _OPENMP
		_ASSERTE(m_pool.size() <= INT_MAX);
		const int n = int(m_pool.size());
		#pragma omp parallel for default(shared) if(n >= n_parallel_thresh)
		for(int i = 0; i < n; ++ i)
			op(m_pool[i]); // todo - write pool::parallel_for_each that would simplify it's pointer arithmetic
#else // _OPENMP
		std::for_each(m_pool.begin(), m_pool.end(), op);
#endif // _OPENMP
	}

	/**
	 *	@copydoc CMultiPool::For_Each(size_t,size_t,COp)
	 */
	template <class COp>
	COp For_Each(size_t n_first, size_t n_last, COp op)
	{
		_ASSERTE(n_first <= n_last);
		_ASSERTE(n_last <= m_pool.size());
		return std::for_each(m_pool.begin() + n_first, m_pool.begin() + n_last, op);
	}

	/**
	 *	@copydoc CMultiPool::For_Each_Parallel(size_t,size_t,COp,const int)
	 */
	template <class COp>
	void For_Each_Parallel(size_t n_first, size_t n_last, COp op, const int n_parallel_thresh = 50)
	{
		_ASSERTE(n_first <= n_last);
		_ASSERTE(n_last <= m_pool.size());
#ifdef _OPENMP
		_ASSERTE(n_last <= INT_MAX);
		const int n = int(n_last);
		#pragma omp parallel for default(shared) if(n - int(n_first) >= n_parallel_thresh)
		for(int i = int(n_first); i < n; ++ i)
			op(m_pool[i]); // todo - write pool::parallel_for_each that would simplify it's pointer arithmetic
#else // _OPENMP
		std::for_each(m_pool.begin() + n_first, m_pool.begin() + n_last, op);
#endif // _OPENMP
	}

	/*template <class COp>
	void For_Each(COp op) const
	{
		std::for_each(m_pool.begin(), m_pool.end(), op);
	}

	template <class COp>
	void For_Each(size_t n_first, size_t n_last, COp op) const
	{
		_ASSERTE(n_first <= n_last);
		_ASSERTE(n_last <= m_pool.size());
		std::for_each(m_pool.begin() + n_first, m_pool.begin() + n_last, op);
	}*/
};

} // ~__multipool

namespace plot_quality {

/**
 *	@brief plot quality profile names
 */
enum EQualityProfile {
	plot_Draft, /**< @brief fast plots for debugging */
	plot_Printing, /**< @brief nice plots for papers or web */
	plot_Printing_LandmarkTicksOnly, /**< @brief same as plot_Printing but only landmarks have ticks (not poses) */
	plot_Printing_NoTicks /**< @brief same as plot_Printing but there are no ticks for poses or landmarks */
};

} // ~plot_quality

/**
 *	@brief optimization system, customized for work with different primitives
 *
 *	@tparam CBaseVertex is base vertex type (all vertex types must be derived from it)
 *	@tparam CVertexTypelist is list of vertex types permitted in the system
 *	@tparam CBaseEdge is base edge type (all edge types must be derived from it)
 *	@tparam CEdgeTypelist is list of edge types permitted in the system
 *	@tparam CUnaryFactorFactory is class, responsible for the initialization of unary factor
 *	@tparam b_allow_fixed_vertices is fixed vertex enable flag (fixed vertices need special processing)
 *	@tparam n_pool_page_size is edge or vertex pool page size, in elements
 */
template <class CBaseVertex, class CVertexTypelist, /*template <class> class CVertexTypeTraits,*/ // unused
	class CBaseEdge, class CEdgeTypelist, class CUnaryFactorFactory = CBasicUnaryFactorFactory,
	const bool b_allow_fixed_vertices = false, const int n_pool_page_size = 1024>
class CFlatSystem {
protected:
	/**
	 *	@brief extracts the first jacobian matrix block type from a binary edge
	 *	@tparam CEdgeType is an edge type name
	 */
	template <class CEdgeType>
	class CEdgeTypeToFirstAJacobianType {
	public:
		typedef typename CEdgeType::_TyMatrix0 _TyResult; /**< @brief Eigen dense matrix block type */
	};

	/**
	 *	@brief extracts the second jacobian matrix block type from a binary edge
	 *	@tparam CEdgeType is an edge type name
	 */
	template <class CEdgeType>
	class CEdgeTypeToSecondAJacobianType {
	public:
		typedef typename CEdgeType::_TyMatrix1 _TyResult; /**< @brief Eigen dense matrix block type */
	};

	/**
	 *	@brief extracts the first unary factor matrix block type from a binary edge
	 *	@tparam CEdgeType is an edge type name
	 */
	template <class CEdgeType>
	class CEdgeTypeToFirstUnaryFactorType { // needed if the measurement dimension is different from either vertex dimensions
	public:
		typedef Eigen::Matrix<double, CEdgeType::n_vertex0_dimension,
			CEdgeType::n_vertex0_dimension> _TyResult; /**< @brief Eigen dense matrix block type */
	};

	/**
	 *	@brief extracts the second unary factor matrix block type from a binary edge
	 *	@tparam CEdgeType is an edge type name
	 */
	template <class CEdgeType>
	class CEdgeTypeToSecondUnaryFactorType { // needed if the measurement dimension is different from either vertex dimensions
	public:
		typedef Eigen::Matrix<double, CEdgeType::n_vertex1_dimension,
			CEdgeType::n_vertex1_dimension> _TyResult; /**< @brief Eigen dense matrix block type */
	};

	typedef typename CTransformTypelist<CEdgeTypelist,
		CEdgeTypeToFirstAJacobianType>::_TyResult TFirstAJacobianList; /**< @brief list of jacobian matrix block types per all the first vertices in all edges (some of the blocks in A) */
	typedef typename CTransformTypelist<CEdgeTypelist,
		CEdgeTypeToSecondAJacobianType>::_TyResult TSecondAJacobianList; /**< @brief list of jacobian matrix block types per all the second vertices in all edges (some of the blocks in A) */
	typedef typename CTransformTypelist<CEdgeTypelist,
		CEdgeTypeToFirstUnaryFactorType>::_TyResult TFirstAUnaryFactorList; /**< @brief list of unary factor matrix block types per all the first vertices in all edges (some of the blocks in A) */
	typedef typename CTransformTypelist<CEdgeTypelist,
		CEdgeTypeToSecondUnaryFactorType>::_TyResult TSecondAUnaryFactorList; /**< @brief list of unary factor matrix block types per all the second vertices in all edges (some of the blocks in A) */

	/**
	 *	@brief void pool type
	 *	@tparam _TyBase is base class of elements stored in this pool
	 *	@note This is implemented only to reduce the overhead of implementing all
	 *		the code that deals with the pools twice. Some of it can be unified like this.
	 */
	template <class _TyBase>
	struct TVoidPool {
		/**
		 *	@copydoc __multipool::CMultiPool::n_Size()
		 */
		inline size_t n_Size() const
		{
			return 0;
		}

		/**
		 *	@copydoc __multipool::CMultiPool::operator []()
		 */
		inline _TyBase &operator [](size_t UNUSED(n_index)) const
		{
			throw std::runtime_error("attempted element access on void pool");
			return *((_TyBase*)0);
		}
	};

public:
	typedef size_t _TyId; /**< @brief data type used as the index in the flat structures */

	/**
	 *	@brief configuration, stored as enum
	 */
	enum {
		pool_PageSize = n_pool_page_size, /**< @brief data storage page size */
		allow_FixedVertices = (b_allow_fixed_vertices)? 1 : 0 /**< @brief allows / disallows fixed vertices */
	};

	typedef CBaseVertex _TyBaseVertex; /**< @brief the data type for storing vertices */
	typedef CVertexTypelist _TyVertexTypelist; /**< @brief list of vertex types */
	typedef CBaseEdge _TyBaseEdge; /**< @brief the data type for storing measurements */
	typedef CEdgeTypelist _TyEdgeTypelist; /**< @brief list of edge types */

	typedef typename CUniqueTypelist<typename CConcatTypelist<typename
		CConcatTypelist<TFirstAJacobianList, TSecondAJacobianList>::_TyResult, typename
		CConcatTypelist<TFirstAUnaryFactorList, TSecondAUnaryFactorList>::_TyResult>::_TyResult>::_TyResult
		_TyJacobianMatrixBlockList; /**< @brief list of jacobian and UF matrix block types (blocks in A) */
	typedef typename __fbs_ut::CBlockSizesAfterPreMultiplyWithSelfTranspose<
		_TyJacobianMatrixBlockList>::_TyResult _TyHessianMatrixBlockList; /**< @brief list of hessian and UF matrix block types (blocks in lambda and L) */

	typedef __multipool::CMultiPool<_TyBaseVertex, _TyVertexTypelist, pool_PageSize> _TyVertexMultiPool; /**< @brief vertex multipool type */
	typedef typename CTypelistItemAt<CTypelist<TVoidPool<_TyBaseVertex>, CTypelist<_TyVertexMultiPool,
		CTypelistEnd> >, allow_FixedVertices>::_TyResult _TyFixVertexMultiPool; /**< @brief fixed vertex multipool type (same as vertex multipool if enabled, otherwise a void type) */
	typedef __multipool::CMultiPool<_TyBaseEdge, _TyEdgeTypelist, pool_PageSize> _TyEdgeMultiPool; /**< @brief edge multipool type */

protected:
	typedef typename CTypelistItemAt<CTypelist<TVoidPool<_TyBaseVertex*>, CTypelist<std::vector<_TyBaseVertex*>,
		CTypelistEnd> >, allow_FixedVertices>::_TyResult _TyVertexLookup; /**< @brief lookup table for the vertices (only used if fixed vertices are enabled) */

	/**
	 *	@brief code for vertex insertion, based on whether the fixed vertices are allowed
	 *	@tparam b_fixed_vertices_enabled is fixed vertex enable flag (they need special processing)
	 */
	template <class _CBaseVertex, class _CVertexTypelist, class _CBaseEdge,
		class _CEdgeTypelist, class _CUnaryFactorFactory,
		const bool b_fixed_vertices_enabled, const int _n_pool_page_size>
	class CGetVertexImpl {
	public:
		/**
		 *	@brief finds vertex by id, or create a new one in case there is no vertex with such id
		 *
		 *	@tparam CVertexType is the vertex type name
		 *	@tparam CInitializer is the lazy initializer type name
		 *
		 *	@param[in] n_id is the id of the vertex required (must be id of an existing vertex,
		 *		or one larger than id of the last vertex)
		 *	@param[in] init is the initializer functor; it is only called (converted to vertex)
		 *		in case new vertex is created
		 *	@param[in,out] r_vertex_pool is vertex pool for optimized vertices
		 *	@param[in] r_fix_vertex_pool is vertex pool for fixed vertices (unused)
		 *	@param[in] r_vertex_lookup is lookup table that contains both fixed and optimized vertices (unused)
		 *	@param[in,out] r_n_vertex_element_num is number of elements of all the
		 *		optimized vertices so far
		 *
		 *	@return Returns reference to the vertex, associated with id n_id.
		 *
		 *	@note This function throws std::bad_alloc.
		 *	@note This function throws std::runtime_error (in case vertex ids aren't contiguous).
		 */
		template <class CVertexType, class CInitializer>
		static CVertexType &r_Get_Vertex(_TyId n_id, CInitializer init,
			_TyVertexMultiPool &r_vertex_pool, _TyFixVertexMultiPool &UNUSED(r_fix_vertex_pool),
			_TyVertexLookup &UNUSED(r_vertex_lookup), size_t &r_n_vertex_element_num) // throw(std::bad_alloc)
		{
			if(n_id == r_vertex_pool.n_Size()) {
				CVertexType &r_t_new_vert = r_vertex_pool.r_Add_Element((CVertexType)(init));
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
				r_t_new_vert.Set_Id(n_id);
#endif // __BASE_TYPES_USE_ID_ADDRESSING
				r_t_new_vert.Set_Order(r_n_vertex_element_num);
				r_n_vertex_element_num += r_t_new_vert.n_Dimension();
				return r_t_new_vert;
			} else {
				if(n_id < r_vertex_pool.n_Size())
					return *(CVertexType*)&r_vertex_pool[n_id]; // a bit of dirty casting
				throw std::runtime_error("vertices must be accessed in incremental manner");

				// if we want to support fixed vertices, we need vertex map for this
				// but otherwise the allocation of the matrices is taken care of by assigning id / order
				// (need to specify that id is *column* id, *not* vertex / edge id)

				// alternately we could keep fixed vertices in the same pool, which will probably break stuff

				// finally, we could have pointer lists for fixed / not fixed vertices
				// that could almost allow changing the fixed-ness on the fly
			}
		}
	};

	/**
	 *	@brief code for vertex insertion, based on whether the fixed vertices are allowed
	 *		(specialization for the case when the fixed vertices are allowed)
	 */
	template <class _CBaseVertex, class _CVertexTypelist, class _CBaseEdge,
		class _CEdgeTypelist, class _CUnaryFactorFactory, const int _n_pool_page_size>
	class CGetVertexImpl<_CBaseVertex, _CVertexTypelist, _CBaseEdge,
		_CEdgeTypelist, _CUnaryFactorFactory, true, _n_pool_page_size> {
	public:
		/**
		 *	@brief finds vertex by id, or create a new one in case there is no vertex with such id
		 *
		 *	@tparam CVertexType is the vertex type name
		 *	@tparam CInitializer is the lazy initializer type name
		 *
		 *	@param[in] n_id is the id of the vertex required (must be id of an existing vertex,
		 *		or one larger than id of the last vertex)
		 *	@param[in] init is the initializer functor; it is only called (converted to vertex)
		 *		in case new vertex is created
		 *	@param[in,out] r_vertex_pool is vertex pool for optimized vertices
		 *	@param[in,out] r_fix_vertex_pool is vertex pool for fixed vertices
		 *	@param[in,out] r_vertex_lookup is lookup table that contains both fixed and optimized vertices
		 *	@param[in,out] r_n_vertex_element_num is number of elements of all the
		 *		optimized vertices so far
		 *
		 *	@return Returns reference to the vertex, associated with id n_id.
		 *
		 *	@note This function throws std::bad_alloc.
		 *	@note This function throws std::runtime_error (in case vertex ids aren't contiguous).
		 */
		template <class CVertexType, class CInitializer>
		static CVertexType &r_Get_Vertex(_TyId n_id, CInitializer init,
			_TyVertexMultiPool &r_vertex_pool, _TyFixVertexMultiPool &r_fix_vertex_pool,
			_TyVertexLookup &r_vertex_lookup, size_t &r_n_vertex_element_num) // throw(std::bad_alloc)
		{
			if(n_id == r_vertex_lookup.n_Size()) {
				CVertexType t_new_vert = (CVertexType)(init);
				// make a new vertex

				bool b_fixed = t_new_vert.b_IsFixed();
				// determine whether the vertex is fixed

				CVertexType &r_t_new_vert = ((b_fixed)? r_fix_vertex_pool :
					r_vertex_pool).r_Add_Element(t_new_vert);
				// add to the appropriate pool

				if(b_fixed) {
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
					r_t_new_vert.Set_Id(-1); // no id; not to be placed in the system matrix
#endif // __BASE_TYPES_USE_ID_ADDRESSING
					r_t_new_vert.Set_Order(-1); // no order; not to be placed in the system matrix
				} else {
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
					r_t_new_vert.Set_Id(r_vertex_pool.n_Size() - 1); // id in vertices that are not fixed
#endif // __BASE_TYPES_USE_ID_ADDRESSING
					r_t_new_vert.Set_Order(r_n_vertex_element_num); // order only of vertices that are not fixed
					r_n_vertex_element_num += r_t_new_vert.n_Dimension();
				}
				// take care of ordering

				r_vertex_lookup.push_back(&r_t_new_vert);
				// add pointer to the lookup table

				return r_t_new_vert;
			} else {
				if(n_id < r_vertex_lookup.n_Size())
					return *(CVertexType*)&r_vertex_lookup[n_id]; // a bit of dirty casting
				throw std::runtime_error("vertices must be accessed in incremental manner");
			}
		}
	};

	typedef CGetVertexImpl<CBaseVertex, CVertexTypelist, CBaseEdge, CEdgeTypelist,
		CUnaryFactorFactory, b_allow_fixed_vertices, n_pool_page_size> _TyGetVertexImpl; /**< @brief type with specialized implementation of GetVertex function */

protected:
	_TyVertexLookup m_vertex_lookup; /**< @brief lookup table for vertices (by id; only used if allow_FixedVertices is set) */
	_TyVertexMultiPool m_vertex_pool; /**< @brief vertex multipool */
	_TyEdgeMultiPool m_edge_pool; /**< @brief edge multipool */
	size_t m_n_vertex_element_num; /**< @brief sum of the numbers of all vertex elements, (also the size of permutation vector) */
	size_t m_n_edge_element_num; /**< @brief sum of all permutation vector components + unary factor rank (also size of the A matrix) */

	_TyFixVertexMultiPool m_fixed_vertex_pool; /**< @brief fixed vertex multipool (only used if allow_FixedVertices is set) */

	CUnaryFactorFactory m_unary_factor_factory; /**< @brief unary factor initializer */
	Eigen::MatrixXd m_t_unary_factor; /**< @brief the unary factor matrix */
	Eigen::VectorXd m_v_unary_error; /**< @brief the error vector associated with the first vertex */

public:
	/**
	 *	@brief default constructor; has no effect
	 *	@param[in] unary_factor_factory is unary factor initializer
	 */
	inline CFlatSystem(CUnaryFactorFactory unary_factor_factory = CUnaryFactorFactory())
		:m_n_vertex_element_num(0), m_n_edge_element_num(0), m_unary_factor_factory(unary_factor_factory)
	{}

	/**
	 *	@brief gets the unary factor matrix
	 *	@return Returns reference to the unary factor matrix.
	 */
	Eigen::MatrixXd &r_t_Unary_Factor()
	{
		return m_t_unary_factor;
	}

	/**
	 *	@brief gets the unary factor matrix
	 *	@return Returns const reference to the unary factor matrix.
	 */
	const Eigen::MatrixXd &r_t_Unary_Factor() const
	{
		return m_t_unary_factor;
	}

	/**
	 *	@brief gets the error vector, associated with unary factor
	 *	@return Returns reference to the error vector, associated with unary factor.
	 */
	Eigen::VectorXd &r_v_Unary_Error()
	{
		return m_v_unary_error;
	}

	/**
	 *	@brief gets the error vector, associated with unary factor
	 *	@return Returns const reference to the error vector, associated with unary factor.
	 */
	const Eigen::VectorXd &r_v_Unary_Error() const
	{
		return m_v_unary_error;
	}

	/**
	 *	@brief gets number of edges in the system
	 *	@return Returns number of edges in the system.
	 */
	inline size_t n_Edge_Num() const
	{
		return m_edge_pool.n_Size();
	}

	/**
	 *	@brief gets number of vertices in the system
	 *	@return Returns number of vertices in the system.
	 *	@note This is only number of vertices that are not fixed,
	 *		and that are potentially included in the system matrix.
	 */
	inline size_t n_Vertex_Num() const
	{
		return m_vertex_pool.n_Size();
	}

	/**
	 *	@brief gets number of fixed vertices in the system
	 *	@return Returns number of fixed vertices in the system.
	 */
	inline size_t n_FixedVertex_Num() const
	{
		return m_fixed_vertex_pool.n_Size();
	}

	/**
	 *	@brief gets the number of all permutation vector components + unary factor rank
	 *	@return Returns the sum of all permutation vector components + unary factor rank (also the size of the A matrix).
	 */
	inline size_t n_EdgeElement_Num() const
	{
		return m_n_edge_element_num;
	}

	/**
	 *	@brief gets the size of the permutation vector
	 *	@return Returns the sum of the numbers of all vertex elements, (also the size of permutation vector).
	 *	@note This is only sum of numbers of elements of vertices that are not fixed,
	 *		and that are potentially included in the system matrix.
	 */
	inline size_t n_VertexElement_Num() const
	{
		return m_n_vertex_element_num;
	}

	/**
	 *	@brief gets reference to edge multipool
	 *	@return Returns reference to edge multipool.
	 */
	inline _TyEdgeMultiPool &r_Edge_Pool()
	{
		return m_edge_pool;
	}

	/**
	 *	@brief gets reference to edge multipool
	 *	@return Returns const reference to edge multipool.
	 */
	inline const _TyEdgeMultiPool &r_Edge_Pool() const
	{
		return m_edge_pool;
	}

	/**
	 *	@brief gets reference to vertex multipool
	 *	@return Returns reference to vertex multipool.
	 *	@note This pool contains only vertices that are not fixed,
	 *		and that are potentially included in the system matrix.
	 */
	inline _TyVertexMultiPool &r_Vertex_Pool()
	{
		return m_vertex_pool;
	}

	/**
	 *	@brief gets reference to vertex multipool
	 *	@return Returns const reference to vertex multipool.
	 *	@note This pool contains only vertices that are not fixed,
	 *		and that are potentially included in the system matrix.
	 */
	inline const _TyVertexMultiPool &r_Vertex_Pool() const
	{
		return m_vertex_pool;
	}

	/**
	 *	@brief gets reference to fixed vertex multipool
	 *	@return Returns reference to fixed vertex multipool.
	 *	@note This is only available if allow_FixedVertices is set.
	 */
	inline _TyFixVertexMultiPool &r_FixedVertex_Pool()
	{
		return m_fixed_vertex_pool;
	}

	/**
	 *	@brief gets reference to fixed vertex multipool
	 *	@return Returns const reference to fixed vertex multipool.
	 *	@note This is only available if allow_FixedVertices is set.
	 */
	inline const _TyFixVertexMultiPool &r_FixedVertex_Pool() const
	{
		return m_fixed_vertex_pool;
	}

	/**
	 *	@brief finds vertex by id, or create a new one in case there is no vertex with such id
	 *
	 *	@tparam CVertexType is the vertex type name
	 *	@tparam CInitializer is the lazy initializer type name
	 *
	 *	@param[in] n_id is the id of the vertex required (must be id of an existing vertex,
	 *		or one larger than id of the last vertex)
	 *	@param[in] init is the initializer functor; it is only called (converted to vertex) in case new vertex is created
	 *
	 *	@return Returns reference to the vertex, associated with id n_id.
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This function throws std::runtime_error (in case vertex ids aren't contiguous).
	 */
	template <class CVertexType, class CInitializer>
	inline CVertexType &r_Get_Vertex(_TyId n_id, CInitializer init) // throw(std::bad_alloc)
	{
		return _TyGetVertexImpl::template r_Get_Vertex<CVertexType>(n_id, init,
			m_vertex_pool, m_fixed_vertex_pool, m_vertex_lookup, m_n_vertex_element_num);
		// use specialized implementation, based on whether fixed vertices are enabled or not

		/*if(n_id == m_vertex_pool.n_Size()) {
			CVertexType &r_t_new_vert = m_vertex_pool.r_Add_Element((CVertexType)(init));
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
			r_t_new_vert.Set_Id(n_id);
#endif // __BASE_TYPES_USE_ID_ADDRESSING
			r_t_new_vert.Set_Order(m_n_vertex_element_num);
			m_n_vertex_element_num += r_t_new_vert.n_Dimension();
			return r_t_new_vert;
		} else {
			if(n_id < m_vertex_pool.n_Size())
				return *(CVertexType*)&m_vertex_pool[n_id]; // a bit of dirty casting
			throw std::runtime_error("vertices must be accessed in incremental manner");
		}*/
		// old implementation that does not allow fixed vertices
	}

	/**
	 *	@brief adds a new edge to the system
	 *
	 *	@param[in] t_edge is the initial edge value
	 *
	 *	@return Returns reference to the new edge.
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This sets edge order, do not overwrite it!
	 */
	template <class _TyEdge>
	_TyEdge &r_Add_Edge(_TyEdge t_edge) // throw(std::bad_alloc)
	{
		bool b_was_empty;
		if((b_was_empty = m_edge_pool.b_Empty())) {
			m_unary_factor_factory(m_t_unary_factor, m_v_unary_error, t_edge);
			m_n_edge_element_num += m_t_unary_factor.rows();
		}
		// in case it is the first edge, gets the unary factor and the error associated with it

		_TyEdge &r_t_new_edge = m_edge_pool.r_Add_Element(t_edge);
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
		r_t_new_edge.Set_Id(m_edge_pool.n_Size() - 1);
#endif // __BASE_TYPES_USE_ID_ADDRESSING
		r_t_new_edge.Set_Order(m_n_edge_element_num);

		m_n_edge_element_num += r_t_new_edge.n_Dimension();
		// the first edge is preceded by the unary factor

		return r_t_new_edge;
	}

	/**
	 *	@brief plots the system as a .tga (targa) image
	 *
	 *	@param[in] p_s_filename is the output file name
	 *	@param[in] n_quality_profile is plot profile (one of plot_*)
	 *
	 *	@return Returns true on success, false on failure.
	 */
	bool Plot2D(const char *p_s_filename, plot_quality::EQualityProfile n_quality_profile) // todo - do the same for 3D
	{
		switch(n_quality_profile) {
		case plot_quality::plot_Draft:
			return Plot2D(p_s_filename, 1024, 8192, 2, 1, 1, 1, false, true, 32, false);
		case plot_quality::plot_Printing:
			return Plot2D(p_s_filename, 2048, 2048, 10, 3, 7, 1, true, false, 10);
		case plot_quality::plot_Printing_LandmarkTicksOnly:
			return Plot2D(p_s_filename, 2048, 2048, 10, 3, 7, 1, true, false, 10, true);
		case plot_quality::plot_Printing_NoTicks:
			return Plot2D(p_s_filename, 2048, 2048, 0, 0, 7, 1, true, false, 4);
		};
		return false;
	}

	/**
	 *	@brief plots the system as a .tga (targa) image
	 *
	 *	@param[in] p_s_filename is the output file name
	 *	@param[in] n_smaller_side_resoluion is the resolution of the shorter side of the image. in pixels
	 *	@param[in] n_max_resolution is the maximal resolution limit for either dimension, in pixels
	 *		(in case the aspect ratio is too high, the longer side is set to n_max_resolution
	 *		and the shorter side is scaled appropriately)
	 *	@param[in] n_cross_size_px is vertex cross size, in pixels (advanced settings)
	 *	@param[in] n_cross_line_width is vertex cross line width, in pixels (advanced settings)
	 *	@param[in] n_edge_line_width is edge line width, in pixels (advanced settings)
	 *	@param[in] n_landmark_edge_line_width is pose-landmark edge line width, in pixels (advanced settings)
	 *	@param[in] b_dark_landmark_edges is pose-landmark color selector; if set,
	 *		it is black, otherwise it is light gray (advanced settings)
	 *	@param[in] b_draw_frame is flag for drawing a frame arround the image
	 *		(the frame is fit tightly arround the image)
	 *	@param[in] n_padding is padding arround the image (arround the frame, if rendered)
	 *	@param[in] b_landmark_ticks_only is vertex rendering flag (if set, only landmarks are rendered,
	 *		otherwise all vertices are rendered)
	 *
	 *	@return Returns true on success, false on failure.
	 */
	bool Plot2D(const char *p_s_filename, int n_smaller_side_resoluion = 1024,
		int n_max_resolution = 8192, int n_cross_size_px = 2, int n_cross_line_width = 1,
		int n_edge_line_width = 1, int n_landmark_edge_line_width = 1,
		bool b_dark_landmark_edges = false, bool b_draw_frame = true, int n_padding = 32,
		bool b_landmark_ticks_only = false)
	{
		Eigen::Vector2d v_min(-1, -1), v_max(1, 1);
		if(!m_vertex_pool.b_Empty()) {
			const Eigen::VectorXd &r_v_state = m_vertex_pool[0].v_State();
			_ASSERTE(r_v_state.rows() >= 2);
			v_min(0) = v_max(0) = r_v_state(0);
			v_min(1) = v_max(1) = r_v_state(1);
			for(size_t i = 0, n = m_vertex_pool.n_Size(); i < n; ++ i) {
				const _TyBaseVertex &r_t_vert = m_vertex_pool[i];
				const Eigen::VectorXd &r_v_state = r_t_vert.v_State();
				for(int j = 0; j < 2; ++ j) {
					v_min(j) = std::min(v_min(j), r_v_state(j));
					v_max(j) = std::max(v_max(j), r_v_state(j));
				}
			}
		}
		// find minima / maxima

		Eigen::Vector2d v_size = v_max - v_min;
		double f_short_side = std::min(v_size(0), v_size(1));
		double f_long_side = std::max(v_size(0), v_size(1));
		double f_scale = n_smaller_side_resoluion / f_short_side;
		if(f_long_side * f_scale > n_max_resolution)
			f_scale = n_max_resolution / f_long_side;

		int n_width = std::max(1, int(v_size(0) * f_scale));
		int n_height = std::max(1, int(v_size(1) * f_scale));
		// calculate image size

		TBmp *p_image;
		if(!(p_image = TBmp::p_Alloc(n_width + 2 * n_padding, n_height + 2 * n_padding)))
			return false;

		p_image->Clear(0xffffffU);
		// white background

		if(b_draw_frame) {
			n_padding -= n_cross_size_px;
			n_width += 2 * n_cross_size_px;
			n_height += 2 * n_cross_size_px;
			p_image->DrawRect(n_padding, n_padding, n_width + n_padding, n_height + n_padding, 0xff000000U);
			/*p_image->DrawLine(n_padding, n_padding, n_width + n_padding, n_padding, 0xff000000U);
			p_image->DrawLine(n_padding, n_height + n_padding, n_width + n_padding, n_height + n_padding, 0xff000000U);
			p_image->DrawLine(n_padding, n_padding, n_padding, n_height + n_padding, 0xff000000U);
			p_image->DrawLine(n_width + n_padding, n_padding, n_width + n_padding, n_height + n_padding, 0xff000000U);*/
			n_padding += n_cross_size_px;
			n_width -= 2 * n_cross_size_px;
			n_height -= 2 * n_cross_size_px;
		}
		// black borders

		for(size_t i = 0, n = m_edge_pool.n_Size(); i < n; ++ i) {
			const _TyBaseEdge &r_edge = m_edge_pool[i];

			if(r_edge.n_Dimension() != 2)
				continue;
			// only landmarks

			const Eigen::VectorXd &r_t_vert0 = m_vertex_pool[r_edge.n_Vertex_Id(0)].v_State();
			const Eigen::VectorXd &r_t_vert1 = m_vertex_pool[r_edge.n_Vertex_Id(1)].v_State();
			_ASSERTE(r_t_vert0.rows() >= 2 && r_t_vert1.rows() >= 2);

			float f_x0 = float((r_t_vert0(0) - v_min(0)) * f_scale);
			float f_y0 = float((r_t_vert0(1) - v_min(1)) * f_scale);
			f_x0 = f_x0 + n_padding;
			f_y0 = n_height - f_y0 + n_padding; // the origin is at the top left corner, flip y

			float f_x1 = float((r_t_vert1(0) - v_min(0)) * f_scale);
			float f_y1 = float((r_t_vert1(1) - v_min(1)) * f_scale);
			f_x1 = f_x1 + n_padding;
			f_y1 = n_height - f_y1 + n_padding; // the origin is at the top left corner, flip y

			uint32_t n_color = (r_edge.n_Dimension() == 2)? 0xffeeeeeeU : 0xff0000ffU;

			if(b_dark_landmark_edges)
				n_color = 0xffaaaaaaU;

			p_image->DrawLine_AA(f_x0, f_y0, f_x1, f_y1, n_color, n_landmark_edge_line_width);
		}
		// draw edges

		for(size_t i = 0, n = m_vertex_pool.n_Size(); i < n; ++ i) {
			const _TyBaseVertex &r_t_vert = m_vertex_pool[i];
			const Eigen::VectorXd &r_v_state = r_t_vert.v_State();

			if(b_landmark_ticks_only && r_t_vert.n_Dimension() != 2)
				continue;

			float f_x = float((r_v_state(0) - v_min(0)) * f_scale);
			float f_y = float((r_v_state(1) - v_min(1)) * f_scale);
			f_x = f_x + n_padding;
			f_y = n_height - f_y + n_padding; // the origin is at the top left corner, flip y

			p_image->DrawLine_SP(f_x - n_cross_size_px, f_y,
				f_x + n_cross_size_px, f_y, 0xffff0000U, n_cross_line_width);
			p_image->DrawLine_SP(f_x, f_y - n_cross_size_px,
				f_x, f_y + n_cross_size_px, 0xffff0000U, n_cross_line_width);
		}
		// draw vertices

		// todo - support fixed vertices

		for(size_t i = 0, n = m_edge_pool.n_Size(); i < n; ++ i) {
			const _TyBaseEdge &r_edge = m_edge_pool[i];

			if(r_edge.n_Dimension() == 2)
				continue;
			// only edges

			const Eigen::VectorXd &r_t_vert0 = m_vertex_pool[r_edge.n_Vertex_Id(0)].v_State();
			const Eigen::VectorXd &r_t_vert1 = m_vertex_pool[r_edge.n_Vertex_Id(1)].v_State();
			_ASSERTE(r_t_vert0.rows() >= 2 && r_t_vert1.rows() >= 2);

			float f_x0 = float((r_t_vert0(0) - v_min(0)) * f_scale);
			float f_y0 = float((r_t_vert0(1) - v_min(1)) * f_scale);
			f_x0 = f_x0 + n_padding;
			f_y0 = n_height - f_y0 + n_padding; // the origin is at the top left corner, flip y

			float f_x1 = float((r_t_vert1(0) - v_min(0)) * f_scale);
			float f_y1 = float((r_t_vert1(1) - v_min(1)) * f_scale);
			f_x1 = f_x1 + n_padding;
			f_y1 = n_height - f_y1 + n_padding; // the origin is at the top left corner, flip y

			uint32_t n_color = (r_edge.n_Dimension() == 2)? 0xffeeeeeeU : 0xff0000ffU;

			p_image->DrawLine_AA(f_x0, f_y0, f_x1, f_y1, n_color, n_edge_line_width);
		}
		// draw edges

		bool b_result = CTgaCodec::Save_TGA(p_s_filename, *p_image, false);
		p_image->Delete();

		return b_result;
	}

	/**
	 *	@brief plots the system as a .tga (targa) image
	 *
	 *	This version assumes a 3D system, it plots it as its projection to the XZ plane.
	 *
	 *	@param[in] p_s_filename is the output file name
	 *	@param[in] n_smaller_side_resoluion is the resolution of the shorter side of the image. in pixels
	 *	@param[in] n_max_resolution is the maximal resolution limit for either dimension, in pixels
	 *		(in case the aspect ratio is too high, the longer side is set to n_max_resolution
	 *		and the shorter side is scaled appropriately)
	 *	@param[in] n_cross_size_px is vertex cross size, in pixels (advanced settings)
	 *	@param[in] n_cross_line_width is vertex cross line width, in pixels (advanced settings)
	 *	@param[in] n_edge_line_width is edge line width, in pixels (advanced settings)
	 *	@param[in] n_landmark_edge_line_width is pose-landmark edge line width, in pixels (advanced settings)
	 *	@param[in] b_dark_landmark_edges is pose-landmark color selector; if set,
	 *		it is black, otherwise it is light gray (advanced settings)
	 *	@param[in] b_draw_frame is flag for drawing a frame arround the image
	 *		(the frame is fit tightly arround the image)
	 *	@param[in] n_padding is padding arround the image (arround the frame, if rendered)
	 *	@param[in] b_landmark_ticks_only is vertex rendering flag (if set, only landmarks are rendered,
	 *		otherwise all vertices are rendered)
	 *
	 *	@return Returns true on success, false on failure.
	 */
	bool Plot3D(const char *p_s_filename, int n_smaller_side_resoluion = 1024,
		int n_max_resolution = 8192, int n_cross_size_px = 2, int n_cross_line_width = 1,
		int n_edge_line_width = 1, int n_landmark_edge_line_width = 1,
		bool b_dark_landmark_edges = false, bool b_draw_frame = true, int n_padding = 32,
		bool b_landmark_ticks_only = false)
	{
		Eigen::Vector2d v_min(-1, -1), v_max(1, 1);
		if(!m_vertex_pool.b_Empty()) {
			const Eigen::VectorXd &r_v_state = m_vertex_pool[0].v_State();
			_ASSERTE(r_v_state.rows() >= 3);
			v_min(0) = v_max(0) = r_v_state(0);
			v_min(1) = v_max(1) = r_v_state(2);
			for(size_t i = 0, n = m_vertex_pool.n_Size(); i < n; ++ i) {
				const _TyBaseVertex &r_t_vert = m_vertex_pool[i];
				const Eigen::VectorXd &r_v_state = r_t_vert.v_State();
				for(int j = 0; j < 2; ++ j) {
					v_min(j) = std::min(v_min(j), r_v_state((j)? 2 : 0));
					v_max(j) = std::max(v_max(j), r_v_state((j)? 2 : 0));
				}
			}
		}
		// find minima / maxima

		Eigen::Vector2d v_size = v_max - v_min;
		double f_short_side = std::min(v_size(0), v_size(1));
		double f_long_side = std::max(v_size(0), v_size(1));
		double f_scale = n_smaller_side_resoluion / f_short_side;
		if(f_long_side * f_scale > n_max_resolution)
			f_scale = n_max_resolution / f_long_side;

		int n_width = std::max(1, int(v_size(0) * f_scale));
		int n_height = std::max(1, int(v_size(1) * f_scale));
		// calculate image size

		TBmp *p_image;
		if(!(p_image = TBmp::p_Alloc(n_width + 2 * n_padding, n_height + 2 * n_padding)))
			return false;

		p_image->Clear(0xffffffU);
		// white background

		if(b_draw_frame) {
			n_padding -= n_cross_size_px;
			n_width += 2 * n_cross_size_px;
			n_height += 2 * n_cross_size_px;
			p_image->DrawRect(n_padding, n_padding, n_width + n_padding, n_height + n_padding, 0xff000000U);
			/*p_image->DrawLine(n_padding, n_padding, n_width + n_padding, n_padding, 0xff000000U);
			p_image->DrawLine(n_padding, n_height + n_padding, n_width + n_padding, n_height + n_padding, 0xff000000U);
			p_image->DrawLine(n_padding, n_padding, n_padding, n_height + n_padding, 0xff000000U);
			p_image->DrawLine(n_width + n_padding, n_padding, n_width + n_padding, n_height + n_padding, 0xff000000U);*/
			n_padding += n_cross_size_px;
			n_width -= 2 * n_cross_size_px;
			n_height -= 2 * n_cross_size_px;
		}
		// black borders

		/*for(size_t i = 0, n = m_edge_pool.n_Size(); i < n; ++ i) {
			const _TyBaseEdge &r_edge = m_edge_pool[i];

			if(r_edge.n_Dimension() != 2)
				continue;
			// only landmarks

			const Eigen::VectorXd &r_t_vert0 = m_vertex_pool[r_edge.n_Vertex_Id(0)].v_State();
			const Eigen::VectorXd &r_t_vert1 = m_vertex_pool[r_edge.n_Vertex_Id(1)].v_State();
			_ASSERTE(r_t_vert0.rows() >= 2 && r_t_vert1.rows() >= 2);

			float f_x0 = float((r_t_vert0(0) - v_min(0)) * f_scale);
			float f_y0 = float((r_t_vert0(1) - v_min(1)) * f_scale);
			f_x0 = f_x0 + n_padding;
			f_y0 = n_height - f_y0 + n_padding; // the origin is at the top left corner, flip y

			float f_x1 = float((r_t_vert1(0) - v_min(0)) * f_scale);
			float f_y1 = float((r_t_vert1(1) - v_min(1)) * f_scale);
			f_x1 = f_x1 + n_padding;
			f_y1 = n_height - f_y1 + n_padding; // the origin is at the top left corner, flip y

			uint32_t n_color = (r_edge.n_Dimension() == 2)? 0xffeeeeeeU : 0xff0000ffU;

			if(b_dark_landmark_edges)
				n_color = 0xffaaaaaaU;

			p_image->DrawLine_AA(f_x0, f_y0, f_x1, f_y1, n_color, n_landmark_edge_line_width);
		}*/
		// draw edges

		for(size_t i = 0, n = m_vertex_pool.n_Size(); i < n; ++ i) {
			const _TyBaseVertex &r_t_vert = m_vertex_pool[i];
			const Eigen::VectorXd &r_v_state = r_t_vert.v_State();

			if(b_landmark_ticks_only && r_t_vert.n_Dimension() != 2)
				continue;

			float f_x = float((r_v_state(0) - v_min(0)) * f_scale);
			float f_y = float((r_v_state(2) - v_min(1)) * f_scale);
			f_x = f_x + n_padding;
			f_y = n_height - f_y + n_padding; // the origin is at the top left corner, flip y

			p_image->DrawLine_SP(f_x - n_cross_size_px, f_y,
				f_x + n_cross_size_px, f_y, 0xffff0000U, n_cross_line_width);
			p_image->DrawLine_SP(f_x, f_y - n_cross_size_px,
				f_x, f_y + n_cross_size_px, 0xffff0000U, n_cross_line_width);
		}
		// draw vertices

		// todo - support fixed vertices

		for(size_t i = 0, n = m_edge_pool.n_Size(); i < n; ++ i) {
			const _TyBaseEdge &r_edge = m_edge_pool[i];

			if(r_edge.n_Dimension() == 2)
				continue;
			// only edges

			const Eigen::VectorXd &r_t_vert0 = m_vertex_pool[r_edge.n_Vertex_Id(0)].v_State();
			const Eigen::VectorXd &r_t_vert1 = m_vertex_pool[r_edge.n_Vertex_Id(1)].v_State();
			_ASSERTE(r_t_vert0.rows() >= 2 && r_t_vert1.rows() >= 2);

			float f_x0 = float((r_t_vert0(0) - v_min(0)) * f_scale);
			float f_y0 = float((r_t_vert0(2) - v_min(1)) * f_scale);
			f_x0 = f_x0 + n_padding;
			f_y0 = n_height - f_y0 + n_padding; // the origin is at the top left corner, flip y

			float f_x1 = float((r_t_vert1(0) - v_min(0)) * f_scale);
			float f_y1 = float((r_t_vert1(2) - v_min(1)) * f_scale);
			f_x1 = f_x1 + n_padding;
			f_y1 = n_height - f_y1 + n_padding; // the origin is at the top left corner, flip y

			uint32_t n_color = (r_edge.n_Dimension() == 2)? 0xffeeeeeeU : 0xff0000ffU;

			p_image->DrawLine_AA(f_x0, f_y0, f_x1, f_y1, n_color, n_edge_line_width);
		}
		// draw edges

		bool b_result = CTgaCodec::Save_TGA(p_s_filename, *p_image, false);
		p_image->Delete();

		return b_result;
	}

	/**
	 *	@brief saves (2D) vertex positions into a text file
	 *	@param[in] p_s_filename is the output file name
	 *	@return Returns true on success, false on failure.
	 */
	bool Dump(const char *p_s_filename) const
	{
		FILE *p_fw;
		if(!(p_fw = fopen(p_s_filename, "w"))) {
			//fprintf(stderr, "error: failed to open \'%\' for writing\n", p_s_filename);
			return false;
		}

		for(size_t i = 0, n = m_vertex_pool.n_Size(); i < n; ++ i) {
			const _TyBaseVertex &r_vert = m_vertex_pool[i];
			const Eigen::VectorXd &r_v_state = r_vert.v_State();
			for(size_t j = 0, m = r_v_state.rows(); j < m; ++ j)
				fprintf(p_fw, (j)? ((j + 1 == m)? " %f\n" : " %f") : ((j + 1 == m)? "%f\n" : "%f"), r_v_state(j)); // t_odo - support any dimensionality
		}

		// todo - support const vertices

		if(ferror(p_fw)) {
			//fprintf(stderr, "error: error writing \'%\' (disk full?)\n", p_s_filename);
			fclose(p_fw);
			return false;
		}
		fclose(p_fw);

		return true;
	}
};

#endif // __FLAT_SYSTEM_TEMPLATE_INCLUDED
