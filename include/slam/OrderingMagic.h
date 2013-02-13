/*
								+-----------------------------------+
								|                                   |
								| *** Matrix ordering utilities *** |
								|                                   |
								|   Copyright  © -tHE SWINe- 2013   |
								|                                   |
								|          OrderingMagic.h          |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __MATRIX_ORDERING_UTILS_INCLUDED
#define __MATRIX_ORDERING_UTILS_INCLUDED

/**
 *	@file include/slam/OrderingMagic.h
 *	@brief matrix ordering utilities
 *	@author -tHE SWINe-
 *	@date 2013-02-07
 */

/**
 *	@def __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT
 *	@brief if enabled, two-level ordering constraint is applied to lambda
 *		in CLastElementOrderingConstraint
 */
#define __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT

#include <vector>
#include <algorithm>
#include "csparse/cs.hpp"
#include "slam/BlockMatrix.h"

/**
 *	@brief ordering constraint for C(COL)AMD libraries
 *
 *	This maintans a constraint vector which forces the last element
 *	to be the last even after symbolic ordering.
 */
class CLastElementOrderingConstraint {
protected:
	std::vector<size_t> m_constraints; /**< @brief storage for the constraint vector */

public:
	/**
	 *	@brief gets the constraint vector of a specified size
	 *	@param[in] n_size is the size of the constraint vector (in elements)
	 *	@return Returns const pointer to the constraint vector
	 *		of the specified size (not to be deleted).
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_Get(size_t n_size); // throws(std::bad_alloc)
};

/**
 *	@brief ordering constraint for C(COL)AMD libraries
 *
 *	This maintans a constraint vector which forces the first and the last elements
 *	to be the first and the last, respecitvely, even after symbolic ordering.
 */
class CFirstLastElementOrderingConstraint {
protected:
	std::vector<size_t> m_constraints; /**< @brief storage for the constraint vector */

public:
	/**
	 *	@brief gets the constraint vector of a specified size
	 *	@param[in] n_size is the size of the constraint vector (in elements)
	 *	@return Returns const pointer to the constraint vector
	 *		of the specified size (not to be deleted).
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_Get(size_t n_size); // throws(std::bad_alloc)
};

/**
 *	@brief ordering constraint for C(COL)AMD libraries
 *
 *	This maintans a constraint vector which forces the first and the last elements
 *	to be the first and the last, respecitvely, even after symbolic ordering.
 */
class CNFirst1LastElementOrderingConstraint {
protected:
	std::vector<size_t> m_constraints; /**< @brief storage for the constraint vector */

public:
	/**
	 *	@brief gets the constraint vector of a specified size
	 *
	 *	@param[in] n_size is the size of the constraint vector (in elements)
	 *	@param[in] n_first_constraint_num is number of the constraints of the first columns
	 *
	 *	@return Returns const pointer to the constraint vector
	 *		of the specified size (not to be deleted).
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_Get(size_t n_size, size_t n_first_constraint_num); // throws(std::bad_alloc)
};

/**
 *	@brief matrix ordering calculator (CAMD wrapper)
 */
class CMatrixOrdering { // t_odo - fill throws, and in related classes
protected:
	std::vector<size_t> m_ordering; /**< @brief storage for the blockwise ordering vector */
	std::vector<size_t> m_ordering_expand; /**< @brief storage for the elementwise ordering vector */
	std::vector<size_t> m_ordering_invert; /**< @brief storage for the inverse elementwise ordering vector */
	cs *m_p_block; /**< @brief matrix block structure (reuses memory storage) */

	std::vector<size_t> m_camd_workspace; /**< @brief workspace for camd_aat() */
	std::vector<size_t> m_camd_workspace1; /**< @brief workspace for camd_1() */

public:
	/**
	 *	@brief default constructor; has no effect
	 */
	inline CMatrixOrdering()
		:m_p_block(0)
	{}

	/**
	 *	@brief destructor; deletes allocated memory
	 */
	~CMatrixOrdering();

	/**
	 *	@brief calculates inverse ordering while maintaining and reusing storage
	 *
	 *	@param[in] p_ordering is an ordering to be inverted
	 *	@param[in] n_ordering_size is size of the given ordering
	 *
	 *	@return Returns const pointer to the inverse of the given ordering (not to be deleted).
	 *
	 *	@note The buffer for inverse ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of CMatrixOrdering
	 *		if multiple inverse orderings need to exist at the same time).
	 *	@note Due to the buffer reuse, this can not invert ordering, inverted
	 *		by the same object (calling the same function on the result of the previous
	 *		call, on the same object). This doesn't apply to the results of the other functions.
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_InvertOrdering(const size_t *p_ordering, size_t n_ordering_size); // throws(std::bad_alloc)

	/**
	 *	@brief extremely intricate function for progressive reordering
	 *
	 *	@param[in] n_order_min is zero-based index of the first element of the ordering to be modified
	 *	@param[in] p_sub_block_ordering is pointer to incremental ordering
	 *	@param[in] n_sub_block_size is size of the incremental ordering; in case it
	 *		is longer than the current ordering, it is extended with identity
	 *
	 *	@return Returns const pointer to the extended ordering (not to be deleted).
	 *
	 *	@note This function invalidates the inverse ordering; needs to be recalculated.
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 */
	const size_t *p_ExtendBlockOrdering_with_SubOrdering(size_t n_order_min,
		const size_t *p_sub_block_ordering, size_t n_sub_block_size); // throws(std::bad_alloc)

	/**
	 *	@brief calculates blockwise ordering in a matrix, using the CAMD library
	 *
	 *	@param[in] r_A is the block matrix (must be symmetric)
	 *	@param[in] p_constraints is the ordering constraint vector (can be 0)
	 *	@param[in] n_constraints_size is size og the ordering constraint vector
	 *		(must match the number of r_A block columns, except it is ignored
	 *		if p_constraints is 0)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of CMatrixOrdering
	 *		if multiple block orderings need to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_BlockOrdering(const CUberBlockMatrix &r_A, const size_t *p_constraints,
		size_t UNUSED(n_constraints_size)); // throws(std::bad_alloc)

	/**
	 *	@brief calculates blockwise ordering in a matrix, using the CAMD library
	 *
	 *	@param[in] r_A is the block matrix (must be symmetric)
	 *	@param[in] n_off_diagonal_num is the number of off-diagonal elements
	 *		(equals number of edges in binary graphs or the sum of edge ranks
	 *		minus the number of vertices in hypergraphs)
	 *	@param[in] p_constraints is the ordering constraint vector (can be 0)
	 *	@param[in] n_constraints_size is size og the ordering constraint vector
	 *		(must match the number of r_A block columns, except it is ignored
	 *		if p_constraints is 0)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of CMatrixOrdering
	 *		if multiple block orderings need to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_HybridBlockOrdering(const CUberBlockMatrix &r_A, size_t n_off_diagonal_num,
		const size_t *p_constraints, size_t UNUSED(n_constraints_size)); // throws(std::bad_alloc)

	/**
	 *	@brief calculates elementwise ordering by expanding blockwise ordering on a block matrix
	 *
	 *	@param[in] r_A is the block matrix (must be symmetric)
	 *	@param[in] b_inverse is inverse ordering flag which chooses the source
	 *		ordering to expand; if set, it is the last ordering produced by a call
	 *		to p_InvertOrdering(), if not set it is the last ordering produced
	 *		by a call to p_BlockOrdering().
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The b_inverse flag must not be set (there is a bug in expanding inverse ordering directly).
	 *	@note The buffer for elementwise ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of CMatrixOrdering
	 *		if multiple elementwise orderings need to exist at the same time).
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_ExpandBlockOrdering(const CUberBlockMatrix &r_A, bool b_inverse); // throws(std::bad_alloc)

	/**
	 *	@brief calculates ordering for an elementwise sparse matrix using CAMD
	 *
	 *	@param[in] p_A is the sparse matrix (must be symmetric)
	 *	@param[in] p_constraints is the ordering constraint vector (can be 0)
	 *	@param[in] n_constraints_size is size og the ordering constraint vector
	 *		(must match the number of p_A columns, except it is ignored
	 *		if p_constraints is 0)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of CMatrixOrdering
	 *		if multiple block orderings need to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note In case the matrix is no longer needed, it is better to call
	 *		p_DestructiveOrdering() as this function needs to make a copy of the matrix.
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_Ordering(const cs *p_A, const size_t *p_constraints,
		size_t UNUSED(n_constraints_size)); // throws(std::bad_alloc)

	/**
	 *	@brief calculates ordering for an elementwise sparse matrix using CAMD
	 *
	 *	@param[in] p_A is the sparse matrix (must be symmetric, the contents will
	 *		be dammaged by this call)
	 *	@param[in] p_constraints is the ordering constraint vector (can be 0)
	 *	@param[in] n_constraints_size is size og the ordering constraint vector
	 *		(must match the number of p_A columns, except it is ignored
	 *		if p_constraints is 0)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of CMatrixOrdering
	 *		if multiple block orderings need to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note In case the matrix is needed to remain untouched, it is needed to call
	 *		p_Ordering() as this function dammages the contents of the matrix.
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_DestructiveOrdering(cs *p_A, const size_t *p_constraints,
		size_t UNUSED(n_constraints_size)); // throws(std::bad_alloc)

	/**
	 *	@brief calculates ordering for an elementwise sparse matrix using CAMD
	 *
	 *	@param[in] p_A is the sparse matrix (must be symmetric, the contents will
	 *		be dammaged by this call)
	 *	@param[in] n_off_diagonal_num is the number of off-diagonal elements
	 *		(equals number of edges in binary graphs or the sum of edge ranks
	 *		minus the number of vertices in hypergraphs)
	 *	@param[in] p_constraints is the ordering constraint vector (can be 0)
	 *	@param[in] n_constraints_size is size og the ordering constraint vector
	 *		(must match the number of p_A columns, except it is ignored
	 *		if p_constraints is 0)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of CMatrixOrdering
	 *		if multiple block orderings need to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note In case the matrix is needed to remain untouched, it is needed to call
	 *		p_Ordering() as this function dammages the contents of the matrix.
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_HybridDestructiveOrdering(cs *p_A, size_t n_off_diagonal_num,
		const size_t *p_constraints, size_t UNUSED(n_constraints_size)); // throws(std::bad_alloc)

	/**
	 *	@brief extends an existing block ordering with identity ordering
	 *
	 *	@param[in] n_new_size is required size of the new ordering (must be greater than old size)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of CMatrixOrdering
	 *		if multiple block orderings need to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note This function throws std::bad_alloc.
	 */
	const size_t *p_ExtendBlockOrdering_with_Identity(size_t n_new_size); // throws(std::bad_alloc)
};

/**
 *	@brief constrained matrix ordering calculator (CAMD wrapper)
 *
 *	This is a convenience wrapper for CLastElementOrderingConstraint and CMatrixOrdering,
 *	calculating the constraint automatically for the caller. The constraint forces the
 *	last matrix block to still be the last, even after the ordering.
 */
class CLastElementConstrainedMatrixOrdering {
protected:
	CLastElementOrderingConstraint m_c; /**< @brief ordering constriant object */
	CMatrixOrdering m_o; /**< @brief matrix ordering (CAMD wrapper) */

public:
	/**
	 *	@brief calculates blockwise ordering in a matrix, using the CAMD library
	 *
	 *	@param[in] r_A is the block matrix (must be symmetric)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of
	 *		CLastElementConstrainedMatrixOrdering if multiple block orderings need
	 *		to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note This function throws std::bad_alloc.
	 */
	inline const size_t *p_BlockOrdering(const CUberBlockMatrix &r_A) // throws(std::bad_alloc)
	{
		size_t n_block_num = r_A.n_BlockColumn_Num();
		return m_o.p_BlockOrdering(r_A, m_c.p_Get(n_block_num), n_block_num);
	}

	/**
	 *	@brief calculates blockwise ordering in a matrix, using the CAMD library
	 *
	 *	@param[in] r_A is the block matrix (must be symmetric)
	 *	@param[in] n_off_diagonal_num is the number of off-diagonal elements
	 *		(equals number of edges in binary graphs or the sum of edge ranks
	 *		minus the number of vertices in hypergraphs)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of
	 *		CLastElementConstrainedMatrixOrdering if multiple block orderings need
	 *		to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note This function throws std::bad_alloc.
	 */
	inline const size_t *p_HybridBlockOrdering(const CUberBlockMatrix &r_A, size_t n_off_diagonal_num) // throws(std::bad_alloc)
	{
		size_t n_block_num = r_A.n_BlockColumn_Num();
		return m_o.p_HybridBlockOrdering(r_A, n_off_diagonal_num, m_c.p_Get(n_block_num), n_block_num);
	}

	/**
	 *	@brief extends an existing block ordering with identity ordering
	 *
	 *	@param[in] n_new_size is required size of the new ordering (must be greater than old size)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for block ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of
	 *		CLastElementConstrainedMatrixOrdering if multiple block orderings need
	 *		to exist at the same time).
	 *	@note The same buffer is used by p_BlockOrdering(), p_HybridBlockOrdering,
	 *		p_Ordering(), p_DestructiveOrdering() p_HybridDestructiveOrdering()
	 *		and p_ExtendBlockOrdering_with_Identity(). A call to any of those
	 *		functions invalidates previous results of all of those functions.
	 *	@note This function throws std::bad_alloc.
	 */
	inline const size_t *p_ExtendBlockOrdering_with_Identity(size_t n_new_size) // throws(std::bad_alloc)
	{
		return m_o.p_ExtendBlockOrdering_with_Identity(n_new_size);
	}

	/**
	 *	@brief calculates elementwise ordering by expanding blockwise ordering on a block matrix
	 *
	 *	This takes blockwise ordering calculated by the last call to p_BlockOrdering()
	 *	or to p_ExtendBlockOrdering_with_Identity() as input. Calling this function before
	 *	calling one of those is not permitted.
	 *
	 *	@param[in] r_A is the block matrix (must be symmetric)
	 *
	 *	@return Returns const pointer to the ordering vector (not to be deleted).
	 *
	 *	@note The buffer for elementwise ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of
	 *		CLastElementConstrainedMatrixOrdering if multiple block orderings need
	 *		to exist at the same time).
	 *	@note This function throws std::bad_alloc.
	 */
	inline const size_t *p_ExpandBlockOrdering(const CUberBlockMatrix &r_A) // throws(std::bad_alloc)
	{
		return m_o.p_ExpandBlockOrdering(r_A, false);
	}

	/**
	 *	@brief calculates inverse ordering while maintaining and reusing storage
	 *
	 *	@param[in] p_ordering is an ordering to be inverted
	 *	@param[in] n_ordering_size is size of the given ordering
	 *
	 *	@return Returns const pointer to the inverse of the given ordering (not to be deleted).
	 *
	 *	@note The buffer for inverse ordering is reused in the next function call,
	 *		invalidating the previous result (create more instances of
	 *		CLastElementConstrainedMatrixOrdering if multiple block orderings need
	 *		to exist at the same time).
	 *	@note Due to the buffer reuse, this can not invert ordering, inverted
	 *		by the same object (calling the same function on the result of the previous
	 *		call, on the same object). This doesn't apply to the results of the other functions.
	 *	@note This function throws std::bad_alloc.
	 */
	inline const size_t *p_InvertOrdering(const size_t *p_ordering, size_t n_ordering_size) // throws(std::bad_alloc)
	{
		return m_o.p_InvertOrdering(p_ordering, n_ordering_size);
	}
};

#endif // __MATRIX_ORDERING_UTILS_INCLUDED
