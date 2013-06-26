/*
								+-----------------------------------+
								|                                   |
								| *** Matrix ordering utilities *** |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2013  |
								|                                   |
								|         OrderingMagic.cpp         |
								|                                   |
								+-----------------------------------+
*/

/**
 *	@file src/slam/OrderingMagic.cpp
 *	@brief matrix ordering utilities
 *	@author -tHE SWINe-
 *	@date 2013-02-07
 */

#include <math.h>
#include <stdio.h>
#include "cholmod/AMD/amd.h"
#include "cholmod/CAMD/camd.h" // a bit unfortunate, now there are two camd.h, one in the global namespace and the other in __camd_internal
#include "slam/Timer.h"

#include "slam/OrderingMagic.h"

#ifdef __MATRIX_ORDERING_USE_AMD1

/**
 *	@brief hides functions from camd_internal.h which are not normally user-callable
 */
namespace __camd_internal {

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
#ifndef DLONG
#define DLONG
#endif // DLONG
#ifdef NLONG
#undef NLONG
#endif // NLONG
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
// makes sure that DLONG is defined in x64 (for CAMD)

extern "C" {
#include "cholmod/CAMD/camd_internal.h" // required by some other blah

#if defined(NDEBUG) && defined(_DEBUG)
#undef NDEBUG
#endif // NDEBUG && _DEBUG
// camd_internal.h defines NDEBUG if it is not defined alteady; might confuse programs

}
// t_odo - put this inside another file, this is inconvenient

}

/**
 *	@brief hides functions from amd_internal.h which are not normally user-callable
 */
namespace __amd_internal {

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
#ifndef DLONG
#define DLONG
#endif // DLONG
#ifdef NLONG
#undef NLONG
#endif // NLONG
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
// makes sure that DLONG is defined in x64 (for AMD)

extern "C" {

#ifdef PRINTF
#undef PRINTF
#endif // !PRINTF
// avoid macro redefinition in amd_internal.h

#include "cholmod/AMD/amd_internal.h" // required by some other blah

#if defined(NDEBUG) && defined(_DEBUG)
#undef NDEBUG
#endif // NDEBUG && _DEBUG
// amd_internal.h defines NDEBUG if it is not defined alteady; might confuse programs

}

}

#endif // __MATRIX_ORDERING_USE_AMD1

/*
 *								=== global ===
 */

/**
 *	@brief debugging predicate
 *	@tparam _Ty is an integral type
 *	@param[in] n_x is a number
 *	@return Returns true if the given number is nonzero, otherwise returns false.
 */
template <class _Ty>
static inline bool b_IsNonzero(_Ty n_x) // just for debugging
{
	return n_x != 0;
}

/*
 *								=== ~global ===
 */

/*
 *								=== CLastElementOrderingConstraint ===
 */

const size_t *CLastElementOrderingConstraint::p_Get(size_t n_size) // throw(std::bad_alloc)
{
#if defined(__MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT)
	const size_t n_groups_size = 14; // todo move this to the template arg or something
#endif // __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT

	if(!m_constraints.empty()) {
		if(m_constraints.capacity() < n_size)
			m_constraints.clear(); // memset() cheaper than memcpy()
		else {
#if defined(__MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT)
			_ASSERTE(m_constraints.back() == ((m_constraints.size() == 1)? 0 :
				((m_constraints.size() > n_groups_size)? 2 : 1))); // t_odo copy to L as well
			size_t n = m_constraints.size();
			for(size_t i = n - std::min(n, n_groups_size); i < n; ++ i)
				m_constraints[i] = 0;
#else // __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT
			_ASSERTE(m_constraints.back() == 1);
			m_constraints.back() = 0;
#endif // __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT
			// there is a 'one' at the end,
		}
	}
	m_constraints.resize(n_size, 0); // resize to contain all zeroes
	_ASSERTE(std::find_if(m_constraints.begin(), m_constraints.end(),
		b_IsNonzero<size_t>) == m_constraints.end()); // makes sure it's all zeroes
	// maintain vector of zeroes (and do it efficiently)

#if defined(__MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT)
	if(n_size <= 2)
		m_constraints[n_size - 1] = (n_size > 1)? 1 : 0; // either "0" or "0, 1"
	else if(n_size > n_groups_size) {
		for(size_t i = n_size - n_groups_size; i < n_size; ++ i)
			m_constraints[i] = 1;
		m_constraints[n_size - 1] = 2; // "0, 0, 0, ... , 0, 1, 1, 1, ... , 1, 2"
	} else {
		for(size_t i = 0; i < n_size; ++ i)
			m_constraints[i] = 0;
		m_constraints[n_size - 1] = 1; // "0, 0, 0, ... , 0, 1"
	}
	// make two groups, one forcing the last pose to be the last and the other forcing several other recent poses to still be at the end
#else // __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT
	m_constraints[n_size - 1] = 1; // impose no constraint on ordering // t_odo remove me
	// put a one at the end
#endif // __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT

	return &m_constraints[0];
	// return pointer for use by camd / ccolamd
}

/*
 *								=== ~CLastElementOrderingConstraint ===
 */

/*
 *								=== CFirstLastElementOrderingConstraint ===
 */

const size_t *CFirstLastElementOrderingConstraint::p_Get(size_t n_size) // throw(std::bad_alloc)
{
	const size_t n_groups_size = 14; // todo move this to the template arg or something

	if(m_constraints.capacity() < n_size) {
		m_constraints.clear();
		m_constraints.reserve(std::max(n_size, 2 * m_constraints.capacity()));
	}
	m_constraints.resize(n_size);

	if(n_size == 1)
		m_constraints[0] = 0;
	else if(n_size == 2) {
		m_constraints[0] = 0;
		m_constraints[1] = 1;
	} else if(n_size <= 2 + n_groups_size) {
		m_constraints[0] = 0;
		_ASSERTE(n_size - 1 > 1); // fill at least one
		for(size_t i = 1; i < n_size - 1; ++ i)
			m_constraints[i] = 1;
		m_constraints[n_size - 1] = 2;
	} else {
		m_constraints[0] = 0;
		_ASSERTE(n_size - (n_groups_size + 1) > 1); // fill at least one
		for(size_t i = 1; i < n_size - (n_groups_size + 1); ++ i)
			m_constraints[i] = 1;
		_ASSERTE(n_size - (n_groups_size + 1) < n_size - 1); // fill at least one
		for(size_t i = n_size - (n_groups_size + 1); i < n_size - 1; ++ i)
			m_constraints[i] = 2;
		m_constraints[n_size - 1] = 3;
	}
	// fixes the first elem as well; very complicated (don't even attempting to do lazy updates on that)

	return &m_constraints[0];
	// return pointer for use by camd / ccolamd
}

/*
 *								=== ~CFirstLastElementOrderingConstraint ===
 */

/*
 *								=== CNFirst1LastElementOrderingConstraint ===
 */

const size_t *CNFirst1LastElementOrderingConstraint::p_Get(size_t n_size, size_t n_first_constraint_num) // throw(std::bad_alloc)
{
	const size_t n_groups_size = 14; // 14 was the best on VP, molson35 and kill

	if(m_constraints.capacity() < n_size) {
		m_constraints.clear();
		m_constraints.reserve(std::max(n_size, 2 * m_constraints.capacity()));
	}
	m_constraints.resize(n_size);

	_ASSERTE(n_size > 0 && n_first_constraint_num > 0);
	if(n_size == 1)
		m_constraints[0] = 0;
	else if(n_size == 2) {
		m_constraints[0] = 0;
		m_constraints[1] = 1;
	} else if(n_size <= n_first_constraint_num + 1) { // onle the constraints and the last one
		_ASSERTE(n_size - 1 > 0); // fill at least one
		for(size_t i = 0; i < n_size - 1; ++ i)
			m_constraints[i] = i; // constrain all
		m_constraints[n_size - 1] = n_size - 1;
	} else if(n_size <= n_first_constraint_num + 1 + n_groups_size) { // only the constraints, the last group and the last one
		_ASSERTE(n_first_constraint_num > 0); // fill at least one
		for(size_t i = 0; i < n_first_constraint_num; ++ i)
			m_constraints[i] = i; // constrain all
		_ASSERTE(n_first_constraint_num < n_size - 1); // fill at least one
		_ASSERTE(n_first_constraint_num + n_groups_size >= n_size - 1); // fill at most n_groups_size
		for(size_t i = n_first_constraint_num; i < n_size - 1; ++ i)
			m_constraints[i] = n_first_constraint_num; // the first last group
		m_constraints[n_size - 1] = n_first_constraint_num + 1; // the second last group
	} else { // the constraints, free elements, last group, last one
		_ASSERTE(n_first_constraint_num > 0); // fill at least one
		for(size_t i = 0; i < n_first_constraint_num; ++ i)
			m_constraints[i] = i; // constrain all
		_ASSERTE(n_first_constraint_num < n_size - (n_groups_size + 1)); // fill at least one
		for(size_t i = n_first_constraint_num; i < n_size - (n_groups_size + 1); ++ i)
			m_constraints[i] = n_first_constraint_num; // free elements
		_ASSERTE(n_size - (n_groups_size + 1) < n_size - 1); // fill at least one
		for(size_t i = n_size - (n_groups_size + 1); i < n_size - 1; ++ i)
			m_constraints[i] = n_first_constraint_num + 1; // last group
		m_constraints[n_size - 1] = n_first_constraint_num + 2; // last one
	}
	// fixes the first elem as well; very complicated (don't even attempting to do lazy updates on that)

	return &m_constraints[0];
	// return pointer for use by camd / ccolamd
}

/*
 *								=== ~CNFirst1LastElementOrderingConstraint ===
 */

/*
 *								=== CMatrixOrdering ===
 */

CMatrixOrdering::~CMatrixOrdering()
{
	if(m_p_block)
		cs_spfree(m_p_block);
}

const size_t *CMatrixOrdering::p_InvertOrdering(const size_t *p_ordering, size_t n_ordering_size) // throw(std::bad_alloc)
{
	_ASSERTE(m_ordering_invert.empty() || p_ordering != &m_ordering_invert[0]); // can't invert twice
	if(m_ordering_invert.capacity() < n_ordering_size) {
		m_ordering_invert.clear();
		m_ordering_invert.reserve(std::max(n_ordering_size, 2 * m_ordering_invert.capacity()));
	}
	m_ordering_invert.resize(n_ordering_size);
	for(size_t i = 0; i < n_ordering_size; ++ i)
		m_ordering_invert[p_ordering[i]] = i;
	return &m_ordering_invert[0];
}

const size_t *CMatrixOrdering::p_ExtendBlockOrdering_with_SubOrdering(size_t n_order_min,
	const size_t *p_sub_block_ordering, size_t n_sub_block_size) // throw(std::bad_alloc)
{
	_ASSERTE(n_order_min + n_sub_block_size >= m_ordering.size());
	_ASSERTE(n_order_min <= m_ordering.size());
	if(n_order_min < m_ordering.size()) { // probably all cases in SLAM (otherwise would have a disconnected graph)
#if 1
		// this is working well and it works on direct ordering

		if(m_ordering_invert.capacity() < n_sub_block_size) {
			m_ordering_invert.clear();
			m_ordering_invert.reserve(std::max(n_sub_block_size, 2 * m_ordering_invert.capacity()));
		}
		m_ordering_invert.resize(n_sub_block_size);
		size_t *p_temp = &m_ordering_invert[0];
		// will use the inverse ordering as temporary ordering (can't do that inplace)

		size_t n_old_size = m_ordering.size();
		for(size_t i = 0; i < n_sub_block_size; ++ i) {
			size_t n_order = p_sub_block_ordering[i] + n_order_min;
			_ASSERTE(n_order >= n_order_min);
			p_temp[i] = (n_order < n_old_size)? m_ordering[n_order] : n_order;
		}
		// permute the end of the old ordering using the new ordering, extend with identity where needed

		m_ordering.erase(m_ordering.begin() + n_order_min, m_ordering.end()); // !!
		m_ordering.insert(m_ordering.end(), p_temp, p_temp + n_sub_block_size);
		// append with the new ordering
#elif 0
		// buggy

		size_t *p_inv_orig = &m_ordering_invert[n_sub_block_size];
		// will use the inverse ordering as temporary ordering (can't do that inplace)

		size_t n_old_size = m_ordering.size();
		for(size_t i = 0; i < n_order_min + n_sub_block_size; ++ i) {
			size_t n_order = (i < n_old_size)? m_ordering[i] : i;
			// extend the new ordering with identity

			if(n_order >= n_order_min) {
				_ASSERTE(n_order < n_order_min + n_sub_block_size);
				p_inv_orig[n_order - n_order_min] = i;
			}
		}
		// invert the original ordering to a temp array

		for(size_t i = 0; i < n_sub_block_size; ++ i) {
			size_t n_order = p_sub_block_ordering[i] + n_order_min;
			_ASSERTE(n_order >= n_order_min);
			p_temp[i] = p_inv_orig[n_order - n_order_min];//(n_order < n_old_size)? m_ordering[n_order] : n_order; // permute in inverse
		}
		// permute the end of the old ordering using the new ordering, extend with identity where needed

		m_ordering.erase(m_ordering.begin() + n_order_min, m_ordering.end()); // !!
		m_ordering.insert(m_ordering.end(), p_temp, p_temp + n_sub_block_size);
#elif 0
		// buggy

		_ASSERTE(m_ordering_invert.size() >= m_ordering.size()); // make sure the inverse is up-to-date
		// this works with inverse since we are ordering on an already ordered matrix

		std::vector<size_t> p_temp(n_sub_block_size);
		size_t n_old_size = m_ordering.size();
		for(size_t i = 0; i < n_sub_block_size; ++ i) {
			size_t n_order = p_sub_block_ordering[i] + n_order_min;
			_ASSERTE(n_order >= n_order_min);
			p_temp[i] = (n_order < n_old_size)? m_ordering_invert[n_order] : n_order;
		}
		// permute the end of inverse of the old ordering using the new ordering, extend with identity where needed

		m_ordering_invert.erase(m_ordering_invert.begin() + n_order_min, m_ordering_invert.end()); // !!
		m_ordering_invert.insert(m_ordering_invert.end(), p_temp.begin(), p_temp.begin() + n_sub_block_size);
		// append with the new ordering

		size_t n_new_size = m_ordering_invert.size();
		if(m_ordering.capacity() < n_new_size) {
			m_ordering.clear();
			m_ordering.resize(std::max(n_new_size, 2 * m_ordering.capacity()));
		}
		for(size_t i = 0; i < n_new_size; ++ i)
			m_ordering[m_ordering_invert[i]] = i;
		// inverse back
#else
		// this is working well and it works on inverse ordering

		size_t n_old_size = m_ordering.size();
		size_t n_new_size = n_order_min + n_sub_block_size;
		m_ordering.resize(n_new_size);
		for(size_t i = n_old_size; i < n_new_size; ++ i)
			m_ordering[i] = i;
		// extend with identity

		p_InvertOrdering(&m_ordering[0], m_ordering.size());

		for(size_t i = 0; i < n_new_size; ++ i) {
			size_t n_order = m_ordering_invert[i];
			m_ordering[(n_order < n_order_min)? n_order :
				p_sub_block_ordering[n_order - n_order_min] + n_order_min] = i;
		}
		// reorder the order
#endif
	} else {
		_ASSERTE(n_order_min == m_ordering.size());
		m_ordering.resize(n_order_min + n_sub_block_size);
		for(size_t i = 0; i < n_sub_block_size; ++ i)
			m_ordering[n_order_min + i] = p_sub_block_ordering[i] + n_order_min;
		// just resize and append the ordering with the new numbers
	}

	return &m_ordering[0];
}

const size_t *CMatrixOrdering::p_BlockOrdering(const CUberBlockMatrix &r_A, const size_t *p_constraints,
	size_t UNUSED(n_constraints_size), bool b_need_inverse) // throw(std::bad_alloc)
{
	_ASSERTE(r_A.b_BlockSquare());
	_ASSERTE(!p_constraints || n_constraints_size == r_A.n_BlockColumn_Num());

	if(!(m_p_block = r_A.p_BlockStructure_to_Sparse(m_p_block)))
		throw std::bad_alloc(); // rethrow
	// get block structure of A

	return p_DestructiveOrdering(m_p_block, p_constraints,
		n_constraints_size, b_need_inverse);
	// calculate ordering by block, can destroy the block matrix in the process
}

const size_t *CMatrixOrdering::p_BlockOrdering_MiniSkirt(const CUberBlockMatrix &r_A,
	size_t n_matrix_cut, size_t n_matrix_diag, const size_t *p_constraints,
	size_t UNUSED(n_constraints_size), bool b_need_inverse) // throw(std::bad_alloc)
{
	_ASSERTE(r_A.b_BlockSquare());
	_ASSERTE(!p_constraints || n_constraints_size == r_A.n_BlockColumn_Num() - n_matrix_cut);

	if(!(m_p_block = r_A.p_BlockStructure_to_Sparse_Apron(n_matrix_cut, n_matrix_diag, m_p_block)))
		throw std::bad_alloc(); // rethrow
	// get block structure of A

	return p_DestructiveOrdering(m_p_block, p_constraints,
		n_constraints_size, b_need_inverse);
	// calculate ordering by block, can destroy the block matrix in the process
}

/**
 *	@brief CRC calculation template
 *
 *	@param TScalar is CRC scalar type (uint16_t for crc16, uint32_t for crc32)
 *	@param n_polynom is CRC polynom
 *	@param n_start is CRC starting value
 *	@param n_final_xor is final reflection value (usually all 0's or all 1's)
 *
 *	@note It is possible to use predefined specializations CCrc_16 and CCrc_32,
 *		implementing standard CRC16 and CRC32 algorithms, respectively.
 */
template<class TScalar, const TScalar n_polynom, const TScalar n_start, const TScalar n_final_xor>
class CCrc {
public:
	/**
	 *	@brief gets CRC starting value
	 *	@return Returns starting value of CRC.
	 */
	static inline TScalar n_Start()
	{
		return n_start;
	}

	/**
	 *	@brief calculates final xor
	 *
	 *	@param[in] n_prev_crc is either n_Start() (CRC of empty buffer),
	 *		or result of previous call to n_Crc() (when calculating CRC of one or more
	 *		concatenated buffers)
	 *
	 *	@return Returns final value of CRC.
	 */
	static inline TScalar n_Finalize(TScalar n_prev_crc)
	{
		return n_prev_crc ^ n_final_xor;
	}

	/**
	 *	@brief calculates CRC on a stream of data
	 *
	 *	Calculates CRC of n_size bytes from buffer p_data.
	 *
	 *	@param[in] n_size is size of input data, in bytes
	 *	@param[in] p_data is input data buffer, contains n_size bytes
	 *	@param[in] n_prev_crc is either n_Start() (starting new CRC calculation),
	 *		or result of previous call to n_Crc() (when calculating CRC of more
	 *		concatenated buffers)
	 *
	 *	@return Returns CRC value, no data reflection, no final xor.
	 */
	static TScalar n_Crc(size_t n_size, const void *p_data, TScalar n_prev_crc)
	{
		const uint8_t *_p_data = (const uint8_t*)p_data;

		static TScalar p_crc_table[256];
		static bool b_crc_table_ready = false;
		if(!b_crc_table_ready) {
			b_crc_table_ready = true;
			for(int i = 0; i < 256; ++ i) {
				TScalar n_root = TScalar(i);
				for(int j = 0; j < 8; ++ j) {
					if(n_root & 1)
						n_root = (n_root >> 1) ^ n_polynom;
					else
						n_root >>= 1;
				}
				p_crc_table[i] = n_root;
			}
		}
		// prepare the table

		TScalar n_crc = n_prev_crc;
		for(const uint8_t *p_end = _p_data + n_size; _p_data != p_end; ++ _p_data)
			n_crc = p_crc_table[(n_crc ^ *_p_data) & 0xff] ^ (n_crc >> 8);
		// calculate CRC

		return n_crc;
	}

	/**
	 *	@brief calculates CRC
	 *
	 *	Calculates CRC of n_size bytes from buffer p_data.
	 *
	 *	@param[in] n_size is size of input data, in bytes
	 *	@param[in] p_data is input data buffer, contains n_size bytes
	 *
	 *	@return Returns the final CRC value.
	 */
	static inline TScalar n_Crc(size_t n_size, const void *p_data)
	{
		return n_Finalize(n_Crc(n_size, p_data, n_Start()));
	}
};

typedef CCrc<uint16_t, 0xA001, 0xffff, 0xffff> CCrc_16; /**< CRC-16, as defined by ANSI (X3.28) */
typedef CCrc<uint32_t, 0xedb88320U, 0xffffffffU, 0xffffffffU> CCrc_32; /**< CRC-32, as defined by POSIX */

const size_t *CMatrixOrdering::p_BlockOrdering(const CUberBlockMatrix &r_A, bool b_need_inverse) // throw(std::bad_alloc)
{
	_ASSERTE(r_A.b_BlockSquare());

	/*CTimer timer;
	double f_start = timer.f_Time();*/

	if(!(m_p_block = r_A.p_BlockStructure_to_Sparse(m_p_block)))
		throw std::bad_alloc(); // rethrow
	// get block structure of A

	_ASSERTE(r_A.b_SymmetricLayout()); // well, it should be
	// it should also be upper-triangular

	//const size_t *p_result = p_Ordering(m_p_block); // AMD, not CAMD
	// can do with any kind of matrix

	const size_t *p_result = p_DestructiveOrdering(m_p_block, b_need_inverse); // AMD, not CAMD
	// only upper-diagonal, only well sorted (guaranteed by p_BlockStructure_to_Sparse())

	/*double f_time = timer.f_Time() - f_start;

#if 1
	uint32_t n_ordering_crc32 = CCrc_32::n_Finalize(CCrc_32::n_Crc(r_A.n_BlockColumn_Num() *
		sizeof(size_t), p_result));
	printf("calculated ordering of size " PRIsize ", crc32: 0x%08x\n",
		r_A.n_BlockColumn_Num(), n_ordering_crc32);
	printf("it took %f sec\n", f_time);
#endif // 1*/
	// debug the orderings generated

	return p_result;
	// calculate ordering by block, can destroy the block matrix in the process
}

#if 0
const size_t *CMatrixOrdering::p_HybridBlockOrdering(const CUberBlockMatrix &r_A, size_t n_off_diagonal_num,
	const size_t *p_constraints, size_t UNUSED(n_constraints_size)) // throw(std::bad_alloc)
{
	_ASSERTE(r_A.b_BlockSquare());
	_ASSERTE(!p_constraints || n_constraints_size == r_A.n_BlockColumn_Num());

	if(!(m_p_block = r_A.p_BlockStructure_to_Sparse(m_p_block)))
		throw std::bad_alloc(); // rethrow
	// get block structure of A

	return p_HybridDestructiveOrdering(m_p_block, n_off_diagonal_num, p_constraints, n_constraints_size);
	// calculate ordering by block, can destroy the block matrix in the process
}
#endif // 0

const size_t *CMatrixOrdering::p_ExpandBlockOrdering(const CUberBlockMatrix &r_A, bool b_inverse) // throw(std::bad_alloc)
{
	_ASSERTE(r_A.b_BlockSquare());
	const size_t n_column_num = r_A.n_Column_Num();
	const size_t n_column_block_num = r_A.n_BlockColumn_Num();

	_ASSERTE(m_ordering.size() == n_column_block_num);
	// make sure the ordering can be possibly applied to the matrix

	_ASSERTE(!b_inverse);
	// can't expand inverse ordering, the r_A.n_BlockColumn_Base(n_order)
	// corresponds to a different block in the inverse ordering, than it did in the original ordering
	// r_A would have to be the permutation of r_A, but that makes this function dangerous
	// just assume we don't need the inverse ordering at the moment

	if(m_ordering_expand.size() < n_column_num) {
		m_ordering_expand.clear();
		m_ordering_expand.resize(std::max(n_column_num, 2 * m_ordering_expand.capacity()));
	}
	size_t n_scalar_offset = 0;
	const size_t *p_order_ptr = (b_inverse)? &m_ordering_invert[0] : &m_ordering[0];
	for(size_t i = 0; i < n_column_block_num; ++ i, ++ p_order_ptr) {
		const size_t n_order = *p_order_ptr;
		size_t n_block_base = r_A.n_BlockColumn_Base(n_order); // wonder what would happen if it was m_ordering[n_order] or something like that ... would that get a correct column in the inverse ordering?
		size_t n_block_width = r_A.n_BlockColumn_Column_Num(n_order);
		for(size_t j = 0; j < n_block_width; ++ j, ++ n_scalar_offset, ++ n_block_base)
			m_ordering_expand[n_scalar_offset] = n_block_base;
	}
	_ASSERTE(n_scalar_offset == n_column_num);
	// blow up the permutation from block level to scalar level

	return &m_ordering_expand[0];
}

const size_t *CMatrixOrdering::p_Ordering(const cs *p_A, const size_t *p_constraints,
	size_t UNUSED(n_constraints_size)) // throw(std::bad_alloc)
{
	_ASSERTE(!p_constraints || n_constraints_size == p_A->n);

	if(m_ordering.capacity() < size_t(p_A->n)) {
		m_ordering.clear(); // avoid copying data in resizing
		m_ordering.reserve(std::max(size_t(p_A->n), 2 * m_ordering.capacity()));
	}
	m_ordering.resize(p_A->n);
	// maintain an ordering array

	_ASSERTE(sizeof(SuiteSparse_long) == sizeof(size_t));
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	switch(camd_l_order(p_A->n, p_A->p, p_A->i, (SuiteSparse_long*)&m_ordering[0],
	   0, 0, (SuiteSparse_long*)p_constraints)) {
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	switch(camd_order(p_A->n, p_A->p, p_A->i, (int*)&m_ordering[0], 0, 0, (int*)p_constraints)) {
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	case CAMD_OK:
	case CAMD_OK_BUT_JUMBLED:
		break;
	case CAMD_OUT_OF_MEMORY:
		throw std::bad_alloc();
	case CAMD_INVALID:
		return 0;
	}
	// call camd

	return &m_ordering[0];
	// return ordering
}

const size_t *CMatrixOrdering::p_Ordering(const cs *p_A) // throw(std::bad_alloc)
{
	if(m_ordering.capacity() < size_t(p_A->n)) {
		m_ordering.clear(); // avoid copying data in resizing
		m_ordering.reserve(std::max(size_t(p_A->n), 2 * m_ordering.capacity()));
	}
	m_ordering.resize(p_A->n);
	// maintain an ordering array

	_ASSERTE(sizeof(SuiteSparse_long) == sizeof(size_t));
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	switch(amd_l_order(p_A->n, p_A->p, p_A->i, (SuiteSparse_long*)&m_ordering[0], 0, 0)) {
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	switch(amd_order(p_A->n, p_A->p, p_A->i, (int*)&m_ordering[0], 0, 0)) {
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	case AMD_OK:
	case AMD_OK_BUT_JUMBLED:
		break;
	case AMD_OUT_OF_MEMORY:
		throw std::bad_alloc();
	case AMD_INVALID:
		return 0;
	}
	// call amd

	return &m_ordering[0];
	// return ordering
}

const size_t *CMatrixOrdering::p_DestructiveOrdering(cs *p_A, bool b_need_inverse) // throw(std::bad_alloc)
{
	_ASSERTE(p_A->n >= 0 && p_A->n <= SIZE_MAX);
	if(m_ordering.capacity() < size_t(p_A->n)) {
		m_ordering.clear(); // avoid copying data in resizing
		m_ordering.reserve(std::max(size_t(p_A->n), 2 * m_ordering.capacity()));
	}
	m_ordering.resize(p_A->n);
	// maintain an ordering array

	if(b_need_inverse) {
		if(m_ordering_invert.capacity() < size_t(p_A->n)) {
			m_ordering_invert.clear(); // avoid copying data in resizing
			m_ordering_invert.reserve(std::max(size_t(p_A->n), 2 * m_ordering_invert.capacity()));
		}
		m_ordering_invert.resize(p_A->n);
	}
	// allocate an extra inverse ordering array, if needed

	const size_t n = p_A->n, nz = p_A->p[p_A->n];
	if(n >= SIZE_MAX / sizeof(size_t) ||
	   nz >= SIZE_MAX / sizeof(size_t))
		throw std::bad_alloc(); // would get AMD_OUT_OF_MEMORY

#ifdef __MATRIX_ORDERING_CACHE_ALIGN
	const size_t n_al = n_Align_Up_POT(n, __MATRIX_ORDERING_CACHE_ALIGNMENT_BYTES / sizeof(size_t));
#else // __MATRIX_ORDERING_CACHE_ALIGN
	const size_t n_al = n;
#endif // __MATRIX_ORDERING_CACHE_ALIGN
	// to align, or not to align ...

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	_ASSERTE(amd_l_valid(n, n, p_A->p, p_A->i) == AMD_OK); // note that this fails on "jumbled" matrices (different behavior than p_Ordering())
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	_ASSERTE(amd_valid(n, n, p_A->p, p_A->i) == AMD_OK); // note that this fails on "jumbled" matrices (different behavior than p_Ordering())
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	// check if it's valid (debug only)

	_ASSERTE(sizeof(SuiteSparse_long) == sizeof(size_t));
	_ASSERTE(sizeof(size_t) > 1); // note 2 * n does not overflow if sizeof(size_t) > 1
	if(m_camd_workspace.capacity() < ((b_need_inverse)? n : n + n_al)) {
		m_camd_workspace.clear(); // avoid copying data if it needs to realloc
		m_camd_workspace.reserve(std::max((b_need_inverse)? n : n + n_al, m_camd_workspace.capacity() * 2));
	}
	m_camd_workspace.resize((b_need_inverse)? n : n + n_al);
	size_t *Len = &m_camd_workspace[0];
	size_t *Pinv = (b_need_inverse)? &m_ordering_invert[0] : &m_camd_workspace[n_al]; // this is actually inverse ordering, we need it most of the time, it is a waste to throw it away
	// alloc workspace

#ifdef __MATRIX_ORDERING_USE_AMD_AAT
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	size_t nzaat = __amd_internal::amd_l_aat(n, p_A->p, p_A->i,
		(SuiteSparse_long*)Len, (SuiteSparse_long*)&m_ordering[0], 0);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	size_t nzaat = __amd_internal::amd_aat(n, p_A->p, p_A->i,
		(int*)Len, (int*)&m_ordering[0], 0);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	_ASSERTE(nzaat == 2 * (nz - n)); // A is upper-triangular, nz of A+A^T (excluding diagonal) should equal 2 * (nz - n)
#else // __MATRIX_ORDERING_USE_AMD_AAT
	{
		const ptrdiff_t *Ap = p_A->p, *Ai = p_A->i; // antialiass
		/*memset(Len, 0, p_A->n * sizeof(size_t)); // memset Len to 0
		for(size_t i = 0, p1 = Ap[0]; i < n; ++ i) {
			size_t p2 = Ap[i + 1];
			_ASSERTE(p1 < p2); // no empty columns!
			_ASSERTE(Ai[p2 - 1] == i); // diagonal is here (skip)
			-- p2; // skip the diagonal
			_ASSERTE(p2 >= p1); // makes sure the next line will not overflow
			Len[i] += p2 - p1;
			for(; p1 < p2; ++ p1) {
				size_t j = Ai[p1];
				_ASSERTE(j < i); // upper-triangular and we skipped the diagonal before the loop
				++ Len[j]; // off-diagonal
			}
			++ p1;
			_ASSERTE(p1 == Ap[i + 1]); // this will soo fail on jumbled matrices (_DEBUG checks for them)
		}*/
		for(size_t i = 0, p1 = Ap[0]; i < n; ++ i) {
			size_t p2 = Ap[i + 1];
			_ASSERTE(p1 < p2); // no empty columns!
			_ASSERTE(Ai[p2 - 1] == i); // diagonal is here (skip)
			_ASSERTE(p2 > p1); // makes sure the next line will not underflow
			Len[i] = p2 - p1 - 2; // could actually underflow by 1, the diag entry will fix it
			p1 = p2;
			_ASSERTE(p1 == Ap[i + 1]); // this will soo fail on jumbled matrices (_DEBUG checks for them)
		}
		for(size_t i = 0; i < nz; ++ i)
			++ Len[Ai[i]];
		// this could be vectorized, essentially len = Ap+1 + Ap + Ai - 2, at least the first two could do with integer sse
		// alternately, memset could be ommited by splitting to one loop over Ap (assigning) and one over Ai (incrementing)
	}
	size_t nzaat = 2 * (nz - n); // we know this
#endif // __MATRIX_ORDERING_USE_AMD_AAT
	_ASSERTE((std::max(nz, n) - n <= nzaat) && (nzaat <= 2 * (size_t)nz));
	//printf("nzaat = " PRIsize "\n", nzaat); // 805348 for 100k, my custom code, amd_aat() gives the same
	// determine the symmetry and count off-diagonal nonzeros in A+A'

	size_t slen = nzaat;
	if(slen > SIZE_MAX - nzaat / 5) // check overflow
		throw std::bad_alloc(); // would get AMD_OUT_OF_MEMORY
	slen += nzaat / 5;
	if(n > SIZE_MAX / 8 || (n == SIZE_MAX && SIZE_MAX % 8) || slen > SIZE_MAX - 8 * n) // check overflow
		throw std::bad_alloc(); // would get AMD_OUT_OF_MEMORY
	slen += 8 * n_al;
	if(m_camd_workspace1.size() < slen) {
		m_camd_workspace1.clear(); // avoid copying data if it needs to realloc
		m_camd_workspace1.resize(std::max(slen, m_camd_workspace1.size() * 2));
	}
	_ASSERTE(m_camd_workspace1.size() >= slen);
	slen = m_camd_workspace1.size();
	// allocate workspace for matrix, elbow room (size n), and 7 (n + 1) vectors

	// todo finish with aligning, compare the actual orderings

#ifdef __MATRIX_ORDERING_USE_AMD1
	size_t *S = &m_camd_workspace1[0];
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	__amd_internal::amd_l1(n, p_A->p, p_A->i, (SuiteSparse_long*)&m_ordering[0],
		(SuiteSparse_long*)Pinv, (SuiteSparse_long*)Len, slen,
		(SuiteSparse_long*)S, 0, 0);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	__amd_internal::amd_1(n, p_A->p, p_A->i, (int*)&m_ordering[0], (int*)Pinv,
		(int*)Len, slen, (int*)S, 0, 0);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	// order the matrix
#else // __MATRIX_ORDERING_USE_AMD1
	{
		const ptrdiff_t *Ap = p_A->p, *Ai = p_A->i; // antialiass

		size_t iwlen = slen - 6 * n_al;
		size_t *Pe, *Nv, *Head, *Elen, *Degree, *W, *Iw;
		{
			size_t *S = &m_camd_workspace1[0];
			Pe = S;		S += n_al;
			Nv = S;		S += n_al;
			Head = S;	S += n_al;
			Elen = S;	S += n_al;
			Degree = S;	S += n_al;
			W = S;
			Iw = S + n_al;
		}
		// one array actually stores many smaller ones // t_odo play with cache here? could align n to 64B offset so that each array fills different cache entry (if 6bit cache off)

		size_t pfree = 0; // this is needed in amd_2()
		{
			*Pe = 0; // not really calculated, needs to be written here
			size_t *p_column_dest = Pe + 1; // calculate Pe inplace
#ifdef _DEBUG
			size_t *p_column_end = Head; // use Head as workspace (do not use Nv, the first item gets overwritten)
#endif // _DEBUG
			for(size_t j = 0; j < n; ++ j) {
				p_column_dest[j] = pfree; // pointer to the next entry to write in column j
				pfree += Len[j];
#ifdef _DEBUG
				p_column_end[j] = pfree; // (j+1)-th column pointer of A+A^T
#endif // _DEBUG
			} // t_odo could only accum p_column_dest, and employ that p_column_dest[i] == Pe[i + 1] at the end. just need to off by 1 and add the starting 0
			_ASSERTE(iwlen >= pfree + n);
			// initialize column pointers

			for(size_t k = 0; k < n; ++ k) {
				size_t p1 = Ap[k];
				size_t p2 = Ap[k + 1];

				_ASSERTE(p1 < p2); // no empty columns!
				_ASSERTE(Ai[p2 - 1] == k); // diagonal is here (skip)
				-- p2; // skip the diagonal
				_ASSERTE(p2 >= p1); // makes sure the next line will not overflow

				for(size_t p = p1; p < p2; ++ p) {
					size_t j = Ai[p];
					_ASSERTE(j >= 0 && j < n);
					_ASSERTE(j < k); // only above-diagonal elements
					// scan the upper triangular part of A

					_ASSERTE(p_column_dest[j] < p_column_end[j]);
					_ASSERTE(p_column_dest[k] < p_column_end[k]);
					Iw[p_column_dest[j] ++] = k;
					Iw[p_column_dest[k] ++] = j;
					// entry A(j, k) in the strictly upper triangular part
				}
				// construct one column of A+A^T
			}
			// builds A+A^T from a symmetric matrix (todo - make this a part of CUberBlockMatrix, avoid forming A, build A+A^T directly)

#ifdef _DEBUG
			for(size_t j = 0, n_cumsum = 0; j < n; ++ j) {
				_ASSERTE(Pe[j] == n_cumsum); // points at the beginning of the column
				n_cumsum += Len[j];
				_ASSERTE(p_column_dest[j] == n_cumsum); // points at the end of the column
			}
			// make sure all the columns are filled, according to Len
#endif // _DEBUG
		}

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
		amd_l2(n, (SuiteSparse_long*)Pe, (SuiteSparse_long*)Iw,
			(SuiteSparse_long*)Len, iwlen, pfree, (SuiteSparse_long*)Nv,
			(SuiteSparse_long*)Pinv, (SuiteSparse_long*)&m_ordering[0],
			(SuiteSparse_long*)Head, (SuiteSparse_long*)Elen,
			(SuiteSparse_long*)Degree, (SuiteSparse_long*)W, 0, 0);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
		amd_2(n, (int*)Pe, (int*)Iw, (int*)Len, iwlen, pfree, (int*)Nv, (int*)Pinv,
			(int*)&m_ordering[0], (int*)Head, (int*)Elen, (int*)Degree, (int*)W, 0, 0);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	}
	// original code taken from AMD by Tim Davis
#endif // __MATRIX_ORDERING_USE_AMD1

	if(b_need_inverse) {
		size_t *p_ordering = &m_ordering[0];
		for(size_t i = 0; i < n; ++ i)
			Pinv[p_ordering[i]] = i;
	}
	// AMD does not output the inverse ordering on output, at least not in this version // todo - or does it?

	// note that amd_1 constructs A^T | A (A is a binary matrix) first,
	// it could be optimized / calculated inplace / implemented directly in block matrix // t_odo

	return &m_ordering[0];
	// return ordering
}

const size_t *CMatrixOrdering::p_DestructiveOrdering(cs *p_A, const size_t *p_constraints,
	size_t UNUSED(n_constraints_size), bool b_need_inverse) // throw(std::bad_alloc)
{
	_ASSERTE(!p_constraints || n_constraints_size == p_A->n);

	_ASSERTE(p_A->n >= 0 && p_A->n <= SIZE_MAX);
	if(m_ordering.capacity() < size_t(p_A->n)) {
		m_ordering.clear(); // avoid copying data in resizing
		m_ordering.reserve(std::max(size_t(p_A->n), 2 * m_ordering.capacity()));
	}
	m_ordering.resize(p_A->n);
	// maintain an ordering array

	if(b_need_inverse) {
		if(m_ordering_invert.capacity() < size_t(p_A->n)) {
			m_ordering_invert.clear(); // avoid copying data in resizing
			m_ordering_invert.reserve(std::max(size_t(p_A->n), 2 * m_ordering_invert.capacity()));
		}
		m_ordering_invert.resize(p_A->n);
	}
	// allocate an extra inverse ordering array, if needed

	const size_t n = p_A->n, nz = p_A->p[p_A->n];
	if(n >= SIZE_MAX / sizeof(size_t) ||
	   nz >= SIZE_MAX / sizeof(size_t))
		throw std::bad_alloc(); // would get CAMD_OUT_OF_MEMORY

#ifdef __MATRIX_ORDERING_CACHE_ALIGN
	const size_t n_al = n_Align_Up_POT(n, __MATRIX_ORDERING_CACHE_ALIGNMENT_BYTES / sizeof(size_t));
#else // __MATRIX_ORDERING_CACHE_ALIGN
	const size_t n_al = n;
#endif // __MATRIX_ORDERING_CACHE_ALIGN
	// to align, or not to align ...

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	_ASSERTE(camd_l_valid(n, n, p_A->p, p_A->i) == CAMD_OK); // note that this fails on "jumbled" matrices (different behavior than p_Ordering())
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	_ASSERTE(camd_valid(n, n, p_A->p, p_A->i) == CAMD_OK); // note that this fails on "jumbled" matrices (different behavior than p_Ordering())
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	// check if it's valid (debug only)

	_ASSERTE(sizeof(SuiteSparse_long) == sizeof(size_t));
	_ASSERTE(sizeof(size_t) > 1); // note 2 * n does not overflow if sizeof(size_t) > 1
	if(m_camd_workspace.capacity() < ((b_need_inverse)? n : n + n_al)) {
		m_camd_workspace.clear(); // avoid copying data if it needs to realloc
		m_camd_workspace.reserve(std::max((b_need_inverse)? n : n + n_al, m_camd_workspace.capacity() * 2));
	}
	m_camd_workspace.resize((b_need_inverse)? n : n + n_al);
	size_t *Len = &m_camd_workspace[0];
	size_t *Pinv = (b_need_inverse)? &m_ordering_invert[0] : &m_camd_workspace[n_al]; // this is actually inverse ordering, we need it most of the time, it is a waste to throw it away
	// alloc workspace

#ifdef __MATRIX_ORDERING_USE_AMD_AAT
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	size_t nzaat = __camd_internal::camd_l_aat(n, p_A->p, p_A->i,
		(SuiteSparse_long*)Len, (SuiteSparse_long*)&m_ordering[0], 0);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	size_t nzaat = __camd_internal::camd_aat(n, p_A->p, p_A->i,
		(int*)Len, (int*)&m_ordering[0], 0);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	_ASSERTE(nzaat == 2 * (nz - n)); // A is upper-triangular, nz of A+A^T (excluding diagonal) should equal 2 * (nz - n)
#else // __MATRIX_ORDERING_USE_AMD_AAT
	{
		const ptrdiff_t *Ap = p_A->p, *Ai = p_A->i; // antialiass
		/*memset(Len, 0, p_A->n * sizeof(size_t)); // memset Len to 0
		for(size_t i = 0, p1 = Ap[0]; i < n; ++ i) {
			size_t p2 = Ap[i + 1];
			_ASSERTE(p1 < p2); // no empty columns!
			_ASSERTE(Ai[p2 - 1] == i); // diagonal is here (skip)
			-- p2; // skip the diagonal
			_ASSERTE(p2 >= p1); // makes sure the next line will not overflow
			Len[i] += p2 - p1;
			for(; p1 < p2; ++ p1) {
				size_t j = Ai[p1];
				_ASSERTE(j < i); // upper-triangular and we skipped the diagonal before the loop
				++ Len[j]; // off-diagonal
			}
			++ p1;
			_ASSERTE(p1 == Ap[i + 1]); // this will soo fail on jumbled matrices (_DEBUG checks for them)
		}*/
		for(size_t i = 0, p1 = Ap[0]; i < n; ++ i) {
			size_t p2 = Ap[i + 1];
			_ASSERTE(p1 < p2); // no empty columns!
			_ASSERTE(Ai[p2 - 1] == i); // diagonal is here (skip)
			_ASSERTE(p2 > p1); // makes sure the next line will not underflow
			Len[i] = p2 - p1 - 2; // could actually underflow by 1, the diag entry will fix it
			p1 = p2;
			_ASSERTE(p1 == Ap[i + 1]); // this will soo fail on jumbled matrices (_DEBUG checks for them)
		}
		for(size_t i = 0; i < nz; ++ i)
			++ Len[Ai[i]];
	}
	size_t nzaat = 2 * (nz - n); // we know this
#endif // __MATRIX_ORDERING_USE_AMD_AAT
	_ASSERTE((std::max(nz, n) - n <= nzaat) && (nzaat <= 2 * (size_t)nz));
	// determine the symmetry and count off-diagonal nonzeros in A+A'

	size_t slen = nzaat;
	if(slen > SIZE_MAX - nzaat / 5) // check overflow
		throw std::bad_alloc(); // would get CAMD_OUT_OF_MEMORY
	slen += nzaat / 5;
	size_t n1 = n_al + 1; // for sure does not overflow
	if(n1 > SIZE_MAX / 9 || (n1 == SIZE_MAX && SIZE_MAX % 9) || slen > SIZE_MAX - 9 * n1) // check overflow
		throw std::bad_alloc(); // would get CAMD_OUT_OF_MEMORY
	slen += 9 * n1; // note it is 1 bigger than it has to be // t_odo - check if we can save memory (probably not, CAMD requires 1 more than AMD, that will be it) // not a good idea, more memory runs faster
	if(m_camd_workspace1.size() < slen) {
		m_camd_workspace1.clear(); // avoid copying data if it needs to realloc
		m_camd_workspace1.resize(std::max(slen, m_camd_workspace1.size() * 2));
	}
	_ASSERTE(m_camd_workspace1.size() >= slen);
	slen = m_camd_workspace1.size(); // give it all we've got
	// allocate workspace for matrix, elbow room (size n), and 7 (n + 1) vectors

#ifdef __MATRIX_ORDERING_USE_AMD1
	size_t *S = &m_camd_workspace1[0];
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	__camd_internal::camd_l1(n, p_A->p, p_A->i, (SuiteSparse_long*)&m_ordering[0],
		(SuiteSparse_long*)Pinv, (SuiteSparse_long*)Len, slen,
		(SuiteSparse_long*)S, 0, 0, (SuiteSparse_long*)p_constraints);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	__camd_internal::camd_1(n, p_A->p, p_A->i, (int*)&m_ordering[0], (int*)Pinv,
		(int*)Len, slen, (int*)S, 0, 0, (int*)p_constraints);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	// order the matrix
#else // __MATRIX_ORDERING_USE_AMD1
	{
		const ptrdiff_t *Ap = p_A->p, *Ai = p_A->i; // antialiass

		size_t iwlen = slen - 7 * n_al - 2;
		size_t *Pe, *Nv, *Head, *Elen, *Degree, *W, *Iw, *BucketSet;
		{
			size_t *S = &m_camd_workspace1[0];
			Pe = S;		S += n_al;
			Nv = S;		S += n_al;
			Head = S;	S += n_al + 1; // note: 1 more than AMD
			Elen = S;	S += n_al;
			Degree = S;	S += n_al;
			W = S;		S += n_al + 1; // note: 1 more than AMD
			BucketSet = S;
			Iw = S + n_al;
		}
		// one array actually stores many smaller ones

		size_t pfree = 0; // this is needed in amd_2()
		{
			*Pe = 0; // not really calculated, needs to be written here
			size_t *p_column_dest = Pe + 1; // calculate Pe inplace
#ifdef _DEBUG
			size_t *p_column_end = Head; // use Head as workspace (do not use Nv, the first item gets overwritten)
#endif // _DEBUG
			for(size_t j = 0; j < n; ++ j) {
				p_column_dest[j] = pfree; // pointer to the next entry to write in column j
				pfree += Len[j];
#ifdef _DEBUG
				p_column_end[j] = pfree; // (j+1)-th column pointer of A+A^T
#endif // _DEBUG
			} // t_odo could only accum p_column_dest, and employ that p_column_dest[i] == Pe[i + 1] at the end. just need to off by 1 and add the starting 0
			_ASSERTE(iwlen >= pfree + n);
			// initialize column pointers

			for(size_t k = 0; k < n; ++ k) {
				size_t p1 = Ap[k];
				size_t p2 = Ap[k + 1];

				_ASSERTE(p1 < p2); // no empty columns!
				_ASSERTE(Ai[p2 - 1] == k); // diagonal is here (skip)
				-- p2; // skip the diagonal
				_ASSERTE(p2 >= p1); // makes sure the next line will not overflow

				for(size_t p = p1; p < p2; ++ p) {
					size_t j = Ai[p];
					_ASSERTE(j >= 0 && j < n);
					_ASSERTE(j < k); // only above-diagonal elements
					// scan the upper triangular part of A

					_ASSERTE(p_column_dest[j] < p_column_end[j]);
					_ASSERTE(p_column_dest[k] < p_column_end[k]);
					Iw[p_column_dest[j] ++] = k;
					Iw[p_column_dest[k] ++] = j;
					// entry A(j, k) in the strictly upper triangular part
				}
				// construct one column of A+A^T
			}
			// builds A+A^T from a symmetric matrix (todo - make this a part of CUberBlockMatrix, avoid forming A, build A+A^T directly)

#ifdef _DEBUG
			for(size_t j = 0, n_cumsum = 0; j < n; ++ j) {
				_ASSERTE(Pe[j] == n_cumsum); // points at the beginning of the column
				n_cumsum += Len[j];
				_ASSERTE(p_column_dest[j] == n_cumsum); // points at the end of the column
			}
			// make sure all the columns are filled, according to Len
#endif // _DEBUG
		}

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
		camd_l2(n, (SuiteSparse_long*)Pe, (SuiteSparse_long*)Iw,
			(SuiteSparse_long*)Len, iwlen, pfree, (SuiteSparse_long*)Nv,
			(SuiteSparse_long*)Pinv, (SuiteSparse_long*)&m_ordering[0],
			(SuiteSparse_long*)Head, (SuiteSparse_long*)Elen,
			(SuiteSparse_long*)Degree, (SuiteSparse_long*)W, 0, 0,
			(SuiteSparse_long*)p_constraints, (SuiteSparse_long*)BucketSet);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
		camd_2(n, (int*)Pe, (int*)Iw, (int*)Len, iwlen, pfree, (int*)Nv, (int*)Pinv,
			(int*)&m_ordering[0], (int*)Head, (int*)Elen, (int*)Degree, (int*)W, 0, 0,
			(int*)p_constraints, (int*)BucketSet);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	}
	// original code taken from AMD by Tim Davis
#endif // __MATRIX_ORDERING_USE_AMD1

	if(b_need_inverse) {
		size_t *p_ordering = &m_ordering[0];
		for(size_t i = 0; i < n; ++ i)
			Pinv[p_ordering[i]] = i;
	}
	// CAMD does not output the inverse ordering on output, at least not in this version

	// note that camd_1 constructs A^T | A (A is a binary matrix) first,
	// it could be optimized / calculated inplace / implemented directly in block matrix // t_odo

	return &m_ordering[0];
	// return ordering
}

#if 0
const size_t *CMatrixOrdering::p_HybridDestructiveOrdering(cs *p_A, size_t n_off_diagonal_num,
	const size_t *p_constraints, size_t UNUSED(n_constraints_size)) // throw(std::bad_alloc)
{
	_ASSERTE(!p_constraints || n_constraints_size == p_A->n);

	_ASSERTE(p_A->n >= 0 && p_A->n <= SIZE_MAX);
	if(m_ordering.capacity() < size_t(p_A->n))
		m_ordering.clear(); // avoid copying data in resizing
	m_ordering.resize(p_A->n);
	// maintain an ordering array

	const size_t n = p_A->n, nz = p_A->p[p_A->n];
	if(n >= SIZE_MAX / sizeof(size_t) ||
	   nz >= SIZE_MAX / sizeof(size_t))
		throw std::bad_alloc(); // would get CAMD_OUT_OF_MEMORY
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	_ASSERTE(camd_l_valid(n, n, p_A->p, p_A->i) == CAMD_OK); // note that this fails on "jumbled" matrices (different behavior than p_Ordering())
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	_ASSERTE(camd_valid(n, n, p_A->p, p_A->i) == CAMD_OK); // note that this fails on "jumbled" matrices (different behavior than p_Ordering())
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	// check if it's valid (debug only)

	_ASSERTE(sizeof(SuiteSparse_long) == sizeof(size_t));
	_ASSERTE(sizeof(size_t) > 1); // note 2 * n does not overflow if sizeof(size_t) > 1
	if(m_camd_workspace.capacity() < 2 * n) {
		m_camd_workspace.clear(); // avoid copying data if it needs to realloc
		m_camd_workspace.reserve(std::max(2 * n, m_camd_workspace.capacity() * 2));
	}
	m_camd_workspace.resize(2 * n);
	size_t *Len = &m_camd_workspace[0];
	size_t *Pinv = &m_camd_workspace[n];
	// alloc workspace

#if 1
	size_t nzaat = /*p_A->n +*/ 2 * n_off_diagonal_num; // calculate A+A^T nnz elements (A is upper triangular) // excluding diagonal!
	memset(Len, 0, p_A->n * sizeof(size_t)); // memset Len to 0
/*#if 0 //def _OPENMP
	if(p_A->p[p_A->n] > 1000) {
		const ptrdiff_t *p_i = p_A->i;
		int nnz = int(p_A->p[p_A->n]);
		#pragma omp parallel for default(none) shared(Len, p_i, inz) schedule(static)
		for(int i = 0; i < inz; ++ i)
			++ Len[p_i[i]]; // access violation, need to use atomics, slow
	} else
#endif // _OPENMP
	{
		const ptrdiff_t *p_i = p_A->i;
		size_t nnz = p_A->p[n];
		for(size_t i = 0; i < nnz; ++ i, ++ p_i)
			++ Len[*p_i]; // off-diagonal or discounted diagonal
	}
	for(size_t i = 0, n = p_A->n; i < n; ++ i) {
		size_t p1 = p_A->p[i], p2 = p_A->p[i + 1];
		_ASSERTE(p1 < p2); // no empty columns!
		_ASSERTE(p_A->i[p2 - 1] == i); // diagonal is here (skip)
		//-- p2; // skip the diagonal
		_ASSERTE(Len[i] > 0); // makes sure the next line will not underflow
		//-- Len[i]; // counted the diagonal as well
		_ASSERTE(p2 > p1); // makes sure the next line will not overflow
		Len[i] += p2 - p1 - 2; // accounts for diagonal, twice
	}*/
	const ptrdiff_t *Ap = p_A->p, *Ai = p_A->i; // antialiass
	for(size_t i = 0, p1 = Ap[0]; i < n; ++ i) {
		size_t p2 = Ap[i + 1];
		_ASSERTE(p1 < p2); // no empty columns!
		_ASSERTE(Ai[p2 - 1] == i); // diagonal is here (skip)
		-- p2; // skip the diagonal
		_ASSERTE(p2 >= p1); // makes sure the next line will not overflow
		Len[i] += p2 - p1;
		for(; p1 < p2; ++ p1) {
			size_t j = Ai[p1];
			_ASSERTE(j < i); // upper-triangular and we skipped the diagonal before the loop
			++ Len[j]; // off-diagonal
		}
		++ p1;
		_ASSERTE(p1 == Ap[i + 1]); // this will soo fail on jumbled matrices (_DEBUG checks for them)
	}
	/*//size_t nzaat2 = 0;
	for(size_t i = 0, n = p_A->n; i < n; ++ i) {
		size_t p1 = p_A->p[i], p2 = p_A->p[i + 1];
		_ASSERTE(p1 < p2); // no empty columns!
		_ASSERTE(p_A->i[p2 - 1] == i); // diagonal is here (skip)
		-- p2; // skip the diagonal
		_ASSERTE(p2 >= p1); // makes sure the next line will not overflow
		Len[i] += p2 - p1;
		for(; p1 < p2; ++ p1) {
			size_t j = p_A->i[p1];
			_ASSERTE(j <= i); // upper-triangular
			_ASSERTE(i != j); // we skipped the diagonal before the loop
			//if(i != j) {
			//++ Len[i]; // excluding diagonal! // moved above
			//++ nzaat2;
				++ Len[j]; // off-diagonal
				//++ nzaat2;
			//}
		}
	}
	//_ASSERTE(nzaat == nzaat2);*/
	// simpler implementation of camd_aat() for upper-triangular matrices
#else // 1
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	size_t nzaat = __camd_internal::camd_l_aat(n, p_A->p, p_A->i,
		(SuiteSparse_long*)Len, (SuiteSparse_long*)&m_ordering[0], 0);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	size_t nzaat = __camd_internal::camd_aat(n, p_A->p, p_A->i,
		(int*)Len, (int*)&m_ordering[0], 0);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
#endif // 1
	_ASSERTE((std::max(nz, n) - n <= nzaat) && (nzaat <= 2 * (size_t)nz));
	// determine the symmetry and count off-diagonal nonzeros in A+A'

	size_t slen = nzaat;
	if(slen > SIZE_MAX - nzaat / 5) // check overflow
		throw std::bad_alloc(); // would get CAMD_OUT_OF_MEMORY
	slen += nzaat / 5;
	size_t n1 = n + 1; // for sure does not overflow
	if(n1 > SIZE_MAX / 8 || (n1 == SIZE_MAX && SIZE_MAX % 8) || slen > SIZE_MAX - 8 * n1) // check overflow
		throw std::bad_alloc(); // would get CAMD_OUT_OF_MEMORY
	slen += 8 * n1; // note it is 1 bigger than it has to be // todo - check if we can save memory (probably not, CAMD requires 1 more than AMD, that will be it)
	if(m_camd_workspace1.capacity() < slen)
		m_camd_workspace1.clear(); // avoid copying data if it needs to realloc
	m_camd_workspace1.resize(slen);
	size_t *S = &m_camd_workspace1[0];
	// allocate workspace for matrix, elbow room (size n), and 7 (n + 1) vectors

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	__camd_internal::camd_l1(n, p_A->p, p_A->i, (SuiteSparse_long*)&m_ordering[0],
		(SuiteSparse_long*)Pinv, (SuiteSparse_long*)Len, slen,
		(SuiteSparse_long*)S, 0, 0, (SuiteSparse_long*)p_constraints);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	__camd_internal::camd_1(n, p_A->p, p_A->i, (int*)&m_ordering[0], (int*)Pinv,
		(int*)Len, slen, (int*)S, 0, 0, (int*)p_constraints);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	// order the matrix

	// note that camd_1 constructs A^T | A (A is a binary matrix) first,
	// it could be optimized / calculated inplace / implemented directly in block matrix // todo

	return &m_ordering[0];
	// return ordering
}
#endif // 0

const size_t *CMatrixOrdering::p_ExtendBlockOrdering_with_Identity(size_t n_new_size) // throw(std::bad_alloc)
{
	size_t n_old_size = m_ordering.size();
	_ASSERTE(n_new_size >= n_old_size);
	m_ordering.resize(n_new_size);
	for(size_t i = n_old_size; i < n_new_size; ++ i)
		m_ordering[i] = i;
	// extend ordering with identity ordering

	return &m_ordering[0];
	// return ordering
}

/*
 *								=== ~CMatrixOrdering ===
 */
