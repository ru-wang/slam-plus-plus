/*
								+--------------------------------+
								|                                |
								| *** Portable integer types *** |
								|                                |
								|  Copyright © -tHE SWINe- 2008  |
								|                                |
								|            Integer.h           |
								|                                |
								+--------------------------------+
*/

#pragma once
#ifndef __INTEGER_TYPES_INCLUDED
#define __INTEGER_TYPES_INCLUDED

/**
 *	@file include/slam/Integer.h
 *	@author -tHE SWINe-
 *	@brief standard integer types
 *
 *	@note This was never tested on compilers other than MSVC or g++, in case you are using another
 *		compiler, make sure either inttypes.h is included, or that the integer types really have
 *		size and signedness they are supposed to.
 *
 *	@date 2008-08-02
 *
 *	added some template-based code for determining signedness and maximal value of integer types
 *
 *	@date 2008-08-08
 *
 *	added \#ifdef for windows 64
 *
 *	@date 2009-11-12
 *
 *	added macro SIZE_MAX for maximal size of size_t type (assumes it is unsigned)
 *
 *	@date 2010-10-18
 *
 *	Added some common "bit hack" functions: b_Is_POT(), n_Make_POT(), n_Make_Lower_POT(),
 *	n_Align_Up(), n_Align_Up_POT(), n_SetBit_Num() and n_Log2().
 *
 *	@date 2010-10-29
 *
 *	Unified windows detection macro to "#if defined(_WIN32) || defined(_WIN64)".
 *
 *	@date 2010-11-25
 *
 *	Added compile-time versions of the common "bit hack" functions: b_Is_POT_Static(), n_Make_POT_Static(),
 *	n_Make_Lower_POT_Static(), n_Align_Up_Static(), n_Align_Up_POT_Static(), n_SetBit_Num_Static()
 *	and n_Log2_Static()
 *
 *	@date 2010-11-28
 *
 *	Note for linux builds - the size_t type is not declared in this header, and must be declared
 *	before including this header. Including "new" or "sys/types.h" should do.
 *
 *	@date 2011-07-18
 *
 *	Changed definitions of signed integer types int8_t, int16_t and int32_t (added explicit 'signed').
 *
 *	@date 2011-11-01
 *
 *	Changed some 64-bit integer constants to two-part 32-bit constants to avoid some g++ complaints.
 *
 *	@date 2013-01-13
 *
 *	Added the PRIsize and PRIdiff macros for printing size_t and ptrdiff_t integer types.
 *
 */

#include "slam/Unused.h"

#if !defined(_WIN32) && !defined(_WIN64)
#include <sys/types.h>
#endif // !_WIN32 && !_WIN64
// for linux builds

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
	#include <inttypes.h>
#elif defined(__sun__)
	#include <inttypes.h>
	#if defined(__STDC__)
		typedef char int8_t; /**< 8 bit signed integer */
		typedef unsigned char uint8_t; /**< 8 bit unsigned integer */
		typedef short int16_t; /**< 16 bit signed integer */
		typedef unsigned short uint16_t; /**< 16 bit unsigned integer */
		typedef int int32_t; /**< 32 bit signed integer */
		typedef unsigned int uint32_t; /**< 32 bit unsigned integer */
		#if defined(__arch64__)
			typedef long int64_t; /**< 64 bit signed integer */
			typedef unsigned long uint64_t; /**< 64 bit unsigned integer */
		#else
			typedef long long int int64_t; /**< 64 bit signed integer */
			typedef unsigned long long int uint64_t; /**< 64 bit unsigned integer */
		#endif /* __arch64__ */
	#endif /* __STDC__ */
#elif defined(__VMS)
	#include <inttypes.h>
#elif defined(__SCO__) || defined(__USLC__)
	#include <stdint.h>
#elif defined(__UNIXOS2__) || defined(__SOL64__)
	typedef char int8_t; /**< 8 bit signed integer */
	typedef unsigned char uint8_t; /**< 8 bit unsigned integer */
	typedef short int16_t; /**< 16 bit signed integer */
	typedef unsigned short uint16_t; /**< 16 bit unsigned integer */
	typedef long int int32_t; /**< 32 bit signed integer */
	typedef unsigned long int uint32_t; /**< 32 bit unsigned integer */
	typedef long long int int64_t; /**< 64 bit signed integer */
	typedef unsigned long long int uint64_t; /**< 64 bit unsigned integer */
#elif defined(_WIN32) && defined(__GNUC__) // MinGW
	#include <stdint.h>
#elif defined(_WIN32) || defined (_WIN64) // MSVC
	typedef signed __int8 int8_t; /**< 8 bit signed integer */
	typedef signed __int16 int16_t; /**< 8 bit unsigned integer */
	typedef signed __int32 int32_t; /**< 16 bit signed integer */
	typedef signed __int64 int64_t; /**< 16 bit unsigned integer */
	typedef unsigned __int8 uint8_t; /**< 32 bit signed integer */
	typedef unsigned __int16 uint16_t; /**< 32 bit unsigned integer */
	typedef unsigned __int32 uint32_t; /**< 64 bit signed integer */
	typedef unsigned __int64 uint64_t; /**< 64 bit unsigned integer */
#else // fallback
	#include <inttypes.h>
#endif
// define int types

#ifndef UINT8_MAX
/**
 *	@brief maximal value of 8 bit unsigned integer (255)
 */
#define UINT8_MAX  (uint8_t(-1))
#endif

#ifndef UINT16_MAX
/**
 *	@brief maximal value of 16 bit unsigned integer (65535)
 */
#define UINT16_MAX (uint16_t(-1))
#endif

#ifndef UINT32_MAX
/**
 *	@brief maximal value of 32 bit unsigned integer (4294967295)
 */
#define UINT32_MAX (uint32_t(-1))
#endif

#ifndef UINT64_MAX
/**
 *	@brief maximal value of 64 bit unsigned integer (18446744073709551615)
 */
#define UINT64_MAX (uint64_t(-1))
#endif

#ifndef INT8_MAX
/**
 *	@brief maximal value of 8 bit signed integer (127)
 */
#define INT8_MAX  (int8_t(UINT8_MAX / 2))
#endif

#ifndef INT16_MAX
/**
 *	@brief maximal value of 16 bit signed integer (32767)
 */
#define INT16_MAX (int16_t(UINT16_MAX / 2))
#endif

#ifndef INT32_MAX
/**
 *	@brief maximal value of 32 bit signed integer (2147483647)
 */
#define INT32_MAX (int32_t(UINT32_MAX / 2))
#endif

#ifndef INT64_MAX
/**
 *	@brief maximal value of 64 bit signed integer (9223372036854775807)
 */
#define INT64_MAX (int64_t(UINT64_MAX / 2))
#endif

#ifndef INT8_MIN
/**
 *	@brief minimal value of 8 bit signed integer (-128)
 */
#define INT8_MIN  (int8_t(-INT8_MAX - 1))
#endif

#ifndef INT16_MIN
/**
 *	@brief minimal value of 16 bit signed integer (-32768)
 */
#define INT16_MIN (int16_t(-INT16_MAX - 1))
#endif

#ifndef INT32_MIN
/**
 *	@brief minimal value of 32 bit signed integer (-2147483648)
 */
#define INT32_MIN (int32_t(-INT32_MAX - 1))
#endif

#ifndef INT64_MIN
/**
 *	@brief minimal value of 64 bit signed integer (-9223372036854775808)
 */
#define INT64_MIN (int64_t(-INT64_MAX - 1))
#endif

#ifndef SIZE_MAX
/**
 *	@brief maximal value of type size_t (same as UINT32_MAX or UINT64_MAX)
 */
#define SIZE_MAX  (size_t(-1))
#endif
// calc limits

/**
 *	@def PRIu64
 *	@brief printf format strings for unsigned 64-bit integer
 */
/**
 *	@def PRI64
 *	@brief printf format strings for signed 64-bit integer
 */
/**
 *	@def _PRIu64
 *	@brief printf format strings for unsigned 64-bit integer, without the percent character
 */
/**
 *	@def _PRI64
 *	@brief printf format strings for signed 64-bit integer, without the percent character
 */
#ifndef PRIu64
#if defined(_MSC_VER) && !defined(__MWERKS__)
#define PRIu64 "%I64u"
#define PRI64 "%I64d"
#define _PRIu64 "I64u"
#define _PRI64 "I64d"
#else // MSVC
#define PRIu64 "%llu"
#define PRI64 "%ll"
#define _PRIu64 "llu"
#define _PRI64 "ll"
#endif // MSVC
#endif // !PRI64U
// printf format strings for 64-bit integers

/**
 *	@def PRIsize
 *	@brief printf format strings for size_t
 */
/**
 *	@def _PRIsize
 *	@brief printf format strings for size_t, without the percent character
 */
#ifndef PRIsize
#if defined(_MSC_VER) && !defined(__MWERKS__)
#define PRIsize "%u"
#define _PRIsize "u"
#else // MSVC
#define PRIsize "%zd"
#define _PRIsize "zd"
#endif // MSVC
#endif // !PRIsize

/**
 *	@def PRIdiff
 *	@brief printf format strings for ptrdiff_t
 */
/**
 *	@def _PRIdiff
 *	@brief printf format strings for ptrdiff_t, without the percent character
 */
#ifndef PRIdiff
#if defined(_MSC_VER) && !defined(__MWERKS__)
#define PRIdiff "%ld"
#define _PRIdiff "ld"
#else // MSVC
#define PRIdiff "%td"
#define _PRIdiff "td"
#endif // MSVC
#endif // !PRIsize

/*
 *								=== templates for maximal value of integer types ===
 */

/**
 *	@brief data type signedness predicate
 *
 *	Decides wheter a type is signed, or not. Value of CIsSigned<_Ty>::result
 *		is true if type _Ty is signed, or false if it's unsigned.
 *
 *	@tparam _Ty is type tested for being signed
 */
template <class _Ty>
struct CIsSigned {
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		result = _Ty(-1) < 1 /**< contains true if type _Ty is signed, or false if it's unsigned */
	};
	// should be '_Ty(-1) < 0', using 1 to avoid g++ warning 'unsigned comparison < 0 is always false'
};

/**
 *	@brief integer type predicate
 *
 *	Decides wheter a type is integer, or floating point. Value of CIsInteger<_Ty>::result
 *		is set to true uf _Ty is integer type, otherwise it's false.
 *
 *	@tparam _Ty is type tested for being integer
 */
template <class _Ty>
struct CIsInteger {
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		result = _Ty(.5) == 0 /**< is set to true uf _Ty is integer type, otherwise it's false */
	};
};

/**
 *	@brief maximal integer value template
 *
 *	Determines maximal value of integer type _Ty, result is returned by CMaxIntValue<_Ty>::result().
 *
 *	@tparam _Ty is type to determine maximal value of
 *
 *	@note _Ty must be integer type, it won't work with float or double.
 */
template <class _Ty>
struct CMaxIntValue {
private:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		n_bit_num = 8 * sizeof(_Ty),
		b_unsigned = (CIsSigned<_Ty>::result)? 0 : 1
	};

/**
 *	@brief this calculates maximal value of _Ty as 64-bit unsigned integer
 *	@note this may not be used outside CMaxIntValue
 */
#define _Ty_max (((uint64_t(1) << (n_bit_num - 1)) - 1) + b_unsigned * (uint64_t(1) << (n_bit_num - 1)))

#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400 // MSVC 6.0
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		n_result_hi = uint32_t(_Ty_max >> 32),
		n_result_lo = uint32_t(_Ty_max & 0xffffffff)
	};
	// msvc doesn't allow static const initialization so we keep it as hi-lo
	// enum values (note this relies on enum values being 32-bit integers)

public:
	/**
	 *	@brief gets maximal value of type _Ty
	 *
	 *	@return Returns maximal value of type _Ty.
	 *
	 *	@note Result is uint64_t, in case larger types are implemented by your compiler,
	 *		this will not work for them.
	 *	@note This doesn't carry-out any calculations, it is just not possible to store
	 *		64-bit number using enum.
	 */
	static inline uint64_t result()
	{
		_ASSERTE(CIsInteger<_Ty>::result); // only works for integers
		return (uint64_t(n_result_hi) << 32) | uint32_t(n_result_lo);
	}
#else // MSVC 6.0
	static const uint64_t n_result = _Ty_max; // g++ allows this (preferable)

public:
	/**
	 *	@brief gets maximal value of type _Ty
	 *
	 *	@return Returns maximal value of type _Ty.
	 *
	 *	@note Result is uint64_t, in case larger types are implemented by your compiler,
	 *		this will not work for them.
	 *	@note This doesn't carry-out any calculations, it is just not possible to store
	 *		64-bit number using enum.
	 */
	static inline uint64_t result()
	{
		_ASSERTE(CIsInteger<_Ty>::result); // only works for integers
		return n_result;
	}
#endif // MSVC

#undef _Ty_max
};

/**
 *	@brief gets maximal value of type _Ty
 *	@tparam _Ty is integer data type
 *	@param[in] x is parameter of type _Ty (the value is ignored here)
 *	@return Returns maximal value of type _Ty.
 *	@note Result is uint64_t, in case larger types are implemented by your compiler,
 *		this will not work for them.
 */
template <class _Ty>
inline uint64_t n_MaxIntValue(_Ty UNUSED(x))
{
	return CMaxIntValue<_Ty>::result();
}

/*
 *								=== ~templates for maximal value of integer types ===
 */

/*
 *								=== bit hacks ===
 */

#ifndef _ASSERTE
#include <assert.h>
#define _ASSERTE(x) assert(x)
#endif // _ASSERTE

/**
 *	@brief hides inside implementation of helper classes for the fast bit hacks
 */
namespace ul_bithacks {

#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
/**
 *	@brief simple template for calculating power of two
 *
 *	MSVC60 doesn't handle partial template specializations, so we need to hide _Ty parameter to CMakePOT
 *	in order to be able to declare full specialization of CShifter (otherwise gives infinite template recursion error).
 *
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 */
template <class _Ty>
class CMakePOT {
public:
	/**
	 *	@brief shift series expansion template
	 *	@tparam n_shift is amount to shift the number
	 */
	template <int n_shift>
	class CShifter {
	public:
		/**
		 *	@brief calculates shift series
		 *	@param[in] n_x is input number
		 *	@return Returns n_x or-ed together with all power-of-two fractions of n_shift down to zero.
		 */
		static inline _Ty n_Or_Shift(_Ty n_x)
		{
			n_x = CShifter<n_shift / 2>::n_Or_Shift(n_x);
			return n_x | (n_x >> n_shift);
		}
	};

	/**
	 *	@brief shift series template specialization (recursion terminator)
	 */
	template <>
	class CShifter<0> {
	public:
		/**
		 *	@brief identity function
		 *	@param[in] n_x is input number
		 *	@return Returns n_x.
		 */
		static inline _Ty n_Or_Shift(_Ty n_x)
		{
			return n_x;
		}
	};

	/**
	 *	@brief calculates power of two greater or equal to the argument
	 *
	 *	@param[in] n_x is input number
	 *
	 *	@return Returns power of two greater or equal to n_x. In case such number is not representable by given type, returns null.
	 *
	 *	@note In case _Ty is signed and n_x is greater than the largest power of two, representable by this type (can be set to zero by masking-out the sign bit).
	 */
	static inline _Ty n_Make_POT(_Ty n_x)
	{
		return CShifter<sizeof(_Ty) * 8 / 2>::n_Or_Shift(n_x - 1) + 1;
	}

	/**
	 *	@brief fill ones right of the highest one that is set
	 *
	 *	@param[in] n_x is input number
	 *
	 *	@return Returns n_x which is right-filled with ones.
	 */
	static inline _Ty n_RightFill_Ones(_Ty n_x)
	{
		return CShifter<sizeof(_Ty) * 8 / 2>::n_Or_Shift(n_x);
	}

	/**
	 *	@brief calculates power of two below or equal to the argument
	 *
	 *	@param[in] n_x is input number
	 *
	 *	@return Returns power of two below or equal to n_x.
	 */
	static inline _Ty n_Make_Lower_POT(_Ty n_x)
	{
		return CShifter<sizeof(_Ty) * 8 / 2>::n_Or_Shift(n_x >> 1) + 1; // this is very similar to n_Make_POT(), we just divide the number by two and add one in order to get lower power of two
	}
};
#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
// this is the right way to go

/**
 *	@brief shift series expansion template
 *
 *	@tparam _Ty is integer data type
 *	@tparam n_shift is amount to shift the number
 */
template <class _Ty, int n_shift>
class CShifter {
public:
	/**
	 *	@brief calculates shift series
	 *	@param[in] n_x is input number
	 *	@return Returns n_x or-ed together with all power-of-two fractions of n_shift down to zero.
	 */
	static inline _Ty n_Or_Shift(_Ty n_x)
	{
		n_x = CShifter<_Ty, n_shift / 2>::n_Or_Shift(n_x);
		return n_x | (n_x >> n_shift);
	}
};

/**
 *	@brief shift series template specialization (recursion terminator)
 *	@tparam _Ty is integer data type
 */
template <class _Ty>
class CShifter<_Ty, 0> {
public:
	/**
	 *	@brief identity function
	 *	@param[in] n_x is input number
	 *	@return Returns n_x.
	 */
	static inline _Ty n_Or_Shift(_Ty n_x)
	{
		return n_x;
	}
};
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400

#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
/**
 *	@brief static shift series expansion template
 *
 *	MSVC60 doesn't handle partial template specializations, so we need to hide the n_x parameter to CMakePOT
 *	in order to be able to declare full specialization of CShifter_Static.
 *
 *	@tparam n_x is the input value to shift
 */
template <const size_t n_x>
class CMakePOT_Static {
public:
	/**
	 *	@brief static shift series expansion template
	 *	@tparam n_shift is amount to shift n_x
	 */
	template <const int n_shift>
	class CShifter_Static {
	public:
		/**
		 *	@brief result, stored as enum
		 */
		enum {
			/**
			 *	@brief equals n_x or-ed together with all power-of-two fractions of n_shift down to zero
			 */
			result = CShifter_Static<n_shift / 2>::result | (CShifter_Static<n_shift / 2>::result >> n_shift)
		};
	};

	/**
	 *	@brief static shift series template specialization (recursion terminator)
	 */
	template <>
	class CShifter_Static<0> {
	public:
		/**
		 *	@brief result, stored as enum
		 */
		enum {
			/**
			 *	@brief equals n_x
			 */
			result = n_x
		};
	};
};
#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
// this is the right way to go

/**
 *	@brief static shift series expansion template
 *	@tparam n_x is the input value to shift
 *	@tparam n_shift is amount to shift n_x
 */
template <const size_t n_x, const int n_shift>
class CShifter_Static {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		/**
		 *	@brief equals n_x or-ed together with all power-of-two fractions of n_shift down to zero
		 */
		result = CShifter_Static<n_x, n_shift / 2>::result | (CShifter_Static<n_x, n_shift / 2>::result >> n_shift)
	};
};

/**
 *	@brief static shift series template specialization (recursion terminator)
 *	@tparam n_x is the input value to shift
 */
template <const size_t n_x>
class CShifter_Static<n_x, 0> {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		/**
		 *	@brief equals n_x
		 */
		result = n_x
	};
};
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400

/**
 *	@brief simple template for counting number of bits set in a number
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 *	@note Due to the limit of enum values to 32 bits, this only works for integers
 *		up to 32 bits. Specialized versions are required to handle 64 bit integers.
 */
template <class _Ty>
class CBitCounter {
public:
	/**
	 *	@brief template for calculating repeating byte patterns of variable lengths
	 *	@tparam n_byte is value of the pattern to be repeated
	 */
	template <const int n_byte>
	class CByteRepeater {
	public:
		/**
		 *	@brief result, stored as enum
		 */
		enum {
			n_value = n_byte |
				((sizeof(_Ty) >= 4)? n_byte << 24 : 0) |
				((sizeof(_Ty) >= 3)? n_byte << 16 : 0) |
				((sizeof(_Ty) >= 2)? n_byte << 8 : 0) /**< the result; contains n_byte, shifted and or-ed to fill entire _Ty */
		};
	};

	/**
	 *	@brief calculates number of bits set, in parallel
	 *
	 *	@param[in] n_x is input value, which must not be negative
	 *
	 *	@return Returns number of bits set in n_x.
	 */
	static inline unsigned int n_SetBit_Num(_Ty n_x)
	{
		_ASSERTE(sizeof(_Ty) <= 4);
		n_x -= (n_x >> 1) & _Ty(CByteRepeater<0x55>::n_value); // sums even + odd (wow), all bits are still significant
		n_x = (n_x & _Ty(CByteRepeater<0x33>::n_value)) + ((n_x >> 2) & _Ty(CByteRepeater<0x33>::n_value)); // sums pairs, only three bits in each nibble are now significant - can add without masking
		n_x = (n_x + (n_x >> 4)) & _Ty(CByteRepeater<0x0f>::n_value); // sums nibbles, using common mask; only every other nibble is now significant
		return (unsigned int)(_Ty(n_x * _Ty(CByteRepeater<0x01>::n_value)) >> (sizeof(_Ty) * 8 - 8)); // can sum everything using multiplication, result is guaranteed not to overflow a single byte
	}
};

/**
 *	@brief simple template for counting number of bits set in an int64_t number
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 */
template <>
class CBitCounter<int64_t> {
public:
	/**
	 *	@brief calculates number of bits set, in parallel
	 *	@param[in] n_x is input value, which must not be negative
	 *	@return Returns number of bits set in n_x.
	 */
	static inline unsigned int n_SetBit_Num(int64_t n_x)
	{
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
		n_x -= (n_x >> 1) & 0x5555555555555555; // sums even + odd (wow), all bits are still significant
		n_x = (n_x & 0x3333333333333333) + ((n_x >> 2) & 0x3333333333333333); // sums pairs, only three bits in each nibble are now significant - can add without masking
		n_x = (n_x + (n_x >> 4)) & 0x0f0f0f0f0f0f0f0f; // sums nibbles, using common mask; only every other nibble is now significant
		return (unsigned int)((n_x * 0x0101010101010101) >> 56); // can sum everything using multiplication, result is guaranteed not to overflow a single byte
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		const int64_t n_fives = 0x55555555 | (int64_t(0x55555555) << 32);
		const int64_t n_threes = 0x33333333 | (int64_t(0x33333333) << 32);
		const int64_t n_zerofs = 0x0f0f0f0f | (int64_t(0x0f0f0f0f) << 32);
		const int64_t n_zeroones = 0x01010101 | (int64_t(0x01010101) << 32);
		// it would be better to use LL suffixes, but g++ is yappin about it

		n_x -= (n_x >> 1) & n_fives; // sums even + odd (wow), all bits are still significant
		n_x = (n_x & n_threes) + ((n_x >> 2) & n_threes); // sums pairs, only three bits in each nibble are now significant - can add without masking
		n_x = (n_x + (n_x >> 4)) & n_zerofs; // sums nibbles, using common mask; only every other nibble is now significant
		return (unsigned int)((n_x * n_zeroones) >> 56); // can sum everything using multiplication, result is guaranteed not to overflow a single byte
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	}
};

/**
 *	@brief simple template for counting number of bits set in an uint64_t number
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 */
template <>
class CBitCounter<uint64_t> {
public:
	/**
	 *	@brief calculates number of bits set, in parallel
	 *	@param[in] n_x is input value
	 *	@return Returns number of bits set in n_x.
	 */
	static inline unsigned int n_SetBit_Num(uint64_t n_x)
	{
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
		n_x -= (n_x >> 1) & 0x5555555555555555U; // sums even + odd (wow), all bits are still significant
		n_x = (n_x & 0x3333333333333333U) + ((n_x >> 2) & 0x3333333333333333U); // sums pairs, only three bits in each nibble are now significant - can add without masking
		n_x = (n_x + (n_x >> 4)) & 0x0f0f0f0f0f0f0f0fU; // sums nibbles, using common mask; only every other nibble is now significant
		return (unsigned int)((n_x * 0x0101010101010101U) >> 56); // can sum everything using multiplication, result is guaranteed not to overflow a single byte
#else // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		const uint64_t n_fives = 0x55555555U | (int64_t(0x55555555U) << 32);
		const uint64_t n_threes = 0x33333333U | (int64_t(0x33333333U) << 32);
		const uint64_t n_zerofs = 0x0f0f0f0fU | (int64_t(0x0f0f0f0fU) << 32);
		const uint64_t n_zeroones = 0x01010101U | (int64_t(0x01010101U) << 32);
		// it would be better to use ULL suffixes, but g++ is yappin about it

		n_x -= (n_x >> 1) & n_fives; // sums even + odd (wow), all bits are still significant
		n_x = (n_x & n_threes) + ((n_x >> 2) & n_threes); // sums pairs, only three bits in each nibble are now significant - can add without masking
		n_x = (n_x + (n_x >> 4)) & n_zerofs; // sums nibbles, using common mask; only every other nibble is now significant
		return (unsigned int)((n_x * n_zeroones) >> 56); // can sum everything using multiplication, result is guaranteed not to overflow a single byte
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
	}
};

#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400

/**
 *	@brief simple compile-time bit counter template
 *	@tparam x is an input value where the bits should be counted
 */
template <const size_t x>
class CBitCounter_Static {
public:
	/**
	 *	@brief simple compile-time bit counter template
	 *	@tparam n is the number of remaining bits to count in x
	 */
	template <const size_t n>
	class CBitCounter2_Static {
	public:
		/**
		 *	@brief result, stored as enum
		 */
		enum {
			result = (x & 1) + CBitCounter_Static<x >> 1>::CBitCounter2_Static<n - 1>::result /**< @brief contains the sum of bits set in x */
		};
	};

	/**
	 *	@brief simple compile-time bit counter template
	 */
	template <>
	class CBitCounter2_Static<0> {
	public:
		/**
		 *	@brief result, stored as enum
		 */
		enum {
			result = 0; /**< @brief contains the sum of bits set in x */
		};
	};
};

#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
// this is the right way to go

/**
 *	@brief simple compile-time bit counter template
 *
 *	@tparam x is an input value where the bits should be counted
 *	@tparam n is the number of remaining bits to count in x
 */
template <const size_t x, const size_t n>
class CBitCounter_Static {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		result = (x & 1) + CBitCounter_Static<x / 2, n - 1>::result /**< @brief contains the sum of bits set in x */
	};
};

/**
 *	@brief simple compile-time bit counter template
 *	@tparam x is an input value where the bits should be counted
 */
template <const size_t x>
class CBitCounter_Static<x, 0> {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		result = 0 /**< @brief contains the sum of bits set in x */
	};
};

#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400

/**
 *	@brief simple compile-time base 2 logarithm template
 *	@tparam n is the number of remaining bits to count in x
 */
template <const size_t n>
class CLog2_Static {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		result = CLog2_Static<n / 2>::result + 1 /**< @brief result = log<sub>2</sub>(n) */
	};
};

/**
 *	@brief simple compile-time base 2 logarithm template
 */
template <>
class CLog2_Static<1> {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		result = 0 /**< @brief result = log<sub>2</sub>(n) */
	};
};

/**
 *	@brief simple compile-time base 2 logarithm template
 */
template <>
class CLog2_Static<0> {
public:
	/**
	 *	@brief result, stored as enum
	 */
	enum {
		result = 0 /**< @brief result = log<sub>2</sub>(n) */
	};
};

} // ul_bithacks

/**
 *	@brief determines whether a number is power of two, or not
 *	@tparam _Ty is integer data type
 *	@param[in] n_x is number being tested. note it must be positive
 *	@return Returns true if n_x is power of two, otherwise returns false.
 */
template <class _Ty>
inline bool b_Is_POT(_Ty n_x)
{
	return !(n_x & (n_x - 1));
}

/**
 *	@brief determines whether a number is power of two, or not, can be evaluated at compile-time
 *	@param[in] n_x is number being tested. note it must be positive
 *	@return Returns true if n_x is power of two, otherwise returns false.
 */
#define b_Is_POT_Static(n_x) (!((n_x) & ((n_x) - 1)))

/**
 *	@brief calculates power of two greater or equal to the argument
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 *	@param[in] n_x is input number, which must not be negative
 *	@return Returns power of two greater or equal to n_x. In case such number is not representable by given type, returns null.
 *	@note In case _Ty is signed and n_x is greater than the largest power of two, representable by this type (can be set to zero by masking-out the sign bit).
 */
template <class _Ty>
inline _Ty n_Make_POT(_Ty n_x)
{
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
	return ul_bithacks::CMakePOT<_Ty>::n_Make_POT(n_x);
#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
	return ul_bithacks::CShifter<_Ty, sizeof(_Ty) * 8 / 2>::n_Or_Shift(n_x - 1) + 1;
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
}

/**
 *	@brief calculates power of two greater or equal to the argument, at compile-time
 *	@param[in] n_x is input number, which must not be negative
 *	@return Returns power of two below or equal to n_x.
 *	@note As this is evaluated at compile-time and the speed is of little interest, this just calls
 *		n_Align_Up_Static() macro and the power-of-two requirement is therefore relaxed.
 */
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
#define n_Make_POT_Static(n_x) ul_bithacks::CMakePOT_Static<(n_x) - 1>::CShifter_Static<sizeof(n_x) * 8 / 2>::result + 1
#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
#define n_Make_POT_Static(n_x) ul_bithacks::CShifter_Static<(n_x) - 1, sizeof(n_x) * 8 / 2>::result + 1
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400

/**
 *	@brief calculates power of two below or equal to the argument
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 *	@param[in] n_x is input number, which must not be negative
 *	@return Returns power of two below or equal to n_x.
 */
template <class _Ty>
inline _Ty n_Make_Lower_POT(_Ty n_x)
{
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
	return ul_bithacks::CMakePOT<_Ty>::n_Make_Lower_POT(n_x);
#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
	return ul_bithacks::CShifter<_Ty, sizeof(_Ty) * 8 / 2>::n_Or_Shift(n_x >> 1) + 1; // this is very similar to n_Make_POT(), we just divide the number by two and add one in order to get lower power of two
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
}

/**
 *	@brief calculates power of two below or equal to the argument, at compile-time
 *	@param[in] n_x is input number, which must not be negative
 *	@return Returns power of two below or equal to n_x.
 *	@note As this is evaluated at compile-time and the speed is of little interest, this just calls
 *		n_Align_Up_Static() macro and the power-of-two requirement is therefore relaxed.
 */
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
#define n_Make_Lower_POT_Static(n_x) ul_bithacks::CMakePOT_Static<(n_x) >> 1>::CShifter_Static<sizeof(n_x) * 8 / 2>::result + 1
#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
#define n_Make_Lower_POT_Static(n_x) ul_bithacks::CShifter_Static<(n_x) / 2, sizeof(n_x) * 8 / 2>::result + 1
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400

/**
 *	@brief aligns number up to the next multiple of given alignment
 *
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 *
 *	@param[in] n_x is value to be aligned
 *	@param[in] n_alignment is alignment
 *
 *	@return Returns such y, so that y >= n_x and y < n_x + n_alignment and y % n_alignment = 0.
 *
 *	@note This function uses the modulo operator; for power of two alignments,
 *		it might be more effective to use n_Align_Up_POT().
 */
template <class _Ty>
inline _Ty n_Align_Up(_Ty n_x, _Ty n_alignment)
{
	n_x += n_alignment - 1;
	return n_x - n_x % n_alignment;
}

/**
 *	@brief aligns number up to the next multiple of given alignment, can be evaluated at compile-time
 *
 *	@param[in] n_x is value to be aligned
 *	@param[in] n_alignment is alignment
 *
 *	@return Returns such y, so that y >= n_x and y < n_x + n_alignment and y % n_alignment = 0.
 */
#define n_Align_Up_Static(n_x,n_alignment) (((n_x) + (n_alignment) - 1) - ((n_x) + (n_alignment) - 1) % (n_alignment))

/**
 *	@brief aligns number up to the next multiple of given power-of-two alignment
 *
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 *
 *	@param[in] n_x is value to be aligned
 *	@param[in] n_pot_alignment is alignment, it must be power of two
 *
 *	@return Returns such y, so that y >= n_x and y < n_x + n_alignment and y % n_alignment = 0.
 */
template <class _Ty>
inline _Ty n_Align_Up_POT(_Ty n_x, _Ty n_pot_alignment)
{
	_ASSERTE(b_Is_POT(n_pot_alignment)); // alignment must be power of two
	-- n_pot_alignment;
	return (n_x + n_pot_alignment) & ~n_pot_alignment;
}

/**
 *	@brief aligns number up to the next multiple of given (power-of-two) alignment, can be evaluated at compile-time
 *
 *	@param[in] n_x is value to be aligned
 *	@param[in] n_alignment is alignment
 *
 *	@return Returns such y, so that y >= n_x and y < n_x + n_alignment and y % n_alignment = 0.
 *
 *	@note As this is evaluated at compile-time and the speed is of little interest, this just calls
 *		n_Align_Up_Static() macro and the power-of-two requirement is therefore relaxed.
 */
#define n_Align_Up_POT_Static(n_x,n_alignment) (n_Align_Up_Static((n_x), (n_alignment)))

/**
 *	@brief calculates number of bits set, in parallel
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 *	@param[in] n_x is input value, which must not be negative
 *	@return Returns number of bits set in n_x.
 */
template <class _Ty>
inline unsigned int n_SetBit_Num(_Ty n_x)
{
	return ul_bithacks::CBitCounter<_Ty>::n_SetBit_Num(n_x);
}

/**
 *	@brief calculates number of bits set in a given number at compile-time
 *	@param[in] n_x is input value, which must not be negative
 *	@return Returns number of bits set in n_x.
 *	@note This uses templates and enums, and may therefore have limited range of input
 *		values (some compilers limit maximum size of enum to 31 or 32 bits).
 */
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
#define n_SetBit_Num_Static(n_x) (ul_bithacks::CBitCounter_Static<(n_x)>::CBitCounter2_Static<sizeof(n_x) * 8>::result)
#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
#define n_SetBit_Num_Static(n_x) (ul_bithacks::CBitCounter_Static<(n_x), sizeof(n_x) * 8>::result)
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400

/**
 *	@brief calculates base two logarithm of a given number (round down)
 *	@tparam _Ty is integer data type; this works for both signed and unsigned types
 *	@param[in] n_x is input value, which must not be negative
 *	@return Returns floor(log2(n_x)).
 */
template <class _Ty>
inline _Ty n_Log2(_Ty n_x)
{
	if(!n_x)
		return 0;
	n_x = n_Make_Lower_POT(n_x) - 1; // set all bits left of highest set bit in n_x
	return n_SetBit_Num(n_x); // count bits set
}

/**
 *	@brief calculates base two logarithm of a given number (round down) at compile-time
 *	@param[in] n_x is input value, which must not be negative
 *	@return Returns floor(log2(n_x)).
 *	@note This uses templates and enums, and may therefore have limited range of input
 *		values (some compilers limit maximum size of enum to 31 or 32 bits).
 */
#define n_Log2_Static(n_x) (ul_bithacks::CLog2_Static<(n_x)>::result)

/**
 *	@brief fill ones right of the highest one that is set
 *	@param[in] n_x is input number
 *	@return Returns n_x which is right-filled with ones.
 */
template <class _Ty>
inline _Ty n_RightFill_Ones(_Ty n_x)
{
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER < 1400
	return ul_bithacks::CMakePOT<_Ty>::n_RightFill_Ones(n_x);
#else // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
	return ul_bithacks::CShifter<_Ty, sizeof(_Ty) * 8 / 2>::n_Or_Shift(n_x);
#endif // _MSC_VER && !__MWERKS__ && _MSC_VER < 1400
}

/**
 *	@brief fill ones right of the highest one that is set at compile-time
 *	@param[in] n_x is input number
 *	@return Returns n_x which is right-filled with ones.
 *	@note This uses templates and enums, and may therefore have limited range of input
 *		values (some compilers limit maximum size of enum to 31 or 32 bits).
 */
#define n_RightFill_Ones_Static(n_x) (n_Make_POT_Static((n_x) + 1) - 1)

/**
 *	@brief calculates mask with a specified right-most bits set
 *	@tparam _Ty is integer data type
 *	@param[in] n_bit_num is number of bits to set (must be less than or equal to 8 * sizeof(_Ty))
 *	@return Returns mask with n_bit_num right-most bits set
 */
template <class _Ty>
inline _Ty n_Mask(_Ty n_bit_num)
{
	return (n_bit_num > 0)? (_Ty(1) << (n_bit_num - 1)) |
		((n_bit_num > 1)? ((_Ty(1) << (n_bit_num - 1)) - 1) : 0) : 0;
}

/**
 *	@brief calculates mask with a specified right-most bits set at compile-time
 *	@param[in] n_bit_num is number of bits to set (must be less than or equal to 8 * sizeof(_Ty))
 *	@return Returns mask with n_bit_num right-most bits set
 *	@note This uses templates and enums, and may therefore have limited range of input
 *		values (some compilers limit maximum size of enum to 31 or 32 bits).
 */
#define n_Mask_Static(n_bit_num) \
	((n_bit_num > 0)? ((1 << ((n_bit_num) - 1)) | \
	((n_bit_num > 1)? ((1 << ((n_bit_num) - 1)) - 1) : 0)) : 0)

/**
 *	@brief calculates number of leading zeros
 *	@tparam _Ty is integer data type
 *	@param[in] n_x is the input number
 *	@return Returns number of leading zeros in n_x.
 */
template <class _Ty>
inline _Ty n_LeadingZero_Num(_Ty n_x)
{
	return n_SetBit_Num(~n_RightFill_Ones(n_x));
}

/**
 *	@brief calculates number of leading zeros at compile-time
 *	@param[in] n_x is the input number
 *	@return Returns number of leading zeros in n_x.
 *	@note This uses templates and enums, and may therefore have limited range of input
 *		values (some compilers limit maximum size of enum to 31 or 32 bits).
 */
#define n_LeadingZero_Num_Static(n_x) (n_SetBit_Num_Static(~n_RightFill_Ones_Static(n_x)))

/*
 *								=== ~bit hacks ===
 */

#endif // __INTEGER_TYPES_INCLUDED
