/*
								+-----------------------------------+
								|                                   |
								| ***  Unused macro definition  *** |
								|                                   |
								|   Copyright  © -tHE SWINe- 2010   |
								|                                   |
								|             Unused.h              |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __UNUSED_INCLUDED
#define __UNUSED_INCLUDED

/**
 *	@file include/slam/Unused.h
 *	@brief The UNUSED() macro definition
 *	@author -tHE SWINe-
 *	@date 2010-11-25
 */

#ifdef UNUSED
#elif defined(__GNUC__)
/**
 *	@def UNUSED
 *	@brief marks function argument / variable as deliberately unused
 *	@param x is parameter to be marked as unused
 *	@note This is especially useful for template programming or defining common interface classes
 *		where functions having unused parameters are pretty common cause of numerous g++ warnings.
 */
#define UNUSED(x) x __attribute__((unused))
#elif defined(__LCLINT__)
/**
 *	@def UNUSED
 *	@brief marks function argument / variable as deliberately unused
 *	@param x is parameter to be marked as unused
 *	@note This is especially useful for template programming or defining common interface classes
 *		where functions having unused parameters are pretty common cause of numerous g++ warnings.
 */
#define UNUSED(x) /*@unused@*/ x
#else
/**
 *	@def UNUSED
 *	@brief marks function argument / variable as deliberately unused
 *	@param x is parameter to be marked as unused
 *	@note This is especially useful for template programming or defining common interface classes
 *		where functions having unused parameters are pretty common cause of numerous g++ warnings.
 */
#define UNUSED(x) x
#endif

#endif // __UNUSED_INCLUDED
