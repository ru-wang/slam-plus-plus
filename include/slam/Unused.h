/*
								+-----------------------------------+
								|                                   |
								| ***  Unused macro definition  *** |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2010  |
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

#ifndef UNUSED
/**
 *	@def UNUSED
 *	@brief marks function argument	 / variable as deliberately unused
 *	@param x is parameter to be marked as unused
 *	@note This is especially useful for template programming or defining common interface classes
 *		where functions having unused parameters are pretty common cause of numerous g++ warnings.
 */
#if defined(__GNUC__)
#define UNUSED(x) x __attribute__((unused))
#elif defined(__LCLINT__)
#define UNUSED(x) /*@unused@*/ x
#else // __GNUC__
#define UNUSED(x) x
#endif // __GNUC__
#endif // UNUSED

#endif // __UNUSED_INCLUDED
