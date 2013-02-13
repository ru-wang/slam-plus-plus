/*
								+----------------------------------+
								|                                  |
								| ***   Multi-platform timer   *** |
								|                                  |
								|   Copyright © -tHE SWINe- 2006   |
								|                                  |
								|             Timer.h              |
								|                                  |
								+----------------------------------+
*/

#pragma once
#ifndef __TIMER_INCLUDED
#define __TIMER_INCLUDED

/**
 *	@file include/slam/Timer.h
 *	@brief multi-platform timer
 *	@author -tHE SWINe-
 *
 *	@date 2006-08-23
 *
 *	passed code revision
 *
 *	removed some unnecessary \#defines for linux
 *	added higher precision timer for linux (instead of clock(), gettimeofday() can be used)
 *
 *	@date 2007-05-10
 *
 *	added conditional \#include <windows.h> to make it easier to include
 *
 *	@date 2008-04-09
 *
 *	added GetTickCount method, redesigned timer code. should cope nicely with counter overflow
 *	(previous version of timer with QueryPerformanceCounter had some trouble on core duo cpu's)
 *
 *	@date 2008-08-02
 *
 *	cleared up code arround n_MaxIntValue() a bit, removed MAX_VALUE_FUNCTION and MAX_VALUE_CONST
 *	and moved all the stuff to Integer.h
 *
 *	@date 2008-08-08
 *
 *	added \#ifdef for windows 64
 *
 *	@date 2009-05-04
 *
 *	fixed mixed windows / linux line endings
 *
 *	@date 2010-01-08
 *
 *	added the PRItime and PRItimeparams macros for easy printfing of time values
 *
 *	@date 2010-10-29
 *
 *	Unified windows detection macro to "#if defined(_WIN32) || defined(_WIN64)".
 *
 *	Added CTimer::Reset() function in place of CTimer::ResetTimer() (the old function
 *	is kept for backward compatibility).
 *
 *	@date 2013-01-13
 *
 *	Added new and more precise TIMER_USE_READREALTIME for linux.
 *
 */

/**
 *	@def PRItime
 *
 *	@brief printf macro for printing time in the "hh:mm:ss.ss" format
 *
 *	This macro is used similar to standard PRI64 macro, it just defines
 *	proper formatting string. It should be used in conjunction with the
 *	PRItimeparams macro, as follows:
 *	@code
 *	double f_time = 123456;
 *	printf("elapsed time: " PRItime "\n", PRItimeparams(f_time));
 *	@endcode
 */
#define PRItime "%s%02d:%02d:%05.2f"

/**
 *	@def PRItimeprecise
 *
 *	@brief printf macro for printing time in the "hh:mm:ss.sssss" format
 *
 *	This macro is used similar to standard PRI64 macro, it just defines
 *	proper formatting string. It should be used in conjunction with the
 *	PRItimeparams macro, as follows:
 *	@code
 *	double f_time = 123456;
 *	printf("elapsed time: " PRItimeprecise "\n", PRItimeparams(f_time));
 *	@endcode
 */
#define PRItimeprecise "%s%02d:%02d:%08.5f"

/**
 *	@def PRItimeparams
 *
 *	@brief splits it's argument to integer hours, integer minutes
 *		and floating-point seconds
 *
 *	Calculates sign and three values in hours, minutes and seconds from given time.
 *	These values are separated using colons, and are intended to be used as
 *	function parameters. It could be used as:
 *
 *	@code
 *	void PrintTime3(const char *sign, int hh, int mm, float ss)
 *	{
 *		printf("time is: %s%02d:%02d:%05.2f ...\n", sign, hh, mm, ss);
 *	}
 *
 *	void PrintTime(float seconds)
 *	{
 *		PrintTime3(PRItimeparams(seconds));
 *	}
 *	@endcode
 *
 *	It's however intended to be used along with PRItime. Refer to PRItime
 *	macro documentation for more details.
 *
 *	@param[in] f is time in seconds to be split into hours, minutes, seconds
 */
#define PRItimeparams(f) (((f) >= 0)? "" : "-"), int(f) / 3600, \
	(((f) >= 0)? 1 : -1) * (int(f) / 60 - 60 * (int(f) / 3600)), \
	(((f) >= 0)? 1 : -1) * ((f) - 60 * (int(f) / 60))

/**
 *	@def TIMER_ALLOW_QPC
 *	@brief this macro enables use of QueryPerformanceCounter() function (win32-only)
 */
#define TIMER_ALLOW_QPC

/**
 *	@def TIMER_ALLOW_GETTICKCOUNT
 *	@brief this macro enables use of GetTickCount() function (win32-only)
 */
#define TIMER_ALLOW_GETTICKCOUNT

/**
 *	@def TIMER_ALLOW_READREALTIME
 *	@brief this macro enables use of read_real_time() function (unix-only)
 *	@note This is only available in AIX.
 */
//#define TIMER_ALLOW_READREALTIME

/**
 *	@def TIMER_ALLOW_GETTIMEOFDAY
 *	@brief this macro enables use of gettimeofday() function (unix-only)
 */
#define TIMER_ALLOW_GETTIMEOFDAY

/**
 *	@def TIMER_METHOD_OVERRIDE
 *	@brief this macro disables automatic timer method selection
 *
 *	One of following macros must be defined as well to choose timer method:
 *		TIMER_USE_QPC			(win32-only, resolution less than 1 usec, depends on cpu)
 *		TIMER_USE_GETTICKCOUNT	(win32-only, resolution 1 msec)
 *		TIMER_USE_READREALTIME	(unix-only, resolution less than 1 usec, depends on cpu)
 *		TIMER_USE_GETTIMEOFDAY	(unix-only, resolution 1 usec)
 *		TIMER_USE_CLOCK			(default fallback, resolution 1 msec or less, depends on os)
 */
//#define TIMER_METHOD_OVERRIDE

#ifndef TIMER_METHOD_OVERRIDE
#if defined(_WIN32) || defined (_WIN64)
#if defined(TIMER_ALLOW_QPC)
/**
 *	@def TIMER_USE_QPC
 *	@brief selected timer method
 */
#define TIMER_USE_QPC
#elif defined(TIMER_ALLOW_GETTICKCOUNT)
/**
 *	@def TIMER_USE_GETTICKCOUNT
 *	@brief selected timer method
 */
#define TIMER_USE_GETTICKCOUNT
#else
/**
 *	@def TIMER_USE_CLOCK
 *	@brief selected timer method
 */
#define TIMER_USE_CLOCK
#endif
#else // _WIN32, _WIN64
#if defined(TIMER_ALLOW_READREALTIME)
/**
 *	@def TIMER_USE_READREALTIME
 *	@brief selected timer method
 */
#define TIMER_USE_READREALTIME
#elif defined(TIMER_ALLOW_GETTIMEOFDAY)
/**
 *	@def TIMER_USE_GETTIMEOFDAY
 *	@brief selected timer method
 */
#define TIMER_USE_GETTIMEOFDAY
#else
/**
 *	@def TIMER_USE_CLOCK
 *	@brief selected timer method
 */
#define TIMER_USE_CLOCK
#endif
#endif // _WIN32, _WIN64
#endif // TIMER_METHOD_OVERRIDE

#include "slam/Integer.h"

#ifdef TIMER_USE_READREALTIME
#include <sys/time.h>
//#include <sys/systemcfg.h> // timebasestruct_t
#endif // TIMER_USE_READREALTIME

/**
 *	@brief a simple timer class
 */
class CTimer {
	int64_t m_n_freq;
#ifdef TIMER_USE_READREALTIME
	mutable timebasestruct_t m_t_time; // feeling a bit paranoid; todo - use uint64_t as well, test on multicore platforms
#else // TIMER_USE_READREALTIME
	mutable int64_t m_n_time;
#endif // TIMER_USE_READREALTIME
	mutable double m_f_time;

public:
	/**
	 *	@brief default constructor
	 */
	CTimer();

	/**
	 *	@brief resets timer (sets time to zero)
	 */
	void Reset();

	/**
	 *	@brief resets timer (sets time to zero)
	 *	@deprecated This function is deprecated in favor of Reset().
	 */
	inline void ResetTimer()
	{
		Reset();
	}

	/**
	 *	@brief gets time in seconds
	 *	@return Returns time since creation of this object or since
	 *		the last call to ResetTimer(), in seconds.
	 *	@note This should cope nicely with counter overflows.
	 */
	double f_Time() const;

	/**
	 *	@brief gets timer frequency
	 *	@return Returns timer frequency (inverse of the smallest time step).
	 */
	int64_t n_Frequency() const;

protected:
	/**
	 *	@brief gets raw timer sample
	 *	@return Returns raw timer sample, semantic of the value depends on selected timer method.
	 */
	static inline int64_t n_SampleTimer();
};

#endif // __TIMER_INCLUDED
