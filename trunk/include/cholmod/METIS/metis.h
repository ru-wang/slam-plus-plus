/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.1 1998/11/27 17:59:21 karypis Exp $
 */


#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#if !defined(_WIN32) && !defined(_WIN64)
#include <strings.h> // not in MSVC
#endif // !_WIN32 && !_WIN64
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include <defs.h>
#include <struct.h>
#include <macros.h>
#include <rename.h>
#include <proto.h>

