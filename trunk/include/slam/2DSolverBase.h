/*
								+----------------------------------+
								|                                  |
								| *** Base class for 2D solver *** |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2012  |
								|                                  |
								|          2DSolverBase.h          |
								|                                  |
								+----------------------------------+
*/

#pragma once
#ifndef __2D_SOLVER_BASE_INCLUDED
#define __2D_SOLVER_BASE_INCLUDED

/**
 *	@file include/slam/2DSolverBase.h
 *	@brief a simple base class for 2D solver
 *	@author -tHE SWINe-
 *	@date 2012-04-05
 */

#ifndef _USE_MATH_DEFINES
/**
 *	@def _USE_MATH_DEFINES
 *	@brief enables math defines such as M_PI in MSVC
 */
#define _USE_MATH_DEFINES
#endif // _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
//#include "slam/BlockMatrix.h"
#include "eigen/Eigen/Cholesky"

#if !defined(_WIN32) && !defined(_WIN64)
/**
 *	@def _finite
 *	@brief determines if an input value is a finite number
 *	@param[in] x is the input value
 *	@return Returns true if the input value is a finite number, otherwise returns false.
 *	@note This is only compiled on linux where standard library function _finite()
 *		is not available. Always returns true (even for inf or nan).
 */
#define _finite(x) (true)

/**
 *	@def _isnan
 *	@brief determines if an input value is not a number
 *	@param[in] x is the input value
 *	@return Returns true if the input value is not a number, otherwise returns false.
 *	@note This is only compiled on linux where standard library function _isnan()
 *		is not available. Always returns false (even for inf or nan).
 */
#define _isnan(x) (false)
#endif // !_WIN32 && !_WIN64
// we are using some advanced math functions to catch evil numbers; ignore those on linux

/**
 *	@brief implementation of Jacobian calculations, required by 2D solvers
 */
class C2DJacobians {
public:
	/**
	 *	@brief modifies angle so it ends up in the [-2Pi, 2Pi] interval
	 *	@param[in] f_angle is the angle
	 *	@return Returns modulo of the angle inside the [-2Pi, 2Pi] interval.
	 */
	static double f_ClampAngle_2Pi(double f_angle)
	{
#if 0
		if(f_angle > M_PI * 100 || f_angle < -M_PI * 100 || _isnan(f_angle)) // could simply use if(_finite(f_angle)) ... that checks for both inf and nan
			return 0;
		// handle inf, nan

#ifdef _DEBUG
		double f_original = f_angle;
#endif // _DEBUG
		while(f_angle >= M_PI * 2)
			f_angle -= M_PI * 2;
		while(f_angle <= -M_PI * 2)
			f_angle += M_PI * 2;
#ifdef _DEBUG
		_ASSERTE(fabs(fmod(f_original, M_PI * 2) - f_angle) < 1e-5); // t_odo - use modulo instead
#endif // _DEBUG
		// fmod

		return f_angle;
#else // 0
		return (_finite(f_angle))? fmod(f_angle, M_PI * 2) : 0; // that's all there is to it
#endif // 0
	}

	/**
	 *	@brief finds minimum of absolute values
	 *
	 *	@param[in] f_a is the first value
	 *	@param[in] f_b is the second value
	 *	@param[in] f_c is the third value
	 *
	 *	@return Returns f_a if |f_a| < |f_b| and |f_a| < |f_c|,
	 *		returns f_b if |f_b| < |f_a| and |f_b| < |f_c|, otherwise returns f_c.
	 */
	static double f_MinimumAbsolute_3(double f_a, double f_b, double f_c)
	{
		double f_min_abs_a_b = (fabs(f_a) < fabs(f_b))? f_a : f_b;
		return (fabs(f_min_abs_a_b) < fabs(f_c))? f_min_abs_a_b : f_c;
	}

	//! fixup the error so the absolute value is the lowest possible (considering it is modulo 2pi)
	static double f_ClampAngularError_2Pi(double f_error)
	{
		f_error = f_ClampAngle_2Pi(f_error);
		return f_MinimumAbsolute_3(f_error, f_error - 2 * M_PI, f_error + 2 * M_PI);
	}

	/**
	 *	@brief converts xyt coordinates from relative measurement to absolute measurement
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second vertex, relative to the first one
	 *	@param[out] r_t_dest is filled with absolute coordinates of the second vertex
	 */
	static void Relative_to_Absolute(const Eigen::Vector3d &r_t_vertex1,
		const Eigen::Vector3d &r_t_vertex2, Eigen::Vector3d &r_t_dest)
	{
		double p1e = r_t_vertex1(0);
		double p1n = r_t_vertex1(1);
		double p1a = r_t_vertex1(2);

		double p2e = r_t_vertex2(0);
		double p2n = r_t_vertex2(1);
		double p2a = r_t_vertex2(2);

		//double pre, prn, o, co, so, p3e, p3n, p3a;

		/**
		syms p1e p1n p1a p2e p2n p2a real
		 */

		double o = p1a;
		double co = cos(o);
		double so = sin(o);

		double pre = co * p2e - so * p2n;
		double prn = so * p2e + co * p2n;

		double p3e = p1e + pre;
		double p3n = p1n + prn;
		double p3a = p1a + p2a;

		/**
		p1 = [p1e,p1n,p1a]';
		p2 = [p2e,p2n,p2a]';
		p3 = [p3e,p3n,p3a]';
		simplify(jacobian(p3,p1))
		simplify(jacobian(p3,p2))

		ans =
		[ 1, 0, - p2n*cos(p1a) - p2e*sin(p1a)]
		[ 0, 1,   p2e*cos(p1a) - p2n*sin(p1a)]
		[ 0, 0,                             1]
		ans =
		[ cos(p1a), -sin(p1a), 0]
		[ sin(p1a),  cos(p1a), 0]
		[        0,         0, 1]
		 */

		p3a = f_ClampAngle_2Pi(p3a);

		r_t_dest(0) = p3e;
		r_t_dest(1) = p3n;
		r_t_dest(2) = p3a;
		// write the result
	}

	/**
	 *	@brief converts xyt coordinates from absolute measurement to relative measurement
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 */
	template <class _TyDestVector> // want to be able to call this with differnent dest types (sometimes generic VectorXd, sometimes with Vector3d)
	static void Absolute_to_Relative(const Eigen::Vector3d &r_t_vertex1,
		const Eigen::Vector3d &r_t_vertex2, _TyDestVector &r_t_dest)
	{
		double p1e = r_t_vertex1(0);
		double p1n = r_t_vertex1(1);
		double p1a = r_t_vertex1(2);

		double p2e = r_t_vertex2(0);
		double p2n = r_t_vertex2(1);
		double p2a = r_t_vertex2(2);

		/*
		function Pr=Absolute2Relative(p1,p2)
		% computes the relative position of p2 in coordinates of p1.
		%p1,p2 must be column vectors and Pr is column vector

		d(1:2,1) = p2(1:2) - p1(1:2);
		d(3,1) = pi2pi(p2(3) - p1(3));
		%d=p2-p1;
		o=-p1(3);
		R=[[ cos(o) -sin(o) 0];
		   [ sin(o) cos(o) 0];
		   [ 0 0 1]];
		Pr=R*d;
		*/

		double de = p2e - p1e;
		double dn = p2n - p1n;
		double da = p2a - p1a;

		double o = -p1a;
		double co = cos(o);
		double so = sin(o);

		double prf = co * de - so * dn;
		double prl = so * de + co * dn;
		double pra = da;

		/**
		p1 = [p1e,p1n,p1a]';
		p2 = [p2e,p2n,p2a]';
		pr = [prf,prl,pra]';
		simplify(jacobian(pr,p1))
		simplify(jacobian(pr,p2))

		ans =
		[ -cos(p1a), -sin(p1a), sin(p1a)*(p1e - p2e) - cos(p1a)*(p1n - p2n)]
		[  sin(p1a), -cos(p1a), cos(p1a)*(p1e - p2e) + sin(p1a)*(p1n - p2n)]
		[         0,         0,                                          -1]
		ans =
		[  cos(p1a), sin(p1a), 0]
		[ -sin(p1a), cos(p1a), 0]
		[         0,        0, 1]
		 */

		pra = f_ClampAngle_2Pi(pra);

		r_t_dest(0) = prf;
		r_t_dest(1) = prl;
		r_t_dest(2) = pra;
	}

	/**
	 *	@brief converts xyt coordinates from absolute measurement to relative measurement
	 *		and calculates the jacobians
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 *	@param[out] r_t_pose3_pose1 is filled with the first jacobian
	 *	@param[out] r_t_pose3_pose2 is filled with the second jacobian
	 */
	template <class _TyDestVector, class _TyDestMatrix0, class _TyDestMatrix1> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d / Matrix3d)
	static void Absolute_to_Relative(const Eigen::Vector3d &r_t_vertex1,
		const Eigen::Vector3d &r_t_vertex2, _TyDestVector &r_t_dest,
		_TyDestMatrix0 &r_t_pose3_pose1, _TyDestMatrix1 &r_t_pose3_pose2)
	{
		double p1e = r_t_vertex1(0);
		double p1n = r_t_vertex1(1);
		double p1a = r_t_vertex1(2);

		double p2e = r_t_vertex2(0);
		double p2n = r_t_vertex2(1);
		double p2a = r_t_vertex2(2);

		double de = p2e - p1e;
		double dn = p2n - p1n;
		double da = p2a - p1a;

		double o = -p1a;
		double co = cos(o);
		double so = sin(o);

		double prf = co * de - so * dn;
		double prl = so * de + co * dn;
		double pra = da;

		pra = f_ClampAngle_2Pi(pra);

		r_t_dest(0) = prf;
		r_t_dest(1) = prl;
		r_t_dest(2) = pra;

		double cp1a = cos(p1a);
		double sp1a = sin(p1a);

		{
			_TyDestMatrix0 &M = r_t_pose3_pose1;
			M(0, 0) = -cp1a;	M(0, 1) = -sp1a;	M(0, 2) = sp1a * (p1e - p2e) - cp1a * (p1n - p2n);
			M(1, 0) =  sp1a;	M(1, 1) = -cp1a;	M(1, 2) = cp1a * (p1e - p2e) + sp1a * (p1n - p2n);
			M(2, 0) =  0;		M(2, 1) =  0;		M(2, 2) = -1;
		}
		{
			_TyDestMatrix1 &M = r_t_pose3_pose2;
			M(0, 0) =  cp1a;	M(0, 1) = sp1a;		M(0, 2) = 0;
			M(1, 0) = -sp1a;	M(1, 1) = cp1a;		M(1, 2) = 0;
			M(2, 0) =  0;		M(2, 1) = 0;		M(2, 2) = 1;
		}
	}

	/**
	 *	@brief converts range bearing coordinates from absolute measurement
	 *		to relative measurement and calculates the jacobians
	 *
	 *	@param[in] r_v_pose is the first vertex (pose), in absolute coordinates
	 *	@param[in] r_v_landmark is the second vertex (landmark), also in absolute coordinates
	 *	@param[out] r_v_observation is filled with relative coordinates of the second vertex (the landmark)
	 *	@param[out] r_t_observation_pose is filled with the first jacobian
	 *	@param[out] r_t_observation_landmark is filled with the second jacobian
	 */
	template <class _TyDestVector, class _TyDestMatrix0, class _TyDestMatrix1> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d, Vector2d / Matrix2x3d, Matrix2x2d)
	static void Observation2D_RangeBearing(const Eigen::Vector3d &r_v_pose, const Eigen::Vector2d &r_v_landmark,
		_TyDestVector &r_v_observation, _TyDestMatrix0 &r_t_observation_pose, _TyDestMatrix1 &r_t_observation_landmark) // t_odo - mind the convention
	{
		double pe = r_v_pose(0);
		double pn = r_v_pose(1);
		double pa = r_v_pose(2);

		double le = r_v_landmark(0);
		double ln = r_v_landmark(1);

		double de = le - pe;
		double dn = ln - pn;

		double obsr = sqrt(de * de + dn * dn);
		double obsb = f_ClampAngle_2Pi(atan2(dn, de) - pa);

		/**
		  v=pt-P1(1:2);
		  d=norm(v);
		  H1=[-v(1)/d     -v(2)/d      0;
			   v(2)/d^2   -v(1)/d^2   -1 ];
		  H2=[ v(1)/d    v(2)/d;
			  -v(2)/d^2  v(1)/d^2 ];
		 */

		if(fabs(obsr) < 1e-5)
			obsr = 1e-5;
		// @todo - what the hell? are we having numerical issues?

		r_v_observation(0) = obsr;
		r_v_observation(1) = obsb;

		double v1 = de;
		double v2 = dn;
		double d  = obsr;
		double d2 = d * d;
		{
			_TyDestMatrix0 &M = r_t_observation_pose;
			M(0, 0) = -v1 / d;		M(0, 1) = -v2 / d;		M(0, 2) = 0;
			M(1, 0) =  v2 / d2;		M(1, 1) = -v1 / d2;		M(1, 2) = -1; // hell, it *is* 3 by 2
		}
		{
			_TyDestMatrix1 &M = r_t_observation_landmark;
			M(0, 0) =  v1 / d;		M(0, 1) = v2 / d;		/*M(0, 2) = 0;*/
			M(1, 0) = -v2 / d2;		M(1, 1) = v1 / d2;		/*M(1, 2) = 0;*/ // t_odo - this is just a fixup; remove it (should be 2x2)
		}
	}
};

#endif // __2D_SOLVER_BASE_INCLUDED