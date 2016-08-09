/*
								+----------------------------------+
								|                                  |
								| *** Base class for 3D solver *** |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2013  |
								|                                  |
								|          3DSolverBase.h          |
								|                                  |
								+----------------------------------+
*/

#pragma once
#ifndef __3D_SOLVER_BASE_INCLUDED
#define __3D_SOLVER_BASE_INCLUDED

/**
 *	@file include/slam/3DSolverBase.h
 *	@brief a simple base class for 3D solver, made according to 2D solver
 *	@author soso
 *	@date 2013-01-14
 *
 *	@date 2013-01-28
 *
 *	Surrounded some old stuff with ifdefs to speed up the compilation
 *	(basically vertex and edge structures based on the legacy code,
 *	which would be unused in the new code).
 *
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
#include "slam/BlockMatrix.h"
#include "eigen/Eigen/Cholesky"
#include "eigen/Eigen/Geometry" // quaternions

/**
 *	@brief calculates 3D jacobians of difference (u [-] v) using code exported from Matlab
 *
 *	@param[out] Ju is jacobian for the first vertex
 *	@param[out] Jv is jacobian for the second vertex
 *	@param[in] theta is rotation angle of the first vertex
 *	@param[in] omega is rotation angle of the second vertex
 *	@param[in] ux is position of the first vertex
 *	@param[in] uy is position of the first vertex
 *	@param[in] uz is position of the first vertex
 *	@param[in] ua is axis-angle rotation of the first vertex
 *	@param[in] ub is axis-angle rotation of the first vertex
 *	@param[in] uc is axis-angle rotation of the first vertex
 *	@param[in] vx is position of the second vertex
 *	@param[in] vy is position of the second vertex
 *	@param[in] vz is position of the second vertex
 *	@param[in] va is axis-angle rotation of the second vertex
 *	@param[in] vb is axis-angle rotation of the second vertex
 *	@param[in] vc is axis-angle rotation of the second vertex
 *
 *	@note The code is in 3DJacs.cpp, which may not be a part of the distribution.
 */
void Calculate3DJacobians(double Ju[6][6], double Jv[6][6], double theta, double omega,
						  double ux, double uy, double uz, double ua, double ub, double uc,
						  double vx, double vy, double vz, double va, double vb, double vc);

/**
 *	@brief calculates 3D jacobian of unitary plus (u [+] epsilon) using code exported from Matlab
 *
 *	@param[out] Ju is jacobian for the first vertex
 *	@param[in] ux is position of the first vertex
 *	@param[in] uy is position of the first vertex
 *	@param[in] uz is position of the first vertex
 *	@param[in] ua is axis-angle rotation of the first vertex
 *	@param[in] ub is axis-angle rotation of the first vertex
 *	@param[in] uc is axis-angle rotation of the first vertex
 *
 *	@note The code is in 3DJacs.cpp, which may not be a part of the distribution.
 */
void Calculate3DJacobians_Plus(double Ju[6][6], double ux, double uy, double uz, double ua, double ub, double uc);

/**
 *	@brief calculates difference using code exported from Matlab
 *
 *	@param[out] d is difference between the two poses (u [-] v)
 *	@param[in] theta is rotation angle of the first vertex
 *	@param[in] omega is rotation angle of the second vertex
 *	@param[in] ux is position of the first vertex
 *	@param[in] uy is position of the first vertex
 *	@param[in] uz is position of the first vertex
 *	@param[in] ua is axis-angle rotation of the first vertex
 *	@param[in] ub is axis-angle rotation of the first vertex
 *	@param[in] uc is axis-angle rotation of the first vertex
 *	@param[in] vx is position of the second vertex
 *	@param[in] vy is position of the second vertex
 *	@param[in] vz is position of the second vertex
 *	@param[in] va is axis-angle rotation of the second vertex
 *	@param[in] vb is axis-angle rotation of the second vertex
 *	@param[in] vc is axis-angle rotation of the second vertex
 *
 *	@note The code is in 3DJacs.cpp, which may not be a part of the distribution.
 */
void CalculateD(double d[6][1], double theta, double omega,
				double ux, double uy, double uz, double ua, double ub, double uc,
				double vx, double vy, double vz, double va, double vb, double vc);

/**
 *	@brief implementation of Jacobian calculations, required by 3D solvers
 */
class C3DJacobians {
public:
	typedef Eigen::Matrix<double, 6, 1> Vector6d; /**< @brief 6D vector type */
	typedef Eigen::Matrix<double, 6, 6> Matrix6d; /**< @brief 6x6 matrix type */

public:
	/**
	 *	@brief converts from axis angle rep to rotation matrix
	 *	@param[in] v_vec is axis-angle rotation (angle is encoded as the magnitude of the axis)
	 *	@return Returns rotation matrix, corresponding to the input.
	 *	@deprecated This function is deprecated in favor of C3DJacobians::t_AxisAngle_to_RotMatrix(). Please, do not use it.
	 */
	static Eigen::Matrix3d Operator_rot(Eigen::Vector3d v_vec)
	{
		return t_AxisAngle_to_RotMatrix(v_vec); // no need for 2nd implementation
	}

	/**
	 *	@brief converts from axis angle representation to rotation matrix
	 *	@tparam Derived is Eigen derived matrix type for the first matrix argument
	 *	@param[in] v_vec is axis-angle rotation (angle is encoded as the magnitude of the axis)
	 *	@return Returns rotation matrix, corresponding to the input.
	 *	@note Consider whether converting to a quaternion would not be faster: converting to a quaternion
	 *		and transforming a single point will be faster than converting to a rotation matrix. Only when
	 *		transforming two or more points, rotation matrix will be faster.
	 */
	template <class Derived>
	static Eigen::Matrix3d t_AxisAngle_to_RotMatrix(const Eigen::MatrixBase<Derived> &v_vec)
	{
		DimensionCheck<Eigen::Vector3d>(v_vec);

#if 1 // reuse conversion to quaternions, to avoid repeating code. recently improved numerical stability there.
		Eigen::Quaternion<double> rot;
		AxisAngle_to_Quat(v_vec, rot);
		return rot.toRotationMatrix();
#else // 1
		double f_angle = v_vec.norm();//sqrt(x*x + y*y + z*z); // SSE
		// this is really angle in radians

		if(f_angle > 0) {
			v_vec /= f_angle; // SSE
			//x = x / f_angle;
			//y = y / f_angle;
			//z = z / f_angle;
		} else
			return Eigen::Matrix3d::Identity(); // no rotation, save some hassle below
		// normalize the axis

#if 1
		double f_half_angle = f_angle * .5;
		v_vec *= sin(f_half_angle); // SSE
		Eigen::Quaternion<double> rot(cos(f_half_angle), v_vec(0), v_vec(1), v_vec(2));
		return rot.toRotationMatrix();
		// use quaternion, Eigen might provide SSE accelerated conversion to matrix
#else // 1
		double s = sin(f_angle);
		double c = cos(f_angle);
		double v = 1 - c;

		double x = v_vec(0), y = v_vec(1), z = v_vec(2);
		double xyv = x * y * v;
		double yzv = y * z * v;
		double xzv = x * z * v;

		Eigen::Matrix3d Q;
		Q << x * x * v + c, xyv - z * s, xzv + y * s,
			 xyv + z * s, y * y * v + c, yzv - x * s,
			 xzv - y * s, yzv + x * s, z * z * v + c;

		return Q;
#endif // 1
#endif // 1
	}

	/**
	 *	@brief converts from rotation matrix rep to axis angle
	 *	@param[in] Q is the rotation matrix
	 *	@return Returns axis-amgle rotation, where the angle is encoded as magnitude of the axis.
	 *	@deprecated This function is deprecated in favor of C3DJacobians::v_RotMatrix_to_AxisAngle(). Please, do not use it.
	 */
	static Eigen::Vector3d Operator_arot(const Eigen::Matrix3d &Q) // note that Eigen probably implements this somewhere
	{
		return v_RotMatrix_to_AxisAngle(Q); // no need for 2nd implementation
	}

	/**
	 *	@brief converts from a rotation matrix to axis anglerepresentation
	 *	@tparam Derived is Eigen derived matrix type for the first matrix argument
	 *	@param[in] Q is the rotation matrix
	 *	@return Returns axis-amgle rotation, where the angle is encoded as magnitude of the axis.
	 */
	template <class Derived>
	static Eigen::Vector3d v_RotMatrix_to_AxisAngle(const Eigen::MatrixBase<Derived> &Q)
	{
		DimensionCheck<Eigen::Matrix3d>(Q);

#if 1 // reuse conversion to quaternions, to avoid repeating code. recently improved numerical stability there.
		Eigen::Quaternion<double> quat(Q); // converts rotation matrix to quaternion (hopefully SSE)
		Eigen::Vector3d v_aang;
		Quat_to_AxisAngle(quat, v_aang); // reuse existing code
		return v_aang;
#else // 1
#if 1
		Eigen::Quaternion<double> quat(Q); // converts rotation matrix to quaternion (hopefully SSE)
		double f_half_angle = acos(quat.w());
		if(f_half_angle == 0)
			return Eigen::Vector3d(0, 0, 0); // handle zero angle rotation
		return quat.vec() * (2 * f_half_angle / sin(f_half_angle)); // SSE
		// need to divide by sine of half angle to get the normalized rotation axis,
		// then multiply by angle to get axis-angle
#else // 1
		const double epsilon = 1e-12;

		Eigen::Vector4d r;
		Eigen::Vector3d axis;
		double matrace = Q.trace();

		if(fabs(matrace - 3.0) <= epsilon)
			r << 0, 1, 0, 0;
		else if(fabs(matrace + 1.0) <= epsilon) {
			printf("fabs(matrace + 1.0) <= epsilon\n");

			for(int a = 0; a < 3; a ++) {
				double f_axis_a = (Q(a, a) + 1) / 2;
				axis(a) = sqrt((f_axis_a > 0)? f_axis_a : 0);
				axis(a) = (axis(a) > epsilon)? axis(a) : 0; // clamp values below epsilon
			}

			//flipping
			Eigen::Vector3d m_upper(Q(1,2), Q(0,2), Q(0,1)), signs;

			for(int a = 0; a < 3; a ++)
				signs(a) = ((m_upper(a) > epsilon)? 1 : ((m_upper(a) < -epsilon)? -1 : 0)); // same amount of comparisons, one multiplication fewer
				//signs(a) = ((m_upper(a) > 0)? 1 : ((m_upper(a) == 0)? 0 : -1)) * (fabs(m_upper(a)) > epsilon);

			Eigen::Vector3d flip; // delay creation of new objects to the latest possible moment
			if(signs.sum() >= 0) // signs(0) + signs(1) + signs(2) >= 0 // SSE
				flip << 1, 1, 1;
			else if(signs(0) != 0 && signs(1) != 0 && signs(2) != 0)
				flip = -signs;
			else {
				Eigen::Vector3d shifted(signs(2), signs(0), signs(1)); // limit the scope of the objects where possible
				for(int a = 0; a < 3; a ++)
					flip(a) = shifted(a) + (shifted(a) == 0);
			}

			//flip the axis
			//for(int a = 0; a < 3; a ++)
			//	axis(a) = axis(a) * flip(a);
			//r << axis(0), axis(1), axis(2), M_PI;
			r.head<3>() = axis.cwiseProduct(flip); // elementwise operation using .array() or .cwiseProduct(), accelerated with SSE; also no need to modify axis, it is about to be overwritten
			r(3) = M_PI;
		} else {
			double phi = acos((matrace - 1) / 2.0);
			double den = 2.0 * sin(phi);
			//axis << (Q(2, 1) - Q(1, 2)) / den, (Q(0, 2) - Q(2, 0)) / den, (Q(1, 0) - Q(0, 1)) / den;
			//r << axis(0), axis(1), axis(2), phi;
			r.head<3>() = (Eigen::Vector3d(Q(2, 1), Q(0, 2), Q(1, 0)) -
				Eigen::Vector3d(Q(1, 2), Q(2, 0), Q(0, 1))) / den; // vector ops, SSE again
			r(3) = phi;
		}

		double sum = r.head<3>().norm(); // sqrt(r(0) * r(0) + r(1) * r(1) + r(2) * r(2));
		axis = r.head<3>() * (r(3) / sum);
		//for(int a = 0; a < 3; a ++)
		//	axis(a) = r(a) * (r(3) / sum);
		// scale axis by angle

		return axis;
#endif // 1
#endif // 1
	}

	/**
	 *	@brief converts axis-angle representation to quaternion
	 *
	 *	@tparam Derived is Eigen derived matrix type for the first matrix argument
	 *
	 *	@param[in] r_axis_angle is a vector where unit direction gives axis, and magnitude gives angle in radians
	 *	@param[out] r_quat is quaternion representing the same rotation as r_axis_angle
	 */
	template <class Derived>
	static void AxisAngle_to_Quat(const Eigen::MatrixBase<Derived> &r_axis_angle, Eigen::Quaterniond &r_quat)
	{
		DimensionCheck<Eigen::Vector3d>(r_axis_angle);

		f_AxisAngle_to_Quat(r_axis_angle, r_quat); // ignore the return value, probably optimized away
		/*double f_angle = r_axis_angle.norm();
		if(f_angle < 1e-12)
			r_quat = Eigen::Quaterniond(1, 0, 0, 0); // cos(0) = 1
		else {
			//_ASSERTE(f_angle <= M_PI); // sometimes broken
			f_angle = fmod(f_angle, M_PI * 2);
			double q = (sin(f_angle * .5) / f_angle);
			r_quat = Eigen::Quaterniond(cos(f_angle * .5), r_axis_angle(0) * q,
				r_axis_angle(1) * q, r_axis_angle(2) * q);
			r_quat.normalize();
		}*/
	}

	/**
	 *	@brief converts axis-angle representation to quaternion
	 *
	 *	@tparam Derived is Eigen derived matrix type for the first matrix argument
	 *
	 *	@param[in] r_axis_angle is a vector where unit direction gives axis, and magnitude gives angle in radians
	 *	@param[out] r_quat is quaternion representing the same rotation as r_axis_angle
	 *
	 *	@return Returns rotation angle in radians (after flushing small angles to zero).
	 */
	template <class Derived>
	static double f_AxisAngle_to_Quat(const Eigen::MatrixBase<Derived> &r_axis_angle, Eigen::Quaterniond &r_quat)
	{
		DimensionCheck<Eigen::Vector3d>(r_axis_angle);

		double f_angle = r_axis_angle.norm();
		if(f_angle < 1e-12) { // increasing this does not help
			r_quat = Eigen::Quaterniond(1, 0, 0, 0); // cos(0) = 1
			return 0;//M_PI * 2;
		} else {
			//_ASSERTE(f_angle <= M_PI); // sometimes broken
			double f_half_angle = fmod(f_angle, M_PI * 2) * .5;
			double q = sin(f_half_angle) / f_angle; // sin is [0, 1], angle is about [-2pi, 2pi] but could be slightly more / less due to optimization incrementing / decrementing it
			r_quat = Eigen::Quaterniond(cos(f_half_angle), r_axis_angle(0) * q,
				r_axis_angle(1) * q, r_axis_angle(2) * q);
			r_quat.normalize(); // should already be normalized unles there was a great roundoff in calculation of q; could detect that.
		}
		return f_angle;
	}

	/**
	 *	@brief converts quaternion to axis-angle representation
	 *
	 *	@tparam Derived is Eigen derived matrix type for the second matrix argument
	 *
	 *	@param[in] r_quat is quaternion representing the same rotation as r_axis_angle
	 *	@param[out] r_axis_angle is a vector where unit direction gives axis, and magnitude gives angle in radians
	 */
	template <class Derived>
	static void Quat_to_AxisAngle(const Eigen::Quaterniond &r_quat, Eigen::MatrixBase<Derived> &r_axis_angle)
	{
		DimensionCheck<Eigen::Vector3d>(r_axis_angle);

		f_Quat_to_AxisAngle(r_quat, r_axis_angle); // ignore the return value, probably optimized away
		////double f_half_angle = /*(r_quat.w() <= 0)? asin(r_quat.vec().norm()) :*/ acos(r_quat.w()); // 0 .. pi
		//double f_half_angle = (r_quat.w() == 1.0)? asin(r_quat.vec().norm()) : acos(r_quat.w()); // 0 .. pi // more numerically robust. normalization of quaternions sometimes yields (eps eps eps 1.0), the above line would calculate angle of 0.0, this line calculates angle of ~eps
		//_ASSERTE(f_half_angle >= 0);
		//if(f_half_angle < 1e-12)
		//	r_axis_angle = Eigen::Vector3d(0, 0, 0); // lim(sin(x) / x) for x->0 equals 1, we're therefore multiplying a null vector by 1
		//else {
		//	double f_angle = 2 * ((r_quat.w() <= 0)? f_half_angle - M_PI : f_half_angle);
		//	r_axis_angle = r_quat.vec() * (f_angle / sin(f_half_angle));
		//}
	}

	/**
	 *	@brief converts quaternion to axis-angle representation
	 *
	 *	@tparam Derived is Eigen derived matrix type for the second matrix argument
	 *
	 *	@param[in] r_quat is quaternion representing the same rotation as r_axis_angle
	 *	@param[out] r_axis_angle is a vector where unit direction gives axis, and magnitude gives angle in radians
	 *
	 *	@return Returns rotation angle in radians (after flushing small angles to zero).
	 */
	template <class Derived>
	static double f_Quat_to_AxisAngle(const Eigen::Quaterniond &r_quat, Eigen::MatrixBase<Derived> &r_axis_angle)
	{
		DimensionCheck<Eigen::Vector3d>(r_axis_angle);

		//double f_half_angle = /*(r_quat.w() <= 0)? asin(r_quat.vec().norm()) :*/ acos(r_quat.w()); // 0 .. pi

		double f_half_angle = (/*fabs*/(r_quat.w()) >= 1.0)? asin(r_quat.vec().norm()) : acos(r_quat.w()); // 0 .. pi
		// more numerically robust. normalization of quaternions sometimes yields (eps eps eps 1.0)
		// product of normalized quats sometimes yields (eps eps eps (1.0 + 1 ULP)),
		// the above line would calculate angle of 0.0, this line calculates angle of ~eps
		// note that r_quat.w() doesn't ever seem to be negative one, otherwise enable fabs()

		_ASSERTE(f_half_angle >= 0); // detects NaNs
		if(f_half_angle < 1e-12) {
			r_axis_angle = Eigen::Vector3d(0, 0, 0); // lim(sin(x) / x) for x->0 equals 1, we're therefore multiplying a null vector by 1
			return 0;
		} else {
			double f_angle = 2 * ((r_quat.w() <= 0)? f_half_angle - M_PI : f_half_angle);
			r_axis_angle = r_quat.vec() * (f_angle / sin(f_half_angle));
			return f_angle;
		}
	}

	/**
	 *	@brief operator plus for 3D vertices
	 *
	 *	@param[in] r_t_vertex1 is the vertex to be modified
	 *	@param[in] r_t_vertex2 is the delta vector
	 *	@param[out] r_t_dest is the result of the operation
	 *
	 *	@deprecated This is deprecated in favor of C3DJacobians::Relative_to_Absolute(). Please, do not use this function.
	 */
	static void Smart_Plus(const Vector6d &r_t_vertex1, const Vector6d &r_t_vertex2, Vector6d &r_t_dest)
	{
		Relative_to_Absolute(r_t_vertex1, r_t_vertex2, r_t_dest);
		// 2nd implementation not required
	}

	/**
	 *	@brief composes a pair of relative transformations to yield an aboslute transformation
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first transformation (in absolute coordinates)
	 *	@param[in] r_t_vertex2 is the second transformation (relative to the first one)
	 *	@param[out] r_t_dest is filled with the composed transformation
	 */
	template <class Derived0, class Derived1, class Derived2>
	static inline void Relative_to_Absolute(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector6d>(r_t_vertex2);
		DimensionCheck<Vector6d>(r_t_dest);

#if 0
		/* TODO: make more efficient */

		//r_t_dest.resize(6, 1); // no need to resize, it is 6D vector at compile time
		//r_t_dest(0) = r_t_vertex1(0) + r_t_vertex2(0);
		//r_t_dest(1) = r_t_vertex1(1) + r_t_vertex2(1);
		//r_t_dest(2) = r_t_vertex1(2) + r_t_vertex2(2);
		r_t_dest.head<3>() = r_t_vertex1.head<3>() + r_t_vertex2.head<3>(); // accelerated using SSE

		//sum the rotations
		Eigen::Matrix3d pQ = t_AxisAngle_to_RotMatrix(r_t_vertex1.tail<3>());
		Eigen::Matrix3d dQ = t_AxisAngle_to_RotMatrix(r_t_vertex2.tail<3>());

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // multiplication without intermediate storage
		//Eigen::Vector3d axis = v_RotMatrix_to_AxisAngle(QQ);
		r_t_dest.tail<3>() = v_RotMatrix_to_AxisAngle(QQ);

		//r_t_dest(3) = axis(0);
		//r_t_dest(4) = axis(1);
		//r_t_dest(5) = axis(2);
#else
		//Eigen::Vector3d p1 = r_t_vertex1.template head<3>();
		Eigen::Quaterniond q1;
		AxisAngle_to_Quat(r_t_vertex1.template tail<3>(), q1);

		//Eigen::Vector3d p2 = r_t_vertex2.template head<3>();
		Eigen::Quaterniond q2;
		AxisAngle_to_Quat(r_t_vertex2.template tail<3>(), q2);

		r_t_dest.template head<3>() = r_t_vertex1.template head<3>() +
			q1._transformVector(r_t_vertex2.template head<3>()); // p2 rotates!
		Eigen::Vector3d v_aang;
		Quat_to_AxisAngle((q1 * q2)/*.normalized()*/, v_aang); // product doesn't denormalize
		r_t_dest.template tail<3>() = v_aang;
#endif

#if 0
#if 0
		/* TODO: make more efficient */
		//r_t_dest.resize(6, 1); // no need to resize fixed size expressions
		Eigen::Vector3d p = r_t_vertex1.head<3>();
		Eigen::Matrix3d pQ = t_AxisAngle_to_RotMatrix(r_t_vertex1.tail<3>());
		Eigen::Vector3d d = r_t_vertex2.head<3>();
		Eigen::Matrix3d dQ = t_AxisAngle_to_RotMatrix(r_t_vertex2.tail<3>());
		r_t_dest.head<3>() = p + pQ * d;

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // multiplication without intermediate storage
		r_t_dest.tail<3>() = v_RotMatrix_to_AxisAngle(QQ);
#else // 0
		Eigen::Vector3d p1 = r_t_vertex1.head<3>();
		Eigen::Quaterniond q1;
		AxisAngle_to_Quat(r_t_vertex1.tail<3>(), q1);

		Eigen::Vector3d p2 = r_t_vertex2.head<3>();
		Eigen::Quaterniond q2;
		AxisAngle_to_Quat(r_t_vertex2.tail<3>(), q2);

		r_t_dest.head<3>() = p1 + q1._transformVector(p2);
		Eigen::Vector3d v_aang;
		Quat_to_AxisAngle((q1 * q2)/*.normalized()*/, v_aang); // product doesn't denormalize
		r_t_dest.tail<3>() = v_aang;
#endif // 0
#endif // 0
	}

	/**
	 *	@brief converts xyz axis angle coordinates from absolute measurement to relative measurement
	 *		and calculates the jacobians
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 */
	template <class Derived0, class Derived1, class Derived2>
	static void Absolute_to_Relative(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector6d>(r_t_vertex2);
		DimensionCheck<Vector6d>(r_t_dest);

#if 0
		/* TODO: make more efficient */
		Eigen::Vector3d p1 = r_t_vertex1.template head<3>();
		Eigen::Matrix3d pQ1 = t_AxisAngle_to_RotMatrix(r_t_vertex1.template tail<3>());

		Eigen::Vector3d p2 = r_t_vertex2.template head<3>();
		Eigen::Matrix3d pQ2 = t_AxisAngle_to_RotMatrix(r_t_vertex2.template tail<3>());

		Eigen::Matrix3d pQ1_inv = pQ1.inverse();
		//Eigen::Matrix3d pQ2_inv = pQ2.inverse();

		//r_t_dest.resize(6, 1); // no need to resize fixed size expressions
		r_t_dest.head<3>() = pQ1_inv * (p2 - p1);

		Eigen::Matrix3d QQ = pQ1_inv * pQ2;
		r_t_dest.tail<3>() = v_RotMatrix_to_AxisAngle(QQ);	//TODO: is this right??
#else
		//Eigen::Vector3d p1 = r_t_vertex1.template head<3>();
		Eigen::Quaterniond q1;
		AxisAngle_to_Quat(r_t_vertex1.template tail<3>(), q1);

		//Eigen::Matrix3d pQ1 = t_AxisAngle_to_RotMatrix(r_t_vertex1.tail<3>());
		//Eigen::Matrix3d pQ1_inv = pQ1.inverse();
		// matrix not required after all

		//Eigen::Vector3d p2 = r_t_vertex2.template head<3>();
		Eigen::Quaterniond q2;
		AxisAngle_to_Quat(r_t_vertex2.template tail<3>(), q2);

		Eigen::Quaterniond q1_inv = q1.conjugate(); // the inverse rotation (also have .inverse() but that is not needed)

		r_t_dest.template head<3>() = q1_inv._transformVector(r_t_vertex2.template head<3>() -
			r_t_vertex1.template head<3>()); // this is precise enough
		//r_t_dest.head<3>() = pQ1_inv * (p2 - p1); // or can use matrix to transform position

		//Eigen::Quaterniond prod = (q1_inv * q2)/*.normalized()*/; // quaternion product does not denormalize
		//if(1 || (prod.w()) >= 0) { //printf("bloop %f\n", prod.w());
			Eigen::Vector3d v_aang;
			Quat_to_AxisAngle(q1_inv * q2, v_aang); // use quaternion to calculate the rotation
			r_t_dest.template tail<3>() = v_aang;
		//} else { //printf("bleep %f\n", prod.w());
		//	//r_t_dest.tail<3>() = v_RotMatrix_to_AxisAngle(pQ1_inv * q2); // use quaternion to calculate the rotation
		//	r_t_dest.tail<3>() = v_RotMatrix_to_AxisAngle((q1_inv * q2).toRotationMatrix()); // use quaternion to calculate the rotation
		//}
		// arot is needed, but toRotationMatrix() immediately followed by arot() can likely be optimized
#endif
	}

	/**
	 *	@brief square
	 *	@param[in] x is input argument
	 *	@return Returns x squared.
	 *	@note This is utility function to support Matlab's ccode() directly.
	 */
	static inline double sqr(double x)
	{
		return x * x;
	}

	/**
	 *	@brief conjugate
	 *	@param[in] x is input argument
	 *	@return Returns x (it is real).
	 *	@note This is utility function to support Matlab's ccode() directly.
	 */
	static inline double conjugate(double x)
	{
		return x;
	}

	/**
	 *	@brief magnitude
	 *
	 *	@param[in] x is input argument
	 *	@param[in] y is input argument
	 *
	 *	@return Returns length of the (x, y) vector.
	 *
	 *	@note This is utility function to support Matlab's ccode() directly.
	 */
	static inline double abs(double x, double y) // complex abs
	{
		return sqrt(x * x + y * y);
	}

	/**
	 *	@brief calculates skew-symmetric matrix for a given vector
	 *	@param[in] v is input vector
	 *	@return Returns skew-symmetric matrix of v.
	 */
	static Eigen::Matrix3d v_SkewSym(Eigen::Vector3d v)
	{
		Eigen::Matrix3d m;
		m << 0, -v(2), v(1),
			v(2), 0, -v(0),
			-v(1), v(0), 0;
		return m;
	}

	/**
	 *	@brief converts xyz axis angle coordinates from absolute measurement to relative measurement
	 *		and calculates the jacobians
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *	@tparam Derived3 is Eigen derived matrix type for the fourth matrix argument
	 *	@tparam Derived4 is Eigen derived matrix type for the fifth matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 *	@param[out] r_t_pose3_pose1 is filled with the first jacobian
	 *	@param[out] r_t_pose3_pose2 is filled with the second jacobian
	 */
	template <class Derived0, class Derived1, class Derived2, class Derived3, class Derived4> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d / Matrix3d)
	static void Absolute_to_Relative(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest,
		Eigen::MatrixBase<Derived3> &r_t_pose3_pose1, Eigen::MatrixBase<Derived4> &r_t_pose3_pose2)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector6d>(r_t_vertex2);
		DimensionCheck<Vector6d>(r_t_dest);
		DimensionCheck<Matrix6d>(r_t_pose3_pose1);
		DimensionCheck<Matrix6d>(r_t_pose3_pose2);

		typedef Eigen::MatrixBase<Derived3> _TyDestMatrix0;
		typedef Eigen::MatrixBase<Derived4> _TyDestMatrix1;
		// used on multiple places below

#if 0
		// use the magical analytical jacobians (magobians)

		Eigen::Vector3d p1 = r_t_vertex1.head<3>();
		Eigen::Quaterniond q1;
		double f_theta = f_AxisAngle_to_Quat(r_t_vertex1.tail<3>(), q1);
		_ASSERTE(f_theta >= 0);

		Eigen::Vector3d p2 = r_t_vertex2.head<3>();
		Eigen::Quaterniond q2;
		double f_omega = f_AxisAngle_to_Quat(r_t_vertex2.tail<3>(), q2);
		_ASSERTE(f_omega >= 0);

		Eigen::Quaterniond q1_inv = q1.conjugate(); // the inverse rotation (also have .inverse() but that is not needed)

		r_t_dest.template head<3>() = q1_inv._transformVector(p2 - p1); // this is precise enough

		Eigen::Vector3d v_aang;
		Eigen::Quaterniond qdest = (q1_inv * q2)/*.normalized()*/; // quaternion product does not denormalize
		f_Quat_to_AxisAngle(qdest, v_aang); // use quaternion to calculate the rotation
		r_t_dest.template tail<3>() = v_aang;

		double ux = r_t_vertex1(0);
		double uy = r_t_vertex1(1);
		double uz = r_t_vertex1(2);
		double ua = r_t_vertex1(3);
		double ub = r_t_vertex1(4);
		double uc = r_t_vertex1(5);
		double vx = r_t_vertex2(0);
		double vy = r_t_vertex2(1);
		double vz = r_t_vertex2(2);
		double va = r_t_vertex2(3);
		double vb = r_t_vertex2(4);
		double vc = r_t_vertex2(5);
		// just rename the variables

#if 0
		double thetaMLVar = f_theta, sintheta = f_sin_theta, costheta = f_cos_theta,
			omega = f_omega, cosomega = f_cos_omega, sinomega = f_sin_omega;
		if(omega < 0.1 || thetaMLVar < 0.1 || /*fabs(omega - M_PI) < 0.1 || fabs(thetaMLVar - M_PI) < 0.1 ||*/ omega >= 2 * M_PI - 0.1 || thetaMLVar >= 2 * M_PI - 0.1) {
			// use the numerical jacobians
			/* TODO: make more efficient */

			//lets try it according to g2o
			const double delta = 1e-9;
			const double scalar = 1.0 / (delta);


			Matrix6d Eps;// = delta * Eigen::MatrixXd::Identity(6, 6); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
			Eps = Matrix6d::Identity() * delta; // faster, all memory on stack

			_TyDestMatrix0 &H1 = r_t_pose3_pose1;
			_TyDestMatrix1 &H2 = r_t_pose3_pose2;
			_ASSERTE(H1.rows() == 6 && H1.cols() == 6 && H2.rows() == 6 && H2.cols() == 6); // just make sure the shape is right
			// can actually work inplace

			//Absolute_to_Relative(r_t_vertex1, r_t_vertex2, r_t_dest);
			for(int j = 0; j < 6; ++ j) {
				Vector6d d1, p_delta;
				Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
				Absolute_to_Relative(p_delta, r_t_vertex2, d1);
				H1.col(j) = (d1 - r_t_dest) * scalar;

				Vector6d d2;
				Relative_to_Absolute(r_t_vertex2, Eps.col(j), p_delta);
				Absolute_to_Relative(r_t_vertex1, p_delta, d2);
				H2.col(j) = (d2 - r_t_dest) * scalar;
			}
			// return jacobs

			return;
		}
		// handle low angles using numerical jacobians (avoid division by zero in the analytical equations)
#endif // 0

		double Ju[6][6] = {0}, Jv[6][6] = {0};
		Calculate3DJacobians(Ju, Jv, f_theta, f_omega, 
			ux, uy, uz, ua, ub, uc, vx, vy, vz, va, vb, vc);
		// todo - simplify plus, actually only interested in

		_TyDestMatrix0 &ju = r_t_pose3_pose1;
		_TyDestMatrix1 &jv = r_t_pose3_pose2;
		for(int i = 0; i < 6; ++ i) {
			for(int j = 0; j < 6; ++ j) {
				ju(i, j) = Ju[i][j];
				jv(i, j) = Jv[i][j];
			}
		}
		// get the data out

		/*FILE *p_fw = fopen("jacobs3d_anal_manifold.m", "a");
		if(p_fw) {
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, r_t_vertex1, "u = ", ";\n");
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, r_t_vertex2, "v = ", ";\n");
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, r_t_dest, "d_gt = ", "; % difference between u and v\n");
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, ju, "H1_gt = ", "; % the first jacobian\n");
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, jv, "H2_gt = ", "; % the second jacobian\n");
			fprintf(p_fw, "%s", "\n% ---\n\n");
			fclose(p_fw);
		}*/
		// get some "ground truth" for analytical jacobian debugging

#ifdef _DEBUG
		Matrix6d num_H1;
		Matrix6d num_H2;
		Matrix6d num_H1_skew;
		Matrix6d num_H2_skew;
		// can actually work inplace

		{
			const double delta = 1e-9;
			const double scalar = 1.0 / (delta);

			Matrix6d Eps;// = delta * Eigen::MatrixXd::Identity(6, 6); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
			Eps = Matrix6d::Identity() * delta; // faster, all memory on stack

			//Absolute_to_Relative(r_t_vertex1, r_t_vertex2, r_t_dest); // already did that
			for(int j = 0; j < 6; ++ j) {
				{
					Vector6d d1, p_delta;
					Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
					//p_delta = r_t_vertex1 + Eps.col(j);
					Absolute_to_Relative(p_delta, r_t_vertex2, d1);
					num_H1.col(j) = (d1 - r_t_dest) * scalar;

					Vector6d d2;
					Relative_to_Absolute(r_t_vertex2, Eps.col(j), p_delta);
					//p_delta = r_t_vertex2 + Eps.col(j);
					Absolute_to_Relative(r_t_vertex1, p_delta, d2);
					num_H2.col(j) = (d2 - r_t_dest) * scalar;
				}
				// proper jacobians

				/*{
					Vector6d d1, p_delta;
					//Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
					p_delta = r_t_vertex1 + Eps.col(j);
					Absolute_to_Relative(p_delta, r_t_vertex2, d1);
					num_H1_skew.col(j) = (d1 - r_t_dest) * scalar;

					Vector6d d2;
					//Relative_to_Absolute(r_t_vertex2, Eps.col(j), p_delta);
					p_delta = r_t_vertex2 + Eps.col(j);
					Absolute_to_Relative(r_t_vertex1, p_delta, d2);
					num_H2_skew.col(j) = (d2 - r_t_dest) * scalar;
				}*/
				// jacobians missing [+] at the input (skewed)
			}
			// return jacobs
		}

		double d[6][1] = {0};
		CalculateD(d, thetaMLVar, omega, 
			ux, uy, uz, ua, ub, uc,
			vx, vy, vz, va, vb, vc);
		Vector6d dv;
		for(int i = 0; i < 6; ++ i)
			dv(i, 0) = d[i][0];
		// convert to a vector

		double f_norm_vertex = (dv - r_t_dest).norm();
		_ASSERTE(f_norm_vertex < 1e-10);

		double f_norm1 = (ju - num_H1/*_skew*/).norm();
		double f_norm2 = (jv - num_H2/*_skew*/).norm();
		double f_norm1p = (ju * jup - num_H1).norm();
		double f_norm2p = (jv * jvp - num_H2).norm();
		double f_normp1 = (jup * ju - num_H1).norm();
		double f_normp2 = (jvp * jv - num_H2).norm();
		/*double f_norm1t = (ju.transpose() - num_H1).norm();
		double f_norm2t = (jv.transpose() - num_H2).norm();*/
		/*double f_norm11 = (ju.block<3, 3>(0, 0) - num_H1.block<3, 3>(0, 0)).norm();
		double f_norm12 = (ju.block<3, 3>(3, 0) - num_H1.block<3, 3>(3, 0)).norm();
		double f_norm13 = (ju.block<3, 3>(0, 3) - num_H1.block<3, 3>(0, 3)).norm();
		double f_norm14 = (ju.block<3, 3>(3, 3) - num_H1.block<3, 3>(3, 3)).norm();
		double f_norm21 = (jv.block<3, 3>(0, 0) - num_H2.block<3, 3>(0, 0)).norm();
		double f_norm22 = (jv.block<3, 3>(3, 0) - num_H2.block<3, 3>(3, 0)).norm();
		double f_norm23 = (jv.block<3, 3>(0, 3) - num_H2.block<3, 3>(0, 3)).norm();
		double f_norm24 = (jv.block<3, 3>(3, 3) - num_H2.block<3, 3>(3, 3)).norm();*/
		/*double f_norm11t = (ju.transpose().block<3, 3>(0, 0) - num_H1.block<3, 3>(0, 0)).norm();
		double f_norm12t = (ju.transpose().block<3, 3>(3, 0) - num_H1.block<3, 3>(3, 0)).norm();
		double f_norm13t = (ju.transpose().block<3, 3>(0, 3) - num_H1.block<3, 3>(0, 3)).norm();
		double f_norm14t = (ju.transpose().block<3, 3>(3, 3) - num_H1.block<3, 3>(3, 3)).norm();
		double f_norm21t = (jv.transpose().block<3, 3>(0, 0) - num_H2.block<3, 3>(0, 0)).norm();
		double f_norm22t = (jv.transpose().block<3, 3>(3, 0) - num_H2.block<3, 3>(3, 0)).norm();
		double f_norm23t = (jv.transpose().block<3, 3>(0, 3) - num_H2.block<3, 3>(0, 3)).norm();
		double f_norm24t = (jv.transpose().block<3, 3>(3, 3) - num_H2.block<3, 3>(3, 3)).norm();*/
		double f_worse = std::max(f_norm1, f_norm2); // useless, just let the debugger break after all is calculated
#endif // _DEBUG
#elif 0
		// use the new Hauke jacobians (not working at the moment)

		Eigen::Vector3d tu = r_t_vertex1.template head<3>();
		Eigen::Matrix3d Ru = t_AxisAngle_to_RotMatrix(r_t_vertex1.template tail<3>());
		Eigen::Vector3d tv = r_t_vertex2.template head<3>();
		Eigen::Matrix3d Rv = t_AxisAngle_to_RotMatrix(r_t_vertex2.template tail<3>());
		// call RotMat

		Eigen::Vector3d td = Ru.transpose() * (tv - tu);
		Eigen::Matrix3d Rd = Ru.transpose() * Rv;

		/*Matrix6d Mu, Mv, Md;

		Mu.block<3, 3>(0, 0) = -v_SkewSym(r_t_vertex1.tail<3>());
		Mu.block<3, 3>(0, 3) = -v_SkewSym(r_t_vertex1.head<3>());
		Mu.block<3, 3>(3, 0).setZero();
		Mu.block<3, 3>(3, 3) = -v_SkewSym(r_t_vertex1.tail<3>());

		Mv.block<3, 3>(0, 0) = v_SkewSym(r_t_vertex2.tail<3>());
		Mv.block<3, 3>(0, 3) = v_SkewSym(r_t_vertex2.head<3>());
		Mv.block<3, 3>(3, 0).setZero();
		Mv.block<3, 3>(3, 3) = v_SkewSym(r_t_vertex2.tail<3>());

		Md.block<3, 3>(0, 0) = Rd;
		Md.block<3, 3>(0, 3) = v_SkewSym(td) * Rd;
		Md.block<3, 3>(3, 0).setZero();
		Md.block<3, 3>(3, 3) = Rd;

		if(0) {
			std::cout << "Mu =" << Mu << std::endl;
			std::cout << "Mv =" << Mv << std::endl;
			std::cout << "Md =" << Md << std::endl;
		}
		// debug

		Matrix6d I6;
		I6.setIdentity();

		_TyDestMatrix0 &H1 = r_t_pose3_pose1;
		_TyDestMatrix1 &H2 = r_t_pose3_pose2;

		H1 = (I6 + Mu * .5) * Md;
		H2 = -(I6 + Mv * .5);*/
		// eval jacobians

		Matrix6d Mu, Mv;

		r_t_dest.template head<3>() = td;
		r_t_dest.template tail<3>() = v_RotMatrix_to_AxisAngle(Rd);
		// calculate distance

		_TyDestVector &d = r_t_dest; // rename

		Mu.template block<3, 3>(0, 0) = -v_SkewSym(d.template tail<3>());
		Mu.template block<3, 3>(0, 3) = -v_SkewSym(d.template head<3>());
		Mu.template block<3, 3>(3, 0).setZero();
		Mu.template block<3, 3>(3, 3) = -v_SkewSym(d.template tail<3>());

		Mv.template block<3, 3>(0, 0) = v_SkewSym(d.template tail<3>());
		Mv.template block<3, 3>(0, 3) = v_SkewSym(d.template head<3>());
		Mv.template block<3, 3>(3, 0).setZero();
		Mv.template block<3, 3>(3, 3) = v_SkewSym(d.template tail<3>());

		if(0) {
			std::cout << "Mu =" << Mu << std::endl;
			std::cout << "Mv =" << Mv << std::endl;
		}
		// debug

		Matrix6d I6;
		I6.setIdentity();

		_TyDestMatrix0 &H1 = r_t_pose3_pose1;
		_TyDestMatrix1 &H2 = r_t_pose3_pose2;

		H1 = -(I6 + Mu * .5);
		H2 = (I6 + Mv * .5);
		// eval jacobians
#else // 0
		// use the numerical jacobians
		/* TODO: make more efficient */

		//lets try it according to g2o
		const double delta = 1e-9;
		const double scalar = 1.0 / (delta);

		//jacobians have to be computed here
		//double eps = 1.e-9; //TODO: what is this???
		//double eps_ = sqrt(eps);
		//double eps_2 = 7.45e-9;

		Matrix6d Eps;// = delta * Eigen::MatrixXd::Identity(6, 6); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
		//Matrix6d H1;// = Eigen::MatrixXd::Zero(6, 6);
		//Matrix6d H2;// = Eigen::MatrixXd::Zero(6, 6);
		Eps = Matrix6d::Identity() * delta; // faster, all memory on stack
		//H1.setZero();
		//H2.setZero(); // faster, no extra memory; actually not needed at all, the below code overwrites both H1 and H2

		_TyDestMatrix0 &H1 = r_t_pose3_pose1;
		_TyDestMatrix1 &H2 = r_t_pose3_pose2;
		_ASSERTE(H1.rows() == 6 && H1.cols() == 6 && H2.rows() == 6 && H2.cols() == 6); // just make sure the shape is right
		// can actually work inplace

		//Vector6d d;
		Absolute_to_Relative(r_t_vertex1, r_t_vertex2, r_t_dest);
		//r_t_dest = d; // possibly an unnecessary copy

		for(int j = 0; j < 6; ++ j) {
			Vector6d d1, p_delta;
			Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
			Absolute_to_Relative(p_delta, r_t_vertex2, d1);
			H1.col(j) = (d1 - r_t_dest) * scalar;

			Vector6d d2;
			Relative_to_Absolute(r_t_vertex2, Eps.col(j), p_delta);
			Absolute_to_Relative(r_t_vertex1, p_delta, d2);
			H2.col(j) = (d2 - r_t_dest) * scalar;
		}
		// return jacobs

#if 0
		Eigen::Vector4d dist_4D;
		dist_4D.head<3>() = r_t_dest.head<3>(); // dx, dy, dz
		dist_4D(3) = r_t_dest.tail<3>().norm(); // a total angle of rotation
		// the 4D distance

		Eigen::Matrix<double, 4, 6> H1_4, H2_4;
		for(int j = 0; j < 6; ++ j) {
			Vector6d d1, p_delta;
			Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
			Absolute_to_Relative(p_delta, r_t_vertex2, d1);
			Eigen::Vector4d d1_4;
			d1_4.head<3>() = d1.head<3>();
			d1_4(3) = d1.tail<3>().norm();
			H1_4.col(j) = (d1_4 - dist_4D) * scalar;

			Vector6d d2;
			Relative_to_Absolute(r_t_vertex2, Eps.col(j), p_delta);
			Absolute_to_Relative(r_t_vertex1, p_delta, d2);
			Eigen::Vector4d d2_4;
			d2_4.head<3>() = d2.head<3>();
			d2_4(3) = d2.tail<3>().norm();
			H2_4.col(j) = (d2_4 - dist_4D) * scalar;
		}
		// the 4D jacobs

		Eigen::Matrix<double, 4, 6> H1_transform = /*-*/Eigen::Matrix<double, 4, 6>::Identity(),
			H2_transform = Eigen::Matrix<double, 4, 6>::Identity(); // start with identity
		{
			/*const double Aa = r_t_vertex1(3), Ab = r_t_vertex1(4), Ac = r_t_vertex1(5);
			/const double Ba = r_t_vertex2(3), Bb = r_t_vertex2(4), Bc = r_t_vertex2(5);*/
			const double invD = 1 / dist_4D(3);
			/*H1_transform(3, 3) = 1 / 2 * invD * (2*Aa-2*Ba);
			H1_transform(3, 4) = 1 / 2 * invD * (2*Ab-2*Bb);
			H1_transform(3, 5) = 1 / 2 * invD * (2*Ac-2*Bc);
			H2_transform(3, 3) = 1 / 2 * invD * (-2*Aa+2*Ba);
			H2_transform(3, 4) = 1 / 2 * invD * (-2*Ab+2*Bb);
			H2_transform(3, 5) = 1 / 2 * invD * (-2*Ac+2*Bc);*/
			const double Aa = r_t_dest(3), Ab = r_t_dest(4), Ac = r_t_dest(5);
			H1_transform(3, 3) = invD * Aa;
			H1_transform(3, 4) = invD * Ab;
			H1_transform(3, 5) = invD * Ac;
			// matlab code
		}
		// get the projection function to go from 6D jacobians to 4D jacobians // this seems to work extremely well

		double f_err_H1 = (H1_4 - H1_transform * H1).norm();
		double f_err_H2 = (H2_4 - H1_transform/*H2_transform*/ * H2).norm();
		f_err_H1 /= std::max(1.0, H1.norm());
		f_err_H2 /= std::max(1.0, H2.norm());
		fprintf(stderr, "error in jacobians is %g %g\n", f_err_H1, f_err_H2);
		// evaluate precision
#endif // 0
		// debug - projection of the jacobian to a 4D quantity

		/*FILE *p_fw = fopen("jacobs3d.m", "a");
		if(p_fw) {
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, r_t_vertex1, "u = ", ";\n");
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, r_t_vertex2, "v = ", ";\n");
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, r_t_dest, "d_gt = ", "; % difference between u and v\n");
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, H1, "H1_gt = ", "; % the first jacobian\n");
			CDebug::Print_DenseMatrix_in_MatlabFormat(p_fw, H2, "H2_gt = ", "; % the second jacobian\n");
			fprintf(p_fw, "%s", "\n% ---\n\n");
			fclose(p_fw);
		}*/
		// get some "ground truth" for analytical jacobian debugging
#endif // 0
	}

	/**
	 *	@brief converts xyz axis angle coordinates from relative measurement to absolute measurement
	 *		and calculates the jacobians
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *	@tparam Derived3 is Eigen derived matrix type for the fourth matrix argument
	 *	@tparam Derived4 is Eigen derived matrix type for the fifth matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 *	@param[out] r_t_pose3_pose1 is filled with the first jacobian
	 *	@param[out] r_t_pose3_pose2 is filled with the second jacobian
	 */
	template <class Derived0, class Derived1, class Derived2, class Derived3, class Derived4> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d / Matrix3d)
	static void Relative_to_Absolute(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest,
		Eigen::MatrixBase<Derived3> &r_t_pose3_pose1, Eigen::MatrixBase<Derived4> &r_t_pose3_pose2)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector6d>(r_t_vertex2);
		DimensionCheck<Vector6d>(r_t_dest);
		DimensionCheck<Matrix6d>(r_t_pose3_pose1);
		DimensionCheck<Matrix6d>(r_t_pose3_pose2);

		//lets try it according to g2o
		const double delta = 1e-9;
		const double scalar = 1.0 / (delta);

		//jacobians have to be computed here
		//double eps = 1.e-9; //TODO: what is this???
		//double eps_ = sqrt(eps);
		//double eps_2 = 7.45e-9;

		Matrix6d Eps;// = delta * Eigen::MatrixXd::Identity(6, 6); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
		//Matrix6d H1;// = Eigen::MatrixXd::Zero(6, 6);
		//Matrix6d H2;// = Eigen::MatrixXd::Zero(6, 6);
		Eps = Matrix6d::Identity() * delta; // faster, all memory on stack
		//H1.setZero();
		//H2.setZero(); // faster, no extra memory; actually not needed at all, the below code overwrites both H1 and H2

		Eigen::MatrixBase<Derived3> &H1 = r_t_pose3_pose1;
		Eigen::MatrixBase<Derived4> &H2 = r_t_pose3_pose2;
		_ASSERTE(H1.rows() == 6 && H1.cols() == 6 && H2.rows() == 6 && H2.cols() == 6); // just make sure the shape is right
		// can actually work inplace

		//Vector6d d;
		Relative_to_Absolute(r_t_vertex1, r_t_vertex2, r_t_dest);
		//Absolute_to_Relative(r_t_vertex1, r_t_vertex2, r_t_dest);
		//r_t_dest = d; // possibly an unnecessary copy

		for(int j = 0; j < 6; ++ j) {
			Vector6d d1, p_delta;
			Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
			Relative_to_Absolute(p_delta, r_t_vertex2, d1);
			H1.col(j) = (d1 - r_t_dest) * scalar;

			Vector6d d2;
			Relative_to_Absolute(r_t_vertex2, Eps.col(j), p_delta);
			Relative_to_Absolute(r_t_vertex1, p_delta, d2);
			H2.col(j) = (d2 - r_t_dest) * scalar;
		}
		// return jacobs
	}

	/**
	 *	@brief converts xyz axis angle coordinates from absolute measurement to relative measurement
	 *		and calculates the jacobians
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first 6DOF vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second 3DOF vertex, also in absolute coordinates
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 */
	template <class Derived0, class Derived1, class Derived2>
	static void Absolute_to_Relative_Landmark(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector3d>(r_t_dest);

		Eigen::Vector3d p1 = r_t_vertex1.template head<3>();
		Eigen::Quaterniond q1;
		AxisAngle_to_Quat(r_t_vertex1.template tail<3>(), q1);

		Eigen::Vector3d p2 = r_t_vertex2.template head<3>();

		Eigen::Quaterniond q1_inv = q1.conjugate(); // the inverse rotation (also have .inverse() but that is not needed)

		r_t_dest = q1_inv._transformVector(p2 - p1); // this is precise enough
		//r_t_dest.head<3>() = pQ1_inv * (p2 - p1); // or can use matrix to transform position
	}

	/**
	 *	@brief composes a relative transformation and a landmark position to yield
	 *		an aboslute landmark position
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first 6DOF vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second 3DOF vertex, in relative coordinates
	 *	@param[out] r_t_dest is filled with absolute coordinates of the composed vertex
	 */
	template <class Derived0, class Derived1, class Derived2>
	static void Relative_to_Absolute_Landmark(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector3d>(r_t_dest);

		Eigen::Vector3d p1 = r_t_vertex1.template head<3>();
		Eigen::Quaterniond q1;
		AxisAngle_to_Quat(r_t_vertex1.template tail<3>(), q1);
		Eigen::Vector3d p2 = r_t_vertex2.template head<3>();
		r_t_dest = p1 + q1._transformVector(p2);
	}

	/**
	 *	@brief converts xyz axis angle coordinates from absolute measurement to relative measurement
	 *		and calculates the jacobians
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *	@tparam Derived3 is Eigen derived matrix type for the fourth matrix argument
	 *	@tparam Derived4 is Eigen derived matrix type for the fifth matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 *	@param[out] r_t_pose3_pose1 is filled with the first jacobian
	 *	@param[out] r_t_pose3_pose2 is filled with the second jacobian
	 */
	template <class Derived0, class Derived1, class Derived2, class Derived3, class Derived4> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d / Matrix3d)
	static void Absolute_to_Relative_Landmark(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest,
		Eigen::MatrixBase<Derived3> &r_t_pose3_pose1, Eigen::MatrixBase<Derived4> &r_t_pose3_pose2)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector3d>(r_t_dest);
		DimensionCheck<Eigen::Matrix<double, 3, 6> >(r_t_pose3_pose1);
		DimensionCheck<Eigen::Matrix3d>(r_t_pose3_pose2);

		//lets try it according to g2o
		const double delta = 1e-9;
		const double scalar = 1.0 / (delta);

		//jacobians have to be computed here
		//double eps = 1.e-9; //TODO: what is this???
		//double eps_ = sqrt(eps);
		//double eps_2 = 7.45e-9;

		Matrix6d Eps;// = delta * Eigen::MatrixXd::Identity(6, 6); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
		//Matrix6d H1;// = Eigen::MatrixXd::Zero(6, 6);
		//Matrix6d H2;// = Eigen::MatrixXd::Zero(6, 6);
		Eps = Matrix6d::Identity() * delta; // faster, all memory on stack
		//H1.setZero();
		//H2.setZero(); // faster, no extra memory; actually not needed at all, the below code overwrites both H1 and H2

		Eigen::MatrixBase<Derived3> &H1 = r_t_pose3_pose1;
		Eigen::MatrixBase<Derived4> &H2 = r_t_pose3_pose2;
		// can actually work inplace

		Absolute_to_Relative_Landmark(r_t_vertex1, r_t_vertex2, r_t_dest);
		//r_t_dest = d; // possibly an unnecessary copy

		for(int j = 0; j < 6; ++ j) {
			Vector6d p_delta;
			Eigen::Vector3d d1;
			Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
			Absolute_to_Relative_Landmark(p_delta, r_t_vertex2, d1);
			H1.col(j) = (d1 - r_t_dest) * scalar;
		}

		for(int j = 0; j < 3; ++ j) {
			Eigen::Vector3d p_delta, d1;
			p_delta = r_t_vertex2 + Eps.col(j).template head<3>();//col(j);
			Absolute_to_Relative_Landmark(r_t_vertex1, p_delta, d1);
			H2.col(j) = (d1 - r_t_dest) * scalar;
		}
	}
};

#endif // !__3D_SOLVER_BASE_INCLUDED
