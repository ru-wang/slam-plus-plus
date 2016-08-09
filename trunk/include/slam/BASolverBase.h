/*
								+----------------------------------+
								|                                  |
								| *** Base class for 3D solver *** |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2013  |
								|                                  |
								|          BASolverBase.h          |
								|                                  |
								+----------------------------------+
*/

#pragma once
#ifndef __BA_SOLVER_BASE_INCLUDED
#define __BA_SOLVER_BASE_INCLUDED

//#define DATA_UPSIDE_DOWN

/**
 *	@file include/slam/BASolverBase.h
 *	@brief a simple base class for BA solver, made according to 2D solver
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
#include <iostream>
#include "slam/BlockMatrix.h"
#include "eigen/Eigen/Cholesky"
#include "slam/3DSolverBase.h" // conversions between axis angles and rotation matrices

/**
 *	@brief implementation of Jacobian calculations, required by 3D solvers
 */
class CBAJacobians {
public:
	typedef Eigen::Matrix<double, 5, 1> Vector5d; /**< @brief 5D vector type */
	typedef Eigen::Matrix<double, 6, 1> Vector6d; /**< @brief 6D vector type */
	typedef Eigen::Matrix<double, 6, 6> Matrix6d; /**< @brief 6x6 matrix type */

public:

	/**
	 *	@brief converts from axis angle rep to a rotation matrix
	 *	@param[in] v_vec is axis-angle rotation (angle is encoded as the magnitude of the axis)
	 *	@return Returns rotation matrix, corresponding to the input.
	 *	@deprecated This function is deprecated in favor of C3DJacobians::t_AxisAngle_to_RotMatrix(). Please, do not use it.
	 */
	static Eigen::Matrix3d Operator_rot(Eigen::Vector3d v_vec)
	{
		return C3DJacobians::t_AxisAngle_to_RotMatrix(v_vec);
		/*double f_angle = v_vec.norm();//sqrt(x*x + y*y + z*z); // SSE
		// this is really angle in radians

		if(f_angle > 0) {
			v_vec /= f_angle; // SSE
			//x = x / f_angle;
			//y = y / f_angle;
			//z = z / f_angle;
		} else
			return Eigen::Matrix3d::Identity(); // no rotation, save some hassle below
		// normalize the axis

		double f_half_angle = f_angle * .5;
		v_vec *= sin(f_half_angle); // SSE
		Eigen::Quaternion<double> rot(cos(f_half_angle), v_vec(0), v_vec(1), v_vec(2));
		return rot.toRotationMatrix();
		// use quaternion, Eigen might provide SSE accelerated conversion to matrix*/
	}

	/**
	 *	@brief converts from a rotation matrix rep to axis angle
	 *	@param[in] Q is the rotation matrix
	 *	@return Returns axis-amgle rotation, where the angle is encoded as magnitude of the axis.
	 *	@deprecated This function is deprecated in favor of C3DJacobians::v_RotMatrix_to_AxisAngle(). Please, do not use it.
	 */
	static Eigen::Vector3d Operator_arot(const Eigen::Matrix3d &Q) // note that Eigen probably implements this somewhere
	{
		return C3DJacobians::v_RotMatrix_to_AxisAngle(Q);
		/*Eigen::Quaternion<double> quat(Q); // converts rotation matrix to quaternion (hopefully SSE)
		double f_half_angle = acos(quat.w());
		if(f_half_angle == 0)
			return Eigen::Vector3d(0, 0, 0); // handle zero angle rotation
		return quat.vec() * (2 * f_half_angle / sin(f_half_angle)); // SSE
		// need to divide by sine of half angle to get the normalized rotation axis,
		// then multiply by angle to get axis-angle*/
	}

	/**
	 *	@brief operator plus for Cam vertices
	 *
	 *	@param[in] r_t_vertex1 is the vertex to be modified
	 *	@param[in] r_t_vertex2 is the delta vector
	 *	@param[out] r_t_dest is the result of the operation
	 *
	 *	@deprecated This function is deprecated in favor of C3DJacobians::Relative_to_Absolute(). Please, do not use it.
	 */
	template <class _TyDestVector>
	static void Smart_Plus_Cam(const Vector6d &r_t_vertex1,
		const Vector6d &r_t_vertex2, _TyDestVector &r_t_dest)
	{
		C3DJacobians::Relative_to_Absolute(r_t_vertex1, r_t_vertex2, r_t_dest);

		/*r_t_dest.template head<3>() = r_t_vertex1.head<3>() + r_t_vertex2.head<3>(); // accelerated using SSE
		//SOSO: ignore fx, fy, cx, cy, we have them fixed
		//r_t_dest.tail<4>() = r_t_vertex1.tail<4>();// + r_t_vertex2.tail<4>(); // accelerated using SSE

		//sum the rotations
		Eigen::Matrix3d pQ = C3DJacobians::t_AxisAngle_to_RotMatrix(r_t_vertex1.tail<3>());
		Eigen::Matrix3d dQ = C3DJacobians::t_AxisAngle_to_RotMatrix(r_t_vertex2.tail<3>());

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // multiplication without intermediate storage
		_ASSERTE(r_t_dest.rows() == 6);
		r_t_dest.template tail<3>() = C3DJacobians::v_RotMatrix_to_AxisAngle(QQ);*/
	}

	/**
	 *	@brief operator plus for XYZ vertices
	 *
	 *	@param[in] r_t_vertex1 is the vertex to be modified
	 *	@param[in] r_t_vertex2 is the delta vector
	 *	@param[out] r_t_dest is the result of the operation
	 *
	 *	@deprecated This function is deprecated in favor of CBAJacobians::Relative_to_Absolute_XYZ(). Please, do not use it.
	 */
	static void Smart_Plus_XYZ(const Eigen::Vector3d &r_t_vertex1,
		const Eigen::Vector3d &r_t_vertex2, Eigen::Vector3d &r_t_dest)
	{
		Relative_to_Absolute_XYZ(r_t_vertex1, r_t_vertex2, r_t_dest);
	}

	/**
	 *	@brief composition for XYZ vertices
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the vertex to be modified
	 *	@param[in] r_t_vertex2 is the delta vector
	 *	@param[out] r_t_dest is the result of the composition
	 */
	template <class Derived0, class Derived1, class Derived2>
	static void Relative_to_Absolute_XYZ(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest)
	{
		DimensionCheck<Eigen::Vector3d>(r_t_vertex1);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector3d>(r_t_dest);

		r_t_dest = r_t_vertex1 + r_t_vertex2; // accelerated using SSE
	}

	/**
	 *	@brief operator plus for intrinsics vertex
	 *
	 *	@param[in] r_t_vertex1 is the vertex to be modified
	 *	@param[in] r_t_vertex2 is the delta vector
	 *	@param[out] r_t_dest is the result of the operation
	 *
	 *	@deprecated This function is deprecated in favor of CBAJacobians::Relative_to_Absolute_Intrinsics(). Please, do not use it.
	 */
	static void Smart_Plus_Intrinsics(const Vector5d &r_t_vertex1,
		const Vector5d &r_t_vertex2, Vector5d &r_t_dest)
	{
		Relative_to_Absolute_Intrinsics(r_t_vertex1, r_t_vertex2, r_t_dest);
	}

	/**
	 *	@brief composition for intrinsic parameters vertices
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the vertex to be modified
	 *	@param[in] r_t_vertex2 is the delta vector
	 *	@param[out] r_t_dest is the result of the composition
	 */
	template <class Derived0, class Derived1, class Derived2>
	static void Relative_to_Absolute_Intrinsics(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_vertex2, Eigen::MatrixBase<Derived2> &r_t_dest)
	{
		DimensionCheck<Vector5d>(r_t_vertex1);
		DimensionCheck<Vector5d>(r_t_vertex2);
		DimensionCheck<Vector5d>(r_t_dest);

		r_t_dest = r_t_vertex1 + r_t_vertex2; // accelerated using SSE
	}

	/**
	 *	@brief undistorts 2d point
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *
	 *	@param[in] pt is the point to be undistorted
	 *	@param[in] r_t_intrinsics is instrinsic parameters of the camera
	 *	@param[out] r_t_dest is the result of the operation
	 */
	template <class Derived0, class Derived1, class Derived2>
	static void Undistort_Point2D(const Eigen::MatrixBase<Derived0> &pt,
		const Eigen::MatrixBase<Derived1> &r_t_intrinsics, Eigen::MatrixBase<Derived2> &r_t_dest)
	{
		DimensionCheck<Eigen::Vector2d>(pt);
		DimensionCheck<Vector5d>(r_t_intrinsics);
		DimensionCheck<Eigen::Vector2d>(r_t_dest);

		double fx = r_t_intrinsics(0);
		double fy = r_t_intrinsics(1);
		double cx = r_t_intrinsics(2);
		double cy = r_t_intrinsics(3);
		double k = r_t_intrinsics(4) / (.5 * (fx + fy));

		Eigen::Vector2d c(cx, cy);
		double r2 = (pt - c).squaredNorm();

		//result
		r_t_dest = c + (1 + r2 * k) * (pt - c);
	}

	/**
	 *	@brief converts projects point from 3D to 2D camera image
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *	@tparam Derived3 is Eigen derived matrix type for the fourth matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex - camera
	 *	@param[in] r_t_intrinsics is instrinsic parameters of the camera
	 *	@param[in] r_t_vertex2 is the second vertex - point
	 *	@param[out] r_t_dest is filled with absolute coordinates of result
	 */
	template <class Derived0, class Derived1, class Derived2, class Derived3>
	static void Project_P2C(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_intrinsics,
		const Eigen::MatrixBase<Derived2> &r_t_vertex2,
		Eigen::MatrixBase<Derived3> &r_t_dest)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector5d>(r_t_intrinsics);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector2d>(r_t_dest);

		//construct camera r_t_intrinsics matrix
		double fx = r_t_intrinsics(0);
		double fy = r_t_intrinsics(1);
		//double ap = 1.0029;
		double cx = r_t_intrinsics(2);
		double cy = r_t_intrinsics(3);
		double k = r_t_intrinsics(4) / (.5 * (fx + fy))/*sqrt(fx * fy)*/; // rescale intrinsics// since post-ICRA2015 release the distortion is scaled by focal length in the internal representation

		Eigen::Matrix3d A;
		/*A(0, 0) = fx; 	A(0, 1) = 0; 	A(0, 2) = cx;
		A(1, 0) = 0; 	A(1, 1) = fy; 	A(1, 2) = cy;
		A(2, 0) = 0; 	A(2, 1) = 0; 	A(2, 2) = 1;*/
		A << fx,  0, cx,
			  0, fy, cy,
			  0,  0,  1;
		Eigen::Vector2d c(cx, cy); // the ?optical? center? is there a better letter describing it?

		//construct [R | t]
		Eigen::Matrix<double, 3, 4> Rt;
		Rt.template block<3, 3>(0, 0) = C3DJacobians::t_AxisAngle_to_RotMatrix(r_t_vertex1.template tail<3>());
		Rt.col(3) = r_t_vertex1.template head<3>();

		//construct point vector
		Eigen::Vector4d X;
		X.template head<3>() = r_t_vertex2;
		/*X(0) = r_t_vertex2(0);
		X(1) = r_t_vertex2(1);
		X(2) = r_t_vertex2(2);*/
		X(3) = 1;

		//project world to camera
		Eigen::Vector3d x = Rt * X;

		//project camera to image
		Eigen::Vector3d uv = A * x;
		//normalize
		uv /= uv(2);

		//apply radial distortion
		double r = (uv.head<2>() - c).norm();//sqrt((uv(0)-A(0, 2))*(uv(0)-A(0, 2)) + (uv(1)-A(1, 2))*(uv(1)-A(1, 2)));
		/*uv(0) =  A(0, 2) + (1 + r*k) * (uv(0) - A(0, 2));
		uv(1) =  A(1, 2) + (1 + r*k) * (uv(1) - A(1, 2));*/
		uv.template head<2>() = c + (1 + r*r * k) * (uv.template head<2>() - c);

		//lets try 04J09 paper distortion
		/*int sign = (uv(0) > 0) ? 1 : ((uv(0) < 0) ? -1 : 0);
		double K = uv(1) / uv(0);	//K = y/x
		std::cout << "r: " << sign * uv(0) * uv(0) << " f: " <<  k * sqrt(1 + K*K) << std::endl;
		uv(0) = uv(0) - k * sqrt(1 + K*K) * sign * uv(0) * uv(0);	//we dont use k2 so rest is 0
		uv(1) = K * uv(0);*/

		r_t_dest(0) = uv(0);
#ifdef DATA_UPSIDE_DOWN
		r_t_dest(1) = -uv(1);
#else
		r_t_dest(1) = uv(1);
#endif
	}

	/**
	 *	@brief converts projects point from 3D to 2D camera image
	 *
	 *	@tparam Derived0 is Eigen derived matrix type for the first matrix argument
	 *	@tparam Derived1 is Eigen derived matrix type for the second matrix argument
	 *	@tparam Derived2 is Eigen derived matrix type for the third matrix argument
	 *	@tparam Derived3 is Eigen derived matrix type for the fourth matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex - camera
	 *	@param[in] r_t_intrinsics is instrinsic parameters of the camera
	 *	@param[in] r_t_vertex2 is the second vertex - point
	 *	@param[out] r_t_dest is filled with absolute coordinates of result
	 */
	template <class Derived0, class Derived1, class Derived2, class Derived3>
	static void Project_P2SC(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_intrinsics,
		const Eigen::MatrixBase<Derived2> &r_t_vertex2,
		Eigen::MatrixBase<Derived3> &r_t_dest)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector6d>(r_t_intrinsics);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector3d>(r_t_dest);

		//construct camera r_t_intrinsics matrix
		double fx = r_t_intrinsics(0);
		double fy = r_t_intrinsics(1);
		//double ap = 1.0029;
		double cx = r_t_intrinsics(2);
		double cy = r_t_intrinsics(3);
		double k = r_t_intrinsics(4) / (.5 * (fx + fy))/*sqrt(fx * fy)*/; // rescale intrinsics// since post-ICRA2015 release the distortion is scaled by focal length in the internal representation
		double b = r_t_intrinsics(5);

		Eigen::Matrix3d A;
		/*A(0, 0) = fx; 	A(0, 1) = 0; 				A(0, 2) = cx;
		A(1, 0) = 0; 				A(1, 1) = fy; 	A(1, 2) = cy;
		A(2, 0) = 0; 				A(2, 1) = 0; 				A(2, 2) = 1;*/
		A << fx,  0, cx,
			  0, fy, cy,
			  0,  0,  1;
		Eigen::Vector2d c(cx, cy); // the ?optical? center? is there a better letter describing it?

		//construct [R | t]
		Eigen::Matrix<double, 3, 4> Rt;
		Rt.template block<3, 3>(0, 0) = C3DJacobians::t_AxisAngle_to_RotMatrix(r_t_vertex1.template tail<3>());
		Rt.col(3) = r_t_vertex1.template head<3>();

		//construct point vector
		Eigen::Vector4d X;
		X.template head<3>() = r_t_vertex2;
		/*X(0) = r_t_vertex2(0);
		X(1) = r_t_vertex2(1);
		X(2) = r_t_vertex2(2);*/
		X(3) = 1;

		//project world to left camera
		Eigen::Vector3d x = Rt * X;

		//project camera to image
		Eigen::Vector3d uv = A * x;
		//normalize
		uv = uv / uv(2);

		//apply radial distortion
		/*double r = sqrt((uv(0)-A(0, 2))*(uv(0)-A(0, 2)) + (uv(1)-A(1, 2))*(uv(1)-A(1, 2)));
		uv(0) =  A(0, 2) + (1 + r*k) * (uv(0) - A(0, 2));
		uv(1) =  A(1, 2) + (1 + r*k) * (uv(1) - A(1, 2));*/
		{
			double r = (uv.head<2>() - c).norm();
			uv.template head<2>() = c + (1 + r * k) * (uv.template head<2>() - c);
		}

		//project point to right cam - instead of moving cam to left, we move point to right
		//X(0) = X(0) - b;
		X.head<3>() -= b * Rt.row(0).template head<3>(); // the first row of the rotation matrix is the direction of the x axis
		x = Rt * X;
		Eigen::Vector3d uv2 = A * x;
		uv2 = uv2 / uv2(2);
		/*r = sqrt((uv2(0)-A(0, 2))*(uv2(0)-A(0, 2)) + (uv2(1)-A(1, 2))*(uv2(1)-A(1, 2)));
		uv2(0) =  A(0, 2) + (1 + r*k) * (uv2(0) - A(0, 2));
		uv2(1) =  A(1, 2) + (1 + r*k) * (uv2(1) - A(1, 2));*/
		{
			double r = (uv2.template head<2>() - c).norm();
			uv2.template head<2>() = c + (1 + r * k) * (uv2.template head<2>() - c);
		}

		r_t_dest(0) = uv(0);
		r_t_dest(1) = uv(1);
		r_t_dest(2) = uv2(0);
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
	 *	@tparam Derived5 is Eigen derived matrix type for the sixth matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates - CAM
	 *	@param[in] r_t_intrinsics is instrinsic parameters of the camera
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates - XYZ
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 *	@param[out] r_t_pose3_pose1 is filled with the first jacobian
	 *	@param[out] r_t_pose3_pose2 is filled with the second jacobian
	 */
	template <class Derived0, class Derived1, class Derived2, class Derived3, class Derived4, class Derived5> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d / Matrix3d)
	static void Project_P2C(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_intrinsics,
		const Eigen::MatrixBase<Derived2> &r_t_vertex2, Eigen::MatrixBase<Derived3> &r_t_dest,
		Eigen::MatrixBase<Derived4> &r_t_pose3_pose1, Eigen::MatrixBase<Derived5> &r_t_pose3_pose2)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector5d>(r_t_intrinsics);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector2d>(r_t_dest);
		DimensionCheck<Eigen::Matrix<double, 2, 6> >(r_t_pose3_pose1);
		DimensionCheck<Eigen::Matrix<double, 2, 3> >(r_t_pose3_pose2);

		/* TODO: make more efficient */
		//lets try it according to g2o
		const double delta = 1e-9;
		const double scalar = 1.0 / (delta);
		//const double delta_pixel = 1e-6;
		//const double scalar_pixel = 1.0 / (delta_pixel);

		Matrix6d Eps;// = delta * Eigen::MatrixXd::Identity(10, 10); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
		Eps = Matrix6d::Identity() * delta; // faster, all memory on stack
		//Eps(6, 6) = delta_pixel; Eps(7, 7) = delta_pixel;
		//Eps(9, 9) = delta_pixel; Eps(8, 8) = delta_pixel;

		Eigen::MatrixBase<Derived4> &H1 = r_t_pose3_pose1;
		Eigen::MatrixBase<Derived5> &H2 = r_t_pose3_pose2;
		// can actually work inplace

		Project_P2C(r_t_vertex1, r_t_intrinsics, r_t_vertex2, r_t_dest);
		//r_t_dest = d; // possibly an unnecessary copy

		Eigen::Vector2d d1;
		Vector6d p_delta;
		Eigen::Vector3d p_delta2;

		//for XYZ and RPY
		for(int j = 0; j < 6; ++ j) {
			C3DJacobians::Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
			Project_P2C(p_delta, r_t_intrinsics, r_t_vertex2, d1);
			//if(j < 6)
				H1.col(j) = (d1 - r_t_dest) * scalar;
			//else {
				//put just zeros in here
			//	H1.block(0, j, 2, 1) = Eigen::Matrix<double, 2, 1>::Zero(2, 1)/*(d1 - r_t_dest) * scalar_pixel*/;
			//}
			//if(j < 4)
			//	H1.block(0, 6+j, 2, 1) = Eigen::Matrix<double, 2, 1>::Zero(2, 1);
		}
		for(int j = 0; j < 3; ++ j) {
			//if(j < 3) {
				Relative_to_Absolute_XYZ(r_t_vertex2, Eps.col(j).template head<3>(), p_delta2);
				Project_P2C(r_t_vertex1, r_t_intrinsics, p_delta2, d1);
				H2.col(j) = (d1 - r_t_dest) * scalar;
			//}
		}

		/*fprintf(stderr, "H2 - num\n");
		for(int a = 0; a < 2; a++) {
			for(int b = 0; b < 3; b++) {
				fprintf(stderr, "%f ", H2(a,b));
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");*/

		//lets do analytical H2
		/*float fx = r_t_intrinsics(0); float fy = r_t_intrinsics(1);
		float tx = r_t_vertex1(0); float ty = r_t_vertex1(1); float tz = r_t_vertex1(2);
		float x = r_t_vertex2(0); float y = r_t_vertex2(1); float z = r_t_vertex2(2);
		Eigen::Matrix3d R = C3DJacobians::t_AxisAngle_to_RotMatrix(r_t_vertex1.tail<3>());
		float r11 = R(0, 0); float r12 = R(0, 1); float r13 = R(0, 2);
		float r21 = R(1, 0); float r22 = R(1, 1); float r23 = R(1, 2);
		float r31 = R(2, 0); float r32 = R(2, 1); float r33 = R(2, 2);

		H2(0, 0) = (fx*(r11*tz + r11*r32*y + r11*r33*z) - fx*r31*(tx + r12*y + r13*z))/pow((tz + r31*x + r32*y + r33*z),2);
		H2(0, 1) = (fx*(r12*tz + r12*r31*x + r12*r33*z) - fx*r32*(tx + r11*x + r13*z))/pow((tz + r31*x + r32*y + r33*z),2);
		H2(0, 2) = -(fx*(r33*tx - r13*tz) + fx*x*(r11*r33 - r13*r31) + fx*y*(r12*r33 - r13*r32))/pow((tz + r31*x + r32*y + r33*z),2);

		H2(1, 0) = (fy*(r21*tz + r21*r32*y + r21*r33*z) - fy*r31*(ty + r22*y + r23*z))/pow((tz + r31*x + r32*y + r33*z),2);
		H2(1, 1) = (fy*(r22*tz + r22*r31*x + r22*r33*z) - fy*r32*(ty + r21*x + r23*z))/pow((tz + r31*x + r32*y + r33*z),2);
		H2(1, 2) = -(fy*(r33*ty - r23*tz) + fy*x*(r21*r33 - r23*r31) + fy*y*(r22*r33 - r23*r32))/pow((tz + r31*x + r32*y + r33*z),2);*/

		/*fprintf(stderr, "-------------------\n\n");
		fprintf(stderr, "H1\n");
		for(int a = 0; a < 2; a++) {
			for(int b = 0; b < 10; b++) {
				fprintf(stderr, "%f ", H1(a,b));
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
		fprintf(stderr, "H2 - anal\n");
		for(int a = 0; a < 2; a++) {
			for(int b = 0; b < 3; b++) {
				fprintf(stderr, "%f ", H2(a,b));
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");*/
		// return jacobs
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
	 *	@tparam Derived5 is Eigen derived matrix type for the sixth matrix argument
	 *	@tparam Derived6 is Eigen derived matrix type for the seventh matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates - CAM
	 *	@param[in] r_t_intrinsics is instrinsic parameters of the camera
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates - XYZ
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 *	@param[out] r_t_cam_measurement is filled with the first jacobian
	 *	@param[out] r_t_xyz_measurement is filled with the second jacobian
	 *	@param[out] r_t_intrinsics_measurement is filled with the third jacobian
	 */
	template <class Derived0, class Derived1, class Derived2,
		class Derived3, class Derived4, class Derived5, class Derived6>
	static void Project_P2CI(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_intrinsics,
		const Eigen::MatrixBase<Derived2> &r_t_vertex2, Eigen::MatrixBase<Derived3> &r_t_dest,
		Eigen::MatrixBase<Derived4> &r_t_cam_measurement, Eigen::MatrixBase<Derived5> &r_t_xyz_measurement,
		Eigen::MatrixBase<Derived6> &r_t_intrinsics_measurement)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector5d>(r_t_intrinsics);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector2d>(r_t_dest);
		DimensionCheck<Eigen::Matrix<double, 2, 6> >(r_t_cam_measurement);
		DimensionCheck<Eigen::Matrix<double, 2, 3> >(r_t_xyz_measurement);
		DimensionCheck<Eigen::Matrix<double, 2, 5> >(r_t_intrinsics_measurement);

		/* TODO: make more efficient */
		//lets try it according to g2o
		const double delta = 1e-9;
		const double scalar = 1.0 / (delta);
		/*const double delta_pixel = 1e-6;
		const double scalar_pixel = 1.0 / (delta_pixel);
		const double delta_distortion = 1e-12;
		const double scalar_distortion = 1.0 / (delta_pixel);*/

		Matrix6d Eps;// = delta * Eigen::MatrixXd::Identity(10, 10); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
		Eps = Matrix6d::Identity() * delta; // faster, all memory on stack

		//compute camera jacobian
		Eigen::MatrixBase<Derived4> &J_cam = r_t_cam_measurement;
		Eigen::MatrixBase<Derived5> &J_xyz = r_t_xyz_measurement;
		Eigen::MatrixBase<Derived6> &J_intrinsics = r_t_intrinsics_measurement;

		Project_P2C(r_t_vertex1, r_t_intrinsics, r_t_vertex2, r_t_dest);

		Eigen::Vector2d d1;
		Vector6d p_delta_cam;
		Eigen::Vector3d p_delta_xyz;
		Vector5d p_delta_intrinsics;

		//for camera
		for(int j = 0; j < 6; ++ j) {
			C3DJacobians::Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta_cam);
			Project_P2C(p_delta_cam, r_t_intrinsics, r_t_vertex2, d1);
			J_cam.col(j) = (d1 - r_t_dest) * scalar;
		}
		//for xyz
		for(int j = 0; j < 3; ++ j) {
			Relative_to_Absolute_XYZ(r_t_vertex2, Eps.col(j).template head<3>(), p_delta_xyz);
			Project_P2C(r_t_vertex1, r_t_intrinsics, p_delta_xyz, d1);
			J_xyz.col(j) = (d1 - r_t_dest) * scalar;
		}
		//for intrinsics
		//Eps(4, 4) = 1e-12; // can't!
		for(int j = 0; j < 5; ++ j) {
			Relative_to_Absolute_Intrinsics(r_t_intrinsics, Eps.col(j).template head<5>(), p_delta_intrinsics);
			Project_P2C(r_t_vertex1, p_delta_intrinsics, r_t_vertex2, d1);
			J_intrinsics.col(j) = (d1 - r_t_dest) * scalar;
		}

		//send jacobians out
		/*r_t_cam_measurement = J_cam;
		r_t_xyz_measurement = J_xyz;
		r_t_intrinsics_measurement = J_intrinsics;*/
		// avoid unnecessary copy: they are all references
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
	 *	@tparam Derived5 is Eigen derived matrix type for the sixth matrix argument
	 *
	 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates - CAM
	 *	@param[in] r_t_intrinsics is instrinsic parameters of the camera
	 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates - XYZ
	 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
	 *	@param[out] r_t_pose3_pose1 is filled with the first jacobian
	 *	@param[out] r_t_pose3_pose2 is filled with the second jacobian
	 */
	template <class Derived0, class Derived1, class Derived2, class Derived3, class Derived4, class Derived5> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d / Matrix3d)
	static void Project_P2SC(const Eigen::MatrixBase<Derived0> &r_t_vertex1,
		const Eigen::MatrixBase<Derived1> &r_t_intrinsics,
		const Eigen::MatrixBase<Derived2> &r_t_vertex2, Eigen::MatrixBase<Derived3> &r_t_dest,
		Eigen::MatrixBase<Derived4> &r_t_pose3_pose1, Eigen::MatrixBase<Derived5> &r_t_pose3_pose2)
	{
		DimensionCheck<Vector6d>(r_t_vertex1);
		DimensionCheck<Vector6d>(r_t_intrinsics);
		DimensionCheck<Eigen::Vector3d>(r_t_vertex2);
		DimensionCheck<Eigen::Vector3d>(r_t_dest);
		DimensionCheck<Eigen::Matrix<double, 3, 6> >(r_t_pose3_pose1);
		DimensionCheck<Eigen::Matrix3d>(r_t_pose3_pose2);

		/* TODO: make more efficient */
		//lets try it according to g2o
		const double delta = 1e-9;
		const double scalar = 1.0 / (delta);
		//const double delta_pixel = 1e-6;
		//const double scalar_pixel = 1.0 / (delta_pixel);

		Matrix6d Eps;// = delta * Eigen::MatrixXd::Identity(10, 10); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
		Eps = Matrix6d::Identity() * delta; // faster, all memory on stack
		//Eps(6, 6) = delta_pixel; Eps(7, 7) = delta_pixel;
		//Eps(9, 9) = delta_pixel; Eps(8, 8) = delta_pixel;

		Eigen::MatrixBase<Derived4> &H1 = r_t_pose3_pose1;
		Eigen::MatrixBase<Derived5> &H2 = r_t_pose3_pose2;
		// can actually work inplace

		Project_P2SC(r_t_vertex1, r_t_intrinsics, r_t_vertex2, r_t_dest);
		//r_t_dest = d; // possibly an unnecessary copy

		Eigen::Vector3d d1;
		Vector6d p_delta;
		Eigen::Vector3d p_delta2;
		//for XYZ and RPY
		for(int j = 0; j < 6; ++ j) {

			C3DJacobians::Relative_to_Absolute(r_t_vertex1, Eps.col(j), p_delta);
			Project_P2SC(p_delta, r_t_intrinsics, r_t_vertex2, d1);
			//if(j < 6)
				H1.col(j) = (d1 - r_t_dest) * scalar;
			//else {
				//put just zeros in here
			//	H1.block(0, j, 2, 1) = Eigen::Matrix<double, 2, 1>::Zero(2, 1)/*(d1 - r_t_dest) * scalar_pixel*/;
			//}
			//if(j < 4)
			//	H1.block(0, 6+j, 2, 1) = Eigen::Matrix<double, 2, 1>::Zero(2, 1);

			if(j < 3) {
				Relative_to_Absolute_XYZ(r_t_vertex2, Eps.col(j).template head<3>(), p_delta2);
				Project_P2SC(r_t_vertex1, r_t_intrinsics, p_delta2, d1);
				H2.col(j) = (d1 - r_t_dest) * scalar;
			}
		}
	}
};

#endif // !__BA_SOLVER_BASE_INCLUDED
