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
//#include "slam/BlockMatrix.h"
#include "eigen/Eigen/Cholesky"

/**
 *	@brief base class for a 3D solver
 */
class CBaseBASolver /*: public COptimizationMethod<TVertex3D, TEdge3D>*/ {
public:
	/**
	 *	@brief implementation of Jacobian calculations, required by 3D solvers
	 */
	class CBAJacobians {
	public:
		/**
		 *	@brief converts from axis angle rep to RPY
		 *	@param[in] v_vec is axis-angle rotation (angle is encoded as the magnitude of the axis)
		 *	@return Returns rotation matrix, corresponding to the input.
		 */
		static Eigen::Matrix3d Operator_rot(Eigen::Vector3d v_vec)
		{
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

			double f_half_angle = f_angle * .5;
			v_vec *= sin(f_half_angle); // SSE
			Eigen::Quaternion<double> rot(cos(f_half_angle), v_vec(0), v_vec(1), v_vec(2));
			return rot.toRotationMatrix();
			// use quaternion, Eigen might provide SSE accelerated conversion to matrix
		}

		/**
		 *	@brief converts from RYP rep to axis angle
		 *	@param[in] Q is the rotation matrix
		 *	@return Returns axis-amgle rotation, where the angle is encoded as magnitude of the axis.
		 */
		static Eigen::Vector3d Operator_arot(const Eigen::Matrix3d &Q) // note that Eigen probably implements this somewhere
		{
			Eigen::Quaternion<double> quat(Q); // converts rotation matrix to quaternion (hopefully SSE)
			double f_half_angle = acos(quat.w());
			if(f_half_angle == 0)
				return Eigen::Vector3d(0, 0, 0); // handle zero angle rotation
			return quat.vec() * (2 * f_half_angle / sin(f_half_angle)); // SSE
			// need to divide by sine of half angle to get the normalized rotation axis,
			// then multiply by angle to get axis-angle
		}

		/**
		 *	@brief operator plus for Cam vertices
		 *
		 *	@param[in] r_t_vertex1 is the vertex to be modified
		 *	@param[in] r_t_vertex2 is the delta vector
		 *	@param[out] r_t_dest is the result of the operation
		 */
		static void Smart_Plus_Cam(const Eigen::Matrix<double, 6, 1> &r_t_vertex1,
			const Eigen::Matrix<double, 6, 1> &r_t_vertex2, Eigen::Matrix<double, 6, 1> &r_t_dest)
		{
			_ASSERTE(r_t_vertex1.rows() == 6 && r_t_vertex2.rows() == 6);

			r_t_dest.head<3>() = r_t_vertex1.head<3>() + r_t_vertex2.head<3>(); // accelerated using SSE
			//SOSO: ignore fx, fy, cx, cy, we ahve them fixed
			//r_t_dest.tail<4>() = r_t_vertex1.tail<4>();// + r_t_vertex2.tail<4>(); // accelerated using SSE

			//sum the rotations
			Eigen::Matrix3d pQ = Operator_rot(r_t_vertex1.segment<3>(3));
			Eigen::Matrix3d dQ = Operator_rot(r_t_vertex2.segment<3>(3));

			Eigen::Matrix3d QQ;// = pQ * dQ;
			QQ.noalias() = pQ * dQ; // multiplication without intermediate storage
			r_t_dest.segment<3>(3) = Operator_arot(QQ);
		}

		/**
		 *	@brief operator plus for XYZ vertices
		 *
		 *	@param[in] r_t_vertex1 is the vertex to be modified
		 *	@param[in] r_t_vertex2 is the delta vector
		 *	@param[out] r_t_dest is the result of the operation
		 */
		static void Smart_Plus_XYZ(const Eigen::Matrix<double, 3, 1> &r_t_vertex1,
			const Eigen::Matrix<double, 3, 1> &r_t_vertex2, Eigen::Matrix<double, 3, 1> &r_t_dest)
		{
			_ASSERTE(r_t_vertex1.rows() == 3 && r_t_vertex2.rows() == 3);

			r_t_dest = r_t_vertex1 + r_t_vertex2; // accelerated using SSE
		}

		/**
		 *	@brief converts projects point from 3D to 2D camera image
		 *
		 *	@param[in] r_t_vertex1 is the first vertex - point
		 *	@param[in] r_t_vertex2 is the second vertex - camera
		 *	@param[out] r_t_dest is filled with absolute coordinates of result
		 */
		static void Project_P2C(const Eigen::Matrix<double, 6, 1> &r_t_vertex1, const Eigen::Matrix<double, 5, 1> &intrinsics,
			const Eigen::Matrix<double, 3, 1> &r_t_vertex2, Eigen::Matrix<double, 2, 1> &r_t_dest)
		{
			//construct camera intrinsics matrix
			double fx = intrinsics(0);
			double fy = intrinsics(1);
			//double ap = 1.0029;
			double cx = intrinsics(2);
			double cy = intrinsics(3);
			double k = intrinsics(4);

			Eigen::Matrix<double, 3, 3> A;
			A(0, 0) = fx; 	A(0, 1) = 0; 				A(0, 2) = cx;
			A(1, 0) = 0; 				A(1, 1) = fy; 	A(1, 2) = cy;
			A(2, 0) = 0; 				A(2, 1) = 0; 				A(2, 2) = 1;

			//construct [R | t]
			Eigen::Matrix<double, 3, 4> Rt;
			Rt.block<3, 3>(0, 0) = Operator_rot(r_t_vertex1.segment<3>(3));
			Rt.col(3) = r_t_vertex1.head<3>();

			//construct point vector
			Eigen::Matrix<double, 4, 1> X;
			X.head<3>() = r_t_vertex2;
			/*X(0) = r_t_vertex2(0);
			X(1) = r_t_vertex2(1);
			X(2) = r_t_vertex2(2);*/
			X(3) = 1;

			//project world to camera
			Eigen::Matrix<double, 3, 1> x = Rt * X;

			//project camera to image
			Eigen::Matrix<double, 3, 1> uv = A * x;
			//normalize
			uv = uv / uv(2);

			//apply radial distortion
			double r = sqrt((uv(0)-A(0, 2))*(uv(0)-A(0, 2)) + (uv(1)-A(1, 2))*(uv(1)-A(1, 2)));
			uv(0) =  A(0, 2) + (1 + r*k) * (uv(0) - A(0, 2));
			uv(1) =  A(1, 2) + (1 + r*k) * (uv(1) - A(1, 2));

			r_t_dest(0) = uv(0);
#ifdef DATA_UPSIDE_DOWN
			r_t_dest(1) = -uv(1);
#else
			r_t_dest(1) = uv(1);
#endif
		}

		/**
		 *	@brief converts xyz axis angle coordinates from absolute measurement to relative measurement
		 *		and calculates the jacobians
		 *
		 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates - CAM
		 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates - XYZ
		 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
		 *	@param[out] r_t_pose3_pose1 is filled with the first jacobian
		 *	@param[out] r_t_pose3_pose2 is filled with the second jacobian
		 */
		template <class _TyDestVector, class _TyDestMatrix0, class _TyDestMatrix1> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d / Matrix3d)
		static void Project_P2C(const Eigen::Matrix<double, 6, 1> &r_t_vertex1, const Eigen::Matrix<double, 5, 1> &intrinsics,
			const Eigen::Matrix<double, 3, 1> &r_t_vertex2, _TyDestVector &r_t_dest,
			_TyDestMatrix0 &r_t_pose3_pose1, _TyDestMatrix1 &r_t_pose3_pose2)
		{
			/* TODO: make more efficient */
			//lets try it according to g2o
			const double delta = 1e-9;
			const double scalar = 1.0 / (delta);
			//const double delta_pixel = 1e-6;
			//const double scalar_pixel = 1.0 / (delta_pixel);

			Eigen::Matrix<double, 6, 6> Eps;// = delta * Eigen::MatrixXd::Identity(10, 10); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
			Eps = Eigen::Matrix<double, 6, 6>::Identity() * delta; // faster, all memory on stack
			//Eps(6, 6) = delta_pixel; Eps(7, 7) = delta_pixel;
			//Eps(9, 9) = delta_pixel; Eps(8, 8) = delta_pixel;

			_TyDestMatrix0 &H1 = r_t_pose3_pose1;
			_TyDestMatrix1 &H2 = r_t_pose3_pose2;
			_ASSERTE(H1.rows() == 2 && H1.cols() == 6 && H2.rows() == 2 && H2.cols() == 3); // just make sure the shape is right
			// can actually work inplace

			Project_P2C(r_t_vertex1, intrinsics, r_t_vertex2, r_t_dest);
			//r_t_dest = d; // possibly an unnecessary copy

			Eigen::Matrix<double, 2, 1> d1;
			Eigen::Matrix<double, 6, 1> p_delta;
			Eigen::Matrix<double, 3, 1> p_delta2;
			//for XYZ and RPY
			for(int j = 0; j < 6; ++ j) {

				Smart_Plus_Cam(r_t_vertex1, Eps.block(0, j, 6, 1), p_delta);
				Project_P2C(p_delta, intrinsics, r_t_vertex2, d1);
				//if(j < 6)
					H1.block(0, j, 2, 1) = (d1 - r_t_dest) * scalar;
				//else {
					//put just zeros in here
				//	H1.block(0, j, 2, 1) = Eigen::Matrix<double, 2, 1>::Zero(2, 1)/*(d1 - r_t_dest) * scalar_pixel*/;
				//}
				//if(j < 4)
				//	H1.block(0, 6+j, 2, 1) = Eigen::Matrix<double, 2, 1>::Zero(2, 1);

				if(j < 3) {
					Smart_Plus_XYZ(r_t_vertex2, Eps.block(0, j, 3, 1), p_delta2);
					Project_P2C(r_t_vertex1, intrinsics, p_delta2, d1);
					H2.block(0, j, 2, 1) = (d1 - r_t_dest) * scalar;
				}
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
			/*float fx = intrinsics(0); float fy = intrinsics(1);
			float tx = r_t_vertex1(0); float ty = r_t_vertex1(1); float tz = r_t_vertex1(2);
			float x = r_t_vertex2(0); float y = r_t_vertex2(1); float z = r_t_vertex2(2);
			Eigen::Matrix3d R = Operator_rot(r_t_vertex1.segment<3>(3));
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
	};
};

#endif // __BA_SOLVER_BASE_INCLUDED
