/*
								+----------------------------------+
								|                                  |
								| *** Base class for 3D solver *** |
								|                                  |
								|   Copyright © -tHE SWINe- 2013   |
								|                                  |
								|          3DSolverBase.h          |
								|                                  |
								+----------------------------------+
*/

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
//#include "slam/BlockMatrix.h"
#include "eigen/Eigen/Cholesky"
#ifdef __3D_SOLVER_BASE_COMPILE_LEGACY_CODE
#include "slam/System.h"
#endif // __3D_SOLVER_BASE_COMPILE_LEGACY_CODE

#ifdef __3D_SOLVER_BASE_COMPILE_LEGACY_CODE

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
 *	@brief a 2D vertex (a vertex for 2D problems)
 */
struct TVertex3D { friend class CBase3DSolver; // CBase2DSolver can access _v_state directly, it knows what it is doing (derived classes such as IncrementalCholesky.* can't)
protected:
	Eigen::VectorXd _v_state; /**< @brief the state vector (actually stores 2D or 3D vector, that's why it is not supposed to be accessed directly) */ // t_odo - make this an array of doubles. that should make things much easier // no need, superseded by the flat system

public:
	int n_vertex_dimension; /**< @brief there are *two* kinds of vertices, vertices (3D) and landmarks (2D); both are stored using Vector3d to avoid dynamic memory allocation, but this is the actual dimension */

	size_t n_x_vector_first_elem_index; /**< @brief permutation function */
	//int n_x_vector_block_index; // t_odo - will i ever need this? (with the block matrix probably, yes ...) // no

	/**
	 *	@brief state vector accessor
	 *	@param[in] n_index is zero-based index of the element to be accessed
	 *	@return Returns the value of the selected element of state vector.
	 */
	inline double operator [](int n_index) const
	{
		_ASSERTE(n_index >= 0 && n_index < n_vertex_dimension);
		return _v_state(n_index);
	}

	void (*p_ordering_function)(TVertex3D &r_t_vertex, size_t &r_n_x_index); /**< @brief calculates indices into X vector */
	void (*p_operator_plus)(TVertex3D &r_t_vertex, const double *p_global_deltas); /**< @brief should add v_state += p_global_deltas[X_vector_indices] */
	void (*p_operator_minus)(TVertex3D &r_t_vertex, const double *p_global_deltas); /**< @brief should subtract v_state -= p_global_deltas[X_vector_indices] */

	// t_odo - convert these "virtual" functions to members; there's always (AFAIK) just a single kind of vertex in the system so it may as well directly implement those
	// there are *two* kinds of vertices, vertices (3D) and landmarks (2D)
	// t_odo - maybe it would be possible to use virtual functions if the derived classes had the same size // it is

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 *	@brief a 3D edge (an edge for 3D problems)
 */
struct TEdge3D {
	size_t p_vertex_id[2]; /**< @brief ids of referenced vertices */
	TVertex3D *p_vertex[2]; /**< @brief pointers to the referenced vertices */
	VectorXd_constrained v_measurement; /**< @brief the measurement */
	MatrixXd_constrained t_sigma_inv; /**< @brief information matrix */
	MatrixXd_constrained t_square_root_sigma_inv_upper; /**< @brief the R matrix (upper diagonal) = transpose(chol(t_sigma_inv)) */

	// the below vectors / matrices are allocated/filled on the run, can be swapped for reference matrices
	VectorXd_constrained v_error; /**< @brief error vector (needs to be filled explicitly by calling p_error_function) */
	VectorXd_constrained v_expectation; /**< @brief expectation vector (needs to be filled explicitly by calling p_error_function) */
	MatrixXd_constrained p_jacobian[2]; /**< @brief the jacobians per either vertex (Ai and Bi in the french code) */ // t_odo - make them an array[2]

	double *p_RH[2]; /**< @brief pointers to jacobians inside A (for block slam) */

	void (*p_error_function)(TEdge3D &r_t_edge); /**< @brief calculates jacobians and error; r_t_edge is reference to this edge */

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // __3D_SOLVER_BASE_COMPILE_LEGACY_CODE

/**
 *	@brief base class for a 3D solver
 */
class CBase3DSolver /*: public COptimizationMethod<TVertex3D, TEdge3D>*/ {
public:
	/**
	 *	@brief implementation of Jacobian calculations, required by 3D solvers
	 */
	class C3DJacobians {
	public:
#if 0
		static Eigen::Matrix3d Operator_rot(double x, double y, double z) // it is a vector, treat it like that
		{
			double angle = sqrt(x*x + y*y + z*z);

			if(angle > 0) {
				x = x / angle;
				y = y / angle;
				z = z / angle;
			}

			double s = sin(angle);
			double c = cos(angle);
			double v = 1 - c;

			double xyv = x * y * v;
			double yzv = y * z * v;
			double xzv = x * z * v;

			Eigen::Matrix3d Q;
			Q << x * x * v + c, xyv - z * s, xzv + y * s,
				 xyv + z * s, y * y * v + c, yzv - x * s,
				 xzv - y * s, yzv + x * s, z * z * v + c;

			return Q;
		}
#endif // 0

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
		}

		/*static Eigen::Vector3d Operator_arot(Eigen::Matrix3d  &Q)
		{
			double c = (Q.trace() - 1.0) / 2.0;
			if(c < -1)
				c = -1;
			else if(c > 1)
				c = 1;

			double angle = acos(c);

			double factor;
			if(angle != 0)
				factor = angle / (2 * sin(angle));
			else
				factor = 1 / 2.0;

			Eigen::Vector3d axis(factor * (Q(2,1) - Q(1,2)), factor * (Q(0,2) - Q(2,0)), factor * (Q(1,0) - Q(0,1)));

			//fprintf(stderr,"%f %f %f\n", axis(0),axis(1),axis(2));

			return axis;
		}*/

		/**
		 *	@brief converts from RYP rep to axis angle
		 *	@param[in] Q is the rotation matrix
		 *	@return Returns axis-amgle rotation, where the angle is encoded as magnitude of the axis.
		 */
		static Eigen::Vector3d Operator_arot(const Eigen::Matrix3d &Q) // note that Eigen probably implements this somewhere
		{
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
		}

#ifdef __3D_SOLVER_BASE_COMPILE_LEGACY_CODE

		/**
		 *	@brief operator plus for 3D vertices
		 *	@param[in,out] r_t_vertex is the vertex to be modified
		 *	@param[in] p_global_deltas is the array of deltas, indexed by the permutation function
		 */
		static void Vertex3D_Plus(TVertex3D &r_t_vertex, const double *p_global_deltas)
		{
			_ASSERTE(r_t_vertex.n_vertex_dimension == 6);

			/* TODO: make more efficient */

			p_global_deltas += r_t_vertex.n_x_vector_first_elem_index;

			r_t_vertex._v_state(0) += p_global_deltas[0/*r_t_vertex.X_vector_indices[0]*/];
			r_t_vertex._v_state(1) += p_global_deltas[1/*r_t_vertex.X_vector_indices[1]*/];
			r_t_vertex._v_state(2) += p_global_deltas[2/*r_t_vertex.X_vector_indices[1]*/];

			//sum the rotations
			Eigen::Matrix3d pQ = Operator_rot(r_t_vertex._v_state(3), r_t_vertex._v_state(4), r_t_vertex._v_state(5));
			Eigen::Matrix3d dQ = Operator_rot(p_global_deltas[3], p_global_deltas[4], p_global_deltas[5]);

			Eigen::Matrix3d QQ = (pQ * dQ);
			Eigen::Vector3d axis = Operator_arot(QQ);

			r_t_vertex._v_state(3) = axis(0);
			r_t_vertex._v_state(4) = axis(1);
			r_t_vertex._v_state(5) = axis(2);
		}

		/**
		 *	@brief operator minus for 3D vertices
		 *	@param[in,out] r_t_vertex is the vertex to be modified
		 *	@param[in] p_global_deltas is the array of deltas, indexed by the permutation function
		 */
		static void Vertex3D_Minus(TVertex3D &r_t_vertex, const double *p_global_deltas)
		{
			_ASSERTE(r_t_vertex.n_vertex_dimension == 6);

			/* TODO: make more efficient, this is probable very inefficient */

			p_global_deltas += r_t_vertex.n_x_vector_first_elem_index;

			r_t_vertex._v_state(0) -= p_global_deltas[0/*r_t_vertex.X_vector_indices[0]*/];
			r_t_vertex._v_state(1) -= p_global_deltas[1/*r_t_vertex.X_vector_indices[1]*/];
			r_t_vertex._v_state(2) -= p_global_deltas[2/*r_t_vertex.X_vector_indices[1]*/];

			//sum the rotations
			Eigen::Matrix3d pQ = Operator_rot(r_t_vertex._v_state(3), r_t_vertex._v_state(4), r_t_vertex._v_state(5));
			Eigen::Matrix3d dQ = Operator_rot(p_global_deltas[3], p_global_deltas[4], p_global_deltas[5]);
			Eigen::Matrix3d dQ_inv = dQ.inverse();

			Eigen::Matrix3d QQ = pQ * dQ_inv;
			Eigen::Vector3d axis = Operator_arot(QQ);

			r_t_vertex._v_state(3) = axis(0);
			r_t_vertex._v_state(4) = axis(1);
			r_t_vertex._v_state(5) = axis(2);
		}

#endif // __3D_SOLVER_BASE_COMPILE_LEGACY_CODE

		/**
		 *	@brief operator plus for 3D vertices
		 *
		 *	@param[in] r_t_vertex1 is the vertex to be modified
		 *	@param[in] r_t_vertex2 is the delta vector
		 *	@param[out] r_t_dest is the result of the operation
		 */
		static void Smart_Plus(const Eigen::Matrix<double, 6, 1> &r_t_vertex1,
			const Eigen::Matrix<double, 6, 1> &r_t_vertex2, Eigen::Matrix<double, 6, 1> &r_t_dest)
		{
			_ASSERTE(r_t_vertex1.rows() == 6 && r_t_vertex2.rows() == 6);

			/* TODO: make more efficient */

			//r_t_dest.resize(6, 1); // no need to resize, it is 6D vector at compile time
			//r_t_dest(0) = r_t_vertex1(0) + r_t_vertex2(0);
			//r_t_dest(1) = r_t_vertex1(1) + r_t_vertex2(1);
			//r_t_dest(2) = r_t_vertex1(2) + r_t_vertex2(2);
			r_t_dest.head<3>() = r_t_vertex1.head<3>() + r_t_vertex2.head<3>(); // accelerated using SSE

			//sum the rotations
			Eigen::Matrix3d pQ = Operator_rot(r_t_vertex1.tail<3>());
			Eigen::Matrix3d dQ = Operator_rot(r_t_vertex2.tail<3>());

			Eigen::Matrix3d QQ;// = pQ * dQ;
			QQ.noalias() = pQ * dQ; // multiplication without intermediate storage
			//Eigen::Vector3d axis = Operator_arot(QQ);
			r_t_dest.tail<3>() = Operator_arot(QQ);

			//r_t_dest(3) = axis(0);
			//r_t_dest(4) = axis(1);
			//r_t_dest(5) = axis(2);
		}

		/**
		 *	@brief converts xyzaxis angles coordinates from relative measurement to absolute measurement
		 *
		 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
		 *	@param[in] r_t_vertex2 is the second vertex, relative to the first one
		 *	@param[out] r_t_dest is filled with absolute coordinates of the second vertex
		 */
		static void Relative_to_Absolute(const Eigen::Matrix<double, 6, 1> &r_t_vertex1,
			const Eigen::Matrix<double, 6, 1> &r_t_vertex2, Eigen::Matrix<double, 6, 1> &r_t_dest)
		{
			/* TODO: make more efficient */
			//r_t_dest.resize(6, 1); // no need to resize fixed size expressions
			Eigen::Vector3d p = r_t_vertex1.head(3);
			Eigen::Matrix3d pQ = Operator_rot(r_t_vertex1.tail<3>());
			Eigen::Vector3d d = r_t_vertex2.head(3);
			Eigen::Matrix3d dQ = Operator_rot(r_t_vertex2.tail<3>());
			r_t_dest.head<3>() = p + pQ * d;

			Eigen::Matrix3d QQ;// = pQ * dQ;
			QQ.noalias() = pQ * dQ; // multiplication without intermediate storage
			r_t_dest.tail<3>() = Operator_arot(QQ);
		}

		/**
		 *	@brief converts xyz axis angle coordinates from absolute measurement to relative measurement
		 *		and calculates the jacobians
		 *
		 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
		 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates
		 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
		 */
		static void Absolute_to_Relative(const Eigen::Matrix<double, 6, 1> &r_t_vertex1,
			const Eigen::Matrix<double, 6, 1> &r_t_vertex2, Eigen::Matrix<double, 6, 1> &r_t_dest)
		{
			/* TODO: make more efficient */
			Eigen::Vector3d p1 = r_t_vertex1.head(3);
			Eigen::Matrix3d pQ1 = Operator_rot(r_t_vertex1.tail<3>());

			Eigen::Vector3d p2 = r_t_vertex2.head(3);
			Eigen::Matrix3d pQ2 = Operator_rot(r_t_vertex2.tail<3>());

			Eigen::Matrix3d pQ1_inv = pQ1.inverse();
			//Eigen::Matrix3d pQ2_inv = pQ2.inverse();

			//r_t_dest.resize(6, 1); // no need to resize fixed size expressions
			r_t_dest.head(3) = pQ1_inv * (p2 - p1);

			Eigen::Matrix3d QQ = pQ1_inv * pQ2;
			r_t_dest.tail(3) = Operator_arot(QQ);	//TODO: is this right??
		}

		/**
		 *	@brief converts xyz axis angle coordinates from absolute measurement to relative measurement
		 *		and calculates the jacobians
		 *
		 *	@param[in] r_t_vertex1 is the first vertex, in absolute coordinates
		 *	@param[in] r_t_vertex2 is the second vertex, also in absolute coordinates
		 *	@param[out] r_t_dest is filled with relative coordinates of the second vertex
		 *	@param[out] r_t_pose3_pose1 is filled with the first jacobian
		 *	@param[out] r_t_pose3_pose2 is filled with the second jacobian
		 */
		template <class _TyDestVector, class _TyDestMatrix0, class _TyDestMatrix1> // want to be able to call this with differnent dest types (sometimes generic VectorXd / MatrixXd, sometimes with Vector3d / Matrix3d)
		static void Absolute_to_Relative(const Eigen::Matrix<double, 6, 1> &r_t_vertex1,
			const Eigen::Matrix<double, 6, 1> &r_t_vertex2, _TyDestVector &r_t_dest,
			_TyDestMatrix0 &r_t_pose3_pose1, _TyDestMatrix1 &r_t_pose3_pose2)
		{
			/* TODO: make more efficient */

			//lets try it according to g2o
			const double delta = 1e-9;
			const double scalar = 1.0 / (delta);

			//jacobians have to be computed here
			//double eps = 1.e-9; //TODO: what is this???
			//double eps_ = sqrt(eps);
			//double eps_2 = 7.45e-9;

			Eigen::Matrix<double, 6, 6> Eps;// = delta * Eigen::MatrixXd::Identity(6, 6); // MatrixXd needs to allocate storage on heap ... many milliseconds lost
			//Eigen::Matrix<double, 6, 6> H1;// = Eigen::MatrixXd::Zero(6, 6);
			//Eigen::Matrix<double, 6, 6> H2;// = Eigen::MatrixXd::Zero(6, 6);
			Eps = Eigen::Matrix<double, 6, 6>::Identity() * delta; // faster, all memory on stack
			//H1.setZero();
			//H2.setZero(); // faster, no extra memory; actually not needed at all, the below code overwrites both H1 and H2

			_TyDestMatrix0 &H1 = r_t_pose3_pose1;
			_TyDestMatrix1 &H2 = r_t_pose3_pose2;
			_ASSERTE(H1.rows() == 6 && H1.cols() == 6 && H2.rows() == 6 && H2.cols() == 6); // just make sure the shape is right
			// can actually work inplace

			//Eigen::Matrix<double, 6, 1> d;
			Absolute_to_Relative(r_t_vertex1, r_t_vertex2, r_t_dest);
			//r_t_dest = d; // possibly an unnecessary copy

			for(int j = 0; j < 6; ++ j) {
				Eigen::Matrix<double, 6, 1> d1, p_delta;
				Smart_Plus(r_t_vertex1, Eps.block(0, j, 6, 1), p_delta);
				Absolute_to_Relative(p_delta, r_t_vertex2, d1);
				H1.block(0, j, 6, 1) = (d1 - r_t_dest) * scalar;

				Eigen::Matrix<double, 6, 1> d2;
				Smart_Plus(r_t_vertex2, Eps.block(0, j, 6, 1), p_delta);
				Absolute_to_Relative(r_t_vertex1, p_delta, d2);
				H2.block(0, j, 6, 1) = (d2 - r_t_dest) * scalar;
			}
			// return jacobs
		}

#ifdef __3D_SOLVER_BASE_COMPILE_LEGACY_CODE

		/**
		 *	@brief vertex initialization functor
		 *	Calculates vertex position from the first vertex and an XYT edge.
		 */
		class CRelative_to_Absolute_XYZ_Initializer {
		protected:
			const Eigen::VectorXd &m_r_v_pose1; /**< @brief the first vertex */
			const CParser::CParseEntity_XYZ_Edge_3D &m_r_edge; /**< @brief the edge, shared by r_v_vertex1 and the vertex being initialized */

		public:
			/**
			 *	@brief default constructor
			 *	@param[in] r_v_vertex1 is the first vertex
			 *	@param[in] r_edge is the edge, shared by r_v_vertex1 and the vertex being initialized
			 */
			inline CRelative_to_Absolute_XYZ_Initializer(const Eigen::VectorXd &r_v_vertex1,
				const CParser::CParseEntity_XYZ_Edge_3D &r_edge)
				:m_r_v_pose1(r_v_vertex1), m_r_edge(r_edge)
			{}

			/**
			 *	@brief function operator
			 *	@return Returns the value of the vertex being initialized.
			 */
			inline TVertex3D operator ()() const
			{
				TVertex3D t_pose2;
				t_pose2.n_vertex_dimension = 6;
				C3DJacobians::Relative_to_Absolute(m_r_v_pose1, m_r_edge.m_v_delta, t_pose2._v_state);

				t_pose2.p_ordering_function = NULL;
				t_pose2.p_operator_plus = &Vertex3D_Plus;
				t_pose2.p_operator_minus = &Vertex3D_Minus;
				return t_pose2;
			}
		};

		/**
		 *	@brief a simple vertex initialization functor (initializes state vector to null)
		 *	@note Note that this could be a static function, but whatevs ...
		 */
		class C3DVertex_Null_Initializer {
		public:
			/**
			 *	@brief function operator
			 *	@return Returns the value of the vertex being initialized.
			 */
			inline TVertex3D operator ()() const
			{
				TVertex3D t_vertex;
				t_vertex.n_vertex_dimension = 6;
				for(int i = 0; i < 6; ++ i)
					t_vertex._v_state(i) = 0;
				t_vertex.p_ordering_function = NULL;
				t_vertex.p_operator_plus = &Vertex3D_Plus;
				t_vertex.p_operator_minus = &Vertex3D_Minus;
				return t_vertex;
			}
		};

		/**
		 *	@brief calculates jacobians and updates error vector for an XYZ edge
		 *	@param[in,out] r_t_edge is the edge to be updated
		 */
		static void Calculate_Jacobians_Error_XYZ(TEdge3D &r_t_edge)
		{
			_ASSERTE(r_t_edge.p_vertex[0] && r_t_edge.p_vertex[1]);
			_ASSERTE(r_t_edge.p_vertex[0]->n_vertex_dimension == 6);
			_ASSERTE(r_t_edge.p_vertex[1]->n_vertex_dimension == 6);
			const Eigen::VectorXd &p1 = r_t_edge.p_vertex[0]->_v_state;
			const Eigen::VectorXd &p2 = r_t_edge.p_vertex[1]->_v_state;

			Absolute_to_Relative(p1, p2, r_t_edge.v_expectation, r_t_edge.p_jacobian[0], r_t_edge.p_jacobian[1]); // calculates "h"
			// calculates the expectation and the jacobians

			//r_t_edge.v_error.head<3> = r_t_edge.v_measurement.head<3> - r_t_edge.v_expectation.head<3>; // r_t_edge.v_measurement is z

			r_t_edge.v_error(0) = p1(0) - p2(0);
			r_t_edge.v_error(1) = p1(1) - p2(1);
			r_t_edge.v_error(2) = p1(2) - p2(2);

			//sum the rotations
			Eigen::Matrix3d pQ = Operator_rot(p1(3), p1(4), p1(5));
			Eigen::Matrix3d dQ = Operator_rot(p2(3), p2(4), p2(5));
			Eigen::Matrix3d dQ_inv = dQ.inverse();

			Eigen::Matrix3d QQ = pQ * dQ_inv;
			Eigen::Vector3d axis = Operator_arot(QQ);

			r_t_edge.v_error(3) = axis(0);
			r_t_edge.v_error(4) = axis(1);
			r_t_edge.v_error(5) = axis(2);
			//TODO: fix the angle somehow? this cannot be right

			// calculates error
		}

#endif // __3D_SOLVER_BASE_COMPILE_LEGACY_CODE
	};

public:
#ifdef __3D_SOLVER_BASE_COMPILE_LEGACY_CODE

	/**
	 *	@copydoc COptimizationMethod<TVertex3D, TEdge3D>::COptimizationMethod
	 */
	CBase3DSolver(COptimizationSystem<TVertex3D, TEdge3D> &r_system)
		:COptimizationMethod<TVertex3D, TEdge3D>(r_system)
	{}

	/**
	 *	@copydoc CParser::CParserAdaptor::AppendSystem(const CParser::TEdge3D&)
	 *	@return Returns reference to the new edge.
	 */
	TEdge3D &r_AppendSystemWith_XYZ_Edge(const CParser::CParseEntity_XYZ_Edge_3D &r_t_edge) // throws(std::bad_alloc)
	{
		TVertex3D &r_t_vertex_0 = m_r_system.r_Get_Vertex(r_t_edge.m_n_node_0,
			C3DJacobians::C3DVertex_Null_Initializer());
		Eigen::VectorXd &r_v_pose1 = r_t_vertex_0._v_state;
		TVertex3D &r_t_vertex_1 = m_r_system.r_Get_Vertex(r_t_edge.m_n_node_1,
			C3DJacobians::CRelative_to_Absolute_XYZ_Initializer(r_v_pose1, r_t_edge));
		// gets the vertices, calculates the coordinates of the second
		// one based on the measurement (assumes the vertices come in order)

		TEdge3D &r_new_edge = m_r_system.r_Add_Edge(3 * ((m_r_system.r_Edge_Map().empty())? 2 : 1)); // the first edge will bear the unary factor
		r_new_edge.p_vertex_id[0] = r_t_edge.m_n_node_0;
		r_new_edge.p_vertex_id[1] = r_t_edge.m_n_node_1;
		r_new_edge.p_vertex[0] = &r_t_vertex_0;
		r_new_edge.p_vertex[1] = &r_t_vertex_1;
		r_new_edge.v_measurement = r_t_edge.m_v_delta; // t_odo - maybe this doesn't assign (v_measurement is dynamically sized)
		r_new_edge.v_error.resize(6,1);
		r_new_edge.v_error << 0, 0, 0, 0, 0, 0;
		r_new_edge.v_expectation.resize(6,1);
		r_new_edge.v_expectation << 0, 0, 0, 0, 0, 0;
		r_new_edge.t_sigma_inv = r_t_edge.m_t_inv_sigma; // t_odo - it might be transpose // no it's not, it's a (major) diagonal matrix
		r_new_edge.t_square_root_sigma_inv_upper = r_t_edge.m_t_inv_sigma.llt().matrixU(); // eigen calculates upper cholesky (U), no need to transpose // t_odo - fixme
		r_new_edge.p_jacobian[0].resize(6, 6); // Ai in the french code
		r_new_edge.p_jacobian[1].resize(6, 6); // Bi in the french code
		r_new_edge.p_error_function = &C3DJacobians::Calculate_Jacobians_Error_XYZ;
		// add measurement

		// note the width of the jacobians is supposed to be equal to measurement vector dimensions,
		// and the height should equal vertex dimensions (based on addEntriesInSparseSystem)

		return r_new_edge;
	}

#endif // __3D_SOLVER_BASE_COMPILE_LEGACY_CODE
};

#endif // __3D_SOLVER_BASE_INCLUDED
