/*
								+----------------------------------+
								|                                  |
								| *** Base class for 2D solver *** |
								|                                  |
								|   Copyright © -tHE SWINe- 2012   |
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
#ifdef __2D_SOLVER_BASE_COMPILE_LEGACY_CODE
#include "slam/System.h"
#endif // __2D_SOLVER_BASE_COMPILE_LEGACY_CODE

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

#ifdef __2D_SOLVER_BASE_COMPILE_LEGACY_CODE

/**
 *	@brief a 2D vertex (a vertex for 2D problems)
 */
struct TVertex2D { friend class CBase2DSolver; // CBase2DSolver can access _v_state directly, it knows what it is doing (derived classes such as IncrementalCholesky.* can't)
protected:
	Eigen::Vector3d _v_state; /**< @brief the state vector (actually stores 2D or 3D vector, that's why it is not supposed to be accessed directly) */ // t_odo - make this an array of doubles. that should make things much easier // no need, superseded by the flat system

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

	void (*p_ordering_function)(TVertex2D &r_t_vertex, size_t &r_n_x_index); /**< @brief calculates indices into X vector */
	void (*p_operator_plus)(TVertex2D &r_t_vertex, const double *p_global_deltas); /**< @brief should add v_state += p_global_deltas[X_vector_indices] */
	void (*p_operator_minus)(TVertex2D &r_t_vertex, const double *p_global_deltas); /**< @brief should subtract v_state -= p_global_deltas[X_vector_indices] */

	// t_odo - convert these "virtual" functions to members; there's always (AFAIK) just a single kind of vertex in the system so it may as well directly implement those
	// there are *two* kinds of vertices, vertices (3D) and landmarks (2D)
	// t_odo - maybe it would be possible to use virtual functions if the derived classes had the same size // it is

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 *	@brief a 2D edge (an edge for 2D problems)
 */
struct TEdge2D {
	size_t p_vertex_id[2]; /**< @brief ids of referenced vertices */
	TVertex2D *p_vertex[2]; /**< @brief pointers to the referenced vertices */
	VectorXd_constrained v_measurement; /**< @brief the measurement (2D for range, 3D for range bearing) */
	MatrixXd_constrained t_sigma_inv; /**< @brief information matrix (2x2 or 3x3) */
	MatrixXd_constrained t_square_root_sigma_inv_upper; /**< @brief the R matrix (upper diagonal) = transpose(chol(t_sigma_inv)) */

	// the below vectors / matrices are allocated/filled on the run, can be swapped for reference matrices
	VectorXd_constrained v_error; /**< @brief error vector (needs to be filled explicitly by calling p_error_function) */
	VectorXd_constrained v_expectation; /**< @brief expectation vector (needs to be filled explicitly by calling p_error_function) */
	MatrixXd_constrained p_jacobian[2]; /**< @brief the jacobians per either vertex (Ai and Bi in the french code) */ // t_odo - make them an array[2]

	double *p_RH[2]; /**< @brief pointers to jacobians inside A (for block slam) */ // t_odo - write CBase2DSolver_Blocky and TEdge2D_Blocky and all that stuff // no. superseded by flat system

	void (*p_error_function)(TEdge2D &r_t_edge); /**< @brief calculates jacobians and error; r_t_edge is reference to this edge */

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// t_odo - decide how to cope with 3D problems: template the solvers (+speed, -flexibility) or create more generic edges and vertices to be able to hold 3D data (-speed, +flexibility) // will use templates + flat system

#endif // __2D_SOLVER_BASE_COMPILE_LEGACY_CODE

/**
 *	@brief base class for a 2D solver
 */
class CBase2DSolver
#ifdef __2D_SOLVER_BASE_COMPILE_LEGACY_CODE
	: public COptimizationMethod<TVertex2D, TEdge2D>
#endif // __2D_SOLVER_BASE_COMPILE_LEGACY_CODE
	{
public:
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

#ifdef __2D_SOLVER_BASE_COMPILE_LEGACY_CODE

		/**
		 *	@brief operator plus for 2D vertices
		 *	@param[in,out] r_t_vertex is the vertex to be modified
		 *	@param[in] p_global_deltas is the array of deltas, indexed by the permutation function
		 */
		static void Vertex2D_Plus(TVertex2D &r_t_vertex, const double *p_global_deltas)
		{
			_ASSERTE(r_t_vertex.n_vertex_dimension == 3);
			p_global_deltas += r_t_vertex.n_x_vector_first_elem_index;
			r_t_vertex._v_state(0) += p_global_deltas[0/*r_t_vertex.X_vector_indices[0]*/];
			r_t_vertex._v_state(1) += p_global_deltas[1/*r_t_vertex.X_vector_indices[1]*/];
			r_t_vertex._v_state(2) = f_ClampAngle_2Pi(r_t_vertex._v_state(2) + p_global_deltas[2/*r_t_vertex.X_vector_indices[2]*/]);
		}

		/**
		 *	@brief operator minus for 2D vertices
		 *	@param[in,out] r_t_vertex is the vertex to be modified
		 *	@param[in] p_global_deltas is the array of deltas, indexed by the permutation function
		 */
		static void Vertex2D_Minus(TVertex2D &r_t_vertex, const double *p_global_deltas)
		{
			_ASSERTE(r_t_vertex.n_vertex_dimension == 3);
			p_global_deltas += r_t_vertex.n_x_vector_first_elem_index;
			r_t_vertex._v_state(0) -= p_global_deltas[0/*r_t_vertex.X_vector_indices[0]*/];
			r_t_vertex._v_state(1) -= p_global_deltas[1/*r_t_vertex.X_vector_indices[1]*/];
			r_t_vertex._v_state(2) = f_ClampAngle_2Pi(r_t_vertex._v_state(2) - p_global_deltas[2/*r_t_vertex.X_vector_indices[2]*/]);
		}

		/**
		 *	@brief operator plus for 2D landmarks
		 *	@param[in,out] r_t_vertex is the landmark to be modified
		 *	@param[in] p_global_deltas is the array of deltas, indexed by the permutation function
		 */
		static void Landmark2D_Plus(TVertex2D &r_t_vertex, const double *p_global_deltas)
		{
			_ASSERTE(r_t_vertex.n_vertex_dimension == 2);
			p_global_deltas += r_t_vertex.n_x_vector_first_elem_index;
			r_t_vertex._v_state(0) += p_global_deltas[0/*r_t_vertex.X_vector_indices[0]*/];
			r_t_vertex._v_state(1) += p_global_deltas[1/*r_t_vertex.X_vector_indices[1]*/];
		}

		/**
		 *	@brief operator minus for 2D landmarks
		 *	@param[in,out] r_t_vertex is the landmark to be modified
		 *	@param[in] p_global_deltas is the array of deltas, indexed by the permutation function
		 */
		static void Landmark2D_Minus(TVertex2D &r_t_vertex, const double *p_global_deltas)
		{
			_ASSERTE(r_t_vertex.n_vertex_dimension == 2);
			p_global_deltas += r_t_vertex.n_x_vector_first_elem_index;
			r_t_vertex._v_state(0) -= p_global_deltas[0/*r_t_vertex.X_vector_indices[0]*/];
			r_t_vertex._v_state(1) -= p_global_deltas[1/*r_t_vertex.X_vector_indices[1]*/];
		}

#endif // __2D_SOLVER_BASE_COMPILE_LEGACY_CODE

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

#ifdef __2D_SOLVER_BASE_COMPILE_LEGACY_CODE

		/**
		 *	@brief updates permutation function by adding a new vertex
		 *
		 *	@param[in,out] r_t_vertex is the new vertex
		 *	@param[in,out] r_n_x_index is index of the next free slot in the matrix
		 */
		static void Vertex2D_OrderingFunction(TVertex2D &r_t_vertex, size_t &r_n_x_index) // t_odo - add element / block index // no need
		{
			if(1/*!r_t_vertex.fixedValue &&*/) {
				_ASSERTE(r_t_vertex.n_vertex_dimension == 3);
				//int sizeXi = r_t_vertex.n_vertex_dimension;
				r_t_vertex.n_x_vector_first_elem_index = r_n_x_index;
				r_n_x_index += 3/*sizeXi*/;
			}
		}

		/**
		 *	@brief updates permutation function by adding a new landmark
		 *
		 *	@param[in,out] r_t_vertex is the new landmark
		 *	@param[in,out] r_n_x_index is index of the next free slot in the matrix
		 */
		static void Landmark2D_OrderingFunction(TVertex2D &r_t_vertex, size_t &r_n_x_index) // t_odo - add element / block index // no need
		{
			if(1/*!r_t_vertex.fixedValue &&*/) {
				_ASSERTE(r_t_vertex.n_vertex_dimension == 2);
				//int sizeXi = r_t_vertex.n_vertex_dimension;
				r_t_vertex.n_x_vector_first_elem_index = r_n_x_index;
				r_n_x_index += 2/*sizeXi*/;
			}
		}

		/**
		 *	@brief vertex initialization functor
		 *	Calculates vertex position from the first vertex and a range-bearing edge.
		 */
		class CRelative_to_Absolute_RangeBearing_Initializer { // t_odo - optimize for z=0 // can't
		protected:
			const Eigen::Vector3d &m_r_v_pose1; /**< @brief the first vertex */
			const CParser::CParseEntity_RangeBearing_Edge_2D &m_r_edge; /**< @brief the edge, shared by r_v_vertex1 and the vertex being initialized */

		public:
			/**
			 *	@brief default constructor
			 *	@param[in] r_v_vertex1 is the first vertex
			 *	@param[in] r_edge is the edge, shared by r_v_vertex1 and the vertex being initialized
			 */
			inline CRelative_to_Absolute_RangeBearing_Initializer(const Eigen::Vector3d &r_v_vertex1,
				const CParser::CParseEntity_RangeBearing_Edge_2D &r_edge)
				:m_r_v_pose1(r_v_vertex1), m_r_edge(r_edge)
			{}

			/**
			 *	@brief function operator
			 *	@return Returns the value of the vertex being initialized.
			 */
			inline TVertex2D operator ()() const
			{
				Eigen::Vector3d pose1(m_r_v_pose1(0), m_r_v_pose1(1), m_r_v_pose1(2)); // fixme - is this correct? is it 3D or just 2D?
				Eigen::Vector3d relPose(m_r_edge.m_v_delta(0), m_r_edge.m_v_delta(1), 0); // 2D to 3D (append 0)
				TVertex2D lmk3dPose;
				lmk3dPose.n_vertex_dimension = 2;
				C2DJacobians::Relative_to_Absolute(pose1, relPose, lmk3dPose._v_state);
				lmk3dPose._v_state(2) = 0; // debug
				lmk3dPose.p_ordering_function = &Landmark2D_OrderingFunction;
				lmk3dPose.p_operator_plus = &Landmark2D_Plus;
				lmk3dPose.p_operator_minus = &Landmark2D_Minus;
				return lmk3dPose;
			}
		};

		/**
		 *	@brief vertex initialization functor
		 *	Calculates vertex position from the first vertex and an XYT edge.
		 */
		class CRelative_to_Absolute_XYT_Initializer {
		protected:
			const Eigen::Vector3d &m_r_v_pose1; /**< @brief the first vertex */
			const CParser::CParseEntity_XYT_Edge_2D &m_r_edge; /**< @brief the edge, shared by r_v_vertex1 and the vertex being initialized */

		public:
			/**
			 *	@brief default constructor
			 *	@param[in] r_v_vertex1 is the first vertex
			 *	@param[in] r_edge is the edge, shared by r_v_vertex1 and the vertex being initialized
			 */
			inline CRelative_to_Absolute_XYT_Initializer(const Eigen::Vector3d &r_v_vertex1,
				const CParser::CParseEntity_XYT_Edge_2D &r_edge)
				:m_r_v_pose1(r_v_vertex1), m_r_edge(r_edge)
			{}

			/**
			 *	@brief function operator
			 *	@return Returns the value of the vertex being initialized.
			 */
			inline TVertex2D operator ()() const
			{
				TVertex2D t_pose2;
				t_pose2.n_vertex_dimension = 3;
				C2DJacobians::Relative_to_Absolute(m_r_v_pose1, m_r_edge.m_v_delta, t_pose2._v_state);
				t_pose2.p_ordering_function = &Vertex2D_OrderingFunction;
				t_pose2.p_operator_plus = &Vertex2D_Plus;
				t_pose2.p_operator_minus = &Vertex2D_Minus;
				return t_pose2;
			}
		};

		/**
		 *	@brief a simple vertex initialization functor (initializes state vector to null)
		 *	@note Note that this could be a static function, but whatevs ...
		 */
		class C2DVertex_Null_Initializer {
		public:
			/**
			 *	@brief function operator
			 *	@return Returns the value of the vertex being initialized.
			 */
			inline TVertex2D operator ()() const
			{
				TVertex2D t_vertex;
				t_vertex.n_vertex_dimension = 3;
				for(int i = 0; i < 3; ++ i)
					t_vertex._v_state(i) = 0;
				t_vertex.p_ordering_function = &Vertex2D_OrderingFunction;
				t_vertex.p_operator_plus = &Vertex2D_Plus;
				t_vertex.p_operator_minus = &Vertex2D_Minus;
				return t_vertex;
			}
		};

		/**
		 *	@brief calculates jacobians and updates error vector for an XYT edge
		 *	@param[in,out] r_t_edge is the edge to be updated
		 */
		static void Calculate_Jacobians_Error_XYT(TEdge2D &r_t_edge)
		{
			_ASSERTE(r_t_edge.p_vertex[0] && r_t_edge.p_vertex[1]);
			_ASSERTE(r_t_edge.p_vertex[0]->n_vertex_dimension == 3);
			_ASSERTE(r_t_edge.p_vertex[1]->n_vertex_dimension == 3);
			const Eigen::Vector3d &p1 = r_t_edge.p_vertex[0]->_v_state;
			const Eigen::Vector3d &p2 = r_t_edge.p_vertex[1]->_v_state;

			Absolute_to_Relative(p1, p2, r_t_edge.v_expectation, r_t_edge.p_jacobian[0], r_t_edge.p_jacobian[1]); // calculates "h"
			// calculates the expectation and the jacobians

			r_t_edge.v_error = r_t_edge.v_measurement - r_t_edge.v_expectation; // r_t_edge.v_measurement is z
			r_t_edge.v_error(2) = f_ClampAngle_2Pi(r_t_edge.v_error(2));
			// calculates error

			double res = f_MinimumAbsolute_3(r_t_edge.v_error(2), r_t_edge.v_error(2) - 2 * M_PI, r_t_edge.v_error(2) + 2 * M_PI); // t_odo - what the hell is going on here?
			//_ASSERTE(res == f_ClampAngle_2Pi(res)); // t_odo - if this passes, remove the last clamp below
			r_t_edge.v_error(2) = res;//f_ClampAngle_2Pi(res); // presumably no-op
			// fixup the error so the absolute value is the lowest possible (considering it is modulo 2pi)
		}

		/**
		 *	@brief calculates jacobians and updates error vector for a range-bearing edge
		 *	@param[in,out] r_t_edge is the edge to be updated
		 */
		static void Calculate_Jacobians_Error_RangeBearing(TEdge2D &r_t_edge)
		{
			_ASSERTE(r_t_edge.p_vertex[0] && r_t_edge.p_vertex[1]);
			_ASSERTE(r_t_edge.p_vertex[0]->n_vertex_dimension == 3);
			_ASSERTE(r_t_edge.p_vertex[1]->n_vertex_dimension == 2);
			const Eigen::Vector3d &r_v_pose = r_t_edge.p_vertex[0]->_v_state;
			const Eigen::Vector2d v_landmark(r_t_edge.p_vertex[1]->_v_state(0),
				r_t_edge.p_vertex[1]->_v_state(1)); // just 2D

			Observation2D_RangeBearing(r_v_pose, v_landmark, r_t_edge.v_expectation,
				r_t_edge.p_jacobian[0], r_t_edge.p_jacobian[1]);

			r_t_edge.v_error = r_t_edge.v_measurement - r_t_edge.v_expectation;
			r_t_edge.v_error(1) = f_ClampAngle_2Pi(r_t_edge.v_error(1));
			// calculates error

			double res = f_MinimumAbsolute_3(r_t_edge.v_error(1), r_t_edge.v_error(1) - 2 * M_PI, r_t_edge.v_error(1) + 2 * M_PI); // t_odo - what the hell is going on here?
			//_ASSERTE(res == f_ClampAngle_2Pi(res)); // t_odo - if this passes, remove the last clamp below
			r_t_edge.v_error(1) = res;//f_ClampAngle_2Pi(res); // presumably no-op
			// fixup the error so the absolute value is the lowest possible (considering it is modulo 2pi)
		}

#endif // __2D_SOLVER_BASE_COMPILE_LEGACY_CODE

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

public:
#ifdef __2D_SOLVER_BASE_COMPILE_LEGACY_CODE

	/**
	 *	@copydoc COptimizationMethod<TVertex2D, TEdge2D>::COptimizationMethod
	 */
	CBase2DSolver(COptimizationSystem<TVertex2D, TEdge2D> &r_system)
		:COptimizationMethod<TVertex2D, TEdge2D>(r_system)
	{}

	/**
	 *	@copydoc CParser::CParserAdaptor::AppendSystem(const CParser::TEdge2D&)
	 *	@return Returns reference to the new edge.
	 */
	TEdge2D &r_AppendSystemWith_RangeBearing_Edge(const CParser::CParseEntity_RangeBearing_Edge_2D &r_t_edge) // throws(std::bad_alloc)
	{
		TVertex2D &r_t_vertex_0 = m_r_system.r_Get_Vertex(r_t_edge.m_n_node_0,
			C2DJacobians::C2DVertex_Null_Initializer());
		Eigen::Vector3d &r_v_pose1 = r_t_vertex_0._v_state;
		TVertex2D &r_t_vertex_1 = m_r_system.r_Get_Vertex(r_t_edge.m_n_node_1,
			C2DJacobians::CRelative_to_Absolute_RangeBearing_Initializer(r_v_pose1, r_t_edge));
		// gets the vertices, calculates the coordinates of the second
		// one based on the measurement (assumes the vertices come in order)

		Eigen::Matrix2d lazerInvCov;
		for(int i = 0; i < 2; ++ i) {
			for(int j = 0; j < 2; ++ j)
				lazerInvCov(i, j) = (i == j)? 1 : 0;
		}
		// t_odo - ask ela why do this? why not use r_t_edge.m_t_inv_sigma instead? // because don't know how to convert covariance xy to rb

		double dx = r_t_edge.m_v_delta(0);
		double dy = r_t_edge.m_v_delta(1);
		Eigen::Vector2d v_laser_data(sqrt(dx * dx + dy * dy), C2DJacobians::f_ClampAngle_2Pi(atan2(dy, dx)));
		// ...

		TEdge2D &r_new_edge = m_r_system.r_Add_Edge(2 * ((m_r_system.r_Edge_Map().empty())? 2 : 1)); // the first edge will bear the unary factor
		r_new_edge.p_vertex_id[0] = r_t_edge.m_n_node_0;
		r_new_edge.p_vertex_id[1] = r_t_edge.m_n_node_1;
		r_new_edge.p_vertex[0] = &r_t_vertex_0;
		r_new_edge.p_vertex[1] = &r_t_vertex_1;
		r_new_edge.v_measurement = v_laser_data; // t_odo - maybe this doesn't assign (v_measurement is dynamically sized)
		r_new_edge.v_error = r_new_edge.v_expectation = Eigen::Vector2d(0, 0); // error and expectation have both the same dimensionality as measurement
		r_new_edge.t_sigma_inv = lazerInvCov;
		r_new_edge.t_square_root_sigma_inv_upper = lazerInvCov.llt().matrixU(); // eigen calculates upper cholesky (U), no need to transpose // t_odo - fixme // t_odo - this is constant. optimize away! // no need, really, this is superseded by flat system
		r_new_edge.p_jacobian[0].resize(2, 3); // Ai in the french code
		r_new_edge.p_jacobian[1].resize(2, 2); // Bi in the french code // t_odo - this is just a fixup (should be 2x2)
		r_new_edge.p_error_function = &C2DJacobians::Calculate_Jacobians_Error_RangeBearing;
		// add measurement

		// note the width of the jacobians is supposed to be equal to measurement vector dimensions,
		// and the height should equal vertex dimensions (based on addEntriesInSparseSystem)

		return r_new_edge;
	}

	/**
	 *	@copydoc CParser::CParserAdaptor::AppendSystem(const CParser::TEdge2D&)
	 *	@return Returns reference to the new edge.
	 */
	TEdge2D &r_AppendSystemWith_XYT_Edge(const CParser::CParseEntity_XYT_Edge_2D &r_t_edge) // throws(std::bad_alloc)
	{
		TVertex2D &r_t_vertex_0 = m_r_system.r_Get_Vertex(r_t_edge.m_n_node_0,
			C2DJacobians::C2DVertex_Null_Initializer());
		Eigen::Vector3d &r_v_pose1 = r_t_vertex_0._v_state;
		TVertex2D &r_t_vertex_1 = m_r_system.r_Get_Vertex(r_t_edge.m_n_node_1,
			C2DJacobians::CRelative_to_Absolute_XYT_Initializer(r_v_pose1, r_t_edge));
		// gets the vertices, calculates the coordinates of the second
		// one based on the measurement (assumes the vertices come in order)

		TEdge2D &r_new_edge = m_r_system.r_Add_Edge(3 * ((m_r_system.r_Edge_Map().empty())? 2 : 1)); // the first edge will bear the unary factor
		r_new_edge.p_vertex_id[0] = r_t_edge.m_n_node_0;
		r_new_edge.p_vertex_id[1] = r_t_edge.m_n_node_1;
		r_new_edge.p_vertex[0] = &r_t_vertex_0;
		r_new_edge.p_vertex[1] = &r_t_vertex_1;
		r_new_edge.v_measurement = r_t_edge.m_v_delta; // t_odo - maybe this doesn't assign (v_measurement is dynamically sized)
		r_new_edge.v_error = r_new_edge.v_expectation = Eigen::Vector3d(0, 0, 0); // error and expectation have both the same dimensionality as measurement
		r_new_edge.t_sigma_inv = r_t_edge.m_t_inv_sigma; // t_odo - it might be transpose // no it's not, it's a (major) diagonal matrix
		r_new_edge.t_square_root_sigma_inv_upper = r_t_edge.m_t_inv_sigma.llt().matrixU(); // eigen calculates upper cholesky (U), no need to transpose // t_odo - fixme
		r_new_edge.p_jacobian[0].resize(3, 3); // Ai in the french code
		r_new_edge.p_jacobian[1].resize(3, 3); // Bi in the french code
		r_new_edge.p_error_function = &C2DJacobians::Calculate_Jacobians_Error_XYT;
		// add measurement

		// note the width of the jacobians is supposed to be equal to measurement vector dimensions,
		// and the height should equal vertex dimensions (based on addEntriesInSparseSystem)

		return r_new_edge;
	}

#endif // __2D_SOLVER_BASE_COMPILE_LEGACY_CODE
};

#endif // __2D_SOLVER_BASE_INCLUDED
