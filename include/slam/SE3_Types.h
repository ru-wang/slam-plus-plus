/**
 *	@file include/slam/SE3_Types.h
 *	@author Soso
 *	@date 2013
 *	@brief projection type for SE3
 */

#pragma once
#ifndef __SE3_TYPES_INCLUDED
#define __SE3_TYPES_INCLUDED

#include "slam/BaseTypes.h"
#include "slam/SE2_Types.h"
//#include "slam/2DSolverBase.h"
#include "slam/3DSolverBase.h"

/**
 *	@brief SE(3) pose vertex type
 */
class CVertexPose3D : public CSEBaseVertexImpl<CVertexPose3D, 6> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
public:
	__SE3_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CVertexPose3D() // copy this
	{}

	/**
	 *	@brief constructor; initializes state vector
	 *	@param[in] r_v_state is state vector initializer
	 */
	inline CVertexPose3D(const Eigen::Matrix<double, 6, 1> &r_v_state) // copy this, change the dimension of the vector to appropriate
		:CSEBaseVertexImpl<CVertexPose3D, 6>(r_v_state) // change the dimension here as well
	{}

	/**
	 *	@copydoc CSEBaseVertex::Operator_Plus()
	 */
	inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
	{
		// todo - use the code from 3DSolverBase.h !!! duplicates code, causes bugs !!!

		//sum the rotations
		Eigen::Matrix3d pQ = C3DJacobians::Operator_rot(m_v_state.tail<3>());
		Eigen::Matrix3d dQ = C3DJacobians::Operator_rot(r_v_delta.segment<3>(m_n_order + 3));

		//m_v_state(0) += r_v_delta(m_n_order); // pick part of the delta vector, belonging to this vertex, apply +
		//m_v_state(1) += r_v_delta(m_n_order+1);
		//m_v_state(2) += r_v_delta(m_n_order+2);
		m_v_state.head<3>() += pQ * r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // noalias says that the destination is not the same as the multiplicands and that the multiplication can be carried out without allocating intermediate storage
		//Eigen::Vector3d axis = C3DJacobians::Operator_arot(QQ);
		//m_v_state.tail(3) = axis;
		//m_v_state.tail<3>() = axis; // we have the information that it is 3 dimensional at compile time, if we put it in the <> brackets, the compiler will avoid allocating dynamic storage for the intermediate
		m_v_state.tail<3>() = C3DJacobians::Operator_arot(QQ); // also, no need for the intermediate object
	}

	/**
	 *	@copydoc CSEBaseVertex::Operator_Minus()
	 */
	inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
	{
		Operator_Plus(-r_v_delta); // call plus with negative delta, that should do the trick
	}
};

/**
 *	@brief SE(3) pose-pose edge
 */
class CEdgePose3D : public CSEBaseEdgeImpl<CEdgePose3D, CVertexPose3D, CVertexPose3D, 6> { // again, this tells that base implementation for base edge for type that will be called CEdgePose2D, and it will be an edge between two vertices of type CVertexPose2D, and the measurement will have 3 dimensions (in order to have hyperedges, you need to provide your own version of CSEBaseEdgeImpl, that is an advanced topic)
public:
	/**
	 *	@brief vertex initialization functor
	 *	Calculates vertex position from the first vertex and an XYT edge.
	 */
	class CRelative_to_Absolute_XYZ_Initializer { // this is an object which is used to lazy initialize vertices (copy it)
	protected:
		const Eigen::Matrix<double, 6, 1> &m_r_v_pose1; /**< @brief the first vertex */
		const Eigen::Matrix<double, 6, 1> &m_r_v_delta; /**< @brief the edge, shared by r_v_vertex1 and the vertex being initialized */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in] r_v_vertex1 is the first vertex
		 *	@param[in] r_v_delta is the edge, shared by r_v_vertex1 and the vertex being initialized
		 */
		inline CRelative_to_Absolute_XYZ_Initializer(const Eigen::Matrix<double, 6, 1> &r_v_vertex1,
			const Eigen::Matrix<double, 6, 1> &r_v_delta) // just change the types, same as above
			:m_r_v_pose1(r_v_vertex1), m_r_v_delta(r_v_delta)
		{}

		/**
		 *	@brief function operator
		 *	@return Returns the value of the vertex being initialized.
		 */
		inline operator CVertexPose3D() const // this function calculates initial prior from the state of the first vertex m_r_v_pose1 and from the edge measurement m_r_edge
		{
			Eigen::Matrix<double, 6, 1> v_pose2;
			C3DJacobians::Relative_to_Absolute(m_r_v_pose1, m_r_v_delta, v_pose2); // implement your own equation here
			return CVertexPose3D(v_pose2);
		}
	};

public:
	__SE3_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, just copy this

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CEdgePose3D() // copy this
	{}

	/**
	 *	@brief constructor; converts parsed edge to edge representation
	 *
	 *	@tparam CSystem is type of system where this edge is being stored
	 *
	 *	@param[in] r_t_edge is parsed edge
	 *	@param[in,out] r_system is reference to system (used to query edge vertices)
	 */
	template <class CSystem>
	CEdgePose3D(const CParserBase::TEdge3D &r_t_edge, CSystem &r_system) // this is edge constructor how it is called in the parse loop; you need to change type of the edge
		:CSEBaseEdgeImpl<CEdgePose3D, CVertexPose3D, CVertexPose3D, 6>(r_t_edge.m_n_node_0,
		r_t_edge.m_n_node_1, r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		//fprintf(stderr, "%f %f %f\n", r_t_edge.m_v_delta(0), r_t_edge.m_v_delta(1), r_t_edge.m_v_delta(2));

		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexPose3D>(r_t_edge.m_n_node_0, CInitializeNullVertex());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexPose3D>(r_t_edge.m_n_node_1,
			CRelative_to_Absolute_XYZ_Initializer(m_p_vertex0->v_State(), r_t_edge.m_v_delta)); // rename your initializer if required

		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(((CSEBaseVertex*)m_p_vertex0)->n_Dimension() == 6);
		_ASSERTE(((CSEBaseVertex*)m_p_vertex1)->n_Dimension() == 6); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief constructor; converts parsed edge to edge representation
	 *
	 *	@tparam CSystem is type of system where this edge is being stored
	 *
	 *	@param[in] n_node0 is (zero-based) index of the first (origin) node
	 *	@param[in] n_node1 is (zero-based) index of the second (endpoint) node
	 *	@param[in] r_v_delta is measurement vector
	 *	@param[in] r_t_inv_sigma is the information matrix
	 *	@param[in,out] r_system is reference to system (used to query edge vertices)
	 */
	template <class CSystem>
	CEdgePose3D(int n_node0, int n_node1, Eigen::Matrix<double, 6, 1> &r_v_delta,
		Eigen::Matrix<double, 6, 6> &r_t_inv_sigma, CSystem &r_system) // this is edge constructor how it is called in the parse loop; you need to change type of the edge
		:CSEBaseEdgeImpl<CEdgePose3D, CVertexPose3D, CVertexPose3D, 6>(n_node0,
		n_node1, r_v_delta, r_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		//fprintf(stderr, "%f %f %f\n", r_t_edge.m_v_delta(0), r_t_edge.m_v_delta(1), r_t_edge.m_v_delta(2));

		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexPose3D>(n_node0, CInitializeNullVertex());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexPose3D>(n_node1,
			CRelative_to_Absolute_XYZ_Initializer(m_p_vertex0->v_State(), r_v_delta)); // rename your initializer if required

		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(((CSEBaseVertex*)m_p_vertex0)->n_Dimension() == 6);
		_ASSERTE(((CSEBaseVertex*)m_p_vertex1)->n_Dimension() == 6); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief updates the edge with a new measurement
	 *
	 *	@param[in] r_t_edge is parsed edge
	 */
	inline void Update(const CParserBase::TEdge3D &r_t_edge)
	{
		CSEBaseEdgeImpl<CEdgePose3D, CVertexPose3D, CVertexPose3D,
			6>::Update(r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma);
	}

	/**
	 *	@brief calculates jacobians, expectation and error
	 *
	 *	@param[out] r_t_jacobian0 is jacobian, associated with the first vertex
	 *	@param[out] r_t_jacobian1 is jacobian, associated with the second vertex
	 *	@param[out] r_v_expectation is expecation vector
	 *	@param[out] r_v_error is error vector
	 */
	inline void Calculate_Jacobians_Expectation_Error(Eigen::Matrix<double, 6, 6> &r_t_jacobian0,
		Eigen::Matrix<double, 6, 6> &r_t_jacobian1, Eigen::Matrix<double, 6, 1> &r_v_expectation,
		Eigen::Matrix<double, 6, 1> &r_v_error) const // change dimensionality of eigen types, if required
	{
		C3DJacobians::Absolute_to_Relative(m_p_vertex0->v_State(),
			m_p_vertex1->v_State(), r_v_expectation, r_t_jacobian0, r_t_jacobian1); // write your jacobian calculation code here (vertex state vectors are inputs, the rest are the outputs that you need to fill)
		// calculates the expectation and the jacobians

		//const Eigen::Matrix<double, 6, 1> &p1 = m_v_measurement;
		//const Eigen::Matrix<double, 6, 1> &p2 = r_v_expectation;

		//r_v_error(0) = p1(0) - p2(0);
		//r_v_error(1) = p1(1) - p2(1);
		//r_v_error(2) = p1(2) - p2(2);
		r_v_error.head<3>() = m_v_measurement.head<3>() - r_v_expectation.head<3>();

		//sum the rotations
		Eigen::Matrix3d pQ = C3DJacobians::Operator_rot(m_v_measurement.tail<3>());
		Eigen::Matrix3d dQ = C3DJacobians::Operator_rot(r_v_expectation.tail<3>());
		//Eigen::Matrix3d dQ_inv = dQ.inverse();

		Eigen::Matrix3d QQ;// = pQ * dQ_inv;
		QQ.noalias() = pQ * dQ.inverse(); // eigen likes long expressions; multiplication needs .noalias() so that eigen knows it can do it without allocating intermediate storage
		//Eigen::Vector3d axis = C3DJacobians::Operator_arot(QQ);

		//r_v_error(3) = axis(0);
		//r_v_error(4) = axis(1);
		//r_v_error(5) = axis(2);
		r_v_error.tail<3>() = C3DJacobians::Operator_arot(QQ); // avoid copying arround
		//TODO: fix the angle somehow? this cannot be right
		// calculates error (possibly re-calculates, if running A-SLAM)
	}

	/**
	 *	@brief calculates chi-square error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Chi_Squared_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		Eigen::Matrix<double, 6, 6> p_jacobi[2];
		Eigen::Matrix<double, 6, 1> v_expectation, v_error;
		Calculate_Jacobians_Expectation_Error(p_jacobi[0], p_jacobi[1], v_expectation, v_error);
		// calculates the expectation, error and the jacobians

		return (v_error.transpose() * m_t_sigma_inv).dot(v_error); // ||h_i(O_i) - z_i||^2 lambda_i
	}
};

/**
 *	@brief edge traits for SE(3) solver
 */
template <class CParsedStructure>
class CSE3OnlyPoseEdgeTraits {
public:
	typedef CFailOnEdgeType _TyEdge; /**< @brief it should fail on unknown edge types */

	/**
	 *	@brief gets reason for error
	 *	@return Returns const null-terminated string, containing
	 *		description of the error (human readable).
	 */
	static const char *p_s_Reason()
	{
		return "unknown edge type occured";
	}
};

/**
 *	@brief edge traits for SE(3) solver (specialized for CParser::TEdge3D)
 */
template <>
class CSE3OnlyPoseEdgeTraits<CParserBase::TEdge3D> {
public:
	typedef CEdgePose3D _TyEdge; /**< @brief the edge type to construct from the parsed type */
};


/**
 *	@brief edge traits for SE(3) solver (specialized for CParser::TEdge3D)
 */
template <>
class CSE3OnlyPoseEdgeTraits<CParserBase::TVertex3D> {
public:
	typedef CIgnoreEdgeType _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

#endif // __SE3_TYPES_INCLUDED
