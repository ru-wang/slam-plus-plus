/**
 *	@file include/slam/BA_Types.h
 *	@author Soso
 *	@date 2013
 *	@brief projection type for BA
 */

#pragma once
#ifndef __BA_TYPES_INCLUDED
#define __BA_TYPES_INCLUDED

#include "slam/BaseTypes.h"
#include "slam/SE2_Types.h"
#include "slam/BASolverBase.h"
#include <iostream> // SOSO

/*class CVertexIntrinsics : public CSEBaseVertexImpl<CVertexIntrinsics, 4> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

	/ **
	 *	@brief default constructor; has no effect
	 * /
    inline CVertexPoseBA() // copy this
    {}

	/ **
	 *	@brief constructor from vector type
	 *	@param[in] r_v_state is state of the vertex
	 * /
    inline CVertexIntrinsics(const Eigen::VectorXd &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexIntrinsics, 4>(r_v_state) // change the dimension here as well
    {}

	/ **
	 *	@copydoc CSEBaseVertex::Operator_Plus()
	 * /
    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
    	m_v_state.segment<4>() += r_v_delta.segment<4>(m_n_order);
    }
};*/

/**
 *	@brief bundle adjustment camera pose
 */
class CVertexCam : public CSEBaseVertexImpl<CVertexCam, 6> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
protected:
    Eigen::Matrix<double, 5, 1> m_v_intrinsics; /**< @brief vertex cam should hold own camera params (these are not being optimized) */

public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

	/**
	 *	@brief default constructor; has no effect
	 */
    inline CVertexCam() // copy this
    {}

	/**
	 *	@brief constructor from vector type
	 *	@param[in] r_v_state is state of the vertex
	 */
    inline CVertexCam(const Eigen::VectorXd &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexCam, 6>(r_v_state.head<6>()), // change the dimension here as well
		m_v_intrinsics(r_v_state.tail<5>())
    {}

	/**
	 *	@brief constructor from parsed type
	 *	@param[in] r_t_vertex is state of the vertex
	 */
    inline CVertexCam(const CParserBase::TVertexCam3D &r_t_vertex) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexCam, 6>(r_t_vertex.m_v_position.head<6>()), // change the dimension here as well
		m_v_intrinsics(r_t_vertex.m_v_position.tail<5>())
    {}

	/**
	 *	@brief gets intrinsic parameters of the camera
	 *	@return Returns const reference to the parameters of the camera.
	 */
	inline const Eigen::Matrix<double, 5, 1> &v_Intrinsics() const
	{
		return m_v_intrinsics;
	}

	/**
	 *	@copydoc CSEBaseVertex::Operator_Plus()
	 */
    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
        //focal lengths and principal points are added normally
    	//SOSO: ignore adding of fx and cx
    	//m_v_state.segment<4>(6) += r_v_delta.segment<4>(m_n_order + 6);

    	//pose is added as SE3
    	m_v_state.head<3>() += r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated

		//sum the rotations
		Eigen::Matrix3d pQ = C3DJacobians::Operator_rot(m_v_state.segment<3>(3));
		Eigen::Matrix3d dQ = C3DJacobians::Operator_rot(r_v_delta.segment<3>(m_n_order + 3));

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // noalias says that the destination is not the same as the multiplicands and that the multiplication can be carried out without allocating intermediate storage
		m_v_state.segment<3>(3) = CBAJacobians::Operator_arot(QQ); // also, no need for the intermediate object
    }

	/**
	 *	@copydoc CSEBaseVertex::Operator_Minus()
	 */
    inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
	{
		Operator_Plus(-r_v_delta); // ...
	}
};

/**
 *	@brief bundle adjustment (observed) keypoint vertex
 */
class CVertexXYZ : public CSEBaseVertexImpl<CVertexXYZ, 3> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

	/**
	 *	@brief default constructor; has no effect
	 */
    inline CVertexXYZ() // copy this
    {}

	/**
	 *	@brief constructor from vector type
	 *	@param[in] r_v_state is state of the vertex
	 */
    inline CVertexXYZ(const Eigen::Vector3d &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexXYZ, 3>(r_v_state) // change the dimension here as well
    {}

	/**
	 *	@brief constructor from parsed type
	 *	@param[in] r_t_vertex is state of the vertex
	 */
    inline CVertexXYZ(const CParserBase::TVertexXYZ &r_t_vertex) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexXYZ, 3>(r_t_vertex.m_v_position) // change the dimension here as well
    {}

	/**
	 *	@copydoc CSEBaseVertex::Operator_Plus()
	 */
    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
    	//pose is added as SE3
    	m_v_state += r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated
    }

	/**
	 *	@copydoc CSEBaseVertex::Operator_Minus()
	 */
    inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" plus
	{
		//pose is added as SE3
		m_v_state -= r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated
	}
};

/**
 *	@brief BA vex-cam edge
 */
class CEdgeP2C3D : public CSEBaseEdgeImpl<CEdgeP2C3D, CVertexCam, CVertexXYZ, 2> { // again, this tells that base implementation for base edge for type that will be called CEdgePose2D, and it will be an edge between two vertices of type CVertexPose2D, and the measurement will have 3 dimensions (in order to have hyperedges, you need to provide your own version of CSEBaseEdgeImpl, that is an advanced topic)
public:
	__SE3_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, just copy this

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CEdgeP2C3D() // copy this
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
	CEdgeP2C3D(const CParserBase::TEdgeP2C3D &r_t_edge, CSystem &r_system) // this is edge constructor how it is called in the parse loop; you need to change type of the edge
		:CSEBaseEdgeImpl<CEdgeP2C3D, CVertexCam, CVertexXYZ, 2>(r_t_edge.m_n_node_1,
		r_t_edge.m_n_node_0, r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexCam>(r_t_edge.m_n_node_1, CInitializeNullVertex());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexXYZ>(r_t_edge.m_n_node_0, CInitializeNullVertex1()); // rename your initializer if required

		//std::swap(r_t_edge.m_n_node_1, r_t_edge.m_n_node_0);

		//fprintf(stderr, "verts: %d %d \n", r_t_edge.m_n_node_0, r_t_edge.m_n_node_1);
		//fprintf(stderr, ">>%f %f %f %f %f %f %f %f %f %f \n", m_p_vertex0->r_v_State()(0), m_p_vertex0->r_v_State()(1), m_p_vertex0->r_v_State()(2), m_p_vertex0->r_v_State()(3),
		//		m_p_vertex0->r_v_State()(4), m_p_vertex0->r_v_State()(5), m_p_vertex0->r_v_State()(6), m_p_vertex0->r_v_State()(7), m_p_vertex0->r_v_State()(8), m_p_vertex0->r_v_State()(9));

		//fprintf(stderr, "dims: %d %d \n", ((CSEBaseVertex*)m_p_vertex0)->n_Dimension(), ((CSEBaseVertex*)m_p_vertex1)->n_Dimension());

		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(((CSEBaseVertex*)m_p_vertex0)->n_Dimension() == 6);
		_ASSERTE(((CSEBaseVertex*)m_p_vertex1)->n_Dimension() == 3); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief updates the edge with a new measurement
	 *
	 *	@param[in] r_t_edge is parsed edge
	 */
	inline void Update(const CParserBase::TEdgeP2C3D &r_t_edge)
	{
		CSEBaseEdgeImpl<CEdgeP2C3D, CVertexCam, CVertexXYZ,
			2>::Update(r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma);
	}

	/**
	 *	@brief calculates jacobians, expectation and error
	 *
	 *	@param[out] r_t_jacobian0 is jacobian, associated with the first vertex
	 *	@param[out] r_t_jacobian1 is jacobian, associated with the second vertex
	 *	@param[out] r_v_expectation is expecation vector
	 *	@param[out] r_v_error is error vector
	 */
	inline void Calculate_Jacobians_Expectation_Error(Eigen::Matrix<double, 2, 6> &r_t_jacobian0,
		Eigen::Matrix<double, 2, 3> &r_t_jacobian1, Eigen::Matrix<double, 2, 1> &r_v_expectation,
		Eigen::Matrix<double, 2, 1> &r_v_error) const // change dimensionality of eigen types, if required
	{
		CBAJacobians::Project_P2C(m_p_vertex0->v_State() /* Cam */, m_p_vertex0->v_Intrinsics() /* camera intrinsic params */,
			m_p_vertex1->v_State() /* XYZ */, r_v_expectation, r_t_jacobian0, r_t_jacobian1); // write your jacobian calculation code here (vertex state vectors are inputs, the rest are the outputs that you need to fill)
		// calculates the expectation and the jacobians

		//std::cout << "me: " << m_v_measurement[0] << " " << m_v_measurement[1] << " ----  " << r_v_expectation[0] << " " << r_v_expectation[1] << std::endl;
		//std::cout << "J0:" << r_t_jacobian0 << std::endl;
		//std::cout << "J1:" << r_t_jacobian1 << std::endl;

		r_v_error = m_v_measurement - r_v_expectation;
		// calculates error (possibly re-calculates, if running A-SLAM)*/
	}

	/**
	 *	@brief calculates chi-square error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Chi_Squared_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		Eigen::Matrix<double, 2, 6> p_jacobi0;
		Eigen::Matrix<double, 2, 3> p_jacobi1;
		Eigen::Matrix<double, 2, 1> v_expectation, v_error;
		Calculate_Jacobians_Expectation_Error(p_jacobi0, p_jacobi1, v_expectation, v_error);
		// calculates the expectation, error and the jacobians
		//std::cerr << v_error[0] << " " << v_error[1] << " ----  " << (v_error.transpose() * m_t_sigma_inv).dot(v_error) << std::endl;
		//std::cout << sqrt( v_error[0]*v_error[0] + v_error[1]*v_error[1] ) << std::endl;

		return (v_error.transpose() * m_t_sigma_inv).dot(v_error); // ||h_i(O_i) - z_i||^2 lambda_i*/
	}

	/**
	 *	@brief calculates reprojection error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Reprojection_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		Eigen::Matrix<double, 2, 6> p_jacobi0;
		Eigen::Matrix<double, 2, 3> p_jacobi1;
		Eigen::Matrix<double, 2, 1> v_expectation, v_error;
		Calculate_Jacobians_Expectation_Error(p_jacobi0, p_jacobi1, v_expectation, v_error);
		// calculates the expectation, error and the jacobians
		//std::cerr << v_error[0] << " " << v_error[1] << " ----  " << (v_error.transpose() * m_t_sigma_inv).dot(v_error) << std::endl;

		return sqrt( (v_error[0] - v_expectation[0])*(v_error[0] - v_expectation[0]) +
				(v_error[1] - v_expectation[1])*(v_error[1] - v_expectation[1]) );
	}
};

/**
 *	@brief edge traits for BA solver
 */
template <class CParsedStructure>
class CBAEdgeTraits {
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
 *	@brief edge traits for BA solver
 */
template <>
class CBAEdgeTraits<CParserBase::TEdgeP2C3D> {
public:
	typedef CEdgeP2C3D _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBAEdgeTraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBAEdgeTraits<CParserBase::TVertexCam3D> {
public:
	typedef CVertexCam _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for BA solver
 */
template <class CParsedStructure>
class CBAVertexTraits {
public:
	typedef CFailOnVertexType _TyVertex; /**< @brief it should fail on unknown vertex types */

	/**
	 *	@brief gets reason for error
	 *	@return Returns const null-terminated string, containing
	 *		description of the error (human readable).
	 */
	static const char *p_s_Reason()
	{
		return "unknown vertex type occured";
	}
};

/**
 *	@brief vertex traits for BA solver
 */
template <>
class CBAVertexTraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for BA solver
 */
template <>
class CBAVertexTraits<CParserBase::TVertexCam3D> {
public:
	typedef CVertexCam _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

#endif
