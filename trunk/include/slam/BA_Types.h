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
//#include "slam/SE3_Types.h" // not this
#include "slam/BASolverBase.h"
#include "slam/3DSolverBase.h" // that. want rot / arot

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
	 *	@copydoc base_iface::CVertexFacade::Operator_Plus()
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
    Eigen::Matrix<double, 5, 1, Eigen::DontAlign> m_v_intrinsics; /**< @brief vertex cam should hold own camera params (these are not being optimized) */

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
    inline CVertexCam(const Eigen::Matrix<double, 11, 1> &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexCam, 6>(r_v_state.head<6>()), // change the dimension here as well
		m_v_intrinsics(r_v_state.tail<5>())
    {}

	/**
	 *	@brief constructor from optimized vector type (because of CSEBaseVertexImpl<CVertexCam, 6>)
	 *	@param[in] r_v_state is state of the vertex
	 */
    inline CVertexCam(const Eigen::Matrix<double, 6, 1> &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexCam, 6>(r_v_state.head<6>()) // change the dimension here as well
    {
		m_v_intrinsics.setZero(); // we did not get this information
	}

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
	inline const Eigen::Matrix<double, 5, 1, Eigen::DontAlign> &v_Intrinsics() const
	{
		return m_v_intrinsics;
	}

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Plus()
	 */
    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
		C3DJacobians::Relative_to_Absolute(m_v_state, r_v_delta.segment<6>(m_n_order), m_v_state);
		//CBAJacobians::Smart_Plus_Cam(m_v_state, r_v_delta.segment<6>(m_n_order), m_v_state);

        /*//focal lengths and principal points are added normally
    	//SOSO: ignore adding of fx and cx
    	//m_v_state.segment<4>(6) += r_v_delta.segment<4>(m_n_order + 6);

    	//pose is added as SE3
    	m_v_state.head<3>() += r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated

		//sum the rotations
		Eigen::Matrix3d pQ = C3DJacobians::t_AxisAngle_to_RotMatrix(m_v_state.tail<3>());
		Eigen::Matrix3d dQ = C3DJacobians::t_AxisAngle_to_RotMatrix(r_v_delta.segment<3>(m_n_order + 3));

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // noalias says that the destination is not the same as the multiplicands and that the multiplication can be carried out without allocating intermediate storage
		m_v_state.tail<3>() = C3DJacobians::v_RotMatrix_to_AxisAngle(QQ); // also, no need for the intermediate object*/
    }

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Minus()
	 */
    inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
	{
		//Operator_Plus(-r_v_delta); // ...
		C3DJacobians::Relative_to_Absolute(m_v_state, -r_v_delta.segment<6>(m_n_order), m_v_state); // avoid calculating negative of the whole r_v_delta
	}
};

/**
 *	@brief bundle adjustment camera pose
 */
class CVertexIntrinsics : public CSEBaseVertexImpl<CVertexIntrinsics, 5> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions

public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

	/**
	 *	@brief default constructor; has no effect
	 */
    inline CVertexIntrinsics() // copy this
    {}

	/**
	 *	@brief constructor from vector type
	 *	@param[in] r_v_state is state of the vertex
	 */
    inline CVertexIntrinsics(const Eigen::Matrix<double, 5, 1> &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexIntrinsics, 5>(r_v_state)
    {}

	/**
	 *	@brief constructor from parsed type
	 *	@param[in] r_t_vertex is state of the vertex
	 */
    inline CVertexIntrinsics(const CParserBase::TVertexIntrinsics &r_t_vertex) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexIntrinsics, 5>(r_t_vertex.m_v_position)
    {}

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Plus()
	 */
    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
    	m_v_state = m_v_state + r_v_delta.segment<5>(m_n_order);
    }

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Minus()
	 */
    inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
	{
    	m_v_state = m_v_state - r_v_delta.segment<5>(m_n_order);
	}
};

/**
 *	@brief bundle adjustment camera pose
 */
class CVertexSCam : public CSEBaseVertexImpl<CVertexSCam, 6> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
protected:
    Eigen::Matrix<double, 6, 1, Eigen::DontAlign> m_v_intrinsics; /**< @brief vertex cam should hold own camera params (these are not being optimized) */

public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

	/**
	 *	@brief default constructor; has no effect
	 */
    inline CVertexSCam() // copy this
    {}

	/**
	 *	@brief constructor from vector type
	 *	@param[in] r_v_state is state of the vertex
	 */
    inline CVertexSCam(const Eigen::Matrix<double, 12, 1> &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexSCam, 6>(r_v_state.head<6>()), // change the dimension here as well
		m_v_intrinsics(r_v_state.tail<6>())
    {}

	/**
	 *	@brief constructor from optimized vector type (because of CSEBaseVertexImpl<CVertexCam, 6>)
	 *	@param[in] r_v_state is state of the vertex
	 */
    inline CVertexSCam(const Eigen::Matrix<double, 6, 1> &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexSCam, 6>(r_v_state.head<6>()) // change the dimension here as well
    {
		m_v_intrinsics.setZero(); // we did not get this information
	}

	/**
	 *	@brief constructor from parsed type
	 *	@param[in] r_t_vertex is state of the vertex
	 */
    inline CVertexSCam(const CParserBase::TVertexSCam3D &r_t_vertex) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexSCam, 6>(r_t_vertex.m_v_position.head<6>()), // change the dimension here as well
		m_v_intrinsics(r_t_vertex.m_v_position.tail<6>())
    {}

	/**
	 *	@brief gets intrinsic parameters of the camera
	 *	@return Returns const reference to the parameters of the camera.
	 */
	inline const Eigen::Matrix<double, 6, 1, Eigen::DontAlign> &v_Intrinsics() const
	{
		return m_v_intrinsics;
	}

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Plus()
	 */
    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
		C3DJacobians::Relative_to_Absolute(m_v_state, r_v_delta.segment<6>(m_n_order), m_v_state);

        /*//focal lengths and principal points are added normally
    	//SOSO: ignore adding of fx and cx
    	//m_v_state.segment<4>(6) += r_v_delta.segment<4>(m_n_order + 6);

    	//pose is added as SE3
    	m_v_state.head<3>() += r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated

		//sum the rotations
		Eigen::Matrix3d pQ = C3DJacobians::t_AxisAngle_to_RotMatrix(m_v_state.tail<3>());
		Eigen::Matrix3d dQ = C3DJacobians::t_AxisAngle_to_RotMatrix(r_v_delta.segment<3>(m_n_order + 3));

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // noalias says that the destination is not the same as the multiplicands and that the multiplication can be carried out without allocating intermediate storage
		m_v_state.tail<3>() = C3DJacobians::v_RotMatrix_to_AxisAngle(QQ); // also, no need for the intermediate object*/
    }

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Minus()
	 */
    inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
	{
		//Operator_Plus(-r_v_delta); // ...
		C3DJacobians::Relative_to_Absolute(m_v_state, -r_v_delta.segment<6>(m_n_order), m_v_state);
	}
};

/**
 *	@brief bundle adjustment camera pose
 */
class CVertexSpheron : public CSEBaseVertexImpl<CVertexSpheron, 6> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

	/**
	 *	@brief default constructor; has no effect
	 */
    inline CVertexSpheron() // copy this
    {}

	/**
	 *	@brief constructor from vector type
	 *	@param[in] r_v_state is state of the vertex
	 */
    inline CVertexSpheron(const Eigen::Matrix<double, 6, 1> &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexSpheron, 6>(r_v_state)
    {}

	/**
	 *	@brief constructor from parsed type
	 *	@param[in] r_t_vertex is state of the vertex
	 */
    inline CVertexSpheron(const CParserBase::TVertexSpheron &r_t_vertex) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexSpheron, 6>(r_t_vertex.m_v_position)
    {}

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Plus()
	 */
    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
		C3DJacobians::Relative_to_Absolute(m_v_state, r_v_delta.segment<6>(m_n_order), m_v_state);

       	/*//pose is added as SE3
    	m_v_state.head<3>() += r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated

		//sum the rotations
		Eigen::Matrix3d pQ = C3DJacobians::t_AxisAngle_to_RotMatrix(m_v_state.tail<3>());
		Eigen::Matrix3d dQ = C3DJacobians::t_AxisAngle_to_RotMatrix(r_v_delta.segment<3>(m_n_order + 3));

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // noalias says that the destination is not the same as the multiplicands and that the multiplication can be carried out without allocating intermediate storage
		m_v_state.tail<3>() = C3DJacobians::v_RotMatrix_to_AxisAngle(QQ); // also, no need for the intermediate object*/
    }

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Minus()
	 */
    inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
	{
		//Operator_Plus(-r_v_delta); // ...
		C3DJacobians::Relative_to_Absolute(m_v_state, -r_v_delta.segment<6>(m_n_order), m_v_state);
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
	 *	@copydoc base_iface::CVertexFacade::Operator_Plus()
	 */
    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
    	//pose is added as SE3
    	m_v_state += r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated
    }

	/**
	 *	@copydoc base_iface::CVertexFacade::Operator_Minus()
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
class CEdgeP2C3D : public CBaseEdgeImpl<CEdgeP2C3D, MakeTypelist(CVertexCam, CVertexXYZ), 2> { // again, this tells that base implementation for base edge for type that will be called CEdgePose2D, and it will be an edge between two vertices of type CVertexPose2D, and the measurement will have 3 dimensions (in order to have hyperedges, you need to provide your own version of CBaseEdgeImpl, that is an advanced topic)
public:
	__BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, just copy this

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
		:CBaseEdgeImpl<CEdgeP2C3D, MakeTypelist(CVertexCam, CVertexXYZ), 2>(r_t_edge.m_n_node_1,
		r_t_edge.m_n_node_0, r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexCam>(r_t_edge.m_n_node_1, CInitializeNullVertex<>());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexXYZ>(r_t_edge.m_n_node_0, CInitializeNullVertex<CVertexXYZ>()); // rename your initializer if required

		//std::swap(r_t_edge.m_n_node_1, r_t_edge.m_n_node_0);

		//fprintf(stderr, "verts: %d %d \n", r_t_edge.m_n_node_0, r_t_edge.m_n_node_1);
		//fprintf(stderr, ">>%f %f %f %f %f %f %f %f %f %f \n", m_p_vertex0->r_v_State()(0), m_p_vertex0->r_v_State()(1), m_p_vertex0->r_v_State()(2), m_p_vertex0->r_v_State()(3),
		//		m_p_vertex0->r_v_State()(4), m_p_vertex0->r_v_State()(5), m_p_vertex0->r_v_State()(6), m_p_vertex0->r_v_State()(7), m_p_vertex0->r_v_State()(8), m_p_vertex0->r_v_State()(9));

		//fprintf(stderr, "dims: %d %d \n", (m_p_vertex0)->n_Dimension(), (m_p_vertex1)->n_Dimension());

		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_1].n_Dimension() == 6); // get the vertices from the vertex pool to ensure a correct type is used, do not use m_p_vertex0 / m_p_vertex1 for this
		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_0].n_Dimension() == 3); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief constructor; initializes edge with data
	 *
	 *	@tparam CSystem is type of system where this edge is being stored
	 *
	 *	@param[in] n_node0 is (zero-based) index of the first (origin) node (camera)
	 *	@param[in] n_node1 is (zero-based) index of the second (endpoint) node (vertex)
	 *	@param[in] v_delta is vector of delta x position, delta y-position
	 *	@param[in] r_t_inv_sigma is the information matrix
	 *	@param[in,out] r_system is reference to system (used to query edge vertices)
	 */
	template <class CSystem>
	CEdgeP2C3D(size_t n_node1, size_t n_node0, const Eigen::Vector2d &v_delta, // won't align non-reference arg
		const Eigen::Matrix2d &r_t_inv_sigma, CSystem &r_system) // respect const-ness!
		:CBaseEdgeImpl<CEdgeP2C3D, MakeTypelist(CVertexCam, CVertexXYZ), 2>(n_node0,
		n_node1, v_delta, r_t_inv_sigma)
	{
		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexCam>(n_node0, CInitializeNullVertex<>());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexXYZ>(n_node1, CInitializeNullVertex<CVertexXYZ>());
		// get vertices (initialize if required)
		// "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(r_system.r_Vertex_Pool()[n_node0].n_Dimension() == 6); // get the vertices from the vertex pool to ensure a correct type is used, do not use m_p_vertex0 / m_p_vertex1 for this
		_ASSERTE(r_system.r_Vertex_Pool()[n_node1].n_Dimension() == 3);
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief updates the edge with a new measurement
	 *
	 *	@param[in] r_t_edge is parsed edge
	 */
	inline void Update(const CParserBase::TEdgeP2C3D &r_t_edge)
	{
		CBaseEdgeImpl<CEdgeP2C3D, MakeTypelist(CVertexCam, CVertexXYZ),
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
		CBAJacobians::Project_P2C(m_p_vertex0->r_v_State() /* Cam */, m_p_vertex0->v_Intrinsics() /* camera intrinsic params */,
			m_p_vertex1->r_v_State() /* XYZ */, r_v_expectation, r_t_jacobian0, r_t_jacobian1); // write your jacobian calculation code here (vertex state vectors are inputs, the rest are the outputs that you need to fill)
		// calculates the expectation and the jacobians

		r_v_error = m_v_measurement - r_v_expectation;
		// calculates error (possibly re-calculates, if running A-SLAM)*/
	}

	/**
	 *	@brief calculates chi-square error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Chi_Squared_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		/*Eigen::Matrix<double, 2, 6> p_jacobi0;
		Eigen::Matrix<double, 2, 3> p_jacobi1;
		Eigen::Matrix<double, 2, 1> v_expectation, v_error;
		Calculate_Jacobians_Expectation_Error(p_jacobi0, p_jacobi1, v_expectation, v_error);
		// calculates the expectation, error and the jacobians
		//std::cerr << v_error[0] << " " << v_error[1] << " ----  " << (v_error.transpose() * m_t_sigma_inv).dot(v_error) << std::endl;
		//std::cout << sqrt( v_error[0]*v_error[0] + v_error[1]*v_error[1] ) << std::endl;*/

		Eigen::Vector2d v_error;
		CBAJacobians::Project_P2C(m_p_vertex0->r_v_State() /* Cam */, m_p_vertex0->v_Intrinsics(),
			m_p_vertex1->r_v_State() /* XYZ */, v_error);
		//std::cout << v_error.transpose() << " | " << m_v_measurement.transpose() << std::endl;

		v_error -= m_v_measurement; // this is actually negative error, but it is squared below so the chi2 is the same
		// this should be faster, as it does not calculate the jacobians
		// (in BA, this is actually timed as it is used in solver step validation)

		return (v_error.transpose() * m_t_sigma_inv).dot(v_error); // ||z_i - h_i(O_i)||^2 lambda_i*/
	}

	/**
	 *	@brief calculates reprojection error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Reprojection_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		/*Eigen::Matrix<double, 2, 6> p_jacobi0;
		Eigen::Matrix<double, 2, 3> p_jacobi1;
		Eigen::Matrix<double, 2, 1> v_expectation, v_error;
		Calculate_Jacobians_Expectation_Error(p_jacobi0, p_jacobi1, v_expectation, v_error);*/
		// calculates the expectation, error and the jacobians
		//std::cerr << v_error[0] << " " << v_error[1] << " ----  " << (v_error.transpose() * m_t_sigma_inv).dot(v_error) << std::endl;

		Eigen::Vector2d v_error;
		CBAJacobians::Project_P2C(m_p_vertex0->r_v_State() /* Cam */, m_p_vertex0->v_Intrinsics(),
			m_p_vertex1->r_v_State() /* XYZ */, v_error);
		v_error -= m_v_measurement; // this is actually negative error, but we only need the norm so sign does not matter
		// this should be faster, as it does not calculate the jacobians
		// (in BA, this is actually timed as it is used in solver step validation)

		return v_error.norm();
		/*return sqrt( (v_error[0] - v_expectation[0])*(v_error[0] - v_expectation[0]) +
				(v_error[1] - v_expectation[1])*(v_error[1] - v_expectation[1]) );*/
	}
};

/**
 *	@brief BA vex-cam edge
 */
class CEdgeP2CI3D : public CBaseEdgeImpl<CEdgeP2CI3D, MakeTypelist(CVertexCam, CVertexXYZ, CVertexIntrinsics), 2> { // again, this tells that base implementation for base edge for type that will be called CEdgePose2D, and it will be an edge between two vertices of type CVertexPose2D, and the measurement will have 3 dimensions (in order to have hyperedges, you need to provide your own version of CBaseEdgeImpl, that is an advanced topic)
public:
	typedef CBaseEdgeImpl<CEdgeP2CI3D, MakeTypelist(CVertexCam, CVertexXYZ, CVertexIntrinsics), 2> _TyBase; /**< @brief base edge type */

public:
	__BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, just copy this

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CEdgeP2CI3D() // copy this
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
	CEdgeP2CI3D(const CParserBase::TEdgeP2CI3D &r_t_edge, CSystem &r_system) // this is edge constructor how it is called in the parse loop; you need to change type of the edge
		:_TyBase(typename _TyBase::_TyVertexIndexTuple(r_t_edge.m_n_node_1,
				r_t_edge.m_n_node_0, r_t_edge.m_n_node_2), r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		m_vertex_ptr.Get<0>() = &r_system.template r_Get_Vertex<CVertexCam>(r_t_edge.m_n_node_1, CInitializeNullVertex<>());
		m_vertex_ptr.Get<1>() = &r_system.template r_Get_Vertex<CVertexXYZ>(r_t_edge.m_n_node_0, CInitializeNullVertex<CVertexXYZ>());
		m_vertex_ptr.Get<2>() = &r_system.template r_Get_Vertex<CVertexIntrinsics>(r_t_edge.m_n_node_2, CInitializeNullVertex<CVertexIntrinsics>()); // rename your initializer if required
		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_0].n_Dimension() == 3); // get the vertices from the vertex pool to ensure a correct type is used, do not use m_p_vertex0 / m_p_vertex1 for this
		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_1].n_Dimension() == 6);
		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_2].n_Dimension() == 5); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief constructor; initializes edge with data
	 *
	 *	@tparam CSystem is type of system where this edge is being stored
	 *
	 *	@param[in] n_node_0 is (zero-based) index of the first (origin) node (camera)
	 *	@param[in] n_node_1 is (zero-based) index of the second (endpoint) node (xyz vertex)
	 *	@param[in] n_node_2 is (zero-based) index of the third (endpoint) node (intrinsics vertex)
	 *	@param[in] v_delta is vector of delta x position, delta y-position
	 *	@param[in] r_t_inv_sigma is the information matrix
	 *	@param[in,out] r_system is reference to system (used to query edge vertices)
	 */
	template <class CSystem>
	CEdgeP2CI3D(size_t n_node_0, size_t n_node_1, size_t n_node_2, const Eigen::Vector2d &v_delta, // won't align non-reference arg
		const Eigen::Matrix2d &r_t_inv_sigma, CSystem &r_system) // respect const-ness!
		:_TyBase(typename _TyBase::_TyVertexIndexTuple(n_node_0,
				n_node_1, n_node_2), v_delta, r_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		m_vertex_ptr.Get<0>() = &r_system.template r_Get_Vertex<CVertexCam>(n_node_0, CInitializeNullVertex<>());
		m_vertex_ptr.Get<1>() = &r_system.template r_Get_Vertex<CVertexXYZ>(n_node_1, CInitializeNullVertex<CVertexXYZ>());
		m_vertex_ptr.Get<2>() = &r_system.template r_Get_Vertex<CVertexIntrinsics>(n_node_2, CInitializeNullVertex<CVertexIntrinsics>()); // rename your initializer if required
		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(r_system.r_Vertex_Pool()[n_node_1].n_Dimension() == 3); // get the vertices from the vertex pool to ensure a correct type is used, do not use m_p_vertex0 / m_p_vertex1 for this
		_ASSERTE(r_system.r_Vertex_Pool()[n_node_0].n_Dimension() == 6);
		_ASSERTE(r_system.r_Vertex_Pool()[n_node_2].n_Dimension() == 5); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief updates the edge with a new measurement
	 *
	 *	@param[in] r_t_edge is parsed edge
	 */
	inline void Update(const CParserBase::TEdgeP2CI3D &r_t_edge)
	{
		CBaseEdgeImpl<CEdgeP2CI3D, MakeTypelist(CVertexCam, CVertexXYZ, CVertexIntrinsics),
			2>::Update(r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma);
	}

	/**
	 *	@brief calculates jacobians, expectation and error
	 *
	 *	@param[out] r_t_jacobian_tuple structure holding all jacobians
	 *	@param[out] r_v_expectation is expecation vector
	 *	@param[out] r_v_error is error vector
	 */
	inline void Calculate_Jacobians_Expectation_Error(_TyBase::_TyJacobianTuple &r_t_jacobian_tuple,
		Eigen::Matrix<double, 2, 1> &r_v_expectation, Eigen::Matrix<double, 2, 1> &r_v_error) const // change dimensionality of eigen types, if required
	{
		CBAJacobians::Project_P2CI(m_vertex_ptr.Get<0>()->v_State(), m_vertex_ptr.Get<2>()->v_State(),
			m_vertex_ptr.Get<1>()->v_State(), r_v_expectation,
			r_t_jacobian_tuple.Get<0>(), r_t_jacobian_tuple.Get<1>(),
			r_t_jacobian_tuple.Get<2>()); // write your jacobian calculation code here (vertex state vectors are inputs, the rest are the outputs that you need to fill)
		// calculates the expectation and the jacobians

		//std::cout << r_v_expectation.transpose() << " | " << m_v_measurement.transpose() << std::endl;

		r_v_error = m_v_measurement - r_v_expectation;
		// calculates error (possibly re-calculates, if running A-SLAM)*/
	}

	/**
	 *	@brief calculates chi-square error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Chi_Squared_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		Eigen::Vector2d v_error;
		CBAJacobians::Project_P2C(m_vertex_ptr.Get<0>()->v_State() /* Cam */, m_vertex_ptr.Get<2>()->v_State(),
			m_vertex_ptr.Get<1>()->v_State() /* XYZ */, v_error);

		v_error -= m_v_measurement; // this is actually negative error, but it is squared below so the chi2 is the same
		// this should be faster, as it does not calculate the jacobians
		// (in BA, this is actually timed as it is used in solver step validation)

		return (v_error.transpose() * m_t_sigma_inv).dot(v_error); // ||z_i - h_i(O_i)||^2 lambda_i*/
	}

	/**
	 *	@brief calculates reprojection error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Reprojection_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		Eigen::Vector2d v_expectation;
		CBAJacobians::Project_P2C(m_vertex_ptr.Get<0>()->v_State() /* Cam */, m_vertex_ptr.Get<2>()->v_State(),
			m_vertex_ptr.Get<1>()->v_State() /* XYZ */, v_expectation);

		Eigen::Vector2d v_error = m_v_measurement - v_expectation;
		// this should be faster, as it does not calculate the jacobians
		// (in BA, this is actually timed as it is used in solver step validation)

		return v_error.norm();
	}
};

/**
 *	@brief BA vex-cam edge
 */
class CEdgeP2SC3D : public CBaseEdgeImpl<CEdgeP2SC3D, MakeTypelist(CVertexSCam, CVertexXYZ), 3> { // again, this tells that base implementation for base edge for type that will be called CEdgePose2D, and it will be an edge between two vertices of type CVertexPose2D, and the measurement will have 3 dimensions (in order to have hyperedges, you need to provide your own version of CBaseEdgeImpl, that is an advanced topic)
public:
	__BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, just copy this

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CEdgeP2SC3D() // copy this
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
	CEdgeP2SC3D(const CParserBase::TEdgeP2SC3D &r_t_edge, CSystem &r_system) // this is edge constructor how it is called in the parse loop; you need to change type of the edge
		:CBaseEdgeImpl<CEdgeP2SC3D, MakeTypelist(CVertexSCam, CVertexXYZ), 3>(r_t_edge.m_n_node_1,
		r_t_edge.m_n_node_0, r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexSCam>(r_t_edge.m_n_node_1, CInitializeNullVertex<>());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexXYZ>(r_t_edge.m_n_node_0, CInitializeNullVertex<CVertexXYZ>()); // rename your initializer if required

		//std::swap(r_t_edge.m_n_node_1, r_t_edge.m_n_node_0);

		//fprintf(stderr, "verts: %d %d \n", r_t_edge.m_n_node_0, r_t_edge.m_n_node_1);
		//fprintf(stderr, ">>%f %f %f %f %f %f %f %f %f %f \n", m_p_vertex0->r_v_State()(0), m_p_vertex0->r_v_State()(1), m_p_vertex0->r_v_State()(2), m_p_vertex0->r_v_State()(3),
		//		m_p_vertex0->r_v_State()(4), m_p_vertex0->r_v_State()(5), m_p_vertex0->r_v_State()(6), m_p_vertex0->r_v_State()(7), m_p_vertex0->r_v_State()(8), m_p_vertex0->r_v_State()(9));

		//fprintf(stderr, "dims: %d %d \n", (m_p_vertex0)->n_Dimension(), (m_p_vertex1)->n_Dimension());

		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_1].n_Dimension() == 6); // get the vertices from the vertex pool to ensure a correct type is used, do not use m_p_vertex0 / m_p_vertex1 for this
		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_0].n_Dimension() == 3); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief constructor; initializes edge with data
	 *
	 *	@tparam CSystem is type of system where this edge is being stored
	 *
	 *	@param[in] n_node0 is (zero-based) index of the first (origin) node (camera)
	 *	@param[in] n_node1 is (zero-based) index of the second (endpoint) node (vertex)
	 *	@param[in] v_delta is vector of delta x position, delta y-position
	 *	@param[in] r_t_inv_sigma is the information matrix
	 *	@param[in,out] r_system is reference to system (used to query edge vertices)
	 */
	template <class CSystem>
	CEdgeP2SC3D(size_t n_node1, size_t n_node0, const Eigen::Vector3d &v_delta, // won't align non-reference arg
		const Eigen::Matrix3d &r_t_inv_sigma, CSystem &r_system) // respect const-ness!
		:CBaseEdgeImpl<CEdgeP2SC3D, MakeTypelist(CVertexSCam, CVertexXYZ), 3>(n_node0,
		n_node1, v_delta, r_t_inv_sigma)
	{
		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexSCam>(n_node0, CInitializeNullVertex<>());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexXYZ>(n_node1, CInitializeNullVertex<CVertexXYZ>());
		// get vertices (initialize if required)
		// "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(r_system.r_Vertex_Pool()[n_node0].n_Dimension() == 6); // get the vertices from the vertex pool to ensure a correct type is used, do not use m_p_vertex0 / m_p_vertex1 for this
		_ASSERTE(r_system.r_Vertex_Pool()[n_node1].n_Dimension() == 3);
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief updates the edge with a new measurement
	 *
	 *	@param[in] r_t_edge is parsed edge
	 */
	inline void Update(const CParserBase::TEdgeP2SC3D &r_t_edge)
	{
		CBaseEdgeImpl<CEdgeP2SC3D, MakeTypelist(CVertexSCam, CVertexXYZ),
			3>::Update(r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma);
	}

	/**
	 *	@brief calculates jacobians, expectation and error
	 *
	 *	@param[out] r_t_jacobian0 is jacobian, associated with the first vertex
	 *	@param[out] r_t_jacobian1 is jacobian, associated with the second vertex
	 *	@param[out] r_v_expectation is expecation vector
	 *	@param[out] r_v_error is error vector
	 */
	inline void Calculate_Jacobians_Expectation_Error(Eigen::Matrix<double, 3, 6> &r_t_jacobian0,
		Eigen::Matrix<double, 3, 3> &r_t_jacobian1, Eigen::Matrix<double, 3, 1> &r_v_expectation,
		Eigen::Matrix<double, 3, 1> &r_v_error) const // change dimensionality of eigen types, if required
	{
		CBAJacobians::Project_P2SC(m_p_vertex0->r_v_State() /* Cam */, m_p_vertex0->v_Intrinsics() /* camera intrinsic params */,
			m_p_vertex1->r_v_State() /* XYZ */, r_v_expectation, r_t_jacobian0, r_t_jacobian1); // write your jacobian calculation code here (vertex state vectors are inputs, the rest are the outputs that you need to fill)
		// calculates the expectation and the jacobians

		//std::cout << "me: " << m_v_measurement[0] << " " << m_v_measurement[1] << " " << m_v_measurement[2] << " ----  " <<
		//		r_v_expectation[0] << " " << r_v_expectation[1] << " " << r_v_expectation[2] << std::endl;
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
		Eigen::Matrix<double, 3, 6> p_jacobi0;
		Eigen::Matrix<double, 3, 3> p_jacobi1;
		Eigen::Matrix<double, 3, 1> v_error;
		//Calculate_Jacobians_Expectation_Error(p_jacobi0, p_jacobi1, v_expectation, v_error);
		CBAJacobians::Project_P2SC(m_p_vertex0->r_v_State() /* Cam */, m_p_vertex0->v_Intrinsics(),
				m_p_vertex1->r_v_State() /* XYZ */, v_error);
		// calculates the expectation, error and the jacobians
		v_error -= m_v_measurement; // this is actually negative error, but it is squared below so the chi2 is the same
		// this should be faster, as it does not calculate the jacobians

		return (v_error.transpose() * m_t_sigma_inv).dot(v_error); // ||z_i - h_i(O_i)||^2 lambda_i*/
	}

	/**
	 *	@brief calculates reprojection error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Reprojection_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		/*Eigen::Matrix<double, 2, 6> p_jacobi0;
		Eigen::Matrix<double, 2, 3> p_jacobi1;
		Eigen::Matrix<double, 3, 1> v_expectation, v_error;
		Calculate_Jacobians_Expectation_Error(p_jacobi0, p_jacobi1, v_expectation, v_error);
		// calculates the expectation, error and the jacobians
		//std::cerr << v_error[0] << " " << v_error[1] << " ----  " << (v_error.transpose() * m_t_sigma_inv).dot(v_error) << std::endl;

		return sqrt( (v_error[0] - v_expectation[0])*(v_error[0] - v_expectation[0]) +
				(v_error[1] - v_expectation[1])*(v_error[1] - v_expectation[1]) );*/

		return 0;
	}
};

/**
 *	@brief edge for 3D landmark observation by a spherical camera
 *	@note This is formally a SLAM edge, not BA.
 */
class CEdgeSpheronXYZ : public CBaseEdgeImpl<CEdgeSpheronXYZ, MakeTypelist(CVertexSpheron, CVertexXYZ), 3> {
public:
	/**
	 *	@brief vertex initialization functor
	 *	Calculates vertex position from the first vertex and an XYT edge.
	 */
	class CRelative_to_Absolute_XYZ_Initializer { // this is an object which is used to lazy initialize vertices (copy it)
	protected:
		const Eigen::Matrix<double, 6, 1> &m_r_v_pose1; /**< @brief the first vertex */
		const Eigen::Matrix<double, 3, 1> &m_r_v_delta; /**< @brief the edge, shared by r_v_vertex1 and the vertex being initialized */

	public:
		/**
		 *	@brief default constructor
		 *
		 *	@param[in] r_v_vertex1 is the first vertex
		 *	@param[in] r_v_delta is the edge, shared by r_v_vertex1 and the vertex being initialized
		 */
		inline CRelative_to_Absolute_XYZ_Initializer(const Eigen::Matrix<double, 6, 1> &r_v_vertex1,
			const Eigen::Matrix<double, 3, 1> &r_v_delta) // just change the types, same as above
			:m_r_v_pose1(r_v_vertex1), m_r_v_delta(r_v_delta)
		{}

		/**
		 *	@brief function operator
		 *	@return Returns the value of the vertex being initialized.
		 */
		inline operator CVertexXYZ() const // this function calculates initial prior from the state of the first vertex m_r_v_pose1 and from the edge measurement m_r_edge
		{
			Eigen::Matrix<double, 6, 1> v_pose2;
			Eigen::Matrix<double, 6, 1> relative_pose;
			relative_pose << m_r_v_delta(0), m_r_v_delta(1), m_r_v_delta(2), 0, 0, 0;	//dummy rotation
			/*std::cout << "p1: " << m_r_v_pose1(0) << " " << m_r_v_pose1(1) << " " << m_r_v_pose1(2) << " " << m_r_v_pose1(3) << " " <<
					m_r_v_pose1(4) << " " << m_r_v_pose1(5) << std::endl;
			std::cout << "p2: " << relative_pose(0) << " " << relative_pose(1) << " " << relative_pose(2) << " " << relative_pose(3) << " " <<
					relative_pose(4) << " " << relative_pose(5) << std::endl;*/
			C3DJacobians::Relative_to_Absolute(m_r_v_pose1, relative_pose, v_pose2); // implement your own equation here
			/*std::cout << "res: " << v_pose2(0) << " " << v_pose2(1) << " " << v_pose2(2) << " " << v_pose2(3) << " " <<
					v_pose2(4) << " " << v_pose2(5) << std::endl;*/
			return CVertexXYZ(v_pose2.segment<3>(0));
		}
	};

public:
	__BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, just copy this

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CEdgeSpheronXYZ() // copy this
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
	CEdgeSpheronXYZ(const CParserBase::TEdgeSpheronXYZ &r_t_edge, CSystem &r_system) // this is edge constructor how it is called in the parse loop; you need to change type of the edge
		:CBaseEdgeImpl<CEdgeSpheronXYZ, MakeTypelist(CVertexSpheron, CVertexXYZ), 3>(r_t_edge.m_n_node_0,
		r_t_edge.m_n_node_1, r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		//fprintf(stderr, "%f %f %f\n", r_t_edge.m_v_delta(0), r_t_edge.m_v_delta(1), r_t_edge.m_v_delta(2));

		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexSpheron>(r_t_edge.m_n_node_0, CInitializeNullVertex<>());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexXYZ>(r_t_edge.m_n_node_1,
			CRelative_to_Absolute_XYZ_Initializer(m_p_vertex0->r_v_State(), r_t_edge.m_v_delta)); // rename your initializer if required

		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		//print me edge
		//std::cout << r_t_edge.m_n_node_0 << " to " << r_t_edge.m_n_node_1 << std::endl;
		//std::cout << r_t_edge.m_v_delta << std::endl;

		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_0].n_Dimension() == 6);
		_ASSERTE(r_system.r_Vertex_Pool()[r_t_edge.m_n_node_1].n_Dimension() == 3); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
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
	CEdgeSpheronXYZ(size_t n_node0, size_t n_node1, const Eigen::Matrix<double, 3, 1> &r_v_delta,
		const Eigen::Matrix<double, 3, 3> &r_t_inv_sigma, CSystem &r_system) // this is edge constructor how it is called in the parse loop; you need to change type of the edge
		:CBaseEdgeImpl<CEdgeSpheronXYZ, MakeTypelist(CVertexSpheron, CVertexXYZ), 3>(n_node0,
		n_node1, r_v_delta, r_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
	{
		//fprintf(stderr, "%f %f %f\n", r_t_edge.m_v_delta(0), r_t_edge.m_v_delta(1), r_t_edge.m_v_delta(2));

		m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexSpheron>(n_node0, CInitializeNullVertex<>());
		m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexXYZ>(n_node1,
			CRelative_to_Absolute_XYZ_Initializer(m_p_vertex0->r_v_State(), r_v_delta)); // rename your initializer if required

		// get vertices (initialize if required)
		// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"

		_ASSERTE(r_system.r_Vertex_Pool()[n_node0].n_Dimension() == 6);
		_ASSERTE(r_system.r_Vertex_Pool()[n_node1].n_Dimension() == 3); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
		// make sure the dimensionality is correct (might not be)
	}

	/**
	 *	@brief updates the edge with a new measurement
	 *
	 *	@param[in] r_t_edge is parsed edge
	 */
	inline void Update(const CParserBase::TEdgeSpheronXYZ &r_t_edge)
	{
		CBaseEdgeImpl<CEdgeSpheronXYZ, MakeTypelist(CVertexSpheron, CVertexXYZ),
			3>::Update(r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma);
	}

	/**
	 *	@brief calculates jacobians, expectation and error
	 *
	 *	@param[out] r_t_jacobian0 is jacobian, associated with the first vertex
	 *	@param[out] r_t_jacobian1 is jacobian, associated with the second vertex
	 *	@param[out] r_v_expectation is expecation vector
	 *	@param[out] r_v_error is error vector
	 */
	inline void Calculate_Jacobians_Expectation_Error(Eigen::Matrix<double, 3, 6> &r_t_jacobian0,
		Eigen::Matrix<double, 3, 3> &r_t_jacobian1, Eigen::Matrix<double, 3, 1> &r_v_expectation,
		Eigen::Matrix<double, 3, 1> &r_v_error) const // change dimensionality of eigen types, if required
	{
		C3DJacobians::Absolute_to_Relative_Landmark(m_p_vertex0->r_v_State(),
			m_p_vertex1->r_v_State(), r_v_expectation, r_t_jacobian0, r_t_jacobian1); // write your jacobian calculation code here (vertex state vectors are inputs, the rest are the outputs that you need to fill)
		// calculates the expectation and the jacobians

		/*std::cout << "v1: " << m_p_vertex0->r_v_State()(0) << " " << m_p_vertex0->r_v_State()(1) << " " <<
				m_p_vertex0->r_v_State()(2) << " " << m_p_vertex0->r_v_State()(3) << " " << m_p_vertex0->r_v_State()(4) << " " <<
				m_p_vertex0->r_v_State()(5) << std::endl;
		std::cout << "v2: " << m_p_vertex1->r_v_State()(0) << " " << m_p_vertex1->r_v_State()(1) << " " <<
						m_p_vertex1->r_v_State()(2) << std::endl;*/
		//std::cout << "expect: " << r_v_expectation(0) << " " << r_v_expectation(1) << " " << r_v_expectation(2) << " -- ";
		///std::cout << m_v_measurement(0) << " " << m_v_measurement(1) << " " << m_v_measurement(2) << std::endl;

		//const Eigen::Matrix<double, 6, 1> &p1 = m_v_measurement;
		//const Eigen::Matrix<double, 6, 1> &p2 = r_v_expectation;

		//r_v_error(0) = p1(0) - p2(0);
		//r_v_error(1) = p1(1) - p2(1);
		//r_v_error(2) = p1(2) - p2(2);
		r_v_error = m_v_measurement - r_v_expectation;
	}

	/**
	 *	@brief calculates chi-square error
	 *	@return Returns (unweighted) chi-square error for this edge.
	 */
	inline double f_Chi_Squared_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
	{
		/*Eigen::Matrix<double, 3, 6> p_jacobi0;
		Eigen::Matrix<double, 3, 3> p_jacobi1;
		Eigen::Matrix<double, 3, 1> v_expectation, v_error;
		Calculate_Jacobians_Expectation_Error(p_jacobi0, p_jacobi1, v_expectation, v_error);*/
		// calculates the expectation, error and the jacobians

		//std::cout << (v_error.transpose() * m_t_sigma_inv).dot(v_error) << std::endl;

		Eigen::Matrix<double, 3, 1> v_error;
		C3DJacobians::Absolute_to_Relative_Landmark(m_p_vertex0->r_v_State(), m_p_vertex1->r_v_State(), v_error);

		v_error -= m_v_measurement; // this is actually negative error, but it is squared below so the chi2 is the same
		// this should be faster, as it does not calculate the jacobians

		return (v_error.transpose() * m_t_sigma_inv).dot(v_error); // ||z_i - h_i(O_i)||^2 lambda_i*/
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
		return "unknown edge type occured 1";
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

/**
 *	@brief edge traits for BA solver
 */
template <class CParsedStructure>
class CBASEdgeTraits {
public:
	typedef CFailOnEdgeType _TyEdge; /**< @brief it should fail on unknown edge types */

	/**
	 *	@brief gets reason for error
	 *	@return Returns const null-terminated string, containing
	 *		description of the error (human readable).
	 */
	static const char *p_s_Reason()
	{
		return "unknown edge type occured 2";
	}
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBASEdgeTraits<CParserBase::TEdgeP2SC3D> {
public:
	typedef CEdgeP2SC3D _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBASEdgeTraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBASEdgeTraits<CParserBase::TVertexSCam3D> {
public:
	typedef CVertexSCam _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for BA solver
 */
template <class CParsedStructure>
class CBASVertexTraits {
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
class CBASVertexTraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for BA solver
 */
template <>
class CBASVertexTraits<CParserBase::TVertexSCam3D> {
public:
	typedef CVertexSCam _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <class CParsedStructure>
class CSpheronEdgeTraits {
public:
	typedef CFailOnEdgeType _TyEdge; /**< @brief it should fail on unknown edge types */

	/**
	 *	@brief gets reason for error
	 *	@return Returns const null-terminated string, containing
	 *		description of the error (human readable).
	 */
	static const char *p_s_Reason()
	{
		return "unknown edge type occured 2";
	}
};

/**
 *	@brief edge traits for Spheron solver
 */
template <>
class CSpheronEdgeTraits<CParserBase::TEdgeSpheronXYZ> {
public:
	typedef CEdgeSpheronXYZ _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for Spheron solver
 */
template <>
class CSpheronEdgeTraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for Spheron solver
 */
template <>
class CSpheronEdgeTraits<CParserBase::TVertexSpheron> {
public:
	typedef CVertexSpheron _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for Spheron solver
 */
template <class CParsedStructure>
class CSpheronVertexTraits {
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
 *	@brief vertex traits for Spheron solver
 */
template <>
class CSpheronVertexTraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for Spheron solver
 */
template <>
class CSpheronVertexTraits<CParserBase::TVertexSpheron> {
public:
	typedef CVertexSpheron _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <class CParsedStructure>
class CBAIntrinsicsEdgeTraits {
public:
	typedef CFailOnEdgeType _TyEdge; /**< @brief it should fail on unknown edge types */

	/**
	 *	@brief gets reason for error
	 *	@return Returns const null-terminated string, containing
	 *		description of the error (human readable).
	 */
	static const char *p_s_Reason()
	{
		return "unknown edge type occured 1";
	}
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBAIntrinsicsEdgeTraits<CParserBase::TEdgeP2CI3D> {
public:
	typedef CEdgeP2CI3D _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBAIntrinsicsEdgeTraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBAIntrinsicsEdgeTraits<CParserBase::TVertexCam3D> {
public:
	typedef CVertexCam _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBAIntrinsicsEdgeTraits<CParserBase::TVertexIntrinsics> {
public:
	typedef CVertexIntrinsics _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for BA solver
 */
template <class CParsedStructure>
class CBAIntrinsicsVertexTraits {
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
class CBAIntrinsicsVertexTraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for BA solver
 */
template <>
class CBAIntrinsicsVertexTraits<CParserBase::TVertexCam3D> {
public:
	typedef CVertexCam _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief vertex traits for BA solver
 */
template <>
class CBAIntrinsicsVertexTraits<CParserBase::TVertexIntrinsics> {
public:
	typedef CVertexIntrinsics _TyVertex; /**< @brief the edge type to construct from the parsed type */
};

#endif // __BA_TYPES_INCLUDED
