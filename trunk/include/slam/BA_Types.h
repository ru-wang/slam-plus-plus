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
//#include "slam/2DSolverBase.h"
//#include "slam/3DSolverBase.h"
#include "slam/BASolverBase.h"
#include <iostream> //SOSO

/*class CVertexIntrinsics : public CSEBaseVertexImpl<CVertexIntrinsics, 4> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this
    inline CVertexPoseBA() // copy this
    {}

    inline CVertexIntrinsics(const Eigen::VectorXd &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexIntrinsics, 4>(r_v_state) // change the dimension here as well
    {}

    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
    	m_v_state.segment<4>() += r_v_delta.segment<4>(m_n_order);
    }
};*/

class CVertexCam : public CSEBaseVertexImpl<CVertexCam, 6> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

    //vertex cam should hold own camera params
    Eigen::Matrix<double, 5, 1> intrinsics;

    inline CVertexCam() // copy this
    {}

    inline CVertexCam(const Eigen::VectorXd &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexCam, 6>(r_v_state.head<6>()) // change the dimension here as well
    {
    	intrinsics = r_v_state.tail<5>();
    }

	/*inline void Alloc_HessianBlocks(CUberBlockMatrix &r_lambda)
	{
		m_p_HtSiH = r_lambda.p_FindBlock(m_n_order, m_n_order, 2, n_dimension, true, false);
		// find a block for hessian on the diagonal (edges can't do it, would have conflicts)
	}*/ // no!

    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
        //focal lengths and principal points are added normally
    	//SOSO: ignore adding of fx and cx
    	//m_v_state.segment<4>(6) += r_v_delta.segment<4>(m_n_order + 6);

    	//pose is added as SE3
    	m_v_state.head<3>() += r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated

		//sum the rotations
		Eigen::Matrix3d pQ = CBase3DSolver::C3DJacobians::Operator_rot(m_v_state.segment<3>(3));
		Eigen::Matrix3d dQ = CBase3DSolver::C3DJacobians::Operator_rot(r_v_delta.segment<3>(m_n_order + 3));

		Eigen::Matrix3d QQ;// = pQ * dQ;
		QQ.noalias() = pQ * dQ; // noalias says that the destination is not the same as the multiplicands and that the multiplication can be carried out without allocating intermediate storage
		m_v_state.segment<3>(3) = CBaseBASolver::CBAJacobians::Operator_arot(QQ); // also, no need for the intermediate object
    }
};

class CVertexXYZ : public CSEBaseVertexImpl<CVertexXYZ, 3> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
public:
    __BA_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this

    inline CVertexXYZ() // copy this
    {}

    inline CVertexXYZ(const Eigen::Vector3d &r_v_state) // copy this, change the dimension of the vector to appropriate
        :CSEBaseVertexImpl<CVertexXYZ, 3>(r_v_state) // change the dimension here as well
    {}

	/*inline void Alloc_HessianBlocks(CUberBlockMatrix &r_lambda)
	{
		m_p_HtSiH = r_lambda.p_FindBlock(m_n_order, m_n_order, 2, n_dimension, true, false);
		// find a block for hessian on the diagonal (edges can't do it, would have conflicts)
	}*/ // no!

    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
    {
    	//pose is added as SE3
    	m_v_state += r_v_delta.segment<3>(m_n_order); // can select 3D segment starting at m_n_order; SSE accelerated
    }

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
		CBaseBASolver::CBAJacobians::Project_P2C(m_p_vertex0->v_State() /* Cam */, m_p_vertex0->intrinsics /* camera intrinsic params */,
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
class CBATraits {
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
class CBATraits<CParserBase::TEdgeP2C3D> {
public:
	typedef CEdgeP2C3D _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBATraits<CParserBase::TVertexXYZ> {
public:
	typedef CVertexXYZ _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief edge traits for BA solver
 */
template <>
class CBATraits<CParserBase::TVertexCam3D> {
public:
	typedef CVertexCam _TyEdge; /**< @brief the edge type to construct from the parsed type */
};

/**
 *	@brief parser loop consumer, working with flat system
 *
 *	@tparam CSystem is flat system type
 *	@tparam CNonlinearSolver is nonlinear solver type
 *	@tparam CEdgeTraits is edge traits template (performs
 *		lookup of edge representation type by parser type)
 *
 *	@note This is currently limited to binary edges. A modification to edge map would be needed
 *		(or removal of the edge map, which would work if there are no duplicate edges).
 *	@todo This is redundant, need to make vertex traits and use only a single parse loop.
 */
template <class CSystem, class CNonlinearSolver, template <class> class CEdgeTraits>
class CBAParseLoop /*: public CParser::CParserAdaptor*/ {
protected:
	CSystem &m_r_system; /**< @brief reference to the system being parsed into */
	CNonlinearSolver &m_r_solver; /**< @brief reference to the solver (for incremental solving) */
	std::map<std::pair<size_t, size_t>, size_t> m_edge_map; /**< @brief map of edges by vertices (currently handles binary edges only) */

public:
	/**
	 *	@brief default constructor; sets the parse loop up
	 *
	 *	@param[in] r_system is reference to the system being parsed into
	 *	@param[in] r_solver is reference to the solver (for incremental solving)
	 */
	inline CBAParseLoop(CSystem &r_system, CNonlinearSolver &r_solver)
		:m_r_system(r_system), m_r_solver(r_solver)
	{}

	/**
	 *	@brief processes an edge, based on it's type and the edge traits
	 *	@tparam CParsedEdge is parsed edge type
	 *	@param[in] r_edge is reference to the parsed edge
	 *	@note This function throws std::bad_alloc, and also
	 *		std::runtime_error as a means of error reporting.
	 */
	template <class CParsedEdge>
	void AppendSystem(const CParsedEdge &r_edge) // throws(std::bad_alloc, std::runtime_error)
	{
		CProcessEdge<CParsedEdge, typename CEdgeTraits<CParsedEdge>::_TyEdge>::Do(m_r_system,
			m_r_solver, r_edge, m_edge_map);
	}

	class CMakeVertexXYZ {
		const Eigen::Matrix<double, 3, 1> &m_r_init;

	public:
		inline CMakeVertexXYZ(const Eigen::Matrix<double, 3, 1> &r_init)
			:m_r_init(r_init)
		{}

		inline operator CVertexXYZ() const
		{
			return CVertexXYZ(m_r_init);
		}
	};

	void AppendSystem(const CParserBase::TVertexXYZ &r_vertex) // throws(std::bad_alloc, std::runtime_error)
	{
		m_r_system.template r_Get_Vertex<CVertexXYZ>(r_vertex.m_n_id, CMakeVertexXYZ(r_vertex.m_v_position));
	}

	class CMakeVertexCam3D {
		const Eigen::Matrix<double, 11, 1> &m_r_init;

	public:
		inline CMakeVertexCam3D(const Eigen::Matrix<double, 11, 1> &r_init)
			:m_r_init(r_init)
		{}

		inline operator CVertexCam() const
		{
			return CVertexCam(m_r_init);
		}
	};

	void AppendSystem(const CParserBase::TVertexCam3D &r_vertex) // throws(std::bad_alloc, std::runtime_error)
	{
		m_r_system.template r_Get_Vertex<CVertexCam>(r_vertex.m_n_id, CMakeVertexCam3D(r_vertex.m_v_position));
	}

protected:
	/**
	 *	@brief edge processing functor
	 *
	 *	@tparam CParsedEdge is type of edge as parsed
	 *	@tparam CRepresentation is type of edge as represented in the system
	 */
	template <class CParsedEdge, class CRepresentation>
	class CProcessEdge {
	public:
		/**
		 *	@brief edge processing function
		 *
		 *	@param[in] r_system is reference to the system being parsed into
		 *	@param[in] r_solver is reference to the solver (for incremental solving)
		 *	@param[in] r_edge is reference to the parsed edge
		 *	@param[in,out] r_edge_map is map of edge indices by edge vertices
		 *
		 *	@note This function throws std::bad_alloc.
		 */
		static inline void Do(CSystem &r_system, CNonlinearSolver &r_solver,
			const CParsedEdge &r_edge, std::map<std::pair<size_t, size_t>, size_t> &r_edge_map) // throws(std::bad_alloc)
		{
			std::pair<size_t, size_t> edge_verts(r_edge.m_n_node_0, r_edge.m_n_node_1);
			// this needs to be changed for multi-edge datasets

			std::map<std::pair<size_t, size_t>, size_t>::const_iterator p_edge_id_it;
			if((p_edge_id_it = r_edge_map.find(edge_verts)) != r_edge_map.end()) {
				size_t n_edge_id = (*p_edge_id_it).second;
				// gets repeated edge id from the map

				CRepresentation &r_rep_edge = *(CRepresentation*)&r_system.r_Edge_Pool()[n_edge_id];
				// gets the edge itself (kind of relies on the type correctness but that should be ok)

				r_rep_edge.Update(r_edge); // t_odo - not implemented for SE2 or SE3
				// update the reading

				//printf("handling dup edge ...\n");
				r_solver.Incremental_Step(r_rep_edge);
				// notify the solver of change
			} else {
				r_edge_map[edge_verts] = r_system.r_Edge_Pool().n_Size();
				// record edge id

				CRepresentation &r_rep_edge = r_system.r_Add_Edge(CRepresentation(r_edge, r_system));
				// add the edge to the system (convert parsed edge to internal representation)

				r_solver.Incremental_Step(r_rep_edge);
				// call solver with the new edge
			}
		}
	};

	/**
	 *	@brief edge processing functor (specialization for ignored edge types)
	 *
	 *	@tparam CParsedEdge is type of edge as parsed
	 */
	template <class CParsedEdge>
	class CProcessEdge<CParsedEdge, CIgnoreEdgeType> {
	public:
		/**
		 *	@brief edge processing function (for edges to be ignored)
		 *
		 *	@param[in] r_system is reference to the system being parsed into (unused)
		 *	@param[in] r_solver is reference to the solver (unused)
		 *	@param[in] r_edge is reference to the parsed edge (unused)
		 *	@param[in] r_edge_map is map of edge indices by edge vertices (unused)
		 */
		static inline void Do(CSystem &UNUSED(r_system),
			CNonlinearSolver &UNUSED(r_solver), const CParsedEdge &UNUSED(r_edge),
			std::map<std::pair<size_t, size_t>, size_t> &UNUSED(r_edge_map))
		{}
	};

	/**
	 *	@brief edge processing functor (specialization for edge
	 *		types that cause the parse loop to fail)
	 *
	 *	@tparam CParsedEdge is type of edge as parsed
	 */
	template <class CParsedEdge>
	class CProcessEdge<CParsedEdge, CFailOnEdgeType> {
	public:
		/**
		 *	@brief edge processing function (for edges that cause the parse loop to fail)
		 *
		 *	@param[in] r_system is reference to the system being parsed into (unused)
		 *	@param[in] r_solver is reference to the solver (unused)
		 *	@param[in] r_edge is reference to the parsed edge (unused)
		 *	@param[in] r_edge_map is map of edge indices by edge vertices (unused)
		 *
		 *	@note This function throws std::runtime_error as a means of error reporting.
		 */
		static inline void Do(CSystem &UNUSED(r_system), CNonlinearSolver &UNUSED(r_solver),
			const CParsedEdge &UNUSED(r_edge), std::map<std::pair<size_t, size_t>, size_t> &UNUSED(r_edge_map)) // throws(std::runtime_error)
		{
			typedef CEdgeTraits<CParsedEdge> TEdgeType; // g++ doesn't like 'typename' here
			throw std::runtime_error(TEdgeType::p_s_Reason());
			// "CParseLoop encountered edge type that is not permitted by the configuration"
		}
	};
};

#endif
