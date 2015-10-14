/*
								+----------------------------------+
								|                                  |
								|   ***  Parse loop adaptor  ***   |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2012  |
								|                                  |
								|            ParseLoop.h           |
								|                                  |
								+----------------------------------+
*/

#pragma once
#ifndef __PARSE_LOOP_INCLUDED
#define __PARSE_LOOP_INCLUDED

/**
 *	@file include/slam/ParseLoop.h
 *	@brief parse loop adaptor
 *	@author -tHE SWINe-
 *	@date 2012-09-03
 */

//#include "slam/Parser.h"
#include "slam/FlatSystem.h"
#include <map>

/**
 *	@brief set this in edge traits for the CParseLoop to ignore this edge
 */
class CIgnoreEdgeType {};

/**
 *	@brief set this in edge traits for the CParseLoop to throw error on this edge
 *	@note The trait must implement p_s_Reason() function, that will give error description
 */
class CFailOnEdgeType {};

/**
 *	@brief set this in vertex traits for the CParseLoop to ignore this vertex
 */
typedef CIgnoreEdgeType CIgnoreVertexType;

/**
 *	@brief set this in vertex traits for the CParseLoop to throw error on this vertex
 *	@note The trait must implement p_s_Reason() function, that will give error description
 */
typedef CFailOnEdgeType CFailOnVertexType;

/**
 *	@brief vertex traits for solvers that ignore vertices
 */
template <class CParsedStructure>
class CIgnoreAllVertexTraits {
public:
	typedef CIgnoreVertexType _TyVertex; /**< @brief it should quietly ignore unknown vertex types */
};

/**
 *	@brief parser loop consumer, working with flat system
 *
 *	@tparam CSystem is flat system type
 *	@tparam CNonlinearSolver is nonlinear solver type
 *	@tparam CEdgeTraits is edge traits template (performs
 *		lookup of edge representation type by parser type)
 *	@tparam CVertexTraits is vertex traits template (performs
 *		lookup of vertex representation type by parser type)
 *
 *	@note This is currently limited to binary edges. A modification to edge map would be needed
 *		(or removal of the edge map, which would work if there are no duplicate edges).
 */
template <class CSystem, class CNonlinearSolver, template <class> class CEdgeTraits,
	template <class> class CVertexTraits = CIgnoreAllVertexTraits>
class CParseLoop {
protected:
	CSystem &m_r_system; /**< @brief reference to the system being parsed into */
	CNonlinearSolver &m_r_solver; /**< @brief reference to the solver (for incremental solving) */
	//std::map<std::pair<size_t, size_t>, size_t> m_edge_map; /**< @brief map of edges by vertices (currently handles binary edges only) */

public:
	/**
	 *	@brief default constructor; sets the parse loop up
	 *
	 *	@param[in] r_system is reference to the system being parsed into
	 *	@param[in] r_solver is reference to the solver (for incremental solving)
	 */
	inline CParseLoop(CSystem &r_system, CNonlinearSolver &r_solver)
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
	void AppendSystem(const CParsedEdge &r_edge) // throw(std::bad_alloc, std::runtime_error)
	{
		CProcessEdge<CParsedEdge, typename CEdgeTraits<CParsedEdge>::_TyEdge>::Do(m_r_system,
			m_r_solver, r_edge/*, m_edge_map*/);
	}

	/**
	 *	@brief processes a vertex, based on it's type and the vertex traits
	 *	@tparam CParsedVertex is parsed vertex type
	 *	@param[in] r_vertex is reference to the parsed vertex
	 *	@note This function throws std::bad_alloc, and also
	 *		std::runtime_error as a means of error reporting.
	 */
	template <class CParsedVertex>
	void InitializeVertex(const CParsedVertex &r_vertex) // throw(std::bad_alloc, std::runtime_error)
	{
		CProcessVertex<CParsedVertex, typename CVertexTraits<CParsedVertex>::_TyVertex>::Do(m_r_system, r_vertex);
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
			const CParsedEdge &r_edge/*, std::map<std::pair<size_t, size_t>, size_t> &r_edge_map*/) // throw(std::bad_alloc)
		{
#if 0
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
#else // 0
			{ // duplicate edges are handled, can be put in the system
#endif // 0
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
			CNonlinearSolver &UNUSED(r_solver), const CParsedEdge &UNUSED(r_edge)/*,
			std::map<std::pair<size_t, size_t>, size_t> &UNUSED(r_edge_map)*/)
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
			const CParsedEdge &UNUSED(r_edge)/*, std::map<std::pair<size_t, size_t>, size_t> &UNUSED(r_edge_map)*/) // throw(std::runtime_error)
		{
			typedef CEdgeTraits<CParsedEdge> TEdgeType; // g++ doesn't like 'typename' here
			throw std::runtime_error(TEdgeType::p_s_Reason());
			// "CParseLoop encountered edge type that is not permitted by the configuration"
		}
	};

	/**
	 *	@brief vertex processing functor
	 *
	 *	@tparam CParsedVertex is type of vertex as parsed
	 *	@tparam CRepresentation is type of vertex as represented in the system
	 */
	template <class CParsedVertex, class CRepresentation>
	class CProcessVertex {
	public:
		/**
		 *	@brief vertex processing function
		 *
		 *	@param[in] r_system is reference to the system being parsed into
		 *	@param[in] r_vertex is reference to the parsed vertex
		 *
		 *	@note This function throws std::bad_alloc.
		 */
		static inline void Do(CSystem &r_system, const CParsedVertex &r_vertex) // throw(std::bad_alloc)
		{
			r_system.template r_Get_Vertex<CRepresentation>(r_vertex.m_n_id, CRepresentation(r_vertex));
			// add the vertex to the system (convert parsed vertex to internal representation)
		}
	};

	/**
	 *	@brief vertex processing functor (specialization for ignored vertex types)
	 *
	 *	@tparam CParsedVertex is type of vertex as parsed
	 */
	template <class CParsedVertex>
	class CProcessVertex<CParsedVertex, CIgnoreVertexType> {
	public:
		/**
		 *	@brief vertex processing function (for vertexs to be ignored)
		 *
		 *	@param[in] r_system is reference to the system being parsed into (unused)
		 *	@param[in] r_vertex is reference to the parsed vertex (unused)
		 */
		static inline void Do(CSystem &UNUSED(r_system), const CParsedVertex &UNUSED(r_vertex))
		{}
	};

	/**
	 *	@brief vertex processing functor (specialization for vertex
	 *		types that cause the parse loop to fail)
	 *
	 *	@tparam CParsedVertex is type of vertex as parsed
	 */
	template <class CParsedVertex>
	class CProcessVertex<CParsedVertex, CFailOnVertexType> {
	public:
		/**
		 *	@brief vertex processing function (for vertexs that cause the parse loop to fail)
		 *
		 *	@param[in] r_system is reference to the system being parsed into (unused)
		 *	@param[in] r_vertex is reference to the parsed vertex (unused)
		 *
		 *	@note This function throws std::runtime_error as a means of error reporting.
		 */
		static inline void Do(CSystem &UNUSED(r_system), const CParsedVertex &UNUSED(r_vertex)) // throw(std::runtime_error)
		{
			typedef CVertexTraits<CParsedVertex> TVertexType; // g++ doesn't like 'typename' here
			throw std::runtime_error(TVertexType::p_s_Reason());
			// "CParseLoop encountered vertex type that is not permitted by the configuration"
		}
	};
};

#endif // __PARSE_LOOP_INCLUDED