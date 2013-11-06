/*
								+-----------------------------------+
								|                                   |
								|  ***  Solvers config traits  ***  |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2013  |
								|                                   |
								|          SolverTraits.h           |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __SOLVER_TRAITS_INCLUDED
#define __SOLVER_TRAITS_INCLUDED

/**
 *	@file include/slam/SolverTraits.h
 *	@author -tHE SWINe-
 *	@brief a file, containing helper traits for ConfigSolvers.h (to keep it clean and short)
 *	@date 2013
 */

/**
 *	@brief token, used as a placeholder for solver templates,
 *		that were not included in the build
 *
 *	@tparam CSystem is the system type (unused)
 *	@tparam CLinearSolver is a linear solver (unused)
 *	@tparam CBlockSizes is a list of matrix block sizes (unused)
 */
template <class CSystem, class CLinearSolver, class CBlockSizes>
class CSolverNotIncluded {};

/**
 *	@brief token, used as a placeholder for solver templates,
 *		that are not supported by SE types
 *
 *	@tparam CSystem is the system type (unused)
 *	@tparam CLinearSolver is a linear solver (unused)
 *	@tparam CBlockSizes is a list of matrix block sizes (unused)
 */
template <class CSystem, class CLinearSolver, class CBlockSizes>
class CSolverNotSupported {};

/**
 *	@brief pair of nonlinear solver id and the type, along with traits
 *
 *	@tparam n_solver_type_id is nonlinear solver id
 *	@tparam CNonlinearSolverType is nonlinear solver template name
 */
template <const int n_solver_type_id,
	template <class, class, class> class CNonlinearSolverType>
class CSolverTypeIdPair {
public:
	/**
	 *	@brief configuration, stored as enum
	 */
	enum {
		solver_type_Id = n_solver_type_id /**< @brief nonlinear solver id */
	};

	/**
	 *	@brief runs the application, using this solver
	 *
	 *	@tparam CSystemType is system type (derived from CFlatSystem)
	 *	@tparam CEdgeTraitsType is edge traits template name
	 *	@tparam CVertexTraitsType is vertex traits template name
	 *	@tparam CParseLoopType is parse loop template name
	 *	@tparam CRunEnvironment is type of environtment the solver is supposed to run in
	 *
	 *	@param[in] t_env is environtment the solver is supposed to run in
	 *
	 *	@return Returns true on success, false on failure.
	 */
	template <class CSystemType, template <class> class CEdgeTraitsType,
		template <class> class CVertexTraitsType,
		template <class, class, template <class> class vcneedsnamehere,
		template <class> class vcneedsnamehereaswell>
		class CParseLoopType, class CRunEnvironment>
	static inline bool Run_MainApp(CRunEnvironment t_env) // throw(std::runtime_error, std::bad_alloc)
	{
		return t_env.template Run<CSystemType, CNonlinearSolverType,
			CEdgeTraitsType, CVertexTraitsType, CParseLoopType>();
		// run with parameters
	}
};

/**
 *	@brief pair of nonlinear solver id and the type,
 *		along with traits (specialization for solvers that were not included)
 *	@tparam n_solver_type_id is nonlinear solver id
 */
template <const int n_solver_type_id>
class CSolverTypeIdPair<n_solver_type_id, CSolverNotIncluded> {
public:
	/**
	 *	@brief configuration, stored as enum
	 */
	enum {
		solver_type_Id = n_solver_type_id /**< @brief nonlinear solver id */
	};

	/**
	 *	@brief runs the application, using this solver
	 *
	 *	@tparam CSystemType is system type (derived from CFlatSystem)
	 *	@tparam CEdgeTraitsType is edge traits template name
	 *	@tparam CVertexTraitsType is vertex traits template name
	 *	@tparam CParseLoopType is parse loop template name
	 *	@tparam CRunEnvironment is type of environtment the solver is supposed to run in
	 *
	 *	@param[in] t_env is environtment the solver is supposed to run in
	 *
	 *	@return Returns true on success, false on failure.
	 */
	template <class CSystemType, template <class> class CEdgeTraitsType,
		template <class> class CVertexTraitsType,
		template <class, class, template <class> class vcneedsnamehere,
		template <class> class vcneedsnamehereaswell>
		class CParseLoopType, class CRunEnvironment>
	static inline bool Run_MainApp(CRunEnvironment UNUSED(t_env))
	{
		fprintf(stderr, "error: the selected solver was not included\n");
		return false;
	}
};

/**
 *	@brief pair of nonlinear solver id and the type,
 *		along with traits (specialization for unsupported solvers)
 *	@tparam n_solver_type_id is nonlinear solver id
 */
template <const int n_solver_type_id>
class CSolverTypeIdPair<n_solver_type_id, CSolverNotSupported> {
public:
	/**
	 *	@brief configuration, stored as enum
	 */
	enum {
		solver_type_Id = n_solver_type_id /**< @brief nonlinear solver id */
	};

	/**
	 *	@brief runs the application, using this solver
	 *
	 *	@tparam CSystemType is system type (derived from CFlatSystem)
	 *	@tparam CEdgeTraitsType is edge traits template name
	 *	@tparam CVertexTraitsType is vertex traits template name
	 *	@tparam CParseLoopType is parse loop template name
	 *	@tparam CRunEnvironment is type of environtment the solver is supposed to run in
	 *
	 *	@param[in] t_env is environtment the solver is supposed to run in
	 *
	 *	@return Returns true on success, false on failure.
	 */
	template <class CSystemType, template <class> class CEdgeTraitsType,
		template <class> class CVertexTraitsType,
		template <class, class, template <class> class vcneedsnamehere,
		template <class> class vcneedsnamehereaswell>
		class CParseLoopType, class CRunEnvironment>
	static inline bool Run_MainApp(CRunEnvironment UNUSED(t_env))
	{
		fprintf(stderr, "error: the selected solver is not supported by the SE types\n");
		return false;
	}
};

#endif // __SOLVER_TRAITS_INCLUDED
