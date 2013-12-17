/*
								+-----------------------------------+
								|                                   |
								|      ***  Base SE types  ***      |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2012  |
								|                                   |
								|            BaseTypes.h            |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __BASE_SE_PRIMITIVE_TYPES_INCLUDED
#define __BASE_SE_PRIMITIVE_TYPES_INCLUDED

/**
 *	@file include/slam/BaseTypes.h
 *	@brief base SE primitive types
 *	@author -tHE SWINe-
 *	@date 2012-09-21
 */

#include "slam/FlatSystem.h"
#include "slam/ParseLoop.h"

// looking for #define __BASE_TYPES_USE_ID_ADDRESSING? it is in slam/FlatSystem.h ...

/*#if !defined(_WIN32) && !defined(_WIN64)
#define __SE2_TYPES_USE_ALIGNED_VECTORS
#define __SE2_TYPES_ALIGN_OPERATOR_NEW EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#else // !_WIN32 && !_WIN64
#define __SE2_TYPES_ALIGN_OPERATOR_NEW
#endif // !_WIN32 && !_WIN64*/ // does not work // todo - investigate (low priority, only a small fraction of time is spent there)
#define __SE2_TYPES_ALIGN_OPERATOR_NEW
#define __SE3_TYPES_ALIGN_OPERATOR_NEW
#define __BA_TYPES_ALIGN_OPERATOR_NEW
// decide whether to align types // todo - document it and add manual overrides

class CSEBaseVertex; // forward declaration

/**
 *	@brief base edge for SE types
 */
class CSEBaseEdge { // need to have this before vertices because they use it's interface. messy. separate .h and .cpp
protected:
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
	size_t m_n_id; /**< @brief edge id (block row in a matrix where the edge block is inseted) */
#endif // __BASE_TYPES_USE_ID_ADDRESSING
	size_t m_p_vertex_id[2]; /**< @brief ids of referenced vertices */
	size_t m_n_order; /**< @brief edge order (row in a matrix where the edge block is inseted) */

public:
	__SE2_TYPES_ALIGN_OPERATOR_NEW

	/**
	 *	@brief gets edge dimension
	 *	@return Returns edge dimension.
	 */
	inline size_t n_Dimension() const
	{
		return n_Dimension_Derived();
	}

	/**
	 *	@brief gets id of edge vertex
	 *	@param[in] n_vertex is zero-based vertex index
	 *	@return Returns id of the selected edge vertex.
	 */
	inline size_t n_Vertex_Id(size_t n_vertex) const
	{
		_ASSERTE(n_vertex < 2);
		return m_p_vertex_id[n_vertex];
	}

	/**
	 *	@brief gets edge order
	 *	@return Returns edge order (row in a matrix where the edge block is inseted).
	 */
	inline size_t n_Order() const
	{
		return m_n_order;
	}

	/**
	 *	@brief sets edge order
	 *	@param[in] n_first_element_index is edge order
	 *		(row in a matrix where the edge block is inseted).
	 */
	inline void Set_Order(size_t n_first_element_index)
	{
		m_n_order = n_first_element_index;
	}

#ifdef __BASE_TYPES_USE_ID_ADDRESSING

	/**
	 *	@brief gets edge id
	 *	@return Returns edge id (block row in a matrix where the edge block is inseted).
	 */
	inline size_t n_Id() const
	{
		return m_n_id;
	}

	/**
	 *	@brief sets edge order
	 *	@param[in] n_id is edge order
	 *		(block row in a matrix where the edge block is inseted).
	 */
	inline void Set_Id(size_t n_id)
	{
		m_n_id = n_id;
	}

#endif // __BASE_TYPES_USE_ID_ADDRESSING

	/**
	 *	@brief calculates chi-square error
	 *	@return Returns (unweighted) chi-square error.
	 *	@note This is supposed to be divided by (m - n), where m is number of edges
	 *		and n is number of degrees of freedom (equals to n_Dimension()).
	 */
	inline double f_Chi_Squared_Error() const
	{
		return f_Chi_Squared_Error_Derived();
	}

#ifdef __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific functions ---

	/**
	 *	@brief allocates hessian block matrices
	 *	@param[out] r_A is the target matrix where the blocks are stored
	 *	@note This function is required for CNonlinearSolver_A.
	 */
	inline void Alloc_JacobianBlocks(CUberBlockMatrix &r_A)
	{
		Alloc_JacobianBlocks_Derived(r_A);
	}

	/**
	 *	@brief calculates hessians
	 *	@note This function is required for CNonlinearSolver_A.
	 */
	inline void Calculate_Jacobians()
	{
		Calculate_Jacobians_Derived();
	}

	/**
	 *	@brief calculates error, multiplied by the R matrix
	 *		(upper diagonal Choleski of inverse sigma matrix)
	 *	@param[out] r_v_dest is the error vector
	 *	@note This function is required for CNonlinearSolver_A.
	 */
	inline void Get_R_Error(Eigen::VectorXd &r_v_dest) const
	{
		Get_R_Error_Derived(r_v_dest);
	}

	/**
	 *	@brief calculates bare error
	 *	@param[out] r_v_dest is the error vector
	 *	@note This function is required for CNonlinearSolver_A.
	 */
	inline void Get_Error(Eigen::VectorXd &r_v_dest) const
	{
		Get_Error_Derived(r_v_dest);
	}

#endif // __SE_TYPES_SUPPORT_A_SOLVERS // --- ~A-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

	/**
	 *	@brief notifies an edge of a conflict (duplicate edge)
	 *	@return Returns true if this is the original owner of the block.
	 */
	inline bool Notify_HessianBlock_Conflict(double *p_block, CUberBlockMatrix &r_lambda)
	{
		return Notify_HessianBlock_Conflict_Derived(p_block, r_lambda);
	}

	/**
	 *	@brief allocates hessian block matrices
	 *	@param[out] r_lambda is the target matrix where the blocks are stored
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline void Alloc_HessianBlocks(CUberBlockMatrix &r_lambda)
	{
		Alloc_HessianBlocks_Derived(r_lambda);
	}

	/**
	 *	@brief calculates off-diagonal lambda hessian blocks
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline void Calculate_Hessians()
	{
		Calculate_Hessians_Derived();
	}

	virtual double f_Max_VertexHessianDiagValue() const = 0; // todo surround with levenberg support ifdef

	/**
	 *	@brief gets edge contributions for lambda and eta blocks
	 *	@param[in] p_which is vertex, associated with queried contributions
	 *	@return Returns a pair of pointers to memory blocks where lambda contribution
	 *		and eta contribution are stored.
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline std::pair<const double*, const double*>
		t_Get_LambdaEta_Contribution(const CSEBaseVertex *p_which)
	{
		return t_Get_LambdaEta_Contribution_Derived(p_which);
	}

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	/**
	 *	@brief allocates L factor block matrices
	 *	@param[out] r_L is the target matrix where the blocks are stored
	 *	@note This function is required for CNonlinearSolver_L.
	 */
	inline void Alloc_LBlocks(CUberBlockMatrix &r_L) const
	{
		Alloc_LBlocks_Derived(r_L);
	}

	/**
	 *	@brief calculates just a part of the lambda matrix, called omega
	 *
	 *	This is used in CNonlinearSolver_L for incremental updates to the factor.
	 *	It calculates omega from scratch, uses serial execution, as opposed
	 *	to Calculate_Hessians(). While it may seem slower, it can't likely be parallelized
	 *	due to very small size of omega (can be just a single edge).
	 *
	 *	@param[out] r_omega is the omega matrix to be filled (must be initially empty)
	 *	@param[in] n_min_vertex_order is order offset, in elements, used to position
	 *		the hessian blocks in the upper-left corner of omega
	 */
	inline void Calculate_Omega(CUberBlockMatrix &r_omega, size_t n_min_vertex_order) const
	{
		Calculate_Omega_Derived(r_omega, n_min_vertex_order);
	}

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---

protected:
	/**
	 *	@copydoc CSEBaseEdge::n_Dimension()
	 */
	virtual size_t n_Dimension_Derived() const = 0;

	/**
	 *	@copydoc f_Chi_Squared_Error()
	 */
	virtual double f_Chi_Squared_Error_Derived() const = 0;

#ifdef __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific functions ---

	/**
	 *	@copydoc Alloc_JacobianBlocks()
	 */
	virtual void Alloc_JacobianBlocks_Derived(CUberBlockMatrix &r_A) = 0;

	/**
	 *	@copydoc Calculate_Jacobians()
	 */
	virtual void Calculate_Jacobians_Derived() = 0;

	/**
	 *	@copydoc Get_R_Error()
	 */
	virtual void Get_R_Error_Derived(Eigen::VectorXd &r_v_dest) const = 0;

	/**
	 *	@copydoc Get_Error()
	 */
	virtual void Get_Error_Derived(Eigen::VectorXd &r_v_dest) const = 0;

#endif // __SE_TYPES_SUPPORT_A_SOLVERS // --- ~A-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

	/**
	 *	@copydoc Alloc_HessianBlocks()
	 */
	virtual void Alloc_HessianBlocks_Derived(CUberBlockMatrix &r_lambda) = 0;

	/**
	 *	@copydoc Notify_HessianBlock_Conflict_Derived()
	 */
	virtual bool Notify_HessianBlock_Conflict_Derived(double *p_block, CUberBlockMatrix &r_lambda) = 0;

	/**
	 *	@copydoc Calculate_Hessians()
	 */
	virtual void Calculate_Hessians_Derived() = 0;

	/**
	 *	@copydoc t_Get_LambdaEta_Contribution()
	 */
	virtual std::pair<const double*, const double*>
		t_Get_LambdaEta_Contribution_Derived(const CSEBaseVertex *p_which) = 0;

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	/**
	 *	@copydoc Alloc_LBlocks()
	 */
	virtual void Alloc_LBlocks_Derived(CUberBlockMatrix &r_L) const = 0;

	/**
	 *	@copydoc Calculate_Omega()
	 */
	virtual void Calculate_Omega_Derived(CUberBlockMatrix &r_omega, size_t n_min_vertex_order) const = 0;

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---

	// virtual functions call inline implementations in derived types (note that this is just
	// a performance optimization; virtual functions can be used without any limitations)
};

/**
 *	@brief base vertex for SE types
 */
class CSEBaseVertex {
public:
	typedef CSEBaseEdge _TyBaseEdge; /**< @brief base edge type, assumed for this vertex */ 

protected:
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
	size_t m_n_id; /**< @brief vertex id (block column in a matrix where the vertex block is inseted) */
#endif // __BASE_TYPES_USE_ID_ADDRESSING
	size_t m_n_order; /**< @brief vertex order (column in a matrix where the vertex block is inseted) */

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific members ---

	std::vector<_TyBaseEdge*> m_edge_list; /**< @brief a list of edges, referencing this vertex (only used by Lambda solver) */

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific members ---

public:
	__SE2_TYPES_ALIGN_OPERATOR_NEW

	/**
	 *	@brief gets vertex state vector
	 *	@return Returns vertex state vector.
	 *	@note This probably requires vector data allocation and is slow.
	 */
	inline Eigen::VectorXd v_State() const
	{
		return v_State_Derived();
	}

	/**
	 *	@brief gets vertex dimension
	 *	@return Returns vertex dimension.
	 */
	inline size_t n_Dimension() const
	{
		return n_Dimension_Derived();
	}

	/**
	 *	@brief gets vertex order
	 *	@return Returns vertex order (column in a matrix where the vertex block is inseted).
	 */
	inline size_t n_Order() const
	{
		return m_n_order;
	}

	/**
	 *	@brief sets vertex order
	 *	@param[in] n_first_element_index is vertex order
	 *		(column in a matrix where the vertex block is inseted).
	 */
	inline void Set_Order(size_t n_first_element_index)
	{
		m_n_order = n_first_element_index;
	}

#ifdef __BASE_TYPES_USE_ID_ADDRESSING

	/**
	 *	@brief gets vertex id
	 *	@return Returns vertex order (block column in a matrix where the vertex block is inseted).
	 */
	inline size_t n_Id() const
	{
		return m_n_id;
	}

	/**
	 *	@brief sets vertex id
	 *	@param[in] n_id is vertex id
	 *		(block column in a matrix where the vertex block is inseted).
	 */
	inline void Set_Id(size_t n_id)
	{
		m_n_id = n_id;
	}

#endif // __BASE_TYPES_USE_ID_ADDRESSING

	/**
	 *	@brief swaps solution vector with vertex state and minds the manifold space
	 *	@param[in,out] r_v_x is solution vector
	 *	@note This is required by the experimental SPCG solver, no need to implement
	 *		it everywhere (now it is just in SE(2)).
	 */
	inline void SwapState(Eigen::VectorXd &r_v_x)
	{
		SwapState_Derived(r_v_x);
	}

	/**
	 *	@brief saves the vertex state in a vector
	 *	@param[out] r_v_x is the state vector to copy the state to
	 */
	inline void SaveState(Eigen::VectorXd &r_v_x) const
	{
		SaveState_Derived(r_v_x);
	}

	/**
	 *	@brief restores the vertex state from a vector
	 *	@param[out] r_v_x is the state vector to copy the state from
	 */
	inline void LoadState(const Eigen::VectorXd &r_v_x)
	{
		LoadState_Derived(r_v_x);
	}

	/**
	 *	@brief adds delta vector to vertex state and minds the manifold space ("smart" plus)
	 *	@param[in] r_v_delta is delta vector
	 */
	inline void Operator_Plus(const Eigen::VectorXd &r_v_delta)
	{
		Operator_Plus_Derived(r_v_delta);
	}

	/**
	 *	@brief subtracts delta vector to vertex state and minds the manifold space ("smart" minus)
	 *	@param[in] r_v_delta is delta vector
	 */
	inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
	{
		Operator_Minus_Derived(r_v_delta);
	}

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

	/**
	 *	@brief adds an edge that references this vertex
	 *
	 *	@param[in] p_edge is an edge that references this vertex (must be called
	 *		only once with each edge to avoid duplicates)
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline void Add_ReferencingEdge(_TyBaseEdge *p_edge) // throw(std::bad_alloc)
	{
		_ASSERTE(!b_IsReferencingEdge(p_edge));
		m_edge_list.push_back(p_edge);
	}

	/**
	 *	@brief checks if an edge references this vertex
	 *	@param[in] p_edge is pointer to an edge
	 *	@return Returns true if the given edge is in the list
	 *		of referencing edges, otherwise returns false.
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline bool b_IsReferencingEdge(const _TyBaseEdge *p_edge)
	{
		return std::find(m_edge_list.begin(), m_edge_list.end(), p_edge) != m_edge_list.end();
	}

	/**
	 *	@brief gets list of referencing edges
	 *	@return Returns const reference to the list of referencing edges.
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline const std::vector<_TyBaseEdge*> &r_ReferencingEdge_List() const
	{
		return m_edge_list;
	}

	/**
	 *	@brief allocates hessian block matrices
	 *	@param[out] r_lambda is the target matrix where the blocks are stored
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline void Alloc_HessianBlocks(CUberBlockMatrix &r_lambda)
	{
		Alloc_HessianBlocks_Derived(r_lambda);
	}

	/**
	 *	@brief calculates diagonal lambda hessian blocks
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline void Calculate_Hessians()
	{
		Calculate_Hessians_Derived();
	}

	/**
	 *	@brief gets a portion of right-hand side vector, associated with this vertex
	 *	@param[out] r_v_eta is the right-hand side vector
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	inline void Get_RightHandSide_Vector(Eigen::VectorXd &r_v_eta)
	{
		Get_RightHandSide_Vector_Derived(r_v_eta);
	}

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	/**
	 *	@brief allocates L factor block matrices
	 *	@param[out] r_L is the target matrix where the blocks are stored
	 *	@note This function is required for CNonlinearSolver_L.
	 */
	inline void Alloc_LBlocks(CUberBlockMatrix &r_L) const
	{
		Alloc_LBlocks_Derived(r_L);
	}

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---

protected:
	/**
	 *	@copydoc v_State()
	 */
	virtual Eigen::VectorXd v_State_Derived() const = 0;

	/**
	 *	@copydoc SwapState()
	 */
	virtual void SwapState_Derived(Eigen::VectorXd &r_v_x) = 0;

	/**
	 *	@copydoc SaveState()
	 */
	virtual void SaveState_Derived(Eigen::VectorXd &r_v_x) const = 0;

	/**
	 *	@copydoc LoadState()
	 */
	virtual void LoadState_Derived(const Eigen::VectorXd &r_v_x) = 0;

	/**
	 *	@copydoc Operator_Plus()
	 */
	virtual void Operator_Plus_Derived(const Eigen::VectorXd &r_v_delta) = 0;

	/**
	 *	@copydoc Operator_Minus()
	 */
	virtual void Operator_Minus_Derived(const Eigen::VectorXd &r_v_delta) = 0;

	/**
	 *	@copydoc CSEBaseVertex::n_Dimension()
	 */
	virtual size_t n_Dimension_Derived() const = 0;

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

	/**
	 *	@copydoc Alloc_HessianBlocks()
	 */
	virtual void Alloc_HessianBlocks_Derived(CUberBlockMatrix &r_lambda) = 0;

	/**
	 *	@copydoc Calculate_Hessians()
	 */
	virtual void Calculate_Hessians_Derived() = 0;

	/**
	 *	@copydoc Get_RightHandSide_Vector()
	 */
	virtual void Get_RightHandSide_Vector_Derived(Eigen::VectorXd &r_v_eta) = 0;

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	/**
	 *	@copydoc Alloc_LBlocks()
	 */
	virtual void Alloc_LBlocks_Derived(CUberBlockMatrix &r_L) const = 0;

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---

	// virtual functions call inline implementations in derived types (note that this is just
	// a performance optimization; virtual functions can be used without any limitations)
	// also note that MSVC can inline calls to even virtual functions, if resolvable (can it?)
};

/**
 *	@brief implementation of solver required functions for a generic SE vertex type
 *
 *	@tparam CDerivedVertex is name of the derived vertex class
 *	@tparam _n_dimension is vertex dimension (number of DOFs)
 */
template <class CDerivedVertex, int _n_dimension>
class CSEBaseVertexImpl : public CSEBaseVertex {
public:
	/**
	 *	@brief copy of template parameters
	 */
	enum {
		n_dimension = _n_dimension /**< @brief vertex dimension (number of DOFs) */
	};

	typedef Eigen::Matrix<double, n_dimension, 1> _TyVector; /**< @brief state vector type */
	typedef Eigen::Matrix<double, n_dimension, n_dimension> _TyMatrix; /**< @brief a compatible matrix type */

#ifdef __SE2_TYPES_USE_ALIGNED_VECTORS
	typedef Eigen::Matrix<double, n_dimension, 1> _TyVectorAlign; /**< @brief state vector storage type */
	typedef Eigen::Matrix<double, n_dimension, n_dimension> _TyMatrixAlign; /**< @brief a compatible matrix type */
#else // __SE2_TYPES_USE_ALIGNED_VECTORS
	typedef Eigen::Matrix<double, n_dimension, 1, Eigen::DontAlign> _TyVectorAlign; /**< @brief state vector storage type (have to use unaligned version due to the use in pools) */
	typedef Eigen::Matrix<double, n_dimension, n_dimension, Eigen::DontAlign> _TyMatrixAlign; /**< @brief a compatible matrix type */
#endif // __SE2_TYPES_USE_ALIGNED_VECTORS

protected:
	_TyVectorAlign m_v_state; /**< @brief state vector */

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific members ---

	double *m_p_HtSiH; /**< @brief pointer to diagonal block in the hessian matrix */
	_TyVectorAlign m_v_right_hand_side; /**< @brief right-hand side vector */

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific members ---

public:
	__SE2_TYPES_ALIGN_OPERATOR_NEW

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CSEBaseVertexImpl()
	{}

	/**
	 *	@brief constructor; initializes state vector
	 *	@param[in] r_v_state is state vector initializer
	 */
	inline CSEBaseVertexImpl(const _TyVector &r_v_state)
		:m_v_state(r_v_state)
	{}

	/**
	 *	@brief gets vertex state vector
	 *	@return Returns reference to vertex state vector.
	 */
	inline _TyVectorAlign &r_v_State()
	{
		return m_v_state;
	}

	/**
	 *	@brief gets vertex state vector
	 *	@return Returns const reference to vertex state vector.
	 */
	const inline _TyVectorAlign &v_State() const
	{
		return m_v_state;
	}

	/**
	 *	@copydoc CSEBaseVertex::n_Dimension()
	 */
	inline size_t n_Dimension() const
	{
		return n_dimension;
	}

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

	inline void Alloc_HessianBlocks(CUberBlockMatrix &r_lambda)
	{
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
		m_p_HtSiH = r_lambda.p_GetBlock_Log(m_n_id, m_n_id, n_dimension, n_dimension, true, false);
#else // __BASE_TYPES_USE_ID_ADDRESSING
		m_p_HtSiH = r_lambda.p_FindBlock(m_n_order, m_n_order, n_dimension, n_dimension, true, false);
#endif // __BASE_TYPES_USE_ID_ADDRESSING
		// find a block for hessian on the diagonal (edges can't do it, would have conflicts)

		_ASSERTE(m_p_HtSiH);
		// if this triggers, most likely __BASE_TYPES_USE_ID_ADDRESSING is enabled (see FlatSystem.h)
		// and the edges come in some random order
	}

	inline void Calculate_Hessians()
	{
		Eigen::Map<_TyMatrix, CUberBlockMatrix::map_Alignment> t_HtSiH(m_p_HtSiH); // t_odo - support aligned if the uberblockmatrix is aligned!
		t_HtSiH.setZero(); // !!
		m_v_right_hand_side.setZero(); // !!
		_ASSERTE(!m_edge_list.empty());
		for(size_t i = 0, n = m_edge_list.size(); i < n; ++ i) {
			std::pair<const double*, const double*> t_contrib = m_edge_list[i]->t_Get_LambdaEta_Contribution(this);
			Eigen::Map<_TyMatrix> t_lambda_contrib((double*)t_contrib.first); // todo - support aligned if the member matrices are aligned!
			t_HtSiH += t_lambda_contrib;
			Eigen::Map<_TyVector> v_eta_contrib((double*)t_contrib.second);  // todo - support aligned if the member matrices are aligned!
			m_v_right_hand_side += v_eta_contrib;
		}
		// calculate running sum (can't parallelize, in 10k or 100k there is
		// only about 6 edges / vertex, and that is the top density you get)
	}

	inline void Get_RightHandSide_Vector(Eigen::VectorXd &r_v_eta)
	{
		r_v_eta.segment<n_dimension>(m_n_order) = m_v_right_hand_side;
	}

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	inline void Alloc_LBlocks(CUberBlockMatrix &r_L) const
	{
		double *UNUSED(p_block_addr);
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
		p_block_addr = r_L.p_GetBlock_Log(m_n_id, m_n_id, n_dimension, n_dimension, true, true); // don't care about the pointer, it is handled differently
#else // __BASE_TYPES_USE_ID_ADDRESSING
		p_block_addr = r_L.p_FindBlock(m_n_order, m_n_order, n_dimension, n_dimension, true, true); // don't care about the pointer, it is handled differently
#endif // __BASE_TYPES_USE_ID_ADDRESSING
		// find a block for hessian on the diagonal (edges can't do it, would have conflicts) // also clear it to zero! L is updated additively

		_ASSERTE(p_block_addr);
		// if this triggers, most likely __BASE_TYPES_USE_ID_ADDRESSING is enabled (see FlatSystem.h)
		// and the edges come in some random order
	}

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---

	/**
	 *	@copydoc CSEBaseVertex::SwapState()
	 */
	inline void SwapState(Eigen::VectorXd &r_v_x)
	{
		m_v_state.swap(r_v_x.segment<n_dimension>(m_n_order)); // swap the numbers
	}

	/**
	 *	@copydoc CSEBaseVertex::SaveState()
	 */
	inline void SaveState(Eigen::VectorXd &r_v_x) const
	{
		r_v_x.segment<n_dimension>(m_n_order) = m_v_state; // save the state vector
	}

	/**
	 *	@copydoc CSEBaseVertex::LoadState()
	 */
	inline void LoadState(const Eigen::VectorXd &r_v_x)
	{
		m_v_state = r_v_x.segment<n_dimension>(m_n_order); // restore the state vector
	}

protected:
	virtual Eigen::VectorXd v_State_Derived() const
	{
		return CSEBaseVertexImpl<CDerivedVertex, _n_dimension>::v_State();
	}

	virtual void SwapState_Derived(Eigen::VectorXd &r_v_x)
	{
		((CDerivedVertex*)this)->SwapState(r_v_x);
	}

	virtual void SaveState_Derived(Eigen::VectorXd &r_v_x) const
	{
		((CDerivedVertex*)this)->SaveState(r_v_x);
	}

	virtual void LoadState_Derived(const Eigen::VectorXd &r_v_x)
	{
		((CDerivedVertex*)this)->LoadState(r_v_x);
	}

	virtual void Operator_Plus_Derived(const Eigen::VectorXd &r_v_delta)
	{
		((CDerivedVertex*)this)->Operator_Plus(r_v_delta);
	}

	virtual void Operator_Minus_Derived(const Eigen::VectorXd &r_v_delta)
	{
		((CDerivedVertex*)this)->Operator_Minus(r_v_delta);
	}

	virtual size_t n_Dimension_Derived() const
	{
		return CSEBaseVertexImpl<CDerivedVertex, _n_dimension>::n_Dimension();
	}

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

	virtual void Alloc_HessianBlocks_Derived(CUberBlockMatrix &r_lambda)
	{
		CSEBaseVertexImpl<CDerivedVertex, _n_dimension>::Alloc_HessianBlocks(r_lambda);
	}

	virtual void Calculate_Hessians_Derived()
	{
		CSEBaseVertexImpl<CDerivedVertex, _n_dimension>::Calculate_Hessians();
	}

	virtual void Get_RightHandSide_Vector_Derived(Eigen::VectorXd &r_v_eta)
	{
		CSEBaseVertexImpl<CDerivedVertex, _n_dimension>::Get_RightHandSide_Vector(r_v_eta);
	}

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	virtual void Alloc_LBlocks_Derived(CUberBlockMatrix &r_L) const
	{
		CSEBaseVertexImpl<CDerivedVertex, _n_dimension>::Alloc_LBlocks(r_L);
	}

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---
};

/**
 *	@brief SE base edge implementation
 *
 *	@tparam CDerivedEdge is the name of derived edge class
 *	@tparam CVertex0 is the name of the first edge vertex class
 *	@tparam CVertex1 is the name of the second edge vertex class
 *	@tparam _n_measurement_dimension is measurement vector dimension
 */
template <class CDerivedEdge, class CVertex0, class CVertex1, int _n_measurement_dimension>
class CSEBaseEdgeImpl : public CSEBaseEdge {
public:
	/**
	 *	@brief copy of template parameters
	 */
	enum {
		n_vertex0_dimension = CVertex0::n_dimension, /**< @brief the first vertex dimension */
		n_vertex1_dimension = CVertex1::n_dimension, /**< @brief the second vertex dimension */
		n_measurement_dimension = _n_measurement_dimension /**< @brief measurement vector dimension (edge dimension) */
	};

	typedef CVertex0 _TyVertex0; /**< @brief name of the first edge vertex class */
	typedef CVertex1 _TyVertex1; /**< @brief name of the second edge vertex class */

	/**
	 *	@brief vertex initialization functor
	 *	Just clears the vertex to null.
	 */
	class CInitializeNullVertex {
	public:
		/**
		 *	@brief function operator
		 *	@return Returns the value of the vertex being initialized.
		 */
		inline operator CVertex0() const
		{
			typename CVertex0::_TyVector v_null;
			v_null.setZero();
			return CVertex0(v_null);
		}
	};

	/**
	 *	@brief vertex initialization functor
	 *	Just clears the vertex to null.
	 */
	class CInitializeNullVertex1 {
	public:
		/**
		 *	@brief function operator
		 *	@return Returns the value of the vertex being initialized.
		 */
		inline operator CVertex1() const
		{
			typename CVertex1::_TyVector v_null;
			v_null.setZero();
			return CVertex1(v_null);
		}
	};

	typedef Eigen::Matrix<double, n_measurement_dimension, 1> _TyVector; /**< @brief edge dimension vector type */
	typedef Eigen::Matrix<double, n_measurement_dimension, n_measurement_dimension> _TyMatrix; /**< @brief edge dimension matrix type */
	typedef Eigen::Matrix<double, n_measurement_dimension, n_vertex0_dimension> _TyMatrix0; /**< @brief mixed edge / the first vertex dimension matrix type */
	typedef Eigen::Matrix<double, n_measurement_dimension, n_vertex1_dimension> _TyMatrix1; /**< @brief mixed edge / the second vertex dimension matrix type */

protected:
#ifdef __SE2_TYPES_USE_ALIGNED_VECTORS
	typedef Eigen::Matrix<double, n_measurement_dimension, 1> _TyVectorAlign; /**< @brief edge dimension vector type with member alignment */
	typedef Eigen::Matrix<double, n_measurement_dimension, n_measurement_dimension> _TyMatrixAlign; /**< @brief edge dimension matrix type with member alignment */ // note this is dirty, it is still called _Unalign, but it is aligned
#else // __SE2_TYPES_USE_ALIGNED_VECTORS
	typedef Eigen::Matrix<double, n_measurement_dimension, 1, Eigen::DontAlign> _TyVectorAlign; /**< @brief edge dimension vector type with member alignment */
	typedef Eigen::Matrix<double, n_measurement_dimension, n_measurement_dimension, Eigen::DontAlign> _TyMatrixAlign; /**< @brief edge dimension matrix type with member alignment */
#endif // __SE2_TYPES_USE_ALIGNED_VECTORS

	CVertex0 *m_p_vertex0; /**< @brief pointers to the referenced vertices */
	CVertex1 *m_p_vertex1; /**< @brief pointers to the referenced vertices */
	_TyVectorAlign m_v_measurement; /**< @brief the measurement */
	_TyMatrixAlign m_t_sigma_inv; /**< @brief information matrix */

#ifdef __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific members ---

	_TyMatrixAlign m_t_square_root_sigma_inv_upper; /**< @brief the R matrix (upper diagonal) = transpose(chol(t_sigma_inv)) */
	_TyVectorAlign m_v_error; /**< @brief error vector (needs to be filled explicitly by calling p_error_function) */
	//_TyVectorAlign m_v_expectation; /**< @brief expectation vector (needs to be filled explicitly by calling p_error_function) */ // no need to store this
	double *m_p_RH[2]; /**< @brief blocks of memory where the R * H matrices are stored */

#endif // __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific members ---

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific members ---

	double *m_p_HtSiH; /**< @brief block of memory where Ht * inv(Sigma) * H matrix is stored (the one above diagonal) */
	double *m_p_HtSiH_reduce; /**< @brief block of memory where Ht * inv(Sigma) * H matrix is stored in the matrix (should there be conflicts between duplicate edges, this is the target of reduction) */
	bool m_b_first_reductor; /**< @brief reduction flag */
	typename _TyVertex0::_TyMatrixAlign m_t_HtSiH_vertex0; /**< @brief block of memory where Ht * inv(Sigma) * H matrix is stored (the vertex contribution) */
	typename _TyVertex1::_TyMatrixAlign m_t_HtSiH_vertex1; /**< @brief block of memory where Ht * inv(Sigma) * H matrix is stored (the vertex contribution) */
	typename _TyVertex0::_TyVectorAlign m_t_right_hand_vertex0; /**< @brief block of memory where right-hand side vector is stored (the vertex contribution) */
	typename _TyVertex1::_TyVectorAlign m_t_right_hand_vertex1; /**< @brief block of memory where right-hand side vector is stored (the vertex contribution) */

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific members ---

public:
	__SE2_TYPES_ALIGN_OPERATOR_NEW

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CSEBaseEdgeImpl()
	{}

	/**
	 *	@brief constructor
	 *
	 *	@param[in] n_vertex0_id is id of the first vertex
	 *	@param[in] n_vertex1_id is id of the second vertex
	 *	@param[in] r_v_measurement is the measurement vecot
	 *	@param[in] r_t_sigma_inv is inverse sigma matrix
	 *
	 *	@note This fills the structure, except for the vertex pointers.
	 *	@note Any derived class should check for vertex dimensions
	 *		being really what it expects (not checked elsewhere).
	 */
	inline CSEBaseEdgeImpl(size_t n_vertex0_id, size_t n_vertex1_id,
		const _TyVector &r_v_measurement, const _TyMatrix &r_t_sigma_inv)
		:m_v_measurement(r_v_measurement), m_t_sigma_inv(r_t_sigma_inv)
#ifdef __SE_TYPES_SUPPORT_A_SOLVERS
		, m_t_square_root_sigma_inv_upper(r_t_sigma_inv.llt().matrixU()) // calculate R
#endif // __SE_TYPES_SUPPORT_A_SOLVERS
	{
		m_p_vertex_id[0] = n_vertex0_id;
		m_p_vertex_id[1] = n_vertex1_id;
		// get references to the vertices, initialize the vertices, if neccessary

		// note that errors, expectations and jacobian matrices are not cleared
	}

	/**
	 *	@brief updates the edge with new measurement
	 *
	 *	@param[in] r_v_measurement is the measurement vecot
	 *	@param[in] r_t_sigma_inv is inverse sigma matrix
	 */
	inline void Update(const _TyVector &r_v_measurement, const _TyMatrix &r_t_sigma_inv)
	{
		m_v_measurement = r_v_measurement;
		m_t_sigma_inv = r_t_sigma_inv;
#ifdef __SE_TYPES_SUPPORT_A_SOLVERS
		m_t_square_root_sigma_inv_upper = r_t_sigma_inv.llt().matrixU(); // calculate R
#endif // __SE_TYPES_SUPPORT_A_SOLVERS
	}

	/**
	 *	@copydoc CSEBaseEdge::n_Dimension()
	 */
	inline size_t n_Dimension() const
	{
		return n_measurement_dimension;
	}

	/**
	 *	@brief gets inverse sigma of the measurement
	 *	@return Returns const reference to the inverse sigma matrix.
	 */
	inline const _TyMatrixAlign &t_Sigma_Inv() const
	{
		return m_t_sigma_inv;
	}

#ifdef __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific functions ---

	inline void Alloc_JacobianBlocks(CUberBlockMatrix &r_A)
	{
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
		m_p_RH[0] = r_A.p_GetBlock_Log(m_n_id + 1, m_p_vertex_id[0],
			n_measurement_dimension, n_vertex0_dimension, true, false);
		m_p_RH[1] = r_A.p_GetBlock_Log(m_n_id + 1, m_p_vertex_id[1],
			n_measurement_dimension, n_vertex1_dimension, true, false);
#else // __BASE_TYPES_USE_ID_ADDRESSING
		m_p_RH[0] = r_A.p_FindBlock(m_n_order, m_p_vertex0->n_Order(),
			n_measurement_dimension, n_vertex0_dimension, true, false);
		m_p_RH[1] = r_A.p_FindBlock(m_n_order, m_p_vertex1->n_Order(),
			n_measurement_dimension, n_vertex1_dimension, true, false);
#endif // __BASE_TYPES_USE_ID_ADDRESSING
		// just place blocks at (edge order, vertex order)

		_ASSERTE(m_p_RH[0]);
		_ASSERTE(m_p_RH[1]);
		// if this triggers, most likely __BASE_TYPES_USE_ID_ADDRESSING is enabled (see FlatSystem.h)
		// and the edges come in some random order
	}

	inline void Calculate_Jacobians()
	{
		_TyMatrix0 t_jacobian0;
		_TyMatrix1 t_jacobian1;
		_TyVector v_expectation, v_error;
		((CDerivedEdge*)this)->Calculate_Jacobians_Expectation_Error(t_jacobian0,
			t_jacobian1, v_expectation, v_error);

		m_v_error = v_error;
		// calculates the expectation, error and the jacobians (implemented by the edge)

		Eigen::Map<_TyMatrix0, CUberBlockMatrix::map_Alignment> t_RH0(m_p_RH[0]);
		Eigen::Map<_TyMatrix1, CUberBlockMatrix::map_Alignment> t_RH1(m_p_RH[1]); // t_odo - support aligned if the uberblockmatrix is aligned!
		// map hessian matrices

		t_RH0 = m_t_square_root_sigma_inv_upper * t_jacobian0;
		t_RH1 = m_t_square_root_sigma_inv_upper * t_jacobian1;
		// recalculate RH (transpose cholesky of sigma times the jacobian)
		// note this references the A block matrix
	}

	inline void Get_R_Error(Eigen::VectorXd &r_v_dest) const
	{
		r_v_dest.segment<n_measurement_dimension>(m_n_order) = m_t_square_root_sigma_inv_upper * m_v_error;
	}

	inline void Get_Error(Eigen::VectorXd &r_v_dest) const
	{
		r_v_dest.segment<n_measurement_dimension>(m_n_order) = m_v_error;
	}

#endif // __SE_TYPES_SUPPORT_A_SOLVERS // --- ~A-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- lambda-SLAM specific functions ---

	virtual inline bool Notify_HessianBlock_Conflict(double *p_block, CUberBlockMatrix &r_lambda)
	{
		if(m_p_HtSiH == p_block) {
			if(!m_p_HtSiH_reduce) {
				m_p_HtSiH = r_lambda.p_Get_DenseStorage(n_vertex0_dimension * n_vertex1_dimension);
				// get a temporary block

				memcpy(m_p_HtSiH, p_block, n_vertex0_dimension * n_vertex1_dimension * sizeof(double));
				// this edge's hessian might already have significant value, need to store it
			}
			// first time arround this is called, alloc a different block to store this edge's hessian

			m_p_HtSiH_reduce = p_block;
			m_b_first_reductor = true; // i was the first

			return true; // acknowledge notify
		}
		return false; // not me
	}

	inline void Alloc_HessianBlocks(CUberBlockMatrix &r_lambda)
	{
		m_p_vertex0->Add_ReferencingEdge(this);
		m_p_vertex1->Add_ReferencingEdge(this);
		// the vertices need to know this edge in order to calculate sum of jacobians on the diagonal of lambda // t_odo - should the lambda solver call this instead? might avoid duplicate records in the future // lambda solver is calling this only once in it's lifetime, no duplicates will occur

		_ASSERTE((m_p_vertex_id[0] > m_p_vertex_id[1]) ==
			(m_p_vertex0->n_Order() > m_p_vertex1->n_Order()));
		// if this triggers, then the edge has the vertices assigned in a different
		// order than the ids (vertex[0] is id[1] and vice versa) consequently, the hessians
		// will have correct shape in the matrix, but the data will be transposed
		// and you will get either not pos def / rubbish solutions

		bool b_uninit_block;
		size_t n_dimension0 = n_vertex0_dimension;
		size_t n_dimension1 = n_vertex1_dimension;
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
		size_t n_id_0 = m_p_vertex_id[0];
		size_t n_id_1 = m_p_vertex_id[1];
		if(n_id_0 > n_id_1) {
			std::swap(n_id_0, n_id_1);
			std::swap(n_dimension0, n_dimension1);
		}
		// make sure the order is sorted (if swapping, will have to transpose the result,
		// but we will deal with that laters)

		_ASSERTE(n_id_0 != n_id_1); // will otherwise overwrite blocks with vertex blocks (malformed system)
		_ASSERTE(n_id_0 < n_id_1);
		// they dont care about the man that ends up under the stage, they only care about the one above (column > row)

		m_p_HtSiH = r_lambda.p_GetBlock_Log_Alloc(n_id_0, n_id_1, n_dimension0, n_dimension1, b_uninit_block);
#else // __BASE_TYPES_USE_ID_ADDRESSING
		size_t n_order_0 = m_p_vertex0->n_Order();
		size_t n_order_1 = m_p_vertex1->n_Order();
		if(n_order_0 > n_order_1) {
			std::swap(n_order_0, n_order_1);
			std::swap(n_dimension0, n_dimension1);
		}
		// make sure the order is sorted (if swapping, will have to transpose the result,
		// but we will deal with that laters)

		_ASSERTE(n_order_0 != n_order_1); // will otherwise overwrite blocks with vertex blocks (malformed system)
		_ASSERTE(n_order_0 < n_order_1);
		// they dont care about the man that ends up under the stage, they only care about the one above (column > row)

		m_p_HtSiH = r_lambda.p_FindBlock_Alloc(n_order_0, n_order_1, n_dimension0, n_dimension1, b_uninit_block);
#endif // __BASE_TYPES_USE_ID_ADDRESSING
		// find a block for hessian above the diagonal, and with the right shape

		// t_odo - need a new function that finds a block, and says if it was there before, or not.
		// that could break, if the lambda was already allocated. does that ever happen?
		// note that m_p_vertex[]->Add_ReferencingEdge() is called here as well, and it is guaranteed
		// to only be called once. that brings some guarantees ...

		_ASSERTE(m_p_HtSiH);
		// if this triggers, most likely __BASE_TYPES_USE_ID_ADDRESSING is enabled (see FlatSystem.h)
		// and the edges come in some random order

		bool b_duplicate_block = !b_uninit_block; // was it there already?
		if(b_duplicate_block) {
			const std::vector<typename _TyVertex1::_TyBaseEdge*> &r_list =
				m_p_vertex1->r_ReferencingEdge_List(); // the second vertex has presumably less references, making the search faster; could select the smaller one manually
			for(size_t i = 0, n = r_list.size(); i < n; ++ i) {
				if(r_list[i]->Notify_HessianBlock_Conflict(m_p_HtSiH, r_lambda))
					break;
				_ASSERTE(i + 1 != n); // make sure we found it
			}
			// find the first reductor!

			m_b_first_reductor = false; // this is duplicate, there already was another
			m_p_HtSiH_reduce = m_p_HtSiH; // where to reduce
			m_p_HtSiH = r_lambda.p_Get_DenseStorage(n_dimension0 * n_dimension1); // a temporary
			// this edge is a duplicate edge; reduction of diagonal blocks is taken care of in the vertices.
			// reduction of off-diagonal blocks needs to be taken care of explicitly

			// need a new storage for the block data, so that Calculate_Hessians() has where to put it
			// could put it in the lambda matrix pool, without putting it in any column (a bit dirty, but ok)

			// finally, the reduction needs to be carried out somehow

			// each duplicate edge between two vertices can be processed by either of the vertices,
			// in CSEBaseVertexImpl::Calculate_Hessians(). there are some problems with load balancing, but
			// it should generally work. there is a problem with matrix size, as the vertices do not have
			// to be the same size, and each vertex does not know about the other vertex' dimension

			// it needs to be done in the edge, but potentially invoked by the vertex
			// to avoid conflicts, the reduction should be done by vertex with e.g. the larger id
			// this will require one more pointer in edge (for the orig block) and one more list in the vertex
			// otherwise, could have one more field in vertex, which would say where to reduce. that saves space
			// as those lists would be normally empty

			// what will happen in N-ary edges?
		} else {
			//m_b_first_reductor = false; // leave uninitialized, don't care
			m_p_HtSiH_reduce = 0; // do not reduce
		}
	}

	inline void Calculate_Hessians()
	{
		_ASSERTE(m_p_vertex0->b_IsReferencingEdge(this));
		_ASSERTE(m_p_vertex1->b_IsReferencingEdge(this));
		// makes sure this is amongst the referencing edges

		_TyMatrix0 t_jacobian0;
		_TyMatrix1 t_jacobian1;
		_TyVector v_expectation, v_error;
		((CDerivedEdge*)this)->Calculate_Jacobians_Expectation_Error(t_jacobian0,
			t_jacobian1, v_expectation, v_error);
		// calculates the expectation and the jacobians

		bool b_transpose_hessian;
		size_t n_dimension0 = n_vertex0_dimension;
		size_t n_dimension1 = n_vertex1_dimension;
#if 1
		size_t n_id_0 = m_p_vertex_id[0];
		size_t n_id_1 = m_p_vertex_id[1]; // this is closer in cache
		_ASSERTE((n_id_0 > n_id_1) == (m_p_vertex0->n_Order() > m_p_vertex1->n_Order())); // if this triggers, then the edge has the vertices assigned in different order than the ids (vertex[0] is id[1] and vice versa)
		if((b_transpose_hessian = (n_id_0 > n_id_1))) {
			std::swap(n_id_0, n_id_1);
			std::swap(n_dimension0, n_dimension1);
		}
		// make sure the order is sorted (if swapping, will have to transpose the result,
		// but we will deal with that laters)

		_ASSERTE(n_id_0 != n_id_1); // will otherwise overwrite blocks with vertex blocks (malformed system)
		_ASSERTE(n_id_0 < n_id_1);
		// they dont care about the man that ends up under the stage, they only care about the one above (column > row)
#else // 1
		size_t n_order_0 = m_p_vertex0->n_Order();
		size_t n_order_1 = m_p_vertex1->n_Order();
		if((b_transpose_hessian = (n_order_0 > n_order_1))) {
			std::swap(n_order_0, n_order_1);
			std::swap(n_dimension0, n_dimension1);
		}
		// make sure the order is sorted (if swapping, will have to transpose the result,
		// but we will deal with that laters)

		_ASSERTE(n_order_0 != n_order_1); // will otherwise overwrite blocks with vertex blocks (malformed system)
		_ASSERTE(n_order_0 < n_order_1);
#endif // 1
		// they dont care about the man that ends up under the stage, they only care about the one above (column > row)

		Eigen::Matrix<double, n_vertex0_dimension, n_measurement_dimension> t_H0_sigma_inv =
			t_jacobian0.transpose() * m_t_sigma_inv;

		if(b_transpose_hessian) {
			Eigen::Map<Eigen::Matrix<double, n_vertex1_dimension, n_vertex0_dimension>,
				CUberBlockMatrix::map_Alignment> t_HtSiH(m_p_HtSiH); // t_odo - support aligned if the uberblockmatrix is aligned!
			// map the matrix above diagonal

			t_HtSiH.noalias() = t_jacobian1.transpose() * t_H0_sigma_inv.transpose();
		} else {
			Eigen::Map<Eigen::Matrix<double, n_vertex0_dimension, n_vertex1_dimension>,
				CUberBlockMatrix::map_Alignment> t_HtSiH(m_p_HtSiH); // t_odo - support aligned if the uberblockmatrix is aligned!
			// map the matrix above diagonal

			t_HtSiH.noalias() = t_H0_sigma_inv * t_jacobian1;
		}
		// calculate the matrix above diagonal

		m_t_HtSiH_vertex0.noalias() = t_H0_sigma_inv * t_jacobian0;
		m_t_HtSiH_vertex1.noalias() = t_jacobian1.transpose() * m_t_sigma_inv * t_jacobian1;
		// calculate vertex hessian contributions

		m_t_right_hand_vertex0.noalias() = t_jacobian0.transpose() * (m_t_sigma_inv * v_error);
		m_t_right_hand_vertex1.noalias() = t_jacobian1.transpose() * (m_t_sigma_inv * v_error);
		// calculate right hand side vector contributions
	}

	virtual double f_Max_VertexHessianDiagValue() const // todo surround with levenberg support ifdef
	{
		double f_max = 0;
		for(size_t i = 0; i < n_vertex0_dimension; ++ i)
			f_max = std::max(f_max, m_t_HtSiH_vertex0(i, i));
		for(size_t i = 0; i < n_vertex1_dimension; ++ i)
			f_max = std::max(f_max, m_t_HtSiH_vertex1(i, i));
		return f_max;
	}

	inline std::pair<const double*, const double*> t_Get_LambdaEta_Contribution(const CSEBaseVertex *p_which)
	{
		if(m_p_HtSiH_reduce) {
			if(m_b_first_reductor) {
				Eigen::Map<Eigen::Matrix<double, n_vertex0_dimension, n_vertex1_dimension>,
					CUberBlockMatrix::map_Alignment> t_HtSiH_reduce(m_p_HtSiH_reduce);
				t_HtSiH_reduce.setZero();
			}
			// !!

			size_t n_id_0 = m_p_vertex_id[0], n_id_1 = m_p_vertex_id[1];
			const CSEBaseVertex *p_max_vert = (n_id_0 > n_id_1)?
				(const CSEBaseVertex*)m_p_vertex0 : (const CSEBaseVertex*)m_p_vertex1;
			if(p_max_vert == p_which) {
				/*if(n_id_0 > n_id_1) { // note that it does not matter, it is just elementwise add
					Eigen::Map<Eigen::Matrix<double, n_vertex1_dimension, n_vertex0_dimension>,
						CUberBlockMatrix::map_Alignment> t_HtSiH(m_p_HtSiH), t_HtSiH_reduce(m_p_HtSiH_reduce);
					t_HtSiH_reduce.noalias() += t_HtSiH;
				} else*/ {
					Eigen::Map<Eigen::Matrix<double, n_vertex0_dimension, n_vertex1_dimension>,
						CUberBlockMatrix::map_Alignment> t_HtSiH(m_p_HtSiH), t_HtSiH_reduce(m_p_HtSiH_reduce);
					t_HtSiH_reduce.noalias() += t_HtSiH;
				}
			}
		}
		// take care of hessian reduction, kind of silently, but it is guaranteed
		// that this is called for all the updated vertices

		_ASSERTE(p_which == m_p_vertex0 || p_which == m_p_vertex1);
		return (p_which == m_p_vertex0)?
			std::make_pair(m_t_HtSiH_vertex0.data(), m_t_right_hand_vertex0.data()) :
			std::make_pair(m_t_HtSiH_vertex1.data(), m_t_right_hand_vertex1.data());
		// return pointer to the matrix data with hessian contribution
	}

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	inline void Alloc_LBlocks(CUberBlockMatrix &UNUSED(r_L)) const
	{
		// no vertex refs, that is already done in Alloc_HessianBlocks()
#if 0
		size_t n_dimension0 = n_vertex0_dimension;
		size_t n_dimension1 = n_vertex1_dimension;
		size_t n_order_0 = m_p_vertex0->n_Order();
		size_t n_order_1 = m_p_vertex1->n_Order();
		if(n_order_0 > n_order_1) {
			std::swap(n_order_0, n_order_1);
			std::swap(n_dimension0, n_dimension1);
		}
		// make sure the order is sorted (if swapping, will have to transpose the result,
		// but we will deal with that laters)

		_ASSERTE(n_order_0 != n_order_1); // will otherwise overwrite blocks with vertex blocks (malformed system)
		_ASSERTE(n_order_0 < n_order_1);
		// they dont care about the man that ends up under the stage, they only care about the one above (column > row)

		r_L.p_FindBlock(n_order_0, n_order_1, n_dimension0, n_dimension1, true, true); // pointer not remembered, will handle it differently
		// find a block for hessian above the diagonal, and with the right shape // also set it to zero, L is updated additively
#endif // 0
		// edges should not allocate off-diagonal blocks, it will just screw the reordering and increase fill-in
	}

	/**
	 *	@brief calculates just a part of the lambda matrix, called omega
	 *
	 *	This is used in CNonlinearSolver_L for incremental updates to the factor.
	 *	It calculates omega from scratch, uses serial execution, as opposed
	 *	to Calculate_Hessians(). While it may seem slower, it can't likely be parallelized
	 *	due to very small size of omega (can be just a single edge).
	 *
	 *	@param[out] r_omega is the omega matrix to be filled (must be initially empty)
	 *	@param[in] n_min_vertex_order is order offset, in elements, used to position
	 *		the hessian blocks in the upper-left corner of omega
	 */
	inline void Calculate_Omega(CUberBlockMatrix &r_omega, size_t n_min_vertex_order) const // todo - can reuse the calculations done for lambda, the hessians are all already there, ready to be used (add a flag?)
	{
		//_ASSERTE(r_omega.b_Empty()); // should be empty // can't assert that here, other edges might have added something already

		bool b_transpose_hessian;
		size_t n_dimension0 = n_vertex0_dimension;
		size_t n_dimension1 = n_vertex1_dimension;

#if 0 // this does not work, omega tends to be very sparse, might need to create some empty columns / rows and p_GetBlock_Log() does not support that
		size_t n_id_0 = m_p_vertex_id[0];
		size_t n_id_1 = m_p_vertex_id[1];
		if((b_transpose_hessian = (n_id_0 > n_id_1))) {
			std::swap(n_id_0, n_id_1);
			std::swap(n_dimension0, n_dimension1);
		}
		// make sure the order is sorted (if swapping, will have to transpose the result,
		// but we will deal with that laters)

		_ASSERTE(n_id_0 != n_id_1); // will otherwise overwrite blocks with vertex blocks (malformed system)
		_ASSERTE(n_id_0 < n_id_1);
		// they dont care about the man that ends up under the stage, they only care about the one above (column > row)

		if(n_id_0 < n_min_vertex_id) // note this might be superficial // todo - make this an assert
			return;
		// this edge doesn't have any hessians inside omega

		double *p_v0 = r_omega.p_GetBlock_Log(n_id_0 - n_min_vertex_id,
			n_id_0 - n_min_vertex_id, n_dimension0, n_dimension0, true, true);
		double *p_v1 = r_omega.p_GetBlock_Log(n_id_1 - n_min_vertex_id,
			n_id_1 - n_min_vertex_id, n_dimension1, n_dimension1, true, true);

		double *p_edge = r_omega.p_GetBlock_Log(n_id_0 - n_min_vertex_id,
			n_id_1 - n_min_vertex_id, n_dimension0, n_dimension1, true, false); // doesn't need to be initialized, there's only one
		// must come last, n_id_1 might not exist in r_omega before the vertices are added
		// this does not work anyway, omega tends to be very sparse, might need to create some empty columns / rows and p_GetBlock_Log() does not support that
#else // 0
		size_t n_order_0 = m_p_vertex0->n_Order();
		size_t n_order_1 = m_p_vertex1->n_Order();
		if((b_transpose_hessian = (n_order_0 > n_order_1))) {
			std::swap(n_order_0, n_order_1);
			std::swap(n_dimension0, n_dimension1);
		}
		// make sure the order is sorted (if swapping, will have to transpose the result,
		// but we will deal with that laters)

		_ASSERTE(n_order_0 != n_order_1); // will otherwise overwrite blocks with vertex blocks (malformed system)
		_ASSERTE(n_order_0 < n_order_1);
		// they dont care about the man that ends up under the stage, they only care about the one above (column > row)

		if(n_order_0 < n_min_vertex_order) // note this might be superficial // todo - make this an assert
			return;
		// this edge doesn't have any hessians inside omega

		double *p_edge = r_omega.p_FindBlock(n_order_0 - n_min_vertex_order,
			n_order_1 - n_min_vertex_order, n_dimension0, n_dimension1, true, false); // doesn't need to be initialized, there's only one
		double *p_v0 = r_omega.p_FindBlock(n_order_0 - n_min_vertex_order,
			n_order_0 - n_min_vertex_order, n_dimension0, n_dimension0, true, true);
		double *p_v1 = r_omega.p_FindBlock(n_order_1 - n_min_vertex_order,
			n_order_1 - n_min_vertex_order, n_dimension1, n_dimension1, true, true);
#endif // 0
		// alloc and initialize / find existing blocks for all the hessians, above the diagonal only

		bool b_recalculate = false;
		// recalculate from scratch, or use values already calculated by the lambda solver?

		_TyMatrix0 t_jacobian0;
		_TyMatrix1 t_jacobian1;
		if(b_recalculate) {
			_TyVector v_expectation, v_error;
			((CDerivedEdge*)this)->Calculate_Jacobians_Expectation_Error(t_jacobian0,
				t_jacobian1, v_expectation, v_error);
		}
		// calculates the expectation and the jacobians

		Eigen::Matrix<double, n_vertex0_dimension, n_measurement_dimension> t_H0_sigma_inv;
		if(b_recalculate)
			t_H0_sigma_inv = t_jacobian0.transpose() * m_t_sigma_inv;

		if(b_transpose_hessian) {
			Eigen::Map<Eigen::Matrix<double, n_vertex1_dimension, n_vertex0_dimension>,
				CUberBlockMatrix::map_Alignment> t_HtSiH(p_edge); // t_odo - support aligned if the uberblockmatrix is aligned!
			// map the matrix above diagonal

			if(b_recalculate)
				t_HtSiH.noalias() = t_jacobian1.transpose() * t_H0_sigma_inv.transpose();
			else {
				Eigen::Map<Eigen::Matrix<double, n_vertex1_dimension, n_vertex0_dimension>,
					CUberBlockMatrix::map_Alignment> t_HtSiH_lambda(m_p_HtSiH);
				t_HtSiH = t_HtSiH_lambda;
			}
		} else {
			Eigen::Map<Eigen::Matrix<double, n_vertex0_dimension, n_vertex1_dimension>,
				CUberBlockMatrix::map_Alignment> t_HtSiH(p_edge); // t_odo - support aligned if the uberblockmatrix is aligned!
			// map the matrix above diagonal

			if(b_recalculate)
				t_HtSiH.noalias() = t_H0_sigma_inv * t_jacobian1;
			else {
				Eigen::Map<Eigen::Matrix<double, n_vertex0_dimension, n_vertex1_dimension>,
					CUberBlockMatrix::map_Alignment> t_HtSiH_lambda(m_p_HtSiH);
				t_HtSiH = t_HtSiH_lambda;
			}
		}
		// calculate the matrix above diagonal

		//if(n_vertex0_dimension != n_vertex1_dimension)
		//	b_transpose_hessian = !b_transpose_hessian;
		// cause a bug deliberately to test the block bounds checks in block matrix

		Eigen::Map<Eigen::Matrix<double, n_vertex0_dimension, n_vertex0_dimension>,
			CUberBlockMatrix::map_Alignment> t_HtSiH_v0((b_transpose_hessian)? p_v1 : p_v0);
		Eigen::Map<Eigen::Matrix<double, n_vertex1_dimension, n_vertex1_dimension>,
			CUberBlockMatrix::map_Alignment> t_HtSiH_v1((b_transpose_hessian)? p_v0 : p_v1);
		if(b_recalculate) {
			t_HtSiH_v0.noalias() += t_H0_sigma_inv * t_jacobian0;
			t_HtSiH_v1.noalias() += t_jacobian1.transpose() * m_t_sigma_inv * t_jacobian1;
		} else {
			t_HtSiH_v0 += m_t_HtSiH_vertex0;
			t_HtSiH_v1 += m_t_HtSiH_vertex1;
		}
		// calculate vertex hessian contributions (note the increments)
	}

	// todo - will need error to incrementally update right hand side vector

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---

protected:
	virtual size_t n_Dimension_Derived() const
	{
		return CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1, _n_measurement_dimension>::n_Dimension();
	}

	virtual double f_Chi_Squared_Error_Derived() const
	{
		return ((const CDerivedEdge*)this)->f_Chi_Squared_Error();
	}

#ifdef __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific functions ---

	virtual void Alloc_JacobianBlocks_Derived(CUberBlockMatrix &r_A)
	{
		CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::Alloc_JacobianBlocks(r_A);
	}

	virtual void Calculate_Jacobians_Derived()
	{
		CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::Calculate_Jacobians();
	}

	virtual void Get_R_Error_Derived(Eigen::VectorXd &r_v_dest) const
	{
		CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::Get_R_Error(r_v_dest);
	}

	virtual void Get_Error_Derived(Eigen::VectorXd &r_v_dest) const
	{
		CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::Get_Error(r_v_dest);
	}

#endif // __SE_TYPES_SUPPORT_A_SOLVERS // --- ~A-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

	virtual bool Notify_HessianBlock_Conflict_Derived(double *p_block, CUberBlockMatrix &r_lambda)
	{
		return Notify_HessianBlock_Conflict(p_block, r_lambda);
	}

	virtual void Alloc_HessianBlocks_Derived(CUberBlockMatrix &r_lambda)
	{
		CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::Alloc_HessianBlocks(r_lambda);
	}

	virtual void Calculate_Hessians_Derived()
	{
		CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::Calculate_Hessians();
	}

	virtual std::pair<const double*, const double*>
		t_Get_LambdaEta_Contribution_Derived(const CSEBaseVertex *p_which)
	{
		return CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::t_Get_LambdaEta_Contribution(p_which);
	}

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	virtual void Alloc_LBlocks_Derived(CUberBlockMatrix &r_L) const
	{
		CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::Alloc_LBlocks(r_L);
	}

	virtual void Calculate_Omega_Derived(CUberBlockMatrix &r_omega, size_t n_min_vertex_order) const
	{
		CSEBaseEdgeImpl<CDerivedEdge, CVertex0, CVertex1,
			_n_measurement_dimension>::Calculate_Omega(r_omega, n_min_vertex_order);
	}

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---
};

#endif // __BASE_SE_PRIMITIVE_TYPES_INCLUDED
