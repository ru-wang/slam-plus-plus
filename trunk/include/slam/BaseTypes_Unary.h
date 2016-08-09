/*
								+-----------------------------------+
								|                                   |
								|      ***  Base SE types  ***      |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2014  |
								|                                   |
								|         BaseTypes_Unary.h         |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __BASE_SE_PRIMITIVE_TYPES_UNARY_EDGE_IMPL_SPECIALIZATION_INCLUDED
#define __BASE_SE_PRIMITIVE_TYPES_UNARY_EDGE_IMPL_SPECIALIZATION_INCLUDED

/**
 *	@file include/slam/BaseTypes_Unary.h
 *	@brief base SE primitive types
 *	@author -tHE SWINe-
 *	@date 2014-11-05
 *	@note This file is not supposed to be included by itself, it is included from include/slam/BaseTypes.h.
 */

#include "slam/BaseTypes.h" // for doxygen

/**
 *	@brief SE base edge implementation (specialization for unary edges)
 *
 *	@tparam CDerivedEdge is the name of derived edge class
 *	@tparam CVertexTypeList is list of types of the vertices
 *	@tparam _n_residual_dimension is residual vector dimension
 *	@tparam _n_storage_dimension is state vector storage dimension (or -1 if the same as _n_residual_dimension)
 *
 *	@todo Write a generic code, permitting any number of edges and remove this.
 */
template <class CDerivedEdge, class CVertexType, int _n_residual_dimension, int _n_storage_dimension /*= -1*/>
class CBaseEdgeImpl<CDerivedEdge, CTypelist<CVertexType, CTypelistEnd>,
	_n_residual_dimension, _n_storage_dimension> : public CBaseEdge {
public:
	typedef CTypelist<CVertexType, CTypelistEnd> _TyVertices; /**< @brief list of vertex types */
	typedef typename CTypelistItemAt<_TyVertices, 0>::_TyResult _TyVertex; /**< @brief name of the vertex class */

	/**
	 *	@brief copy of template parameters
	 */
	enum {
		n_vertex_num = 1, /**< @brief number of vertices */
		n_measurement_dimension = _n_residual_dimension, /**< @brief measurement vector dimension (edge dimension) @deprecated This is slightly unclear, use n_residual_dimension and n_storage_dimension instead. */
		n_residual_dimension = _n_residual_dimension, /**< @brief residual vector dimension */
		n_storage_dimension = (_n_storage_dimension == -1)? _n_residual_dimension : _n_storage_dimension, /**< @brief edge state storage dimension */
#ifdef __SE2_TYPES_USE_ALIGNED_VECTORS
		n_matrix_alignment = Eigen::AutoAlign | Eigen::ColMajor
#else // __SE2_TYPES_USE_ALIGNED_VECTORS
		n_matrix_alignment = Eigen::DontAlign | Eigen::ColMajor
#endif // __SE2_TYPES_USE_ALIGNED_VECTORS
	};

	/**
	 *	@brief local vertex traits
	 *	@tparam n_index is zero-based index of a vertex in CVertexTypeList
	 */
	template <const int n_index>
	class CVertexTraits {
	public:
		typedef typename CTypelistItemAt<_TyVertices, n_index>::_TyResult _TyVertex; /**< @brief type of vertex */

		/**
		 *	@brief copy of vertex parameters
		 */
		enum {
			n_dimension = _TyVertex::n_dimension, /**< @brief vertex dimension */
		};

		typedef fbs_ut::CCTSize2D<n_dimension, n_dimension> _TyMatrixSize; /**< @brief size of the matrix as fbs_ut::CCTSize2D */
		typedef typename _TyVertex::_TyVector _TyVector; /**< @brief state vector type */
		typedef typename _TyVertex::_TyMatrix _TyMatrix; /**< @brief a compatible matrix type */
		typedef typename _TyVertex::_TyVectorAlign _TyVectorAlign; /**< @brief aligned state vector type */
		typedef typename _TyVertex::_TyMatrixAlign _TyMatrixAlign; /**< @brief a compatible aligned matrix type */
		typedef Eigen::Matrix<double, n_residual_dimension, n_dimension> _TyJacobianMatrix; /**< @brief mixed edge / vertex dimension matrix type */
	};

	/**
	 *	@brief vertex initialization functor
	 *	Just clears the vertex to null.
	 */
	template <class CVertex = _TyVertex>
	class CInitializeNullVertex {
	public:
		/**
		 *	@brief function operator
		 *	@return Returns the value of the vertex being initialized.
		 */
		inline operator CVertex() const
		{
			typename CVertex::_TyVector v_null;
			v_null.setZero();
			return CVertex(v_null);
		}
	};

private:
	typedef Eigen::Matrix<double, n_residual_dimension, 1> _TyVector; /**< @brief edge dimension vector type */
	typedef Eigen::Matrix<double, n_storage_dimension, 1> _TyStorageVector; /**< @brief storage dimension vector type */
	typedef Eigen::Matrix<double, n_residual_dimension, n_residual_dimension> _TyMatrix; /**< @brief edge dimension matrix type */

	typedef Eigen::Matrix<double, n_storage_dimension, 1, n_matrix_alignment> _TyStorageVectorAlign; /**< @brief edge dimension storage vector type with member alignment */
	typedef Eigen::Matrix<double, n_residual_dimension, 1, n_matrix_alignment> _TyVectorAlign; /**< @brief edge dimension vector type with member alignment */
	typedef Eigen::Matrix<double, n_residual_dimension, n_residual_dimension, n_matrix_alignment> _TyMatrixAlign; /**< @brief edge dimension matrix type with member alignment */

protected:
	size_t m_p_vertex_id[n_vertex_num]; /**< @brief ids of referenced vertices */
	_TyVertex *m_p_vertex0; /**< @brief pointers to the referenced vertices */
	_TyStorageVectorAlign m_v_measurement; /**< @brief the measurement */
	_TyMatrixAlign m_t_sigma_inv; /**< @brief information matrix */

private:
#ifdef __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific members ---

	_TyMatrixAlign m_t_square_root_sigma_inv_upper; /**< @brief the R matrix (upper diagonal) = transpose(chol(t_sigma_inv)) */
	_TyVectorAlign m_v_error; /**< @brief error vector (needs to be filled explicitly by calling p_error_function) */
	double *m_p_RH[n_vertex_num]; /**< @brief blocks of memory where the R * H matrices are stored */

#endif // __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific members ---

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific members ---
#ifdef __LAMBDA_USE_V2_REDUCTION_PLAN

	double *m_p_HtSiH_vert[n_vertex_num];
	double *m_p_RHS_vert[n_vertex_num];
	// this is used only by the new lambda reduction

#else // __LAMBDA_USE_V2_REDUCTION_PLAN

	typename CVertexTraits<0>::_TyMatrixAlign m_t_HtSiH_vertex0; /**< @brief block of memory where Ht * inv(Sigma) * H matrix is stored (the vertex contribution) */
	typename CVertexTraits<0>::_TyVectorAlign m_t_right_hand_vertex0; /**< @brief block of memory where right-hand side vector is stored (the vertex contribution) */
	// t_odo - this also needs to be passed through the reductor, did not think of this (it is a slightly different
	// philosophy, it is supposed to reduce to a dense vector, but it will be reduced to a sparse blah. maybe do
	// the reduction using an offset table only)
	// this is only used by the oldschool lambda reduction

#endif // __LAMBDA_USE_V2_REDUCTION_PLAN
#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific members ---

public:
	__SE2_TYPES_ALIGN_OPERATOR_NEW

	/**
	 *	@brief default constructor; has no effect
	 */
	inline CBaseEdgeImpl()
	{}

	/**
	 *	@brief constructor
	 *
	 *	@param[in] n_vertex_id is id of the vertex
	 *	@param[in] r_v_measurement is the measurement vecot
	 *	@param[in] r_t_sigma_inv is inverse sigma matrix
	 *
	 *	@note This fills the structure, except for the vertex pointers.
	 *	@note Any derived class should check for vertex dimensions
	 *		being really what it expects (not checked elsewhere).
	 */
	inline CBaseEdgeImpl(size_t n_vertex_id, const _TyStorageVector &r_v_measurement,
		const _TyMatrix &r_t_sigma_inv)
		:m_v_measurement(r_v_measurement), m_t_sigma_inv(r_t_sigma_inv)
#ifdef __SE_TYPES_SUPPORT_A_SOLVERS
		, m_t_square_root_sigma_inv_upper(r_t_sigma_inv.llt().matrixU()) // calculate R
#endif // __SE_TYPES_SUPPORT_A_SOLVERS
	{
		m_p_vertex_id[0] = n_vertex_id;
		// get references to the vertices, initialize the vertices, if neccessary

		// note that errors, expectations and jacobian matrices are not cleared
	}

	/**
	 *	@brief updates the edge with new measurement
	 *
	 *	@param[in] r_v_measurement is the measurement vector
	 *	@param[in] r_t_sigma_inv is inverse sigma matrix
	 */
	inline void Update(const _TyStorageVector &r_v_measurement, const _TyMatrix &r_t_sigma_inv)
	{
		m_v_measurement = r_v_measurement;
		m_t_sigma_inv = r_t_sigma_inv;
#ifdef __SE_TYPES_SUPPORT_A_SOLVERS
		m_t_square_root_sigma_inv_upper = r_t_sigma_inv.llt().matrixU(); // calculate R
#endif // __SE_TYPES_SUPPORT_A_SOLVERS
	}

	/**
	 *	@copydoc base_iface::CConstEdgeFacade::n_Dimension()
	 */
	inline size_t n_Dimension() const
	{
		return n_residual_dimension;
	}

	/**
	 *	@copydoc base_iface::CConstEdgeFacade::n_Vertex_Num
	 */
	inline size_t n_Vertex_Num() const
	{
		return n_vertex_num;
	}

	/**
	 *	@copydoc base_iface::CConstEdgeFacade::n_Vertex_Id
	 */
	inline size_t n_Vertex_Id(size_t n_vertex) const
	{
		_ASSERTE(n_vertex < n_vertex_num);
		return m_p_vertex_id[n_vertex];
	}

	/**
	 *	@brief gets the measurement
	 *	@return Returns const reference to the measurement.
	 */
	inline const _TyStorageVectorAlign &v_Measurement() const
	{
		return m_v_measurement;
	}

	/**
	 *	@brief gets inverse sigma of the measurement
	 *	@return Returns const reference to the inverse sigma matrix.
	 */
	inline const _TyMatrixAlign &t_Sigma_Inv() const
	{
		return m_t_sigma_inv;
	}

	/**
	 *	@brief calculates Jacobian, based on a state of an arbitrary vertex
	 *
	 *	@param[out] r_t_jacobian0 is jacobian, corresponding to the vertex
	 *	@param[out] r_v_expectation is expectation
	 *	@param[out] r_v_error is error
	 *	@param[in] r_vertex0 is reference to the vertex
	 *	@param[in] r_v_measurement is measurement
	 *
	 *	@note In case the derived edge has extended state, it should implement
	 *		a custom version of this function.
	 */
	static inline void Calculate_Jacobians(
		typename CVertexTraits<0>::_TyJacobianMatrix &r_t_jacobian0,
		_TyVector &r_v_expectation, _TyVector &r_v_error,
		const typename CVertexTraits<0>::_TyVertex &r_vertex0,
		const _TyStorageVector &r_v_measurement)
	{
		_ASSERTE(n_vertex_num == 1); // this is always true, this is a specialization for binary edges
		// t_odo - support n-ary edges (need a different interface for that)

		CDerivedEdge dummy; // relies on having default constructor (most edges should)
		dummy.m_p_vertex0 = &const_cast<typename CVertexTraits<0>::_TyVertex&>(r_vertex0);
		dummy.m_v_measurement = r_v_measurement;
		((const CDerivedEdge&)dummy).Calculate_Jacobian_Expectation_Error(r_t_jacobian0, r_v_expectation, r_v_error);
		// Calculate_Jacobians_Expectation_Error in CDerivedEdge should have been
		// static and it would be simple, but now we can'T force the users to rewrite
		// all their code. meh.
	}

#ifdef __SE_TYPES_SUPPORT_A_SOLVERS // --- A-SLAM specific functions ---

	/**
	 *	@copydoc base_iface::CEdgeFacade::Alloc_JacobianBlocks
	 */
	inline void Alloc_JacobianBlocks(CUberBlockMatrix &r_A)
	{
		// easily implemented using a loop
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
		m_p_RH[0] = r_A.p_GetBlock_Log(m_n_id + 1, m_p_vertex_id[0],
			n_residual_dimension, CVertexTraits<0>::n_dimension, true, false);
#else // __BASE_TYPES_USE_ID_ADDRESSING
		m_p_RH[0] = r_A.p_FindBlock(m_n_order, m_p_vertex0->n_Order(),
			n_residual_dimension, CVertexTraits<0>::n_dimension, true, false);
#endif // __BASE_TYPES_USE_ID_ADDRESSING
		// just place blocks at (edge order, vertex order)

		_ASSERTE(m_p_RH[0]);
		// if this triggers, most likely __BASE_TYPES_USE_ID_ADDRESSING is enabled (see FlatSystem.h)
		// and the edges come in some random order
	}

	/**
	 *	@copydoc base_iface::CEdgeFacade::Calculate_Jacobians
	 */
	inline void Calculate_Jacobians()
	{
		// easily implemented using a loop
		typename CVertexTraits<0>::_TyJacobianMatrix t_jacobian0;
		_TyVector v_expectation, v_error;
		((CDerivedEdge*)this)->Calculate_Jacobian_Expectation_Error(t_jacobian0, v_expectation, v_error);

		m_v_error = v_error;
		// calculates the expectation, error and the jacobians (implemented by the edge)

		Eigen::Map<typename CVertexTraits<0>::_TyJacobianMatrix, CUberBlockMatrix::map_Alignment> t_RH0(m_p_RH[0]);
		// map hessian matrices

		t_RH0 = m_t_square_root_sigma_inv_upper * t_jacobian0;
		// recalculate RH (transpose cholesky of sigma times the jacobian)
		// note this references the A block matrix
	}

	/**
	 *	@copydoc base_iface::CConstEdgeFacade::Get_R_Error
	 */
	inline void Get_R_Error(Eigen::VectorXd &r_v_dest) const
	{
		r_v_dest.segment<n_residual_dimension>(m_n_order) = m_t_square_root_sigma_inv_upper * m_v_error;
	}

	/**
	 *	@copydoc base_iface::CConstEdgeFacade::Get_Error
	 */
	inline void Get_Error(Eigen::VectorXd &r_v_dest) const
	{
		r_v_dest.segment<n_residual_dimension>(m_n_order) = m_v_error;
	}

#endif // __SE_TYPES_SUPPORT_A_SOLVERS // --- ~A-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- lambda-SLAM specific functions ---
#ifdef __LAMBDA_USE_V2_REDUCTION_PLAN

	/**
	 *	@brief allocates hessian block matrices
	 *
	 *	@tparam CReductionPlan is lambda reduction plan ?nstantiation
	 *
	 *	@param[out] r_lambda is the target matrix where the blocks are stored
	 *	@param[in,out] r_rp is lambda reduction plan
	 *
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 */
	template <class CReductionPlan>
	inline void Alloc_HessianBlocks_v2(CUberBlockMatrix &r_lambda, CReductionPlan &r_rp)
	{
#ifdef __BASE_TYPES_USE_ID_ADDRESSING
		size_t n_id_0 = m_p_vertex_id[0];
		m_p_HtSiH_vert[0] = r_lambda.p_GetBlock_Log(n_id_0, n_id_0,
			CVertexTraits<0>::n_dimension, CVertexTraits<0>::n_dimension, true, false);
#else // __BASE_TYPES_USE_ID_ADDRESSING
		const size_t n_order_0 = m_p_vertex0->n_Order();
		m_p_HtSiH_vert[0] = r_lambda.p_FindBlock(n_order_0, n_order_0,
			CVertexTraits<0>::n_dimension, CVertexTraits<0>::n_dimension, true, false);
#endif // __BASE_TYPES_USE_ID_ADDRESSING
		// find a block for vertices' hessian on the diagonal (do not use the potentially swapped id / order)
		// note that if this is added after the edge hessian, it will be likely
		// added at the end of the block list in the matrix

		{
			typename CReductionPlan::CLambdaReductor &rp = r_rp.r_Lambda_ReductionPlan();

			size_t n_id_0 = m_p_vertex_id[0];
			m_p_HtSiH_vert[0] = rp.template p_Diagonal_GetTempBlock<typename
				CVertexTraits<0>::_TyMatrixSize>(n_id_0, n_id_0, m_p_HtSiH_vert[0]);
		}
		// make space for the vertex in the reduction scheme

		{
			typename CReductionPlan::CRHSReductor &rp = r_rp.r_RHS_ReductionPlan();

			m_p_RHS_vert[0] = rp.template p_Get_ReductionBlock<CVertexTraits<0>::n_dimension>(m_p_vertex0->n_Order());
		}
		// alloc the RHS reduction temporary

		// in virtual call scenario, either system or solver would need to specialize a virtual interface
		// on reduction plan, which would in turn need to export non-template functions with runtime seatch
		// of appropriate bin for a given block size
	}

	/**
	 *	@brief sums up edge hessian block contributions to get values of the blocks in lambda
	 *
	 *	@tparam CReductionPlan is lambda reduction plan ?nstantiation
	 *
	 *	@param[in] r_lambda is the target matrix where the blocks are stored
	 *	@param[in,out] r_rp is lambda reduction plan
	 *
	 *	@note This function is required for CNonlinearSolver_Lambda.
	 *	@note This only reduces the hessians, Calculate_Hessians_v2() must be called
	 *		before to have stuff to calculate.
	 */
	template <class CReductionPlan>
	inline void Reduce_Hessians_v2(const CUberBlockMatrix &r_lambda, CReductionPlan &r_rp) // todo - doc
	{
		typename CReductionPlan::CLambdaReductor &rp = r_rp.r_Lambda_ReductionPlan();
		const size_t n_id_0 = m_p_vertex_id[0];//, n_id_1 = m_p_vertex_id[1];

		typename CReductionPlan::CLambdaReductor::TKey p_key_vert[n_vertex_num];
		if(CReductionPlan::CLambdaReductor::b_use_block_coord_keys) // compile-time const
			p_key_vert[0] = CReductionPlan::CLambdaReductor::t_MakeKey(n_id_0, n_id_0, 0); // pointer ignored
		else {
			const double *p_HtSiH_vert[n_vertex_num];
			p_HtSiH_vert[0] = r_lambda.p_GetBlock_Log(n_id_0, n_id_0,
				CVertexTraits<0>::n_dimension, CVertexTraits<0>::n_dimension);
			_ASSERTE(p_HtSiH_vert[0]);
			// need the address of the unique *original* block, not of one of the many reduction temporaries
			// t_odo - perhaps it would be better to use (row, col) as unique index to the reduction, this lookup would then be avoided
			p_key_vert[0] = CReductionPlan::CLambdaReductor::t_MakeKey(n_id_0, n_id_0, p_HtSiH_vert[0]);
		}
		rp.template ReduceSingle<typename CVertexTraits<0>::_TyMatrixSize>(p_key_vert[0]);
		// reduce edge hessian and both vertex hessians
	}

	/**
	 *	@brief calculates hessian contributions
	 *
	 *	@note This only calculates the hessians, either Reduce_Hessians_v2() or the
	 *		reduction plan needs to be used to fill lambda with values.
	 */
	inline void Calculate_Hessians_v2()
	{
		typename CVertexTraits<0>::_TyJacobianMatrix t_jacobian0;
		_TyVector v_expectation, v_error;
		((CDerivedEdge*)this)->Calculate_Jacobian_Expectation_Error(t_jacobian0, v_expectation, v_error);
		// calculates the expectation and the jacobian

		Eigen::Map<typename CVertexTraits<0>::_TyMatrixAlign, CUberBlockMatrix::map_Alignment>
			m_t_HtSiH_vertex0(m_p_HtSiH_vert[0]);
		m_t_HtSiH_vertex0.noalias() = t_jacobian0.transpose() * m_t_sigma_inv * t_jacobian0;
		// calculate vertex hessian contributions

		Eigen::Map<typename CVertexTraits<0>::_TyVectorAlign, CUberBlockMatrix::map_Alignment>
			m_t_right_hand_vertex0(m_p_RHS_vert[0]);
		m_t_right_hand_vertex0.noalias() = t_jacobian0.transpose() * (m_t_sigma_inv * v_error);
		// calculate right hand side vector contributions
	}

#else // __LAMBDA_USE_V2_REDUCTION_PLAN

	inline bool Notify_HessianBlock_Conflict(double *p_block, CUberBlockMatrix &r_lambda)
	{
		// can't happen with unary edges
		return false; // not me
	}

	inline void Alloc_HessianBlocks(CUberBlockMatrix &r_lambda)
	{
		m_p_vertex0->Add_ReferencingEdge(this);
		// the vertices need to know this edge in order to calculate sum of jacobians on the diagonal of lambda // t_odo - should the lambda solver call this instead? might avoid duplicate records in the future // lambda solver is calling this only once in it's lifetime, no duplicates will occur

		// alloc nothing, the vertex will allocate the diagonal block
	}

	inline void Calculate_Hessians()
	{
		_ASSERTE(m_p_vertex0->b_IsReferencingEdge(this));
		// makes sure this is amongst the referencing edges

		typename CVertexTraits<0>::_TyJacobianMatrix t_jacobian0;
		_TyVector v_expectation, v_error;
		((CDerivedEdge*)this)->Calculate_Jacobian_Expectation_Error(t_jacobian0, v_expectation, v_error);
		// calculates the expectation and the jacobian

		m_t_HtSiH_vertex0.noalias() = t_jacobian0.transpose() * m_t_sigma_inv * t_jacobian0;
		// calculate vertex hessian contributions

		m_t_right_hand_vertex0.noalias() = t_jacobian0.transpose() * (m_t_sigma_inv * v_error);
		// calculate right hand side vector contributions
	}

	inline std::pair<const double*, const double*> t_Get_LambdaEta_Contribution(const CBaseVertex *p_which)
	{
		// no redfuction here, it all takes place in the vertex

		_ASSERTE(p_which == m_p_vertex0);
		return std::make_pair(m_t_HtSiH_vertex0.data(), m_t_right_hand_vertex0.data());
		// return pointer to the matrix data with hessian contribution
	}

#endif // __LAMBDA_USE_V2_REDUCTION_PLAN

	/**
	 *	@copydoc base_iface::CConstEdgeFacade::f_Max_VertexHessianDiagValue
	 */
	double f_Max_VertexHessianDiagValue() const // todo surround with levenberg support ifdef
	{
		double f_max = 0;

#ifdef __LAMBDA_USE_V2_REDUCTION_PLAN
		Eigen::Map<typename CVertexTraits<0>::_TyMatrix,
			CUberBlockMatrix::map_Alignment> m_t_HtSiH_vertex0(m_p_HtSiH_vert[0]);
		// need to map them first (note that storing a map instead of a pointer would save some hurdles)
#endif // __LAMBDA_USE_V2_REDUCTION_PLAN

		for(size_t i = 0; i < CVertexTraits<0>::n_dimension; ++ i)
			f_max = std::max(f_max, m_t_HtSiH_vertex0(i, i));
		return f_max;
	}

#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS // --- ~lambda-SLAM specific functions ---

#ifdef __SE_TYPES_SUPPORT_L_SOLVERS // --- L-SLAM specific functions ---

	/**
	 *	@copydoc base_iface::CConstEdgeFacade::Alloc_LBlocks
	 */
	inline void Alloc_LBlocks(CUberBlockMatrix &UNUSED(r_L)) const
	{
		// no offdiag blocks
	}

#endif // __SE_TYPES_SUPPORT_L_SOLVERS // --- ~L-SLAM specific functions ---

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
		const size_t n_dimension0 = CVertexTraits<0>::n_dimension;
		const size_t n_order_0 = m_p_vertex0->n_Order();
		// make sure the order is sorted (if swapping, will have to transpose the result,
		// but we will deal with that laters)

		_ASSERTE(n_order_0 >= n_min_vertex_order);
		double *p_v0 = r_omega.p_FindBlock(n_order_0 - n_min_vertex_order,
			n_order_0 - n_min_vertex_order, n_dimension0, n_dimension0, true, true);
		// alloc and initialize / find existing blocks for all the hessians, above the diagonal only

		bool b_recalculate = false;
		// recalculate from scratch, or use values already calculated by the lambda solver?

		typename CVertexTraits<0>::_TyJacobianMatrix t_jacobian0;
		if(b_recalculate) {
			_TyVector v_expectation, v_error;
			((CDerivedEdge*)this)->Calculate_Jacobian_Expectation_Error(t_jacobian0,
				v_expectation, v_error);
		}
		// calculates the expectation and the jacobians

		Eigen::Matrix<double, CVertexTraits<0>::n_dimension, n_residual_dimension> t_H0_sigma_inv;
		if(b_recalculate)
			t_H0_sigma_inv = t_jacobian0.transpose() * m_t_sigma_inv;

		Eigen::Map<typename CVertexTraits<0>::_TyMatrix,
			CUberBlockMatrix::map_Alignment> t_HtSiH_v0(p_v0);
		if(b_recalculate) {
			t_HtSiH_v0.noalias() += t_H0_sigma_inv * t_jacobian0;
		} else {
#ifdef __LAMBDA_USE_V2_REDUCTION_PLAN
			Eigen::Map<typename CVertexTraits<0>::_TyMatrix,
				CUberBlockMatrix::map_Alignment> t_HtSiH_vertex0(m_p_HtSiH_vert[0]);
			t_HtSiH_v0 += t_HtSiH_vertex0;
#else // __LAMBDA_USE_V2_REDUCTION_PLAN
			t_HtSiH_v0 += m_t_HtSiH_vertex0;
#endif // __LAMBDA_USE_V2_REDUCTION_PLAN
		}
		// calculate vertex hessian contributions (note the increments)
	}
};

#endif // !__BASE_SE_PRIMITIVE_TYPES_UNARY_EDGE_IMPL_SPECIALIZATION_INCLUDED
