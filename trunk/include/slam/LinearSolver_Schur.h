/*
								+-----------------------------------+
								|                                   |
								|   ***  Schur linear solver  ***   |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2013  |
								|                                   |
								|       LinearSolver_Schur.h        |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __LINEAR_SOLVER_SCHUR_INCLUDED
#define __LINEAR_SOLVER_SCHUR_INCLUDED

/**
 *	@file include/slam/LinearSolver_Schur.h
 *	@brief blocky linear solver wrapper that enables Schur complement optimization
 *	@author -tHE SWINe-
 *	@date 2013-01-28
 */

#include "slam/LinearSolverTags.h"
#include "csparse/cs.hpp"
#include "slam/BlockMatrix.h" // includes Eigen as well

/**
 *	@brief implementation of matrix orderings for Schur complement
 */
class CSchurOrdering {
public:
	/**
	 *	@brief a simple IOTA functor
	 */
	class CIOTA {
	protected:
		size_t m_n_counter; /**< @brief counter value */

	public:
		/**
		 *	@brief default constructor; initializes the counter
		 */
		inline CIOTA()
			:m_n_counter(-1)
		{}

		/**
		 *	@brief counter invokation
		 *	@return Returns the next number in sequence, starting with zero.
		 */
		inline size_t operator ()()
		{
			return ++ m_n_counter;
		}
	};

	/**
	 *	@brief a simple function object for vertex degree comparison
	 */
	class CCompareVertexDegree {
	protected:
		const cs *m_p_graph; /**< @brief pointer to the graph */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] p_graph is const pointer to the graph
		 */
		inline CCompareVertexDegree(const cs *p_graph)
			:m_p_graph(p_graph)
		{}

		/**
		 *	@brief vertex degree comparison operator
		 *
		 *	@param[in] n_vertex_a is zero-based index of the first vertex, in the given graph
		 *	@param[in] n_vertex_b is zero-based index of the second vertex, in the given graph
		 *
		 *	@return Returns true if the first vertex has lower degree than the second vertex.
		 */
		inline bool operator ()(size_t n_vertex_a, size_t n_vertex_b)
		{
			_ASSERTE(n_vertex_a < size_t(m_p_graph->n) && n_vertex_b < size_t(m_p_graph->n));
			return m_p_graph->p[n_vertex_a + 1] - m_p_graph->p[n_vertex_a] <
				m_p_graph->p[n_vertex_b + 1] - m_p_graph->p[n_vertex_b];
		}
	};

	enum {
		max_Recurse_Size = 128,
		thread_Num = 1 << 8, // make a bit more to allow for better scheduling (at each other step, the graph is only cut a little, need to have enough steps to make all cores busy with the big graphs)
		thread_Num_Log2 = n_Log2_Static(thread_Num)
	};

	/**
	 *	@brief explicit stack frame, for t_MIS_ExStack() and t_MIS_Parallel()
	 */
	struct TMIS_StackFrame {
		std::vector<size_t> complement_set; /**< @brief storage for the vertex complement set */
		size_t n_intermediate; /**< @brief intermediate value (used to store the pivot) */
		std::vector<size_t> int_result; /**< @brief intermediate result */
		const cs *p_graph; /**< @brief pointer to the original graph (not to be deleted) */
		cs *p_subgraph; /**< @brief pointer to a subgraph (to be deleted when this stack frame returns) */
		int n_phase; /**< @brief index of the phase of the algorithm, in which the stack frame was saved */
	};

public:
	/**
	 *	@brief calculates ordering for Schur
	 *
	 *	@param[out] r_ordering is filled with the ordering (does not need to be allocated)
	 *	@param[in] r_lambda is the matrix to be ordered on
	 *
	 *	@return Returns number of connected vertices (the size of the Schur factor).
	 *
	 *	@note There is some workspace, which could be reused.
	 */
	static size_t n_Calculate_Ordering(std::vector<size_t> &r_ordering,
		const CUberBlockMatrix &r_lambda) // throw(std::bad_alloc)
	{
		const size_t n = r_lambda.n_BlockColumn_Num();
		if(r_ordering.capacity() < n) {
			r_ordering.clear();
			r_ordering.reserve(std::max(2 * r_ordering.capacity(), n));
		}
		r_ordering.resize(n);
		// allocate space

		cs *p_lambda, *p_lambda_t, *p_lambda_ata;
		if(!(p_lambda = r_lambda.p_BlockStructure_to_Sparse()))
			throw std::bad_alloc(); // rethrow
		// get block structure as a sparse CSC

		if(!(p_lambda_t = cs_transpose(p_lambda, false))) {
			cs_spfree(p_lambda);
			throw std::bad_alloc(); // rethrow
		}
		if(!(p_lambda_ata = cs_add(p_lambda, p_lambda_t, 1, 1))) {
			cs_spfree(p_lambda);
			cs_spfree(p_lambda_t);
			throw std::bad_alloc(); // rethrow
		}
		cs_spfree(p_lambda);
		cs_spfree(p_lambda_t);
		p_lambda = p_lambda_ata;
		// modify it to have symmetric structure // todo - implement this directly in CUberBlockMatrix

		try {
			std::vector<size_t> mis = CSchurOrdering::t_MIS_FirstFit(p_lambda);
			// calculate MIS

			CSchurOrdering::Complement_VertexSet(r_ordering, mis, n);
			// collect the rest of the vertices not in MIS

			r_ordering.insert(r_ordering.end(), mis.begin(), mis.end());
			// append with MIS to gain the final ordering

			cs_spfree(p_lambda);
			// cleanup

			return n - mis.size();
			// return number of connected vertices
		} catch(std::bad_alloc &r_exc) {
			cs_spfree(p_lambda);
			throw r_exc; // rethrow
		}
	}

	/**
	 *	@brief calculates guided ordering for Schur
	 *
	 *	@param[out] r_ordering is filled with the ordering (does not need to be allocated)
	 *	@param[in] n_pose_vertex_dimension is dimension of pose (e.g. 3 for 2D SLAM)
	 *	@param[in] n_landmark_vertex_dimension is dimension of landmark (e.g. 2 for 2D SLAM)
	 *	@param[in] r_lambda is the matrix to be ordered on
	 *
	 *	@return Returns number of poses (the size of the Schur factor).
	 *
	 *	@note In case this is used, simpler linear solver can be used
	 *		as the block sizes are known beforehand, and are constant
	 *		for each section of the factorized matrix.
	 */
	static size_t n_Calculate_GuidedOrdering(std::vector<size_t> &r_ordering,
		size_t n_pose_vertex_dimension, size_t n_landmark_vertex_dimension,
		const CUberBlockMatrix &r_lambda) // throw(std::bad_alloc)
	{
		const size_t n = r_lambda.n_BlockColumn_Num();
		if(r_ordering.capacity() < n) {
			r_ordering.clear();
			r_ordering.reserve(std::max(2 * r_ordering.capacity(), n));
		}
		r_ordering.resize(n);
		// allocate space

		size_t n_pose_num = 0, n_landmark_num = 0;
		{
			size_t i = 0;

			for(; i < n; ++ i) {
				if(r_lambda.n_BlockColumn_Column_Num(i) == n_pose_vertex_dimension) {
					r_ordering[n_pose_num] = i;
					++ n_pose_num;
				} else
					break;
			}
			// look for all the poses

			n_landmark_num = n_pose_num; // offset the destination index
			for(; i < n; ++ i) {
				if(r_lambda.n_BlockColumn_Column_Num(i) == n_landmark_vertex_dimension) {
					r_ordering[n_landmark_num] = i;
					++ n_landmark_num;
				} else
					break;
			}
			n_landmark_num -= n_pose_num; // offset back
			// look for all the landmarks (assume landmarks are smaller than poses)

			if(i < n) {
				std::vector<size_t> &r_poses = r_ordering; // keep only poses in the destination ordering
				std::vector<size_t> landmarks(n - n_pose_num); // allocate space for the remaining landmarks
				// get memory

				std::copy(r_ordering.begin() + n_pose_num, r_ordering.begin() +
					(n_pose_num + n_landmark_num), landmarks.begin());
				// copy the landmarks away to a second array

				for(; i < n; ++ i) {
					if(r_lambda.n_BlockColumn_Column_Num(i) == n_pose_vertex_dimension) {
						r_poses[n_pose_num] = i;
						++ n_pose_num;
					} else {
						_ASSERTE(r_lambda.n_BlockColumn_Column_Num(i) == n_landmark_vertex_dimension);
						landmarks[n_landmark_num] = i;
						++ n_landmark_num;
					}
				}
				// loop through the rest of the vertices

				std::copy(landmarks.begin(), landmarks.begin() + n_landmark_num,
					r_ordering.begin() + n_pose_num);
				// copy the landmarks back
			}
		}
		_ASSERTE(n_pose_num + n_landmark_num == r_lambda.n_BlockColumn_Num());
		// calculate the simple guided ordering (assumes that landmarks have smaller
		// dimension than poses, and that landmarks are not connected among each other)

		return n_pose_num;
	}

	/**
	 *	@brief calculates MIS using sorted first-fit, followed by vertex swap improvement
	 *
	 *	@param[in] p_graph is A^T+A graph (the diagonal may or may not be present)
	 *
	 *	@return Returns (sorted) maximum independent set of vertices.
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This function is only approximate, and therefore, it is fast. Usually,
	 *		it comes within 10% of the perfect solution (at least on graph SLAM problems).
	 *		To calculate perfect MIS, use t_MIS_Parallel() or t_MIS_ExStack().
	 */
	static std::vector<size_t> t_MIS_FirstFit(const cs *p_graph) // throw(std::bad_alloc)
	{
		const size_t n = p_graph->n;

		std::vector<size_t> perm(n);
		std::generate(perm.begin(), perm.end(), CIOTA());
		// generate an identity permutation

		std::sort(perm.begin(), perm.end(), CCompareVertexDegree(p_graph));
		// make ordering, starting with smallest degree (turns out to be better
		// than starting with the greatest degree)
		// note that ordering with vertex folding could yield even better results

		std::vector<size_t> mis;
		mis.reserve(n / 2);
		// for most of the pose graphs, MIS of under half of the graph is resonable

		std::vector<bool> vertex_coverage(n, false); // need small auxiliary space
		for(size_t j = 0; j < n; ++ j) {
			size_t i = perm[j];
			// pick a vertex from permutation

			if(vertex_coverage[i])
				continue;
			// skip neighbours of vertices that are already in MIS

			mis.push_back(i);
			// add uncovered vertex to MIS

			for(size_t p0 = p_graph->p[i], p1 = p_graph->p[i + 1]; p0 < p1; ++ p0)
				vertex_coverage[p_graph->i[p0]] = true;
			// set all neighbours (including self) to covered
		}
		// a simple greedy approach to MIS

		std::vector<bool> &mis_membership = vertex_coverage;
		// note that vertex_coverage is not needed below, could reuse as mis_membership

		std::fill(mis_membership.begin(), mis_membership.end(), false);
		for(size_t j = 0, m = mis.size(); j < m; ++ j)
			mis_membership[mis[j]] = true;
		// need a fast lookup of what is in MIS

		bool b_swap;
		do {
			b_swap = false;

			for(size_t j = 0, m = mis.size(); j < m; ++ j) {
				size_t i = mis[j];
				// get a vertex that belongs to MIS

				mis_membership[i] = false;
				// assume we remove it

				std::vector<size_t> mis_add;
				for(size_t p0 = p_graph->p[i], p1 = p_graph->p[i + 1]; p0 < p1; ++ p0) {
					size_t v = p_graph->i[p0];
					if(v == i)
						continue;
					// get a vertex, that is a neighbor of i ...

					_ASSERTE(!mis_membership[v]);
					// for sure is not a member of MIS (as it is a neighbor of i), and it should be covered

					bool b_covered = false;
					for(size_t r0 = p_graph->p[v], r1 = p_graph->p[v + 1]; r0 < r1; ++ r0) {
						if(mis_membership[p_graph->i[r0]]) {
							b_covered = true;
							break;
						}
					}
					// see if this vertex is covered

					if(!b_covered) {
						mis_membership[v] = true; // assume temporary membership
						mis_add.push_back(v);
					}
					// this vertex could be traded for i
				}
				// go through neighbors of i

				if(mis_add.size() > 1) {
					size_t n_v0 = mis_add.front();
					mis[j] = n_v0;
					mis_membership[n_v0] = true;
					// add the first vertex in place of the old vertex

					mis.reserve(m = mis.size() + mis_add.size() - 1); // also update m
					for(size_t k = 1, o = mis_add.size(); k < o; ++ k) {
						size_t n_vk = mis_add[k];
						mis.push_back(n_vk);
						mis_membership[n_vk] = true;
					}
					// add the next vertices at the end

					-- j; // can try again on the first one (after all, the adding is greedy)
					b_swap = true;
					// we just traded some vertices
				} else if(!mis_add.empty() && CCompareVertexDegree(p_graph)(i, mis_add.front())) {
					_ASSERTE(mis_add.size() == 1); // not empty and less than 2
					size_t n_v0 = mis_add.front();
					mis[j] = n_v0;
					mis_membership[n_v0] = true;
					// put the vertex in place of the old vertex

					-- j; // can try again on this one (will not swap back, degree won't allow it)
					b_swap = true;
					// we just traded i for a vertex with smaller degree
					// (MIS can theoretically contain more vertices of smaller degree)
				} else {
					for(size_t k = 0, o = mis_add.size(); k < o; ++ k)
						mis_membership[mis_add[k]] = false; // remove the temporary memberships
					mis_membership[i] = true; // in the end, we do not remove it
				}
				// trade the vertices

#if 0
				_ASSERTE(m == mis.size());
				for(size_t k = 0; k < m; ++ k) {
					size_t i = mis[k];
					_ASSERTE(mis_membership[i]);
					for(size_t p0 = p_graph->p[i], p1 = p_graph->p[i + 1]; p0 < p1; ++ p0)
						_ASSERTE(p_graph->i[p0] == i || !mis_membership[p_graph->i[p0]]);
				}
				// debug checking - makes it really slow
#endif // 0
			}
		} while(b_swap);
		// improve MIS by trading membership of single vertices for more than one vertices

		if(n < 64) {
			// below algorithm is O(n^2), does not improve the solution much (maybe disable completely)
			// but sometimes it helps to find the perfect solution in small graphs

			bool b_swap_outer;
			do {
				b_swap_outer = false;

				for(size_t l = 0, m = mis.size(); l < m; ++ l) {
					const size_t v0 = mis[l];
					// get a vertex that belongs to MIS

					mis_membership[v0] = false;
					// assume we remove it

					std::vector<size_t> mis_add;
					for(size_t p0 = p_graph->p[v0], p1 = p_graph->p[v0 + 1]; p0 < p1; ++ p0) {
						size_t v = p_graph->i[p0];
						if(v == v0)
							continue;
						// get a vertex, that is a neighbor of v0 ...

						_ASSERTE(!mis_membership[v]);
						// for sure is not a member of MIS (as it is a neighbor of i), and it should be covered

						bool b_covered = false;
						for(size_t r0 = p_graph->p[v], r1 = p_graph->p[v + 1]; r0 < r1; ++ r0) {
							if(mis_membership[p_graph->i[r0]]) {
								b_covered = true;
								break;
							}
						}
						// see if this vertex is covered

						if(!b_covered) {
							mis_membership[v] = true; // assume temporary membership
							mis_add.push_back(v);
						}
						// this vertex could be traded for v0
					}
					// go through the neighbors of v0

					if(mis_add.empty()) {
						mis_membership[v0] = true; // !!
						continue;
					}
					_ASSERTE(mis_add.size() == 1); // the others already eliminated above or in the inner loop

					const size_t v1 = mis_add.front();
					_ASSERTE(mis_membership[v1] == true); // already set
					mis[l] = v1;
					// assume the swap helps something

					_ASSERTE(m == mis.size());
					size_t n_before_add = m;

					do {
						b_swap = false;

						for(size_t j = 0; j < m; ++ j) {
							size_t i = mis[j];
							// get a vertex that belongs to MIS

							mis_membership[i] = false;
							// assume we remove it

							std::vector<size_t> mis_add;
							for(size_t p0 = p_graph->p[i], p1 = p_graph->p[i + 1]; p0 < p1; ++ p0) {
								size_t v = p_graph->i[p0];
								if(v == i)
									continue;
								// get a vertex, that is a neighbor of i ...

								_ASSERTE(!mis_membership[v]);
								// for sure is not a member of MIS (as it is a neighbor of i), and it should be covered

								bool b_covered = false;
								for(size_t r0 = p_graph->p[v], r1 = p_graph->p[v + 1]; r0 < r1; ++ r0) {
									if(mis_membership[p_graph->i[r0]]) {
										b_covered = true;
										break;
									}
								}
								// see if this vertex is covered

								if(!b_covered) {
									mis_membership[v] = true; // assume temporary membership
									mis_add.push_back(v);
								}
								// this vertex could be traded for i
							}
							// go through the neighbors of i

							if(mis_add.size() > 1) {
								size_t n_v0 = mis_add.front();
								mis[j] = n_v0;
								mis_membership[n_v0] = true;
								// add the first vertex in place of the old vertex

								mis.reserve(m = mis.size() + mis_add.size() - 1); // also update m
								for(size_t k = 1, o = mis_add.size(); k < o; ++ k) {
									size_t n_vk = mis_add[k];
									mis.push_back(n_vk);
									mis_membership[n_vk] = true;
								}
								// add the next vertices at the end

								-- j; // can try again on the first one (after all, the adding is greedy)
								b_swap = true;
								// we just traded some vertices
							} else {
								for(size_t k = 0, o = mis_add.size(); k < o; ++ k)
									mis_membership[mis_add[k]] = false; // remove the temporary memberships
								mis_membership[i] = true; // in the end, we do not remove it
							}
							// trade the vertices
						}
					} while(b_swap);

					_ASSERTE(m >= n_before_add); // certainly not smaller
					if(m == n_before_add) { // no change
						mis[l] = v0;
						mis_membership[v0] = true;
						mis_membership[v1] = false;
						// the swap does not help anything, put it back
					} else
						b_swap_outer = true; // we did it, we managed a chain swap
					// see if the swap was proficient

#if 0
					_ASSERTE(m == mis.size());
					for(size_t k = 0; k < m; ++ k) {
						size_t i = mis[k];
						_ASSERTE(mis_membership[i]);
						for(size_t p0 = p_graph->p[i], p1 = p_graph->p[i + 1]; p0 < p1; ++ p0)
							_ASSERTE(p_graph->i[p0] == i || !mis_membership[p_graph->i[p0]]);
					}
					// debug checking - makes it really slow
#endif // 0
				}
			} while(b_swap_outer);
		}
		// try to swap chains of two vertices

		std::sort(mis.begin(), mis.end()); // "postorder" .. ?
		// note that mis is implicitly sorted, if there is no permutation

		// t_odo - approaches with sorting of the vertices based on the degree,
		// and swapping of vertices

		return mis;
	}

	/**
	 *	@brief calculates maximum independent set
	 *
	 *	@param[in] p_graph is A^T+A graph (the diagonal may or may not be present)
	 *
	 *	@return Returns the size of the maximum independent set.
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This function has exponential time complexity,
	 *		and is unfeasible to call it on sizable graphs.
	 */
	static size_t n_MIS(const cs *p_graph) // throw(std::bad_alloc)
	{
		if(!b_IsConnected(p_graph)) {
			std::vector<size_t> C;
			Get_MinConnectedSet(p_graph, C);
			// get minimal connected set

			cs *p_GC = p_Subgraph(p_graph, C.begin(), C.end());
			size_t n_result = n_MIS(p_GC);
			// recurse with the rest of the graph

			cs_spfree(p_GC);
			if(C.size() <= 2)
				n_result += 1; // one of vertices in C used (no matter which one)
			else {
				std::vector<size_t> inv_C;
				Complement_VertexSet(inv_C, C, p_graph->n);
				cs *p_C = p_Subgraph(p_graph, inv_C.begin(), inv_C.end());
				n_result += n_MIS(p_C);
				cs_spfree(p_C);
				// biggest subset; need to recurse
			}

			return n_result;
		}
		// handle disconnected graphs

		if(p_graph->n <= 1)
			return p_graph->n;
		// handle trivial graphs

		size_t n_greatest = 0;
		size_t n_degree = p_graph->p[1] - p_graph->p[0];
		const size_t n = p_graph->n;
		for(size_t i = 1; i < n; ++ i) {
			size_t n_deg = p_graph->p[i + 1] - p_graph->p[i];
			if(n_deg > n_degree) {
				n_degree = n_deg;
				n_greatest = i;
			}
		}
		// find a vertex with maximum degree (note that this is approximate as there is only the upper triangular, some vertices may have bigger degrees; would need n workspace to calculate degrees in O(nnz))

		cs *p_GB = p_Subgraph(p_graph, &n_greatest, &n_greatest + 1);
		size_t n_not_used = n_MIS(p_GB);
		cs_spfree(p_GB);
		// case when the vertex (B) is not used

		_ASSERTE(sizeof(csi) == sizeof(size_t));
		cs *p_GNB = p_Subgraph(p_graph, (size_t*)p_graph->i + p_graph->p[n_greatest],
			(size_t*)p_graph->i + p_graph->p[n_greatest + 1]);
		size_t n_used = n_MIS(p_GNB);
		cs_spfree(p_GNB);
		// case when the vertex (B) is used

		return std::max(n_used, n_not_used);
	}

	/**
	 *	@brief calculates maximum independent set
	 *
	 *	@param[in] p_graph is A^T+A graph (the diagonal may or may not be present)
	 *
	 *	@return Returns list of (zero-based) indices of vertices in the maximum independent set.
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This function has exponential time complexity,
	 *		and is unfeasible to call it on sizable graphs.
	 */
	static std::vector<size_t> t_MIS(const cs *p_graph) // throw(std::bad_alloc)
	{
		if(p_graph->n <= 1) {
			std::vector<size_t> result;
			if(p_graph->n)
				result.push_back(0);
			return result;
		}
		// handle trivial graphs

		if(!b_IsConnected(p_graph)) {
			std::vector<size_t> C;
			Get_MinConnectedSet(p_graph, C);
			// get minimal connected set

			cs *p_GC = p_Subgraph(p_graph, C.begin(), C.end());
			std::vector<size_t> result = t_MIS(p_GC); // t_odo - transform

			// recursion unroll save point 0

			Transform_Indices_Complement(result, C);
			// recurse with the rest of the graph

			cs_spfree(p_GC);
			if(C.size() <= 2) {
				size_t n_vertex = C.front(); // one of vertices in C used (no matter which one)
				result.insert(std::lower_bound(result.begin(), result.end(), n_vertex), n_vertex);
			} else {
				std::vector<size_t> inv_C;
				Complement_VertexSet(inv_C, C, p_graph->n);
				cs *p_C = p_Subgraph(p_graph, inv_C.begin(), inv_C.end());
				std::vector<size_t> result_C = t_MIS(p_C); // t_odo - transform

				// recursion unroll save point 1

				Transform_Indices_Complement(result_C, inv_C);
				result.insert(result.end(), result_C.begin(), result_C.end());
				cs_spfree(p_C);
				// biggest subset; need to recurse

				std::sort(result.begin(), result.end());
				// !!
			}

			return result;
		}
		// handle disconnected graphs

		size_t n_greatest = 0;
		size_t n_degree = p_graph->p[1] - p_graph->p[0];
		const size_t n = p_graph->n;
		for(size_t i = 1; i < n; ++ i) {
			size_t n_deg = p_graph->p[i + 1] - p_graph->p[i];
			if(n_deg > n_degree) {
				n_degree = n_deg;
				n_greatest = i;
			}
		}
		// find a vertex with maximum degree

		cs *p_GB = p_Subgraph(p_graph, &n_greatest, &n_greatest + 1);
		std::vector<size_t> not_used = t_MIS(p_GB); // t_odo - transform
		cs_spfree(p_GB);
		// case when the vertex (B) is not used

		// recursion unroll save point 2

		_ASSERTE(sizeof(csi) == sizeof(size_t));
		cs *p_GNB = p_Subgraph(p_graph, (size_t*)p_graph->i + p_graph->p[n_greatest],
			(size_t*)p_graph->i + p_graph->p[n_greatest + 1]);
		std::vector<size_t> used = t_MIS(p_GNB); // t_odo - transform
		cs_spfree(p_GNB);
		// case when the vertex (B) is used

		// recursion unroll save point 3

		if(not_used.size() > used.size() + 1) {
			Transform_Indices_Complement(not_used, &n_greatest, &n_greatest + 1);
			return not_used;
		} else {
			Transform_Indices_Complement(used, p_graph->i + p_graph->p[n_greatest],
				p_graph->i + p_graph->p[n_greatest + 1]);
			used.insert(std::lower_bound(used.begin(), used.end(), n_greatest), n_greatest);
			return used;
		}
		// transform indices and return
	}

	/**
	 *	@brief calculates maximum independent set, using explicit stack
	 *
	 *	@param[in] p_graph is A^T+A graph (the diagonal may or may not be present)
	 *
	 *	@return Returns list of (zero-based) indices of vertices in the maximum independent set.
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This function has exponential time complexity,
	 *		and is unfeasible to call it on sizable graphs.
	 *	@note This function uses explicit stack, so it can work on graphs of
	 *		arbitrary size (unline t_MIS()), but it is slightly slower.
	 */
	static std::vector<size_t> t_MIS_ExStack(const cs *p_graph) // throw(std::bad_alloc)
	{
		std::vector<size_t> recursion_result;

		TMIS_StackFrame t_base;
		t_base.p_graph = p_graph;
		t_base.p_subgraph = 0;
		t_base.n_phase = 0;
		std::vector<TMIS_StackFrame> stack;
		stack.push_back(t_base);
		while(!stack.empty()) {
			TMIS_StackFrame &t_fr = stack.back();

			p_graph = t_fr.p_graph;
			// assume current stack frame graph

			if(t_fr.n_phase == 0) {
				recursion_result.clear();
				// always begin with empty result

				if(p_graph->n <= 1) {
					//recursion_result.clear(); // no need, done just above
					if(p_graph->n)
						recursion_result.push_back(0);
					// leave result in recursion_result

					stack.erase(stack.end() - 1); // note that the reference to t_fr is invalidated
					continue;
					// pop stack frame, return to the one above
				}
				// handle trivial graphs

				if(!b_IsConnected(p_graph)) {
					std::vector<size_t> &C = t_fr.complement_set; // store in this frame
					Get_MinConnectedSet(p_graph, C);
					// get minimal connected set

					cs *p_GC = p_Subgraph(p_graph, C.begin(), C.end());
					t_fr.p_subgraph = p_GC; // save to delete later
					++ t_fr.n_phase; // jump to the next phase

					if(p_GC->n < max_Recurse_Size)
						recursion_result = t_MIS(p_GC); // the graph is small, can use implicit stack (faster)
					else {
						TMIS_StackFrame t_recurse;
						t_recurse.n_phase = 0;
						t_recurse.p_graph = p_GC;
						t_recurse.p_subgraph = 0;
						stack.push_back(t_recurse); // note that the reference to t_fr is invalidated
						continue;
						// push stack frame, recurse
					}
				} else
					t_fr.n_phase = 4; // jump to phase 4 (skip phase 1 - 3)
			}
			// handle the first phase of the algorithm

			if(t_fr.n_phase == 1) {
				_ASSERTE(p_graph->n > 1 && !b_IsConnected(p_graph));
				// we were inside the disconnected graph branch

				{
					std::vector<size_t> &C = t_fr.complement_set;
					cs *p_GC = t_fr.p_subgraph; // stored in this frame
					t_fr.p_subgraph = 0; // ...
					std::vector<size_t> &result = recursion_result; // got result from the recursion

					Transform_Indices_Complement(result, C);
					cs_spfree(p_GC);
					if(C.size() <= 2) {
						result.push_back(C.front()); // one of vertices in C used (no matter which one)
						t_fr.int_result.swap(result); // save result
						t_fr.n_phase = 3; // go to phase 3
					} else {
						{
							std::vector<size_t> inv_C;
							Complement_VertexSet(inv_C, C, p_graph->n);
							t_fr.complement_set.swap(inv_C);
						}
						std::vector<size_t> &inv_C = t_fr.complement_set;
						cs *p_C = p_Subgraph(p_graph, inv_C.begin(), inv_C.end());

						t_fr.p_subgraph = p_C; // save to delete later
						t_fr.int_result.swap(result); // save result
						++ t_fr.n_phase; // jump to the next phase

						if(p_C->n < max_Recurse_Size)
							recursion_result = t_MIS(p_C); // the graph is small, can use implicit stack (faster)
						else {
							TMIS_StackFrame t_recurse;
							t_recurse.n_phase = 0;
							t_recurse.p_graph = p_C;
							t_recurse.p_subgraph = 0;
							stack.push_back(t_recurse); // note that the reference to t_fr is invalidated
							continue;
							// push stack frame, recurse
						}
					}
				}
			}

			if(t_fr.n_phase == 2) {
				_ASSERTE(p_graph->n > 1 && !b_IsConnected(p_graph));
				// we were inside the disconnected graph branch, C was bigger than 2

				std::vector<size_t> &inv_C = t_fr.complement_set; // size smaller than or equal to n - 2
				cs *p_C = t_fr.p_subgraph;
				t_fr.p_subgraph = 0; // ...
				std::vector<size_t> &result = t_fr.int_result, &result_C = recursion_result;
				// resture locals

				{
					Transform_Indices_Complement(result_C, inv_C);
					result.insert(result.end(), result_C.begin(), result_C.end());
					cs_spfree(p_C);
					// biggest subset; need to recurse

					++ t_fr.n_phase; // move to the next phase
				}
			}

			if(t_fr.n_phase == 3) {
				_ASSERTE(p_graph->n > 1 && !b_IsConnected(p_graph));
				// we were inside the disconnected graph branch, C was bigger than 2

				std::vector<size_t> &result = t_fr.int_result;
				// resture locals

				std::sort(result.begin(), result.end()); // note that in the upper branch, sorting could be replaced by insertion at lower_bound
				// !!

				recursion_result.swap(result); // leave result in recursion_result
				stack.erase(stack.end() - 1); // note that the reference to t_fr is invalidated
				continue;
				// pop stack frame, return to the one above
			}

			if(t_fr.n_phase == 4) {
				_ASSERTE(p_graph->n > 1 && b_IsConnected(p_graph));
				// we are below the disconnected graph branch

				size_t n_greatest = 0;
				size_t n_degree = p_graph->p[1] - p_graph->p[0];
				const size_t n = p_graph->n;
				for(size_t i = 1; i < n; ++ i) {
					size_t n_deg = p_graph->p[i + 1] - p_graph->p[i];
					if(n_deg > n_degree) {
						n_degree = n_deg;
						n_greatest = i;
					}
				}
				// find a vertex with maximum degree

				cs *p_GB = p_Subgraph(p_graph, &n_greatest, &n_greatest + 1);

				t_fr.n_intermediate = n_greatest;
				t_fr.p_subgraph = p_GB;
				++ t_fr.n_phase;

				if(p_GB->n < max_Recurse_Size)
					recursion_result = t_MIS(p_GB); // the graph is small, can use implicit stack (faster)
				else {
					TMIS_StackFrame t_recurse;
					t_recurse.n_phase = 0;
					t_recurse.p_graph = p_GB;
					t_recurse.p_subgraph = 0;
					stack.push_back(t_recurse); // note that the reference to t_fr is invalidated
					continue;
					// push stack frame, recurse
				}
			}

			if(t_fr.n_phase == 5) {
				_ASSERTE(p_graph->n > 1 && b_IsConnected(p_graph));
				// we are below the disconnected graph branch

				cs *p_GB = t_fr.p_subgraph;
				t_fr.p_subgraph = 0; // ...
				size_t n_greatest = t_fr.n_intermediate;
				// restore locals

				cs_spfree(p_GB);
				// case when the vertex (B) is not used

				_ASSERTE(sizeof(csi) == sizeof(size_t));
				cs *p_GNB = p_Subgraph(p_graph, (size_t*)p_graph->i + p_graph->p[n_greatest],
					(size_t*)p_graph->i + p_graph->p[n_greatest + 1]);

				t_fr.int_result.swap(recursion_result); // save this
				t_fr.p_subgraph = p_GNB;
				++ t_fr.n_phase;

				if(p_GNB->n < max_Recurse_Size)
					recursion_result = t_MIS(p_GNB); // the graph is small, can use implicit stack (faster)
				else {
					TMIS_StackFrame t_recurse;
					t_recurse.n_phase = 0;
					t_recurse.p_graph = p_GNB;
					t_recurse.p_subgraph = 0;
					stack.push_back(t_recurse); // note that the reference to t_fr is invalidated
					continue;
					// push stack frame, recurse
				}
			}

			if(t_fr.n_phase == 6) {
				_ASSERTE(p_graph->n > 1 && b_IsConnected(p_graph));
				// we are below the disconnected graph branch

				cs *p_GNB = t_fr.p_subgraph;
				t_fr.p_subgraph = 0; // ...
				size_t n_greatest = t_fr.n_intermediate;
				std::vector<size_t> &not_used = t_fr.int_result;
				std::vector<size_t> &used = recursion_result;
				// restore locals

				cs_spfree(p_GNB);
				// case when the vertex (B) is used

				if(not_used.size() > used.size() + 1) {
					Transform_Indices_Complement(not_used, &n_greatest, &n_greatest + 1);

					recursion_result.swap(not_used); // leave result in recursion_result
					stack.erase(stack.end() - 1); // note that the reference to t_fr is invalidated
					continue;
					// pop stack frame, return to the one above
				} else {
					Transform_Indices_Complement(used, p_graph->i + p_graph->p[n_greatest],
						p_graph->i + p_graph->p[n_greatest + 1]);
					used.insert(std::lower_bound(used.begin(), used.end(), n_greatest), n_greatest);

					//recursion_result.swap(used); // already there
					stack.erase(stack.end() - 1); // note that the reference to t_fr is invalidated
					continue;
					// pop stack frame, return to the one above
				}
				// transform indices and return
			}
		}

		return recursion_result;
	}

	/**
	 *	@brief calculates maximum independent set, in parallel
	 *
	 *	@param[in] p_graph is A^T+A graph (the diagonal may or may not be present)
	 *
	 *	@return Returns list of (zero-based) indices of vertices in the maximum independent set.
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This function has exponential time complexity,
	 *		and is unfeasible to call it on sizable graphs.
	 *	@note This function uses explicit stack, so it can work on graphs of
	 *		arbitrary size (unline t_MIS()), but it is slightly slower.
	 *	@note This function uses OpenMP to run the calculation in parallel;
	 *		the scheduling is, however, only static - and thus imperfect.
	 */
	static std::vector<size_t> t_MIS_Parallel(const cs *p_graph) // throw(std::bad_alloc)
	{
		if(p_graph->n <= thread_Num_Log2)
			return t_MIS_ExStack(p_graph);
		// don't parallelize vary small graphs

		const int n_max_level = std::min(int(thread_Num_Log2), int(n_Log2(p_graph->n) - 1));
		const int n_thread_num = 1 << n_max_level;
		//printf("max_level: %d\n", n_max_level); // verbose
		// decide maxmium generator recursion

#ifdef _OPENMP
		int n_CPU_num;
		#pragma omp parallel
		{
			#pragma omp master
			{
				n_CPU_num = omp_get_num_threads();
			}
		}
		omp_set_num_threads(std::max(1, n_CPU_num - 2));
#endif // _OPENMP
		// don't use all the threads

		TMIS_StackFrame search_tree[thread_Num_Log2 + 1][thread_Num];
		search_tree[0][0].p_graph = p_graph;
		for(int i = 0; i < n_max_level; ++ i) {
			int n_stage_size = 1 << i;

			#pragma omp parallel for schedule(static, 1)
			for(int j = 0; j < n_stage_size; ++ j) {
				const cs *p_subgraph = search_tree[i][j].p_graph;
				// get the graph

				size_t n_greatest = 0;
				size_t n_degree = -1;
				const size_t n = p_subgraph->n;
				for(size_t k = 0; k < n; ++ k) {
					size_t n_deg = p_subgraph->p[k + 1] - p_subgraph->p[k];
					if(n_deg > n_degree) {
						n_degree = n_deg;
						n_greatest = k;
					}
				}
				// find a vertex with maximum degree

				search_tree[i][j].n_intermediate = n_greatest;

				_ASSERTE(sizeof(csi) == sizeof(size_t));
				cs *p_GB = p_Subgraph(p_subgraph, &n_greatest, &n_greatest + 1);
				cs *p_GNB = p_Subgraph(p_subgraph, (size_t*)p_subgraph->i + p_subgraph->p[n_greatest],
					(size_t*)p_subgraph->i + p_subgraph->p[n_greatest + 1]);
				// generate subgraphs

				search_tree[i + 1][j * 2 + 0].p_graph = p_GB;
				search_tree[i + 1][j * 2 + 1].p_graph = p_GNB;
				// generate two branches of the recursion tree
			}

			//printf(".");
		}
		//printf("\n");
		// generate a set of subgraphs to be solved in parallel

		#pragma omp parallel for schedule(dynamic, 1)
		for(int i = 0; i < n_thread_num; ++ i) {
			search_tree[n_max_level][i].int_result =
				t_MIS_ExStack(search_tree[n_max_level][i].p_graph);
			//printf("%d\n", int(search_tree[n_max_level][i].int_result.size()));
			//printf("*");
		}
		//printf("\n");
		// calculate the intermediate result for each subgraph (can be done in parallel)

		for(int i = n_max_level; i > 0;) {
			-- i;
			int n_stage_size = 1 << i;

			#pragma omp parallel for schedule(static, 1)
			for(int j = 0; j < n_stage_size; ++ j) {
				p_graph = search_tree[i][j].p_graph;
				// get the graph

				std::vector<size_t> &not_used = search_tree[i + 1][j * 2 + 0].int_result;
				std::vector<size_t> &used = search_tree[i + 1][j * 2 + 1].int_result;
				// get the results of the two subgraphs

				size_t n_greatest = search_tree[i][j].n_intermediate;

				if(not_used.size() > used.size() + 1) {
					Transform_Indices_Complement(not_used, &n_greatest, &n_greatest + 1);
					search_tree[i][j].int_result.swap(not_used);
				} else {
					Transform_Indices_Complement(used, p_graph->i + p_graph->p[n_greatest],
						p_graph->i + p_graph->p[n_greatest + 1]);
					used.insert(std::lower_bound(used.begin(), used.end(), n_greatest), n_greatest);
					search_tree[i][j].int_result.swap(used);
				}
				// choose the bigger MIS

				cs_spfree((cs*)search_tree[i + 1][j * 2 + 0].p_graph);
				cs_spfree((cs*)search_tree[i + 1][j * 2 + 1].p_graph);
				// delete the graphs
			}

			//printf(".");
		}
		//printf("\n");
		// reduce the search tree results

		return search_tree[0][0].int_result;
		// the result is in the root
	}

	/**
	 *	@brief determines whether a graph is connected
	 *	@param[in] p_graph is pointer to a graph symmetric matrix (the diagonal is not accessed)
	 *	@return Returns true if the given graph is connected, otherwise returns false.
	 *	@note This function throws std::bad_alloc.
	 */
	static bool b_IsConnected(const cs *p_graph) // throw(std::bad_alloc)
	{
		const size_t n = p_graph->n;
		if(!n)
			return true;

		size_t n_connected_num = 0;
		// calculates the number of connected vertices with vertex 0

		std::vector<bool> close(n, false);
		std::vector<size_t> open;
		open.push_back(0);
		while(!open.empty()) {
			size_t v = open.back();
			open.erase(open.end() - 1);
			if(close[v])
				continue; // already closed
			close[v] = true;
			++ n_connected_num;
			// get a vertex

			for(size_t p0 = p_graph->p[v], p1 = p_graph->p[v + 1]; p0 < p1; ++ p0) {
				size_t v1 = p_graph->i[p0];
				if(!close[v1] && std::find(open.begin(), open.end(), v1) == open.end())
					open.push_back(v1);
			}
			// add all the neighbours to open
		}
		// recurse the graph in DFS fashion; maybe unnecessarily complicated

		return n_connected_num == n;
	}

	/*void Get_MaxConnectedSet(const cs *p_graph, std::vector<size_t> &r_max_conn) // not required
	{
		const size_t n = p_graph->n;

		r_max_conn.clear();

		size_t n_connected_so_far = 0;
		size_t n_max_connected_size = 0;

		std::vector<size_t> open;
		std::vector<bool> close(n, false);
		std::vector<size_t> connected_set;
		for(;;) {
			size_t n_not_closed = 0;
			for(; n_not_closed < n; ++ n_not_closed) {
				if(!close[n_not_closed])
					break;
			}
			if(n_not_closed == n)
				break; // all vertices closed
			// find a vertex that is not closed

			connected_set.clear();
			open.push_back(n_not_closed);
			while(!open.empty()) {
				size_t v = open.back();
				open.erase(open.end() - 1);
				if(close[v])
					continue; // already closed
				close[v] = true;
				connected_set.push_back(v);
				// get a vertex

				for(size_t p0 = p_graph->p[v], p1 = p_graph->p[v + 1]; p0 < p1; ++ p0) {
					size_t v1 = p_graph->i[p0];
					if(!close[v1] && std::find(open.begin(), open.end(), v1) == open.end())
						open.push_back(v1);
				}
				// add all the neighbours to open
			}
			// get connected set with n_not_closed

			n_connected_so_far += connected_set.size();
			if(connected_set.size() > n_max_connected_size) {
				n_max_connected_size = connected_set.size();
				r_max_conn.swap(connected_set);
			}
			// look for the greatest connected set

			if(n - n_connected_so_far < n_max_connected_size)
				break; // early termination
			// even if all the remaining vertices were connected, it wouldn't be enough
		}

		_ASSERTE(!n || !r_max_conn.empty());
		// there should be at least one vertex in the maximum
		// connected set, in case the graph is not empty
	}*/

	/**
	 *	@brief gets minimum connected set
	 *
	 *	@param[in] p_graph is pointer to a graph symmetric matrix (the diagonal is not accessed)
	 *	@param[out] r_min_conn is filled with minimum connected set
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	static void Get_MinConnectedSet(const cs *p_graph, std::vector<size_t> &r_min_conn)
	{
		const size_t n = p_graph->n;

		r_min_conn.clear();
		if(!n)
			return;

		size_t n_min_connected_size = n + 1;

		std::vector<size_t> open;
		std::vector<bool> close(n, false);
		std::vector<size_t> connected_set;
		for(;;) {
			size_t n_not_closed = 0;
			for(; n_not_closed < n; ++ n_not_closed) {
				if(!close[n_not_closed])
					break;
			}
			if(n_not_closed == n)
				break; // all vertices closed
			// find a vertex that is not closed

			connected_set.clear();
			open.push_back(n_not_closed);
			while(!open.empty()) {
				size_t v = open.back();
				open.erase(open.end() - 1);
				if(close[v])
					continue; // already closed
				close[v] = true;
				connected_set.push_back(v);
				// get a vertex

				for(size_t p0 = p_graph->p[v], p1 = p_graph->p[v + 1]; p0 < p1; ++ p0) {
					size_t v1 = p_graph->i[p0];
					if(!close[v1] && std::find(open.begin(), open.end(), v1) == open.end())
						open.push_back(v1);
				}
				// add all the neighbours to open
			}
			// get connected set with n_not_closed

			if(connected_set.size() < n_min_connected_size) {
				n_min_connected_size = connected_set.size();
				r_min_conn.swap(connected_set);

				if(n_min_connected_size == 1)
					break;
				// won't get any smaller now
			}
			// look for the smallest connected set
		}

		_ASSERTE(!n || !r_min_conn.empty());
		// there should be at least one vertex in the maximum
		// connected set, in case the graph is not empty

		std::sort(r_min_conn.begin(), r_min_conn.end());
		// p_Subgraph() needs a sorted set
	}

	/**
	 *	@brief determines whether a set of vertex indices is sorted
	 *
	 *	@tparam _TyConstIter is vertex index const iterator type
	 *
	 *	@param[in] p_vertex_begin is iterator, pointing to the first vertex index
	 *	@param[in] p_vertex_end is iterator, pointing to one past the last vertex index
	 *
	 *	@return Returns true if the set is sorted in ascending
	 *		order, and does not contain duplicate elements.
	 */
	template <class _TyConstIter>
	static bool b_IsSortedSet(_TyConstIter p_vertex_begin, _TyConstIter p_vertex_end)
	{
		_ASSERTE(p_vertex_end >= p_vertex_begin);
		if(p_vertex_begin == p_vertex_end)
			return true;
		// empty set is sorted

		size_t n_prev = *p_vertex_begin;
		for(++ p_vertex_begin; p_vertex_begin != p_vertex_end; ++ p_vertex_begin) {
			size_t n_cur = *p_vertex_begin;
			if(n_prev >= n_cur)
				return false; // not sorted, or contains repeating elements
			n_prev = n_cur;
		}
		return true;
	}

	/**
	 *	@brief determines whether a set of vertex indices is sorted
	 *	@param[in] r_set is vector, containing a set of vertex indices
	 *	@return Returns true if the set is sorted in ascending
	 *		order, and does not contain duplicate elements.
	 */
	static inline bool b_IsSortedSet(const std::vector<size_t> &r_set)
	{
		return b_IsSortedSet(r_set.begin(), r_set.end());
	}

	/**
	 *	@brief builds a subgraph of a graph
	 *
	 *	@tparam _TyConstIter is vertex index const iterator type
	 *
	 *	@param[in] p_graph is pointer to a graph symmetric matrix (the diagonal is not accessed)
	 *	@param[in] p_vertex_begin is iterator, pointing to the first vertex index
	 *	@param[in] p_vertex_end is iterator, pointing to one past the last vertex index
	 *
	 *	@return Returns pointer to the specified subgraph
	 *		(the caller is responsible for freeing it).
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	template <class _TyConstIter>
	static cs *p_Subgraph(const cs *p_graph,
		_TyConstIter p_vertex_begin, _TyConstIter p_vertex_end) // throw(std::bad_alloc)
	{
		_ASSERTE(b_IsSortedSet(p_vertex_begin, p_vertex_end));

		const size_t n = p_graph->n;
		_ASSERTE(p_vertex_end >= p_vertex_begin);
		size_t n_remove_verts = p_vertex_end - p_vertex_begin;
		_ASSERTE(n >= n_remove_verts);

		cs *p_new_graph;
		if(!(p_new_graph = cs_spalloc(n - n_remove_verts,
		   n - n_remove_verts, p_graph->p[n], 0, 1)))
			throw std::bad_alloc(); // rehtrow
		// alloc workspace for a graph with up to last graph nnz entries (conservative estimate)

		size_t n_dest = 0;
		for(size_t v = 0; v < n; ++ v) {
			_TyConstIter p_lb_it = std::lower_bound(p_vertex_begin, p_vertex_end, v);
			if(p_lb_it != p_vertex_end && *p_lb_it == v)
				continue;
			// find out if the vertex is removed

			size_t n_col = v - (p_lb_it - p_vertex_begin);
			// calculate destination column (column minus number of vertices skipped)

			for(size_t p0 = p_graph->p[v], p1 = p_graph->p[v + 1]; p0 < p1; ++ p0) {
				size_t v1 = p_graph->i[p0];
				_TyConstIter p_lb2_it = std::lower_bound(p_vertex_begin, p_vertex_end, v1);
				if(p_lb2_it != p_vertex_end && *p_lb2_it == v1)
					continue;
				// find out if the (referenced) vertex is removed

				size_t n_row = v1 - (p_lb2_it - p_vertex_begin);
				// calculate destination row (row minus number of vertices skipped)

				p_new_graph->p[n_dest] = n_col;
				p_new_graph->i[n_dest] = n_row;
				++ n_dest;
			}
		}
		// assemble the graph in triplet form

		p_new_graph->nz = n_dest;
		cs *p_new_graph_c;
		if(!(p_new_graph_c = cs_compress(p_new_graph)))
			throw std::bad_alloc(); // rethrow
		cs_spfree(p_new_graph);
		// compress and free

		return p_new_graph_c;
	}

	/**
	 *	@brief calculates a complement of a set of vertex indices
	 *
	 *	@param[out] r_complement is vector, containing a complement of a set of vertex indices
	 *	@param[in] r_set is vector, containing a set of vertex indices
	 *	@param[in] n_graph_size is number of vertices in the whole graph
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	static void Complement_VertexSet(std::vector<size_t> &r_complement,
		const std::vector<size_t> &r_set, size_t n_graph_size) // throw(std::bad_alloc)
	{
		_ASSERTE(b_IsSortedSet(r_set));

		const size_t n = n_graph_size; // rename
		r_complement.clear();
		if(!r_set.empty()) {
			if(r_set.size() == n)
				return;
			// complement of everything is nothing

			for(size_t i = 0, m = r_set.size(); i < m; ++ i) {
				if(!i) {
					for(size_t j = 0, e = r_set[i]; j < e; ++ j)
						r_complement.push_back(j);
					// include all the vertices before the first one in the set
				} else {
					for(size_t j = r_set[i - 1] + 1, e = r_set[i]; j < e; ++ j)
						r_complement.push_back(j);
					// include all the vertices between two set elements
				}
			}
			for(size_t j = r_set.back() + 1; j < n; ++ j)
				r_complement.push_back(j);
			// include all the vertices after the last one in the set
		} else {
			for(size_t i = 0; i < n; ++ i)
				r_complement.push_back(i);
			// complement of nothing is everything
		}

		_ASSERTE(b_IsSortedSet(r_complement));
	}

	/**
	 *	@brief transforms vertex indices of a subgraph to vertex indices of the original graph,
	 *		using complement to a set of vertices of the subgraph
	 *
	 *	@tparam _TyConstIter is vertex index const iterator type
	 *
	 *	@param[in,out] r_indices is a list of indices to be transformed (works inplace)
	 *	@param[in] r_generating_set_complement is complement vertex set
	 */
	static inline void Transform_Indices_Complement(std::vector<size_t> &r_indices,
		std::vector<size_t> &r_generating_set_complement)
	{
		Transform_Indices_Complement(r_indices, r_generating_set_complement.begin(),
			r_generating_set_complement.end()); 
	}

	/**
	 *	@brief transforms vertex indices of a subgraph to vertex indices of the original graph,
	 *		using complement to a set of vertices of the subgraph
	 *
	 *	@tparam _TyConstIter is vertex index const iterator type
	 *
	 *	@param[in,out] r_indices is a list of indices to be transformed (works inplace)
	 *	@param[in] p_gsc_begin is iterator, pointing to the first complement vertex index
	 *	@param[in] p_gsc_end is iterator, pointing to one past the last complement vertex index
	 */
	template <class _TyConstIter>
	static void Transform_Indices_Complement(std::vector<size_t> &r_indices,
		_TyConstIter p_gsc_begin, _TyConstIter p_gsc_end)
	{
		_ASSERTE(p_gsc_begin <= p_gsc_end);
		_ASSERTE(b_IsSortedSet(p_gsc_begin, p_gsc_end));
		size_t n_off = 0;
		for(size_t i = 0, m = r_indices.size(); i < m; ++ i) {
			size_t v = r_indices[i];
			while(p_gsc_begin != p_gsc_end && size_t(*p_gsc_begin) <= v + n_off) {
				++ p_gsc_begin;
				++ n_off;
			}
			r_indices[i] = v + n_off;
		}
		_ASSERTE(b_IsSortedSet(r_indices));
		// transform indices back, based on the generating set complement
	}
};

/**
 *	@brief linear solver model based on Schur complement (a wrapper for another linear solver)
 *
 *	@tparam CBaseSolver is a type of basic linear solver,
 *		to be used to solve the Schur complement subproblem
 *	@tparam CAMatrixBlockSizes is is typelist, containing Eigen
 *		matrices with known compile-time sizes
 */
template <class CBaseSolver, class CAMatrixBlockSizes>
class CLinearSolver_Schur {
public:
	typedef CBlockwiseLinearSolverTag _Tag; /**< @brief solver type tag */
	typedef CBaseSolver _TyBaseSolver; /**< @brief name of the base linear solver */

	typedef typename CUniqueTypelist<CAMatrixBlockSizes>::_TyResult _TyAMatrixBlockSizes; /**< @brief list of block matrix sizes */
	typedef typename __fbs_ut::CBlockMatrixTypesAfterPreMultiplyWithSelfTranspose<
		_TyAMatrixBlockSizes>::_TyResult _TyLambdaMatrixBlockSizes; /**< @brief possible block matrices, found in lambda and L */
	//typedef typename _TyBaseSolver::_TySystem::_TyHessianMatrixBlockList _TyLambdaMatrixBlockSizes; // could use that instead, but there might be other prior knowledge of block sizes

	typedef typename CTransformTypelist<_TyAMatrixBlockSizes,
		__fbs_ut::CEigenToDimension>::_TyResult CDimsList; /**< @brief list of block sizes as CCTSize2D */
	typedef typename CUniqueTypelist<typename CTransformTypelist<CDimsList,
		__fbs_ut::CTransformDimensionColumnsToSize>::_TyResult>::_TyResult CWidthList; /**< @brief list of A block widths as CCTSize */
	typedef typename CSortTypelist<CWidthList,
		__fbs_ut::CCompareScalar_Less>::_TyResult _TyVertexSizeList; /**< @brief sorted list of all possible vertex sizes */

	/**
	 *	@brief configuration, stored as enum
	 */
	enum {
		vertex_Size_Num = CTypelistLength<_TyVertexSizeList>::n_result, /**< @brief number of different vertex dimensions */

		vertex_Size_SmallIdx = 0, /**< @brief index of "small" vertex size */
		vertex_Size_BigIdx = (vertex_Size_Num > 1)? 1 : 0, /**< @brief index of "big" vertex size */
		// choose index that is in the list (to avoid bounds violation)

		vertex_Size_Small = CTypelistItemAt<typename CConcatTypelist<_TyVertexSizeList,
			__fbs_ut::CCTSize<-1> >::_TyResult, vertex_Size_SmallIdx>::_TyResult::n_size, /**< @brief value of "small" vertex size */
		vertex_Size_Big = CTypelistItemAt<typename CConcatTypelist<_TyVertexSizeList,
			__fbs_ut::CCTSize<-1> >::_TyResult, vertex_Size_BigIdx>::_TyResult::n_size /**< @brief value of "big" vertex size */
		// padd the list with a single -1 item to allow compilation even with empty lists
		// (maybe not very useful, but at least will not give confusing complicated error messages)
	};

	typedef typename _TyBaseSolver::_Tag _TyBaseSolverTag; /**< @brief linear solver tag */

protected:
	_TyBaseSolver m_linear_solver; /**< @brief linear solver */

	typedef CLinearSolverWrapper<_TyBaseSolver, _TyBaseSolverTag> _TyLinearSolverWrapper; /**< @brief wrapper for the base linear solvers (shields solver capability to solve blockwise) */

	std::vector<double> m_double_workspace; /**< @brief temporary workspace, used to reorder the vectors (not used outside Schur_Solve()) */
	std::vector<size_t> m_order_workspace; /**< @brief temporary workspace, used to invert the ordering (not used outside n_Schur_Ordering()) */
	std::vector<size_t> m_order; /**< @brief ordering that separates matrix into diagonal part and dense parts */
	size_t m_n_matrix_cut; /**< @brief separation between the diagonal part and the dense parts (in blocks) */
	bool m_b_base_solver_reorder; /**< @brief base solver reordering flag */

public:
	/**
	 *	@brief default constructor (does nothing)
	 */
	inline CLinearSolver_Schur(_TyBaseSolver &r_solver)
		:m_linear_solver(r_solver), m_b_base_solver_reorder(true)
	{}

	/**
	 *	@brief copy constructor (has no effect)
	 *	@param[in] r_other is the solver to copy from
	 */
	inline CLinearSolver_Schur(const CLinearSolver_Schur &r_other)
		:m_linear_solver(r_other.m_linear_solver), m_b_base_solver_reorder(true)
	{}

	/**
	 *	@brief copy operator (has no effect; memory for lambda not copied)
	 *	@param[in] r_other is the solver to copy from
	 *	@return Returns reference to this.
	 */
	CLinearSolver_Schur &operator =(const CLinearSolver_Schur &r_other)
	{
		m_linear_solver = r_other.m_linear_solver;
		return *this;
	}

	/**
	 *	@brief solves linear system given by positive-definite matrix using Schur complement
	 *
	 *	@param[in] r_lambda is positive-definite matrix
	 *	@param[in,out] r_v_eta is the right-side vector, and is overwritten with the solution
	 *
	 *	@return Returns true on success, false on failure.
	 *
	 *	@note This function throws std::bad_alloc.
	 */
	bool Solve_PosDef(const CUberBlockMatrix &r_lambda, Eigen::VectorXd &r_v_eta) // throw(std::bad_alloc)
	{
		if(m_order.capacity() < r_lambda.n_BlockColumn_Num()) {
			m_order.clear(); // prevent copying
			m_order.reserve(std::max(2 * m_order.capacity(), r_lambda.n_BlockColumn_Num()));
		}
		m_order.resize(r_lambda.n_BlockColumn_Num());
		// allocate the ordering array

		size_t n_matrix_cut = n_Schur_Ordering(r_lambda, &m_order[0], m_order.size()); // calculate ordering
		// calculate the ordering

		return Schur_Solve(r_lambda, r_v_eta, r_v_eta, &m_order[0], m_order.size(), n_matrix_cut); // solve
	}

	/**
	 *	@brief deletes symbolic decomposition, if calculated (forces a symbolic
	 *		decomposition update in the next call to Solve_PosDef_Blocky())
	 */
	inline void Clear_SymbolicDecomposition()
	{
		m_order.clear();
		m_b_base_solver_reorder = true;
		m_n_matrix_cut = size_t(-1);
	}

	/**
	 *	@brief calculates symbolic decomposition of a block
	 *		matrix for later (re)use in solving it
	 *	@param[in] r_lambda is positive-definite matrix
	 *	@return Returns true on success, false on failure.
	 */
	bool SymbolicDecomposition_Blocky(const CUberBlockMatrix &r_lambda) // throw(std::bad_alloc)
	{
		if(m_order.capacity() < r_lambda.n_BlockColumn_Num()) {
			m_order.clear(); // prevent copying
			m_order.reserve(std::max(2 * m_order.capacity(), r_lambda.n_BlockColumn_Num()));
		}
		m_order.resize(r_lambda.n_BlockColumn_Num());
		// allocate the ordering array

		m_n_matrix_cut = n_Schur_Ordering(r_lambda, &m_order[0], m_order.size()); // calculate ordering
		// calculate the ordering

		m_b_base_solver_reorder = true;
		// will need to reorder later on

		return true;
	}

	/**
	 *	@brief solves linear system given by positive-definite matrix
	 *
	 *	@param[in] r_lambda is positive-definite matrix
	 *	@param[in,out] r_v_eta is the right-side vector, and is overwritten with the solution
	 *
	 *	@return Returns true on success, false on failure.
	 *
	 *	@note This function throws std::bad_alloc.
	 *	@note This enables a reuse of previously calculated symbolic decomposition,
	 *		it can be either calculated in SymbolicDecomposition_Blocky(), or it is
	 *		calculated automatically after the first call to this function,
	 *		or after Clear_SymbolicDecomposition() was called (preferred).
	 */
	bool Solve_PosDef_Blocky(const CUberBlockMatrix &r_lambda, Eigen::VectorXd &r_v_eta) // throw(std::bad_alloc)
	{
		if(r_lambda.n_BlockColumn_Num() != m_order.size())
			SymbolicDecomposition_Blocky(r_lambda);
		// in case the ordering is nonconforming, calculate a new one

		bool b_keep_ordering = !m_b_base_solver_reorder;
		m_b_base_solver_reorder = false; // only once, or until the ordering changes
		// in case we need to reorder the base linear solver as well

		/*Eigen::VectorXd sol_copy = r_v_eta;
		_TyBaseSolver solver;
		solver.Solve_PosDef(r_lambda, sol_copy);*/
		// calculate reference "good" solution

		bool b_result = Schur_Solve(r_lambda, r_v_eta, r_v_eta, &m_order[0],
			m_order.size(), m_n_matrix_cut, b_keep_ordering); // solve

		/*sol_copy -= r_v_eta;
		printf("error in solution is %g (head %g, tail %g)\n",
			sol_copy.norm(), sol_copy.head(m_n_matrix_cut * 3).norm(),
			sol_copy.tail((r_lambda.n_BlockColumn_Num() - m_n_matrix_cut) * 2).norm());*/
		// calculate difference in the solution

		return b_result;
	}

	/**
	 *	@brief computes solution for linear system using schur complement
	 */
	bool Schur_Solve(const CUberBlockMatrix &r_lambda,
		const Eigen::VectorXd &r_v_eta, Eigen::VectorXd &r_v_schur_sol,
		const size_t *p_ordering, size_t UNUSED(n_ordering_size), size_t n_matrix_cut,
		bool b_reuse_block_structure = false)
	{
		_ASSERTE(r_lambda.b_SymmetricLayout()); // pos-def is supposed to be symmetric
		_ASSERTE(r_v_eta.rows() == r_lambda.n_Row_Num()); // make sure the vector has correct dimension

		// note that &r_v_eta == &r_v_schur_sol is permissible (can solve inplace)
		_ASSERTE(r_lambda.b_SymmetricLayout());
		_ASSERTE(r_lambda.n_BlockColumn_Num() == n_ordering_size);
		_ASSERTE(n_matrix_cut > 0 && n_matrix_cut + 1 < n_ordering_size);
		_ASSERTE(r_v_schur_sol.rows() == r_lambda.n_Column_Num());
		_ASSERTE(r_v_eta.rows() == r_lambda.n_Column_Num());

		CUberBlockMatrix lambda_perm;
		r_lambda.Permute_UpperTriangular_To(lambda_perm,
			p_ordering, n_ordering_size, true);
		// reorder the matrix

		const size_t n = lambda_perm.n_BlockColumn_Num();

		CUberBlockMatrix A, U, C, V;
		lambda_perm.SliceTo(A, 0, n_matrix_cut, 0, n_matrix_cut, true);
		lambda_perm.SliceTo(U, 0, n_matrix_cut, n_matrix_cut, n, true);
		lambda_perm.SliceTo(C, n_matrix_cut, n, n_matrix_cut, n, true);
		U.TransposeTo(V); // because lower-triangular of lambda is not calculated
		/*CUberBlockMatrix &A = lambda_perm;
		A.SliceTo(A, n_matrix_cut, n_matrix_cut, true);*/ // can't, would free data that U, C and V references
		// cut Lambda matrix into pieces
		// \lambda = | A U |
		//           | V C |

		CUberBlockMatrix C_inv;
		C_inv.InverseOf_Symmteric_FBS<_TyLambdaMatrixBlockSizes>(C); // C is block diagonal (should also be symmetric)
		// inverse of C

		/*CUberBlockMatrix unity, u_ref;
		unity.ProductOf(C, C_inv);
		unity.CopyLayoutTo(u_ref);
		u_ref.SetIdentity();
		unity.AddTo(u_ref, -1);
		double f_error = u_ref.f_Norm();
		fprintf(stderr, "error of matrix inverse is: %g\n", f_error);*/
		// check inverse

		CUberBlockMatrix minus_U_Cinv;
		U.MultiplyToWith_FBS<_TyLambdaMatrixBlockSizes,
			_TyLambdaMatrixBlockSizes>(minus_U_Cinv, C_inv);	// U*(C^-1) // U is not symmetric, it is rather dense, tmp is again dense
		minus_U_Cinv.Scale(-1.0);	// -U*(C^-1)
		CUberBlockMatrix schur_compl; // not needed afterwards
		minus_U_Cinv.MultiplyToWith_FBS<_TyLambdaMatrixBlockSizes,
			_TyLambdaMatrixBlockSizes>(schur_compl, V, true); // -U*(C^-1)V // UV is symmetric, the whole product should be symmetric, calculate only the upper triangular part
		A.AddTo_FBS<_TyLambdaMatrixBlockSizes>(schur_compl); // -U*(C^-1)V + A // A is symmetric, if schur_compl is symmetric, the difference also is
		// compute left-hand side A - U(C^-1)V

		// note that the sum and difference of two symmetric matrices is again symmetric,
		// but this is not always true for the product

		/*lambda_perm.Save_MatrixMarket("lambda_perm.mtx", "lambda_perm.bla");
		A.Save_MatrixMarket("lambda_perm00.mtx", "lambda_perm00.bla");
		U.Save_MatrixMarket("lambda_perm01.mtx", "lambda_perm01.bla");
		V.Save_MatrixMarket("lambda_perm10.mtx", "lambda_perm10.bla");
		C.Save_MatrixMarket("lambda_perm11.mtx", "lambda_perm11.bla");
		C_inv.Save_MatrixMarket("lambda_perm11_inv.mtx", "lambda_perm11_inv.bla");
		schur_compl.Save_MatrixMarket("schur.mtx", "schur.bla");*/
		/*lambda_perm.Rasterize("schur0_lambda_perm.tga", 3);
		A.Rasterize("schur1_lambda_perm00.tga", 3);
		U.Rasterize("schur2_lambda_perm01.tga", 3);
		V.Rasterize("schur3_lambda_perm10.tga", 3);
		C.Rasterize("schur4_lambda_perm11.tga", 3);
		schur_compl.Rasterize("schur5_A-(UC-1V).tga", 3);
		C_inv.Rasterize("schur6_lambda_perm11_inv.tga", 3);*/
		// debug

		size_t n_rhs_vector_size = r_lambda.n_Column_Num();
		size_t n_pose_vector_size = A.n_Column_Num(); // 6 * n_matrix_cut;
		size_t n_landmark_vector_size = U.n_Column_Num(); // 3 * (n - n_matrix_cut);
		// not block columns! element ones

		if(m_double_workspace.capacity() < n_rhs_vector_size) {
			m_double_workspace.clear(); // avoid data copying
			m_double_workspace.reserve(std::max(2 * m_double_workspace.capacity(), n_rhs_vector_size));
		}
		m_double_workspace.resize(n_rhs_vector_size);
		double *p_double_workspace = &m_double_workspace[0];
		// alloc workspace

		lambda_perm.InversePermute_RightHandSide_Vector(p_double_workspace,
			&r_v_eta(0), n_rhs_vector_size, p_ordering, n_ordering_size);
		// need to permute the vector !!

		Eigen::VectorXd v_x = Eigen::Map<Eigen::VectorXd>(p_double_workspace, n_pose_vector_size); // don't really need a copy, but need Eigen::VectorXd for _TyLinearSolverWrapper::Solve()
		Eigen::VectorXd v_l = Eigen::Map<Eigen::VectorXd>(p_double_workspace +
			n_pose_vector_size, n_landmark_vector_size); // need a copy, need one vector of workspace
		// get eta and cut it into pieces
		// \eta = | x |
		//        | l |

		// we are now solving:
		// \lambda          \eta
		// | A U | | dx | = | x |
		// | V C | | dl |   | l |

		minus_U_Cinv.PreMultiply_Add_FBS<_TyLambdaMatrixBlockSizes>(&v_x(0),
			n_pose_vector_size, &v_l(0), n_landmark_vector_size); // x - U(C^-1)l
		// compute right-hand side x - U(C^-1)l

		if(!b_reuse_block_structure)
			_TyLinearSolverWrapper::FinalBlockStructure(m_linear_solver, schur_compl); // the ordering on schur_compl will not change, can calculate it only in the first pass and then reuse
		bool b_result = _TyLinearSolverWrapper::Solve(m_linear_solver, schur_compl, v_x);
		Eigen::VectorXd &v_dx = v_x; // rename, solves inplace
		// solve for dx = A - U(C^-1)V / x

		// note that schur_compl only contains pose-sized blocks when guided ordering is used! could optimize for that
		// also note that schur_compl is not completely dense if it is not many times smaller than C

		Eigen::Map<Eigen::VectorXd>(p_double_workspace, n_pose_vector_size) = v_dx; // an unnecessary copy, maybe could work arround
		// obtained a first part of the solution

		Eigen::Map<Eigen::VectorXd> v_dl(p_double_workspace +
			n_pose_vector_size, n_landmark_vector_size); // calculated inplace
		v_dl.setZero();
		V.PreMultiply_Add_FBS<_TyLambdaMatrixBlockSizes>(&v_dl(0),
			n_landmark_vector_size, &v_dx(0), n_pose_vector_size); // V * dx
		v_l -= v_dl; // (l - dl)
		v_dl.setZero();
		C_inv.PreMultiply_Add_FBS<_TyLambdaMatrixBlockSizes>(&v_dl(0),
			n_landmark_vector_size, &v_l(0), n_landmark_vector_size); // (C^-1)(l - V * dx)
		// solve for dl = (C^-1)(l - V * dx)
		// the second part of the solution is calculated inplace in the dest vector

		lambda_perm.Permute_RightHandSide_Vector(&r_v_schur_sol(0), p_double_workspace,
			n_rhs_vector_size, p_ordering, n_ordering_size);
		// permute back!

		// todo - do some profiling

		return b_result;
	}

	/**
	 *	@brief computes ordering on lambda that spearate poses from landmarks
	 */
	size_t n_Schur_Ordering(const CUberBlockMatrix &r_lambda,
		size_t *p_ordering, size_t UNUSED(n_ordering_size)) // throw(std::bad_alloc)
	{
		_ASSERTE(p_ordering);
		_ASSERTE(n_ordering_size == r_lambda.n_BlockColumn_Num());

		size_t n_pose_num;
		if(vertex_Size_Num == 2) { // count of 2 also implies that the sizes are different
			n_pose_num = CSchurOrdering::n_Calculate_GuidedOrdering(m_order_workspace,
				vertex_Size_Big, vertex_Size_Small, r_lambda);
			// this is the "guided ordering", need to have two types of vertices with a different dimensionality
			// note this may give inferior results, i.e. on landmark datasets like Victoria park
		}
		if(vertex_Size_Num != 2 || n_pose_num >= r_lambda.n_BlockColumn_Num() / 2) { // pose-only without -po also tries to use guided ordering (but there are no landmarks - wasted time)
			n_pose_num = CSchurOrdering::n_Calculate_Ordering(m_order_workspace, r_lambda);
			// this is full ordering, using the graph structure, always gives best results
		}

		/*printf("debug: Schur inverse size: " PRIsize " vertices (out of " PRIsize ")\n",
			r_lambda.n_BlockColumn_Num() - n_pose_num, r_lambda.n_BlockColumn_Num());*/
		// see how successful the ordering was

		for(size_t i = 0; i < n_ordering_size; ++ i)
			p_ordering[m_order_workspace[i]] = i;
		//	p_ordering[i] = m_order_workspace[i]; // bad ordering - to test inverse of not purely diagonal matrix
		// have to use inverse ordering for matrix reordering

		return n_pose_num;
		// calculate ordering

		// t_odo - calculate proper schur ordering based on the graph dissection
	}
};

#endif // __LINEAR_SOLVER_SCHUR_INCLUDED
