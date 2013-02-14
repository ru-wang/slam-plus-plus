/*
								+-----------------------------------+
								|                                   |
								|      ***  .graph Parser  ***      |
								|                                   |
								|   Copyright  � -tHE SWINe- 2012   |
								|                                   |
								|             Parser.h              |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __GRAPH_PARSER_INCLUDED
#define __GRAPH_PARSER_INCLUDED

/**
 *	@file include/slam/Parser.h
 *	@brief a simple .graph file parser
 *	@author -tHE SWINe-
 *	@date 2012-03-30
 *
 *	@date 2012-04-02
 *
 *	Added the CParseEntity_Vertex2D and CParseEntity_Edge2D classes.
 *	These contain edge and vertex data for 2D case. Note that these should
 *	not be used directly, use one of the derived classes instead.
 *
 *	Added the explicit "2D" to some type names to be clear that they
 *	are to be used for the 2D SLAM case.
 *
 *	Changed the wording from "pose" to "node", as suggested by Ela.
 *
 *	The matrices are stored as "row-major upper triangular and diagonal".
 *	That is, the sequece of elements (0, 1, 2, 3, 4, 5) make up the following matrix:
 *	@code
 *	|0 1 2|
 *	|1 3 4|
 *	|2 4 5|
 *	@endcode
 *	Which corresponds with the matrix ordering, used in the "french code":
 *	@code
 *	|0    |
 *	|1 3  |
 *	|2 4 5|
 *	@endcode
 *	That is hereby deemed correct (and the small slam example is incorrect).
 *
 *	@date 2012-04-20
 *
 *	Clarified format in which the data come from the parser:
 *
 *	Edges contain the vertices in origin - endpoint order. Sigma matrices actually
 *	contain component-wise square roots of inverse covariance.
 *
 *	And so it shall be (watch out).
 */

#include <string>
#include <vector>
#include <stdio.h>
#include <ctype.h>
#include <map>
#include <set>
#include "slam/TypeList.h"
#include "slam/Integer.h"
#include "eigen/Eigen/Dense"
#include "slam/2DSolverBase.h" // CBase2DSolver::Absolute_to_Relative() and such
#include "slam/3DSolverBase.h" // CBase3DSolver::Absolute_to_Relative() and such

#if 0 // legacy code

/**
 *	@def __SLAM_CONSTRAIN_EIGEN_TYPE_SIZE
 *	@brief constrains size of eigen matrices and vectors, effectively avoiding
 *		dynamic allocation (but sometimes wasting memory)
 *	@note Use with caution (and _DEBUG first).
 */
#define __SLAM_CONSTRAIN_EIGEN_TYPE_SIZE
// saves .5 sec on 10k (16%)
// saves .9 sec on vp (first 1500) (0%)
// is supposed to save 23% (according to profiler)

//typedef Eigen::Matrix<double, 2, 2> Matrix2d; /**< @brief 2 by 2 matrix */
//typedef Eigen::Matrix<double, 2, 3> Matrix2x3d; /**< @brief 2 by 3 matrix */
//typedef Eigen::Matrix<double, 3, 3> Matrix3d; /**< @brief 3 by 3 matrix */
//typedef Eigen::Matrix<double, 2, 1> Vector2d; /**< @brief 2D vector */
//typedef Eigen::Matrix<double, 3, 1> Vector3d; /**< @brief 3D vector */
//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd; /**< @brief dynamically sized vector */
//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd; /**< @brief dynamically sized matrix */
// t_odo - typedef matrix and vector type here, reuse throughout the sources

#ifndef __SLAM_CONSTRAIN_EIGEN_TYPE_SIZE
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_constrained; /**< @brief dynamically sized constrained matrix */
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd_constrained; /**< @brief dynamically sized constrained vector */
#else // !__SLAM_CONSTRAIN_EIGEN_TYPE_SIZE
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign, 3, 3> MatrixXd_constrained; /**< @brief dynamically sized constrained matrix */
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::AutoAlign, 3, 1> VectorXd_constrained; /**< @brief dynamically sized constrained vector */
#endif // !__SLAM_CONSTRAIN_EIGEN_TYPE_SIZE
// on 10k, it is possible to save up to .5 sec (out of 3) if matrix and vector size is constrained (static allocation). npot doesn't seem to matter, better make it small.

#endif // 0

/**
 *	@brief a simple .graph file parser
 */
class CParserBase {
public:
	/**
	 *	@brief base parsed entity class
	 */
	class CParseEntity {};

	/**
	 *	@brief range-bearing measurement base class
	 */
	struct TLandmark2D : public CParseEntity {
		int m_n_node_0; /**< @brief (zero-based) index of the first node (origin) */
		int m_n_node_1; /**< @brief (zero-based) index of the second node (endpoint) */
		Eigen::Vector2d m_v_delta; /**< @brief dealte measurement (also called "z") */
		Eigen::Matrix2d m_t_inv_sigma; /**< @brief inverse sigma matrix, elements are square roots (also called "Sz") */

		/**
		 *	@brief default constructor
		 *
		 *	@param[in] n_node_0 is (zero-based) index of the first (origin) node
		 *	@param[in] n_node_1 is (zero-based) index of the second (endpoint) node
		 *	@param[in] f_delta_x is delta x position
		 *	@param[in] f_delta_y is delta y position
		 *	@param[in] p_upper_matrix_2x2 is row-major upper triangular and diagonal 2x2 matrix
		 *		containing transformation between the nodes, elements are square roots
		 *
		 *	The matrix is stored row by row from top to bottom,
		 *	with left-to-right column order. Example:
		 *
		 *	@code
		 *	|0 1|
		 *	|  2|
		 *	@endcode
		 */
		inline TLandmark2D(int n_node_0, int n_node_1,
			double f_delta_x, double f_delta_y, const double *p_upper_matrix_2x2)
			:m_n_node_0(n_node_0), m_n_node_1(n_node_1), m_v_delta(f_delta_x, f_delta_y)
		{
			m_t_inv_sigma <<
				p_upper_matrix_2x2[0], p_upper_matrix_2x2[1],
				p_upper_matrix_2x2[1], p_upper_matrix_2x2[2];
			// fill the matrix
		}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	/**
	 *	@brief XYT measurement base class
	 */
	struct TEdge2D : public CParseEntity {
		int m_n_node_0; /**< @brief (zero-based) index of the first (origin) node */
		int m_n_node_1; /**< @brief (zero-based) index of the second (endpoint) node */
		Eigen::Vector3d m_v_delta; /**< @brief dealte measurement (also called "z") */
		Eigen::Matrix3d m_t_inv_sigma; /**< @brief inverse sigma matrix, elements are square roots (also called "Sz") */

		/**
		 *	@brief default constructor
		 *
		 *	@param[in] n_node_0 is (zero-based) index of the first (origin) node
		 *	@param[in] n_node_1 is (zero-based) index of the second (endpoint) node
		 *	@param[in] f_delta_x is delta x position
		 *	@param[in] f_delta_y is delta y position
		 *	@param[in] f_delta_theta is delta theta
		 *	@param[in] p_upper_matrix_3x3 is row-major upper triangular and diagonal 3x3 matrix
		 *		containing transformation between the nodes, elements are square roots
		 *
		 *	The matrix is stored row by row from top to bottom,
		 *	with left to right column order. Example:
		 *	@code
		 *	|0 1 2|
		 *	|  3 4|
		 *	|    5|
		 *	@endcode
		 */
		inline TEdge2D(int n_node_0, int n_node_1,
			double f_delta_x, double f_delta_y, double f_delta_theta,
			const double *p_upper_matrix_3x3)
			:m_n_node_0(n_node_0), m_n_node_1(n_node_1),
			m_v_delta(f_delta_x, f_delta_y, f_delta_theta)
		{
			m_t_inv_sigma <<
				p_upper_matrix_3x3[0], p_upper_matrix_3x3[1], p_upper_matrix_3x3[2],
				p_upper_matrix_3x3[1], p_upper_matrix_3x3[3], p_upper_matrix_3x3[4],
				p_upper_matrix_3x3[2], p_upper_matrix_3x3[4], p_upper_matrix_3x3[5];
			// fill the matrix
		}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	/**
	 *	@brief node measurement class ("VERTEX" or "VERTEX2" in the datafile)
	 *	@note Vertices are not processed by the french code.
	 */
	struct TVertex2D : public CParseEntity {
		int m_n_id; /**< @brief vertex id */
		Eigen::Vector3d m_v_position; /**< @brief vertex position (the state vector) */

		/**
		 *	@brief default constructor
		 *
		 *	@param[in] n_node_id is (zero-based) index of the node
		 *	@param[in] f_x is x position
		 *	@param[in] f_y is y position
		 *	@param[in] f_theta is theta
		 *
		 *	@note The vertices are only used as ground truth. Those are mostly not processed.
		 */
		inline TVertex2D(int n_node_id, double f_x, double f_y, double f_theta)
			:m_n_id(n_node_id), m_v_position(f_x, f_y, f_theta)
		{}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	/**
	 *	@brief XYZRPY measurement base class
	 */
	struct TEdge3D : public CParseEntity {
		int m_n_node_0; /**< @brief (zero-based) index of the first (origin) node */
		int m_n_node_1; /**< @brief (zero-based) index of the second (endpoint) node */
		Eigen::Matrix<double, 6, 1> m_v_delta; /**< @brief dealte measurement (also called "z") */
		Eigen::Matrix<double, 6, 6> m_t_inv_sigma; /**< @brief inverse sigma matrix, elements are square roots (also called "Sz") */

		/**
		 *	@brief default constructor
		 *
		 *	@param[in] n_node_0 is (zero-based) index of the first (origin) node
		 *	@param[in] n_node_1 is (zero-based) index of the second (endpoint) node
		 *	@param[in] f_delta_x is delta x position
		 *	@param[in] f_delta_y is delta y position
		 *	@param[in] f_delta_z is delta z position
		 *	@param[in] f_delta_roll is rotation around X axis
		 *	@param[in] f_delta_pitch is rotation around Y axis
		 *	@param[in] f_delta_yaw is is rotation around Z axis
		 *	@param[in] p_upper_matrix_6x6 is row-major upper triangular and diagonal 6x6 matrix
		 *		containing transformation between the nodes, elements are square roots
		 *
		 *	The matrix is stored row by row from top to bottom,
		 *	with left to right column order. Example:
		 *	@code
		 *	|0 1  2  3  4  5|
		 *	|  6  7  8  9 10|
		 *	|    11 12 13 14|
		 *	|       15 16 17|
		 *	|          18 19|
		 *	|             20|
		 *	@endcode
		 */
		inline TEdge3D(int n_node_0, int n_node_1,
			double f_delta_x, double f_delta_y, double f_delta_z,
			double f_delta_roll, double f_delta_pitch, double f_delta_yaw,
			const double *p_upper_matrix_6x6)
			:m_n_node_0(n_node_0), m_n_node_1(n_node_1)
		{
			m_v_delta << f_delta_x, f_delta_y, f_delta_z,
				f_delta_roll, f_delta_pitch, f_delta_yaw;
			// no constructor for 6-valued vector

			const double *p_u = p_upper_matrix_6x6;
			m_t_inv_sigma <<
				p_u[0],  p_u[1],  p_u[2],  p_u[3],  p_u[4],  p_u[5],
				p_u[1],  p_u[6],  p_u[7],  p_u[8],  p_u[9], p_u[10],
				p_u[2],  p_u[7], p_u[11], p_u[12], p_u[13], p_u[14],
				p_u[3],  p_u[8], p_u[12], p_u[15], p_u[16], p_u[17],
				p_u[4],  p_u[9], p_u[13], p_u[16], p_u[18], p_u[19],
				p_u[5], p_u[10], p_u[14], p_u[17], p_u[19], p_u[20];
			// fill the matrix
		}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	/**
	 *	@brief Point projection measurement base class
	 */
	struct TProjection3D : public CParseEntity {
		int m_n_node_0; /**< @brief (zero-based) index of the first (origin) node */
	 	int m_n_node_1; /**< @brief (zero-based) index of the second (endpoint) node */
	 	Eigen::Vector2d m_v_delta; /**< @brief dealte measurement (also called "z") */
	 	Eigen::Matrix2d m_t_inv_sigma; /**< @brief inverse sigma matrix, elements are square roots (also called "Sz") */

	 	/**
		 *	@brief default constructor
		 *
		 *	@param[in] n_node_0 is (zero-based) index of the 3D point
		 *	@param[in] n_node_1 is (zero-based) index of the view the point is seen in
		 *	@param[in] f_delta_x is delta x position
		 *	@param[in] f_delta_y is delta y position
		 *	@param[in] p_upper_matrix_2x2 is row-major upper triangular and diagonal 2x2 matrix
		 *		containing transformation between the nodes, elements are square roots
		 *
		 *	The matrix is stored row by row from top to bottom,
		 *	with left to right column order. Example:
		 *	@code
		 *	|0 1|
		 *	|  2|
		 *	@endcode
		 */
	 	inline TProjection3D(int n_node_0, int n_node_1,
		 	double f_delta_x, double f_delta_y, const double *p_upper_matrix_2x2) // you can change the parameters as you wish
 	    	:m_n_node_0(n_node_0), m_n_node_1(n_node_1), // fill the indices of edge vertices (these should be zero-based)
 	        m_v_delta(f_delta_x, f_delta_y) // fill the measurement vector
 	    {
 	        m_t_inv_sigma <<
 	            p_upper_matrix_2x2[0], p_upper_matrix_2x2[1],
 	            p_upper_matrix_2x2[1], p_upper_matrix_2x2[2]; // fill the inverse sigma matrix
 	    }

	 	EIGEN_MAKE_ALIGNED_OPERATOR_NEW // this is imposed by the use of eigen, copy it
	};

	/**
	 *	@brief node measurement class ("VERTEX3" in the datafile)
	 */
	struct TVertex3D : public CParseEntity {
		int m_n_id; /**< @brief vertex id */
		Eigen::Matrix<double, 6, 1> m_v_position; /**< @brief vertex position (the state vector) */

		/**
		 *	@brief default constructor
		 *
		 *	@param[in] n_node_id is (zero-based) index of the node
		 *	@param[in] f_x is x position
		 *	@param[in] f_y is y position
		 *	@param[in] f_z is z position
		 *	@param[in] ax is x part of axis vector
		 *	@param[in] ay is y part of axis vector
		 *	@param[in] az is z part of axis vector
		 *
		 *	@note The vertices are only used as ground truth. Those are mostly not processed.
		 */
		inline TVertex3D(int n_node_id, double f_x, double f_y, double f_z, double ax, double ay, double az)
			:m_n_id(n_node_id)
		{
			m_v_position << f_x, f_y, f_z, ax, ay, az;
			// no constructor for 6-valued vector
		}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	/**
	 *	@brief a simple callback class, to be used by the parser
	 *	@note This has the disadvantage of having to be modified after adding a new
	 *		parsed types, and having to rewrite all implementations of this. This is
	 *		now only used by CDatasetPeeker.
	 */
	class CParserAdaptor {
	public:
		/**
		 *	@brief appends the system with an odometry measurement
		 *	@param[in] r_t_edge is the measurement to be appended
		 */
		virtual void AppendSystem(const TEdge2D &r_t_edge) = 0;

		/**
		 *	@brief appends the system with a landmark measurement
		 *	@param[in] r_t_landmark is the measurement to be appended
		 */
		virtual void AppendSystem(const TLandmark2D &r_t_landmark) = 0;

		/**
		 *	@brief appends the system with vertex position
		 *	@param[in] r_t_vertex is the vertex to be appended
		 *	@note The vertices can be ignored in most of the solvers.
		 */
		virtual void AppendSystem(const TVertex2D &r_t_vertex) = 0;

		/**
		 *	@brief appends the system with an odometry measurement
		 *	@param[in] r_t_edge is the measurement to be appended
		 */
		virtual void AppendSystem(const TEdge3D &r_t_edge) = 0;

		/**
		 *	@brief appends the system with a 3D projection measurement
		 *	@param[in] r_t_edge is the measurement to be appended
		 */
		virtual void AppendSystem(const TProjection3D &r_t_edge) = 0;

		/**
		 *	@brief appends the system with vertex position
		 *	@param[in] r_t_vertex is the vertex to be appended
		 *	@note The vertices can be ignored in most of the solvers.
		 */
		virtual void AppendSystem(const TVertex3D &r_t_vertex) = 0;
	};

protected:
	/**
	 *	@brief a simple function object which adds token names to a lookup map
	 *	@note This is actually used by CParserTemplate, but it is not dependent on its parameters.
	 */
	class CAddToMap {
	protected:
		int m_n_index; /**< @brief index of the current parse primitive */
		std::map<std::string, int> &m_r_token_map; /**< @brief reference to the token map */

	public:
		/**
		 *	@brief default constructor
		 *	@param[in] r_token_map is a reference to the token map (output)
		 */
		inline CAddToMap(std::map<std::string, int> &r_token_map)
			:m_n_index(0), m_r_token_map(r_token_map)
		{}

		/**
		 *	@brief function operator
		 *	@tparam CParsePrimitive is a parse primitive type to be added to the map
		 *	@note This function throws std::bad_alloc.
		 */
		template <class CParsePrimitive>
		inline void operator ()() // throws(std::bad_alloc)
		{
			CParsePrimitive::EnumerateTokens(m_r_token_map, m_n_index);
			++ m_n_index; // the original version had this one at compile time; a bit more elegant
			// add all tokens for this primitive type
		}
	};

	static std::map<std::string, std::pair<std::set<std::string>,
		std::pair<bool, bool> > > m_per_file_warned_token_set; /**< @brief contains persistent warned tokens per file, in order to not warn twice (once per dataset peeker, once per main loop) */

public:
	/**
	 *	@brief removes whitespace from the beginning and from the end of the string
	 *	@param[in,out] r_s_string is the string to remove whitespace from
	 */
	static void TrimSpace(std::string &r_s_string);

	/**
	 *	@brief reads line form a file
	 *
	 *	@param[out] r_s_line is output string, containing one line read from a file
	 *	@param[in] p_fr is pointer to a file
	 *
	 *	@return Returns true on success, false on failure (not enough memory / input error).
	 *
	 *	@note In case file is at it's end, output lines are empty, but the function still succeeds.
	 *	@note Output lines may contain carriage-return character(s), for example if the file
	 *		is opened for binary reading. Line-feed character marks end of line and is never
	 *		included.
	 */
	static bool ReadLine(std::string &r_s_line, FILE *p_fr);

	/**
	 *	@brief splits a string by a separator
	 *
	 *	@param[out] r_s_dest is destination vector for substrings
	 *	@param[in] r_s_string is the input string
	 *	@param[in] p_s_separator is the separator
	 *	@param[in] n_thresh is minimal output string threshold
	 *		(only strings longer than threshold are stored in r_s_dest)
	 *
	 *	@return Returns true on success, false on failure (not enough memory).
	 */
	static bool Split(std::vector<std::string> &r_s_dest, const std::string &r_s_string,
		const char *p_s_separator, int n_thresh);
};

/**
 *	@brief template of a parser implementation with a list of parsed primitives
 *
 *	@tparam CParseLoopType is a type of parse loop (the sink); it can be CParserBase::CParserAdaptor,
 *		a specialization of CParseLoop or something completely different
 *	@tparam CParsePrimitiveTypelist is a list of parse primitive handlers
 */
template <class CParseLoopType, class CParsePrimitiveTypelist>
class CParserTemplate : public CParserBase {
protected:
	typedef CParseLoopType _TyParseLoop; /**< @brief parse loop type */
	typedef typename CUniqueTypelist<CParsePrimitiveTypelist>::_TyResult
		_TyParsePrimitiveTypelist; /**< @brief a list of parse primitive handlers */

	/**
	 *	@brief a simple function object that passes a line of text to
	 *		a parse primitive and dispatches the parsed item to parse loop
	 */
	class CPrimitiveFire {
	protected:
		size_t m_n_line_no; /**< @brief current line number (zero-based, for error reporting) */
		const std::string &m_r_s_line; /**< @brief the current line without the token */
		const std::string &m_r_s_token; /**< @brief the name of the token */
		_TyParseLoop &m_r_parse_loop; /**< @brief parse loop to pass the parsed item to */
		bool m_b_result; /**< @brief the result of the operation */

	public:
		/**
		 *	@brief default constructor; just copies the arguments
		 *
		 *	@param[in] n_line_no is current line number (zero-based, for error reporting)
		 *	@param[in] r_s_line is the current line without the token
		 *	@param[in] r_s_token is the name of the token
		 *	@param[in,out] r_parse_loop is the parse loop to pass the parsed item to
		 */
		inline CPrimitiveFire(size_t n_line_no, const std::string &r_s_line,
			const std::string &r_s_token, _TyParseLoop &r_parse_loop)
			:m_n_line_no(n_line_no), m_r_s_line(r_s_line),
			m_r_s_token(r_s_token), m_r_parse_loop(r_parse_loop), m_b_result(false)
		{}

		/**
		 *	@brief performs the parsing using the selected parse primitive
		 *	@tparam CParsePrimitive is the selected parse primitive for the current token
		 */
		template <class CParsePrimitive>
		inline void operator ()()
		{
			m_b_result = CParsePrimitive::Parse_and_Dispatch(m_n_line_no,
				m_r_s_line, m_r_s_token, m_r_parse_loop);
			// parse this using the primitive
		}

		/**
		 *	@brief gets the result of the parsing
		 *	@return Returns true on success, false on failure.
		 */
		operator bool() const
		{
			return m_b_result;
		}
	};

public:
	/**
	 *	@brief parses a .graph file, emitting parsed entities in the process
	 *
	 *	@param[in] p_s_filename is a null-terminated string containing input file name
	 *	@param[in] r_callback is a reference to the callback object that will receive parsed entities
	 *	@param[in] n_line_limit is the limit of the number of processed lines (use e.g. 1500 for debugging)
	 *	@param[in] b_ignore_case is case sensitivity flag; if set, token case is ignored
	 *
	 *	@return Returns true on success, false on failure.
	 */
	bool Parse(const char *p_s_filename, CParseLoopType &r_callback,
		size_t n_line_limit = 0, bool b_ignore_case = false)
	{
		FILE *p_fr;
#if defined(_MSC_VER) && !defined(__MWERKS__) && _MSC_VER >= 1400
		if(fopen_s(&p_fr, p_s_filename, "r"))
#else //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
		if(!(p_fr = fopen(p_s_filename, "r")))
#endif //_MSC_VER && !__MWERKS__ && _MSC_VER >= 1400
			return false;
		// open the input file

		try {
			std::map<std::string, int> token_name_map;
			CTypelistForEach<_TyParsePrimitiveTypelist, CAddToMap>::Run(CAddToMap(token_name_map));
			// fill the token map, based on the parsed primitives

			bool b_new_file = m_per_file_warned_token_set.find(std::string(p_s_filename)) ==
				m_per_file_warned_token_set.end();
			bool &r_b_warned_about_duplicates =
				m_per_file_warned_token_set[std::string(p_s_filename)].second.first;
			bool &r_b_warned_about_reverse_order =
				m_per_file_warned_token_set[std::string(p_s_filename)].second.second;
			if(b_new_file) {
				r_b_warned_about_duplicates = false;
				r_b_warned_about_reverse_order = false;
			}
			std::set<std::string> &r_bad_token_name_set =
				m_per_file_warned_token_set[std::string(p_s_filename)].first;
			// set of unknown token names (remember them for application runtime)

			std::string s_line, s_token;
			for(size_t n_line_no = 0; !feof(p_fr) && (!n_line_limit || n_line_no < n_line_limit); ++ n_line_no) {
				if(!ReadLine(s_line, p_fr)) {
					fclose(p_fr);
					return false;
				}
				TrimSpace(s_line);
				if(s_line.empty())
					continue;
				// read the file line by line, skip the empty lines

				s_token.erase();
				size_t n_pos;
				if((n_pos = s_line.find_first_of(" \t")) != std::string::npos) {
					s_token.insert(s_token.begin(), s_line.begin(), s_line.begin() + n_pos);
					s_line.erase(0, n_pos + 1);
					TrimSpace(s_line); // ...
				} else {
					s_token.insert(s_token.begin(), s_line.begin(), s_line.end());
					s_line.erase();
				}
				// read the token (separated by space; otherwise the entire line is a token)
				// and remove the token from the line

				if(b_ignore_case)
					std::for_each(s_token.begin(), s_token.end(), toupper);
				// map is case sensitive, if there is dataset with lowercase tokens, it needs to be fixed

				std::map<std::string, int>::const_iterator p_tok_it;
				if((p_tok_it = token_name_map.find(s_token)) == token_name_map.end()) {
					if(r_bad_token_name_set.find(s_token) == r_bad_token_name_set.end()) {
						fprintf(stderr, "warning: unknown token: \'%s\' (ignored)\n", s_token.c_str());
						r_bad_token_name_set.insert(s_token);
					}
					// t_odo - add set of unknown token names so each warning is only printed once

					continue;
				}
				// identify the token, skip unknown tokens

				size_t n_token_type = (*p_tok_it).second;
				CPrimitiveFire dispatcher(n_line_no, s_line, s_token, r_callback);
				if(!CTypelistItemSelect<_TyParsePrimitiveTypelist,
				   CPrimitiveFire>::Select(n_token_type, dispatcher)) {
					fclose(p_fr);
					return false;
				}
				// pass the token to the appropriate implementation, parse it and dispatch it to the callback
			}
		} catch(std::bad_alloc&) {
			fclose(p_fr);
			return false;
		}

		return true;
	}
};

#include "slam/ParsePrimitives.h"

#endif // __GRAPH_PARSER_INCLUDED