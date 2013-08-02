/*
								+-----------------------------------+
								|                                   |
								|        ***  SLAM Main  ***        |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2013  |
								|                                   |
								|              Main.h               |
								|                                   |
								+-----------------------------------+
*/

#pragma once
#ifndef __SLAMPP_MAIN_INCLUDED
#define __SLAMPP_MAIN_INCLUDED

/**
 *	@file include/slam/Main.h
 *	@brief contains some of the common classes required in main() and the documentation pages
 *	@author -tHE SWINe-
 *	@date 2013-06-14
 */

#include "slam/Config.h"

/**
 *	@mainpage SLAM ++
 *
 *  @section intro_sec Introduction
 *
 *	SLAM ++ is a fast nonlinear-optimization package. At the moment, it is not fully
 *	featured yet, but it is very fast in both batch and incremental modes.
 *
 *	@section install_sec Building
 *
 *	There is a CMakeFile. To be able to change the code and commit back
 *	to the svn, do an out-of-source build, like this:
 *
 *	@code{.sh}
 *	$ mkdir slam
 *	$ cd slam
 *	$ svn checkout svn://svn.code.sf.net/p/slam-plus-plus/code/trunk .
 *	$ cd build
 *	$ cmake -i ..
 *	@endcode
 *
 *	This will configure the project (user input is required: press enter many times, real quick).
 *	One interesting option is to specify default linear solver. In Windows CHOLMOD is the fastest,
 *	in Linux CSparse is even faster.
 *
 *	@code{.sh} $ make @endcode
 *
 *	Will finish the building process. You can also use fast parallel build, like this:
 *
 *	@code{.sh}
 *	$ make -j 8
 *	@endcode
 *
 *	You should now be able to run by typing:
 *
 *	@code{.sh} $ ../bin/SLAM_plus_plus --help @endcode
 *
 *	@section use_sec Usage
 *
 *	The simplest use case is:
 *
 *	@code{.sh} $ ../bin/SLAM_plus_plus -i ../data/manhattanOlson3500.txt --pose-only @endcode
 *
 *	Which should reveal something like this (also depending on the actual configuration):
 *
 *	@code{.txt}
 *	General use:
 *		./SLAM_plus_plus -i <filename> --no-detailed-timing
 *
 *	To run the pose-only datasets more quickly:
 *		./SLAM_plus_plus -i <filename> --pose-only --no-detailed-timing
 *
 *	To run incrementally:
 *		./SLAM_plus_plus -lsp <optimize-each-N-steps> -i <filename> --no-detailed-timing
 *
 *	This generates initial.txt and initial.tga, a description and image of the
 *	system before the final optimization, and solution.txt and solution.tga, a
 *	description and image of the final optimized system (unless --no-bitmaps
 *	is specified).
 *
 *	--help|-h         displays this help screen
 *	--verbose|-v      displays verbose output while running (may slow down,
 *					  especially in windows and if running incrementally)
 *	--silent|-s       suppresses displaying verbose output
 *	--no-show|-ns     doesn't show output image (windows only)
 *	--no-commandline|-nc    doesn't echo command line
 *	--no-flags|-nf    doesn't show compiler flags
 *	--no-detailed-timing    doesn't show detailed timing breakup (use this, you'll
 *					  get confused)
 *	--no-bitmaps      doesn't write bitmaps initial.tga and solution.tga (neither
 *					  the text files)
 *	--pose-only|-po   enables optimisation for pose-only slam (will warn and ignore
 *					  on datasets with landmarks (only the first 1000 lines checked
 *					  in case there are landmarks later, it would segfault))
 *	--use-old-code|-uogc    uses the old CSparse code (no block matrices in it)
 *	--a-slam|-A       uses A-SLAM
 *	--lambda|-,\      uses lambda-SLAM (default, preferred batch solver)
 *	--l-slam|-L       uses L-SLAM
 *	--fast-l-slam|-fL uses the new fast L-SLAM solver (preferred incremental solver)
 *	--use-schur|-us   uses Schur complement to accelerate linear solving
 *	--infile|-i <filename>    specifies input file <filename>; it can cope with
 *					  many file types and conventions
 *	--parse-lines-limit|-pll <N>    sets limit of lines read from the input file
 *					  (handy for testing), note this does not set limit of vertices
 *					  nor edges!
 *	--linear-solve-period|-lsp <N>    sets period for incrementally running linear
 *					  solver (default 0: disabled)
 *	--nonlinear-solve-period|-nsp <N>    sets period for incrementally running
 *					  non-linear solver (default 0: disabled)
 *	--max-nonlinear-solve-iters|-mnsi <N>    sets maximal number of nonlinear
 *					  solver iterations (default 10)
 *	--nonlinear-solve-error-thresh|-nset <f>    sets nonlinear solve error threshold
 *					  (default 20)
 *	--max-final-nonlinear-solve-iters|-mfnsi <N>    sets maximal number of final
 *					  optimization iterations (default 5)
 *	--final-nonlinear-solve-error-thresh|-fnset <f>    sets final nonlinear solve
 *					  error threshold (default .01)
 *	--run-matrix-benchmarks|-rmb <benchmark-name> <benchmark-type>    runs block
 *					  matrix benchmarks (benchmark-name is name of a folder with
 *					  UFLSMC benchmark, benchmark-type is one of alloc, factor, all)
 *	--run-matrix-unit-tests|-rmut    runs block matrix unit tests
 *	@endcode
 *
 *	This seemed to work on Ubuntu (x64 or x86 alike), and on Mac.
 *
 *	In case you are working on Windows, you can conveniently use CMakeGUI and
 *	generate solution for Visual Studio (Graph\@FIT members can also get pre-configured
 *	workspace from Matylda).
 *
 *	@section adv_sec Advanced
 *
 *	Some advanced topics are covered in the following pages:
 *	* \ref ownsolvers shows how to implement solver for a custom problem
 */

/**
 *	@page ownsolvers Implementing own solvers
 *
 *	In order to implement solver for your problem of desired dimensionality and
 *	with appropriate jacobians, one needs to execute the following steps:
 *
 *	* implement a new parser code if required (the parser now only processes 2D and 3D types)
 *		* note that this might not be required if you are going to use your own means of passing data to the solver
 *	* create a new header file for the new code
 *	* implement new vertex types (contain "smart" plus code)
 *	* implement new edge types (contain jacobians code)
 *	* implement edge type traits to connect to parse loop if using the builtin parser, or implement your own parse loop
 *	* write specialization of a system that would use your new types
 *	* add a new branch in Main.cpp to call your new code, or call the optimizer yourself
 *
 *	@section sec0 Implement a new parser code
 *
 *	In case the edge types can not be parsed using the existing code, you need to implement parser for
 *	them. Note that this might not be required if you are going to use your own means of passing data
 *	to the solver (such as e.g. from a video stream if implementing visual SLAM).
 *
 *	Nevertheless if you are going to use the built-in parser, you need to implement storage structure
 *	for your new edge, let's start by copying the CParseEntity_XYT_Edge_2D type that can be found in Parser.h
 *	in the class CParserBase (copy it just below it):
 *
 *	@code{.cpp}
 *	struct CParseEntity_XYT_Edge_2D : public CParseEntity { // change name here (and also below)
 *	    int m_n_node_0;
 *	    int m_n_node_1;
 *	    Eigen::Vector3d m_v_delta;
 *	    Eigen::Matrix3d m_t_inv_sigma; // don't change names of variables unless you must, just change the types from 3d to whatever you need
 *
 *	    inline CParseEntity_XYT_Edge_2D(int n_node_0, int n_node_1,
 *	        double f_delta_x, double f_delta_y, double f_delta_theta,
 *	        const double *p_upper_matrix_3x3) // you can change the parameters as you wish
 *	        :m_n_node_0(n_node_0), m_n_node_1(n_node_1), // fill the indices of edge vertices (these should be zero-based)
 *	        m_v_delta(f_delta_x, f_delta_y, f_delta_theta) // fill the measurement vector
 *	    {
 *	        m_t_inv_sigma <<
 *	            p_upper_matrix_3x3[0], p_upper_matrix_3x3[1], p_upper_matrix_3x3[2],
 *	            p_upper_matrix_3x3[1], p_upper_matrix_3x3[3], p_upper_matrix_3x3[4],
 *	            p_upper_matrix_3x3[2], p_upper_matrix_3x3[4], p_upper_matrix_3x3[5]; // fill the inverse sigma matrix
 *	    }
 *
 *	    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // this is imposed by the use of eigen, copy it
 *	};
 *	@endcode
 *
 *	Once you have parsed type, you need to implement code to parse it. Start by making a copy
 *	of CIgnoreParsePrimitive in ParsePrimitives.h:
 *
 *	@code{.cpp}
 *	class MyNewEdgeParsePrimitive {
 *	public:
 *		static void EnumerateTokens(std::map<std::string, int> &r_token_name_map,
 *			int n_assigned_id) // throw(std::bad_alloc)
 *		{
 *			r_token_name_map["MY-NEW-TOKEN-NAME-IN-UPPERCASE"] = n_assigned_id;
 *			r_token_name_map["MY-OTHER-NEW-TOKEN-NAME-IN-UPPERCASE"] = n_assigned_id;
 *			// associate all the tokens you want to handle in this parse primitive
 *		}
 *
 *		template <class _TyParseLoop>
 *		static bool Parse_and_Dispatch(size_t n_line_no, const std::string &r_s_line,
 *			const std::string &UNUSED(r_s_token), _TyParseLoop &r_parse_loop)
 *		{
 *			// here, a primitive of type r_s_token should be parsed from r_s_line
 *			// and if successful, passed to r_parse_loop
 *
 *			return true;
 *		}
 *	};
 *	@endcode
 *
 *	You will need to associate the token names inside EnumerateTokens(). Note that the names should
 *	all point to the same primitive (e.g. it is ok to put "EDGE:SE2" and "EDGE2", but not so ok to
 *	put "EDGE" and "LANDMARK" together as these carry different data). If handling of more types
 *	is required, more parse primitives need to be written.
 *
 *	In the second function, you need to implement the actual parsing and passing of the parsed
 *	structure to the parse loop for processing. Add your custom code for the parsing itself:
 *
 *	@code{.cpp}
 *	template <class _TyParseLoop>
 *	static bool Parse_and_Dispatch(size_t n_line_no, const std::string &r_s_line,
 *		const std::string &UNUSED(r_s_token), _TyParseLoop &r_parse_loop)
 *	{
 *	    int p_pose_idx[2];
 *	    double p_measurement[3];
 *	    double p_matrix[6];
 *	    // data that are stored on the line (change the measurement and matrix dimensions as required)
 *
 *	    if(sscanf(r_s_line.c_str(), "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", // modify the numbers of %lf to reflect number of numbers you need to read from the line
 *	       p_pose_idx, p_pose_idx + 1, p_measurement, p_measurement + 1, p_measurement + 2,
 *	       p_matrix, p_matrix + 1, p_matrix + 2, p_matrix + 3, p_matrix + 4, // put all the matrix and measurement elements, in the order they come
 *	       p_matrix + 5) != 2 + 3 + 6) { // modify the condition
 *	        _ASSERTE(n_line_no < SIZE_MAX);
 *	        fprintf(stderr, "error: line " PRIsize ": line is truncated\n", n_line_no + 1);
 *	        return false;
 *	    }
 *	    // read the individual numbers
 *
 *	    CParseEntity_XYT_Edge_2D edge(p_pose_idx[0], p_pose_idx[1],
 *	        p_measurement[0], p_measurement[1], p_measurement[2], p_matrix); // this creates an object of type that you introduced to Parser.h (rename it)
 *	    // process the measurement
 *
 *	    r_parse_loop.AppendSystem(edge);
 *	    // append the measurement to the system, or something
 *
 *		return true;
 *	}
 *	@endcode
 *
 *	This adds handling of the new edge type. If you want to extend the code and make your new type
 *	one of the default handled types, you need to add your new type to the list of standard types
 *	at the end of ParsePrimitives.h:
 *
 *	@code{.cpp}
 *	typedef MakeTypelist_Safe((CEdge2DParsePrimitive, CLandmark2DParsePrimitive,
 *		CEdge3DParsePrimitive, CVertex2DParsePrimitive, CVertex3DParsePrimitive,
 *		CIgnoreParsePrimitive, MyNewEdgeParsePrimitive)) CStandardParsedPrimitives;
 *	@endcode
 *
 *	If you did that, you also need to change the CParserBase::CParserAdaptor class,
 *	back in Parser.h to include your new type(s). Create a variant of the CParserBase::CParserAdaptor::AppendSystem() function
 *	for every new type added (just copy-paste and rename). This also needs to be done in
 *	TDatasetPeeker in Main.cpp since it inherits from parser adaptor.
 *
 *	This is, however, not always neccessary. If you just want to solve your specific problem,
 *	you can leave CStandardParsedPrimitives and CParserBase::CParserAdaptor as is.
 *
 *	@section sec1 Create a new header file for the new code
 *
 *	To accomplish the second step, create a new file, it should be in the include/slam folder
 *	and it should be named Something_Types.h, where "Something" is replaced by something that makes
 *	more sense (e.g. SE2 types are stored in SE2_Types.h). The file should be blank and it should contain
 *	the following include directives:
 *
 *	@code{.cpp}
 *	//   @file include/slam/Something_Types.h // change this to your name
 *	//   @author <your name here>
 *	//   @date 2012
 *	//   @brief <brief description of what types are being defined here and what problem it can solve>
 *
 *	#ifndef __SOMETHING_TYPES_INCLUDED
 *	#define __SOMETHING_TYPES_INCLUDED // change this to your name
 *
 *	#include "slam/BaseTypes.h" // this needs to be included for base type implementations, required by the solvers
 *	#include "slam/SE2_Types.h" // include this as well to have some common classes (type traits)
 *
 *	// you can start implementing stuff here
 *
 *	#endif // __SOMETHING_TYPES_INCLUDED
 *	@endcode
 *
 *	Note that modifications to CMakeFiles.txt are not required since it is only a header file.
 *
 *	Also note that for SE2 and SE3 types there was the convention of making "xDSolverBase.h" file (such as 2DSolverBase.h and 3DSolverBase.h)
 *	with some common functions. This is not completely neccessary as it was only required by the legacy code
 *	and this file is not really needed for anything anymore.
 *
 *	@section sec2 Implement new vertex types
 *
 *	In order to start implementing your types, it is best to begin by exploring of what
 *	is in slam/SE2_Types.h. Note the CVertexPose2D class:
 *
 *	@code{.cpp}
 *	class CVertexPose2D : public CSEBaseVertexImpl<CVertexPose2D, 3> { // this says to use base vertex implementation for class with name CVertexPose2D, while the vertex has 3 dimensions; this will generate member variable m_v_state ("member vector" state), which will be Eigen dense column vector with the given number of dimensions
 *	public:
 *	    __SE2_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, copy this
 *
 *	    inline CVertexPose2D() // copy this
 *	    {}
 *
 *	    inline CVertexPose2D(const Eigen::Vector3d &r_v_state) // copy this, change the dimension of the vector to appropriate
 *	        :CSEBaseVertexImpl<CVertexPose2D, 3>(r_v_state) // change the dimension here as well
 *	    {}
 *
 *	    inline void Operator_Plus(const Eigen::VectorXd &r_v_delta) // "smart" plus
 *	    {
 *	        m_v_state += r_v_delta.segment<3>(m_n_order); // pick part of the delta vector, belonging to this vertex, apply +
 *	        m_v_state(2) = CBase2DSolver::C2DJacobians::f_ClampAngle_2Pi(m_v_state(2)); // clamp angle
 *	    }
 *
 *	    inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
 *	    {
 *	        m_v_state -= r_v_delta.segment<3>(m_n_order); // pick part of the delta vector, belonging to this vertex, apply -
 *	        m_v_state(2) = CBase2DSolver::C2DJacobians::f_ClampAngle_2Pi(m_v_state(2)); // clamp angle
 *	    }
 *	};
 *	@endcode
 *
 *	Basically what you need to do is to change the name of the class (don't delete this one,
 *	copy it to your Something_Types.h and edit it there), set the number of dimensions and
 *	implement Operator_Plus(). You can leave operator minus empty, it is not really used,
 *	or at least do this:
 *
 *	@code{.cpp}
 *	    inline void Operator_Minus(const Eigen::VectorXd &r_v_delta) // "smart" minus
 *	    {
 *	        Operator_Plus(-r_v_delta); // call plus with negative delta, that should do the trick
 *	    }
 *	@endcode
 *
 *	That gives you your new shiny vertex type. Once you have that, you need to implement
 *	the edge type(s).
 *
 *	@section sec3 Implement new edge types
 *
 *	Now since the parser is ready, you can implement your edge type; let's begin by copying
 *	CEdgePose2D from SE2_Types.h to your Something_Types.h:
 *
 *	@code{.cpp}
 *	class CEdgePose2D : public CSEBaseEdgeImpl<CEdgePose2D, CVertexPose2D, CVertexPose2D, 3> { // again, this tells that base implementation for base edge for type that will be called
 *	CEdgePose2D, and it will be an edge between two vertices of type CVertexPose2D, and the measurement will have 3 dimensions (in order to have hyperedges, you need to provide your own
 *	version of CSEBaseEdgeImpl, that is an advanced topic)
 *	public:
 *	    class CRelative_to_Absolute_XYT_Initializer { // this is an object which is used to lazy initialize vertices (copy it)
 *	    protected:
 *	        const Eigen::Vector3d &m_r_v_pose1; // this is a reference to the first vertex state (change the dimensions if required)
 *	        const CParserBase::CParseEntity_XYT_Edge_2D &m_r_edge; // this is a reference to parser entity, containing the current edge (change the type to your type)
 *
 *	    public:
 *	        inline CRelative_to_Absolute_XYT_Initializer(const Eigen::Vector3d &r_v_vertex1,
 *	            const CParserBase::CParseEntity_XYT_Edge_2D &r_edge) // just change the types, same as above
 *	            :m_r_v_pose1(r_v_vertex1), m_r_edge(r_edge)
 *	        {}
 *
 *	        inline operator CVertexPose2D() const // this function calculates initial prior from the state of the first vertex m_r_v_pose1 and from the edge measurement m_r_edge
 *	        {
 *	            Eigen::Vector3d v_pose2;
 *	            CBase2DSolver::C2DJacobians::Relative_to_Absolute(m_r_v_pose1, m_r_edge.m_v_delta, v_pose2); // implement your own equation here
 *	            return CVertexPose2D(v_pose2);
 *	        }
 *	    };
 *
 *	public:
 *	    __SE2_TYPES_ALIGN_OPERATOR_NEW // imposed by the use of eigen, just copy this
 *
 *	    inline CEdgePose2D() // copy this
 *	    {}
 *
 *	    template <class CSystem>
 *	    CEdgePose2D(const CParserBase::CParseEntity_XYT_Edge_2D &r_t_edge, CSystem &r_system) // this is edge constructor how it is called in the parse loop; you need to change type of the edge
 *	        :CSEBaseEdgeImpl<CEdgePose2D, CVertexPose2D, CVertexPose2D, 3>(r_t_edge.m_n_node_0,
 *	        r_t_edge.m_n_node_1, r_t_edge.m_v_delta, r_t_edge.m_t_inv_sigma) // this just calls the base edge implementation, you need to change types and maintain the parameters if
 *	        required (these are: index of first and of second vertex, the measurement vector and the inverse sigma matrix)
 *	    {
 *	        m_p_vertex0 = &r_system.template r_Get_Vertex<CVertexPose2D>(r_t_edge.m_n_node_0, CInitializeNullVertex());
 *	        m_p_vertex1 = &r_system.template r_Get_Vertex<CVertexPose2D>(r_t_edge.m_n_node_1, CRelative_to_Absolute_XYT_Initializer(m_p_vertex0->v_State(), r_t_edge)); // rename your initializer if required
 *	        // get vertices (initialize if required)
 *			// the strange syntax with "template" is required by g++, otherwise gives "expected primary-expression before '>' token"
 *
 *	        _ASSERTE(((CSEBaseVertex*)m_p_vertex0)->n_Dimension() == 3);
 *	        _ASSERTE(((CSEBaseVertex*)m_p_vertex1)->n_Dimension() == 3); // change the dimensionality here, if required (full type control currently not possible to make sure that the indices in the edge really point to vertices that are of the expected type. e.g. the edge might expect both vertices to be poses, but if the dataset is not "nice", one of the vertices can very well be a landmark)
 *	        // make sure the dimensionality is correct (might not be)
 *	    }
 *
 *	    inline void Calculate_Jacobians_Expectation_Error(Eigen::Matrix3d &r_t_jacobian0,
 *	        Eigen::Matrix3d &r_t_jacobian1, Eigen::Vector3d &r_v_expectation,
 *	        Eigen::Vector3d &r_v_error) const // change dimensionality of eigen types, if required
 *	    {
 *	        CBase2DSolver::C2DJacobians::Absolute_to_Relative(m_p_vertex0->v_State(),
 *	            m_p_vertex1->v_State(), r_v_expectation, r_t_jacobian0, r_t_jacobian1); // write your jacobian calculation code here (vertex state vectors are inputs, the rest are the outputs
 *	            that you need to fill)
 *	        // calculates the expectation and the jacobians
 *
 *	        r_v_error = m_v_measurement - r_v_expectation;
 *	        r_v_error(2) = CBase2DSolver::C2DJacobians::f_ClampAngularError_2Pi(r_v_error(2)); // write your error calculation code here
 *	        // calculates error (possibly re-calculates, if running A-SLAM)
 *	    }
 *
 *	    inline double f_Chi_Squared_Error() const // this function should mostly work as is, you just need to change dimensions of the vectors and matrices
 *	    {
 *	        Eigen::Matrix3d p_jacobi[2];
 *	        Eigen::Vector3d v_expectation, v_error;
 *	        Calculate_Jacobians_Expectation_Error(p_jacobi[0], p_jacobi[1], v_expectation, v_error);
 *	        // calculates the expectation, error and the jacobians
 *
 *	        return (v_error.transpose() * m_t_sigma_inv).dot(v_error); // ||h_i(O_i) - z_i||^2 lambda_i
 *	    }
 *	};
 *	@endcode
 *
 *	@section sec4 Implement edge type traits to connect to parse loop
 *
 *	This chapter only applies if using the built-in parser (otherwise see the \ref sec5 "next chapter").
 *
 *	And that is it, you have an edge! At this point, it is still a little bit to do. In case
 *	you will be using the built-in parser, you need to hook up the edges with the parse loop.
 *	For every edge type you've created, you need to write trait class specialization. These
 *	traits help the parse loop decide which edge type to create from which parsed structure
 *	(or whether to ignore / fail):
 *
 *	@code{.cpp}
 *	template <class CParsedStructure>
 *	class CMyNewEdgeTraits { // replace "MyNew" with something (e.g. SE2 types have CSE2EdgeTraits)
 *	    typedef CFailOnEdgeType _TyEdge; // it should fail on unknown edge types
 *	    static const char *p_s_Reason()
 *	    { return "unknown edge type occured"; }
 *	};
 *
 *	template <>
 *	class CMyNewEdgeTraits<CParserBase::TOdometry2D> { // this is a trait for *parsed* type CParserBase::TOdometry2D
 *	public:
 *	    typedef CEdgePose2D _TyEdge; // here we say to create CEdgePose2D type
 *	};
 *
 *	template <>
 *	class CSE2OnlyPoseEdgeTraits<CParserBase::TLandmark2D> { // an example, don't write "CSE2OnlyPoseEdgeTraits", that is already defined
 *	public:
 *	    typedef CFailOnEdgeType _TyEdge; // it should fail on this edge type
 *
 *	    static const char *p_s_Reason()
 *	    {
 *	        return "landmark edges not permitted in pose-only solver"; // tell us why
 *	    }
 *	};
 *	@endcode
 *
 *	You need to fill traits with all the types you want to permit in the system. You can add
 *	also the edges which you don't want to handle to provide more meaningful error messages
 *	than "unknown edge type occured", but it is not mandatory.
 *
 *	With all that, the parser knows how to parse your new edge types, the parse loop knows
 *	how to redirect those to the system and how to call optimizer.
 *
 *	@section sec5 Implement your own parse loop
 *
 *	If you don't want to use
 *	the built-in parser, look inside ParseLoop.h and see how the edges are passed to the
 *	system for optimization. You need to write your own code that does the same. Basically
 *	you need to:
 *
 *	@code{.cpp}
 *	while(have more edges) {
 *	    TParsedType edge_in;
 *	    // fill this with the new measurement somehow
 *
 *	    CYourEdgeType &r_internal_edge_rep = system.r_Add_Edge(CYourEdgeType(edge_in, system));
 *	    // add the edge to the system ("convert" parsed edge to internal representation)
 *
 *	    solver.Incremental_Step(r_internal_edge_rep);
 *	    // call solver with the new edge
 *	}
 *	@endcode
 *
 *	Note that at this point you don't have <tt>system</tt> and you don't have <tt>solver</tt>. The next chapters
 *	contain information on where to get those.
 *
 *	@section sec6 Write specialization of a system
 *
 *	Once you have that, finally we need to make a solver to solve using your jacobians and
 *	system to contain your edges and vertices. You write:
 *
 *	@code{.cpp}
 *	typedef MakeTypelist_Safe((CVertexPose2D, CVertexLandmark2D)) TVertexTypelist; // just put your vertex types in the list (note the double parentheses are required)
 *	typedef MakeTypelist_Safe((CEdgePose2D, CEdgePoseLandmark2D)) TEdgeTypelist; // just put your edge types in the list (note the double parentheses are required)
 *	typedef CFlatSystem<CSEBaseVertex, TVertexTypelist,
 *		CSEBaseEdge, TEdgeTypelist> CSystemType; // the Base types can be your types in case you are not using CBaseSEVertexImpl / CBaseSEEdgeImpl
 *	// make a system permitting SE(2) vertex and edge types
 *	@endcode
 *
 *	@section sec7 Calling the new system
 *
 *	@code{.cpp}
 *	CSystemType system;
 *	// declare the system
 *
 *	typedef CLinearSolver_CSparse CLinearSolverType; // use CSparse as linear solver (or any other)
 *	CLinearSolverType linear_solver;
 *	// you need a linear solver
 *
 *	typedef CNonlinearSolver_Lambda CNonlinearSolverType; // or use A, these two are rather stable and we're confident they will give correct results
 *	typedef CNonlinearSolverType<CSystemType, CLinearSolverType> CSpecializedNonlinearSolverType;
 *	CSpecializedNonlinearSolverType nonlinear_solver(system, n_linear_solve_each_n_steps,
 *		n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
 *		f_nonlinear_solve_error_threshold, b_verbose, linear_solver); // initialize nonlinear solver
 *	// prepare nonlinear solver
 *
 *	typedef CMyNewEdgeTraits CEdgeTraitsType; // use any edge traits appropriate
 *	typedef CParseLoop<CSystemType, CSpecializedNonlinearSolverType,
 *	    CEdgeTraitsType> CSpecializedParseLoopType;
 *	CSpecializedParseLoopType parse_loop(system, nonlinear_solver);
 *	// prepare parse loop
 *
 *	typedef MakeTypelist_Safe((MyNewEdgeParsePrimitive)) CMyParsedPrimitives; // list of the new parsed primitives
 *	typedef typename CConcatTypelist<CMyParsedPrimitives,
 *		CStandardParsedPrimitives>::_TyResult CParsePrimitivesList; // list of parsed primitives
 *	// this is only required if you did not modify CStandardParsedPrimitives already
 *	// also, for handling your custom problems, CMyParsedPrimitives might be just enough (no concatenation required)
 *
 *	CParserTemplate<CSpecializedParseLoopType, CParsePrimitivesList> p;
 *	if(!p.Parse(p_s_input_file, &parse_loop, n_max_lines_to_process)) {
 *	    fprintf(stderr, "error: failed to parse input file\n");
 *	    exit(-1);
 *	}
 *	// run the parser
 *
 *	system.Plot2D("initial.tga");
 *	system.Dump("initial.txt");
 *	// save the solution
 *
 *	nonlinear_solver.Optimize(n_max_final_optimization_iteration_num, f_final_optimization_threshold);
 *	// perform the final optimization
 *
 *	system.Plot2D("solution.tga");
 *	system.Dump("solution.txt");
 *	// save the solution
 *	@endcode
 *
 *	And as they say, you have optimized.
 *
 *	@section sec8 Modifying Main.cpp to use your new code via commandline
 *
 *	In case you want to update the SLAM_plus_plus program with your new types then you don't have
 *	to write all this, you just find the variable b_10k_opts and follow it arround, mirroring
 *	everything that is done for it - that should also include declaring new system type with your
 *	list of vertices and edges per every solver type you want to support (you might e.g. want to
 *	write b_bundle_adjustment).
 *
 *	@section sec9 Final Thoughts
 *
 *	There is kind of a lot of duplicate editing in vertex and edge types, some of it might be handled
 *	using some \#defines but i prefer to keep the code clean. It is not so much extra work, just use
 *	"find and replace all", making sure to tick "Case Sensitive" and "Whole Words Only".
 *
 *	If you run in any problems while implementing your SLAM, let me know (to help you and to improve
 *	this tutorial / the code).
 *
 */

/**
 *	@brief peeks at the dataset to detect what's inside
 */
struct TDatasetPeeker : public CParserBase::CParserAdaptor {
	bool b_has_odometry; /**< @brief set, if 2D odometry token was found in the file */
	bool b_has_landmark; /**< @brief set, if 2D landmark token was found in the file */
	bool b_has_edge2d; /**< @brief set, if 2D edge (odometry) token was found in the file */
	bool b_has_vertex; /**< @brief set, if 2D vertex token was found in the file */
	bool b_has_edge3d; /**< @brief set, if 3D edge (odometry) token was found in the file */
	bool b_has_vertex3d; /**< @brief set, if 3D vertex token was found in the file */
	bool b_has_ba; /**< @brief set, if BA token was found in the file */

	/**
	 *	@brief peeks at the dataset
	 */
	TDatasetPeeker(const char *p_s_input_file, size_t n_max_lines_to_process = 0)
		:b_has_odometry(false), b_has_landmark(false), b_has_edge2d(false),
		b_has_vertex(false), b_has_edge3d(false), b_has_vertex3d(false), b_has_ba(false)
	{
		CStandardParser p;
		if(!p.Parse(p_s_input_file, *this, (n_max_lines_to_process > 0)?
		   n_max_lines_to_process : 1000)) {
			fprintf(stderr, "warning: failed to peek-parse input file\n");
			/*b_has_odometry = true;
			b_has_landmark = true;
			b_has_edge2d = true;
			b_has_vertex = true;*/
			// don't do that, it will fail anyway and this only makes up for the --pose-only error
			// which will confuse everyone in this case
		}
	}

	/**
	 *	@brief appends the system with an odometry measurement
	 *	@param[in] r_t_odo is the measurement to be appended
	 */
	virtual void AppendSystem(const CParserBase::TEdge2D &UNUSED(r_t_odo))
	{
		b_has_odometry = true;
	}

	/**
	 *	@brief appends the system with a landmark measurement
	 *	@param[in] r_t_landmark is the measurement to be appended
	 */
	virtual void AppendSystem(const CParserBase::TLandmark2D &UNUSED(r_t_landmark))
	{
		b_has_landmark = true;
	}

	/**
	 *	@brief appends the system with vertex position
	 *	@param[in] r_t_vertex is the vertex to be appended
	 *	@note The vertices can be ignored in most of the solvers.
	 */
	virtual void AppendSystem(const CParserBase::TVertex2D &UNUSED(r_t_vertex))
	{
		b_has_vertex = true;
	}

	/**
	 *	@brief appends the system with an odometry measurement
	 *	@param[in] r_t_edge is the measurement to be appended
	 */
	virtual void AppendSystem(const CParserBase::TEdge3D &UNUSED(r_t_edge))
	{
		b_has_edge3d = true;
	}

	/**
	 *	@brief appends the system with vertex position
	 *	@param[in] r_t_vertex is the vertex to be appended
	 *	@note The vertices can be ignored in most of the solvers.
	 */
	virtual void AppendSystem(const CParserBase::TVertex3D &UNUSED(r_t_vertex))
	{
		b_has_vertex3d = true;
	}

	/**
	 *	@brief appends the system with vertex position
	 *	@param[in] r_t_vertex is the vertex to be appended
	 *	@note The vertices can be ignored in most of the solvers.
	 */
	virtual void AppendSystem(const CParserBase::TVertexXYZ &UNUSED(r_t_vertex))
	{
		b_has_ba = true;
	}

	/**
	 *	@brief appends the system with camera position and parameters
	 *	@param[in] r_t_vertex is the vertex to be appended
	 *	@note The vertices can be ignored in most of the solvers.
	 */
	virtual void AppendSystem(const CParserBase::TVertexCam3D &UNUSED(r_t_vertex))
	{
		b_has_ba = true;
	}

	/**
	 *	@brief appends the system with a projection measurement
	 *	@param[in] r_t_edge is the measurement to be appended
	 */
	virtual void AppendSystem(const CParserBase::TEdgeP2C3D &UNUSED(r_t_edge))
	{
		b_has_ba = true;
	}
};

/**
 *	@brief a helper template for running tests with given type configuration (to avoid repeating the same code)
 *
 *	@tparam CSystemType is system type (derived from CFlatSystem)
 *	@tparam CNonlinearSolverType is nonlinear solver template name
 *	@tparam CEdgeTraitsType is edge traits template name
 *	@tparam CParseLoopType is parse loop template name
 */
template <class CSystemType,
	template <class, class, class> class CNonlinearSolverType,
	template <class> class CEdgeTraitsType,
	template <class, class, template <class> class> class CParseLoopType>
class CTester {
public:
#ifdef __USE_NATIVE_CHOLESKY
	typedef typename CSystemType::_TyJacobianMatrixBlockList _TyAMatrixBlockSizes; /** @brief list of jacobian block sizes */
	typedef typename __fbs_ut::CBlockSizesAfterPreMultiplyWithSelfTranspose<
		_TyAMatrixBlockSizes>::_TyResult _TyLambdaMatrixBlockSizes; /**< @brief possible block matrices, found in lambda and L */

	typedef CLinearSolver_UberBlock<_TyLambdaMatrixBlockSizes> CLinearSolverType; /**< @brief linear solver type (native) */
#elif defined(__USE_CHOLMOD)
	typedef CLinearSolver_CholMod CLinearSolverType; /**< @brief linear solver type (CHOLMOD) */
#else // __USE_NATIVE_CHOLESKY
#ifdef __USE_CXSPARSE
	typedef CLinearSolver_CXSparse CLinearSolverType; /**< @brief linear solver type (CXSparse) */
#else // __USE_CXSPARSE
	typedef CLinearSolver_CSparse CLinearSolverType; /**< @brief linear solver type (CSparse) */
#endif // __USE_CXSPARSE
#endif // __USE_NATIVE_CHOLESKY

public:
	/**
	 *	@brief runs optimization of system stored in a file, measures timing
	 *
	 *	@param[in] p_s_input_file is the name of input file
	 *	@param[in] n_max_lines_to_process is limit of parsed lines (0 means unlimited)
	 *	@param[in] n_linear_solve_each_n_steps is period for calculating
	 *		linear solution to the system (0 means never)
	 *	@param[in] n_nonlinear_solve_each_n_steps is period for calculating
	 *		nonlinear solution to the system (0 means never)
	 *	@param[in] n_max_nonlinear_solve_iteration_num is maximal number
	 *		of iterations spent in calculating nonlinear solution to the system
	 *	@param[in] f_nonlinear_solve_error_threshold is error threshold for early stopping
	 *		the incremenatl calculation of nonlinear solution to the system
	 *	@param[in] n_max_final_optimization_iteration_num is maximal number of iterations
	 *		spent in calculating nonlinear solution to the system at the end
	 *	@param[in] f_final_optimization_threshold is error threshold for early stopping
	 *		the final calculation of nonlinear solution to the system
	 *	@param[in] b_verbose is verbosity flag
	 *	@param[in] b_use_schur is Schur complement flag
	 *	@param[in] b_show_detailed_timing is detailed timing flag
	 *	@param[in] b_write_bitmaps is flag for writing (2D) bitmaps of initial and optimized system
	 *
	 *	@return Returns true on success, false on failure.
	 */
	static inline bool Run_and_Shout(const char *p_s_input_file,
		size_t n_max_lines_to_process, size_t n_linear_solve_each_n_steps,
		size_t n_nonlinear_solve_each_n_steps, size_t n_max_nonlinear_solve_iteration_num,
		double f_nonlinear_solve_error_threshold, size_t n_max_final_optimization_iteration_num,
		double f_final_optimization_threshold, bool b_verbose, bool b_use_schur,
		bool b_show_detailed_timing, bool b_write_bitmaps)
	{
		CSystemType system;

		CLinearSolverType linear_solver;
		// prepare a linear solver

		typedef typename CSystemType::_TyJacobianMatrixBlockList CMatrixTypelist; // use the default from the system
		typedef CNonlinearSolverType<CSystemType, CLinearSolverType, CMatrixTypelist> // have to supply CMatrixTypelist because in CTester template there is no default param for CNonlinearSolverType :(
			CSpecializedNonlinearSolverType;
		CSpecializedNonlinearSolverType nonlinear_solver(system, n_linear_solve_each_n_steps,
			n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
			f_nonlinear_solve_error_threshold, b_verbose, linear_solver, b_use_schur);
		// prepare nonlinear solver

		typedef CParseLoopType<CSystemType, CSpecializedNonlinearSolverType,
			CEdgeTraitsType> CSpecializedParseLoop;
		CSpecializedParseLoop parse_loop(system, nonlinear_solver);
		// prepare parse loop

		if(b_verbose && b_use_schur)
			printf("using Schur complement\n");
		// verbose

		CTimer t;
		t.ResetTimer();
		// start meassuring time

		//printf("n_max_lines_to_process = %d\n", int(n_max_lines_to_process)); // debug
		//CParser p;
		CParserTemplate<CSpecializedParseLoop, CStandardParsedPrimitives> p;
		if(!p.Parse(p_s_input_file, parse_loop, n_max_lines_to_process)) {
			fprintf(stderr, "error: failed to parse input file\n");
			return false;
		}
		// run the parser, solver incremental function is called

		if(!n_linear_solve_each_n_steps && !n_nonlinear_solve_each_n_steps) {
			if(b_verbose)
				fprintf(stderr, "warning: running in batch mode. ignoring time spent in parser\n");
			t.ResetTimer();
		}
		// in case we're running in batch mode, don't count the time
		// spent in parser (because no computation was done)

		double f_time_initial_save_start = t.f_Time();
		if(b_write_bitmaps) {
			system.Plot2D("initial.tga");
			system.Dump("initial.txt");
			// save the initial configuration to a file

			/*double f_error = nonlinear_solver.f_Chi_Squared_Error_Denorm();
			printf("denormalized initial chi2 error: %.2lf\n", f_error);*/
			// need to calculate jacobians & errors first, don't do it ...
		}
		double f_time_initial_save_end = t.f_Time();

		nonlinear_solver.Optimize(n_max_final_optimization_iteration_num, f_final_optimization_threshold);
		// perform the final optimization

		double f_time = t.f_Time() - (f_time_initial_save_end - f_time_initial_save_start); // don't count saving of the initial system state (rasterizing the image takes some time)
		printf("\ndone. it took " PRItimeprecise " (%f sec)\n", PRItimeparams(f_time), f_time);
		//printf("time / 1.6 = " PRItimeprecise "\n", PRItimeparams(f_time / 1.6));
		// display time it took

		if(b_show_detailed_timing) {
			nonlinear_solver.Dump(f_time);
			double f_error = nonlinear_solver.f_Chi_Squared_Error_Denorm();
			printf("denormalized chi2 error: %.2f\n", f_error);
		}

		if(b_write_bitmaps) {
			/*if(!nonlinear_solver.Dump_SystemMatrix("g:\\system_matrix.tga", 5) &&
			   !nonlinear_solver.Dump_SystemMatrix("g:\\system_matrix.tga", 4) &&
			   !nonlinear_solver.Dump_SystemMatrix("g:\\system_matrix.tga", 3) &&
			   !nonlinear_solver.Dump_SystemMatrix("g:\\system_matrix.tga", 2))
				fprintf(stderr, "error: failed to dump system matrix\n");*/
			system.Plot2D("solution.tga");
			system.Plot2D("solution_print.tga", 2048, 2048, 10, 3, 7, 1, true, false, 10); // print size images
			system.Plot2D("solution_print_landmarks.tga", 2048, 2048, 10, 3, 7, 1, true, false, 10, true); // print size images
			system.Plot2D("solution_print_noticks.tga", 2048, 2048, 0, 0, 7, 1, true, false, 4); // print size images
			//system.Plot3D("solution_print3d.tga", 2048, 2048, 10, 3, 7, 1, true, false, 10); // print size images
			//system.Plot3D("solution_print_landmarks3d.tga", 2048, 2048, 10, 3, 7, 1, true, false, 10, true); // print size images
			//system.Plot3D("solution_print_noticks3d.tga", 2048, 2048, 0, 0, 7, 1, true, false, 4); // print size images
			//nonlinear_solver.Save_SystemMatrix_MM("system_matrix.mtx");
			system.Dump("solution.txt");
		}
		// save the solution to a file

		return true;
	}
};

extern int n_dummy_param;

#ifdef __UBER_BLOCK_MATRIX_BENCHMARK_INCLUDED

/**
 *	@brief matrix benchmark runner
 *	@tparam n_block_size is matrix block size for the given benchmark
 */
template <const int n_block_size>
class CBlockBenchRunner {
public:
	/**
	 *	@brief runs a benchmark
	 *
	 *	@param[in] n_iter_num is number of iterations to take
	 *	@param[in] p_s_bench_name is name of a folder with the benchmark
	 *	@param[in] p_s_bench_type is type of benchmark to be run ("alloc", "chol" ro "all")
	 *
	 *	@return Returns true on success, false on failure (benchmark
	 *		failed or benchmars were not compiled).
	 */
	static bool Run(int n_iter_num, const char *p_s_bench_name, const char *p_s_bench_type)
	{
		if(!n_iter_num)
			n_iter_num = 1;

		CBlockMatrixBenchmark bmb(n_iter_num, CBlockMatrixBenchmark::order_RectOutline); // more realistic type of fill as incremental SLAM

		char p_s_infile[256];
#if defined(_WIN32) || defined(_WIN64)
		strcpy(p_s_infile, "G:\\uflsmc\\");
#else // _WIN32 || _WIN64
		strcpy(p_s_infile, "../data/"); // ubuntu
#endif // _WIN32 || _WIN64
		strcat(p_s_infile, p_s_bench_name);
		if(!strcmp(p_s_bench_type, "alloc") || !strcmp(p_s_bench_type, "all")) {
			if(!bmb.template Run_AllocationBenchmark<n_block_size>(p_s_infile)) {
				fprintf(stderr, "error: benchmark failed\n");
				bmb.ShowResults(); // show at least partial results
				return false;
			}
		}
		if(!strcmp(p_s_bench_type, "factor") || !strcmp(p_s_bench_type, "all")) {
			if(!bmb.template Run_FactorizationBenchmark<n_block_size>(p_s_infile)) {
				fprintf(stderr, "error: benchmark failed\n");
				bmb.ShowResults(); // show at least partial results
				return false;
			}
		}
		bmb.ShowResults();
		char p_s_outfile[256];
#ifdef __BLOCK_BENCH_CHOLESKY_USE_AMD
		sprintf(p_s_outfile, "%s_%d_AMD.csv", p_s_bench_name, n_block_size);
#else // __BLOCK_BENCH_CHOLESKY_USE_AMD
		sprintf(p_s_outfile, "%s_%d.csv", p_s_bench_name, n_block_size);
#endif // __BLOCK_BENCH_CHOLESKY_USE_AMD
		if(!bmb.Save_ResultSheet(p_s_outfile))
			fprintf(stderr, "error: i/o error while writing results\n");
#ifdef __BLOCK_BENCH_CHOLESKY_USE_AMD
		sprintf(p_s_outfile, "%s_%d_avg_AMD.csv", p_s_bench_name, n_block_size);
#else // __BLOCK_BENCH_CHOLESKY_USE_AMD
		sprintf(p_s_outfile, "%s_%d_avg.csv", p_s_bench_name, n_block_size);
#endif // __BLOCK_BENCH_CHOLESKY_USE_AMD
		if(!bmb.Save_TestBased_ResultSheet(p_s_outfile))
			fprintf(stderr, "error: i/o error while writing results\n");
		// benchmarks on saved matrices (see if they save properly)

		return true;
	}
};

#endif // __UBER_BLOCK_MATRIX_BENCHMARK_INCLUDED

/**
 *	@brief runs a block matrix benchmark
 *
 *	@param[in] n_block_size is matrix block size for the given benchmark
 *		(one of 1, 2, 3, 4, 5, 6, 8, 10, 15, 16, 20, 25 or 30)
 *	@param[in] p_s_bench_name is benchmark name (not path)
 *	@param[in] p_s_bench_type is benchmark type, one of "factor", "alloc" or "all"
 *
 *	@note This handles block sizes 1 - 30 with certain step, the implementation
 *		is split to BlockBenchImpl0.cpp and BlockBenchImpl1.cpp in order to
 *		make parallel building faster (each file instantiates test templates
 *		with different block sizes).
 *	@note The benchmarks mostly use four block sizes, the specified n x n
 *		and (n + 1) x n, n x (n + 1) and (n + 1) x (n + 1), as happens
 *		in landmark datasets.
 */
int n_Run_BlockBenchmark(int n_block_size,
	const char *p_s_bench_name, const char *p_s_bench_type);

/**
 *	@copydoc n_Run_BlockBenchmark()
 *	@note This handles block sizes 8 - 30.
 */
int n_Run_BlockBenchmark1(int n_block_size,
	const char *p_s_bench_name, const char *p_s_bench_type);

/**
 *	@brief prints commandline help
 */
void PrintHelp();

/**
 *	@brief nonlinear solver type
 */
enum ENonlinearSolverType {
	nlsolver_A, /**< @brief nonlinear solver A */
	nlsolver_Lambda, /**< @brief nonlinear solver lambda */
	nlsolver_LambdaLM, /**< @brief nonlinear solver lambda with Levenberg-Marquardt */
	nlsolver_L, /**< @brief nonlinear solver L */
	nlsolver_FastL, /**< @brief nonlinear progressively reordering solver L */
	nlsolver_SPCG /**< @brief nonlinear solver SPCG */
};

/**
 *	@brief structure, containing values of all the commandline arguments
 */
struct TCommandLineArgs {
	ENonlinearSolverType n_solver_choice; /**< @brief nonlinear solver selector */
	bool b_write_bitmaps;
	bool b_no_show;
	bool b_show_commandline;
	bool b_show_flags;
	bool b_show_detailed_timing;
	bool b_verbose; /**< @brief verbosity flag; true means verbose, false means silent */
	bool b_use_schur;
	bool b_run_matrix_benchmarks;
	bool b_run_matrix_unit_tests;
	bool b_use_old_system;
	bool b_10k_opts;
	bool b_use_SE3; // note this is not overriden in commandline but detected in peek-parsing
	bool b_use_BA; // note this is not overriden in commandline but detected in peek-parsing
	const char *p_s_input_file; /** <@brief path to the data file */
	int n_max_lines_to_process; /** <@brief maximal number of lines to process */
	size_t n_linear_solve_each_n_steps; /**< @brief linear solve period, in steps (0 means disabled) */
	size_t n_nonlinear_solve_each_n_steps; /**< @brief nonlinear solve period, in steps (0 means disabled) */
	size_t n_max_nonlinear_solve_iteration_num; /**< @brief maximal number of iterations in nonlinear solve step */
	double f_nonlinear_solve_error_threshold; /**< @brief error threshold for nonlinear solve */
	size_t n_max_final_optimization_iteration_num; // as many other solvers
	double f_final_optimization_threshold;
	const char *p_s_bench_name;
	const char *p_s_bench_type;

	/**
	 *	@brief selects default values for commandline args
	 */
	void Defaults()
	{
		n_solver_choice = nlsolver_Lambda; /**< @brief nonlinear solver selector */
		// solver selection

		b_write_bitmaps = true;
		b_no_show = false;
		b_show_commandline = true;
		b_show_flags = true;
		b_show_detailed_timing = true;
		b_verbose = true;
		// verbosity

		b_use_schur = false;

		b_run_matrix_benchmarks = false;
		b_run_matrix_unit_tests = false;
		b_use_old_system = false; // t_odo - make this commandline
		b_10k_opts = false;
		b_use_SE3 = false; // note this is not overriden in commandline but detected in peek-parsing
		b_use_BA = false; // note this is not overriden in commandline but detected in peek-parsing

		p_s_input_file = 0; /** <@brief path to the data file */
		n_max_lines_to_process = 0; /** <@brief maximal number of lines to process */

		n_linear_solve_each_n_steps = 0; /**< @brief linear solve period, in steps (0 means disabled) */
		n_nonlinear_solve_each_n_steps = 0; /**< @brief nonlinear solve period, in steps (0 means disabled) */
		n_max_nonlinear_solve_iteration_num = 10; /**< @brief maximal number of iterations in nonlinear solve step */
		f_nonlinear_solve_error_threshold = 20; /**< @brief error threshold for nonlinear solve */
		n_max_final_optimization_iteration_num = 5; // as many other solvers
		f_final_optimization_threshold = .01;
		// optimization mode for slam

		p_s_bench_name = 0;
		p_s_bench_type = "all";
	}

	/**
	 *	@brief parse commandline arguments
	 *
	 *	@param[in] n_arg_num is number of commandline arguments
	 *	@param[in] p_arg_list is the list of commandline arguments
	 *
	 *	@return Returns true on success, false on failure.
	 */
	bool Parse(int n_arg_num, const char **p_arg_list)
	{
		for(int i = 1; i < n_arg_num; ++ i) {
			if(!strcmp(p_arg_list[i], "--help") || !strcmp(p_arg_list[i], "-h")) {
				PrintHelp();
				//fprintf(stderr, "no help for you! mwuhahaha! (please read Main.cpp)\n"); // t_odo
				return false; // quit
			} else if(!strcmp(p_arg_list[i], "--verbose") || !strcmp(p_arg_list[i], "-v"))
				b_verbose = true;
			else if(!strcmp(p_arg_list[i], "--silent") || !strcmp(p_arg_list[i], "-s"))
				b_verbose = false;
			else if(!strcmp(p_arg_list[i], "--use-schur") || !strcmp(p_arg_list[i], "-us"))
				b_use_schur = true;
			else if(!strcmp(p_arg_list[i], "--no-show") || !strcmp(p_arg_list[i], "-ns"))
				b_no_show = true;
			else if(!strcmp(p_arg_list[i], "--no-commandline") || !strcmp(p_arg_list[i], "-nc"))
				b_show_commandline = false;
			else if(!strcmp(p_arg_list[i], "--lambda") || !strcmp(p_arg_list[i], "-,\\"))
				n_solver_choice = nlsolver_Lambda;
			else if(!strcmp(p_arg_list[i], "--precond") || !strcmp(p_arg_list[i], "-spcg"))
				n_solver_choice = nlsolver_SPCG;
			else if(!strcmp(p_arg_list[i], "--no-flags") || !strcmp(p_arg_list[i], "-nf"))
				b_show_flags = false;
			else if(!strcmp(p_arg_list[i], "--run-matrix-unit-tests") || !strcmp(p_arg_list[i], "-rmut"))
				b_run_matrix_unit_tests = true;
			else if(!strcmp(p_arg_list[i], "--no-detailed-timing") || !strcmp(p_arg_list[i], "-ndt"))
				b_show_detailed_timing = false;
			else if(!strcmp(p_arg_list[i], "--use-old-code") || !strcmp(p_arg_list[i], "-uogc"))
				b_use_old_system = true;
			else if(!strcmp(p_arg_list[i], "--no-bitmaps") || !strcmp(p_arg_list[i], "-nb")) {
				b_write_bitmaps = false;
				b_no_show = true; // no bitmaps ... what can it show?
			} else if(!strcmp(p_arg_list[i], "--pose-only") || !strcmp(p_arg_list[i], "-po"))
				b_10k_opts = true;
			else if(!strcmp(p_arg_list[i], "--a-slam") || !strcmp(p_arg_list[i], "-A"))
				n_solver_choice = nlsolver_A;
			else if(!strcmp(p_arg_list[i], "--l-slam") || !strcmp(p_arg_list[i], "-L"))
				n_solver_choice = nlsolver_L;
			else if(!strcmp(p_arg_list[i], "--fast-l-slam") || !strcmp(p_arg_list[i], "-fL"))
				n_solver_choice = nlsolver_FastL;
			else if(i + 1 == n_arg_num) {
				fprintf(stderr, "error: argument \'%s\': missing value or an unknown argument\n", p_arg_list[i]);
				return false;
			} else if(!strcmp(p_arg_list[i], "--infile") || !strcmp(p_arg_list[i], "-i"))
				p_s_input_file = p_arg_list[++ i];
			else if(!strcmp(p_arg_list[i], "--parse-lines-limit") || !strcmp(p_arg_list[i], "-pll"))
				n_max_lines_to_process = atol(p_arg_list[++ i]);
			else if(!strcmp(p_arg_list[i], "--linear-solve-period") || !strcmp(p_arg_list[i], "-lsp"))
				n_linear_solve_each_n_steps = atol(p_arg_list[++ i]);
			else if(!strcmp(p_arg_list[i], "--nonlinear-solve-period") || !strcmp(p_arg_list[i], "-nsp"))
				n_nonlinear_solve_each_n_steps = atol(p_arg_list[++ i]);
			else if(!strcmp(p_arg_list[i], "--max-nonlinear-solve-iters") || !strcmp(p_arg_list[i], "-mnsi"))
				n_max_nonlinear_solve_iteration_num = atol(p_arg_list[++ i]);
			else if(!strcmp(p_arg_list[i], "--nonlinear-solve-error-thresh") || !strcmp(p_arg_list[i], "-nset"))
				f_nonlinear_solve_error_threshold = atof(p_arg_list[++ i]);
			else if(!strcmp(p_arg_list[i], "--max-final-nonlinear-solve-iters") || !strcmp(p_arg_list[i], "-mfnsi"))
				n_max_final_optimization_iteration_num = atol(p_arg_list[++ i]);
			else if(!strcmp(p_arg_list[i], "--final-nonlinear-solve-error-thresh") || !strcmp(p_arg_list[i], "-fnset"))
				f_final_optimization_threshold = atof(p_arg_list[++ i]);
			else if(!strcmp(p_arg_list[i], "--run-matrix-benchmarks") || !strcmp(p_arg_list[i], "-rmb")) {
				if(i + 2 >= n_arg_num) {
					fprintf(stderr, "error: argument \'%s\': missing the second value\n", p_arg_list[i]);
					return false;
				}
				b_run_matrix_benchmarks = true;
				p_s_bench_name = p_arg_list[++ i];
				p_s_bench_type = p_arg_list[++ i];
				if(strcmp(p_s_bench_type, "alloc") &&
				   strcmp(p_s_bench_type, "factor") &&
				   strcmp(p_s_bench_type, "all")) {
					fprintf(stderr, "error: argument \'%s\': unknown benchmark type\n", p_arg_list[i]);
					return false;
				}
			} else if(!strcmp(p_arg_list[i], "--dummy-param") || !strcmp(p_arg_list[i], "-dp"))
				n_dummy_param = atol(p_arg_list[++ i]);
			else {
				fprintf(stderr, "error: argument \'%s\': an unknown argument\n", p_arg_list[i]);
				return false;
			}
		}
		// "parse" cmdline

		return true;
	}
};

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
	enum {
		solver_type_Id = n_solver_type_id /**< @brief nonlinear solver id */
	};

	/**
	 *	@brief runs the application, using this solver
	 *
	 *	@tparam CSystemType is system type (derived from CFlatSystem)
	 *	@tparam CEdgeTraitsType is edge traits template name
	 *	@tparam CParseLoopType is parse loop template name
	 *
	 *	@param[in] t_args is parsed commandline arguments value
	 *
	 *	@return Returns true on success, false on failure.
	 */
	template <class CSystemType, template <class> class CEdgeTraitsType,
		template <class, class, template <class> class vcneedsnamehere> class CParseLoopType>
	static inline bool Run_MainApp(TCommandLineArgs t_args) // throw(std::runtime_error, std::bad_alloc)
	{
		return CTester<CSystemType, CNonlinearSolverType, CEdgeTraitsType, CParseLoopType>::Run_and_Shout(
			t_args.p_s_input_file, t_args.n_max_lines_to_process, t_args.n_linear_solve_each_n_steps,
			t_args.n_nonlinear_solve_each_n_steps, t_args.n_max_nonlinear_solve_iteration_num,
			t_args.f_nonlinear_solve_error_threshold, t_args.n_max_final_optimization_iteration_num,
			t_args.f_final_optimization_threshold, t_args.b_verbose, t_args.b_use_schur,
			t_args.b_show_detailed_timing, t_args.b_write_bitmaps);
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
	enum {
		solver_type_Id = n_solver_type_id /**< @brief nonlinear solver id */
	};

	/**
	 *	@brief runs the application, using this solver
	 *
	 *	@tparam CSystemType is system type (derived from CFlatSystem)
	 *	@tparam CEdgeTraitsType is edge traits template name
	 *	@tparam CParseLoopType is parse loop template name
	 *
	 *	@param[in] t_args is parsed commandline arguments value
	 *
	 *	@return Returns true on success, false on failure.
	 */
	template <class CSystemType, template <class> class CEdgeTraitsType,
		template <class, class, template <class> class vcneedsnamehere> class CParseLoopType>
	static inline bool Run_MainApp(TCommandLineArgs UNUSED(t_args))
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
	enum {
		solver_type_Id = n_solver_type_id /**< @brief nonlinear solver id */
	};

	/**
	 *	@brief runs the application, using this solver
	 *
	 *	@tparam CSystemType is system type (derived from CFlatSystem)
	 *	@tparam CEdgeTraitsType is edge traits template name
	 *	@tparam CParseLoopType is parse loop template name
	 *
	 *	@param[in] t_args is parsed commandline arguments value
	 *
	 *	@return Returns true on success, false on failure.
	 */
	template <class CSystemType, template <class> class CEdgeTraitsType,
		template <class, class, template <class> class vcneedsnamehere> class CParseLoopType>
	static inline bool Run_MainApp(TCommandLineArgs UNUSED(t_args))
	{
		fprintf(stderr, "error: the selected solver is not supported by the SE types\n");
		return false;
	}
};

/**
 *	@brief functor for CTypelistForEach; selects requested solver and runs the application
 *
 *	@tparam CSystemType is system type (derived from CFlatSystem)
 *	@tparam CEdgeTraitsType is edge traits template name
 *	@tparam CParseLoopType is parse loop template name
 */
template <class CSystemType, template <class> class CEdgeTraitsType,
	template <class, class, template <class> class> class CParseLoopType>
class CSolverCaller {
protected:
	TCommandLineArgs m_t_args; /**< @brief copy of parsed commandline args */
	int m_n_result; /**< @brief final result of the application */

public:
	/**
	 *	@brief default constructor
	 *	@param[in] _t_args is a copy of parsed commandline args
	 */
	inline CSolverCaller(TCommandLineArgs t_args)
		:m_t_args(t_args), m_n_result(-1)
	{}

	/**
	 *	@brief function operator; checks for suitable solver and runs
	 *	@tparam CSolverType is a pair of solver id and type
	 */
	template <class CSolverType>
	inline void operator ()() // throw(std::runtime_error, std::bad_alloc)
	{
		if(m_t_args.n_solver_choice == int(CSolverType::solver_type_Id)) {
			m_n_result = (CSolverType::template Run_MainApp<CSystemType, CEdgeTraitsType,
				CParseLoopType>(m_t_args))? 0 : -1;
			// call the traits
		}
		// find the solver in the list
	}

	inline int n_Result() const
	{
		return m_n_result;
	}
};

/**
 *	@brief runs bundle adjustment with a specified solver
 *	@param[in] t_args is a copy of parsed commandline arguments
 *	@return Returns 0 on success, -1 on failure.
 */
int n_Run_BA_Solver(TCommandLineArgs t_args); // throw(std::runtime_error, std::bad_alloc)

/**
 *	@brief runs (pose-only) 3D SLAM with a specified solver
 *	@param[in] t_args is a copy of parsed commandline arguments
 *	@return Returns 0 on success, -1 on failure.
 */
int n_Run_SE3_Solver(TCommandLineArgs t_args); // throw(std::runtime_error, std::bad_alloc)

/**
 *	@brief runs 2D SLAM with a specified solver
 *	@param[in] t_args is a copy of parsed commandline arguments
 *	@return Returns 0 on success, -1 on failure.
 */
int n_Run_SE2_Solver(TCommandLineArgs t_args); // throw(std::runtime_error, std::bad_alloc)

/**
 *	@brief runs pose-only 2D SLAM with a specified solver
 *	@param[in] t_args is a copy of parsed commandline arguments
 *	@return Returns 0 on success, -1 on failure.
 */
int n_Run_SE2PoseOnly_Solver(TCommandLineArgs t_args); // throw(std::runtime_error, std::bad_alloc)

#endif // __SLAMPP_MAIN_INCLUDED
