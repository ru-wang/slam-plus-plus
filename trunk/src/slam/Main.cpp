/*
								+-----------------------------------+
								|                                   |
								|        ***  SLAM Main  ***        |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2012  |
								|                                   |
								|             Main.cpp              |
								|                                   |
								+-----------------------------------+
*/

/**
 *	@file src/slam/Main.cpp
 *	@brief contains the main() function of the SLAM program
 *	@author -tHE SWINe-
 *	@date 2012-03-30
 */

int n_dummy_param = 0;
// a dummy parameter, used as a convenient commandline input, intended for debugging / testing

#include "slam/Main.h"

/**
 *	@brief prints all the important compiler / optimization switches this app was built with
 */
static void DisplaySwitches()
{
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
	printf("SLAM++ version x64 (compiled at %s)\nbuilt with the following flags:\n", __DATE__);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
	printf("SLAM++ version x86 (compiled at %s)\nbuilt with the following flags:\n", __DATE__);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64

#ifdef _DEBUG
	printf("%s\n", "_DEBUG");
#endif // _DEBUG
#ifdef _OPENMP
#pragma omp parallel
#pragma omp master
	{
		printf("%s (%d threads)\n", "_OPENMP", omp_get_num_threads());
	}
#endif // _OPENMP

#ifdef __USE_NATIVE_CHOLESKY
	printf("%s\n", "__USE_NATIVE_CHOLESKY");
#elif defined(__USE_CHOLMOD)
	printf("%s\n", "__USE_CHOLMOD");
#ifdef __CHOLMOD_SHORT
	printf("%s\n", "__CHOLMOD_SHORT");
#endif // __CHOLMOD_SHORT
#else // __USE_CHOLMOD
#ifdef __USE_CXSPARSE
	printf("%s\n", "__USE_CXSPARSE");
#ifdef __CXSPARSE_SHORT
	printf("%s\n", "__CXSPARSE_SHORT");
#endif // __CXSPARSE_SHORT
#else // __USE_CXSPARSE
	printf("%s\n", "__USE_CSPARSE");
#endif // __USE_CXSPARSE
#endif // __USE_CHOLMOD

#ifdef GPU_BLAS
	printf("%s\n", "GPU_BLAS");
#endif // GPU_BLAS
#ifdef EIGEN_VECTORIZE
	printf("%s\n", "EIGEN_VECTORIZE");
#endif // EIGEN_VECTORIZE
#ifdef EIGEN_VECTORIZE_SSE
	printf("%s\n", "EIGEN_VECTORIZE_SSE");
#endif // EIGEN_VECTORIZE_SSE
#ifdef EIGEN_VECTORIZE_SSE2
	printf("%s\n", "EIGEN_VECTORIZE_SSE2");
#endif // EIGEN_VECTORIZE_SSE2
#ifdef EIGEN_VECTORIZE_SSE3
	printf("%s\n", "EIGEN_VECTORIZE_SSE3");
#endif // EIGEN_VECTORIZE_SSE3

#ifdef __BLOCK_BENCH_DUMP_MATRIX_IMAGES
	printf("%s\n", "__BLOCK_BENCH_DUMP_MATRIX_IMAGES");
#endif // __BLOCK_BENCH_DUMP_MATRIX_IMAGES
#ifdef __BLOCK_BENCH_BLOCK_TYPE_A
	printf("%s\n", "__BLOCK_BENCH_BLOCK_TYPE_A");
#else // __BLOCK_BENCH_BLOCK_TYPE_A
	printf("%s\n", "__BLOCK_BENCH_BLOCK_TYPE_B"); // nonexistent, but want to know
#endif // __BLOCK_BENCH_BLOCK_TYPE_A
#ifdef __BLOCK_BENCH_CHOLESKY_USE_AMD
	printf("%s\n", "__BLOCK_BENCH_CHOLESKY_USE_AMD");
#endif // __BLOCK_BENCH_CHOLESKY_USE_AMD
#ifdef __BLOCK_BENCH_DUMP_MATRICES_IN_MATLAB_FORMAT
	printf("%s\n", "__BLOCK_BENCH_DUMP_MATRICES_IN_MATLAB_FORMAT");
#endif // __BLOCK_BENCH_DUMP_MATRICES_IN_MATLAB_FORMAT
#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
	printf("%s\n", "__BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
	printf("%s\n", "__BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

#if 0 // this is not so interesting at this point (tuning the solvers, the matrix is stable)
#ifdef __UBER_BLOCK_MATRIX_SUPRESS_FBS
	printf("%s\n", "__UBER_BLOCK_MATRIX_SUPRESS_FBS");
#endif // __UBER_BLOCK_MATRIX_SUPRESS_FBS
#ifdef __UBER_BLOCK_MATRIX_LAZY_PRODUCT
	printf("%s\n", "__UBER_BLOCK_MATRIX_LAZY_PRODUCT");
#endif // __UBER_BLOCK_MATRIX_LAZY_PRODUCT
#ifdef __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
	printf("%s\n", "__UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT");
#endif // __UBER_BLOCK_MATRIX_FBS_LAZY_PRODUCT
#ifdef __UBER_BLOCK_MATRIX_PERFCOUNTERS
	printf("%s\n", "__UBER_BLOCK_MATRIX_PERFCOUNTERS");
#endif // __UBER_BLOCK_MATRIX_PERFCOUNTERS
#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_MEMORY_DEBUGGING
	printf("%s\n", "__UBER_BLOCK_MATRIX_MULTIPLICATION_MEMORY_DEBUGGING");
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_MEMORY_DEBUGGING
#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
	printf("%s\n", "__UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR");
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_LINEAR
#ifdef __UBER_BLOCK_MATRIX_HYBRID_AT_A
	printf("%s\n", "__UBER_BLOCK_MATRIX_HYBRID_AT_A");
#endif // __UBER_BLOCK_MATRIX_HYBRID_AT_A
#ifdef __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
	printf("%s\n", "__UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS");
#endif // __UBER_BLOCK_MATRIX_MULTIPLICATION_PREALLOCATES_BLOCK_LISTS
#ifdef __UBER_BLOCK_MATRIX_ALIGN_BLOCK_MEMORY
	printf("%s\n", "__UBER_BLOCK_MATRIX_ALIGN_BLOCK_MEMORY");
#endif // __UBER_BLOCK_MATRIX_ALIGN_BLOCK_MEMORY
#endif // 0

#ifdef __SLAM_COUNT_ITERATIONS_AS_VERTICES
	printf("%s\n", "__SLAM_COUNT_ITERATIONS_AS_VERTICES");
#endif // __SLAM_COUNT_ITERATIONS_AS_VERTICES
#ifdef __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT
	printf("%s\n", "__MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT");
#endif // __MATRIX_ORDERING_TWO_LEVEL_CONSTRAINT
#ifdef __MATRIX_ORDERING_CACHE_ALIGN
	printf("%s\n", "__MATRIX_ORDERING_CACHE_ALIGN");
#endif // __MATRIX_ORDERING_CACHE_ALIGN
#ifdef __MATRIX_ORDERING_USE_AMD1
	printf("%s\n", "__MATRIX_ORDERING_USE_AMD1");
#endif // __MATRIX_ORDERING_USE_AMD1
#ifdef __MATRIX_ORDERING_USE_AMD_AAT
	printf("%s\n", "__MATRIX_ORDERING_USE_AMD_AAT");
#endif // __MATRIX_ORDERING_USE_AMD_AAT

#ifdef __SEGREGATED_MAKE_CHECKED_ITERATORS
	printf("%s\n", "__SEGREGATED_MAKE_CHECKED_ITERATORS");
#endif // __SEGREGATED_MAKE_CHECKED_ITERATORS
#ifdef __SEGREGATED_COMPILE_FAP_TESTS
	printf("%s\n", "__SEGREGATED_COMPILE_FAP_TESTS");
#endif // __SEGREGATED_COMPILE_FAP_TESTS

#ifdef __SE_TYPES_SUPPORT_A_SOLVERS
	printf("%s\n", "__SE_TYPES_SUPPORT_A_SOLVERS");
#endif // __SE_TYPES_SUPPORT_A_SOLVERS
#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
	printf("%s\n", "__SE_TYPES_SUPPORT_LAMBDA_SOLVERS");
#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
#ifdef __SE_TYPES_SUPPORT_L_SOLVERS
	printf("%s\n", "__SE_TYPES_SUPPORT_L_SOLVERS");
#endif // __SE_TYPES_SUPPORT_L_SOLVERS

#ifdef __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2
	printf("%s\n", "__NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2");
#endif // __NONLINEAR_SOLVER_LAMBDA_DUMP_CHI2

#ifdef __NONLINEAR_SOLVER_L_DUMP_TIMESTEPS
	printf("%s\n", "__NONLINEAR_SOLVER_L_DUMP_TIMESTEPS");
#endif // __NONLINEAR_SOLVER_L_DUMP_TIMESTEPS
#ifdef __NONLINEAR_SOLVER_L_DUMP_CHI2
	printf("%s\n", "__NONLINEAR_SOLVER_L_DUMP_CHI2");
#ifdef __NONLINEAR_SOLVER_L_DUMP_CHI2_AT_LAST_EDGE
	printf("%s\n", "__NONLINEAR_SOLVER_L_DUMP_CHI2_AT_LAST_EDGE");
#endif // __NONLINEAR_SOLVER_L_DUMP_CHI2_AT_LAST_EDGE
#endif // __NONLINEAR_SOLVER_L_DUMP_CHI2
#ifdef __NONLINEAR_SOLVER_L_DUMP_DENSITY
	printf("%s\n", "__NONLINEAR_SOLVER_L_DUMP_DENSITY");
#endif // __NONLINEAR_SOLVER_L_DUMP_DENSITY
#ifdef __NONLINEAR_SOLVER_L_DUMP_L_UPDATE_VARIANT_TIMES
	printf("%s\n", "__NONLINEAR_SOLVER_L_DUMP_L_UPDATE_VARIANT_TIMES");
#endif // __NONLINEAR_SOLVER_L_DUMP_L_UPDATE_VARIANT_TIMES
#ifdef __NONLINEAR_SOLVER_L_DETAILED_TIMING
	printf("%s\n", "__NONLINEAR_SOLVER_L_DETAILED_TIMING");
#endif // __NONLINEAR_SOLVER_L_DETAILED_TIMING
#ifdef __NONLINEAR_SOLVER_L_USE_SPARSE_BACKSUBST
	printf("%s\n", "__NONLINEAR_SOLVER_L_USE_SPARSE_BACKSUBST");
#endif // __NONLINEAR_SOLVER_L_USE_SPARSE_BACKSUBST
#ifdef __NONLINEAR_SOLVER_L_USE_RESUMED_BACKSUBST
	printf("%s\n", "__NONLINEAR_SOLVER_L_USE_RESUMED_BACKSUBST");
#endif // __NONLINEAR_SOLVER_L_USE_RESUMED_BACKSUBST
#ifdef __NONLINEAR_SOLVER_L_ENABLE_DENSE_CHOLESKY
	printf("%s\n", "__NONLINEAR_SOLVER_L_ENABLE_DENSE_CHOLESKY");
#endif // __NONLINEAR_SOLVER_L_ENABLE_DENSE_CHOLESKY

#ifdef __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1");
#endif // __NONLINEAR_SOLVER_FAST_L_BACKSUBSTITUTE_EACH_1
#ifdef __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE");
#endif // __NONLINEAR_SOLVER_FAST_L_ALWAYS_L_UPDATE
#ifdef __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING");
#endif // __NONLINEAR_SOLVER_FAST_L_VERIFY_PERM_FOLDING
#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS");
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_TIMESTEPS
#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_DUMP_CHI2");
#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE");
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2_AT_LAST_EDGE
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_CHI2
#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY");
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_DENSITY
#ifdef __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES");
#endif // __NONLINEAR_SOLVER_FAST_L_DUMP_L_UPDATE_VARIANT_TIMES
#ifdef __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING");
#endif // __NONLINEAR_SOLVER_FAST_L_DETAILED_TIMING
#ifdef __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
	printf("%s\n", "__NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY");
#endif // __NONLINEAR_SOLVER_FAST_L_ENABLE_DENSE_CHOLESKY
	printf("\n");
}

void PrintHelp()
{
	printf("General use:\n"
		"    ./SLAM_plus_plus -i <filename> --no-detailed-timing\n"
		"\n"
		"To run the pose-only datasets more quickly:\n"
		"    ./SLAM_plus_plus -i <filename> --pose-only --no-detailed-timing\n"
		"\n"
		"To run incrementally:\n"
		"    ./SLAM_plus_plus -lsp <optimize-each-N-steps> -i <filename> --no-detailed-timing\n"
		"\n"
		"This generates initial.txt and initial.tga, a description and image of the\n"
		"system before the final optimization, and solution.txt and solution.tga, a\n"
		"description and image of the final optimized system (unless --no-bitmaps\n"
		"is specified).\n"
		"\n"
		"--help|-h         displays this help screen\n"
		"--verbose|-v      displays verbose output while running (may slow down,\n"
		"                  especially in windows and if running incrementally)\n"
		"--silent|-s       suppresses displaying verbose output\n"
		"--no-show|-ns     doesn't show output image (windows only)\n"
		"--no-commandline|-nc    doesn't echo command line\n"
		"--no-flags|-nf    doesn't show compiler flags\n"
		"--no-detailed-timing    doesn't show detailed timing breakup (use this, you'll\n"
		"                  get confused)\n"
		"--no-bitmaps      doesn't write bitmaps initial.tga and solution.tga (neither\n"
		"                  the text files)\n"
		"--pose-only|-po   enables optimisation for pose-only slam (will warn and ignore\n"
		"                  on datasets with landmarks (only the first 1000 lines checked\n"
		"                  in case there are landmarks later, it would segfault))\n"
		"--use-old-code|-uogc    uses the old CSparse code (no block matrices in it)\n"
		"--a-slam|-A       uses A-SLAM\n"
		"--lambda|-,\\      uses lambda-SLAM (default, preferred batch solver)\n"
		"--l-slam|-L       uses L-SLAM\n"
		"--fast-l-slam|-fL uses the new fast L-SLAM solver (preferred incremental solver)\n"
		"--infile|-i <filename>    specifies input file <filename>; it can cope with\n"
		"                  many file types and conventions\n"
		"--parse-lines-limit|-pll <N>    sets limit of lines read from the input file\n"
		"                  (handy for testing), note this does not set limit of vertices\n"
		"                  nor edges!\n"
		"--linear-solve-period|-lsp <N>    sets period for incrementally running linear\n"
		"                  solver (default 0: disabled)\n"
		"--nonlinear-solve-period|-nsp <N>    sets period for incrementally running\n"
		"                  non-linear solver (default 0: disabled)\n"
		"--max-nonlinear-solve-iters|-mnsi <N>    sets maximal number of nonlinear\n"
		"                  solver iterations (default 10)\n"
		"--nonlinear-solve-error-thresh|-nset <f>    sets nonlinear solve error threshold\n"
		"                  (default 20)\n"
		"--max-final-nonlinear-solve-iters|-mfnsi <N>    sets maximal number of final\n"
		"                  optimization iterations (default 5)\n"
		"--final-nonlinear-solve-error-thresh|-fnset <f>    sets final nonlinear solve\n"
		"                  error threshold (default .01)\n"
		"--run-matrix-benchmarks|-rmb <benchmark-name> <benchmark-type>    runs block\n"
		"                  matrix benchmarks (benchmark-name is name of a folder with\n"
		"                  UFLSMC benchmark, benchmark-type is one of alloc, factor, all)\n"
		"--run-matrix-unit-tests|-rmut    runs block matrix unit tests\n");
}

/**
 *	@def p_s_Number_Suffix
 *	@brief gets a suffix for a number
 *	@param[in] n_x is a number for which the suffix is to be returned, must be positive
 *	@return Returns suffix for the given number.
 */
#define p_s_Number_Suffix(n_x) (((n_x + 9) % 10 > 2 || (n_x > 10 && n_x < 20))? "th" : \
	((n_x) % 10 == 1)? "st" : ((n_x) % 10 == 2)? "nd" : "rd")

/**
 *	@brief main
 *
 *	@param[in] n_arg_num is number of commandline arguments
 *	@param[in] p_arg_list is the list of commandline arguments
 *
 *	@return Returns 0 on success, -1 on failure.
 */
int main(int n_arg_num, const char **p_arg_list)
{
	//tetris_main(n_arg_num, p_arg_list);

#if (defined(_WIN32) || defined(_WIN64)) && defined(__PIMP_MY_WINDOW)
	{
		system("title SLAM ++\n");
		HWND h_console = FindWindow(0, "SLAM ++");
		int n_all_mons_width = GetSystemMetrics(SM_CXVIRTUALSCREEN);
		int n_all_mons_height = GetSystemMetrics(SM_CYVIRTUALSCREEN);
		int n_pri_mon_width = GetSystemMetrics(SM_CXSCREEN);
		int n_pri_mon_height = GetSystemMetrics(SM_CYSCREEN);
		if(n_all_mons_width > n_pri_mon_width + n_pri_mon_width / 2) {
			RECT t_rect;
			GetWindowRect(h_console, &t_rect);
			if(t_rect.left < n_pri_mon_width) // running batch job? don't float away
				SetWindowPos(h_console, HWND_TOP, t_rect.left + n_pri_mon_width, t_rect.top, 0, 0, SWP_NOSIZE);
		} else if(n_all_mons_height > n_pri_mon_height + n_pri_mon_height / 2) {
			RECT t_rect;
			GetWindowRect(h_console, &t_rect);
			if(t_rect.top < n_pri_mon_height) // running batch job? don't float away
				SetWindowPos(h_console, HWND_TOP, t_rect.left, t_rect.top + n_pri_mon_height, 0, 0, SWP_NOSIZE);
		}
	}
	// windows hack - make the console window appear on the secondary monitor
#endif // (_WIN32 || _WIN64) && __PIMP_MY_WINDOW

	TCommandLineArgs t_cmd_args;
	t_cmd_args.Defaults(); // set defaults
	if(!t_cmd_args.Parse(n_arg_num, p_arg_list))
		return -1;
	// parse commandline

	if(t_cmd_args.b_run_matrix_unit_tests) {
#ifdef __UBER_BLOCK_MATRIX_UNIT_TESTS_INCLUDED
		CBlockMatrixUnitTests::RunAll();
		return 0;
#else // __UBER_BLOCK_MATRIX_UNIT_TESTS_INCLUDED
		fprintf(stderr, "warning: unit tests not included: skipping\n");
		return -1;
#endif // __UBER_BLOCK_MATRIX_UNIT_TESTS_INCLUDED
	}
	// run unit tests

	if(t_cmd_args.b_run_matrix_benchmarks)
		return n_Run_BlockBenchmark(n_dummy_param, t_cmd_args.p_s_bench_name, t_cmd_args.p_s_bench_type);
	// run benchmarks (allow without input file)

	if(!t_cmd_args.p_s_input_file) {
		fprintf(stderr, "error: no input file specified; find help below:\n");
		PrintHelp();
		return -1;
	}
	if(t_cmd_args.n_linear_solve_each_n_steps && t_cmd_args.n_nonlinear_solve_each_n_steps &&
	   t_cmd_args.n_nonlinear_solve_each_n_steps < t_cmd_args.n_linear_solve_each_n_steps) {
		fprintf(stderr, "warning: nonlinear solver is scheduled to run every %d%s step, "
			"linear solver runs every %d%s step. the nonlinear solver will never run.\n",
			int(t_cmd_args.n_linear_solve_each_n_steps),
			p_s_Number_Suffix(t_cmd_args.n_linear_solve_each_n_steps),
			int(t_cmd_args.n_nonlinear_solve_each_n_steps),
			p_s_Number_Suffix(t_cmd_args.n_nonlinear_solve_each_n_steps));
	}
	// check inputs

	if(t_cmd_args.b_show_commandline) {
		printf("> ./SLAM_plus_plus");
		for(int i = 1; i < n_arg_num; ++ i)
			printf(" %s", p_arg_list[i]);
		printf("\n");
	}
	// display commandline

	if(t_cmd_args.b_show_flags)
		DisplaySwitches();
	// display switches

	TDatasetPeeker peek(t_cmd_args.p_s_input_file, t_cmd_args.n_max_lines_to_process);
	if(t_cmd_args.b_10k_opts && peek.b_has_landmark) {
		fprintf(stderr, "error: --pose-only flag detected, but the system has landmarks (flag cleared)\n");
		t_cmd_args.b_10k_opts = false;
	} else if(!t_cmd_args.b_10k_opts && !peek.b_has_landmark) {
		fprintf(stderr, "warning: the system doesn't seem to have landmarks"
			" and --pose-only flag not present: could run faster\n");
	}
	if(peek.b_has_edge3d || peek.b_has_vertex3d) {
		t_cmd_args.b_use_SE3 = true;
		if(t_cmd_args.b_verbose)
			fprintf(stderr, "detected SE3\n");
	}
	if(peek.b_has_ba) {
		t_cmd_args.b_use_BA = true;
		if(t_cmd_args.b_verbose)
			fprintf(stderr, "detected BA\n");
	}
	// detect landmarks, clear b_10k_opts where it would cause the code to crash

	//Eigen::initParallel();
	//Eigen::setNbThreads(0);
	//printf("Eigen is running in %d threads\n", Eigen::nbThreads(n));
	// ...

	/*Test_FBS();
	return 0;*/

	if(t_cmd_args.b_use_old_system) {
#if 0
		CTimer t;
		t.ResetTimer();

		COptimizationSystem<TVertex2D, TEdge2D> system;
		// system to be optimized (2D case)

		bool b_no_blocks = true; // todo - make this commandline
		if(!b_no_blocks && b_L_slam)
			fprintf(stderr, "warning: L-SLAM requested but using the experimental blocky solver\n");
		else if(!b_no_blocks && !b_L_slam)
			fprintf(stderr, "warning: A-SLAM requested but using the experimental blocky solver\n");
		else if(b_no_blocks)
			fprintf(stderr, "warning: *not* using the experimental blocky solver\n");
		COptimizationMethod_Interface *p_solver = (b_no_blocks)? ((b_L_slam)?
			(COptimizationMethod_Interface*) new(std::nothrow) CIncrementalCholeskyLSLAMSolver(system, n_linear_solve_each_n_steps) :
			(COptimizationMethod_Interface*) new(std::nothrow) CIncrementalCholeskyASLAMSolver(system,
			n_linear_solve_each_n_steps, n_nonlinear_solve_each_n_steps,
			n_max_nonlinear_solve_iteration_num, f_nonlinear_solve_error_threshold, b_verbose)) :
			(COptimizationMethod_Interface*) new(std::nothrow)  CIncrementalCholeskyASLAMSolver_Blocky(system,
			n_linear_solve_each_n_steps, n_nonlinear_solve_each_n_steps,
			n_max_nonlinear_solve_iteration_num, f_nonlinear_solve_error_threshold, b_verbose);
		if(!p_solver) {
			fprintf(stderr, "error: not enough memory\n");
			return -1;
		}
		// the solver to optimize the system

		CParser p;
		CParserBase::CParserAdaptor *p_sink = p_solver; // solvers must implement parser adaptor
		if(!p.Parse(p_s_input_file, p_sink, n_max_lines_to_process)) {
			fprintf(stderr, "error: failed to parse input file\n");
			delete p_solver;
			return -1;
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
		}
		double f_time_initial_save_end = t.f_Time();

		p_solver->Optimize(n_max_final_optimization_iteration_num, f_final_optimization_threshold);
		// perform the final optimization

		double f_time = t.f_Time() - (f_time_initial_save_end - f_time_initial_save_start); // don't count saving of the initial system state (rasterizing the image takes some time)
		printf("\ndone. it took " PRItimeprecise " (%lf sec)\n", PRItimeparams(f_time), f_time);
		//printf("time / 1.6 = " PRItimeprecise "\n", PRItimeparams(f_time / 1.6));
		// display time it took

		if(!b_no_blocks)
			((CIncrementalCholeskyASLAMSolver_Blocky*)p_solver)->Dump();

		if(b_write_bitmaps) {
			system.Plot2D("solution.tga");
			system.Dump("solution.txt");
		}
		// save the solution to a file

		delete p_solver;
#else // 0
		fprintf(stderr, "error: the legacy solver was not compiled\n");
		return -1;
#endif // 0
	} else {
		try {
			if(t_cmd_args.b_use_BA) {
				if(n_Run_BA_Solver(t_cmd_args))
					return -1;
			} else if(t_cmd_args.b_use_SE3) {
				if(n_Run_SE3_Solver(t_cmd_args))
					return -1;
			} else if(t_cmd_args.b_10k_opts) {
				if(n_Run_SE2PoseOnly_Solver(t_cmd_args))
					return -1;
			} else {
				if(n_Run_SE2_Solver(t_cmd_args))
					return -1;
			}
			// a simple selection of solver and problem type, without all the ifdef mess
#if 0
			if(t_cmd_args.n_solver_choice == nlsolver_SPCG) {
#ifdef __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_SPCG
#ifdef __NONLINEAR_BLOCKY_SOLVER_SPCG_INCLUDED
				if(b_use_SE3) {
					fprintf(stderr, "error: SE(3) solver not instanced for SPCG\n");
					return -1;
				} else if(b_10k_opts) {
					typedef MakeTypelist1(CVertexPose2D) TVertexTypelist;
					typedef MakeTypelist1(CEdgePose2D) TEdgeTypelist;
					typedef CFlatSystem<CVertexPose2D, TVertexTypelist,
						/*CVertexTypeTraits,*/ CEdgePose2D, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) pose vertex and edge types (only poses, not landmarks)

					//typedef MakeTypelist_Safe((Eigen::Matrix3d)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_SPCG, CSE2OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else {
					typedef MakeTypelist2(CVertexPose2D, CVertexLandmark2D) TVertexTypelist;
					typedef MakeTypelist2(CEdgePose2D, CEdgePoseLandmark2D) TEdgeTypelist;
					typedef CFlatSystem<CSEBaseVertex, TVertexTypelist,
						/*CVertexTypeTraits,*/ CSEBaseEdge, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) vertex and edge types

					//typedef MakeTypelist_Safe((Eigen::Matrix3d, Eigen::Matrix2d, Eigen::Matrix<double, 2, 3>)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_SPCG, CSE2EdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				}
#else // __NONLINEAR_BLOCKY_SOLVER_SPCG_INCLUDED
				fprintf(stderr, "error: CNonlinearSolver_SPCG was not included\n"); // not very precise, but tells the story well
				return -1;
#endif // __NONLINEAR_BLOCKY_SOLVER_SPCG_INCLUDED
#else // __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_SPCG
				fprintf(stderr, "error: CNonlinearSolver_SPCG was not compiled\n"); // not very precise, but tells the story well
				return -1;
#endif // __SE_TYPES_SUPPORT_NONLINEAR_SOLVER_SPCG
			} else if(t_cmd_args.n_solver_choice == nlsolver_Lambda) {
#ifdef __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
				if(b_use_BA) {	//BA
					typedef MakeTypelist_Safe((CVertexCam, CVertexXYZ)) TVertexTypelist_SE3; // just put your vertex types in the list (note the double parentheses are required)
					typedef MakeTypelist_Safe((CEdgeP2C3D)) TEdgeTypelist_SE3; // just put your edge types in the list (note the double parentheses are required)
					typedef CFlatSystem<CSEBaseVertex, TVertexTypelist_SE3,
						CSEBaseEdge, TEdgeTypelist_SE3> CSystemType; // the Base types can be your types in case you are not using CBaseSEVertexImpl / CBaseSEEdgeImpl

					//typedef MakeTypelist_Safe((Eigen::Matrix3d)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_Lambda_LM, CBATraits, CBAParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				}
#else // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
				if(b_use_BA) {
					fprintf(stderr, "error: CNonlinearSolver_Lambda_LM was not included\n");
					return -1;
				}
#endif // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_LM_INCLUDED
#ifdef __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
				if(b_use_SE3) { // SE3
					typedef MakeTypelist_Safe((CVertexPose3D)) TVertexTypelist_SE3; // just put your vertex types in the list (note the double parentheses are required)
					typedef MakeTypelist_Safe((CEdgePose3D)) TEdgeTypelist_SE3; // just put your edge types in the list (note the double parentheses are required)
					typedef CFlatSystem<CSEBaseVertex, TVertexTypelist_SE3,
					 CSEBaseEdge, TEdgeTypelist_SE3> CSystemType; // the Base types can be your types in case you are not using CBaseSEVertexImpl / CBaseSEEdgeImpl

					//typedef MakeTypelist_Safe((Eigen::Matrix3d)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_Lambda, CSE3OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else if(t_cmd_args.b_10k_opts) {
					typedef MakeTypelist1(CVertexPose2D) TVertexTypelist;
					typedef MakeTypelist1(CEdgePose2D) TEdgeTypelist;
					typedef CFlatSystem<CVertexPose2D, TVertexTypelist,
						/*CVertexTypeTraits,*/ CEdgePose2D, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) pose vertex and edge types (only poses, not landmarks)

					//typedef MakeTypelist_Safe((Eigen::Matrix3d)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_Lambda, CSE2OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else {
					typedef MakeTypelist2(CVertexPose2D, CVertexLandmark2D) TVertexTypelist;
					typedef MakeTypelist2(CEdgePose2D, CEdgePoseLandmark2D) TEdgeTypelist;
					typedef CFlatSystem<CSEBaseVertex, TVertexTypelist,
						/*CVertexTypeTraits,*/ CSEBaseEdge, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) vertex and edge types

					//typedef MakeTypelist_Safe((Eigen::Matrix3d, Eigen::Matrix2d, Eigen::Matrix<double, 2, 3>)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_Lambda, CSE2EdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				}
#else // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
				if(!b_use_BA) {
					fprintf(stderr, "error: CNonlinearSolver_Lambda was not included\n"); // not very precise, but tells the story well
					return -1;
				}
#endif // __NONLINEAR_BLOCKY_SOLVER_LAMBDA_INCLUDED
#else // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
				fprintf(stderr, "error: SE types do not support lambda solvers\n");
				return -1;
#endif // __SE_TYPES_SUPPORT_LAMBDA_SOLVERS
			} else if(t_cmd_args.n_solver_choice == nlsolver_L) {
#ifdef __SE_TYPES_SUPPORT_L_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_L_INCLUDED
				if(b_use_SE3) { // SE3
					typedef MakeTypelist_Safe((CVertexPose3D)) TVertexTypelist_SE3; // just put your vertex types in the list (note the double parentheses are required)
					typedef MakeTypelist_Safe((CEdgePose3D)) TEdgeTypelist_SE3; // just put your edge types in the list (note the double parentheses are required)
					typedef CFlatSystem<CVertexPose3D, TVertexTypelist_SE3,
							CEdgePose3D, TEdgeTypelist_SE3> CSystemType; // the Base types can be your types in case you are not using CBaseSEVertexImpl / CBaseSEEdgeImpl

					if(!CTester<CSystemType, CNonlinearSolver_L, CSE3OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else if(b_10k_opts) {
					typedef MakeTypelist1(CVertexPose2D) TVertexTypelist;
					typedef MakeTypelist1(CEdgePose2D) TEdgeTypelist;
					typedef CFlatSystem<CVertexPose2D, TVertexTypelist,
						/*CVertexTypeTraits,*/ CEdgePose2D, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) pose vertex and edge types (only poses, not landmarks)

					//typedef MakeTypelist_Safe((Eigen::Matrix3d)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_L, CSE2OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else {
					typedef MakeTypelist2(CVertexPose2D, CVertexLandmark2D) TVertexTypelist;
					typedef MakeTypelist2(CEdgePose2D, CEdgePoseLandmark2D) TEdgeTypelist;
					typedef CFlatSystem<CSEBaseVertex, TVertexTypelist,
						/*CVertexTypeTraits,*/ CSEBaseEdge, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) vertex and edge types

					//typedef MakeTypelist_Safe((Eigen::Matrix3d, Eigen::Matrix2d, Eigen::Matrix<double, 2, 3>)) TSE2MatrixTypes;
					//typedef MakeTypelist_Safe((Eigen::Matrix3d, Eigen::Matrix2d,
					//	Eigen::Matrix<double, 2, 3>, Eigen::Matrix<double, 3, 2>)) TSE2MatrixTypes; // use combinations as well (lazy me, but it is needed for most of FBS operations inside anyway) // t_odo - write the bloody template to calculate the carthesian product for me
					if(!CTester<CSystemType, CNonlinearSolver_L, CSE2EdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				}
#else // __NONLINEAR_BLOCKY_SOLVER_L_INCLUDED
				fprintf(stderr, "error: CNonlinearSolver_L was not included\n"); // not very precise, but tells the story well
				return -1;
#endif // __NONLINEAR_BLOCKY_SOLVER_L_INCLUDED
#else // __SE_TYPES_SUPPORT_L_SOLVERS
				fprintf(stderr, "error: CNonlinearSolver_L was not compiled\n"); // not very precise, but tells the story well
				return -1;
#endif // __SE_TYPES_SUPPORT_L_SOLVERS
			} else if(t_cmd_args.n_solver_choice == nlsolver_FastL) {
#ifdef __SE_TYPES_SUPPORT_L_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED
				if(b_use_SE3) {
					typedef MakeTypelist_Safe((CVertexPose3D)) TVertexTypelist_SE3; // just put your vertex types in the list (note the double parentheses are required)
					typedef MakeTypelist_Safe((CEdgePose3D)) TEdgeTypelist_SE3; // just put your edge types in the list (note the double parentheses are required)
					typedef CFlatSystem<CVertexPose3D, TVertexTypelist_SE3,
							CEdgePose3D, TEdgeTypelist_SE3> CSystemType; // the Base types can be your types in case you are not using CBaseSEVertexImpl / CBaseSEEdgeImpl

					if(!CTester<CSystemType, CNonlinearSolver_FastL, CSE3OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else if(b_10k_opts) {
					typedef MakeTypelist1(CVertexPose2D) TVertexTypelist;
					typedef MakeTypelist1(CEdgePose2D) TEdgeTypelist;
					typedef CFlatSystem<CVertexPose2D, TVertexTypelist,
						/*CVertexTypeTraits,*/ CEdgePose2D, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) pose vertex and edge types (only poses, not landmarks)

					//typedef MakeTypelist_Safe((Eigen::Matrix3d)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_FastL, CSE2OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else {
					typedef MakeTypelist2(CVertexPose2D, CVertexLandmark2D) TVertexTypelist;
					typedef MakeTypelist2(CEdgePose2D, CEdgePoseLandmark2D) TEdgeTypelist;
					typedef CFlatSystem<CSEBaseVertex, TVertexTypelist,
						/*CVertexTypeTraits,*/ CSEBaseEdge, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) vertex and edge types

					//typedef MakeTypelist_Safe((Eigen::Matrix3d, Eigen::Matrix2d, Eigen::Matrix<double, 2, 3>)) TSE2MatrixTypes;
					//typedef MakeTypelist_Safe((Eigen::Matrix3d, Eigen::Matrix2d,
					//	Eigen::Matrix<double, 2, 3>, Eigen::Matrix<double, 3, 2>)) TSE2MatrixTypes; // use combinations as well (lazy me, but it is needed for most of FBS operations inside anyway) // t_odo - write the bloody template to calculate the carthesian product for me
					if(!CTester<CSystemType, CNonlinearSolver_FastL, CSE2EdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				}
#else // __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED
				fprintf(stderr, "error: CNonlinearSolver_FastL was not included\n"); // not very precise, but tells the story well
				return -1;
#endif // __NONLINEAR_BLOCKY_SOLVER_FAST_L_INCLUDED
#else // __SE_TYPES_SUPPORT_L_SOLVERS
				fprintf(stderr, "error: CNonlinearSolver_FastL was not compiled\n"); // not very precise, but tells the story well
				return -1;
#endif // __SE_TYPES_SUPPORT_L_SOLVERS
			} else {
				_ASSERTE(t_cmd_args.n_solver_choice == nlsolver_A);
#ifdef __SE_TYPES_SUPPORT_A_SOLVERS
#ifdef __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED
				if(b_use_BA) {	//BA
					typedef MakeTypelist_Safe((CVertexCam, CVertexXYZ)) TVertexTypelist_SE3; // just put your vertex types in the list (note the double parentheses are required)
					typedef MakeTypelist_Safe((CEdgeP2C3D)) TEdgeTypelist_SE3; // just put your edge types in the list (note the double parentheses are required)
					typedef CFlatSystem<CSEBaseVertex, TVertexTypelist_SE3,
						CSEBaseEdge, TEdgeTypelist_SE3> CSystemType; // the Base types can be your types in case you are not using CBaseSEVertexImpl / CBaseSEEdgeImpl

					//typedef MakeTypelist_Safe((Eigen::Matrix3d)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_A, CBATraits, CBAParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else if(b_use_SE3) { // SE3
					typedef MakeTypelist_Safe((CVertexPose3D)) TVertexTypelist_SE3; // just put your vertex types in the list (note the double parentheses are required)
					typedef MakeTypelist_Safe((CEdgePose3D)) TEdgeTypelist_SE3; // just put your edge types in the list (note the double parentheses are required)
					typedef CFlatSystem<CVertexPose3D, TVertexTypelist_SE3,
							CEdgePose3D, TEdgeTypelist_SE3> CSystemType; // the Base types can be your types in case you are not using CBaseSEVertexImpl / CBaseSEEdgeImpl

					if(!CTester<CSystemType, CNonlinearSolver_A, CSE3OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else if(b_10k_opts) {
					typedef MakeTypelist1(CVertexPose2D) TVertexTypelist;
					typedef MakeTypelist1(CEdgePose2D) TEdgeTypelist;
					typedef CFlatSystem<CVertexPose2D, TVertexTypelist,
						/*CVertexTypeTraits,*/ CEdgePose2D, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) pose vertex and edge types (only poses, not landmarks)

					//typedef MakeTypelist_Safe((Eigen::Matrix3d)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_A, CSE2OnlyPoseEdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				} else {
					typedef MakeTypelist2(CVertexPose2D, CVertexLandmark2D) TVertexTypelist;
					typedef MakeTypelist2(CEdgePose2D, CEdgePoseLandmark2D) TEdgeTypelist;
					typedef CFlatSystem<CSEBaseVertex, TVertexTypelist,
						/*CVertexTypeTraits,*/ CSEBaseEdge, TEdgeTypelist> CSystemType;
					// make a system permitting SE(2) vertex and edge types

					//typedef MakeTypelist_Safe((Eigen::Matrix3d, Eigen::Matrix2d, Eigen::Matrix<double, 2, 3>)) TSE2MatrixTypes;
					if(!CTester<CSystemType, CNonlinearSolver_A, CSE2EdgeTraits, CParseLoop>::Run_and_Shout(
					   p_s_input_file, n_max_lines_to_process, n_linear_solve_each_n_steps,
					   n_nonlinear_solve_each_n_steps, n_max_nonlinear_solve_iteration_num,
					   f_nonlinear_solve_error_threshold, n_max_final_optimization_iteration_num,
					   f_final_optimization_threshold, b_verbose, b_show_detailed_timing,
					   b_write_bitmaps))
						return -1;
					// run with parameters
				}
#else // __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED
				fprintf(stderr, "error: CNonlinearSolver_A was not included\n"); // not very precise, but tells the story well
				return -1;
#endif // __NONLINEAR_BLOCKY_SOLVER_A_INCLUDED
#else // __SE_TYPES_SUPPORT_A_SOLVERS
				fprintf(stderr, "error: CNonlinearSolver_A was not compiled\n"); // not very precise, but tells the story well
				return -1;
#endif // __SE_TYPES_SUPPORT_A_SOLVERS
			}
#endif // 0
		} catch(std::runtime_error &r_exc) {
			fprintf(stderr, "error: uncaught exception: \'%s\'\n", r_exc.what());
			return -1;
		}
	}

	if(!t_cmd_args.b_no_show) { // to get rid of the "set but not used" warning
#if defined(_WIN32) || defined(_WIN64)
		/*char p_s_class_name[256];
		GetClassName(FindWindow(0, "XnView - [solution.tga]"), p_s_class_name, 255);
		printf("classname is: \'%s\'\n", p_s_class_name);*/
		// this is how you get classname of your image browser window

		if(!FindWindow("XmainClass", 0))
			ShellExecute(0, "open", "solution.tga", 0, 0, SW_SHOW);
		else
			fprintf(stderr, "warning: XnView already running: new window not opened\n");
		// detect if image browser is running
#endif // _WIN32 || _WIN64
	}
	// maybe a bad idea ...

	return 0;
}

#ifdef __COMPILE_LIBRARY_TESTS

#include <iostream>

/**
 *	@brief makes sure that CSparse is linked correctly (and is working)
 */
static void CSparse_Test()
{
	cs *test_my_matrix = cs_spalloc(4, 4, 16, 1, 1);
	for(int i = 0; i < 4; ++ i)
		cs_entry(test_my_matrix, i, i, 1);
	cs *test_my_matrix2 = cs_compress(test_my_matrix);
	cs_print(test_my_matrix, 0);
	cs_print(test_my_matrix2, 0);
	cs_spfree(test_my_matrix);
	cs_spfree(test_my_matrix2);
	// just make sure that CSparse is working
}

/**
 *	@brief makes sure that eigen is linked correctly (and is working)
 */
static void Eigen3_Test()
{
	Eigen::Matrix<double, 3, 3> C;
	const double alpha = 30 / 180.0 * M_PI;

	double cosa = cos(alpha),
		   sina = sin(alpha);
	C << cosa, -sina,   0.0,
		 sina,  cosa,   0.0,
		 0.0,   0.0,   1.0;
	// Create the principal rotation matrix C_3(alpha)

	Eigen::Matrix<double, 3, 1> a(1, 2, 3);
	Eigen::Matrix<double, 3, 1> b = C * a;

	std::cout << "C = " << C << std::endl;
	std::cout << "a = " << a << std::endl;
	std::cout << "C * a = " << b << std::endl;
}

#if 0

/**
 *	@brief tests reference matrix concept (not a library test, actually)
 */
static void Test_RefMatrix() // t_odo - see assembly to see what's going on // don't have ReferencingMatrixXd anymore, using Eigen::Map<> now.
{
	double p_data[16] = {
		0, 1, 2, 3,
		4, 5, 6, 7,
		8, 9, 10, 11,
		12, 13, 14, 15
	};

	Eigen::ReferencingMatrixXd matrix(4, 4, p_data); // no copy of the data is made (at least i hope :))

	Eigen::Matrix4d b, c; // ordinary matrices

	b << 1, 0, 0, 0,
		 0, 1, 0, 0,
		 0, 0, 1, 0,
		 0, 0, 0, 1;

	c = matrix * b; // see what happens

	b = matrix * matrix;
	// seems legit - result of multiplying two reference matrices is matrixXd // t_odo - check it

	printf("b = {\n");
	for(int i = 0; i < 4; ++ i) {
		printf("    ");
		for(int j = 0; j < 4; ++ j)
			printf("%3g ", b(i, j));
		printf("\n");
	}
	printf("}\n");
	printf("c = {\n");
	for(int i = 0; i < 4; ++ i) {
		printf("    ");
		for(int j = 0; j < 4; ++ j)
			printf("%2g ", c(i, j));
		printf("\n");
	}
	printf("}\n");
}

#endif // 0

/**
 *	@brief tests forward allocated pool (not a library test, actually)
 */
static void Test_FAP()
{
	{
		forward_allocated_pool<int, 0> fap100(100);
		forward_allocated_pool<int, 0> fap200(200);
		fap100.swap(fap200);
	}
	{
		forward_allocated_pool<int, 0> fap100(100);
		forward_allocated_pool<int, 0> fap200(100);
		fap100.swap(fap200);
	}
	{
		forward_allocated_pool<int, 100> fap100;
		forward_allocated_pool<int, 0> fap200(100);
		fap100.swap(fap200);
	}
	{
		forward_allocated_pool<int, 100> fap100;
		forward_allocated_pool<int, 200> fap200;
		fap100.swap(fap200);
	}
	{
		forward_allocated_pool<int, 100> fap100;
		forward_allocated_pool<int, 100> fap200;
		fap100.swap(fap200);
	}

	{
		forward_allocated_pool<int, 100> fap100;
		_ASSERTE(fap100.page_size() == 100);
		forward_allocated_pool<int, 0> fap(100);
		_ASSERTE(fap.page_size() == 100);
		forward_allocated_pool<int> fap2(100);
		_ASSERTE(fap2.page_size() == 100);
	}

	if(!CFAP_Test<100, false>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");
	if(!CFAP_Test<100, true>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");

	if(!CFAP_Test<10, false>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");
	if(!CFAP_Test<10, true>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");

	if(!CFAP_Test<13, false>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");
	if(!CFAP_Test<13, true>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");

	if(!CFAP_Test<128, false>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");
	if(!CFAP_Test<128, true>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");

	if(!CFAP_Test<1, false>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");
	if(!CFAP_Test<1, true>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");

	if(!CFAP_Test<2, false>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");
	if(!CFAP_Test<2, true>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");

	if(!CFAP_Test<4, false>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");
	if(!CFAP_Test<4, true>::FAP_Test(1000))
		fprintf(stderr, "error: tests failed\n");
}

/**
 *	@brief tests block matrix allocator using tetris (not a library test, actually)
 */
static void Test_Tetris(int n_arg_num, const char **p_arg_list)
{
	tetris_main(n_arg_num, p_arg_list);
}

#endif // __COMPILE_LIBRARY_TESTS

/*
 *	end-of-file
 */
