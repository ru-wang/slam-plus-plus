/*
								+-----------------------------------+
								|                                   |
								|   ***  Schur linear solver  ***   |
								|                                   |
								|  Copyright  (c) -tHE SWINe- 2015  |
								|                                   |
								|    LinearSolver_Schur_GPU.cpp     |
								|                                   |
								+-----------------------------------+
*/

/**
 *	@file src/slam/LinearSolver_Schur_GPU.cpp
 *	@brief GPU routines for Schur complement solver
 *	@author -tHE SWINe-
 *	@date 2015-03-09
 */

#include "slam/LinearSolver_Schur.h"
#include "slam/MemUsage.h" // PRIsizeB

//#define __SCHUR_PROFILING // defined in slam/LinearSolver_Schur.h
// want to see profile info?

#if !defined(__DISABLE_GPU)

#include <cula.hpp>
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <signal.h>

#if (defined(_WIN32) || defined(_WIN64)) && defined(__GPU_EXPLICIT_LINKAGE)

//#pragma comment(lib, "cula.lib") // r3.1
//#pragma comment(lib, "cula_lapack.lib") // r17
#pragma comment(lib, "cula_core.lib") // r16a
#pragma comment(lib, "cula_lapack.lib") // r16a
// link with CULA (windows only)

#pragma comment(lib, "cuda.lib") // cuda
#pragma comment(lib, "cublas.lib") // cublas
#pragma comment(lib, "cusparse.lib") // cusparse

#endif // (_WIN32 || _WIN64) && __GPU_EXPLICIT_LINKAGE

/**
 *	@def __GPU_SCHUR_VERIFY_RESULT
 *	@brief calculates Schur on both GPU and CPU and compares them (makes it slower)
 */
//#define __GPU_SCHUR_VERIFY_RESULT

/**
 *	@def __GPU_SCHUR_EASY_PROD_ONLY
 *	@brief this only runs U*(C^-1) on the GPU, leaving U*(C^-1)*V to the CPU
 *
 *	Both branches are functional, but the full product on GPU will likely
 *	be slower than the SSE optimized implementation unless running
 *	on a high-end GPU.
 */
#define __GPU_SCHUR_EASY_PROD_ONLY

/**
 *	@def MAKE_FILE_LINE_DEBUG_STR1
 *
 *	@brief concatenates filename and line to a single "C" string
 *
 *	@param[in] a is the opening bracket
 *	@param[in] b is file name
 *	@param[in] c is comma separator
 *	@param[in] d is line number
 *	@param[in] e is the closing bracket
 */
#define MAKE_FILE_LINE_DEBUG_STR1(a,b,c,d,e) a b c #d e

/**
 *	@def MAKE_FILE_LINE_DEBUG_STR0
 *
 *	@brief concatenates filename and line to a single "C" string
 *
 *	@param[in] a is the opening bracket
 *	@param[in] b is file name
 *	@param[in] c is comma separator
 *	@param[in] d is line number
 *	@param[in] e is the closing bracket
 */
#define MAKE_FILE_LINE_DEBUG_STR0(a,b,c,d,e) MAKE_FILE_LINE_DEBUG_STR1(a, b, c, d, e)

/**
 *	@def MAKE_FILE_LINE_DEBUG_STR
 *
 *	@brief concatenates filename and line to a single "C" string
 *
 *	@param[in] f is file name
 *	@param[in] l is line number
 */
#define MAKE_FILE_LINE_DEBUG_STR(f,l) MAKE_FILE_LINE_DEBUG_STR0("(", f, ", line ", l, ")")

/**
 *	@def FILE_LINE
 *	@brief expands to a string, containing parenthesized file and line as a string literal
 */
#define FILE_LINE MAKE_FILE_LINE_DEBUG_STR(__FILE__, __LINE__)

/*
 *								=== CGPUGuard ===
 */

/**
 *	@brief a simple helper that helps to avoid crashes caused by interrupting the program while calling GPU kernels
 *
 *	At the beginning, call CGPUGuard::Register_SignalHandlers() to register signal handlers.
 *	In your program, place __GPU_Function; at the first line of each function / at the beginning
 *	of each block that is using GPU. That will mark the regions where the program should not
 *	be terminated. If e.g. SIGINT is caught, the program quits as soon as it exits all the regions
 *	that are using GPU.
 *
 *	@note This is OpenMP aware and should work even if multiple threads are using the GPU.
 *	@note The __GPU_Function presents a small overhead, and should be placed reasonably. If, on the
 *		other hand, there is a single function looping on data and using GPU, it might be a good
 *		idea to make the scope the loop iteration instead of the whole function. It is possible
 *		to use __GPU_ScopeCheck in short functions to avoid __GPU_Function overhead but to still
 *		make sure that the caller used __GPU_Function.
 */
class CGPUGuard {
protected:
	static bool b_quit; /**< @brief quit flag; set once interrupted in a GPU scope (when not possible to quit right away) */
	static size_t n_using_GPU; /**< @brief number of function recursions / blocks that are marked as using GPU */

public:
	/**
	 *	@brief default constructor; marks the beginning of a GPU scope
	 */
	CGPUGuard()
	{
		if(b_quit)
			SafeExit();
		Set_UsingGPU(true);
	}

	/**
	 *	@brief destructor; marks the end of a GPU scope
	 */
	~CGPUGuard()
	{
		Set_UsingGPU(false);
		if(b_quit)
			SafeExit();
	}

	/**
	 *	@brief GPU-aware safe exit function; quits when the GPU is not in use
	 *	@note If the GPU is in use right now, it schedules exit at the end of the current GPU scope.
	 */
	static void SafeExit()
	{
		Set_UsingGPU(false, true); // the first arg is ignored in this case
		b_quit = true; // were not able to quit now? nevermind, we'll get there
	}

	/**
	 *	@brief gets GPU scope flag
	 *	@return Returns true if the GPU is currently in use, otherwise returns false.
	 *	@note This is not thread-safe: the value might have changed since this function returned.
	 */
	static bool b_Using_GPU()
	{
		return n_using_GPU != 0;
	}

	/**
	 *	@brief registers the common program termination signal handlers
	 */
	static void Register_SignalHandlers()
	{
		signal(SIGINT, &SigInt_Handler); // ctrl+c
		signal(SIGTERM, &SigInt_Handler);
#if defined(_WIN32) || defined(_WIN64)
		signal(SIGBREAK, &SigInt_Handler); // on closing console window
#else // _WIN32 || _WIN64
		signal(SIGQUIT, &SigInt_Handler); // ctrl+\, ctrl+d
		signal(SIGSTOP, &SigInt_Handler); // ctrl+z
		signal(SIGKILL, &SigInt_Handler); // ctrl+break
#endif // _WIN32 || _WIN64
		// avoid letting the user quit the program in the middle
		// of some computation, as that sometimes crashes the GPU
	}

protected:
	/**
	 *	@brief interrupt signal handler; quits if not using GPU, or schedules a quit for when done
	 *	@param[in] n_signal is signal number (unused)
	 */
	static void SigInt_Handler(int UNUSED(n_signal))
	{
		SafeExit(); // see if we are inside a GPU function right now - if not, may as well just quit
		printf("\ninterrupted: will quit as soon as GPU processing finishes\n");
		b_quit = true;
	}

	/**
	 *	@brief the only function that accesses the GPU use counter
	 *
	 *	This is a slightly awkward design of a single function that does two things.
	 *	This is dictated by the use of a critical section instead of a mutex (all the code
	 *	that reads / writes needs to be inside this section). This makes the code
	 *	slightly simpler.
	 *
	 *	@param[in] b_using_GPU is GPU use flag (if set, the counter is incremented, otherwise decremented)
	 *	@param[in] b_quit_if_not_using is quit flag (if set, b_using_GPU is *ignored* and the program
	 *		exits if the counter is zero)
	 *
	 *	@note This does not read (or write) the b_quit flag.
	 */
	static void Set_UsingGPU(bool b_using_GPU, bool b_quit_if_not_using = false)
	{
#ifdef _OPENMP
		#pragma omp critical
#endif // _OPENMP
		{
			if(b_quit_if_not_using) {
				if(!n_using_GPU)
					exit(-1);
				// just see if GPU is in use, don't increment / decrement
			} else {
				n_using_GPU += (b_using_GPU)? 1 : -1;
				_ASSERTE(n_using_GPU != size_t(-1)); // make sure it did not underflow
			}
		}
	}
};

bool CGPUGuard::b_quit = false;
size_t CGPUGuard::n_using_GPU = 0;

/**
 *	@def GPU_FUNCGUARD_PRETTYNAME_CAT30
 *	@brief helper macro for making the names of the guard variables more random in order to avoid name colissions
 */
#define GPU_FUNCGUARD_PRETTYNAME_CAT30(a,b,c) a##b##c

/**
 *	@def GPU_FUNCGUARD_PRETTYNAME_CAT3
 *	@brief helper macro for making the names of the guard variables more random in order to avoid name colissions
 */
#define GPU_FUNCGUARD_PRETTYNAME_CAT3(a,b,c) GPU_FUNCGUARD_PRETTYNAME_CAT30(a,b,c)

/**
 *	@def __GPU_Function
 *	@brief GPU usage scope marker
 */
#define __GPU_Function CGPUGuard GPU_FUNCGUARD_PRETTYNAME_CAT3(__gpu_function_guard_, __LINE__, __)

/**
 *	@def __GPU_ScopeCheck
 *	@brief GPU usage scope checker: asserts that it is in GPU scope, only in debug
 */
#define __GPU_ScopeCheck do { _ASSERTE(CGPUGuard::b_Using_GPU()); } while(0)

/*
 *								=== ~CGPUGuard ===
 */

/*
 *								=== TSharedGPUData ===
 */

/**
 *	@brief GPU data shared by multiple instances of the solvers
 */
struct TSharedGPUData { // basically a member storage of CLinearSolver_DenseGPU, except that I don't want to modify CLinearSolver_DenseGPU itself to avoid recompiling, and also it is practically a singleton
	int n_device; /**< @brief chosen CUDA device index */
	CUdevice t_device; /**< @brief chosen CUDA device handle */
	CUcontext t_context; /**< @brief CUDA context, using the chosen device */ // simplify management by using a single context: only one context is active in SLAM++ and no context switching needs to take place
	// needed for both CUDA and CULA

	/*CUdeviceptr dense_p_lambda_storage, dense_p_rhs_storage;
	size_t dense_n_lambda_size; // affects p_lambda_storage and p_rhs_storage*/
	// CULA does not have a context like CUDA / CUBLAS / CUsparse do

	/*cusparseHandle_t t_cusparse;
	cublasHandle_t t_cublas;
	cusparseMatDescr_t t_matrix_descr, t_sym_matrix_descr;
	cs *p_A, *p_B;
	int *csrRowPtrD, csrRowPtrD_size, nnzD;*/
	// moved to CLinearSolver_Schur_GPUBase members

	/**
	 *	@brief gets an instance of the object
	 *	@return Returns a reference to the instance of the object.
	 */
	static TSharedGPUData &r_GetInstance()
	{
		static TSharedGPUData instance; // singleton
		return instance;
	}

	/**
	 *	@brief inizialize CUDA and choose a device
	 *	@note This function throws std::runtime_error.
	 */
	void CUInit() // throw(std::runtime_error)
	{
		__GPU_Function;

		if(n_device != -1)
			return;
		// already inizialized

		int n_dev_num;
		int result, dn_result;
		if((result = cuInit(0)) != CUDA_SUCCESS || (dn_result =
		   cuDeviceGetCount(&n_dev_num)) != CUDA_SUCCESS || !n_dev_num) {
			if(result != CUDA_SUCCESS)
				fprintf(stderr, "error: cuInit() returns %d\n", result);
			else if(dn_result != CUDA_SUCCESS)
				fprintf(stderr, "error: cuDeviceGetCount() returns %d\n", dn_result);
			else
				fprintf(stderr, "error: cuDeviceGetCount() found no devices\n");
			throw std::runtime_error("cuInit() failed");
		}

		size_t n_max_mem;
		for(int i = 0; i < n_dev_num; ++ i) {
			if(cuDeviceGet(&t_device, i) != CUDA_SUCCESS)
				throw std::runtime_error("cuDeviceGet() failed");
			size_t n_mem;
			if(cuDeviceTotalMem(&n_mem, t_device) != CUDA_SUCCESS)
				throw std::runtime_error("cuDeviceTotalMem() failed");
			if(!i || n_mem > n_max_mem) {
				n_max_mem = n_mem;
				n_device = i;
			}
		}
		// choose the best device, based on memory size (tends to prefer Tesla devices before GTX ones)

		if(cuDeviceGet(&t_device, n_device) != CUDA_SUCCESS)
			throw std::runtime_error("cuDeviceGet() failed");
		// get handle to the chosen device

		char p_s_dev_name[1024];
		if(cuDeviceGetName(p_s_dev_name, sizeof(p_s_dev_name) / sizeof(char) - 1, t_device) != CUDA_SUCCESS)
			throw std::runtime_error("cuDeviceGetName() failed");
		printf("debug: chosen GPU device: \'%s\' (" PRIsizeB "B RAM)\n", p_s_dev_name, PRIsizeBparams(n_max_mem));
		// print the chosen device name

		if(cuCtxCreate(&t_context, CU_CTX_SCHED_AUTO, t_device) != CUDA_SUCCESS ||
		   cuCtxSetCurrent(t_context) != CUDA_SUCCESS)
			throw std::runtime_error("cuCtxCreate() failed");
		if(cuCtxSetSharedMemConfig(CU_SHARED_MEM_CONFIG_EIGHT_BYTE_BANK_SIZE) != CUDA_SUCCESS)
			throw std::runtime_error("cuCtxSetSharedMemConfig() failed");
		// create context and request 8-byte memory banks: will work mostly with doubles
	}

protected:
	/**
	 *	@brief default constructor; initializes empry data
	 */
	TSharedGPUData()
		:n_device(-1), t_device(0), t_context(0)
	{}

	/**
	 *	@brief destructor; frees GPU handles
	 */
	~TSharedGPUData()
	{
		__GPU_Function;

		if(t_context)
			cuCtxDestroy(t_context);
		// delete cuda context
	}

	/**
	 *	@brief copy-constructor (not implemented)
	 *	@param[in] r_other is other instance to copy from
	 */
	TSharedGPUData(const TSharedGPUData &r_other); // no-copy

	/**
	 *	@brief copy operator (not implemented)
	 *	@param[in] r_other is other instance to copy from
	 *	@return Returns reference to this.
	 */
	TSharedGPUData &operator =(const TSharedGPUData &r_other); // no-copy
};

static TSharedGPUData &gpu = TSharedGPUData::r_GetInstance(); /**< @brief GPU data shared by dense and Schur GPU solvers */

/*
 *								=== ~TSharedGPUData ===
 */

/*
 *								=== CLinearSolver_DenseGPU ===
 */

bool CLinearSolver_DenseGPU::b_cula_initialized = false;
size_t CLinearSolver_DenseGPU::n_instance_num = 0;

CLinearSolver_DenseGPU::CLinearSolver_DenseGPU() // throw(std::runtime_error)
{
	__GPU_Function;

	++ n_instance_num;
	if(!b_cula_initialized) {
		CGPUGuard::Register_SignalHandlers();
		// avoid letting the user quit the program in the middle
		// of some computation, as that sometimes crashes the GPU
		// now have to look at b_quit

		{
			gpu.CUInit();
			// initialize CUDA

			if(culaSelectDevice(gpu.n_device) != culaNoError)
				throw std::runtime_error("culaSelectDevice() failed");
			// "To bind without error, this function must be called before culaInitialize."
		}

		//printf("debug: culaInitialize()\n"); // seems to do that only once, when needed
		culaStatus s;
		switch(s = culaInitialize()) {
		case culaNoError:
			b_cula_initialized = true;
			break;
		case culaInsufficientRuntime:
			throw std::runtime_error("failed to initialize CULA: no compatible driver found");
		case culaInsufficientComputeCapability:
			throw std::runtime_error("failed to initialize CULA: no hardware with sufficient comp. cap. found");
		case culaNoHardware:
			throw std::runtime_error("failed to initialize CULA: no compatible hardware found");
		default:
			fprintf(stderr, "error: CULA error: \'%s\'\n", culaGetStatusString(s));
			{
				culaInfo i = culaGetErrorInfo();
				char p_s_error[1024];
				culaGetErrorInfoString(s, i, p_s_error, sizeof(p_s_error) / sizeof(char) - 1);
				fprintf(stderr, "error: CULA error: %d, \'%s\'\n", int(s), culaGetStatusString(s));
				fprintf(stderr, "error: CULA error info: %d, \'%s\'\n", int(i), p_s_error);
			}
			throw std::runtime_error("failed to initialize CULA: unspecified error");
		}
	}
	// the first one is in charge of initializing CULA
	// don't want to do lazy initialization here, just initialize it
	// when creating the solver so it is ready for when it is used

	{
		_ASSERTE(sizeof(CUdeviceptr) == sizeof(m_gpu.p_lambda_storage));
		_ASSERTE(sizeof(CUdeviceptr) == sizeof(m_gpu.p_rhs_storage));
		m_gpu.p_lambda_storage = 0;
		m_gpu.p_rhs_storage = 0;
		m_gpu.n_lambda_size = 0;
	}
	// in?tialize per-instance data (this is lazily allocated)
}

CLinearSolver_DenseGPU::~CLinearSolver_DenseGPU()
{
	__GPU_Function;

	{
		if(m_gpu.p_lambda_storage)
			cuMemFree((CUdeviceptr)m_gpu.p_lambda_storage);
		if(m_gpu.p_rhs_storage)
			cuMemFree((CUdeviceptr)m_gpu.p_rhs_storage);
	}
	// delete per instance data

	_ASSERTE(n_instance_num > 0); // this exists
	if(!(-- n_instance_num)) {
		if(b_cula_initialized) {
			b_cula_initialized = false;
			//printf("debug: culaShutdown()\n"); // seems to do that only once, when needed
			culaShutdown();
		}
	}
	// the last one takes care of the cleanup
}

void CLinearSolver_DenseGPU::Free_Memory()
{
	__GPU_Function;

	{
		Eigen::MatrixXd empty;
		std::swap(empty, m_t_lambda);
	}
	// free host memory

	{
		if(m_gpu.p_lambda_storage)
			cuMemFree((CUdeviceptr)m_gpu.p_lambda_storage);
		if(m_gpu.p_rhs_storage)
			cuMemFree((CUdeviceptr)m_gpu.p_rhs_storage);
		m_gpu.p_lambda_storage = 0;
		m_gpu.p_rhs_storage = 0;
		m_gpu.n_lambda_size = 0;
	}
	// free GPU memory
}

bool CLinearSolver_DenseGPU::Solve_PosDef(const CUberBlockMatrix &r_lambda, Eigen::VectorXd &r_eta) // throw(std::bad_alloc)
{
	__GPU_Function;

	r_lambda.Convert_to_Dense(m_t_lambda);

#if 0 && defined(_DEBUG) // seems to work
	Eigen::VectorXd ref;
	{
		typedef Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> _TyDecomposition; /**< @brief matrix factorization type (Cholesky) */
		_TyDecomposition m_t_factorization;

		m_t_factorization.compute(m_t_lambda);
		if(m_t_factorization.info() != Eigen::Success)
			return false;
		// Cholesky

		ref = m_t_factorization.solve(r_eta);
		// solve
	}
	// just a mockup
#endif // 0 && _DEBUG

	// on Venice:
	// Tesla K40 GPU:
	// device Cholesky		0.261068 sec / iteration, chi2 after 5 iters 234013913.66 (quite precise, I guess there is better IEEE-754)
	// Cholesky				0.297229 sec / iteration, chi2 after 5 iters 234013913.66

	// GTX 680 GPU:
	// device Cholesky		0.876127 sec / iteration, chi2 after 5 iters 234013920.81
	// Cholesky				1.018086 sec / iteration, chi2 after 5 iters 234013920.81
	// LU					2.977428 sec / iteration, chi2 after 5 iters 234084067.61 234013915.10 (for some reason very imprecise) | unable to reproduce
	// Chol + sw backsubst	1.389326 sec / iteration, chi2 after 5 iters 234013915.55
	// LU + sw backsubst	2.956637 sec / iteration, chi2 after 5 iters 234013920.97

	// CPU:
	// sparse				20 sec / iteration, chi2 after 5 iters 234013918.175338
	// LLT					9 sec / iteration, chi2 after 5 iters 234013924.73
	// LDLT					70 sec / iteration
	// PartialPivLU			19 sec / iteration // without parallelization
	// PartialPivLU			10 sec / iteration // with parallelization (8 cores), still slower than serial LLT, chi2 234013919.58
	// FullPivLU			486 sec / iteration // very slightly lower chi2
	// HouseholderQR		64 sec / iteration // pivotless

	size_t n = m_t_lambda.cols();
	_ASSERTE(m_t_lambda.rows() == n && r_eta.rows() == n);
	if(n > INT_MAX)
		throw std::runtime_error("matrix too large for CULA"); // likely wouldn't fit in the memory at this point
#if 1
	if(m_gpu.n_lambda_size < n) {
		if(m_gpu.p_lambda_storage)
			cuMemFree((CUdeviceptr)m_gpu.p_lambda_storage);
		if(m_gpu.p_rhs_storage)
			cuMemFree((CUdeviceptr)m_gpu.p_rhs_storage);
		if(cuMemAlloc((CUdeviceptr*)&m_gpu.p_lambda_storage, n * n * sizeof(double)) != CUDA_SUCCESS ||
		   cuMemAlloc((CUdeviceptr*)&m_gpu.p_rhs_storage, n * sizeof(double)) != CUDA_SUCCESS) {
			m_gpu.n_lambda_size = 0; // !!
			char p_s_error[256];
			sprintf(p_s_error, "CLinearSolver_DenseGPU::Solve_PosDef failed to allocate"
				" memory (%d x %d, " PRIsize ")", n, n, r_lambda.n_NonZero_Num());
			throw std::runtime_error(p_s_error);
		}
		m_gpu.n_lambda_size = n;
	}
	// (re)allocate storage for lambda and the RHS

	if(cuMemcpyHtoDAsync((CUdeviceptr)m_gpu.p_lambda_storage, m_t_lambda.data(), n * n * sizeof(double), 0) != CUDA_SUCCESS ||
	   cuMemcpyHtoDAsync((CUdeviceptr)m_gpu.p_rhs_storage, r_eta.data(), n * sizeof(double), 0) != CUDA_SUCCESS ||
	   cuCtxSynchronize() != CUDA_SUCCESS)
		throw std::runtime_error("CLinearSolver_DenseGPU::Solve_PosDef failed to upload data");
	// copy the data to the GPU

	culaStatus n_retval = culaDevicePosv('U', int(n), 1, (culaDeviceDouble*)m_gpu.p_lambda_storage,
		int(n), (culaDeviceDouble*)m_gpu.p_rhs_storage, int(n)); // Cholesky solve
	if(n_retval) {
		culaInfo n_info = culaGetErrorInfo();
		char p_s_error[1024];
		culaGetErrorInfoString(n_retval, n_info, p_s_error, sizeof(p_s_error) / sizeof(char) - 1);
		fprintf(stderr, "error: GPU solution returns %d (culaGetErrorInfo() returns %d, says \'%s\')\n",
			n_retval, n_info, p_s_error);
		throw std::runtime_error("CLinearSolver_DenseGPU::Solve_PosDef failed");
	}
	// factorize and solve on GPU

	if(cuMemcpyDtoH(r_eta.data(), (CUdeviceptr)m_gpu.p_rhs_storage, n * sizeof(double)) != CUDA_SUCCESS)
		throw std::runtime_error("CLinearSolver_DenseGPU::Solve_PosDef failed to download data");
	// copy only the RHS back, save time that would be spent transferring the factorization back - we don't need it
#elif 0
	culaStatus n_retval = culaPosv('U', int(n), 1, m_t_lambda.data(), int(n), r_eta.data(), int(n)); // Cholesky solve
	if(n_retval) {
		culaInfo n_info = culaGetErrorInfo();
		char p_s_error[1024];
		culaGetErrorInfoString(n_retval, n_info, p_s_error, sizeof(p_s_error) / sizeof(char) - 1);
		fprintf(stderr, "error: GPU solution returns %d (culaGetErrorInfo() returns %d, says \'%s\')\n",
			n_retval, n_info, p_s_error);
		throw std::runtime_error("CLinearSolver_DenseGPU::Solve_PosDef failed");
	}
	// factorize and solve on GPU
#elif 1
	m_t_lambda.triangularView<Eigen::StrictlyLower>() =
		m_t_lambda.triangularView<Eigen::StrictlyUpper>().transpose();
	// make the matrix symmetric for GPU LU

	std::vector<int> pivot(n);
	culaStatus n_retval = culaGesv(int(n), 1, m_t_lambda.data(), int(n), &pivot[0], r_eta.data(), int(n)); // LU solve
	if(n_retval) {
		culaInfo n_info = culaGetErrorInfo();
		char p_s_error[1024];
		culaGetErrorInfoString(n_retval, n_info, p_s_error, sizeof(p_s_error) / sizeof(char) - 1);
		fprintf(stderr, "error: GPU solution returns %d (culaGetErrorInfo() returns %d, says \'%s\')\n",
			n_retval, n_info, p_s_error);
		throw std::runtime_error("CLinearSolver_DenseGPU::Solve_PosDef failed");
	}
	// factorize and solve on GPU
#elif 0
	m_t_lambda.triangularView<Eigen::StrictlyLower>() =
		m_t_lambda.triangularView<Eigen::StrictlyUpper>().transpose();
	// make the matrix symmetric for GPU LU

	std::vector<int> pivot_order(n);
	culaStatus n_retval = culaGetrf(int(n), int(n), m_t_lambda.data(), int(n), &pivot_order[0]); // LU
	if(n_retval) {
		culaInfo n_info = culaGetErrorInfo();
		char p_s_error[1024];
		culaGetErrorInfoString(n_retval, n_info, p_s_error, sizeof(p_s_error) / sizeof(char) - 1);
		fprintf(stderr, "error: GPU solution returns %d (culaGetErrorInfo() returns %d, says \'%s\')\n",
			n_retval, n_info, p_s_error);
		throw std::runtime_error("CLinearSolver_DenseGPU::Solve_PosDef failed");
	}
	// factorize on GPU

	for(size_t i = 0; i < n; ++ i)
		std::swap(r_eta(i), r_eta(pivot_order[i] - 1));
	// permute inplace using row swaps

	r_eta = m_t_lambda.triangularView<Eigen::UnitLower>().solve(r_eta); // backsubstitution
	r_eta = m_t_lambda.triangularView<Eigen::Upper>().solve(r_eta); // backsubstitution
	// solve

	// culaGetrf() performs only row pivoting, no unpermuting required at this point
#else // 1
	int n_retval = culaPotrf('U', int(n), m_t_lambda.data(), int(n)); // Cholesky
	if(n_retval == 8) {
		fprintf(stderr, "error: GPU solution returns %d (not pos def; culaGetErrorInfo() returns %d)\n", n_retval, culaGetErrorInfo());
		throw std::runtime_error("CLinearSolver_DenseGPU::Solve_PosDef: not pos def");
	} else if(n_retval) {
		fprintf(stderr, "error: GPU solution returns %d (culaGetErrorInfo() returns %d)\n", n_retval, culaGetErrorInfo());
		throw std::runtime_error("CLinearSolver_DenseGPU::Solve_PosDef failed");
	}
	r_eta = m_t_lambda.triangularView<Eigen::Upper>().transpose().solve(r_eta); // backsubstitution
	r_eta = m_t_lambda.triangularView<Eigen::Upper>().solve(r_eta); // backsubstitution
#endif // 1
	// do it on a GPU?

#if 0 && defined(_DEBUG) // seems to work
	printf("GPU solution error: %g\n", (ref - r_eta).norm());
#endif // 0 && _DEBUG

	return true;
}

/*
 *								=== ~CLinearSolver_DenseGPU ===
 */

/*
 *								=== CLinearSolver_Schur_GPUBase ===
 */

bool CLinearSolver_Schur_GPUBase::b_cuda_initialized = false;
size_t CLinearSolver_Schur_GPUBase::n_instance_num = 0;

CLinearSolver_Schur_GPUBase::CLinearSolver_Schur_GPUBase()
{
	__GPU_Function;

	++ n_instance_num;
	if(!b_cuda_initialized) {
		CGPUGuard::Register_SignalHandlers();
		// avoid letting the user quit the program in the middle
		// of some computation, as that sometimes crashes the GPU
		// now have to look at b_quit

		{
			gpu.CUInit();
			// initialize CUDA
		}
	}

	{
		if(cusparseCreate((cusparseHandle_t*)&m_gpu.t_cusparse) != CUSPARSE_STATUS_SUCCESS)
			throw std::runtime_error("cusparseCreate() failed");
		// create a cusparse context

		if(cusparseCreateMatDescr((cusparseMatDescr_t*)&m_gpu.t_matrix_descr) != CUSPARSE_STATUS_SUCCESS ||
		   cusparseSetMatType((cusparseMatDescr_t)m_gpu.t_matrix_descr, CUSPARSE_MATRIX_TYPE_GENERAL) != CUSPARSE_STATUS_SUCCESS ||
		   cusparseSetMatIndexBase((cusparseMatDescr_t)m_gpu.t_matrix_descr, CUSPARSE_INDEX_BASE_ZERO) != CUSPARSE_STATUS_SUCCESS)
			throw std::runtime_error("cusparseCreateMatDescr() failed");
		if(cusparseCreateMatDescr((cusparseMatDescr_t*)&m_gpu.t_sym_matrix_descr) != CUSPARSE_STATUS_SUCCESS ||
		   cusparseSetMatType((cusparseMatDescr_t)m_gpu.t_sym_matrix_descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC) != CUSPARSE_STATUS_SUCCESS ||
		   cusparseSetMatIndexBase((cusparseMatDescr_t)m_gpu.t_sym_matrix_descr, CUSPARSE_INDEX_BASE_ZERO) != CUSPARSE_STATUS_SUCCESS)
			throw std::runtime_error("cusparseCreateMatDescr() failed");
		// create matrix descriptors (will recycle that)

		if(cusparseSetPointerMode((cusparseHandle_t)m_gpu.t_cusparse, CUSPARSE_POINTER_MODE_HOST) != CUSPARSE_STATUS_SUCCESS)
			throw std::runtime_error("cusparseSetPointerMode() failed");
		// cusparse misc settings

		if(cublasCreate((cublasHandle_t*)&m_gpu.t_cublas) != CUBLAS_STATUS_SUCCESS)
			throw std::runtime_error("cublasCreate() failed");
		// initialize cublas

		if(cublasSetPointerMode((cublasHandle_t)m_gpu.t_cublas, CUBLAS_POINTER_MODE_HOST) != CUBLAS_STATUS_SUCCESS)
			throw std::runtime_error("cublasSetPointerMode() failed");
		// cublas misc settings
	}
	// the rest of the state is per instance
}

CLinearSolver_Schur_GPUBase::~CLinearSolver_Schur_GPUBase()
{
	__GPU_Function;

	{
		if(m_gpu.p_A)
			cs_spfree(m_gpu.p_A);
		if(m_gpu.p_B)
			cs_spfree(m_gpu.p_B);
		if(m_gpu.csrRowPtrD)
			cuMemFree((CUdeviceptr)m_gpu.csrRowPtrD);
		if(m_gpu.t_matrix_descr)
			cusparseDestroyMatDescr((cusparseMatDescr_t)m_gpu.t_matrix_descr);
		if(m_gpu.t_sym_matrix_descr)
			cusparseDestroyMatDescr((cusparseMatDescr_t)m_gpu.t_sym_matrix_descr);
		// delete GPU resources

		if(m_gpu.t_cublas)
			cublasDestroy((cublasHandle_t)m_gpu.t_cublas);
		if(m_gpu.t_cusparse)
			cusparseDestroy((cusparseHandle_t)m_gpu.t_cusparse);
		// delete context
	}
	// delete per instance data

	_ASSERTE(n_instance_num > 0); // this exists
	if(!(-- n_instance_num)) {
		// nothing really happens, CUDA remains initialized
	}
	// the last one takes care of the cleanup
}

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
cs *cs_spalloc32(csi m, csi n, csi nzmax, csi values, csi triplet);
// in BlockMatrix.cpp, needed for x64 builds as cusparse is always using 32-bit indices
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64

bool CLinearSolver_Schur_GPUBase::GPUSolve(const CUberBlockMatrix &r_lambda,
	Eigen::VectorXd &r_v_eta, bool b_keep_ordering, size_t n_landmark_dimension,
	std::vector<double> &m_double_workspace, const std::vector<size_t> &m_order,
	const size_t m_n_matrix_cut, _TyBaseSolver &m_linear_solver)
{
	__GPU_Function;

	_ASSERTE(r_lambda.b_SymmetricLayout()); // pos-def is supposed to be symmetric
	_ASSERTE(r_v_eta.rows() == r_lambda.n_Row_Num()); // make sure the vector has correct dimension

	_ASSERTE(r_lambda.b_SymmetricLayout());
	_ASSERTE(r_lambda.n_BlockColumn_Num() == m_order.size());
	_ASSERTE((m_order.empty() && !m_n_matrix_cut) || (!m_order.empty() &&
		m_n_matrix_cut > 0 && m_n_matrix_cut < SIZE_MAX && m_n_matrix_cut + 1 < m_order.size()));
	_ASSERTE(r_v_eta.rows() == r_lambda.n_Column_Num());

#ifdef __SCHUR_PROFILING
	CDeltaTimer dt;
#endif // __SCHUR_PROFILING

	CUberBlockMatrix lambda_perm;
	r_lambda.Permute_UpperTriangular_To(lambda_perm, &m_order[0], m_order.size(), true);
	// reorder the matrix

#ifdef __SCHUR_PROFILING
	double f_reperm_time = dt.f_Time();
#endif // __SCHUR_PROFILING

	const size_t n = lambda_perm.n_BlockColumn_Num();

	CUberBlockMatrix A, U, C;
	const size_t n_matrix_cut = m_n_matrix_cut; // antialiass
	lambda_perm.SliceTo(A, 0, n_matrix_cut, 0, n_matrix_cut, true);
	lambda_perm.SliceTo(U, 0, n_matrix_cut, n_matrix_cut, n, true);
	lambda_perm.SliceTo(C, n_matrix_cut, n, n_matrix_cut, n, true);
	// cut Lambda matrix into pieces
	// \lambda = | A U |
	//           | V C |

	const size_t n_rhs_vector_size = r_lambda.n_Column_Num();
	const size_t n_pose_vector_size = A.n_Column_Num(); // 6 * m_n_matrix_cut;
	const size_t n_landmark_vector_size = U.n_Column_Num(); // 3 * (n - m_n_matrix_cut);
	// not block columns! element ones

#ifdef __SCHUR_PROFILING
	double f_slice_time = dt.f_Time();
	double f_transpose_time = 0;//dt.f_Time();
#endif // __SCHUR_PROFILING

#if 0
	CUberBlockMatrix C_inv;
	if(n_landmark_dimension == 3) // this can be anticipated
		C_inv.InverseOf_Symmteric_FBS<MakeTypelist_Safe((fbs_ut::CCTSize2D<3, 3>))>(C); // C is block diagonal (should also be symmetric)
	else
		C_inv.InverseOf_Symmteric(C);
#else // 0
	_ASSERTE(C.b_BlockDiagonal()); // it is, unless the ordering is bad
	CUberBlockMatrix &C_inv = C; // can do it inplace
	if(n_landmark_dimension == 3) // this can be anticipated in BA
		C_inv.InverseOf_BlockDiag_FBS_Parallel<MakeTypelist_Safe((fbs_ut::CCTSize2D<3, 3>))>(C); // faster version, slightly less general
	else
		C_inv.InverseOf_Symmteric(C); // there is no non-FBS InverseOf_BlockDiag(), just use this one instead 
#endif // 0
	// inverse of C

#ifdef __SCHUR_PROFILING
	double f_inverse_time = dt.f_Time();
#endif // __SCHUR_PROFILING

	/*CUberBlockMatrix unity, u_ref;
	unity.ProductOf(C, C_inv);
	unity.CopyLayoutTo(u_ref);
	u_ref.SetIdentity();
	unity.AddTo(u_ref, -1);
	double f_error = u_ref.f_Norm();
	fprintf(stderr, "error of matrix inverse is: %g\n", f_error);*/
	// check inverse

#ifdef __GPU_SCHUR_VERIFY_RESULT
	CUberBlockMatrix U_Cinv_ref, minus_schur_ref;
	{
		U.MultiplyToWith(U_Cinv_ref, C_inv);
#ifndef __GPU_SCHUR_EASY_PROD_ONLY
		CUberBlockMatrix V;
		U.TransposeTo(V); // because lower-triangular of lambda is not calculated
		U_Cinv_ref.MultiplyToWith(minus_schur_ref, V);
#endif // !__GPU_SCHUR_EASY_PROD_ONLY
	}
	// debug - calculate product on the CPU
#endif // __GPU_SCHUR_VERIFY_RESULT

#ifdef __SCHUR_PROFILING
	double f_mul0_time; // inside
#endif // __SCHUR_PROFILING
	CUberBlockMatrix minus_U_Cinv, schur_compl;

#if 0 // all prods on CPU (overrides __GPU_SCHUR_NO_PRODS in LinearSolver_Schur.h without having to rebuild all; note that it is non-FBS and therefore slower)
	C_inv.Scale(-1);
	//U.MultiplyToWith(minus_U_Cinv, C_inv);
	U.MultiplyToWith_FBS<MakeTypelist_Safe((Eigen::Matrix<double, 6, 3>)),
		MakeTypelist_Safe((Eigen::Matrix<double, 3, 3>))>(minus_U_Cinv, C_inv); // use FBS here (guess the block sizes)
#ifdef __SCHUR_PROFILING
	f_mul0_time = dt.f_Time();
#endif // __SCHUR_PROFILING
	printf("UBlock GEMM1 time: %f\n", f_mul0_time);
	{
		CUberBlockMatrix V;
		U.TransposeTo(V); // because lower-triangular of lambda is not calculated
		//minus_U_Cinv.MultiplyToWith_FBS<>(schur_compl, V);
		minus_U_Cinv.MultiplyToWith_FBS<MakeTypelist_Safe((Eigen::Matrix<double, 6, 3>)),
			MakeTypelist_Safe((Eigen::Matrix<double, 3, 6>))>(schur_compl, V);
#ifdef __SCHUR_PROFILING
		printf("UBlock GEMM2 time: %f (full matrix, not only upper half, " PRIsize " NNZ blocks)\n", dt.f_Time(),
			schur_compl.n_Block_Num());
#endif // __SCHUR_PROFILING
		// time with perfect ordering

		/*{
			cs *p_blayout = V.p_BlockStructure_to_Sparse(cs_spalloc(V.n_BlockRow_Num(),
				V.n_BlockColumn_Num(), V.n_Block_Num(), 1, 0));
			if(!p_blayout)
				throw std::runtime_error("V.p_BlockStructure_to_Sparse() failed");
			std::vector<size_t> row_cumsums(p_blayout->m), col_cumsums(p_blayout->n), workspace;
			for(size_t i = 0, n = row_cumsums.size(); i < n; ++ i)
				row_cumsums[i] = i + 1;
			for(size_t i = 0, n = col_cumsums.size(); i < n; ++ i)
				col_cumsums[i] = i + 1;
			CUberBlockMatrix V_bs(row_cumsums.begin(), row_cumsums.end(),
				col_cumsums.begin(), col_cumsums.end());
			if(!V_bs.From_Sparse(0, 0, p_blayout, false, workspace))
				throw std::runtime_error("V_bs.From_Sparse() failed");
			cs_spfree(p_blayout);
			V_bs.Save_MatrixMarket("G:\\uflsmc\\sparse\\mat\\Venice\\V_struct.mtx");
		}
		// save V's structure as .mtx

		V.Save_MatrixMarket("G:\\uflsmc\\sparse\\mat\\Venice\\V_full.mtx");
		// save full V as .mtx

		{
			cs *p_blayout = minus_U_Cinv.p_BlockStructure_to_Sparse(cs_spalloc(minus_U_Cinv.n_BlockRow_Num(),
				minus_U_Cinv.n_BlockColumn_Num(), minus_U_Cinv.n_Block_Num(), 1, 0));
			if(!p_blayout)
				throw std::runtime_error("UCinv.p_BlockStructure_to_Sparse() failed");
			std::vector<size_t> row_cumsums(p_blayout->m), col_cumsums(p_blayout->n), workspace;
			for(size_t i = 0, n = row_cumsums.size(); i < n; ++ i)
				row_cumsums[i] = i + 1;
			for(size_t i = 0, n = col_cumsums.size(); i < n; ++ i)
				col_cumsums[i] = i + 1;
			CUberBlockMatrix V_bs(row_cumsums.begin(), row_cumsums.end(),
				col_cumsums.begin(), col_cumsums.end());
			if(!V_bs.From_Sparse(0, 0, p_blayout, false, workspace))
				throw std::runtime_error("UCinv_bs.From_Sparse() failed");
			cs_spfree(p_blayout);
			V_bs.Save_MatrixMarket("G:\\uflsmc\\sparse\\mat\\Venice\\UCinv_struct.mtx");
		}
		// save minus_U_Cinv's structure as .mtx

		minus_U_Cinv.Save_MatrixMarket("G:\\uflsmc\\sparse\\mat\\Venice\\UCinv_full.mtx");
		// save full minus_U_Cinv as .mtx

		for(size_t n_slice_num = 4; n_slice_num < 64; n_slice_num *= 2) {
			for(int i = 0; i < n_slice_num; ++ i) {
				const size_t n = minus_U_Cinv.n_BlockColumn_Num();
				const size_t b = (n * i) / n_slice_num;
				const size_t e = (i + 1 < n_slice_num)? (n * (i + 1)) / n_slice_num : n;

				CUberBlockMatrix slice;
				minus_U_Cinv.SliceTo(slice, 0, minus_U_Cinv.n_BlockRow_Num(), b, e, true);
				char p_s_name[1024];
				sprintf(p_s_name, "G:\\uflsmc\\sparse\\mat\\Venice\\UCinv_slice_%dof%d.mtx",
					int(i + 1), int(n_slice_num));
				slice.Save_MatrixMarket(p_s_name); 
			}
		}*/
		// save different granularity slices of minus_U_Cinv as .mtx

		/*
#ifdef __SCHUR_PROFILING
		dt.f_Time(); // don't count the stuff before
#endif // __SCHUR_PROFILING

		cs *p_U = U.p_Convert_to_Sparse();
		cs *p_C_inv = C_inv.p_Convert_to_Sparse();
		cs *p_U_Cinv = cs_multiply(p_U, p_C_inv);

#ifdef __SCHUR_PROFILING
		printf("CSparse GEMM1 takes %f sec\n", dt.f_Time());
#endif // __SCHUR_PROFILING

		cs_spfree(p_U);
		cs_spfree(p_C_inv);
		cs_spfree(p_U_Cinv);

#ifdef __SCHUR_PROFILING
		dt.f_Time(); // don't count spfree()
#endif // __SCHUR_PROFILING

		cs *p_minus_U_Cinv = minus_U_Cinv.p_Convert_to_Sparse();
		cs *p_V = V.p_Convert_to_Sparse();
		cs *p_SC = cs_multiply(p_minus_U_Cinv, p_V);

#ifdef __SCHUR_PROFILING
		printf("CSparse GEMM2 takes %f sec\n", dt.f_Time());
#endif // __SCHUR_PROFILING

		cs_spfree(p_minus_U_Cinv);
		cs_spfree(p_V);
		cs_spfree(p_SC);

#ifdef __SCHUR_PROFILING
		dt.f_Time(); // don't count spfree()
#endif // __SCHUR_PROFILING
		*/
		// time CSparse
	}
	C_inv.Scale(-1);
#else // 0
	// GPU gemm
	{
		cs *&p_A = m_gpu.p_A, *&p_B = m_gpu.p_B;
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
		p_A = C_inv.p_Convert_to_Sparse32(p_A);
		p_B = U.p_Convert_to_Sparse32(p_B);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
		p_A = C_inv.p_Convert_to_Sparse(p_A);
		p_B = U.p_Convert_to_Sparse(p_B);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
		// swap order! it does it in transpose and transposes the result

		const size_t nnzA = ((uint32_t*)p_A->p)[p_A->n], nnzB = ((uint32_t*)p_B->p)[p_B->n];
		const size_t m = p_A->n; // A is CUSPARSE_OPERATION_NON_TRANSPOSE
		const size_t n = p_B->m; // B is CUSPARSE_OPERATION_NON_TRANSPOSE
		const size_t k = p_A->m; // A is CUSPARSE_OPERATION_NON_TRANSPOSE
		// m, n are also dimensions of C

		/*printf("A is " PRIsize " x " PRIsize "\n", p_A->m, p_A->n);
		printf("B is " PRIsize " x " PRIsize "\n", p_B->m, p_B->n);
		printf("entered m, n, k: " PRIsize ", " PRIsize ", " PRIsize "\n", m, n, k);
		printf("A has " PRIsize " nnz, B has " PRIsize " nnz\n", U.n_NonZero_Num(), C_inv.n_NonZero_Num());
		U.MultiplyToWith(minus_U_Cinv, C_inv);
		printf("A * B has " PRIsize " nnz (int max is %d; %d)\n", minus_U_Cinv.n_NonZero_Num(),
			INT_MAX, minus_U_Cinv.n_NonZero_Num() < INT_MAX);
		printf("A * B is " PRIsize " x " PRIsize "\n", minus_U_Cinv.n_Row_Num(), minus_U_Cinv.n_Column_Num());
		printf("A has " PRIsize " nnz < INT_MAX, B has " PRIsize " nnz < INT_MAX\n",
			U.n_NonZero_Num() < INT_MAX, C_inv.n_NonZero_Num() < INT_MAX);
		CUberBlockMatrix minus_U_Cinv_ref;
		minus_U_Cinv_ref.Swap(minus_U_Cinv);
		minus_U_Cinv.Clear();
		// debug - matrix sizes
		printf("A has " PRIsize " nnz, B has " PRIsize " nnz\n", nnzA, nnzB);*/
		// debug

		int *csrRowPtrC;
		/*{
			int result = cuMemAlloc((CUdeviceptr*)&csrRowPtrC, sizeof(int) * (m + 1));
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("failed to allocate product matrix rowptr on GPU");
		}*/
		// the same as csrRowPtrB

		int *csrRowPtrA, *csrColIndA;
		double *csrValA;
		{
			int result = cuMemAlloc((CUdeviceptr*)&csrRowPtrA, sizeof(int) * (p_A->n + 1));
			result |= cuMemAlloc((CUdeviceptr*)&csrColIndA, sizeof(int) * nnzA);
			result |= cuMemAlloc((CUdeviceptr*)&csrValA, sizeof(double) * nnzA);
			if(result)
				fprintf(stderr, "error: cuMemAlloc() failed (%d, " PRIsize ")\n", result, nnzA);
			result |= cuMemcpyHtoDAsync((CUdeviceptr)csrRowPtrA, p_A->p, sizeof(int) * (p_A->n + 1), 0);
			result |= cuMemcpyHtoDAsync((CUdeviceptr)csrColIndA, p_A->i, sizeof(int) * nnzA, 0);
			result |= cuMemcpyHtoDAsync((CUdeviceptr)csrValA, p_A->x, sizeof(double) * nnzA, 0);
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("failed to send the A matrix to GPU");
		}
		int *csrRowPtrB, *csrColIndB;
		double *csrValB;
		{
			int result = cuMemAlloc((CUdeviceptr*)&csrRowPtrB, sizeof(int) * (p_B->n + 1));
			result |= cuMemAlloc((CUdeviceptr*)&csrColIndB, sizeof(int) * nnzB);
			result |= cuMemAlloc((CUdeviceptr*)&csrValB, sizeof(double) * nnzB);
			if(result)
				fprintf(stderr, "error: cuMemAlloc() failed (%d, " PRIsize ")\n", result, nnzB);
			result |= cuMemcpyHtoDAsync((CUdeviceptr)csrRowPtrB, p_B->p, sizeof(int) * (p_B->n + 1), 0);
			result |= cuMemcpyHtoDAsync((CUdeviceptr)csrColIndB, p_B->i, sizeof(int) * nnzB, 0);
			result |= cuMemcpyHtoDAsync((CUdeviceptr)csrValB, p_B->x, sizeof(double) * nnzB, 0);
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("failed to send the B matrix to GPU");
		}
		// B is symmetric - the same in SCR and CSC
		// todo - see which combination of transpose flags yields faster results

		if(cuCtxSynchronize() != CUDA_SUCCESS)
			throw std::runtime_error("cuCtxSynchronize() failed " FILE_LINE);

		//printf("debug: have matrices in GPU\n"); // debug

		// p_A, p_B are not needed anymore // todo - reuse for the second product

		{
			const double f_x = -1;
			_ASSERTE(nnzA < INT_MAX);
			cublasDscal((cublasHandle_t)m_gpu.t_cublas, int(nnzA), &f_x, csrValA, 1);
		}
		// use cublasDscal() to flip signs on A (A = C^-1)

		size_t nnzC/*, nnzCstat*/;
#if 1
		nnzC = nnzB;
		/*{
			int result = cuMemcpyDtoD((CUdeviceptr)csrRowPtrC, (CUdeviceptr)csrRowPtrB, (p_B->n + 1) * sizeof(int));
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("cuMemcpyDtoD() failed");
			// this is the same thing, the result has the same structure as B
		}*/
		csrRowPtrC = csrRowPtrB; // in fact, just use the same array, it will not get overwritten
#else // 1
		{
			_ASSERTE(m <= INT_MAX && n <= INT_MAX && k <= INT_MAX);
			_ASSERTE(nnzA <= INT_MAX && nnzB <= INT_MAX);
			int nnzCi;
			int status = cusparseXcsrgemmNnz((cusparseHandle_t)m_gpu.t_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
				CUSPARSE_OPERATION_NON_TRANSPOSE, int(m), int(n), int(k), (cusparseMatDescr_t)m_gpu.t_matrix_descr,
				int(nnzA), csrRowPtrA, csrColIndA, (cusparseMatDescr_t)m_gpu.t_matrix_descr, int(nnzB), csrRowPtrB,
				csrColIndB, (cusparseMatDescr_t)m_gpu.t_matrix_descr, csrRowPtrC, &nnzCi);
			nnzC = nnzCi; // !!
			// A is symmetric (A is C^-1)

			if(status != CUSPARSE_STATUS_SUCCESS)
				throw std::runtime_error("cusparseXcsrgemmNnz() failed");
			//nnzCstat = status; // debug
		}
		// calculate the amount of NNZ in C
		// t_odo - recaculate this only if !b_keep_ordering, otherwise reuse
		// t_odo - here we are multiplying a block matrix by a block diagonal matrix
		//		  the sparsity is therefore always the same, as is the nnz

		if(cuCtxSynchronize() != CUDA_SUCCESS)
			throw std::runtime_error("cuCtxSynchronize() failed " FILE_LINE);
		// about to read the NNZs back for cuMemAlloc()
#endif // 1

		//printf("debug: the product has %d nnz (%d)\n", nnzC, nnzCstat); // debug

		int *csrColIndC;
		double *csrValC;
		{
			int result = CUDA_SUCCESS;
#ifdef _DEBUG
			int baseC, _nnzC;
			result |= cuMemcpyDtoH(&_nnzC, (CUdeviceptr)(csrRowPtrC + m), sizeof(int)); // read the nnz
			_ASSERTE(_nnzC == nnzC);
			result |= cuMemcpyDtoH(&baseC, (CUdeviceptr)csrRowPtrC, sizeof(int)); // read the base
			_ASSERTE(!baseC); // debug
#endif // _DEBUG
			result |= cuMemAlloc((CUdeviceptr*)&csrColIndC, sizeof(int) * nnzC);
			result |= cuMemAlloc((CUdeviceptr*)&csrValC, sizeof(double) * nnzC);
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("failed to allocate product matrix on GPU");
		}
		// alloc C

		//printf("debug: allocated product storage\n"); // debug

		{
			_ASSERTE(m <= INT_MAX && n <= INT_MAX && k <= INT_MAX);
			_ASSERTE(nnzA <= INT_MAX && nnzB <= INT_MAX);
			int status = cusparseDcsrgemm((cusparseHandle_t)m_gpu.t_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
				CUSPARSE_OPERATION_NON_TRANSPOSE, int(m), int(n), int(k), (cusparseMatDescr_t)m_gpu.t_matrix_descr,
				int(nnzA), csrValA, csrRowPtrA, csrColIndA, (cusparseMatDescr_t)m_gpu.t_matrix_descr, int(nnzB), csrValB,
				csrRowPtrB, csrColIndB, (cusparseMatDescr_t)m_gpu.t_matrix_descr, csrValC, csrRowPtrC, csrColIndC);
			// A is symmetric (A is C^-1)

			if(status != CUSPARSE_STATUS_SUCCESS)
				throw std::runtime_error("cusparseDcsrgemm() failed");
		}
		// gemm

		//printf("debug: GEMM finished\n"); // debug

		if(cuCtxSynchronize() != CUDA_SUCCESS)
			throw std::runtime_error("cuCtxSynchronize() failed " FILE_LINE);

		{
			int result = cuMemFree((CUdeviceptr)csrRowPtrA);
			result |= cuMemFree((CUdeviceptr)csrColIndA);
			result |= cuMemFree((CUdeviceptr)csrValA);
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("failed to free matrix storage on GPU");
		}
		// can delete A at this point

#ifndef __GPU_SCHUR_EASY_PROD_ONLY

		const size_t m1 = /*p_B->m*/U.n_Row_Num(); // B is CUSPARSE_OPERATION_TRANSPOSE
		const size_t n1 = /*p_C->*/n; // C is CUSPARSE_OPERATION_NON_TRANSPOSE
		const size_t k1 = /*p_B->n*/U.n_Column_Num(); // B is CUSPARSE_OPERATION_TRANSPOSE

		int *&csrRowPtrD = m_gpu.csrRowPtrD, &nnzD = m_gpu.nnzD;
		if(!csrRowPtrD || !b_keep_ordering) {
			if(!csrRowPtrD || m_gpu.csrRowPtrD_size != (n + 1)) {
				if(csrRowPtrD)
					cuMemFree((CUdeviceptr)csrRowPtrD);
				m_gpu.csrRowPtrD_size = 0;
				int result = cuMemAlloc((CUdeviceptr*)&csrRowPtrD, sizeof(int) * (n + 1));
				if(result != CUDA_SUCCESS)
					throw std::runtime_error("failed to allocate product matrix rowptr on GPU");
				m_gpu.csrRowPtrD_size = n + 1;
			}
			// alloc

			_ASSERTE(m1 <= INT_MAX && n1 <= INT_MAX && k1 <= INT_MAX);
			_ASSERTE(nnzB <= INT_MAX && nnzC <= INT_MAX);
			int status = cusparseXcsrgemmNnz((cusparseHandle_t)m_gpu.t_cusparse, CUSPARSE_OPERATION_TRANSPOSE,
				CUSPARSE_OPERATION_NON_TRANSPOSE, int(m1), int(n1), int(k1), (cusparseMatDescr_t)m_gpu.t_matrix_descr,
				int(nnzB), csrRowPtrB, csrColIndB, (cusparseMatDescr_t)m_gpu.t_matrix_descr, int(nnzC), csrRowPtrC, csrColIndC,
				(cusparseMatDescr_t)m_gpu.t_matrix_descr, csrRowPtrD, &nnzD); // t_odo - D = U(C^-1)V is symmetric, try what it does // fails
			if(status != CUSPARSE_STATUS_SUCCESS)
				throw std::runtime_error("cusparseXcsrgemmNnz() failed");
			// fill
		}
		// recalculate this only if ordering changes

		if(cuCtxSynchronize() != CUDA_SUCCESS)
			throw std::runtime_error("cuCtxSynchronize() failed " FILE_LINE);
		// about to read the NNZs back for cuMemAlloc()

		/*printf("the product of all 3 matrices will be " PRIsize " x " PRIsize " and will have " PRIsize " nnz\n",
			m1, n1, nnzD);*/

		int *csrColIndD; // todo - could cache this nut it probably chanes often
		double *csrValD;
		{
			int result = CUDA_SUCCESS;
#ifdef _DEBUG
			int baseD, _nnzD;
			result |= cuMemcpyDtoH(&_nnzD, (CUdeviceptr)(csrRowPtrD + m1), sizeof(int)); // read the nnz
			_ASSERTE(_nnzD == nnzD);
			result |= cuMemcpyDtoH(&baseD, (CUdeviceptr)csrRowPtrD, sizeof(int)); // read the base
			_ASSERTE(!baseD); // debug
#endif // _DEBUG
			result |= cuMemAlloc((CUdeviceptr*)&csrColIndD, sizeof(int) * nnzD);
			result |= cuMemAlloc((CUdeviceptr*)&csrValD, sizeof(double) * nnzD);
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("failed to allocate product matrix on GPU");
		}
		// alloc D

		//printf("debug: allocated product storage\n"); // debug

		{
			int status = cusparseDcsrgemm((cusparseHandle_t)m_gpu.t_cusparse, CUSPARSE_OPERATION_TRANSPOSE,
				CUSPARSE_OPERATION_NON_TRANSPOSE, int(m1), int(n1), int(k1), (cusparseMatDescr_t)m_gpu.t_matrix_descr,
				int(nnzB), csrValB, csrRowPtrB, csrColIndB, (cusparseMatDescr_t)m_gpu.t_matrix_descr, int(nnzC), csrValC,
				csrRowPtrC, csrColIndC, (cusparseMatDescr_t)m_gpu.t_matrix_descr, csrValD, csrRowPtrD, csrColIndD);
			// todo - D = U(C^-1)V is symmetric, try what it does

			if(status != CUSPARSE_STATUS_SUCCESS)
				throw std::runtime_error("cusparseDcsrgemm() failed");
		}
		// gemm

#endif // !__GPU_SCHUR_EASY_PROD_ONLY

		//printf("debug: GEMM finished\n"); // debug

		{
			//U.CopyLayoutTo(minus_U_Cinv);
			minus_U_Cinv = U; // also preallocates the required blocks
			//U.Swap(minus_U_Cinv); // U is needed below
			// prepare block layout

			if(p_B->n != m || p_B->m != n || p_B->nzmax < 0 || size_t(p_B->nzmax) < nnzC)
				throw std::runtime_error("can't reuse p_B");
			cs *p_C = p_B;//cs_spalloc32(n, m, nnzC, 1, 0); // can reuse the storage
			{
				int result = cuMemcpyDtoHAsync(p_C->p, (CUdeviceptr)csrRowPtrC, sizeof(int) * (m + 1), 0);
				result |= cuMemcpyDtoHAsync(p_C->i, (CUdeviceptr)csrColIndC, sizeof(int) * nnzC, 0);
				result |= cuMemcpyDtoHAsync(p_C->x, (CUdeviceptr)csrValC, sizeof(double) * nnzC, 0);
				if(result != CUDA_SUCCESS)
					throw std::runtime_error("failed to download GEMM result from the GPU");
			}
			// copy data but don't free the product yet, the last GEMM is still referencing it

			if(cuCtxSynchronize() != CUDA_SUCCESS)
				throw std::runtime_error("cuCtxSynchronize() failed " FILE_LINE);

			std::vector<size_t> workspace0;
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
			if(!minus_U_Cinv.From_Sparse32_Parallel(0, 0, p_C, false, workspace0)) {
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
			if(!minus_U_Cinv.From_Sparse_Parallel(0, 0, p_C, false, workspace0)) {
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
				fprintf(stderr, "error: C is " PRIsize " x " PRIsize ", UC--1 expects " PRIsize " x " PRIsize "\n",
					p_C->m, p_C->n, minus_U_Cinv.n_Row_Num(), minus_U_Cinv.n_Column_Num());
				throw std::runtime_error("From_Sparse_Parallel() failed");
			}
			// make it back into a block matrix

			//cs_spfree(p_C);
			// not needed anymore
		}
		// overlap this with the GPU computation
		// contains synchronization

		{
			int result = cuMemFree((CUdeviceptr)csrRowPtrB);
			result |= cuMemFree((CUdeviceptr)csrColIndB);
			result |= cuMemFree((CUdeviceptr)csrValB);
			//result |= cuMemFree((CUdeviceptr)csrRowPtrC); // reused as csrRowPtrC
			result |= cuMemFree((CUdeviceptr)csrColIndC);
			result |= cuMemFree((CUdeviceptr)csrValC);
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("failed to free matrix storage on GPU");
		}
		// free the operands

		//printf("debug: operands deleted\n"); // debug

#ifdef __SCHUR_PROFILING
		f_mul0_time = dt.f_Time();
#endif // __SCHUR_PROFILING

#ifndef __GPU_SCHUR_EASY_PROD_ONLY

#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
		cs *p_D = cs_spalloc32(n1, m1, nnzD, 1, 0);
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
		cs *p_D = cs_spalloc(n1, m1, nnzD, 1, 0);
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
		{
			int result = cuMemcpyDtoHAsync(p_D->p, (CUdeviceptr)csrRowPtrD, sizeof(int) * (m1 + 1), 0);
			result |= cuMemcpyDtoHAsync(p_D->i, (CUdeviceptr)csrColIndD, sizeof(int) * nnzD, 0);
			result |= cuMemcpyDtoHAsync(p_D->x, (CUdeviceptr)csrValD, sizeof(double) * nnzD, 0);
			//result |= cuMemFree((CUdeviceptr)csrRowPtrD); // don't, it is cached
			result |= cuMemFree((CUdeviceptr)csrColIndD);
			result |= cuMemFree((CUdeviceptr)csrValD);
			if(result != CUDA_SUCCESS)
				throw std::runtime_error("failed to download GEMM result from the GPU");
		}
		// copy, free the product

		//printf("debug: product deleted\n"); // debug

		A.CopyLayoutTo(schur_compl);
		// prepare block layout

		if(cuCtxSynchronize() != CUDA_SUCCESS)
			throw std::runtime_error("cuCtxSynchronize() failed " FILE_LINE);
		// wait for the copy to finish

		//printf("debug: calling From_Sparse32_Parallel()\n"); // debug

		std::vector<size_t> workspace;
#if defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) || defined(__x86_64) || defined(__amd64) || defined(__ia64)
		if(!schur_compl.From_Sparse32_Parallel(0, 0, p_D, false, workspace)) {
#else // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
		if(!schur_compl.From_Sparse_Parallel(0, 0, p_D, false, workspace)) {
#endif // _M_X64 || _M_AMD64 || _M_IA64 || __x86_64 || __amd64 || __ia64
			fprintf(stderr, "error: D is " PRIsize " x " PRIsize ", UC--1 expects " PRIsize " x " PRIsize "\n",
				p_D->m, p_D->n, schur_compl.n_Row_Num(), schur_compl.n_Column_Num());
			throw std::runtime_error("From_Sparse_Parallel() failed");
		}
		// make it back into a block matrix

		cs_spfree(p_D);
		// not needed anymore

#endif // !__GPU_SCHUR_EASY_PROD_ONLY

		/*printf("debug: done\n");*/

#ifdef __GPU_SCHUR_VERIFY_RESULT
		{
			double f_correct_norm = U_Cinv_ref.f_Norm();
			minus_U_Cinv.AddTo(U_Cinv_ref);
			printf("GPU UC^-1 precise to: %g (%g rel)\n", U_Cinv_ref.f_Norm(),
				U_Cinv_ref.f_Norm() / f_correct_norm);
		}
#ifndef __GPU_SCHUR_EASY_PROD_ONLY
		{
			double f_correct_norm = minus_schur_ref.f_Norm();
			schur_compl.AddTo(minus_schur_ref);
			printf("GPU UC^-1V precise to: %g (%g rel)\n", minus_schur_ref.f_Norm(),
				minus_schur_ref.f_Norm() / f_correct_norm);
		}
#endif // !__GPU_SCHUR_EASY_PROD_ONLY
		// debug - verify the GPU results

		/*U_Cinv_ref.Scale(-1);
		minus_U_Cinv.Swap(U_Cinv_ref);
#ifndef __GPU_SCHUR_EASY_PROD_ONLY
		minus_schur_ref.Scale(-1);
		schur_compl.Swap(minus_schur_ref);
#endif // !__GPU_SCHUR_EASY_PROD_ONLY*/
		// debug - replace GPU results by CPU
#endif // __GPU_SCHUR_VERIFY_RESULT
	}

#ifdef __GPU_SCHUR_EASY_PROD_ONLY
	{
		CUberBlockMatrix V;
		U.TransposeTo(V); // because lower-triangular of lambda is not calculated
		// if the other half of the product runs on CPU, need V; U gets destroyed in the process
		minus_U_Cinv.MultiplyToWith/*_FBS<_TyLambdaMatrixBlockSizes,
			_TyLambdaMatrixBlockSizes>*/(schur_compl, V, true); // -U*(C^-1)V // UV is symmetric, the whole product should be symmetric, calculate only the upper triangular part
	}
#endif // __GPU_SCHUR_EASY_PROD_ONLY
#endif // 0

#ifdef __SCHUR_PROFILING
	double f_mul1_time = dt.f_Time();
#endif // __SCHUR_PROFILING

	/*minus_U_Cinv.Scale(-1.0);	// -U*(C^-1)

	double f_scale_time = dt.f_Time();*/

	//C_inv.Scale(-1.0);	// -(C^-1)

#ifdef __SCHUR_PROFILING
	double f_scale_time = 0;//dt.f_Time();
#endif // __SCHUR_PROFILING

	//CUberBlockMatrix schur_compl; // not needed afterwards
	//minus_U_Cinv.MultiplyToWith/*_FBS<_TyLambdaMatrixBlockSizes,
	//	_TyLambdaMatrixBlockSizes>*/(schur_compl, V, true); // -U*(C^-1)V // UV is symmetric, the whole product should be symmetric, calculate only the upper triangular part

	//printf("the product of all 3 matrices is " PRIsize " x " PRIsize " and will have " PRIsize " nnz\n",
	//	schur_compl.n_Row_Num(), schur_compl.n_Column_Num(), schur_compl.n_NonZero_Num());
	// debug

#ifdef __SCHUR_PROFILING
	//double f_mul1_time = dt.f_Time();
#endif // __SCHUR_PROFILING

	A.AddTo/*_FBS<_TyLambdaMatrixBlockSizes>*/(schur_compl); // -U*(C^-1)V + A // A is symmetric, if schur_compl is symmetric, the difference also is
	// compute left-hand side A - U(C^-1)V
	// todo - need multiplication with transpose matrices (in this case schur * U^T)

#ifdef __SCHUR_PROFILING
	double f_add_time = dt.f_Time();
#endif // __SCHUR_PROFILING

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

	if(m_double_workspace.capacity() < n_rhs_vector_size) {
		m_double_workspace.clear(); // avoid data copying
		m_double_workspace.reserve(std::max(2 * m_double_workspace.capacity(), n_rhs_vector_size));
	}
	m_double_workspace.resize(n_rhs_vector_size);
	double *p_double_workspace = &m_double_workspace[0];
	// alloc workspace

	lambda_perm.InversePermute_RightHandSide_Vector(p_double_workspace,
		&r_v_eta(0), n_rhs_vector_size, &m_order[0], m_order.size());
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

	minus_U_Cinv.PreMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_x(0),
		n_pose_vector_size, &v_l(0), n_landmark_vector_size); // x - U(C^-1)l
	// compute right-hand side x - U(C^-1)l

#ifdef __SCHUR_PROFILING
	double f_RHS_time = dt.f_Time();
#endif // __SCHUR_PROFILING

	if(!b_keep_ordering)
		_TyLinearSolverWrapper::FinalBlockStructure(m_linear_solver, schur_compl); // the ordering on schur_compl will not change, can calculate it only in the first pass and then reuse
	bool b_result = _TyLinearSolverWrapper::Solve(m_linear_solver, schur_compl, v_x);
	Eigen::VectorXd &v_dx = v_x; // rename, solves inplace
	// solve for dx = A - U(C^-1)V / x

#ifdef __SCHUR_PROFILING
	double f_linsolve0_time = dt.f_Time();
#endif // __SCHUR_PROFILING

	// note that schur_compl only contains pose-sized blocks when guided ordering is used! could optimize for that
	// also note that schur_compl is not completely dense if it is not many times smaller than C

	Eigen::Map<Eigen::VectorXd>(p_double_workspace, n_pose_vector_size) = v_dx; // an unnecessary copy, maybe could work around
	// obtained a first part of the solution

#if 1
	Eigen::Map<Eigen::VectorXd> v_dl(p_double_workspace +
		n_pose_vector_size, n_landmark_vector_size); // calculated inplace
	v_dl.setZero();
	//V.PreMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_dl(0),
	//	n_landmark_vector_size, &v_dx(0), n_pose_vector_size); // V * dx
	U.PostMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_dl(0),
		n_landmark_vector_size, &v_dx(0), n_pose_vector_size); // dx^T * U^T = (V * dx)^T and the vector transposes are ignored
	v_l -= v_dl; // (l - dl)
	v_dl.setZero();
	C_inv.PreMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_dl(0),
		n_landmark_vector_size, &v_l(0), n_landmark_vector_size); // (C^-1)(l - V * dx)
	// solve for dl = (C^-1)(l - V * dx)
	// the second part of the solution is calculated inplace in the dest vector
#else // 1
	Eigen::Map<Eigen::VectorXd> v_dl(p_double_workspace +
		n_pose_vector_size, n_landmark_vector_size); // calculated inplace
	v_l = -v_l; // trades setZero() for sign negation and a smaller sign negation above
	V.PreMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_l(0),
		n_landmark_vector_size, &v_dx(0), n_pose_vector_size); // V * dx
	// l = V * dx - l
	v_dl.setZero();
	C_inv.PreMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_dl(0),
		n_landmark_vector_size, &v_l(0), n_landmark_vector_size); // (C^-1)(l - V * dx)
	// solve for dl = (C^-1)(l - V * dx) = -(C^-1)(V * dx - l)
	// the second part of the solution is calculated inplace in the dest vector
	// t_odo - carry this modification to the original schur as well
	// todo - free U from memory after calculating U(C^-1) in the original schur as well - or do it inplace using a specialized diagonal multiply kernel
#endif // 1

	lambda_perm.Permute_RightHandSide_Vector(&r_v_eta(0), p_double_workspace,
		n_rhs_vector_size, &m_order[0], m_order.size());
	// permute back!

#ifdef __SCHUR_PROFILING
	double f_linsolve1_time = dt.f_Time();
	double f_totel_time = f_reperm_time + f_slice_time + f_transpose_time +
		f_inverse_time + f_mul0_time + f_scale_time + f_mul1_time +
		f_add_time + f_RHS_time + f_linsolve0_time + f_linsolve1_time;

	printf("Schur took %f sec, out of which:\n", f_totel_time);
	printf("   reperm: %f\n", f_reperm_time);
	printf("    slice: %f\n", f_slice_time);
	printf("transpose: %f\n", f_transpose_time);
	printf("  inverse: %f\n", f_inverse_time);
	printf(" multiply: %f, out of which:\n", f_mul0_time + f_scale_time + f_mul1_time);
	printf("\tdiag gemm: %f\n", f_mul0_time);
	printf("\t    scale: %f\n", f_scale_time);
	printf("\t     gemm: %f\n", f_mul1_time);
	printf("      add: %f\n", f_add_time);
	printf(" RHS prep: %f\n", f_RHS_time);
	printf("  cholsol: %f (" PRIsize " x " PRIsize ", " PRIsize " nnz (%.2f %%))\n",
		f_linsolve0_time, schur_compl.n_Row_Num(), schur_compl.n_Column_Num(),
		schur_compl.n_NonZero_Num(), 100 * float(schur_compl.n_NonZero_Num()) /
		(schur_compl.n_Row_Num() * schur_compl.n_Column_Num()));
	printf(" dy solve: %f\n", f_linsolve1_time);
	// debug - do some profiling
#endif // __SCHUR_PROFILING

	/*static size_t n_iter = 0;
	if(!n_iter) {
		std::string s_name;
		{
			char p_s_it_nr[256];
			sprintf(p_s_it_nr, "schur/lambda_perm_" PRIsize "_", n_iter);
			s_name = p_s_it_nr;
			++ n_iter;
		}
		A.Save_MatrixMarket((s_name + "00.mtx").c_str(), (s_name + "00.bla").c_str());
		U.Save_MatrixMarket((s_name + "01.mtx").c_str(), (s_name + "01.bla").c_str());
		V.Save_MatrixMarket((s_name + "10.mtx").c_str(), (s_name + "10.bla").c_str());
		C.Save_MatrixMarket((s_name + "11.mtx").c_str(), (s_name + "11.bla").c_str());
		C_inv.Save_MatrixMarket((s_name + "11_inv.mtx").c_str(), (s_name + "11_inv.bla").c_str());
		schur_compl.Save_MatrixMarket((s_name + "schur.mtx").c_str(), (s_name + "schur.bla").c_str());
	}*/
	// debug - dump the matrices

	return b_result;
}

/*
 *								=== ~CLinearSolver_Schur_GPUBase ===
 */

#else // /*(_WIN32 || _WIN64) &&*/ !__DISABLE_GPU

// GPU disabled, fallback to reference CPU implementation

/*
 *								=== CLinearSolver_DenseGPU ===
 */

bool CLinearSolver_DenseGPU::b_cula_initialized = false;
size_t CLinearSolver_DenseGPU::n_instance_num = 0;

CLinearSolver_DenseGPU::CLinearSolver_DenseGPU() // throw(std::runtime_error)
{
	if(!n_instance_num) {
		fprintf(stderr, "warning: built without GPU support: fallback to software\n");
		//fprintf(stderr, "warning: this software implementation is slower than SLAM++ default\n"); // not true, dense Eigen would be used anyway and the non-fbs path won't be called
		++ n_instance_num;
	}
}

CLinearSolver_DenseGPU::~CLinearSolver_DenseGPU()
{}

void CLinearSolver_DenseGPU::Free_Memory()
{
	{
		Eigen::MatrixXd empty;
		std::swap(empty, m_t_lambda);
	}
}

bool CLinearSolver_DenseGPU::Solve_PosDef(const CUberBlockMatrix &r_lambda, Eigen::VectorXd &r_eta) // throw(std::bad_alloc)
{
	r_lambda.Convert_to_Dense(m_t_lambda);

	typedef Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> _TyDecomposition; /**< @brief matrix factorization type (Cholesky) */
	_TyDecomposition m_t_factorization;

	m_t_factorization.compute(m_t_lambda);
	if(m_t_factorization.info() != Eigen::Success)
		return false;
	// Cholesky

	r_eta = m_t_factorization.solve(r_eta);
	// solve

	return true;
}

/*
 *								=== ~CLinearSolver_DenseGPU ===
 */

/*
 *								=== CLinearSolver_Schur_GPUBase ===
 */

bool CLinearSolver_Schur_GPUBase::GPUSolve(const CUberBlockMatrix &r_lambda,
	Eigen::VectorXd &r_v_eta, bool b_keep_ordering, size_t n_landmark_dimension,
	std::vector<double> &m_double_workspace, const std::vector<size_t> &m_order,
	const size_t m_n_matrix_cut, _TyBaseSolver &m_linear_solver)
{
#if 1 // don't repeat code (though the code below works just fine)
	throw std::runtime_error("no GPU");
	return false;
#else // 1
	_ASSERTE(r_lambda.b_SymmetricLayout()); // pos-def is supposed to be symmetric
	_ASSERTE(r_v_eta.rows() == r_lambda.n_Row_Num()); // make sure the vector has correct dimension

	_ASSERTE(r_lambda.b_SymmetricLayout());
	_ASSERTE(r_lambda.n_BlockColumn_Num() == m_order.size());
	_ASSERTE((m_order.empty() && !m_n_matrix_cut) || (!m_order.empty() &&
		m_n_matrix_cut > 0 && m_n_matrix_cut < SIZE_MAX && m_n_matrix_cut + 1 < m_order.size()));
	_ASSERTE(r_v_eta.rows() == r_lambda.n_Column_Num());

	CDeltaTimer dt;

	CUberBlockMatrix lambda_perm;
	r_lambda.Permute_UpperTriangular_To(lambda_perm, &m_order[0], m_order.size(), true);
	// reorder the matrix

	double f_reperm_time = dt.f_Time();

	const size_t n = lambda_perm.n_BlockColumn_Num();

	CUberBlockMatrix A, U, C, V;
	const size_t n_matrix_cut = m_n_matrix_cut; // antialiass
	lambda_perm.SliceTo(A, 0, n_matrix_cut, 0, n_matrix_cut, true);
	lambda_perm.SliceTo(U, 0, n_matrix_cut, n_matrix_cut, n, true);
	lambda_perm.SliceTo(C, n_matrix_cut, n, n_matrix_cut, n, true);
	// cut Lambda matrix into pieces
	// \lambda = | A U |
	//           | V C |

	const size_t n_rhs_vector_size = r_lambda.n_Column_Num();
	const size_t n_pose_vector_size = A.n_Column_Num(); // 6 * m_n_matrix_cut;
	const size_t n_landmark_vector_size = U.n_Column_Num(); // 3 * (n - m_n_matrix_cut);
	// not block columns! element ones

	double f_slice_time = dt.f_Time();

	U.TransposeTo(V);

	double f_transpose_time = dt.f_Time();

	_ASSERTE(C.b_BlockDiagonal()); // it is, unless the ordering is bad
	CUberBlockMatrix &C_inv = C; // can do it inplace
	if(n_landmark_dimension == 3) // this can be anticipated in BA
		C_inv.InverseOf_BlockDiag_FBS_Parallel<MakeTypelist_Safe((fbs_ut::CCTSize2D<3, 3>))>(C); // faster version, slightly less general
	else
		C_inv.InverseOf_Symmteric(C); // there is no non-FBS InverseOf_BlockDiag(), just use this one instead 
	// inverse of C

	double f_inverse_time = dt.f_Time();

	C_inv.Scale(-1.0);	// -U*(C^-1)

	double f_scale_time = dt.f_Time();

	CUberBlockMatrix minus_U_Cinv;
	U.MultiplyToWith/*_FBS<_TyLambdaMatrixBlockSizes,
		_TyLambdaMatrixBlockSizes>*/(minus_U_Cinv, C_inv);

	double f_mul0_time = dt.f_Time();

	CUberBlockMatrix schur_compl; // not needed afterwards
	minus_U_Cinv.MultiplyToWith/*_FBS<_TyLambdaMatrixBlockSizes,
		_TyLambdaMatrixBlockSizes>*/(schur_compl, V, true); // -U*(C^-1)V // UV is symmetric, the whole product should be symmetric, calculate only the upper triangular part

	double f_mul1_time = dt.f_Time();

	A.AddTo/*_FBS<_TyLambdaMatrixBlockSizes>*/(schur_compl); // -U*(C^-1)V + A // A is symmetric, if schur_compl is symmetric, the difference also is
	// compute left-hand side A - U(C^-1)V

	double f_add_time = dt.f_Time();

	// note that the sum and difference of two symmetric matrices is again symmetric,
	// but this is not always true for the product

	if(m_double_workspace.capacity() < n_rhs_vector_size) {
		m_double_workspace.clear(); // avoid data copying
		m_double_workspace.reserve(std::max(2 * m_double_workspace.capacity(), n_rhs_vector_size));
	}
	m_double_workspace.resize(n_rhs_vector_size);
	double *p_double_workspace = &m_double_workspace[0];
	// alloc workspace

	lambda_perm.InversePermute_RightHandSide_Vector(p_double_workspace,
		&r_v_eta(0), n_rhs_vector_size, &m_order[0], m_order.size());
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

	minus_U_Cinv.PreMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_x(0),
		n_pose_vector_size, &v_l(0), n_landmark_vector_size); // x - U(C^-1)l
	// compute right-hand side x - U(C^-1)l

	double f_RHS_time = dt.f_Time();

	if(!b_keep_ordering)
		_TyLinearSolverWrapper::FinalBlockStructure(m_linear_solver, schur_compl); // the ordering on schur_compl will not change, can calculate it only in the first pass and then reuse
	bool b_result = _TyLinearSolverWrapper::Solve(m_linear_solver, schur_compl, v_x);
	Eigen::VectorXd &v_dx = v_x; // rename, solves inplace
	// solve for dx = A - U(C^-1)V / x

	double f_linsolve0_time = dt.f_Time();

	// note that schur_compl only contains pose-sized blocks when guided ordering is used! could optimize for that
	// also note that schur_compl is not completely dense if it is not many times smaller than C

	Eigen::Map<Eigen::VectorXd>(p_double_workspace, n_pose_vector_size) = v_dx; // an unnecessary copy, maybe could work around
	// obtained a first part of the solution

	Eigen::Map<Eigen::VectorXd> v_dl(p_double_workspace +
		n_pose_vector_size, n_landmark_vector_size); // calculated inplace
	v_l = -v_l; // trades setZero() for sign negation and a smaller sign negation above
	Eigen::VectorXd v_l_copy = v_l;

	//V.PreMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_l(0),
	//	n_landmark_vector_size, &v_dx(0), n_pose_vector_size); // V * dx
	U.PostMultiply_Add_Parallel/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_l_copy(0),
		n_landmark_vector_size, &v_dx(0), n_pose_vector_size); // dx^T * U^T = (V * dx)^T and the vector transposes are ignored
	//printf("post-mad error %g (rel %g)\n", (v_l_copy - v_l).norm(), (v_l_copy - v_l).norm() / v_l.norm());

	// l = V * dx - l
	v_dl.setZero();
	C_inv.PreMultiply_Add/*_FBS<_TyLambdaMatrixBlockSizes>*/(&v_dl(0),
		n_landmark_vector_size, &v_l(0), n_landmark_vector_size); // (C^-1)(l - V * dx)
	// solve for dl = (C^-1)(l - V * dx) = -(C^-1)(V * dx - l)
	// the second part of the solution is calculated inplace in the dest vector
	// t_odo - carry this modification to the original schur as well
	// todo - free U from memory after calculating U(C^-1) in the original schur as well - or do it inplace using a specialized diagonal multiply kernel

	lambda_perm.Permute_RightHandSide_Vector(&r_v_eta(0), p_double_workspace,
		n_rhs_vector_size, &m_order[0], m_order.size());
	// permute back!

#if defined(_DEBUG) || defined(__SCHUR_PROFILING)
	double f_linsolve1_time = dt.f_Time();
	double f_totel_time = f_reperm_time + f_slice_time + f_transpose_time +
		f_inverse_time + f_mul0_time + f_scale_time + f_mul1_time +
		f_add_time + f_RHS_time + f_linsolve0_time + f_linsolve1_time;

	printf("Schur took %f sec, out of which:\n", f_totel_time);
	printf("   reperm: %f\n", f_reperm_time);
	printf("    slice: %f\n", f_slice_time);
	printf("transpose: %f\n", f_transpose_time);
	printf("  inverse: %f\n", f_inverse_time);
	printf(" multiply: %f, out of which:\n", f_mul0_time + f_scale_time + f_mul1_time);
	printf("\tdiag gemm: %f\n", f_mul0_time);
	printf("\t    scale: %f\n", f_scale_time);
	printf("\t     gemm: %f\n", f_mul1_time);
	printf("      add: %f\n", f_add_time);
	printf(" RHS prep: %f\n", f_RHS_time);
	printf("  cholsol: %f (" PRIsize " x " PRIsize ", " PRIsize " nnz (%.2f %%))\n",
		f_linsolve0_time, schur_compl.n_Row_Num(), schur_compl.n_Column_Num(),
		schur_compl.n_NonZero_Num(), 100 * float(schur_compl.n_NonZero_Num()) /
		(schur_compl.n_Row_Num() * schur_compl.n_Column_Num()));
	printf(" dy solve: %f\n", f_linsolve1_time);
	// debug - do some profiling
#endif // _DEBUG || __SCHUR_PROFILING

	return b_result;
#endif // 1
}

/*
 *								=== ~CLinearSolver_Schur_GPUBase ===
 */

#endif // !__DISABLE_GPU
