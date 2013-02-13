/*
								+---------------------------------+
								|                                 |
								| *** Block matrix unit tests *** |
								|                                 |
								|  Copyright  © -tHE SWINe- 2012  |
								|                                 |
								|          BlockUnit.cpp          |
								|                                 |
								+---------------------------------+
*/

/**
 *	@file src/slam/BlockUnit.cpp
 *	@date 2012
 *	@author -tHE SWINe-
 *	@brief block matrix unit tests
 */

#include <stdio.h>
#include "slam/BlockUnit.h"
#include "slam/Parser.h" // uses CParserBase::ReadLine() and CParserBase::TrimSpace()

/*
 *								=== CBlockMatrixUnitTests ===
 */

bool CBlockMatrixUnitTests::RunAll()
{
	try {
		MatrixMultiplication_UnitTest();
		MatrixAddition_UnitTest();
	} catch(std::runtime_error &r_exc) {
		fprintf(stderr, "error: block matrix unit tests failed: \'%s\'\n", r_exc.what());
		return false;
	} catch(std::bad_alloc&) {
		fprintf(stderr, "error: caught std::bad_alloc while running matrix unit tests\n");
		return false;
	}
	return true;
}

bool CBlockMatrixUnitTests::ReadLine(FILE *p_fr, std::string &r_s_line, int &r_n_cur_line)
{
	while(!feof(p_fr)) {
		if(!CParserBase::ReadLine(r_s_line, p_fr))
			return false;
		// read line

		++ r_n_cur_line;
		// line counter for file debugging

		if(r_s_line.find('%') != std::string::npos)
			r_s_line.erase(r_s_line.find('%'));
		// throw away line comment

		CParserBase::TrimSpace(r_s_line);
		// throw away begin / end whitespace

		if(r_s_line.empty())
			continue;
		// skip empty lines

		return true;
	}

	return false;
}

/**
 *	@brief unit test for block matrix multiplication functions (A * B, A^T * A)
 *	@note This also tests the _FBS A^T * A function.
 *	@todo Add rectangular block tests.
 */
void CBlockMatrixUnitTests::MatrixMultiplication_UnitTest() // throws(std::bad_alloc, std::runtime_error)
{
	_TyBlockMatrix t_block_A, t_block_B, t_block_C;
	for(int i = 0, k = 0; i < block_Size; ++ i) {
		for(int j = 0; j < block_Size; ++ j, ++ k) {
			t_block_A(i, j) = k;
			t_block_B(i, j) = (i == j || block_Size - 1 - i == j)? 1 : 0;
			t_block_C(i, j) = ((i + j) % 3 == 0)? 3 : 0;
		}
	}
	// get some matrix blocks

	printf("test mul (rudimentary test cases) ...\n");

	{
		CUberBlockMatrix A, B, C, R;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_B, 0, 0);
		b_add &= C.Append_Block(_TyBlockMatrix(t_block_A * t_block_B), 0, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.ProductOf(A, B)) // R = A * B
			throw std::runtime_error("ProductOf() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		R.Rasterize("3_AxB.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.b_Equal(C))
			throw std::runtime_error("MatrixMultiplication_UnitTest() failed");
	}
	// a simple case where the matrices only contain a single block

	{
		CUberBlockMatrix A, B, C, R;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_B, 0, 0);
		b_add &= A.Append_Block(t_block_C, block_Size, block_Size);
		b_add &= B.Append_Block(t_block_C, block_Size, block_Size);
		b_add &= C.Append_Block(_TyBlockMatrix(t_block_A * t_block_B), 0, 0);
		b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * t_block_C), block_Size, block_Size);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.ProductOf(A, B)) // R = A * B
			throw std::runtime_error("ProductOf() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		R.Rasterize("3_AxB.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.b_Equal(C))
			throw std::runtime_error("MatrixMultiplication_UnitTest() failed");
	}
	// a simple case where the matrices only contain blocks on diagonal

	{
		CUberBlockMatrix A, B, R;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= A.Append_Block(t_block_B, block_Size, block_Size);
		b_add &= B.Append_Block(t_block_B, 0, 0);
		b_add &= B.Append_Block(t_block_C, 0, block_Size);
		b_add &= B.Append_Block(t_block_C, block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

		cs *p_A = A.p_Convert_to_Sparse(), *p_B = B.p_Convert_to_Sparse(), *p_C = cs_multiply(p_A, p_B);
		cs_spfree(p_A); cs_spfree(p_B);
		// calculate reference using CSparse

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		cs *p_C_dup = cs_add(p_C, p_C, 0, 1);
		for(int i = 0; i < p_C_dup->nzmax; ++ i)
			p_C_dup->x[i] += 1e-3f; // to be nonzero
		CDebug::Dump_SparseMatrix("2_C.tga", p_C, p_C_dup);
		cs_spfree(p_C_dup);
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.ProductOf(A, B)) // R = A * B
			throw std::runtime_error("ProductOf() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		R.Rasterize("3_AxB.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		cs *p_R = R.p_Convert_to_Sparse(), *p_diff = cs_add(p_C, p_R, 1, -1);
		csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
		cs_spfree(p_R); cs_spfree(p_C); cs_spfree(p_diff);
		if(n_nonzero_difference_num)
			throw std::runtime_error("MatrixMultiplication_UnitTest() failed");
		// check results using CSparse
	}
	// a simple case where the matrices contain some blocks, are square and there are no dummy cols / rows

	{
		CUberBlockMatrix A, B, R;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= A.Append_Block(t_block_C, block_Size, 2 * block_Size);
		b_add &= A.Append_Block(t_block_B, 2 * block_Size, 3 * block_Size);
		b_add &= B.Append_Block(t_block_B, 0, 0);
		b_add &= B.Append_Block(t_block_C, 0, 2 * block_Size);
		b_add &= B.Append_Block(t_block_C, 3 * block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

		cs *p_A = A.p_Convert_to_Sparse(), *p_B = B.p_Convert_to_Sparse(), *p_C = cs_multiply(p_A, p_B);
		cs_spfree(p_A); cs_spfree(p_B);
		// calculate reference using CSparse

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		cs *p_C_dup = cs_add(p_C, p_C, 0, 1);
		for(int i = 0; i < p_C_dup->nzmax; ++ i)
			p_C_dup->x[i] += 1e-3f; // to be nonzero
		CDebug::Dump_SparseMatrix("2_C.tga", p_C, p_C_dup);
		cs_spfree(p_C_dup);
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.ProductOf(A, B)) // R = A * B
			throw std::runtime_error("ProductOf() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		R.Rasterize("3_AxB.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		cs *p_R = R.p_Convert_to_Sparse(), *p_diff = cs_add(p_C, p_R, 1, -1);
		csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
		cs_spfree(p_R); cs_spfree(p_C); cs_spfree(p_diff);
		if(n_nonzero_difference_num)
			throw std::runtime_error("MatrixMultiplication_UnitTest() failed");
		// check results using CSparse
	}
	// a case where the matrices are not square and there are some dummy blocks

	{
		CUberBlockMatrix A, B, R;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_C, 2 * block_Size, block_Size);
		b_add &= A.Append_Block(t_block_B, 2 * block_Size, 3 * block_Size);
		b_add &= B.Append_Block(t_block_B, 0, 0);
		b_add &= B.Append_Block(t_block_C, 0, 2 * block_Size);
		b_add &= B.Append_Block(t_block_C, 3 * block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

		cs *p_A = A.p_Convert_to_Sparse(), *p_B = B.p_Convert_to_Sparse(), *p_C = cs_multiply(p_A, p_B);
		cs_spfree(p_A); cs_spfree(p_B);
		// calculate reference using CSparse

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		cs *p_C_dup = cs_add(p_C, p_C, 0, 1);
		for(int i = 0; i < p_C_dup->nzmax; ++ i)
			p_C_dup->x[i] += 1e-3f; // to be nonzero
		CDebug::Dump_SparseMatrix("2_C.tga", p_C, p_C_dup);
		cs_spfree(p_C_dup);
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.ProductOf(A, B)) // R = A * B
			throw std::runtime_error("ProductOf() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		R.Rasterize("3_AxB.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		cs *p_R = R.p_Convert_to_Sparse(), *p_diff = cs_add(p_C, p_R, 1, -1);
		csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
		cs_spfree(p_R); cs_spfree(p_C); cs_spfree(p_diff);
		if(n_nonzero_difference_num)
			throw std::runtime_error("MatrixMultiplication_UnitTest() failed");
		// check results using CSparse
	}
	// a case where the matrices are not square and there are some dummy blocks

	{
		CUberBlockMatrix A, B, R;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= A.Append_Block(t_block_C, block_Size, 2 * block_Size);
		b_add &= A.Append_Block(t_block_B, 2 * block_Size, 3 * block_Size);
		b_add &= A.Append_Block(t_block_B, 5 * block_Size, 3 * block_Size);
		b_add &= B.Append_Block(t_block_B, 0, 0);
		b_add &= B.Append_Block(t_block_C, 0, 2 * block_Size);
		b_add &= B.Append_Block(t_block_C, 3 * block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

		cs *p_A = A.p_Convert_to_Sparse(), *p_B = B.p_Convert_to_Sparse(), *p_C = cs_multiply(p_A, p_B);
		cs_spfree(p_A); cs_spfree(p_B);
		// calculate reference using CSparse

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		cs *p_C_dup = cs_add(p_C, p_C, 0, 1);
		for(int i = 0; i < p_C_dup->nzmax; ++ i)
			p_C_dup->x[i] += 1e-3f; // to be nonzero
		CDebug::Dump_SparseMatrix("2_C.tga", p_C, p_C_dup);
		cs_spfree(p_C_dup);
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.ProductOf(A, B)) // R = A * B
			throw std::runtime_error("ProductOf() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		R.Rasterize("3_AxB.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		cs *p_R = R.p_Convert_to_Sparse(), *p_diff = cs_add(p_C, p_R, 1, -1);
		csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
		cs_spfree(p_R); cs_spfree(p_C); cs_spfree(p_diff);
		if(n_nonzero_difference_num)
			throw std::runtime_error("MatrixMultiplication_UnitTest() failed");
		// check results using CSparse
	}
	// a case where the matrices are not square and there are some dummy blocks, the result is not square

	{
		CUberBlockMatrix A, B, R;

		bool b_add = true;
		b_add &= B.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_C, 2 * block_Size, block_Size);
		b_add &= B.Append_Block(t_block_B, 3 * block_Size, 2 * block_Size);
		b_add &= B.Append_Block(t_block_B, 3 * block_Size, 5 * block_Size);
		b_add &= A.Append_Block(t_block_B, 0, 0);
		b_add &= A.Append_Block(t_block_C, 2 * block_Size, 0);
		b_add &= A.Append_Block(t_block_C, 0, 3 * block_Size);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

		cs *p_A = A.p_Convert_to_Sparse(), *p_B = B.p_Convert_to_Sparse(), *p_C = cs_multiply(p_A, p_B);
		cs_spfree(p_A); cs_spfree(p_B);
		// calculate reference using CSparse

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		cs *p_C_dup = cs_add(p_C, p_C, 0, 1);
		for(int i = 0; i < p_C_dup->nzmax; ++ i)
			p_C_dup->x[i] += 1e-3f; // to be nonzero
		CDebug::Dump_SparseMatrix("2_C.tga", p_C, p_C_dup);
		cs_spfree(p_C_dup);
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!R.ProductOf(A, B)) // R = A * B
			throw std::runtime_error("ProductOf() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
		R.Rasterize("3_AxB.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

		cs *p_R = R.p_Convert_to_Sparse(), *p_diff = cs_add(p_C, p_R, 1, -1);
		csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
		cs_spfree(p_R); cs_spfree(p_C); cs_spfree(p_diff);
		if(n_nonzero_difference_num)
			throw std::runtime_error("MatrixMultiplication_UnitTest() failed");
		// check results using CSparse
	}
	// a case where the matrices are not square and there are some dummy blocks, the result is not square

	for(int n_test = 0; n_test < 5; ++ n_test) {
		const struct {
			size_t n_A_rows, n_A_cols_B_rows, n_B_cols;
			const char *p_s_test_desc;
		} p_test_cfg[5] = {
			{10, 10, 10, "10x10 matrices, result is 10x10"}, // 10x10 matrices, result is 10x10
			{10, 7, 10, "10x7, 7x10 matrices, result is 10x10"}, // 10x7, 7x10 matrices, result is 10x10
			{10, 13, 10, "10x13, 13x10 matrices, result is 10x10"}, // 10x13, 13x10 matrices, result is 10x10
			{10, 13, 7, "10x13, 13x7 matrices, result is 10x7"}, // 10x13, 13x7 matrices, result is 10x7
			{7, 13, 10, "7x13, 13x10 matrices, result is 7x10"} // 7x13, 13x10 matrices, result is 7x10
		};
		// different test configs

		printf("testing mul and AtA on random %s ...\n", p_test_cfg[n_test].p_s_test_desc);
		// verbose

		for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
			const size_t n_block_rows_A = p_test_cfg[n_test].n_A_rows;
			const size_t n_block_cols_A = p_test_cfg[n_test].n_A_cols_B_rows;
			const size_t n_block_rows_B = p_test_cfg[n_test].n_A_cols_B_rows;
			const size_t n_block_cols_B = p_test_cfg[n_test].n_B_cols;
			const size_t n_max_blocks = n_block_rows_A * n_block_cols_B;
			// size of matrices in blocks

			CUberBlockMatrix A, B, R, S, Q, U, V, W, X, J, K, L, M;

			bool b_add = true;
			if(rand() & 1) {
				b_add &= A.Append_Block(t_block_A, 0, 0);
				b_add &= A.Append_Block(t_block_C, block_Size * (n_block_rows_A - 1), block_Size * (n_block_cols_A - 1));
				b_add &= B.Append_Block(t_block_B, 0, block_Size * (n_block_cols_B - 1));
				b_add &= B.Append_Block(t_block_C, block_Size * (n_block_rows_B - 1), 0);
			} else {
				b_add &= B.Append_Block(t_block_A, 0, 0);
				b_add &= B.Append_Block(t_block_C, block_Size * (n_block_rows_B - 1), block_Size * (n_block_cols_B - 1));
				b_add &= A.Append_Block(t_block_B, 0, block_Size * (n_block_cols_A - 1));
				b_add &= A.Append_Block(t_block_C, block_Size * (n_block_rows_A - 1), 0);
			}
			// append blocks to make the matrices of equal size

			for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
				int r = rand() % n_block_rows_A;
				int c = rand() % n_block_cols_A;
				const _TyBlockMatrix *p_block;
				{
					int b = rand() % 3;
					p_block = (b == 0)? &t_block_A : (b == 1)? &t_block_B : &t_block_C;
				}
				// pick row / col / block at random

				A.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
				// add block to the A matrix
			}
			for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
				int r = rand() % n_block_rows_B;
				int c = rand() % n_block_cols_B;
				const _TyBlockMatrix *p_block;
				{
					int b = rand() % 3;
					p_block = (b == 0)? &t_block_A : (b == 1)? &t_block_B : &t_block_C;
				}
				// pick row / col / block at random

				B.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
				// add block to the B matrix
			}
			if(!b_add)
				throw std::runtime_error("Append_Block() failed");
			// add random blocks to A and to B

			cs *p_A = A.p_Convert_to_Sparse(), *p_B = B.p_Convert_to_Sparse(), *p_C = cs_multiply(p_A, p_B);
			cs *p_At = cs_transpose(p_A, 1), *p_Bt = cs_transpose(p_B, 1);
			cs *p_AtA = cs_multiply(p_At, p_A), *p_BtB = cs_multiply(p_Bt, p_B);
			cs_spfree(p_A); cs_spfree(p_B); cs_spfree(p_At); cs_spfree(p_Bt);
			// calculate reference using CSparse

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("0_A.tga");
			B.Rasterize("1_B.tga");
			cs *p_C_dup = cs_add(p_C, p_C, 0, 1);
			for(int i = 0; i < p_C_dup->nzmax; ++ i)
				p_C_dup->x[i] = fabs(p_C_dup->x[i]) + 1e-3f; // to be nonzero
			CDebug::Dump_SparseMatrix("2_C.tga", p_C, p_C_dup);
			cs_spfree(p_C_dup);
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!R.ProductOf(A, B)) // R = A * B
				throw std::runtime_error("ProductOf() failed");
			A.PreMultiplyWithSelfTransposeTo(S); // S = A^T * A
			B.PreMultiplyWithSelfTransposeTo(Q); // Q = B^T * B

			typedef MakeTypelist_Safe((_TyBlockMatrix)) one_type_list;
			typedef MakeTypelist_Safe((_TyBlockMatrix,
				Eigen::Matrix<double, block_Size + 1, block_Size + 1>,
				Eigen::Matrix<double, block_Size, block_Size + 1>,
				Eigen::Matrix<double, block_Size + 1, block_Size>)) two_type_list; // that's correct, instances
			// block size lists with one and two (four actually, two per dimension) different block sizes

			A.PreMultiplyWithSelfTransposeTo_FBS<one_type_list>(U); // U = A^T * A
			B.PreMultiplyWithSelfTransposeTo_FBS<one_type_list>(V); // V = B^T * B
			A.PreMultiplyWithSelfTransposeTo_FBS<two_type_list>(W); // W = A^T * A
			B.PreMultiplyWithSelfTransposeTo_FBS<two_type_list>(X); // X = B^T * B

			CUberBlockMatrix At, Bt;
			A.TransposeTo(At);
			B.TransposeTo(Bt);
			if(!J.ProductOf_FBS<one_type_list, one_type_list>(A, B) ||
			   !K.ProductOf_FBS<two_type_list, two_type_list>(A, B) ||
			   !L.ProductOf_FBS<two_type_list, two_type_list>(At, A) ||
			   !M.ProductOf_FBS<two_type_list, two_type_list>(Bt, B))
				throw std::runtime_error("ProductOf_FBS() failed");

#ifdef __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES
			R.Rasterize("3_AxB.tga");
			S.Rasterize("4_AtA.tga");
			T.Rasterize("5_BtB.tga");
#endif // __BLOCK_MUL_UNIT_TEST_DUMP_MATRIX_IMAGES

			//R.t_FindBlock(0, 0, block_Size, block_Size) += t_block_A;
			// cause error deliberately

			/*size_t n_R_nz = R.n_NonZero_Num();
			csi _n_C_nz = p_C->p[p_C->n];
			csi n_C_nz = cs_droptol(p_C, -1);*/ // debug

			{
				cs *p_R = R.p_Convert_to_Sparse(), *p_diff = cs_add(p_R, p_C, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_R); /*cs_spfree(p_C);*/ cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (AB)");
			}
			{
				cs *p_J = J.p_Convert_to_Sparse(), *p_diff = cs_add(p_J, p_C, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_J); /*cs_spfree(p_C);*/ cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (AB)");
			}
			{
				cs *p_K = K.p_Convert_to_Sparse(), *p_diff = cs_add(p_K, p_C, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_K); cs_spfree(p_C); cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (AB)");
			}
			{
				cs *p_L = L.p_Convert_to_Sparse(), *p_diff = cs_add(p_L, p_AtA, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_L); /*cs_spfree(p_AtA);*/ cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (AtA FBS single)");
			}
			{
				cs *p_M = M.p_Convert_to_Sparse(), *p_diff = cs_add(p_M, p_BtB, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_M); /*cs_spfree(p_BtB);*/ cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (BtB FBS single)");
			}
			{
				cs *p_U = U.p_Convert_to_Sparse(), *p_diff = cs_add(p_U, p_AtA, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_U); /*cs_spfree(p_AtA);*/ cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (AtA FBS single)");
			}
			{
				cs *p_V = V.p_Convert_to_Sparse(), *p_diff = cs_add(p_V, p_BtB, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_V); /*cs_spfree(p_BtB);*/ cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (BtB FBS single)");
			}
			{
				cs *p_W = W.p_Convert_to_Sparse(), *p_diff = cs_add(p_W, p_AtA, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_W); /*cs_spfree(p_AtA);*/ cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (AtA FBS two)");
			}
			{
				cs *p_X = X.p_Convert_to_Sparse(), *p_diff = cs_add(p_X, p_BtB, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_X); /*cs_spfree(p_BtB);*/ cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (BtB FBS two)");
			}
			{
				cs *p_S = S.p_Convert_to_Sparse(), *p_diff = cs_add(p_S, p_AtA, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_S); cs_spfree(p_AtA); cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (AtA)");
			}
			{
				cs *p_Q = Q.p_Convert_to_Sparse(), *p_diff = cs_add(p_Q, p_BtB, -1, 1);
				csi n_nonzero_difference_num = cs_droptol(p_diff, 1e-10);
				cs_spfree(p_Q); cs_spfree(p_BtB); cs_spfree(p_diff);
				if(n_nonzero_difference_num)
					throw std::runtime_error("MatrixMultiplication_UnitTest() failed (BtB)");
			}
			// check results using CSparse
		}
	}

	printf("matrix multiplication unit tests passed\n");
}

void CBlockMatrixUnitTests::MatrixAddition_UnitTest() // throws(std::bad_alloc, std::runtime_error)
{
	_TyBlockMatrix t_block_A, t_block_B, t_block_C;
	for(int i = 0, k = 0; i < block_Size; ++ i) {
		for(int j = 0; j < block_Size; ++ j, ++ k) {
			t_block_A(i, j) = k;
			t_block_B(i, j) = (i == j || block_Size - 1 - i == j)? 1 : 0;
			t_block_C(i, j) = ((i + j) % 3 == 0)? 3 : 0;
		}
	}
	// get some matrix blocks

	_TyBlockMatrix2 t_block_D, t_block_E, t_block_F;
	for(int i = 0, k = 0; i < block_Size2; ++ i) {
		for(int j = 0; j < block_Size2; ++ j, ++ k) {
			t_block_D(i, j) = k;
			t_block_E(i, j) = (i == j || block_Size2 - 1 - i == j)? 1 : 0;
			t_block_F(i, j) = ((i + j) % 3 == 0)? 3 : 0;
		}
	}
	// get some more matrix blocks

	const int n_max_block_size = 10; // not size_t, Eigen would get in trouble
	Eigen::MatrixXd t_block_G(n_max_block_size, n_max_block_size),
		t_block_H(n_max_block_size, n_max_block_size),
		t_block_I(n_max_block_size, n_max_block_size);
	for(int i = 0, k = 0; i < n_max_block_size; ++ i) {
		for(int j = 0; j < n_max_block_size; ++ j, ++ k) {
			t_block_G(i, j) = k;
			t_block_H(i, j) = (i == j || block_Size2 - 1 - i == j)? 1 : 0;
			t_block_I(i, j) = ((i + j) % 3 == 0)? 3 : 0;
		}
	}
	// get some more large matrix blocks

#if 0
	printf("test add (rudimentary test cases) ...\n");

	{
		CUberBlockMatrix A, B, C;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_B, 0, 0);
		b_add &= C.Append_Block(_TyBlockMatrix(t_block_A + t_block_B), 0, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where the matrices only contain a single block

	{
		CUberBlockMatrix A, B, C;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_B, 0, 0);
		b_add &= A.Append_Block(t_block_C, block_Size, block_Size);
		b_add &= B.Append_Block(t_block_B, block_Size, block_Size);
		b_add &= A.Append_Block(t_block_C, 0, block_Size);
		b_add &= B.Append_Block(t_block_A, 0, block_Size);
		b_add &= C.Append_Block(_TyBlockMatrix(t_block_A + t_block_B), 0, 0);
		b_add &= C.Append_Block(_TyBlockMatrix(t_block_C + t_block_B), block_Size, block_Size);
		b_add &= C.Append_Block(_TyBlockMatrix(t_block_C + t_block_A), 0, block_Size);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where the matrices contain blocks at the same locations

	{
		CUberBlockMatrix A, B, C;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_B, block_Size, block_Size);
		b_add &= A.Append_Block(t_block_C, 0, block_Size);
		b_add &= A.Append_Block(t_block_C, block_Size, 0);
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(t_block_B, block_Size, block_Size);
		b_add &= C.Append_Block(t_block_C, 0, block_Size);
		b_add &= C.Append_Block(t_block_C, block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where the matrices contain blocks at different locations

	{
		CUberBlockMatrix A, B, C;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= A.Append_Block(t_block_C, block_Size, block_Size);
		b_add &= B.Append_Block(t_block_B, block_Size, block_Size);
		b_add &= A.Append_Block(t_block_C, 0, block_Size);
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(_TyBlockMatrix(t_block_C + t_block_B), block_Size, block_Size);
		b_add &= C.Append_Block(t_block_C, 0, block_Size);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where the matrices contain blocks at the same or different locations

	{
		CUberBlockMatrix A, B, C;

		/*   0369
		A = |a   |0
			|    |block_Size
			|   c|2 * block_Size
		B = |   c|0
			|  b |block_Size
			|c   |2 * block_Size
		C = |a  c|0
			|  b |block_Size
			|c  c|2 * block_Size
		*/

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= A.Append_Block(t_block_C, 2 * block_Size, 3 * block_Size);
		b_add &= B.Append_Block(t_block_B, block_Size, 2 * block_Size);
		b_add &= B.Append_Block(t_block_C, 0, 3 * block_Size);
		b_add &= B.Append_Block(t_block_C, 2 * block_Size, 0);
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(t_block_B, block_Size, 2 * block_Size);
		b_add &= C.Append_Block(t_block_C, 2 * block_Size, 3 * block_Size);
		b_add &= C.Append_Block(t_block_C, 0, 3 * block_Size);
		b_add &= C.Append_Block(t_block_C, 2 * block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where the columns need to get subdivided, rows are the same

	{
		CUberBlockMatrix A, B, C;

		bool b_add = true;
		b_add &= B.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_C, 2 * block_Size, 3 * block_Size);
		b_add &= A.Append_Block(t_block_B, block_Size, 2 * block_Size);
		b_add &= A.Append_Block(t_block_C, 0, 3 * block_Size);
		b_add &= A.Append_Block(t_block_C, 2 * block_Size, 0);
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(t_block_B, block_Size, 2 * block_Size);
		b_add &= C.Append_Block(t_block_C, 2 * block_Size, 3 * block_Size);
		b_add &= C.Append_Block(t_block_C, 0, 3 * block_Size);
		b_add &= C.Append_Block(t_block_C, 2 * block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where the columns need to get subdivided, rows are the same (swapped A and B)

	{
		CUberBlockMatrix A, B, C;

		bool b_add = true;
		b_add &= B.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_C, 3 * block_Size, 2 * block_Size);
		b_add &= A.Append_Block(t_block_B, 2 * block_Size, block_Size);
		b_add &= A.Append_Block(t_block_C, 3 * block_Size, 0);
		b_add &= A.Append_Block(t_block_C, 0, 2 * block_Size);
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(t_block_B, 2 * block_Size, block_Size);
		b_add &= C.Append_Block(t_block_C, 3 * block_Size, 2 * block_Size);
		b_add &= C.Append_Block(t_block_C, 3 * block_Size, 0);
		b_add &= C.Append_Block(t_block_C, 0, 2 * block_Size);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where the columns need to get subdivided, rows are the same (swapped A and B, transposed)

	{
		CUberBlockMatrix A, B, C;

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= A.Append_Block(t_block_C, 3 * block_Size, 2 * block_Size);
		b_add &= B.Append_Block(t_block_B, 2 * block_Size, block_Size);
		b_add &= B.Append_Block(t_block_C, 3 * block_Size, 0);
		b_add &= B.Append_Block(t_block_C, 0, 2 * block_Size);
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(t_block_B, 2 * block_Size, block_Size);
		b_add &= C.Append_Block(t_block_C, 3 * block_Size, 2 * block_Size);
		b_add &= C.Append_Block(t_block_C, 3 * block_Size, 0);
		b_add &= C.Append_Block(t_block_C, 0, 2 * block_Size);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where the columns need to get subdivided, rows are the same (transposed)

	{
		CUberBlockMatrix A, B, C;

		/*   0369
		A = |a   |0
			|    |block_Size
			|    |2 * block_Size
			|   c|3 * block_Size
		B = |   c|0
			|    |block_Size
			|  b |2 * block_Size
			|c   |3 * block_Size
		C = |a  c|0
			|    |block_Size
			|  b |2 * block_Size
			|c  c|3 * block_Size
		*/

		bool b_add = true;
		b_add &= A.Append_Block(t_block_A, 0, 0);
		b_add &= A.Append_Block(t_block_C, 3 * block_Size, 3 * block_Size);
		b_add &= B.Append_Block(t_block_B, 2 * block_Size, 2 * block_Size);
		b_add &= B.Append_Block(t_block_C, 0, 3 * block_Size);
		b_add &= B.Append_Block(t_block_C, 3 * block_Size, 0);
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(t_block_B, 2 * block_Size, 2 * block_Size);
		b_add &= C.Append_Block(t_block_C, 3 * block_Size, 3 * block_Size);
		b_add &= C.Append_Block(t_block_C, 0, 3 * block_Size);
		b_add &= C.Append_Block(t_block_C, 3 * block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where both columns and rows need to get subdivided

	{
		CUberBlockMatrix A, B, C;

		/*   0369
		A = |a   |0
			|    |block_Size
			|    |2 * block_Size
			|   c|3 * block_Size
		B = |   c|0
			|    |block_Size
			|  b |2 * block_Size
			|c   |3 * block_Size
		C = |a  c|0
			|    |block_Size
			|  b |2 * block_Size
			|c  c|3 * block_Size
		*/

		bool b_add = true;
		b_add &= B.Append_Block(t_block_A, 0, 0);
		b_add &= B.Append_Block(t_block_C, 3 * block_Size, 3 * block_Size);
		b_add &= A.Append_Block(t_block_B, 2 * block_Size, 2 * block_Size);
		b_add &= A.Append_Block(t_block_C, 0, 3 * block_Size);
		b_add &= A.Append_Block(t_block_C, 3 * block_Size, 0);
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(t_block_B, 2 * block_Size, 2 * block_Size);
		b_add &= C.Append_Block(t_block_C, 3 * block_Size, 3 * block_Size);
		b_add &= C.Append_Block(t_block_C, 0, 3 * block_Size);
		b_add &= C.Append_Block(t_block_C, 3 * block_Size, 0);
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!B.AddTo(A))
			throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(!A.b_Equal(C))
			throw std::runtime_error("MatrixAddition_UnitTest() failed");
	}
	// a simple case where both columns and rows need to get subdivided (matrices swapped)

	printf("test random add on bigger matrices ...\n");

	for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 10;
		// size of matrices in blocks

		CUberBlockMatrix A, B, C;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_A, 0, 0);
			b_add &= A.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1));
			b_add &= B.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1));
			b_add &= B.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0);
		} else {
			b_add &= B.Append_Block(t_block_A, 0, 0);
			b_add &= B.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1));
			b_add &= A.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1));
			b_add &= A.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0);
		}
		b_add &= C.Append_Block(t_block_A, 0, 0);
		b_add &= C.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1));
		b_add &= C.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1));
		b_add &= C.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0);
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_A : (b == 1)? &t_block_B : &t_block_C;
			}
			// pick row / col / block at random

			A.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			C.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_A : (b == 1)? &t_block_B : &t_block_C;
			}
			// pick row / col / block at random

			B.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			C.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(rand() & 1) {
			if(!B.AddTo(A))
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo(B))
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test simple add

	printf("test add with src factor ...\n");

	for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 10;
		const double f_factor = 1000; // simple factor that makes it impossible to accidentally commutate results
		// size of matrices in blocks

		CUberBlockMatrix A, B, C;

		bool b_add_to_A = (rand() & 1) != 0;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_A, 0, 0);
			b_add &= A.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1));
			b_add &= B.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1));
			b_add &= B.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0);

			if(b_add_to_A) {
				b_add &= C.Append_Block(t_block_A, 0, 0); // A
				b_add &= C.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1)); // A
				b_add &= C.Append_Block(_TyBlockMatrix(t_block_B * f_factor), 0, block_Size * (n_block_cols - 1)); // B
				b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * f_factor), block_Size * (n_block_rows - 1), 0); // B
			} else {
				b_add &= C.Append_Block(_TyBlockMatrix(t_block_A * f_factor), 0, 0); // A
				b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * f_factor), block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1)); // A
				b_add &= C.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1)); // B
				b_add &= C.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0); // B
			}
		} else {
			b_add &= B.Append_Block(t_block_A, 0, 0);
			b_add &= B.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1));
			b_add &= A.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1));
			b_add &= A.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0);

			if(b_add_to_A) {
				b_add &= C.Append_Block(_TyBlockMatrix(t_block_A * f_factor), 0, 0); // B
				b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * f_factor), block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1)); // B
				b_add &= C.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1)); // A
				b_add &= C.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0); // A
			} else {
				b_add &= C.Append_Block(t_block_A, 0, 0); // B
				b_add &= C.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1)); // B
				b_add &= C.Append_Block(_TyBlockMatrix(t_block_B * f_factor), 0, block_Size * (n_block_cols - 1)); // A
				b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * f_factor), block_Size * (n_block_rows - 1), 0); // A
			}
		}
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_A : (b == 1)? &t_block_B : &t_block_C;
			}
			// pick row / col / block at random

			A.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			if(b_add_to_A)
				C.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block; // b_add_to_A => A += B * f_factor (A unmodified)
			else
				C.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block * f_factor;
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_A : (b == 1)? &t_block_B : &t_block_C;
			}
			// pick row / col / block at random

			B.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			if(b_add_to_A)
				C.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block * f_factor; // b_add_to_A => A += B * f_factor (B modified)
			else
				C.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(b_add_to_A) {
			if(!B.AddTo(A, f_factor)) // A += B * f_factor
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo(B, f_factor))
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test add with src factor

	printf("test add with both src and dest factors ...\n");

	for(int n_pass = 0; n_pass < 10000; ++ n_pass) { // this requires more iterations to test border cases
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 100; // this requires more blocks to test border cases
		const double f_factor_A = 19;
		const double f_factor_B = 23; // simple factors that makes it impossible to accidentally commutate results (primes)
		// size of matrices in blocks

		CUberBlockMatrix A, B, C;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_A, 0, 0);
			b_add &= A.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1));
			b_add &= B.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1));
			b_add &= B.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0);

			b_add &= C.Append_Block(_TyBlockMatrix(t_block_A * f_factor_A), 0, 0); // A
			b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * f_factor_A), block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1)); // A
			b_add &= C.Append_Block(_TyBlockMatrix(t_block_B * f_factor_B), 0, block_Size * (n_block_cols - 1)); // B
			b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * f_factor_B), block_Size * (n_block_rows - 1), 0); // B
		} else {
			b_add &= B.Append_Block(t_block_A, 0, 0);
			b_add &= B.Append_Block(t_block_C, block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1));
			b_add &= A.Append_Block(t_block_B, 0, block_Size * (n_block_cols - 1));
			b_add &= A.Append_Block(t_block_C, block_Size * (n_block_rows - 1), 0);

			b_add &= C.Append_Block(_TyBlockMatrix(t_block_A * f_factor_B), 0, 0); // B
			b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * f_factor_B), block_Size * (n_block_rows - 1), block_Size * (n_block_cols - 1)); // B
			b_add &= C.Append_Block(_TyBlockMatrix(t_block_B * f_factor_A), 0, block_Size * (n_block_cols - 1)); // A
			b_add &= C.Append_Block(_TyBlockMatrix(t_block_C * f_factor_A), block_Size * (n_block_rows - 1), 0); // A
		}
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_A : (b == 1)? &t_block_B : &t_block_C;
			}
			// pick row / col / block at random

			A.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			C.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block * f_factor_A;
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_A : (b == 1)? &t_block_B : &t_block_C;
			}
			// pick row / col / block at random

			B.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block;
			C.t_FindBlock(r * block_Size, c * block_Size, block_Size, block_Size) += *p_block * f_factor_B;
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(rand() & 1) {
			if(!B.AddTo(A, f_factor_A, f_factor_B)) // A = A * f_factor_A + B * f_factor_B
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo(B, f_factor_B, f_factor_A)) // B = B * f_factor_B + A * f_factor_A
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test add with both src and dest factor

	printf("test add on bigger random rectangular matrices ...\n");

	for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 10;
		// size of matrices in blocks

		std::vector<size_t> row_sizes, col_sizes;
		std::vector<size_t> row_cumsums, col_cumsums;
		row_cumsums.push_back(0);
		for(size_t n_row = 0; n_row < n_block_rows; ++ n_row) {
			row_sizes.push_back(((rand() & 1)? block_Size : block_Size2));
			row_cumsums.push_back(row_cumsums.back() + row_sizes.back());
		}
		col_cumsums.push_back(0);
		for(size_t n_col = 0; n_col < n_block_cols; ++ n_col) {
			col_sizes.push_back(((rand() & 1)? block_Size : block_Size2));
			col_cumsums.push_back(col_cumsums.back() + col_sizes.back());
		}
		// generate row / column cumsums

		_ASSERTE(block_Size <= block_Size2);
		// make sure that matrix of block_Size2 is larger than of block_Size
		// (will use its upper-left submatrices as block seeds)

		CUberBlockMatrix A, B, C;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= A.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);
		} else {
			b_add &= B.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= B.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);
		}
		b_add &= C.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
		b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
			row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
		b_add &= C.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
			0, col_cumsums[n_block_cols - 1]);
		b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
			row_cumsums[n_block_rows - 1], 0);
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix2 *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_D : (b == 1)? &t_block_E : &t_block_F;
			}
			// pick row / col / block at random

			A.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix2 *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_D : (b == 1)? &t_block_E : &t_block_F;
			}
			// pick row / col / block at random

			B.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

		typedef MakeTypelist_Safe((Eigen::Matrix<double, block_Size, block_Size>,
			Eigen::Matrix<double, block_Size2, block_Size>,
			Eigen::Matrix<double, block_Size, block_Size2>,
			Eigen::Matrix<double, block_Size2, block_Size2>)) two_typelist;

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(rand() & 1) {
			if(!B.AddTo_FBS<two_typelist>(A))
				throw std::runtime_error("AddTo_FBS() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo_FBS<two_typelist>(B))
				throw std::runtime_error("AddTo_FBS() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test simple add with different tile sizes (FBS-able)

	printf("test add on bigger random unruly rectangular matrices ...\n");

	for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 10;
		const int n_max_block_size = 10; // not size_t, Eigen would get in trouble
		// size of matrices in blocks

		Eigen::MatrixXd t_block_G(n_max_block_size, n_max_block_size),
			t_block_H(n_max_block_size, n_max_block_size),
			t_block_I(n_max_block_size, n_max_block_size);
		for(int i = 0, k = 0; i < n_max_block_size; ++ i) {
			for(int j = 0; j < n_max_block_size; ++ j, ++ k) {
				t_block_G(i, j) = k;
				t_block_H(i, j) = (i == j || block_Size2 - 1 - i == j)? 1 : 0;
				t_block_I(i, j) = ((i + j) % 3 == 0)? 3 : 0;
			}
		}
		// get some more matrix blocks

		std::vector<size_t> row_sizes, col_sizes;
		std::vector<size_t> row_cumsums, col_cumsums;
		row_cumsums.push_back(0);
		for(size_t n_row = 0; n_row < n_block_rows; ++ n_row) {
			row_sizes.push_back(rand() % n_max_block_size + 1);
			row_cumsums.push_back(row_cumsums.back() + row_sizes.back());
		}
		col_cumsums.push_back(0);
		for(size_t n_col = 0; n_col < n_block_cols; ++ n_col) {
			col_sizes.push_back(rand() % n_max_block_size + 1);
			col_cumsums.push_back(col_cumsums.back() + col_sizes.back());
		}
		// generate row / column cumsums

		CUberBlockMatrix A, B, C;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= A.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);
		} else {
			b_add &= B.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= B.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);
		}
		b_add &= C.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
		b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
			row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
		b_add &= C.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
			0, col_cumsums[n_block_cols - 1]);
		b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
			row_cumsums[n_block_rows - 1], 0);
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const Eigen::MatrixXd *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_G : (b == 1)? &t_block_H : &t_block_I;
			}
			// pick row / col / block at random

			A.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const Eigen::MatrixXd *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_G : (b == 1)? &t_block_H : &t_block_I;
			}
			// pick row / col / block at random

			B.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(rand() & 1) {
			if(!B.AddTo(A))
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo(B))
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test simple add with different tile sizes (not FBS-able)

	printf("test add with src factor on bigger random rectangular matrices ...\n");

	for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 10;
		const double f_factor = 1000; // simple factor that makes it impossible to accidentally commutate results
		// size of matrices in blocks

		bool b_add_to_A = (rand() & 1) != 0;

		std::vector<size_t> row_sizes, col_sizes;
		std::vector<size_t> row_cumsums, col_cumsums;
		row_cumsums.push_back(0);
		for(size_t n_row = 0; n_row < n_block_rows; ++ n_row) {
			row_sizes.push_back(((rand() & 1)? block_Size : block_Size2));
			row_cumsums.push_back(row_cumsums.back() + row_sizes.back());
		}
		col_cumsums.push_back(0);
		for(size_t n_col = 0; n_col < n_block_cols; ++ n_col) {
			col_sizes.push_back(((rand() & 1)? block_Size : block_Size2));
			col_cumsums.push_back(col_cumsums.back() + col_sizes.back());
		}
		// generate row / column cumsums

		_ASSERTE(block_Size <= block_Size2);
		// make sure that matrix of block_Size2 is larger than of block_Size
		// (will use its upper-left submatrices as block seeds)

		CUberBlockMatrix A, B, C;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= A.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);

			if(b_add_to_A) {
				b_add &= C.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0); // A
				b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
					row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // A
				b_add &= C.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]) * f_factor,
					0, col_cumsums[n_block_cols - 1]); // B
				b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]) * f_factor,
					row_cumsums[n_block_rows - 1], 0); // B
			} else {
				b_add &= C.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]) * f_factor, 0, 0);
				b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]) * f_factor,
					row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
				b_add &= C.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
					0, col_cumsums[n_block_cols - 1]);
				b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
					row_cumsums[n_block_rows - 1], 0);
			}
		} else {
			b_add &= B.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= B.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);

			if(b_add_to_A) {
				b_add &= C.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]) * f_factor, 0, 0); // B
				b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]) * f_factor,
					row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // B
				b_add &= C.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
					0, col_cumsums[n_block_cols - 1]); // A
				b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
					row_cumsums[n_block_rows - 1], 0); // A
			} else {
				b_add &= C.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0); // B
				b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
					row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // B
				b_add &= C.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]) * f_factor,
					0, col_cumsums[n_block_cols - 1]); // A
				b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]) * f_factor,
					row_cumsums[n_block_rows - 1], 0); // A
			}
		}
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix2 *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_D : (b == 1)? &t_block_E : &t_block_F;
			}
			// pick row / col / block at random

			A.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			if(b_add_to_A) {
				C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
					p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			} else {
				C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
					p_block->block(0, 0, row_sizes[r], col_sizes[c]) * f_factor;
			}
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix2 *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_D : (b == 1)? &t_block_E : &t_block_F;
			}
			// pick row / col / block at random

			B.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			if(b_add_to_A) {
				C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
					p_block->block(0, 0, row_sizes[r], col_sizes[c]) * f_factor;
			} else {
				C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
					p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			}
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

		typedef MakeTypelist_Safe((Eigen::Matrix<double, block_Size, block_Size>,
			Eigen::Matrix<double, block_Size2, block_Size>,
			Eigen::Matrix<double, block_Size, block_Size2>,
			Eigen::Matrix<double, block_Size2, block_Size2>)) two_typelist;

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(b_add_to_A) {
			if(!B.AddTo_FBS<two_typelist>(A, f_factor))
				throw std::runtime_error("AddTo_FBS() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo_FBS<two_typelist>(B, f_factor))
				throw std::runtime_error("AddTo_FBS() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test simple add with different tile sizes (FBS-able)

	printf("test add with src factor on bigger random unruly rectangular matrices ...\n");

	for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 10;
		const int n_max_block_size = 10; // not size_t, Eigen would get in trouble
		const double f_factor = 1000; // simple factor that makes it impossible to accidentally commutate results
		// size of matrices in blocks

		bool b_add_to_A = (rand() & 1) != 0;

		Eigen::MatrixXd t_block_G(n_max_block_size, n_max_block_size),
			t_block_H(n_max_block_size, n_max_block_size),
			t_block_I(n_max_block_size, n_max_block_size);
		for(int i = 0, k = 0; i < n_max_block_size; ++ i) {
			for(int j = 0; j < n_max_block_size; ++ j, ++ k) {
				t_block_G(i, j) = k;
				t_block_H(i, j) = (i == j || block_Size2 - 1 - i == j)? 1 : 0;
				t_block_I(i, j) = ((i + j) % 3 == 0)? 3 : 0;
			}
		}
		// get some more matrix blocks

		std::vector<size_t> row_sizes, col_sizes;
		std::vector<size_t> row_cumsums, col_cumsums;
		row_cumsums.push_back(0);
		for(size_t n_row = 0; n_row < n_block_rows; ++ n_row) {
			row_sizes.push_back(rand() % n_max_block_size + 1);
			row_cumsums.push_back(row_cumsums.back() + row_sizes.back());
		}
		col_cumsums.push_back(0);
		for(size_t n_col = 0; n_col < n_block_cols; ++ n_col) {
			col_sizes.push_back(rand() % n_max_block_size + 1);
			col_cumsums.push_back(col_cumsums.back() + col_sizes.back());
		}
		// generate row / column cumsums

		CUberBlockMatrix A, B, C;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= A.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);

			if(b_add_to_A) {
				b_add &= C.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0); // A
				b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
					row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // A
				b_add &= C.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]) * f_factor,
					0, col_cumsums[n_block_cols - 1]); // B
				b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]) * f_factor,
					row_cumsums[n_block_rows - 1], 0); // B
			} else {
				b_add &= C.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]) * f_factor, 0, 0);
				b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]) * f_factor,
					row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
				b_add &= C.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
					0, col_cumsums[n_block_cols - 1]);
				b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
					row_cumsums[n_block_rows - 1], 0);
			}
		} else {
			b_add &= B.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= B.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);

			if(b_add_to_A) {
				b_add &= C.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]) * f_factor, 0, 0); // B
				b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]) * f_factor,
					row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // B
				b_add &= C.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
					0, col_cumsums[n_block_cols - 1]); // A
				b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
					row_cumsums[n_block_rows - 1], 0); // A
			} else {
				b_add &= C.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0); // B
				b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
					row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // B
				b_add &= C.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]) * f_factor,
					0, col_cumsums[n_block_cols - 1]); // A
				b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]) * f_factor,
					row_cumsums[n_block_rows - 1], 0); // A
			}
		}
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const Eigen::MatrixXd *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_G : (b == 1)? &t_block_H : &t_block_I;
			}
			// pick row / col / block at random

			A.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			if(b_add_to_A) {
				C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
					p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			} else {
				C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
					p_block->block(0, 0, row_sizes[r], col_sizes[c]) * f_factor;
			}
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const Eigen::MatrixXd *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_G : (b == 1)? &t_block_H : &t_block_I;
			}
			// pick row / col / block at random

			B.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			if(b_add_to_A) {
				C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
					p_block->block(0, 0, row_sizes[r], col_sizes[c]) * f_factor;
			} else {
				C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
					p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			}
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(b_add_to_A) {
			if(!B.AddTo(A, f_factor))
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo(B, f_factor))
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test simple add with different tile sizes (not FBS-able)

	printf("test add with src and dest factor on bigger random rectangular matrices ...\n");

	for(int n_pass = 0; n_pass < 10000; ++ n_pass) {
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 10;
		const double f_factor_A = 19;
		const double f_factor_B = 23; // simple factors that makes it impossible to accidentally commutate results (primes)
		// size of matrices in blocks

		std::vector<size_t> row_sizes, col_sizes;
		std::vector<size_t> row_cumsums, col_cumsums;
		row_cumsums.push_back(0);
		for(size_t n_row = 0; n_row < n_block_rows; ++ n_row) {
			row_sizes.push_back(((rand() & 1)? block_Size : block_Size2));
			row_cumsums.push_back(row_cumsums.back() + row_sizes.back());
		}
		col_cumsums.push_back(0);
		for(size_t n_col = 0; n_col < n_block_cols; ++ n_col) {
			col_sizes.push_back(((rand() & 1)? block_Size : block_Size2));
			col_cumsums.push_back(col_cumsums.back() + col_sizes.back());
		}
		// generate row / column cumsums

		_ASSERTE(block_Size <= block_Size2);
		// make sure that matrix of block_Size2 is larger than of block_Size
		// (will use its upper-left submatrices as block seeds)

		CUberBlockMatrix A, B, C;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= A.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);

			b_add &= C.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]) * f_factor_A, 0, 0); // A
			b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]) * f_factor_A,
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // A
			b_add &= C.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]) * f_factor_B,
				0, col_cumsums[n_block_cols - 1]); // B
			b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]) * f_factor_B,
				row_cumsums[n_block_rows - 1], 0); // B
		} else {
			b_add &= B.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= B.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);

			b_add &= C.Append_Block(t_block_D.block(0, 0, row_sizes[0], col_sizes[0]) * f_factor_B, 0, 0); // B
			b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]) * f_factor_B,
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // B
			b_add &= C.Append_Block(t_block_E.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]) * f_factor_A,
				0, col_cumsums[n_block_cols - 1]); // A
			b_add &= C.Append_Block(t_block_F.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]) * f_factor_A,
				row_cumsums[n_block_rows - 1], 0); // A
		}
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix2 *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_D : (b == 1)? &t_block_E : &t_block_F;
			}
			// pick row / col / block at random

			A.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]) * f_factor_A;
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const _TyBlockMatrix2 *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_D : (b == 1)? &t_block_E : &t_block_F;
			}
			// pick row / col / block at random

			B.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]) * f_factor_B;
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

		typedef MakeTypelist_Safe((Eigen::Matrix<double, block_Size, block_Size>,
			Eigen::Matrix<double, block_Size2, block_Size>,
			Eigen::Matrix<double, block_Size, block_Size2>,
			Eigen::Matrix<double, block_Size2, block_Size2>)) two_typelist;

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(rand() & 1) {
			if(!B.AddTo(A, f_factor_A, f_factor_B)) // A = A * f_factor_A + B * f_factor_B
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo(B, f_factor_B, f_factor_A)) // B = B * f_factor_B + A * f_factor_A
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test simple add with different tile sizes (FBS-able)

	printf("test add with src and dest factor on bigger random unruly rectangular matrices ...\n");

	for(int n_pass = 0; n_pass < 10000; ++ n_pass) {
		const size_t n_block_cols = 10;
		const size_t n_block_rows = 10;
		const size_t n_max_blocks = 10;
		const int n_max_block_size = 10; // not size_t, Eigen would get in trouble
		const double f_factor_A = 19;
		const double f_factor_B = 23; // simple factors that makes it impossible to accidentally commutate results (primes)
		// size of matrices in blocks

		Eigen::MatrixXd t_block_G(n_max_block_size, n_max_block_size),
			t_block_H(n_max_block_size, n_max_block_size),
			t_block_I(n_max_block_size, n_max_block_size);
		for(int i = 0, k = 0; i < n_max_block_size; ++ i) {
			for(int j = 0; j < n_max_block_size; ++ j, ++ k) {
				t_block_G(i, j) = k;
				t_block_H(i, j) = (i == j || block_Size2 - 1 - i == j)? 1 : 0;
				t_block_I(i, j) = ((i + j) % 3 == 0)? 3 : 0;
			}
		}
		// get some more matrix blocks

		std::vector<size_t> row_sizes, col_sizes;
		std::vector<size_t> row_cumsums, col_cumsums;
		row_cumsums.push_back(0);
		for(size_t n_row = 0; n_row < n_block_rows; ++ n_row) {
			row_sizes.push_back(rand() % n_max_block_size + 1);
			row_cumsums.push_back(row_cumsums.back() + row_sizes.back());
		}
		col_cumsums.push_back(0);
		for(size_t n_col = 0; n_col < n_block_cols; ++ n_col) {
			col_sizes.push_back(rand() % n_max_block_size + 1);
			col_cumsums.push_back(col_cumsums.back() + col_sizes.back());
		}
		// generate row / column cumsums

		CUberBlockMatrix A, B, C;

		bool b_add = true;
		if(rand() & 1) {
			b_add &= A.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= A.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= B.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);

			b_add &= C.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]) * f_factor_A, 0, 0); // A
			b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]) * f_factor_A,
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // A
			b_add &= C.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]) * f_factor_B,
				0, col_cumsums[n_block_cols - 1]); // B
			b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]) * f_factor_B,
				row_cumsums[n_block_rows - 1], 0); // B
		} else {
			b_add &= B.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]), 0, 0);
			b_add &= B.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]),
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]),
				0, col_cumsums[n_block_cols - 1]);
			b_add &= A.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]),
				row_cumsums[n_block_rows - 1], 0);

			b_add &= C.Append_Block(t_block_G.block(0, 0, row_sizes[0], col_sizes[0]) * f_factor_B, 0, 0); // B
			b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[n_block_cols - 1]) * f_factor_B,
				row_cumsums[n_block_rows - 1], col_cumsums[n_block_cols - 1]); // B
			b_add &= C.Append_Block(t_block_H.block(0, 0, row_sizes[0], col_sizes[n_block_cols - 1]) * f_factor_A,
				0, col_cumsums[n_block_cols - 1]); // A
			b_add &= C.Append_Block(t_block_I.block(0, 0, row_sizes[n_block_rows - 1], col_sizes[0]) * f_factor_A,
				row_cumsums[n_block_rows - 1], 0); // A
		}
		// append blocks to make the matrices of equal size

		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const Eigen::MatrixXd *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_G : (b == 1)? &t_block_H : &t_block_I;
			}
			// pick row / col / block at random

			A.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]) * f_factor_A;
			// add block to the A matrix, and to the C matrix
		}
		for(int i = 0, n = rand() % n_max_blocks; i < n; ++ i) {
			int r = rand() % n_block_rows;
			int c = rand() % n_block_cols;
			const Eigen::MatrixXd *p_block;
			{
				int b = rand() % 3;
				p_block = (b == 0)? &t_block_G : (b == 1)? &t_block_H : &t_block_I;
			}
			// pick row / col / block at random

			B.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]);
			C.t_FindBlock(row_cumsums[r], col_cumsums[c], row_sizes[r], col_sizes[c]) +=
				p_block->block(0, 0, row_sizes[r], col_sizes[c]) * f_factor_B;
			// add block to the B matrix, and to the C matrix
		}
		if(!b_add)
			throw std::runtime_error("Append_Block() failed");
		// add random blocks to A, C and to B, C (C now contains result, A and B contains parts)

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		A.Rasterize("0_A.tga");
		B.Rasterize("1_B.tga");
		C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

		if(rand() & 1) {
			if(!B.AddTo(A, f_factor_A, f_factor_B)) // A = A * f_factor_A + B * f_factor_B
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!A.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		} else {
			if(!A.AddTo(B, f_factor_B, f_factor_A)) // B = B * f_factor_B + A * f_factor_A
				throw std::runtime_error("AddTo() failed");

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
			B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

			if(!B.b_Equal(C))
				throw std::runtime_error("MatrixAddition_UnitTest() failed");
		}
		// add A += B or B += A, check result corectness
	}
	// test simple add with different tile sizes (not FBS-able)
#endif // 0

	for(int n_size = 0; n_size < 5; ++ n_size) {
		int p_block_rows[] = {10, 7, 12, 13, 5};
		int p_block_cols[] = {10, 12, 7, 5, 13};
		int p_block_fill[] = {10, 10, 10, 10, 10};
		// different size and density configurations

		for(int n_method = 0; n_method < 6; ++ n_method) {
			printf("test add (%d x %d blocks %s block size rectangular matrices%s%s) ...\n",
				p_block_rows[n_size], p_block_cols[n_size],
				(n_method < 3)? "free" : "fixed", (n_method % 3 == 0)? "" :
				(n_method % 3 == 1)? " with src factor" : " with src and dest factor",
				(n_method < 3)? "" : ", FBS");
			// verbose ...

			for(int n_polarity = 0; n_polarity < 2; ++ n_polarity) {
				double p_factor_A[] = {1.0, (!n_polarity)? 1.0 : 1000.0, 17.0};
				double p_factor_B[] = {1.0, (!n_polarity)? 1000.0 : 1.0, 17.0};
				// factors (the second method can only use one factor, based on polarity) ...

				for(int n_generator = 0; n_generator < ((n_method < 3)? 2 : 1); ++ n_generator) {
					for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
						CUberBlockMatrix A, B, C;

						if(!n_generator) {
							Generate_Addable_Matrices(A, B, p_factor_A[n_method % 3],
								p_factor_B[n_method % 3], C, t_block_D, t_block_E, t_block_F, n_method < 3,
								p_block_rows[n_size], p_block_cols[n_size], p_block_fill[n_size]);
							// stick to block_Size and block_Size2
						} else {
							Generate_Addable_Matrices(A, B, p_factor_A[n_method % 3],
								p_factor_B[n_method % 3], C, t_block_G, t_block_H, t_block_I, true,
								p_block_rows[n_size], p_block_cols[n_size], p_block_fill[n_size]);
							// use really big blocks to cause bigger layout variance (only non-fbs tests)
						}
						// generate some pretty random matrices

						Test_Add(A, B, p_factor_A[n_method % 3],
							p_factor_B[n_method % 3], C, n_polarity != 0, n_method);
						// test
					}
					// test add failures with different tile sizes
				}
			}
		}
	}

	printf("test add failures (rudimentary test cases) ...\n");

	for(int n_generator = 0; n_generator < 2; ++ n_generator) {
		for(int n_polarity = 0; n_polarity < 2; ++ n_polarity) {
			for(int n_method = 0; n_method < 6; ++ n_method) {
				CUberBlockMatrix A, B;

				_ASSERTE(block_Size != block_Size2);

				bool b_add = true;
				if(!n_generator) {
					b_add &= A.Append_Block(t_block_A, 0, 0);
					b_add &= B.Append_Block(t_block_D, 0, 0);
					_ASSERTE(A.n_Column_Num() != B.n_Column_Num() &&
						A.n_Row_Num() != B.n_Row_Num());
					// size of the matrices is different
				} else {
					b_add &= A.Append_Block(t_block_A, block_Size2, block_Size2);
					b_add &= B.Append_Block(t_block_D, block_Size, block_Size);
					_ASSERTE(A.n_Column_Num() == B.n_Column_Num() &&
						A.n_Row_Num() == B.n_Row_Num());
					// size of the matrices is the same
				}
				if(!b_add)
					throw std::runtime_error("Append_Block() failed");
				// generate

				Test_AddFailure(A, B, !n_polarity, n_method);
				// test
			}
		}
	}
	// a simple case where the matrices only contain a single block of different size
	// (size of the matrix is the same in one test and different in the other)

	for(int n_size = 0; n_size < 5; ++ n_size) {
		int p_block_rows[] = {10, 7, 12, 13, 5};
		int p_block_cols[] = {10, 12, 7, 5, 13};
		int p_block_fill[] = {10, 10, 10, 10, 10};
		// different size and density configurations

		for(int n_method = 0; n_method < 6; ++ n_method) {
			printf("test add failures (%d x %d blocks %s block size rectangular matrices%s%s) ...\n",
				p_block_rows[n_size], p_block_cols[n_size],
				(n_method < 3)? "free" : "fixed", (n_method % 3 == 0)? "" :
				(n_method % 3 == 1)? " with src factor" : " with src and dest factor",
				(n_method < 3)? "" : ", FBS");
			// verbose ...

			for(int n_pass = 0; n_pass < 1000; ++ n_pass) {
				CUberBlockMatrix A, B;

				Generate_UnAddable_Matrices(A, B, t_block_D, t_block_E, t_block_F, n_method < 3,
					p_block_rows[n_size], p_block_cols[n_size], p_block_fill[n_size]);
				// generate some pretty random matrices

				Test_AddFailure(A, B, !(rand() & 1), n_method);
				// test
			}
			// test add failures with different tile sizes (FBS-able)
		}
	}

	printf("matrix add unit tests passed\n");
}

void CBlockMatrixUnitTests::Test_Add(CUberBlockMatrix &r_A, CUberBlockMatrix &r_B,
	double f_factor_A, double f_factor_B, const CUberBlockMatrix &UNUSED(r_C),
	bool b_polarity, int n_method) // throws(std::bad_alloc, std::runtime_error)
{
#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
	r_A.Rasterize("0_A.tga");
	r_B.Rasterize("1_B.tga");
	r_C.Rasterize("2_C.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

	typedef MakeTypelist_Safe((Eigen::Matrix<double, block_Size, block_Size>,
		Eigen::Matrix<double, block_Size2, block_Size>,
		Eigen::Matrix<double, block_Size, block_Size2>,
		Eigen::Matrix<double, block_Size2, block_Size2>)) two_typelist;

	if(b_polarity) {
		switch(n_method) {
		case 0:
			if(!r_B.AddTo(r_A))
				throw std::runtime_error("AddTo() failed");
			break;
		case 1:
			_ASSERTE(f_factor_B == 1.0); // note this will compare as it is initialized with the exactly same literal
			if(!r_B.AddTo(r_A, f_factor_A))
				throw std::runtime_error("AddTo() failed");
			break;
		case 2:
			if(!r_B.AddTo(r_A, f_factor_A, f_factor_B))
				throw std::runtime_error("AddTo() failed");
			break;
		case 3:
			if(!r_B.AddTo_FBS<two_typelist>(r_A))
				throw std::runtime_error("AddTo_FBS() failed");
			break;
		case 4:
			_ASSERTE(f_factor_B == 1.0); // note this will compare as it is initialized with the exactly same literal
			if(!r_B.AddTo_FBS<two_typelist>(r_A, f_factor_A))
				throw std::runtime_error("AddTo_FBS() failed");
			break;
		default:
		case 5:
			if(!r_B.AddTo_FBS<two_typelist>(r_A, f_factor_A, f_factor_B))
				throw std::runtime_error("AddTo_FBS() failed");
			break;
		}

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		r_B.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
	} else {
		switch(n_method) {
		case 0:
			if(!r_A.AddTo(r_B))
				throw std::runtime_error("AddTo() failed");
			break;
		case 1:
			_ASSERTE(f_factor_A == 1.0); // note this will compare as it is initialized with the exactly same literal
			if(!r_A.AddTo(r_B, f_factor_B))
				throw std::runtime_error("AddTo() failed");
			break;
		case 2:
			if(!r_A.AddTo(r_B, f_factor_B, f_factor_A))
				throw std::runtime_error("AddTo() failed");
			break;
		case 3:
			if(!r_A.AddTo_FBS<two_typelist>(r_B))
				throw std::runtime_error("AddTo_FBS() failed");
			break;
		case 4:
			_ASSERTE(f_factor_A == 1.0); // note this will compare as it is initialized with the exactly same literal
			if(!r_A.AddTo_FBS<two_typelist>(r_B, f_factor_B))
				throw std::runtime_error("AddTo_FBS() failed");
			break;
		default:
		case 5:
			if(!r_A.AddTo_FBS<two_typelist>(r_B, f_factor_B, f_factor_A))
				throw std::runtime_error("AddTo_FBS() failed");
			break;
		}

#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
		r_A.Rasterize("3_A+B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
	}
	// add A += B or B += A
}

void CBlockMatrixUnitTests::Test_AddFailure(CUberBlockMatrix &r_A, CUberBlockMatrix &r_B,
	bool b_polarity, int n_method) // throws(std::bad_alloc, std::runtime_error)
{
#ifdef __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES
	A.Rasterize("0_A.tga");
	B.Rasterize("1_B.tga");
#endif // __BLOCK_ADD_UNIT_TEST_DUMP_MATRIX_IMAGES

	typedef MakeTypelist_Safe((Eigen::Matrix<double, block_Size, block_Size>,
		Eigen::Matrix<double, block_Size2, block_Size>,
		Eigen::Matrix<double, block_Size, block_Size2>,
		Eigen::Matrix<double, block_Size2, block_Size2>)) two_typelist;

	if(b_polarity) {
		switch(n_method) {
		case 0:
			if(r_B.AddTo(r_A))
				throw std::runtime_error("AddTo() succeeded (should have failed)");
			break;
		case 1:
			if(r_B.AddTo(r_A, 1000.0))
				throw std::runtime_error("AddTo() succeeded (should have failed)");
			break;
		case 2:
			if(r_B.AddTo(r_A, 17.0, 19.0))
				throw std::runtime_error("AddTo() succeeded (should have failed)");
			break;
		case 3:
			if(r_B.AddTo_FBS<two_typelist>(r_A))
				throw std::runtime_error("AddTo_FBS() succeeded (should have failed)");
			break;
		case 4:
			if(r_B.AddTo_FBS<two_typelist>(r_A, 1000.0))
				throw std::runtime_error("AddTo_FBS() succeeded (should have failed)");
			break;
		default:
		case 5:
			if(r_B.AddTo_FBS<two_typelist>(r_A, 17.0, 19.0))
				throw std::runtime_error("AddTo_FBS() succeeded (should have failed)");
			break;
		}
	} else {
		switch(n_method) {
		case 0:
			if(r_A.AddTo(r_B))
				throw std::runtime_error("AddTo() succeeded (should have failed)");
			break;
		case 1:
			if(r_A.AddTo(r_B, 1000.0))
				throw std::runtime_error("AddTo() succeeded (should have failed)");
			break;
		case 2:
			if(r_A.AddTo(r_B, 17.0, 19.0))
				throw std::runtime_error("AddTo() succeeded (should have failed)");
			break;
		case 3:
			if(r_A.AddTo_FBS<two_typelist>(r_B))
				throw std::runtime_error("AddTo_FBS() succeeded (should have failed)");
			break;
		case 4:
			if(r_A.AddTo_FBS<two_typelist>(r_B, 1000.0))
				throw std::runtime_error("AddTo_FBS() succeeded (should have failed)");
			break;
		default:
		case 5:
			if(r_A.AddTo_FBS<two_typelist>(r_B, 17.0, 19.0))
				throw std::runtime_error("AddTo_FBS() succeeded (should have failed)");
			break;
		}
	}
	// add A += B or B += A
}

/*
 *								=== ~CBlockMatrixUnitTests ===
 */
