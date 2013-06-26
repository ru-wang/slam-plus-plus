/*
								+----------------------------------+
								|                                  |
								|         ***  COVADO  ***         |
								|                                  |
								|  Copyright (c) -tHE SWINe- 2013  |
								|                                  |
								|            Covado.cpp            |
								|                                  |
								+----------------------------------+
*/

/**
 *	@file src/slam/Covado.cpp
 *	@brief covariance conversion utility
 *	@author -tHE SWINe-
 *	@date 2013-06-14
 *
 *	This is a simple "utility" that we used to convert different covariance
 *	representation in datasets, e.g. to convert covariance of quaternion measurements
 *	to covariance of axis-angle measurements. Note that this conversion is not
 *	fully automated and requires editing this source code (it is not hard, though).
 *	Also, the resulting covariances are no longer diagonal matrices.
 */

/**
 *	@def __RUN_COVADO
 *	@brief define if in need of covariance conversion (hijacks the application)
 */
//#define __RUN_COVADO

/**
 *	@def __RUN_COVARIANCE_CONVERSION
 *	@brief define if in need of covariance conversion, otherwise runs cholesky
 */
//#define __RUN_COVARIANCE_CONVERSION

/**
 *	@def __COVADO_CONVERT_QUATCOV_TO_AXISANGLECOV
 *	@brief chooses quaternion covariance to axis-angle covariance (otherwise quaternion to RPY)
 */
//#define __COVADO_CONVERT_QUATCOV_TO_AXISANGLECOV

#ifdef __RUN_COVADO

#include "slam/Config.h"
#include "slam/SE2_Types.h"
#include "slam/SE3_Types.h"
#include "slam/BA_Types.h"

#include <iostream> // want to print Eigen matrices in easy way

/**
 *	@brief covariance conversion for dataset operation
 */
class CCovado {
public:
	/**
	 *	@brief default constructor; runs COVADO
	 */
	CCovado()
	{
#ifndef __RUN_COVARIANCE_CONVERSION
		printf("running cholesky\n");

/*

         100  2.15114e-05 -3.46442e-05    0.0101222    -0.128098    -0.342183
 2.15114e-05          100 -2.66435e-05  -0.00611334    -0.160674    -0.101645
-3.46442e-05 -2.66435e-05          100   -0.0230755     0.354741     0.117736
   0.0101222  -0.00611334   -0.0230755      156.068  -0.00913407    0.0357455
   -0.128098    -0.160674     0.354741  -0.00913407      2500.13     0.011439
   -0.342183    -0.101645     0.117736    0.0357455     0.011439      2499.07

          10            0            0            0            0            0
           0           10            0            0            0            0
           0            0           10            0            0            0
           0            0            0      12.4927            0            0
           0            0            0 -0.000731152      50.0013            0
           0            0            0   0.00286131  0.000228816      49.9907

10 0 0 0 0 0 10 0 0 0 0 10 0 0 0 12.4927 -0.000731152 0.00286131 50.0013 0.000228816 49.9907

*/ // YPR

		//const double p_u[] = {100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 2499.07, 0.011439, 0.0357455, 2500.13, -0.00913407, 156.068};
		//const double p_u[] = {100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 2498.43, 0.000198889, 0.0397551, 2498.43, -0.0116887, 156.089};
		const double p_u[] = {100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 156.068, -0.00913407, 0.0357455, 2500.13, 0.011439, 2499.07};

		Eigen::Matrix<double, 6, 6> information;
		information <<
			p_u[0],  p_u[1],  p_u[2],  p_u[3],  p_u[4],  p_u[5],
			p_u[1],  p_u[6],  p_u[7],  p_u[8],  p_u[9], p_u[10],
			p_u[2],  p_u[7], p_u[11], p_u[12], p_u[13], p_u[14],
			p_u[3],  p_u[8], p_u[12], p_u[15], p_u[16], p_u[17],
			p_u[4],  p_u[9], p_u[13], p_u[16], p_u[18], p_u[19],
			p_u[5], p_u[10], p_u[14], p_u[17], p_u[19], p_u[20];

		Eigen::LLT<Eigen::Matrix<double, 6, 6>, Eigen::Lower> chol_information(information);

		Eigen::Matrix<double, 6, 6> information_sqrt = chol_information.matrixL();

		std::cout << "information" << std::endl << information << std::endl;
		std::cout << "sqrt information" << std::endl << information_sqrt << std::endl;

		exit(-1);
#else // __RUN_COVARIANCE_CONVERSION
		printf("running covado\n");

		const double p_u[] = {100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 10000, 0, 0, 10000, 0, 625};
		//		      100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 10000, 0, 0, 10000, 0, 625

		Eigen::Matrix<double, 6, 6> information_g2o;
		information_g2o <<
			p_u[0],  p_u[1],  p_u[2],  p_u[3],  p_u[4],  p_u[5],
			p_u[1],  p_u[6],  p_u[7],  p_u[8],  p_u[9], p_u[10],
			p_u[2],  p_u[7], p_u[11], p_u[12], p_u[13], p_u[14],
			p_u[3],  p_u[8], p_u[12], p_u[15], p_u[16], p_u[17],
			p_u[4],  p_u[9], p_u[13], p_u[16], p_u[18], p_u[19],
			p_u[5], p_u[10], p_u[14], p_u[17], p_u[19], p_u[20];


		Eigen::Matrix<double, 6, 6> cov_g2o = information_g2o.inverse();

		std::cout << "g2o inf" << std::endl << information_g2o << std::endl;
		std::cout << "g2o cov" << std::endl << cov_g2o << std::endl;

		std::vector<Eigen::Matrix<double, 6, 1> > samples;
		for(int i = 0; i < 100000000; ++ i) {
			Eigen::Matrix<double, 6, 1> v;
			v << f_Rand(), f_Rand(), f_Rand(), f_Rand(), f_Rand(), f_Rand();
			// make uniform distribution 6D vector

			samples.push_back(v);
		}

		Eigen::Matrix<double, 6, 6> covariance;
		{
			Eigen::Matrix<double, 6, 1> v_average;
			v_average.setZero();
			for(size_t i = 0, n = samples.size(); i < n; ++ i)
				v_average += samples[i];
			v_average /= samples.size();
			// calculate mean

			std::cout << "mean" << std::endl << v_average << std::endl;

			covariance.setZero();
			for(size_t i = 0, n = samples.size(); i < n; ++ i)
				covariance += (samples[i] - v_average) * (samples[i] - v_average).transpose();
			covariance /= (samples.size() - 1);
		}
		// calculate covariance of the samples

		std::cout << "cov" << std::endl << covariance << std::endl;

		Eigen::LLT<Eigen::Matrix<double, 6, 6>, Eigen::Lower> cholCov(covariance.inverse());
		Eigen::LLT<Eigen::Matrix<double, 6, 6>, Eigen::Lower> cholg2oCov(cov_g2o);
		Eigen::Matrix<double, 6, 6> transform = Eigen::Matrix<double, 6, 6>(cholCov.matrixL()) * Eigen::Matrix<double, 6, 6>(cholg2oCov.matrixL());
		// calculate transform

		for(size_t i = 0, n = samples.size(); i < n; ++ i)
			samples[i] = transform * samples[i];

		{
			Eigen::Matrix<double, 6, 1> v_average;
			v_average.setZero();
			for(size_t i = 0, n = samples.size(); i < n; ++ i)
				v_average += samples[i];
			v_average /= samples.size();
			// calculate mean

			std::cout << "mean" << std::endl << v_average << std::endl;

			Eigen::Matrix<double, 6, 6> covariance2;
			covariance2.setZero();
			for(size_t i = 0, n = samples.size(); i < n; ++ i)
				covariance2 += (samples[i] - v_average) * (samples[i] - v_average).transpose();
			covariance2 /= (samples.size() - 1);
			// calculate covariance of the samples

			std::cout << "cov after transform" << std::endl << covariance2 << std::endl;
		}

#ifdef __COVADO_CONVERT_QUATCOV_TO_AXISANGLECOV
		for(size_t i = 0, n = samples.size(); i < n; ++ i) {
			Eigen::Vector3d v_pos = samples[i].segment<3>(0);
			Eigen::Vector3d v_rot = samples[i].segment<3>(3);

			Eigen::Quaterniond v_quat;
			v_quat.w() = sqrt(1 - v_rot.norm() * v_rot.norm()); // sin2 + cos2 = 1
			v_quat.vec() = v_rot;
			// transform the rotation to a (normalized) quaternion

			Eigen::Vector3d v_axis_angle;
			CBase3DSolver::C3DJacobians::Quat_to_AxisAngle(v_quat, v_axis_angle);

			samples[i].segment<3>(3) = v_axis_angle;
			// the rotation part changes
		}
#else // __COVADO_CONVERT_QUATCOV_TO_AXISANGLECOV
		for(size_t i = 0, n = samples.size(); i < n; ++ i) {
			Eigen::Vector3d v_pos = samples[i].segment<3>(0);
			Eigen::Vector3d v_rot = samples[i].segment<3>(3);

			Eigen::Quaterniond v_quat;
			v_quat.w() = sqrt(1 - v_rot.norm() * v_rot.norm()); // sin2 + cos2 = 1
			v_quat.vec() = v_rot;
			// transform the rotation to a (normalized) quaternion

			Eigen::Vector3d v_RPY;
			double q0 = v_quat.w(), q1 = v_quat.x(), q2 = v_quat.y(), q3 = v_quat.z(); // wikipedia notation
#if 0
			v_RPY.x() = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2));
			v_RPY.y() = asin(2 * (q0 * q2 - q3 * q1));
			v_RPY.z() = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));
			// RPY
#else
			v_RPY.z() = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2));
			v_RPY.y() = asin(2 * (q0 * q2 - q3 * q1));
			v_RPY.x() = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));
			// YPR
#endif
			// transform that quaternion to roll pitch yaw (blame [http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles])

			samples[i].segment<3>(3) = v_RPY;
			// the rotation part changes
		}
#endif // __COVADO_CONVERT_QUATCOV_TO_AXISANGLECOV

		{
			Eigen::Matrix<double, 6, 1> v_average;
			v_average.setZero();
			for(size_t i = 0, n = samples.size(); i < n; ++ i)
				v_average += samples[i];
			v_average /= samples.size();
			// calculate mean

			std::cout << "mean" << std::endl << v_average << std::endl;

			Eigen::Matrix<double, 6, 6> covariance2;
			covariance2.setZero();
			for(size_t i = 0, n = samples.size(); i < n; ++ i)
				covariance2 += (samples[i] - v_average) * (samples[i] - v_average).transpose();
			covariance2 /= (samples.size() - 1);
			// calculate covariance of the samples

			std::cout << "cov after space conversion" << std::endl << covariance2 << std::endl;

			std::cout << "information after space conversion" << std::endl << covariance2.inverse() << std::endl;
		}

		exit(-1);
#endif // __RUN_COVARIANCE_CONVERSION
	}

	static double f_Rand() // uniform distribution in [-1, 1]
	{
		return double(rand()) / RAND_MAX * 2 - 1;
	}
} covado;

/*

running covado
g2o inf
  100     0     0     0     0     0
    0   100     0     0     0     0
    0     0   100     0     0     0
    0     0     0 10000     0     0
    0     0     0     0 10000     0
    0     0     0     0     0   625
g2o cov
  0.01      0      0      0      0      0
     0   0.01      0      0      0      0
     0      0   0.01      0      0      0
     0      0      0 0.0001      0      0
     0      0      0      0 0.0001      0
     0      0      0      0      0 0.0016
mean
6.12858e-05
2.58306e-05
-4.77033e-05
0.000116611
3.85139e-05
-1.02417e-05
cov
     0.33334   4.8468e-05  5.63561e-06 -2.55678e-05 -9.45738e-06  1.80132e-05
  4.8468e-05     0.333349 -1.02368e-05 -7.46834e-06  -1.1752e-05 -1.08595e-05
 5.63561e-06 -1.02368e-05     0.333318  8.86376e-06  2.63536e-05 -4.10111e-05
-2.55678e-05 -7.46834e-06  8.86376e-06     0.333302 -1.49279e-06 -2.83143e-05
-9.45738e-06  -1.1752e-05  2.63536e-05 -1.49279e-06     0.333311  8.28715e-06
 1.80132e-05 -1.08595e-05 -4.10111e-05 -2.83143e-05  8.28715e-06     0.333277
mean
1.06149e-05
4.47234e-06
-8.26269e-06
2.02099e-06
6.68224e-07
-7.10914e-07
cov after transform
        0.01  2.41864e-11 -1.14956e-10  6.90354e-07  2.55384e-07  -3.2426e-07
 2.41864e-11         0.01  -3.9222e-12  2.01706e-07  3.17323e-07  1.95496e-07
-1.14956e-10  -3.9222e-12         0.01 -2.39243e-07 -7.11627e-07  7.38345e-07
 6.90354e-07  2.01706e-07 -2.39243e-07       0.0001  4.03325e-11 -2.55194e-08
 2.55384e-07  3.17323e-07 -7.11627e-07  4.03325e-11       0.0001  7.40846e-09
 -3.2426e-07  1.95496e-07  7.38345e-07 -2.55194e-08  7.40846e-09       0.0016
mean
1.06149e-05
4.47234e-06
-8.26269e-06
4.04273e-06
1.33738e-06
-1.42425e-06
cov after space conversion
        0.01  2.41864e-11 -1.14956e-10  1.38107e-06  5.10993e-07 -6.48934e-07
 2.41864e-11         0.01  -3.9222e-12  4.03601e-07  6.34849e-07  3.91598e-07
-1.14956e-10  -3.9222e-12         0.01 -4.78611e-07 -1.42374e-06  1.47798e-06
 1.38107e-06  4.03601e-07 -4.78611e-07  0.000400251  1.31997e-10 -1.02086e-07
 5.10993e-07  6.34849e-07 -1.42374e-06  1.31997e-10  0.000400251  2.97538e-08
-6.48934e-07  3.91598e-07  1.47798e-06 -1.02086e-07  2.97538e-08   0.00640658
information after space conversion
         100  2.13929e-05  -3.5038e-05    -0.345049    -0.127669    0.0101243
 2.13929e-05          100 -2.64658e-05    -0.100839    -0.158612  -0.00611329
 -3.5038e-05 -2.64658e-05          100     0.119572     0.355714   -0.0230695
   -0.345049    -0.100839     0.119572      2498.43  0.000198889    0.0397551
   -0.127669    -0.158612     0.355714  0.000198889      2498.43   -0.0116887
   0.0101243  -0.00611329   -0.0230695    0.0397551   -0.0116887      156.089

so the corresponding information matrix for axis-angle representations should be:
100 0 0 0 0 0 100 0 0 0 0 100 0 0 0 2498.43 0.000198889 0.0397551 2498.43 -0.0116887 156.089

running cholesky
information
        100           0           0           0           0           0
          0         100           0           0           0           0
          0           0         100           0           0           0
          0           0           0     2498.43 0.000198889   0.0397551
          0           0           0 0.000198889     2498.43  -0.0116887
          0           0           0   0.0397551  -0.0116887     156.089
sqrt information
          10            0            0            0            0            0
           0           10            0            0            0            0
           0            0           10            0            0            0
           0            0            0      49.9843            0            0
           0            0            0  3.97903e-06      49.9843            0
           0            0            0  0.000795352 -0.000233848      12.4936

so the corresponding information matrix for axis-angle square-rooted representations should be:
10 0 0 0 0 0 10 0 0 0 0 10 0 0 0 49.9843 3.97903e-06 0.000795352 49.9843 -0.000233848 12.4936

*/

/*

running covado
g2o inf
  100     0     0     0     0     0
    0   100     0     0     0     0
    0     0   100     0     0     0
    0     0     0 10000     0     0
    0     0     0     0 10000     0
    0     0     0     0     0   625
g2o cov
  0.01      0      0      0      0      0
     0   0.01      0      0      0      0
     0      0   0.01      0      0      0
     0      0      0 0.0001      0      0
     0      0      0      0 0.0001      0
     0      0      0      0      0 0.0016
mean
6.12858e-05
2.58306e-05
-4.77033e-05
0.000116611
3.85139e-05
-1.02417e-05
cov
     0.33334   4.8468e-05  5.63561e-06 -2.55678e-05 -9.45738e-06  1.80132e-05
  4.8468e-05     0.333349 -1.02368e-05 -7.46834e-06  -1.1752e-05 -1.08595e-05
 5.63561e-06 -1.02368e-05     0.333318  8.86376e-06  2.63536e-05 -4.10111e-05
-2.55678e-05 -7.46834e-06  8.86376e-06     0.333302 -1.49279e-06 -2.83143e-05
-9.45738e-06  -1.1752e-05  2.63536e-05 -1.49279e-06     0.333311  8.28715e-06
 1.80132e-05 -1.08595e-05 -4.10111e-05 -2.83143e-05  8.28715e-06     0.333277
mean
1.06149e-05
4.47234e-06
-8.26269e-06
2.02099e-06
6.68224e-07
-7.10914e-07
cov after transform
        0.01  2.41864e-11 -1.14956e-10  6.90354e-07  2.55384e-07  -3.2426e-07
 2.41864e-11         0.01  -3.9222e-12  2.01706e-07  3.17323e-07  1.95496e-07
-1.14956e-10  -3.9222e-12         0.01 -2.39243e-07 -7.11627e-07  7.38345e-07
 6.90354e-07  2.01706e-07 -2.39243e-07       0.0001  4.03325e-11 -2.55194e-08
 2.55384e-07  3.17323e-07 -7.11627e-07  4.03325e-11       0.0001  7.40846e-09
 -3.2426e-07  1.95496e-07  7.38345e-07 -2.55194e-08  7.40846e-09       0.0016
mean
1.06149e-05
4.47234e-06
-8.26269e-06
4.05567e-06
1.385e-06
-1.4247e-06
cov after space conversion
        0.01  2.41864e-11 -1.14956e-10  1.36925e-06  5.12359e-07 -6.48861e-07
 2.41864e-11         0.01  -3.9222e-12  4.06722e-07  6.42664e-07  3.91656e-07
-1.14956e-10  -3.9222e-12         0.01 -4.71132e-07 -1.41889e-06  1.47858e-06
 1.36925e-06  4.06722e-07 -4.71132e-07  0.000400148 -1.66802e-09 -9.17922e-08
 5.12359e-07  6.42664e-07 -1.41889e-06 -1.66802e-09   0.00039998  2.31919e-08
-6.48861e-07  3.91656e-07  1.47858e-06 -9.17922e-08  2.31919e-08   0.00640748
information after space conversion
         100  2.15114e-05 -3.46442e-05    -0.342183    -0.128098    0.0101222
 2.15114e-05          100 -2.66435e-05    -0.101645    -0.160674  -0.00611334
-3.46442e-05 -2.66435e-05          100     0.117736     0.354741   -0.0230755
   -0.342183    -0.101645     0.117736      2499.07     0.011439    0.0357455
   -0.128098    -0.160674     0.354741     0.011439      2500.13  -0.00913407
   0.0101222  -0.00611334   -0.0230755    0.0357455  -0.00913407      156.068

so the corresponding information matrix for RPY representations should be:
100 0 0 0 0 0 100 0 0 0 0 100 0 0 0 2499.07 0.011439 0.0357455 2500.13 -0.00913407 156.068

running cholesky
information
        100           0           0           0           0           0
          0         100           0           0           0           0
          0           0         100           0           0           0
          0           0           0     2499.07    0.011439   0.0357455
          0           0           0    0.011439     2500.13 -0.00913407
          0           0           0   0.0357455 -0.00913407     156.068
sqrt information
         10           0           0           0           0           0
          0          10           0           0           0           0
          0           0          10           0           0           0
          0           0           0     49.9907           0           0
          0           0           0 0.000228823     50.0013           0
          0           0           0 0.000715043 -0.00018268     12.4927


so the corresponding information matrix for RPY square-rooted representations should be:
10 0 0 0 0 0 10 0 0 0 0 10 0 0 0 49.9907 0.000228823 0.000715043 50.0013 -0.00018268 12.4927

*/

#endif // __RUN_COVADO

//#define __DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
#ifdef __DUMP_RSS2013_PRESENTATION_ANIMATION_DATA

#include "slam/Config.h"
class CRSSImagesCompositor {
protected:
	TBmp *p_lambda_add;
	TBmp *p_lambda_hot;
	TBmp *p_lambda_cold;
	TBmp *p_chol_hot;
	TBmp *p_chol_cold;
	TBmp *p_chol_add;
	TBmp *p_rhs_add;
	TBmp *p_rhs_hot;
	TBmp *p_rhs_cold;
	TBmp *p_lhs_hot;
	TBmp *p_lhs_cold;

public:
	CRSSImagesCompositor()
		:p_lambda_add(0), p_lambda_hot(0), p_lambda_cold(0), p_chol_add(0),
		p_chol_hot(0), p_chol_cold(0), p_rhs_hot(0), p_lhs_hot(0),
		p_rhs_cold(0), p_lhs_cold(0), p_rhs_add(0)
	{
		int n_image_num = 17;
		int n_padding = 10; // px

		LoadImages(n_image_num);

		int n_mat_w = p_lambda_add->n_width;
		int n_mat_h = p_lambda_add->n_height;
		int n_vec_w = p_rhs_hot->n_width;
		int n_frame_width = 2 * n_mat_w + n_vec_w + 1 * n_padding + n_mat_w / 2;
		int n_frame_height = n_mat_w;
		TBmp *p_frame = TBmp::p_Alloc(n_frame_width, n_frame_height);

		TBmp *p_prev_chol = TBmp::p_Alloc(1, 1);
		p_prev_chol->Clear(-1);
		TBmp *p_prev_dx = TBmp::p_Alloc(1, 1);
		p_prev_dx->Clear(-1);

#if 0
		for(int i = 1; i <= n_image_num; ++ i) {
			LoadImages(i);
			// load all the images

			for(int n_pass = 0; n_pass < 4; ++ n_pass) {
				{
					p_frame->Clear(0xffffffffU);
					// clear the frame

					Blit(p_frame, (n_pass == 0)? p_lambda_add : (n_pass == 2)? p_lambda_hot : p_lambda_cold, 0, 0);
					// the first is lambda

					Blit(p_frame, (n_pass == 0)? p_rhs_add : (n_pass == 2)? p_rhs_hot : p_rhs_cold,
						/*2 **/ (n_mat_w + n_padding), 0);
					// rhs

					//if(n_pass == 1 /*|| n_pass == 2*/ || (i == n_image_num && n_pass > 1)) {
						Blit(p_frame, (n_pass == 0)? p_prev_chol : (n_pass == 1)? p_chol_hot : p_chol_cold,
							n_mat_w + n_vec_w + n_padding + n_mat_w / 2, 0);
						//Blit(p_frame, (n_pass == 0)? p_prev_chol : /*(n_pass == 1)? p_chol_hot :*/
						//	p_chol_cold, n_mat_w + n_vec_w + n_padding + n_mat_w / 2, 0);
						// then cholesky

						//Blit(p_frame, (n_pass == 0)? p_prev_dx : (n_pass == 1)? p_lhs_hot : p_lhs_cold,
						//	n_mat_w + n_padding, n_mat_w + n_padding); // no dX
						// lhs
					//}
				}

				char p_s_filename[256];
				sprintf(p_s_filename, "rss2013/anim_%05d_%d.tga", i, n_pass);
				CTgaCodec::Save_TGA(p_s_filename, *p_frame, true, true);
				// save as an image
			}

			p_prev_chol->Delete();
			p_prev_chol = p_chol_cold->p_Clone();
			p_prev_dx->Delete();
			p_prev_dx = p_lhs_cold->p_Clone();
		}
		// make RSS-style very flashing animation
#else // 0
		for(int i = 1; i <= n_image_num; ++ i) {
			LoadImages(i);
			// load all the images

			for(int n_pass = 0; n_pass < 2; ++ n_pass) {
				p_frame->Clear(0xffffffffU);
				// clear the frame

				Blit(p_frame, p_lambda_add, 0, 0);
				// the first is lambda

				Blit(p_frame, p_rhs_add,
					/*2 **/ (n_mat_w + n_padding), 0);
				// rhs

				//Blit(p_frame, p_chol_add, n_mat_w + n_vec_w + n_padding + n_mat_w / 2, 0);
				if(n_pass)
					Blit(p_frame, p_chol_cold, n_mat_w + n_vec_w + n_padding + n_mat_w / 2, 0);

				char p_s_filename[256];
				sprintf(p_s_filename, "rss2013/anim_%05d_%d.tga", i, n_pass);
				CTgaCodec::Save_TGA(p_s_filename, *p_frame, true, true);
				// save as an image
			}

			p_prev_chol->Delete();
			p_prev_chol = p_chol_cold->p_Clone();
			p_prev_dx->Delete();
			p_prev_dx = p_lhs_cold->p_Clone();
		}
		// make IAV-style incremental animation
#endif // 0

		p_frame->Delete();
		p_prev_chol->Delete();
		p_prev_dx->Delete();

		int n_anim_frame_num = 4 * n_image_num; // number of frames in the sequence
		float n_anim_FPS = 25;
		float f_anim_wait_start = 5;//10;
		float f_anim_wait_end = 10; // a bit over-time but at least it won't finish in the slide
		float p_anim_wait_1stloop_stages[] = {10, 6, 7, 5};
		float f_anim_wait_1stloop_stages_sum =
			p_anim_wait_1stloop_stages[0] + p_anim_wait_1stloop_stages[1] +
			p_anim_wait_1stloop_stages[2] + p_anim_wait_1stloop_stages[3];
		float f_anim_length = 55;//60; // seconds

		int n_anim_wait_start_frame_num = int(ceil(n_anim_FPS * f_anim_wait_start));
		int n_anim_wait_end_frame_num = int(ceil(n_anim_FPS * f_anim_wait_end));
		int p_anim_wait_1stloop_stages_frame_num[4];
		for(int i = 0; i < 4; ++ i)
			p_anim_wait_1stloop_stages_frame_num[i] = int(ceil(n_anim_FPS * p_anim_wait_1stloop_stages[i]));
		int n_anim_wait_next_stages_frame_num = int(ceil(n_anim_FPS * (f_anim_length -
			f_anim_wait_start - f_anim_wait_end - f_anim_wait_1stloop_stages_sum) / (n_anim_frame_num - 5))); // except the first ones and the last one

		int n_wait_slide_id = 10;

		FILE *p_fw;
		if((p_fw = fopen("rss2013\\list.txt", "w"))) {
			{
				int i = 0/*n_image_num*/, n_pass = 3;
				char p_s_filename[256];
				sprintf(p_s_filename, "anim_%05d_%d.png", i, n_pass); // start with the /*last*/ zero-th, manually edited frame

				for(int j = 0; j < n_anim_wait_start_frame_num; ++ j)
					fprintf(p_fw, "%s\n", p_s_filename);
			}
			// wait at the beginning

			for(int i = 1; i <= n_image_num; ++ i) {
				for(int n_pass = 0; n_pass < ((i == n_image_num)? 3 : 4); ++ n_pass) { // the last one displayed explicitly
					char p_s_filename[256];
					sprintf(p_s_filename, "anim_%05d_%d.png", i, n_pass);
					int n_frame_num = (i == n_wait_slide_id)? p_anim_wait_1stloop_stages_frame_num[n_pass] :
						n_anim_wait_next_stages_frame_num;
					for(int j = 0; j < n_frame_num; ++ j)
						fprintf(p_fw, "%s\n", p_s_filename);
				}
			}
			// put all the other animation frames

			{
				int i = n_image_num, n_pass = 3;
				char p_s_filename[256];
				sprintf(p_s_filename, "anim_%05d_%d.png", i, n_pass); // end with the last frame

				for(int j = 0; j < n_anim_wait_end_frame_num; ++ j)
					fprintf(p_fw, "%s\n", p_s_filename);
			}
			// wait at the end

			fclose(p_fw);
		}

		if((p_fw = fopen("rss2013\\list2.txt", "w"))) {
			{
				for(int j = 0; j < n_anim_wait_start_frame_num; ++ j)
					fprintf(p_fw, "%s\n", "arrow_blank.png"); // no arrow
			}
			// wait at the beginning

			for(int i = 1; i <= n_image_num; ++ i) {
				for(int n_pass = 0; n_pass < ((i == n_image_num)? 3 : 4); ++ n_pass) { // the last one displayed explicitly
					char p_s_filename[256];
					sprintf(p_s_filename, "arrow_step%d.png", n_pass);
					int n_frame_num = (i == n_wait_slide_id)? p_anim_wait_1stloop_stages_frame_num[n_pass] :
						n_anim_wait_next_stages_frame_num;
					for(int j = 0; j < n_frame_num; ++ j)
						fprintf(p_fw, "%s\n", p_s_filename);
				}
			}
			// put all the other animation frames

			{
				for(int j = 0; j < n_anim_wait_start_frame_num; ++ j)
					fprintf(p_fw, "%s\n", "arrow_step3.png"); // no arrow
			}
			// wait at the end

			fclose(p_fw);
		}

		if((p_fw = fopen("rss2013\\list3.txt", "w"))) {
			{
				for(int j = 0; j < n_anim_wait_start_frame_num; ++ j)
					fprintf(p_fw, "%s\n", "arrow_blank.png"); // no arrow
			}
			// wait at the beginning

			for(int i = 1; i <= n_image_num; ++ i) {
				for(int n_pass = 0; n_pass < ((i == n_image_num)? 3 : 4); ++ n_pass) { // the last one displayed explicitly
					char p_s_filename[256];
					sprintf(p_s_filename, "arrow_step%d.png", n_pass);
					int n_frame_num = (i == n_wait_slide_id)? p_anim_wait_1stloop_stages_frame_num[n_pass] :
						n_anim_wait_next_stages_frame_num;
					for(int j = 0; j < n_frame_num; ++ j)
						fprintf(p_fw, "%s\n", (i != n_wait_slide_id)? "arrow_blank.png" : p_s_filename);
				}
			}
			// put all the other animation frames

			{
				for(int j = 0; j < n_anim_wait_start_frame_num; ++ j)
					fprintf(p_fw, "%s\n", "arrow_step3.png"); // no arrow
			}
			// wait at the end

			fclose(p_fw);
		}

		exit(-3);
	}

	void Blit(TBmp *p_dest, const TBmp *p_src, int dx, int dy)
	{
		const int dw = p_dest->n_width, dh = p_dest->n_height;
		for(int y = 0, w = p_src->n_width, h = p_src->n_height; y < h; ++ y) {
			for(int x = 0; x < w; ++ x) {
				if(x + dx >= 0 && x + dx < dw && y + dy >= 0 && y + dy < dh)
					p_dest->p_buffer[x + dx + (y + dy) * dw] = p_src->p_buffer[x + y * w];
			}
		}
	}

	~CRSSImagesCompositor()
	{
		DeleteImages(); // !!
	}

	void LoadImages(int n)
	{
		DeleteImages(); // !!

		char p_s_filename[256];
		sprintf(p_s_filename, "rss2013/%05d_0_lambda.tga", n);
		p_lambda_cold = CTgaCodec::p_Load_TGA(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_1_lambda2.tga", n);
		p_lambda_add = CTgaCodec::p_Load_TGA(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_2_lambda3.tga", n);
		p_lambda_hot = CTgaCodec::p_Load_TGA(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_a_Lnoord_red.tga", n);
		p_chol_hot = CTgaCodec::p_Load_TGA(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_b_Lnoord_inc.tga", n);
		p_chol_add = CTgaCodec::p_Load_TGA(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_3_Lnoord.tga", n);
		p_chol_cold = CTgaCodec::p_Load_TGA(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_6_rhs_red.tga", n);
		p_rhs_hot = CTgaCodec::p_Load_TGA(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_7_rhs_nnz.tga", n);
		p_rhs_cold = CTgaCodec::p_Load_TGA(p_s_filename);
		sprintf(p_s_filename, "rss2013/%05d_9_rhs_add.tga", n);
		p_rhs_add = CTgaCodec::p_Load_TGA(p_s_filename);

		p_lhs_cold = p_rhs_cold->p_Clone(true);
		p_lhs_hot = p_rhs_hot->p_Clone(true);
		std::swap(p_lhs_cold->n_width, p_lhs_cold->n_height);
		std::swap(p_lhs_hot->n_width, p_lhs_hot->n_height);
		for(int y = 0, w = p_lhs_cold->n_width, h = p_lhs_cold->n_height; y < h; ++ y) {
			for(int x = 0; x < w; ++ x) {
				p_lhs_hot->p_buffer[x + w * y] = p_rhs_hot->p_buffer[y + h * x];
				p_lhs_cold->p_buffer[x + w * y] = p_rhs_cold->p_buffer[y + h * x];
			}
		}
		// transpose here
	}

	void DeleteImages()
	{
		if(p_lambda_add) {
			p_lambda_add->Delete();
			p_lambda_add = 0;
		}
		if(p_lambda_cold) {
			p_lambda_cold->Delete();
			p_lambda_cold = 0;
		}
		if(p_lambda_hot) {
			p_lambda_hot->Delete();
			p_lambda_hot = 0;
		}
		if(p_chol_add) {
			p_chol_add->Delete();
			p_chol_add = 0;
		}
		if(p_chol_hot) {
			p_chol_hot->Delete();
			p_chol_hot = 0;
		}
		if(p_rhs_hot) {
			p_rhs_hot->Delete();
			p_rhs_hot = 0;
		}
		if(p_rhs_cold) {
			p_rhs_cold->Delete();
			p_rhs_cold = 0;
		}
		if(p_rhs_add) {
			p_rhs_add->Delete();
			p_rhs_add = 0;
		}
		if(p_lhs_hot) {
			p_lhs_hot->Delete();
			p_lhs_hot = 0;
		}
		if(p_lhs_cold) {
			p_lhs_cold->Delete();
			p_lhs_cold = 0;
		}
	}
} compo;
// composite rss-generated images

#endif // __DUMP_RSS2013_PRESENTATION_ANIMATION_DATA
