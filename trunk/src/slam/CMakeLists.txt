
# Set required variables.
SET(SLAM_P_P_BINARY_NAME "slam_plus_plus")

FILE(GLOB_RECURSE SLAM_P_P_HEADERS ${SLAM_P_P_SOURCE_DIR}/include/slam/*.h)

ADD_DEFINITIONS(-DEIGEN_MPL2_ONLY)

IF(SLAM_P_P_SHORT_CXSPARSE)
	ADD_DEFINITIONS(-D__CXSPARSE_SHORT)
ENDIF(SLAM_P_P_SHORT_CXSPARSE)

IF(SLAM_P_P_SHORT_CHOLMOD)
	ADD_DEFINITIONS(-D__CHOLMOD_SHORT)
ENDIF(SLAM_P_P_SHORT_CHOLMOD)

# support for GPU BLAS
IF(SLAM_P_P_GPU_BLAS)
	ADD_DEFINITIONS(
		-DGPU_BLAS
		-D__CHANNEL_DESCRIPTOR_H__
		-D__CUDA_RUNTIME_H__)
ENDIF(SLAM_P_P_GPU_BLAS)

ADD_DEFINITIONS(-DCERES_NO_PROTOCOL_BUFFERS)

SET(SLAM_P_P_SOURCES
#  ${PROJECT_SOURCE_DIR}/external/Ceres/internal/ceres/block_sparse_matrix.cc
#  ${PROJECT_SOURCE_DIR}/external/Ceres/internal/ceres/block_structure.cc
#  ${PROJECT_SOURCE_DIR}/external/Ceres/internal/ceres/sparse_matrix.cc
#  ${PROJECT_SOURCE_DIR}/external/Ceres/internal/ceres/triplet_sparse_matrix.cc
#  ${PROJECT_SOURCE_DIR}/external/Ceres/internal/ceres/linear_operator.cc
#  ${PROJECT_SOURCE_DIR}/external/Ceres/internal/ceres/miniglog/glog/logging.cc
  BlockBench.cpp
  BlockMatrix.cpp
  BlockUnit.cpp
  LinearSolver_CholMod.cpp
  LinearSolver_CSparse.cpp
  LinearSolver_CXSparse.cpp
  OrderingMagic.cpp
  Main.cpp
  Parser.cpp
  Tga.cpp
  Timer.cpp)

IF(SLAM_P_P_LINEAR_SOLVER_TYPE)
	ADD_DEFINITIONS("-D__LINEAR_SOLVER_OVERRIDE=${SLAM_P_P_LINEAR_SOLVER_TYPE}")
ENDIF(SLAM_P_P_LINEAR_SOLVER_TYPE)

# Build and install SLAM_plus_plus executable.
ADD_EXECUTABLE(${SLAM_P_P_BINARY_NAME} ${SLAM_P_P_SOURCES} ${SLAM_P_P_HEADERS})
TARGET_LINK_LIBRARIES(${SLAM_P_P_BINARY_NAME} ${SLAM_P_P_LIBRARIES})
INSTALL(TARGETS ${SLAM_P_P_BINARY_NAME}
  RUNTIME DESTINATION bin)
