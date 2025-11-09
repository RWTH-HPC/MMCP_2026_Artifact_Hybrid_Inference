// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


#include "maiamath.h"

#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wduplicated-branches"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-dtor"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-result"
#pragma GCC diagnostic ignored "-Wnull-dereference" // Required for solveSparseMatrixIterative
#endif
#ifdef MAIA_CLANG_COMPILER
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wextra-semi-stmt"
#pragma clang diagnostic ignored "-Wfloat-equal"
#endif
#include <Eigen/Eigen/Eigen>
#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic pop
#endif
#ifdef MAIA_CLANG_COMPILER
#pragma clang diagnostic pop
#endif

namespace maia {
namespace math {

/// all assume dense matrix calculations
using MFloatMatrix2d = Eigen::Matrix<MFloat, 2, 2>;
using MFloatMatrix3d = Eigen::Matrix<MFloat, 3, 3>;
using MFloatMatrixXd = Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic>;
template <MInt nDim>
using MFloatMatrix = Eigen::Matrix<MFloat, nDim, nDim>;

using MFloatVector2d = Eigen::Vector<MFloat, 2>;
using MFloatVector3d = Eigen::Vector<MFloat, 3>;
using MFloatVectorXd = Eigen::Vector<MFloat, Eigen::Dynamic>;
template <MInt nDim>
using MFloatVector = Eigen::Vector<MFloat, nDim>;
using MTriplet = Eigen::Triplet<MFloat>;


/// convert array to eigen matrix
template <typename T, std::size_t N>
Eigen::MatrixXd ConvertToEigenMatrix(std::array<std::array<T, N>, N> data) {
  Eigen::Matrix<T, N, N> m;
  for(std::size_t i = 0; i < N; ++i) {
    m.row(i) = Eigen::VectorXd::Map(&data[i][0], N);
  }
  return m;
}


/// Row of adjoint of a matrix
/// \tparam T Type
/// \tparam N Dimension of the matrix
/// \param m Matrix (one dimensional array containing the matrix elements)
/// \param r Row of the adjuncte
/// \return Given row of the adjuncte of the matrix m in A
template <typename T, std::size_t N>
void adjointRow(std::array<std::array<T, N>, N>& m, std::array<T, N>& A, const MInt r) {
  Eigen::Matrix<T, N, N> m2 = ConvertToEigenMatrix(m);
  Eigen::Matrix<T, N, N> adjMatrix(m2.adjoint());
  for(std::size_t i = 0; i < N; i++) {
    A[i] = adjMatrix.row(r)[i];
  }
}
template void adjointRow<MFloat, 3>(std::array<std::array<MFloat, 3>, 3>& m, std::array<MFloat, 3>& A, const MInt);


/// First row of adjoint of a 4 x 4 matrix
/// \tparam T Type
/// \param m Matrix (two dimensional array containing the matrix elements)
/// \return First row of the adjuncte of the 4 x 4 matrix m in A
template <typename T>
void adjoint1stRow4x4(std::array<std::array<T, 4>, 4>& m, std::array<T, 4>& A) {
  A[0] = m[1][1] * m[2][2] * m[3][3] + m[1][2] * m[2][3] * m[3][1] + m[1][3] * m[2][1] * m[3][2]
         - m[1][3] * m[2][2] * m[3][1] - m[1][2] * m[2][1] * m[3][3] - m[1][1] * m[2][3] * m[3][2];
  A[1] = -m[0][1] * m[2][2] * m[3][3] - m[0][2] * m[2][3] * m[3][1] - m[0][3] * m[2][1] * m[3][2]
         + m[0][3] * m[2][2] * m[3][1] + m[0][2] * m[2][1] * m[3][3] + m[0][1] * m[2][3] * m[3][2];
  A[2] = m[0][1] * m[1][2] * m[3][3] + m[0][2] * m[1][3] * m[3][1] + m[0][3] * m[1][1] * m[3][2]
         - m[0][3] * m[1][2] * m[3][1] - m[0][2] * m[1][1] * m[3][3] - m[0][1] * m[1][3] * m[3][2];
  A[3] = -m[0][1] * m[1][2] * m[2][3] - m[0][2] * m[1][3] * m[2][1] - m[0][3] * m[1][1] * m[2][2]
         + m[0][3] * m[1][2] * m[2][1] + m[0][2] * m[1][1] * m[2][3] + m[0][1] * m[1][3] * m[2][2];
}
template void adjoint1stRow4x4<MFloat>(std::array<std::array<MFloat, 4>, 4>& m, std::array<MFloat, 4>& A);

/// First row of adjoint of a matrix
/// \tparam T Type
/// \param m Matrix (two dimensional array containing the matrix elements)
/// \return First row of the adjuncte of the matrix m in A
template <typename T, std::size_t N>
void adjoint1stRow(std::array<std::array<T, N>, N>& m, std::array<T, N>& A) {
  if constexpr(N == 4) {
    adjoint1stRow4x4(m, A);
  } else {
    adjointRow(m, A, 0);
  }
}
template void adjoint1stRow<MFloat, 4>(std::array<std::array<MFloat, 4>, 4>& m, std::array<MFloat, 4>& A);
template void adjoint1stRow<MFloat, 3>(std::array<std::array<MFloat, 3>, 3>& m, std::array<MFloat, 3>& A);

/// Determinant of a matrix
/// \tparam T Type
/// \tparam N Size of the matrix (rows * columns)
/// \param m Matrix (one dimensional array containing the matrix elements)
/// \return Determinant of the matrix
template <typename T, std::size_t N>
MFloat determinant(std::array<T, N>& m) {
  constexpr MInt dim = N == 2 ? 1 : 2;
  Eigen::Matrix<T, dim, dim> matrix(m.data());
  return matrix.determinant();
}
template MFloat determinant<MFloat, 2ul>(std::array<MFloat, 2ul>&);
template MFloat determinant<MFloat, 3ul>(std::array<MFloat, 3ul>&);
template MFloat determinant<MFloat, 4ul>(std::array<MFloat, 4ul>&);

/// Determinant of a matrix
/// \tparam T Type
/// \tparam N Dimension of the square matrix
/// \param m Matrix (two dimensional array containing the matrix elements)
/// \return Determinant of the matrix
template <typename T, std::size_t N>
MFloat determinant(std::array<std::array<T, N>, N>& m) {
  Eigen::Matrix<T, N, N> matrix(&m[0][0]);
  return matrix.determinant();
}
template MFloat determinant<MFloat, 2ul>(std::array<std::array<MFloat, 2ul>, 2ul>&);
template MFloat determinant<MFloat, 3ul>(std::array<std::array<MFloat, 3ul>, 3ul>&);
template MFloat determinant<MFloat, 4ul>(std::array<std::array<MFloat, 4ul>, 4ul>&);

// inline void jacobiSVD(MFloat* A, MFloat* U, MFloat* S, const MInt m, const MInt n) {
//   Eigen::Map<Eigen::MatrixXd> mA(A, m, n);
//   Eigen::Map<Eigen::MatrixXd> mU(U, m, n);
//   Eigen::Map<Eigen::VectorXd> vS(S, std::min(m, n));
// #ifdef useEigenOld
//   Eigen::JacobiSVD<Eigen::MatrixXd> svd(mA, Eigen::ComputeThinU | Eigen::ComputeThinV);
// #else
//   Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mA);
// #endif
//   vS = svd.singularValues();
//   mU = svd.matrixU();
//   mA = svd.matrixV();
// }

// inline void jacobiSVD(MFloatTensor& A, MFloatTensor& U, MFloatTensor& S) {
//   jacobiSVD(&A(0, 0), &U(0, 0), &S(0), A.dim0(), A.dim1());
// }

void invert(MFloat* A, const MInt m, const MInt n) {
  Eigen::Map<Eigen::MatrixXd> mA(A, m, n);
  mA = mA.inverse();
}

/// Note: don't use this for new code!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// Use instead:
/// A.colPivHouseholderQR.solve(b) (general case)
/// A.LLT.solve(b) (if the matrix is symmetric positive definite) (fastest)
/// A.HouseHolderQR.solve(b) (need for speed) (less accurate)
///
///  Get the pseudo-inverse of A by using QR
/// \param A system matrix (m x n)
/// \param AInv pseudo-inverse
/// \param m number of columns
/// \param n number of rows
/// \author Sven Berger
template <class T>
void invert(T& A, T& AInv, const MInt m, const MInt n) {
  if(n == 0 || m == 0) {
    std::cerr << globalDomainId() << " Warning: empty eq. sys. (" << m << "x" << n << "), skip." << std::endl;
    return;
  }
  Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic> B(m, n);
  Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic> inv(n, m);
  for(MInt i = 0; i < m; i++) {
    for(MInt j = 0; j < n; j++) {
      B(i, j) = A(i, j);
    }
  }
  //  if (n != m) {
  // system is overdetermined/underdetermined
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> orth(B);
  inv = orth.pseudoInverse();
  // todo labels:totest SSE optimization leads to false-positive during sanitizer run
  //  }
  //  else {
  //    inv = B.inverse();
  //  }
  for(MInt i = 0; i < m; i++) {
    for(MInt j = 0; j < n; j++) {
      AInv(i, j) = inv(i, j);
    }
  }
}
template void invert<ScratchSpace<MFloat>>(ScratchSpace<MFloat>&, ScratchSpace<MFloat>&, const MInt, const MInt);
template void invert<maia::tensor::Tensor<MFloat>>(maia::tensor::Tensor<MFloat>&, maia::tensor::Tensor<MFloat>&,
                                                   const MInt, const MInt);


/// Note: don't use this for new code!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// Use instead:
/// A.colPivHouseholderQR.solve(b) (general case)
/// A.LLT.solve(b) (if the matrix is symmetric positive definite) (fastest)
/// A.HouseHolderQR.solve(b) (need for speed) (less accurate)
///
///  Get the pseudo-inverse of A by using QR
/// \param A system matrix (m x n)
/// \param AInv pseudo-inverse
/// \param weights weights applied to A
/// \param m number of columns
/// \param n number of rows
/// \author Sven Berger
template <class T>
void invert(T& A, T& weights, T& AInv, const MInt m, const MInt n) {
  if(n == 0 || m == 0) {
    std::cerr << globalDomainId() << " Warning: empty eq. sys. (" << m << "x" << n << "), skip." << std::endl;
    return;
  }
  Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic> B(m, n);
  Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic> inv(n, m);
  for(MInt i = 0; i < m; i++) {
    for(MInt j = 0; j < n; j++) {
      B(i, j) = weights[i] * A(i, j);
    }
  }
  if(n != m) {
    // system is overdetermined/underdetermined
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> orth(B);
    inv = orth.pseudoInverse();
  } else {
    inv = B.inverse();
  }
  for(MInt k = 0; k < n; ++k) {
    for(MInt i = 0; i < m; ++i) {
      AInv(k, i) = inv(k, i) * weights(i);
    }
  }
}
template void invert<ScratchSpace<MFloat>>(ScratchSpace<MFloat>&, ScratchSpace<MFloat>&, ScratchSpace<MFloat>&,
                                           const MInt, const MInt);
template void invert<maia::tensor::Tensor<MFloat>>(maia::tensor::Tensor<MFloat>&, maia::tensor::Tensor<MFloat>&,
                                                   maia::tensor::Tensor<MFloat>&, const MInt, const MInt);

/// Note: don't use this for new code!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// Use instead:
/// A.colPivHouseholderQR.solve(b) (general case)
/// A.LLT.solve(b) (if the matrix is symmetric positive definite) (fastest)
/// A.HouseHolderQR.solve(b) (need for speed) (less accurate)
///
///  Get the pseudo-inverse of A by using QR
/// \param A system matrix (m x n)
/// \param AInv pseudo-inverse
/// \param weights weights applied to A
/// \param m number of columns
/// \param n number of rows
/// \return rank of the matrix A
/// \author Sven Berger
template <class T>
MInt invertR(T& A, T& weights, T& AInv, const MInt m, const MInt n) {
  if(n == 0 || m == 0) {
    std::cerr << globalDomainId() << " Warning: empty eq. sys. (" << m << "x" << n << "), skip." << std::endl;
    return -1;
  }
  Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic> B(m, n);
  Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic> inv(n, m);
  for(MInt i = 0; i < m; i++) {
    for(MInt j = 0; j < n; j++) {
      B(i, j) = weights[i] * A(i, j);
    }
  }

  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> orth(B);
  inv = orth.pseudoInverse();
  const MInt rank = orth.rank();

  for(MInt k = 0; k < n; ++k) {
    for(MInt i = 0; i < m; ++i) {
      AInv(k, i) = inv(k, i) * weights(i);
    }
  }
  return rank;
}
template MInt invertR<ScratchSpace<MFloat>>(ScratchSpace<MFloat>&, ScratchSpace<MFloat>&, ScratchSpace<MFloat>&,
                                            const MInt, const MInt);

/// Solve b = A * x for sparse matrix A using LU decomposition
/// \param A_coeff non-zero coefficients of the system matrix
/// \param pos holds the position of the coefficients in the matrix A
/// \param n number of non-zero coefficients
/// \param m size of system matrix A
/// \param b known vector to be calculated by the product A * x
/// \param x unknown vector to be calculated by A^(-1) * b
/// \author Moritz Waldmann
void solveSparseMatrix(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b, MFloat* x) {
  std::vector<MTriplet> coefficients;
  Eigen::Map<MFloatVectorXd> B(b, m);
  Eigen::Map<MFloatVectorXd> X(x, m);

  coefficients.reserve(n);
  for(MInt i = 0; i < n; i++) {
    coefficients.emplace_back(MTriplet(pos[i][0], pos[i][1], A_coeff[i]));
  }

  // set matrix from coefficients
  Eigen::SparseMatrix<MFloat, Eigen::ColMajor> A(m, m);
  A.setFromTriplets(coefficients.begin(), coefficients.end());

  Eigen::SparseLU<Eigen::SparseMatrix<MFloat, Eigen::ColMajor>, Eigen::COLAMDOrdering<MInt>> slu;
  slu.compute(A);

  X = slu.solve(B);
  MBool a_solution_exists = (A * X).isApprox(B, 1e-10);

  MFloat resultA = F1B2 * X.transpose() * A * X;
  MFloat resultB = X.transpose() * B;

  MFloat result = resultA - resultB;

  std::cout << "Energie norm is " << resultA << " " << resultB << " " << result << std::endl;

  if(a_solution_exists) {
    std::cerr << "SOLUTION IS CORRECT" << std::endl;
  } else {
    std::cerr << "SOLUTION IS NOT CORRECT" << std::endl;
  }
}

/// Solve b = A^(-1) * x for sparse matrix A using BiCGSTAB
/// \param A_coeff non-zero coefficients of the system matrix
/// \param pos holds the position of the coefficients in the matrix A
/// \param n number of non-zero coefficients
/// \param m size of system matrix A
/// \param b known vector to be calculated by the product A * x
/// \param x unknown vector to be calculated by A^(-1) * b
/// \author Moritz Waldmann
void solveSparseMatrixIterative(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b, MFloat* x) {
  std::vector<MTriplet> coefficients;
  Eigen::Map<MFloatVectorXd> B(b, m);
  Eigen::Map<MFloatVectorXd> X(x, m);

  coefficients.reserve(n);
  for(MInt i = 0; i < n; i++) {
    coefficients.emplace_back(MTriplet(pos[i][0], pos[i][1], A_coeff[i]));
  }

  // set matrix from coefficients
  Eigen::SparseMatrix<MFloat, Eigen::ColMajor> A(m, m);
  A.setFromTriplets(coefficients.begin(), coefficients.end());

  Eigen::BiCGSTAB<Eigen::SparseMatrix<MFloat, Eigen::ColMajor>> BCGST;
  BCGST.compute(A);

  X = BCGST.solve(B);

  MBool a_solution_exists = (A * X).isApprox(B, 1e-10);

  MFloat resultA = F1B2 * X.transpose() * A * X;
  MFloat resultB = X.transpose() * B;

  MFloat result = resultA - resultB;

  std::cout << "Energie norm is " << resultA << " " << resultB << " " << result << std::endl;

  if(a_solution_exists) {
    std::cout << "SOLUTION IS CORRECT" << std::endl;
  } else {
    std::cout << "SOLUTION IS NOT CORRECT" << std::endl;
  }
}

/// Solve x = A^(-1) * b for dense matrix A using LU decomposition
/// \param A_coeff non-zero coefficients of the system matrix
/// \param pos holds the position of the coefficients in the matrix A
/// \param n number of non-zero coefficients
/// \param m size of system matrix A
/// \param b known vector to be calculated by the product A * x
/// \param x unknown vector to be calculated by A^(-1) * b
/// \author Moritz Waldmann
void solveDenseMatrix(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b, MFloat* x) {
  Eigen::Map<MFloatVectorXd> B(b, m);
  Eigen::Map<MFloatVectorXd> X(x, m);

  MFloatMatrixXd A = MFloatMatrixXd::Zero(m, m);

  for(MInt i = 0; i < n; i++) {
    A(pos[i][0], pos[i][1]) = A_coeff[i];
  }

  Eigen::FullPivLU<MFloatMatrixXd> lu;
  lu.compute(A);

  X = lu.solve(B);
  MBool a_solution_exists = (A * X).isApprox(B, 1e-10);

  if(a_solution_exists) {
    std::cout << "SOLUTION IS CORRECT" << std::endl;
  } else {
    std::cout << "SOLUTION IS NOT CORRECT" << std::endl;
  }
}

//// Solve b = A * x for dense matrix A using Cholesky or LU
/// \param A_coeff non-zero coefficients of the system matrix
/// \param pos holds the position of the coefficients in the matrix A
/// \param n number of non-zero coefficients
/// \param m size of system matrix A
/// \param b unknown vector to be calculated by the product A * x
/// \param x known vector
/// \author Moritz Waldmann
void multiplySparseMatrixVector(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b, MFloat* x) {
  Eigen::Map<MFloatVectorXd> B(b, m);
  Eigen::Map<MFloatVectorXd> X(x, m);

  MFloatMatrixXd A = MFloatMatrixXd::Zero(m, m);

  for(MInt i = 0; i < n; i++) {
    A(pos[i][0], pos[i][1]) = A_coeff[i];
  }
  B = A * X;
}

///  Solve linear system using QR (square matrix)
/// \param A system matrix
/// \param b rhs
/// \tparam nDim dimension of the matrix
/// \author Sven Berger
template <MInt nDim>
void solveQR(std::array<std::array<MFloat, nDim>, nDim>& A_, std::array<MFloat, nDim>& b_) {
  TRACE();

  Eigen::Matrix<MFloat, nDim, nDim> A(&A_[0][0]);
  Eigen::Vector<MFloat, nDim> b(&b_[0]);
  b = A.colPivHouseholderQr().solve(b);
  std::memcpy(&b_[0], b.data(), nDim * sizeof(MFloat));
}
template void solveQR<2>(std::array<std::array<MFloat, 2>, 2>&, std::array<MFloat, 2>&);
template void solveQR<3>(std::array<std::array<MFloat, 3>, 3>&, std::array<MFloat, 3>&);

// ///  Solve linear system using QR (square matrix)
// /// \param A system matrix
// /// \param b rhs
// /// \tparam nDim dimension of the matrix
// /// \author Sven Berger
// inline void solveQRXd(MFloatMatrixXd& A, MFloatVectorXd& b, MFloatVectorXd& x) {
//   TRACE();
//   x = A.colPivHouseholderQr().solve(b);
// }


///  Get the Eigenvalues of matrix A (symmetrical)
/// \param A system matrix
/// \param w Eigenvalues
/// \author Sven Berger
void calcEigenValues(MFloat A[3][3], MFloat w[3]) {
  Eigen::Map<Eigen::Matrix<MFloat, 3, 3>> mA(&A[0][0]);
  Eigen::Map<Eigen::Matrix<MFloat, 3, 1>> _w(w);

  // Assume matrix is symmetrical
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<MFloat, 3, 3>> es(mA, Eigen::EigenvaluesOnly);
  _w = es.eigenvalues().real();
}

///  Get the Eigenvalues of matrix A_in (symmetrical)
/// \param A_in system matrix
/// \param lambda_in vector of eigen values
/// \param m size of square matrix A
/// \author Moritz Waldmann
void calcEigenValues(MFloat** A_in, MFloat* lambda_in, const MInt m) {
  Eigen::Map<Eigen::MatrixXd> A(&A_in[0][0], m, m);
  Eigen::Map<Eigen::VectorXd> lambda(lambda_in, m);

  // Assume matrix is symmetrical
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Eigen::EigenvaluesOnly);
  lambda = es.eigenvalues().real();
}

///  Get the Eigenvalues and Eigenvectors of matrix A (symmetrical)
/// \param A system matrix
/// \param w Eigenvalues
/// \param Q Eigenvectors
/// \author Sven Berger
void calcEigenVectors(MFloat A[3][3], MFloat Q[3][3], MFloat w[3]) {
  Eigen::Map<Eigen::Matrix<MFloat, 3, 3>> mA(&A[0][0]);
  Eigen::Map<Eigen::Matrix<MFloat, 3, 3>> mQ(&Q[0][0]);
  Eigen::Map<Eigen::Matrix<MFloat, 3, 1>> _w(w);

  // Assume matrix is symmetrical
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<MFloat, 3, 3>> es(mA);
  mQ = es.eigenvectors().real();
  _w = es.eigenvalues().real();
}

std::vector<MFloat> svd(MFloat* const A, MFloat* const b, const MInt m, const MInt n, MFloat* const x) {
  Eigen::Map<Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A_eigen(A, m, n);

  Eigen::Map<MFloatVectorXd> b_eigen(b, m);
  //      Eigen::Map<maia::math::MFloatVector<Eigen::Dynamic>> S_eigen(&S[0], std::min(m, n));

#ifdef useEigenOld
  Eigen::BDCSVD<Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> SVD =
      A_eigen.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
#else
  Eigen::BDCSVD<Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>,
                Eigen::ComputeThinU | Eigen::ComputeThinV>
      SVD = A_eigen.bdcSvd<Eigen::ComputeThinU | Eigen::ComputeThinV>();
#endif

  // TODO labels:CONTROLLER this should work as well
  MFloatVectorXd x_eigen = SVD.solve(b_eigen);
  //      Eigen::Matrix<MFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> x_eigen =
  //          SVD.solve(b_eigen);
  MFloatVectorXd S = SVD.singularValues();

  for(MInt i = 0; i < n; i++) {
    x[i] = x_eigen(i);
  }

  std::vector<MFloat> singularValues(S.size());
  for(MInt i = 0; i < S.size(); i++) {
    singularValues[i] = S[i];
  }
  return singularValues;
}

} // namespace math
} // namespace maia
