// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SBPINTERPOLATION_H_
#define SBPINTERPOLATION_H_

#include <algorithm>
#include <iomanip>

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/tensor.h"
#include "enums.h"
#include "sbpcartesianoperators.h"

namespace maia {
namespace dg {

/**
 * \brief Holds helper functions for the construction
 *        of operators and interpolation in SBP mode.
 *
 * \author Julian Vorspohl<j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 */
namespace interpolation {


/**
 * \brief Generates an equidistant node distribution.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 * \param[in] noNodes Number of nodes
 * \param[out] nodes Node distribution
 */
inline void calcSBPNodes(const MInt noNodes, MFloat* nodes) {
  TRACE();

  ASSERT(noNodes > 0, "Number of nodes must be greater than zero!");

  // Create equidistant nodes on [-1,1]
  const MFloat dx = F2 / (noNodes - 1);
  std::fill_n(&nodes[0], noNodes, F0);
  for(MInt i = 0; i < noNodes; i++) {
    nodes[i] = -F1 + i * dx;
  }
}

/**
 * \brief Calulates the diagonal weight matrix (H) entries.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 * \param[in] noNodes Number of nodes
 * \param[in] sbpP Coefficients for the boundary nodes
 * \param[out] wInt Integration matrix (1D diagonal)
 */
inline void calcSBPWeights(const MInt noNodes, MFloatVector sbpP, MFloat* wInt) {
  TRACE();

  ASSERT(noNodes > 0, "Number of nodes must be greater than zero!");

  const MInt lenP = sbpP.dim0();
  const MFloat dx = F2 / (noNodes - 1);

  // Create (diagonal) H Matrix as Vector from 'p' Vector
  std::fill_n(&wInt[0], noNodes, dx);
  for(MInt i = 0; i < lenP; i++) {
    wInt[i] = sbpP(i) * dx;
    wInt[noNodes - i - 1] = sbpP(i) * dx;
  }
}

/**
 * \brief Helper function for checking the SBP property of a given operator.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 * \param[in] Q (hopefully) SBP operator
 */
inline MBool checkSBPProp(const MFloatMatrix& Q) {
  const MInt noNodes = Q.dim(0);

  MFloatMatrix B(noNodes, noNodes);
  B.set(F0);
  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < noNodes; j++) {
      B(i, j) = Q(i, j) + Q(j, i);
    }
  }

  MBool sbpProp = true;
  if(!approx(B(0, 0), -F1, MFloatEps) || !approx(B(noNodes - 1, noNodes - 1), F1, MFloatEps)) {
    sbpProp = false;
  }

  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < noNodes; j++) {
      if(!approx(B(i, j), 0.0, MFloatEps) && i != j) {
        sbpProp = false;
        break;
      }
      if(!sbpProp) break;
    }
  }
  std::cout << "SBP Propery " << sbpProp << std::endl;
  return sbpProp;
}

/**
 * \brief Reads .csv-file on rank 0 and broadcasts it to all ranks
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2020-01-01
 *
 * \param[in] path Absolute path to .csv-file
 * \param[out] data Read data
 */
inline void readCSV(std::string path, std::vector<std::vector<MFloat>>& data) {
  std::ifstream stream(path);

  MInt rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!rank) {
    // Read data on rank 0
    MInt l = 0;
    while(stream) {
      l++;
      std::string s;
      if(!getline(stream, s)) break;
      if(s[0] != '#') {
        std::istringstream ss(s);
        std::vector<MFloat> record;
        while(ss) {
          std::string line;
          if(!getline(ss, line, ',')) break;
          record.push_back(std::stof(line));
        }
        data.push_back(record);
      }
    }

    if(!stream.eof()) {
      std::cerr << "Could not read file " << path << "\n";
      std::__throw_invalid_argument("File not found.");
    }

    MInt dim0 = data.size();
    MInt dim1 = data[0].size();

    // transfer read data to sendBuffer
    MFloatVector sendBuffer;
    sendBuffer.resize(dim0 * dim1);
    MInt it = 0;
    for(auto& line : data) {
      for(auto& coef : line) {
        sendBuffer[it] = coef;
        it++;
      }
    }

    // Send data
    MPI_Bcast(&dim0, 1, MPI_INT, 0, MPI_COMM_WORLD, AT_, "dim0");
    MPI_Bcast(&dim1, 1, MPI_INT, 0, MPI_COMM_WORLD, AT_, "dim1");
    MPI_Bcast(&sendBuffer[0], dim0 * dim1, MPI_DOUBLE, 0, MPI_COMM_WORLD, AT_, "sendBuffer");

  } else {
    // Receive data
    MInt dim0;
    MInt dim1;
    MPI_Bcast(&dim0, 1, MPI_INT, 0, MPI_COMM_WORLD, AT_, "dim0");
    MPI_Bcast(&dim1, 1, MPI_INT, 0, MPI_COMM_WORLD, AT_, "dim1");

    MFloatVector recvBuffer;
    recvBuffer.resize(dim0 * dim1);
    MPI_Bcast(&recvBuffer[0], dim0 * dim1, MPI_DOUBLE, 0, MPI_COMM_WORLD, AT_, "recvBuffer");

    for(MInt i = 0; i < dim0; i++) {
      std::vector<MFloat> line;
      for(MInt j = 0; j < dim1; j++) {
        line.push_back(recvBuffer(dim1 * i + j));
      }
      data.push_back(line);
    }
  }
}

/**
 * \brief Reads DDRP operator and nodal distribution from file.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 * \param[out] nodesVector Nodal distribution from file
 * \param[out] wIntVector Integration matrix (1D diagonal)
 * \param[out] dHatMatrix DDRP operator (already scaled by wInt)
 */
inline void readDDRP(MFloat* nodesVector, MFloat* wIntVector, MFloat* dhatMatrix) {
  std::string path = "./operatorsSBP/go4/DDRP307_N32/";
  std::vector<std::vector<MFloat>> dData;
  std::vector<std::vector<MFloat>> wData;
  std::vector<std::vector<MFloat>> xData;

  readCSV(path + "diff_matrix.csv", dData);
  readCSV(path + "weights_int.csv", wData);
  readCSV(path + "nodes.csv", xData);

  const MInt nNodes = xData.size();

  MFloatMatrix Dhat(&dhatMatrix[0], nNodes, nNodes);
  MFloatVector wInt(&wIntVector[0], nNodes);
  MFloatVector nodes(&nodesVector[0], nNodes);

  // Copy the 2D Vector to the corresponding data structures
  // Transpose and normalize differentiation matrix D
  for(MInt i = 0; i < nNodes; i++) {
    for(MInt j = 0; j < nNodes; j++) {
      Dhat(i, j) = dData[j][i] * wData[j][0] / wData[i][0];
    }
    wInt(i) = wData[i][0];
    nodes(i) = xData[i][0];
  }
}

/**
 * \brief Calculates SBP operator from coefficients
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 * \param[in] noNodes Number of nodes
 * \param[in] sbpA Coefficients for the inner stencil
 * \param[in] sbpQ Coefficients for the boundary solver
 * \param[in] wInt Integration weights
 * \param[out] dHatMatrix SBP operator (already scaled by wInt)
 */
inline void calcDhatSBP(const MInt noNodes, MFloatVector sbpA, MFloatVector sbpQ, const MFloat* wInt,
                        MFloat* dhatMatrix) {
  TRACE();

  if(noNodes <= 0) {
    TERM(-1);
  }

  MInt lenA = sbpA.dim0();
  MInt lenQ = sbpQ.dim0();
  MInt lenR = (1 + sqrt(1 + 8 * lenQ)) / 2;

  // Create Q Matrix from 'q' and 'a' Vector
  MFloatMatrix Q = MFloatMatrix(noNodes, noNodes);
  Q.set(F0);

  // Write upper rxr solver
  Q(0, 0) = -0.5;
  MInt k = 0;
  for(MInt i = 0; i < lenR; i++) {
    for(MInt j = i; j < lenR - 1; j++) {
      Q(i, j + 1) = sbpQ(k);
      Q(j + 1, i) = -sbpQ(k);
      k++;
    }
  }

  // Write interior points
  for(MInt i = lenR; i < noNodes - lenR; i++) {
    k = lenA - 1;
    for(MInt j = i - lenA; j < i; j++) {
      Q(i, j) = -sbpA(k);
      k--;
    }
    k = 0;
    for(MInt j = i + 1; j < i + lenA + 1; j++) {
      Q(i, j) = sbpA(k);
      k++;
    }
  }

  // Write 'a' points around the upper rxr solver
  for(MInt i = 0; i < lenA; i++) {
    for(MInt j = i; j < lenA; j++) {
      Q(lenR - i - 1, lenR - i + j) = sbpA(j);
    }
  }

  // Mirror upper rxr solver onto lower half
  for(MInt i = 0; i < lenR; i++) {
    for(MInt j = 0; j < lenR + lenA; j++) {
      Q(noNodes - i - 1, noNodes - j - 1) = -Q(i, j);
    }
  }

  Q.transpose();

  // checkSBPProp(Q);

  MFloatMatrix Dhat(&dhatMatrix[0], noNodes, noNodes);
  // Add normalization using integration weights
  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < noNodes; j++) {
      Dhat(i, j) = -Q(i, j) / wInt[i];
    }
  }
}

/**
 * \brief Reads all coefficient files to contruct SBP operator
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 * \param[in] noNodes1D Number of nodes
 * \param[in] opName Name of the SBP operator
 * \param[out] sbpA Coefficients for the inner stencil
 * \param[out] sbpP Coefficients for the integration matrix
 * \param[out] sbpQ Coefficients for the boundary solver
 */
inline void readSbpOperator(MInt noNodes1D, MString opName, MFloatVector& sbpA, MFloatVector& sbpP,
                            MFloatVector& sbpQ) {
  MString path = "./operatorsSBP/" + opName;
  std::vector<std::vector<MFloat>> aData;
  std::vector<std::vector<MFloat>> pData;
  std::vector<std::vector<MFloat>> qData;

  readCSV(path + "/a.csv", aData);
  readCSV(path + "/p.csv", pData);
  readCSV(path + "/q.csv", qData);

  MInt lenA = aData.size();
  MInt lenP = pData.size();
  MInt lenR = lenP;
  MInt lenQ = lenP * (lenP - 1) / 2;

  ASSERT(noNodes1D >= 2 * lenR,
         "Not enough nodes for the " + opName + " operator! " << noNodes1D << " set but " << 2 * lenR << " required.");

  sbpA = MFloatVector(lenA);
  MInt i = 0;
  for(auto& a : aData) {
    sbpA(i) = a[0];
    i++;
  }

  sbpQ = MFloatVector(lenQ);
  i = 0;
  for(auto& q : qData) {
    sbpQ(i) = q[0];
    i++;
  }

  sbpP = MFloatVector(lenP);
  i = 0;
  for(auto& p : pData) {
    sbpP(i) = p[0];
    i++;
  }
}

/**
 * \brief Gets SBP Coefficients from corresponding header if existent
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2020-01-01
 *
 * \param[in] noNodes1D Number of nodes
 * \param[in] opName Name of SBP operator
 * \param[out] sbpA Coefficients for the inner stencil
 * \param[out] sbpP Coefficients for the integration matrix
 * \param[out] sbpQ Coefficients for the boundary solver
 *
 * \returns true if operator exists in header
 */
inline MBool getSBPOperator(const MInt noNodes1D, const MString opName, MFloatVector& sbpA, MFloatVector& sbpP,
                            MFloatVector& sbpQ) {
  SBPOperator* op = nullptr;
  MBool exists = false;

  if(opName == "go4/s306") {
    op = &s306;
    exists = true;
  } else {
    return exists = false;
  }

  MInt lenA = op->a.size();
  MInt lenP = op->p.size();
  MInt lenR = lenP;
  MInt lenQ = lenP * (lenP - 1) / 2;

  ASSERT(noNodes1D >= 2 * lenR,
         "Not enough nodes for the " + opName + " operator!" << noNodes1D << " set but " << 2 * lenR << " required.");

  sbpA = MFloatVector(lenA);
  MInt i = 0;
  for(auto& a : op->a) {
    sbpA(i) = a;
    i++;
  }

  sbpQ = MFloatVector(lenQ);
  i = 0;
  for(auto& q : op->q) {
    sbpQ(i) = q;
    i++;
  }

  sbpP = MFloatVector(lenP);
  i = 0;
  for(auto& p : op->p) {
    sbpP(i) = p;
    i++;
  }

  return exists;
}

/**
 * \brief Finds the neighboring indices to point in one dimension.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 * \note Even if point is collocating with a node, two different indices are returned.
 *
 * \param[in] point 1D Point
 * \param[in] nodes Nodal distribution
 * \param[in] noNodes Number of nodes
 * \param[out] idx1 First index (more negativer)
 */
inline void findNodeIndicesAtPoint(const MFloat point, const MFloat* nodes, const MInt noNodes, MInt& idx1,
                                   MInt& idx2) {
  // TODO labels:DG Change to bisection/Catch collocating points
  for(MInt i = 0; i < noNodes; i++) {
    if(nodes[i] > point || approx(nodes[i], point, MFloatEps)) {
      idx1 = i - 1;
      idx2 = i;
      break;
    }
  }
  if(idx1 == -1) {
    idx1++;
    idx2++;
  }
}

/**
 * \brief Calculates linear interpolation base for point (1D).
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-11-05
 *
 * \note f(x) = nodalValues * polynomials
 *
 * \param[in] x Point at which the solution is interpolated
 * \param[in] noNodes Number of nodes
 * \param[in] nodes Nodal distribution
 * \param[out] polynomials Coefficients for the linear interpolation
 */
inline void calcLinearInterpolationBase(const MFloat x, const MInt noNodes, const MFloat* nodes,
                                        MFloat* const polynomials) {
  std::fill_n(&polynomials[0], noNodes, F0);

  // Framing indices of point on grid
  MInt idx[2] = {0, 0};
  findNodeIndicesAtPoint(x, nodes, noNodes, idx[0], idx[1]);

  const MFloat dx = 2.0 / (noNodes - 1);
  polynomials[idx[0]] = fabs(nodes[idx[1]] - x) / dx;
  polynomials[idx[1]] = fabs(nodes[idx[0]] - x) / dx;
}

/**
 * \brief Calculates the linear interpolation matrix (Vandermonde) to
 *        interpolate from one set of nodes to another.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-05-14
 *
 * \tparam T,U,V Any container type that overloads operator[] for element
 *                 access (including pointers).
 * \param[in] noNodesIn Number of incoming nodes
 * \param[in] nodesIn Incoming nodal distribution
 * \param[in] noNodesOut Number of outgoing nodes
 * \param[in] nodesOut Outgoing nodal distribution
 * \param[out] vandermonde Interpolation matrix
 */
template <class T, class U, class V>
void calcLinearInterpolationMatrix(const MInt noNodesIn, const T nodesIn, const MInt noNodesOut, const U nodesOut,
                                   V vandermonde) {
  TRACE();

  MFloatMatrix vdm(&vandermonde[0], noNodesOut, noNodesIn);

  for(MInt k = 0; k < noNodesOut; k++) {
    MBool rowHasMatch = false;
    for(MInt j = 0; j < noNodesIn; j++) {
      vdm(k, j) = F0;
      if(approx(nodesOut[k], nodesIn[j], MFloatEps)) {
        rowHasMatch = true;
        vdm(k, j) = F1;
      }
    }

    if(rowHasMatch == false) {
      calcLinearInterpolationBase(nodesOut[k], noNodesIn, nodesIn, &vdm(k, 0));
    }
  }
}

/**
 * \brief Calculates bilinear interpolation at given point (2D)
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-05-14
 *
 * \param[in] point Point at wich interpolation is evaluated
 * \param[in] nodes Nodal distribution
 * \param[in] noNodes Number of nodes 1D
 * \param[in] noVars Number of variables per node
 * \param[in] u Solution field at all nodes
 * \param[out] out Interpolated solution at point
 */

inline void calcBilinearInterpolation(MFloat* point, const MFloat* nodes, const MInt noNodes, const MInt noVars,
                                      const MFloat* u, MFloat* const out) {
  // Framing indices of point on grid
  // nDim x 2
  MInt idx[2][2] = {{0, 0}, {0, 0}};
  findNodeIndicesAtPoint(point[0], nodes, noNodes, idx[0][0], idx[0][1]);
  findNodeIndicesAtPoint(point[1], nodes, noNodes, idx[1][0], idx[1][1]);

  MFloatTensor U(const_cast<MFloat*>(u), noNodes, noNodes, 1, noVars);
  const MFloat dx = 2.0 / (noNodes - 1);
  const MFloat dx2 = pow(dx, 2);
  // Not sure if necessary
  std::fill_n(out, noVars, 0.0);

  // aij is the corresponding weight to the value Uij
  // aij is the fraction of the diagonally opposite subarea

  for(MInt i = 0; i < 2; i++) {
    for(MInt j = 0; j < 2; j++) {
      const MFloat a = fabs((nodes[idx[0][1 - i]] - point[0]) * (nodes[idx[1][1 - j]] - point[1])) / dx2;
      for(MInt v = 0; v < noVars; v++) {
        out[v] += a * U(idx[0][i], idx[1][j], 0, v);
      }
    }
  }
}

/**
 * \brief Calculates trilinear interpolation at given point (3D)
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-05-14
 *
 * \param[in] point Point at wich interpolation is evaluated
 * \param[in] nodes Nodal distribution
 * \param[in] noNodes Number of nodes 1D
 * \param[in] noVars Number of variables per node
 * \param[in] u Solution field at all nodes
 * \param[out] out Interpolated solution at point
 */
inline void calcTrilinearInterpolation(MFloat* point, const MFloat* nodes, const MInt noNodes, const MInt noVars,
                                       const MFloat* u, MFloat* const out) {
  // Framing indices of point on grid
  // nDim x 2
  MInt idx[3][2] = {{0, 0}, {0, 0}, {0, 0}};
  findNodeIndicesAtPoint(point[0], nodes, noNodes, idx[0][0], idx[0][1]);
  findNodeIndicesAtPoint(point[1], nodes, noNodes, idx[1][0], idx[1][1]);
  findNodeIndicesAtPoint(point[2], nodes, noNodes, idx[2][0], idx[2][1]);

  MFloatTensor U(const_cast<MFloat*>(u), noNodes, noNodes, noNodes, noVars);
  const MFloat dx = 2.0 / (noNodes - 1);
  const MFloat dx3 = pow(dx, 3);
  // Not sure if necessary
  std::fill_n(out, noVars, 0.0);

  // aijk is the corresponding weight to the value Uijk
  // aijk is the fraction of the diagonally opposite subvolume

  for(MInt i = 0; i < 2; i++) {
    for(MInt j = 0; j < 2; j++) {
      for(MInt k = 0; k < 2; k++) {
        const MFloat a = fabs((nodes[idx[0][1 - i]] - point[0]) * (nodes[idx[1][1 - j]] - point[1])
                              * (nodes[idx[2][1 - k]] - point[2]))
                         / dx3;

        for(MInt v = 0; v < noVars; v++) {
          out[v] += a * U(idx[0][i], idx[1][j], idx[2][k], v);
        }
      }
    }
  }
}

// ~jv Rework
/*inline void calcBicubicSplineInterpolation(spline* spline, MFloat *point,
                                           const MFloat* nodes, const MInt noNodes,
                                           const MInt noVars, const MFloat* u,
                                           MFloat* const out){

  spline->set_boundary(spline->second_deriv,0,spline->second_deriv,0,false);
  std::vector<MFloat> X(noNodes), Utmp(noNodes), Uinter1D(noNodes);
  for (MInt i = 0; i < noNodes; ++i) {
    X[i] = nodes[i];
  }
  MFloatTensor U(const_cast<MFloat*>(u),noNodes,noNodes,1,noVars);

  for (MInt v = 0; v < noVars; ++v) {
    auto debug = false;
    if (v==2 && U(0,0,0,0)>1E30){
      debug=true;
      std::cout << "BICUBIC OUTPUT for p at " << point[0] << " " << point[1] << std::endl;}
    Uinter1D.clear();
    for (MInt i = 0; i < noNodes; ++i) {
      Utmp.clear();
      for (MInt j = 0; j < noNodes; ++j) {
        Utmp.push_back(U(i,j,0,v));
      }
      if (debug){
      std::cout << "TMP" << " " << point[1] << " " << X[i] << std::endl;
      //for (auto en : X){std::cout << en << std::endl;}
      //for (auto en : Utmp){std::cout << en << std::endl;}
      for (MInt k = 0; k < noNodes; ++k) {
        std::cout << X[k] << " : " << Utmp[k] << std::endl;
      }}
      //mTerm(1,AT_,"DEBUG");
      spline->set_points(X,Utmp);
      Uinter1D.push_back(spline->operator()(point[1]));
      if(debug)std::cout << "TMP RES " << Uinter1D[i] << std::endl;
    }
    spline->set_points(X,Uinter1D);
    out[v] = spline->operator()(point[0]);
    if(debug){
    std::cout << "END " << point[0] << std::endl;
    for (MInt k = 0; k < noNodes; ++k) {
      std::cout << X[k] << " : " << Uinter1D[k] << std::endl;
    }
    std::cout << "END RES " << out[v] << std::endl;}
    //if (out[v] <= -0.004){mTerm(1,AT_,"DEBUG TERM");}
  }

  //calcBilinearInterpolation(point,nodes,noNodes,noVars,u,out);
}
*/

} // namespace interpolation
} // namespace dg
} // namespace maia

#endif /* SBPINTERPOLATION_H_ */
