// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SBPMORTAR_H_
#define SBPMORTAR_H_

#include "INCLUDE/maiatypes.h"
#include "UTIL/tensor.h"
#include "sbpcartesianinterpolation.h"

#include <string>

namespace maia {
namespace sbp {
namespace mortar {


const MBool forward = true;
const MBool reverse = false;
const MInt lower = 0;
const MInt upper = 1;

/**
 * \brief Reads csv from path to MFloatMatrix
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-05-11
 *
 * \param[in] path Path to csv
 * \param[out] mat Matrix of read coefficients
 *
 * \returns true if read successfully
 */
inline MBool readCSV(std::string path, MFloat* mat) {
  std::vector<std::vector<MFloat>> data;
  std::ifstream stream(path);

  if(stream.fail()) {
    std::cerr << "File doesn't exist: " << path << "\n";
    return false;
  }

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
    return false;
  }

  const MInt m = data.size();
  const MInt n = data[0].size();

  MFloatMatrix matrix(mat, m, n);

  for(MInt i = 0; i < m; i++) {
    for(MInt j = 0; j < n; j++) {
      matrix(i, j) = data[i][j];
    }
  }

  return true;
}

/**
 * \brief Reads coeffficients from csv files for Carpenter operators
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-05-11
 *
 * \param[in] path Path to coefficients
 * \param[out] projBn Matrix of boundary coefficients
 * \param[out] projIn Vector of inner stencil
 */
inline void readProjectionCoeffs(std::string path, MFloatMatrix& projBn, MFloatVector& projIn) {
  projIn = MFloatVector(2);
  projIn.set(0.2);

  std::string bnFileName = path + "bn.csv";
  std::string inFileName = path + "in.csv";

  std::vector<std::vector<MFloat>> bnData;
  std::ifstream bnFile(bnFileName);
  MInt l = 0;

  while(bnFile) {
    l++;
    std::string s;
    if(!getline(bnFile, s)) break;
    if(s[0] != '#') {
      std::istringstream ss(s);
      std::vector<MFloat> record;

      while(ss) {
        std::string line;
        if(!getline(ss, line, ',')) break;
        record.push_back(std::stof(line));
      }
      bnData.push_back(record);
    }
  }

  if(!bnFile.eof()) {
    std::cerr << "Could not read file " << bnFileName << "\n";
    std::__throw_invalid_argument("File not found.");
  }

  std::vector<std::vector<MFloat>> inData;
  std::ifstream inFile(inFileName);
  l = 0;

  while(inFile) {
    l++;
    std::string s;
    if(!getline(inFile, s)) break;
    if(s[0] != '#') {
      std::istringstream ss(s);
      std::vector<MFloat> record;

      while(ss) {
        std::string line;
        if(!getline(ss, line, ',')) break;
        record.push_back(std::stof(line));
      }
      inData.push_back(record);
    }
  }

  if(!inFile.eof()) {
    std::cerr << "Could not read file " << inFileName << "\n";
    std::__throw_invalid_argument("File not found.");
  }

  const MInt q = bnData.size();
  const MInt r = bnData[0].size();

  projBn = MFloatMatrix(q, r);

  for(MInt i = 0; i < q; i++) {
    for(MInt j = 0; j < r; j++) {
      projBn(i, j) = bnData[i][j];
    }
  }

  const MInt numInternal = inData.size();

  projIn = MFloatVector(numInternal);

  for(MInt i = 0; i < numInternal; i++) {
    projIn(i) = inData[i][0];
  }
}

/**
 * \brief Reads projection coefficients and stores them
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-05-11
 *
 * \param[in] sourceOp SBP operator of source side
 * \param[in] targetOp SBP operator of target side
 * \param[in] sourceNoNodes Number of nodes on source side
 * \param[in] targetNoNodes Number of nodes on target side
 * \param[out] f forward projection operator
 * \param[out] b backward projection operator
 */
inline void calcMortarProjectionMatrixPSBP(const MString sourceOp,
                                           const MString targetOp,
                                           const MInt sourceNoNodes,
                                           const MInt targetNoNodes,
                                           MFloat* const f,
                                           MFloat* const b) {
  std::string dir = "./operatorsSBP/projection/p/";
  std::string name = sourceOp.substr(sourceOp.find('/') + 1) + targetOp.substr(targetOp.find('/') + 1);

  std::string pathF = dir + name + "/" + std::to_string(sourceNoNodes) + "_" + std::to_string(targetNoNodes) + "_f.csv";
  std::string pathB = dir + name + "/" + std::to_string(sourceNoNodes) + "_" + std::to_string(targetNoNodes) + "_b.csv";

  auto fRead = readCSV(pathF, f);
  auto bRead = readCSV(pathB, b);

  if(!fRead || !bRead) {
    std::cerr << "P projection matrix for  " << sourceOp << " with " << sourceNoNodes << " nodes to " << targetOp
              << " with " << targetNoNodes << " nodes NOT CALCULATED" << std::endl;
    MFloatMatrix fMatrix(f, targetNoNodes, sourceNoNodes);
    MFloatMatrix bMatrix(b, sourceNoNodes, targetNoNodes);
    fMatrix.set(0);
    bMatrix.set(0);
  } else {
    MFloatMatrix fMatrix(f, targetNoNodes, sourceNoNodes);
    MFloatMatrix bMatrix(b, sourceNoNodes, targetNoNodes);
  }
}

/**
 * \brief Constructs h-refinement projection operator according to Carpenter
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-05-11
 *
 * \param[in] noNodes Number of nodes 1D
 * \param[in] sourceIntWeights Integration weights of used SBP operator
 * \param[out] forwardLower Projection operator
 * \param[out] forwardUpper Projection operator
 * \param[out] backwardLower Projection operator
 * \param[out] backwardUpper Projection operator
 */
inline void calcMortarProjectionMatrixHSBPCarpenter(const MInt noNodes,
                                                    const MFloat* const sourceIntWeights,
                                                    MFloat* const forwardLower,
                                                    MFloat* const forwardUpper,
                                                    MFloat* const backwardLower,
                                                    MFloat* const backwardUpper) {
  // Create targetIntWeight by stacking up sourceIntWeight
  // Since wInt is unscaled we don't need to apply any factor here
  MFloatVector targetIntWeights(2 * noNodes);

  for(MInt i = 0; i < noNodes; i++) {
    targetIntWeights(i) = sourceIntWeights[i];
    targetIntWeights(i + noNodes) = sourceIntWeights[i];
  }

  MFloatMatrix stenBn;
  MFloatVector coeffsIn;
  std::string path = "./operatorsSBP/projection/s204/";
  readProjectionCoeffs(path, stenBn, coeffsIn);

  // TODO labels:DG s = #coeffs, carpenter s + 1 = #coeffs
  const MInt s = coeffsIn.dim0();
  const MInt q = stenBn.dim0();
  const MInt r = stenBn.dim1();

  // Create inner stencil
  // Coeffs: 1 2 3
  // Stencil 3 2 1 2 3
  MFloatVector stenIn(2 * s - 1);
  for(MInt i = 0; i < s; i++) {
    stenIn(s - 1 + i) = coeffsIn(i);
    stenIn(s - 1 - i) = coeffsIn(i);
  }

  MFloatMatrix f(2 * noNodes, noNodes);
  MFloatMatrix b(noNodes, 2 * noNodes);

  b.set(0.0);
  f.set(0.0);

  if(r > 2 * noNodes || 2 * q > noNodes || 2 * s - 1 > 2 * noNodes) {
    return;
  }

  // Set upper boundary solver
  for(MInt i = 0; i < q; i++) {
    for(MInt j = 0; j < r; j++) {
      b(i, j) = stenBn(i, j);
    }
  }

  // Set inner stencil
  MInt k = 2 * q - s + 1;
  for(MInt i = q; i < noNodes / 2; i++) {
    for(MInt j = 0; j < 2 * s - 1; j++) {
      b(i, j + k) = stenIn(j);
    }
    k += 2;
  }

  // Mirror
  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < 2 * noNodes; j++) {
      b(noNodes - i - 1, 2 * noNodes - j - 1) = b(i, j);
    }
  }

  // Create foward operator by applying the H-Compatibility condition eq.
  // Apply dx scale factor since wInt is unscaled.
  // For the 1:2 href case it's always 2
  // TODO labels:DG scale beforehand
  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < 2 * noNodes; j++) {
      f(j, i) = b(i, j) * sourceIntWeights[i] / targetIntWeights(j) * 2;
    }
  }

  // Split operators into upper and lower part
  MFloatMatrix fl(forwardLower, noNodes, noNodes);
  MFloatMatrix fu(forwardUpper, noNodes, noNodes);
  MFloatMatrix bl(backwardLower, noNodes, noNodes);
  MFloatMatrix bu(backwardUpper, noNodes, noNodes);

  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < noNodes; j++) {
      fl(i, j) = f(i, j);
      fu(i, j) = f(i + noNodes, j);
      bl(i, j) = b(i, j);
      bu(i, j) = b(i, j + noNodes);
    }
  }
}

/**
 * \brief Reads and constructs h-refinement projection operator
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 2019-05-11
 *
 * \param[in] noNodes Number of nodes 1D
 * \param[in] sbpOperator Used SBP operator
 * \param[out] forwardLower Projection operator
 * \param[out] forwardUpper Projection operator
 * \param[out] backwardLower Projection operator
 * \param[out] backwardUpper Projection operator
 */
inline void calcMortarProjectionMatrixHSBP(const MInt noNodes, const MString sbpOperator, MFloat* const forwardLower,
                                           MFloat* const forwardUpper, MFloat* const backwardLower,
                                           MFloat* const backwardUpper) {
  std::string path = "./operatorsSBP/projection/h/";

  std::string pathF = path + sbpOperator + "/" + std::to_string(noNodes) + "_f.csv";
  std::string pathB = path + sbpOperator + "/" + std::to_string(noNodes) + "_b.csv";

  MFloatMatrix f(2 * noNodes, noNodes);
  MFloatMatrix b(noNodes, 2 * noNodes);

  auto fRead = readCSV(pathF, &f[0]);
  auto bRead = readCSV(pathB, &b[0]);

  if(!fRead || !bRead) {
    std::cerr << "H projection matrix for " << noNodes << " nodes NOT CALCULATED" << std::endl;
    f.set(0);
    b.set(0);
  }

  // Split operators into upper and lower part
  MFloatMatrix fl(forwardLower, noNodes, noNodes);
  MFloatMatrix fu(forwardUpper, noNodes, noNodes);
  MFloatMatrix bl(backwardLower, noNodes, noNodes);
  MFloatMatrix bu(backwardUpper, noNodes, noNodes);

  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < noNodes; j++) {
      fl(i, j) = f(i, j);
      fu(i, j) = f(i + noNodes, j);
      bl(i, j) = b(i, j);
      bu(i, j) = b(i, j + noNodes);
    }
  }
}

} // namespace mortar
} // namespace sbp
} // namespace maia

#endif // define DGMORTAR_H_
