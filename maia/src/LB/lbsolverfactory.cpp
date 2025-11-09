// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbsolverfactory.h"

#include <sstream>
#include "IO/context.h"
#include "UTIL/functions.h"
#include "lbsolverdxqy.cpp" // yes, .cpp is correct here
#include "lbsyseqn.h"

namespace maia::lb {

template <MInt nDim>
std::unique_ptr<Solver> LbSolverFactory<nDim>::create(const MInt solverId, maia::grid::Proxy<nDim>& gridProxy,
                                                      Geometry<nDim>& geometry, const MPI_Comm comm) {
  //---read properties----------------------------------------------------------
  /*! \page propertyPage1
    \section noDistributions
    <code>MInt LbSolver::m_noDistributions</code>\n
    default = <code>9</code>\n \n
    Set the number of distributions to use in the Lattice Boltzmann solver in 3D, i.e. the
    discretization model. \n \n
    Possible values are:
    <ul>
      <li>See the possible values of enum "MAIALbModelType"</li>
    </ul>
    Keywords: <i>LATTICE_BOLTZMANN, DISTRIBUTIONS, DISCRETIZATION</i>
  */
  MInt noDistributions = 9;
  noDistributions = Context::getSolverProperty<MInt>("noDistributions", solverId, AT_, &noDistributions);

  // << temporary
  const MString solverMethod = Context::getSolverProperty<MString>("solverMethod", solverId, AT_);
  std::set<MInt> compressibleMethods = {
      MAIA_LATTICE_CUMULANT,        MAIA_LATTICE_BGK_THERMAL,           MAIA_LATTICE_BGK_INNERENERGY,
      MAIA_LATTICE_BGK_TOTALENERGY, MAIA_LATTICE_BGK_THERMAL_TRANSPORT, MAIA_LATTICE_BGKC};
  MBool isCompressible = compressibleMethods.count(string2enum(solverMethod));
  // temporary >>
  isCompressible = Context::getSolverProperty<MBool>("isCompressible", solverId, AT_, &isCompressible);

  //---create solver and return-------------------------------------------------
  if(isCompressible) {
    //
    if constexpr(nDim == 2) {
      switch(noDistributions) {
        case 9: {
          return std::make_unique<LbSolverDxQy<2, 9, LbSysEqnCompressible<2, 9>>>(solverId, 9, gridProxy, geometry,
                                                                                  gridProxy.mpiComm());
          break;
        }
        default: {
          std::stringstream ss;
          ss << "Error: D" << nDim << "Q" << noDistributions << " is unsupported !" << std::endl;
          mTerm(1, ss.str());
        }
      }
    } else if constexpr(nDim == 3) {
      switch(noDistributions) {
        case 19: {
          return std::make_unique<LbSolverDxQy<3, 19, LbSysEqnCompressible<3, 19>>>(solverId, 19, gridProxy, geometry,
                                                                                    comm);
          break;
        }
        case 27: {
          return std::make_unique<LbSolverDxQy<3, 27, LbSysEqnCompressible<3, 27>>>(solverId, 27, gridProxy, geometry,
                                                                                    comm);
          break;
        }
        default: {
          std::stringstream ss;
          ss << "Error: D" << nDim << "Q" << noDistributions << " is unsupported !" << std::endl;
          mTerm(1, ss.str());
        }
      }
    }
  } else {
    if constexpr(nDim == 2) {
      switch(noDistributions) {
        case 9: {
          return std::make_unique<LbSolverDxQy<2, 9, LbSysEqnIncompressible<2, 9>>>(solverId, 9, gridProxy, geometry,
                                                                                    gridProxy.mpiComm());
        }
        default: {
          std::stringstream ss;
          ss << "Error: D" << nDim << "Q" << noDistributions << " is unsupported !" << std::endl;
          mTerm(1, ss.str());
          break;
        }
      }
    } else if constexpr(nDim == 3) {
      switch(noDistributions) {
        case 19: {
          return std::make_unique<LbSolverDxQy<nDim, 19, LbSysEqnIncompressible<3, 19>>>(solverId, 19, gridProxy,
                                                                                         geometry, gridProxy.mpiComm());
          break;
        }
        case 27: {
          return std::make_unique<LbSolverDxQy<nDim, 27, LbSysEqnIncompressible<3, 27>>>(solverId, 27, gridProxy,
                                                                                         geometry, gridProxy.mpiComm());
          break;
        }
        default: {
          std::stringstream ss;
          ss << "Error: D" << nDim << "Q" << noDistributions << " is unsupported !" << std::endl;
          mTerm(1, ss.str());
          break;
        }
      }
    }
  }
}

// Explicit instantiation
template class LbSolverFactory<2>;
template class LbSolverFactory<3>;

} // namespace maia::lb
