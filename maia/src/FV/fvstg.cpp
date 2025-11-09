// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstg.h"

using namespace std;

template <class base_iterator, class SolverType>
typename nDim_iterator_t<base_iterator, SolverType>::value_type
nDim_iterator_t<base_iterator, SolverType>::getNghbr(MInt dir) const {
  return p->m_solver->a_reconstructionNeighborId(p->m_stgToCellId[*(*this)], dir);
}


template <class base_iterator, class SolverType>
typename nDim_iterator_t<base_iterator, SolverType>::value_type
nDim_iterator_t<base_iterator, SolverType>::getCellId() const {
  return p->m_stgToCellId[*(*this)];
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
MSTG<nDim, SolverTypeR, SolverTypeL>::MSTG(FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>* solver,
                                           MInt bcId,
                                           const MPI_Comm commStg,
                                           MInt* sortedBndryCellIds,
                                           MInt noStgLCells,
                                           MInt stgFaceNormalDir,
                                           MInt stgDir,
                                           MInt stgWallNormalDir,
                                           MInt wallDir,
                                           MBool cutOff)

  : m_solver(solver),
    m_solverId(solver->m_solverId),
    m_bcId(bcId),
    m_sutherlandConstant(solver->m_sutherlandConstant),
    m_sutherlandPlusOne(solver->m_sutherlandPlusOne),
    m_commStg(commStg),
    m_noStgLCells(noStgLCells),
    m_stgFaceNormalDir(stgFaceNormalDir),
    m_stgWallNormalDir(stgWallNormalDir),
    m_stgDir(stgDir),
    m_wallDir(wallDir),
    m_cutOff(cutOff),
    a(new Accessor(sortedBndryCellIds, noStgLCells, m_solver, commStg, cutOff)) {
  TRACE();

  IF_CONSTEXPR(nDim == 2) mTerm(1, AT_, "Only implemented for 3D!");

  for(MInt i = 0; i < nDim; i++) {
    if(i != m_stgDir && i != m_wallDir) {
      m_periodicDir = i;
    }
  }

  m_stgLocal = (noStgLCells > 0);

  setSTGProperties();

  // Allocate m_stgVariables
  mAlloc(m_stgLVariables, STG::noStgVars, a->sizeStg(), "m_stgLVariables", F0, AT_);

  if(m_newStgMethod) {
    mAlloc(m_stgEddieCoverage, 2, a->sizeStg(), "m_stgEddieCoverage", F0, AT_);
  }
}


/*
 *    random number generator
 *    @author Marian Albers
 *    @date 01.01.1010
 *
 */
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
MFloat MSTG<nDim, SolverTypeR, SolverTypeL>::generate_rand() {
  std::uniform_real_distribution<MFloat> randRealDistNumber(0, 1);
  return randRealDistNumber(randomEddyPosGenerator());
  //  return rand() / double(RAND_MAX);
}

/**
 *  Generate weighted random numbers, for a distribution with
 *  a higher probability towards the higher bound pass a number < 1.0
 *  for a higher probability towards lower bound pass a number > 1.0
 *  (see tfs)
 * @author Marian Albers
 * @date 01.01.1010
 *
 */
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
MFloat MSTG<nDim, SolverTypeR, SolverTypeL>::generate_rand_weighted() {
  std::uniform_real_distribution<MFloat> randRealDistNumber(0, 1);
  return pow(randRealDistNumber(randomEddyPosGenerator()), m_stgEddieDistribution);
  // return pow((rand() / double(RAND_MAX)), m_stgEddieDistribution);
}

template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
MFloat MSTG<nDim, SolverTypeR, SolverTypeL>::generate_rand_normalDist() {
  return 15.0 * pow((rand() / double(RAND_MAX)) - 0.35, F3) + F1;
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
MFloat MSTG<nDim, SolverTypeR, SolverTypeL>::get_angle(MFloat y, MFloat z) {
  MFloat angle = atan(z / y);
  if(y < 0) {
    if(z >= 0) {
      angle += PI;
    } else {
      angle -= PI;
    }
  }
  return angle;
}

template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::setSTGProperties() {
  TRACE();
  m_stgInitialStartup = false;
  m_stgRootRank = false;

  if(m_stgLocal) {
    MPI_Comm_rank(m_commStg, &m_commStgMyRank);
  } else {
    m_commStgMyRank = -1;
  }

  //
  std::string prop_name;

  // default error message
  std::stringstream errorMsg;
  errorMsg << ": This property is not allowed to exist without an assignment to a BC" << std::endl;

  // STG seed
  prop_name = "bc" + std::to_string(m_bcId) + "STGSeed";
  if(Context::propertyExists(prop_name, m_solverId)) {
    m_randomEddySeed = Context::getSolverProperty<MInt>(prop_name, m_solverId, AT_, &m_randomEddySeed);
    m_PRNGEddy.seed(m_randomEddySeed);
  }

  m_zonal = m_solver->m_zonal;

  prop_name = "bc" + std::to_string(m_bcId) + "preliminary";
  m_preliminary = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_preliminary);

  prop_name = "bc" + std::to_string(m_bcId) + "preliminaryRans2D";
  m_preliminaryRans2D = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_preliminary);

  if(m_preliminaryRans2D) m_preliminary = true;

  prop_name = "bc" + std::to_string(m_bcId) + "newStgMethod";
  m_newStgMethod = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_newStgMethod);

  prop_name = "bc" + std::to_string(m_bcId) + "stgCylinderTransformation";
  m_cylinderTransformation = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_cylinderTransformation);

  // isotropicTurbulence
  if(Context::propertyExists("isotropicTurbulence", m_solverId)) {
    m_isotropicTurbulence =
        Context::getSolverProperty<MBool>("isotropicTurbulence", m_solverId, AT_, &m_isotropicTurbulence);
  }

  // freeStreamTurbulence
  m_uuFS = F0;
  m_vvFS = F0;
  m_wwFS = F0;
  m_SijSijFS = F0;
  m_BLEddieFraction = F1;
  m_freeStreamTurbulence = false;
  if(Context::propertyExists("freeStreamTurbulence", m_solverId)) {
    m_freeStreamTurbulence =
        Context::getSolverProperty<MBool>("freeStreamTurbulence", m_solverId, AT_, &m_freeStreamTurbulence);
    m_uuFS = Context::getSolverProperty<MFloat>("uuFS", m_solverId, AT_, &m_uuFS);
    m_vvFS = Context::getSolverProperty<MFloat>("vvFS", m_solverId, AT_, &m_vvFS);
    m_wwFS = Context::getSolverProperty<MFloat>("wwFS", m_solverId, AT_, &m_wwFS);
    m_SijSijFS = Context::getSolverProperty<MFloat>("SijSijFS", m_solverId, AT_, &m_SijSijFS);
    m_BLEddieFraction = Context::getSolverProperty<MFloat>("BLEddieFraction", m_solverId, AT_, &m_BLEddieFraction);
  }

  /*! \page propertyPage1
    \section stgSubSup
    <code>MFloat FvStructuredSolver::m_stgSubSup </code>\n
    default = <code>0</code>\n \n
    Use mixed subsonics/subsonic formulation\n
    of the STG boundary.\n
    possible values are:
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>STG, STRUCTURED</i>
  */
  if(Context::propertyExists("stgSubSup", m_solverId)) mTerm(1, AT_, "stgSubSup" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgSubSup";
  m_stgSubSup = false;
  if(Context::propertyExists(prop_name, m_solverId)) {
    m_stgSubSup = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_stgSubSup);
  }

  /*! \page propertyPage1
    \section stgSupersonic
    <code>MFloat FvStructuredSolver::m_stgSupersonic </code>\n
    default = <code>0</code>\n \n
    Use supersonic STG boundary formulation.\n
    possible values are:
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>STG, STRUCTURED</i>
  */
  if(Context::propertyExists("stgSupersonic", m_solverId)) mTerm(1, AT_, "stgSupersonic" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgSupersonic";
  m_stgSupersonic = false;
  if(Context::propertyExists(prop_name, m_solverId)) {
    m_stgSupersonic = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_stgSupersonic);

    if(m_stgSupersonic && m_stgSubSup) {
      m_stgSubSup = false;
      if(m_solver->domainId() == 0) {
        std::cout << "WARNING: You activated conflicting properties stgSubSup "
                  << "and stgSupersonic, thus only the pure supersonic formulation will be used. "
                  << "Switch off stgSupersonic to get the mixed formulation" << std::endl;
      }
    }
  }

  /*! \page propertyPage1
    \section BLT1
    <code>MFloat FvStructuredSolver::m_stgBLT1 </code>\n
    default = <code>1.0</code>\n \n
    Defines the size of the STG virtual box\n
    in the x-direction as a fraction of the\n
    delta0 specified.\n
    possible values are:
    <ul>
    <li>Floating point > 0.0</li>
    </ul>
    Keywords: <i>STG, STRUCTURED</i>
  */
  if(Context::propertyExists("BLT1", m_solverId)) mTerm(1, AT_, "BLT1" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "BLT1";
  m_stgBLT1 = 1.0;
  m_stgBLT1 = Context::getSolverProperty<MFloat>(prop_name, m_solverId, AT_, &m_stgBLT1);

  /*! \page propertyPage1
    \section BLT2
    <code>MFloat FvStructuredSolver::m_stgBLT2 </code>\n
    default = <code>2</code>\n \n
    Defines the size of the STG virtual box\n
    in the y-direction as a fraction of the\n
    delta0 specified.\n
    possible values are:
    <ul>
    <li>Floating point > 0.0</li>
    </ul>
    Keywords: <i>STG, STRUCTURED</i>
  */
  if(Context::propertyExists("BLT2", m_solverId)) mTerm(1, AT_, "BLT2" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "BLT2";
  m_stgBLT2 = 1.1;
  m_stgBLT2 = Context::getSolverProperty<MFloat>(prop_name, m_solverId, AT_, &m_stgBLT2);

  /*! \page propertyPage1
    \section BLT3
    <code>MFloat FvStructuredSolver::m_stgBLT3 </code>\n
    default = <code>1.1</code>\n \n
    Defines the size of the STG virtual box\n
    in the z-direction as a fraction of the\n
    delta0 specified.\n
    possible values are:
    <ul>
    <li>Floating point > 0.0</li>
    </ul>
    Keywords: <i>STG, STRUCTURED</i>
  */
  if(Context::propertyExists("BLT3", m_solverId)) mTerm(1, AT_, "BLT3" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "BLT3";
  m_stgBLT3 = 1.0;
  m_stgBLT3 = Context::getSolverProperty<MFloat>(prop_name, m_solverId, AT_, &m_stgBLT3);

  mAlloc(m_stgLengthFactors, 3, "m_stgLengthFactors", F0, AT_);
  m_stgLengthFactors[0] = 1.0;
  m_stgLengthFactors[1] = 0.6;
  m_stgLengthFactors[2] = 1.5;

  /*! \page propertyPage1
    \section stgLengthFactors
    <code>MFloat FvStructuredSolver::m_stgLengthFactors </code>\n
    default = <code>1.0, 0.6, 1.5</code>\n \n
    The factor to scale the length scales\n
    in each coordinate direction with. For higher\n
    Reynolds number the values [1.0, 0.5, 1.4] \n
    produce better results.\n
    possible values are:
    <ul>
    <li>Float<3> > 0.0</li>
    </ul>
    Keywords: <i>STG, STRUCTURED</i>
  */
  if(Context::propertyExists("stgLengthFactors", m_solverId)) mTerm(1, AT_, "stgLengthFactors" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgLengthFactors";
  if(Context::propertyExists(prop_name, m_solverId)) {
    for(MInt i = 0; i < 3; i++) {
      m_stgLengthFactors[i] = Context::getSolverProperty<MFloat>(prop_name, m_solverId, AT_, &m_stgLengthFactors[i], i);
    }
  }

  mAlloc(m_stgRSTFactors, 3, "m_stgRSTFactors", F0, AT_);
  m_stgRSTFactors[0] = 1.0;
  m_stgRSTFactors[1] = 0.4;
  m_stgRSTFactors[2] = 0.5;

  /*! \page propertyPage1
    \section stgRSTFactors
    <code>MFloat FvStructuredSolver::m_stgRSTFactors </code>\n
    default = <code>1.0, 0.4, 0.5</code>\n \n
    The factor to scale the length scales\n
    in each coordinate direction with. For higher\n
    Reynolds number the values [1.0, 0.4, 0.5] \n
    produce better results.\n
    possible values are:
    <ul>
    <li>Float<3> > 0.0</li>
    </ul>
    Keywords: <i>STG, STRUCTURED</i>
  */
  if(Context::propertyExists("stgRSTFactors", m_solverId)) mTerm(1, AT_, "stgRSTFactors" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgRSTFactors";
  if(Context::propertyExists(prop_name, m_solverId)) {
    for(MInt i = 0; i < 3; i++) {
      m_stgRSTFactors[i] = Context::getSolverProperty<MFloat>(prop_name, m_solverId, AT_, &m_stgRSTFactors[i], i);
    }
  }

  /*! \page propertyPage1
    \section stgMaxNoEddies
    <code>MFloat FvStructuredSolver::m_stgMaxNoEddies </code>\n
    default = <code>200</code>\n \n
    Number of Eddies in the STG virtual box.\n
    possible values are:
    <ul>
    <li>Integer > 0</li>
    </ul>
    Keywords: <i>STG, STRUCTURED</i>
  */
  if(Context::propertyExists("stgMaxNoEddies", m_solverId)) mTerm(1, AT_, "stgMaxNoEddies" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgMaxNoEddies";
  m_stgMaxNoEddies = 200;
  m_stgMaxNoEddies = Context::getSolverProperty<MInt>(prop_name, m_solverId, AT_, &m_stgMaxNoEddies);

  /*! \page propertyPage1
    \section stgExple
    <code>MFloat FvStructuredSolver::m_stgExple </code>\n
    default = <code>0.5</code>\n \n
    Exponent of the STG LengthScale law.\n
    possible values are:
    <ul>
    <li>Floating point > 0.0</li>
    </ul>
Keywords: <i>STG, STRUCTURED</i>
*/
  if(Context::propertyExists("stgExple", m_solverId)) mTerm(1, AT_, "stgExple" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgExple";
  m_stgExple = 0.5;
  m_stgExple = Context::getSolverProperty<MFloat>(prop_name, m_solverId, AT_, &m_stgExple);

  /*! \page propertyPage1
    \section stgEddieDistribution
    <code>MFloat FvStructuredSolver::m_stgEddieDistribution </code>\n
    default = <code>1.0</code>\n \n
    Shift die eddie distribution more to the wall\n
    or boundary layer edge.\n
    possible values are:
    <ul>
    <li>Floating point > 0.0</li>
    </ul>
Keywords: <i>STG, STRUCTURED</i>
*/
  if(Context::propertyExists("stgEddieDistribution", m_solverId))
    mTerm(1, AT_, "stgEddieDistribution" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgEddieDistribution";
  m_stgEddieDistribution = 1.0;
  if(Context::propertyExists(prop_name, m_solverId)) {
    m_stgEddieDistribution = Context::getSolverProperty<MFloat>(prop_name, m_solverId, AT_, &m_stgEddieDistribution);
  }

  /*! \page propertyPage1
    \section stgCreateNewEddies
    <code>MFloat FvStructuredSolver::m_stgCreateNewEddies </code>\n
    default = <code>0</code>\n \n
    Enforces the creation of all new eddies in STG virtual box\n
    or boundary layer edge.\n
    possible values are:
    <ul>
    <li>true/false</li>
    </ul>
Keywords: <i>STG, STRUCTURED</i>
*/
  if(Context::propertyExists("stgCreateNewEddies", m_solverId)) mTerm(1, AT_, "stgCreateNewEddies" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgCreateNewEddies";
  m_stgCreateNewEddies = false;
  if(Context::propertyExists(prop_name, m_solverId)) {
    m_stgCreateNewEddies = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_stgCreateNewEddies);
  }

  /*! \page propertyPage1
    \section stgInitialStartup
    <code>MFloat FvStructuredSolver::m_stgInitialStartup </code>\n
    default = <code>0</code>\n \n
    Initialize STG Method at Startup\n
    possible values are:
    <ul>
    <li>true/false</li>
    </ul>
Keywords: <i>STG, STRUCTURED</i>
*/
  if(Context::propertyExists("stgInitialStartup", m_solverId)) mTerm(1, AT_, "stgInitialStartup" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgInitialStartup";
  m_stgInitialStartup = false;
  m_stgInitialStartup = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_stgInitialStartup);

  /*! \page propertyPage1
    \section stgEddieLengthScales
    <code>MFloat FvStructuredSolver::m_stgEddieLengthScales </code>\n
    default = <code>0</code>\n \n
    Connect length scales to eddies, not cells.\n
    possible values are:
    <ul>
    <li>true/false</li>
    </ul>
Keywords: <i>STG, STRUCTURED</i>
*/
  if(Context::propertyExists("stgEddieLengthScales", m_solverId))
    mTerm(1, AT_, "stgEddieLengthScales" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgEddieLengthScales";
  m_stgEddieLengthScales = false;
  if(Context::propertyExists(prop_name, m_solverId)) {
    m_stgEddieLengthScales = Context::getSolverProperty<MBool>(prop_name, m_solverId, AT_, &m_stgEddieLengthScales);
  }

  /*! \page propertyPage1
    \section stgShapeFunction
    <code>MFloat FvStructuredSolver::m_stgFunction </code>\n
    default = <code>4</code>\n \n
    Shape function to be used in STG method.\n
    possible values are:
    <ul>
    <li>Integer >= 0</li>
    </ul>
Keywords: <i>STG, STRUCTURED</i>
*/
  if(Context::propertyExists("stgShapeFunction", m_solverId)) mTerm(1, AT_, "stgShapeFunction" + errorMsg.str());
  prop_name = "bc" + std::to_string(m_bcId) + "stgShapeFunction";
  m_stgShapeFunction = 4;
  if(Context::propertyExists(prop_name, m_solverId)) {
    m_stgShapeFunction = Context::getSolverProperty<MInt>(prop_name, m_solverId, AT_, &m_stgShapeFunction);
  }

  // TODO
  /*  if(m_stgInitialStartup) {
      //activate the nut FQ field
      FQ->neededFQVariables[FQ->NU_T] = 1;
    }*/

  m_stgNoEddieProperties = 6;
  if(m_stgEddieLengthScales) {
    m_stgNoEddieProperties = 9;
  }
  mAlloc(m_stgEddies, m_stgMaxNoEddies, m_stgNoEddieProperties, "m_solver->m_stgEddies", -F1, AT_);

  // JANNIK:new method
  mAlloc(m_stgEddyStrength, m_stgMaxNoEddies, "m_solver->m_stgEddyStrength", F1, AT_);

  m_log << "===========================================================" << std::endl
        << "                    STG PROPERTIES " << std::endl
        << "===========================================================" << std::endl
        << "Initial Start: " << m_stgInitialStartup << std::endl
        << "SubSup (Mixed subsonic/supersonic bc): " << m_stgSubSup << std::endl
        << "Supersonic BC: " << m_stgSupersonic << std::endl
        << "BLT 1,2,3: " << m_stgBLT1 << ", " << m_stgBLT2 << ", " << m_stgBLT3 << std::endl
        << "Length factors: " << m_stgLengthFactors[0] << ", " << m_stgLengthFactors[1] << ", " << m_stgLengthFactors[2]
        << std::endl
        << "Number of eddies: " << m_stgMaxNoEddies << std::endl
        << "Length scale exponent: " << m_stgExple << std::endl
        << "Eddie distribution: " << m_stgEddieDistribution << std::endl
        << "Create new eddies: " << m_stgCreateNewEddies << std::endl
        << "Eddie lengthscales: " << m_stgEddieLengthScales << std::endl
        << "Shape function: " << m_stgShapeFunction << std::endl
        << "Number of eddie properties: " << m_stgNoEddieProperties << std::endl
        << "Number of stg variables: " << STG::noStgVars << std::endl
        << "===========================================================" << std::endl;
}


/**
 *      Synthetic Turbulence Generation
 *
 */
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::init(MInt commStgRoot) {
  TRACE();

  MPI_Comm_rank(m_commStg, &m_stgMyRank);
  m_commStgRoot = commStgRoot;

  mAlloc(m_stgMaxVel, m_nDim, "m_stgMaxVel", 0.0, AT_);
  mAlloc(m_stgVbStart, m_nDim, "m_stgVbStart", 0.0, AT_);
  mAlloc(m_stgVbEnd, m_nDim, "m_stgVbEnd", 0.0, AT_);

#ifndef NDEBUG
  if(m_stgMyRank == m_commStgRoot) {
    std::cout << globalTimeStep << " Initializing BC " + std::to_string(m_bcId) << std::endl;
  }
#endif

  readRANSProfileStg();

  if(m_solver->m_resetInitialCondition) {
    m_stgDelta99Inflow = 0.0;
    m_initialRange = true;
  } else {
    getBoundaryLayerThickness();
  }

  // Compute size of the Box, take inflow as middle x-coordinate
  mAlloc(m_inflowStart, m_nDim, "m_inflowStart", AT_);
  mAlloc(m_inflowEnd, m_nDim, "m_inflowEnd", AT_);
  getInflowStartEnd(m_inflowStart, m_inflowEnd);
  setVb(m_inflowStart, m_inflowEnd);

  if(m_newStgMethod) {
    // determinePeriodicCells();
  }

  // TODO: check if this is the right place
  if((globalTimeStep == m_solver->m_restartTimeStep && globalTimeStep > m_solver->m_stgStartTimeStep
      && !m_stgCreateNewEddies)
     || (m_solver->m_wasAdapted) || (m_solver->m_wasBalancedZonal)) {
    loadStg();
  } else {
    // // create new eddies if it's an initial start or the createNewEddie flag is set
    MFloatScratchSpace bcast_eddies(m_stgMaxNoEddies * 6, AT_, "bcast_eddies");
    MFloatScratchSpace bcast_eddyStrength(m_stgMaxNoEddies, AT_, "bcast_eddyStrength");

    if(m_stgMyRank == m_commStgRoot) {
#ifndef NDEBUG
      std::cerr << "-------Creating new Eddies inside Virtual Box------" << std::endl;
#endif
      MFloat xk[m_nDim];
      MFloat epsi[m_nDim];
      for(MInt n = 0; n < m_stgMaxNoEddies; n++) {
        generateNewEddies(xk, epsi);

        bcast_eddies[n + m_stgMaxNoEddies * 0] = xk[0];
        bcast_eddies[n + m_stgMaxNoEddies * 1] = xk[1];
        bcast_eddies[n + m_stgMaxNoEddies * 2] = xk[2];
        bcast_eddies[n + m_stgMaxNoEddies * 3] = epsi[0];
        bcast_eddies[n + m_stgMaxNoEddies * 4] = epsi[1];
        bcast_eddies[n + m_stgMaxNoEddies * 5] = epsi[2];

        bcast_eddyStrength[n] = F1;
        if(m_newStgMethod) {
          bcast_eddyStrength[n] = generate_rand_normalDist();
        }
      }
    }

    // Broadcast the new/updated eddies to all relevant processes
    MPI_Bcast(bcast_eddies.begin(), 6 * m_stgMaxNoEddies, MPI_DOUBLE, m_commStgRoot, m_commStg, AT_,
              "bcast_eddies.begin()");
    MPI_Bcast(bcast_eddyStrength.begin(), m_stgMaxNoEddies, MPI_DOUBLE, m_commStgRoot, m_commStg, AT_,
              "bcast_eddyStrength.begin()");
    // Copy data into m_FQeddie vector
    for(MInt n = 0; n < m_stgMaxNoEddies; n++) {
      m_stgEddies[n][0] = bcast_eddies[n + m_stgMaxNoEddies * 0];
      m_stgEddies[n][1] = bcast_eddies[n + m_stgMaxNoEddies * 1];
      m_stgEddies[n][2] = bcast_eddies[n + m_stgMaxNoEddies * 2];
      m_stgEddies[n][3] = bcast_eddies[n + m_stgMaxNoEddies * 3];
      m_stgEddies[n][4] = bcast_eddies[n + m_stgMaxNoEddies * 4];
      m_stgEddies[n][5] = bcast_eddies[n + m_stgMaxNoEddies * 5];

      m_stgEddyStrength[n] = bcast_eddyStrength[n];
    }
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::determinePeriodicCells() {
  TRACE();

  MInt noStgBcCells = 0;
  MInt noBcStgLocations = 0;
  MFloat periodicL = 0;
  MBool first = true;
  m_stgWallNormalLocations.clear();

  // ======================================================
  // 1) Determine local locations in wall-normal direction for spanwise average
  // ======================================================
  for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
    const MInt cellId = a->getCellId(it);

    MFloat halfCellLength = m_solver->grid().halfCellLength(cellId);
    if(first) {
      periodicL = m_solver->a_coordinate(cellId, m_periodicDir) - F1B3 * halfCellLength;
      first = false;
    }
    if(abs(m_solver->a_coordinate(cellId, m_periodicDir) + halfCellLength - eps - periodicL) < halfCellLength) {
      m_stgWallNormalLocations.push_back(m_solver->a_coordinate(cellId, m_wallDir));
      noBcStgLocations++;
    }
    noStgBcCells++;
  }

  // ======================================================
  // 2) Create communicator
  // ======================================================
  MInt comm_size;
  MPI_Comm_size(m_commStg, &comm_size);

  // ======================================================
  // 3) Determine global locations in wall-normal direction for spanwise average
  // ======================================================
  MInt globalNoBcStgLocations = 0;
  MPI_Allreduce(&noBcStgLocations, &globalNoBcStgLocations, 1, MPI_INT, MPI_SUM, m_commStg, AT_, "noBcStgLocations",
                "globalNoBcStgLocations");

  ScratchSpace<MFloat> globalBcStgLocations(globalNoBcStgLocations, "globalBcStgLocations", FUN_);
  ScratchSpace<MInt> recvbuf(comm_size, "recvbuf", FUN_);
  ScratchSpace<MInt> displs(comm_size, "displspos", FUN_);

  recvbuf.fill(0);

  MPI_Gather(&noBcStgLocations, 1, MPI_INT, &recvbuf[0], 1, MPI_INT, 0, m_commStg, AT_, "noBcStgLocations", "recvbuf");

  if(m_stgMyRank == m_commStgRoot) {
    MInt offset = 0;
    for(MInt dom = 0; dom < comm_size; dom++) {
      displs[dom] = offset;
      offset += recvbuf[dom];
    }
  }

  MPI_Gatherv(&m_stgWallNormalLocations[0], noBcStgLocations, MPI_DOUBLE, &globalBcStgLocations[0],
              &recvbuf[m_stgMyRank], &displs[m_stgMyRank], MPI_DOUBLE, 0, m_commStg, AT_, "m_StgwallNormalLocations",
              "globalBcStgLocations");

  MPI_Bcast(&globalBcStgLocations[0], globalNoBcStgLocations, MPI_DOUBLE, 0, m_commStg, AT_, "globalBcStgLocations");

  m_stgGlobalWallNormalLocations.clear();

  for(MInt i = 0; i < globalNoBcStgLocations; i++) {
    MFloat L = globalBcStgLocations[i];
    if(std::find(m_stgGlobalWallNormalLocations.begin(), m_stgGlobalWallNormalLocations.end(), L)
       == m_stgGlobalWallNormalLocations.end()) {
      m_stgGlobalWallNormalLocations.push_back(L);
    }
  }


  m_stgGlobalNoWallNormalLocations = (MInt)m_stgGlobalWallNormalLocations.size();

  std::sort(m_stgGlobalWallNormalLocations.begin(), m_stgGlobalWallNormalLocations.end());

  // ======================================================
  // 4) Determine local cells in periodic locations of global wall-normal locations
  // ======================================================
  mAlloc(m_stgPeriodicCellId, m_stgGlobalNoWallNormalLocations, "m_stgPeriodicCellIndex", FUN_);
  mAlloc(m_stgGlobalNoPeriodicLocations, m_stgGlobalNoWallNormalLocations, "m_stgGlobalNoPeriodicLocations", 0, FUN_);

  for(MInt i = 0; i < m_stgGlobalNoWallNormalLocations; i++) {
    m_stgPeriodicCellId[i].clear();
  }

  mAlloc(m_stgPeriodicIndex, noStgBcCells, "m_stgPeriodicIndex", -1, FUN_);

  vector<MInt> noPeriodicLocations(m_stgGlobalNoWallNormalLocations, F0);

  for(MInt i = 0; i < m_stgGlobalNoWallNormalLocations; i++) {
    for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
      const MInt cellId = a->getCellId(it);
      const MInt IBC = a->getStgId(it);
      if(abs(m_solver->a_coordinate(cellId, m_wallDir) - m_stgGlobalWallNormalLocations[i]) < eps) {
        m_stgPeriodicIndex[IBC] = i;
        if(!m_solver->a_isHalo(cellId)) {
          m_stgPeriodicCellId[i].push_back(IBC);
          ++noPeriodicLocations[i];
        }
      }
    }
  }

  MPI_Allreduce(&noPeriodicLocations[0], &m_stgGlobalNoPeriodicLocations[0], m_stgGlobalNoWallNormalLocations, MPI_INT,
                MPI_SUM, m_commStg, AT_, "noPeriodicLocations", "m_stgGlobalNoPeriodicLocations");
}

template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::saveStg() {
  TRACE();

  using namespace maia::parallel_io;
  std::stringstream filename;
  filename << m_solver->outputDir() << "stg" << std::to_string(m_bcId) << "RestartNew_" << globalTimeStep
           << ParallelIo::fileExt();

  m_log << "Writing restart file " << filename.str() << " ..." << std::endl;
#ifndef NDEBUG
  if(m_stgMyRank == m_commStgRoot) {
    cerr << "Writing restart file " << filename.str() << " ..." << std::endl;
  }
#endif

  std::stringstream stgPrefix_;
  stgPrefix_ << "stgVar" << m_bcId << "_";
  MString stgPrefix = stgPrefix_.str();

  ParallelIo parallelIo((filename.str()).c_str(), PIO_REPLACE, m_commStg);
  parallelIo.setAttributes(&(m_stgMaxNoEddies), (stgPrefix + "stgMaxNoEddies").c_str(), 1);

  parallelIo.defineArray(PIO_FLOAT, (stgPrefix + "FQeddies").c_str(), m_stgMaxNoEddies * m_stgNoEddieProperties);

  parallelIo.defineArray(PIO_FLOAT, (stgPrefix + "FQeddyStrength").c_str(), m_stgMaxNoEddies);

  parallelIo.setAttribute("FQeddies", "name", (stgPrefix + "FQeddies").c_str());
  parallelIo.setOffset(m_stgMaxNoEddies * m_stgNoEddieProperties, 0);
  parallelIo.writeArray(&(m_stgEddies[0][0]), (stgPrefix + "FQeddies").c_str());

  parallelIo.setAttribute("FQeddyStrength", "name", (stgPrefix + "FQeddyStrength").c_str());
  parallelIo.setOffset(m_stgMaxNoEddies, 0);
  parallelIo.writeArray(&(m_stgEddyStrength[0]), (stgPrefix + "FQeddyStrength").c_str());
}

template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::loadStg() {
  TRACE();

  using namespace maia::parallel_io;
  std::stringstream filename;
  filename << m_solver->restartDir() << "stg" << std::to_string(m_bcId) << "RestartNew_" << globalTimeStep
           << ParallelIo::fileExt();

  m_log << "Loading restart file " << filename.str() << " ..." << std::endl;
#ifndef NDEBUG
  if(m_stgMyRank == m_commStgRoot) {
    cerr << "Loading restart file " << filename.str() << " ..." << std::endl;
  }
#endif

  std::stringstream stgPrefix_;
  stgPrefix_ << "stgVar" << m_bcId << "_";
  MString stgPrefix = stgPrefix_.str();

  // All stg ranks read m_stgMaxNoEddies & m_stgEddies
  {
    ParallelIo parallelIo((filename.str()).c_str(), PIO_READ, m_commStg);
    MInt stgMaxNoEddies;
    parallelIo.getAttribute(&stgMaxNoEddies, (stgPrefix + "stgMaxNoEddies").c_str());
    if(stgMaxNoEddies != m_stgMaxNoEddies && !m_stgCreateNewEddies)
      TERMM(-1, "bc" + std::to_string(m_bcId) + ": Number of eddies has changed!");
    ParallelIo::size_type stgEddies_ = parallelIo.getArraySize((stgPrefix + "FQeddies").c_str());
    if(stgEddies_ != stgMaxNoEddies * m_stgNoEddieProperties) TERMM(-1, "");
    ParallelIo::size_type stgEddyStrength_ = parallelIo.getArraySize((stgPrefix + "FQeddyStrength").c_str());
    if(stgEddyStrength_ != stgMaxNoEddies) TERMM(-1, "");

    parallelIo.setOffset(stgEddies_, 0);
    parallelIo.readArray(&(m_stgEddies[0][0]), (stgPrefix + "FQeddies").c_str()); // TODO: check this

    parallelIo.setOffset(stgEddyStrength_, 0);
    parallelIo.readArray(&(m_stgEddyStrength[0]), (stgPrefix + "FQeddyStrength").c_str()); // TODO: check this
  }
}

/** \brief Reformulated Synthetic Turbulence Generation
 * Synthetic Turbulence Generation Method
 * Computes fluctuations from a given averaged
 * velocity and nu_t profile and adds them at
 * the inflow of the domain.
 * Same computations as STG in TFS by Benedikt Roidl
 * /author Jannik Borgelt
 */
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::bc7909() {
  if((globalTimeStep != m_solver->m_restartTimeStep) && !m_solver->m_wasAdapted
     && globalTimeStep % m_solver->m_restartInterval == 0 && m_solver->m_RKStep == 0) {
    saveStg();
  }

  if(m_solver->m_RKStep == 0) {
    if(m_zonal
       && (globalTimeStep % m_solver->m_zonalTransferInterval == 0 || globalTimeStep == m_solver->m_restartTimeStep
           || m_solver->m_wasAdapted)) {
      readRANSProfileStg();

      getBoundaryLayerThickness();
      getInflowStartEnd(m_inflowStart, m_inflowEnd);
      setVb(m_inflowStart, m_inflowEnd);

      calcStrainRate();
      calcReynoldsStressConvVelLengthScale();
    }

    if(!m_zonal && (globalTimeStep == m_solver->m_restartTimeStep || globalTimeStep <= 1)) {
      calcStrainRate();
      calcReynoldsStressConvVelLengthScale();
    }

    /*********************** J *****************/
    // The virtual box part - executed by Master Solver at the inflow
    /*********************** J *****************/

    if(!m_solver->m_wasAdapted) {
      advanceEddies();
    }

#ifndef NDEBUG
    // Summary of synth turb parameters
    // printSTGSummary();
#endif
    /*********************** L *****************/
    // Calculation of the fluctuation induced by all eddies on each cell
    /*********************** L *****************/
    calcTotalFluctuationCholesky();

    if(m_newStgMethod) {
      // calcEddieCoverage();
    }

  } // RKStep end if

  /////////////////////////////////////////////////////////
  ////////////// APPLY TO BC //////////////////////////////
  /////////////////////////////////////////////////////////
  // Now comes the BC stuff that we need to do every RK step
  apply();
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::getBoundaryLayerThickness() {
  TRACE();

  std::stringstream errorMsg;
  std::string prop_name;
  errorMsg << ": This property is not allowed to exist without an assignment to a BC" << std::endl;

  m_stgDelta99Inflow = -1.0;

  MInt factor = (m_stgWallNormalDir % 2 == 0) ? -1 : 1;

  if(!m_zonal) {
    prop_name = "bc" + std::to_string(m_bcId) + "deltaIn";
    if(!Context::propertyExists(prop_name, m_solverId)) {
      mTerm(1, AT_, "deltaIn" + errorMsg.str());
    }
    m_stgDelta99Inflow = Context::getSolverProperty<MFloat>(prop_name, m_solverId, AT_, &m_stgDelta99Inflow);
  } else {
    MFloat maxVel = F0;
    MFloat maxVelPos = std::numeric_limits<MFloat>::max();
    MFloat minPos = std::numeric_limits<MFloat>::max();
    MFloat maxPos = std::numeric_limits<MFloat>::lowest();
    MFloat halfCellLength = F0;
    for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
      const MInt IBC = a->getStgId(it);
      const MInt cellId = a->getCellId(it);
      MFloat y = LES->a_coordinate(cellId, m_wallDir);
      MFloat z = LES->a_coordinate(cellId, m_periodicDir);
      MFloat r = sqrt(y * y + z * z);
      if(m_stgLVariables[m_stgDir][IBC] >= maxVel) {
        maxVel = m_stgLVariables[m_stgDir][IBC];
        if(m_cylinderTransformation) {
          maxVelPos = r;
        } else {
          maxVelPos = m_solver->a_coordinate(cellId, m_wallDir);
        }
      }
      if(m_cylinderTransformation) {
        if(r < minPos) {
          minPos = r;
        }
        if(r > maxPos) {
          maxPos = r;
        }
      } else {
        if(m_solver->a_coordinate(cellId, m_wallDir) < minPos) {
          minPos = m_solver->a_coordinate(cellId, m_wallDir);
          halfCellLength = m_solver->grid().halfCellLength(cellId);
        }
        if(m_solver->a_coordinate(cellId, m_wallDir) > maxPos) {
          maxPos = m_solver->a_coordinate(cellId, m_wallDir);
        }
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &maxVel, 1, MPI_DOUBLE, MPI_MAX, m_commStg, AT_, "maxVel", "1");
    MPI_Allreduce(MPI_IN_PLACE, &maxVelPos, 1, MPI_DOUBLE, MPI_MIN, m_commStg, AT_, "maxVelPos", "1");
    MPI_Allreduce(MPI_IN_PLACE, &minPos, 1, MPI_DOUBLE, MPI_MIN, m_commStg, AT_, "minPos", "1");
    MPI_Allreduce(MPI_IN_PLACE, &maxPos, 1, MPI_DOUBLE, MPI_MAX, m_commStg, AT_, "maxPos", "1");

    if(maxVel < 0.9 * m_solver->m_UInfinity) {
      m_initialRange = true;
      m_stgDelta99Inflow = 2 * halfCellLength;
    }

    if(!m_initialRange) {
      MFloat delta0 = std::numeric_limits<MFloat>::max();
      MFloat delta0_1 = minPos;
      MFloat delta0_2 = maxPos;
      MFloat delta0_u1 = std::numeric_limits<MFloat>::max();
      MFloat delta0_u2 = std::numeric_limits<MFloat>::lowest();
      for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
        const MInt cellId = a->getCellId(it);
        const MInt IBC = a->getStgId(it);
        MFloat pos = m_solver->a_coordinate(cellId, m_wallDir);
        if(m_cylinderTransformation) {
          MFloat y = LES->a_coordinate(cellId, m_wallDir);
          MFloat z = LES->a_coordinate(cellId, m_periodicDir);
          pos = sqrt(y * y + z * z);
        }

        if(m_solver->a_isHalo(cellId)) continue;

        if(factor * m_stgLVariables[m_stgDir][IBC] > factor * 0.99 * maxVel && pos < delta0_2 && pos < maxVelPos) {
          delta0_2 = pos;
          delta0_u2 = m_stgLVariables[m_stgDir][IBC];
        }

        if(factor * m_stgLVariables[m_stgDir][IBC] < factor * 0.99 * maxVel && pos > delta0_1 && pos < maxVelPos) {
          delta0_1 = pos;
          delta0_u1 = m_stgLVariables[m_stgDir][IBC];
        }
      }

      MFloat deltaOffset = (factor == 1) ? minPos : maxPos;
      delta0 =
          factor
          * (delta0_1 + (delta0_2 - delta0_1) / (delta0_u2 - delta0_u1) * (0.99 * maxVel - delta0_u1) - deltaOffset);

      if(delta0 < 1e-16) delta0 = std::numeric_limits<MFloat>::max();

      MPI_Allreduce(&delta0, &m_stgDelta99Inflow, 1, MPI_DOUBLE, MPI_MIN, m_commStg, AT_, "m_stgDelta99Inflow",
                    "delta0");
    }
  }

#ifndef NDEBUG
  if(m_stgMyRank == m_commStgRoot) {
    std::cerr << "---------- Boundary layer thickness (" << globalTimeStep << "," << m_bcId << "): ----------"
              << std::endl;
    std::cerr << m_stgDelta99Inflow << std::endl;
  }
#endif
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::getBoundaryLayerThickness(/*SolverTraits<nDim, MAIA_STRUCTURED>**/) {
  MInt bndryLayerIndex = 0;
  MFloat u99_0 = F0, u99_1 = F0;
  MInt start_[] = {0, a->start(1) + a->m_noGhostLayers, 2};
  MInt end_[] = {1, a->end(1) - a->m_noGhostLayers, 3};
  for(typename Accessor::nDim_citerator it = a->nDim_citerator_begin(start_, end_); it != a->nDim_citerator_end();
      ++it) {
    const MInt IBC = a->getStgId(it);
    const MInt IJMK = a->getNghbrStg(it, 2);
    if(m_stgLVariables[PV->U][IBC] >= 0.99 * LES->UInfinity()
       && m_stgLVariables[PV->U][IJMK] < 0.99 * LES->UInfinity()) {
      u99_0 = m_stgLVariables[PV->U][IJMK];
      u99_1 = m_stgLVariables[PV->U][IBC];
      bndryLayerIndex = it->getijk(1) - 1;
      break;
    }
  }
  MFloat bndryLayerThicknessLocal = 0.0;
  MFloat bndryLayerThicknessGlobal = 0.0;

  if(bndryLayerIndex > 0) {
    bndryLayerThicknessLocal = LES->a_coordinate(a->cellIndex(0, bndryLayerIndex, 2), 1)
                               + (0.99 * LES->UInfinity() - u99_0) / (u99_1 - u99_0)
                                     * (LES->a_coordinate(a->cellIndex(0, bndryLayerIndex + 1, 2), 1)
                                        - LES->a_coordinate(a->cellIndex(0, bndryLayerIndex, 2), 1));
  }
  MPI_Allreduce(&bndryLayerThicknessLocal, &bndryLayerThicknessGlobal, 1, MPI_DOUBLE, MPI_MAX, m_commStg, AT_,
                "bndryLayerThicknessLocal", "bndryLayerThicknessGlobal");
  if(m_stgMyRank == m_commStgRoot) {
    std::cout << "Boundary Layer thickness delta99: " << bndryLayerThicknessGlobal << std::endl;
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::getInflowStartEnd(MFloat* inflowStart, MFloat* inflowEnd) {
  TRACE();

  MFloat inflowStartLocal[] = {std::numeric_limits<MFloat>::max(), std::numeric_limits<MFloat>::max(),
                               std::numeric_limits<MFloat>::max()};
  MFloat inflowEndlocal[] = {-std::numeric_limits<MFloat>::max(), -std::numeric_limits<MFloat>::max(),
                             -std::numeric_limits<MFloat>::max()};
  for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
    MInt cellId = a->getCellId(it);
    if(m_solver->a_isHalo(cellId)) continue;

    if(m_cylinderTransformation) {
      MFloat y = LES->a_coordinate(cellId, m_wallDir);
      MFloat z = LES->a_coordinate(cellId, m_periodicDir);
      MFloat alpha = get_angle(y, z);
      MFloat r = sqrt(y * y + z * z);
      if(inflowStartLocal[0] > LES->a_coordinate(cellId, 0)) inflowStartLocal[0] = LES->a_coordinate(cellId, 0);
      if(inflowEndlocal[0] < LES->a_coordinate(cellId, 0)) inflowEndlocal[0] = LES->a_coordinate(cellId, 0);
      if(inflowStartLocal[1] > r) inflowStartLocal[1] = r;
      if(inflowEndlocal[1] < r) inflowEndlocal[1] = r;
      if(inflowStartLocal[2] > alpha) inflowStartLocal[2] = alpha;
      if(inflowEndlocal[2] < alpha) inflowEndlocal[2] = alpha;

    } else {
      for(MInt d = 0; d < nDim; ++d) {
        if(inflowStartLocal[d] > LES->a_coordinate(cellId, d)) inflowStartLocal[d] = LES->a_coordinate(cellId, d);
        if(inflowEndlocal[d] < LES->a_coordinate(cellId, d)) inflowEndlocal[d] = LES->a_coordinate(cellId, d);
      }
    }
  }

  std::fill(inflowStart, inflowStart + nDim, 9999999);
  std::fill(inflowEnd, inflowEnd + nDim, -9999999);
  MPI_Allreduce(inflowStartLocal, inflowStart, nDim, MPI_DOUBLE, MPI_MIN, m_commStg, AT_, "inflowStartLocal",
                "inflowStart");
  MPI_Allreduce(inflowEndlocal, inflowEnd, nDim, MPI_DOUBLE, MPI_MAX, m_commStg, AT_, "inflowEndlocal", "inflowEnd");

  if(m_stgWallNormalDir % 2 == 0) {
    MFloat temp = inflowStart[m_wallDir];
    inflowStart[m_wallDir] = inflowEnd[m_wallDir];
    inflowEnd[m_wallDir] = temp;
  }

#ifndef NDEBUG
  if(m_stgMyRank == m_commStgRoot) {
    std::cerr << "inflowStartEnd(" << m_bcId << "): " << inflowStart[0] << " " << inflowStart[1] << " "
              << inflowStart[2] << " " << inflowEnd[0] << " " << inflowEnd[1] << " " << inflowEnd[2] << std::endl;
  }
#endif
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::getInflowStartEnd(MFloat* inflowStart, MFloat* inflowEnd) {
  MInt minIndex = getPointIdFromCell(a->m_noGhostLayers, a->m_noGhostLayers, a->m_noGhostLayers);
  MInt maxIndex =
      getPointIdFromCell(a->m_noGhostLayers, a->m_noGhostLayers, a->m_noGhostLayers + m_solver->m_nActiveCells[0]);
  MFloat inflowStartLocal[3] = {LES->a_coordinates(minIndex, 0), LES->a_coordinates(minIndex, 1),
                                LES->a_coordinates(minIndex, 2)};
  MFloat inflowEndlocal[3] = {LES->a_coordinates(maxIndex, 0), LES->a_coordinates(maxIndex, 1),
                              LES->a_coordinates(maxIndex, 2)};
  std::fill(inflowStart, inflowStart + nDim, 99999.9);
  std::fill(inflowEnd, inflowEnd + nDim, F0);
  MPI_Allreduce(inflowStartLocal, inflowStart, nDim, MPI_DOUBLE, MPI_MIN, m_commStg, AT_, "inflowStartLocal",
                "inflowStart");
  MPI_Allreduce(inflowEndlocal, inflowEnd, nDim, MPI_DOUBLE, MPI_MAX, m_commStg, AT_, "inflowEndlocal", "inflowEnd");
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::setVb(MFloat* inflowStart, MFloat* inflowEnd) {
  MFloatScratchSpace bcast_vb(6, AT_, "bcast_vb");


  const MFloat vbDepth = (inflowEnd[m_periodicDir] - inflowStart[m_periodicDir]) * m_stgBLT3;
  const MFloat offset = F1B2 * (inflowEnd[m_periodicDir] - inflowStart[m_periodicDir]) * (m_stgBLT3 - F1);

  if(m_stgMyRank == m_commStgRoot) {
    m_stgRootRank = true;

    // Get the coordinate of the inflow
    // Points in stgDir
    MInt factor = (m_stgFaceNormalDir % 2 == 0) ? 1 : -1;
    bcast_vb[m_stgDir] = inflowStart[m_stgDir] - factor * m_stgBLT1 * m_stgDelta99Inflow * F1B2;
    bcast_vb[m_stgDir + 3] = inflowStart[m_stgDir] + factor * m_stgBLT1 * m_stgDelta99Inflow * F1B2;

    // Points in wallDir
    factor = (m_stgWallNormalDir % 2 != 0) ? 1 : -1;
    bcast_vb[m_wallDir] = inflowStart[m_wallDir];
    bcast_vb[m_wallDir + 3] = inflowStart[m_wallDir] + factor * m_stgBLT2 * m_stgDelta99Inflow;

    if(m_freeStreamTurbulence) {
      m_stgWallEnd = inflowStart[m_wallDir] + factor * m_stgBLT2 * m_stgDelta99Inflow;
      bcast_vb[m_wallDir + 3] = inflowEnd[m_wallDir];
    }

    if(m_isotropicTurbulence) {
      m_stgWallEnd = inflowStart[m_wallDir];
      bcast_vb[m_wallDir + 3] = inflowEnd[m_wallDir];
    }

    // Points in periodicDir
    bcast_vb[m_periodicDir] = inflowStart[m_periodicDir] - offset;
    bcast_vb[m_periodicDir + 3] = inflowStart[m_periodicDir] - offset + vbDepth;
  }

  MPI_Bcast(bcast_vb.begin(), 6, MPI_DOUBLE, m_commStgRoot, m_commStg, AT_, "bcast_vb.begin()");

  if(m_freeStreamTurbulence) {
    MPI_Bcast(&m_stgWallEnd, 1, MPI_DOUBLE, m_commStgRoot, m_commStg, AT_, "m_stgWallEnd");
  }

  m_stgVbStart[0] = bcast_vb[0];
  m_stgVbStart[1] = bcast_vb[1];
  m_stgVbStart[2] = bcast_vb[2];
  m_stgVbEnd[0] = bcast_vb[3];
  m_stgVbEnd[1] = bcast_vb[4];
  m_stgVbEnd[2] = bcast_vb[5];

#ifndef NDEBUG
  if(m_stgMyRank == m_commStgRoot) {
    std::cout << "setVbStart(" << m_bcId << "): " << m_stgVbStart[0] << " " << m_stgVbStart[1] << " " << m_stgVbStart[2]
              << std::endl;
    std::cout << "setVbEnd(" << m_bcId << "): " << m_stgVbEnd[0] << " " << m_stgVbEnd[1] << " " << m_stgVbEnd[2]
              << std::endl;
    if(m_freeStreamTurbulence) {
      std::cout << "setVbEndFS(" << m_bcId << "): " << m_stgWallEnd << endl;
    }
  }
#endif
}


/** \brief write RANSValues from solver to stgLVariables (fully coupled)
 * /author Jannik Borgelt
 */
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeR == MAIA_FINITE_VOLUME, _*>,
          std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::readRANSProfileStg() {
  TRACE();

  if(m_zonal && !m_preliminary) {
    for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
      const MInt IBC = a->getStgId(it);
      const MInt cellId = a->m_stgToCellId[IBC];

      // Debugging
      for(MInt varId = 0; varId < (MInt)m_solver->m_noRANSVariables; varId++) {
        ASSERT(cellId < (MInt)m_solver->m_RANSValues[varId].size(),
               "Trying to access data [" + std::to_string(cellId) + "] in RANSValues with length "
                   + std::to_string(m_solver->m_RANSValues[varId].size())
                   + ", domainId: " + std::to_string(m_solver->domainId()));
      }

      m_stgLVariables[STG::AVG_RHO][IBC] = m_solver->m_RANSValues[PV->RHO][cellId];
      m_stgLVariables[STG::AVG_U][IBC] = m_solver->m_RANSValues[PV->U][cellId];
      m_stgLVariables[STG::AVG_V][IBC] = m_solver->m_RANSValues[PV->V][cellId];
      m_stgLVariables[STG::AVG_W][IBC] = m_solver->m_RANSValues[PV->W][cellId];
      m_stgLVariables[STG::AVG_P][IBC] = m_solver->m_RANSValues[PV->P][cellId];
      // JANNIK: generalize this for other RANS models
      m_stgLVariables[STG::NU_T][IBC] = m_solver->m_RANSValues[PV->N][cellId];

      for(MInt d = 0; d < m_nDim; d++) {
        m_stgLVariables[STG::AVG_UU[d]][IBC] = m_solver->m_RANSValues[PV->VV[d]][cellId];
      }

      // reconstruct fluc values since not saved anymore
      // (for timsm function in calcTotalFluctuationCholesky)
      if(globalTimeStep == m_solver->m_restartTimeStep) {
        m_stgLVariables[STG::FLUC_U][IBC] =
            m_solver->a_pvariable(cellId, PV->U) - m_solver->m_RANSValues[PV->U][cellId];
        m_stgLVariables[STG::FLUC_V][IBC] =
            m_solver->a_pvariable(cellId, PV->V) - m_solver->m_RANSValues[PV->V][cellId];
        m_stgLVariables[STG::FLUC_W][IBC] =
            m_solver->a_pvariable(cellId, PV->W) - m_solver->m_RANSValues[PV->W][cellId];
      }
    }
  } else {
    // read in text-file with
    if(m_preliminaryRans2D) {
      stringstream fn;
      fn.clear();
      MString prop_name = "bc" + std::to_string(m_bcId);
      fn << prop_name << "preliminaryDataRans2D.txt";
      MString fname = fn.str();
      if(m_stgMyRank == m_commStgRoot) cerr << "loading STG preliminary data from " << fname << "...";

      ifstream preliminaryData;
      preliminaryData.open(fname);

      vector<MFloat> data;
      MFloat num;
      string line;

      while(preliminaryData >> num) {
        data.push_back(num);
      }

      // MInt count = 0;
      prop_name = "bc" + std::to_string(m_bcId) + "preliminaryRansDataCount";
      MInt preliminaryDataVarCount =
          Context::getSolverProperty<MInt>(prop_name, m_solverId, AT_, &preliminaryDataVarCount);

      MInt dataCount = data.size() / preliminaryDataVarCount;

      vector<vector<MFloat>> preData(dataCount, vector<MFloat>(preliminaryDataVarCount, 0));

      MInt index = 0;
      for(MInt d = 0; d < dataCount; d++) {
        for(MInt dd = 0; dd < preliminaryDataVarCount; dd++) {
          index = preliminaryDataVarCount * d + dd;
          preData[d][dd] = data[index];
        }
      }

      preliminaryData.close();

      for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
        const MInt IBC = a->getStgId(it);
        const MInt cellId = a->m_stgToCellId[IBC];
        MFloat pos = m_solver->a_coordinate(cellId, m_wallDir);
        // interpolate data
        for(MInt d = 0; d < dataCount - 1; d++) {
          if(pos >= preData[d][0] && pos < preData[d + 1][0]) {
            MInt var = 1;
            m_stgLVariables[STG::AVG_U][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            var = 2;
            m_stgLVariables[STG::AVG_V][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            var = 3;
            m_stgLVariables[STG::AVG_RHO][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            var = 4;
            m_stgLVariables[STG::AVG_P][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            var = 5;
            m_stgLVariables[STG::NU_T][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            var = 6;
            m_stgLVariables[STG::S11][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            var = 7;
            m_stgLVariables[STG::S12][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            var = 8;
            m_stgLVariables[STG::S22][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            var = 9;
            m_stgLVariables[STG::SIJSIJ][IBC] =
                preData[d][var]
                + (preData[d + 1][var] - preData[d][var]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            m_stgLVariables[STG::AVG_W][IBC] = F0;
            m_stgLVariables[STG::S13][IBC] = F0;
            m_stgLVariables[STG::S23][IBC] = F0;
            m_stgLVariables[STG::S33][IBC] = F0;
          }
        }
      }
    }

    if(m_preliminary && !m_preliminaryRans2D) {
      stringstream fn;
      fn.clear();
      MString prop_name = "bc" + std::to_string(m_bcId);
      fn << prop_name << "preliminaryData.txt";

      MString fname = fn.str();

      if(m_stgMyRank == m_commStgRoot) cerr << "loading STG preliminary data from " << fname << "...";

      ifstream preliminaryData;
      preliminaryData.open(fname);

      vector<MFloat> data;
      MFloat num;

      while(preliminaryData >> num) {
        data.push_back(num);
      }

      // MInt count = 0;
      MInt preliminaryDataVarCount = 1 + m_solver->noVariables() + 8; // y + primVars + rms + p'p' + SijSij
      MInt dataCount = data.size() / preliminaryDataVarCount;

      vector<vector<MFloat>> preData(dataCount, vector<MFloat>(preliminaryDataVarCount, 0));

      MInt index = 0;
      for(MInt d = 0; d < dataCount; d++) {
        for(MInt dd = 0; dd < preliminaryDataVarCount; dd++) {
          index = preliminaryDataVarCount * d + dd;
          preData[d][dd] = data[index];
        }
      }

      preliminaryData.close();

      if(m_stgMyRank == m_commStgRoot) cerr << "ok" << endl;

      for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
        const MInt IBC = a->getStgId(it);
        const MInt cellId = a->m_stgToCellId[IBC];

        MFloat pos = m_solver->a_coordinate(cellId, m_wallDir);
        if(m_cylinderTransformation) {
          MFloat y = m_solver->a_coordinate(cellId, m_wallDir);
          MFloat z = m_solver->a_coordinate(cellId, m_periodicDir);
          pos = sqrt(y * y + z * z);
        }

        for(MInt d = 0; d < dataCount; d++) {
          if(abs(pos - preData[d][0]) / pos < 0.001) {
            // no interpolation needed
            for(MInt var = 0; var < nDim + 2; var++) {
              m_stgLVariables[var][IBC] = preData[d][var + 1];
              if(var < nDim) {
                m_stgLVariables[STG::AVG_UU[var]][IBC] = m_stgLVariables[var][IBC];
              }
            }
            MInt v = 6;
            for(MInt fluc = STG::FLUC_UU; fluc < STG::FLUC_WW + 1; fluc++) {
              m_stgLVariables[fluc][IBC] = preData[d][v];
              v++;
            }
            MInt index_ = preliminaryDataVarCount - 1;
            m_stgLVariables[STG::SIJSIJ][IBC] = preData[d][index_];
          } else {
            if(d == dataCount - 1) continue; // this is necessary for interpolation
            // interpolation
            if(pos >= preData[d][0] && pos < preData[d + 1][0]) {
              for(MInt var = 0; var < nDim + 2; var++) {
                m_stgLVariables[var][IBC] = preData[d][var + 1]
                                            + (preData[d + 1][var + 1] - preData[d][var + 1]) * (pos - preData[d][0])
                                                  / (preData[d + 1][0] - preData[d][0]);

                if(var < nDim) {
                  m_stgLVariables[STG::AVG_UU[var]][IBC] = m_stgLVariables[var][IBC];
                }
              }
              MInt v = 6;
              for(MInt fluc = STG::FLUC_UU; fluc < STG::FLUC_WW + 1; fluc++) {
                m_stgLVariables[fluc][IBC] =
                    preData[d][v]
                    + (preData[d + 1][v] - preData[d][v]) * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);

                v++;
              }

              MInt index_ = preliminaryDataVarCount - 1;
              m_stgLVariables[STG::SIJSIJ][IBC] = preData[d][index_]
                                                  + (preData[d + 1][index_] - preData[d][index_])
                                                        * (pos - preData[d][0]) / (preData[d + 1][0] - preData[d][0]);
            }
          }
        }
      }
    }
  }
}


/** \brief write RANSValues from solver to stgLVariables (fully coupled zonal)
 * /author Sutharsan
 */
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeR == MAIA_STRUCTURED, _*>,
          std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::readRANSProfileStg() {
  // Sut: if it is not a zonal run, fill a_pvariable of ghost cells with RANS values
  //     and if it is a m_stgInitialStartup, than also fill m_stgLVariables with RANS values
  if(m_stgInitialStartup) {
    // Reading in RANS Profile for BC
    // Iterate over all stg cells for gradient calculations

    //    ScratchSpace<MFloat> stgCoords(m_nDim, a->sizeStg(), AT_, "stgCoords" );
    std::vector<std::vector<MFloat>> stgCoords(m_nDim);
    for(MInt d = 0; d < m_nDim; ++d)
      stgCoords[d].resize(a->sizeStg());

    // loop over all stg cells
    for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
      const MInt stgId = a->getStgId(it);
      const MInt cellId = a->getCellId(it);
      stgCoords[0][stgId] = LES->a_coordinate(cellId, 0);
      stgCoords[1][stgId] = LES->a_coordinate(cellId, 1);
      stgCoords[2][stgId] = LES->a_coordinate(cellId, 2);
    }

    // hack
    MInt c = 0;
    std::vector<MFloat*> stgCoords_ptr(m_nDim);
    for(auto& s_ptr : stgCoords) {
      stgCoords_ptr[c++] = s_ptr.data();
    }
    //    MInt* temp = new MInt[a->sizeStg()];
    MInt temp2[] = {a->sizeStg(), 1, 1};

    auto* structuredInterpolation = new StructuredInterpolation<3>(m_commStg);
    // TODO: check this later
    //    structuredInterpolation->prepareInterpolation(a->sizeStg(), stgCoords_ptr.data(), temp);
    structuredInterpolation->prepareInterpolationField(&temp2[0], stgCoords_ptr.data());
    std::array<MString, m_nDim + 2> pvariableNames;
    /*    pvariableNames[PV->U] = "u";
        pvariableNames[PV->V] = "v";
        if (nDim==3)
          pvariableNames[PV->W] = "w";
        pvariableNames[PV->P] = "p";
        pvariableNames[PV->RHO] = "rho";*/
    pvariableNames[STG::AVG_U] = "u";
    pvariableNames[STG::AVG_V] = "v";
    IF_CONSTEXPR(nDim == 3) { pvariableNames[STG::AVG_W] = "w"; }
    pvariableNames[STG::AVG_P] = "p";
    pvariableNames[STG::AVG_RHO] = "rho";
    const MInt noVars = nDim + 2;
    for(MInt var = 0; var < noVars; var++) {
      // BE CAREFUL WITH THE ORDER OF THE VARIABLES
      structuredInterpolation->interpolateField(pvariableNames[var], m_stgLVariables[var]);
    }
    // TODO "nu_t" is hardcoded
    structuredInterpolation->interpolateField("nu_t", m_stgLVariables[STG::NU_T]);

    // SANITY CHECK
    MInt cnt = 0;
    for(MInt var = 0; var < noVars; var++) {
      MFloat minVar = std::numeric_limits<MFloat>::max();
      MFloat maxVar = -std::numeric_limits<MFloat>::max();
      for(typename Accessor::nDim_citerator it = a->iterateAll(); it != a->iterateAll_nDim_citerator_end(); ++it) {
        cnt++;
        const MInt stgId = a->getStgId(it);
        MFloat variable = m_stgLVariables[var][stgId];
        if(std::isnan(variable)) TERMM(1, "NAN");
        if(variable > maxVar) maxVar = variable;
        if(variable < minVar) minVar = variable;
      }
      MFloat minVarGlobal;
      MFloat maxVarGlobal;
      MPI_Reduce(&minVar, &minVarGlobal, 1, MPI_DOUBLE, MPI_MIN, m_commStgRoot, m_commStg, AT_, "MPI_IN_PLACE",
                 "minVar");
      MPI_Reduce(&maxVar, &maxVarGlobal, 1, MPI_DOUBLE, MPI_MAX, m_commStgRoot, m_commStg, AT_, "MPI_IN_PLACE",
                 "maxVar");
      if(m_stgMyRank == m_commStgRoot)
        std::cout << " CHECK after readRANSProfileStg of rank=" << m_solver->domainId() << ":var=" << var
                  << std::scientific << ": minVar=" << minVarGlobal << "; maxVar=" << maxVarGlobal << "; " << cnt
                  << std::endl;
    }
    // SANITY CHECK ENDS

  } else {
    // TODO: check if it should be B1 or Ghost cells
    for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
      const MInt cellId = a->getCellId(it);
      const MInt IBC = a->getStgId(it);
      for(int var = 0; var < PV->noVariables; var++) {
        LES->a_pvariable(cellId, var) = m_stgLVariables[var][IBC];
      }
    }
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::calcStrainRate() {
  TRACE();

  if(!m_preliminary) {
    for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
      const MInt cellId = a->getCellId(it);
      const MInt recData = m_solver->a_reconstructionData(cellId);
      MFloat du[nDim][nDim]{{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}};
      const MFloat u[nDim] = {m_stgLVariables[STG::AVG_U][a->getStgId(it)],
                              m_stgLVariables[STG::AVG_V][a->getStgId(it)],
                              m_stgLVariables[STG::AVG_W][a->getStgId(it)]};

      for(MInt nghbr = 0; nghbr < m_solver->a_noReconstructionNeighbors(cellId); nghbr++) {
        const MInt recNghbrId = m_solver->a_reconstructionNeighborId(cellId, nghbr);
        const MInt nghbrStgId = a->getNghbrMapping(a->getStgId(it), nghbr);

        const MFloat recConst_x = m_solver->m_reconstructionConstants[nDim * (recData + nghbr) + 0];
        const MFloat recConst_y = m_solver->m_reconstructionConstants[nDim * (recData + nghbr) + 1];
        const MFloat recConst_z = m_solver->m_reconstructionConstants[nDim * (recData + nghbr) + 2];

        MFloat delta_u = F0;
        for(MInt d = 0; d < nDim; ++d) {
          if(m_cutOff) {
            if(recNghbrId > 0 && recNghbrId < m_solver->c_noCells()) {
              delta_u = m_solver->m_RANSValues[d][recNghbrId] - u[d];
            }
          } else {
            if(nghbrStgId > 0) delta_u = m_stgLVariables[STG::AVG_UU[d]][nghbrStgId] - u[d];
          }
          du[d][0] += recConst_x * delta_u;
          du[d][1] += recConst_y * delta_u;
          du[d][2] += recConst_z * delta_u;
        }
      }
      if(m_solver->m_reConstSVDWeightMode == 3) {
        if(!m_solver->a_hasProperty(cellId, /*SolverCell*/ FvCell::HasCoarseNghbr)) {
          continue;
        }
        std::array<MBool, nDim> dirsJmp = {};
        for(MInt d = 0; d < nDim; ++d) {
          if(m_solver->m_cells.nghbrInterface(cellId, 2 * d) == 3
             || m_solver->m_cells.nghbrInterface(cellId, 2 * d + 1) == 3) {
            dirsJmp[d] = true;
          } else {
            dirsJmp[d] = false;
          }
        }

        for(MUint d = 0; d < nDim; ++d) {
          for(MInt dd = 0; dd < nDim; ++dd) {
            if(dirsJmp[dd]) {
              du[d][dd] = 0.0;
            }
          }
        }

        for(MInt nghbr = 0; nghbr < m_solver->a_noReconstructionNeighbors(cellId); nghbr++) {
          const MInt nghbrStgId = a->getNghbrMapping(a->getStgId(it), nghbr);
          const MInt nghbrId = m_solver->a_reconstructionNeighborId(cellId, nghbr);

          for(MUint d = 0; d < nDim; ++d) {
            MFloat tmp = m_stgLVariables[STG::AVG_UU[d]][nghbrStgId] - u[d];
            for(MInt dd = 0; dd < nDim; ++dd) {
              if(!dirsJmp[dd]) {
                const MFloat dx = LES->a_coordinate(nghbrId, dd) - LES->a_coordinate(cellId, dd);
                tmp -= du[d][dd] * dx;
              }
            }
            for(MInt dd = 0; dd < nDim; ++dd) {
              if(dirsJmp[dd]) {
                du[d][dd] += m_solver->m_reconstructionConstants[nDim * (recData + nghbr) + dd] * tmp;
              }
            }
          }
        }
      }

      const MFloat s11 = 2.0 * du[0][0];
      const MFloat s22 = 2.0 * du[1][1];
      const MFloat s33 = 2.0 * du[2][2];

      const MFloat s12 = du[1][0] + du[0][1];
      const MFloat s13 = du[2][0] + du[0][2];
      const MFloat s23 = du[2][1] + du[1][2];

      const MFloat s21 = s12;
      const MFloat s31 = s13;
      const MFloat s32 = s23;

      // Strain tensor
      const MFloat SijSij =
          F1B4
          * (s11 * s11 + s12 * s12 + s13 * s13 + s21 * s21 + s22 * s22 + s23 * s23 + s31 * s31 + s32 * s32 + s33 * s33);

      m_stgLVariables[STG::S11][a->getStgId(it)] = s11;
      m_stgLVariables[STG::S22][a->getStgId(it)] = s22;
      m_stgLVariables[STG::S33][a->getStgId(it)] = s33;

      m_stgLVariables[STG::S12][a->getStgId(it)] = s12;
      m_stgLVariables[STG::S13][a->getStgId(it)] = s13;
      m_stgLVariables[STG::S23][a->getStgId(it)] = s23;

      m_stgLVariables[STG::SIJSIJ][a->getStgId(it)] = SijSij;
    }
  }
}


#if 0
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::calcStrainRate() {

  const MInt ii = 1;

  for(MInt k = start[2]+1; k < end[2]-1; k++) {
    for(MInt j = start[1]+1; j < end[1]-1; j++) {
      /*********************** C *****************/
      //This is the metrics / strain tensor / max shear part
      /*********************** C *****************/

      MInt IPJK, IMJK, IJPK, IJMK, IJKP, IJKM;

      I = cellIndex(ii,j,k);
      IPJK = cellIndex(ii+1,j,k);
      IMJK = cellIndex(ii-1,j,k);
      IJPK = cellIndex(ii,j+1,k);
      IJMK = cellIndex(ii,j-1,k);
      IJKP = cellIndex(ii,j,k+1);
      IJKM = cellIndex(ii,j,k-1);


      const MFloat dxdi = F1B2 * (m_cells->coordinates[0][IPJK] - m_cells->coordinates[0][IMJK]);
      const MFloat dxdj = F1B2 * (m_cells->coordinates[0][IJPK] - m_cells->coordinates[0][IJMK]);
      const MFloat dxdk = F1B2 * (m_cells->coordinates[0][IJKP] - m_cells->coordinates[0][IJKM]);
      const MFloat dydi = F1B2 * (m_cells->coordinates[1][IPJK] - m_cells->coordinates[1][IMJK]);
      const MFloat dydj = F1B2 * (m_cells->coordinates[1][IJPK] - m_cells->coordinates[1][IJMK]);
      const MFloat dydk = F1B2 * (m_cells->coordinates[1][IJKP] - m_cells->coordinates[1][IJKM]);
      const MFloat dzdi = F1B2 * (m_cells->coordinates[2][IPJK] - m_cells->coordinates[2][IMJK]);
      const MFloat dzdj = F1B2 * (m_cells->coordinates[2][IJPK] - m_cells->coordinates[2][IJMK]);
      const MFloat dzdk = F1B2 * (m_cells->coordinates[2][IJKP] - m_cells->coordinates[2][IJKM]);

      const MFloat dxl = sqrt(dxdi*dxdi + dydi*dydi + dzdi*dzdi);
      const MFloat dyl = sqrt(dxdj*dxdj + dydj*dydj + dzdj*dzdj);
      const MFloat dzl = sqrt(dxdk*dxdk + dydk*dydk + dzdk*dzdk);

      const MFloat dxidx = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][0];
      const MFloat dxidy = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][1];
      const MFloat dxidz = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][2];

      const MFloat detadx = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][3 + 0];
      const MFloat detady = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][3 + 1];
      const MFloat detadz = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][3 + 2];

      const MFloat dzetadx = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][6 + 0];
      const MFloat dzetady = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][6 + 1];
      const MFloat dzetadz = (1./std::max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][6 + 2];

      IBC = cellIndexBC(ii,j,k);
      IPJK = cellIndexBC(ii+1,j,k);
      IMJK = cellIndexBC(ii-1,j,k);
      IJPK = cellIndexBC(ii,j+1,k);
      IJMK = cellIndexBC(ii,j-1,k);
      IJKP = cellIndexBC(ii,j,k+1);
      IJKM = cellIndexBC(ii,j,k-1);

      const MFloat frho = 1.0 / m_stgVariables[STG::AVG_RHO][IBC];

      //dud?
      const MFloat dudxi = F1B2*(m_stgVariables[STG::AVG_U][IPJK] - m_stgVariables[STG::AVG_U][IMJK]);
      const MFloat dudeta = F1B2*(m_stgVariables[STG::AVG_U][IJPK] - m_stgVariables[STG::AVG_U][IJMK]);
      const MFloat dudzeta = F1B2*(m_stgVariables[STG::AVG_U][IJKP] - m_stgVariables[STG::AVG_U][IJKM]);

      //dvd?
      const MFloat dvdxi = F1B2*(m_stgVariables[STG::AVG_V][IPJK] - m_stgVariables[STG::AVG_V][IMJK]);
      const MFloat dvdeta = F1B2*(m_stgVariables[STG::AVG_V][IJPK] - m_stgVariables[STG::AVG_V][IJMK]);
      const MFloat dvdzeta = F1B2*(m_stgVariables[STG::AVG_V][IJKP] - m_stgVariables[STG::AVG_V][IJKM]);

      //dwd?
      const MFloat dwdxi = F1B2*(m_stgVariables[STG::AVG_W][IPJK] - m_stgVariables[STG::AVG_W][IMJK]);
      const MFloat dwdeta = F1B2*(m_stgVariables[STG::AVG_W][IJPK] - m_stgVariables[STG::AVG_W][IJMK]);
      const MFloat dwdzeta = F1B2*(m_stgVariables[STG::AVG_W][IJKP] - m_stgVariables[STG::AVG_W][IJKM]);

      const MFloat dudx = dudxi*dxidx + dudeta*detadx + dudzeta*dzetadx;
      const MFloat dudy = dudxi*dxidy + dudeta*detady + dudzeta*dzetady;
      const MFloat dudz = dudxi*dxidz + dudeta*detadz + dudzeta*dzetadz;

      const MFloat dvdx = dvdxi*dxidx + dvdeta*detadx + dvdzeta*dzetadx;
      const MFloat dvdy = dvdxi*dxidy + dvdeta*detady + dvdzeta*dzetady;
      const MFloat dvdz = dvdxi*dxidz + dvdeta*detadz + dvdzeta*dzetadz;

      const MFloat dwdx = dwdxi*dxidx + dwdeta*detadx + dwdzeta*dzetadx;
      const MFloat dwdy = dwdxi*dxidy + dwdeta*detady + dwdzeta*dzetady;
      const MFloat dwdz = dwdxi*dxidz + dwdeta*detadz + dwdzeta*dzetadz;


      m_stgVariables[STG::S11][cellIndexBC(ii,j,k)] = 2.0*dudx;
      m_stgVariables[STG::S22][cellIndexBC(ii,j,k)] = 2.0*dvdy;
      m_stgVariables[STG::S33][cellIndexBC(ii,j,k)] = 2.0*dwdz;

      m_stgVariables[STG::S12][cellIndexBC(ii,j,k)] = dvdx+dudy;
      m_stgVariables[STG::S13][cellIndexBC(ii,j,k)] = dwdx+dudz;
      m_stgVariables[STG::S23][cellIndexBC(ii,j,k)] = dwdy+dvdz;
    }
  }
}
#endif


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::calcReynoldsStressConvVelLengthScale() {
  // Arrays for MPI operations
  MFloatScratchSpace maxValsLocal(3, AT_, "maxValsLocal");
  MInt size;
  MPI_Comm_size(m_commStg, &size);
  MFloatScratchSpace maxValsGlobal(3 * size, AT_, "maxValsGlobal");
  maxValsLocal.fill(F0);
  maxValsGlobal.fill(F0);

  const MFloat fre = 1.0 / m_solver->sysEqn().m_Re0;

  const MFloat delta_in = m_stgDelta99Inflow;

  // Initialize max values
  MFloat utaumax = F0, umax = F0, vmax = F0, wmax = F0, minLengthLocal = F0;

  for(typename Accessor::nDim_citerator it = a->iterateSlopes(); it != a->iterateSlopes_nDim_citerator_end(); ++it) {
    const MInt id = a->getStgId(it);
    const MInt cellId = a->getCellId(it);

    const MFloat frho = 1.0 / m_stgLVariables[STG::AVG_RHO][id];

    MFloat SijSij = m_stgLVariables[STG::SIJSIJ][id];

    if(!m_preliminary || m_preliminaryRans2D) {
      // Strain tensor
      // const MFloat s11 = m_stgLVariables[STG::S11][id];
      // const MFloat s22 = m_stgLVariables[STG::S22][id];
      // const MFloat s33 = m_stgLVariables[STG::S33][id];
      const MFloat s12 = m_stgLVariables[STG::S12][id];
      const MFloat s23 = m_stgLVariables[STG::S23][id];
      const MFloat s13 = m_stgLVariables[STG::S13][id];
      const MFloat s21 = s12;
      const MFloat s32 = s23;
      const MFloat s31 = s13;

      if(m_freeStreamTurbulence) {
        MFloat SijSij_ = SijSij;
        MFloat tmp = LES->a_coordinate(cellId, m_wallDir);
        if(m_cylinderTransformation) {
          MFloat y = LES->a_coordinate(cellId, m_wallDir);
          MFloat z = LES->a_coordinate(cellId, m_periodicDir);
          tmp = sqrt(y * y + z * z);
        }
        if(tmp > m_stgWallEnd) {
          SijSij_ = m_SijSijFS;
        }
        SijSij = SijSij_;
        m_stgLVariables[STG::SIJSIJ][id] = SijSij;
      }

      //>marian: in TFS code this isn't SQRT
      // m_stgLVariables[STG::SIJSIJ][id] = SijSij;

      // Assume a one-, or two equation turbulence model
      // Assume a simplified directivity:
      // Read from RANS profile

      MFloat nu_t = m_stgLVariables[STG::NU_T][id];

      const MFloat sr1 = (s12 + s21) * (s12 + s21);
      const MFloat sr2 = (s23 + s32) * (s23 + s32);
      const MFloat sr3 = (s13 + s31) * (s13 + s31);
      const MFloat srt = std::max(sqrt(sr1 + sr2 + sr3), epsl);

      const MFloat rr1 = sqrt(sr1) / srt;
      const MFloat rr2 = sqrt(sr2) / srt;
      const MFloat rr3 = sqrt(sr3) / srt;

      const MFloat uv = -sqrt(2.0 * SijSij) * rr1 * nu_t * fre;
      const MFloat vw = -sqrt(2.0 * SijSij) * rr2 * nu_t * fre;
      const MFloat uw = -sqrt(2.0 * SijSij) * rr3 * nu_t * fre;
      const MFloat uu = a1 * std::abs(uv) * m_stgRSTFactors[0];
      const MFloat vv = a1 * std::abs(uv) * m_stgRSTFactors[1];
      const MFloat ww = a1 * std::abs(uv) * m_stgRSTFactors[2];

      m_stgLVariables[STG::FLUC_UU][id] = uu;
      m_stgLVariables[STG::FLUC_UV][id] = uv;
      m_stgLVariables[STG::FLUC_VV][id] = vv;
      m_stgLVariables[STG::FLUC_WW][id] = ww;
      m_stgLVariables[STG::FLUC_VW][id] = vw;
      m_stgLVariables[STG::FLUC_UW][id] = uw;
    } //! m_preliminary
    else {
      MFloat uv = m_stgLVariables[STG::FLUC_UV][id];
      m_stgLVariables[STG::FLUC_UU][id] = a1 * std::abs(uv) * m_stgRSTFactors[0];
      m_stgLVariables[STG::FLUC_VV][id] = a1 * std::abs(uv) * m_stgRSTFactors[1];
      m_stgLVariables[STG::FLUC_WW][id] = a1 * std::abs(uv) * m_stgRSTFactors[2];
      m_stgLVariables[STG::FLUC_VW][id] = 0;
      m_stgLVariables[STG::FLUC_UW][id] = 0;
    }

    if(m_freeStreamTurbulence) {
      m_stgLVariables[STG::FLUC_UU][id] = m_uuFS;
      m_stgLVariables[STG::FLUC_UV][id] = 0;
      m_stgLVariables[STG::FLUC_VV][id] = m_vvFS;
      m_stgLVariables[STG::FLUC_WW][id] = m_wwFS;
      m_stgLVariables[STG::FLUC_VW][id] = 0;
      m_stgLVariables[STG::FLUC_UW][id] = 0;
    }

    // TODO: this is different from the implementation in the Structured solver
    MFloat rhoRANSI1 = m_stgLVariables[PV->RHO][id];
    MFloat pressure1 = m_stgLVariables[PV->P][id];

    const MFloat temp = m_solver->m_gamma * pressure1 / rhoRANSI1;
    const MFloat xmu = SUTHERLANDLAW(temp);

    // Get utau using laminar viscosity
    MFloat utau2_ = sqrt(fre * sqrt(2.0 * SijSij) * xmu * frho);

    const MFloat utau2 = utau2_;

    // TODO: check the following
    const MFloat dV = getCellSize(it);

    // Save values if they are the new maximum
    if(utau2 >= utaumax) {
      minLengthLocal = pow(dV, 0.33);
      utaumax = std::max(utau2, utaumax);
    }

    // In which direction aims the maximum averaged velocity?
    // TODO: think if RANS or LES velocity?
    const MFloat u = m_stgLVariables[STG::AVG_U][id];
    const MFloat v = m_stgLVariables[STG::AVG_V][id];
    const MFloat w = m_stgLVariables[STG::AVG_W][id];

    // TODO: check the following
    umax = (fabs(u) > fabs(umax)) ? u : umax; // std::max(u, umax);
    vmax = (fabs(v) > fabs(vmax)) ? v : vmax; // std::max(v, vmax);
    wmax = (fabs(w) > fabs(wmax)) ? w : wmax; // std::max(w, wmax);

    // We need a global length scale to compare
    m_stgLVariables[STG::LENGTH_SCALE][id] = pow(sqrt(std::max(2.0 * SijSij, epss)), -m_stgExple);

    // TODO: next section only for Structured:
    // Zero gradient extrapolation of boundary values
    // extrapolateToBoundary(it);

    // Save max direction vector and max tau
    maxValsLocal[0] = umax;
    maxValsLocal[1] = vmax;
    maxValsLocal[2] = wmax;

    if(std::isnan(utaumax) || std::isinf(utaumax)) utaumax = F0;
  }

  // Communication: Exchange min and max values
  MFloat minLengthGlobal = 0.0;
  MFloat utaux;
  MPI_Allgather(maxValsLocal.begin(), 3, MPI_DOUBLE, maxValsGlobal.begin(), 3, MPI_DOUBLE, m_commStg, AT_,
                "maxValsLocal.begin()", "maxValsGlobal.begin()");
  MPI_Allreduce(&utaumax, &utaux, 1, MPI_DOUBLE, MPI_MAX, m_commStg, AT_, "utaumax", "utaux");

  MPI_Allreduce(&minLengthLocal, &minLengthGlobal, 1, MPI_DOUBLE, MPI_MIN, m_commStg, AT_, "minLengthLocal",
                "minLengthGlobal");

  // Maximum convection velocities at inflow
  std::fill_n(&m_stgMaxVel[0], 3, 0.0);
  for(MInt d = 0; d < size; ++d) {
    m_stgMaxVel[0] = (fabs(m_stgMaxVel[0]) > fabs(maxValsGlobal[3 * d])) ? m_stgMaxVel[0] : maxValsGlobal[3 * d];
    m_stgMaxVel[1] =
        (fabs(m_stgMaxVel[1]) > fabs(maxValsGlobal[3 * d + 1])) ? m_stgMaxVel[1] : maxValsGlobal[3 * d + 1];
    m_stgMaxVel[2] =
        (fabs(m_stgMaxVel[2]) > fabs(maxValsGlobal[3 * d + 2])) ? m_stgMaxVel[2] : maxValsGlobal[3 * d + 2];
  }

#ifndef NDEBUG
  if(m_stgMyRank == m_commStgRoot) {
    cerr << "m_stgMaxVel: " << utaux << " " << m_stgMaxVel[0] << " " << m_stgMaxVel[1] << " " << m_stgMaxVel[2] << endl;
  }
#endif

  for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
    const MInt IBC = it.getStgId();
    // const MInt cellId = a->getCellId(it);

    MFloat maxLength = delta_in;

    // Length scale in main flow direction
    const MFloat xlength = std::max(
        std::min(
            m_stgLengthFactors[0]
                * std::max(m_stgLVariables[STG::LENGTH_SCALE][IBC] * delta_in * pow(utaux / delta_in, m_stgExple), eps),
            maxLength * 1.0),
        minLengthGlobal);

    // Length scale in the direction of main shear
    const MFloat ylength = std::max(
        std::min(
            m_stgLengthFactors[1]
                * std::max(m_stgLVariables[STG::LENGTH_SCALE][IBC] * delta_in * pow(utaux / delta_in, m_stgExple), eps),
            maxLength * 0.66),
        minLengthGlobal);

    // Length scale in the direction perpendicular of x and y
    const MFloat zlength = std::max(
        std::min(
            m_stgLengthFactors[2]
                * std::max(m_stgLVariables[STG::LENGTH_SCALE][IBC] * delta_in * pow(utaux / delta_in, m_stgExple), eps),
            maxLength * 1.0),
        minLengthGlobal);

    m_stgLVariables[STG::LENGTH_X][IBC] = xlength;
    m_stgLVariables[STG::LENGTH_Y][IBC] = ylength;
    m_stgLVariables[STG::LENGTH_Z][IBC] = zlength;
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*>>
MFloat MSTG<nDim, SolverTypeR, SolverTypeL>::getCellSize(typename Accessor::nDim_citerator it) {
  return m_solver->a_cellVolume(a->getCellId(it));
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*>>
MFloat MSTG<nDim, SolverTypeR, SolverTypeL>::getCellSize(typename Accessor::nDim_citerator it) {
  const MInt IPJK = a->getNghbr(it, 1);
  const MInt IMJK = a->getNghbr(it, 0);
  const MInt IJPK = a->getNghbr(it, 3);
  const MInt IJMK = a->getNghbr(it, 2);
  const MInt IJKP = a->getNghbr(it, 5);
  const MInt IJKM = a->getNghbr(it, 4);

  const MFloat dxdi = F1B2 * (LES->a_coordinate(IPJK, 0) - LES->a_coordinate(IMJK, 0));
  const MFloat dxdj = F1B2 * (LES->a_coordinate(IJPK, 0) - LES->a_coordinate(IJMK, 0));
  const MFloat dxdk = F1B2 * (LES->a_coordinate(IJKP, 0) - LES->a_coordinate(IJKM, 0));
  const MFloat dydi = F1B2 * (LES->a_coordinate(IPJK, 1) - LES->a_coordinate(IMJK, 1));
  const MFloat dydj = F1B2 * (LES->a_coordinate(IJPK, 1) - LES->a_coordinate(IJMK, 1));
  const MFloat dydk = F1B2 * (LES->a_coordinate(IJKP, 1) - LES->a_coordinate(IJKM, 1));
  const MFloat dzdi = F1B2 * (LES->a_coordinate(IPJK, 2) - LES->a_coordinate(IMJK, 2));
  const MFloat dzdj = F1B2 * (LES->a_coordinate(IJPK, 2) - LES->a_coordinate(IJMK, 2));
  const MFloat dzdk = F1B2 * (LES->a_coordinate(IJKP, 2) - LES->a_coordinate(IJKM, 2));

  const MFloat dxl = sqrt(dxdi * dxdi + dydi * dydi + dzdi * dzdi);
  const MFloat dyl = sqrt(dxdj * dxdj + dydj * dydj + dzdj * dzdj);
  const MFloat dzl = sqrt(dxdk * dxdk + dydk * dydk + dzdk * dzdk);

  return dxl * dyl * dzl;
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::generateNewEddies(MFloat* xk, MFloat* epsi) {
  xk[m_stgDir] = m_stgVbStart[m_stgDir] + generate_rand() * (m_stgVbEnd[m_stgDir] - m_stgVbStart[m_stgDir]);
  xk[m_periodicDir] =
      m_stgVbStart[m_periodicDir] + generate_rand() * (m_stgVbEnd[m_periodicDir] - m_stgVbStart[m_periodicDir]);
  xk[m_wallDir] =
      m_stgVbStart[m_wallDir] + generate_rand_weighted() * (m_stgVbEnd[m_wallDir] - m_stgVbStart[m_wallDir]);

  if(m_isotropicTurbulence) {
    xk[m_wallDir] = m_stgVbStart[m_wallDir] + generate_rand() * (m_stgVbEnd[m_wallDir] - m_stgVbStart[m_wallDir]);
  }

  if(m_cylinderTransformation) {
    xk[m_wallDir] = m_stgVbStart[m_wallDir] + generate_rand() * (m_stgVbEnd[m_wallDir] - m_stgVbStart[m_wallDir]);
  }

  if(m_freeStreamTurbulence) {
    if(generate_rand() < m_BLEddieFraction) {
      // generate new eddie in BL
      xk[m_wallDir] = m_stgVbStart[m_wallDir] + generate_rand_weighted() * (m_stgWallEnd - m_stgVbStart[m_wallDir]);
    } else {
      // generate new eddie in FS
      xk[m_wallDir] = m_stgWallEnd + generate_rand() * (m_stgVbEnd[m_wallDir] - m_stgWallEnd);
    }
  }

  epsi[m_stgDir] = 2.0 * generate_rand() - 1.0;
  epsi[m_stgDir] = epsi[m_stgDir] / std::max(std::abs(epsi[m_stgDir]), eps);
  epsi[m_wallDir] = 2.0 * generate_rand() - 1.0;
  epsi[m_wallDir] = epsi[m_wallDir] / std::max(std::abs(epsi[m_wallDir]), eps);
  epsi[m_periodicDir] = 2.0 * generate_rand() - 1.0;
  epsi[m_periodicDir] = epsi[m_periodicDir] / std::max(std::abs(epsi[m_periodicDir]), eps);
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::advanceEddies() {
  // Array for MPI operation
  MFloatScratchSpace eddyBcastBuffer(m_stgMaxNoEddies * m_stgNoEddieProperties, AT_, "eddyBcastBuffer");
  eddyBcastBuffer.fill(F0);

  MFloatScratchSpace eddyBcastStrength(m_stgMaxNoEddies, AT_, "eddyBcastStrength");
  eddyBcastStrength.fill(F0);


  MFloat xk[3];
  MFloat epsi[3];
  if(m_stgMyRank == m_commStgRoot) {
    for(MInt n = 0; n < m_stgMaxNoEddies; n++) {
      xk[0] = m_stgEddies[n][0];
      xk[1] = m_stgEddies[n][1];
      xk[2] = m_stgEddies[n][2];

      eddyBcastStrength[n] = m_stgEddyStrength[n];

      // Check if the eddie has left the Virtual Box
      if(xk[0] > m_stgVbEnd[0] || xk[0] < m_stgVbStart[0] || xk[1] > m_stgVbEnd[1] || xk[1] < m_stgVbStart[1]
         || xk[2] > m_stgVbEnd[2] || xk[2] < m_stgVbStart[2]) {
        generateNewEddies(xk, epsi);

        eddyBcastStrength[n] = F1;
        if(m_newStgMethod) {
          eddyBcastStrength[n] = generate_rand_normalDist();
        }

      } else {
        xk[0] = m_stgMaxVel[0] * m_solver->m_timeStep + m_stgEddies[n][0];
        xk[1] = m_stgMaxVel[1] * m_solver->m_timeStep + m_stgEddies[n][1];
        xk[2] = m_stgMaxVel[2] * m_solver->m_timeStep + m_stgEddies[n][2];

        epsi[0] = m_stgEddies[n][3];
        epsi[1] = m_stgEddies[n][4];
        epsi[2] = m_stgEddies[n][5];
      }

      eddyBcastBuffer[n + m_stgMaxNoEddies * 0] = xk[0];
      eddyBcastBuffer[n + m_stgMaxNoEddies * 1] = xk[1];
      eddyBcastBuffer[n + m_stgMaxNoEddies * 2] = xk[2];
      eddyBcastBuffer[n + m_stgMaxNoEddies * 3] = epsi[0];
      eddyBcastBuffer[n + m_stgMaxNoEddies * 4] = epsi[1];
      eddyBcastBuffer[n + m_stgMaxNoEddies * 5] = epsi[2];
    }
  }

  // Broadcast the new/updated eddies to all relevant processes
  MPI_Bcast(eddyBcastBuffer.begin(), m_stgMaxNoEddies * m_stgNoEddieProperties, MPI_DOUBLE, m_commStgRoot, m_commStg,
            AT_, "eddyBcastBuffer.begin()");

  MPI_Bcast(eddyBcastStrength.begin(), m_stgMaxNoEddies, MPI_DOUBLE, m_commStgRoot, m_commStg, AT_,
            "eddyBcastStrength.begin()");

  // Copy data into m_FQeddie vector
  for(MInt n = 0; n < m_stgMaxNoEddies; n++) {
    for(MInt p = 0; p < m_stgNoEddieProperties; p++) {
      m_stgEddies[n][p] = eddyBcastBuffer[n + m_stgMaxNoEddies * p];
    }
    m_stgEddyStrength[n] = eddyBcastStrength[n];
  }
}

template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::printSTGSummary() {
  const MFloat BLT3 = std::abs(m_stgVbEnd[m_periodicDir] - m_stgVbStart[m_periodicDir]);
  const MFloat Vb = m_stgBLT2 * BLT3 * m_stgBLT1 * POW2(m_stgDelta99Inflow);

  // Summary of synth turb parameters
  if(m_stgMyRank == m_commStgRoot && globalTimeStep == m_solver->m_restartTimeStep) {
    std::cout << "**************************" << std::endl
              << "Synthetic turbulence (bcId=" << m_bcId << "):" << std::endl
              << "zones: 1" << std::endl
              << "Dir: " << m_stgFaceNormalDir << " / " << m_stgWallNormalDir << " / " << m_stgDir << " / " << m_wallDir
              << std::endl
              << "nr. eddies: " << m_stgMaxNoEddies << std::endl
              << "conv. vel: " << sqrt(POW2(m_stgMaxVel[0]) + POW2(m_stgMaxVel[1]) + POW2(m_stgMaxVel[2])) << std::endl
              << "umax = " << m_stgMaxVel[0] << std::endl
              << "vmax = " << m_stgMaxVel[1] << std::endl
              << "wmax = " << m_stgMaxVel[2] << std::endl
              << "virtual box volume: " << Vb << std::endl
              << "Vb/stgMaxNoEddies = " << Vb / m_stgMaxNoEddies << std::endl
              << "**************************" << std::endl;
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::calcTotalFluctuationCholesky() {
  TRACE();

  const MFloat BLT3 = std::abs(m_stgVbEnd[m_periodicDir] - m_stgVbStart[m_periodicDir]);
  MFloat Vb = m_stgBLT1 * m_stgBLT2 * POW2(m_stgDelta99Inflow) * BLT3;
  if(m_freeStreamTurbulence || m_isotropicTurbulence) {
    const MFloat BLT2 = std::abs(m_stgVbEnd[m_wallDir] - m_stgVbStart[m_wallDir]);
    Vb = m_stgBLT1 * m_stgDelta99Inflow * BLT2 * BLT3;
  }
  // const MFloat Vb = m_stgBLT1 * m_stgBLT2 * POW2(m_stgDelta99Inflow) * BLT3;
  const MFloat vbFactor = sqrt(Vb / m_stgMaxNoEddies);

  const MFloat umax = m_stgMaxVel[0];
  const MFloat vmax = m_stgMaxVel[1];
  const MFloat wmax = m_stgMaxVel[2];

  for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
    const MInt cellIdBC = a->getStgId(it);
    const MInt cellIdBCFirst = cellIdBC;
    const MInt cellId = a->getCellId(it);

    MFloat help1 = F0, help2 = F0, help3 = F0, help4 = F0, help5 = F0, help6 = F0;

    const MFloat uu = m_stgLVariables[STG::FLUC_UU][cellIdBCFirst];
    const MFloat vv = m_stgLVariables[STG::FLUC_VV][cellIdBCFirst];
    const MFloat ww = m_stgLVariables[STG::FLUC_WW][cellIdBCFirst];
    const MFloat uv = m_stgLVariables[STG::FLUC_UV][cellIdBCFirst];
    const MFloat vw = m_stgLVariables[STG::FLUC_VW][cellIdBCFirst];
    const MFloat uw = m_stgLVariables[STG::FLUC_UW][cellIdBCFirst];

    // Cholesky decomposition of the Reynolds stress tensor
    const MFloat a11 = sqrt(std::max(uu, epsl));
    const MFloat a21 = uv / a11;
    const MFloat a31 = uw / a11;
    const MFloat a22 = sqrt(std::max((vv - a21 * a21), epsl));
    const MFloat a32 = (vw - a21 * a31) / a22;
    const MFloat a33 = sqrt(std::max((ww - a31 * a31 - a32 * a32), epsl));

    for(MInt n = 0; n < m_stgMaxNoEddies; n++) {
      // Tent function to determine the symmetric function to model the
      // decay of the fluctuations
      const MFloat xk1t = m_stgEddies[n][0];
      const MFloat xk2t = m_stgEddies[n][1];
      const MFloat xk3t = m_stgEddies[n][2];

      MFloat distX = LES->a_coordinate(cellId, 0) - xk1t;
      MFloat distY = LES->a_coordinate(cellId, 1) - xk2t;
      MFloat distZ = LES->a_coordinate(cellId, 2) - xk3t;

      if(m_cylinderTransformation) {
        MFloat y = LES->a_coordinate(cellId, m_wallDir);
        MFloat z = LES->a_coordinate(cellId, m_periodicDir);
        MFloat r = sqrt(y * y + z * z);
        MFloat alpha = get_angle(y, z);
        distY = r - xk2t;
        distZ = std::fmod(r * (alpha - xk3t), 2 * r * PI); // Distance of 2*PI means zero distance
      }

      // if(m_newStgMethod){
      const MFloat xLb1 = m_stgLVariables[STG::LENGTH_X][cellIdBCFirst] * m_stgEddyStrength[n];
      const MFloat xLb2 = m_stgLVariables[STG::LENGTH_Y][cellIdBCFirst] * m_stgEddyStrength[n];
      const MFloat xLb3 = m_stgLVariables[STG::LENGTH_Z][cellIdBCFirst] * m_stgEddyStrength[n];
      //}

      const MFloat fxLb1 = 1.0 / xLb1;
      const MFloat fxLb2 = 1.0 / xLb2;
      const MFloat fxLb3 = 1.0 / xLb3;

      const MFloat fsqrtxLb1 = 1.0 / sqrt(xLb1);
      const MFloat fsqrtxLb2 = 1.0 / sqrt(xLb2);
      const MFloat fsqrtxLb3 = 1.0 / sqrt(xLb3);

      const MFloat fsqrtPixLb1 = 1.0 / sqrt(3.141 * xLb1);
      const MFloat fsqrtPixLb2 = 1.0 / sqrt(3.141 * xLb2);
      const MFloat fsqrtPixLb3 = 1.0 / sqrt(3.141 * xLb3);

      const MFloat aDistX = fabs(distX);
      const MFloat aDistY = fabs(distY);
      const MFloat aDistZ = fabs(distZ);

      // only compute contribution if eddie is in vicinity, i.e.,
      // if it is within 4 eddy lengthscales
      if(aDistX < 4.0 * xLb1 && aDistY < 4.0 * xLb2 && aDistZ < 4.0 * xLb3) {
        const MFloat zacfq1 = distX / aDistX;
        const MFloat rol1H = zacfq1 * std::min(aDistX * fxLb1, 1.0);

        const MFloat zacfq2 = distY / aDistY;
        const MFloat rol2H = zacfq2 * std::min(aDistY * fxLb2, 1.0);

        const MFloat zacfq3 = distZ / aDistZ;
        const MFloat rol3H = zacfq3 * std::min(aDistZ * fxLb3, 1.0);

        const MFloat fl1 = 2.0 * fsqrtPixLb1 * exp(-((distX)*2 * fxLb1) * ((distX)*2 * fxLb1));
        const MFloat fl2 = 2.0 * fsqrtPixLb2 * exp(-((distY)*2 * fxLb2) * ((distY)*2 * fxLb2));
        const MFloat fl3 = 2.0 * fsqrtPixLb3 * exp(-((distZ)*2 * fxLb3) * ((distZ)*2 * fxLb3));

        // Normalization factor cannot be chosen as Pamies did... or we...
        MFloat fH1 = (F1 - cos(2.0 * 3.141 * rol1H)) / (2.0 * 3.141 * rol1H * 0.44);
        MFloat fH2 = (F1 - cos(2.0 * 3.141 * rol2H)) / (2.0 * 3.141 * rol2H * 0.44);
        MFloat fH3 = (F1 - cos(2.0 * 3.141 * rol3H)) / (2.0 * 3.141 * rol3H * 0.44);

        fH1 = aniso * fH1 * fsqrtxLb1 + fabs(aniso - 1.0) * fl1;
        fH2 = aniso * fH2 * fsqrtxLb2 + fabs(aniso - 1.0) * fl2;
        fH3 = aniso * fH3 * fsqrtxLb3 + fabs(aniso - 1.0) * fl3;

        const MFloat epsik1 = m_stgEddies[n][3];
        const MFloat epsik2 = m_stgEddies[n][4];
        const MFloat epsik3 = m_stgEddies[n][5];

        help4 += vbFactor * epsik1 * fl1 * fl2 * fH3;
        help5 += vbFactor * epsik2 * fl1 * fl2 * fH3;
        help6 += vbFactor * epsik3 * fl1 * fH2 * fl3;

        // if(m_newStgMethod){
        //   MFloat coverage = POW2(fl1 * fl2 * fH3);
        //   MFloat T = 2 * m_stgLVariables[STG::LENGTH_X][cellIdBC] / m_stgMaxVel[m_stgDir];
        //   MFloat alpha = m_solver->timeStep() / T;
        //   m_stgEddieCoverage[0][cellIdBC] = alpha * coverage + (1 - alpha) * m_stgEddieCoverage[0][cellIdBC];

        //   //transfer
        //   m_solver->m_stgEddieCoverage[0][cellId] = m_stgEddieCoverage[0][cellIdBC];
        // }
      }
    }

    help1 = help4 * a11;                             // Fluctuation u'
    help2 = help4 * a21 + help5 * a22;               // Fluctuation v'
    help3 = help4 * a31 + help5 * a32 + help6 * a33; // Fluctuation w'

    const MFloat velmax = sqrt(umax * umax + vmax * vmax + wmax * wmax);

    MFloat ufluc = std::min(std::max(help1, -0.3 * velmax), 0.3 * velmax);
    MFloat vfluc = std::min(std::max(help2, -0.3 * velmax), 0.3 * velmax);
    MFloat wfluc = std::min(std::max(help3, -0.3 * velmax), 0.3 * velmax);

    if(m_cylinderTransformation) {
      MFloat y = m_solver->a_coordinate(cellId, m_wallDir);
      MFloat z = m_solver->a_coordinate(cellId, m_periodicDir);
      MFloat alpha = get_angle(y, z);
      MFloat vtemp = vfluc;
      MFloat wtemp = wfluc;
      vfluc = vtemp * cos(alpha) + wtemp * sin(alpha);
      wfluc = vtemp * sin(alpha) - wtemp * cos(alpha);
    }

    m_stgLVariables[STG::FLUC_U][cellIdBC] += timsm * (ufluc - m_stgLVariables[STG::FLUC_U][cellIdBC]);
    m_stgLVariables[STG::FLUC_V][cellIdBC] += timsm * (vfluc - m_stgLVariables[STG::FLUC_V][cellIdBC]);
    m_stgLVariables[STG::FLUC_W][cellIdBC] += timsm * (wfluc - m_stgLVariables[STG::FLUC_W][cellIdBC]);
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
void MSTG<nDim, SolverTypeR, SolverTypeL>::calcEddieCoverage() {
  TRACE();

  for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
    const MInt IBC = a->getStgId(it);

    m_stgEddieCoverage[1][IBC] = F0;
  }

  std::vector<MFloat> eddieCoverage(m_stgGlobalNoWallNormalLocations, F0);

  for(MInt i = 0; i < m_stgGlobalNoWallNormalLocations; i++) {
    for(MInt p = 0; p < (MInt)m_stgPeriodicCellId[i].size(); p++) {
      MInt id = m_stgPeriodicCellId[i][p];
      eddieCoverage[i] += m_stgEddieCoverage[0][id] / m_stgGlobalNoPeriodicLocations[i];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &eddieCoverage[0], m_stgGlobalNoWallNormalLocations, MPI_DOUBLE, MPI_SUM, m_commStg, AT_,
                "MPI_IN_PLACE", "m_LESPeriodicAverage");

  for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
    const MInt IBC = a->getStgId(it);
    const MInt cellId = a->getCellId(it);

    MFloat index = m_stgPeriodicIndex[IBC];
    m_stgEddieCoverage[1][IBC] = eddieCoverage[index];

    m_solver->m_stgEddieCoverage[1][cellId] = m_stgEddieCoverage[1][IBC];
    m_solver->m_stgEddieCoverage[2][cellId] = m_stgEddieCoverage[1][IBC] - m_stgEddieCoverage[0][IBC];
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::extrapolateToBoundary(typename Accessor::nDim_citerator it) {
  const MInt noSTGVariables = STG::noSTGVars - PV->noVariables;

  const MInt k = it.getijk(2);
  if(k == a->start(2) + 1) {
    // 1st layer
    for(MInt var = PV->noVariables; var < noSTGVariables; ++var) {
      m_stgLVariables[var][it.getNghbrStg(4)] = m_stgLVariables[var][it.getStgId()];
    }
  } else if(k == a->end(2) - a->m_noGhostLayers) {
    // 1st layer
    for(MInt var = PV->noVariables; var < noSTGVariables; ++var) {
      m_stgLVariables[var][it.getNghbrStg(5)] = m_stgLVariables[var][it.getStgId()];
    }
  }

  const MInt j = it.getijk(1);
  if(j == a->start(1) + 1) {
    for(MInt var = PV->noVariables; var < noSTGVariables; ++var) {
      m_stgLVariables[var][it.getNghbrStg(2)] = m_stgLVariables[var][it.getStgId()];
    }
  } else if(j == a->end(1) - a->m_noGhostLayers) {
    for(MInt var = PV->noVariables; var < noSTGVariables; ++var) {
      m_stgLVariables[var][it.getNghbrStg(3)] = m_stgLVariables[var][it.getStgId()];
    }
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::extrapolateToGX(typename Accessor::nDim_citerator it) {
  const MInt cellIdG1 = a->getCellId(it);
  const MInt cellIdA1 = a->getNghbr(it, 1); // express this by means of face
  const MInt cellIdG2 = a->getNghbr(it, 0);
  //  const MInt cellIdG2 = cellIndex(0,j,k);
  // extrapolate into second ghost cell
  for(MInt var = 0; var < PV->noVariables; var++) {
    LES->a_pvariable(cellIdG2, var) = F2 * LES->a_pvariable(cellIdG1, var) - LES->a_pvariable(cellIdA1, var);
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_FINITE_VOLUME, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::apply() {
  TRACE();

  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;

  if(m_initialRange || globalTimeStep < m_solver->m_stgStartTimeStep) {
    for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
      const MInt cellIdG1 = a->getCellId(it);
      const MInt IBC = a->getStgId(it);

      if(m_cutOff) {
        LES->a_pvariable(cellIdG1, PV->RHO) = m_stgLVariables[STG::AVG_RHO][IBC];
        LES->a_pvariable(cellIdG1, PV->U) = m_stgLVariables[STG::AVG_U][IBC];
        LES->a_pvariable(cellIdG1, PV->V) = m_stgLVariables[STG::AVG_V][IBC];
        LES->a_pvariable(cellIdG1, PV->W) = m_stgLVariables[STG::AVG_W][IBC];
      } else {
        const MInt ghostCellId = a->getGhostIdFromStgId(IBC);
        // Density
        LES->a_pvariable(ghostCellId, PV->RHO) =
            2.0 * m_stgLVariables[STG::AVG_RHO][IBC] - LES->a_pvariable(cellIdG1, PV->RHO);
        // Velocitiy
        LES->a_pvariable(ghostCellId, PV->U) =
            2.0 * m_stgLVariables[STG::AVG_U][IBC] - LES->a_pvariable(cellIdG1, PV->U);
        LES->a_pvariable(ghostCellId, PV->V) =
            2.0 * m_stgLVariables[STG::AVG_V][IBC] - LES->a_pvariable(cellIdG1, PV->V);
        LES->a_pvariable(ghostCellId, PV->W) =
            2.0 * m_stgLVariables[STG::AVG_W][IBC] - LES->a_pvariable(cellIdG1, PV->W);
        // Pressure
        LES->a_pvariable(ghostCellId, PV->P) =
            2.0 * m_stgLVariables[STG::AVG_P][IBC] - LES->a_pvariable(cellIdG1, PV->P);
      }
    }

  } else {
    for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
      const MInt cellIdG1 = a->getCellId(it);
      const MInt IBC = a->getStgId(it);
      // TODO: CHECK THIS LATER
      const MInt cellIdA1 = cellIdG1;
      // if(m_cutOff) {
      //   const MInt cellIdA1 = a->getNghbr(it, 1); // cellIdG1; // a->getNghbr(it, 1); //express this by means of face
      // }

      // TODO: especially check if this is true if stg is applied at +x
      MFloat dxidx = 1;
      MFloat dxidy = 0;
      MFloat dxidz = 0;

      MFloat gradxi = -F1 / sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

      MFloat dxHelp = dxidx * gradxi;
      MFloat dyHelp = dxidy * gradxi;
      MFloat dzHelp = dxidz * gradxi;

      // compute Mach number perpendicular to inlet
      const MFloat rhoBC = LES->a_pvariable(cellIdG1, PV->RHO);
      const MFloat pBC = LES->a_pvariable(cellIdG1, PV->P);
      const MFloat fRhoBC = F1 / rhoBC;
      const MFloat aBC = sqrt(gamma * pBC * fRhoBC);
      const MFloat uBC = LES->a_pvariable(cellIdG1, PV->U);
      const MFloat vBC = LES->a_pvariable(cellIdG1, PV->V);
      const MFloat wBC = LES->a_pvariable(cellIdG1, PV->W);

      const MFloat maBC = (dxHelp * uBC + dyHelp * vBC + dzHelp * wBC) / aBC;
      /////////

      /**This is where the fluctuations should be added to velocity!**/

      // get mean values from the rans
      const MFloat rhoRANS = m_stgLVariables[STG::AVG_RHO][IBC];
      const MFloat uRANS = m_stgLVariables[STG::AVG_U][IBC];
      const MFloat vRANS = m_stgLVariables[STG::AVG_V][IBC];
      const MFloat wRANS = m_stgLVariables[STG::AVG_W][IBC];
      const MFloat pRANS = m_stgLVariables[STG::AVG_P][IBC];

      // fluctuation values from the STG
      const MFloat u_prime = m_stgLVariables[STG::FLUC_U][IBC];
      const MFloat v_prime = m_stgLVariables[STG::FLUC_V][IBC];
      const MFloat w_prime = m_stgLVariables[STG::FLUC_W][IBC];

      // superpose onto mean RANS variables
      const MFloat uSTG = std::max(uRANS + u_prime, epsl);
      const MFloat vSTG = vRANS + v_prime;
      const MFloat wSTG = wRANS + w_prime;

      // compute correct density
      const MFloat u9a = LES->UInfinity();
      const MFloat u9ff = u_prime;
      const MFloat alok = sqrt(gamma * LES->PInfinity() / LES->rhoInfinity());
      const MFloat flucc =
          u9ff / u9a * POW2((LES->UInfinity() / alok)) * gammaMinusOne * m_stgLVariables[STG::AVG_RHO][IBC];
      const MFloat zdir = flucc / std::max(fabs(flucc), 0.0000001);
      const MFloat rhoSTG = rhoRANS + zdir * std::min(fabs(flucc), 0.1 * rhoRANS);

      // get field values inside the integration domain
      const MFloat pField = LES->a_pvariable(cellIdA1, PV->P);
      const MFloat rhoField = LES->a_pvariable(cellIdA1, PV->RHO);
      const MFloat uField = LES->a_pvariable(cellIdA1, PV->U);
      const MFloat vField = LES->a_pvariable(cellIdA1, PV->V);
      const MFloat wField = LES->a_pvariable(cellIdA1, PV->W);
      const MFloat aField = sqrt(gamma * pField / rhoField);

      /////////////////////////////////////////////////
      //////////// SUBSONIC PART //////////////////////
      /////////////////////////////////////////////////
      const MFloat pSub =
          F1B2
          * (pField + pRANS
             + rhoField * aField * (+dxHelp * (uField - uSTG) + dyHelp * (vField - vSTG) + dzHelp * (wField - wSTG)));
      const MFloat rhoSub = rhoSTG + (pSub - pRANS) / (POW2(aField));
      const MFloat rhoSubHelp = (pSub - pRANS) / (rhoField * aField);

      // Multiply velocities with density
      const MFloat uSub = uSTG + dxHelp * rhoSubHelp;
      const MFloat vSub = vSTG + dyHelp * rhoSubHelp;
      const MFloat wSub = wSTG + dzHelp * rhoSubHelp;

      /////////////////////////////////////////////////
      //////////// SUPERSONIC PART ////////////////////
      /////////////////////////////////////////////////
      const MFloat rhoSup = rhoSTG;
      const MFloat uSup = uSTG;
      const MFloat vSup = vSTG;
      const MFloat wSup = wSTG;
      const MFloat pSup = pRANS;

      //////////////////////////////////////////////////
      /////////// SUB/SUP INTERPOLATION ////////////////
      //////////////////////////////////////////////////

      // by default the subsonic formulation is used
      // switch on "stgSubSup" to get the mixed formulation
      // or "stgSupersonic" to use the pure supersonic formulation
      MFloat xSub = F1;
      MFloat xSup = F0;

      if(m_stgSubSup) {
        const MFloat maBCAbs = fabs(maBC);
        const MFloat alpha = 14.0;
        const MFloat b = 0.95;
        const MFloat count = alpha * (maBCAbs - b);
        const MFloat denom = (F1 - 0.99 * b) * maBCAbs + b;
        const MFloat ratio = count / denom;
        const MFloat wfun = F1B2 * (F1 + tanh(ratio) / tanh(alpha));

        xSub = fabs(wfun - F1);
        xSup = fabs(wfun);
      } else if(m_stgSupersonic) {
        xSub = F0;
        xSup = F1;
      }

      const MFloat rho_target = rhoSub * xSub + rhoSup * xSup;
      const MFloat u_target = uSub * xSub + uSup * xSup;
      const MFloat v_target = vSub * xSub + vSup * xSup;
      const MFloat w_target = wSub * xSub + wSup * xSup;
      const MFloat p_target = pSub * xSub + pSup * xSup;

      if(m_cutOff) {
        LES->a_pvariable(cellIdG1, PV->RHO) = rho_target;
        LES->a_pvariable(cellIdG1, PV->U) = u_target;
        LES->a_pvariable(cellIdG1, PV->V) = v_target;
        LES->a_pvariable(cellIdG1, PV->W) = w_target;
      } else {
        const MInt ghostCellId = a->getGhostIdFromStgId(IBC);
        LES->a_pvariable(ghostCellId, PV->RHO) = 2.0 * rho_target - LES->a_pvariable(cellIdG1, PV->RHO);
        LES->a_pvariable(ghostCellId, PV->U) = 2.0 * u_target - LES->a_pvariable(cellIdG1, PV->U);
        LES->a_pvariable(ghostCellId, PV->V) = 2.0 * v_target - LES->a_pvariable(cellIdG1, PV->V);
        LES->a_pvariable(ghostCellId, PV->W) = 2.0 * w_target - LES->a_pvariable(cellIdG1, PV->W);
        LES->a_pvariable(ghostCellId, PV->P) = 2.0 * p_target - LES->a_pvariable(cellIdG1, PV->P);
      }
    }
  }
}


template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
template <class _, std::enable_if_t<SolverTypeL == MAIA_STRUCTURED, _*>>
void MSTG<nDim, SolverTypeR, SolverTypeL>::apply() {
  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;

  for(typename Accessor::nDim_citerator it = a->iterateB1(); it != a->iterateB1_nDim_citerator_end(); ++it) {
    const MInt cellIdG1 = a->getCellId(it);
    const MInt IBC = a->getStgId(it);
    const MInt cellIdA1 = a->getNghbr(it, 1); // express this by means of face


    MFloat dxidx = m_solver->m_cells->cellMetrics[cellIdA1][0];
    MFloat dxidy = m_solver->m_cells->cellMetrics[cellIdA1][1];
    MFloat dxidz = m_solver->m_cells->cellMetrics[cellIdA1][2];

    MFloat gradxi = -F1 / sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

    MFloat dxHelp = dxidx * gradxi;
    MFloat dyHelp = dxidy * gradxi;
    MFloat dzHelp = dxidz * gradxi;

    // compute Mach number perpendicular to inlet
    const MFloat rhoBC = LES->a_pvariable(cellIdG1, PV->RHO);
    const MFloat pBC = LES->a_pvariable(cellIdG1, PV->P);
    const MFloat fRhoBC = F1 / rhoBC;
    const MFloat aBC = sqrt(gamma * pBC * fRhoBC);
    const MFloat uBC = LES->a_pvariable(cellIdG1, PV->U);
    const MFloat vBC = LES->a_pvariable(cellIdG1, PV->V);
    const MFloat wBC = LES->a_pvariable(cellIdG1, PV->W);

    const MFloat maBC = (dxHelp * uBC + dyHelp * vBC + dzHelp * wBC) / aBC;
    /////////

    /**This is where the fluctuations should be added to velocity!**/

    // get mean values from the rans
    const MFloat rhoRANS = m_stgLVariables[STG::AVG_RHO][IBC];
    const MFloat uRANS = m_stgLVariables[STG::AVG_U][IBC];
    const MFloat vRANS = m_stgLVariables[STG::AVG_V][IBC];
    const MFloat wRANS = m_stgLVariables[STG::AVG_W][IBC];
    const MFloat pRANS = m_stgLVariables[STG::AVG_P][IBC];

    // fluctuation values from the STG
    const MFloat u_prime = m_stgLVariables[STG::FLUC_U][IBC];
    const MFloat v_prime = m_stgLVariables[STG::FLUC_V][IBC];
    const MFloat w_prime = m_stgLVariables[STG::FLUC_W][IBC];

    // superpose onto mean RANS variables
    const MFloat uSTG = std::max(uRANS + u_prime, epsl);
    const MFloat vSTG = vRANS + v_prime;
    const MFloat wSTG = wRANS + w_prime;

    // compute correct density
    const MFloat u9a = LES->UInfinity();
    const MFloat u9ff = u_prime;
    const MFloat alok = sqrt(gamma * LES->PInfinity() / LES->rhoInfinity());
    const MFloat flucc =
        u9ff / u9a * POW2((LES->UInfinity() / alok)) * gammaMinusOne * m_stgLVariables[STG::AVG_RHO][IBC];
    const MFloat zdir = flucc / std::max(fabs(flucc), 0.0000001);
    const MFloat rhoSTG = rhoRANS + zdir * std::min(fabs(flucc), 0.1 * rhoRANS);

    // get field values inside the integration domain
    const MFloat pField = LES->a_pvariable(cellIdA1, PV->P);
    const MFloat rhoField = LES->a_pvariable(cellIdA1, PV->RHO);
    const MFloat uField = LES->a_pvariable(cellIdA1, PV->U);
    const MFloat vField = LES->a_pvariable(cellIdA1, PV->V);
    const MFloat wField = LES->a_pvariable(cellIdA1, PV->W);
    const MFloat aField = sqrt(gamma * pField / rhoField);

    /////////////////////////////////////////////////
    //////////// SUBSONIC PART //////////////////////
    /////////////////////////////////////////////////
    const MFloat pSub =
        F1B2
        * (pField + pRANS
           + rhoField * aField * (+dxHelp * (uField - uSTG) + dyHelp * (vField - vSTG) + dzHelp * (wField - wSTG)));
    const MFloat rhoSub = rhoSTG + (pSub - pRANS) / (POW2(aField));
    const MFloat rhoSubHelp = (pSub - pRANS) / (rhoField * aField);

    // Multiply velocities with density
    const MFloat uSub = uSTG + dxHelp * rhoSubHelp;
    const MFloat vSub = vSTG + dyHelp * rhoSubHelp;
    const MFloat wSub = wSTG + dzHelp * rhoSubHelp;

    /////////////////////////////////////////////////
    //////////// SUPERSONIC PART ////////////////////
    /////////////////////////////////////////////////
    const MFloat rhoSup = rhoSTG;
    const MFloat uSup = uSTG;
    const MFloat vSup = vSTG;
    const MFloat wSup = wSTG;
    const MFloat pSup = pRANS;

    //////////////////////////////////////////////////
    /////////// SUB/SUP INTERPOLATION ////////////////
    //////////////////////////////////////////////////

    // by default the subsonic formulation is used
    // switch on "stgSubSup" to get the mixed formulation
    // or "stgSupersonic" to use the pure supersonic formulation
    MFloat xSub = F1;
    MFloat xSup = F0;

    if(m_stgSubSup) {
      const MFloat maBCAbs = fabs(maBC);
      const MFloat alpha = 14.0;
      const MFloat b = 0.95;
      const MFloat count = alpha * (maBCAbs - b);
      const MFloat denom = (F1 - 0.99 * b) * maBCAbs + b;
      const MFloat ratio = count / denom;
      const MFloat wfun = F1B2 * (F1 + tanh(ratio) / tanh(alpha));

      xSub = fabs(wfun - F1);
      xSup = fabs(wfun);
    } else if(m_stgSupersonic) {
      xSub = F0;
      xSup = F1;
    }

    LES->a_pvariable(cellIdG1, PV->RHO) = rhoSub * xSub + rhoSup * xSup;
    LES->a_pvariable(cellIdG1, PV->U) = uSub * xSub + uSup * xSup;
    LES->a_pvariable(cellIdG1, PV->V) = vSub * xSub + vSup * xSup;
    LES->a_pvariable(cellIdG1, PV->W) = wSub * xSub + wSup * xSup;
    LES->a_pvariable(cellIdG1, PV->P) = pSub * xSub + pSup * xSup;


    //////////////////////////////////////////////////
    //////// EXTRAPOLATE TO SECOND GC ////////////////
    //////////////////////////////////////////////////
    extrapolateToGX(it);
    //    }
  }
}


template class MSTG<3, MAIA_FINITE_VOLUME, MAIA_FINITE_VOLUME>;
template class MSTG<3, MAIA_STRUCTURED, MAIA_FINITE_VOLUME>;
