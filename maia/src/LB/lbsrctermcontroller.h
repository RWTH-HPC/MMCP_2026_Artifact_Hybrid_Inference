// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBSRCTERMCONTROLLER_H
#define LBSRCTERMCONTROLLER_H

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// DO NOT CHANGE THIS FILE UNLESS YOU KNOW WHAT YOU ARE DOING !
//
// If you want to add a new source term you can do so by directly adding it into
// the lbsrcterm.cpp (and nothing else anywhere). Alternatively, you might want
// to create a new file lbsrcterm_xy.cpp. This needs to include
// lbsrctermcontroller.h and lbsrcterm.h.
// Examples can be found in lbsrcterm.cpp.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <functional>
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "UTIL/debug.h"
#include "lbsrcterm.h"

template <MInt nDim, MInt nDist, class SysEqn>
class LbSolverDxQy;

namespace maia::lb {

/** \brief  Front-end to create source term objects
 *  \author Miro Gondrum
 *  \date   01.02.2022
 *
 * This class holds a list (function_reg) where all source terms are
 * automatically registered if they are derived from LbSrcTerm and include
 * in their implementation the line:
 *
 * static LbRegSrcTerm<LbSrcTermNAME> reg("unique_srcTerm_name");
 *
 * To create a new source term, simply copy one of the existing ones and change
 * the name. It will be found automatically, no includes are necessary.
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTermFactory {
  using srcTermConstructor = std::function<LbSrcTerm<nDim, nDist, SysEqn>*(LbSolverDxQy<nDim, nDist, SysEqn>*)>;

 public:
  /**  \brief  Creates a static instance of LbSrcTermFactory
   *   \author Miro Gondrum
   *   \date   01.02.2022
   *   \return Pointer to a static instance of LbSrcTermFactory
   */
  static LbSrcTermFactory* instance();

  /**  \brief  Adds a new LbSrcTerm object to the function registry.
   *   \author Miro Gondrum
   *   \date   01.02.2022
   *   \param[in]  p_name         The name of the source term
   *   \param[in]  fact_function  The source term constructor
   */
  void reg_function(const MString& p_name, srcTermConstructor fact_function) { m_function_reg[p_name] = fact_function; }

  /**  \brief  Create an lb source term object and return a pointer to it
   *   \author Miro Gondrum
   *   \date   01.02.2022
   *   \param[in]  p_name  name of the source term (as stored in the registry)
   *   \return     Pointer to the source term object
   */
  LbSrcTerm<nDim, nDist, SysEqn>* create_srcTerm(const MString& p_name, LbSolverDxQy<nDim, nDist, SysEqn>* p_solver);

 private:
  LbSrcTermFactory(){};
  std::map<std::string, srcTermConstructor> m_function_reg; ///< data structure to associate srcTerm names with their
                                                            /// implementations
};

/** \brief  Class for registering a source term to the factory registry
 *  \author Miro Gondrum
 *  \date   01.02.2022
 *  \tparam  T The type T should be the Class of the source term to be added.
 */
template <class T>
class LbRegSrcTerm {
 public:
  /**  \brief  Constructor, does the registering
   *   \author Miro Gondrum
   *   \date   01.02.2022
   *   \param[in]  p_name  Name of the source term
   */
  LbRegSrcTerm(const MString& p_name) {
    LbSrcTermFactory<T::nDim, T::nDist, typename T::SysEqn>::instance()->reg_function(
        p_name,
        [](LbSolverDxQy<T::nDim, T::nDist, typename T::SysEqn> * p_solver)
            -> LbSrcTerm<T::nDim, T::nDist, typename T::SysEqn>* { return new T(p_solver); });
  }
};

/**  \brief  Front-end to control all source terms in a wrapping manner
 *   \author Miro Gondrum
 *   \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTermController {
 public:
  LbSrcTermController(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : m_solver(p_solver){};

  void init();
  void initSrcTerms();
  void addSrcTerm(const MString& p_name);
  void apply_preCollision();
  void apply_postCollision();
  void apply_postPropagation();

 private:
  LbSolverDxQy<nDim, nDist, SysEqn>* m_solver;
  MInt m_noSrcTerms{};
  std::vector<std::unique_ptr<LbSrcTerm<nDim, nDist, SysEqn>>> m_srcTerms;
};

/**  \brief  Initialize the source term controller
 *   \author Miro Gondrum
 *   \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermController<nDim, nDist, SysEqn>::init() {
  TRACE();
  std::stringstream ss;
  ss << "- Init source terms:" << std::endl;
  m_noSrcTerms = 0;
  if(Context::propertyExists("lbSrcTerms", m_solver->m_solverId)) {
    m_noSrcTerms = Context::propertyLength("lbSrcTerms");
  }
  m_srcTerms.clear();
  ss << " no source terms =" << m_noSrcTerms << std::endl;
  ss << "  ";
  if(m_noSrcTerms > 0) {
    for(MInt i = 0; i < m_noSrcTerms; i++) {
      const MString srcTermName = Context::getBasicProperty<MString>("lbSrcTerms", AT_, &srcTermName, i);
      m_srcTerms.emplace_back(LbSrcTermFactory<nDim, nDist, SysEqn>::instance()->create_srcTerm(srcTermName, m_solver));
      ss << srcTermName << " ";
    }
    ss << std::endl;
    m_log << ss.str();
    if(m_solver->domainId() == 0) std::cout << ss.str();
  }
}

/**  \brief  Initialize the source term controller
 *   \author Julian Vorspohl
 *   \date   01.05.2024
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermController<nDim, nDist, SysEqn>::initSrcTerms() {
  // Call init routine for each source term instance created
  for(auto& srcTerm : m_srcTerms) {
    srcTerm->init();
  }
}

/**  \brief  Add a source term to the controller by its property tag
 *   \author Miro Gondrum
 *   \date   01.02.2022
 *   \param[in]  p_name The name of the source term
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermController<nDim, nDist, SysEqn>::addSrcTerm(const MString& p_name) {
  TRACE();
  m_srcTerms.emplace_back(LbSrcTermFactory<nDim, nDist, SysEqn>::instance()->create_srcTerm(p_name, m_solver));
}

/**  \brief  Call the pre collision routines of all source terms
 *   \author Miro Gondrum
 *   \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermController<nDim, nDist, SysEqn>::apply_preCollision() {
  TRACE();
  for(auto& srcTerm : m_srcTerms) {
    srcTerm->apply_preCollision();
  }
}

/**  \brief  Call the post collision routines of all source terms
 *   \author Miro Gondrum
 *   \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermController<nDim, nDist, SysEqn>::apply_postCollision() {
  TRACE();
  for(auto& srcTerm : m_srcTerms) {
    srcTerm->apply_postCollision();
  }
}

/**  \brief  Call the post collision routines of all source terms
 *   \author Julian Vorspohl
 *   \date   01.05.2024
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermController<nDim, nDist, SysEqn>::apply_postPropagation() {
  TRACE();
  for(auto& srcTerm : m_srcTerms) {
    srcTerm->apply_postPropagation();
  }
}

} // namespace maia::lb

#endif
