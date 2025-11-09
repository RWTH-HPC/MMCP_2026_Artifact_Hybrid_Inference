// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


#ifndef EXECUTIONRECIPE_H_
#define EXECUTIONRECIPE_H_

#include <algorithm>
#include <vector>

#include "IO/context.h"
#include "UTIL/functions.h"

#include "COUPLER/coupling.h"
#include "solver.h"

/// Base recipe provides public interface to Application
class ExecutionRecipe {
 public:
  ExecutionRecipe(std::vector<std::unique_ptr<Solver>>* const solvers,
                  std::vector<std::unique_ptr<Coupling>>* const couplers)
    : m_solvers(solvers), m_couplers(couplers) {
    // This function initializes function pointers to the preTimeStep, solutionStep
    // and postTimeStep functions of each solver
    initFunctionPointers();
  }

  // public functions calle from the run-loop:

  virtual void preTimeStep() final;
  virtual void timeStep();
  virtual void postTimeStep() final;

  virtual void preCouple() final;
  virtual void postCouple() final;

  virtual MBool updateCallOrder() { return true; };

  MInt a_step() const { return m_step; }
  MInt& a_step() { return m_step; }

  MBool callAdaptation() const { return m_allowAdaptation; }

 protected:
  void readCallOrder();
  void initFunctionPointers();

  void setSolverStatus(MInt, MBool);
  void setCouplerStatus(const MInt couplerId, const MBool active);

  MInt noSolvers() const { return m_solvers->size(); }
  MInt noCouplers() const { return m_couplers->size(); }

  void nextStep() { ++m_step; }

  MBool solverOrder(const MInt solverId) const { return m_solverOrder[m_step][solverId]; }

  MBool couplerOrder(const MInt couplerId) const { return m_couplerOrder[m_step][couplerId]; }

  void setAdaptation() { m_allowAdaptation = m_adaptationOrder[m_step]; }

  MInt maxNoSteps() const { return m_maxNoSteps; }

  MBool solutionStep(const MInt solverId) { return m_solutionStep[solverId](); }

  void preSolutionStep(const MInt solverId, const MInt mode) { return m_preSolutionStep[solverId](mode); }

  MBool postSolutionStep(const MInt solverId) { return m_postSolutionStep[solverId](); }

  void subCouple(const MInt couplerId, const MInt step, const MInt solverId, std::vector<MBool>& solverCompleted) {
    m_subCouple[couplerId](step, solverId, solverCompleted);
  }

  const std::vector<std::unique_ptr<Solver>>* a_solvers() const { return m_solvers; }

  const std::vector<std::unique_ptr<Coupling>>* a_couplers() const { return m_couplers; }

  MInt m_maxSolutionIteration = 1;

  std::map<MInt, MInt> m_swapSolverIds{};

  MInt swapedSolverId(const MInt oldSolverId) {
    // switch solver order during execution, this becomes necessary during intraStep-coupling,
    // when however one solver depends on a certain solver initialisation, thus prescrining the
    // solverId-oder, but during the time-step the order needs to be the other way around!

    ASSERT(!m_swapSolverIds.empty(), "");

    MInt solverId = oldSolverId;
    auto it = m_swapSolverIds.find(solverId);
    if(it != m_swapSolverIds.end()) {
      solverId = it->second;
    }

    ASSERT(solverOrder(oldSolverId) == solverOrder(solverId), "");

    return solverId;
  }

  void startLoadTimer(const MInt solverId) { m_solvers->at(solverId)->startLoadTimer(AT_); }

  void stopLoadTimer(const MInt solverId) { m_solvers->at(solverId)->stopLoadTimer(AT_); }

  MBool solverIsActive(const MInt solverId) { return m_solvers->at(solverId)->isActive(); }

 private:
  const std::vector<std::unique_ptr<Solver>>* const m_solvers;
  const std::vector<std::unique_ptr<Coupling>>* const m_couplers;

  std::vector<std::function<void()>> m_preStep;
  std::vector<std::function<void(const MInt)>> m_preSolutionStep;
  std::vector<std::function<MBool()>> m_solutionStep;
  std::vector<std::function<MBool()>> m_postSolutionStep;
  std::vector<std::function<void()>> m_postStep;

  std::vector<std::function<void(const MInt)>> m_preCouple;
  // std::vector< std::function< void(const MInt) > > m_preSolutionCouple;
  std::vector<std::function<void(const MInt, const MInt, std::vector<MBool>&)>> m_subCouple;
  // std::vector< std::function< void(const MInt) > > m_postSolutionCouple;
  std::vector<std::function<void(const MInt)>> m_postCouple;


  MBool m_allowAdaptation = true;

  MInt m_step = -1;
  MInt m_maxNoSteps = 1;
  std::vector<std::vector<MBool>> m_solverOrder{};
  std::vector<std::vector<MBool>> m_couplerOrder{};
  std::vector<MBool> m_adaptationOrder{};
};


/** \brief: Reads the call order of solvers, couplers and adaptation

    For some multisolver applications the solvers can not be executed simultaneously.
    E.g. in level-set moving boundary problems the new level-set solution needs
    to be known before the fvMbSolver can perform its time step.
    Therefore, a callOrder is specified in the propertiesFile for the solvers, couplers
    and the adaptation. According to the callOrder solvers and couplers are set to
    active or empty and the adaptation can be turned on and off.
    The callOrder is updated in each step of the advance time step iteration inside
    the time step loop ( updateCallOrder() ).
    In the first step of level-set moving boundary problems the LS solver is set active
    and the fvMb solver is set to idle. In the subsequent step the LS solver is
    deactivated, while the fvMb solver is activated.

  * \author Thomas Hoesgen
  * \date 01/2020
  */
void ExecutionRecipe::readCallOrder() {
  TRACE();

  a_step() = 0;

  // TODO labels:COUPLER,DOC,toenhance allow a default call order to be defined by the couplers (when there is only one
  // coupler)
  // @ansgar this would make sense e.g. for the multilevel interpolation coupler!

  m_maxNoSteps = 1;
  m_maxNoSteps = Context::getBasicProperty<MInt>("recipeMaxNoSteps", AT_, &m_maxNoSteps);


  m_maxSolutionIteration = 1;
  m_maxSolutionIteration = Context::getBasicProperty<MInt>("maxIterations", AT_, &m_maxSolutionIteration);

  m_solverOrder.resize(m_maxNoSteps, std::vector<MBool>(noSolvers()));
  m_couplerOrder.resize(m_maxNoSteps, std::vector<MBool>(noCouplers()));
  m_adaptationOrder.resize(m_maxNoSteps);

  // Read solverOrder_x properties defining when the solvers should be executed, e.g. solverOrder_1 =
  // [1, 0] to execute the solver #1 in the first step but not in the second step
  for(MInt solver = 0; solver < noSolvers(); solver++) {
    const MString propName = "solverOrder_" + std::to_string(solver);
    Context::assertPropertyLength(propName, m_maxNoSteps);

    for(MInt step = 0; step < m_maxNoSteps; step++) {
      m_solverOrder[step][solver] = (MBool)Context::getBasicProperty<MInt>(propName, AT_, step);
    }
  }

  // Read couplerOrder_x properties defining when the coupler should be executed, e.g. [0, 1] to
  // execute the coupler only in the second step
  for(MInt coupler = 0; coupler < noCouplers(); coupler++) {
    const MString propName = "couplerOrder_" + std::to_string(coupler);
    Context::assertPropertyLength(propName, m_maxNoSteps);

    for(MInt step = 0; step < m_maxNoSteps; step++) {
      m_couplerOrder[step][coupler] = (MBool)Context::getBasicProperty<MInt>(propName, AT_, step);
    }
  }

  for(MInt step = 0; step < m_maxNoSteps; step++) {
    // Property adaptationOrder is an array of length m_maxNoSteps - e.g. [0,1]
    Context::assertPropertyLength("adaptationOrder", m_maxNoSteps);
    m_adaptationOrder[step] = (MBool)Context::getBasicProperty<MInt>("adaptationOrder", AT_, step);
  }

  // Set solver and coupler statuses for recipe step 0
  for(MInt blck = 0; blck < noSolvers(); blck++) {
    setSolverStatus(blck, m_solverOrder[a_step()][blck]);
  }

  for(MInt cpl = 0; cpl < noCouplers(); cpl++) {
    setCouplerStatus(cpl, m_couplerOrder[a_step()][cpl]);
  }

  m_allowAdaptation = m_adaptationOrder[a_step()];

  // Print execution recipe call order to m_log
  m_log << " === Execution recipe: call order with " << m_maxNoSteps << " steps" << std::endl;
  m_log << "step:      ";
  for(MInt step = 0; step < m_maxNoSteps; step++) {
    m_log << " " << step;
  }
  m_log << std::endl;

  for(MInt solver = 0; solver < noSolvers(); solver++) {
    m_log << "solver #" << solver << ":  ";
    for(MInt step = 0; step < m_maxNoSteps; step++) {
      m_log << " " << m_solverOrder[step][solver];
    }
    m_log << std::endl;
  }
  for(MInt cpl = 0; cpl < noCouplers(); cpl++) {
    m_log << "coupler #" << cpl << ":";
    for(MInt step = 0; step < m_maxNoSteps; step++) {
      m_log << " " << m_couplerOrder[step][cpl];
    }
    m_log << std::endl;
  }

  for(MInt blck = 0; blck < noSolvers(); blck++) {
    setSolverStatus(blck, solverOrder(blck));
  }

  for(MInt cpl = 0; cpl < noCouplers(); cpl++) {
    setCouplerStatus(cpl, couplerOrder(cpl));
  }

  setAdaptation();

  m_swapSolverIds.clear();
  if(Context::propertyExists("swapSolverSolutionStep")) {
    MInt size = Context::propertyLength("swapSolverSolutionStep") / 2;
    for(MInt i = 0; i < size; i++) {
      const MInt oldId = Context::getBasicProperty<MInt>("swapSolverSolutionStep", AT_, i);
      const MInt newId = Context::getBasicProperty<MInt>("swapSolverSolutionStep", AT_, i + 1);
      m_swapSolverIds.insert(std::make_pair(oldId, newId));
      m_swapSolverIds.insert(std::make_pair(newId, oldId));
    }
  }
}

/** \brief: Wrapper function to set solvers active or idle
 * \author Thomas Hoesgen
 * \date 01/2020
 */
void ExecutionRecipe::setSolverStatus(MInt solverId, MBool active) {
  TRACE();

  // Tell the solver its current status such that a coupler can ask the solver if its currently active
  // inside of the execution recipe
  m_solvers->at(solverId)->setSolverStatus(active);


  if(active) {
    m_preStep[solverId] = [=]() { m_solvers->at(solverId)->preTimeStep(); };
    m_preSolutionStep[solverId] = [=](MInt mode) { m_solvers->at(solverId)->preSolutionStep(mode); };
    m_solutionStep[solverId] = [=]() { return m_solvers->at(solverId)->solutionStep(); };
    m_postSolutionStep[solverId] = [=]() { return m_solvers->at(solverId)->postSolutionStep(); };
    m_postStep[solverId] = [=]() { m_solvers->at(solverId)->postTimeStep(); };
  } else {
    // set to empty lambda
    m_preStep[solverId] = []() {};
    m_preSolutionStep[solverId] = [](MInt) {};
    m_solutionStep[solverId] = []() { return true; };
    m_postSolutionStep[solverId] = []() { return true; };
    m_postStep[solverId] = []() {};
  }
}

/** \brief: Wrapper function to set couplers active or empty
 * \author Thomas Hoesgen
 * \date 01/2020
 */
void ExecutionRecipe::setCouplerStatus(const MInt couplerId, const MBool active) {
  TRACE();

  if(active) {
    m_preCouple[couplerId] = [=](MInt i) { m_couplers->at(couplerId)->preCouple(i); };
    m_subCouple[couplerId] = [=](const MInt i, const MInt s, std::vector<MBool>& sc) {
      return m_couplers->at(couplerId)->subCouple(i, s, sc);
    };
    m_postCouple[couplerId] = [=](MInt i) { m_couplers->at(couplerId)->postCouple(i); };
  } else {
    // set to empty lambda
    m_preCouple[couplerId] = [](const MInt) {};
    m_subCouple[couplerId] = [](const MInt, const MInt, std::vector<MBool>&) { return true; };
    m_postCouple[couplerId] = [](const MInt) {};
  }
}

/** \brief: Initialize the vector containing function pointers to preTimeStep,
            solutionStep and postTimestep as well as preCouple, subCouple and
            postCouple. All are set active.
  * \author Thomas Hoesgen
  * \date 01/2020
  */
void ExecutionRecipe::initFunctionPointers() {
  TRACE();

  for(MInt i = 0; i < noSolvers(); i++) {
    m_preStep.emplace_back([=]() { m_solvers->at(i)->preTimeStep(); });
    m_preSolutionStep.emplace_back([=](MInt b) { m_solvers->at(i)->preSolutionStep(b); });
    m_solutionStep.emplace_back([=]() { return m_solvers->at(i)->solutionStep(); });
    m_postSolutionStep.emplace_back([=]() { return m_solvers->at(i)->postSolutionStep(); });
    m_postStep.emplace_back([=]() { m_solvers->at(i)->postTimeStep(); });
  }

  for(MInt i = 0; i < noCouplers(); i++) {
    m_preCouple.emplace_back([=](MInt v) { m_couplers->at(i)->preCouple(v); });
    m_subCouple.emplace_back(
        [=](const MInt v, const MInt s, std::vector<MBool>& sc) { return m_couplers->at(i)->subCouple(v, s, sc); });
    m_postCouple.emplace_back([=](MInt v) { m_couplers->at(i)->postCouple(v); });
  }
}

/** \brief: Calls each solvers preTimeStep - might be empty
 * \author Thomas Hoesgen
 * \date 01/2020
 */
inline void ExecutionRecipe::preTimeStep() {
  TRACE();
  for(auto&& solver : *m_solvers) {
    if(!solver->isActive()) {
      continue; // skip inactive solvers
    }

    solver->startLoadTimer(AT_);
    m_preStep[solver->solverId()]();
    solver->stopLoadTimer(AT_);
  }
}

/** \brief: Single solver time step function.
            Calls solutionStep() of the specific solver
  * \author Thomas Hoesgen
  * \date 01/2020
  */
inline void ExecutionRecipe::timeStep() {
  TRACE();

  // For a single solver run, there is obviously only one solver!
  for(auto&& solver : *m_solvers) {
    if(!solver->isActive()) {
      continue; // skip inactive solvers
    }

    // Outer loop over substep iterations
    MBool completed = false;
    while(!completed) {
      // Solution step
      solver->startLoadTimer(AT_);
      completed = m_solutionStep[solver->solverId()]();
      solver->stopLoadTimer(AT_);
    }
  }
}

/** \brief: Calls each solvers postTimeStep - might be empty
 * \author Thomas Hoesgen
 * \date 01/2020
 */
void ExecutionRecipe::postTimeStep() {
  TRACE();
  for(auto&& solver : *m_solvers) {
    const MInt solverId = m_swapSolverIds.empty() ? solver->solverId() : swapedSolverId(solver->solverId());
    if(!solverIsActive(solverId)) {
      continue; // skip inactive solvers
    }
    startLoadTimer(solverId);
    m_postStep[solverId]();
    stopLoadTimer(solverId);
  }
}

/** \brief: Calls each couplers preCouple - might be empty
 * \author Thomas Hoesgen
 * \date 01/2020
 */
inline void ExecutionRecipe::preCouple() {
  TRACE();

  for(auto&& coupler : *m_couplers) {
    coupler->startLoadTimer(AT_);
    m_preCouple[coupler->couplerId()](a_step());
    coupler->stopLoadTimer(AT_);
  }
}

/** \brief: Calls each couplers postCouple - might be empty
 * \author Thomas Hoesgen
 * \date 01/2020
 */
void ExecutionRecipe::postCouple() {
  TRACE();

  for(auto&& coupler : *m_couplers) {
    coupler->startLoadTimer(AT_);
    m_postCouple[coupler->couplerId()](a_step());
    coupler->stopLoadTimer(AT_);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/// Recipe for any number of solvers/couplers, where the couplers (if any) are executed once after each solver's
/// substep.
///
class ExecutionRecipeIntraStepCoupling : public ExecutionRecipe {
 public:
  ExecutionRecipeIntraStepCoupling(std::vector<std::unique_ptr<Solver>>* const solvers,
                                   std::vector<std::unique_ptr<Coupling>>* const couplers)
    : ExecutionRecipe(solvers, couplers) {
    // This function initializes function pointers to the preTimeStep, solutionStep
    // and postTimeStep functions of each solver
    initFunctionPointers();

    // Read the callOrder from the propertiesFile
    readCallOrder();
  }

  // Set solvers and couplers active or idle according to the callOrder read in
  // readCallOrder(). Also sets adaptation active or inactive.
  // The function is called in the advance time step iteration of the time step loop.
  // If the maximum number of iteration steps is reached (all solvers have been
  // executed) advanceTimeStep is set to true.
  MBool updateCallOrder() override {
    MBool advanceTimeStep = false;

    nextStep();

    if(a_step() == maxNoSteps()) {
      a_step() = 0;
      advanceTimeStep = true;
    }


    for(MInt solver = 0; solver < noSolvers(); solver++) {
      setSolverStatus(solver, solverOrder(solver));
    }

    for(MInt cpl = 0; cpl < noCouplers(); cpl++) {
      setCouplerStatus(cpl, couplerOrder(cpl));
    }

    setAdaptation();

    return advanceTimeStep;
  }

  // The time step function.
  void timeStep() override {
    TRACE();

    MBool completed = false;
    // Store solver completed status to allow solvers to have different number of RK steps, i.e.,
    // skip solutionStep of a solver that is already finished
    std::vector<MBool> solverCompleted(noSolvers(), false);

    while(!completed) {
      // Exit once all solvers are completed.
      completed = true;

      // Intermediate loop over all solvers
      for(auto&& solver : *a_solvers()) {
        const MInt solverId = m_swapSolverIds.empty() ? solver->solverId() : swapedSolverId(solver->solverId());

        // Solver is active on this domain (has cells) and has not already completed its current
        // solution step, i.e., all RK stages
        if(solverIsActive(solverId) && !solverCompleted[solverId]) {
          // Call solver specific solutionStep() - solvers might be idle (in this step)
          startLoadTimer(solverId);
          solverCompleted[solverId] = solutionStep(solverId);
          stopLoadTimer(solverId);
        }

        // Inner-most loop over all couplers
        // Each couplers subCouple() is executed after each solvers solutionStep()
        // Couplers can, however, be disabled via the callOrder
        for(auto&& coupler : *a_couplers()) {
          coupler->startLoadTimer(AT_);
          subCouple(coupler->couplerId(), a_step(), solverId, solverCompleted);
          coupler->stopLoadTimer(AT_);
        }
        if(solverIsActive(solverId) && !solverCompleted[solverId]) {
          completed &= solverCompleted[solverId];
        }
      }
    }
  };
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// new recepie class to handle multiple solution loops

class ExecutionRecipeSolutionIteration : public ExecutionRecipe {
 public:
  ExecutionRecipeSolutionIteration(std::vector<std::unique_ptr<Solver>>* const solvers,
                                   std::vector<std::unique_ptr<Coupling>>* const couplers)
    : ExecutionRecipe(solvers, couplers) {
    // This function initializes function pointers to the preTimeStep, solutionStep
    // and postTimeStep functions of each solver
    initFunctionPointers();

    // Read the callOrder from the propertiesFile
    readCallOrder();
  }

  // Set solvers and couplers active or idle according to the callOrder read in
  // readCallOrder(). Also sets adaptation active or inactive.
  // The function is called in the advance time step iteration of the time step loop.
  // If the maximum number of iteration steps is reached (all solvers have been
  // executed) advanceTimeStep is set to true.
  MBool updateCallOrder() override {
    MBool advanceTimeStep = false;

    nextStep();

    if(a_step() == maxNoSteps()) {
      a_step() = 0;
      advanceTimeStep = true;
    }


    for(MInt solver = 0; solver < noSolvers(); solver++) {
      setSolverStatus(solver, solverOrder(solver));
    }

    for(MInt cpl = 0; cpl < noCouplers(); cpl++) {
      setCouplerStatus(cpl, couplerOrder(cpl));
    }

    setAdaptation();

    return advanceTimeStep;
  }

  // NOTE: this timeStep allows for a iteration of multiple solution Steps
  //      (i.e. multiple full RungeKutta-Steps during a single time-step!
  //      however a sub-couple is not implemented yet!
  void timeStep() override {
    TRACE();

    MInt iteration = 0;
    while(iteration < m_maxSolutionIteration) {
      // Intermediate loop over all solvers
      for(auto&& solver : *a_solvers()) {
        const MInt solverId = solver->solverId();

        if(iteration > 0) {
          solver->startLoadTimer(AT_);
          preSolutionStep(solverId, -1);
          solver->stopLoadTimer(AT_);
        }

        // Outer loop over substep iterations
        MBool timeStepCompleted = false;
        // Store solver completed status to allow solvers to have different number of RK steps, i.e.,
        // skip solutionStep of a solver that is already finished
        std::vector<MBool> solverTimeStepCompleted(noSolvers(), false);

        while(!timeStepCompleted) {
          timeStepCompleted = true;

          // Solver is active on this domain (has cells) and has not already completed its current
          // solution step, i.e., all RK stages
          if(solver->isActive() && !solverTimeStepCompleted[solverId]) {
            // Call solver specific solutionStep() - solvers might be idle (in this step)
            solver->startLoadTimer(AT_);
            solverTimeStepCompleted[solverId] = solutionStep(solver->solverId());
            timeStepCompleted &= solverTimeStepCompleted[solverId];
            solver->stopLoadTimer(AT_);
          }
        }

        if(m_maxSolutionIteration) {
          solver->startLoadTimer(AT_);
          MBool iterationConverged = postSolutionStep(solverId);
          if(iterationConverged) iteration = m_maxSolutionIteration;
          solver->stopLoadTimer(AT_);
        }
        iteration++;
      }
    }
  };
};


#endif // ifndef EXECUTIONRECIPE_H_
