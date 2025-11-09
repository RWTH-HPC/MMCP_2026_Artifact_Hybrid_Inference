// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredsolverwindowinfo.h"
#include "COMM/mpioverride.h"
#include "IO/parallelio_hdf5.h"
#include "UTIL/kdtree.h"
#include "UTIL/pointbox.h"
#include "fvstructuredwindowmapping.h"
#include "globals.h"
#include "vector"

using namespace std;


template <MInt nDim>
FvStructuredSolverWindowInfo<nDim>::FvStructuredSolverWindowInfo(StructuredGrid<nDim>* grid_,
                                                                 MPI_Comm structuredCommunicator_,
                                                                 const MInt noDomains_,
                                                                 const MInt domainId_,
                                                                 const MInt solverId_)
  : m_grid(grid_),
    m_noBlocks(grid_->m_noBlocks),
    m_blockId(grid_->m_blockId),
    m_noDomains(noDomains_),
    m_domainId(domainId_),
    m_solverId(solverId_),
    m_noGhostLayers(grid_->m_noGhostLayers),
    m_StructuredComm(structuredCommunicator_) {}

template <MInt nDim>
FvStructuredSolverWindowInfo<nDim>::~FvStructuredSolverWindowInfo() {
  // do nothing
}

template <>
void FvStructuredSolverWindowInfo<2>::multiBlockAssembling() {
  constexpr MInt nDim = 2;
  MInt countConnection = 0; // numConnections;
  unique_ptr<StructuredWindowMap<nDim>> newWindow;
  MBool found, notexisted, labell;
  connectionNode temp1(nDim);
  MInt numWindows[nDim] = {0, 0};

  //  set<connectionNode>::iterator it;

  while(connectionset.size() != 0) {
    auto element = connectionset.begin();
    MInt order[2] = {-1, -1};
    MInt step1[2] = {1, 1};
    MInt step2[2] = {1, 1};
    MInt pos1[3], pos2[3], b1, b2;

    pos1[0] = element->pos1[0];
    pos1[1] = element->pos1[1];
    pos2[0] = element->pos2[0];
    pos2[1] = element->pos2[1];
    b1 = element->blockId1;
    b2 = element->blockId2;

    newWindow = make_unique<StructuredWindowMap<nDim>>();
    mapCreate(b1, pos1, pos1, step1, b2, pos2, pos2, step2, order, element->BC, newWindow);
    newWindow->Nstar = -1;

    for(MInt i = 0; i < 2; ++i) {
      found = false;
      pos1[0] = newWindow->start1[0];
      pos1[1] = newWindow->start1[1];

      if(newWindow->order[i] == -1) {
        // D1+

        pos1[i] = newWindow->start1[i] + 1;
        for(MInt j = 0; j < 2; ++j) {
          // D2
          pos2[0] = newWindow->start2[0];
          pos2[1] = newWindow->start2[1];

          notexisted = true;
          for(MInt k = 0; k < 2; ++k) {
            if(newWindow->order[k] == j) notexisted = false;
          }

          if(notexisted) {
            pos2[j] = newWindow->start2[j] + 1;

            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.Nstar = -1;
            found = findConnection(temp1);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            // D2-
            pos2[j] = newWindow->start2[j] - 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.Nstar = -1;
            found = findConnection(temp1);
            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }

        if(found) {
          continue;
        }

        // D1-
        pos1[i] = newWindow->start1[i] - 1;
        for(MInt j = 0; j < 2; ++j) {
          // D2+
          pos2[0] = newWindow->start2[0];
          pos2[1] = newWindow->start2[1];
          notexisted = true;
          for(MInt k = 0; k < 2; ++k) {
            if(newWindow->order[k] == j) {
              notexisted = false;
            }
          }

          if(notexisted) {
            pos2[j] = newWindow->start2[j] + 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.Nstar = -1;
            found = findConnection(temp1);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            // D2-
            pos2[j] = newWindow->start2[j] - 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.Nstar = -1;
            found = findConnection(temp1);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }
      }
    }

    MInt ordercount = 0;
    MBool facewindow, directionInWindow[nDim];
    for(MInt i = 0; i < nDim; ++i) {
      directionInWindow[i] = false;
      if(newWindow->order[i] != -1) {
        ordercount++;
        directionInWindow[i] = true;
      }
    }

    if(ordercount > 2) {
      cout << "Invalid volume mapping found! Are your blocks overlapping or is the grid epsilon too large?" << endl;
    }

    if(ordercount == 2) {
      facewindow = true;
    } else {
      facewindow = false;
    }

    for(MInt i = 0; i < 2; ++i) {
      MInt j = 0;

      if(newWindow->order[i] == -1) {
        for(j = 0; j < 2; ++j) {
          labell = true;
          for(MInt k = 0; k < 2; ++k) {
            if(newWindow->order[k] == j) {
              labell = false;
            }
          }

          if(labell == true) {
            newWindow->order[i] = j;
            break;
          }
        }

        if(facewindow) {
          if(newWindow->start1[i] == 0) {
            newWindow->dc1 = i + 1;
          } else {
            newWindow->dc1 = -i - 1;
          }

          if(newWindow->start2[j] == 0) {
            newWindow->dc2 = j + 1;
          } else {
            newWindow->dc2 = -j - 1;
          }
        } else {
          newWindow->dc1 = 999;
          newWindow->dc2 = 999;
        }
      }
    }

    MInt start1[2], end1[2], start2[2], end2[2];
    MBool goGo = true;
    MInt countDim, countDim2;
    MInt ii, jj;
    while(goGo) {
      goGo = false;
      for(countDim = 0; countDim < nDim; ++countDim) {
        if(directionInWindow[countDim]) {
          countDim2 = newWindow->order[countDim];

          for(MInt i = 0; i < 2; ++i) {
            start1[i] = newWindow->start1[i];
            end1[i] = newWindow->end1[i];
            start2[i] = newWindow->start2[i];
            end2[i] = newWindow->end2[i];
          }
          end1[countDim] = end1[countDim] + newWindow->step1[countDim];
          end2[countDim2] = end2[countDim2] + newWindow->step2[countDim2];
          start1[countDim] = end1[countDim];
          start2[countDim2] = end2[countDim2];


          pos2[newWindow->order[1]] = start2[newWindow->order[1]];
          jj = start1[1];
          do {
            pos2[newWindow->order[0]] = start2[newWindow->order[0]];
            ii = start1[0];
            do {
              pos1[0] = ii;
              pos1[1] = jj;
              if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
                labell = true;
              } else {
                labell = false;
              }

              connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
              found = findConnection(temp2);

              if(!found && domainId() == 0 && newWindow->BC == 4401 && newWindow->Id1 == 2 && newWindow->Id2 != 2) {
                // temp2.print();
                connectionNode temp3(newWindow->BC, newWindow->Id2, pos2, newWindow->Id1, pos1, labell, nDim);
                found = findConnection(temp3);

                if(found) {
                  cout << "wooooo!! somthing is  wrong and i do not know where and why!!!!!!" << endl;
                }
              }
              if(!found) {
                break;
              }
              pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
              ii = ii + newWindow->step1[0];

            } while(ii >= start1[0] && ii <= end1[0]);
            pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
            jj = jj + newWindow->step1[1];

          } while(jj >= start1[1] && jj <= end1[1]);


          if(found) {
            // all connections have been found
            for(MInt m = 0; m < 2; ++m) {
              newWindow->end1[m] = end1[m];
              newWindow->end2[m] = end2[m];
            }
            goGo = true;
          }

          for(MInt i = 0; i < 2; ++i) {
            start1[i] = newWindow->start1[i];
            end1[i] = newWindow->end1[i];
            start2[i] = newWindow->start2[i];
            end2[i] = newWindow->end2[i];
          }
          start1[countDim] = start1[countDim] - newWindow->step1[countDim];
          start2[countDim2] = start2[countDim2] - newWindow->step2[countDim2];
          end1[countDim] = start1[countDim];
          end2[countDim2] = start2[countDim2];

          pos2[newWindow->order[1]] = start2[newWindow->order[1]];
          jj = start1[1];
          do {
            pos2[newWindow->order[0]] = start2[newWindow->order[0]];
            ii = start1[0];
            do {
              pos1[0] = ii;
              pos1[1] = jj;
              if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
                labell = true;
              } else {
                labell = false;
              }

              connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
              found = findConnection(temp2);

              if(!found && domainId() == 0 && newWindow->BC == 4401 && newWindow->Id1 == 2 && newWindow->Id2 != 2) {
                connectionNode temp3(newWindow->BC, newWindow->Id2, pos2, newWindow->Id1, pos1, labell, nDim);
                found = findConnection(temp3);
                if(found) {
                  cout << "wooooo!! somthing is  wrong and i do not know where and why!!!!!!" << endl;
                }
              }
              if(!found) {
                break;
              }
              pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
              ii = ii + newWindow->step1[0];

            } while(ii >= start1[0] && ii <= end1[0]);
            pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
            jj = jj + newWindow->step1[1];

          } while(jj >= start1[1] && jj <= end1[1]);


          if(found) {
            // all connections have been found
            for(MInt m = 0; m < 2; ++m) {
              newWindow->start1[m] = start1[m];
              newWindow->start2[m] = start2[m];
            }
            goGo = true;
          }
        }
      }
    }

    if(newWindow->BC >= 6000 && newWindow->BC < 6010
       && newWindow->Id2
              < newWindow->Id1) { //(newWindow->BC == 6000 && newWindow->Id2 < newWindow->Id1) { //6000er-change
      mapInvert1(newWindow);
    }
    mapNormalize3(newWindow);

    // delete treated connections
    pos2[newWindow->order[1]] = newWindow->start2[newWindow->order[1]];
    for(MInt j = newWindow->start1[1]; j <= newWindow->end1[1]; j = j + newWindow->step1[1]) {
      pos2[newWindow->order[0]] = newWindow->start2[newWindow->order[0]];
      for(MInt i = newWindow->start1[0]; i <= newWindow->end1[0]; i = i + newWindow->step1[0]) {
        pos1[0] = i;
        pos1[1] = j;

        if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
          labell = true;
        } else {
          labell = false;
        }

        connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
        found = findConnection(temp2);

        if(!found) {
          cout << "error!! can not delete the element!!!" << endl;
          cout << "BC: " << newWindow->BC << " id1:" << newWindow->Id1 << " id1:" << newWindow->Id2 << endl;
          cout << "i: " << i << " iend:" << newWindow->end1[0] << " j: " << j << " jend:" << newWindow->end1[1] << endl;
        }

        if(found) {
          removeConnection(temp2);
          countConnection = countConnection + 1;
        }

        pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
      }
      pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
    }

    MInt facecount;
    if(!mapCheck(newWindow)) {
      cout << "invalid window!" << endl;
    }

    facecount = 0;
    for(MInt i = 0; i < nDim; ++i) {
      if(newWindow->end1[i] - newWindow->start1[i] != 0) {
        facecount++;
      }
    }

    switch(facecount) {
      case 1: {
        window1d.push_back(std::move(newWindow));
        numWindows[1]++;
      } break;

      case 0: {
        window0d.push_back(std::move(newWindow));
        numWindows[0]++;
      } break;

      default:
        cout << "ERROR!!! face dim is wrong!!!" << endl;
    }
  }

  m_log << "Connection Statistics: " << numWindows[1] << " 1D-, and " << numWindows[0] << " 0D-connections found"
        << endl;
}


template <>
void FvStructuredSolverWindowInfo<2>::singularityAssembling() {
  TRACE();
  const MInt nDim = 2;
  // TODO_SS labels:FV,toenhance In 2d singularities are points. The following algorithm is a unnecessary overhead

  // cout<<"number of connections: "<<numConnections<<endl;
  MInt countConnection = 0, numConnections = 0, Nstar;
  unique_ptr<StructuredWindowMap<nDim>> newWindow;
  MBool found, notexisted, labell;
  connectionNode temp1(nDim);
  MInt numWindows[nDim]{};
  numConnections = singularconnectionset.size();

  //  typename set<connectionNode>::iterator it;

  while(singularconnectionset.size() != 0) {
    auto element = singularconnectionset.begin();
    MInt order[2] = {-1, -1};
    MInt step1[2] = {1, 1};
    MInt step2[2] = {1, 1};
    MInt pos1[3], pos2[3], b1, b2;

    pos1[0] = element->pos1[0];
    pos1[1] = element->pos1[1];
    pos2[0] = element->pos2[0];
    pos2[1] = element->pos2[1];
    b1 = element->blockId1;
    b2 = element->blockId2;
    Nstar = element->Nstar;

    newWindow = make_unique<StructuredWindowMap<nDim>>();
    mapCreate(b1, pos1, pos1, step1, b2, pos2, pos2, step2, order, element->BC, newWindow);
    newWindow->Nstar = Nstar;

    for(MInt i = 0; i < nDim; ++i) {
      found = false;
      pos1[0] = newWindow->start1[0];
      pos1[1] = newWindow->start1[1];

      if(newWindow->order[i] == -1) {
        // D1+

        pos1[i] = newWindow->start1[i] + 1;
        for(MInt j = 0; j < nDim; ++j) {
          // D2
          pos2[0] = newWindow->start2[0];
          pos2[1] = newWindow->start2[1];

          notexisted = true;
          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) {
              notexisted = false;
            }
          }

          if(notexisted) {
            pos2[j] = newWindow->start2[j] + 1;

            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.Nstar = Nstar;
            found = findConnection(temp1, Nstar);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            // D2-
            pos2[j] = newWindow->start2[j] - 1;
            // temp1(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2);
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.Nstar = Nstar;
            found = findConnection(temp1, Nstar);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }

        // D1-
        pos1[i] = newWindow->start1[i] - 1;
        for(MInt j = 0; j < nDim; ++j) {
          // D2+
          pos2[0] = newWindow->start2[0];
          pos2[1] = newWindow->start2[1];
          notexisted = true;

          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) notexisted = false;
          }

          if(notexisted) {
            pos2[j] = newWindow->start2[j] + 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.Nstar = Nstar;
            found = findConnection(temp1, Nstar);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }
            // D2-
            pos2[j] = newWindow->start2[j] - 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.Nstar = Nstar;
            found = findConnection(temp1, Nstar);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }
      }
    }

    MInt ordercount = 0;
    MBool facewindow, directionInWindow[nDim];

    for(MInt i = 0; i < nDim; ++i) {
      directionInWindow[i] = false;
      if(newWindow->order[i] != -1) {
        ordercount++;
        directionInWindow[i] = true;
      }
    }

    ASSERT(ordercount == 0, "In 2D only ordercount==0 makes sense!");

    if(ordercount > 2) {
      cout << "Invalid volume mapping found! "
           << "Are your blockIds overlapping or is the grid epsilon too large?" << endl;
    }

    if(ordercount == 2) {
      facewindow = true;
    } else {
      facewindow = false;
    }

    for(MInt i = 0; i < nDim; ++i) {
      MInt j = 0;
      if(newWindow->order[i] == -1) {
        for(j = 0; j < nDim; ++j) {
          labell = true;
          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) {
              labell = false;
            }
          }
          if(labell == true) {
            newWindow->order[i] = j;
            break;
          }
        }

        if(facewindow) {
          if(newWindow->start1[i] == 0) {
            newWindow->dc1 = i + 1;
          } else {
            newWindow->dc1 = -i - 1;
          }

          if(newWindow->start2[j] == 0) {
            newWindow->dc2 = j + 1;
          } else {
            newWindow->dc2 = -j - 1;
          }
        } else {
          newWindow->dc1 = 999;
          newWindow->dc2 = 999;
        }
      }
    }

    MInt start1[2], end1[2], start2[2], end2[2];
    MBool goGo = true;
    MInt countDim, countDim2;
    MInt ii, jj;

    while(goGo) {
      goGo = false;
      for(countDim = 0; countDim < nDim; ++countDim) {
        if(directionInWindow[countDim]) {
          countDim2 = newWindow->order[countDim];
          for(MInt i = 0; i < nDim; ++i) {
            start1[i] = newWindow->start1[i];
            end1[i] = newWindow->end1[i];
            start2[i] = newWindow->start2[i];
            end2[i] = newWindow->end2[i];
          }
          end1[countDim] = end1[countDim] + newWindow->step1[countDim];
          end2[countDim2] = end2[countDim2] + newWindow->step2[countDim2];
          start1[countDim] = end1[countDim];
          start2[countDim2] = end2[countDim2];

          pos2[newWindow->order[1]] = start2[newWindow->order[1]];
          jj = start1[1];
          do {
            pos2[newWindow->order[0]] = start2[newWindow->order[0]];
            ii = start1[0];
            do {
              pos1[0] = ii;
              pos1[1] = jj;
              if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
                labell = true;
              } else {
                labell = false;
              }

              connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
              temp2.Nstar = Nstar;
              found = findConnection(temp2, Nstar);

              if(!found) {
                break;
              }
              pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
              ii = ii + newWindow->step1[0];

            } while(ii >= start1[0] && ii <= end1[0]);
            pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
            jj = jj + newWindow->step1[1];

          } while(jj >= start1[1] && jj <= end1[1]);

          if(found) {
            // all connections have been found
            for(MInt m = 0; m < nDim; ++m) {
              newWindow->end1[m] = end1[m];
              newWindow->end2[m] = end2[m];
            }
            goGo = true;
          }

          for(MInt i = 0; i < nDim; ++i) {
            start1[i] = newWindow->start1[i];
            end1[i] = newWindow->end1[i];
            start2[i] = newWindow->start2[i];
            end2[i] = newWindow->end2[i];
          }

          start1[countDim] = start1[countDim] - newWindow->step1[countDim];
          start2[countDim2] = start2[countDim2] - newWindow->step2[countDim2];
          end1[countDim] = start1[countDim];
          end2[countDim2] = start2[countDim2];

          pos2[newWindow->order[1]] = start2[newWindow->order[1]];
          jj = start1[1];
          do {
            pos2[newWindow->order[0]] = start2[newWindow->order[0]];
            ii = start1[0];
            do {
              pos1[0] = ii;
              pos1[1] = jj;

              if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
                labell = true;
              } else {
                labell = false;
              }

              connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
              temp2.Nstar = Nstar;
              found = findConnection(temp2, Nstar);

              if(!found) {
                break;
              }

              pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
              ii = ii + newWindow->step1[0];

            } while(ii >= start1[0] && ii <= end1[0]);
            pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
            jj = jj + newWindow->step1[1];

          } while(jj >= start1[1] && jj <= end1[1]);

          if(found) {
            // all connections have been found
            for(MInt m = 0; m < nDim; ++m) {
              newWindow->start1[m] = start1[m];
              newWindow->start2[m] = start2[m];
            }

            goGo = true;
          }
        }
      }
    }

    if(newWindow->BC >= 6000 && newWindow->BC < 6010
       && newWindow->Id2
              < newWindow->Id1) { //(newWindow->BC == 6000 && newWindow->Id2 < newWindow->Id1) { //6000er-change
      mapInvert1(newWindow);
    }
    mapNormalize3(newWindow);

    // delete treated connections
    pos2[newWindow->order[1]] = newWindow->start2[newWindow->order[1]];
    for(MInt j = newWindow->start1[1]; j <= newWindow->end1[1]; j = j + newWindow->step1[1]) {
      pos2[newWindow->order[0]] = newWindow->start2[newWindow->order[0]];
      for(MInt i = newWindow->start1[0]; i <= newWindow->end1[0]; i = i + newWindow->step1[0]) {
        pos1[0] = i;
        pos1[1] = j;

        if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
          labell = true;
        } else {
          labell = false;
        }

        connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
        temp2.Nstar = Nstar;
        found = findConnection(temp2, Nstar);

        if(!found) {
          cout << "singular  error!! can not delete the element!!!" << endl;
        }

        removeConnection(temp2, Nstar);
        if(found) {
          countConnection = countConnection + 1;
        }
        pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
      }
      pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
    }

    // set Id for singularmap
    newWindow->SingularId = numWindows[0];
    singularwindow.push_back(std::move(newWindow));
    numWindows[0]++;

    if(domainId() == 0) {
      cout << numWindows[0] << " singular connections found (" << countConnection * 100 / numConnections << "% done)"
           << endl;
    }
  }
}


template <>
void FvStructuredSolverWindowInfo<3>::multiBlockAssembling() {
  const MInt nDim = 3;
  MInt countConnection = 0; // numConnections;
  unique_ptr<StructuredWindowMap<nDim>> newWindow;
  MBool found, notexisted, labell;
  connectionNode temp1(nDim);
  MInt numWindows[3] = {0, 0, 0};

  //  set<connectionNode>::iterator it;
  // numConnections=connectionset.size();

  while(connectionset.size() != 0) {
    auto element = connectionset.begin();
    MInt order[3] = {-1, -1, -1};
    MInt step1[3] = {1, 1, 1};
    MInt step2[3] = {1, 1, 1};
    MInt pos1[3], pos2[3], b1, b2;

    pos1[0] = element->pos1[0];
    pos1[1] = element->pos1[1];
    pos1[2] = element->pos1[2];
    pos2[0] = element->pos2[0];
    pos2[1] = element->pos2[1];
    pos2[2] = element->pos2[2];
    b1 = element->blockId1;
    b2 = element->blockId2;

    newWindow = make_unique<StructuredWindowMap<nDim>>();
    mapCreate(b1, pos1, pos1, step1, b2, pos2, pos2, step2, order, element->BC, newWindow);
    newWindow->Nstar = -1;

    for(MInt i = 0; i < nDim; ++i) {
      found = false;
      pos1[0] = newWindow->start1[0];
      pos1[1] = newWindow->start1[1];
      pos1[2] = newWindow->start1[2];

      if(newWindow->order[i] == -1) {
        // D1+

        pos1[i] = newWindow->start1[i] + 1;
        for(MInt j = 0; j < nDim; ++j) {
          // D2
          pos2[0] = newWindow->start2[0];
          pos2[1] = newWindow->start2[1];
          pos2[2] = newWindow->start2[2];

          notexisted = true;
          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) notexisted = false;
          }

          if(notexisted) {
            pos2[j] = newWindow->start2[j] + 1;

            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos1[2] = pos1[2];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.pos2[2] = pos2[2];
            temp1.Nstar = -1;
            found = findConnection(temp1);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            // D2-
            pos2[j] = newWindow->start2[j] - 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos1[2] = pos1[2];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.pos2[2] = pos2[2];
            temp1.Nstar = -1;
            found = findConnection(temp1);
            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }

        if(found) {
          continue;
        }

        // D1-
        pos1[i] = newWindow->start1[i] - 1;
        for(MInt j = 0; j < nDim; ++j) {
          // D2+
          pos2[0] = newWindow->start2[0];
          pos2[1] = newWindow->start2[1];
          pos2[2] = newWindow->start2[2];
          notexisted = true;
          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) {
              notexisted = false;
            }
          }

          if(notexisted) {
            pos2[j] = newWindow->start2[j] + 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos1[2] = pos1[2];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.pos2[2] = pos2[2];
            temp1.Nstar = -1;
            found = findConnection(temp1);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            // D2-
            pos2[j] = newWindow->start2[j] - 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos1[2] = pos1[2];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.pos2[2] = pos2[2];
            temp1.Nstar = -1;
            found = findConnection(temp1);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }
      }
    }

    MInt ordercount = 0;
    MBool facewindow, directionInWindow[3];
    for(MInt i = 0; i < nDim; ++i) {
      directionInWindow[i] = false;
      if(newWindow->order[i] != -1) {
        ordercount++;
        directionInWindow[i] = true;
      }
    }

    if(ordercount > 2) {
      cout << "Invalid volume mapping found! Are your blockIds overlapping or is the grid epsilon too large?" << endl;
    }

    if(ordercount == 2) {
      facewindow = true;
    } else {
      facewindow = false;
    }

    for(MInt i = 0; i < nDim; ++i) {
      MInt j = 0;

      if(newWindow->order[i] == -1) {
        for(j = 0; j < nDim; ++j) {
          labell = true;
          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) {
              labell = false;
            }
          }

          if(labell == true) {
            newWindow->order[i] = j;
            break;
          }
        }

        if(facewindow) {
          if(newWindow->start1[i] == 0) {
            newWindow->dc1 = i + 1;
          } else {
            newWindow->dc1 = -i - 1;
          }

          if(newWindow->start2[j] == 0) {
            newWindow->dc2 = j + 1;
          } else {
            newWindow->dc2 = -j - 1;
          }
        } else {
          newWindow->dc1 = 999;
          newWindow->dc2 = 999;
        }
      }
    }

    MInt start1[3], end1[3], start2[3], end2[3];
    MBool goGo = true;
    MInt countDim, countDim2;
    MInt ii, jj, kk;
    while(goGo) {
      goGo = false;
      for(countDim = 0; countDim < nDim; ++countDim) {
        if(directionInWindow[countDim]) {
          countDim2 = newWindow->order[countDim];

          for(MInt i = 0; i < nDim; ++i) {
            start1[i] = newWindow->start1[i];
            end1[i] = newWindow->end1[i];
            start2[i] = newWindow->start2[i];
            end2[i] = newWindow->end2[i];
          }
          end1[countDim] = end1[countDim] + newWindow->step1[countDim];
          end2[countDim2] = end2[countDim2] + newWindow->step2[countDim2];
          start1[countDim] = end1[countDim];
          start2[countDim2] = end2[countDim2];

          pos2[newWindow->order[2]] = start2[newWindow->order[2]];
          kk = start1[2];
          do {
            pos2[newWindow->order[1]] = start2[newWindow->order[1]];
            jj = start1[1];
            do {
              pos2[newWindow->order[0]] = start2[newWindow->order[0]];
              ii = start1[0];
              do {
                pos1[0] = ii;
                pos1[1] = jj;
                pos1[2] = kk;
                if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
                  labell = true;
                } else {
                  labell = false;
                }

                connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
                found = findConnection(temp2);

                if(!found && domainId() == 0 && newWindow->BC == 4401 && newWindow->Id1 == 2 && newWindow->Id2 != 2) {
                  // temp2.print();
                  connectionNode temp3(newWindow->BC, newWindow->Id2, pos2, newWindow->Id1, pos1, labell, nDim);
                  found = findConnection(temp3);

                  if(found) {
                    cout << "wooooo!! somthing is  wrong and i do not know where and why!!!!!!" << endl;
                  }
                }
                if(!found) {
                  break;
                }
                pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
                ii = ii + newWindow->step1[0];

              } while(ii >= start1[0] && ii <= end1[0]);
              pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
              jj = jj + newWindow->step1[1];

            } while(jj >= start1[1] && jj <= end1[1]);
            pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
            kk = kk + newWindow->step1[2];

          } while(kk >= start1[2] && kk <= end1[2]);

          if(found) {
            // all connections have been found
            for(MInt m = 0; m < nDim; ++m) {
              newWindow->end1[m] = end1[m];
              newWindow->end2[m] = end2[m];
            }
            goGo = true;
          }

          for(MInt i = 0; i < nDim; ++i) {
            start1[i] = newWindow->start1[i];
            end1[i] = newWindow->end1[i];
            start2[i] = newWindow->start2[i];
            end2[i] = newWindow->end2[i];
          }
          start1[countDim] = start1[countDim] - newWindow->step1[countDim];
          start2[countDim2] = start2[countDim2] - newWindow->step2[countDim2];
          end1[countDim] = start1[countDim];
          end2[countDim2] = start2[countDim2];

          pos2[newWindow->order[2]] = start2[newWindow->order[2]];
          kk = start1[2];
          do {
            pos2[newWindow->order[1]] = start2[newWindow->order[1]];
            jj = start1[1];
            do {
              pos2[newWindow->order[0]] = start2[newWindow->order[0]];
              ii = start1[0];
              do {
                pos1[0] = ii;
                pos1[1] = jj;
                pos1[2] = kk;
                if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
                  labell = true;
                } else {
                  labell = false;
                }

                connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
                found = findConnection(temp2);

                if(!found && domainId() == 0 && newWindow->BC == 4401 && newWindow->Id1 == 2 && newWindow->Id2 != 2) {
                  connectionNode temp3(newWindow->BC, newWindow->Id2, pos2, newWindow->Id1, pos1, labell, nDim);
                  found = findConnection(temp3);
                  if(found) {
                    cout << "wooooo!! somthing is  wrong and i do not know where and why!!!!!!" << endl;
                  }
                }
                if(!found) {
                  break;
                }
                pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
                ii = ii + newWindow->step1[0];

              } while(ii >= start1[0] && ii <= end1[0]);
              pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
              jj = jj + newWindow->step1[1];

            } while(jj >= start1[1] && jj <= end1[1]);
            pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
            kk = kk + newWindow->step1[2];

          } while(kk >= start1[2] && kk <= end1[2]);


          if(found) {
            // all connections have been found
            for(MInt m = 0; m < nDim; ++m) {
              newWindow->start1[m] = start1[m];
              newWindow->start2[m] = start2[m];
            }
            goGo = true;
          }
        }
      }
    }

    if(newWindow->BC >= 6000 && newWindow->BC < 6010
       && newWindow->Id2
              < newWindow->Id1) { //(newWindow->BC == 6000 && newWindow->Id2 < newWindow->Id1) { //6000er-change
      mapInvert1(newWindow);
    }
    mapNormalize3(newWindow);

    // delete treated connections
    pos2[newWindow->order[2]] = newWindow->start2[newWindow->order[2]];
    for(MInt k = newWindow->start1[2]; k <= newWindow->end1[2]; k = k + newWindow->step1[2]) {
      pos2[newWindow->order[1]] = newWindow->start2[newWindow->order[1]];
      for(MInt j = newWindow->start1[1]; j <= newWindow->end1[1]; j = j + newWindow->step1[1]) {
        pos2[newWindow->order[0]] = newWindow->start2[newWindow->order[0]];
        for(MInt i = newWindow->start1[0]; i <= newWindow->end1[0]; i = i + newWindow->step1[0]) {
          pos1[0] = i;
          pos1[1] = j;
          pos1[2] = k;

          if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
            labell = true;
          } else {
            labell = false;
          }

          connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
          found = findConnection(temp2);

          if(!found) {
            cout << "error!! can not delete the element!!!" << endl;
            cout << "BC: " << newWindow->BC << " id1:" << newWindow->Id1 << " id1:" << newWindow->Id2 << endl;
            cout << "i: " << i << " iend:" << newWindow->end1[0] << " j: " << j << " jend:" << newWindow->end1[1]
                 << " k: " << k << " kend:" << newWindow->end1[2] << endl;
          }

          if(found) {
            removeConnection(temp2);
            countConnection = countConnection + 1;
          }

          pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
        }
        pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
      }
      pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
    }

    MInt facecount;
    if(!mapCheck(newWindow)) {
      cout << "invalid window!" << endl;
    }

    facecount = 0;
    for(MInt i = 0; i < nDim; ++i) {
      if(newWindow->end1[i] - newWindow->start1[i] != 0) {
        facecount++;
      }
    }

    switch(facecount) {
      case 2: {
        window2d.push_back(std::move(newWindow));
        numWindows[2]++;
      } break;

      case 1: {
        window1d.push_back(std::move(newWindow));
        numWindows[1]++;
      } break;

      case 0: {
        window0d.push_back(std::move(newWindow));
        numWindows[0]++;
      } break;

      default:
        cout << "ERROR!!! face dim is wrong!!!" << endl;
    }
  }

  m_log << "Connection Statistics: " << numWindows[2] << " 2D-, " << numWindows[1] << " 1D-, and " << numWindows[0]
        << " 0D-connections found" << endl;
}

template <>
void FvStructuredSolverWindowInfo<3>::singularityAssembling() {
  const MInt nDim = 3;
  // cout<<"number of connections: "<<numConnections<<endl;
  MInt countConnection = 0, numConnections = 0, Nstar;
  unique_ptr<StructuredWindowMap<nDim>> newWindow;
  MBool found, notexisted, labell;
  connectionNode temp1(nDim);
  MInt numWindows[3] = {0, 0, 0};
  numConnections = singularconnectionset.size();

  //  typename set<connectionNode>::iterator it;

  while(singularconnectionset.size() != 0) {
    auto element = singularconnectionset.begin();
    MInt order[3] = {-1, -1, -1};
    MInt step1[3] = {1, 1, 1};
    MInt step2[3] = {1, 1, 1};
    MInt pos1[3], pos2[3], b1, b2;

    pos1[0] = element->pos1[0];
    pos1[1] = element->pos1[1];
    pos1[2] = element->pos1[2];
    pos2[0] = element->pos2[0];
    pos2[1] = element->pos2[1];
    pos2[2] = element->pos2[2];
    b1 = element->blockId1;
    b2 = element->blockId2;
    Nstar = element->Nstar;

    newWindow = make_unique<StructuredWindowMap<nDim>>();
    mapCreate(b1, pos1, pos1, step1, b2, pos2, pos2, step2, order, element->BC, newWindow);
    newWindow->Nstar = Nstar;

    for(MInt i = 0; i < nDim; ++i) {
      found = false;
      pos1[0] = newWindow->start1[0];
      pos1[1] = newWindow->start1[1];
      pos1[2] = newWindow->start1[2];

      if(newWindow->order[i] == -1) {
        // D1+

        pos1[i] = newWindow->start1[i] + 1;
        for(MInt j = 0; j < nDim; ++j) {
          // D2
          pos2[0] = newWindow->start2[0];
          pos2[1] = newWindow->start2[1];
          pos2[2] = newWindow->start2[2];

          notexisted = true;
          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) {
              notexisted = false;
            }
          }

          if(notexisted) {
            pos2[j] = newWindow->start2[j] + 1;

            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos1[2] = pos1[2];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.pos2[2] = pos2[2];
            temp1.Nstar = Nstar;
            found = findConnection(temp1, Nstar);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            // D2-
            pos2[j] = newWindow->start2[j] - 1;
            // temp1(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2);
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos1[2] = pos1[2];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.pos2[2] = pos2[2];
            temp1.Nstar = Nstar;
            found = findConnection(temp1, Nstar);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }

        // D1-
        pos1[i] = newWindow->start1[i] - 1;
        for(MInt j = 0; j < nDim; ++j) {
          // D2+
          pos2[0] = newWindow->start2[0];
          pos2[1] = newWindow->start2[1];
          pos2[2] = newWindow->start2[2];
          notexisted = true;

          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) notexisted = false;
          }

          if(notexisted) {
            pos2[j] = newWindow->start2[j] + 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos1[2] = pos1[2];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.pos2[2] = pos2[2];
            temp1.Nstar = Nstar;
            found = findConnection(temp1, Nstar);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }
            // D2-
            pos2[j] = newWindow->start2[j] - 1;
            temp1.BC = newWindow->BC;
            temp1.blockId1 = newWindow->Id1;
            temp1.blockId2 = newWindow->Id2;
            temp1.pos1[0] = pos1[0];
            temp1.pos1[1] = pos1[1];
            temp1.pos1[2] = pos1[2];
            temp1.pos2[0] = pos2[0];
            temp1.pos2[1] = pos2[1];
            temp1.pos2[2] = pos2[2];
            temp1.Nstar = Nstar;
            found = findConnection(temp1, Nstar);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }
      }
    }

    MInt ordercount = 0;
    MBool facewindow, directionInWindow[3];

    for(MInt i = 0; i < nDim; ++i) {
      directionInWindow[i] = false;
      if(newWindow->order[i] != -1) {
        ordercount++;
        directionInWindow[i] = true;
      }
    }

    if(ordercount > 2) {
      cout << "Invalid volume mapping found! "
           << "Are your blockIds overlapping or is the grid epsilon too large?" << endl;
    }

    if(ordercount == 2) {
      facewindow = true;
    } else {
      facewindow = false;
    }

    for(MInt i = 0; i < nDim; ++i) {
      MInt j = 0;
      if(newWindow->order[i] == -1) {
        for(j = 0; j < nDim; ++j) {
          labell = true;
          for(MInt k = 0; k < nDim; ++k) {
            if(newWindow->order[k] == j) {
              labell = false;
            }
          }
          if(labell == true) {
            newWindow->order[i] = j;
            break;
          }
        }

        if(facewindow) {
          if(newWindow->start1[i] == 0) {
            newWindow->dc1 = i + 1;
          } else {
            newWindow->dc1 = -i - 1;
          }

          if(newWindow->start2[j] == 0) {
            newWindow->dc2 = j + 1;
          } else {
            newWindow->dc2 = -j - 1;
          }
        } else {
          newWindow->dc1 = 999;
          newWindow->dc2 = 999;
        }
      }
    }

    MInt start1[3], end1[3], start2[3], end2[3];
    MBool goGo = true;
    MInt countDim, countDim2;
    MInt ii, jj, kk;

    while(goGo) {
      goGo = false;
      for(countDim = 0; countDim < nDim; ++countDim) {
        if(directionInWindow[countDim]) {
          countDim2 = newWindow->order[countDim];
          for(MInt i = 0; i < nDim; ++i) {
            start1[i] = newWindow->start1[i];
            end1[i] = newWindow->end1[i];
            start2[i] = newWindow->start2[i];
            end2[i] = newWindow->end2[i];
          }
          end1[countDim] = end1[countDim] + newWindow->step1[countDim];
          end2[countDim2] = end2[countDim2] + newWindow->step2[countDim2];
          start1[countDim] = end1[countDim];
          start2[countDim2] = end2[countDim2];

          pos2[newWindow->order[2]] = start2[newWindow->order[2]];
          kk = start1[2];

          do {
            pos2[newWindow->order[1]] = start2[newWindow->order[1]];
            jj = start1[1];
            do {
              pos2[newWindow->order[0]] = start2[newWindow->order[0]];
              ii = start1[0];
              do {
                pos1[0] = ii;
                pos1[1] = jj;
                pos1[2] = kk;
                if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
                  labell = true;
                } else {
                  labell = false;
                }

                connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
                temp2.Nstar = Nstar;
                found = findConnection(temp2, Nstar);

                if(!found) {
                  break;
                }
                pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
                ii = ii + newWindow->step1[0];

              } while(ii >= start1[0] && ii <= end1[0]);
              pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
              jj = jj + newWindow->step1[1];

            } while(jj >= start1[1] && jj <= end1[1]);
            pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
            kk = kk + newWindow->step1[2];

          } while(kk >= start1[2] && kk <= end1[2]);

          if(found) {
            // all connections have been found
            for(MInt m = 0; m < nDim; ++m) {
              newWindow->end1[m] = end1[m];
              newWindow->end2[m] = end2[m];
            }
            goGo = true;
          }

          for(MInt i = 0; i < nDim; ++i) {
            start1[i] = newWindow->start1[i];
            end1[i] = newWindow->end1[i];
            start2[i] = newWindow->start2[i];
            end2[i] = newWindow->end2[i];
          }

          start1[countDim] = start1[countDim] - newWindow->step1[countDim];
          start2[countDim2] = start2[countDim2] - newWindow->step2[countDim2];
          end1[countDim] = start1[countDim];
          end2[countDim2] = start2[countDim2];

          pos2[newWindow->order[2]] = start2[newWindow->order[2]];
          kk = start1[2];

          do {
            pos2[newWindow->order[1]] = start2[newWindow->order[1]];
            jj = start1[1];
            do {
              pos2[newWindow->order[0]] = start2[newWindow->order[0]];
              ii = start1[0];
              do {
                pos1[0] = ii;
                pos1[1] = jj;
                pos1[2] = kk;

                if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
                  labell = true;
                } else {
                  labell = false;
                }

                connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
                temp2.Nstar = Nstar;
                found = findConnection(temp2, Nstar);

                if(!found) {
                  break;
                }

                pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
                ii = ii + newWindow->step1[0];

              } while(ii >= start1[0] && ii <= end1[0]);
              pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
              jj = jj + newWindow->step1[1];

            } while(jj >= start1[1] && jj <= end1[1]);
            pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
            kk = kk + newWindow->step1[2];

          } while(kk >= start1[2] && kk <= end1[2]);

          if(found) {
            // all connections have been found
            for(MInt m = 0; m < nDim; ++m) {
              newWindow->start1[m] = start1[m];
              newWindow->start2[m] = start2[m];
            }

            goGo = true;
          }
        }
      }
    }

    if(newWindow->BC >= 6000 && newWindow->BC < 6010
       && newWindow->Id2
              < newWindow->Id1) { //(newWindow->BC == 6000 && newWindow->Id2 < newWindow->Id1) { //6000er-change
      mapInvert1(newWindow);
    }
    mapNormalize3(newWindow);

    // delete treated connections
    pos2[newWindow->order[2]] = newWindow->start2[newWindow->order[2]];
    for(MInt k = newWindow->start1[2]; k <= newWindow->end1[2]; k = k + newWindow->step1[2]) {
      pos2[newWindow->order[1]] = newWindow->start2[newWindow->order[1]];
      for(MInt j = newWindow->start1[1]; j <= newWindow->end1[1]; j = j + newWindow->step1[1]) {
        pos2[newWindow->order[0]] = newWindow->start2[newWindow->order[0]];
        for(MInt i = newWindow->start1[0]; i <= newWindow->end1[0]; i = i + newWindow->step1[0]) {
          pos1[0] = i;
          pos1[1] = j;
          pos1[2] = k;

          if(newWindow->BC >= 6000 && newWindow->BC < 6010) { //(newWindow->BC==6000) { //6000er-change
            labell = true;
          } else {
            labell = false;
          }

          connectionNode temp2(newWindow->BC, newWindow->Id1, pos1, newWindow->Id2, pos2, labell, nDim);
          temp2.Nstar = Nstar;
          found = findConnection(temp2, Nstar);

          if(!found) {
            cout << "singular  error!! can not delete the element!!!" << endl;
          }

          removeConnection(temp2, Nstar);
          if(found) {
            countConnection = countConnection + 1;
          }
          pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
        }
        pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
      }
      pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
    }

    // set Id for singularmap
    newWindow->SingularId = numWindows[0];
    singularwindow.push_back(std::move(newWindow));
    numWindows[0]++;

    if(domainId() == 0) {
      cout << numWindows[0] << " singular connections found (" << countConnection * 100 / numConnections << "% done)"
           << endl;
    }
  }
}

template <>
void FvStructuredSolverWindowInfo<3>::readWindowCoordinates(MFloat* periodicDisplacements) {
  constexpr MInt nDim = 3;

  MInt hasConnectionInfo = 0;
  if(globalDomainId() == 0) {
    ParallelIoHdf5 pio(m_grid->m_gridInputFileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);
    MBool attributeExists = pio.hasAttribute("hasConnectionInfo", "Connectivity");
    if(attributeExists) {
      m_log << "Grid file has connection info!" << endl;
      pio.getAttribute(&hasConnectionInfo, "hasConnectionInfo", "Connectivity");
    }
    MPI_Bcast(&hasConnectionInfo, 1, MPI_INT, 0, m_StructuredComm, AT_, "hasConnectionInfo");
  } else {
    MPI_Bcast(&hasConnectionInfo, 1, MPI_INT, 0, m_StructuredComm, AT_, "hasConnectionInfo");
  }

  MBool readInConnectionInfo = true;

  if(!hasConnectionInfo) {
    readInConnectionInfo = false;
  } else {
    if(Context::propertyExists("readInConnectionInfo", m_solverId)) {
      readInConnectionInfo =
          Context::getSolverProperty<MBool>("readInConnectionInfo", solverId(), AT_, &readInConnectionInfo);
    }
  }

  if(readInConnectionInfo) {
    m_log << "Connection info available, reading connection information from grid file!" << endl;
    readConnectionWindowInformation3D(periodicDisplacements);
  } else {
    if(domainId() == 0) {
      cout << "Starting grid connection info search..." << endl;
    }
    m_log << " Starting grid connection info search..." << endl;
    MInt offset[3], size[3], count = 0;
    // 3) read in the coordinates of the grid points
    // open file for reading the grid data
    // pointType<nDim>* nodeMap;
    MInt totalGridCells = 0;
    MInt** totalGridBlockDim = nullptr;
    MInt* totalGridBlockCells = nullptr;

    mAlloc(totalGridBlockDim, m_noBlocks, 3, "totalGridBlockDim", -1, AT_);
    mAlloc(totalGridBlockCells, m_noBlocks, "totalGridBlockCells", -1, AT_);

    for(MInt i = 0; i < m_noBlocks; ++i) {
      MInt temp = 1;
      for(MInt dim = 0; dim < 3; ++dim) {
        temp *= m_grid->getBlockNoCells(i, dim) + 1;
        totalGridBlockDim[i][dim] = m_grid->getBlockNoCells(i, dim) + 1; // number of points in the grid File
      }
      totalGridBlockCells[i] = temp;
      if(temp != 1) totalGridCells += temp;
    }

    MInt countNode = 0;
    for(MInt i = 0; i < m_noBlocks; i++) {
      countNode = countNode + 2 * totalGridBlockDim[i][1] * totalGridBlockDim[i][2]
                  + 2 * totalGridBlockDim[i][2] * (totalGridBlockDim[i][0] - 2)
                  + 2 * (totalGridBlockDim[i][0] - 2) * (totalGridBlockDim[i][1] - 2);
    }

    MFloatScratchSpace coordinates(3, countNode, AT_, "coordinates");
    coordinates.fill(-1.01010101);
    std::vector<pointType<nDim>> nodeMap;
    nodeMap.resize(countNode);


    // auxiliary lambda functions to determine BC of interface points
    const MInt converter[] = {2, 1, 0};
    auto getCurrentWindow = [this, converter](const MInt blockId, const MInt* const offset_, const MInt* const size_) {
      std::vector<MInt> currentWindows;
      for(MInt windowId = 0; windowId < noInputWindowInformation; ++windowId) {
        if(inputWindows[windowId]->blockId != blockId) continue;
        if(inputWindows[windowId]->BC < 6000 || inputWindows[windowId]->BC >= 6010) continue;
        // Besides the windows of the same face, windows of adjacent faces are also candidates for the corner points
        MBool takeIt = true;
        for(MInt dim = 0; dim < nDim; ++dim) {
          if(inputWindows[windowId]->startindex[dim] >= offset_[converter[dim]] + size_[converter[dim]]
             || inputWindows[windowId]->endindex[dim] < offset_[converter[dim]]) {
            takeIt = false;
            break;
          }
        }
        if(takeIt) currentWindows.push_back(windowId);
      }
      return currentWindows;
    };
    std::vector<MInt> pointsFoundPerWindow(noInputWindowInformation); // for detailed output
    auto getBCFromWindow = [this, &pointsFoundPerWindow](const MInt idx[nDim],
                                                         const std::vector<MInt>& currentWindows) {
      std::vector<MInt> BCset; // CHANGE_SET
      // Use set in order to have a unique list
      //      std::set<MInt> BCset;
      //      MInt BC[3] = {0, 0, 0}; //we can have at most 3 different BCs shared by one point
      //      MInt cnt = 0;
      for(auto windowId : currentWindows) {
        MBool takeIt = true;
        for(MInt dim = 0; dim < nDim; ++dim) {
          if(inputWindows[windowId]->startindex[dim] > idx[dim] || inputWindows[windowId]->endindex[dim] < idx[dim]) {
            takeIt = false;
            break;
          }
        }
        if(takeIt) {
          ++pointsFoundPerWindow[windowId];
          BCset.push_back(inputWindows[windowId]->BC); // CHANGE_SET
          //          BCset.insert(inputWindows[windowId]->BC);
          //          BC[cnt++] = inputWindows[windowId]->BC;
        }
      }
      //      if (cnt>3) mTerm(1, "");
      //      if (cnt>0) return MAX3(BC[0], BC[1], BC[2]); //higher BC has higher priority; in future think of something
      //      more sophisticated return 6000; //default //-1;
      // we can have at most 3 different BCs shared by one point
      if(BCset.size() > 3) mTerm(1, "");
      if(BCset.size() == 0) BCset.push_back(6000); // BCset.insert(6000); //CHANGE_SET
      MInt BC[3] = {0, 0, 0};
      MInt cnt = 0;
      for(auto it = BCset.begin(); it != BCset.end(); ++it)
        BC[cnt++] = *it;
      return std::tuple<MInt, MInt, MInt>(BC[0], BC[1], BC[2]);
    };

    MInt memsize = 0;
    // domain 0 reads the coordinates for all 6 sides of all blocks
    if(domainId() == 0) {
      cout << "Reading in all coordinates of all block faces" << endl;
    }

    ParallelIoHdf5 pio(m_grid->m_gridInputFileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);

    for(MInt i = 0; i < m_noBlocks; ++i) {
      MString bName = "/block";
      stringstream number;
      number << i << "/";
      bName += number.str();

      ParallelIo::size_type ioOffset[3] = {0, 0, 0};
      ParallelIo::size_type ioSize[3] = {0, 0, 0};
      // in grid file 2:i  1:j   0:k
      //////////////////////////////
      ///////// FACE 1 /////////////
      //////////////////////////////
      // face +x  face 1

      if(m_domainId == 0) {
        ioOffset[2] = 0;
        ioSize[2] = 1; // i
        ioOffset[1] = 0;
        ioSize[1] = m_grid->getBlockNoCells(i, 1) + 1; // j
        ioOffset[0] = 0;
        ioSize[0] = m_grid->getBlockNoCells(i, 0) + 1; // k

        pio.readArray(&coordinates(0, memsize), bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(2, memsize), bName, "z", 3, ioOffset, ioSize);
      } else {
        ioOffset[2] = 0;
        ioSize[2] = 0; // i
        ioOffset[1] = 0;
        ioSize[1] = 0; // j
        ioOffset[0] = 0;
        ioSize[0] = 0; // k

        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "z", 3, ioOffset, ioSize);
      }
      for(int dim = 0; dim < 3; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      std::vector<MInt> currentWindows = getCurrentWindow(i, offset, size);

      // write all computational coordinates
      // of the face to nodeMap object
      for(MInt a = offset[2]; a < offset[2] + size[2]; ++a) {
        for(MInt c = offset[0]; c < offset[0] + size[0]; ++c) {
          for(MInt b = offset[1]; b < offset[1] + size[1]; ++b) {
            nodeMap[count].blockId = i;
            nodeMap[count].pos[0] = a;
            nodeMap[count].pos[1] = b;
            nodeMap[count].pos[2] = c;
            std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1], nodeMap[count].BC[2]) =
                getBCFromWindow(nodeMap[count].pos, currentWindows);
            count++;
          }
        }
      }
      memsize += size[0] * size[1] * size[2];

      //////////////////////////////
      ///////// FACE 2 /////////////
      //////////////////////////////
      // face -x  face 2
      if(m_domainId == 0) {
        ioOffset[2] = m_grid->getBlockNoCells(i, 2);
        ioSize[2] = 1; // i
        ioOffset[1] = 0;
        ioSize[1] = m_grid->getBlockNoCells(i, 1) + 1; // j
        ioOffset[0] = 0;
        ioSize[0] = m_grid->getBlockNoCells(i, 0) + 1; // k
        pio.readArray(&coordinates(0, memsize), bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(2, memsize), bName, "z", 3, ioOffset, ioSize);
      } else {
        ioOffset[2] = 0;
        ioSize[2] = 0; // i
        ioOffset[1] = 0;
        ioSize[1] = 0; // j
        ioOffset[0] = 0;
        ioSize[0] = 0; // k
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "z", 3, ioOffset, ioSize);
      }
      for(int dim = 0; dim < 3; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      currentWindows = getCurrentWindow(i, offset, size);

      for(MInt a = offset[2]; a < offset[2] + size[2]; ++a) {
        for(MInt c = offset[0]; c < offset[0] + size[0]; ++c) {
          for(MInt b = offset[1]; b < offset[1] + size[1]; ++b) {
            nodeMap[count].blockId = i;
            nodeMap[count].pos[0] = a;
            nodeMap[count].pos[1] = b;
            nodeMap[count].pos[2] = c;
            std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1], nodeMap[count].BC[2]) =
                getBCFromWindow(nodeMap[count].pos, currentWindows);
            count++;
          }
        }
      }
      memsize += size[0] * size[1] * size[2];

      //////////////////////////////
      ///////// FACE 3 /////////////
      //////////////////////////////
      // face +y  face 3
      if(m_domainId == 0) {
        ioOffset[2] = 1;
        ioSize[2] = m_grid->getBlockNoCells(i, 2) - 1; // i
        ioOffset[1] = 0;
        ioSize[1] = 1; // j
        ioOffset[0] = 0;
        ioSize[0] = m_grid->getBlockNoCells(i, 0) + 1; // k
        pio.readArray(&coordinates(0, memsize), bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(2, memsize), bName, "z", 3, ioOffset, ioSize);
      } else {
        ioOffset[2] = 0;
        ioSize[2] = 0; // i
        ioOffset[1] = 0;
        ioSize[1] = 0; // j
        ioOffset[0] = 0;
        ioSize[0] = 0; // k
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "z", 3, ioOffset, ioSize);
      }
      for(int dim = 0; dim < 3; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      currentWindows = getCurrentWindow(i, offset, size);

      for(MInt b = offset[1]; b < offset[1] + size[1]; ++b) {
        for(MInt c = offset[0]; c < offset[0] + size[0]; ++c) {
          for(MInt a = offset[2]; a < offset[2] + size[2]; ++a) {
            nodeMap[count].blockId = i;
            nodeMap[count].pos[0] = a;
            nodeMap[count].pos[1] = b;
            nodeMap[count].pos[2] = c;
            std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1], nodeMap[count].BC[2]) =
                getBCFromWindow(nodeMap[count].pos, currentWindows);
            count++;
          }
        }
      }
      memsize += size[0] * size[1] * size[2];

      //////////////////////////////
      ///////// FACE 4 /////////////
      //////////////////////////////
      // face -y  face 4
      if(m_domainId == 0) {
        ioOffset[2] = 1;
        ioSize[2] = m_grid->getBlockNoCells(i, 2) - 1; // i
        ioOffset[1] = m_grid->getBlockNoCells(i, 1);
        ioSize[1] = 1; // j
        ioOffset[0] = 0;
        ioSize[0] = m_grid->getBlockNoCells(i, 0) + 1; // k
        pio.readArray(&coordinates(0, memsize), bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(2, memsize), bName, "z", 3, ioOffset, ioSize);
      } else {
        ioOffset[2] = 0;
        ioSize[2] = 0; // i
        ioOffset[1] = 0;
        ioSize[1] = 0; // j
        ioOffset[0] = 0;
        ioSize[0] = 0; // k
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "z", 3, ioOffset, ioSize);
      }
      for(int dim = 0; dim < 3; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      currentWindows = getCurrentWindow(i, offset, size);

      for(MInt b = offset[1]; b < offset[1] + size[1]; ++b) {
        for(MInt c = offset[0]; c < offset[0] + size[0]; ++c) {
          for(MInt a = offset[2]; a < offset[2] + size[2]; ++a) {
            nodeMap[count].blockId = i;
            nodeMap[count].pos[0] = a;
            nodeMap[count].pos[1] = b;
            nodeMap[count].pos[2] = c;
            std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1], nodeMap[count].BC[2]) =
                getBCFromWindow(nodeMap[count].pos, currentWindows);
            count++;
          }
        }
      }
      memsize += size[0] * size[1] * size[2];

      //////////////////////////////
      ///////// FACE 5 /////////////
      //////////////////////////////
      // face +z  face 5
      if(m_domainId == 0) {
        ioOffset[2] = 1;
        ioSize[2] = m_grid->getBlockNoCells(i, 2) - 1; // i
        ioOffset[1] = 1;
        ioSize[1] = m_grid->getBlockNoCells(i, 1) - 1; // j
        ioOffset[0] = 0;
        ioSize[0] = 1; // k
        pio.readArray(&coordinates(0, memsize), bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(2, memsize), bName, "z", 3, ioOffset, ioSize);
      } else {
        ioOffset[2] = 0;
        ioSize[2] = 0; // i
        ioOffset[1] = 0;
        ioSize[1] = 0; // j
        ioOffset[0] = 0;
        ioSize[0] = 0; // k
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "z", 3, ioOffset, ioSize);
      }
      for(int dim = 0; dim < 3; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      currentWindows = getCurrentWindow(i, offset, size);

      for(MInt c = offset[0]; c < offset[0] + size[0]; ++c) {
        for(MInt b = offset[1]; b < offset[1] + size[1]; ++b) {
          for(MInt a = offset[2]; a < offset[2] + size[2]; ++a) {
            nodeMap[count].blockId = i;
            nodeMap[count].pos[0] = a;
            nodeMap[count].pos[1] = b;
            nodeMap[count].pos[2] = c;
            std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1], nodeMap[count].BC[2]) =
                getBCFromWindow(nodeMap[count].pos, currentWindows);
            count++;
          }
        }
      }
      memsize += size[0] * size[1] * size[2];

      //////////////////////////////
      ///////// FACE 6 /////////////
      //////////////////////////////
      // face -z  face 6
      if(m_domainId == 0) {
        ioOffset[2] = 1;
        ioSize[2] = m_grid->getBlockNoCells(i, 2) - 1; // i
        ioOffset[1] = 1;
        ioSize[1] = m_grid->getBlockNoCells(i, 1) - 1; // j
        ioOffset[0] = m_grid->getBlockNoCells(i, 0);
        ioSize[0] = 1; // k
        pio.readArray(&coordinates(0, memsize), bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&coordinates(2, memsize), bName, "z", 3, ioOffset, ioSize);
      } else {
        ioOffset[2] = 0;
        ioSize[2] = 0; // i
        ioOffset[1] = 0;
        ioSize[1] = 0; // j
        ioOffset[0] = 0;
        ioSize[0] = 0; // k
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 3, ioOffset, ioSize);
        pio.readArray(&empty, bName, "z", 3, ioOffset, ioSize);
      }
      for(int dim = 0; dim < 3; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      currentWindows = getCurrentWindow(i, offset, size);

      for(MInt c = offset[0]; c < offset[0] + size[0]; ++c) {
        for(MInt b = offset[1]; b < offset[1] + size[1]; ++b) {
          for(MInt a = offset[2]; a < offset[2] + size[2]; ++a) {
            nodeMap[count].blockId = i;
            nodeMap[count].pos[0] = a;
            nodeMap[count].pos[1] = b;
            nodeMap[count].pos[2] = c;
            std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1], nodeMap[count].BC[2]) =
                getBCFromWindow(nodeMap[count].pos, currentWindows);
            count++;
          }
        }
      }
      memsize += size[0] * size[1] * size[2];
    }

    // Sanity check
    for(MInt windowId = 0; windowId < noInputWindowInformation; ++windowId) {
      if(inputWindows[windowId]->BC < 6000 || inputWindows[windowId]->BC >= 6010) continue;
      MInt noWindowPoints = 1;
      for(MInt dim = 0; dim < nDim; ++dim) {
        noWindowPoints *= inputWindows[windowId]->endindex[dim] - inputWindows[windowId]->startindex[dim];
      }
      // Every point should be found at least once
      if(noWindowPoints > pointsFoundPerWindow[windowId]) mTerm(1, "noWindowPOints > pointsFoundPerWindow");
    }

    // now broadcast the surface coordinates to every processor
    MPI_Bcast(&count, 1, MPI_INT, 0, m_StructuredComm, AT_, "count");
    MPI_Bcast(&memsize, 1, MPI_INT, 0, m_StructuredComm, AT_, "memsize");
    MPI_Bcast(&coordinates(0, 0), memsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "coordinates(0");
    MPI_Bcast(&coordinates(1, 0), memsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "coordinates(1");
    MPI_Bcast(&coordinates(2, 0), memsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "coordinates(2");

    // also broadcast the nodeMap to every processor
    MPI_Datatype matrix;
    MPI_Type_contiguous(2 + 2 * nDim, MPI_INT, &matrix, AT_);
    MPI_Type_commit(&matrix, AT_);
    MPI_Bcast(nodeMap.data(), memsize, matrix, 0, m_StructuredComm, AT_, "nodeMap");

    /////////////////////////////////
    /////// CONNECTIONS /////////////
    /////////////////////////////////
    if(domainId() == 0) {
      cout << "Building up connections for multiblock connection search" << endl;
    }
    vector<Point<3>> pts;
    for(MInt j = 0; j < count; ++j) {
      Point<3> a(coordinates(0, j), coordinates(1, j), coordinates(2, j));
      pts.push_back(a);
      nodeMap[j].found = false;
    }

    MFloat m_gridEps = 0.0000001;

    KDtree<3> tree(pts);
    MBool add;
    MInt numConnections = 0, numSingularConnections = 0, nfound, tempnum = 0, tempcount;
    MInt results[10];
    for(MInt i = 0; i < count; ++i) {
      if(!nodeMap[i].found) {
        Point<3> a(coordinates(0, i), coordinates(1, i), coordinates(2, i));
        nfound = tree.locatenear(a, m_gridEps, results, 10, false); // 0.00001
        for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
          for(MInt countNode3 = countNode2 + 1; countNode3 < nfound; ++countNode3) {
            MBool check = false;
            for(MInt bc2 = 0; bc2 < nDim; ++bc2) {
              if(nodeMap[results[countNode2]].BC[bc2] == 0) break;
              for(MInt bc3 = 0; bc3 < nDim; ++bc3) {
                if(nodeMap[results[countNode2]].BC[bc2] == nodeMap[results[countNode3]].BC[bc3]) {
                  check = true;
                  // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
                  if(addConnection(
                         nodeMap[results[countNode2]].BC[bc2], // Does this make sense? Check! //6000, //6000er-change
                         nodeMap[results[countNode2]].blockId,
                         nodeMap[results[countNode2]].pos,
                         nodeMap[results[countNode3]].blockId,
                         nodeMap[results[countNode3]].pos)) {
                    numConnections = numConnections + 1;
                  }
                }
              }
            }
            if(!check) {
              cout << "DANGER: x|y|z=" << coordinates(0, i) << "|" << coordinates(1, i) << "|" << coordinates(2, i)
                   << " nodeMap2: inpuBlockId=" << nodeMap[results[countNode2]].blockId << " BC =";
              for(MInt d = 0; d < nDim; ++d)
                cout << " " << nodeMap[results[countNode2]].BC[d];
              cout << " nodeMap3: blockId=" << nodeMap[results[countNode3]].blockId << " BC =";
              for(MInt d = 0; d < nDim; ++d)
                cout << " " << nodeMap[results[countNode3]].BC[d];
              cout << endl;
              //              mTerm(1, "");
            }

            /*            if (nodeMap[results[countNode2]].BC!=nodeMap[results[countNode3]].BC) {
                          // By now the higher BC has a higher priority; Check if this makes sense
                          nodeMap[results[countNode2]].BC = mMax(nodeMap[results[countNode2]].BC,
                                                                   nodeMap[results[countNode3]].BC);
                          cout << "WARNING: 2 nodes with different BCs!!!" << endl;
                        }
                        //CHANGE_SET: currently addConnection will always return true, because we are now using a
               multiset if (addConnection( nodeMap[results[countNode2]].BC, // Does this make sense? Check! //6000,
               //6000er-change nodeMap[results[countNode2]].blockId, nodeMap[results[countNode2]].pos,
                                           nodeMap[results[countNode3]].blockId,
                                           nodeMap[results[countNode3]].pos)) {
                          numConnections = numConnections + 1;
                        }*/
          }
          nodeMap[results[countNode2]].found = true;
        }

        std::set<MInt> BCs;
        for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
          for(MInt d = 0; d < nDim; ++d) {
            if(nodeMap[results[countNode2]].BC[d] != 0) {
              BCs.insert(nodeMap[results[countNode2]].BC[d]);
            }
          }
        }

        // if three points share the same
        // coordinate it must be a 3-star
        if(nfound == 3) {
          add = false;
          tempcount = 0;

          // check if it is a singularity 3-star or normal 3-point connection
          for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            for(MInt j = 0; j < m_noBlocks; ++j) {
              if(nodeMap[results[countNode2]].blockId == j) {
                tempnum = j;
                break;
              }
            }
            for(MInt j = 0; j < 3; ++j) {
              // pos[]: ijk  DirLast[]: kji
              if(nodeMap[results[countNode2]].pos[j] == 0
                 || nodeMap[results[countNode2]].pos[j] == m_grid->getBlockNoCells(tempnum, 2 - j)) {
                ++tempcount;
              }
            }
          }

          if(tempcount == 9 || tempcount == 6) {
            add = true;
          }

          if(add) {
            for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
              if(BCs.size() > 1) cout << "CAUTION: This singular point has more than one BC!" << endl;
              for(const auto& BC : BCs) {
                // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
                if(addConnection(BC, // nodeMap[results[countNode2]].BC[0], //6000, //6000er-change
                                 nodeMap[results[countNode2]].blockId,
                                 nodeMap[results[countNode2]].pos,
                                 nodeMap[results[countNode2]].blockId,
                                 nodeMap[results[countNode2]].pos,
                                 3)) {
                  numSingularConnections = numSingularConnections + 1;
                }
              }
            }
          }
        }

        // if three points share the same
        // coordinate it must be a 5-star
        if(nfound == 5) {
          for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            if(BCs.size() > 1) cout << "CAUTION: This singular point has more than one BC!" << endl;
            for(const auto& BC : BCs) {
              // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
              if(addConnection(BC, // nodeMap[results[countNode2]].BC[0], //6000, //6000er-change
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               5)) {
                numSingularConnections = numSingularConnections + 1;
              }
            }
          }
        }
        // coordinate it must be a 6-star
        if(nfound == 6) {
          for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            if(BCs.size() > 1) cout << "CAUTION: This singular point has more than one BC!" << endl;
            for(const auto& BC : BCs) {
              // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
              if(addConnection(BC, // nodeMap[results[countNode2]].BC[0], //6000, //6000er-change
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               6)) {
                numSingularConnections = numSingularConnections + 1;
              }
            }
          }
        }
      }
    }


    MFloat tmppoint[3];
    MInt pcount = 0, plocation = 0, numofpoints[3];

    for(MInt i = 0; i < noInputWindowInformation; i++) {
      if(inputWindows[i]->BC > 4000 && inputWindows[i]->BC < 5000) {
        for(MInt j = 0; j < 3; j++) {
          numofpoints[j] = inputWindows[i]->endindex[j] - inputWindows[i]->startindex[j] + 1;
        }
        pcount += numofpoints[0] * numofpoints[1] * numofpoints[2];
      }
    }

    if(pcount > 0) {
      MFloatScratchSpace periodicCoordinates(3, pcount, AT_, "periodicCoordinates");
      periodicCoordinates.fill(-1.01010101);
      std::vector<pointType<nDim>> nodeMapP;
      nodeMapP.resize(pcount);

      if(domainId() == 0) {
        cout << "Loading periodic face coordinates" << endl;
      }
      MInt pmemsize = 0;
      for(MInt i = 0; i < noInputWindowInformation; i++) {
        if(inputWindows[i]->BC >= 4000 && inputWindows[i]->BC < 5000) {
          MString bName = "/block";
          stringstream number;
          number << inputWindows[i]->blockId << "/";
          bName += number.str();

          ParallelIo::size_type ioOffset[3] = {0, 0, 0};
          ParallelIo::size_type ioSize[3] = {0, 0, 0};
          if(m_domainId == 0) {
            for(MInt j = 0; j < 3; j++) {
              ioOffset[j] = inputWindows[i]->startindex[2 - j];
              ioSize[j] = inputWindows[i]->endindex[2 - j] - inputWindows[i]->startindex[2 - j] + 1;
            }

            pio.readArray(&periodicCoordinates(0, plocation), bName, "x", 3, ioOffset, ioSize);
            pio.readArray(&periodicCoordinates(1, plocation), bName, "y", 3, ioOffset, ioSize);
            pio.readArray(&periodicCoordinates(2, plocation), bName, "z", 3, ioOffset, ioSize);
          } else {
            ioOffset[2] = 0;
            ioSize[2] = 0; // i
            ioOffset[1] = 0;
            ioSize[1] = 0; // j
            ioOffset[0] = 0;
            ioSize[0] = 0; // k
            MFloat empty = 0;
            pio.readArray(&empty, bName, "x", 3, ioOffset, ioSize);
            pio.readArray(&empty, bName, "y", 3, ioOffset, ioSize);
            pio.readArray(&empty, bName, "z", 3, ioOffset, ioSize);
          }
          for(int dim = 0; dim < 3; ++dim) {
            offset[dim] = ioOffset[dim];
            size[dim] = ioSize[dim];
          }


          for(MInt c = offset[0]; c < offset[0] + size[0]; ++c) {
            for(MInt b = offset[1]; b < offset[1] + size[1]; ++b) {
              for(MInt a = offset[2]; a < offset[2] + size[2]; ++a) {
                nodeMapP[plocation].BC[0] = inputWindows[i]->BC;
                nodeMapP[plocation].blockId = inputWindows[i]->blockId;
                nodeMapP[plocation].pos[0] = a;
                nodeMapP[plocation].pos[1] = b;
                nodeMapP[plocation].pos[2] = c;
                plocation++;
              }
            }
          }
          pmemsize += size[0] * size[1] * size[2];
        }
      }

      MPI_Bcast(&plocation, 1, MPI_INT, 0, m_StructuredComm, AT_, "plocation");
      MPI_Bcast(&pmemsize, 1, MPI_INT, 0, m_StructuredComm, AT_, "pmemsize");

      MPI_Bcast(&periodicCoordinates(0, 0), pmemsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "periodicCoordinates(0");
      MPI_Bcast(&periodicCoordinates(1, 0), pmemsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "periodicCoordinates(1");
      MPI_Bcast(&periodicCoordinates(2, 0), pmemsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "periodicCoordinates(2");

      MPI_Bcast(nodeMapP.data(), pmemsize, matrix, 0, m_StructuredComm, AT_, "nodeMapP");

      if(domainId() == 0) {
        cout << "Computing periodic window displacements" << endl;
      }
      m_log << "Computing periodic window displacements" << endl;
      /////////////////////////////////////////////////
      //////// PERIODIC WINDOW DISPLACEMENTS //////////
      /////////////////////////////////////////////////
      MFloatScratchSpace periodicWindowCenter(m_noBlocks, nDim, 2 * nDim, AT_, "periodicWindowCOGS");
      MIntScratchSpace periodicWindowNoNodes(m_noBlocks, 2 * nDim, AT_, "periodicWindowNoNodes");
      cout.precision(10);
      periodicWindowCenter.fill(F0);
      periodicWindowNoNodes.fill(0);

      // compute displacement
      const MInt periodicOffset = 4401;
      for(MInt i = 0; i < pcount; ++i) {
        if(nodeMapP[i].BC[0] < 4401 || nodeMapP[i].BC[0] > 4406) {
          continue;
        }

        // add up all coordinates
        for(MInt dim = 0; dim < nDim; dim++) {
          periodicWindowCenter(nodeMapP[i].blockId, dim, nodeMapP[i].BC[0] - periodicOffset) +=
              periodicCoordinates(dim, i);
        }

        periodicWindowNoNodes(nodeMapP[i].blockId, nodeMapP[i].BC[0] - periodicOffset) += 1;
      }

      // new approach to account also for block distribution which are periodic in spanwise direction
      // and splited in spanwise direction

      // 1) Compute center of the the faces and split even and odd sides
      multimap<MInt, pair<vector<MFloat>, MInt>> oddPeriodicSides;
      multimap<MInt, pair<vector<MFloat>, MInt>> evenPeriodicSides;
      // vector< tuple<MInt, vector<MInt>, MInt>> oddPeriodicSides;
      // vector< tuple<MInt, vector<MInt>, MInt>> evenPeriodicSides;
      // evenPeriodicSides.reserve(m_noBlocks*nDim); //per block we can only have max 3 odd sides
      // oddPeriodicSides.reserve(m_noBlocks*nDim); //per block we can only have mac 3 even sides

      // now compute the center
      //(not weighted, should be basically the same surface in another location)
      // and compute distance between the two centers
      for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
        for(MInt periodicWindowId = 0; periodicWindowId < 2 * nDim; periodicWindowId++) {
          if(periodicWindowNoNodes(blockId, periodicWindowId) <= 0) { // we have no periodic cells
            continue;
          }
          vector<MFloat> centerDimensions(nDim);
          for(MInt dim = 0; dim < nDim; dim++) { // add the periodic cells
            periodicWindowCenter(blockId, dim, periodicWindowId) /=
                (MFloat)periodicWindowNoNodes(blockId, periodicWindowId);
            centerDimensions[dim] =
                periodicWindowCenter(blockId, dim,
                                     periodicWindowId); //(MFloat)periodicWindowNoNodes(blockId,periodicWindowId);
            cout << "periodic Window center " << periodicWindowCenter(blockId, dim, periodicWindowId) << endl;
          }
          if(periodicWindowId % 2 == 0) { // we have an even side
            // evenPeriodicSides.push_back(make_tuple(blockId,centerDimensions, periodicWindowId));
            evenPeriodicSides.insert(pair<MInt, pair<vector<MFloat>, MInt>>(
                blockId, pair<vector<MFloat>, MInt>(centerDimensions, periodicWindowId)));
            cout << "dimensions of centerDimensions (even)" << centerDimensions[0] << " " << centerDimensions[1] << " "
                 << centerDimensions[2] << endl;
          } else { // we have an odd side
            // oddPeriodicSides.push_back(make_tuple(blockId,centerDimensions, periodicWindowId));
            oddPeriodicSides.insert(pair<MInt, pair<vector<MFloat>, MInt>>(
                blockId, pair<vector<MFloat>, MInt>(centerDimensions, periodicWindowId)));
            cout << "dimensions of centerDimensions (odd)" << centerDimensions[0] << " " << centerDimensions[1] << " "
                 << centerDimensions[2] << endl;
          }
        }
      }

      // 2) now set the perodic sides for the blockId
      // a) check for periodic sides in the same blockId

      // get the iterator start and end of the m_block
      pair<multimap<MInt, pair<vector<MFloat>, MInt>>::iterator, multimap<MInt, pair<vector<MFloat>, MInt>>::iterator>
          evenIt;
      pair<multimap<MInt, pair<vector<MFloat>, MInt>>::iterator, multimap<MInt, pair<vector<MFloat>, MInt>>::iterator>
          oddIt;


      MFloatScratchSpace periodicDisplacementsBlockIds(m_noBlocks, nDim * nDim, AT_,
                                                       "periodic Displacements per block");

      for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
        MInt noEvenPerodicWindowsBlockIds = evenPeriodicSides.count(blockId);
        MInt noOddPeriodicWindowsBlockIds = oddPeriodicSides.count(blockId);
        cout << "we have evenSides " << noEvenPerodicWindowsBlockIds << " and we have oddSides "
             << noOddPeriodicWindowsBlockIds << " in BlockId " << blockId << endl;
        if(noEvenPerodicWindowsBlockIds == 0) continue; // there are no periodic connections
        if(noOddPeriodicWindowsBlockIds == 0)
          continue; // there are no periodic connections associated with that blockId
        // else we have a interblock connection (much easier to find
        evenIt = evenPeriodicSides.equal_range(blockId);
        multimap<MInt, pair<vector<MFloat>, MInt>>::iterator eit = evenIt.first;
        while(eit != evenIt.second) {
          MBool tobeErased = false;
          MInt evenSideDim =
              ceil(((MFloat)((*eit).second.second + 1) / 2.0) - 1); // get the dimensional direction of the side
          MInt evenSide = (*eit).second.second;
          MInt oddSide = evenSide + 1;
          oddIt = oddPeriodicSides.equal_range(blockId);
          multimap<MInt, pair<vector<MFloat>, MInt>>::iterator oit = oddIt.first;

          while(oit != oddIt.second) {
            if((*oit).second.second == oddSide) { // yes we found a corresponding odd side
              // compute the periodic connection and put them into the scratchspace
              for(MInt dim = 0; dim < nDim; dim++) {
                cout << "first element is " << (*oit).second.first[dim] << " and we subtract "
                     << (*eit).second.first[dim] << endl;
                periodicDisplacementsBlockIds(blockId, evenSideDim + nDim * dim) =
                    abs((*oit).second.first[dim] - (*eit).second.first[dim]);
              }
              oit = oddPeriodicSides.erase(oit);
              tobeErased = true;
              break;
            } else {
              oit++;
            }
          }
          if(tobeErased) {
            eit = evenPeriodicSides.erase(eit);
          } else {
            eit++;
          }
        }
      }

      MFloat periodicEpsilon = 0.000001;
      // now the list only contains the inner blockId periodic conditions.
      // check for multiblock communications
      if(evenPeriodicSides.size() != 0
         && oddPeriodicSides.size() != 0) { // we did not delete everything in the side container
        for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
          MInt noEvenPerodicWindowsBlockId = evenPeriodicSides.count(blockId);
          if(noEvenPerodicWindowsBlockId == 0) continue; // there are no periodic connections left in the even container
          evenIt = evenPeriodicSides.equal_range(blockId);
          multimap<MInt, pair<vector<MFloat>, MInt>>::iterator eit = evenIt.first;
          while(eit != evenIt.second) {
            MBool tobeErased = false;
            MInt evenSideDim =
                ceil(((MFloat)((*eit).second.second + 1) / 2.0) - 1); // get the dimensional direction of the side
            // MInt evenSide    = (*eit).second.second;
            multimap<MInt, pair<vector<MFloat>, MInt>>::iterator oit = oddPeriodicSides.begin();
            while(oit != oddPeriodicSides.end()) {
              // check the number of conncetions which have the same cooridnates
              MInt equalConnectionDims = 0;
              for(MInt d = 0; d < nDim; d++) {
                if(abs((*oit).second.first[d] - (*eit).second.first[d]) < periodicEpsilon) equalConnectionDims++;
              }
              if(equalConnectionDims == 2) { // yes we found a pair belonging together
                for(MInt dim = 0; dim < nDim; dim++) {
                  periodicDisplacementsBlockIds(blockId, evenSideDim + nDim * dim) =
                      abs((*oit).second.first[dim] - (*eit).second.first[dim]);
                  // now store the odd side also for multiblock communication
                  periodicDisplacementsBlockIds((*oit).first, evenSideDim + nDim * dim) =
                      periodicDisplacementsBlockIds(blockId, evenSideDim + nDim * dim);
                }
                tobeErased = true;
                oit = oddPeriodicSides.erase(oit);
              } else {
                oit++;
              }
            }
            if(tobeErased) {
              eit = evenPeriodicSides.erase(eit);
            } else {
              eit++;
            }
          }
        }
      }


      for(MInt periodicWindowId = 0; periodicWindowId < 2 * nDim; periodicWindowId++) {
        if(periodicWindowId % 2 == 1) {
          for(MInt dim = 0; dim < nDim; dim++) {
            const MInt displacementId = (MFloat)(periodicWindowId + 1) / 2.0 - 1;
            periodicDisplacements[dim * nDim + displacementId] =
                periodicWindowCenter(m_blockId, dim, periodicWindowId)
                - periodicWindowCenter(m_blockId, dim, periodicWindowId - 1);
            m_log << m_blockId << " periodicWindowCenter for window: " << periodicWindowId + periodicOffset
                  << "  dim: " << dim << " displacementId: " << displacementId
                  << " displacement: " << periodicDisplacements[dim * nDim + displacementId] << endl;
          }
        }
      }

      cout << "old approach" << endl;
      for(MInt periodicWindowId = 0; periodicWindowId < 2 * nDim; periodicWindowId++) {
        if(periodicWindowId % 2 == 1) {
          for(MInt dim = 0; dim < nDim; dim++) {
            const MInt displacementId = (MFloat)(periodicWindowId + 1) / 2.0 - 1;
            periodicDisplacements[dim * nDim + displacementId] =
                periodicWindowCenter(m_blockId, dim, periodicWindowId)
                - periodicWindowCenter(m_blockId, dim, periodicWindowId - 1);
            cout << "blockId = " << m_blockId
                 << "periodicWindowCenter for window: " << periodicWindowId + periodicOffset << "  dim: " << dim
                 << " displacementId: " << displacementId
                 << " displacement: " << periodicDisplacements[dim * nDim + displacementId] << endl;
          }
        }
      }
      cout << "new approach " << endl;
      for(MInt periodicWindowId = 0; periodicWindowId < 2 * nDim; periodicWindowId++) {
        if(periodicWindowId % 2 == 1) {
          const MInt displacementId = (MFloat)(periodicWindowId + 1) / 2.0 - 1;
          for(MInt dim = 0; dim < nDim; dim++) {
            cout << "blockId = " << m_blockId
                 << "periodicWindowCenter for window: " << periodicWindowId + periodicOffset << "  dim: " << dim
                 << " displacementId: " << displacementId
                 << " displacement: " << periodicDisplacementsBlockIds(m_blockId, dim * nDim + displacementId) << endl;
          }
        }
      }
      MPI_Barrier(m_StructuredComm, AT_);
      if(domainId() == 0) {
        for(MInt b = 0; b < m_noBlocks; b++) {
          cout << "BlockId " << b << ": " << endl;
          cout << "\t";
          for(MInt i = 0; i < nDim * nDim; i++) {
            cout << periodicDisplacementsBlockIds(b, i) << " ";
          }
          cout << endl;
        }
      }
      MPI_Barrier(m_StructuredComm, AT_);
      for(MInt b = 0; b < m_noBlocks; b++) {
        if(b == m_blockId) {
          for(MInt i = 0; i < nDim * nDim; i++) {
            periodicDisplacements[i] = periodicDisplacementsBlockIds(b, i);
          }
        }
      }
      // mTerm(-1, AT_, "I want to terminate ");

      /////////////////////////////////////////////////
      //////// PERIODIC CONNECTIONS ///////////////////

      /////////////////////////////////////////////////

      for(MInt i = 0; i < pcount; ++i) {
        tmppoint[0] = periodicCoordinates(0, i);
        tmppoint[1] = periodicCoordinates(1, i);
        tmppoint[2] = periodicCoordinates(2, i);
        periodicPointsChange(tmppoint, nodeMapP[i].BC[0], periodicDisplacements);
        Point<3> aaa(tmppoint[0], tmppoint[1], tmppoint[2]);
        nfound = tree.locatenear(aaa, m_gridEps, results, 10, false); // 0.0001
        for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
          // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
          if(addConnection(nodeMapP[i].BC[0],
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           nodeMap[results[countNode2]].blockId,
                           nodeMap[results[countNode2]].pos)) {
            numConnections = numConnections + 1;
          }
        }

        if(nfound == 5) {
          // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
          if(addConnection(nodeMapP[i].BC[0], nodeMapP[i].blockId, nodeMapP[i].pos, nodeMapP[i].blockId,
                           nodeMapP[i].pos, 5)) {
            numSingularConnections = numSingularConnections + 1;
          }
        }
        // 6 star for periodic bc
        if(nfound == 6) {
          // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
          if(addConnection(nodeMapP[i].BC[0], nodeMapP[i].blockId, nodeMapP[i].pos, nodeMapP[i].blockId,
                           nodeMapP[i].pos, 6)) {
            numSingularConnections = numSingularConnections + 1;
          }
        }
      }
    }

    if(domainId() == 0) {
      cout << "Assemble multiblock connections" << endl;
    }
    ///////////////////////////////////////////////////
    /////// CREATE THE WINDOWS FROM CONNECTIONS ///////
    ///////////////////////////////////////////////////
    multiBlockAssembling();
    singularityAssembling();

    // Delete duplicate windows
    deleteDuplicateWindows(window0d, window0d);
    deleteDuplicateWindows(window1d, window1d);
    deleteDuplicateWindows(window2d, window2d);
    deleteDuplicateWindows(window0d, window1d);
    deleteDuplicateWindows(window0d, window2d);
    deleteDuplicateWindows(window1d, window2d);

    //

    for(auto it = singularwindow.cbegin(); it != singularwindow.cend();) {
      if(singularwindow.empty()) {
        break;
      }
      // TODO_SS labels:FV by now don't do anything if BC is perioduc
      if((*it)->BC >= 4400 && (*it)->BC < 4410) continue;
      MInt noSameWindows = 0;
      (*it)->BCsingular[noSameWindows++] = (*it)->BC;
      (*it)->BC = 0;
      for(auto it2 = it + 1; it2 != singularwindow.cend();) {
        if(singularwindow.empty()) {
          break;
        }
        const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
        const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
        const MBool test = mapCompare11(map1, map2);
        if(test) {
          if((*it)->BCsingular[0] == (*it2)->BC || (*it)->Nstar != (*it2)->Nstar) {
            mTerm(1, "There is no hope anymore!");
          }

          if(noSameWindows >= (*it)->Nstar) {
            mTerm(1, "");
          }

          (*it)->BCsingular[noSameWindows++] = (*it2)->BC;
          singularwindow.erase(it2);
        } else {
          ++it2;
        }
      }
      ++it;
    }


    // Create global BC maps for those singularities which don't have one
    /*    for (MInt i = 0; i < (signed)singularwindow.size(); ++i) {
          const auto& mymap = singularwindow[i];
          for (MInt bc = 0; bc < mymap->Nstar; ++bc) {
            if (mymap->BCsingular[bc]!=0 && mymap->BCsingular[bc]!=6000) {
              MBool addMap = true;
              for (MInt j = 0; j < (signed)globalStructuredBndryCndMaps.size(); ++j) {
                if (mymap->Id1==globalStructuredBndryCndMaps[j]->Id1 &&
       mymap->BCsingular[bc]==globalStructuredBndryCndMaps[j]->BC) { MBool addMap_temp = false; for (MInt d = 0; d <
       nDim; ++d) { if (globalStructuredBndryCndMaps[j]->start1[d]>mymap->start1[d] ||
       globalStructuredBndryCndMaps[j]->end1[d]<mymap->end1[d]) { addMap_temp = true;
                    }
                  }
                  addMap = addMap_temp;
                  if (!addMap) break;
                }
              }
              if (addMap) {
                StructuredWindowMap* windowMap=new StructuredWindowMap(nDim);
                MInt order[3]={0,1,2};
                MInt step1[3]={1,1,1};
                mapCreate(mymap->Id1, mymap->start1, mymap->end1, step1,
                          -1 ,  NULL, NULL , NULL, order,mymap->BCsingular[bc], windowMap );
                windowMap->Nstar = mymap->Nstar;
                //add the map to the list!
                if (domainId()==0) {
                  cout << "ADDDING following GlobalBCMap:" << endl;
                  mapPrint(windowMap);
                }
                globalStructuredBndryCndMaps.push_back(std::move(windowMap));
              }
            }
          }
        }*/

    // StructuredWindowMap* windowMap;
    for(MInt i = 0; i < (MInt)window2d.size(); ++i) {
      // To distinguish between communication bcs and normal ones, do the following, negate the communication bcs
      const MInt BC = (window2d[i]->BC >= 6000 && window2d[i]->BC < 6010) ? -window2d[i]->BC : window2d[i]->BC;

      unique_ptr<StructuredWindowMap<nDim>> windowMap1 = make_unique<StructuredWindowMap<nDim>>();
      mapCreate(window2d[i]->Id1, window2d[i]->start1, window2d[i]->end1, window2d[i]->step1, window2d[i]->Id2,
                window2d[i]->start2, window2d[i]->end2, window2d[i]->step2, window2d[i]->order, BC, windowMap1);
      windowMap1->Nstar = -1;
      windowMap1->SingularId = -1;
      windowMap1->dc1 = window2d[i]->dc1;
      windowMap1->dc2 = window2d[i]->dc2;

      // i-surface: change sign of I; same for J and K
      if(windowMap1->dc1 * windowMap1->dc2 > 0) {
        MInt tempdim = abs(windowMap1->dc2) - 1;
        windowMap1->step2[tempdim] = -1;
      }

      // add the map to the list!
      globalStructuredBndryCndMaps.push_back(std::move(windowMap1));


      if(window2d[i]->BC >= 6000 && window2d[i]->BC < 6010) { //( window2d[i]->BC==6000) { //6000er-change
        unique_ptr<StructuredWindowMap<nDim>> windowMap = make_unique<StructuredWindowMap<nDim>>();
        mapCreate(window2d[i]->Id1, window2d[i]->start1, window2d[i]->end1, window2d[i]->step1, window2d[i]->Id2,
                  window2d[i]->start2, window2d[i]->end2, window2d[i]->step2, window2d[i]->order, BC, windowMap);
        windowMap->Nstar = -1;
        windowMap->SingularId = -1;
        windowMap->dc1 = window2d[i]->dc1;
        windowMap->dc2 = window2d[i]->dc2;
        // i-surface: change sign of I; same for J and K
        if(windowMap->dc1 * windowMap->dc2 > 0) {
          MInt tempdim = abs(windowMap->dc2) - 1;
          windowMap->step2[tempdim] = -1;
        }

        mapInvert1(windowMap);
        mapNormalize3(windowMap);
        globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      }
    }


    for(MInt i = 0; i < (MInt)window1d.size(); ++i) {
      // To distinguish between communication bcs and normal ones, do the followin, negate the communication bcs
      const MInt BC = (window1d[i]->BC >= 6000 && window1d[i]->BC < 6010) ? -window1d[i]->BC : window1d[i]->BC;

      unique_ptr<StructuredWindowMap<nDim>> windowMap1 = make_unique<StructuredWindowMap<nDim>>();
      mapCreate(window1d[i]->Id1, window1d[i]->start1, window1d[i]->end1, window1d[i]->step1, window1d[i]->Id2,
                window1d[i]->start2, window1d[i]->end2, window1d[i]->step2, window1d[i]->order, BC, windowMap1);
      windowMap1->Nstar = -1;
      windowMap1->SingularId = -1;
      windowMap1->dc1 = window1d[i]->dc1;
      windowMap1->dc2 = window1d[i]->dc2;

      // i-surface: change sign of I; same for J and K
      for(MInt j = 0; j < 3; ++j) {
        if(windowMap1->start1[j] == windowMap1->end1[j]) {
          if(windowMap1->start1[j] == 0 && windowMap1->start2[windowMap1->order[j]] == 0) {
            windowMap1->step2[windowMap1->order[j]] = -1;
          }
          if(windowMap1->start1[j] > 0 && windowMap1->start2[windowMap1->order[j]] > 0) {
            windowMap1->step2[windowMap1->order[j]] = -1;
          }
        }
      }

      // now for 5 star communicaitons (or ...... if needed)
      // 3 star does not need additional information
      for(MInt j = 0; j < (MInt)singularwindow.size(); ++j) {
        const MInt test = mapCompare11(windowMap1, singularwindow[j]);
        if(test) {
          windowMap1->Nstar = singularwindow[j]->Nstar;
          windowMap1->SingularId = singularwindow[j]->SingularId;
        }
      }

      // add the map to the list!
      globalStructuredBndryCndMaps.push_back(std::move(windowMap1));

      if(window1d[i]->BC >= 6000 && window1d[i]->BC < 6010) { //( window1d[i]->BC==6000) { //6000er-change
        unique_ptr<StructuredWindowMap<nDim>> windowMap = make_unique<StructuredWindowMap<nDim>>();
        mapCreate(window1d[i]->Id1, window1d[i]->start1, window1d[i]->end1, window1d[i]->step1, window1d[i]->Id2,
                  window1d[i]->start2, window1d[i]->end2, window1d[i]->step2, window1d[i]->order, BC, windowMap);
        windowMap->Nstar = -1;
        windowMap->SingularId = -1;
        windowMap->dc1 = window1d[i]->dc1;
        windowMap->dc2 = window1d[i]->dc2;

        // i-surface: change sign of I; same for J and K
        for(MInt j = 0; j < 3; ++j) {
          if(windowMap->start1[j] == windowMap->end1[j]) {
            if(windowMap->start1[j] == 0 && windowMap->start2[windowMap->order[j]] == 0) {
              windowMap->step2[windowMap->order[j]] = -1;
            }
            if(windowMap->start1[j] > 0 && windowMap->start2[windowMap->order[j]] > 0) {
              windowMap->step2[windowMap->order[j]] = -1;
            }
          }
        }

        mapInvert1(windowMap);
        mapNormalize3(windowMap);

        // now for 5 star communicaitons (or ...... if needed)
        // 3 star does not need additional information
        for(MInt j = 0; j < (MInt)singularwindow.size(); ++j) {
          const MInt test = mapCompare11(windowMap, singularwindow[j]);
          if(test) {
            windowMap->Nstar = singularwindow[j]->Nstar;
            windowMap->SingularId = singularwindow[j]->SingularId;
          }
        }

        globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      }
    }


    for(MInt i = 0; i < (MInt)window0d.size(); ++i) {
      // To distinguish between communication bcs and normal ones, do the followin, negate the communication bcs
      const MInt BC = (window0d[i]->BC >= 6000 && window0d[i]->BC < 6010) ? -window0d[i]->BC : window0d[i]->BC;

      // point communicaiton (rare)
      unique_ptr<StructuredWindowMap<nDim>> windowMap1 = make_unique<StructuredWindowMap<nDim>>();
      mapCreate(window0d[i]->Id1, window0d[i]->start1, window0d[i]->end1, window0d[i]->step1, window0d[i]->Id2,
                window0d[i]->start2, window0d[i]->end2, window0d[i]->step2, window0d[i]->order, BC, windowMap1);
      windowMap1->Nstar = -1;
      windowMap1->SingularId = -1;
      windowMap1->dc1 = window0d[i]->dc1;
      windowMap1->dc2 = window0d[i]->dc2;
      // i-surface: change sign of I; same for J and K
      for(MInt j = 0; j < 3; ++j) {
        if(windowMap1->start1[j] == windowMap1->end1[j]) {
          if(windowMap1->start1[j] == 0 && windowMap1->start2[windowMap1->order[j]] == 0) {
            windowMap1->step2[windowMap1->order[j]] = -1;
          }
          if(windowMap1->start1[j] > 0 && windowMap1->start2[windowMap1->order[j]] > 0) {
            windowMap1->step2[windowMap1->order[j]] = -1;
          }
        }
      }

      // now for 5 star communicaitons (or ...... if needed)
      // 3 star does not need additional information
      for(MInt j = 0; j < (MInt)singularwindow.size(); ++j) {
        const MInt test = mapCompare11(windowMap1, singularwindow[j]);
        if(test) {
          windowMap1->Nstar = singularwindow[j]->Nstar;
          windowMap1->SingularId = singularwindow[j]->SingularId;
        }
      }

      // add the map to the list!
      globalStructuredBndryCndMaps.push_back(std::move(windowMap1));

      if(window0d[i]->BC >= 6000 && window0d[i]->BC < 6010) { //( window0d[i]->BC==6000) { //6000er-change
        unique_ptr<StructuredWindowMap<nDim>> windowMap = make_unique<StructuredWindowMap<nDim>>();
        mapCreate(window0d[i]->Id1, window0d[i]->start1, window0d[i]->end1, window0d[i]->step1, window0d[i]->Id2,
                  window0d[i]->start2, window0d[i]->end2, window0d[i]->step2, window0d[i]->order, BC, windowMap);
        windowMap->Nstar = -1;
        windowMap->SingularId = -1;
        windowMap->dc1 = window0d[i]->dc1;
        windowMap->dc2 = window0d[i]->dc2;

        // i-surface: change sign of I; same for J and K
        for(MInt j = 0; j < 3; ++j) {
          if(windowMap->start1[j] == windowMap->end1[j]) {
            if(windowMap->start1[j] == 0 && windowMap->start2[windowMap->order[j]] == 0) {
              windowMap->step2[windowMap->order[j]] = -1;
            }
            if(windowMap->start1[j] > 0 && windowMap->start2[windowMap->order[j]] > 0) {
              windowMap->step2[windowMap->order[j]] = -1;
            }
          }
        }

        mapInvert1(windowMap);
        mapNormalize3(windowMap);

        // now for 5 star communicaitons (or ...... if needed)
        // 3 star does not need additional information
        for(MInt j = 0; j < (MInt)singularwindow.size(); ++j) {
          const MInt test = mapCompare11(windowMap, singularwindow[j]);
          if(test) {
            windowMap->Nstar = singularwindow[j]->Nstar;
            windowMap->SingularId = singularwindow[j]->SingularId;
          }
        }
        globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      }
    }

    // Duplicates are already deleted above. The following call should be unnecessary
    deleteDuplicateCommMaps();

#ifndef NDEBUG
    if(domainId() == 0) {
      cout << "======= GLOBALSTRUCTUREDBNDRYCNDMAPS ======" << endl;
      for(MInt i = 0; i < (signed)globalStructuredBndryCndMaps.size(); ++i) {
        mapPrint(globalStructuredBndryCndMaps[i]);
      }

      cout << " ======= SINGULARWINDOW =============" << endl;
      for(MInt i = 0; i < (signed)singularwindow.size(); ++i) {
        mapPrint(singularwindow[i]);
      }

      cout << " ======== window0d ==============" << endl;
      for(MInt i = 0; i < (signed)window0d.size(); ++i) {
        mapPrint(window0d[i]);
      }

      cout << " ======== window1d ==============" << endl;
      for(MInt i = 0; i < (signed)window1d.size(); ++i) {
        mapPrint(window1d[i]);
      }

      cout << " ======== window2d ==============" << endl;
      for(MInt i = 0; i < (signed)window2d.size(); ++i) {
        mapPrint(window2d[i]);
      }
    }
#endif

    if(!hasConnectionInfo) {
      cout << "Writing connection info" << endl;
      writeConnectionWindowInformation3D(periodicDisplacements);
    } else {
      if(domainId() == 0) {
        cout << "##########################################################################" << endl;
        cout << "WARNING: OLD CONNECTION INFO EXISTS, DELETE BEFORE NEW ONE CAN BE WRITTEN!" << endl;
        cout << "##########################################################################" << endl;
      }
    }

    if(domainId() == 0) {
      cout << "Connection identification and window creation finished! NUM_SINGULAR=" << numSingularConnections << endl;
    }
  }
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::periodicPointsChange(MFloat* pt, MInt type, MFloat* periodicDisplacements) {
  MFloat angle = 0;
  /*
    case 4401 4402: first periodic direction
    case 4403 4404: second periodic direction
    case 4405 4405: third periodic direction
    case 4011 4012:  rotation X axis clockwise and anticlockwise
  */
  MFloat rotationMatrix[3][3]{};
  MFloat tmp[3]{};
  switch(type) {
    case 4401:
    case 4403:
    case 4405: {
      const MInt displacementId = (MFloat)(type - 4400 + 1) / 2.0 - 1;
      for(MInt dim = 0; dim < nDim; ++dim) {
        pt[dim] = pt[dim] + periodicDisplacements[dim * nDim + displacementId];
      }
      break;
    }

    case 4402:
    case 4404:
    case 4406: {
      const MInt displacementId = (MFloat)(type - 4400) / 2.0 - 1;
      for(MInt dim = 0; dim < nDim; ++dim) {
        pt[dim] = pt[dim] - periodicDisplacements[dim * nDim + displacementId];
      }
      break;
    }

    case 4011: {
      ASSERT(nDim == 3, "");
      rotationMatrix[1][1] = cos(angle);
      rotationMatrix[1][2] = -sin(angle);
      rotationMatrix[2][1] = sin(angle);
      rotationMatrix[2][2] = cos(angle);
      tmp[0] = pt[0];
      tmp[1] = pt[1];
      tmp[2] = pt[2];
      for(MInt i = 0; i < nDim; ++i) {
        pt[i] = rotationMatrix[i][0] * tmp[0] + rotationMatrix[i][1] * tmp[1] + rotationMatrix[i][2] * tmp[2];
      }
      break;
    }

    case 4012: {
      ASSERT(nDim == 3, "");
      rotationMatrix[1][1] = cos(-angle);
      rotationMatrix[1][2] = -sin(-angle);
      rotationMatrix[2][1] = sin(-angle);
      rotationMatrix[2][2] = cos(-angle);
      tmp[0] = pt[0];
      tmp[1] = pt[1];
      tmp[2] = pt[2];
      for(MInt i = 0; i < nDim; ++i) {
        pt[i] = rotationMatrix[i][0] * tmp[0] + rotationMatrix[i][1] * tmp[1] + rotationMatrix[i][2] * tmp[2];
      }
      break;
    }

    default: {
      cout << "ERROR!!! periodic type is wrong!!! in windowinfoe" << endl;
    }
  }
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::addConnection(MInt connectiontype, MInt b1, MInt* p1, MInt b2, MInt* p2) {
  // add connections for 6000 and 6000+
  //  pair<multiset<connectionNode>::iterator,MBool> ret; //CHANGE_SET

  if(connectiontype >= 6000 && connectiontype < 6010) { //(connectiontype==6000) { //6000er-change
    connectionNode a(connectiontype, b1, p1, b2, p2, true, nDim);
    /*ret=*/connectionset.insert(a);
  } else {
    connectionNode a(connectiontype, b1, p1, b2, p2, false, nDim);
    /*ret=*/connectionset.insert(a);
  }

  //  if(ret.second==false) { //CHANGE_SET
  //    cout<<"error! same connection can not be added! check the connections!"<<endl;
  //    return false;
  //  }

  return true;
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::addConnection(MInt connectiontype, MInt b1, MInt* p1, MInt b2, MInt* p2,
                                                        MInt Nstar) {
  // add special connections (3 or 5-star)
  //  pair<multiset<connectionNode>::iterator,MBool> ret; //CHANGE_SET

  if(connectiontype >= 6000 && connectiontype < 6010) { //(connectiontype==6000) { //6000er-change
    connectionNode a(connectiontype, b1, p1, b2, p2, true, nDim);
    a.Nstar = Nstar;
    /*ret=*/singularconnectionset.insert(a);
  } else if(connectiontype >= 4000 && connectiontype < 5000) {
    // periodic boundary condition only need to handle 5 star communications
    connectionNode a(connectiontype, b1, p1, b2, p2, false, nDim);
    a.Nstar = Nstar;
    /*ret=*/singularconnectionset.insert(a);
  } else {
    cout << "singular point on BC " << connectiontype << " is not supported!! check it!!!" << endl;
    exit(0);
  }

  //  if(ret.second==false) { //CHANGE_SET
  //    cout<<"error! same connection can not be added! check the connections!"<<endl;
  //    return false;
  //  }

  return true;
}


template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::findConnection(connectionNode a) {
  multiset<connectionNode>::iterator it; // CHANGE_SET
  it = connectionset.find(a);

  if(it != connectionset.end()) {
    return true;
  } else {
    return false;
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::removeConnection(connectionNode a) {
  //  connectionset.erase(a); //CHANGE_SET
  auto it = connectionset.find(a);
  connectionset.erase(it);
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::findConnection(connectionNode a, MInt Nstar) {
  multiset<connectionNode>::iterator it; // CHANGE_SET
  a.Nstar = Nstar;
  it = singularconnectionset.find(a);

  if(it != singularconnectionset.end()) {
    return true;
  } else {
    return false;
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::removeConnection(connectionNode a, MInt Nstar) {
  a.Nstar = Nstar;
  //  singularconnectionset.erase(a); //CHANGE_SET
  auto it = singularconnectionset.find(a);
  singularconnectionset.erase(it);
}


template <>
void FvStructuredSolverWindowInfo<2>::readWindowCoordinates(MFloat* periodicDisplacements) {
  constexpr MInt nDim = 2;

  MInt hasConnectionInfo = 0;
  if(globalDomainId() == 0) {
    ParallelIoHdf5 pio(m_grid->m_gridInputFileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);
    MBool attributeExists = pio.hasAttribute("hasConnectionInfo", "Connectivity");
    if(attributeExists) {
      m_log << "Grid file has connection info!" << endl;
      pio.getAttribute(&hasConnectionInfo, "hasConnectionInfo", "Connectivity");
    }
    MPI_Bcast(&hasConnectionInfo, 1, MPI_INT, 0, m_StructuredComm, AT_, "hasConnectionInfo");
  } else {
    MPI_Bcast(&hasConnectionInfo, 1, MPI_INT, 0, m_StructuredComm, AT_, "hasConnectionInfo");
  }

  MBool readInConnectionInfo = true;

  if(!hasConnectionInfo) {
    readInConnectionInfo = false;
  } else {
    if(Context::propertyExists("readInConnectionInfo", m_solverId)) {
      readInConnectionInfo =
          Context::getSolverProperty<MBool>("readInConnectionInfo", solverId(), AT_, &readInConnectionInfo);
    }
  }

  if(readInConnectionInfo) {
    m_log << "Connection info available, reading connection information from grid file!" << endl;
    readConnectionWindowInformation2D(periodicDisplacements);
  } else {
    if(domainId() == 0) {
      cout << "Starting grid connection info search..." << endl;
    }

    ParallelIoHdf5 pio(m_grid->m_gridInputFileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);
    m_log << " Starting grid connection info search..." << endl;
    MInt offset[2], size[2], count = 0;
    // 3) read in the coordinates of the grid points
    // open file for reading the grid data
    pointType<nDim>* nodeMap;
    MInt totalGridCells = 0;
    MInt** totalGridBlockDim = nullptr;
    MInt* totalGridBlockCells = nullptr;

    mAlloc(totalGridBlockDim, m_noBlocks, 2, "totalGridBlockDim", -1, AT_);
    mAlloc(totalGridBlockCells, m_noBlocks, "totalGridBlockCells", -1, AT_);

    for(MInt i = 0; i < m_noBlocks; ++i) {
      MInt temp = 1;
      for(MInt dim = 0; dim < nDim; ++dim) {
        temp *= m_grid->getBlockNoCells(i, dim) + 1;
        totalGridBlockDim[i][dim] = m_grid->getBlockNoCells(i, dim) + 1; // number of points in the grid File
      }
      totalGridBlockCells[i] = temp;
      if(temp != 1) totalGridCells += temp;
    }

    MInt countNode = 0;
    for(MInt i = 0; i < m_noBlocks; i++) {
      countNode = countNode + 2 * totalGridBlockDim[i][1] + 2 * totalGridBlockDim[i][0] - 2;
    }

    MFloatScratchSpace coordinates(2, countNode, AT_, "coordinates");
    coordinates.fill(-1.01010101);
    nodeMap = new pointType<nDim>[countNode];

    // auxiliary lambda functions to determine BC of interface points
    const MInt converter[] = {1, 0};
    auto getCurrentWindow = [this, converter](const MInt blockId, const MInt* const offset_, const MInt* const size_) {
      std::vector<MInt> currentWindows;
      for(MInt windowId = 0; windowId < noInputWindowInformation; ++windowId) {
        if(inputWindows[windowId]->blockId != blockId) continue;
        if(inputWindows[windowId]->BC < 6000 || inputWindows[windowId]->BC >= 6010) continue;
        // Besides the windows of the same face, windows of adjacent faces are also candidates for the corner points
        MBool takeIt = true;
        for(MInt dim = 0; dim < nDim; ++dim) {
          if(inputWindows[windowId]->startindex[dim] >= offset_[converter[dim]] + size_[converter[dim]]
             || inputWindows[windowId]->endindex[dim] < offset_[converter[dim]]) {
            takeIt = false;
            break;
          }
        }
        if(takeIt) currentWindows.push_back(std::move(windowId));
      }
      return currentWindows;
    };
    std::vector<MInt> pointsFoundPerWindow(noInputWindowInformation); // for detailed output
    auto getBCFromWindow = [this, &pointsFoundPerWindow](const MInt idx[nDim],
                                                         const std::vector<MInt>& currentWindows) {
      std::vector<MInt> BCset; // CHANGE_SET
      // Use set in order to have a unique list
      //      std::set<MInt> BCset;
      //      MInt BC[2] = {0, 0}; //we can have at most 2 different BCs shared by one point
      //      MInt cnt = 0;
      for(auto windowId : currentWindows) {
        MBool takeIt = true;
        for(MInt dim = 0; dim < nDim; ++dim) {
          if(inputWindows[windowId]->startindex[dim] > idx[dim] || inputWindows[windowId]->endindex[dim] < idx[dim]) {
            takeIt = false;
            break;
          }
        }
        if(takeIt) {
          ++pointsFoundPerWindow[windowId];
          //          BCset.insert(inputWindows[windowId]->BC);
          BCset.push_back(inputWindows[windowId]->BC); // CHANGE_SET
          //          BC[cnt++] = inputWindows[windowId]->BC;
        }
      }
      //      if (cnt>2) mTerm(1, "");
      //      if (cnt>0) return mMax(BC[0], BC[1]); //higher BC has higher priority; in future think of something
      //      more sophisticated return 6000; //default //-1;
      // we can have at most 2 different BCs shared by one point
      if(BCset.size() > 2) mTerm(1, "");
      if(BCset.size() == 0) BCset.push_back(6000); // BCset.insert(6000); //CHANGE_SET
      MInt BC[2] = {0, 0};
      MInt cnt = 0;
      for(auto it = BCset.begin(); it != BCset.end(); ++it)
        BC[cnt++] = *it;
      return std::tuple<MInt, MInt>(BC[0], BC[1]);
    };

    MInt memsize = 0;
    // domain 0 reads the coordinates for all 4 sides of all blocks
    if(domainId() == 0) {
      cout << "Reading in all coordinates of all block faces" << endl;
    }
    for(MInt i = 0; i < m_noBlocks; ++i) {
      MString bName = "/block";
      stringstream number;
      number << i << "/";
      bName += number.str();

      ParallelIo::size_type ioOffset[2] = {0, 0};
      ParallelIo::size_type ioSize[2] = {0, 0};
      // in grid file 2:i  1:j   0:k
      //////////////////////////////
      ///////// FACE 1 /////////////
      //////////////////////////////
      // face +x  face 1
      if(m_domainId == 0) {
        ioOffset[1] = 0;
        ioSize[1] = 1; // i
        ioOffset[0] = 0;
        ioSize[0] = m_grid->getBlockNoCells(i, 0) + 1; // j

        cout << "offset[1]=" << ioOffset[1] << " offset[0]=" << ioOffset[0] << endl;
        cout << "size[1]=" << ioSize[1] << " size[0]=" << ioSize[0] << endl;
        pio.readArray(&coordinates(0, memsize), bName, "x", 2, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 2, ioOffset, ioSize);
      } else {
        ioOffset[1] = 0;
        ioSize[1] = 0; // i
        ioOffset[0] = 0;
        ioSize[0] = 0; // j
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 2, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 2, ioOffset, ioSize);
      }
      for(int dim = 0; dim < nDim; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      std::vector<MInt> currentWindows = getCurrentWindow(i, offset, size);

      // write all computational coordinates
      // of the face to nodeMap object
      for(MInt a = offset[1]; a < offset[1] + size[1]; ++a) {
        for(MInt b = offset[0]; b < offset[0] + size[0]; ++b) {
          nodeMap[count].blockId = i;
          nodeMap[count].pos[0] = a;
          nodeMap[count].pos[1] = b;
          std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1]) = getBCFromWindow(nodeMap[count].pos, currentWindows);
          count++;
        }
      }
      memsize += size[0] * size[1];

      //////////////////////////////
      ///////// FACE 2 /////////////
      //////////////////////////////
      // face -x  face 2
      if(m_domainId == 0) {
        ioOffset[1] = m_grid->getBlockNoCells(i, 1);
        ioSize[1] = 1; // i
        ioOffset[0] = 0;
        ioSize[0] = m_grid->getBlockNoCells(i, 0) + 1; // j

        cout << "offset[1]=" << ioOffset[1] << " offset[0]=" << ioOffset[0] << endl;
        cout << "size[1]=" << ioSize[1] << " size[0]=" << ioSize[0] << endl;
        pio.readArray(&coordinates(0, memsize), bName, "x", 2, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 2, ioOffset, ioSize);
      } else {
        ioOffset[1] = 0;
        ioSize[1] = 0; // i
        ioOffset[0] = 0;
        ioSize[0] = 0; // j
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 2, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 2, ioOffset, ioSize);
      }
      for(int dim = 0; dim < nDim; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      currentWindows = getCurrentWindow(i, offset, size);

      for(MInt a = offset[1]; a < offset[1] + size[1]; ++a) {
        for(MInt b = offset[0]; b < offset[0] + size[0]; ++b) {
          nodeMap[count].blockId = i;
          nodeMap[count].pos[0] = a;
          nodeMap[count].pos[1] = b;
          std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1]) = getBCFromWindow(nodeMap[count].pos, currentWindows);
          count++;
        }
      }
      memsize += size[0] * size[1];

      //////////////////////////////
      ///////// FACE 3 /////////////
      //////////////////////////////
      // face +y  face 3
      if(m_domainId == 0) {
        ioOffset[1] = 1;
        ioSize[1] = m_grid->getBlockNoCells(i, 1) - 1; // i
        ioOffset[0] = 0;
        ioSize[0] = 1; // j
        pio.readArray(&coordinates(0, memsize), bName, "x", 2, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 2, ioOffset, ioSize);
      } else {
        ioOffset[1] = 0;
        ioSize[1] = 0; // i
        ioOffset[0] = 0;
        ioSize[0] = 0; // j
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 2, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 2, ioOffset, ioSize);
      }
      for(int dim = 0; dim < nDim; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      currentWindows = getCurrentWindow(i, offset, size);

      for(MInt b = offset[0]; b < offset[0] + size[0]; ++b) {
        for(MInt a = offset[1]; a < offset[1] + size[1]; ++a) {
          nodeMap[count].blockId = i;
          nodeMap[count].pos[0] = a;
          nodeMap[count].pos[1] = b;
          std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1]) = getBCFromWindow(nodeMap[count].pos, currentWindows);
          count++;
        }
      }
      memsize += size[0] * size[1];

      //////////////////////////////
      ///////// FACE 4 /////////////
      //////////////////////////////
      // face -y  face 4
      if(m_domainId == 0) {
        ioOffset[1] = 1;
        ioSize[1] = m_grid->getBlockNoCells(i, 1) - 1; // i
        ioOffset[0] = m_grid->getBlockNoCells(i, 0);
        ioSize[0] = 1; // j
        pio.readArray(&coordinates(0, memsize), bName, "x", 2, ioOffset, ioSize);
        pio.readArray(&coordinates(1, memsize), bName, "y", 2, ioOffset, ioSize);
      } else {
        ioOffset[1] = 0;
        ioSize[1] = 0; // i
        ioOffset[0] = 0;
        ioSize[0] = 0; // j
        MFloat empty = 0;
        pio.readArray(&empty, bName, "x", 2, ioOffset, ioSize);
        pio.readArray(&empty, bName, "y", 2, ioOffset, ioSize);
      }
      for(int dim = 0; dim < nDim; ++dim) {
        offset[dim] = ioOffset[dim];
        size[dim] = ioSize[dim];
      }

      currentWindows = getCurrentWindow(i, offset, size);

      for(MInt b = offset[0]; b < offset[0] + size[0]; ++b) {
        for(MInt a = offset[1]; a < offset[1] + size[1]; ++a) {
          nodeMap[count].blockId = i;
          nodeMap[count].pos[0] = a;
          nodeMap[count].pos[1] = b;
          std::tie(nodeMap[count].BC[0], nodeMap[count].BC[1]) = getBCFromWindow(nodeMap[count].pos, currentWindows);
          count++;
        }
      }
      memsize += size[0] * size[1];
    }

    // Sanity check
    for(MInt windowId = 0; windowId < noInputWindowInformation; ++windowId) {
      if(inputWindows[windowId]->BC < 6000 || inputWindows[windowId]->BC >= 6010) continue;
      MInt noWindowPoints = 1;
      for(MInt dim = 0; dim < nDim; ++dim) {
        noWindowPoints *= inputWindows[windowId]->endindex[dim] - inputWindows[windowId]->startindex[dim];
      }
      // Every point should be found at least once
      if(noWindowPoints > pointsFoundPerWindow[windowId]) mTerm(1, "noWindowPOints > pointFoundPerWindow");
    }

    // now broadcast the surface coordinates to every processor
    MPI_Bcast(&count, 1, MPI_INT, 0, m_StructuredComm, AT_, "count");
    MPI_Bcast(&memsize, 1, MPI_INT, 0, m_StructuredComm, AT_, "memsize");
    MPI_Bcast(&coordinates(0, 0), memsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "coordinates(0");
    MPI_Bcast(&coordinates(1, 0), memsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "coordinates(1");

    // also broadcast the nodeMap to every processor
    MPI_Datatype matrix;
    MPI_Type_contiguous(2 + 2 * nDim, MPI_INT, &matrix, AT_);
    MPI_Type_commit(&matrix, AT_);
    MPI_Bcast(nodeMap, memsize, matrix, 0, m_StructuredComm, AT_, "nodeMap");

    // Periodic
    pointType<nDim>* nodeMapP;
    MFloat tmppoint[nDim];
    MInt pcount = 0, plocation = 0, numofpoints[nDim];

    for(MInt i = 0; i < noInputWindowInformation; i++) {
      if(inputWindows[i]->BC > 4000 && inputWindows[i]->BC < 5000) {
        for(MInt j = 0; j < nDim; j++) {
          numofpoints[j] = inputWindows[i]->endindex[j] - inputWindows[i]->startindex[j] + 1;
        }
        pcount += numofpoints[0] * numofpoints[1];
      }
    }


    MFloatScratchSpace periodicCoordinates(nDim, max(1, pcount), AT_, "periodicCoordinates");
    periodicCoordinates.fill(-1.01010101);
    nodeMapP = new pointType<nDim>[pcount];

    if(domainId() == 0) {
      cout << "Loading periodic face coordinates" << endl;
    }
    MInt pmemsize = 0;
    for(MInt i = 0; i < noInputWindowInformation; i++) {
      if(inputWindows[i]->BC >= 4000 && inputWindows[i]->BC < 5000) {
        MString bName = "/block";
        stringstream number;
        number << inputWindows[i]->blockId << "/";
        bName += number.str();

        ParallelIo::size_type ioOffset[2] = {0, 0};
        ParallelIo::size_type ioSize[2] = {0, 0};
        if(domainId() == 0) {
          for(MInt j = 0; j < nDim; j++) {
            ioOffset[j] = inputWindows[i]->startindex[1 - j];
            ioSize[j] = inputWindows[i]->endindex[1 - j] - inputWindows[i]->startindex[1 - j] + 1;
          }

          pio.readArray(&periodicCoordinates(0, plocation), bName, "x", nDim, ioOffset, ioSize);
          pio.readArray(&periodicCoordinates(1, plocation), bName, "y", nDim, ioOffset, ioSize);
        } else {
          ioOffset[1] = 0;
          ioSize[1] = 0; // i
          ioOffset[0] = 0;
          ioSize[0] = 0; // j
          MFloat empty = 0;
          pio.readArray(&empty, bName, "x", nDim, ioOffset, ioSize);
          pio.readArray(&empty, bName, "y", nDim, ioOffset, ioSize);
        }

        for(int dim = 0; dim < nDim; ++dim) {
          offset[dim] = ioOffset[dim];
          size[dim] = ioSize[dim];
        }

        for(MInt b = offset[0]; b < offset[0] + size[0]; ++b) {
          for(MInt a = offset[1]; a < offset[1] + size[1]; ++a) {
            nodeMapP[plocation].BC[0] = inputWindows[i]->BC;
            nodeMapP[plocation].blockId = inputWindows[i]->blockId;
            nodeMapP[plocation].pos[0] = a;
            nodeMapP[plocation].pos[1] = b;
            plocation++;
          }
        }
        pmemsize += size[0] * size[1];
      }
    }

    MPI_Bcast(&plocation, 1, MPI_INT, 0, m_StructuredComm, AT_, "plocation");
    MPI_Bcast(&pmemsize, 1, MPI_INT, 0, m_StructuredComm, AT_, "pmemsize");

    MPI_Bcast(&periodicCoordinates(0, 0), pmemsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "periodicCoordinates(0");
    MPI_Bcast(&periodicCoordinates(1, 0), pmemsize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "periodicCoordinates(1");

    MPI_Bcast(nodeMapP, pmemsize, matrix, 0, m_StructuredComm, AT_, "nodeMapP");

    /////////////////////////////////
    /////// CONNECTIONS /////////////
    /////////////////////////////////
    if(domainId() == 0) {
      cout << "Building up connections for multiblock connection search" << endl;
    }
    vector<Point<2>> pts;
    for(MInt j = 0; j < count; ++j) {
      Point<2> a(coordinates(0, j), coordinates(1, j));
      pts.push_back(a);
      nodeMap[j].found = false;
    }

    MFloat m_gridEps = 0.0000001;

    KDtree<2> tree(pts);
    MBool add;
    MInt numConnections = 0, numSingularConnections = 0, nfound, tempnum = 0, tempcount;
    MInt results[10];
    // results=new MInt [10];
    for(MInt i = 0; i < count; ++i) {
      if(!nodeMap[i].found) {
        Point<2> a(coordinates(0, i), coordinates(1, i));
        nfound = tree.locatenear(a, m_gridEps, results, 10, false); // 0.00001
        for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
          for(MInt countNode3 = countNode2 + 1; countNode3 < nfound; ++countNode3) {
            MBool check = false;
            for(MInt bc2 = 0; bc2 < nDim; ++bc2) {
              if(nodeMap[results[countNode2]].BC[bc2] == 0) break;
              for(MInt bc3 = 0; bc3 < nDim; ++bc3) {
                if(nodeMap[results[countNode2]].BC[bc2] == nodeMap[results[countNode3]].BC[bc3]) {
                  check = true;
                  // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
                  if(addConnection(
                         nodeMap[results[countNode2]].BC[bc2], // Does this make sense? Check! //6000, //6000er-change
                         nodeMap[results[countNode2]].blockId,
                         nodeMap[results[countNode2]].pos,
                         nodeMap[results[countNode3]].blockId,
                         nodeMap[results[countNode3]].pos)) {
                    numConnections = numConnections + 1;
                  }
                }
              }
            }
            if(!check) {
              cout << "DANGER: x|y=" << coordinates(0, i) << "|" << coordinates(1, i)
                   << " nodeMap2: blockId=" << nodeMap[results[countNode2]].blockId << " BC =";
              for(MInt d = 0; d < nDim; ++d)
                cout << " " << nodeMap[results[countNode2]].BC[d];
              cout << " nodeMap3: blockId=" << nodeMap[results[countNode3]].blockId << " BC =";
              for(MInt d = 0; d < nDim; ++d)
                cout << " " << nodeMap[results[countNode3]].BC[d];
              cout << endl;
              //              mTerm(1, "");
            }

            /*            if (nodeMap[results[countNode2]].BC!=nodeMap[results[countNode3]].BC) {
                          // By now the higher BC has a higher priority; Check if this makes sense
                          nodeMap[results[countNode2]].BC = mMax(nodeMap[results[countNode2]].BC,
                                                                   nodeMap[results[countNode3]].BC);
                          cout << "WARNING: 2 nodes with different BCs!!!" << endl;
                        }
                        //CHANGE_SET: currently addConnection will always return true, because we are now using a
               multiset if (addConnection( nodeMap[results[countNode2]].BC, // Does this make sense? Check! //6000,
               //6000er-change nodeMap[results[countNode2]].blockId, nodeMap[results[countNode2]].pos,
                                           nodeMap[results[countNode3]].blockId,
                                           nodeMap[results[countNode3]].pos)) {
                          numConnections = numConnections + 1;
                        }*/
          }
          nodeMap[results[countNode2]].found = true;
        }

        // TODO_SS labels:FV,cleanup The following algorithm is better for singularity detection then the one below and
        // should replace
        //       the below algorithm; it works in case the singularity is surrounded by 6000er BCs
        MInt Nstar = 1;
        std::set<MInt> BCs;
        for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
          for(MInt bc2 = 0; bc2 < nDim; ++bc2) {
            const MInt BC = nodeMap[results[countNode2]].BC[bc2];
            if(BC < 6000 || BC >= 6010) Nstar = 0;
          }
        }
        if(Nstar) {
          tempcount = 0;
          for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            for(MInt j = 0; j < m_noBlocks; ++j) {
              if(nodeMap[results[countNode2]].blockId == j) {
                tempnum = j;
                break;
              }
            }
            for(MInt j = 0; j < nDim; ++j) {
              // pos[]: ijk  DirLast[]: kji
              if(nodeMap[results[countNode2]].pos[j] == 0
                 || nodeMap[results[countNode2]].pos[j] == m_grid->getBlockNoCells(tempnum, 1 - j)) {
                ++tempcount;
              }
            }
          }
          Nstar = nfound + 2 * nfound - tempcount;
          if(Nstar == 4) Nstar = 0;
          if(Nstar) {
            for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
              for(MInt d = 0; d < nDim; ++d) {
                if(nodeMap[results[countNode2]].BC[d] != 0) {
                  BCs.insert(nodeMap[results[countNode2]].BC[d]);
                }
              }
            }
          }
        }

        // if three points share the same
        // coordinate it must be a 3-star
        if(nfound == 3) {
          add = false;
          tempcount = 0;

          // check if it is a singularity 3-star or normal 3-point connection
          for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            for(MInt j = 0; j < m_noBlocks; ++j) {
              if(nodeMap[results[countNode2]].blockId == j) {
                tempnum = j;
                break;
              }
            }
            for(MInt j = 0; j < nDim; ++j) {
              // pos[]: ijk  DirLast[]: kji
              if(nodeMap[results[countNode2]].pos[j] == 0
                 || nodeMap[results[countNode2]].pos[j] == m_grid->getBlockNoCells(tempnum, 1 - j)) {
                ++tempcount;
              }
            }
          }

          if(tempcount == 6) { //(tempcount==9||tempcount==6) {
            add = true;
          }

          if(add) {
            for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
              if(BCs.size() > 1) cout << "CAUTION: This singular point has more than one BC!" << endl;
              for(const auto& BC : BCs) {
                // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
                if(addConnection(BC, // nodeMap[results[countNode2]].BC[0], //6000, //6000er-change
                                 nodeMap[results[countNode2]].blockId,
                                 nodeMap[results[countNode2]].pos,
                                 nodeMap[results[countNode2]].blockId,
                                 nodeMap[results[countNode2]].pos,
                                 3)) {
                  ASSERT(Nstar == 3, "");
                  Nstar = 0;
                  numSingularConnections = numSingularConnections + 1;
                }
              }
            }
          }
        }

        // if three points share the same
        // coordinate it must be a 5-star
        if(nfound == 5) {
          for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            if(BCs.size() > 1) cout << "CAUTION: This singular point has more than one BC!" << endl;
            for(const auto& BC : BCs) {
              // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
              if(addConnection(BC, // nodeMap[results[countNode2]].BC[0], //6000, //6000er-change
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               5)) {
                ASSERT(Nstar == 5, "");
                Nstar = 0;
                numSingularConnections = numSingularConnections + 1;
              }
            }
          }
        }
        // coordinate it must be a 6-star
        if(nfound == 6) {
          for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            if(BCs.size() > 1) cout << "CAUTION: This singular point has more than one BC!" << endl;
            for(const auto& BC : BCs) {
              // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
              if(addConnection(BC, // nodeMap[results[countNode2]].BC[0], //6000, //6000er-change
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               6)) {
                ASSERT(Nstar == 6, "");
                Nstar = 0;
                numSingularConnections = numSingularConnections + 1;
              }
            }
          }
        }
        if(Nstar != 0) {
          cout << "ERROR: nfound=" << nfound << " Nstar=" << Nstar << endl;
          //          mTerm(1, "Inconsistency in singularity detection!");
        }
      }
    }

    if(domainId() == 0) {
      cout << "Computing periodic window displacements" << endl;
    }
    m_log << "Computing periodic window displacements" << endl;
    /////////////////////////////////////////////////
    //////// PERIODIC WINDOW DISPLACEMENTS //////////
    /////////////////////////////////////////////////
    MFloatScratchSpace periodicWindowCenter(m_noBlocks, nDim, 2 * nDim, AT_, "periodicWindowCOGS");
    MIntScratchSpace periodicWindowNoNodes(m_noBlocks, 2 * nDim, AT_, "periodicWindowNoNodes");
    cout.precision(10);
    periodicWindowCenter.fill(F0);
    periodicWindowNoNodes.fill(0);

    // compute displacement
    const MInt periodicOffset = 4401;
    for(MInt i = 0; i < pcount; ++i) {
      if(nodeMapP[i].BC[0] < 4401 || nodeMapP[i].BC[0] > 4406) {
        continue;
      }

      // add up all coordinates
      for(MInt dim = 0; dim < nDim; dim++) {
        periodicWindowCenter(nodeMapP[i].blockId, dim, nodeMapP[i].BC[0] - periodicOffset) +=
            periodicCoordinates(dim, i);
      }

      periodicWindowNoNodes(nodeMapP[i].blockId, nodeMapP[i].BC[0] - periodicOffset) += 1;
    }

    // new approach to account also for input block distribution which are periodic in spanwise direction
    // and splited in spanwise direction

    // 1) Compute center of the the faces and split even and odd sides
    multimap<MInt, pair<vector<MFloat>, MInt>> oddPeriodicSides;
    multimap<MInt, pair<vector<MFloat>, MInt>> evenPeriodicSides;
    // vector< tuple<MInt, vector<MInt>, MInt>> oddPeriodicSides;
    // vector< tuple<MInt, vector<MInt>, MInt>> evenPeriodicSides;
    // evenPeriodicSides.reserve(m_noBlocks*nDim); //per block we can only have max 3 odd sides
    // oddPeriodicSides.reserve(m_noBlocks*nDim); //per block we can only have mac 3 even sides

    // now compute the center
    //(not weighted, should be basically the same surface in another location)
    // and compute distance between the two centers
    for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
      for(MInt periodicWindowId = 0; periodicWindowId < 2 * nDim; periodicWindowId++) {
        if(periodicWindowNoNodes(blockId, periodicWindowId) <= 0) { // we have no periodic cells
          continue;
        }
        vector<MFloat> centerDimensions(nDim);
        for(MInt dim = 0; dim < nDim; dim++) { // add the periodic cells
          periodicWindowCenter(blockId, dim, periodicWindowId) /=
              (MFloat)periodicWindowNoNodes(blockId, periodicWindowId);
          centerDimensions[dim] =
              periodicWindowCenter(blockId, dim,
                                   periodicWindowId); //(MFloat)periodicWindowNoNodes(blockId,periodicWindowId);
          cout << "periodic Window center " << periodicWindowCenter(blockId, dim, periodicWindowId) << endl;
        }
        if(periodicWindowId % 2 == 0) { // we have an even side
          // evenPeriodicSides.push_back(make_tuple(blockId,centerDimensions, periodicWindowId));
          evenPeriodicSides.insert(pair<MInt, pair<vector<MFloat>, MInt>>(
              blockId, pair<vector<MFloat>, MInt>(centerDimensions, periodicWindowId)));
          cout << "dimensions of centerDimensions (even)" << centerDimensions[0] << " " << centerDimensions[1]
               << endl; // << " " << centerDimensions[2] << endl;
        } else {        // we have an odd side
          // oddPeriodicSides.push_back(make_tuple(blockId,centerDimensions, periodicWindowId));
          oddPeriodicSides.insert(pair<MInt, pair<vector<MFloat>, MInt>>(
              blockId, pair<vector<MFloat>, MInt>(centerDimensions, periodicWindowId)));
          cout << "dimensions of centerDimensions (odd)" << centerDimensions[0] << " " << centerDimensions[1]
               << endl; //" " << centerDimensions[2] << endl;
        }
      }
    }

    // 2) now set the perodic sides for the input blockId
    // a) check for periodic sides in the same blockId

    // get the iterator start and end of the m_block
    pair<multimap<MInt, pair<vector<MFloat>, MInt>>::iterator, multimap<MInt, pair<vector<MFloat>, MInt>>::iterator>
        evenIt;
    pair<multimap<MInt, pair<vector<MFloat>, MInt>>::iterator, multimap<MInt, pair<vector<MFloat>, MInt>>::iterator>
        oddIt;


    MFloatScratchSpace periodicDisplacementsBlockIds(m_noBlocks, nDim * nDim, AT_, "periodic Displacements per block");

    for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
      MInt noEvenPerodicWindowsBlockIds = evenPeriodicSides.count(blockId);
      MInt noOddPeriodicWindowsBlockIds = oddPeriodicSides.count(blockId);
      cout << "we have evenSides " << noEvenPerodicWindowsBlockIds << " and we have oddSides "
           << noOddPeriodicWindowsBlockIds << " in BlockId " << blockId << endl;
      if(noEvenPerodicWindowsBlockIds == 0) continue; // there are no periodic connections
      if(noOddPeriodicWindowsBlockIds == 0) continue; // there are no periodic connections associated with that blockId
      // else we have a interBlock connection (much easier to find
      evenIt = evenPeriodicSides.equal_range(blockId);
      multimap<MInt, pair<vector<MFloat>, MInt>>::iterator eit = evenIt.first;
      while(eit != evenIt.second) {
        MBool tobeErased = false;
        MInt evenSideDim =
            ceil(((MFloat)((*eit).second.second + 1) / 2.0) - 1); // get the dimensional direction of the side
        MInt evenSide = (*eit).second.second;
        MInt oddSide = evenSide + 1;
        oddIt = oddPeriodicSides.equal_range(blockId);
        multimap<MInt, pair<vector<MFloat>, MInt>>::iterator oit = oddIt.first;

        while(oit != oddIt.second) {
          if((*oit).second.second == oddSide) { // yes we found a corresponding odd side
            // compute the periodic connection and put them into the scratchspace
            for(MInt dim = 0; dim < nDim; dim++) {
              cout << "first element is " << (*oit).second.first[dim] << " and we subtract " << (*eit).second.first[dim]
                   << endl;
              periodicDisplacementsBlockIds(blockId, evenSideDim + nDim * dim) =
                  abs((*oit).second.first[dim] - (*eit).second.first[dim]);
            }
            oit = oddPeriodicSides.erase(oit);
            tobeErased = true;
            break;
          } else {
            oit++;
          }
        }
        if(tobeErased) {
          eit = evenPeriodicSides.erase(eit);
        } else {
          eit++;
        }
      }
    }

    MFloat periodicEpsilon = 0.000001;
    // now the list only contains the inner blockId periodic conditions.
    // check for multiblock communications
    if(evenPeriodicSides.size() != 0
       && oddPeriodicSides.size() != 0) { // we did not delete everything in the side container
      for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
        MInt noEvenPerodicWindowsBlockId = evenPeriodicSides.count(blockId);
        if(noEvenPerodicWindowsBlockId == 0) continue; // there are no periodic connections left in the even container
        evenIt = evenPeriodicSides.equal_range(blockId);
        multimap<MInt, pair<vector<MFloat>, MInt>>::iterator eit = evenIt.first;
        while(eit != evenIt.second) {
          MBool tobeErased = false;
          MInt evenSideDim =
              ceil(((MFloat)((*eit).second.second + 1) / 2.0) - 1); // get the dimensional direction of the side
          // MInt evenSide    = (*eit).second.second;
          multimap<MInt, pair<vector<MFloat>, MInt>>::iterator oit = oddPeriodicSides.begin();
          while(oit != oddPeriodicSides.end()) {
            // check the number of conncetions which have the same cooridnates
            MInt equalConnectionDims = 0;
            for(MInt d = 0; d < nDim; d++) {
              if(abs((*oit).second.first[d] - (*eit).second.first[d]) < periodicEpsilon) equalConnectionDims++;
            }
            if(equalConnectionDims == 2) { // yes we found a pair belonging together
              for(MInt dim = 0; dim < nDim; dim++) {
                periodicDisplacementsBlockIds(blockId, evenSideDim + nDim * dim) =
                    abs((*oit).second.first[dim] - (*eit).second.first[dim]);
                // now store the odd side also for multiblock communication
                periodicDisplacementsBlockIds((*oit).first, evenSideDim + nDim * dim) =
                    periodicDisplacementsBlockIds(blockId, evenSideDim + nDim * dim);
              }
              tobeErased = true;
              oit = oddPeriodicSides.erase(oit);
            } else {
              oit++;
            }
          }
          if(tobeErased) {
            eit = evenPeriodicSides.erase(eit);
          } else {
            eit++;
          }
        }
      }
    }


    for(MInt periodicWindowId = 0; periodicWindowId < 2 * nDim; periodicWindowId++) {
      if(periodicWindowId % 2 == 1) {
        for(MInt dim = 0; dim < nDim; dim++) {
          const MInt displacementId = (MFloat)(periodicWindowId + 1) / 2.0 - 1;
          periodicDisplacements[dim * nDim + displacementId] =
              periodicWindowCenter(m_blockId, dim, periodicWindowId)
              - periodicWindowCenter(m_blockId, dim, periodicWindowId - 1);
          m_log << m_blockId << " periodicWindowCenter for window: " << periodicWindowId + periodicOffset
                << "  dim: " << dim << " displacementId: " << displacementId
                << " displacement: " << periodicDisplacements[dim * nDim + displacementId] << endl;
        }
      }
    }

    cout << "old approach" << endl;
    for(MInt periodicWindowId = 0; periodicWindowId < 2 * nDim; periodicWindowId++) {
      if(periodicWindowId % 2 == 1) {
        for(MInt dim = 0; dim < nDim; dim++) {
          const MInt displacementId = (MFloat)(periodicWindowId + 1) / 2.0 - 1;
          periodicDisplacements[dim * nDim + displacementId] =
              periodicWindowCenter(m_blockId, dim, periodicWindowId)
              - periodicWindowCenter(m_blockId, dim, periodicWindowId - 1);
          cout << "blockId = " << m_blockId << "periodicWindowCenter for window: " << periodicWindowId + periodicOffset
               << "  dim: " << dim << " displacementId: " << displacementId
               << " displacement: " << periodicDisplacements[dim * nDim + displacementId] << endl;
        }
      }
    }
    cout << "new approach " << endl;
    for(MInt periodicWindowId = 0; periodicWindowId < 2 * nDim; periodicWindowId++) {
      if(periodicWindowId % 2 == 1) {
        const MInt displacementId = (MFloat)(periodicWindowId + 1) / 2.0 - 1;
        for(MInt dim = 0; dim < nDim; dim++) {
          cout << "blockId = " << m_blockId << "periodicWindowCenter for window: " << periodicWindowId + periodicOffset
               << "  dim: " << dim << " displacementId: " << displacementId
               << " displacement: " << periodicDisplacementsBlockIds(m_blockId, dim * nDim + displacementId) << endl;
        }
      }
    }
    MPI_Barrier(m_StructuredComm, AT_);
    if(domainId() == 0) {
      for(MInt b = 0; b < m_noBlocks; b++) {
        cout << "BlockId " << b << ": " << endl;
        cout << "\t";
        for(MInt i = 0; i < nDim * nDim; i++) {
          cout << periodicDisplacementsBlockIds(b, i) << " ";
        }
        cout << endl;
      }
    }
    MPI_Barrier(m_StructuredComm, AT_);
    for(MInt b = 0; b < m_noBlocks; b++) {
      if(b == m_blockId) {
        for(MInt i = 0; i < nDim * nDim; i++) {
          periodicDisplacements[i] = periodicDisplacementsBlockIds(b, i);
        }
      }
    }
    // mTerm(-1, AT_, "I want to terminate ");

    /////////////////////////////////////////////////
    //////// PERIODIC CONNECTIONS ///////////////////

    /////////////////////////////////////////////////

    for(MInt i = 0; i < pcount; ++i) {
      tmppoint[0] = periodicCoordinates(0, i);
      tmppoint[1] = periodicCoordinates(1, i);
      periodicPointsChange(tmppoint, nodeMapP[i].BC[0], periodicDisplacements);
      Point<2> aaa(tmppoint[0], tmppoint[1]);
      nfound = tree.locatenear(aaa, m_gridEps, results, 10, false); // 0.0001
      for(MInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
        // CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
        if(addConnection(nodeMapP[i].BC[0],
                         nodeMapP[i].blockId,
                         nodeMapP[i].pos,
                         nodeMap[results[countNode2]].blockId,
                         nodeMap[results[countNode2]].pos)) {
          numConnections = numConnections + 1;
        }
      }

      /*
      if(nfound==5) {
        //CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
        if (addConnection( nodeMapP[i].BC[0],
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           5)) {
          numSingularConnections = numSingularConnections + 1;
        }
      }
      //6 star for periodic bc
      if(nfound==6) {
        //CHANGE_SET: currently addConnection will always return true, because we are now using a multiset
        if (addConnection( nodeMapP[i].BC[0],
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           6)) {
          numSingularConnections = numSingularConnections + 1;
        }
      }*/
    }


    if(domainId() == 0) {
      cout << "Assemble multiblock connections" << endl;
    }
    ///////////////////////////////////////////////////
    /////// CREATE THE WINDOWS FROM CONNECTIONS ///////
    ///////////////////////////////////////////////////
    multiBlockAssembling();
    singularityAssembling();

    // Delete duplicate windows
    deleteDuplicateWindows(window0d, window0d);
    deleteDuplicateWindows(window1d, window1d);
    deleteDuplicateWindows(window0d, window1d);

    for(auto it = singularwindow.cbegin(); it != singularwindow.cend(); ++it) {
      // TODO_SS labels:FV by now don't do anything if BC is perioduc
      if((*it)->BC >= 4400 && (*it)->BC < 4410) continue;
      MInt noSameWindows = 0;
      (*it)->BCsingular[noSameWindows++] = (*it)->BC;
      (*it)->BC = 0;
      for(auto it2 = it + 1; it2 != singularwindow.cend();) {
        const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
        const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
        const MBool test = mapCompare11(map1, map2);
        if(test) {
          if((*it)->BCsingular[0] == (*it2)->BC || (*it)->Nstar != (*it2)->Nstar) {
            mTerm(1, "There is no hope anymore!");
          }

          if(noSameWindows >= (*it)->Nstar) {
            mTerm(1, "");
          }

          (*it)->BCsingular[noSameWindows++] = (*it2)->BC;
          singularwindow.erase(it2);
        } else {
          ++it2;
        }
      }
    }


    // Create global BC maps for those singularities which don't have one
    /*    for (MInt i = 0; i < (signed)singularwindow.size(); ++i) {
          const auto& mymap = singularwindow[i];
          for (MInt bc = 0; bc < mymap->Nstar; ++bc) {
            if (mymap->BCsingular[bc]!=0 && mymap->BCsingular[bc]!=6000) {
              MBool addMap = true;
              for (MInt j = 0; j < (signed)globalStructuredBndryCndMaps.size(); ++j) {
                if (mymap->Id1==globalStructuredBndryCndMaps[j]->Id1 &&
       mymap->BCsingular[bc]==globalStructuredBndryCndMaps[j]->BC) { MBool addMap_temp = false; for (MInt d = 0; d <
       nDim; ++d) { if (globalStructuredBndryCndMaps[j]->start1[d]>mymap->start1[d] ||
       globalStructuredBndryCndMaps[j]->end1[d]<mymap->end1[d]) { addMap_temp = true;
                    }
                  }
                  addMap = addMap_temp;
                  if (!addMap) break;
                }
              }
              if (addMap) {
                StructuredWindowMap* windowMap=new StructuredWindowMap(nDim);
                MInt order[3]={0,1,2};
                MInt step1[3]={1,1,1};
                mapCreate(mymap->Id1, mymap->start1, mymap->end1, step1,
                          -1 ,  NULL, NULL , NULL, order,mymap->BCsingular[bc], windowMap );
                windowMap->Nstar = mymap->Nstar;
                //add the map to the list!
                if (domainId()==0) {
                  cout << "ADDDING following GlobalBCMap:" << endl;
                  mapPrint(windowMap);
                }
                globalStructuredBndryCndMaps.push_back(std::move(windowMap));
              }
            }
          }
        }*/

    unique_ptr<StructuredWindowMap<nDim>> windowMap;
    // TODO_SS labels:FV in 2D impossible
    /*    for(MInt i=0;i<(MInt)window2d.size();++i) {

          // To distinguish between communication bcs and normal ones, do the following, negate the communication bcs
          const MInt BC = (window2d[i]->BC>=6000 && window2d[i]->BC<6010) ? -window2d[i]->BC : window2d[i]->BC;

          windowMap=new StructuredWindowMap(nDim);
          mapCreate(window2d[i]->Id1, window2d[i]->start1, window2d[i]->end1, window2d[i]->step1,
                    window2d[i]->Id2, window2d[i]->start2, window2d[i]->end2, window2d[i]->step2,
                    window2d[i]->order, BC, windowMap);
          windowMap->Nstar=-1;windowMap->SingularId=-1;
          windowMap->dc1=window2d[i]->dc1;
          windowMap->dc2=window2d[i]->dc2;

          //i-surface: change sign of I; same for J and K
          if(windowMap->dc1*windowMap->dc2>0) {
            MInt tempdim=abs(windowMap->dc2)-1;
            windowMap->step2[tempdim]=-1;
          }

          //add the map to the list!
          globalStructuredBndryCndMaps.push_back(std::move(windowMap));


          if (window2d[i]->BC>=6000 && window2d[i]->BC<6010) { //( window2d[i]->BC==6000) { //6000er-change
            windowMap=new StructuredWindowMap(nDim);
            mapCreate(window2d[i]->Id1, window2d[i]->start1, window2d[i]->end1, window2d[i]->step1,
                      window2d[i]->Id2, window2d[i]->start2, window2d[i]->end2 ,window2d[i]->step2,
                      window2d[i]->order, BC, windowMap);
            windowMap->Nstar=-1;windowMap->SingularId=-1;
            windowMap->dc1=window2d[i]->dc1;
            windowMap->dc2=window2d[i]->dc2;
            //i-surface: change sign of I; same for J and K
            if(windowMap->dc1*windowMap->dc2>0) {
              MInt tempdim=abs(windowMap->dc2)-1;
              windowMap->step2[tempdim]=-1;
            }

            mapInvert1(windowMap);
            mapNormalize3(windowMap);
            globalStructuredBndryCndMaps.push_back(std::move(windowMap));
          }
        }*/


    for(MInt i = 0; i < (MInt)window1d.size(); ++i) {
      // To distinguish between communication bcs and normal ones, do the following, negate the communication bcs
      const MInt BC = (window1d[i]->BC >= 6000 && window1d[i]->BC < 6010) ? -window1d[i]->BC : window1d[i]->BC;

      windowMap = make_unique<StructuredWindowMap<nDim>>();
      mapCreate(window1d[i]->Id1, window1d[i]->start1, window1d[i]->end1, window1d[i]->step1, window1d[i]->Id2,
                window1d[i]->start2, window1d[i]->end2, window1d[i]->step2, window1d[i]->order, BC, windowMap);
      windowMap->Nstar = -1;
      windowMap->SingularId = -1;
      windowMap->dc1 = window1d[i]->dc1;
      windowMap->dc2 = window1d[i]->dc2;

      // i-surface: change sign of I; same for J and K
      for(MInt j = 0; j < nDim; ++j) {
        if(windowMap->start1[j] == windowMap->end1[j]) {
          if(windowMap->start1[j] == 0 && windowMap->start2[windowMap->order[j]] == 0) {
            windowMap->step2[windowMap->order[j]] = -1;
          }
          if(windowMap->start1[j] > 0 && windowMap->start2[windowMap->order[j]] > 0) {
            windowMap->step2[windowMap->order[j]] = -1;
          }
        }
      }

      // now for 5 star communicaitons (or ...... if needed)
      // 3 star does not need additional information
      for(MInt j = 0; j < (MInt)singularwindow.size(); ++j) {
        const MBool test = mapCompare11(windowMap, singularwindow[j]);
        if(test) {
          mTerm(1, "Impossible behaviour!");
          windowMap->Nstar = singularwindow[j]->Nstar;
          windowMap->SingularId = singularwindow[j]->SingularId;
        }
      }

      // add the map to the list!
      globalStructuredBndryCndMaps.push_back(std::move(windowMap));

      if(window1d[i]->BC >= 6000 && window1d[i]->BC < 6010) { //( window1d[i]->BC==6000) { //6000er-change
        windowMap = make_unique<StructuredWindowMap<nDim>>();
        mapCreate(window1d[i]->Id1, window1d[i]->start1, window1d[i]->end1, window1d[i]->step1, window1d[i]->Id2,
                  window1d[i]->start2, window1d[i]->end2, window1d[i]->step2, window1d[i]->order, BC, windowMap);
        windowMap->Nstar = -1;
        windowMap->SingularId = -1;
        windowMap->dc1 = window1d[i]->dc1;
        windowMap->dc2 = window1d[i]->dc2;

        // i-surface: change sign of I; same for J and K
        for(MInt j = 0; j < nDim; ++j) {
          if(windowMap->start1[j] == windowMap->end1[j]) {
            if(windowMap->start1[j] == 0 && windowMap->start2[windowMap->order[j]] == 0) {
              windowMap->step2[windowMap->order[j]] = -1;
            }
            if(windowMap->start1[j] > 0 && windowMap->start2[windowMap->order[j]] > 0) {
              windowMap->step2[windowMap->order[j]] = -1;
            }
          }
        }

        mapInvert1(windowMap);
        mapNormalize3(windowMap);

        // now for 5 star communicaitons (or ...... if needed)
        // 3 star does not need additional information
        for(MInt j = 0; j < (MInt)singularwindow.size(); ++j) {
          const MBool test = mapCompare11(windowMap, singularwindow[j]);
          if(test) {
            mTerm(1, "Impossible behaviour!");
            windowMap->Nstar = singularwindow[j]->Nstar;
            windowMap->SingularId = singularwindow[j]->SingularId;
          }
        }

        globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      }
    }


    for(MInt i = 0; i < (MInt)window0d.size(); ++i) {
      // To distinguish between communication bcs and normal ones, do the following, negate the communication bcs
      const MInt BC = (window0d[i]->BC >= 6000 && window0d[i]->BC < 6010) ? -window0d[i]->BC : window0d[i]->BC;

      // point communicaiton (rare)
      windowMap = make_unique<StructuredWindowMap<nDim>>();
      mapCreate(window0d[i]->Id1, window0d[i]->start1, window0d[i]->end1, window0d[i]->step1, window0d[i]->Id2,
                window0d[i]->start2, window0d[i]->end2, window0d[i]->step2, window0d[i]->order, BC, windowMap);
      windowMap->Nstar = -1;
      windowMap->SingularId = -1;
      windowMap->dc1 = window0d[i]->dc1;
      windowMap->dc2 = window0d[i]->dc2;
      // i-surface: change sign of I; same for J and K
      for(MInt j = 0; j < nDim; ++j) {
        if(windowMap->start1[j] == windowMap->end1[j]) {
          if(windowMap->start1[j] == 0 && windowMap->start2[windowMap->order[j]] == 0) {
            windowMap->step2[windowMap->order[j]] = -1;
          }
          if(windowMap->start1[j] > 0 && windowMap->start2[windowMap->order[j]] > 0) {
            windowMap->step2[windowMap->order[j]] = -1;
          }
        }
      }

      // now for 5 star communicaitons (or ...... if needed)
      // 3 star does not need additional information
      for(MInt j = 0; j < (MInt)singularwindow.size(); ++j) {
        const MBool test = mapCompare11(windowMap, singularwindow[j]);
        if(test) {
          windowMap->Nstar = singularwindow[j]->Nstar;
          windowMap->SingularId = singularwindow[j]->SingularId;
        }
      }

      // add the map to the list!
      globalStructuredBndryCndMaps.push_back(std::move(windowMap));

      if(window0d[i]->BC >= 6000 && window0d[i]->BC < 6010) { //( window0d[i]->BC==6000) { //6000er-change
        windowMap = make_unique<StructuredWindowMap<nDim>>();
        mapCreate(window0d[i]->Id1, window0d[i]->start1, window0d[i]->end1, window0d[i]->step1, window0d[i]->Id2,
                  window0d[i]->start2, window0d[i]->end2, window0d[i]->step2, window0d[i]->order, BC, windowMap);
        windowMap->Nstar = -1;
        windowMap->SingularId = -1;
        windowMap->dc1 = window0d[i]->dc1;
        windowMap->dc2 = window0d[i]->dc2;

        // i-surface: change sign of I; same for J and K
        for(MInt j = 0; j < nDim; ++j) {
          if(windowMap->start1[j] == windowMap->end1[j]) {
            if(windowMap->start1[j] == 0 && windowMap->start2[windowMap->order[j]] == 0) {
              windowMap->step2[windowMap->order[j]] = -1;
            }
            if(windowMap->start1[j] > 0 && windowMap->start2[windowMap->order[j]] > 0) {
              windowMap->step2[windowMap->order[j]] = -1;
            }
          }
        }

        mapInvert1(windowMap);
        mapNormalize3(windowMap);

        // now for 5 star communicaitons (or ...... if needed)
        // 3 star does not need additional information
        for(MInt j = 0; j < (MInt)singularwindow.size(); ++j) {
          const MBool test = mapCompare11(windowMap, singularwindow[j]);
          if(test) {
            windowMap->Nstar = singularwindow[j]->Nstar;
            windowMap->SingularId = singularwindow[j]->SingularId;
          }
        }
        globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      }
    }

    // Duplicates are already deleted above. The following call should be unnecessary
    deleteDuplicateCommMaps();

#ifndef NDEBUG
    if(domainId() == 0) {
      cout << "======= GLOBALSTRUCTUREDBNDRYCNDMAPS ======" << endl;
      for(MInt i = 0; i < (signed)globalStructuredBndryCndMaps.size(); ++i) {
        mapPrint(globalStructuredBndryCndMaps[i]);
      }

      cout << " ======= SINGULARWINDOW =============" << endl;
      for(MInt i = 0; i < (signed)singularwindow.size(); ++i) {
        mapPrint(singularwindow[i]);
      }

      cout << " ======== window0d ==============" << endl;
      for(MInt i = 0; i < (signed)window0d.size(); ++i) {
        mapPrint(window0d[i]);
      }

      cout << " ======== window1d ==============" << endl;
      for(MInt i = 0; i < (signed)window1d.size(); ++i) {
        mapPrint(window1d[i]);
      }

      cout << " ======== window2d ==============" << endl;
      for(MInt i = 0; i < (signed)window2d.size(); ++i) {
        mapPrint(window2d[i]);
      }
    }
#endif

    if(!hasConnectionInfo) {
      writeConnectionWindowInformation2D(periodicDisplacements);
    } else {
      if(domainId() == 0) {
        cout << "##########################################################################" << endl;
        cout << "WARNING: OLD CONNECTION INFO EXISTS, DELETE BEFORE NEW ONE CAN BE WRITTEN!" << endl;
        cout << "##########################################################################" << endl;
      }
    }

    if(domainId() == 0) {
      cout << "Connection identification and window creation finished!" << endl;
    }
  }
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::setBCsingular(unique_ptr<StructuredWindowMap<nDim>>& s_window,
                                                       unique_ptr<StructuredWindowMap<nDim>>& g_window,
                                                       const MInt singularcount) {
  const MInt s_blockId = s_window->Id1;
  const MInt g_blockId = g_window->Id2;
  stringstream prop_name;
  prop_name << "BCMap_" << std::min(s_blockId, g_blockId) << "_" << std::max(s_blockId, g_blockId);
  const MInt BC_temp = Context::getSolverProperty<MInt>(prop_name.str(), solverId(), AT_);
  if(BC_temp < 6000 || BC_temp >= 6010) mTerm(1, "Your are an idiot!");
  s_window->BCsingular[singularcount + 2] = -BC_temp;
  g_window->BC = -BC_temp;
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::deleteDuplicateWindows(
    std::vector<unique_ptr<StructuredWindowMap<nDim>>>& window1,
    std::vector<unique_ptr<StructuredWindowMap<nDim>>>& window2) {
  // Communication BCs exist in duplicate; the following algorithm should consistently remove both maps
  const MBool isSame = (&window1 == &window2);

  // Loop over all comm Bcs
  for(auto it = window1.cbegin(); it != window1.cend();) {
    if(window1.empty()) {
      break;
    }
    MBool deleted1 = false;
    if((*it)->BC >= 6000 && (*it)->BC < 6010) {
      for(auto it2 = isSame ? it + 1 : window2.cbegin(); it2 != window2.cend();) {
        if(window2.empty()) {
          break;
        }
        MBool deleted2 = false;
        if((*it2)->BC >= 6000 && (*it2)->BC < 6010) {
          // Compare if both BCs describe the communication bewteen same blocks
          // TODO_SS labels:FV,toenhance Isn't it enough to check just for the receive blocks, i.g., to check that
          //      a block receives data for same range of cells multiple times? Afterwards also delete
          //      the inverted maps
          if((*it)->Id1 == (*it2)->Id1 && (*it)->Id2 == (*it2)->Id2) {
            // Determine if one BC is fully included in the other one;
            // Assumption: When on e map is fully included in another map on one side
            //             then it is also included on the opposite side
            MInt diffStart[nDim], diffEnd[nDim];
            for(MInt dim = 0; dim < nDim; ++dim) {
              // Sanity check that map is normalized
              if((*it)->step1[dim] < 0 || (*it2)->step1[dim] < 0) {
                mTerm(1, "");
              }

              diffStart[dim] = (*it)->start1[dim] - (*it2)->start1[dim];
              diffEnd[dim] = (*it)->end1[dim] - (*it2)->end1[dim];
            }

            const MInt minStart = *std::min_element(&diffStart[0], &diffStart[0] + nDim);
            const MInt maxStart = *std::max_element(&diffStart[0], &diffStart[0] + nDim);
            const MInt minEnd = *std::min_element(&diffEnd[0], &diffEnd[0] + nDim);
            const MInt maxEnd = *std::max_element(&diffEnd[0], &diffEnd[0] + nDim);
            if(minStart * maxStart >= 0 && minEnd * maxEnd >= 0
               && ((minStart >= 0 && maxEnd <= 0) || (maxStart <= 0 && minEnd >= 0))) {
              //            if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && (minStart>=maxEnd||maxStart<=minEnd)) {
              //            if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && minStart*minEnd<=0) {
              // TODO_SS labels:FV Previously  if two maps cover exactly same range, nothing was done, except they
              // belong to same BC
              //              const MBool sameBC = (((*it))->BC == (it2)->BC);
              //              if (!(minStart==0 && maxStart==0 && minEnd==0 && maxEnd==0 && !sameBC)) {

              if(minStart == 0 && maxStart == 0 && minEnd == 0 && maxEnd == 0) {
                if((*it)->BC != 6000 || (*it2)->BC != 6000) {
                  if(domainId() == 0) cout << "CAUTION: BC reset to 6000!" << endl;
                  (*it2)->BC = 6000;
                }
              }

              // Check which map is the enclosed map
              //                if (maxStart>0 || minEnd<0) {
              if(minStart >= 0 && maxEnd <= 0) {
#ifndef NDEBUG
                if(domainId() == 0) {
                  cout << "!!! DELETE WINDOW: !!!" << endl;
                  const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
                  mapPrint(map1);
                  cout << "!!! BECAUSE OF: !!!" << endl;
                  const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
                  mapPrint(map2);
                  cout << "!!!!!!!!!!!!!!!!!!!" << endl;
                }
#endif
                it = window1.erase(it);
                deleted1 = true;
                break;
              } else {
#ifndef NDEBUG
                if(domainId() == 0) {
                  cout << "!!! DELETE WINDOW2: !!!" << endl;
                  const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
                  mapPrint(map2);
                  cout << "!!! BECAUSE OF: !!!" << endl;
                  const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
                  mapPrint(map1);
                  cout << "!!!!!!!!!!!!!!!!!!!" << endl;
                }
#endif
                it2 = window2.erase(it2);
                deleted2 = true;
              }

              //              }
            }
          }
        }
        if(!deleted2) ++it2;
      }
    }
    if(!deleted1) ++it;
  }


  // Similar story for periodic BCs: Check if receive maps overlap and delete the one which is
  // completly contained in some other periodic map
  for(auto it = window1.cbegin(); it != window1.cend();) {
    if(window1.empty()) {
      break;
    }
    MBool deleted1 = false;
    if((*it)->BC >= 4401 && (*it)->BC <= 4406) {
      for(auto it2 = isSame ? it + 1 : window2.cbegin(); it2 != window2.cend();) {
        if(window2.empty()) {
          break;
        }
        MBool deleted2 = false;
        if((*it2)->BC >= 4401 && (*it2)->BC <= 4406) {
          // Compare receive blocks
          if((*it)->Id2 == (*it2)->Id2) {
            // Determine if one BC is fully included in the other one;
            // Assumption: When on e map is fully included in another map on one side
            //             then it is also included on the opposite side
            MInt diffStart[nDim], diffEnd[nDim];
            for(MInt dim = 0; dim < nDim; ++dim) {
              auto it_start2 = (*it)->step2[dim] > 0 ? (*it)->start2[dim] : (*it)->end2[dim];
              auto it_end2 = (*it)->step2[dim] > 0 ? (*it)->end2[dim] : (*it)->start2[dim];
              auto it2_start2 = (*it2)->step2[dim] > 0 ? (*it2)->start2[dim] : (*it2)->end2[dim];
              auto it2_end2 = (*it2)->step2[dim] > 0 ? (*it2)->end2[dim] : (*it2)->start2[dim];
              diffStart[dim] = it_start2 - it2_start2;
              diffEnd[dim] = it_end2 - it2_end2;
            }


            const MInt minStart = *std::min_element(&diffStart[0], &diffStart[0] + nDim);
            const MInt maxStart = *std::max_element(&diffStart[0], &diffStart[0] + nDim);
            const MInt minEnd = *std::min_element(&diffEnd[0], &diffEnd[0] + nDim);
            const MInt maxEnd = *std::max_element(&diffEnd[0], &diffEnd[0] + nDim);
            if(minStart * maxStart >= 0 && minEnd * maxEnd >= 0
               && ((minStart >= 0 && maxEnd <= 0) || (maxStart <= 0 && minEnd >= 0))) {
              //            if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && (minStart>=maxEnd||maxStart<=minEnd)) {
              //            if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && minStart*minEnd<=0) {
              // TODO_SS labels:FV By now if two maps cover exactly same range, do nothing
              if(!(minStart == 0 && maxStart == 0 && minEnd == 0 && maxEnd == 0)) {
                // Check which map is the enclosed map
                //                if (maxStart>0 || minEnd<0) {
                if(minStart >= 0 && maxEnd <= 0) {
#ifndef NDEBUG
                  if(domainId() == 0) {
                    cout << "!!! DELETE WINDOW: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
                    mapPrint(map1);
                    cout << minStart << "|" << maxStart << " " << minEnd << "|" << maxEnd << endl;
                    cout << "!!! BECAUSE OF: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
                    mapPrint(map2);
                    cout << "!!!!!!!!!!!!!!!!!!!" << endl;
                  }
#endif
                  it = window1.erase(it);
                  deleted1 = true;
                  break;
                } else {
#ifndef NDEBUG
                  if(domainId() == 0) {
                    cout << "!!! DELETE WINDOW2: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
                    mapPrint(map2);
                    cout << "!!! BECAUSE OF: !!!" << endl;
                    cout << minStart << "|" << maxStart << " " << minEnd << "|" << maxEnd << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
                    mapPrint(map1);
                    cout << "!!!!!!!!!!!!!!!!!!!" << endl;
                  }
#endif
                  it2 = window2.erase(it2);
                  deleted2 = true;
                }
              }
            }
          }
        }
        if(!deleted2) ++it2;
      }
    }
    if(!deleted1) ++it;
  }
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::deleteDuplicateCommMaps() {
  // Communication BCs exist in duplicate; the following algorithm should consistently remove both maps

  // Loop over all comm Bcs
  for(auto it = globalStructuredBndryCndMaps.cbegin(); it != globalStructuredBndryCndMaps.cend();) {
    if(globalStructuredBndryCndMaps.empty()) {
      break;
    }

    MBool deleted1 = false;
    if((*it)->BC <= -6000 && (*it)->BC > -6010) {
      for(auto it2 = it + 1; it2 != globalStructuredBndryCndMaps.cend();) {
        if(globalStructuredBndryCndMaps.empty()) {
          break;
        }

        MBool deleted2 = false;
        if((*it2)->BC <= -6000 && (*it2)->BC > -6010) {
          if((*it)->Nstar != (*it2)->Nstar) {
            ++it2;
            continue;
          }

          // Compare if both BCs describe the communication bewteen same blocks
          // TODO_SS labels:FV,toenhance Isn't it enough to check just for the receive blocks, i.g., to check that
          //      a blocks receives data for same range of cells multiple times? Afterwards also delete
          //      the inverted maps
          if((*it)->Id1 == (*it2)->Id1 && (*it)->Id2 == (*it2)->Id2) {
            // Determine if one BC is fully included in the other one;
            // Assumption: When on e map is fully included in another map on one side
            //             then it is also included on the opposite side
            MInt diffStart[nDim], diffEnd[nDim];
            for(MInt dim = 0; dim < nDim; ++dim) {
              // Sanity check that map is normalized
              if((*it)->step1[dim] < 0 || (*it2)->step1[dim] < 0) mTerm(1, "");
              diffStart[dim] = (*it)->start1[dim] - (*it2)->start1[dim];
              diffEnd[dim] = (*it)->end1[dim] - (*it2)->end1[dim];
            }
            const MInt minStart = *std::min_element(&diffStart[0], &diffStart[0] + nDim);
            const MInt maxStart = *std::max_element(&diffStart[0], &diffStart[0] + nDim);
            const MInt minEnd = *std::min_element(&diffEnd[0], &diffEnd[0] + nDim);
            const MInt maxEnd = *std::max_element(&diffEnd[0], &diffEnd[0] + nDim);
            if(minStart * maxStart >= 0 && minEnd * maxEnd >= 0
               && ((minStart >= 0 && maxEnd <= 0) || (maxStart <= 0 && minEnd >= 0))) {
              //            if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && (minStart>=maxEnd||maxStart<=minEnd)) {
              //            if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && minStart*minEnd<=0) {
              // TODO_SS labels:FV By now if two maps cover exactly same range, do nothing
              if(!(minStart == 0 && maxStart == 0 && minEnd == 0 && maxEnd == 0)) {
                mTerm(1, "Duplicates are supposed to be deleted already in deleteDuplicateWindows!");

                // Check which map is the enclosed map
                //                if (maxStart>0 || minEnd<0) {
                if(minStart >= 0 && maxEnd <= 0) {
#ifndef NDEBUG
                  if(domainId() == 0) {
                    cout << "!!! DELETE MAP: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
                    mapPrint(map1);
                    cout << "!!! BECAUSE OF: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
                    mapPrint(map2);
                    cout << "!!!!!!!!!!!!!!!!!!!" << endl;
                  }
#endif
                  it = globalStructuredBndryCndMaps.erase(it);
                  deleted1 = true;
                  break;
                } else {
#ifndef NDEBUG
                  if(domainId() == 0) {
                    cout << "!!! DELETE MAP2: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
                    mapPrint(map2);
                    cout << "!!! BECAUSE OF: !!!" << endl;
                    cout << minStart << "|" << maxStart << " " << minEnd << "|" << maxEnd << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
                    mapPrint(map1);
                    cout << "!!!!!!!!!!!!!!!!!!!" << endl;
                  }
#endif
                  it2 = globalStructuredBndryCndMaps.erase(it2);
                  deleted2 = true;
                }
              }
            }
          }
        }
        if(!deleted2) ++it2;
      }
    }
    if(!deleted1) ++it;
  }

  // Similar story for periodic BCs: Check if receive maps overlap and delete the one which is
  // completly contained in some other periodic map

  for(auto it = globalStructuredBndryCndMaps.cbegin(); it != globalStructuredBndryCndMaps.cend();) {
    if(globalStructuredBndryCndMaps.empty()) {
      break;
    }

    MBool deleted1 = false;
    if((*it)->BC >= 4401 && (*it)->BC <= 4406) {
      for(auto it2 = it + 1; it2 != globalStructuredBndryCndMaps.cend();) {
        if(globalStructuredBndryCndMaps.empty()) {
          break;
        }

        MBool deleted2 = false;
        if((*it2)->BC >= 4401 && (*it2)->BC <= 4406) {
          // Compare receive blocks
          if((*it)->Id2 == (*it2)->Id2) {
            // Determine if one BC is fully included in the other one;
            // Assumption: When on e map is fully included in another map on one side
            //             then it is also included on the opposite side
            MInt diffStart[nDim], diffEnd[nDim];
            for(MInt dim = 0; dim < nDim; ++dim) {
              auto it_start2 = (*it)->step2[dim] > 0 ? (*it)->start2[dim] : (*it)->end2[dim];
              auto it_end2 = (*it)->step2[dim] > 0 ? (*it)->end2[dim] : (*it)->start2[dim];
              auto it2_start2 = (*it2)->step2[dim] > 0 ? (*it2)->start2[dim] : (*it2)->end2[dim];
              auto it2_end2 = (*it2)->step2[dim] > 0 ? (*it2)->end2[dim] : (*it2)->start2[dim];
              diffStart[dim] = it_start2 - it2_start2;
              diffEnd[dim] = it_end2 - it2_end2;
            }


            const MInt minStart = *std::min_element(&diffStart[0], &diffStart[0] + nDim);
            const MInt maxStart = *std::max_element(&diffStart[0], &diffStart[0] + nDim);
            const MInt minEnd = *std::min_element(&diffEnd[0], &diffEnd[0] + nDim);
            const MInt maxEnd = *std::max_element(&diffEnd[0], &diffEnd[0] + nDim);
            if(minStart * maxStart >= 0 && minEnd * maxEnd >= 0
               && ((minStart >= 0 && maxEnd <= 0) || (maxStart <= 0 && minEnd >= 0))) {
              //            if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && (minStart>=maxEnd||maxStart<=minEnd)) {
              //            if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && minStart*minEnd<=0) {
              // TODO_SS labels:FV By now if two maps cover exactly same range, do nothing
              if(!(minStart == 0 && maxStart == 0 && minEnd == 0 && maxEnd == 0)) {
                mTerm(1, "Duplicates are supposed to be deleted already in deleteDuplicateWindows!");

                // Check which map is the enclosed map
                //                if (maxStart>0 || minEnd<0) {
                if(minStart >= 0 && maxEnd <= 0) {
#ifndef NDEBUG
                  if(domainId() == 0) {
                    cout << "!!! DELETE MAP: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
                    mapPrint(map1);
                    cout << minStart << "|" << maxStart << " " << minEnd << "|" << maxEnd << endl;
                    cout << "!!! BECAUSE OF: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
                    mapPrint(map2);
                    cout << "!!!!!!!!!!!!!!!!!!!" << endl;
                  }
#endif
                  it = globalStructuredBndryCndMaps.erase(it);
                  deleted1 = true;
                  break;
                } else {
#ifndef NDEBUG
                  if(domainId() == 0) {
                    cout << "!!! DELETE MAP2: !!!" << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
                    mapPrint(map2);
                    cout << "!!! BECAUSE OF: !!!" << endl;
                    cout << minStart << "|" << maxStart << " " << minEnd << "|" << maxEnd << endl;
                    const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
                    mapPrint(map1);
                    cout << "!!!!!!!!!!!!!!!!!!!" << endl;
                  }
#endif
                  it2 = globalStructuredBndryCndMaps.erase(it2);
                  deleted2 = true;
                }
              }
            }
          }
        }
        if(!deleted2) ++it2;
      }
    }
    if(!deleted1) ++it;
  }
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::deleteDuplicateBCMaps(
    std::vector<unique_ptr<StructuredWindowMap<nDim>>>& BCMap) {
  for(auto it = BCMap.cbegin(); it != BCMap.cend();) {
    if(BCMap.empty()) {
      break;
    }

    MBool deleted1 = false;
    for(auto it2 = it + 1; it2 != BCMap.cend();) {
      if(BCMap.empty()) {
        break;
      }

      MBool deleted2 = false;

      // It should not matter if I compare "1"-indices or "2"-indices.
      // We could also check here if both maps have same BC, before
      // deciding which one to delete

      // Determine if one BC is fully included in the other one;
      // Assumption: When on e map is fully included in another map on one side
      //             then it is also included on the opposite side
      MInt diffStart[nDim], diffEnd[nDim];
      for(MInt dim = 0; dim < nDim; ++dim) {
        diffStart[dim] = (*it)->start1[dim] - (*it2)->start1[dim];
        diffEnd[dim] = (*it)->end1[dim] - (*it2)->end1[dim];
      }
      const MInt minStart = *std::min_element(&diffStart[0], &diffStart[0] + nDim);
      const MInt maxStart = *std::max_element(&diffStart[0], &diffStart[0] + nDim);
      const MInt minEnd = *std::min_element(&diffEnd[0], &diffEnd[0] + nDim);
      const MInt maxEnd = *std::max_element(&diffEnd[0], &diffEnd[0] + nDim);
      //      if (minStart*maxStart>=0 && minEnd*maxEnd>=0 && (minStart>=maxEnd||maxStart<=minEnd)) {
      if(minStart * maxStart >= 0 && minEnd * maxEnd >= 0
         && ((minStart >= 0 && maxEnd <= 0) || (maxStart <= 0 && minEnd >= 0))) {
        // TODO_SS labels:FV By now if two maps cover exactly same range, do nothing
        if(!(minStart == 0 && maxStart == 0 && minEnd == 0 && maxEnd == 0)) {
          // Check which map is the enclosed map
          //          if (maxStart>0 || minEnd<0) {
          if(minStart >= 0 && maxEnd <= 0) {
#ifndef NDEBUG
            if(domainId() == 0) {
              cout << "!!! DELETE MAP: !!!" << endl;
              const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
              mapPrint(map1);
              cout << "!!! BECAUSE OF: !!!" << endl;
              const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
              mapPrint(map2);
              cout << "!!!!!!!!!!!!!!!!!!!" << endl;
            }
#endif
            it = BCMap.erase(it);
            deleted1 = true;
            break;
          } else {
#ifndef NDEBUG
            if(domainId() == 0) {
              cout << "!!! DELETE MAP2: !!!" << endl;
              const unique_ptr<StructuredWindowMap<nDim>>& map2 = (*it2);
              mapPrint(map2);
              cout << "!!! BECAUSE OF: !!!" << endl;
              const unique_ptr<StructuredWindowMap<nDim>>& map1 = (*it);
              mapPrint(map1);
              cout << "!!!!!!!!!!!!!!!!!!!" << endl;
            }
#endif
            it2 = BCMap.erase(it2);
            deleted2 = true;
          }
        }
      }
      if(!deleted2) ++it2;
    }
    if(!deleted1) ++it;
  }
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::readWindowInfo() {
  m_log << "-->reading windowInformation from file ...(serial->parallel approach)" << endl;
  if(globalDomainId() == 0) {
    ParallelIoHdf5 pio(m_grid->m_gridInputFileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);
    stringstream dummy1;
    // go to the windowInformation folder and get number of windows
    noInputWindowInformation = 0;
    pio.getAttribute(&noInputWindowInformation, "noWindows", "/WindowInformation");
    // pio.getAttribute(&noInputWindowConnections, "noConnec", "Connectivity");
    noInputWindowConnections = 0;
    pio.getAttribute(&noInputBndryCnds, "noBndryCnds", "/BC");

    MIntScratchSpace windowInfo(noInputWindowInformation * 7, AT_, "widnow information");

    MInt indices[MAX_SPACE_DIMENSIONS][2]; // stores the indices of the windowInformation
    // copy the data of the windowInformations into a Scratspace which is send to everyone
    for(MInt i = 0; i < noInputWindowInformation; i++) {
      stringstream windowdummy;
      windowdummy << (i + 1);
      MString window = "window" + windowdummy.str();
      MString window1 = "/WindowInformation/window" + windowdummy.str();
      pio.getAttribute(&(windowInfo[i * 7]), "block", window1);
      ParallelIo::size_type ioSize[2] = {nDim, 2};
      pio.readArray(&(windowInfo[i * 7 + 1]), "WindowInformation", window, 2, ioSize);
    }

    MInt noBcWindows;
    MInt** BCInformation = new MInt*[mMax(noInputBndryCnds, 0)];
    MInt* bcCndId = new MInt[mMax(noInputBndryCnds, 0)];
    MInt* noBcCndWindows = new MInt[mMax(noInputBndryCnds, 0)];

    for(MInt i = 0; i < noInputBndryCnds; i++) {
      stringstream bcDummy;
      stringstream bcDummy1;
      bcDummy << "/BC/BC" << (i + 1);
      bcDummy1 << "BC" << (i + 1);
      char bcIdWindows[80];
      char bcId[80];
      strcpy(bcIdWindows, &bcDummy.str()[0]);
      strcpy(bcId, &bcDummy.str()[0]);
      strcat(bcIdWindows, "/noWindows");
      pio.getAttribute(&noBcWindows, "noWindows", bcDummy.str());
      noBcCndWindows[i] = noBcWindows;
      BCInformation[i] = new MInt[mMax((MInt)noBcWindows, 0)];
      ParallelIo::size_type ioSize = noBcWindows;
      pio.readArray(&BCInformation[i][0], "BC", bcDummy1.str(), 1, &ioSize);
      strcat(bcId, "/BC");
      pio.getAttribute(&bcCndId[i], "BC", bcDummy.str());
    }

    MInt noWindows = 0;
    for(MInt i = 0; i < noInputBndryCnds; i++) {
      noWindows += noBcCndWindows[i];
    }
    // put everything into a vector to send the data
    MInt count = 0;
    MIntScratchSpace bcInfo(noWindows + 3 * noInputBndryCnds, AT_, "send bc info");
    for(MInt i = 0; i < noInputBndryCnds; i++) {
      bcInfo[count] = bcCndId[i];
      count++;
      bcInfo[count] = noBcCndWindows[i];
      count++;
      memcpy(&(bcInfo[count]), &(BCInformation[i][0]), noBcCndWindows[i] * sizeof(MInt));
      count += noBcCndWindows[i];
      bcInfo[count] = bcCndId[i];
      count++;
    }

    MInt countInfo[3] = {noInputWindowInformation, noInputBndryCnds, noWindows};
    MPI_Bcast(countInfo, 3, MPI_INT, 0, m_StructuredComm, AT_, "countInfo");

    // broadcast the information
    MPI_Bcast(windowInfo.getPointer(), noInputWindowInformation * 7, MPI_INT, 0, m_StructuredComm, AT_,
              "windowInfo.getPointer()");
    MPI_Bcast(bcInfo.getPointer(), (noWindows + 3 * noInputBndryCnds), MPI_INT, 0, m_StructuredComm, AT_,
              "bcInfo.getPointer()");

    // put the window information into the right place
    for(MInt i = 0; i < noInputWindowInformation; i++) {
      std::unique_ptr<windowInformation<nDim>> temp = make_unique<windowInformation<nDim>>();
      // fill the windowInformation with life
      // first attribute the windowId to the object
      temp->windowId = i;
      temp->blockId = windowInfo[i * 7];
      memcpy(&indices[0][0], &(windowInfo[i * 7 + 1]), MAX_SPACE_DIMENSIONS * 2 * sizeof(MInt));
      for(MInt dim = 0; dim < nDim; dim++) {
        temp->startindex[dim] = indices[dim][0] - 1;
        temp->endindex[dim] = indices[dim][1] - 1;
      }
      // add the window to the vector
      inputWindows.push_back(std::move(temp));
    }
    // distribute the boundary conditions to the windows!
    for(MInt i = 0; i < noInputBndryCnds; ++i) {
      for(MInt j = 0; j < noBcCndWindows[i]; ++j) {
        MInt windowId = BCInformation[i][j] - 1; // correct the fortran notation
        inputWindows[windowId]->BC = bcCndId[i];
      }
    }

    delete[] bcCndId;
    delete[] noBcCndWindows;
    for(MInt i = 0; i < noInputBndryCnds; i++) {
      delete[] BCInformation[i];
    }
    delete[] BCInformation;
  } else {
    // receive the data first
    MInt countInfo[3] = {0, 0, 0};
    MPI_Bcast(countInfo, 3, MPI_INT, 0, m_StructuredComm, AT_, "countInfo");
    noInputWindowInformation = countInfo[0];
    noInputBndryCnds = countInfo[1];
    MInt noWindows = countInfo[2];
    MIntScratchSpace bcInfo(noWindows + 3 * noInputBndryCnds, AT_, "send bc info");
    MIntScratchSpace windowInfo(noInputWindowInformation * 7, AT_, "widnow information");
    // receive the data from the other processes
    MPI_Bcast(windowInfo.getPointer(), noInputWindowInformation * 7, MPI_INT, 0, m_StructuredComm, AT_,
              "windowInfo.getPointer()");
    MPI_Bcast(bcInfo.getPointer(), (noWindows + 3 * noInputBndryCnds), MPI_INT, 0, m_StructuredComm, AT_,
              "bcInfo.getPointer()");
    MInt indices[MAX_SPACE_DIMENSIONS][2]; // stores the indices of the windowInformation
    // put the window information into the right place
    for(MInt i = 0; i < noInputWindowInformation; i++) {
      std::unique_ptr<windowInformation<nDim>> temp = make_unique<windowInformation<nDim>>();
      // fill the windowInformation with life
      // first attribute the windowId to the object
      temp->windowId = i;
      temp->blockId = windowInfo[i * 7];
      memcpy(&indices[0][0], &(windowInfo[i * 7 + 1]), MAX_SPACE_DIMENSIONS * 2 * sizeof(MInt));
      for(MInt dim = 0; dim < nDim; dim++) {
        temp->startindex[dim] = indices[dim][0] - 1;
        temp->endindex[dim] = indices[dim][1] - 1;
      }
      // add the window to the vector
      inputWindows.push_back(std::move(temp));
    }

    // put the boundary condition info into the right place
    MInt** BCInformation = new MInt*[mMax(noInputBndryCnds, 0)];
    MInt* bcCndId = new MInt[mMax(noInputBndryCnds, 0)];
    MInt* noBcCndWindows = new MInt[mMax(noInputBndryCnds, 0)];
    MInt count = 0;
    for(MInt i = 0; i < noInputBndryCnds; i++) {
      bcCndId[i] = bcInfo[count];
      count++;
      noBcCndWindows[i] = bcInfo[count];
      count++;
      BCInformation[i] = new MInt[mMax(noBcCndWindows[i], 0)];
      memcpy(&(BCInformation[i][0]), &(bcInfo[count]), noBcCndWindows[i] * sizeof(MInt));
      count += noBcCndWindows[i];
      bcCndId[i] = bcInfo[count];
      count++;
    }

    // distribute the boundary conditions to the windows!
    for(MInt i = 0; i < noInputBndryCnds; ++i) {
      for(MInt j = 0; j < noBcCndWindows[i]; ++j) {
        MInt windowId = BCInformation[i][j] - 1; // correct the fortran notation
        inputWindows[windowId]->BC = bcCndId[i];
      }
    }

    delete[] bcCndId;
    delete[] noBcCndWindows;
    for(MInt i = 0; i < noInputBndryCnds; i++) {
      delete[] BCInformation[i];
    }
    delete[] BCInformation;
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::initGlobals() {
  unique_ptr<StructuredWindowMap<nDim>> windowMap;
  // MInt noPeriodicWindows[3] = {0,0,0};
  // function creates a map of all the global boundary conditions
  // global means from Input file !!!

  // first look which periodic direction is not used
  // MBool switchPeriodicIds = true;
  MBool isPeriodicDirUsed[3] = {false, false, false};
  MBool switchPeriodicBC[3] = {true, true, true};

  // check if periodic direction (4401,4403,4405) is already in use (find free slots for channel bc)
  // also check if distinct bc numbers (4401/4402, 4403/4404, 4405/4406) for both surfaces are used
  // or if both periodic surfaces have the same number, in that case swtich
  for(MInt windowId = 0; windowId < noInputWindowInformation; windowId++) {
    switch(inputWindows[windowId]->BC) {
      case 4401:
        isPeriodicDirUsed[0] = true;
        break;
      case 4402:
        switchPeriodicBC[0] = false;
        break;
      case 4403:
        isPeriodicDirUsed[1] = true;
        break;
      case 4404:
        switchPeriodicBC[1] = false;
        break;
      case 4405:
        isPeriodicDirUsed[2] = true;
        break;
      case 4406:
        switchPeriodicBC[2] = false;
        break;
      default:
        // do nothing
        break;
    }
  }

  for(MInt noWin = 0; noWin < noInputWindowInformation; noWin++) {
    windowMap = make_unique<StructuredWindowMap<nDim>>();
    MInt order[3] = {0, 1, 2};
    MInt step1[3] = {1, 1, 1};

    // use the standard mapping to itself
    // pushback all information on the boundary condition into the global container
    // except for periodic and multiblock communication
    if(inputWindows[noWin]->BC != 6000
       && (inputWindows[noWin]->BC < 4000
           || inputWindows[noWin]->BC >= 5000) /*&&inputWindows[noWin]->BC!=2222&&inputWindows[noWin]->BC!=2221*/) {
      mapCreate(inputWindows[noWin]->blockId, inputWindows[noWin]->startindex, inputWindows[noWin]->endindex, step1,
                (inputWindows[noWin]->windowId + 1), nullptr, nullptr, nullptr, order, inputWindows[noWin]->BC,
                windowMap);
      // add the map to the list!
      globalStructuredBndryCndMaps.push_back(std::move(windowMap));
    }

    // Pascal: Why do we have a spcial treatment for that? this is exactely what is done above
    // special treatment for zonal BC 2222
    /*if(inputWindows[noWin]->BC==2222) {
      mapCreate(inputWindows[noWin]->blockId, inputWindows[noWin]->startindex, inputWindows[noWin]->endindex, step1,
      (inputWindows[noWin]->windowId+1) ,  nullptr, nullptr , nullptr, order,inputWindows[noWin]->BC , windowMap );
      //add map for bc2222 to list
      globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      }
    */
    // special treatment for zonal BC 2221
    if(inputWindows[noWin]->BC == 2221) {
      /*mapCreate(inputWindows[noWin]->blockId, inputWindows[noWin]->startindex, inputWindows[noWin]->endindex, step1,
        (inputWindows[noWin]->windowId+1) ,  nullptr, nullptr , nullptr, order,inputWindows[noWin]->BC , windowMap );
        //add map for bc2221 to list
        globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      */
      // now add a copy of the boundary condition of type bc7909 to list (special treatment)
      windowMap = make_unique<StructuredWindowMap<nDim>>();
      mapCreate(inputWindows[noWin]->blockId, inputWindows[noWin]->startindex, inputWindows[noWin]->endindex, step1,
                (inputWindows[noWin]->windowId + 1), nullptr, nullptr, nullptr, order, 7909, windowMap);
      // add the map to the list!
      globalStructuredBndryCndMaps.push_back(std::move(windowMap));
    }


    if(inputWindows[noWin]->BC > 4000 && inputWindows[noWin]->BC < 5000) {
      // add channel flow and rotating periodic BC to physical BC.
      // in order to build the comm groups for channel and rotation.
      //!!!important set numbers for channel and rotation.
      // 4001 and 4002 for rotation
      // 2401 and 2402 for channel

      if(inputWindows[noWin]->BC == 4011) {
        mapCreate(inputWindows[noWin]->blockId, inputWindows[noWin]->startindex, inputWindows[noWin]->endindex, step1,
                  (inputWindows[noWin]->windowId + 1), nullptr, nullptr, nullptr, order, 4001, windowMap);
        globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      }

      if(inputWindows[noWin]->BC == 4012) {
        mapCreate(inputWindows[noWin]->blockId, inputWindows[noWin]->startindex, inputWindows[noWin]->endindex, step1,
                  (inputWindows[noWin]->windowId + 1), nullptr, nullptr, nullptr, order, 4002, windowMap);
        globalStructuredBndryCndMaps.push_back(std::move(windowMap));
      }

      if(inputWindows[noWin]->BC > 4400 && inputWindows[noWin]->BC < 4407) {
        // if same BC number is used for both periodic surfaces
        // switch BC number of the surface at the end of the blocks
        MInt periodicDirection = -1;
        for(MInt dim = 0; dim < nDim; dim++) {
          if(inputWindows[noWin]->startindex[dim] == inputWindows[noWin]->endindex[dim]) {
            periodicDirection = dim;
            break;
          }
        }

        MInt bigBlockId = inputWindows[noWin]->blockId;
        MInt bigBlockSize = m_grid->getBlockNoCells(bigBlockId, nDim - 1 - periodicDirection);
        MInt periodicId = inputWindows[noWin]->BC - 4401;
        MInt periodicIndex = periodicId / 2;

        if(switchPeriodicBC[periodicIndex] && bigBlockSize == inputWindows[noWin]->startindex[periodicDirection]
           && periodicId % 2 == 0) {
          m_log << "Changing BC from " << inputWindows[noWin]->BC << " to " << inputWindows[noWin]->BC + 1 << endl;
          inputWindows[noWin]->BC++;
        }
      }

      // Channel BCs
      // for BC 4005 and 4006 two periodic bcs 4401 and 4402
      // will be created and also physical bcs 2401 and 2402
      // to correct the pressure and density at both sides
      if(inputWindows[noWin]->BC == 4005 || inputWindows[noWin]->BC == 4006) {
        // use the first unused periodic direction
        MInt freePeriodicBC = 4400;
        for(MInt dir = 0; dir < nDim; dir++) {
          if(isPeriodicDirUsed[dir] == false) {
            freePeriodicBC += dir * 2 + 1;
            break;
          }
        }

        m_log << "Using Periodic BC " << freePeriodicBC << " / " << freePeriodicBC + 1 << " for the channel/pipe flow"
              << endl;

        if(inputWindows[noWin]->BC == 4005) {
          m_log << "Identified channel inflow bc 4005, creating periodic map BC " << freePeriodicBC
                << " and physicalMap BC " << 2401 << endl;
          inputWindows[noWin]->BC = freePeriodicBC;
          mapCreate(inputWindows[noWin]->blockId, inputWindows[noWin]->startindex, inputWindows[noWin]->endindex, step1,
                    (inputWindows[noWin]->windowId + 1), nullptr, nullptr, nullptr, order, 2401, windowMap);
          globalStructuredBndryCndMaps.push_back(std::move(windowMap));
        }

        if(inputWindows[noWin]->BC == 4006) {
          m_log << "Identified channel inflow bc 4006, creating periodic map BC " << freePeriodicBC + 1
                << " and physicalMap BC " << 2402 << endl;
          inputWindows[noWin]->BC = freePeriodicBC + 1;
          mapCreate(inputWindows[noWin]->blockId, inputWindows[noWin]->startindex, inputWindows[noWin]->endindex, step1,
                    (inputWindows[noWin]->windowId + 1), nullptr, nullptr, nullptr, order, 2402, windowMap);
          globalStructuredBndryCndMaps.push_back(std::move(windowMap));
        }
      }
    }
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::writeConnectionWindowInformation3D(MFloat* periodicDisplacements) {
  /////////////////////////////////
  /////// NORMAL CONNECTIONS //////
  /////////////////////////////////
  ParallelIoHdf5 pio(m_grid->m_gridInputFileName, maia::parallel_io::PIO_APPEND, m_StructuredComm);
  MInt noConnectionMaps = 0;
  MInt empty = 0;
  for(MInt mapId = 0; mapId < (MInt)globalStructuredBndryCndMaps.size(); mapId++) {
    // The 6000er connection maps have negative BC?!
    // only write out the connection maps
    if((globalStructuredBndryCndMaps[mapId]->BC <= -6000 && globalStructuredBndryCndMaps[mapId]->BC > -6010)
       || // (globalStructuredBndryCndMaps[mapId]->BC==6000) || //6000er-change ; check if this change makes sense!
       (globalStructuredBndryCndMaps[mapId]->BC > 4400 && globalStructuredBndryCndMaps[mapId]->BC < 4407)) {
      stringstream windowName;
      windowName << "map" << noConnectionMaps << endl;
      MString windowNameStr = windowName.str();
      ParallelIo::size_type noWindowInformation = 7 * nDim + 16;
      pio.defineArray(maia::parallel_io::PIO_INT, "Connectivity", windowNameStr, 1, &noWindowInformation);

      MIntScratchSpace windowInformation((MInt)noWindowInformation, AT_, "windowInformation");
      writeMapToArray(globalStructuredBndryCndMaps[mapId], windowInformation.begin());


      ParallelIo::size_type offset = 0;
      if(globalDomainId() == 0) {
        pio.writeArray(windowInformation.begin(), "Connectivity", windowNameStr, 1, &offset, &noWindowInformation);
      } else {
        ParallelIo::size_type zeroWindows = 0;
        pio.writeArray(&empty, "Connectivity", windowNameStr, 1, &offset, &zeroWindows);
      }
      noConnectionMaps++;
    }
  }

  /////////////////////////////////
  /////// SINGULARITIES ///////////
  /////////////////////////////////
  MInt noSingularityMaps = 0;
  for(MInt mapId = 0; mapId < (MInt)singularwindow.size(); mapId++) {
    stringstream windowName;
    windowName << "singularMap" << noSingularityMaps << endl;
    MString windowNameStr = windowName.str();
    ParallelIo::size_type noWindowInformation = 7 * nDim + 16;
    pio.defineArray(maia::parallel_io::PIO_INT, "Connectivity", windowNameStr, 1, &noWindowInformation);

    MIntScratchSpace windowInformation(noWindowInformation, AT_, "windowInformation");
    writeMapToArray(singularwindow[mapId], windowInformation.begin());

    ParallelIo::size_type offset = 0;
    if(domainId() == 0) {
      pio.writeArray(windowInformation.begin(), "Connectivity", windowNameStr, 1, &offset, &noWindowInformation);
    } else {
      ParallelIo::size_type zeroWindows = 0;
      pio.writeArray(&empty, "Connectivity", windowNameStr, 1, &offset, &zeroWindows);
    }
    noSingularityMaps++;
  }

  /////////////////////////////////
  ///// PERIODIC DISPLACEMENTS ////
  /////////////////////////////////
  cout << "Starting write out " << endl;
  MInt noPeriodicDisplacementInfo = nDim * nDim;
  MFloatScratchSpace localPeriodicDisplacements(m_noBlocks * noPeriodicDisplacementInfo, AT_,
                                                "allPeriodicDisplacements");
  MFloatScratchSpace globalPeriodicDisplacements(m_noBlocks * noPeriodicDisplacementInfo, AT_,
                                                 "allPeriodicDisplacements");
  localPeriodicDisplacements.fill(-99999999.0);
  for(MInt periodicId = 0; periodicId < noPeriodicDisplacementInfo; periodicId++) {
    localPeriodicDisplacements(m_blockId * noPeriodicDisplacementInfo + periodicId) = periodicDisplacements[periodicId];
  }

  cout << "Allreduce noBlocks: " << m_noBlocks << " totalSize: " << m_noBlocks * noPeriodicDisplacementInfo << endl;
  MPI_Allreduce(&localPeriodicDisplacements(0), &globalPeriodicDisplacements(0),
                m_noBlocks * noPeriodicDisplacementInfo, MPI_DOUBLE, MPI_MAX, m_StructuredComm, AT_,
                "localPeriodicDisplacements(0)", "globalPeriodicDisplacements(0)");

  cout << "Writing out periodic Window displacements" << endl;

  for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
    stringstream path;
    path << "periodicDisplacements" << blockId;
    MString pathStr = path.str();
    ParallelIo::size_type ioSize = noPeriodicDisplacementInfo;
    if(!pio.hasDataset(pathStr, "Connectivity")) {
      pio.defineArray(maia::parallel_io::PIO_FLOAT, "Connectivity", pathStr, 1, &ioSize);
    }
    ParallelIo::size_type ioOffset = 0;
    if(domainId() == 0) {
      pio.writeArray(&globalPeriodicDisplacements[blockId * noPeriodicDisplacementInfo], "Connectivity", pathStr, 1,
                     &ioOffset, &ioSize);
    } else {
      ParallelIo::size_type zeroWindows = 0;
      MFloat emptyVar = 0;
      pio.writeArray(&emptyVar, "Connectivity", pathStr, 1, &ioOffset, &zeroWindows);
    }
  }

  MInt hasConnectionInfo = 1;
  m_log << "Connection info written to grid file, normalConnections: " << noConnectionMaps
        << " noSingularityMaps: " << noSingularityMaps << endl;

  pio.setAttribute(hasConnectionInfo, "hasConnectionInfo", "Connectivity");
  pio.setAttribute(noConnectionMaps, "noConnections", "Connectivity");
  pio.setAttribute(noSingularityMaps, "noSingularConnections", "Connectivity");
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::readConnectionWindowInformation3D(MFloat* periodicDisplacements) {
  MInt noConnectionMaps = 0;
  MInt noSingularityMaps = 0;
  MInt noPeriodicDisplacementInfo = nDim * nDim;
  MFloatScratchSpace globalPeriodicDisplacements(m_noBlocks * noPeriodicDisplacementInfo, AT_,
                                                 "allPeriodicDisplacements");

  // first read the number of connection information
  // and broadcast it to all partitions
  if(globalDomainId() == 0) {
    ParallelIoHdf5 pio(m_grid->m_gridInputFileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);
    pio.getAttribute(&noConnectionMaps, "noConnections", "Connectivity");
    pio.getAttribute(&noSingularityMaps, "noSingularConnections", "Connectivity");
    MInt dummy[2] = {noConnectionMaps, noSingularityMaps};
    MPI_Bcast(dummy, 2, MPI_INT, 0, m_StructuredComm, AT_, "dummy");
  } else {
    MInt dummy[2] = {0, 0};
    MPI_Bcast(dummy, 2, MPI_INT, 0, m_StructuredComm, AT_, "dummy");
    noConnectionMaps = dummy[0];
    noSingularityMaps = dummy[1];
  }

  // if connections were looked for previously (hasConnectionInfo is set)
  // but there are no connections there is no need to read them
  if(noConnectionMaps == 0 && noSingularityMaps == 0) {
    if(globalDomainId() == 0) {
      cout << "Connections were previously searched for, none found!" << endl
           << "---> hasConnectionInfo = 1" << endl
           << "---> noConnectionMaps = 0" << endl
           << "---> noSingularityMaps = 0" << endl;
    }
    return;
  }

  if(globalDomainId() == 0) {
    ParallelIoHdf5 pio(m_grid->m_gridInputFileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);
    MInt noWindowInformation = 7 * nDim + 16;
    if(noConnectionMaps > 0) {
      MIntScratchSpace windowInformation(noWindowInformation * noConnectionMaps, AT_, "windowInformation");
      for(MInt mapId = 0; mapId < noConnectionMaps; mapId++) {
        stringstream windowName;
        windowName << "map" << mapId << endl;
        MString windowNameStr = windowName.str();
        ParallelIo::size_type ioSize = noWindowInformation;
        pio.readArray(&(windowInformation[mapId * noWindowInformation]), "Connectivity", windowNameStr, 1, &ioSize);
      }
      MPI_Bcast(windowInformation.getPointer(), noWindowInformation * noConnectionMaps, MPI_INT, 0, m_StructuredComm,
                AT_, "windowInformation.getPointer()");

      for(MInt mapId = 0; mapId < noConnectionMaps; mapId++) {
        unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
        readMapFromArray(temp, &windowInformation[mapId * noWindowInformation]);
        globalStructuredBndryCndMaps.push_back(std::move(temp));
      }
    }

    if(noSingularityMaps > 0) {
      MIntScratchSpace singularityInformation(noWindowInformation * noSingularityMaps, AT_, "singularityInformation");
      for(MInt mapId = 0; mapId < noSingularityMaps; mapId++) {
        stringstream windowName;
        windowName << "singularMap" << mapId << endl;
        MString windowNameStr = windowName.str();
        ParallelIo::size_type ioSize = noWindowInformation;
        pio.readArray(&(singularityInformation[mapId * noWindowInformation]), "Connectivity", windowNameStr, 1,
                      &ioSize);
      }
      MPI_Bcast(singularityInformation.getPointer(), noWindowInformation * noSingularityMaps, MPI_INT, 0,
                m_StructuredComm, AT_, "singularityInformation.getPointer()");
      for(MInt mapId = 0; mapId < noSingularityMaps; mapId++) {
        unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
        readMapFromArray(temp, &(singularityInformation[mapId * noWindowInformation]));
        singularwindow.push_back(std::move(temp));
      }
    }

    for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
      stringstream path;
      path << "periodicDisplacements" << blockId;
      MString pathStr = path.str();
      ParallelIo::size_type ioSize = noPeriodicDisplacementInfo;
      if(pio.hasDataset(pathStr, "Connectivity")) {
        // check if path for block individual displacements exists, if not read
        // same dataset for every block
        pio.readArray(&globalPeriodicDisplacements[blockId * noPeriodicDisplacementInfo], "Connectivity", pathStr, 1,
                      &ioSize);
      } else {
        pio.readArray(&globalPeriodicDisplacements[blockId * noPeriodicDisplacementInfo], "Connectivity",
                      "periodicDisplacements", 1, &ioSize);
      }
    }

    MPI_Bcast(&globalPeriodicDisplacements(0), m_noBlocks * noPeriodicDisplacementInfo, MPI_DOUBLE, 0, m_StructuredComm,
              AT_, "globalPeriodicDisplacements(0)");

    for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
      for(MInt i = 0; i < noPeriodicDisplacementInfo; i++) {
        cout << "globalPeriodicDisplacement: " << globalPeriodicDisplacements[blockId * noPeriodicDisplacementInfo + i]
             << endl;
      }
    }
  } else {
    MInt noWindowInformation = 7 * nDim + 16;
    if(noConnectionMaps > 0) {
      MIntScratchSpace windowInformation(noWindowInformation * noConnectionMaps, AT_, "windowInformation");
      MPI_Bcast(windowInformation.getPointer(), noWindowInformation * noConnectionMaps, MPI_INT, 0, m_StructuredComm,
                AT_, "windowInformation.getPointer()");
      for(MInt mapId = 0; mapId < noConnectionMaps; mapId++) {
        unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
        readMapFromArray(temp, &windowInformation[mapId * noWindowInformation]);
        globalStructuredBndryCndMaps.push_back(std::move(temp));
      }
    }
    if(noSingularityMaps > 0) {
      MIntScratchSpace singularityInformation(noWindowInformation * noSingularityMaps, AT_, "singularityInformation");
      MPI_Bcast(singularityInformation.getPointer(), noWindowInformation * noSingularityMaps, MPI_INT, 0,
                m_StructuredComm, AT_, "singularityInformation.getPointer()");
      for(MInt mapId = 0; mapId < noSingularityMaps; mapId++) {
        unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
        readMapFromArray(temp, &(singularityInformation[mapId * noWindowInformation]));
        singularwindow.push_back(std::move(temp));
      }
    }

    MPI_Bcast(&globalPeriodicDisplacements(0), m_noBlocks * noPeriodicDisplacementInfo, MPI_DOUBLE, 0, m_StructuredComm,
              AT_, "globalPeriodicDisplacements(0)");
  }

  for(MInt periodicId = 0; periodicId < noPeriodicDisplacementInfo; periodicId++) {
    periodicDisplacements[periodicId] =
        globalPeriodicDisplacements(m_blockId * noPeriodicDisplacementInfo + periodicId);
  }
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::setSpongeInformation(MInt noSpongeInfo, MFloat* beta, MFloat* sigma,
                                                              MFloat* thickness, MInt* bcInfo, MInt informationType) {
  // we could also use the window information for sponges
  // set through properties.

  // add the sponge Layter type to it too!!.
  if(informationType) {
    for(MUint bcId = 0; bcId < globalStructuredBndryCndMaps.size(); ++bcId) {
      for(MInt i = 0; i < noSpongeInfo; ++i) {
        if(globalStructuredBndryCndMaps[bcId]->BC == bcInfo[i]) {
          unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
          mapCpy(globalStructuredBndryCndMaps[bcId], temp);
          temp->spongeThickness = thickness[i];
          temp->beta = beta[i];
          temp->sigma = sigma[i];
          m_spongeInfoMap.push_back(std::move(temp));
        }
      }
    }
  } else {
    for(MUint bcId = 0; bcId < globalStructuredBndryCndMaps.size(); ++bcId) {
      for(MInt i = 0; i < noSpongeInfo; ++i) {
        if(globalStructuredBndryCndMaps[bcId]->Id2 == bcInfo[i]) {
          unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
          mapCpy(globalStructuredBndryCndMaps[bcId], temp);
          temp->spongeThickness = thickness[i];
          temp->beta = beta[i];
          temp->sigma = sigma[i];
          m_spongeInfoMap.push_back(std::move(temp));
        }
      }
    }
  }
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::setZonalBCInformation() {
  unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC == 2221 || globalStructuredBndryCndMaps[i]->BC == 2222) {
      mapCpy(globalStructuredBndryCndMaps[i], temp);
      m_zonalBCMaps.push_back(std::move(temp));
      temp = make_unique<StructuredWindowMap<nDim>>();
    }
  }
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::checkZonalBCMaps(unique_ptr<StructuredWindowMap<nDim>>& map1,
                                                           unique_ptr<StructuredWindowMap<nDim>>& map2) {
  if(map1->BC == map2->BC && map1->Id1 == map2->Id1) {
    return true;
  }

  return false;
}


// new function to create list of all wall bcs needed for SA RANS
template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::setWallInformation() {
  unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
  for(MUint bcId = 0; bcId < globalStructuredBndryCndMaps.size(); ++bcId) {
    MInt firstDigit = (MInt)(((MFloat)globalStructuredBndryCndMaps[bcId]->BC) / 1000.0);
    if(firstDigit == 1 || globalStructuredBndryCndMaps[bcId]->BC == 2601) {
      mapCpy(globalStructuredBndryCndMaps[bcId], temp);
      m_wallDistInfoMap.push_back(std::move(temp));
      temp = make_unique<StructuredWindowMap<nDim>>();
    }
  }
}


// TODO_SS labels:FV Check if physicalAuxDataMap has the same information except that physicalAuxDataMap is not
// thickened by the # of ghost cells
// TODO_SS labels:FV,toremove The following is not used and can be deleted
// new function to create list of all wall bcs needed for SA RANS
template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::setLocalWallInformation() {
  unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
  for(MUint bcId = 0; bcId < localStructuredBndryCndMaps.size(); ++bcId) {
    MInt firstDigit = (MInt)(((MFloat)localStructuredBndryCndMaps[bcId]->BC) / 1000.0);
    if(firstDigit == 1 || localStructuredBndryCndMaps[bcId]->BC == 2601) {
      mapCpy(localStructuredBndryCndMaps[bcId], temp);
      temp->face = localStructuredBndryCndMaps[bcId]->face;

      // Undo of thickening of the physicalBCMaps at the edges of the domain
      for(MInt dim = 0; dim < nDim; ++dim) {
        if(dim != (MInt)temp->face / 2) {
          if(temp->start1[dim] == m_myMapWithGC->start2[dim]) {
            temp->start1[dim] += m_noGhostLayers;
            temp->start2[dim] += m_noGhostLayers;
          }
          if(temp->end1[dim] == m_myMapWithGC->end2[dim]) {
            temp->end1[dim] -= m_noGhostLayers;
            temp->end2[dim] -= m_noGhostLayers;
          }
        } else { // dim==face/2
          const MInt face = temp->face;
          const MInt n = (face % 2) * 2 - 1; //-1,+1
          temp->start1[dim] += (MInt)(0.5 - (1.5 * (MFloat)n));
          temp->start2[dim] += (MInt)(0.5 - (1.5 * (MFloat)n));
        }
      }

      m_wallDistInfoMap.push_back(std::move(temp));
      temp = make_unique<StructuredWindowMap<nDim>>();
    }
  }
}


template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::createAuxDataMap(const MInt auxDataType,
                                                          const MString blockType,
                                                          const std::vector<MInt>& porousBlockIds,
                                                          const std::vector<MInt>& auxDataWindowIds,
                                                          const MBool force) {
  ///////////////////////////////////////////////////////////
  //////////////// AUX DATA MAP CORRECTION //////////////////
  ///////////////////////////////////////////////////////////
  // For computing the modified wall distance we need to be on the side of the porous block,
  // but for computing the wall forces (cf & cp) we want to be on the fluid side
  //  const MBool considerFPInterfaces = /*(blockType=="porous")
  //                &&*/ (auxDataType==1 || auxDataType==2 || auxDataType==3);

  // create the auxilary data maps
  // needed to write out cf,cp etc ...
  unique_ptr<StructuredWindowMap<nDim>> localMapDummy = make_unique<StructuredWindowMap<nDim>>();
  MInt counter = 0;
  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    //    if (globalStructuredBndryCndMaps[i]->Nstar!=-1) continue;
    const MInt BC = globalStructuredBndryCndMaps[i]->BC;
    MInt firstDigit = (MInt)(((MFloat)BC) / 1000.0);
    // If not solid and not fluid-porous interface continue
    //    if (firstDigit!=1 && (!considerFPInterfaces || BC!=6002)) continue;
    if(firstDigit != 1 && BC != 6002) continue;


    MBool candidate = false;

    if(auxDataType == 3) {
      for(auto windowId : auxDataWindowIds) {
        if(globalStructuredBndryCndMaps[i]->Id2 == windowId) { // && BC!=6000) {
          candidate = true;
          break;
        }
      }
    } else if(auxDataType == 0) { // only solid
      if(firstDigit == 1) candidate = true;
    } else if(auxDataType == 2) { // only fluid-porous
      if(BC == 6002 && porousBlockIds[globalStructuredBndryCndMaps[i]->Id1] == -1 /*blockType=="fluid"*/)
        candidate = true;
    } else if(firstDigit == 1
              || (BC == 6002
                  && porousBlockIds[globalStructuredBndryCndMaps[i]->Id1]
                         == -1 /*blockType=="fluid"*/)) { // fluid-porous interfaces and solids
      candidate = true;
    }

    if(candidate) m_auxDataWindowIds.insert({i, globalStructuredBndryCndMaps[i]->Id2});

    if(candidate || force) { // all the wall should be covered
      mapCombine11(m_myMapWithoutGC, globalStructuredBndryCndMaps[i], localMapDummy);
      MBool test = false;
      test = mapCheck(localMapDummy);
      if(test == true) {
        physicalAuxDataMap.push_back(std::move(localMapDummy));
        for(MInt j = 0; j < nDim; ++j) {
          // correct the global part for the output!!!
          physicalAuxDataMap[counter]->start2[j] -= globalStructuredBndryCndMaps[i]->start2[j];
          physicalAuxDataMap[counter]->end2[j] -= globalStructuredBndryCndMaps[i]->start2[j];
        }
        if(globalStructuredBndryCndMaps[i]->BC == 6002 && blockType == "porous")
          physicalAuxDataMap[counter]->isFluidPorousInterface = true;
        counter++;
        localMapDummy = make_unique<StructuredWindowMap<nDim>>();
      }
    }
  }

  deleteDuplicateBCMaps(physicalAuxDataMap);

  for(MInt i = 0; i < (MInt)physicalAuxDataMap.size(); i++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      // check if the index is the same
      // MPI_Barrier(m_StructuredComm, AT_ );
      if(physicalAuxDataMap[i]->start1[dim] == physicalAuxDataMap[i]->end1[dim]) {
        // check which face it is and save this information
        switch(dim) {
          case 0: { // face 0 or 1
            // check if face 0 or 1
            if(physicalAuxDataMap[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // test of start2 instead of start1.
              physicalAuxDataMap[i]->face = 0;
            } else {
              physicalAuxDataMap[i]->face = 1;
            }
            break;
          }
          case 1: { // face 2 or 3
            // check if face 2 or 3
            if(physicalAuxDataMap[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              physicalAuxDataMap[i]->face = 2;
            } else {
              physicalAuxDataMap[i]->face = 3;
            }
            break;
          }
          case 2: { // face 4 or 5
            // check if face 4 or 5
            if(physicalAuxDataMap[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              physicalAuxDataMap[i]->face = 4;
            } else {
              physicalAuxDataMap[i]->face = 5;
            }
            break;
          }
          default: {
            cerr << "error no side could be attributed" << endl;
            exit(1);
          }
        }
        continue;
      }
    }
  }
}

//#include <chrono>
//#include <thread>
template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::createWindowMapping(
    MPI_Comm* channelIn, MPI_Comm* channelOut, MPI_Comm* channelWorld, MInt* channelRoots, MPI_Comm* commStg,
    MInt* commStgRoot, MInt* commStgRootGlobal, MPI_Comm* commBC2600, MInt* commBC2600Root, MInt* commBC2600RootGlobal,
    MPI_Comm* rescalingCommGrComm, MInt* rescalingCommGrRoot, MInt* rescalingCommGrRootGlobal, MPI_Comm* commPerRotOne,
    MPI_Comm* commPerRotTwo, MPI_Comm* commPerRotWorld, MInt* rotationRoots, MInt& perRotGroup,
    SingularInformation* singularity, MInt* hasSingularity, MPI_Comm* plenumComm, MInt* plenumRoot) {
  // create a map of the partition and check whether it fits any of the globals
  // boundary Condition Maps
  // thus this functions relates the partition to the actual window and there-
  // fore also to the boundary condition.

  MInt start1[3] = {0, 0, 0};
  MInt end1[3] = {0, 0, 0};
  MInt step1[3] = {1, 1, 1};
  MInt order[3] = {0, 1, 2};
  MInt start2[3] = {0, 0, 0};
  MInt end2[3] = {0, 0, 0};

  MInt blockId = m_grid->getBlockId(domainId());

  /////////////////////////////////////////////////////////////
  //////////////// CREATE OWN/PARTITION MAPS //////////////////
  /////////////////////////////////////////////////////////////
  // this is the map of the active cells (no ghostcells) in my own partition
  // shifted by the no of ghost-cells
  for(MInt i = 0; i < nDim; i++) {
    start1[i] = m_grid->getMyOffset(nDim - 1 - i) + m_noGhostLayers; // shifted by the number of ghost layers
    start2[i] = m_noGhostLayers;   // shifted by the number of ghost layers
    end1[i] = start1[i] + m_grid->getMyActivePoints(nDim - 1 - i) - 1;
    end2[i] = start2[i] + m_grid->getMyActivePoints(nDim - 1 - i) - 1;
  }

  m_myMapWithoutGC = make_unique<StructuredWindowMap<nDim>>();
  mapCreate(blockId, start1, end1, step1, domainId(), start2, end2, step1, order, -1, m_myMapWithoutGC);

  // this is the map of all the cells  (including ghost-cells) in my own partition
  for(MInt i = 0; i < nDim; i++) {
    if(m_grid->getMyOffset(nDim - 1 - i) != 0) {
      start1[i] = m_grid->getMyOffset(nDim - 1 - i);
      end1[i] = start1[i] + m_grid->getMyActivePoints(nDim - 1 - i) - 1 + 2 * m_noGhostLayers;
    } else {
      start1[i] = 0;
      end1[i] = start1[i] + m_grid->getMyActivePoints(nDim - 1 - i) - 1 + 2 * m_noGhostLayers;
    }
    start2[i] = 0;
    end2[i] = m_grid->getMyActivePoints(nDim - 1 - i) - 1 + 2 * m_noGhostLayers;
  }
  m_myMapWithGC = make_unique<StructuredWindowMap<nDim>>();
  mapCreate(blockId, start1, end1, step1, domainId(), start2, end2, step1, order, -1, m_myMapWithGC);


  // maps of all partitions WITHOUT the ghostcells
  unique_ptr<StructuredWindowMap<nDim>> localMapDummy = nullptr;
  for(MInt j = 0; j < noDomains(); j++) {
    localMapDummy = make_unique<StructuredWindowMap<nDim>>();

    blockId = m_grid->getBlockId(j);
    for(MInt i = 0; i < nDim; i++) {
      start1[i] = m_grid->getOffset(j, nDim - 1 - i) + m_noGhostLayers;
      start2[i] = 0;
      end1[i] = start1[i] + m_grid->getActivePoints(j, nDim - 1 - i) - 1;
      end2[i] = m_grid->getActivePoints(j, nDim - 1 - i) - 1;
    }
    mapCreate(blockId, start1, end1, step1, j, start2, end2, step1, order, -6000, localMapDummy);
    m_partitionMapsWithoutGC.push_back(std::move(localMapDummy));
  }

  // maps of all partitions WITH the ghostcells
  for(MInt j = 0; j < noDomains(); j++) {
    localMapDummy = make_unique<StructuredWindowMap<nDim>>();

    blockId = m_grid->getBlockId(j);

    for(MInt i = 0; i < nDim; i++) {
      if(m_grid->getOffset(j, nDim - 1 - i) != 0) {
        start1[i] = m_grid->getOffset(j, nDim - 1 - i);
        end1[i] = start1[i] + m_grid->getActivePoints(j, nDim - 1 - i) - 1 + 2 * m_noGhostLayers;
      } else {
        start1[i] = 0;
        end1[i] = start1[i] + m_grid->getActivePoints(j, nDim - 1 - i) - 1 + 2 * m_noGhostLayers;
      }

      start2[i] = 0;
      end2[i] = m_grid->getActivePoints(j, nDim - 1 - i) - 1 + 4;
    }
    mapCreate(blockId, start1, end1, step1, j, start2, end2, step1, order, -6000, localMapDummy);
    m_partitionMapsWithGC.push_back(std::move(localMapDummy));
  }

  //////////////////////////////////////////////////////////////////////////////////
  ////////////////// USE MAPS TO CHECK FOR OVERLAPPING PARTS ///////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  // check overlapping of own partition with other partition and put into SND maps
  localMapDummy = make_unique<StructuredWindowMap<nDim>>();
  for(MInt i = 0; i < noDomains(); i++) {
    // skip if comparison to itself
    if(domainId() == i) {
      continue;
    }
    mapCombine11(m_myMapWithoutGC, m_partitionMapsWithGC[i], localMapDummy);
    MBool test = false;
    test = mapCheck(localMapDummy);
    if(test == true) {
      localMapDummy->BC = -6000;
      sndMap.push_back(std::move(localMapDummy));
      localMapDummy = make_unique<StructuredWindowMap<nDim>>();
    }
  }

  // check overlapping of own partition with other partition and put into RCV maps
  for(MInt i = 0; i < noDomains(); i++) {
    // skip if comparison to itself
    if(domainId() == i) {
      continue;
    }
    mapCombine11(m_myMapWithGC, m_partitionMapsWithoutGC[i], localMapDummy);
    MBool test = false;
    test = mapCheck(localMapDummy);
    if(test == true) {
      localMapDummy->BC = -6000;
      rcvMap.push_back(std::move(localMapDummy));
      localMapDummy = make_unique<StructuredWindowMap<nDim>>();
    }
  }

  // shift globalBC maps by ghost-cells
  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      globalStructuredBndryCndMaps[i]->start1[dim] += m_noGhostLayers;
      globalStructuredBndryCndMaps[i]->end1[dim] += m_noGhostLayers;
      if((globalStructuredBndryCndMaps[i]->BC <= -6000 && globalStructuredBndryCndMaps[i]->BC > -6010)
         || (globalStructuredBndryCndMaps[i]->BC >= 4000 && globalStructuredBndryCndMaps[i]->BC < 5000)
         || globalStructuredBndryCndMaps[i]->BC == 6011) {
        //      if(globalStructuredBndryCndMaps[i]->BC==6000 ||
        //      (globalStructuredBndryCndMaps[i]->BC>=4000&&globalStructuredBndryCndMaps[i]->BC<5000) ||
        //      globalStructuredBndryCndMaps[i]->BC==6001 ) { //6000er-change
        // for these special cases we also need to change the opposite side
        globalStructuredBndryCndMaps[i]->start2[dim] += m_noGhostLayers;
        globalStructuredBndryCndMaps[i]->end2[dim] += m_noGhostLayers;
      }
    }
  }

  // go over all the global boundary maps and check if there is any overlapping part with the local partition and store
  // it in localStructuredBndryCndMaps
  //==>find all the local boundaries of the partition (including the periodic part)

  // check for local boundary maps
  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    mapCombine11(m_myMapWithoutGC, globalStructuredBndryCndMaps[i], localMapDummy);
    MBool test = false;
    test = mapCheck(localMapDummy);

    if(test == true) {
      localMapDummy->BC = globalStructuredBndryCndMaps[i]->BC;
      localMapDummy->Nstar = globalStructuredBndryCndMaps[i]->Nstar;
      localMapDummy->SingularId = globalStructuredBndryCndMaps[i]->SingularId;
      localMapDummy->dc1 = globalStructuredBndryCndMaps[i]->dc1;
      localMapDummy->dc2 = globalStructuredBndryCndMaps[i]->dc2;
      localMapDummy->Id1 = globalStructuredBndryCndMaps[i]->Id1;
      localStructuredBndryCndMaps.push_back(std::move(localMapDummy));
      localMapDummy = make_unique<StructuredWindowMap<nDim>>();
    }
  }

  for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
    for(MInt i = 0; i < nDim; i++) {
      m_partitionMapsWithoutGC[j]->start2[i] += m_noGhostLayers;
      m_partitionMapsWithoutGC[j]->end2[i] += m_noGhostLayers;
    }
  }

  // check overlapping with bc 6000 (multiblock)
  vector<unique_ptr<StructuredWindowMap<nDim>>> SndMapBC6000;

  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC <= -6000
       && globalStructuredBndryCndMaps[i]->BC > -6010) { //(globalStructuredBndryCndMaps[i]->BC==6000) { //6000er-change
      for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
        // skip own domain, but only if 6000 does not connect the same domain
        if(globalStructuredBndryCndMaps[i]->Id1 != globalStructuredBndryCndMaps[i]->Id2) {
          if(domainId() == m_partitionMapsWithoutGC[j]->Id2) {
            continue;
          }
        }

        MInt myBlockId = m_grid->getBlockId(domainId());
        if(globalStructuredBndryCndMaps[i]->Id2 == myBlockId) {
          mapCombine11(globalStructuredBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
          MBool test = false;
          test = mapCheck(localMapDummy);

          if(test == true) {
            localMapDummy->BC = globalStructuredBndryCndMaps[i]->BC;
            localMapDummy->Nstar = globalStructuredBndryCndMaps[i]->Nstar;
            localMapDummy->SingularId = globalStructuredBndryCndMaps[i]->SingularId;
            localMapDummy->dc1 = globalStructuredBndryCndMaps[i]->dc1;
            localMapDummy->dc2 = globalStructuredBndryCndMaps[i]->dc2;

            SndMapBC6000.push_back(std::move(localMapDummy));
            localMapDummy = make_unique<StructuredWindowMap<nDim>>();
          }
        }
      }
    }
  }


  // now check overlapping with bc 4xxx (periodic bc)
  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC >= 4000 && globalStructuredBndryCndMaps[i]->BC < 5000) {
      for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
        // self communication is allowed
        mapCombine11(globalStructuredBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        MBool test = false;
        test = mapCheck(localMapDummy);

        if(test == true) {
          localMapDummy->BC = globalStructuredBndryCndMaps[i]->BC;
          localMapDummy->Nstar = globalStructuredBndryCndMaps[i]->Nstar;
          localMapDummy->SingularId = globalStructuredBndryCndMaps[i]->SingularId;
          localMapDummy->dc1 = globalStructuredBndryCndMaps[i]->dc1;
          localMapDummy->dc2 = globalStructuredBndryCndMaps[i]->dc2;
          SndMapBC6000.push_back(std::move(localMapDummy));
          localMapDummy = make_unique<StructuredWindowMap<nDim>>();
        }
      }
    }
  }

  // now 6000 and 4xxx communication bc are both in SndMapBC6000

  ///////////////////////////////////////////////
  /////// SINGULARITY COMMUNICATION  ////////////
  ///////////////////////////////////////////////
  // special treatment of singularities
  MInt singularPoint = (MInt)singularwindow.size();
  for(MInt i = 0; i < singularPoint; ++i) {
    for(MInt dim = 0; dim < nDim; ++dim) {
      singularwindow[i]->start1[dim] += m_noGhostLayers;
      singularwindow[i]->end1[dim] += m_noGhostLayers;
      // need also to change the opposite side
      singularwindow[i]->start2[dim] += m_noGhostLayers;
      singularwindow[i]->end2[dim] += m_noGhostLayers;
    }
  }


  for(MInt i = 0; i < singularPoint; ++i) {
    MInt singularcount = 0, limit = singularwindow[i]->Nstar - 3, cnt = 0;

    std::set<MInt> BCs;
    if(singularwindow[i]->BC < 4400 || singularwindow[i]->BC >= 4410) {
      // all 6000er singularities have BC=-6000
      singularwindow[i]->BC = -6000;
      for(MInt t = 0; t < singularwindow[i]->Nstar; ++t) {
        if(singularwindow[i]->BCsingular[t] != 0) BCs.insert(-singularwindow[i]->BCsingular[t]);
      }
      singularwindow[i]->BCsingular[0] = 0;
      singularwindow[i]->BCsingular[1] = 0;
    }

    for(MInt j = 0; j < (MInt)globalStructuredBndryCndMaps.size(); j++) {
      if((globalStructuredBndryCndMaps[j]->Nstar == 5 || globalStructuredBndryCndMaps[j]->Nstar == 6)
         && (globalStructuredBndryCndMaps[j]->BC <= -6000 && globalStructuredBndryCndMaps[j]->BC > -6010)
         && globalStructuredBndryCndMaps[j]->Nstar == singularwindow[i]->Nstar) {
        //      if((globalStructuredBndryCndMaps[j]->Nstar==5||globalStructuredBndryCndMaps[j]->Nstar==6)&&globalStructuredBndryCndMaps[j]->BC==6000&&globalStructuredBndryCndMaps[j]->Nstar==singularwindow[i]->Nstar)
        //      { //6000er-change

        // TODO_SS labels:FV I am not sure if the following is always true
        if(singularwindow[i]->BC >= 4400 && singularwindow[i]->BC < 4410) mTerm(1, "ERROR 0");

        MInt labell = 0;
        labell = mapCompare11(singularwindow[i], globalStructuredBndryCndMaps[j]);
        if(labell) {
          singularwindow[i]->SingularBlockId[singularcount] = globalStructuredBndryCndMaps[j]->Id2;
          if(BCs.size() == 1) {
            if(*BCs.cbegin() != globalStructuredBndryCndMaps[j]->BC) mTerm(1, "ERROR 1");
            singularwindow[i]->BCsingular[singularcount + 2] = *BCs.cbegin();
          } else {
            setBCsingular(singularwindow[i], globalStructuredBndryCndMaps[j], singularcount);
          }
          singularcount++;
        }
      }
      if(singularcount > limit)
        cout << "Important!! " << globalStructuredBndryCndMaps[j]->Nstar
             << " star check error!!!!!!!!!!!!!!!!!! check it !!!!!" << endl;

      // set BCsingular[0] and BCsingular[1] in arbritary order
      if(globalStructuredBndryCndMaps[j]->Nstar == -1
         && (globalStructuredBndryCndMaps[j]->BC <= -6000 && globalStructuredBndryCndMaps[j]->BC > -6010)) {
        localMapDummy = make_unique<StructuredWindowMap<nDim>>();
        mapCombine11(singularwindow[i], globalStructuredBndryCndMaps[j], localMapDummy);
        MBool test = false;
        test = mapCheck(localMapDummy);
        if(test) {
          singularwindow[i]->BCsingular[cnt++] = globalStructuredBndryCndMaps[j]->BC;
        }
      }
    }

    if(cnt >= 3) mTerm(1, "");
  }

  // if(domainId()==0)
  //   {
  //     for(MInt i=0; i<singularPoint; ++i) {
  //     if (singularwindow[i]->BC>=6000 && singularwindow[i]->BC<6010) //(singularwindow[i]->BC==6000) //6000er-change
  //       cout<< "singularmaps  Id1: " << singularwindow[i]->Id1 << " Id2: " << singularwindow[i]->Id2 << " BC: " <<
  //       singularwindow[i]->BC<<"  "<<  singularwindow[i]->start1[0] <<"-" <<  singularwindow[i]->end1[0]<< " "<<
  //       singularwindow[i]->start1[1] <<"-" <<  singularwindow[i]->end1[1]<< " "<<  singularwindow[i]->start1[2] <<"-"
  //       <<  singularwindow[i]->end1[2]<< " ID: " <<singularwindow[i]->SingularBlockId[0]<<"
  //       "<<singularwindow[i]->SingularBlockId[1]<< endl;

  //     }
  //   }

  /*  if (domainId()==0) {
      cout << "MY TESTGLOBAL" << endl;
      for (MInt i = 0; i < (signed)globalStructuredBndryCndMaps.size(); ++i) {
        if (abs(globalStructuredBndryCndMaps[i]->BC)>=6000 && abs(globalStructuredBndryCndMaps[i]->BC<6010))
          mapPrint(globalStructuredBndryCndMaps[i]);
      }
      cout << "MY TEST SINGULAR" << endl;
      for (MInt i = 0; i < (signed)singularwindow.size(); ++i) {
        mapPrint(singularwindow[i]);
      }
    }*/

  for(MInt i = 0; i < singularPoint; ++i) {
    mapCombine11(m_myMapWithoutGC, singularwindow[i], localMapDummy);
    MBool test = false;
    test = mapCheck(localMapDummy);
    localMapDummy->Nstar = singularwindow[i]->Nstar;
    localMapDummy->SingularId = singularwindow[i]->SingularId;
    for(MInt n = 0; n < 6; ++n) {
      localMapDummy->BCsingular[n] = singularwindow[i]->BCsingular[n];
    }

    //     // To distinguish between communication bcs and normal ones, do the following, negate the communication bcs
    //    const MInt BC = (singularwindow[i]->BC>=6000 && singularwindow[i]->BC<6010) ? -singularwindow[i]->BC :
    //    singularwindow[i]->BC;
    localMapDummy->BC = singularwindow[i]->BC;

    localMapDummy->SingularBlockId[0] = singularwindow[i]->SingularBlockId[0];
    localMapDummy->SingularBlockId[1] = singularwindow[i]->SingularBlockId[1];
    localMapDummy->SingularBlockId[2] = singularwindow[i]->SingularBlockId[2];
    localMapDummy->SingularBlockId[3] = singularwindow[i]->SingularBlockId[3];

    if(test == true) {
      localSingularMap.push_back(std::move(localMapDummy));
      localMapDummy = make_unique<StructuredWindowMap<nDim>>();
    }
  }

  *hasSingularity = (MInt)localSingularMap.size();
  singularPoint = (MInt)localSingularMap.size();

  for(MInt i = 0; i < (MInt)localSingularMap.size(); ++i) {
    memcpy(singularity[i].start, localSingularMap[i]->start1, nDim * sizeof(MInt));
    memcpy(singularity[i].end, localSingularMap[i]->end1, nDim * sizeof(MInt));
    singularity[i].Nstar = localSingularMap[i]->Nstar;
    singularity[i].count = localSingularMap[i]->SingularId;
    singularity[i].BC = localSingularMap[i]->BC;
    for(MInt n = 0; n < 6; ++n) {
      singularity[i].BCsingular[n] = localSingularMap[i]->BCsingular[n];
    }
    singularity[i].SingularBlockId[0] = localSingularMap[i]->SingularBlockId[0];
    singularity[i].SingularBlockId[1] = localSingularMap[i]->SingularBlockId[1];
    singularity[i].SingularBlockId[2] = localSingularMap[i]->SingularBlockId[2];
    singularity[i].SingularBlockId[3] = localSingularMap[i]->SingularBlockId[3];

    //!!!!nimportant!!!!
    // for local singular map
    // BC6000: line only (no point and face) BC4000: point only.

    if(localSingularMap[i]->BC >= 4400 && localSingularMap[i]->BC < 4410) {
      localSingularMap[i]->face = localSingularMap[i]->BC - 4401;
    } else {
      for(MInt k = 0; k < nDim; ++k) {
        if(localSingularMap[i]->start1[k] != localSingularMap[i]->end1[k]) {
          localSingularMap[i]->face = k * 2;
        }
      }
    }

    MInt displace[nDim], dimcount = 0;
    MInt dimN = -1;
    MInt dimT1, dimT2;
    for(MInt k = 0; k < nDim; ++k) {
      if(singularity[i].start[k] == singularity[i].end[k]) {
        if(singularity[i].start[k] == 2) {
          displace[k] = -1;
        } else {
          displace[k] = 1;
        }
        dimcount++;
      } else {
        displace[k] = 0;
        dimN = k;
      }
    }

    IF_CONSTEXPR(nDim == 2) {
      ASSERT(dimcount == 2, "Impossible!");
      dimN = 0;
    }

    if(dimcount == 3) {
      if(localSingularMap[i]->face != -1) {
        dimN = localSingularMap[i]->face / 2;
        displace[dimN] = 0;
        // TODO_SS labels:FV take care of the next if-clause
        if(singularity[i].BC >= 6000
           && singularity[i].BC < 6010) { //(singularity[i].BC==6000) { //6000er-change ; what is the next line supposed
                                          // to do? 6003 is nowhere used!
          singularity[i].BC = 6003;
        }
      } else {
        cout << "ERROR!!! point singular communication cannot decide the direction!! check the singularity!!!" << endl;
        exit(0);
      }
    }

    if(dimN == 0 && nDim != 2) {
      dimT1 = 1;
      dimT2 = 2;
    } else if(dimN == 1) {
      dimT1 = 0;
      dimT2 = 2;
    } else {
      dimT1 = 0;
      dimT2 = 1;
    } // this else clause is always true in 2D

    switch(singularity[i].Nstar) {
      case 3: {
        // labels:FV bugs possibly exist
        singularity[i].displacement[0][dimN] = 0;
        singularity[i].displacement[0][dimT1] = displace[dimT1];
        singularity[i].displacement[0][dimT2] = 0;

        singularity[i].displacement[1][dimN] = 0;
        singularity[i].displacement[1][dimT1] = 0;
        singularity[i].displacement[1][dimT2] = displace[dimT2];

        break;
      }
      case 5: {
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // labels:FV bugs  exist need to be fixed
        singularity[i].displacement[0][dimN] = 0;
        singularity[i].displacement[0][dimT1] = displace[dimT1];
        singularity[i].displacement[0][dimT2] = 0;

        singularity[i].displacement[1][dimN] = 0;
        singularity[i].displacement[1][dimT1] = 0;
        singularity[i].displacement[1][dimT2] = displace[dimT2];

        singularity[i].displacement[2][dimN] = 0;
        singularity[i].displacement[2][dimT1] = displace[dimT1];
        singularity[i].displacement[2][dimT2] = displace[dimT2];

        singularity[i].displacement[3][dimN] = 0;
        singularity[i].displacement[3][dimT1] = 2 * displace[dimT1];
        singularity[i].displacement[3][dimT2] = 2 * displace[dimT2];


        // communication map change start from 2 to 3
        singularity[i].count = 2;

        break;
      }
      case 6: {
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // labels:FV bugs  exist need to be fixed
        singularity[i].displacement[0][dimN] = 0;
        singularity[i].displacement[0][dimT1] = displace[dimT1];
        singularity[i].displacement[0][dimT2] = 0;

        singularity[i].displacement[1][dimN] = 0;
        singularity[i].displacement[1][dimT1] = 0;
        singularity[i].displacement[1][dimT2] = displace[dimT2];

        singularity[i].displacement[2][dimN] = 0;
        singularity[i].displacement[2][dimT1] = displace[dimT1];
        singularity[i].displacement[2][dimT2] = displace[dimT2];

        singularity[i].displacement[3][dimN] = 0;
        singularity[i].displacement[3][dimT1] = 2 * displace[dimT1];
        singularity[i].displacement[3][dimT2] = 2 * displace[dimT2];

        singularity[i].displacement[4][dimN] = 0;
        singularity[i].displacement[4][dimT1] = displace[dimT1];
        singularity[i].displacement[4][dimT2] = 2 * displace[dimT2];


        // communication map change start from 2 to 3
        singularity[i].count = 2;

        break;
      }
      default: {
        break;
      }
    }
  }

  // expand the singularmap
  for(MInt i = 0; i < (MInt)localSingularMap.size(); ++i) {
    if(localSingularMap[i]->BC <= -6000
       && localSingularMap[i]->BC > -6010) //(localSingularMap[i]->BC==6000) //6000er-change
    {
      for(MInt k = 0; k < nDim; ++k) {
        if(localSingularMap[i]->start1[k] != localSingularMap[i]->end1[k]) {
          // if ( m_myMapWithoutGC->start2[k]==localSingularMap[i]->start1[k])
          {
            localSingularMap[i]->start1[k] -= 2;
            localSingularMap[i]->start2[k] -= 2;
          }
          // if ( m_myMapWithoutGC->end2[k]==localSingularMap[i]->end1[k])
          {
            localSingularMap[i]->end1[k] += 2;
            localSingularMap[i]->end2[k] += 2;
          }
        } // if !=
      }   // for k
    }     // if 6000
  }

  /*  for (MInt i = 0; i < (MInt)localSingularMap.size(); ++i) {
      cout << "SINGULARITY CHECK : d=" << domainId() << " Nstar=" << singularity[i].Nstar
           << " singularBlockId=" << singularity[i].SingularBlockId[0] << "|" <<
    singularity[i].SingularBlockId[1] << "|" << singularity[i].SingularBlockId[2] << "|" <<
    singularity[i].SingularBlockId[3]
           << " BC=" << singularity[i].BC << " BC_singularity=" << singularity[i].BCsingular[0] << "|" <<
    singularity[i].BCsingular[1] << "|" << singularity[i].BCsingular[2] << "|" << singularity[i].BCsingular[3]
           << " start=" << singularity[i].start[0] << "|" << singularity[i].start[1]
           << " end=" << singularity[i].end[0] << "|" << singularity[i].end[1] << endl;
    }*/

  // first put multiblock communication maps into addCommBC,
  // afterwards also add the periodic communcation maps

  vector<unique_ptr<StructuredWindowMap<nDim>>> periodicBC;
  vector<unique_ptr<StructuredWindowMap<nDim>>> addCommBC;
  MPI_Barrier(m_StructuredComm, AT_);

  // auto& localMap = begin(localStructuredBndryCndMaps);
  for(auto& localMap : localStructuredBndryCndMaps) {
    const MInt BC = (localMap->BC <= -6000 && localMap->BC > -6010) ? 6000 : localMap->BC; // 6000er-change
    switch(BC) {
      case 4011: // periodic rotational 1
      case 4012: // periodic rotational 2
      case 4401: // periodic surface 1
      case 4402: // periodic surface 2
      case 4403:
      case 4404:
      case 4405:
      case 4406: {
        break;
      }
      case 6000: // additional block communication (interblock communication)
      {
        // save the additional Communication
        unique_ptr<StructuredWindowMap<nDim>> copy = make_unique<StructuredWindowMap<nDim>>();
        mapCpy(localMap, copy);
        addCommBC.push_back(std::move(copy));
        // localMap = localStructuredBndryCndMaps.erase(localMap);
        break;
      }
      case 2097: { // this is for the plenum inlet where we need the surface informations
        unique_ptr<StructuredWindowMap<nDim>> copy1 = make_unique<StructuredWindowMap<nDim>>();
        mapCpy(localMap, copy1);
        physicalBCMap.push_back(std::move(copy1));

        unique_ptr<StructuredWindowMap<nDim>> copy2 = make_unique<StructuredWindowMap<nDim>>();
        mapCpy(localMap, copy2);
        plenumInletSurfaceIndices.push_back(std::move(copy2)); // needed for identifiaction of areas etc.

        // localMap = localStructuredBndryCndMaps.erase(localMap);
        break;
      }
      case 2401:
      case 2402: // this is for the channel flow
      {
        // we need a copy of the map to find the channel In/Outflow parti-
        // cipating processors.
        unique_ptr<StructuredWindowMap<nDim>> copy1 = make_unique<StructuredWindowMap<nDim>>();
        mapCpy(localMap, copy1);
        physicalBCMap.push_back(std::move(copy1));

        unique_ptr<StructuredWindowMap<nDim>> copy2 = make_unique<StructuredWindowMap<nDim>>();
        mapCpy(localMap, copy2);
        channelSurfaceIndices.push_back(std::move(copy2)); // needed for identifiaction of aereas etc.

        // localMap = localStructuredBndryCndMaps.erase(localMap);
        break;
      }
      default: {
        unique_ptr<StructuredWindowMap<nDim>> copy = make_unique<StructuredWindowMap<nDim>>();
        mapCpy(localMap, copy);
        physicalBCMap.push_back(std::move(copy));

        // localMap = localStructuredBndryCndMaps.erase(localMap);
        break;
      }
    }
  }

  // localMap = begin(localStructuredBndryCndMaps);
  // while (localMap != end(localStructuredBndryCndMaps)) {
  for(auto& localMap : localStructuredBndryCndMaps) {
    switch(localMap->BC) {
      case 4011: // periodic rotational 1
      case 4012: // periodic rotational 2
      case 4401: // periodic surface 1
      case 4402: // periodic surface 2
      case 4403:
      case 4404:
      case 4405:
      case 4406: {
        // save the periodic boundary condition to addcommBC, same position
        unique_ptr<StructuredWindowMap<nDim>> copy = make_unique<StructuredWindowMap<nDim>>();
        mapCpy(localMap, copy);
        addCommBC.push_back(std::move(copy));
        break;
      }
      default: {
        // do nothing
        break;
      }
    }
  }

  ///////////////////////////////////////////
  //////// PLENUM BC 2097 ///////////////////
  ///////////////////////////////////////////
  vector<MInt> plenumPartitions;
  MInt count = 1;
  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC == 2097) {
      // found a face containing the plenum inlet
      count++;
      // loop over all domains without ghost cells
      for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
        // check if boundary condition is contained
        mapCombine11(globalStructuredBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        MBool test = false;
        test = mapCheck(localMapDummy);                                 // true if a match is found
        if(test) {                                                      // match is found
          plenumPartitions.push_back(m_partitionMapsWithoutGC[j]->Id2); // store the domainId
        }
      }
    }
  }

  // get the ranks for building up a communicator
  vector<MInt> plenumInflow;
  if(plenumPartitions.size() != 0) {
    for(auto& plenumPartition : plenumPartitions) {
      // for(MInt i=0; i< plenumPartitions.size(); i++){
      for(MInt j = 0; j < noDomains(); j++) {
        if(plenumPartition == j) {
          plenumInflow.push_back(j);
        }
      }
    }

    m_log << "Number of input partitions plenum Inlet " << plenumInflow.size() << endl;
    MPI_Group groupIn = MPI_GROUP_NULL;
    MPI_Group* newgroupIn = new MPI_Group;
    MPI_Comm_group(m_StructuredComm, &groupIn, AT_, "groupIn");
    MPI_Group_incl(groupIn, (MInt)plenumInflow.size(), plenumInflow.data(), newgroupIn, AT_);
    MPI_Comm_create(m_StructuredComm, newgroupIn[0], plenumComm, AT_, "plenumComm");

    if(domainId() == plenumInflow[0]) {
      MPI_Comm_rank(plenumComm[0], &plenumRoot[0]);
    }
    MPI_Bcast(&plenumRoot[0], 1, MPI_INT, plenumInflow[0], m_StructuredComm, AT_, "plenumRoots[0]");
  }

  ///////////////////////////////////////////
  //////// CHANNEL BC ///////////////////////
  ///////////////////////////////////////////

  // creating the communicators
  vector<MInt> channelface;
  vector<MInt> channelpartitions;
  MInt anzahl = 1;
  MInt counterIn = 0;
  MInt counterOut = 0;
  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC == 2401 || globalStructuredBndryCndMaps[i]->BC == 2402) {
      // found a face containing the channel in/out flow
      // now check the partition it belongs to!!
      // cout << "found " << anzahl << endl;
      anzahl++;
      for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11(globalStructuredBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        MBool test = false;
        test = mapCheck(localMapDummy);
        if(test == true) {
          channelpartitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
          // partition contains part of the boundary
          // check the face

          // now check the face
          if(globalStructuredBndryCndMaps[i]->BC == 2401) {
            channelface.push_back(globalStructuredBndryCndMaps[i]->BC);
            counterIn++;
          } else {
            channelface.push_back(globalStructuredBndryCndMaps[i]->BC);
            counterOut++;
          }
        }
      }
    }
  }

  // now assign the ranks
  MIntScratchSpace channelInflow(counterIn, AT_, "channelInflow");
  MIntScratchSpace channelOutflow(counterOut, AT_, "channelInflow");

  if(channelpartitions.size() != 0) {
    MInt counterI = 0, counterO = 0;
    for(MInt i = 0; i < (MInt)channelpartitions.size(); i++) {
      for(MInt j = 0; j < noDomains(); j++) {
        if(channelpartitions[i] == j) {
          if(channelface[i] == 2401) {
            channelInflow[counterI] = j;
            counterI++;
          }
          if(channelface[i] == 2402) {
            channelOutflow[counterO] = j;
            counterO++;
          }
        }
      }
    }


    m_log << "ChannelIn: " << counterI << " ChannelOut: " << counterOut << endl;

    // create a communicator for the in- and the outflow
    // averaging process
    MPI_Group groupIn = MPI_GROUP_NULL;
    MPI_Group newgroupIn = MPI_GROUP_NULL;
    MPI_Group groupOut = MPI_GROUP_NULL;
    MPI_Group newgroupOut = MPI_GROUP_NULL;

    MPI_Comm_group(m_StructuredComm, &groupIn, AT_, "groupIn");
    MPI_Comm_group(m_StructuredComm, &groupOut, AT_, "groupOut");

    MPI_Group_incl(groupIn, counterIn, channelInflow.begin(), &newgroupIn, AT_);
    MPI_Comm_create(m_StructuredComm, newgroupIn, channelIn, AT_, "channelIn");

    MPI_Group_incl(groupOut, counterOut, channelOutflow.begin(), &newgroupOut, AT_);
    MPI_Comm_create(m_StructuredComm, newgroupOut, channelOut, AT_, "channelOut");

    // finally we need to create a overall channel surface communicator to distribute
    // the averaged values p & T

    // create group of all in- and outflow participants
    // but pay attention to duplicates (one domain shares in- and outflow)
    MInt counterAll = 0;
    MInt* channelAll = new MInt[counterOut + counterIn];
    for(MInt i = 0; i < counterIn; i++) {
      channelAll[i] = channelInflow[i];
      counterAll++;
    }
    for(MInt i = 0; i < counterOut; i++) {
      MBool check = false;
      for(MInt j = 0; j < counterIn; ++j) {
        if(channelAll[j] == channelOutflow[i]) {
          check = true;
          cout << "ATTENTION: Skipping duplicate generation in channelAll" << endl;
          mTerm(1, "In some cases this causes undefined behavior");
          break;
        }
      }
      // if this domain is already in channelAll we must not insert it again
      if(!check) {
        channelAll[counterIn + i] = channelOutflow[i];
        counterAll++;
      }
    }
    MPI_Group groupAll = MPI_GROUP_NULL;
    MPI_Group newgroupAll = MPI_GROUP_NULL;
    MPI_Comm_group(m_StructuredComm, &groupAll, AT_, "groupAll");
    MPI_Group_incl(groupAll, (MInt)(counterAll), channelAll, &newgroupAll, AT_);

    MPI_Comm_create(m_StructuredComm, newgroupAll, channelWorld, AT_, "channelWorld");

    // for global exchange we need a root process
    if(domainId() == channelInflow[0]) {
      MPI_Comm_rank(*channelIn, &channelRoots[0]);
      MPI_Comm_rank(*channelWorld, &channelRoots[2]);
    }
    MPI_Barrier(m_StructuredComm, AT_);
    MPI_Bcast(&channelRoots[0], 1, MPI_INT, channelInflow[0], m_StructuredComm, AT_, "channelRoots[0]");
    MPI_Bcast(&channelRoots[2], 1, MPI_INT, channelInflow[0], m_StructuredComm, AT_, "channelRoots[2]");

    if(domainId() == channelOutflow[0]) {
      MPI_Comm_rank(*channelOut, &channelRoots[1]);
      MPI_Comm_rank(*channelWorld, &channelRoots[3]);
    }
    MPI_Barrier(m_StructuredComm, AT_);
    MPI_Bcast(&channelRoots[1], 1, MPI_INT, channelOutflow[0], m_StructuredComm, AT_, "channelRoots[1]");
    MPI_Bcast(&channelRoots[3], 1, MPI_INT, channelOutflow[0], m_StructuredComm, AT_, "channelRoots[3]");
  }

  //////////////////////////////////////////
  //////// PERIODIC ROTATION BC ////////////
  //////////////////////////////////////////

  vector<MInt> rotatingface;
  vector<MInt> rotatingpartitions;

  MInt noRotPartitions = 0;
  counterIn = 0;
  counterOut = 0;

  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC == 4001
       || globalStructuredBndryCndMaps[i]->BC == 4002) { //== containing rotation info
      // found a face containing the per. rotation BC ==> which partition?
      noRotPartitions++;
      for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11(globalStructuredBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        MBool test = false;
        test = mapCheck(localMapDummy);
        if(test == true) {
          rotatingpartitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
          // partition contains part of the boundary
          if(globalStructuredBndryCndMaps[i]->BC == 4001) {
            rotatingface.push_back(globalStructuredBndryCndMaps[i]->BC);
            counterIn++;
          } else {
            rotatingface.push_back(globalStructuredBndryCndMaps[i]->BC);
            counterOut++;
          }
        }
      }
    }
  }

  MInt* rotationOne = nullptr;
  MInt* rotationTwo = nullptr;
  if(rotatingpartitions.size() != 0) {
    MInt counterI = 0, counterO = 0;
    rotationOne = new MInt[counterIn];
    rotationTwo = new MInt[counterOut];
    for(MInt i = 0; i < (MInt)rotatingpartitions.size(); i++) {
      for(MInt j = 0; j < noDomains(); j++) {
        if(rotatingpartitions[i] == j) {
          if(rotatingface[i] == 4001) {
            rotationOne[counterI] = j;
            ++counterI;
          } else {
            rotationTwo[counterO] = j;
            counterO++;
          }
        }
      }
    }
    // check which group the processor is in
    perRotGroup = 0;
    for(MInt i = 0; i < counterI; ++i) {
      if(domainId() == rotationOne[i]) {
        perRotGroup = 1;
        break;
      }
    }
    for(MInt i = 0; i < counterO; ++i) {
      if(domainId() == rotationTwo[i]) {
        // Maybe domain is in both groups, then set to 3
        if(perRotGroup == 1)
          perRotGroup = 3;
        else
          perRotGroup = 2;

        break;
      }
    }

    // create a communicator for the two sides of the rotation BC
    MPI_Group groupOne = MPI_GROUP_NULL;
    MPI_Group* newgroupOne = new MPI_Group;
    MPI_Group groupTwo = MPI_GROUP_NULL;
    MPI_Group* newgroupTwo = new MPI_Group;

    MPI_Comm_group(m_StructuredComm, &groupOne, AT_, "groupOne");
    MPI_Comm_group(m_StructuredComm, &groupTwo, AT_, "groupTwo");

    MPI_Group_incl(groupOne, (MInt)counterIn, rotationOne, newgroupOne, AT_);

    MPI_Comm_create(m_StructuredComm, *newgroupOne, commPerRotOne, AT_, "commPerRotOne");

    MPI_Group_incl(groupTwo, (MInt)counterOut, rotationTwo, newgroupTwo, AT_);

    MPI_Comm_create(m_StructuredComm, *newgroupTwo, commPerRotTwo, AT_, "commPerRotTwo");

    // we also need one cross communicator containing all processors involved
    // in the periodic rotation boundary condition

    MInt counterAll = 0;
    MInt* rotationAll = new MInt[counterOut + counterIn];
    for(MInt i = 0; i < counterIn; ++i) {
      rotationAll[i] = rotationOne[i];
      counterAll++;
    }

    for(MInt i = 0; i < counterOut; ++i) {
      // check if this domain is already in rotationAll
      MBool check = false;
      for(MInt j = 0; j < counterIn; ++j) {
        if(rotationAll[j] == rotationTwo[i]) {
          check = true;
          cout << "ATTENTION: Skipping duplicate generation in rotationAll" << endl;
          break;
        }
      }

      // if this domain is already in rotationAll we must not insert it again
      if(!check) {
        rotationAll[counterIn + i] = rotationTwo[i];
        counterAll++;
      }
    }

    MPI_Group groupAll = MPI_GROUP_NULL;
    MPI_Group* newgroupAll = new MPI_Group;

    MPI_Comm_group(m_StructuredComm, &groupAll, AT_, "groupAll");
    MPI_Group_incl(groupAll, (MInt)(counterAll), rotationAll, newgroupAll, AT_);

    MPI_Comm_create(m_StructuredComm, *newgroupAll, commPerRotWorld, AT_, "commPerRotWorld");

    // fix the root processes for the communication
    if(domainId() == rotationOne[0]) {
      MPI_Comm_rank(commPerRotOne[0], &rotationRoots[0]);
      MPI_Comm_rank(commPerRotWorld[0], &rotationRoots[2]);
    }
    MPI_Barrier(m_StructuredComm, AT_);
    MPI_Bcast(&rotationRoots[0], 1, MPI_INT, rotationOne[0], m_StructuredComm, AT_, "rotationRoots[0]");
    MPI_Bcast(&rotationRoots[2], 1, MPI_INT, rotationOne[0], m_StructuredComm, AT_, "rotationRoots[2]");

    if(domainId() == rotationTwo[0]) {
      MPI_Comm_rank(commPerRotTwo[0], &rotationRoots[1]);
      MPI_Comm_rank(commPerRotWorld[0], &rotationRoots[3]);
    }
    MPI_Barrier(m_StructuredComm, AT_);
    MPI_Bcast(&rotationRoots[1], 1, MPI_INT, rotationTwo[0], m_StructuredComm, AT_, "rotationRoots[1]");
    MPI_Bcast(&rotationRoots[3], 1, MPI_INT, rotationTwo[0], m_StructuredComm, AT_, "rotationRoots[3]");
  }


  /////////////////////////////////////////////////////////
  ///////////////// STG BC Communication //////////////////
  /////////////////////////////////////////////////////////

  /*
    Create the MPI Group for the STG methods that only includes those
    domains that are located at the STG inflow.
  */

  vector<MInt> stgpartitions;

  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC == 7909) {
      for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11(globalStructuredBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        MBool test = false;
        test = mapCheck(localMapDummy);
        if(test == true) stgpartitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
      }
    }
  }

  MPI_Barrier(m_StructuredComm, AT_);

  if(stgpartitions.size() != 0) {
    m_log << "STG domains: ";
    MInt counterstg = 0;
    MIntScratchSpace stgranks(stgpartitions.size(), "stgranks", AT_);

    for(MInt i = 0; i < (MInt)stgpartitions.size(); i++) {
      for(MInt j = 0; j < noDomains(); j++) {
        if(stgpartitions[i] == j) {
          stgranks[counterstg] = j;
          counterstg++;
          m_log << j << ", ";
        }
      }
    }
    m_log << endl;
    m_log << "Total number of STG domains: " << counterstg << endl;

    MPI_Barrier(m_StructuredComm, AT_);

    MPI_Group groupStg, newgroupStg;
    MInt stgcommsize = (MInt)(stgpartitions.size());

    MPI_Comm_group(m_StructuredComm, &groupStg, AT_, "groupStg");
    MPI_Group_incl(groupStg, stgcommsize, stgranks.begin(), &newgroupStg, AT_);

    MPI_Comm_create(m_StructuredComm, newgroupStg, commStg, AT_, "commStg");

    if(domainId() == stgranks[0]) {
      MPI_Comm_rank(*commStg, commStgRoot);
      MPI_Comm_rank(m_StructuredComm, commStgRootGlobal);
    }

    MPI_Barrier(m_StructuredComm, AT_);
    MPI_Bcast(commStgRoot, 1, MPI_INT, stgranks[0], m_StructuredComm, AT_, "commStgRoot");
    MPI_Bcast(commStgRootGlobal, 1, MPI_INT, stgranks[0], m_StructuredComm, AT_, "commStgRootGlobal");
  }


  //////////////////////////////////////////////////////
  ////////////////BC 2600 Communication////////////////
  /////////////////////////////////////////////////////
  vector<MInt> bc2600partitions;

  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC == 2600) {
      for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11(globalStructuredBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        MBool test = false;
        test = mapCheck(localMapDummy);
        if(test == true) bc2600partitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
      }
    }
  }

  MPI_Barrier(m_StructuredComm, AT_);

  MInt* bc2600ranks;
  if(bc2600partitions.size() != 0) {
    m_log << "bc2600 domains: ";
    MInt counterBC2600 = 0;
    bc2600ranks = new MInt[mMax(MInt(bc2600partitions.size()), 0)];
    for(MInt i = 0; i < (MInt)bc2600partitions.size(); i++) {
      for(MInt j = 0; j < noDomains(); j++) {
        if(bc2600partitions[i] == j) {
          bc2600ranks[counterBC2600] = j;
          counterBC2600++;
          m_log << j << ", ";
        }
      }
    }
    m_log << endl;
    m_log << "Total number of BC2600 domains: " << counterBC2600 << endl;

    MPI_Barrier(m_StructuredComm, AT_);

    MPI_Group groupBC2600, newgroupBC2600;
    MInt bc2600commsize = (MInt)(bc2600partitions.size());

    MPI_Comm_group(m_StructuredComm, &groupBC2600, AT_, "groupBC2600");
    MPI_Group_incl(groupBC2600, bc2600commsize, bc2600ranks, &newgroupBC2600, AT_);
    MPI_Comm_create(m_StructuredComm, newgroupBC2600, commBC2600, AT_, "commBC2600");

    if(domainId() == bc2600ranks[0]) {
      MPI_Comm_rank(commBC2600[0], &commBC2600Root[0]);
      MPI_Comm_rank(m_StructuredComm, &commBC2600RootGlobal[0]);
    }

    MPI_Barrier(m_StructuredComm, AT_);
    MPI_Bcast(commBC2600Root, 1, MPI_INT, bc2600ranks[0], m_StructuredComm, AT_, "commBC2600Root");
    MPI_Bcast(commBC2600RootGlobal, 1, MPI_INT, bc2600ranks[0], m_StructuredComm, AT_, "commBC2600RootGlobal");
  }


  /////////////////////////////////////////////////////////
  ///////////////// Rescaling BC Communication ////////////
  /////////////////////////////////////////////////////////
  vector<MInt> rescalingCommGrPartitions;
  for(MInt i = 0; i < (MInt)globalStructuredBndryCndMaps.size(); i++) {
    if(globalStructuredBndryCndMaps[i]->BC == 2500 || globalStructuredBndryCndMaps[i]->BC == 2501) {
      for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11(globalStructuredBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        MBool test = false;
        test = mapCheck(localMapDummy);
        if(test == true) {
          rescalingCommGrPartitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
        }
      }
    }
  }

  MPI_Barrier(m_StructuredComm, AT_);

  MInt* rescalingCommGrRanks;

  if(rescalingCommGrPartitions.size() != 0) {
    MInt counterCommGrRescaling = 0;
    rescalingCommGrRanks = new MInt[mMax(MInt(rescalingCommGrPartitions.size()), 0)];
    for(MInt i = 0; i < (MInt)rescalingCommGrPartitions.size(); i++) {
      for(MInt j = 0; j < noDomains(); j++) {
        if(rescalingCommGrPartitions[i] == j) {
          rescalingCommGrRanks[counterCommGrRescaling] = j;
          counterCommGrRescaling++;
        }
      }
    }

    MPI_Barrier(m_StructuredComm, AT_);
    MPI_Group groupRescalingCommGr, newgroupRescalingCommGr;
    MInt rescalingCommGrCommsize = (MInt)(rescalingCommGrPartitions.size());

    MPI_Comm_group(m_StructuredComm, &groupRescalingCommGr, AT_, "groupRescalingCommGr");
    MPI_Group_incl(groupRescalingCommGr, rescalingCommGrCommsize, rescalingCommGrRanks, &newgroupRescalingCommGr, AT_);

    MPI_Comm_create(m_StructuredComm, newgroupRescalingCommGr, rescalingCommGrComm, AT_, "rescalingCommGrComm");

    if(domainId() == rescalingCommGrRanks[0]) {
      MPI_Comm_rank(rescalingCommGrComm[0], &rescalingCommGrRoot[0]);
      MPI_Comm_rank(m_StructuredComm, &rescalingCommGrRootGlobal[0]);
    }

    MPI_Barrier(m_StructuredComm, AT_);
    MPI_Bcast(rescalingCommGrRoot, 1, MPI_INT, rescalingCommGrRanks[0], m_StructuredComm, AT_, "rescalingCommGrRoot");
    MPI_Bcast(rescalingCommGrRootGlobal, 1, MPI_INT, rescalingCommGrRanks[0], m_StructuredComm, AT_,
              "rescalingCommGrRootGlobal");
  }

  MPI_Barrier(m_StructuredComm, AT_);

  //////////////////////////////////////////////////////////////////
  ///////////////// MULTIBLOCK/PERIODIC COMMUNICATION  /////////////
  //////////////////////////////////////////////////////////////////

  vector<unique_ptr<StructuredWindowMap<nDim>>> addComm6000Recv; // store the maps for the communication
  vector<unique_ptr<StructuredWindowMap<nDim>>> addComm6000Snd;
  vector<MInt> adjacentPartitionBC6000;
  // we do know from the addCommBC which side we possess but we need to search for the opposite side
  // and the partition containing the opposite side!!!

  // correct indices only for local multiblock bc comm maps
  for(MInt i = 0; i < (MInt)addCommBC.size(); i++) {
    unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
    MInt dimcount = 0;
    for(MInt dim = 0; dim < nDim; dim++) {
      if(addCommBC[i]->start1[dim] != addCommBC[i]->end1[dim]
         && (addCommBC[i]->BC <= -6000 && addCommBC[i]->BC > -6010)) { // 6000er-change
        if(addCommBC[i]->step2[addCommBC[i]->order[dim]] > 0) {
          // TODO_SS labels:FV the next if-clause and the following ones, are a first attempt to fix the bug, but
          //      consistent treatment of corners etc. to have e.g. solutions which are independent of
          //      particular block topologies is still necessary
          if(addCommBC[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
            addCommBC[i]->start1[dim] -= m_noGhostLayers;
            addCommBC[i]->start2[addCommBC[i]->order[dim]] -= m_noGhostLayers;
          }
          if(addCommBC[i]->end1[dim] == m_myMapWithoutGC->end2[dim]) {
            addCommBC[i]->end1[dim] += m_noGhostLayers;
            addCommBC[i]->end2[addCommBC[i]->order[dim]] += m_noGhostLayers;
          }
        } else {
          if(addCommBC[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
            addCommBC[i]->start1[dim] -= m_noGhostLayers;
            addCommBC[i]->start2[addCommBC[i]->order[dim]] += m_noGhostLayers;
          }
          if(addCommBC[i]->end1[dim] == m_myMapWithoutGC->end2[dim]) {
            addCommBC[i]->end1[dim] += m_noGhostLayers;
            addCommBC[i]->end2[addCommBC[i]->order[dim]] -= m_noGhostLayers;
          }
        }
      }
      if(addCommBC[i]->start1[dim] == addCommBC[i]->end1[dim]
         && (addCommBC[i]->BC <= -6000 && addCommBC[i]->BC > -6010)) // 6000er-change
        dimcount++;
    }
    // TODO labels:FV How is it in 2D?
    if(dimcount == 3) {
      cout << "error!!!!!! 0d surface found!!  check it .... " << endl;
    }

    // now compare the multiblock/periodic bc maps to all partitions
    for(MInt j = 0; j < (MInt)m_partitionMapsWithoutGC.size(); j++) {
      mapCombine21(addCommBC[i], m_partitionMapsWithoutGC[j], temp);
      temp->BC = addCommBC[i]->BC;
      temp->Nstar = addCommBC[i]->Nstar;
      temp->SingularId = addCommBC[i]->SingularId;
      temp->dc1 = addCommBC[i]->dc1;
      temp->dc2 = addCommBC[i]->dc2;
      MBool test = false;

      if(addCommBC[i]->BC <= -6000 && addCommBC[i]->BC > -6010) //(addCommBC[i]->BC==6000) //6000er-change
      {
        if(dimcount == 1)
          test = mapCheck2d(temp);
        else if(dimcount == 2)
          test = mapCheck1d(temp);
        else if(dimcount == 3)
          test = mapCheck0d(temp);
      } else {
        test = mapCheck(temp); // map does exist
      }
      if(test == true) {
        addComm6000Recv.push_back(std::move(temp));
        adjacentPartitionBC6000.push_back(j);
        temp = make_unique<StructuredWindowMap<nDim>>();
      }
    }
  }


  // first correct the SndMapBC6000 only for multiblock by GC layers
  for(MInt i = 0; i < (MInt)SndMapBC6000.size(); i++) {
    MInt recvDom = -1;
    for(MInt j = 0; j < noDomains(); ++j) {
      if(SndMapBC6000[i]->Id2 == j) {
        recvDom = j;
        break;
      }
    }
    if(recvDom == -1) mTerm(1, "");

    MInt dimcount = 0;
    for(MInt dim = 0; dim < nDim; dim++) {
      if(SndMapBC6000[i]->start1[dim] != SndMapBC6000[i]->end1[dim]
         && (SndMapBC6000[i]->BC <= -6000 && SndMapBC6000[i]->BC > -6010)) { // 6000er-change
        if(SndMapBC6000[i]->step2[SndMapBC6000[i]->order[dim]] > 0) {
          // TODO_SS labels:FV the next if-clause and the following ones, are a first attempt to fix the bug, but
          //      consistent treatment of corners etc. to have e.g. solutions which are independent of
          //      particular block topologies is still necessary
          if(SndMapBC6000[i]->start2[SndMapBC6000[i]->order[dim]]
             == m_partitionMapsWithoutGC[recvDom]->start2[SndMapBC6000[i]->order[dim]]) {
            SndMapBC6000[i]->start1[dim] -= m_noGhostLayers;
            SndMapBC6000[i]->start2[SndMapBC6000[i]->order[dim]] -= m_noGhostLayers;
          }
          if(SndMapBC6000[i]->end2[SndMapBC6000[i]->order[dim]]
             == m_partitionMapsWithoutGC[recvDom]->end2[SndMapBC6000[i]->order[dim]]) {
            SndMapBC6000[i]->end1[dim] += m_noGhostLayers;
            SndMapBC6000[i]->end2[SndMapBC6000[i]->order[dim]] += m_noGhostLayers;
          }
        } else {
          if(SndMapBC6000[i]->start2[SndMapBC6000[i]->order[dim]]
             == m_partitionMapsWithoutGC[recvDom]->end2[SndMapBC6000[i]->order[dim]]) {
            SndMapBC6000[i]->start1[dim] -= m_noGhostLayers;
            SndMapBC6000[i]->start2[SndMapBC6000[i]->order[dim]] += m_noGhostLayers;
          }
          if(SndMapBC6000[i]->end2[SndMapBC6000[i]->order[dim]]
             == m_partitionMapsWithoutGC[recvDom]->start2[SndMapBC6000[i]->order[dim]]) {
            SndMapBC6000[i]->end1[dim] += m_noGhostLayers;
            SndMapBC6000[i]->end2[SndMapBC6000[i]->order[dim]] -= m_noGhostLayers;
          }
        }
      }
      if(SndMapBC6000[i]->start1[dim] == SndMapBC6000[i]->end1[dim]
         && (SndMapBC6000[i]->BC <= -6000 && SndMapBC6000[i]->BC > -6010)) // 6000er-change
        dimcount++;
    }
    // TODO labels:FV How is it in 2D?
    if(dimcount == 3) {
      cout << "error!!!!!! 0d surface found!!  check it .... " << endl;
    }

    // now compare multiblock/periodic bc map (with ghost-cells) with own map (w/o ghost-cells)
    // for(MInt i=0; i<(MInt)SndMapBC6000.size(); i++) {
    unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();

    mapCombine11(m_myMapWithoutGC, SndMapBC6000[i], temp);
    temp->BC = SndMapBC6000[i]->BC;
    temp->Nstar = SndMapBC6000[i]->Nstar;
    temp->SingularId = SndMapBC6000[i]->SingularId;
    temp->dc1 = SndMapBC6000[i]->dc1;
    temp->dc2 = SndMapBC6000[i]->dc2;

    MBool test = false;
    if(SndMapBC6000[i]->BC <= -6000 && SndMapBC6000[i]->BC > -6010) //(SndMapBC6000[i]->BC==6000) //6000er-change
    {
      if(dimcount == 1)
        test = mapCheck2d(temp);
      else if(dimcount == 2)
        test = mapCheck1d(temp);
      else if(dimcount == 3)
        test = mapCheck0d(temp);
    } else {
      test = mapCheck(temp); // map does exist
    }
    if(test == true) {
      addComm6000Snd.push_back(std::move(temp));
      temp = make_unique<StructuredWindowMap<nDim>>();
    }
  }

  // special treatment for certain cases where number of SND and RCV maps is unequal
  MInt count1[4] = {0, 0, 0, 0}, count2[4] = {0, 0, 0, 0};
  for(MInt i = 0; i < (MInt)addComm6000Snd.size(); ++i) {
    if((addComm6000Snd[i]->BC <= -6000 && addComm6000Snd[i]->BC > -6010)
       && addComm6000Snd[i]->Nstar == -1) // 6000er-change
      count1[0]++;
    else if((addComm6000Snd[i]->BC <= -6000 && addComm6000Snd[i]->BC > -6010)
            && addComm6000Snd[i]->Nstar != -1) // 6000er-change
      count1[1]++;
    else if(addComm6000Snd[i]->BC >= 4000 && addComm6000Snd[i]->BC < 5000 && addComm6000Snd[i]->Nstar == -1)
      count1[2]++;
    else if(addComm6000Snd[i]->BC >= 4000 && addComm6000Snd[i]->BC < 5000 && addComm6000Snd[i]->Nstar != -1)
      count1[3]++;
  }
  for(MInt i = 0; i < (MInt)addComm6000Recv.size(); ++i) {
    if((addComm6000Recv[i]->BC <= -6000 && addComm6000Recv[i]->BC > -6010)
       && addComm6000Recv[i]->Nstar == -1) // 6000er-change
      count2[0]++;
    else if((addComm6000Recv[i]->BC <= -6000 && addComm6000Recv[i]->BC > -6010)
            && addComm6000Recv[i]->Nstar != -1) // 6000er-change
      count2[1]++;
    else if(addComm6000Recv[i]->BC >= 4000 && addComm6000Recv[i]->BC < 5000 && addComm6000Recv[i]->Nstar == -1)
      count2[2]++;
    else if(addComm6000Recv[i]->BC >= 4000 && addComm6000Recv[i]->BC < 5000 && addComm6000Recv[i]->Nstar != -1)
      count2[3]++;
  }

  MPI_Barrier(m_StructuredComm, AT_);


  //////////////////////////////////////////////////////////////////
  ///////////////// CORRECTIONS FOR BC 6000 (MULTIBLOCK) ///////////
  //////////////////////////////////////////////////////////////////

  // determine faces of receiving maps
  // and make maps three-dimensional
  // RCV maps, no singularity
  for(MInt i = 0; i < (MInt)addComm6000Recv.size(); ++i) {
    if(addComm6000Recv[i]->Nstar == -1
       && (addComm6000Recv[i]->BC <= -6000 && addComm6000Recv[i]->BC > -6010)) { // 6000er-change
      unique_ptr<StructuredWindowMap<nDim>> tempRCV = make_unique<StructuredWindowMap<nDim>>();
      mapCpy(addComm6000Recv[i], tempRCV);

      // IMPORTANT:
      // the faces are important for the unique send and rcv tags
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->step2[dim] < 0) {
          MInt haha = addComm6000Recv[i]->start2[dim];
          addComm6000Recv[i]->start2[dim] = addComm6000Recv[i]->end2[dim];
          addComm6000Recv[i]->end2[dim] = haha;
        }
      }

      // check the side for the receiving parts and save the face
      for(MInt part = 0; part < (MInt)m_partitionMapsWithoutGC.size(); ++part) {
        if(addComm6000Recv[i]->Id2 == m_partitionMapsWithoutGC[part]->Id2) {
          for(MInt dim = 0; dim < nDim; dim++) {
            if(addComm6000Recv[i]->start2[dim] == addComm6000Recv[i]->end2[dim]) {
              switch(dim) {
                case 0: { // face 0 or 1
                          // check if face 0 or 1
                  if(addComm6000Recv[i]->start2[dim] == m_partitionMapsWithoutGC[part]->start2[dim]) {
                    // test of start2 instead of start1.
                    tempRCV->face = 0;
                  } else {
                    tempRCV->face = 1;
                  }
                  break;
                }
                case 1: { // face 2 or 3
                          // check if face 2 or 3
                  if(addComm6000Recv[i]->start2[dim] == m_partitionMapsWithoutGC[part]->start2[dim]) {
                    tempRCV->face = 2;
                  } else {
                    tempRCV->face = 3;
                  }
                  break;
                }
                case 2: { // face 4 or 5
                          // check if face 4 or 5
                  if(addComm6000Recv[i]->start2[dim] == m_partitionMapsWithoutGC[part]->start2[dim]) {
                    tempRCV->face = 4;
                  } else {
                    tempRCV->face = 5;
                  }
                  break;
                }
                default: {
                  cerr << "error no side could be attributed" << endl;
                  exit(1);
                }
              }
            }
          }
        }
      }

      // //check how many dimensions the map has
      // MInt mapDimension = 0;
      // for(MInt dim=0; dim<nDim; dim++) {
      //        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
      //          mapDimension++;
      //        }
      // }
      // //this is a 2d face with one zero dimension
      // //in that case extend in both finite dimensions
      // if(mapDimension==1||mapDimension==2) {
      //        for(MInt dim=0; dim<nDim; dim++) {
      //          if(addComm6000Recv[i]->start1[dim]!=addComm6000Recv[i]->end1[dim]) {
      //              tempRCV->start1[dim]-=m_noGhostLayers;
      //              tempRCV->end1[dim]+=m_noGhostLayers;
      //          }
      //        }
      // }

      // make the RCV maps three-dimensional
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
          // values are the same:
          // change the send and receive values
          if(addComm6000Recv[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
            // we are on the start side of the domain
            tempRCV->start1[dim] -= m_noGhostLayers;
          } else {
            // we are at the end side of the domain
            tempRCV->end1[dim] += m_noGhostLayers;
          }
        }
      }
      tempRCV->BC = addComm6000Recv[i]->BC;
      rcvMap.push_back(std::move(tempRCV));
    }
  }


  // determine faces of receiving maps
  // and make maps three-dimensional
  // RCV maps, with singularity
  for(MInt i = 0; i < (MInt)addComm6000Recv.size(); ++i) {
    if(addComm6000Recv[i]->Nstar != -1
       && (addComm6000Recv[i]->BC <= -6000 && addComm6000Recv[i]->BC > -6010)) { // 6000er-change
      // unique_ptr<StructuredWindowMap<nDim>> tempSND= make_unique<StructuredWindowMap<nDim>>();
      unique_ptr<StructuredWindowMap<nDim>> tempRCV = make_unique<StructuredWindowMap<nDim>>();
      mapCpy(addComm6000Recv[i], tempRCV);

      // IMPORTANT:
      // the faces are important for the unique send and rcv tags
      MInt dimcount = 0;
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
          dimcount++;
        }
      }

      if(dimcount == 2) {
        for(MInt dim = 0; dim < nDim; dim++) {
          if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
            // values are the same:
            // change the send and receive values
            if(addComm6000Recv[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // corner
              // rcv Maps:
              tempRCV->end1[dim] += 1;
            } else {
              // rcv Maps:
              tempRCV->start1[dim] -= 1;
            }
          } else {
            // line
            // rcv Maps:
            // tempRCV->start1[dim]-=m_noGhostLayers;
            // tempRCV->end1[dim]  +=m_noGhostLayers;
          }
        }
      } else if(dimcount == 3) {
        MInt dimN;
        if(addComm6000Recv[i]->BC >= 4400 && addComm6000Recv[i]->BC < 4410) {
          addComm6000Recv[i]->face = addComm6000Recv[i]->BC - 4401;
        }
        if(addComm6000Recv[i]->face != -1) {
          dimN = addComm6000Recv[i]->face / 2;
        } else {
          cout << "ERROR!!! point singular communication cannot decide the direction!! check the singularity!!!"
               << endl;
          exit(0);
        }

        for(MInt dim = 0; dim < nDim; dim++) {
          if(dim != dimN) {
            // values are the same:
            // change the send and receive values
            if(addComm6000Recv[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // corner
              // rcv Maps:
              tempRCV->end1[dim] += 1;

            } else {
              // rcv Maps:
              tempRCV->start1[dim] -= 1;
            }

          } else {
            // rcv Maps:
            if(addComm6000Recv[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // rcv Maps:
              tempRCV->start1[dim] -= 2;

            } else {
              // rcv Maps:
              tempRCV->end1[dim] += 2;
            }
          }
        }
      }

      for(MInt j = 0; j < singularPoint; ++j) {
        unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
        /*      MBool b = false;
                for (MInt n = 0; n < 6; ++n) {
                  if (localSingularMap[j]->BCsingular[n]==addComm6000Recv[i]->BC)
                    b = true;
                }*/
        if(localSingularMap[j]->BC == addComm6000Recv[i]->BC) {
          //        if (b) {
          // TODO_SS labels:FV,toenhance temporary fix. Check if there is a better way
          const MInt temporary = localSingularMap[j]->Id1;
          localSingularMap[j]->Id1 = m_grid->getBlockId(domainId());
          mapCombine11(localSingularMap[j], addComm6000Recv[i], temp);
          localSingularMap[j]->Id1 = temporary;
        }
        MBool test = false;
        test = mapCheck(temp); // map does exist
        if(test == true) {
          MInt Recvblock = m_grid->getBlockId(addComm6000Recv[i]->Id2);
          MInt singnumber = -1;
          for(MInt k = 0; k < addComm6000Recv[i]->Nstar - 3; k++) {
            if(Recvblock == singularity[j].SingularBlockId[k]) {
              singnumber = k;
              break;
            }
          }

          // if (Recvblock== singularity[j].SingularBlockId[0] )
          //   singnumber=0;
          // else  if (Recvblock== singularity[j].SingularBlockId[1] )
          //   singnumber=1;
          if(singnumber == -1) {
            cout << "ERROR!!!!!!!!can not find the correct the displacement!!!!!!" << endl;
            cout << "recvid2:" << addComm6000Recv[i]->Id2 << " recvblock:" << Recvblock
                 << " id:" << singularity[j].SingularBlockId[0] << " " << singularity[j].SingularBlockId[1] << endl;
          }

          // MInt singnumber=addComm6000Recv[i]->Nstar;
          for(MInt dim = 0; dim < nDim; dim++) {
            tempRCV->start1[dim] += singularity[j].displacement[singnumber + 2][dim];
            tempRCV->end1[dim] += singularity[j].displacement[singnumber + 2][dim];
          }
          tempRCV->BCsingular[0] = singularity[j].BCsingular[singnumber + 2];
          //  singularity[j].count++;
          tempRCV->SingularId = j;
          break;
        }
      }

      tempRCV->face = 100;

      if(addComm6000Recv[i]->BC >= 4000 && addComm6000Recv[i]->BC < 5000) {
        tempRCV->BC = addComm6000Recv[i]->BC;
      } else {
        tempRCV->BC = 6333;
      }
      rcvMap.push_back(std::move(tempRCV));
    }
  }

  // determine faces of receiving maps
  // and make maps three-dimensional
  // SND maps, no singularity
  for(MInt i = 0; i < (MInt)addComm6000Snd.size(); i++) {
    if(addComm6000Snd[i]->Nstar == -1
       && (addComm6000Snd[i]->BC <= -6000 && addComm6000Snd[i]->BC > -6010)) { // 6000er-change
      unique_ptr<StructuredWindowMap<nDim>> tempSND = make_unique<StructuredWindowMap<nDim>>();
      mapCpy(addComm6000Snd[i], tempSND);
      // IMPORTANT:
      // the faces are important for the unique send and rcv tags
      // first check the side for the sending parts
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim] == addComm6000Snd[i]->end1[dim]) {
          switch(dim) {
            case 0: { // face 0 or 1
                      // check if face 0 or 1
              if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
                // test of start2 instead of start1.
                tempSND->face = 0;
              } else {
                tempSND->face = 1;
              }
              break;
            }
            case 1: { // face 2 or 3
                      // check if face 2 or 3
              if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
                tempSND->face = 2;
              } else {
                tempSND->face = 3;
              }
              break;
            }
            case 2: { // face 4 or 5
                      // check if face 4 or 5
              if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
                tempSND->face = 4;
              } else {
                tempSND->face = 5;
              }
              break;
            }
            default: {
              cerr << "error no side could be attributed" << endl;
              exit(1);
            }
          }
        }
      }

      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim] == addComm6000Snd[i]->end1[dim]) {
          // values are the same:
          // change the send and receive values
          if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
            // we are on the start side of the domain
            tempSND->end1[dim] += m_noGhostLayers;
          } else {
            // we are at the end side of the domain
            tempSND->start1[dim] -= m_noGhostLayers;
          }
        }
      }

      tempSND->BC = addComm6000Snd[i]->BC;
      sndMap.push_back(std::move(tempSND));
    }
  }

  // determine faces of receiving maps
  // and make maps three-dimensional
  // SND maps, with singularity
  for(MInt i = 0; i < (MInt)addComm6000Snd.size(); i++) {
    if(addComm6000Snd[i]->Nstar != -1
       && (addComm6000Snd[i]->BC <= -6000 && addComm6000Snd[i]->BC > -6010)) { // 6000er-change
      unique_ptr<StructuredWindowMap<nDim>> tempSND = make_unique<StructuredWindowMap<nDim>>();
      mapCpy(addComm6000Snd[i], tempSND);

      MInt dimcount = 0;
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
          dimcount++;
        }
      }

      if(dimcount == 2) {
        for(MInt dim = 0; dim < nDim; dim++) {
          if(addComm6000Snd[i]->start1[dim] == addComm6000Snd[i]->end1[dim]) {
            // values are the same:
            // change the send and receive values
            // corner
            if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // snd Maps:
              tempSND->end1[dim] += 1;
            } else {
              // snd Maps:
              tempSND->start1[dim] -= 1;
            }
          } else {
            // line
            // snd Maps:
            // tempSND->start1[dim]-=m_noGhostLayers;
            // tempSND->end1[dim]+=m_noGhostLayers;
          }
        }
      } else if(dimcount == 3) {
        MInt dimN;
        if(addComm6000Snd[i]->BC >= 4400 && addComm6000Snd[i]->BC < 4410) {
          addComm6000Snd[i]->face = addComm6000Snd[i]->BC - 4401;
        }
        if(addComm6000Snd[i]->face != -1) {
          dimN = addComm6000Snd[i]->face / 2;
        } else {
          cout << "erroor!!! point singular communication cannot decide the direction!! check the singylairy!!!!!!"
               << endl;
          exit(0);
        }

        for(MInt dim = 0; dim < nDim; dim++) {
          if(dim != dimN) {
            // values are the same:
            // change the send and receive values
            // corner
            if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // rcv Maps:
              tempSND->end1[dim] += 1;
            } else {
              // rcv Maps:
              tempSND->start1[dim] -= 1;
            }
          } else {
            // line
            // rcv Maps:
            if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // rcv Maps:
              tempSND->end1[dim] += 2;

            } else {
              // rcv Maps:
              tempSND->start1[dim] -= 2;
            }
          }
        }
      }

      for(MInt j = 0; j < singularPoint; ++j) {
        unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
        if(localSingularMap[j]->BC == addComm6000Snd[i]->BC) {
          // TODO_SS labels:FV,toenhance temporary fix. Check if there is a better way
          const MInt temporary = localSingularMap[j]->Id1;
          localSingularMap[j]->Id1 = m_grid->getBlockId(domainId());
          mapCombine11(localSingularMap[j], addComm6000Snd[i], temp);
          localSingularMap[j]->Id1 = temporary;
        }
        MBool test = false;
        test = mapCheck(temp); // map does exist
        if(test == true) {
          MInt Recvblock = m_grid->getBlockId(addComm6000Snd[i]->Id2);
          MInt singnumber = -1;
          for(MInt k = 0; k < addComm6000Snd[i]->Nstar - 3; k++) {
            if(Recvblock == singularity[j].SingularBlockId[k]) {
              singnumber = k;
              break;
            }
          }

          // if (Recvblock== singularity[j].SingularBlockId[0] )
          //   singnumber=0;
          // else  if (Recvblock== singularity[j].SingularBlockId[1] )
          //   singnumber=1;
          if(singnumber == -1) {
            cout << "ERROR!!!!!!!!can not find the correct the displacement!!!!!!" << endl;
            cout << "recvid2:" << addComm6000Snd[i]->Id2 << " recvblock:" << Recvblock
                 << " id:" << singularity[j].SingularBlockId[0] << " " << singularity[j].SingularBlockId[1] << endl;
          }

          tempSND->BCsingular[0] = singularity[j].BCsingular[singnumber + 2];
          //  singularity[j].count++;
          tempSND->SingularId = j;
          break;
        }
      }

      tempSND->face = 100;

      if(addComm6000Snd[i]->BC >= 4000 && addComm6000Snd[i]->BC < 5000) {
        tempSND->BC = addComm6000Snd[i]->BC;
      } else {
        tempSND->BC = 6333;
      }

      sndMap.push_back(std::move(tempSND));
    }
  }

  //////////////////////////////////////////////////////////////////
  ///////////////// CORRECTIONS FOR BC 4xxx (PERIODIC) /////////////
  //////////////////////////////////////////////////////////////////

  // determine faces of receiving maps
  // and make maps three-dimensional
  // RCV maps, no singularity
  for(MInt i = 0; i < (MInt)addComm6000Recv.size(); ++i) {
    if(addComm6000Recv[i]->Nstar == -1
       && (addComm6000Recv[i]->BC > -6000 || addComm6000Recv[i]->BC <= -6010)) { // 6000er-change
      unique_ptr<StructuredWindowMap<nDim>> tempRCV = make_unique<StructuredWindowMap<nDim>>();
      mapCpy(addComm6000Recv[i], tempRCV);

      // IMPORTANT:
      // the faces are important for the unique send and rcv tags
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->step2[dim] < 0) {
          MInt tmp = addComm6000Recv[i]->start2[dim];
          addComm6000Recv[i]->start2[dim] = addComm6000Recv[i]->end2[dim];
          addComm6000Recv[i]->end2[dim] = tmp;
        }
      }
      // check the side for the receiving parts
      for(MInt part = 0; part < (MInt)m_partitionMapsWithoutGC.size(); ++part) {
        if(addComm6000Recv[i]->Id2 == m_partitionMapsWithoutGC[part]->Id2) {
          for(MInt dim = 0; dim < nDim; dim++) {
            if(addComm6000Recv[i]->start2[dim] == addComm6000Recv[i]->end2[dim]) {
              switch(dim) {
                case 0: { // face 0 or 1
                          // check if face 0 or 1
                  if(addComm6000Recv[i]->start2[dim] == m_partitionMapsWithoutGC[part]->start2[dim]) {
                    // test of start2 instead of start1.
                    tempRCV->face = 0;
                  } else {
                    tempRCV->face = 1;
                  }
                  break;
                }
                case 1: { // face 2 or 3
                          // check if face 2 or 3
                  if(addComm6000Recv[i]->start2[dim] == m_partitionMapsWithoutGC[part]->start2[dim]) {
                    tempRCV->face = 2;
                  } else {
                    tempRCV->face = 3;
                  }
                  break;
                }
                case 2: { // face 4 or 5
                          // check if face 4 or 5
                  if(addComm6000Recv[i]->start2[dim] == m_partitionMapsWithoutGC[part]->start2[dim]) {
                    tempRCV->face = 4;
                  } else {
                    tempRCV->face = 5;
                  }
                  break;
                }
                default: {
                  cerr << "error no side could be attributed" << endl;
                  exit(1);
                }
              }
            }
          }
        }
      }

      // check how many dimensions the map has
      MInt mapDimension = 0;
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
          mapDimension++;
        }
      }

      // this is a 2d face with one zero dimension
      // in that case extend in both finite dimensions
      if(mapDimension == 1) {
        for(MInt dim = 0; dim < nDim; dim++) {
          if(addComm6000Recv[i]->start1[dim] != addComm6000Recv[i]->end1[dim]) {
            tempRCV->start1[dim] -= m_noGhostLayers;
            tempRCV->end1[dim] += m_noGhostLayers;
          }
        }
      }

      // now thicken the map:
      // point: thicken in all three dimensions
      // line: thicken in two zero dimensions
      // face: thicken in the only zero dimension
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
          // values are the same:
          // change the send and receive values
          if(addComm6000Recv[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
            // we are on the start side of the domain
            // rcv Maps:
            tempRCV->start1[dim] -= m_noGhostLayers;
          } else {
            // we are at the end side of the domain
            // rcv Maps:
            tempRCV->end1[dim] += m_noGhostLayers;
          }
        }
      }

      tempRCV->BC = addComm6000Recv[i]->BC;
      rcvMap.push_back(std::move(tempRCV));
    }
  }

  // determine faces of receiving maps
  // and make maps three-dimensional
  // RCV maps, with singularity
  for(MInt i = 0; i < (MInt)addComm6000Recv.size(); ++i) {
    if(addComm6000Recv[i]->Nstar != -1
       && (addComm6000Recv[i]->BC > -6000 || addComm6000Recv[i]->BC <= -6010)) { // 6000er-change
      unique_ptr<StructuredWindowMap<nDim>> tempRCV = make_unique<StructuredWindowMap<nDim>>();
      mapCpy(addComm6000Recv[i], tempRCV);

      // IMPORTANT:
      // the faces are important for the unique send and rcv tags
      MInt dimcount = 0;
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
          dimcount++;
        }
      }

      if(dimcount == 2) {
        for(MInt dim = 0; dim < nDim; dim++) {
          if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
            // values are the same:
            // change the send and receive values
            // corner
            if(addComm6000Recv[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // rcv Maps:
              tempRCV->end1[dim] += 1;
            } else {
              // rcv Maps:
              tempRCV->start1[dim] -= 1;
            }
          } else {
            // line
            // rcv Maps:
            tempRCV->start1[dim] -= m_noGhostLayers;
            tempRCV->end1[dim] += m_noGhostLayers;
          }
        }
      } else if(dimcount == 3) {
        MInt dimN;
        if(addComm6000Recv[i]->BC >= 4400 && addComm6000Recv[i]->BC < 4410) {
          addComm6000Recv[i]->face = addComm6000Recv[i]->BC - 4401;
        }
        if(addComm6000Recv[i]->face != -1) {
          dimN = addComm6000Recv[i]->face / 2;
        } else {
          cout << "erroor!!! point singular communication cannot decide the direction!! check the singylairy!!!!!!"
               << endl;
          exit(0);
        }

        for(MInt dim = 0; dim < nDim; dim++) {
          if(dim != dimN) {
            // values are the same:
            // change the send and receive values
            // corner
            if(addComm6000Recv[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // rcv Maps:
              tempRCV->end1[dim] += 1;
            } else {
              // rcv Maps:
              tempRCV->start1[dim] -= 1;
            }
          } else {
            // line
            // rcv Maps:
            if(addComm6000Recv[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // corner
              // rcv Maps:
              tempRCV->start1[dim] -= 2;
            } else {
              // rcv Maps:
              tempRCV->end1[dim] += 2;
            }
          }
        }
      }

      for(MInt j = 0; j < singularPoint; ++j) {
        unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
        if(localSingularMap[j]->BC == addComm6000Recv[i]->BC)
          mapCombine11(localSingularMap[j], addComm6000Recv[i], temp);
        MBool test = false;
        test = mapCheck(temp); // map does exist
        if(test == true) {
          // communication map change start from 2 to 3
          MInt singnumber = singularity[j].count;
          // MInt singnumber=addComm6000Recv[i]->Nstar;
          for(MInt dim = 0; dim < nDim; dim++) {
            tempRCV->start1[dim] += singularity[j].displacement[singnumber][dim];
            tempRCV->end1[dim] += singularity[j].displacement[singnumber][dim];
          }
          singularity[j].count++;
          tempRCV->SingularId = j;
          break;
        }
      }

      tempRCV->face = 100;

      if(addComm6000Recv[i]->BC >= 4000 && addComm6000Recv[i]->BC < 5000) {
        tempRCV->BC = addComm6000Recv[i]->BC;
      } else {
        tempRCV->BC = 6333;
      }

      rcvMap.push_back(std::move(tempRCV));
    }
  }

  // determine faces of receiving maps
  // and make maps three-dimensional
  // SND maps, no singularity
  for(MInt i = 0; i < (MInt)addComm6000Snd.size(); i++) {
    if(addComm6000Snd[i]->Nstar == -1
       && (addComm6000Snd[i]->BC > -6000 || addComm6000Snd[i]->BC <= -6010)) { // 6000er-change
      unique_ptr<StructuredWindowMap<nDim>> tempSND = make_unique<StructuredWindowMap<nDim>>();
      mapCpy(addComm6000Snd[i], tempSND);
      // IMPORTANT:
      // the faces are important for the unique send and rcv tags

      // first check the side for the sending parts
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim] == addComm6000Snd[i]->end1[dim]) {
          switch(dim) {
            case 0: { // face 0 or 1
                      // check if face 0 or 1
              if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
                // test of start2 instead of start1.
                tempSND->face = 0;
              } else {
                tempSND->face = 1;
              }
              break;
            }
            case 1: { // face 2 or 3
                      // check if face 2 or 3
              if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
                tempSND->face = 2;
              } else {
                tempSND->face = 3;
              }
              break;
            }
            case 2: { // face 4 or 5
                      // check if face 4 or 5
              if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
                tempSND->face = 4;
              } else {
                tempSND->face = 5;
              }
              break;
            }
            default: {
              cerr << "error no side could be attributed" << endl;
              exit(1);
            }
          }
        }
      }

      // check how many dimensions the map has
      MInt mapDimension = 0;
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim] == addComm6000Snd[i]->end1[dim]) {
          mapDimension++;
        }
      }

      // this is a 2d face with one zero dimension
      // in that case extend in both finite dimensions
      if(mapDimension == 1) {
        for(MInt dim = 0; dim < nDim; dim++) {
          if(addComm6000Snd[i]->start1[dim] != addComm6000Snd[i]->end1[dim]) {
            // if(tempSND->end1[dim]-tempSND->start1[dim]!=m_noGhostLayers) {
            tempSND->start1[dim] -= m_noGhostLayers;
            tempSND->end1[dim] += m_noGhostLayers;
            // }
          }
        }
      }

      // now thicken the map:
      // point: thicken in all three dimensions
      // line: thicken in two zero dimensions
      // face: thicken in the only zero dimension
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim] == addComm6000Snd[i]->end1[dim]) {
          // values are the same:
          // change the send and receive values
          if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
            // we are on the start side of the domain
            // snd Maps:
            tempSND->end1[dim] += m_noGhostLayers;
          } else {
            // we are at the end side of the domain
            // snd Maps:
            tempSND->start1[dim] -= m_noGhostLayers;
          }
        }
      }

      tempSND->BC = addComm6000Snd[i]->BC;
      sndMap.push_back(std::move(tempSND));
    }
  }

  // determine faces of receiving maps
  // and make maps three-dimensional
  // SND maps, with singularity
  for(MInt i = 0; i < (MInt)addComm6000Snd.size(); i++) {
    if(addComm6000Snd[i]->Nstar != -1
       && (addComm6000Snd[i]->BC > -6000 || addComm6000Snd[i]->BC <= -6010)) { // 6000er-change
      unique_ptr<StructuredWindowMap<nDim>> tempSND = make_unique<StructuredWindowMap<nDim>>();
      mapCpy(addComm6000Snd[i], tempSND);

      MInt dimcount = 0;
      for(MInt dim = 0; dim < nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim] == addComm6000Recv[i]->end1[dim]) {
          dimcount++;
        }
      }

      if(dimcount == 2) {
        for(MInt dim = 0; dim < nDim; dim++) {
          if(addComm6000Snd[i]->start1[dim] == addComm6000Snd[i]->end1[dim]) {
            // values are the same:
            // change the send and receive values
            // corner
            if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // snd Maps:
              tempSND->end1[dim] += 1;
            } else {
              // snd Maps:
              tempSND->start1[dim] -= 1;
            }
          } else {
            // line
            // snd Maps:
            tempSND->start1[dim] -= m_noGhostLayers;
            tempSND->end1[dim] += m_noGhostLayers;
          }
        }
      } else if(dimcount == 3) {
        MInt dimN;
        if(addComm6000Snd[i]->BC >= 4400 && addComm6000Snd[i]->BC < 4410) {
          addComm6000Snd[i]->face = addComm6000Snd[i]->BC - 4401;
        }
        if(addComm6000Snd[i]->face != -1) {
          dimN = addComm6000Snd[i]->face / 2;
        } else {
          cout << "erroor!!! point singular communication cannot decide the direction!! check the singylairy!!!!!!"
               << endl;
          exit(0);
        }

        for(MInt dim = 0; dim < nDim; dim++) {
          if(dim != dimN) {
            // values are the same:
            // change the send and receive values
            // corner
            if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // rcv Maps:
              tempSND->end1[dim] += 1;
            } else {
              // rcv Maps:
              tempSND->start1[dim] -= 1;
            }
          } else {
            // line
            // rcv Maps:
            if(addComm6000Snd[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // corner
              // rcv Maps:
              tempSND->end1[dim] += 2;
            } else {
              // rcv Maps:
              tempSND->start1[dim] -= 2;
            }
          }
        }
      }

      tempSND->face = 100;

      if(addComm6000Snd[i]->BC >= 4000 && addComm6000Snd[i]->BC < 5000) {
        tempSND->BC = addComm6000Snd[i]->BC;
      } else {
        tempSND->BC = 6333;
      }

      sndMap.push_back(std::move(tempSND));
    }
  }

  ////////////////////////////////////////////////////
  /////////// SINGULARITY CORRECTION  ////////////////
  ////////////////////////////////////////////////////

  for(MInt i = 0; i < singularPoint; ++i) {
    for(MInt dim = 0; dim < nDim; dim++) {
      singularity[i].Viscous[dim] = 0;

      if(singularity[i].start[dim] == singularity[i].end[dim]) {
        // values are the same:
        if(singularity[i].start[dim] == m_myMapWithoutGC->start2[dim]) {
          // corner
          singularity[i].end[dim] += 1;
          singularity[i].Viscous[dim] = -1;
        } else {
          singularity[i].start[dim] -= 1;
        }
      } else {
        singularity[i].start[dim] -= m_noGhostLayers;
        singularity[i].end[dim] += m_noGhostLayers;
      }
    }
  }

  //=============> now the snd and rcv maps contain the following order
  // 1) all the communication because of partitioning
  // 2) communication to other partition/block because of multiblock 6000 bc or periodic boundary condition 4xxx

  ///////////////////////////////////////////////////////////
  //////////////// PHYSICAL BC MAP CORRECTION ///////////////
  ///////////////////////////////////////////////////////////
  // correct all the start2/end2 indices on the physical maps
  // to get correct results from mapCombine
  for(MInt i = 0; i < (MInt)physicalBCMap.size(); i++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      physicalBCMap[i]->start2[dim] += m_noGhostLayers;
      physicalBCMap[i]->end2[dim] += m_noGhostLayers;
    }
  }

  // until now, all physical BCs are still 2D or 1D,
  // now thicken them to 3D and extend them if necessary
  for(MInt i = 0; i < (MInt)physicalBCMap.size(); i++) {
    /*    if (physicalBCMap[i]->Nstar!=-1) {
                unique_ptr<StructuredWindowMap<nDim>> map = physicalBCMap[i];
          cout << "TEST00 BC_" << i << " BC=" << map->BC << "; start1=" << map->start1[0] << "|" << map->start1[1]// <<
       "|" << map->start1[2]
               << "; end1=" << map->end1[0] << "|" << map->end1[1] << endl;// << "|" << map->end1[2] << endl
          physicalBCMap[i]->face = -1;
          physicalBCMap[i]->originShape = 0;
          continue;
        }*/

    MInt addDimensionCount = 0;

    for(MInt dim = 0; dim < nDim; dim++) {
      if(physicalBCMap[i]->start1[dim] == physicalBCMap[i]->end1[dim]) {
        // check which face it is and save this information
        // also thicken in the face normal direction by two ghost-cells
        switch(dim) {
          case 0: { // face 0 or 1
            // check if face 0 or 1
            if(physicalBCMap[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              physicalBCMap[i]->face = 0;
              physicalBCMap[i]->start1[dim] -= m_noGhostLayers;
              physicalBCMap[i]->start2[dim] -= m_noGhostLayers;
            } else {
              physicalBCMap[i]->face = 1;
              physicalBCMap[i]->end1[dim] += m_noGhostLayers;
              physicalBCMap[i]->end2[dim] += m_noGhostLayers;
            }
            break;
          }
          case 1: { // face 2 or 3
            if(physicalBCMap[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // check if face 2 or 3
              physicalBCMap[i]->face = 2;
              physicalBCMap[i]->start1[dim] -= m_noGhostLayers;
              physicalBCMap[i]->start2[dim] -= m_noGhostLayers;
            } else {
              physicalBCMap[i]->face = 3;
              physicalBCMap[i]->end1[dim] += m_noGhostLayers;
              physicalBCMap[i]->end2[dim] += m_noGhostLayers;
            }
            break;
          }
          case 2: { // face 4 or 5
            if(physicalBCMap[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
              // check if face 4 or 5
              physicalBCMap[i]->face = 4;
              physicalBCMap[i]->start1[dim] -= m_noGhostLayers;
              physicalBCMap[i]->start2[dim] -= m_noGhostLayers;
            } else {
              physicalBCMap[i]->face = 5;
              physicalBCMap[i]->end1[dim] += m_noGhostLayers;
              physicalBCMap[i]->end2[dim] += m_noGhostLayers;
            }
            break;
          }
          default: {
            cerr << "error no side could be attributed" << endl;
            exit(1);
          }
        }

        addDimensionCount++;
        continue;
      }

      // if the map starts at the begin of the domain, extend in negative direction
      if(physicalBCMap[i]->start1[dim] == m_myMapWithoutGC->start2[dim]) {
        physicalBCMap[i]->start1[dim] -= m_noGhostLayers;
        physicalBCMap[i]->start2[dim] -= m_noGhostLayers;
      }

      // if the map end at the end of the domain, extend in positive direction
      if(physicalBCMap[i]->end1[dim] == m_myMapWithoutGC->end2[dim]) {
        physicalBCMap[i]->end1[dim] += m_noGhostLayers;
        physicalBCMap[i]->end2[dim] += m_noGhostLayers;
      }
    }

    physicalBCMap[i]->originShape = addDimensionCount;
  }

  // I guess this call would have been also possible before the "PHYSICAL BC MAP CORRECTION"?!
  deleteDuplicateBCMaps(physicalBCMap);

  ///////////////////////////////////////////////////////////
  ///////// EXTENSION CORRECTION ////////////////////////////
  ///////////////////////////////////////////////////////////

  // remove the extension in the cases where a partition boundary
  // coincides with two different BCs
  //(if we don't remove it, there might be two overlapping BCs with
  // uncertainty which one is applied first and second, this
  // can cause oscillations and NANs eventually)

  // go over all physical maps
  for(MInt i = 0; i < (MInt)physicalBCMap.size(); i++) {
    // first loop: only maps that were originally surfaces
    if(physicalBCMap[i]->originShape == 2) {
      continue;
    }
    for(MInt j = 0; j < (MInt)physicalBCMap.size(); j++) {
      //      // I don't know for sure why I need to skip this, but ...
      //      if (physicalBCMap[j]->Nstar!=-1) continue;
      // second loop: only maps that were originally lines
      if(physicalBCMap[j]->originShape != 2) {
        continue;
      }
      // skip if BC number is equal (we can leave the extension then)
      if(physicalBCMap[i]->BC == physicalBCMap[j]->BC) {
        continue;
      }

      mapCombine11(physicalBCMap[i], physicalBCMap[j], localMapDummy);
      MBool is3dMap = mapCheck3d(localMapDummy);

      if(is3dMap) {
        MInt faceDim = physicalBCMap[i]->face / 2;

        for(MInt dim = 0; dim < nDim; dim++) {
          // skip if is the face normal direction
          if(dim == faceDim) {
            continue;
          }

          // start of domain
          if(physicalBCMap[j]->start1[dim] == physicalBCMap[i]->start1[dim]
             && physicalBCMap[j]->end1[dim] == m_myMapWithoutGC->start2[dim]) {
            physicalBCMap[i]->start1[dim] += m_noGhostLayers;
            cout << "DomainID: " << domainId()
                 << " cutting extension at start of partition because of partition/bc conflict with BCs "
                 << physicalBCMap[i]->BC << " and " << physicalBCMap[j]->BC << endl;
          }

          // end of domain
          if(physicalBCMap[j]->start1[dim] == m_myMapWithoutGC->end2[dim]
             && physicalBCMap[j]->end1[dim] == physicalBCMap[i]->end1[dim]) {
            physicalBCMap[i]->end1[dim] -= m_noGhostLayers;
            cout << "DomainID: " << domainId()
                 << " cutting extension at end of partition because of partition/bc conflict with BCs "
                 << physicalBCMap[i]->BC << " and " << physicalBCMap[j]->BC << endl;
          }
        }
      }
    }
  }


  // correct all the start2/end2 indices on the rcv maps
  for(MInt j = 0; j < (MInt)rcvMap.size(); j++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      rcvMap[j]->start2[dim] = rcvMap[j]->start1[dim] + m_grid->getMyOffset(nDim - 1 - dim);
      rcvMap[j]->end2[dim] = rcvMap[j]->end1[dim] + m_grid->getMyOffset(nDim - 1 - dim);
    }
  }

  // Create singularity global BCs from rcv map
  for(MInt j = 0; j < (MInt)rcvMap.size(); j++) {
    if(rcvMap[j]->Nstar > 3 && rcvMap[j]->BCsingular[0] < -6000 && rcvMap[j]->BCsingular[0] > -6010) {
      unique_ptr<StructuredWindowMap<nDim>> windowMap = make_unique<StructuredWindowMap<nDim>>();
      mapCreate(rcvMap[j]->Id1, rcvMap[j]->start1, rcvMap[j]->end1, rcvMap[j]->step1, -1, NULL, NULL, NULL,
                rcvMap[j]->order, -rcvMap[j]->BCsingular[0], windowMap);
      windowMap->Nstar = rcvMap[j]->Nstar;
      windowMap->SingularId = rcvMap[j]->SingularId;
      physicalBCMap.push_back(std::move(windowMap));
    }
  }

  // Check that 6000er BCs only cover the range, which is also covered by -6000er communication BCs
  // Therefore check the edges/extensions; This istuation arises at singularities; we don't want to have
  // the extension in case of Nstar==5 or Nstar==6
  for(MInt i = 0; i < (MInt)physicalBCMap.size(); i++) {
    if(physicalBCMap[i]->BC >= 6000 && physicalBCMap[i]->BC < 6010 && physicalBCMap[i]->Nstar == -1) {
      for(MInt d = 0; d < nDim; ++d) {
        if(physicalBCMap[i]->step1[d] < 0) mTerm(1, "");
      }
      const MInt faceDim = physicalBCMap[i]->face / 2;

      MInt start_begin[nDim], end_begin[nDim];
      MInt start_end[nDim], end_end[nDim];
      std::copy(&physicalBCMap[i]->start1[0], &physicalBCMap[i]->start1[0] + nDim, &start_begin[0]);
      std::copy(&physicalBCMap[i]->start1[0], &physicalBCMap[i]->start1[0] + nDim, &end_begin[0]);
      std::copy(&physicalBCMap[i]->end1[0], &physicalBCMap[i]->end1[0] + nDim, &end_end[0]);
      std::copy(&physicalBCMap[i]->end1[0], &physicalBCMap[i]->end1[0] + nDim, &start_end[0]);
      for(MInt d = 0; d < nDim; ++d) {
        if(physicalBCMap[i]->start1[d] + m_noGhostLayers == m_myMapWithoutGC->start2[d] || d == faceDim)
          end_begin[d] += m_noGhostLayers;
        if(physicalBCMap[i]->end1[d] - m_noGhostLayers == m_myMapWithoutGC->end2[d] || d == faceDim)
          start_end[d] -= m_noGhostLayers;
        if(faceDim != -1 && d != faceDim && end_begin[d] > start_end[d]) {
          //          std::this_thread::sleep_for(std::chrono::milliseconds(2000)*domainId());
          cout << "STOP " << d << "|" << faceDim << " " << start_end[d] << "|" << end_begin[d] << endl;
          mapPrint(physicalBCMap[i]);
          mTerm(1, "stop");
        }
      }

      MBool overlap1 = false;
      MBool overlap2 = false;
      for(MInt j = 0; j < (MInt)rcvMap.size(); j++) {
        //        if (rcvMap[j]->BC>-6000 || rcvMap[j]->BC<=-6010) continue;
        MBool overlap1_temp = true;
        MBool overlap2_temp = true;
        for(MInt d = 0; d < nDim; ++d) {
          if(rcvMap[j]->step1[d] < 0) mTerm(1, "");
        }

        //
        for(MInt d = 0; d < nDim; ++d) {
          if(rcvMap[j]->start1[d] > start_begin[d] || rcvMap[j]->end1[d] < end_begin[d]) overlap1_temp = false;
          if(rcvMap[j]->start1[d] > start_end[d] || rcvMap[j]->end1[d] < end_end[d]) overlap2_temp = false;
        }

        if(overlap1_temp) overlap1 = overlap1_temp;
        if(overlap2_temp) overlap2 = overlap2_temp;
      }

      if(!overlap1) {
        cout << "DomainID: " << domainId() << " cutting extension at start of partition! " << endl;
        for(MInt d = 0; d < nDim; ++d) {
          if(d == faceDim) continue;

          physicalBCMap[i]->start1[d] += m_noGhostLayers;
        }
      }
      if(!overlap2) {
        cout << "DomainID: " << domainId() << " cutting extension at end of partition! " << endl;
        for(MInt d = 0; d < nDim; ++d) {
          if(d == faceDim) continue;

          physicalBCMap[i]->end1[d] -= m_noGhostLayers;
        }
      }
    }
  }


  for(MInt i = 0; i < (MInt)physicalBCMap.size(); i++) {
    if(physicalBCMap[i]->BC == 2501) {
      continue;
    }
    for(MInt j = 0; j < (MInt)rcvMap.size(); j++) {
      mapCombine11(rcvMap[j], physicalBCMap[i], localMapDummy);
      MBool is3dMap = mapCheck3d(localMapDummy);
      if(is3dMap && ((rcvMap[j]->BC <= -6000 && rcvMap[j]->BC > -6010) || rcvMap[j]->BC == 6333)) { // 6000er-change
        cout << "############ ATTENTION: BC MAP IS OVERLAPPING WITH 6000 RCV MAP!!! #############" << endl;
        cout << "receiveMap: " << endl;
        mapPrint(rcvMap[j]);
        cout << "physicalBCMap: " << endl;
        mapPrint(physicalBCMap[i]);
        cout << "combined21: " << endl;
        mapPrint(localMapDummy);
      }
    }
  }


  /*  MPI_Barrier(m_StructuredComm, AT_);
    std::this_thread::sleep_for(std::chrono::milliseconds(2000*domainId()));
    cout << "MY DOMAINID=" << domainId() << "  blockId="
    << m_grid->getBlockId(domainId()) << endl; for (MInt i = 0; i <
    (signed)physicalBCMap.size(); ++i) { cout << "GLOBAL BCs" << endl; mapPrint(physicalBCMap[i]);
    }

    for (MInt i = 0; i < (signed)rcvMap.size(); ++i) {
      cout << "RCV MAP " << endl;
      mapPrint(rcvMap[i]);
    }

    for (MInt i = 0; i < (signed)sndMap.size(); ++i) {
      cout << "SND MAP " << endl;
      mapPrint(sndMap[i]);
    }

    for (MInt i = 0; i < singularPoint; ++i) {
      cout << "SINGULARITY: d=" << domainId() << " Nstar=" << singularity[i].Nstar
        << " singularBlockId=" << singularity[i].SingularBlockId[0] << "|" <<
    singularity[i].SingularBlockId[1]
        << "|" << singularity[i].SingularBlockId[2] << "|" << singularity[i].SingularBlockId[3]
        << " BC=" << singularity[i].BC << " BC_singularity=" << singularity[i].BCsingular[0]
        << "|" << singularity[i].BCsingular[1] << "|" << singularity[i].BCsingular[2] << "|" <<
    singularity[i].BCsingular[3]
        << " start=" << singularity[i].start[0] << "|" << singularity[i].start[1]
        << " end=" << singularity[i].end[0] << "|" << singularity[i].end[1] << endl;
    }
    MPI_Barrier(m_StructuredComm, AT_);*/
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::createWaveWindowMapping(MInt waveZeroPos) {
  if(nDim != 3) return;
  waveRcvMap.clear();
  waveSndMap.clear();

  MInt start1[3] = {0, 0, 0};
  MInt end1[3] = {0, 0, 0};
  MInt step1[3] = {1, 1, 1};
  MInt order[3] = {0, 1, 2};
  MInt start2[3] = {0, 0, 0};
  MInt end2[3] = {0, 0, 0};
  MInt offsetCells[3] = {0, 0, 0};
  MInt activeCells[3] = {0, 0, 0};

  MInt blockId = m_grid->getBlockId(domainId());

  /////////////////////////////////////////////////////////////
  //////////////// CREATE OWN/PARTITION MAPS //////////////////
  /////////////////////////////////////////////////////////////
  // this is the map of the active cells (no ghostcells) in my own partition
  // shifted by the no of ghost-cells
  for(MInt i = 0; i < nDim; i++) {
    start1[i] = m_grid->getMyOffset(nDim - 1 - i);                          // shifted by the number of ghost layers
    start2[i] = 0;                                                          // shifted by the number of ghost layers
    end1[i] = start1[i] + m_grid->getMyActivePoints(nDim - 1 - i) - 1;
    end2[i] = start2[i] + m_grid->getMyActivePoints(nDim - 1 - i) - 1;
  }

  unique_ptr<StructuredWindowMap<nDim>> waveMyMap = make_unique<StructuredWindowMap<nDim>>();
  mapCreate(blockId, start1, end1, step1, domainId(), start2, end2, step1, order, -1, waveMyMap);

  MInt allCells[3] = {0, 0, 0};
  for(MInt dim = 0; dim < nDim; dim++) {
    allCells[dim] = m_grid->getBlockNoCells(0, nDim - 1 - dim);
  }

  vector<unique_ptr<StructuredWindowMap<nDim>>> waveSndPartitionMaps;
  unique_ptr<StructuredWindowMap<nDim>> localMapDummy = nullptr;
  for(MInt j = 0; j < noDomains(); j++) {
    localMapDummy = make_unique<StructuredWindowMap<nDim>>();
    blockId = m_grid->getBlockId(j);

    for(MInt dim = 0; dim < nDim; dim++) {
      offsetCells[dim] = m_grid->getOffset(j, nDim - 1 - dim);
      activeCells[dim] = m_grid->getActivePoints(j, nDim - 1 - dim) - 1;
    }

    for(MInt dim = 0; dim < nDim; dim++) {
      start1[dim] = offsetCells[dim];
      end1[dim] = start1[dim] + activeCells[dim];
    }

    mapCreate(blockId, start1, end1, step1, j, start2, end2, step1, order, 0, localMapDummy);
    waveSndPartitionMaps.push_back(std::move(localMapDummy));
  }

  vector<unique_ptr<StructuredWindowMap<nDim>>> wavePartitionMaps;
  // maps of all the partitions moved to the relative system
  for(MInt j = 0; j < noDomains(); j++) {
    localMapDummy = make_unique<StructuredWindowMap<nDim>>();

    // preparation
    MInt shiftedOffset[3] = {0, 0, 0};
    MInt croppedActiveCells[3] = {0, 0, 0};
    MInt overhangCells[3] = {0, 0, 0};
    MInt overhangOffset[3] = {0, 0, 0};
    blockId = m_grid->getBlockId(j);
    for(MInt dim = 0; dim < nDim; dim++) {
      offsetCells[dim] = m_grid->getOffset(j, nDim - 1 - dim);
      activeCells[dim] = m_grid->getActivePoints(j, nDim - 1 - dim) - 1;
    }

    // first fill with the unmoved values
    // important for 0- and 1-direction
    for(MInt dim = 0; dim < nDim; dim++) {
      shiftedOffset[dim] = offsetCells[dim];
      croppedActiveCells[dim] = activeCells[dim];
      overhangCells[dim] = activeCells[dim];
      overhangOffset[dim] = offsetCells[dim];
    }

    if(offsetCells[2] + activeCells[2] <= waveZeroPos) {
      // domain is before wave zero
      // only shift domain
      shiftedOffset[2] = allCells[2] - waveZeroPos + offsetCells[2];
      for(MInt dim = 0; dim < nDim; dim++) {
        start1[dim] = shiftedOffset[dim];
        end1[dim] = start1[dim] + activeCells[dim];
      }
      mapCreate(blockId, start1, end1, step1, j, start2, end2, step1, order, 0, localMapDummy);
      wavePartitionMaps.push_back(std::move(localMapDummy));
    } else if(offsetCells[2] < waveZeroPos && offsetCells[2] + activeCells[2] > waveZeroPos) {
      // create first map
      croppedActiveCells[2] = waveZeroPos - offsetCells[2];
      shiftedOffset[2] = allCells[2] - croppedActiveCells[2];
      for(MInt dim = 0; dim < nDim; dim++) {
        start1[dim] = shiftedOffset[dim];
        end1[dim] = start1[dim] + croppedActiveCells[dim];
      }
      mapCreate(blockId, start1, end1, step1, j, start2, end2, step1, order, 0, localMapDummy);
      wavePartitionMaps.push_back(std::move(localMapDummy));
      localMapDummy = make_unique<StructuredWindowMap<nDim>>();
      // create second map
      overhangCells[2] = activeCells[2] - croppedActiveCells[2];
      overhangOffset[2] = 0;
      for(MInt dim = 0; dim < nDim; dim++) {
        start1[dim] = overhangOffset[dim];
        end1[dim] = start1[dim] + overhangCells[dim];
      }
      mapCreate(blockId, start1, end1, step1, j, start2, end2, step1, order, 1, localMapDummy);
      wavePartitionMaps.push_back(std::move(localMapDummy));
    } else {
      shiftedOffset[2] = offsetCells[2] - waveZeroPos;
      for(MInt dim = 0; dim < nDim; dim++) {
        start1[dim] = shiftedOffset[dim];
        end1[dim] = start1[dim] + activeCells[dim];
      }
      mapCreate(blockId, start1, end1, step1, j, start2, end2, step1, order, 1, localMapDummy);
      wavePartitionMaps.push_back(std::move(localMapDummy));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////
  ////////////////// USE MAPS TO CHECK FOR OVERLAPPING PARTS ///////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  // check overlapping of own non-gc partition with other gc partitions and put into SND maps
  localMapDummy = make_unique<StructuredWindowMap<nDim>>();
  for(MInt i = 0; i < (MInt)wavePartitionMaps.size(); i++) {
    // only do this for own maps
    if(wavePartitionMaps[i]->Id2 != domainId()) {
      continue;
    }

    for(MInt j = 0; j < (MInt)waveSndPartitionMaps.size(); j++) {
      mapCombineWave(wavePartitionMaps[i], waveSndPartitionMaps[j], localMapDummy);
      MBool test = false;
      test = mapCheckWave(localMapDummy);
      if(test == true) {
        localMapDummy->BC = wavePartitionMaps[i]->BC;
        if(localMapDummy->start1[2] < allCells[2] - waveZeroPos) {
          const MInt translationK = allCells[2] - waveZeroPos - localMapDummy->start1[2];
          const MInt sizeK = localMapDummy->end1[2] - localMapDummy->start1[2];
          localMapDummy->start1[2] = allCells[2] - translationK;
          localMapDummy->end1[2] = localMapDummy->start1[2] + sizeK;
        } else {
          const MInt translationK = allCells[2] - waveZeroPos;
          localMapDummy->start1[2] -= translationK;
          localMapDummy->end1[2] -= translationK;
        }
        waveSndMap.push_back(std::move(localMapDummy));
        localMapDummy = make_unique<StructuredWindowMap<nDim>>();
      }
    }
  }

  // check overlapping of own gc partition with other non-gc partitions and put into RCV maps
  for(MInt i = 0; i < (MInt)wavePartitionMaps.size(); i++) {
    mapCombineWave(waveMyMap, wavePartitionMaps[i], localMapDummy);
    MBool test = false;
    test = mapCheckWave(localMapDummy);
    if(test == true) {
      waveRcvMap.push_back(std::move(localMapDummy));
      localMapDummy = make_unique<StructuredWindowMap<nDim>>();
    }
  }

  for(MInt i = 0; i < (MInt)waveRcvMap.size(); i++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      const MInt offset = m_grid->getMyOffset(nDim - 1 - dim);
      waveRcvMap[i]->start1[dim] = waveRcvMap[i]->start1[dim] - offset + m_noGhostLayers;
      waveRcvMap[i]->end1[dim] = waveRcvMap[i]->end1[dim] - offset + m_noGhostLayers;
    }
  }

  for(MInt i = 0; i < (MInt)waveSndMap.size(); i++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      const MInt offset = m_grid->getMyOffset(nDim - 1 - dim);
      waveSndMap[i]->start1[dim] = waveSndMap[i]->start1[dim] - offset + m_noGhostLayers;
      waveSndMap[i]->end1[dim] = waveSndMap[i]->end1[dim] - offset + m_noGhostLayers;
    }
  }

  waveSndPartitionMaps.clear();
  wavePartitionMaps.clear();
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::createWaveCommunicationExchangeFlags(
    vector<unique_ptr<StructuredComm<nDim>>>& sndComm,
    vector<unique_ptr<StructuredComm<nDim>>>& rcvComm,
    const MInt noVariables) {
  sndComm.clear();
  rcvComm.clear();
  // SND Maps
  for(MInt i = 0; i < (MInt)waveSndMap.size(); i++) {
    // compute the buffersizes for sending and receiving
    MInt noCellsSnd = 1;
    for(MInt j = 0; j < nDim; j++) {
      MInt cellSizes = (waveSndMap[i]->end1[j] - waveSndMap[i]->start1[j]);
      if(cellSizes != 0) {
        noCellsSnd *= (cellSizes);
      }
    }

    StructuredCommType commType = PARTITION_NORMAL;

    unique_ptr<StructuredComm<nDim>> snd =
        make_unique<StructuredComm<nDim>>(noVariables, nullptr, noCellsSnd, 0, commType);

    snd->nghbrId = waveSndMap[i]->Id2;
    for(MInt dim = 0; dim < nDim; ++dim) {
      snd->startInfoCells[dim] = waveSndMap[i]->start1[dim];
      snd->endInfoCells[dim] = waveSndMap[i]->end1[dim];
    }
    snd->tagHelper = waveSndMap[i]->BC;

    sndComm.push_back(std::move(snd));
  }

  // RCV Maps
  for(MInt i = 0; i < (MInt)waveRcvMap.size(); i++) {
    // compute the buffersizes for sending and receiving
    MInt noCellsRcv = 1;
    for(MInt j = 0; j < nDim; j++) {
      MInt cellSizes = (waveRcvMap[i]->end1[j] - waveRcvMap[i]->start1[j]);
      if(cellSizes != 0) {
        noCellsRcv *= (cellSizes);
      }
    }

    StructuredCommType commType = PARTITION_NORMAL;

    unique_ptr<StructuredComm<nDim>> rcv =
        make_unique<StructuredComm<nDim>>(noVariables, nullptr, noCellsRcv, 0, commType);

    rcv->nghbrId = waveRcvMap[i]->Id2;
    for(MInt dim = 0; dim < nDim; ++dim) {
      rcv->startInfoCells[dim] = waveRcvMap[i]->start1[dim];
      rcv->endInfoCells[dim] = waveRcvMap[i]->end1[dim];
    }
    rcv->tagHelper = waveRcvMap[i]->BC;

    rcvComm.push_back(std::move(rcv));
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::createCommunicationExchangeFlags(
    vector<unique_ptr<StructuredComm<nDim>>>& sndComm,
    vector<unique_ptr<StructuredComm<nDim>>>& rcvComm,
    const MInt noVariables,
    MFloat* const* const variables) {
  if(rcvMap.size() != sndMap.size()) {
    MInt count1[4] = {0, 0, 0, 0}, count2[4] = {0, 0, 0, 0};
    for(MInt i = 0; i < (MInt)sndMap.size(); ++i) {
      if((sndMap[i]->BC <= -6000 && sndMap[i]->BC > -6010) && sndMap[i]->Nstar == -1) // 6000er-change
        count1[0]++;
      else if((sndMap[i]->BC <= -6000 && sndMap[i]->BC > -6010) && sndMap[i]->Nstar != -1) // 6000er-change
        count1[1]++;
      else if(sndMap[i]->BC >= 4000 && sndMap[i]->BC < 5000 && sndMap[i]->Nstar == -1)
        count1[2]++;
      else if(sndMap[i]->BC >= 4000 && sndMap[i]->BC < 5000 && sndMap[i]->Nstar != -1)
        count1[3]++;
    }
    for(MInt i = 0; i < (MInt)rcvMap.size(); ++i) {
      if((rcvMap[i]->BC <= -6000 && rcvMap[i]->BC > -6010) && rcvMap[i]->Nstar == -1) // 6000er-change
        count2[0]++;
      else if((rcvMap[i]->BC <= -6000 && rcvMap[i]->BC > -6010) && rcvMap[i]->Nstar != -1) // 6000er-change
        count2[1]++;
      else if(rcvMap[i]->BC >= 4000 && rcvMap[i]->BC < 5000 && rcvMap[i]->Nstar == -1)
        count2[2]++;
      else if(rcvMap[i]->BC >= 4000 && rcvMap[i]->BC < 5000 && rcvMap[i]->Nstar != -1)
        count2[3]++;
    }

    // sending and receiving does not correspond ==> error
    cout << "***********************************WARNING************************************** " << endl;
    cout << "this error may be caused by the partition when you have a Step in the grid, try to change number of "
            "domains "
            "and run again!"
         << endl;
    cout << "******************************************************************************** " << endl;
    mTerm(1, AT_, "number of sending and receiving processes does not match");
  }

  for(MInt i = 0; i < (MInt)rcvMap.size(); i++) {
    StructuredCommType commType = PARTITION_NORMAL;
    const MInt BC = (rcvMap[i]->BC <= -6000 && rcvMap[i]->BC > -6010) ? 6000 : rcvMap[i]->BC; // 6000er-change

    MInt bcIdSnd = sndMap[i]->BC;
    MInt bcIdRcv = rcvMap[i]->BC;

    switch(BC) {
      case 6333: // singular communications
      {
        commType = SINGULAR;
        bcIdSnd = sndMap[i]->BCsingular[0];
        bcIdRcv = rcvMap[i]->BCsingular[0];
        break;
      }
      case 4011:
      case 4012:
      case 4401:
      case 4402:
      case 4403:
      case 4404:
      case 4405:
      case 4406: {
        if(rcvMap[i]->Nstar != -1) {
          commType = PERIODIC_BC_SINGULAR;
        } else {
          commType = PERIODIC_BC;
        }

        break;
      }
      // case 6011: {
      //   commType = CHANNEL;
      //   break;
      // }
      default: {
        commType = PARTITION_NORMAL;
      }
    }

    ////////////////////////////////////////////////////
    ///////////////////// SND COMM /////////////////////
    ////////////////////////////////////////////////////
    // compute the buffersizes for sending and receiving
    MInt noCellsSnd = 1, noPointsSnd = 1;
    for(MInt j = 0; j < nDim; j++) {
      MInt cellSizes = (sndMap[i]->end1[j] - sndMap[i]->start1[j]);
      if(cellSizes != 0) {
        noCellsSnd *= (cellSizes);
        noPointsSnd *= (cellSizes + 1);
      }
    }

    unique_ptr<StructuredComm<nDim>> snd =
        make_unique<StructuredComm<nDim>>(noVariables, variables, noCellsSnd, noPointsSnd, commType);

    snd->nghbrId = sndMap[i]->Id2;

    for(MInt dim = 0; dim < nDim; ++dim) {
      snd->startInfoCells[dim] = sndMap[i]->start1[dim];
      snd->endInfoCells[dim] = sndMap[i]->end1[dim];
      snd->startInfoPoints[dim] = sndMap[i]->start1[dim];
      snd->endInfoPoints[dim] = sndMap[i]->end1[dim];
    }

    snd->tagHelper = sndMap[i]->face;
    if(snd->tagHelper < 0) {
      snd->tagHelper = 0;
    }

    snd->bcId = bcIdSnd;

    sndComm.push_back(std::move(snd));


    ////////////////////////////////////////////////////
    ///////////////////// RCV COMM /////////////////////
    ////////////////////////////////////////////////////
    // compute the buffersizes for sending and receiving
    MInt noCellsRcv = 1, noPointsRcv = 1;
    for(MInt j = 0; j < nDim; j++) {
      MInt cellSizes = (rcvMap[i]->end1[j] - rcvMap[i]->start1[j]);
      if(cellSizes != 0) {
        noCellsRcv *= (cellSizes);
        noPointsRcv *= (cellSizes + 1);
      }
    }

    unique_ptr<StructuredComm<nDim>> rcv =
        make_unique<StructuredComm<nDim>>(noVariables, variables, noCellsRcv, noPointsRcv, commType);

    rcv->nghbrId = rcvMap[i]->Id2;
    for(MInt dim = 0; dim < nDim; ++dim) {
      rcv->startInfoCells[dim] = rcvMap[i]->start1[dim];
      rcv->endInfoCells[dim] = rcvMap[i]->end1[dim];
      rcv->startInfoPoints[dim] = rcvMap[i]->start1[dim];
      rcv->endInfoPoints[dim] = rcvMap[i]->end1[dim];
      rcv->orderInfo[dim] = rcvMap[i]->order[dim];
      rcv->stepInfo[dim] = rcvMap[i]->step2[dim];
    }

    rcv->tagHelper = rcvMap[i]->face;
    if(rcv->tagHelper < 0) {
      rcv->tagHelper = 0;
    }

    rcv->bcId = bcIdRcv;

    rcvComm.push_back(std::move(rcv));
  }
}

template <MInt nDim>
MInt FvStructuredSolverWindowInfo<nDim>::mapCompare(unique_ptr<StructuredWindowMap<nDim>>& map1,
                                                    unique_ptr<StructuredWindowMap<nDim>>& map2) {
  if(map1->Id1 != map2->Id1) return 0;
  if(map1->Id2 != map2->Id2) return 0;
  for(MInt i = 0; i < nDim; i++) {
    if(map1->start1[i] != map2->start1[i]) return 0;
    if(map1->end1[i] != map2->end1[i]) return 0;
    // if(map1->start2[i]!=map2->start2[i]) return 0;
    // if(map1->end2[i]!=map2->end2[i]) return 0;
    if(map1->step1[i] != map2->step1[i]) return 0;
    // if(map1->step2[i]!=map2->step2[i]) return 0;
    if(map1->order[i] != map2->order[i]) return 0;
  }
  return 1;
}

template <MInt nDim>
MInt FvStructuredSolverWindowInfo<nDim>::mapCompare11(const unique_ptr<StructuredWindowMap<nDim>>& map1,
                                                      const unique_ptr<StructuredWindowMap<nDim>>& map2) {
  // only compare the Id1, start1 and end1
  if(map1->Id1 != map2->Id1) return 0;
  for(MInt i = 0; i < nDim; i++) {
    if(map1->start1[i] != map2->start1[i]) return 0;
    if(map1->end1[i] != map2->end1[i]) return 0;
  }
  return 1;
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCreate(MInt Id1, MInt* start1, MInt* end1, MInt* step1, MInt Id2,
                                                   MInt* start2, MInt* end2, MInt* step2, MInt* order, MInt BC,
                                                   unique_ptr<StructuredWindowMap<nDim>>& output) {
  output->Id1 = Id1;
  output->BC = BC;
  output->Nstar = -1;
  output->dc1 = 888;
  output->dc2 = 888;
  output->SingularId = -1;

  if(Id2 != -1) {
    output->Id2 = Id2;
  } else {
    output->Id2 = Id1;
  }

  memcpy(output->start1, start1, nDim * sizeof(MInt));
  if(start2 != nullptr) {
    memcpy(output->start2, start2, nDim * sizeof(MInt));
  } else {
    memcpy(output->start2, start1, nDim * sizeof(MInt));
  }

  memcpy(output->end1, end1, nDim * sizeof(MInt));
  if(end2 != nullptr) {
    memcpy(output->end2, end2, nDim * sizeof(MInt));
  } else {
    memcpy(output->end2, end1, nDim * sizeof(MInt));
  }

  if(step1 != nullptr) {
    memcpy(output->step1, step1, nDim * sizeof(MInt));
  } else {
    for(MInt i = 0; i < nDim; i++)
      output->step1[i] = 1;
  }

  if(step2 != nullptr) {
    memcpy(output->step2, step2, nDim * sizeof(MInt));
  } else {
    for(MInt i = 0; i < nDim; i++)
      output->step2[i] = 1;
  }

  if(order != nullptr) {
    memcpy(output->order, order, nDim * sizeof(MInt));
  } else {
    for(MInt i = 1; i < nDim; i++)
      output->order[i] = i;
  }

  output->SingularId = -1;
  output->Nstar = -1;
  output->dc1 = -999;
  output->dc2 = -999;
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::mapCheck(unique_ptr<StructuredWindowMap<nDim>>& input) {
  MBool test;
  for(MInt i = 0; i < nDim; i++) {
    test = false;
    if(input->step1[i] == 0) break;
    if(input->step2[i] == 0) break;
    if(input->step1[i] > 0) {
      if(input->end1[i] < input->start1[i]) break;
    } else {
      if(input->start1[i] < input->end1[i]) break;
    }
    if(input->step2[i] > 0) {
      if(input->end2[i] < input->start2[i]) break;
    } else {
      if(input->start2[i] < input->end2[i]) break;
    }

    if(((input->end1[i] - input->start1[i]) % (input->step1[i])) != 0) break;
    if(((input->end2[i] - input->start2[i]) % (input->step2[i])) != 0) break;

    if(input->order[i] < 0 || input->order[i] > (nDim - 1)) break;

    if(((input->end1[i] - input->start1[i]) / input->step1[i])
       != ((input->end2[input->order[i]] - input->start2[input->order[i]]) / input->step2[input->order[i]]))
      break;

    test = true;
  }
  // check on Boundary Condition not possible
  return test;
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::mapCheck0d(unique_ptr<StructuredWindowMap<nDim>>& input) {
  MInt dummy = 0;
  MBool test;
  for(MInt i = 0; i < nDim; i++) {
    test = false;
    if(input->step1[i] == 0) break;
    if(input->step2[i] == 0) break;
    if(input->step1[i] > 0) {
      if(input->end1[i] < input->start1[i]) break;
    } else {
      if(input->start1[i] < input->end1[i]) break;
    }
    if(input->step2[i] > 0) {
      if(input->end2[i] < input->start2[i]) break;
    } else {
      if(input->start2[i] < input->end2[i]) break;
    }

    if(((input->end1[i] - input->start1[i]) % (input->step1[i])) != 0) break;
    if(((input->end2[i] - input->start2[i]) % (input->step2[i])) != 0) break;

    if(input->order[i] < 0 || input->order[i] > (nDim - 1)) break;

    if(((input->end1[i] - input->start1[i]) / input->step1[i])
       != ((input->end2[input->order[i]] - input->start2[input->order[i]]) / input->step2[input->order[i]]))
      break;

    // check if the map is a line (3D) or a point (2D-problem)
    // if it is a line then two of the three indicies are the same else for a point
    // only one is
    if(input->start1[i] == input->end1[i]) dummy++;


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // line or point is acceptable!!!!

    test = true;
  }
  if(dummy != 3) test = false;
  // check on Boundary Condition not possible
  return test;
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::mapCheck1d(unique_ptr<StructuredWindowMap<nDim>>& input) {
  // MInt pointorline = 2;
  MInt dummy = 0;
  MBool test;
  for(MInt i = 0; i < nDim; i++) {
    test = false;
    if(input->step1[i] == 0) break;
    if(input->step2[i] == 0) break;
    if(input->step1[i] > 0) {
      if(input->end1[i] < input->start1[i]) break;
    } else {
      if(input->start1[i] < input->end1[i]) break;
    }
    if(input->step2[i] > 0) {
      if(input->end2[i] < input->start2[i]) break;
    } else {
      if(input->start2[i] < input->end2[i]) break;
    }

    if(((input->end1[i] - input->start1[i]) % (input->step1[i])) != 0) break;
    if(((input->end2[i] - input->start2[i]) % (input->step2[i])) != 0) break;

    if(input->order[i] < 0 || input->order[i] > (nDim - 1)) break;

    if(((input->end1[i] - input->start1[i]) / input->step1[i])
       != ((input->end2[input->order[i]] - input->start2[input->order[i]]) / input->step2[input->order[i]]))
      break;

    // check if the map is a line (3D) or a point (2D-problem)
    // if it is a line then two of the three indicies are the same else for a point
    // only one is
    if(input->start1[i] == input->end1[i]) dummy++;


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // line or point is acceptable!!!!

    test = true;
  }
  if(dummy != 2) test = false;
  // check on Boundary Condition not possible
  return test;
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::mapCheck2d(unique_ptr<StructuredWindowMap<nDim>>& input) {
  // MInt pointorline = nDim-1;
  MInt dummy = 0;
  MBool result = false;
  for(MInt i = 0; i < nDim; i++) {
    result = false;
    if(input->step1[i] == 0) break;
    if(input->step2[i] == 0) break;
    if(input->step1[i] > 0) {
      if(input->end1[i] < input->start1[i]) break;
    } else {
      if(input->start1[i] < input->end1[i]) break;
    }
    if(input->step2[i] > 0) {
      if(input->end2[i] < input->start2[i]) break;
    } else {
      if(input->start2[i] < input->end2[i]) break;
    }

    if(((input->end1[i] - input->start1[i]) % (input->step1[i])) != 0) break;
    if(((input->end2[i] - input->start2[i]) % (input->step2[i])) != 0) break;

    if(input->order[i] < 0 || input->order[i] > (nDim - 1)) break;

    if(((input->end1[i] - input->start1[i]) / input->step1[i])
       != ((input->end2[input->order[i]] - input->start2[input->order[i]]) / input->step2[input->order[i]]))
      break;


    // check if the map is a line (3D) or a point (2D-problem)
    // if it is a line then two of the three indicies are the same else for a point
    // only one is
    if(input->start1[i] == input->end1[i]) dummy++;

    result = true;
  }

  if(dummy != 1 || result == false) {
    return false;
  } else {
    return true;
  }
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::mapCheck3d(unique_ptr<StructuredWindowMap<nDim>>& input) {
  // MInt pointorline = nDim-1;
  MInt dummy = 0;
  MBool result = false;
  for(MInt i = 0; i < nDim; i++) {
    result = false;
    if(input->step1[i] == 0) break;
    if(input->step2[i] == 0) break;
    if(input->step1[i] > 0) {
      if(input->end1[i] < input->start1[i]) break;
    } else {
      if(input->start1[i] < input->end1[i]) break;
    }
    if(input->step2[i] > 0) {
      if(input->end2[i] < input->start2[i]) break;
    } else {
      if(input->start2[i] < input->end2[i]) break;
    }

    if(((input->end1[i] - input->start1[i]) % (input->step1[i])) != 0) break;
    if(((input->end2[i] - input->start2[i]) % (input->step2[i])) != 0) break;

    if(input->order[i] < 0 || input->order[i] > (nDim - 1)) break;

    if(((input->end1[i] - input->start1[i]) / input->step1[i])
       != ((input->end2[input->order[i]] - input->start2[input->order[i]]) / input->step2[input->order[i]]))
      break;

    // check if the map is a line (3D) or a point (2D-problem)
    // if it is a line then two of the three indicies are the same else for a point
    // only one is
    if(input->start1[i] == input->end1[i]) dummy++;

    result = true;
  }

  if(dummy != 0 || result == false) {
    return false;
  } else {
    return true;
  }
}

template <MInt nDim>
MBool FvStructuredSolverWindowInfo<nDim>::mapCheckWave(unique_ptr<StructuredWindowMap<nDim>>& input) {
  // MInt pointorline = nDim-1;
  MInt dummy = 0;
  MBool result = false;
  for(MInt i = 0; i < nDim; i++) {
    result = false;
    if(input->step1[i] == 0) break;
    if(input->step1[i] > 0) {
      if(input->end1[i] < input->start1[i]) break;
    } else {
      if(input->start1[i] < input->end1[i]) break;
    }

    if(((input->end1[i] - input->start1[i]) % (input->step1[i])) != 0) break;

    // check if the map is a line (3D) or a point (2D-problem)
    // if it is a line then two of the three indicies are the same else for a point
    // only one is
    if(input->start1[i] == input->end1[i]) dummy++;
    result = true;
  }

  if(dummy != 0 || result == false) {
    return false;
  } else {
    return true;
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapPrint(const unique_ptr<StructuredWindowMap<nDim>>& input) {
  cout << "======== MAP INFO =======" << endl;
  cout << "Id1: " << input->Id1 << endl;
  stringstream start1;
  start1 << "start1: ";
  stringstream start2;
  start2 << "start2: ";
  stringstream end1;
  end1 << "end1: ";
  stringstream end2;
  end2 << "end2: ";
  stringstream step1;
  step1 << "step1: ";
  stringstream step2;
  step2 << "step2: ";
  stringstream order;
  order << "order: ";
  for(MInt i = 0; i < nDim; i++) {
    start1 << input->start1[i] << " ";
    start2 << input->start2[i] << " ";
    end1 << input->end1[i] << " ";
    end2 << input->end2[i] << " ";
    step1 << input->step1[i] << " ";
    step2 << input->step2[i] << " ";
    order << input->order[i] << " ";
  }
  cout << start1.str() << endl;
  cout << end1.str() << endl;
  cout << step1.str() << endl;
  cout << order.str() << endl;
  cout << "Id2: " << input->Id2 << endl;
  cout << start2.str() << endl;
  cout << end2.str() << endl;
  cout << step2.str() << endl;
  cout << "BC: " << input->BC << endl;
  if(input->Nstar != -1)
    cout << "BCsingular: " << input->BCsingular[0] << "|" << input->BCsingular[1] << "|" << input->BCsingular[2] << "|"
         << input->BCsingular[3] << "|" << input->BCsingular[4] << "|" << input->BCsingular[5] << endl;
  cout << "Face: " << input->face << endl;

  cout << "spongInfos:" << endl;
  cout << "hasSponge: " << input->hasSponge << endl;
  cout << "spongeThickness: " << input->spongeThickness << endl;
  cout << "beta: " << input->beta << endl;
  cout << "sigma: " << input->sigma << endl;

  cout << "Nstar: " << input->Nstar << endl;
  cout << "SingularId: " << input->SingularId << endl;
  cout << "DC1: " << input->dc1 << "  DC2: " << input->dc2 << endl;
  cout << "======== MAP INFO END=======" << endl;
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapPrintSimple(unique_ptr<StructuredWindowMap<nDim>>& input) {
  stringstream start1;
  start1 << "start1: ";
  stringstream start2;
  start2 << "start2: ";
  stringstream step2;
  step2 << "step2: ";
  stringstream order;
  order << "order: ";

  for(MInt i = 0; i < nDim; i++) {
    start1 << input->start1[i] << "-" << input->end1[i] << " ";
    start2 << input->start2[i] << "-" << input->end2[i] << " ";
    step2 << input->step2[i] << " ";
    order << input->order[i] << " ";
  }
  cout << "Id1: " << input->Id1 << " Id2: " << input->Id2 << " BC: " << input->BC << "  " << start1.str() << "  "
       << start2.str() << step2.str() << order.str() << " BC:" << input->BC << " Nstar:" << input->Nstar << endl;
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapZero(unique_ptr<StructuredWindowMap<nDim>>& output) {
  // output=make_unique<StructuredWindowMap<nDim>>();
  output->Id1 = 0;
  output->Id2 = 0;
  for(MInt i = 0; i < nDim; i++) {
    output->start1[i] = 0;
    output->end1[i] = -1;
    output->step1[i] = 1;
    output->start2[i] = 0;
    output->end2[i] = -1;
    output->step2[i] = 1;
    output->order[i] = 0;
    output->BC = -1;
  }
  output->hasSponge = false;
  output->spongeThickness = F0;
  output->beta = F0;
  output->sigma = F0;
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapInvert(unique_ptr<StructuredWindowMap<nDim>>& input,
                                                   unique_ptr<StructuredWindowMap<nDim>>& output) {
  // output= make_unique<StructuredWindowMap<nDim>>();
  output->Id1 = input->Id2;
  output->Id2 = input->Id1;

  memcpy(output->start2, input->start1, nDim * sizeof(MInt));

  memcpy(output->end2, input->end1, nDim * sizeof(MInt));
  memcpy(output->step2, input->step1, nDim * sizeof(MInt));
  memcpy(output->start1, input->start2, nDim * sizeof(MInt));
  memcpy(output->end1, input->end2, nDim * sizeof(MInt));
  memcpy(output->step1, input->step2, nDim * sizeof(MInt));

  for(MInt i = 0; i < nDim; i++) {
    output->order[input->order[i]] = i;
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapInvert1(unique_ptr<StructuredWindowMap<nDim>>& output) {
  // output= make_unique<StructuredWindowMap<nDim>>();
  unique_ptr<StructuredWindowMap<nDim>> input;
  input = make_unique<StructuredWindowMap<nDim>>();
  mapCpy(output, input);
  output->Id1 = input->Id2;
  output->Id2 = input->Id1;

  output->dc1 = input->dc2;
  output->dc2 = input->dc1;

  memcpy(output->start2, input->start1, nDim * sizeof(MInt));

  memcpy(output->end2, input->end1, nDim * sizeof(MInt));
  memcpy(output->step2, input->step1, nDim * sizeof(MInt));
  memcpy(output->start1, input->start2, nDim * sizeof(MInt));
  memcpy(output->end1, input->end2, nDim * sizeof(MInt));
  memcpy(output->step1, input->step2, nDim * sizeof(MInt));

  for(MInt i = 0; i < nDim; i++) {
    output->order[input->order[i]] = i;
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCpy(unique_ptr<StructuredWindowMap<nDim>>& input,
                                                unique_ptr<StructuredWindowMap<nDim>>& output) {
  output->Id1 = input->Id1;
  output->Id2 = input->Id2;
  memcpy(output->start1, input->start1, nDim * sizeof(MInt));
  memcpy(output->start2, input->start2, nDim * sizeof(MInt));
  memcpy(output->end1, input->end1, nDim * sizeof(MInt));
  memcpy(output->end2, input->end2, nDim * sizeof(MInt));
  memcpy(output->step1, input->step1, nDim * sizeof(MInt));
  memcpy(output->step2, input->step2, nDim * sizeof(MInt));
  memcpy(output->order, input->order, nDim * sizeof(MInt));
  output->BC = input->BC;
  output->Nstar = input->Nstar;
  output->SingularId = input->SingularId;
  output->dc1 = input->dc1;
  output->dc2 = input->dc2;
  memcpy(output->SingularBlockId, input->SingularBlockId, 4 * sizeof(MInt));
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapNormalize3(unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
  mapCpy(output, temp);
  for(MInt i = 0; i < nDim; i++) {
    if(temp->step1[i] < 0) {
      output->step1[i] = -temp->step1[i];
      output->start1[i] = temp->end1[i];
      output->end1[i] = temp->start1[i];
      output->step2[output->order[i]] = -temp->step2[temp->order[i]];
      output->start2[output->order[i]] = temp->end2[temp->order[i]];
      output->end2[output->order[i]] = temp->start2[temp->order[i]];
    }
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapNormalize1(unique_ptr<StructuredWindowMap<nDim>>& input,
                                                       unique_ptr<StructuredWindowMap<nDim>>& output) {
  mapCpy(input, output);
  for(MInt i = 0; i < nDim; i++) {
    if(input->step1[i] < 0) {
      output->step1[i] = -input->step1[i];
      output->start1[i] = input->end1[i];
      output->end1[i] = input->start1[i];
      output->step2[output->order[i]] = -input->step2[input->order[i]];
      output->start2[output->order[i]] = input->end2[input->order[i]];
      output->end2[output->order[i]] = input->start2[input->order[i]];
    }
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapNormalize2(unique_ptr<StructuredWindowMap<nDim>>& input,
                                                       unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> out = make_unique<StructuredWindowMap<nDim>>();
  unique_ptr<StructuredWindowMap<nDim>> out1 = make_unique<StructuredWindowMap<nDim>>();
  mapInvert(input, out);
  mapNormalize1(out, out1);
  mapInvert(out1, output);
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombine11(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                      unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                      unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  MInt shift, shift2;
  if(input1->Id1 != input2->Id1) {
    mapZero(output);
  } else {
    mapNormalize1(input1, output);
    mapNormalize1(input2, tmp);
    for(MInt i = 0; i < nDim; i++) {
      shift = 0;
      while(output->start1[i] + shift * output->step1[i] <= output->end1[i]) {
        if((output->start1[i] + shift * output->step1[i] >= tmp->start1[i])
           && (((output->start1[i] + shift * output->step1[i]) - (tmp->start1[i])) % (tmp->step1[i])) == 0)
          break;
        shift = shift + 1;
      }
      output->start1[i] = output->start1[i] + shift * output->step1[i];
      output->start2[output->order[i]] = output->start2[output->order[i]] + shift * output->step2[output->order[i]];

      shift = 0;

      while((output->end1[i] - shift * output->step1[i]) >= output->start1[i]) {
        if(((output->end1[i] - shift * output->step1[i]) <= tmp->end1[i])
           && ((output->end1[i] - shift * output->step1[i] - tmp->end1[i]) % tmp->step1[i]) == 0)
          break;
        ++shift;
      }

      output->end1[i] = output->end1[i] - shift * output->step1[i];
      output->end2[output->order[i]] = output->end2[output->order[i]] - shift * output->step2[output->order[i]];

      shift = 0;

      while(tmp->start1[i] + shift * tmp->step1[i] <= tmp->end1[i]) {
        if((tmp->start1[i] + shift * tmp->step1[i] >= output->start1[i])
           && ((tmp->start1[i] + (shift * tmp->step1[i]) - output->start1[i]) % output->step1[i]) == 0)
          break;
        ++shift;
      }

      tmp->start1[i] = tmp->start1[i] + shift * tmp->step1[i];
      tmp->start2[tmp->order[i]] = tmp->start2[tmp->order[i]] + shift * tmp->step2[tmp->order[i]];

      shift = 0;
      while(tmp->end1[i] - shift * tmp->step1[i] >= tmp->start1[i]) {
        if(((tmp->end1[i] - shift * tmp->step1[i]) <= output->end1[i])
           && ((tmp->end1[i] - shift * tmp->step1[i] - output->end1[i]) % output->step1[i]) == 0)
          break;
        ++shift;
      }
      tmp->end1[i] = tmp->end1[i] - shift * tmp->step1[i];
      tmp->end2[tmp->order[i]] = tmp->end2[tmp->order[i]] - shift * tmp->step2[tmp->order[i]];


      shift = 1;
      shift2 = 1;

      while(output->step1[i] * shift != tmp->step1[i] * shift2) {
        if(output->step1[i] * shift < tmp->step1[i] * shift2) {
          ++shift;
        } else {
          ++shift2;
        }
      }
      output->step1[i] = shift * output->step1[i];
      output->step2[output->order[i]] = shift * output->step2[output->order[i]];
      tmp->step1[i] = shift2 * tmp->step1[i];
      tmp->step2[tmp->order[i]] = shift2 * tmp->step2[tmp->order[i]];
    }

    unique_ptr<StructuredWindowMap<nDim>> tmp1 = make_unique<StructuredWindowMap<nDim>>();
    mapCpy(output, tmp1);
    mapInvert(tmp1, output);
    for(MInt i = 0; i < nDim; i++) {
      output->order[i] = tmp->order[output->order[i]];
    }
    output->Id2 = tmp->Id2;
    memcpy(output->start2, tmp->start2, nDim * sizeof(MInt));
    memcpy(output->end2, tmp->end2, nDim * sizeof(MInt));
    memcpy(output->step2, tmp->step2, nDim * sizeof(MInt));

    mapCpy(output, tmp1);
    mapNormalize1(tmp1, output);
    if(input1->BC == -1 && input2->BC > 0) {
      output->BC = input2->BC;
    } else {
      output->BC = input1->BC;
    }
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombine12(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                      unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                      unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  mapInvert(input2, tmp);
  mapCombine11(input1, tmp, output);
  // return tmp;
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombine21(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                      unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                      unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  mapInvert(input1, tmp);
  mapCombine11(tmp, input2, output);
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombine22(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                      unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                      unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  unique_ptr<StructuredWindowMap<nDim>> tmp1 = make_unique<StructuredWindowMap<nDim>>();
  mapInvert(input1, tmp);
  mapInvert(input2, tmp1);
  mapCombine11(tmp, tmp1, output);
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombineCell11(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                          unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                          unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  MInt shift, shift2;
  if(input1->Id1 != input2->Id1) {
    mapZero(output);
  } else {
    mapNormalize1(input1, output);
    mapNormalize1(input2, tmp);
    for(MInt i = 0; i < nDim; i++) {
      if(output->start1[i] == output->end1[i] || tmp->start1[i] == tmp->end1[i]) {
        shift = 0;
        while(output->start1[i] + shift * output->step1[i] <= output->end1[i]) {
          if((output->start1[i] + shift * output->step1[i] >= tmp->start1[i])
             && (((output->start1[i] + shift * output->step1[i]) - (tmp->start1[i])) % (tmp->step1[i])) == 0)
            break;
          shift = shift + 1;
        }
        output->start1[i] = output->start1[i] + shift * output->step1[i];
        output->start2[output->order[i]] = output->start2[output->order[i]] + shift * output->step2[output->order[i]];

        shift = 0;

        while((output->end1[i] - shift * output->step1[i]) >= output->start1[i]) {
          if(((output->end1[i] - shift * output->step1[i]) <= tmp->end1[i])
             && ((output->end1[i] - shift * output->step1[i] - tmp->end1[i]) % tmp->step1[i]) == 0)
            break;
          ++shift;
        }

        output->end1[i] = output->end1[i] - shift * output->step1[i];
        output->end2[output->order[i]] = output->end2[output->order[i]] - shift * output->step2[output->order[i]];

        shift = 0;

        while(tmp->start1[i] + shift * tmp->step1[i] <= tmp->end1[i]) {
          if((tmp->start1[i] + shift * tmp->step1[i] >= output->start1[i])
             && ((tmp->start1[i] + (shift * tmp->step1[i]) - output->start1[i]) % output->step1[i]) == 0)
            break;
          ++shift;
        }

        tmp->start1[i] = tmp->start1[i] + shift * tmp->step1[i];
        tmp->start2[tmp->order[i]] = tmp->start2[tmp->order[i]] + shift * tmp->step2[tmp->order[i]];

        shift = 0;
        while(tmp->end1[i] - shift * tmp->step1[i] >= tmp->start1[i]) {
          if(((tmp->end1[i] - shift * tmp->step1[i]) <= output->end1[i])
             && ((tmp->end1[i] - shift * tmp->step1[i] - output->end1[i]) % output->step1[i]) == 0)
            break;
          ++shift;
        }
        tmp->end1[i] = tmp->end1[i] - shift * tmp->step1[i];
        tmp->end2[tmp->order[i]] = tmp->end2[tmp->order[i]] - shift * tmp->step2[tmp->order[i]];
      }


      else {
        shift = 1;
        output->end1[i] = output->end1[i] - shift * output->step1[i];
        output->end2[output->order[i]] = output->end2[output->order[i]] - shift * output->step2[output->order[i]];
        tmp->end1[i] = tmp->end1[i] - shift * tmp->step1[i];
        tmp->end2[tmp->order[i]] = tmp->end2[tmp->order[i]] - shift * tmp->step2[tmp->order[i]];


        shift = 0;
        while(output->start1[i] + shift * output->step1[i] <= output->end1[i]) {
          if((output->start1[i] + shift * output->step1[i] >= tmp->start1[i])
             && (((output->start1[i] + shift * output->step1[i]) - (tmp->start1[i])) % (tmp->step1[i])) == 0)
            break;
          shift = shift + 1;
        }
        output->start1[i] = output->start1[i] + shift * output->step1[i];
        output->start2[output->order[i]] = output->start2[output->order[i]] + shift * output->step2[output->order[i]];

        shift = 0;

        while((output->end1[i] - shift * output->step1[i]) >= output->start1[i]) {
          if(((output->end1[i] - shift * output->step1[i]) <= tmp->end1[i])
             && ((output->end1[i] - shift * output->step1[i] - tmp->end1[i]) % tmp->step1[i]) == 0)
            break;
          ++shift;
        }

        output->end1[i] = output->end1[i] - shift * output->step1[i];
        output->end2[output->order[i]] = output->end2[output->order[i]] - shift * output->step2[output->order[i]];


        shift = 0;

        while(tmp->start1[i] + shift * tmp->step1[i] <= tmp->end1[i]) {
          if((tmp->start1[i] + shift * tmp->step1[i] >= output->start1[i])
             && ((tmp->start1[i] + (shift * tmp->step1[i]) - output->start1[i]) % output->step1[i]) == 0)
            break;
          ++shift;
        }

        tmp->start1[i] = tmp->start1[i] + shift * tmp->step1[i];
        tmp->start2[tmp->order[i]] = tmp->start2[tmp->order[i]] + shift * tmp->step2[tmp->order[i]];

        shift = 0;
        while(tmp->end1[i] - shift * tmp->step1[i] >= tmp->start1[i]) {
          if(((tmp->end1[i] - shift * tmp->step1[i]) <= output->end1[i])
             && ((tmp->end1[i] - shift * tmp->step1[i] - output->end1[i]) % output->step1[i]) == 0)
            break;
          ++shift;
        }
        tmp->end1[i] = tmp->end1[i] - shift * tmp->step1[i];
        tmp->end2[tmp->order[i]] = tmp->end2[tmp->order[i]] - shift * tmp->step2[tmp->order[i]];


        if(output->start1[i] <= output->end1[i] && tmp->start1[i] <= tmp->end1[i]) {
          shift = -1;
          output->end1[i] = output->end1[i] - shift * output->step1[i];
          output->end2[output->order[i]] = output->end2[output->order[i]] - shift * output->step2[output->order[i]];
          tmp->end1[i] = tmp->end1[i] - shift * tmp->step1[i];
          tmp->end2[tmp->order[i]] = tmp->end2[tmp->order[i]] - shift * tmp->step2[tmp->order[i]];
        }
      }

      shift = 1;
      shift2 = 1;

      while(output->step1[i] * shift != tmp->step1[i] * shift2) {
        if(output->step1[i] * shift < tmp->step1[i] * shift2) {
          ++shift;
        } else {
          ++shift2;
        }
      }
      output->step1[i] = shift * output->step1[i];
      output->step2[output->order[i]] = shift * output->step2[output->order[i]];
      tmp->step1[i] = shift2 * tmp->step1[i];
      tmp->step2[tmp->order[i]] = shift2 * tmp->step2[tmp->order[i]];
    }

    unique_ptr<StructuredWindowMap<nDim>> tmp1 = make_unique<StructuredWindowMap<nDim>>();
    mapCpy(output, tmp1);
    mapInvert(tmp1, output);
    for(MInt i = 0; i < nDim; i++) {
      output->order[i] = tmp->order[output->order[i]];
    }
    output->Id2 = tmp->Id2;
    memcpy(output->start2, tmp->start2, nDim * sizeof(MInt));
    memcpy(output->end2, tmp->end2, nDim * sizeof(MInt));
    memcpy(output->step2, tmp->step2, nDim * sizeof(MInt));

    mapCpy(output, tmp1);
    mapNormalize1(tmp1, output);
    if(input1->BC == -1 && input2->BC > 0) {
      output->BC = input2->BC;
    } else {
      output->BC = input1->BC;
    }
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombineWave(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                        unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                        unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  MInt shift;
  if(input1->Id1 != input2->Id1) {
    mapZero(output);
  } else {
    mapNormalize1(input1, output);
    mapNormalize1(input2, tmp);
    for(MInt i = 0; i < nDim; i++) {
      shift = 0;
      // go from map1.start1 to map1.end1 and find starting point of
      // overlapping part
      while(output->start1[i] + shift * output->step1[i] <= output->end1[i]) {
        // stop incrementing if map1.start1 >= map2.start2

        if((output->start1[i] + shift * output->step1[i] >= tmp->start1[i])
           && (((output->start1[i] + shift * output->step1[i]) - (tmp->start1[i])) % (tmp->step1[i])) == 0)
          break;
        shift = shift + 1;
      }

      output->start1[i] = output->start1[i] + shift * output->step1[i];

      output->start2[output->order[i]] = output->start2[output->order[i]] + shift * output->step2[output->order[i]];

      shift = 0;

      // go from map1.end1 to map1.start1 and find ending point of
      // overlapping part
      while((output->end1[i] - shift * output->step1[i]) >= output->start1[i]) {
        // stop incrementing if map1.end1 <= map2.end2
        if(((output->end1[i] - shift * output->step1[i]) <= tmp->end1[i])
           && ((output->end1[i] - shift * output->step1[i] - tmp->end1[i]) % tmp->step1[i]) == 0)
          break;
        ++shift;
      }

      output->end1[i] = output->end1[i] - shift * output->step1[i];

      output->end2[output->order[i]] = output->end2[output->order[i]] - shift * output->step2[output->order[i]];

      shift = 0;
    }

    output->Id1 = output->Id2;
    output->Id2 = tmp->Id2;
    output->BC = tmp->BC;
  }
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombineCell12(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                          unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                          unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  mapInvert(input2, tmp);
  mapCombineCell11(input1, tmp, output);
  // return tmp;
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombineCell21(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                          unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                          unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  mapInvert(input1, tmp);
  mapCombineCell11(tmp, input2, output);
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::mapCombineCell22(unique_ptr<StructuredWindowMap<nDim>>& input1,
                                                          unique_ptr<StructuredWindowMap<nDim>>& input2,
                                                          unique_ptr<StructuredWindowMap<nDim>>& output) {
  unique_ptr<StructuredWindowMap<nDim>> tmp = make_unique<StructuredWindowMap<nDim>>();
  unique_ptr<StructuredWindowMap<nDim>> tmp1 = make_unique<StructuredWindowMap<nDim>>();
  mapInvert(input1, tmp);
  mapInvert(input2, tmp1);
  mapCombineCell11(tmp, tmp1, output);
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::readMapFromArray(unique_ptr<StructuredWindowMap<nDim>>& map, MInt* dataField) {
  for(MInt dim = 0; dim < nDim; dim++) {
    map->start1[dim] = dataField[0 * nDim + dim];
    map->end1[dim] = dataField[1 * nDim + dim];
    map->step1[dim] = dataField[2 * nDim + dim];
    map->start2[dim] = dataField[3 * nDim + dim];
    map->end2[dim] = dataField[4 * nDim + dim];
    map->step2[dim] = dataField[5 * nDim + dim];
    map->order[dim] = dataField[6 * nDim + dim];
  }

  map->Id1 = dataField[7 * nDim + 0];
  map->Id2 = dataField[7 * nDim + 1];
  // map->nDim = dataField[7 * nDim + 2]; removed
  map->BC = dataField[7 * nDim + 3];
  // map->constIndex = dataField[7 * nDim + 4]; removed
  map->face = dataField[7 * nDim + 5];
  map->dir = dataField[7 * nDim + 6];
  map->dc1 = dataField[7 * nDim + 7];
  map->dc2 = dataField[7 * nDim + 8];
  map->originShape = dataField[7 * nDim + 9];
  map->hasSponge = dataField[7 * nDim + 10];
  map->spongeThickness = dataField[7 * nDim + 11];
  map->beta = dataField[7 * nDim + 12];
  map->sigma = dataField[7 * nDim + 13];
  map->Nstar = dataField[7 * nDim + 14];
  map->SingularId = dataField[7 * nDim + 15];
}

template <MInt nDim>
void FvStructuredSolverWindowInfo<nDim>::writeMapToArray(unique_ptr<StructuredWindowMap<nDim>>& map, MInt* dataField) {
  for(MInt dim = 0; dim < nDim; dim++) {
    dataField[0 * nDim + dim] = map->start1[dim];
    dataField[1 * nDim + dim] = map->end1[dim];
    dataField[2 * nDim + dim] = map->step1[dim];
    dataField[3 * nDim + dim] = map->start2[dim];
    dataField[4 * nDim + dim] = map->end2[dim];
    dataField[5 * nDim + dim] = map->step2[dim];
    dataField[6 * nDim + dim] = map->order[dim];
  }

  dataField[7 * nDim + 0] = map->Id1;
  dataField[7 * nDim + 1] = map->Id2;
  // dataField[7 * nDim + 2] = map->nDim; removed
  dataField[7 * nDim + 3] = map->BC;
  // dataField[7 * nDim + 4] = map->constIndex; removed
  dataField[7 * nDim + 5] = map->face;
  dataField[7 * nDim + 6] = map->dir;
  dataField[7 * nDim + 7] = map->dc1;
  dataField[7 * nDim + 8] = map->dc2;
  dataField[7 * nDim + 9] = map->originShape;
  dataField[7 * nDim + 10] = map->hasSponge;
  dataField[7 * nDim + 11] = map->spongeThickness;
  dataField[7 * nDim + 12] = map->beta;
  dataField[7 * nDim + 13] = map->sigma;
  dataField[7 * nDim + 14] = map->Nstar;
  dataField[7 * nDim + 15] = map->SingularId;
}

// Explicit instantiations for 2D and 3D
template class FvStructuredSolverWindowInfo<2>;
template class FvStructuredSolverWindowInfo<3>;
