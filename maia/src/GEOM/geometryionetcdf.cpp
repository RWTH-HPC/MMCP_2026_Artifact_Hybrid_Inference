// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometryionetcdf.h"

#include <cstring>
#include "COMM/mpioverride.h"
#include "IO/parallelio.h"
#include "UTIL/timer.h"
#include "globals.h"
using namespace std;

//--------------------------------------------------------------------------
/**
 *
 *
 */
MInt GeometryIONetcdf::segmentCount() { return m_noSegments; }

//--------------------------------------------------------------------------
/**
 *
 *
 */
void GeometryIONetcdf::makeProperty(GeometryProperty* p, MString name, ParallelIo* bdFile) {
  //    TRACE();
  using namespace maia::parallel_io;
  ParallelIo::size_type length;
  maiabd_type type;
  MInt noDims;
  type = bdFile->getDatasetType(name);
  if(bdFile->hasDataset(name, 0))
    noDims = 0;
  else
    noDims = bdFile->getDatasetNoDims(name);
  switch(type) {
    case PIO_INT: {
      DEBUG("GeometryIONetcdf::makeProperty found integer property :", MAIA_DEBUG_USER1);
      p->propertyType = MINT;
      if(noDims == 0) {
        p->elements = 1;
        p->intField = new MInt[1];
        bdFile->readScalar(p->intField, name);
        DEBUG("GeometryIONetcdf::makeProperty " << p->intField[0], MAIA_DEBUG_USER1);
      }
      if(noDims == 1) {
        /* look up the length of the array */
        length = bdFile->getArraySize(name);
        p->elements = length;
        p->intField = new MInt[length];
        bdFile->setOffset(length, 0);
        bdFile->readArray(p->intField, name);
        //                 for (MInt i=0; i < length; i++)     // delete this line!
        //                     DEBUG("GeometryIONetcdf::makeProperty " << p->intField[i], MAIA_DEBUG_USER1);
      }
      if(noDims == 2) {
        /*NOT YET IMPLEMENTED !*/
      }
      break;
    }
    case PIO_FLOAT: {
      DEBUG("GeometryIONetcdf::makeProperty found float property :", MAIA_DEBUG_USER1);
      p->propertyType = MFLOAT;
      if(noDims == 0) {
        p->elements = 1;
        p->floatField = new MFloat[1];
        bdFile->readScalar(p->floatField, name);
        DEBUG("GeometryIONetcdf::makeProperty " << p->floatField[0], MAIA_DEBUG_USER1);
      }
      if(noDims == 1) {
        /* look up the length of the array */
        length = bdFile->getArraySize(name);
        p->elements = length;
        p->floatField = new MFloat[length];
        bdFile->setOffset(length, 0);
        bdFile->readArray(p->floatField, name);
        //                 for (MInt i=0; i < length; i++)   //delete this line!
        //                     DEBUG("GeometryIONetcdf::makeProperty " << p->floatField[i], MAIA_DEBUG_USER1);
      }
      if(noDims == 2) {
        /*NOT YET IMPLEMENTED !*/
      }
      break;
    }
    case PIO_STRING: {
      DEBUG("GeometryIONetcdf::makeProperty found char property :", MAIA_DEBUG_USER1);
      p->propertyType = MSTRING;
      if(noDims == 0) {
        p->elements = 1;
        p->stringField = new MString[1];
        MString buf;
        bdFile->readScalar(&buf, name);
        p->stringField->append(buf);
      }
      if(noDims == 1) {
        p->elements = 1;
        p->stringField = new MString[1];
        /* look up the length of the string */
        length = bdFile->getArraySize(name, 0);
        MString buf;
        bdFile->setOffset(length, 0);
        bdFile->readArray(&buf, name);
        p->stringField->append(buf);
        DEBUG("GeometryIONetcdf::makeProperty " << p->stringField[0], MAIA_DEBUG_USER1);
      }
      if(noDims == 2) {
        /* look up the number of strings */
        ParallelIo::size_type noStrings = bdFile->getArraySize(name, 0);
        length = bdFile->getArraySize(name, 1);
        p->elements = length;
        p->stringField = new MString[noStrings];
        ParallelIo::size_type start;
        /* loop over number of strings */
        for(MInt i = 0; i < p->elements; i++) {
          start = i;
          bdFile->setOffset(1, start, 2);
          MString buf;
          bdFile->readArray(&buf, name);
          (p->stringField[i]).append(buf);
          DEBUG("GeometryIONetcdf::makeProperty " << p->stringField[i], MAIA_DEBUG_USER1);
        }
      }
      break;
    }
    default: {
      // SX8 cannot compile with excepetion handling
      // throw(IOError(" Error : Unsupported variable input type encountered! \n"));
    }
  }
  const pair<const MString, GeometryProperty*> mp(p->name, p);
  m_geometryPropertyMap->insert(mp);
  DEBUG("GeometryIONetcdf::makeProperty created default property ", MAIA_DEBUG_USER1);
  DEBUG("GeometryIONetcdf::makeProperty elements = " << p->elements, MAIA_DEBUG_USER1);
  DEBUG("GeometryIONetcdf::makeProperty m_noProperties = " << m_geometryPropertyMap->size(), MAIA_DEBUG_USER1);
}


void GeometryIONetcdf::readBodiesNewIOMethod(ParallelIo* bdFile) {
  TRACE();
  m_log << "      * reading body information" << endl;

  if(bdFile->hasDataset("body_segments.") && bdFile->hasDataset("body_segments_names.")
     && bdFile->hasDataset("body_segments_num.")) {
    MInt len_bsnum = bdFile->getArraySize("body_segments_num.");
    MInt len_bs = bdFile->getArraySize("body_segments.");
    MInt len_bsname = bdFile->getArraySize("body_segments_names.");

    MString bsallnames;
    MIntScratchSpace bs(len_bs, AT_, "bs");
    MStringScratchSpace bsname(len_bsnum, AT_, "bsname");
    MIntScratchSpace bsnum(len_bsnum, AT_, "bsnum");

    bdFile->setOffset(len_bs, 0);
    bdFile->readArray(bs.getPointer(), "body_segments.");

    bdFile->setOffset(len_bsname, 0);
    bdFile->readArray(&bsallnames, "body_segments_names.");
    bdFile->setOffset(len_bsnum, 0);
    bdFile->readArray(bsnum.getPointer(), "body_segments_num.");

    for(MInt i = 0; i < len_bsnum; i++) {
      MInt del_pos = bsallnames.find(",");
      bsname[i] = bsallnames.substr(0, del_pos);
      MString tmp = bsallnames.substr(del_pos + 1, bsallnames.length() - del_pos);
      bsallnames = tmp;
    }

    for(MInt i = 0, j = 0; i < len_bsnum; i++) {
      Body* body = new Body;
      m_noBodies++;

      body->name = bsname[i];
      body->noSegments = bsnum[i];
      body->segments = new MInt[body->noSegments];

      for(MInt k = 0; k < body->noSegments; k++, j++)
        body->segments[k] = bs[j];

      if(m_bodyMap->find(body->name) == m_bodyMap->end()) m_bodyMap->insert(make_pair(body->name, body));
    }
  }

  buildDefaultBody();
}

//--------------------------------------------------------------------------
/**
 *
 *
 */
void GeometryIONetcdf::readBodiesOldIOMethod(ParallelIo* bdFile) {
  TRACE();

  m_log << "      * reading body information" << endl;
  MInt noVariables = 0;
  const char MPropertySeperator = '.'; // Seperator that defines different Properties

  /* read and store the body information and the segment info*/
  const char* du; // tmp variable

  MString varName;
  ParallelIo::size_type length;

  vector<MString> varNames = bdFile->getDatasetNames(1);
  noVariables = varNames.size();

  for(MInt n = 0; n < noVariables; n++) {
    varName = varNames[n];
    if(strstr(varName.c_str(), "body_segments.")) {         // if a body has been identified
      du = strrchr(varName.c_str(), MPropertySeperator);    // get pointer to where the "." is

      /* create body element */
      Body* body;
      body = new Body;
      body->name.append(++du);
      DEBUG("GeometryIONetcdf::readBodies Body found (name: " << body->name << ")", MAIA_DEBUG_USER2);
      m_noBodies++;

      /* get segments for the body */
      length = bdFile->getArraySize(varName);
      DEBUG("GeometryIONetcdf::readBodies The Body has " << length << " solvers.", MAIA_DEBUG_USER2);
      body->noSegments = length;
      body->segments = new MInt[length];
      /* There is an implicit conversion from int to unsigned int between netcdf
         and the body segments */
      for(ParallelIo::size_type i = 0; i < length; i++) {
        MInt dsegments;
        bdFile->setOffset(1, i);
        bdFile->readArray(&dsegments, varName);
        body->segments[i] = dsegments;
      }
      /* insert the body into the bodyMap*/
      if(m_bodyMap->find(body->name) == m_bodyMap->end()) {
        m_bodyMap->insert(make_pair(body->name, body));
      }
    }
  }
  buildDefaultBody();
}


//--------------------------------------------------------------------------
/** This function creates a default body for all segments, that do not appear
 * in a body definition.
 *
 */
void GeometryIONetcdf::buildDefaultBody() {
  TRACE();
  /* create a body named "default" */
  Body* body;
  body = new Body;
  body->name.append("default");

  /* search for all segments that are defined */
  list<MInt> segmentList;
  for(bodyMap::const_iterator it = m_bodyMap->begin(); it != m_bodyMap->end(); it++) {
    for(MInt i = 0; i < it->second->noSegments; i++) {
      segmentList.push_back(it->second->segments[i]);
    }
  }
  segmentList.sort();

  /* The number of undefined segments equals the number of total
     segments minus the number of defined segments ... */
  body->noSegments = (m_noSegments - segmentList.size());
  body->segments = new MInt[body->noSegments];

  /* this adds all segmentId's that aren't defined to the default body*/
  MInt i = 0;
  MInt simpleIterator = 0;
  list<MInt>::const_iterator segmentIt = segmentList.begin();

  if(!segmentList.empty()) {
    while(i < m_noSegments) {
      if(*segmentIt == i) {
        i++;
        segmentIt++;
      } else {
        body->segments[simpleIterator] = i;
        DEBUG("GeometryIONetcdf::buildDefaultBody Added segment " << i << " to the default body. ", MAIA_DEBUG_USER1);
        i++;
        simpleIterator++;
      }
    }
  }

  /* insert the default body in the body map */
  m_bodyMap->insert(make_pair(body->name, body));
  DEBUG("GeometryIONetcdf::buildDefaultBody  Default body has " << body->noSegments << " segments.", MAIA_DEBUG_USER1);
  m_noBodies++;
}

//--------------------------------------------------------------------------
/** \brief check if the geometry property file is of new or old type and calls the according function
 *
 * \author Andreas Lintermann
 * \date 06.01.2016
 *
 * Checks for a global attribute in the geometry property file named "newIOMethod"
 * and calls the according function
 *
 * \param[in] name the name of the file to open
 *
 **/
geometryAssembly* GeometryIONetcdf::readPropertyFile(MString name) {
  TRACE();

  // init member variables
  m_newIOMethod = false;
  m_noSegments = 0;
  m_noBodies = 0;
  m_geometryPropertyMap = new geometryPropertyMap;
  m_bodyMap = new bodyMap;
  m_geometryAssembly = new geometryAssembly;

  // read initial variables required for further processing
  if(domainId() == 0) {
    ParallelIo parallelIo(name.c_str(), maia::parallel_io::PIO_READ, MPI_COMM_SELF);
    m_newIOMethod = parallelIo.hasAttribute("newIOMethod");
    parallelIo.readScalar(&m_noSegments, "noSegments");

    MPI_Bcast(&m_newIOMethod, 1, MPI_INT, 0, mpiComm(), AT_, "m_newIOMethod");
    MPI_Bcast(&m_noSegments, 1, MPI_INT, 0, mpiComm(), AT_, "m_noSegments");
  } else {
    MPI_Bcast(&m_newIOMethod, 1, MPI_INT, 0, mpiComm(), AT_, "m_newIOMethod");
    MPI_Bcast(&m_noSegments, 1, MPI_INT, 0, mpiComm(), AT_, "m_noSegments");
  }

  // do the actual reading and distribute the geometry properties
  if(domainId() == 0) {
    m_log << "    - rank 0 reads data from disk" << endl;

    // open the property file
    ParallelIo parallelIo(name.c_str(), maia::parallel_io::PIO_READ, MPI_COMM_SELF);

    NEW_SUB_TIMER(t_readRest, "rest geometry property", g_t_readGeomFile);
    NEW_SUB_TIMER(t_readBodies, "geometry bodies", g_t_readGeomFile);
    NEW_SUB_TIMER(t_distBodies, "distribute bodies", g_t_readGeomFile);
    NEW_SUB_TIMER(t_distRest, "distribute rest geometry property", g_t_readGeomFile);

    if(m_newIOMethod) {
      RECORD_TIMER_START(t_readBodies);
      readBodiesNewIOMethod(&parallelIo);
      RECORD_TIMER_STOP(t_readBodies);

      RECORD_TIMER_START(t_readRest);
      readPropertyFileNewIOMethod(&parallelIo);
      RECORD_TIMER_STOP(t_readRest);
    } else {
      RECORD_TIMER_START(t_readBodies);
      readBodiesOldIOMethod(&parallelIo);
      RECORD_TIMER_STOP(t_readBodies);

      RECORD_TIMER_START(t_readRest);
      readPropertyFileOldIOMethod(&parallelIo);
      RECORD_TIMER_STOP(t_readRest);
    }

    /*
    for(bodyIterator it = m_bodyMap->begin(); it != m_bodyMap->end(); it++)
{
  Body* b = (*it).second;
  cout << "body: " << b->name << " " << b->noSegments << " ";
  for(MInt i = 0; i < b->noSegments; i++)
    cout << b->segments[i] << " ";
  cout << endl;
}

    cout << endl;
    for(geometryPropertyIterator it = m_geometryPropertyMap->begin(); it != m_geometryPropertyMap->end(); it++)
      {
        GeometryProperty* b = (*it).second;
        cout << "prop: " << b->name << " " << b->propertyType << " " << b->elements << " " << b->segmentId << " - ";
        for(MInt i = 0; i < b->elements; i++)
    {
      if(b->propertyType == MINT)
  cout << b->intField[i] << " ";
      else if(b->propertyType == MFLOAT)
  cout << b->floatField[i] << " ";
      else if (b->propertyType == MSTRING)
  cout << b->stringField[i] << " ";

    }
  cout << endl;
}
    */
    if(noDomains() > 1) {
      RECORD_TIMER_START(t_distBodies);
      distributeBodyProperties();
      RECORD_TIMER_STOP(t_distBodies);

      RECORD_TIMER_START(t_distRest);
      distributeGeometryProperties();
      RECORD_TIMER_STOP(t_distRest);
    }
  } else {
    receiveBodyProperties();
    receiveGeometryProperties();

    /*
    for(bodyIterator it = m_bodyMap->begin(); it != m_bodyMap->end(); it++)
{
  Body* b = (*it).second;
  cout << "body: " << b->name << " " << b->noSegments << " ";
  for(MInt i = 0; i < b->noSegments; i++)
    cout << b->segments[i] << " ";
  cout << endl;
}

    cout << endl;
    for(geometryPropertyIterator it = m_geometryPropertyMap->begin(); it != m_geometryPropertyMap->end(); it++)
      {
        GeometryProperty* b = (*it).second;
        cout << "prop: " << b->name << " " << b->propertyType << " " << b->elements << " " << b->segmentId << " - ";
        for(MInt i = 0; i < b->elements; i++)
    {
      if(b->propertyType == MINT)
  cout << b->intField[i] << " ";
      else if(b->propertyType == MFLOAT)
  cout << b->floatField[i] << " ";
      else if (b->propertyType == MSTRING)
  cout << b->stringField[i] << " ";

    }
  cout << endl;
}
    */
  }

  // fill the body and the property maps
  m_geometryAssembly->geometryProperties = m_geometryPropertyMap;
  m_geometryAssembly->bodies = m_bodyMap;

  m_log << endl;

  return m_geometryAssembly;
}

/** \brief reads in the geomertry property file the old way
 *
 * \author Andreas Lintermann
 * \date 05.01.2016
 *
 * Reads in the geometry property file on rank 0 and distributes the
 * information among the processes.
 *
 * \param[in] parallelIo the parallelIo handle
 *
 **/
void GeometryIONetcdf::readPropertyFileOldIOMethod(ParallelIo* parallelIo) {
  TRACE();

  MInt noVariables = 0;
  const char MPropertySeperator = '.'; // Seperator that defines different Properties

  MString varName;

  GeometryProperty* p;
  m_log << "      * reading rest of the properties" << endl;

  // check the consistency of the body
  if(checkBodyConsistency()) {
    vector<MString> varNames = parallelIo->getDatasetNames();
    noVariables = varNames.size();

    for(MInt id = 0; id < noVariables; id++) {
      varName = varNames[id];
      // if default property
      if(!strstr(varName.c_str(), ".")) {
        // if there's no dot in string (default prop)
        // store default property in last solver,
        // the others to their belonging id's

        DEBUG("GeometryIONetcdf::readNCPropertyFile  default property : " << varName, MAIA_DEBUG_USER1);
        p = new GeometryProperty;
        p->segmentId = m_noSegments; // Insert default segment as last segment
        p->name.append(varName);
        makeProperty(p, varName, parallelIo); // create new Property
      }
      // if a property defined for one or more bodies
      else {
        if(strstr(varName.c_str(), "_bodies.")) {
          DEBUG("GeometryIONetcdf::readPropertyFile normal property: " << varName, MAIA_DEBUG_USER1);
          ParallelIo::size_type noDims;
          noDims = parallelIo->getDatasetNoDims(varName);
          MString* bodies = nullptr;
          MInt noBodies = 1;
          DEBUG("GeometryIONetcdf::readPropertyFile no of dimensions = " << noDims, MAIA_DEBUG_USER1);
          switch(noDims) {
            case 0: {
              // if only one char (and one body)
              MString buf;
              parallelIo->readScalar(&buf, varName);
              bodies = new MString(buf);
              break;
            }
            case 1: {
              // look up the length of the string
              ParallelIo::size_type length = parallelIo->getArraySize(varName, 0);
              MString buf;
              parallelIo->setOffset(length, 0);
              parallelIo->readArray(&buf, varName);
              bodies = new MString(buf);
              break;
            }
            case 2: {
              // look up the number of bodies
              ParallelIo::size_type dbodies;
              dbodies = parallelIo->getArraySize(varName, 0);
              noBodies = (MInt)dbodies;
              bodies = new MString[noBodies];
              ParallelIo::size_type start;
              // loop over number of bodies
              for(MInt i = 0; i < noBodies; i++) {
                start = i;
                parallelIo->setOffset(1, start, 2);
                MString buf;
                parallelIo->readArray(&buf, varName);
                (bodies[i]).append(buf);
                DEBUG("GeometryIONetcdf::readPropertyFile " << bodies[i], MAIA_DEBUG_USER1);
              }
              break;
            }
            default: {
            }
          }
          // find all segments for the property

          list<MInt> segmentList;
          for(MInt i = 0; i != noBodies; i++) {
            DEBUG("GeometryIONetcdf::readPropertyFile definition for body " << bodies[i], MAIA_DEBUG_USER1);
            bodyIterator zI;
            // look for the body in the bodyMap
            zI = m_bodyMap->find(bodies[i]);
            // append all segments of the body to the segmentlist
            for(MInt j = 0; j < zI->second->noSegments; j++)
              segmentList.push_back(zI->second->segments[j]);
          }

          char* du;
          du = strrchr(const_cast<MChar*>(varName.c_str()), MPropertySeperator);
          MString dummy(varName);
          // strip the varname of "_bodies.1"
          dummy.replace(dummy.find("_bodies."), dummy.size(), du);
          DEBUG("GeometryIONetcdf::readPropertyFile  found property : " << dummy, MAIA_DEBUG_USER1);
          // find the property id

          list<MInt>::const_iterator it = segmentList.begin();
          dummy.erase(dummy.find(".")); // erase the dot and the id
          // create the property for every segment
          for(; it != segmentList.end(); it++) {
            p = new GeometryProperty;
            p->segmentId = *it;
            p->name.append(dummy);
            makeProperty(p, dummy, parallelIo); // create new Property
            DEBUG("GeometryIONetcdf::readPropertyFile created property for solver " << *it, MAIA_DEBUG_USER1);
          }
        } // end of if ( strstr ( varName,"_bodies." ) )
        else {
          // Determine the segmentId
          char* du;
          du = strrchr(const_cast<MChar*>(varName.c_str()), MPropertySeperator) + 1;
          MInt singleSegmentId = atoi(du);

          MString dummyName(varName);

          if(singleSegmentId || *du == '0') {
            DEBUG("Found single segment property definition for segment " << singleSegmentId, MAIA_DEBUG_IO);
            p = new GeometryProperty;
            p->segmentId = singleSegmentId;

            MString dummy = dummyName;

            dummyName.erase(dummyName.find(".")); // erase the dot and the id
            p->name.append(dummyName);
            makeProperty(p, dummy, parallelIo); // create new Property
          }
        }
      }
    }
    // Insert here the check for the geometry property consistency
    if(checkGeometryPropertyConsistency()) {
      DEBUG("GeometryIONetcdf::readPropertyFile ** Property check successful \n", MAIA_DEBUG_USER1);
    }
  }
}

/** \brief reads in the geomertry property file the old way
 *
 * \author Andreas Lintermann
 * \date 05.01.2016
 *
 * Reads in the geometry property file on rank 0 and distributes the
 * information among the processes.
 *
 * \param[in] name the name of the file to open
 *
 **/
void GeometryIONetcdf::readPropertyFileNewIOMethod(ParallelIo* parallelIo) {
  TRACE();

  m_log << "      * reading rest of the properties" << endl;

  if(checkBodyConsistency()) {
    vector<MString> varNames = parallelIo->getDatasetNames();
    MInt noVariables = varNames.size();

    for(MInt id = 0; id < noVariables; id++) {
      MString varName = varNames[id];

      if(!strstr(varName.c_str(), ".")) {
        GeometryProperty* p = new GeometryProperty;
        p->segmentId = m_noSegments;
        p->name.append(varName);
        makeProperty(p, varName, parallelIo);
      } else {
        if(strstr(varName.c_str(), "BC.")) {
          // create default
          GeometryProperty* p_def = new GeometryProperty;
          p_def->propertyType = MINT;
          p_def->segmentId = m_noSegments;
          p_def->name.append("BC");
          p_def->elements = 1;
          p_def->intField = new MInt[1];
          p_def->intField[0] = 0;

          const pair<const MString, GeometryProperty*> mp_def(p_def->name, p_def);
          m_geometryPropertyMap->insert(mp_def);

          // read all others
          MInt len = parallelIo->getArraySize(varName.c_str());
          MIntScratchSpace bcs(len, AT_, "bcs");
          parallelIo->setOffset(len, 0);
          parallelIo->readArray(bcs.getPointer(), varName.c_str());

          // create a property for each element with increasing segmentId
          for(MInt i = 0; i < len; i++) {
            GeometryProperty* p = new GeometryProperty;
            p->propertyType = MINT;
            p->segmentId = i;
            p->name.append("BC");
            p->elements = 1;
            p->intField = new MInt[1];
            p->intField[0] = bcs[i];

            const pair<const MString, GeometryProperty*> mp(p->name, p);
            m_geometryPropertyMap->insert(mp);
          }
        } else if(strstr(varName.c_str(), "filename.")) {
          // create default
          GeometryProperty* p_def = new GeometryProperty;
          p_def->propertyType = MSTRING;
          p_def->segmentId = m_noSegments;
          p_def->name.append("filename");
          p_def->elements = 1;
          p_def->stringField = new MString[1];
          p_def->stringField[0] = "";

          const pair<const MString, GeometryProperty*> mp_def(p_def->name, p_def);
          m_geometryPropertyMap->insert(mp_def);

          MString allnames;
          MInt lenAll = parallelIo->getArraySize("filename.");
          parallelIo->setOffset(lenAll, 0);
          parallelIo->readArray(&allnames, "filename.");

          MInt num = count(allnames.begin(), allnames.end(), ',') + 1;
          MStringScratchSpace name(num, AT_, "name");
          for(MInt i = 0; i < num; i++) {
            MInt del_pos = allnames.find(",");
            name[i] = allnames.substr(0, del_pos);
            MString tmp = allnames.substr(del_pos + 1, allnames.length() - del_pos);
            allnames = tmp;
          }

          // create a property for each element with increasing segmentId
          for(MInt i = 0; i < num; i++) {
            GeometryProperty* p = new GeometryProperty;
            p->propertyType = MSTRING;
            p->segmentId = i;
            p->name.append("filename");
            p->elements = 1;
            p->stringField = new MString[1];
            p->stringField[0] = name[i];

            const pair<const MString, GeometryProperty*> mp(p->name, p);
            m_geometryPropertyMap->insert(mp);
          }
        }
      }
    }
  }
}

/** \brief distributes the read geometry properties under all processes
 *
 * \author Andreas Lintermann
 * \date 05.01.2016
 *
 * This algorithm does the following:
 *
 *   1.  communicate the number of properties
 *   2.  run over all properties
 *   2.1 communicate the name
 *   2.2 communicate the int data, i.e., propertyType, elements, segmentId
 *   2.3 communicate the arrays
 *
 **/
void GeometryIONetcdf::distributeGeometryProperties() {
  TRACE();

  m_log << "    - rank 0 distributes the geometry property information" << endl;
  // 1. communicate the number of properties
  MInt propsize = m_geometryPropertyMap->size();
  MPI_Bcast(&propsize, 1, MPI_INT, 0, mpiComm(), AT_, "propsize");

  // 2. run over all properties
  for(geometryPropertyMap::iterator it = m_geometryPropertyMap->begin(); it != m_geometryPropertyMap->end(); ++it) {
    GeometryProperty* prop = (*it).second;

    // 2.1 communicate the name
    MInt name_length = prop->name.length();
    MCharScratchSpace c_name(name_length, AT_, "c_name");
    strcpy(c_name.begin(), prop->name.c_str());

    MPI_Bcast(&name_length, 1, MPI_INT, 0, mpiComm(), AT_, "name_length");
    MPI_Bcast(c_name.begin(), name_length, MPI_CHAR, 0, mpiComm(), AT_, "c_name.begin()");

    // 2.2 communicate the int data, i.e., propertyType, elements, segmentId
    MIntScratchSpace int_data(3, AT_, "int_data");
    int_data[0] = prop->propertyType;
    int_data[1] = prop->elements;
    int_data[2] = prop->segmentId;

    MPI_Bcast(int_data.getPointer(), 3, MPI_INT, 0, mpiComm(), AT_, "int_data.getPointer()");

    // 2.3 communicate the arrays
    switch(prop->propertyType) {
      case MINT: {
        MPI_Bcast(prop->intField, prop->elements, MPI_INT, 0, mpiComm(), AT_, "prop->intField");
        break;
      }
      case MFLOAT: {
        MPI_Bcast(prop->floatField, prop->elements, MPI_DOUBLE, 0, mpiComm(), AT_, "prop->floatField");
        break;
      }
      case MSTRING: {
        MIntScratchSpace off(prop->elements + 1, AT_, "off");
        MInt fld_len = 0;
        for(MInt i = 0; i < prop->elements; i++) {
          off[i] = fld_len;
          fld_len += prop->stringField[i].length();
        }
        off[prop->elements] = fld_len;

        char* buf = new char[fld_len];
        for(MInt i = 0; i < prop->elements; i++)
          strncpy(&buf[off[i]], prop->stringField[i].c_str(), prop->stringField[i].length());

        MPI_Bcast(&fld_len, 1, MPI_INT, 0, mpiComm(), AT_, "fld_len");
        MPI_Bcast(off.getPointer(), prop->elements + 1, MPI_INT, 0, mpiComm(), AT_, "off.getPointer()");
        MPI_Bcast(buf, fld_len, MPI_CHAR, 0, mpiComm(), AT_, "buf");
        delete[] buf;
        break;
      }
      default: {
      }
    }
  }
}


/** \brief receives the geometry properties from rank 0
 *
 * \author Andreas Lintermann
 * \date 05.01.2016
 *
 * This algorithm does the following:
 *
 *   1.  receive the number of properties
 *   2.  run over all properties
 *   2.1 receive the name
 *   2.2 receive the int data, i.e., propertyType, elements, segmentId
 *   2.3 receive the arrays
 *   2.4 inster the property into the property map
 *
 **/
void GeometryIONetcdf::receiveGeometryProperties() {
  TRACE();

  // 1. receive the number of properties
  MInt propsize = 0;
  MPI_Bcast(&propsize, 1, MPI_INT, 0, mpiComm(), AT_, "propsize");

  // 2. run over all properties
  for(MInt p = 0; p < propsize; p++) {
    GeometryProperty* prop = new GeometryProperty;

    // 2.1 receive the name
    MInt name_length = 0;
    MPI_Bcast(&name_length, 1, MPI_INT, 0, mpiComm(), AT_, "name_length");
    char* name = new char[name_length + 1];
    name[name_length] = '\0';
    MPI_Bcast(name, name_length, MPI_CHAR, 0, mpiComm(), AT_, "name");
    prop->name = name;

    // 2.2 receive the int data, i.e., propertyType, elements, segmentId
    MIntScratchSpace int_data(3, AT_, "int_data");
    MPI_Bcast(int_data.getPointer(), 3, MPI_INT, 0, mpiComm(), AT_, "int_data.getPointer()");
    prop->propertyType = (VariableType)int_data[0];
    prop->elements = int_data[1];
    prop->segmentId = int_data[2];

    // 2.3 receive the arrays
    switch(prop->propertyType) {
      case MINT: {
        prop->intField = new MInt[prop->elements];
        MPI_Bcast(prop->intField, prop->elements, MPI_INT, 0, mpiComm(), AT_, "prop->intField");
        break;
      }
      case MFLOAT: {
        prop->floatField = new MFloat[prop->elements];
        MPI_Bcast(prop->floatField, prop->elements, MPI_DOUBLE, 0, mpiComm(), AT_, "prop->floatField");
        break;
      }
      case MSTRING: {
        prop->stringField = new MString[prop->elements];

        MIntScratchSpace off(prop->elements + 1, AT_, "off");
        MInt fld_len = 0;
        MPI_Bcast(&fld_len, 1, MPI_INT, 0, mpiComm(), AT_, "fld_len");
        MPI_Bcast(off.getPointer(), prop->elements + 1, MPI_INT, 0, mpiComm(), AT_, "off.getPointer()");

        char* buf = new char[fld_len];
        MPI_Bcast(buf, fld_len, MPI_CHAR, 0, mpiComm(), AT_, "buf");

        for(MInt i = 0; i < prop->elements; i++) {
          MInt charsize = off[i + 1] - off[i];
          char* tmp = new char[charsize + 1];
          strncpy(tmp, &buf[off[i]], charsize);
          tmp[charsize] = '\0';
          MString st(tmp);
          prop->stringField[i] = st;
          delete[] tmp;
        }
        delete[] buf;
        break;
      }
      default: {
      }
    }

    // 2.4 insert the property into the property map
    const pair<const MString, GeometryProperty*> mp(prop->name, prop);
    m_geometryPropertyMap->insert(mp);
  }
}

/** \brief distributes the body information under all processes
 *
 * \author Andreas Lintermann
 * \date 05.01.2016
 *
 * This algorithm does the following:
 *
 *   1.  communicate the number of bodies
 *   2.  run over all bodies
 *   2.1 communicate the name
 *   2.2 communicate the int data
 *
 **/
void GeometryIONetcdf::distributeBodyProperties() {
  TRACE();

  m_log << "    - rank 0 distributes the body information" << endl;

  // 1. communicate the number of bodies
  MInt bosize = m_bodyMap->size();
  MPI_Bcast(&bosize, 1, MPI_INT, 0, mpiComm(), AT_, "bosize");

  // 2. run over all bodies
  for(bodyMap::iterator it = m_bodyMap->begin(); it != m_bodyMap->end(); ++it) {
    Body* body = (*it).second;

    // 2.1 communicate the name
    MInt name_length = body->name.length();
    MCharScratchSpace c_name(name_length, AT_, "c_name");
    strcpy(c_name.begin(), body->name.c_str());

    MPI_Bcast(&name_length, 1, MPI_INT, 0, mpiComm(), AT_, "name_length");
    MPI_Bcast(c_name.begin(), name_length, MPI_CHAR, 0, mpiComm(), AT_, "c_name.begin()");

    // 2.2 communicate the int data
    MPI_Bcast(&body->noSegments, 1, MPI_INT, 0, mpiComm(), AT_, "body->noSegments");
    MPI_Bcast(body->segments, body->noSegments, MPI_INT, 0, mpiComm(), AT_, "body->segments");
  }
}

/** \brief receives the body information from rank 0
 *
 * \author Andreas Lintermann
 * \date 05.01.2016
 *
 * This algorithm does the following:
 *
 *   1.  receive the number of bodies
 *   2.  run over all bodies
 *   2.1 receive the name
 *   2.2 receive the int data
 *   2.3 insert body into body map
 *
 **/
void GeometryIONetcdf::receiveBodyProperties() {
  TRACE();

  // 1. receive the number of bodies
  MInt bosize = 0;
  MPI_Bcast(&bosize, 1, MPI_INT, 0, mpiComm(), AT_, "bosize");

  // 2. run over all bodies
  for(MInt p = 0; p < bosize; p++) {
    Body* body = new Body;
    m_noBodies++;

    // 2.1 receive the name
    MInt name_length = 0;
    MPI_Bcast(&name_length, 1, MPI_INT, 0, mpiComm(), AT_, "name_length");
    char* name = new char[name_length + 1];
    name[name_length] = '\0';
    MPI_Bcast(name, name_length, MPI_CHAR, 0, mpiComm(), AT_, "name");
    body->name = name;

    // 2.2 receive the int data
    MPI_Bcast(&body->noSegments, 1, MPI_INT, 0, mpiComm(), AT_, "body->noSegments");
    body->segments = new MInt[body->noSegments];
    MPI_Bcast(body->segments, body->noSegments, MPI_INT, 0, mpiComm(), AT_, "body->segments");

    // 2.3 insert body into body map
    m_bodyMap->insert(make_pair(body->name, body));
  }

  buildDefaultBody();
}

//--------------------------------------------------------------------------------
/**
 *
 *
 */
MBool GeometryIONetcdf::checkGeometryPropertyConsistency() {
  TRACE();

  for(geometryPropertyIterator i = m_geometryPropertyMap->begin(); i != m_geometryPropertyMap->end(); i++) {
    /* if default property exists, then take next property */

    if(m_geometryPropertyMap->lower_bound(i->second->name)->second->segmentId == m_noSegments) {
      DEBUG("GeometryIONetcdf::checkPropertyConsistency default property exists for :" << i->second->name,
            MAIA_DEBUG_USER1);
    } else {
      return false;
    }
  }

  return true;
}

//--------------------------------------------------------------------------
/** Checks the consistency of the bodies, by creating a list, that contains
 *  all integer values from 0 to the number of the segments, and comparing it
 *  to a list of all actual registered segmentId's. For a consistent bodyList
 *  both lists must be equal.
 */
MBool GeometryIONetcdf::checkBodyConsistency() {
  TRACE();
  list<MInt> segmentList;
  list<MInt> compareList;
  MInt index = 0;
  for(bodyMap::const_iterator it = m_bodyMap->begin(); it != m_bodyMap->end(); it++) {
    for(MInt i = 0; i < it->second->noSegments; i++) {
      //       cerr << it->second->segments[i] << endl;
      //       cerr << index << endl;
      segmentList.push_back(it->second->segments[i]);
      compareList.push_back(index++);
    }
  }
  segmentList.sort();
  if(segmentList == compareList) {
    return true;
  } else {
    return false;
  }
}

//---------------------------------------------------------------------------
/** \brief Write the properties into a netcdf file
 *
 * This function write all properties in the propertyMap pMap to
 * a specified file.
 */

void GeometryIONetcdf::writeProperties(const MChar* fileName, geometryPropertyMap* pMap) {
  TRACE();

  // WARNING: untested switch from NetCDF/Parallel netCDF to ParallelIo
  // The method previously used direct I/O calls, which were replaced by
  // ParallelIo methods in summer 2015. However, since the method was not
  // used by any of the testcases, this code is still *untested*. Thus,
  // if your code uses this part of the code, please make sure that the
  // I/O still works as expected and then remove this warning as well as
  // the subsequent TERMM().
  TERMM(1, "untested I/O method, please see comment for how to proceed");
  ParallelIo parallelIo(fileName, maia::parallel_io::PIO_REPLACE, MPI_COMM_SELF);
  /*
   *
   * ParallelIo define solver --->
   *
   */
  MString dimName;
  for(geometryPropertyMap::iterator i = pMap->begin(); i != pMap->end(); i++) {
    dimName = i->first;
    dimName.append("Dim");
    switch(i->second->type()) {
      case MINT: {
        parallelIo.defineArray(maia::parallel_io::PIO_INT, i->first, i->second->count());
        break;
      }

      case MFLOAT: {
        parallelIo.defineArray(maia::parallel_io::PIO_FLOAT, i->first, i->second->count());
        break;
      }

      case MSTRING: {
        ParallelIo::size_type totalCount[2];
        totalCount[0] = i->second->count();
        totalCount[1] = NC_MAX_NAME;
        parallelIo.defineArray(maia::parallel_io::PIO_STRING, i->first, 2, totalCount);
        break;
      }

      default: {
        mTerm(1, AT_, "Uknown property type");
      }
    }
  }

  /*
   *
   * <----- ParallelIo define solver
   *
   */
  for(geometryPropertyMap::iterator i = pMap->begin(); i != pMap->end(); i++) {
    switch(i->second->type()) {
      case MINT: {
        parallelIo.setOffset(i->second->count(), 0);
        parallelIo.writeArray(i->second->intField, i->first);
        break;
      }

      case MFLOAT: {
        parallelIo.setOffset(i->second->count(), 0);
        parallelIo.writeArray(i->second->floatField, i->first);
        break;
      }

      case MSTRING: {
        parallelIo.setOffset(i->second->count(), 0, 2);
        parallelIo.writeArray(i->second->asString(), i->first);
        break;
      }

      default: {
        mTerm(1, AT_, "Unknown property type");
      }
    }
  }
}
