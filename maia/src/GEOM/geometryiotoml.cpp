// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometryiotoml.h"

#include <cstring>
#include "COMM/mpioverride.h"
#include "IO/tomlutils.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "UTIL/timer.h"

using namespace std;
using namespace maia::io::toml;


MInt GeometryIOToml::segmentCount() { return m_noSegments; }


void GeometryIOToml::makeProperty(GeometryProperty* p, const Property& prop) {
  //    TRACE();

  switch(prop.type()) {
    case MINT: {
      DEBUG("GeometryIOToml::makeProperty found integer property :", MAIA_DEBUG_USER1);
      p->propertyType = MINT;
      p->elements = prop.size();
      p->intField = new MInt[prop.size()];
      copy(prop.asInt().begin(), prop.asInt().end(), p->intField);
      break;
    }

    case MFLOAT: {
      DEBUG("GeometryIOToml::makeProperty found float property :", MAIA_DEBUG_USER1);
      p->propertyType = MFLOAT;
      p->elements = prop.size();
      p->floatField = new MFloat[prop.size()];
      copy(prop.asFloat().begin(), prop.asFloat().end(), p->floatField);
      break;
    }

    case MSTRING: {
      DEBUG("GeometryIOToml::makeProperty found char property :", MAIA_DEBUG_USER1);
      p->propertyType = MSTRING;
      p->elements = prop.size();
      p->stringField = new MString[prop.size()];
      copy(prop.asString().begin(), prop.asString().end(), p->stringField);
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "GeometryIOToml::makeProperty Error: Unsupported variable input type encountered!" << endl;
      TERMM(1, errorMessage.str());
    }
  }

  // Insert property into property map
  const pair<const MString, GeometryProperty*> mp(p->name, p);
  m_geometryPropertyMap->insert(mp);

  DEBUG("GeometryIOToml::makeProperty created default property ", MAIA_DEBUG_USER1);
  DEBUG("GeometryIOToml::makeProperty elements = " << p->elements, MAIA_DEBUG_USER1);
  DEBUG("GeometryIOToml::makeProperty m_noProperties = " << m_geometryPropertyMap->size(), MAIA_DEBUG_USER1);
}


void GeometryIOToml::readBodiesNewIOMethod(const std::vector<Property>& properties) {
  TRACE();
  m_log << "      * reading body information" << endl;

  Property bodySegments, bodySegmentsNames, bodySegmentsNum;
  for(auto&& prop : properties) {
    if(prop.name() == "body_segments") {
      bodySegments = prop;
    } else if(prop.name() == "body_segments_names") {
      bodySegmentsNames = prop;
    } else if(prop.name() == "body_segments_num") {
      bodySegmentsNum = prop;
    }
  }

  if(bodySegments.valid() && bodySegmentsNames.valid() && bodySegmentsNum.valid()) {
    MInt len_bsnum = bodySegmentsNum.size();
    MInt len_bs = bodySegments.size();

    MString bsallnames;
    MIntScratchSpace bs(len_bs, AT_, "bs");
    MStringScratchSpace bsname(len_bsnum, AT_, "bsname");
    MIntScratchSpace bsnum(len_bsnum, AT_, "bsnum");

    copy(bodySegments.asInt().begin(), bodySegments.asInt().end(), bs.data());
    bsallnames = bodySegmentsNames.asString()[0];
    copy(bodySegmentsNum.asInt().begin(), bodySegmentsNum.asInt().end(), bsnum.data());

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


void GeometryIOToml::readBodiesOldIOMethod(const std::vector<Property>& properties) {
  TRACE();

  m_log << "      * reading body information" << endl;

  // Iterate over body_segments
  for(auto&& prop : properties) {
    // Skip if wrong variable name
    if(!strstr(prop.name().c_str(), "body_segments.")) {
      continue;
    }

    const auto pname = prop.name(); // This is necessary to prevent use-after-free errors
    const char* du = strrchr(pname.c_str(), '.');

    /* create body element */
    Body* body;
    body = new Body;
    body->name.append(++du);
    DEBUG("GeometryIOToml::readBodies Body found (name: " << body->name << ")", MAIA_DEBUG_USER2);
    m_noBodies++;

    /* get segments for the body */
    MInt length = prop.size();
    DEBUG("GeometryIOToml::readBodies The Body has " << length << " solvers.", MAIA_DEBUG_USER2);
    body->noSegments = length;
    body->segments = new MInt[length];
    copy(prop.asInt().begin(), prop.asInt().end(), body->segments);

    /* insert the body into the bodyMap*/
    if(m_bodyMap->find(body->name) == m_bodyMap->end()) m_bodyMap->insert(make_pair(body->name, body));
  }

  buildDefaultBody();
}


//--------------------------------------------------------------------------
/** This function creates a default body for all segments, that do not appear
 * in a body definition.
 *
 */
void GeometryIOToml::buildDefaultBody() {
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
        DEBUG("GeometryIOToml::buildDefaultBody Added segment " << i << " to the default body. ", MAIA_DEBUG_USER1);
        i++;
        simpleIterator++;
      }
    }
  }

  /* insert the default body in the body map */
  m_bodyMap->insert(make_pair(body->name, body));
  DEBUG("GeometryIOToml::buildDefaultBody  Default body has " << body->noSegments << " segments.", MAIA_DEBUG_USER1);
  m_noBodies++;
}


//--------------------------------------------------------------------------
/** \brief check if the geometry property file is of new or old type and calls the according
 *function
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
geometryAssembly* GeometryIOToml::readPropertyFile(MString fileName) {
  TRACE();

  // init member variables
  m_newIOMethod = false;
  m_noSegments = 0;
  m_noBodies = 0;
  m_geometryPropertyMap = new geometryPropertyMap;
  m_bodyMap = new bodyMap;
  m_geometryAssembly = new geometryAssembly;

  if(domainId() == 0) {
    m_log << "    - rank 0 reads data from disk" << endl;

    // Read property file
    ifstream ifs(fileName);
    if(!ifs) {
      TERMM(1, "could not open geometry file '" + fileName + "'");
    }
    stringstream ss;
    ss << ifs.rdbuf();
    if(ss.str() == "") {
      TERMM(1, "geometry file '" + fileName + "' appears to be empty");
    }

    // Determine length of null-terminated character string
    MInt length = ss.str().size() + 1;

    // Convert to char array
    vector<MChar> text(length);
    strcpy(text.data(), ss.str().c_str());

    // Broadcast string length and contents
    MPI_Bcast(&length, 1, MPI_INT, 0, mpiComm(), AT_, "length");
    MPI_Bcast(text.data(), length, MPI_CHAR, 0, mpiComm(), AT_, "text.data()");

    // Store file content in string
    m_rawText = ss.str();
  } else {
    // On all other domains, first receive string length
    MInt length = -1;
    MPI_Bcast(&length, 1, MPI_INT, 0, mpiComm(), AT_, "length");

    // Then receive file contents
    vector<MChar> text(length);
    MPI_Bcast(text.data(), length, MPI_CHAR, 0, mpiComm(), AT_, "text.data()");

    // Store file content in string
    m_rawText = text.data();
  }

  // Create parser for TOML files and parse file contents
  stringstream ss;
  ss << m_rawText;
  vector<maia::io::toml::Property> properties;
  collectProperties(cpptoml::parser{ss}.parse(), properties);

  // Read information on newIOMethod and number of segments
  m_newIOMethod = false;
  m_noSegments = -1;
  for(auto&& prop : properties) {
    if(prop.name() == "newIOMethod") {
      m_newIOMethod = true;
      continue;
    }

    if(prop.name() == "noSegments") {
      m_noSegments = prop.asInt()[0];
      continue;
    }
  }

  // do the actual processing of the geometry properties
  NEW_SUB_TIMER(t_readRest, "rest geometry property", g_t_readGeomFile);
  NEW_SUB_TIMER(t_readBodies, "geometry bodies", g_t_readGeomFile);

  if(m_newIOMethod) {
    RECORD_TIMER_START(t_readBodies);
    readBodiesNewIOMethod(properties);
    RECORD_TIMER_STOP(t_readBodies);

    RECORD_TIMER_START(t_readRest);
    readPropertyFileNewIOMethod(properties);
    RECORD_TIMER_STOP(t_readRest);
  } else {
    RECORD_TIMER_START(t_readBodies);
    readBodiesOldIOMethod(properties);
    RECORD_TIMER_STOP(t_readBodies);

    RECORD_TIMER_START(t_readRest);
    readPropertyFileOldIOMethod(properties);
    RECORD_TIMER_STOP(t_readRest);
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
 **/
void GeometryIOToml::readPropertyFileOldIOMethod(const std::vector<Property>& properties) {
  TRACE();

  GeometryProperty* p;
  m_log << "      * reading rest of the properties" << endl;

  // Check body consistency or die
  if(!checkBodyConsistency()) {
    TERMM(1, "Body consistency check failed");
  }

  for(auto&& prop : properties) {
    // Store name for convenience
    const MString varName = prop.name();

    // if default property
    if(!strstr(varName.c_str(), ".")) {
      // if there's no dot in string (default prop)
      // store default property in last solver,
      // the others to their belonging id's

      DEBUG("GeometryIOToml::readPropertyFile  default property : " << varName, MAIA_DEBUG_USER1);
      p = new GeometryProperty;
      p->segmentId = m_noSegments; // Insert default segment as last segment
      p->name.append(varName);
      makeProperty(p, prop); // create new Property
    }
    // if a property defined for one or more bodies
    else {
      if(strstr(varName.c_str(), "_bodies.")) {
        DEBUG("GeometryIOToml::readPropertyFile normal property: " << varName, MAIA_DEBUG_USER1);
        MInt noBodies = prop.size();
        MString* bodies = new MString[noBodies];
        copy(prop.asString().begin(), prop.asString().end(), bodies);

        // find all segments for the property
        list<MInt> segmentList;
        for(MInt i = 0; i != noBodies; i++) {
          DEBUG("GeometryIOToml::readPropertyFile definition for body " << bodies[i], MAIA_DEBUG_USER1);
          bodyIterator zI;
          // look for the body in the bodyMap
          zI = m_bodyMap->find(bodies[i]);
          // append all segments of the body to the segmentlist
          for(MInt j = 0; j < zI->second->noSegments; j++)
            segmentList.push_back(zI->second->segments[j]);
        }

        char* du;
        du = strrchr(const_cast<MChar*>(varName.c_str()), '.');
        MString dummy(varName);
        // strip the varname of "_bodies.1"
        dummy.replace(dummy.find("_bodies."), dummy.size(), du);
        DEBUG("GeometryIOToml::readPropertyFile  found property : " << dummy, MAIA_DEBUG_USER1);
        // find the property id

        list<MInt>::const_iterator it = segmentList.begin();
        dummy.erase(dummy.find(".")); // erase the dot and the id
        // create the property for every segment
        for(; it != segmentList.end(); it++) {
          p = new GeometryProperty;
          p->segmentId = *it;
          p->name.append(dummy);
          makeProperty(p, prop); // create new Property
          DEBUG("GeometryIOToml::readPropertyFile created property for solver " << *it, MAIA_DEBUG_USER1);
        }
      } // end of if ( strstr ( varName,"_bodies." ) )
      else {
        // Determine the segmentId
        char* du;
        du = strrchr(const_cast<MChar*>(varName.c_str()), '.') + 1;
        MInt singleSegmentId = atoi(du);

        MString dummyName(varName);

        if(singleSegmentId || *du == '0') {
          DEBUG("Found single segment property definition for segment " << singleSegmentId, MAIA_DEBUG_IO);
          p = new GeometryProperty;
          p->segmentId = singleSegmentId;

          MString dummy = dummyName;

          dummyName.erase(dummyName.find(".")); // erase the dot and the id
          p->name.append(dummyName);
          makeProperty(p, prop); // create new Property
        }
      }
    }
  }
  // Insert here the check for the geometry property consistency
  if(checkGeometryPropertyConsistency()) {
    DEBUG("GeometryIOToml::readPropertyFile ** Property check successful \n", MAIA_DEBUG_USER1);
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
void GeometryIOToml::readPropertyFileNewIOMethod(const std::vector<Property>& properties) {
  TRACE();

  m_log << "      * reading rest of the properties" << endl;

  // Check body consistency or die
  if(!checkBodyConsistency()) {
    TERMM(1, "Body consistency check failed");
  }

  for(auto&& prop : properties) {
    // Store name for convenience
    const MString varName = prop.name();

    if(!strstr(varName.c_str(), ".")) {
      GeometryProperty* p = new GeometryProperty;
      p->segmentId = m_noSegments;
      p->name.append(varName);
      makeProperty(p, prop);
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
        MInt len = prop.size();
        MIntScratchSpace bcs(len, AT_, "bcs");
        copy(prop.asInt().begin(), prop.asInt().end(), bcs.data());

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

        MString allnames = prop.asString()[0];

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


MBool GeometryIOToml::checkGeometryPropertyConsistency() {
  TRACE();

  for(geometryPropertyIterator i = m_geometryPropertyMap->begin(); i != m_geometryPropertyMap->end(); i++) {
    /* if default property exists, then take next property */

    if(m_geometryPropertyMap->lower_bound(i->second->name)->second->segmentId == m_noSegments) {
      DEBUG("GeometryIOToml::checkPropertyConsistency default property exists for :" << i->second->name,
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
MBool GeometryIOToml::checkBodyConsistency() {
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
