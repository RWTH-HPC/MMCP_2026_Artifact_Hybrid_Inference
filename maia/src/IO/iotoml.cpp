// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "iotoml.h"

#include <cstring>
#include <list>
#include "COMM/mpioverride.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "property.h"
#include "tomlutils.h"

using namespace std;
using namespace maia::io::toml;

MInt IOToml::solverCount() { return m_noSolvers; }


void IOToml::makeProperty(MProperty* p, const maia::io::toml::Property& prop) {
  //    TRACE();

  switch(prop.type()) {
    case MINT: {
      DEBUG("IOToml::makeProperty: found integer property :", MAIA_DEBUG_IO);
      p->propertyType = MINT;
      p->elements = prop.size();
      p->intField = new MInt[prop.size()];
      copy(prop.asInt().begin(), prop.asInt().end(), p->intField);
      break;
    }

    case MFLOAT: {
      DEBUG("IOToml::makeProperty: found float property :", MAIA_DEBUG_IO);
      p->propertyType = MFLOAT;
      p->elements = prop.size();
      p->floatField = new MFloat[prop.size()];
      copy(prop.asFloat().begin(), prop.asFloat().end(), p->floatField);
      break;
    }

    case MBOOL: {
      DEBUG("IOToml::makeProperty: found bool property :", MAIA_DEBUG_IO);
      p->propertyType = MBOOL;
      p->elements = prop.size();
      p->boolField = new MBool[prop.size()];
      copy(prop.asBool().begin(), prop.asBool().end(), p->boolField);
      break;
    }

    case MSTRING: {
      DEBUG("IOToml::makeProperty: found char property :", MAIA_DEBUG_IO);
      p->propertyType = MSTRING;
      p->elements = prop.size();
      p->stringField = new MString[prop.size()];
      copy(prop.asString().begin(), prop.asString().end(), p->stringField);
      break;
    }

    default: {
      stringstream errorMessage;
      errorMessage << "IOToml::makeProperty Error: Unsupported variable input type encountered!" << endl;
      TERMM(1, errorMessage.str());
    }
  }

  // Insert property into property map
  const pair<const MString, MProperty*> mp(p->name, p);

  // produce the lowercase variant of the name to warn if property exists with different case
  MString nameL = p->name;
  std::transform(nameL.begin(), nameL.end(), nameL.begin(), [](unsigned char c) { return std::tolower(c); });

  const pair<const MString, MProperty*> mpL(nameL, p);

  // check if the lowercase variant of this property already exists but no regular-case variant
  if(!(m_propertyMapLowercase->find(nameL) == m_propertyMapLowercase->end())
     && (m_propertyMap->find(p->name) == m_propertyMap->end())) {
    mTerm(1, AT_,
          "There are multiple occurrences of the property " + p->name
              + " with different cases. This should not happen!");
  }
  m_propertyMap->insert(mp);
  m_propertyMapLowercase->insert(mpL);

  DEBUG("IOToml::makeProperty: created default property ", MAIA_DEBUG_IO);
  DEBUG("IOToml::makeProperty: elements = " << p->elements, MAIA_DEBUG_IO);
  DEBUG("IOToml::makeProperty: m_noProperties = " << m_propertyMap->size(), MAIA_DEBUG_IO);
}


/** This function creates a default zone for all solvers, that do not appear
 * in a zone definition.
 *
 */
void IOToml::buildDefaultZone() {
  TRACE();
  /* create a zone named "default" */
  MZone* zone;
  zone = new MZone;
  zone->name.append("default");

  /* search for all solvers that are defined */
  list<MInt> solverList;
  for(zoneMap::const_iterator it = m_zoneMap->begin(); it != m_zoneMap->end(); it++) {
    for(MInt i = 0; i < it->second->noSolvers; i++) {
      solverList.push_back(it->second->solvers[i]);
    }
  }
  solverList.sort();

  /* The number of undefined solvers equals the number of total
     solvers minus the number of defined solvers ... */
  zone->noSolvers = (m_noSolvers - solverList.size());
  zone->solvers = new MInt[zone->noSolvers];

  /* this adds all solverId's that aren't defined to the default zone*/
  MInt i = 0;
  MInt simpleIterator = 0;
  list<MInt>::const_iterator solverIt = solverList.begin();

  while(i < m_noSolvers) {
    if((!solverList.empty()) && (*solverIt == i)) {
      i++;
      solverIt++;
    } else {
      zone->solvers[simpleIterator] = i;
      DEBUG("IOToml::buildDefaultZone: added solver " << i << " to the default zone. ", MAIA_DEBUG_IO);
      i++;
      simpleIterator++;
    }
  }

  /* insert the default zone in the zone map */
  m_zoneMap->insert(make_pair(zone->name, zone));
  DEBUG("IOToml::buildDefaultZone: default zone has " << zone->noSolvers << " solvers.", MAIA_DEBUG_IO);
}


assembly* IOToml::readPropertyFile(const MString& fileName) {
  TRACE();

  m_noSolvers = 0;
  m_propertyMap = new propertyMap;
  m_propertyMapLowercase = new propertyMap;
  m_zoneMap = new zoneMap;

  // Read in file in serial and distribute it to all domains
  if(globalDomainId() == 0) {
    // Read property file
    ifstream ifs(fileName);
    if(!ifs) {
      TERMM(1, "could not open property file '" + fileName + "'");
    }
    stringstream ss;
    ss << ifs.rdbuf();
    if(ss.str().empty()) {
      TERMM(1, "property file '" + fileName + "' appears to be empty");
    }

    // Determine length of null-terminated character string
    MInt length = ss.str().size() + 1;

    // Convert to char array
    vector<MChar> text(length);
    strcpy(text.data(), ss.str().c_str());

    // Broadcast string length and contents
    MPI_Bcast(&length, 1, MPI_INT, 0, globalMaiaCommWorld(), AT_, "length");
    MPI_Bcast(text.data(), length, MPI_CHAR, 0, globalMaiaCommWorld(), AT_, "text.data()");

    // Store file content in string
    m_rawText = ss.str();
  } else {
    // On all other domains, first receive string length
    MInt length = -1;
    MPI_Bcast(&length, 1, MPI_INT, 0, globalMaiaCommWorld(), AT_, "length");

    // Then receive file contents
    vector<MChar> text(length);
    MPI_Bcast(text.data(), length, MPI_CHAR, 0, globalMaiaCommWorld(), AT_, "text.data()");

    // Store file content in string
    m_rawText = text.data();
  }

  // Create parser for TOML files and parse file contents
  stringstream ss;
  ss << m_rawText;

  // collect solver aliases
  const std::shared_ptr<cpptoml::table> table = cpptoml::parser{ss}.parse();
  map<MString, MInt> solverAliases;
  collectSolverAliases(table, solverAliases);

  // Collect property information from TOML table
  vector<maia::io::toml::Property> properties;
  collectProperties(table, properties, solverAliases);

  // Get the number of solvers
  m_noSolvers = 1;
  for(auto&& prop : properties) {
    if(prop.name() == "noSolvers") {
      m_noSolvers = prop.asInt()[0];
    }
  }

  // Do not read in zones (deprecated), just build default zone
  buildDefaultZone();

  // Check zone consistency or die
  if(!checkZoneConsistency()) {
    stringstream errorMessage;
    errorMessage << "IOToml::readPropertyFile Error: Inconsistent solver List!\n Make sure "
                    "that solvers start with index 0 and are consecutive!"
                 << endl;
    TERMM(1, errorMessage.str());
  }

  // Create properties
  for(auto&& prop : properties) {
    /* if default property */
    if(!strstr(prop.name().c_str(), ".")) {
      // store default property in last solver,
      // the others to their belonging id's
      DEBUG("IOToml::readPropertyFile: default property : " << prop.name(), MAIA_DEBUG_IO);
      MProperty* p = new MProperty;
      p->solverId = m_noSolvers; // Insert default solver as last solver
      p->name.append(prop.name());
      makeProperty(p, prop); // create new Property
    }
    /* if a property defined for one or more solvers */
    else {
      const auto pname = prop.name(); // This is necessary to prevent use-after-free errors
      const char* du;
      du = strrchr(pname.c_str(), '.') + 1;
      const MInt singleSolverId = atoi(du);

      if(singleSolverId || *du == '0') {
        // Determine name without '.' and solver id
        const MString name = prop.name().substr(0, prop.name().find("."));

        DEBUG("IOToml::readPropertyFile: Found single solver property definition for solver " << singleSolverId,
              MAIA_DEBUG_IO);

        // Create new property with specific solver id
        MProperty* p = new MProperty;
        p->solverId = singleSolverId;
        p->name.append(name);
        makeProperty(p, prop); // create new Property
      }
    }
  }

  // Check consistency
  if(checkPropertyConsistency()) {
    DEBUG("IOToml::readPropertyFile: property consistency check succeeded", MAIA_DEBUG_IO);
  } else {
    stringstream errorMessage;
    errorMessage << "IOToml::readPropertyFile Error: some properties are not defined for "
                    "all solvers!"
                 << endl;
    TERMM(1, errorMessage.str());
  }

  m_assembly = new assembly;
  m_assembly->properties = m_propertyMap;
  m_assembly->propertiesLowercase = m_propertyMapLowercase;
  m_assembly->zones = m_zoneMap;

  return m_assembly;
}


MBool IOToml::checkPropertyConsistency() {
  TRACE();

  for(propertyIterator i = m_propertyMap->begin(); i != m_propertyMap->end(); i++) {
    /* if default property exists, then take next property */

    std::string propertyName = i->second->name;

    // name doesnot contain a solverId so is a default value
    std::size_t found = propertyName.find('.');
    if(found == std::string::npos) {
      continue;
    }

    const auto range = m_propertyMap->equal_range(propertyName);

    MInt size = 0;
    for(auto j = range.first; j != range.second; j++) {
      size++;
    }
    if(size == m_noSolvers) {
      continue;
    } else {
      std::cerr << "Did not find property " << propertyName << " for all solvers" << endl;
      return false;
    }
  }
  return true;
}

/** Checks the consistency of the zones, by creating a list, that contains
 *  all integer values from 0 to the number of the solvers, and comparing it
 *  to a list of all actual registered solverId's. For a consistent zoneList
 *  both lists must be equal.
 */
MBool IOToml::checkZoneConsistency() {
  TRACE();
  list<MInt> solverList;
  list<MInt> compareList;
  MInt index = 0;
  for(zoneMap::const_iterator it = m_zoneMap->begin(); it != m_zoneMap->end(); it++) {
    for(MInt i = 0; i < it->second->noSolvers; i++) {
      solverList.push_back(it->second->solvers[i]);
      compareList.push_back(index++);
    }
  }
  solverList.sort();
  if(solverList == compareList) {
    return true;
  } else {
    return false;
  }
}
