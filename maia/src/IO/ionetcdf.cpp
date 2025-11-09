// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ionetcdf.h"
#include "COMM/mpioverride.h"
#include "parallelio.h"

#include <cstring>
#include "globals.h"

using namespace std;

//--------------------------------------------------------------------------
IONetcdf::IONetcdf() {
  DEBUG("IONetcdf::IONetcdf: entry", MAIA_DEBUG_ALLOCATION);
  DEBUG("IONetcdf::IONetcdf: return", MAIA_DEBUG_ALLOCATION);
}

//--------------------------------------------------------------------------
IONetcdf::~IONetcdf() {
  DEBUG("IONetcdf::~IONetcdf: entry", MAIA_DEBUG_ALLOCATION);
  DEBUG("IONetcdf::~IONetcdf: return", MAIA_DEBUG_ALLOCATION);
}

//--------------------------------------------------------------------------
MInt IONetcdf::solverCount() { return m_noSolvers; }

void IONetcdf::makeProperty(MProperty* p, const MString& name, ParallelIo* bdFile) {
  //    TRACE();
  using namespace maia::parallel_io;
  ParallelIo::size_type length, start;
  maiabd_type type;
  MInt noDims;
  type = bdFile->getDatasetType(name);
  noDims = bdFile->getDatasetNoDims(name);
#ifndef MAIA_WINDOWS
  switch(type) {
    case PIO_INT: {
      DEBUG("IONetcdf::makeProperty: found integer property :", MAIA_DEBUG_IO);
      p->propertyType = MINT;
      if(noDims == 0) {
        p->elements = 1;
        p->intField = new MInt[1];
        bdFile->readScalar(p->intField, name);
        if(p->intField[0] == NC_FILL_INT) {
          stringstream errorMessage;
          errorMessage << "IONetcdf::makeProperty: Error" << endl
                       << "Property int " << p->name << " is filled by value " << NC_FILL_INT
                       << ". This happens, because the property was declared in the property file but was never "
                          "initialized there. Go and fix it!";
          mTerm(1, AT_, errorMessage.str());
        }
        DEBUG("IONetcdf::makeProperty: " << p->intField[0], MAIA_DEBUG_IO);
      }
      if(noDims == 1) {
        /* look up the length of the array */
        length = bdFile->getArraySize(name);
        p->elements = length;
        p->intField = new MInt[length];
        bdFile->setOffset(length, 0);
        bdFile->readArray(p->intField, name);
        for(MInt i = 0; i < (MInt)length; i++) {
          if(p->intField[i] == NC_FILL_INT) {
            stringstream errorMessage;
            errorMessage << "Property int " << p->name << " array is filled by the fill value " << NC_FILL_INT
                         << " at the position " << i
                         << ". This happens, because the property was declared in the property file but never "
                            "initialized there. Go and fix it!";
            mTerm(1, AT_, errorMessage.str());
          }
        }
        // for (MInt i=0; i < length; i++)     // delete this line!
        //  DEBUG("IONetcdf::makeProperty: " << p->intField[i], MAIA_DEBUG_IO);
      }
      if(noDims == 2) {
        /*NOT YET IMPLEMENTED !*/
      }
      break;
    }
    case PIO_FLOAT: {
      DEBUG("IONetcdf::makeProperty: found float property :", MAIA_DEBUG_IO);
      p->propertyType = MFLOAT;
      if(noDims == 0) {
        p->elements = 1;
        p->floatField = new MFloat[1];
        bdFile->readScalar(p->floatField, name);
        DEBUG("IONetcdf::makeProperty: " << p->floatField[0], MAIA_DEBUG_IO);
        if(approx(p->floatField[0], NC_FILL_DOUBLE, MFloatEps)) {
          stringstream errorMessage;
          errorMessage << "Property double " << p->name << " is filled by the fill value " << NC_FILL_DOUBLE
                       << ". This happens, because the property was declared in the property file but never "
                          "initialized there. Go and fix it!";
          mTerm(1, AT_, errorMessage.str());
        }
      }
      if(noDims == 1) {
        /* look up the length of the array */
        length = bdFile->getArraySize(name);
        p->elements = length;
        p->floatField = new MFloat[length];
        bdFile->setOffset(length, 0);
        bdFile->readArray(p->floatField, name);
        for(MInt i = 0; i < (MInt)length; i++) {
          if(approx(p->floatField[i], NC_FILL_DOUBLE, MFloatEps)) {
            stringstream errorMessage;
            errorMessage << "Property double " << p->name << " array is filled by the fill value " << NC_FILL_DOUBLE
                         << " at the position " << i
                         << ". This happens, because the property was declared in the property file but never "
                            "initialized there. Go and fix it!";
            mTerm(1, AT_, errorMessage.str());
          }
        }
        // for (MInt i=0; i < length; i++)   //delete this line!
        // DEBUG("IONetcdf::makeProperty: " << p->floatField[i], MAIA_DEBUG_IO);
      }
      if(noDims == 2) {
        /*NOT YET IMPLEMENTED !*/
      }
      break;
    }
    case PIO_STRING: {
      DEBUG("IONetcdf::makeProperty: found char property :", MAIA_DEBUG_IO);
      p->propertyType = MSTRING;
      if(noDims == 0) {
        p->elements = 1;
        p->stringField = new MString[1];
        MString buf;
        bdFile->readScalar(&buf, name);
        if(buf[0] == NC_FILL_CHAR) {
          stringstream errorMessage;
          errorMessage << "Property char " << p->name << " is filled by the fill value \"" << NC_FILL_CHAR
                       << "\". This happens, because the property was declared in the property file but never "
                          "initialized there. Go and fix it!";
          mTerm(1, AT_, errorMessage.str());
        }
        p->stringField->append(buf);
      }
      if(noDims == 1) {
        p->elements = 1;
        p->stringField = new MString[1];
        /* look up the length of the string */
        length = bdFile->getArraySize(name);
        MString buf;
        bdFile->setOffset(length, 0);
        bdFile->readArray(&buf, name);
        if(buf[0] == NC_FILL_CHAR) {
          stringstream errorMessage;
          errorMessage << "Property array of char " << p->name << " is filled by the fill value \"" << NC_FILL_CHAR
                       << "\". This happens, because the property was declared in the property file but never "
                          "initialized there. Go and fix it!";
          mTerm(1, AT_, errorMessage.str());
        }
        p->stringField->append(buf);
        DEBUG("IONetcdf::makeProperty: " << p->stringField[0], MAIA_DEBUG_IO);
      }
      if(noDims == 2) {
        /* look up the number of strings */
        length = bdFile->getArraySize(name, 0);
        p->elements = length;
        p->stringField = new MString[length];
        /* loop over number of strings */
        for(MInt i = 0; i < p->elements; i++) {
          start = i;
          bdFile->setOffset(1, start, 2);
          MString buf;
          bdFile->readArray(&buf, name);
          if(buf[0] == NC_FILL_CHAR) {
            stringstream errorMessage;
            errorMessage << "Property array of array of char " << p->name << " element " << i << " of " << p->elements
                         << " is filled by the fill value \"" << NC_FILL_CHAR
                         << "\". This happens, because the property was declared in the property file but never "
                            "initialized there. Go and fix it!";
            mTerm(1, AT_, errorMessage.str());
          }
          p->stringField[i] = buf;
          DEBUG("IONetcdf::makeProperty: " << p->stringField[i], MAIA_DEBUG_IO);
        }
      }
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "IONetcdf::makeProperty Error: Unsupported variable input type encountered!" << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }
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

#else
  // not supported
  TERM(-1);
#endif

  DEBUG("IONetcdf::makeProperty: created default property ", MAIA_DEBUG_IO);
  DEBUG("IONetcdf::makeProperty: elements = " << p->elements, MAIA_DEBUG_IO);
  DEBUG("IONetcdf::makeProperty: m_noProperties = " << m_propertyMap->size(), MAIA_DEBUG_IO);
}

//--------------------------------------------------------------------------
/**
 *
 *
 */
void IONetcdf::readZones(ParallelIo* bdFile) {
  TRACE();
  MInt noVariables = 0;
  const char MPropertySeperator = '.'; // Seperator that defines different Properties

  /* read and store the zone information and the solver info*/
  const char* du; // tmp variable

  MString varName;
  ParallelIo::size_type length;

  vector<MString> varNames = bdFile->getDatasetNames();
  noVariables = varNames.size();
  for(MInt n = 0; n < noVariables; n++) {
    varName = varNames[n];
    if(strstr(varName.c_str(), "zone_solvers.")) {          // if a zone has been identified
      du = strrchr(varName.c_str(), MPropertySeperator);    // get pointer to where the "." is

      /* create zone element */
      MZone* zone;
      zone = new MZone;
      zone->name.append(++du);
      if(zone->name == "default") {
        mTerm(1, AT_, "IONetcdf::readZones Error: The zone name \"default\" is reserved and cannot be used!");
      }
      DEBUG("IONetcdf::readZones: zone found (name: " << zone->name << ")", MAIA_DEBUG_USER2);

      /* get solvers for the zone */
      length = bdFile->getArraySize(varName);
      DEBUG("IONetcdf::readZones: the Zone has " << length << " solvers.", MAIA_DEBUG_USER2);
      zone->noSolvers = length;
      zone->solvers = new MInt[length];
      /* There is an implicit conversion from int to unsigned int between netcdf
         and the zone solvers */
      for(ParallelIo::size_type i = 0; i < length; i++) {
        MInt dsolvers;
        bdFile->setOffset(1, i);
        bdFile->readArray(&dsolvers, varName);
        zone->solvers[i] = dsolvers;
      }
      /* insert the zone into the zoneMap*/
      if(m_zoneMap->find(zone->name) == m_zoneMap->end()) {
        m_zoneMap->insert(make_pair(zone->name, zone));
      } else {
        // if a zone with the same name exists throw error
        stringstream errorMessage;
        errorMessage << " Error : Double occurence of zone id! " << endl;
        mTerm(1, AT_, errorMessage.str());
      }
    }
  }
  buildDefaultZone();
}

//--------------------------------------------------------------------------
/** This function creates a default zone for all solvers, that do not appear
 * in a zone definition.
 *
 */
void IONetcdf::buildDefaultZone() {
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
      DEBUG("IONetcdf::buildDefaultZone: added solver " << i << " to the default zone. ", MAIA_DEBUG_IO);
      i++;
      simpleIterator++;
    }
  }

  /* insert the default zone in the zone map */
  m_zoneMap->insert(make_pair(zone->name, zone));
  DEBUG("IONetcdf::buildDefaultZone: default zone has " << zone->noSolvers << " solvers.", MAIA_DEBUG_IO);
}

//--------------------------------------------------------------------------
/**
 *
 *
 */
assembly* IONetcdf::readPropertyFile(const MString& name) {
  TRACE();
  MBool readNew = true;
  MInt noVariables = 0;
  m_noSolvers = 0;
  m_propertyMap = new propertyMap;
  m_propertyMapLowercase = new propertyMap;
  m_zoneMap = new zoneMap;

  if(readNew) {
    // split into serial read and distribution is done here
    if(globalDomainId() == 0) { // reading  is only done in serial
      //??? todo labels:IO change MPropertySeperator to MPropertySeparator
      const char MPropertySeperator = '.'; // Separator that defines different Properties

      /*open the property file*/
      ParallelIo parallelIo(name, maia::parallel_io::PIO_READ, MPI_COMM_SELF);

      /* get the number of solvers  */
      MInt dsolvers;
      if(!parallelIo.hasDataset("noSolvers", 0)) {
        dsolvers = 1;
        // cerr << "IONetcdf::readPropertyFile, property \"noSolvers\" is not specified in properties file \"" <<
        // name << "\" . Thats why it is set to 1." << endl;
      } else {
        parallelIo.readScalar(&dsolvers, "noSolvers");
      }
      m_noSolvers = dsolvers;
      /* read and store the zone information and the solver info*/

      readZones(&parallelIo);

      MString varName;

      MProperty* p;
      /* check the consistency of the zone*/
      if(checkZoneConsistency()) {
        vector<MString> varNames = parallelIo.getDatasetNames();
        noVariables = varNames.size();

        for(MInt id = 0; id < noVariables; id++) {
          varName = varNames[id];
          /* if default property */
          if(!strstr(varName.c_str(), ".")) {
            // if there's no dot in string (default prop)
            // store default property in last solver,
            // the others to their belonging id's

            DEBUG("IONetcdf::readNCPropertyFile: default property : " << varName, MAIA_DEBUG_IO);
            p = new MProperty;
            p->solverId = m_noSolvers; // Insert default solver as last solver
            p->name.append(varName);
            // skip since netcdf 2 -> netcdf 4 conversion messes this up
            if("standardTextLength" == varName) {
              continue;
            }
            makeProperty(p, varName, &parallelIo); // create new Property
          }
          /* if a property defined for one or more zones */
          else {
            if(strstr(varName.c_str(), "_zones.")) {
              DEBUG("IONetcdf::readPropertyFile: regular property: " << varName, MAIA_DEBUG_IO);
              ParallelIo::size_type length;
              MInt noDims;
              noDims = parallelIo.getDatasetNoDims(varName);
              MString* zones = nullptr;
              MInt noZones = 1;
              DEBUG("IONetcdf::readPropertyFile: no. of dimensions = " << noDims, MAIA_DEBUG_IO);
              switch(noDims) {
                case 0: {
                  /* if only one char (and one zone) */
                  MString buf;
                  parallelIo.readScalar(&buf, varName);
                  zones = new MString(buf);
                  break;
                }
                case 1: {
                  /* look up the length of the string */
                  length = parallelIo.getArraySize(varName);
                  MString buf;
                  parallelIo.setOffset(length, 0);
                  parallelIo.readArray(&buf, varName);
                  zones = new MString(buf);
                  break;
                }
                case 2: {
                  /* look up the number of zones */
                  ParallelIo::size_type dzones;
                  dzones = parallelIo.getArraySize(varName, 0);
                  noZones = (MInt)dzones;
                  zones = new MString[noZones];
                  /* look up the length of the string */
                  length = parallelIo.getArraySize(varName, 1);
                  ParallelIo::size_type start;
                  /* loop over number of zones */
                  for(MInt i = 0; i < noZones; i++) {
                    start = i;
                    parallelIo.setOffset(1, start, 2);
                    MString buf;
                    parallelIo.readArray(&buf, varName);
                    (zones[i]).append(buf);
                    DEBUG("IONetcdf::readPropertyFile: " << zones[i], MAIA_DEBUG_IO);
                  }
                  break;
                }
                default: {
                  stringstream errorMessage;
                  errorMessage << "IONetcdf::readPropertyFile Error: only one dimensional zone lists allowed!" << endl;
                  mTerm(1, AT_, errorMessage.str());
                }
              }
              /* find all solvers for the property */

              list<MInt> solverList;
              for(MInt i = 0; i != noZones; i++) {
                DEBUG("IONetcdf::readPropertyFile: definition for zone " << zones[i], MAIA_DEBUG_IO);
                zoneIterator zI;
                /* look for the zone in the zoneMap */
                zI = m_zoneMap->find(zones[i]);
                /* append all solvers of the zone to the solverlist */
                for(MInt j = 0; j < zI->second->noSolvers; j++) {
                  solverList.push_back(zI->second->solvers[j]);
                }
              }


              const char* du;
              du = strrchr(varName.c_str(), MPropertySeperator);
              MString dummy(varName);
              /* strip the varname of "_zones.1" */
              dummy.replace(dummy.find("_zones."), dummy.size(), du);
              DEBUG("IONetcdf::readPropertyFile: found property: " << dummy, MAIA_DEBUG_IO);
              MString dName(dummy);

              list<MInt>::const_iterator it = solverList.begin();
              dummy.erase(dummy.find(".")); // erase the dot and the id
              /* create the property for every solver */
              for(; it != solverList.end(); it++) {
                p = new MProperty;
                p->solverId = *it;
                p->name.append(dummy);
                makeProperty(p, dName, &parallelIo); // create new Property
                DEBUG("IONetcdf::readPropertyFile: created property for solver " << *it, MAIA_DEBUG_IO);
              }
            } // end of if ( strstr ( varName,"_zones." ) )
            else {
              // Determine the solverId
              const char* du;
              du = strrchr(varName.c_str(), MPropertySeperator) + 1;
              MInt singleSolverId = atoi(du);

              MString dummyName(varName);

              if(singleSolverId || *du == '0') {
                DEBUG("IONetcdf::readPropertyFile: Found single solver property definition for solver "
                          << singleSolverId,
                      MAIA_DEBUG_IO);
                p = new MProperty;
                p->solverId = singleSolverId;

                MString dName(dummyName);

                dummyName.erase(dummyName.find(".")); // erase the dot and the id
                p->name.append(dummyName);
                if("standardTextLength" != varName) {
                  makeProperty(p, dName, &parallelIo); // create new Property
                }
              }
            }
          }
        }
        /* Insert here the check for the property consistency */
        if(checkPropertyConsistency()) {
          DEBUG("IONetcdf::readPropertyFile: property consistency check succeeded", MAIA_DEBUG_IO);
        } else {
          stringstream errorMessage;
          errorMessage << "IONetcdf::readPropertyFile Error: some properties are not defined for all solvers!" << endl;
          mTerm(1, AT_, errorMessage.str());
        }
      } else {
        stringstream errorMessage;
        errorMessage << "IONetcdf::readPropertyFile Error: Inconsistent solver List!\n Make sure that solvers start "
                        "with index 0 and are consecutive!"
                     << endl;
        mTerm(1, AT_, errorMessage.str());
      }
      //-------------NEW: SEND THE PROPERTIES AND ZONES
      //-> PROPERTIES
      //-->gather all the data
      MInt noProperties = m_propertyMap->size();
      //-->get the number of counts of each property and add it up
      MInt totalCount = 0;
      for(propertyIterator prop = m_propertyMap->begin(); prop != m_propertyMap->end(); prop++) {
        totalCount += prop->second->count();
      }
      MInt sizeofallProp = noProperties * (256 + 4 * sizeof(MInt)) + totalCount * (256 + sizeof(MInt));
      //-->allocate the size of the send array [name, type, count, values(asume to be 256chars)] (converted to chars)
      MCharScratchSpace propMapArray(sizeofallProp, AT_, "propertyMapAsChar");
      propMapArray.fill(0);
      //-->put the variables into the scratchspace
      MInt pCounter = 0;
      for(propertyIterator prop = m_propertyMap->begin(); prop != m_propertyMap->end(); prop++) {
        MInt asize = prop->first.size();
        memcpy(&propMapArray(pCounter), &asize, sizeof(MInt));
        pCounter += sizeof(MInt);
        memcpy(&propMapArray(pCounter), prop->first.c_str(), asize * sizeof(MChar)); // copy the name
        pCounter += asize * sizeof(MChar);
        memcpy(&propMapArray(pCounter), (void*)&(prop->second->solverId), sizeof(MInt));
        pCounter += sizeof(MInt);
        memcpy(&propMapArray(pCounter), (void*)&(prop->second->propertyType), sizeof(MInt));
        pCounter += sizeof(MInt);
        MInt count = prop->second->count();
        memcpy(&propMapArray(pCounter), &count, sizeof(MInt));
        pCounter += sizeof(MInt);
        switch(prop->second->propertyType) {
          case MINT: {
            memcpy((void*)&propMapArray(pCounter), (void*)&(prop->second->intField[0]), count * sizeof(MInt));
            pCounter += count * sizeof(MInt);
            break;
          }
          case MFLOAT: {
            memcpy((void*)&propMapArray(pCounter), (void*)&(prop->second->floatField[0]), count * sizeof(MFloat));
            pCounter += count * sizeof(MFloat);
            break;
          }
          case MSTRING: {
            for(MInt i = 0; i < count; i++) {
              asize = prop->second->stringField[i].size();
              memcpy(&propMapArray(pCounter), &asize, sizeof(MInt));
              pCounter += sizeof(MInt);
              memcpy(&propMapArray(pCounter), prop->second->stringField[i].c_str(),
                     prop->second->stringField[i].size() * sizeof(MChar));
              pCounter += asize * sizeof(MChar);
            }
            break;
          }
          default: {
            mTerm(1, AT_, "no such variable type!!!!");
            break;
          }
        }
      }
      //-> ZONES
      //-->gather all the data
      MInt noZones = m_zoneMap->size();
      //->get the number of counts of each zone and add it up
      MInt zoneCount = 0;
      for(zoneIterator zone = m_zoneMap->begin(); zone != m_zoneMap->end(); zone++) {
        zoneCount += zone->second->noSolvers;
      }
      MInt sizeofZones = noProperties * (2 * (256 + sizeof(MInt)) + sizeof(MInt) * 2) + zoneCount * sizeof(MInt);
      MCharScratchSpace zoneMapArray(sizeofZones, AT_, "zoneMapArrayAsChar");
      zoneMapArray.fill(0);

      pCounter = 0;
      for(zoneIterator zone = m_zoneMap->begin(); zone != m_zoneMap->end(); zone++) {
        MInt asize = zone->first.size();
        memcpy(&zoneMapArray(pCounter), &asize, sizeof(MInt));
        pCounter += sizeof(MInt);
        memcpy(&zoneMapArray(pCounter), zone->first.c_str(), asize * sizeof(MChar)); // copy the name
        pCounter += asize * sizeof(MChar);
        memcpy(&zoneMapArray(pCounter), &(zone->second->id), sizeof(MInt));
        pCounter += sizeof(MInt);
        memcpy(&zoneMapArray(pCounter), &(zone->second->noSolvers), sizeof(MInt));
        pCounter += sizeof(MInt);
        memcpy(&zoneMapArray(pCounter), &(zone->second->solvers[0]), sizeof(MInt) * zone->second->noSolvers);
        pCounter += sizeof(MInt) * zone->second->noSolvers;
      }
      // END OF DATA COLLECTION --> SENDING THE INFORMATION TO OTHER PROCESSES
      MInt dummyInt[5] = {noProperties, totalCount, zoneCount, noZones, m_noSolvers};
      MPI_Bcast(&dummyInt, 5, MPI_INT, 0, globalMaiaCommWorld(), AT_, "dummyInt");
      MPI_Bcast(propMapArray.getPointer(), sizeofallProp, MPI_CHAR, 0, globalMaiaCommWorld(), AT_,
                "propMapArray.getPointer()");
      MPI_Bcast(zoneMapArray.getPointer(), sizeofZones, MPI_CHAR, 0, globalMaiaCommWorld(), AT_, "zoneMapArray.getPointer()");
      //->root has finished
    } else { // ALL OTHER PROCESSES
      //--> receive data
      MInt dummyInt[5] = {0, 0, 0, 0, 0};
      MPI_Bcast(&dummyInt, 5, MPI_INT, 0, globalMaiaCommWorld(), AT_, "dummyInt");
      MInt noProperties = dummyInt[0];
      MInt totalCount = dummyInt[1];
      MInt zoneCount = dummyInt[2];
      MInt noZones = dummyInt[3];
      m_noSolvers = dummyInt[4];
      MInt sizeofallProp = noProperties * (256 + 4 * sizeof(MInt)) + totalCount * (256 + sizeof(MInt));
      MInt sizeofZones = noProperties * (2 * (256 + sizeof(MInt)) + sizeof(MInt) * 2) + zoneCount * sizeof(MInt);
      MCharScratchSpace propMapArray(sizeofallProp, AT_, "proptertydistribution");
      propMapArray.fill(0);
      MCharScratchSpace zoneMapArray(sizeofZones, AT_, "zoneMapArrayAsChar");
      zoneMapArray.fill(0);
      //-> receive data from root
      MPI_Bcast(propMapArray.getPointer(), sizeofallProp, MPI_CHAR, 0, globalMaiaCommWorld(), AT_,
                "propMapArray.getPointer()");
      MPI_Bcast(zoneMapArray.getPointer(), sizeofZones, MPI_CHAR, 0, globalMaiaCommWorld(), AT_, "zoneMapArray.getPointer()");
      //-> create properties again
      MProperty* p;
      MInt pCounter = 0;
      for(MInt i = 0; i < noProperties; i++) {
        // create the property again
        p = new MProperty;
        // analyze the field
        MInt asize = 0;
        //->get the size of the name
        memcpy(&asize, &propMapArray[pCounter], sizeof(MInt));
        pCounter += sizeof(MInt);
        char* a = new char[asize + 1];
        //->get the name
        memcpy(a, &propMapArray[pCounter], asize);
        a[asize] = '\0';
        pCounter += asize * sizeof(MChar);
        const char* b = a;
        MString s(b);
        delete[] a;
        p->name = s;
        memcpy(&(p->solverId), &propMapArray[pCounter], sizeof(MInt));
        pCounter += sizeof(MInt);
        MInt propertyType = -1;
        memcpy(&propertyType, &propMapArray[pCounter], sizeof(MInt));
        pCounter += sizeof(MInt);
        switch(propertyType) {
          case 0: {
            p->propertyType = MINT;
            break;
          }
          case 1: {
            p->propertyType = MFLOAT;
            break;
          }
          case 2: {
            p->propertyType = MSTRING;
            break;
          }
          default: {
            mTerm(1, AT_, "no such variable type");
            break;
          }
        }
        memcpy(&(p->elements), &propMapArray[pCounter], sizeof(MInt));
        pCounter += sizeof(MInt);
        switch(p->propertyType) {
          case MINT: {
            p->intField = new MInt[p->elements];
            memcpy(&(p->intField[0]), &propMapArray[pCounter], p->elements * sizeof(MInt));
            pCounter += sizeof(MInt) * p->elements;
            break;
          }
          case MFLOAT: {
            p->floatField = new MFloat[p->elements];
            memcpy(&(p->floatField[0]), &propMapArray[pCounter], p->elements * sizeof(MFloat));
            pCounter += sizeof(MFloat) * p->elements;
            break;
          }
          case MSTRING: {
            p->stringField = new MString[p->elements];
            for(MInt j = 0; j < p->elements; j++) {
              asize = 0;
              memcpy(&asize, &propMapArray(pCounter), sizeof(MINT));
              pCounter += sizeof(MInt);
              char* c = new char[asize + 1];
              memcpy(c, &propMapArray[pCounter], asize * sizeof(MChar));
              pCounter += asize * sizeof(MChar);
              c[asize] = '\0';
              const char* d = c;
              MString st(d);
              p->stringField[j] = st;
            }
            break;
          }
          default: {
            mTerm(1, AT_, "no such variable type");
            break;
          }
        }
        // put property to property map
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
      }
      //-> zone reconstruction
      MZone* z;
      MInt pCounterZ = 0;
      // zones
      for(MInt zo = 0; zo < noZones; zo++) {
        // zones
        MInt asize = 0;
        z = new MZone;
        memcpy(&asize, &zoneMapArray[pCounterZ], sizeof(MInt));
        pCounterZ += sizeof(MInt);
        char* c = new char[asize + 1];
        memcpy(c, &zoneMapArray[pCounterZ], asize * sizeof(MChar));
        pCounterZ += asize * sizeof(MChar);
        c[asize] = '\0';
        const char* d = c;
        MString st(d);
        z->name = st;
        memcpy(&(z->id), &zoneMapArray[pCounterZ], sizeof(MInt));
        pCounterZ += sizeof(MInt);
        memcpy(&(z->noSolvers), &zoneMapArray[pCounterZ], sizeof(MInt));
        pCounterZ += sizeof(MInt);
        z->solvers = new MInt[z->noSolvers];
        memcpy(&(z->solvers[0]), &zoneMapArray[pCounterZ], z->noSolvers * sizeof(MInt));
        pCounterZ += z->noSolvers * sizeof(MInt);
        m_zoneMap->insert(make_pair(z->name, z));
      }
    }
  } else { // use old approach already better and corrected
    //??? todo labels:IO change MPropertySeperator to MPropertySeparator
    const char MPropertySeperator = '.'; // Separator that defines different Properties
    /*open the property file*/
    ParallelIo parallelIo(name, maia::parallel_io::PIO_READ, globalMaiaCommWorld());
    /* get the number of solvers  */
    MInt dsolvers;
    if(!parallelIo.hasDataset("noSolvers", 0)) {
      dsolvers = 1;
    } else {
      parallelIo.readScalar(&dsolvers, "noSolvers");
    }
    m_noSolvers = dsolvers;
    /* read and store the zone information and the solver info*/
    readZones(&parallelIo);
    MString varName;
    MProperty* p;
    /* check the consistency of the zone*/
    if(checkZoneConsistency()) {
      vector<MString> varNames = parallelIo.getDatasetNames();
      noVariables = varNames.size();
      for(MInt id = 0; id < noVariables; id++) {
        varName = varNames[id];
        /* if default property */
        if(!strstr(varName.c_str(), ".")) {
          // if there's no dot in string (default prop)
          // store default property in last solver,
          // the others to their belonging id's
          DEBUG("IONetcdf::readNCPropertyFile: default property : " << varName, MAIA_DEBUG_IO);
          p = new MProperty;
          p->solverId = m_noSolvers; // Insert default solver as last solver
          p->name.append(varName);
          // skip since netcdf 2 -> netcdf 4 conversion messes this up
          if("standardTextLength" == varName) {
            continue;
          }
          makeProperty(p, varName, &parallelIo); // create new Property
        }
        /* if a property defined for one or more zones */
        else {
          if(strstr(varName.c_str(), "_zones.")) {
            DEBUG("IONetcdf::readPropertyFile: regular property: " << varName, MAIA_DEBUG_IO);
            ParallelIo::size_type length;
            MInt noDims;
            noDims = parallelIo.getDatasetNoDims(varName);
            MString* zones = nullptr;
            MInt noZones = 1;
            DEBUG("IONetcdf::readPropertyFile: no. of dimensions = " << noDims, MAIA_DEBUG_IO);
            switch(noDims) {
              case 0: {
                /* if only one char (and one zone) */
                MString buf;
                parallelIo.readScalar(&buf, varName);
                zones = new MString(buf);
                break;
              }
              case 1: {
                /* look up the length of the string */
                length = parallelIo.getArraySize(varName);
                MString buf;
                parallelIo.setOffset(length, 0);
                parallelIo.readArray(&buf, varName);
                zones = new MString(buf);
                break;
              }
              case 2: {
                /* look up the number of zones */
                ParallelIo::size_type dzones;
                dzones = parallelIo.getArraySize(varName, 0);
                noZones = (MInt)dzones;
                zones = new MString[noZones];
                /* look up the length of the string */
                length = parallelIo.getArraySize(varName, 1);
                ParallelIo::size_type start;
                /* loop over number of zones */
                for(MInt i = 0; i < noZones; i++) {
                  start = i;
                  parallelIo.setOffset(1, start, 2);
                  MString buf;
                  parallelIo.readArray(&buf, varName);
                  (zones[i]).append(buf);
                  DEBUG("IONetcdf::readPropertyFile: " << zones[i], MAIA_DEBUG_IO);
                }
                break;
              }
              default: {
                stringstream errorMessage;
                errorMessage << "IONetcdf::readPropertyFile Error: only one dimensional zone lists allowed!" << endl;
                mTerm(1, AT_, errorMessage.str());
              }
            }
            /* find all solvers for the property */

            list<MInt> solverList;
            for(MInt i = 0; i != noZones; i++) {
              DEBUG("IONetcdf::readPropertyFile: definition for zone " << zones[i], MAIA_DEBUG_IO);
              zoneIterator zI;
              /* look for the zone in the zoneMap */
              zI = m_zoneMap->find(zones[i]);
              /* append all solvers of the zone to the solverlist */
              for(MInt j = 0; j < zI->second->noSolvers; j++) {
                solverList.push_back(zI->second->solvers[j]);
              }
            }
            const char* du;
            du = strrchr(varName.c_str(), MPropertySeperator);
            MString dummy(varName);
            /* strip the varname of "_zones.1" */
            dummy.replace(dummy.find("_zones."), dummy.size(), du);
            DEBUG("IONetcdf::readPropertyFile: found property: " << dummy, MAIA_DEBUG_IO);
            MString dName(dummy);

            list<MInt>::const_iterator it = solverList.begin();
            dummy.erase(dummy.find(".")); // erase the dot and the id
            /* create the property for every solver */
            for(; it != solverList.end(); it++) {
              p = new MProperty;
              p->solverId = *it;
              p->name.append(dummy);
              makeProperty(p, dName, &parallelIo); // create new Property
              DEBUG("IONetcdf::readPropertyFile: created property for solver " << *it, MAIA_DEBUG_IO);
            }
          } // end of if ( strstr ( varName,"_zones." ) )
          else {
            // Determine the solverId
            const char* du;
            du = strrchr(varName.c_str(), MPropertySeperator) + 1;
            MInt singleSolverId = atoi(du);
            MString dummyName(varName);
            if(singleSolverId || *du == '0') {
              DEBUG("IONetcdf::readPropertyFile: Found single solver property definition for solver " << singleSolverId,
                    MAIA_DEBUG_IO);
              p = new MProperty;
              p->solverId = singleSolverId;
              MString dName(dummyName);
              dummyName.erase(dummyName.find(".")); // erase the dot and the id
              p->name.append(dummyName);
              if("standardTextLength" != varName) {
                makeProperty(p, dName, &parallelIo); // create new Property
              }
            }
          }
        }
      }
      /* Insert here the check for the property consistency */
      if(checkPropertyConsistency()) {
        DEBUG("IONetcdf::readPropertyFile: property consistency check succeeded", MAIA_DEBUG_IO);
      } else {
        stringstream errorMessage;
        errorMessage << "IONetcdf::readPropertyFile Error: some properties are not defined for all solvers!" << endl;
        mTerm(1, AT_, errorMessage.str());
      }
    } else {
      stringstream errorMessage;
      errorMessage << "IONetcdf::readPropertyFile Error: Inconsistent solver List!\n Make sure that solvers start "
                      "with index 0 and are consecutive!"
                   << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  } // end of else old approach for reading

  m_assembly = new assembly;
  m_assembly->properties = m_propertyMap;
  m_assembly->propertiesLowercase = m_propertyMapLowercase;
  m_assembly->zones = m_zoneMap;

  return m_assembly;
}

//--------------------------------------------------------------------------------
/**
 *
 *
 */
MBool IONetcdf::checkPropertyConsistency() {
  TRACE();

  for(propertyIterator i = m_propertyMap->begin(); i != m_propertyMap->end(); i++) {
    /* if default property exists, then take next property */

    if(m_propertyMap->lower_bound(i->second->name)->second->solverId == m_noSolvers) {
      DEBUG("IONetcdf::checkPropertyConsistency: default property exists for :" << i->second->name, MAIA_DEBUG_IO);
    } else {
      return false;
    }
  }

  return true;
}

//--------------------------------------------------------------------------
/** Checks the consistency of the zones, by creating a list, that contains
 *  all integer values from 0 to the number of the solvers, and comparing it
 *  to a list of all actual registered solverId's. For a consistent zoneList
 *  both lists must be equal.
 */
MBool IONetcdf::checkZoneConsistency() {
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
