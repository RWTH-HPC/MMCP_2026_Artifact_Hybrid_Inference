// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometrycontext.h"

#include "geometryionetcdf.h"
#include "geometryiotoml.h"
#include "globals.h"

using namespace std;

//--------------------------------------------------------------------------
/** Returns a copy of the bodyMap
 */
bodyMap GeometryContext::getBodies() { return *m_bodyMap; }


//--------------------------------------------------------------------------
/** Creates an object of the specified filetype to open the property file
 */
void GeometryContext::readPropertyFile(FileType fileType, const MChar* fileName) {
  TRACE();

  MBool success = false; // successfully accessed the property file ?
  m_name = fileName;
  while(!success) {
    switch(fileType) {
      case NETCDF:
        m_geometryIOBase = new GeometryIONetcdf(mpiComm());
        break;
      case TOML:
        m_geometryIOBase = new GeometryIOToml(mpiComm());
        break;
        /*		 case HDF:
                 m_geometryIOBase  = new GeometryIOHdf;
                 break;
             case ASCII:
                 m_geometryIOBase  = new GeometryIOAscii;
                 break; */
      default:
        break;
    }
    m_geometryAssembly = m_geometryIOBase->readPropertyFile(m_name);
    m_geometryPropertyMap = m_geometryAssembly->geometryProperties;
    m_bodyMap = m_geometryAssembly->bodies;
    m_noSegments = m_geometryIOBase->segmentCount();
    success = true;
  }
}

//--------------------------------------------------------------------------
/** Returns a pointer to the requested property.
 *
 *
 */
GeometryProperty* GeometryContext::getProperty(const MString& name, MInt segmentId) {
  //     TRACE();
  GeometryProperty* property = nullptr;
  if(m_geometryPropertyMap->find(name) == m_geometryPropertyMap->end()) {
    MString message;
    message.append("GeometryContext::getProperty Property: \"");
    message.append(name);
    message.append("\" not found!");
    // SX8 cannot compile with excepetion handling
    // throw ( MPropertyNotFound ( message ) );
    cerr << message;
    // allocate and return empty property
    property = new GeometryProperty();
    return property;
  }

  m_pair = m_geometryPropertyMap->equal_range(name); // lookup all properties
                                                     // with the requested key
  MBool propertyFound = false;
  for(geometryPropertyIterator i = m_pair.first; i != m_pair.second; i++) { // search for the segmentId
    DEBUG("GeometryContext::getProperty  segmentId:" << i->second->segmentId, MAIA_DEBUG_USER1);
    if(i->second->segmentId == segmentId) { // look for the requested segmentId
      property = i->second;                 // and deliver the value
      propertyFound = true;
    }
  }
  /* if no property found use standard value (remember, a standard value allways has
     the segmentId of the last segmentId +1, i.e. it is equal the no of segments
     expressed in m_noSegments). Since the consistency check guarantees the existence
     of a standard value, this should always work!*/
  if(!propertyFound) {
    for(geometryPropertyIterator i = m_pair.first; i != m_pair.second; i++) {
      if(i->second->segmentId == m_noSegments) { // look for a default value
        property = i->second;                    // and deliver the value
        DEBUG("GeometryContext::getProperty :" << name << " using standard property", MAIA_DEBUG_USER1);
      } else {
        DEBUG("GeometryContext::getProperty :" << name << " no standard value found !", MAIA_DEBUG_USER1);
      }
    }
  }
  return property;
}

//--------------------------------------------------------------------------
/** This method cares for the correct deletion of the property structs.
 * therefore it runs through the propertyMap and deletes every occuring
 * property struct. Afterwards the bodyMap is equally treated. And finally
 * the bodyMap and the propertyMap are deleted.
 */
GeometryContext::~GeometryContext() {
  TRACE();
  for(geometryPropertyIterator i = m_geometryPropertyMap->begin(); i != m_geometryPropertyMap->end(); i++) {
    DEBUG("GeometryProperty::clear()  name : " << i->second->name, MAIA_DEBUG_IO);
    DEBUG("GeometryProperty::clear()  segmentId : " << i->second->segmentId, MAIA_DEBUG_IO);
    i->second->clear(); // remove the pointers of property
  }
  m_geometryPropertyMap->clear();

  for(bodyIterator i = m_bodyMap->begin(); i != m_bodyMap->end(); i++) {
    DEBUG("GeometryProperty::clear()  name : " << i->second->name, MAIA_DEBUG_IO);
    DEBUG("GeometryProperty::clear()  noSegments : " << i->second->noSegments, MAIA_DEBUG_IO);
    delete[] i->second->segments;
  }
}


//---------------------------------------------------------------------------
/** \brief This method adds properties
 *
 *
 *
 *
 */
void GeometryContext::addProperty(GeometryProperty* p) {
  TRACE();

  const pair<const MString, GeometryProperty*> mp(p->name, p);
  m_geometryPropertyMap->insert(mp);
}

//---------------------------------------------------------------------------
/** \brief This intializes the property Map
 *
 * Use ts function if you want to create a new property list. If you want
 * to load properties from file use the function readPropertyFile.
 *
 */
void GeometryContext::init() {
  //  TRACE();

  m_geometryPropertyMap = new geometryPropertyMap;
}

//---------------------------------------------------------------------------
/** \brief This function writes all properties in a property File
 *
 * Use ts function if you want to create a new property list. If you want
 * to load properties from file use the function readPropertyFile.
 *
 */
void GeometryContext::writeProperties(MChar* fileName) {
  TRACE();

  m_geometryIOBase = new GeometryIONetcdf(mpiComm());
  m_geometryIOBase->writeProperties(fileName, m_geometryPropertyMap);
}

//--------------------------------------------------------------------------
/**
 *
 *
 */
MBool GeometryContext::propertyExists(MString name, MInt /*solverId*/) {
  if(m_geometryPropertyMap->find(name) == m_geometryPropertyMap->end()) {
    return false;
  } else {
    return true;
  }
}
