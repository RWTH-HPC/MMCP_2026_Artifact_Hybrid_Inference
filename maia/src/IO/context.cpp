// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

/* Copyright (C) Insitute of Aerodynamics, RWTH Aachen - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Contact: office@aia.rwth-aachen.de
 */

#include "context.h"

#include "COMM/mpioverride.h"
#include "UTIL/binarytree.h"
#include "globals.h"
#include "ionetcdf.h"
#include "iotoml.h"
#include "typetraits.h"

using namespace std;

#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
// Anonymous namespace to prevent accidental access from outside of this translation unit
namespace {
ofstream maia_properties;
}
#endif
#ifdef MAIA_PRINT_PROPERTIES
namespace {
ofstream propsFromFile;
ofstream propsDefault;
} // namespace
#endif
/* definition of static members */
propertyMap* Context::m_propertyMap;
propertyMap* Context::m_propertyMapLowercase;
zoneMap* Context::m_zoneMap;
assembly* Context::m_assembly;
IONetcdf* Context::m_IONetcdf;
IOToml* Context::m_IOToml;
pair<propertyMap::iterator, propertyMap::iterator> Context::m_pair;
MInt Context::m_noSolvers;
MString Context::m_name;
MString Context::m_propertyFileOutputName;
MString Context::m_propertyFileText;
MInt Context::m_fileType;
MBool Context::m_checkingProperties;

//--------------------------------------------------------------------------
/** Creates an object of the specified filetype to open the property file
 *
 *
 */
void Context::readPropertyFile(FileType fileType, const MString& fileName) {
  TRACE();
  m_name = fileName;
  m_checkingProperties = false;
  if(fileType == NETCDF) {
    m_fileType = NETCDF;
    m_IONetcdf = new IONetcdf;
    m_assembly = m_IONetcdf->readPropertyFile(m_name);
    m_propertyMap = m_assembly->properties;
    m_propertyMapLowercase = m_assembly->propertiesLowercase;
    m_zoneMap = m_assembly->zones;
    m_noSolvers = m_IONetcdf->solverCount();
  } else if(fileType == TOML) {
    m_fileType = TOML;
    m_IOToml = new IOToml;
    m_assembly = m_IOToml->readPropertyFile(m_name);
    m_propertyMap = m_assembly->properties;
    m_propertyMapLowercase = m_assembly->propertiesLowercase;
    m_zoneMap = m_assembly->zones;
    m_noSolvers = m_IOToml->solverCount();
    m_propertyFileText = m_IOToml->rawText();
  }


/*! \page propertyPage1
    \section propertyFilename
    <code>MString Context::m_propertyFileOutputName</code> \n
    default = <code>"access_properties"</code> \n \n
    name of the file, in which the Property accesses are logged during the solver run and in which a property table is
   printed at the end. \n Keywords: <i>INPUT_OUTPUT</i>
*/
#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
  MInt logId = 0;

  // avoid overwriting previous log files when running a testcase with canary
  m_propertyFileOutputName.assign("access_properties_domain" + to_string(0) + "_0");
  std::ifstream infile(m_propertyFileOutputName);
  MBool fileExists = infile.good();
  while(fileExists) {
    // check if file for rank 0 exists (testcases with varying number of domains)
    m_propertyFileOutputName.assign("access_properties_domain" + to_string(0) + "_" + to_string(logId));
    std::ifstream infile2(m_propertyFileOutputName);
    fileExists = infile2.good();
    if(fileExists) logId++;
  }

  // use separate file for each rank
  m_propertyFileOutputName.assign("access_properties_domain" + to_string(globalDomainId()) + "_" + to_string(logId));
  MPI_Barrier(MPI_COMM_WORLD, AT_);
#endif
#ifdef MAIA_PRINT_PROPERTIES
  cout << "WARNING: MAIA_PRINT_PROPERTIES is enabled. This can slow down the simulation extremely. Furthermore, this "
          "should only be done for small testcases to reduce the output (every process needs to write out its "
          "properties). If you dont want to print out the properties, disable MAIA_PRINT_PROPERTIES in config.h!!!"
       << endl;
  propsFromFile.open("PropertiesFromFile.dat");
  propsDefault.open("DefaultProperties.dat");
#endif
}

//-----------------Get Properties, New approach (Begin)--------------------//
/** \brief Returns the value of the requested basic level property (implicit).
 * \details This template function searchs for a certain property given as an argument (name).
 * If the property exists, it will be returned. If the property does not exist, the default
 * value will be returned (if one is given). If there is also no defaultValue, the code exits
 * with an error.
 * \author Moritz Waldmann
 * \date 12.2018
 */

template <typename T, typename F>
T Context::getBasicPropertyImplicit(const MString& nameOfCallingFunction, const MString& name, MBool has_default,
                                    const T* defaultValue, MInt position, F&& f) {
  if(propertyExists(name)) {
    // lookup all properties with the requested key
    m_pair = m_propertyMap->equal_range(name);
    for(propertyIterator i = m_pair.first; i != m_pair.second; i++) {
      if(i->second->solverId == m_noSolvers) {
#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
        if(!m_checkingProperties) { // Count number of accesses (if not during property checking)
          i->second->noAccesses++;
        }
#endif

        const T pValue = f(i->second, position);

#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
        // Check for default property values specified in the property file
        if(has_default && !m_checkingProperties) {
          const T dValue = *defaultValue;
          const MBool isDefault = isSameValue(dValue, pValue);
          if(!isDefault) {
            // TODO labels:IO could be used to count properties with/without default values
            /* i->second->noOldAccesses++; */
          } else {
            std::cerr << "Property default value specified: " << name << " solverId=" << -1 << " pos=" << position
                      << " value=" << pValue << " default=" << *defaultValue << " " << nameOfCallingFunction
                      << std::endl;
          }
        }
#endif

        return pValue;
      }
    }
  }
  stringstream errormessage;
  errormessage << " The solver property '" << name << "' is  required for the simulation"
               << ", but was not set in the toml file!!!!";

  if(!has_default) {
    mTerm(1, nameOfCallingFunction, errormessage.str());
  } else {
    return *defaultValue;
  }
}

/** \brief Returns the value of the requested basic property. overloaded
 * \details The following 4 (one for MInt, MFloat, MString, MBool) functions
 * provide a value for the template funtion getBasicProperty.
 * They receive their value from the implicite function getBasicPropertyImplicit.
 * For debugging 'MAIA_PRINT_PROPERTIES' can be defined in config.h. Each of the
 * overloaded functions will than print the property name and corresponding value or default
 * value to a file. This can extend the simulation time.
 * \author Moritz Waldmann
 * \date 12.2018
 */

MString Context::getBasicPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                            MBool has_default, const MString* default_value, MInt pos) {
  MString propertyValue;
  propertyValue = getBasicPropertyImplicit(nameOfCallingFunction, name, has_default, default_value, pos,
                                           [](MProperty* i, MInt p) { return (*(i)->asString(p)); });
#ifdef MAIA_PRINT_PROPERTIES
  if(propertyExists(name))
    propsFromFile << name << " = " << propertyValue << ", basic Property (MString)" << endl;
  else
    propsDefault << name << " = " << propertyValue << ", basic Property (MString)" << endl;
#endif
  return propertyValue;
}

MInt Context::getBasicPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name, MBool has_default,
                                         const MInt* default_value, MInt pos) {
  MInt propertyValue;
  propertyValue = getBasicPropertyImplicit(nameOfCallingFunction, name, has_default, default_value, pos,
                                           [](MProperty* i, MInt p) { return (*(i)->asInt(p)); });
#ifdef MAIA_PRINT_PROPERTIES
  if(propertyExists(name))
    propsFromFile << name << " = " << propertyValue << ", basic Property (MInt)" << endl;
  else
    propsDefault << name << " = " << propertyValue << ", basic Property (MInt)" << endl;
#endif

  return propertyValue;
}

MFloat Context::getBasicPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name, MBool has_default,
                                           const MFloat* default_value, MInt pos) {
  MFloat propertyValue;
  propertyValue = getBasicPropertyImplicit(nameOfCallingFunction, name, has_default, default_value, pos,
                                           [](MProperty* i, MInt p) { return (*(i)->asFloat(p)); });
#ifdef MAIA_PRINT_PROPERTIES
  if(propertyExists(name))
    propsFromFile << name << " = " << propertyValue << ", basic Property (MFloat)" << endl;
  else
    propsDefault << name << " = " << propertyValue << ", basic Property (MFloat)" << endl;
#endif

  return propertyValue;
}

MBool Context::getBasicPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name, MBool has_default,
                                          const MBool* default_value, MInt pos) {
  MBool propertyValue;
  if(m_fileType == TOML) {
    propertyValue = getBasicPropertyImplicit(nameOfCallingFunction, name, has_default, default_value, pos,
                                             [](MProperty* i, MInt p) { return (*(i)->asBool(p)); });
  } else {
    // handle the case that we have a netcdf property file
    MInt tmpDef = 0; // false
    if(has_default) {
      tmpDef = static_cast<MInt>(*default_value);
    }
    MInt tmp = getBasicPropertyImplicit(nameOfCallingFunction, name, has_default, &tmpDef, pos,
                                        [](MProperty* i, MInt p) { return (*(i)->asInt(p)); });
    propertyValue = (tmp == 1);
  }
#ifdef MAIA_PRINT_PROPERTIES
  if(propertyExists(name))
    propsFromFile << name << " = " << propertyValue << ", basic Property (MBool)" << endl;
  else
    propsDefault << name << " = " << propertyValue << ", basic Property (MBool)" << endl;
#endif

  return propertyValue;
}


/** \brief Returns the value of the requested solver property (implicit).
 * \details This template function searchs for a certain property of a certain solver given as
 * an argument (name, solverId). If the property exists, it will be returned. If the property
 * does not exist, the basic property will be returned. If also the basic property does not exist,
 * the default value will be returned (if one is given). If there is also no defaultValue, the
 * code exits with an error.
 * \author Moritz Waldmann
 * \date 12.2018
 */

template <typename T, typename F>
T Context::getSolverPropertyImplicit(const MString& nameOfCallingFunction, const MString& name, MInt solverId,
                                     MBool has_default, const T* defaultValue, MInt position, F&& f) {
  if(solverPropertyExists(name, solverId)) {
    // lookup all properties with the requested key
    m_pair = m_propertyMap->equal_range(name);
    for(propertyIterator i = m_pair.first; i != m_pair.second; i++) {
      if(i->second->solverId == solverId) {
#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
        if(!m_checkingProperties) { // Count number of accesses (if not during property checking)
          i->second->noAccesses++;
        }
#endif

        const T pValue = f(i->second, position);

#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
        // Check for default property values specified in the property file
        if(has_default && !m_checkingProperties) {
          const T dValue = *defaultValue;
          const MBool isDefault = isSameValue(dValue, pValue);
          if(!isDefault) {
            // TODO labels:IO could be used to count properties with/without default values
            /* i->second->noOldAccesses++; */
          } else {
            std::cerr << "Property default value specified: " << name << " solverId=" << solverId << " pos=" << position
                      << " value=" << pValue << " default=" << *defaultValue << " " << nameOfCallingFunction
                      << std::endl;
          }
        }
#endif

        return pValue;
      }
    }
  } else if(propertyExists(name)) {
    return getBasicPropertyImplicit(nameOfCallingFunction, name, has_default, defaultValue, position, f);
  }

  stringstream errormessage;
  errormessage << " The solver property '" << name << "' is  required for the simulation"
               << ", but was not set in the toml file!!!!";
  if(!has_default) {
    mTerm(1, nameOfCallingFunction, errormessage.str());
  } else {
    return *defaultValue;
  }
}

/** \brief Returns the value of the requested solver property. overloaded
 * \details The following 4 (one for MInt, MFloat, MString, MBool) functions
 * provide a value for the template funtion getSolverProperty.
 * They receive their value from the implicite function getSolverPropertyImplicit
 * \author Moritz Waldmann
 * For debugging 'MAIA_PRINT_PROPERTIES' can be defined in config.h. Each of the
 * overloaded functions will than print the property name, the solverId and the corresponding
 * value or default value to a file. This can extend the simulation time.
 * \date 12.2018
 */

MString Context::getSolverPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                             const MInt solverId, MBool has_default, const MString* default_value,
                                             MInt pos) {
  MString propertyValue;
  propertyValue = getSolverPropertyImplicit(nameOfCallingFunction, name, solverId, has_default, default_value, pos,
                                            [](MProperty* i, MInt p) { return (*(i)->asString(p)); });
#ifdef MAIA_PRINT_PROPERTIES
  if(solverId != m_noSolvers) {
    if(propertyExists(name, solverId))
      propsFromFile << name << " = " << propertyValue << ", solver Property (solverId = " << solverId
                    << ", type = MString)" << endl;
    else
      propsDefault << name << " = " << propertyValue << ", solver Property (solverId = " << solverId
                   << ", type = MString)" << endl;
  }
#endif
  return propertyValue;
}

MInt Context::getSolverPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                          const MInt solverId, MBool has_default, const MInt* default_value, MInt pos) {
  MInt propertyValue;
  propertyValue = getSolverPropertyImplicit(nameOfCallingFunction, name, solverId, has_default, default_value, pos,
                                            [](MProperty* i, MInt p) { return (*(i)->asInt(p)); });
#ifdef MAIA_PRINT_PROPERTIES
  if(solverId != m_noSolvers) {
    if(propertyExists(name, solverId))
      propsFromFile << name << " = " << propertyValue << ", solver Property (solverId = " << solverId
                    << ", type = MInt)" << endl;
    else
      propsDefault << name << " = " << propertyValue << ", solver Property (solverId = " << solverId << ", type = MInt)"
                   << endl;
  }
#endif
  return propertyValue;
}

MFloat Context::getSolverPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                            const MInt solverId, MBool has_default, const MFloat* default_value,
                                            MInt pos) {
  MFloat propertyValue;
  propertyValue = getSolverPropertyImplicit(nameOfCallingFunction, name, solverId, has_default, default_value, pos,
                                            [](MProperty* i, MInt p) { return (*(i)->asFloat(p)); });
#ifdef MAIA_PRINT_PROPERTIES
  if(solverId != m_noSolvers) {
    if(propertyExists(name, solverId))
      propsFromFile << name << " = " << propertyValue << ", solver Property (solverId = " << solverId
                    << ", type = MFloat)" << endl;
    else
      propsDefault << name << " = " << propertyValue << ", solver Property (solverId = " << solverId
                   << ", type = MFloat)" << endl;
  }
#endif
  return propertyValue;
}

MBool Context::getSolverPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                           const MInt solverId, MBool has_default, const MBool* default_value,
                                           MInt pos) {
  MBool propertyValue;
  if(m_fileType == TOML) {
    propertyValue = getSolverPropertyImplicit(nameOfCallingFunction, name, solverId, has_default, default_value, pos,
                                              [](MProperty* i, MInt p) { return (*(i)->asBool(p)); });
  } else {
    // handle the case that we have a netcdf property file
    MInt tmpDef = 0; // false
    if(has_default) {
      tmpDef = static_cast<MInt>(*default_value);
    }
    MInt tmp = getBasicPropertyImplicit(nameOfCallingFunction, name, has_default, &tmpDef, pos,
                                        [](MProperty* i, MInt p) { return (*(i)->asInt(p)); });
    propertyValue = (tmp == 1);
  }

#ifdef MAIA_PRINT_PROPERTIES
  if(solverId != m_noSolvers) {
    if(propertyExists(name, solverId))
      propsFromFile << name << " = " << propertyValue << ", solver Property (solverId = " << solverId
                    << ", type = MBool)" << endl;
    else
      propsDefault << name << " = " << propertyValue << ", solver Property (solverId = " << solverId
                   << ", type = MBool)" << endl;
  }
#endif
  return propertyValue;
}
//-----------------Get Properties, New approach (End)---------------------//
//--------------------------------------------------------------------------

/** This method  deletes the property structs.
 * It runs over all elements of the propertyMap and zoneMap and deletes the
 * individual structs.
 * Finally the zoneMap and propertyMap are deleted.
 */
void Context::clear() {
  TRACE();
  for(propertyIterator i = m_propertyMap->begin(); i != m_propertyMap->end(); i++) {
    DEBUG("MProperty::clear:  name " << i->second->name, MAIA_DEBUG_IO);
    DEBUG("MProperty::clear:  solverId " << i->second->solverId, MAIA_DEBUG_IO);
    i->second->clear(); // remove the pointers of property
  }
  m_propertyMap->clear();

  for(propertyIterator i = m_propertyMapLowercase->begin(); i != m_propertyMapLowercase->end(); i++) {
    DEBUG("MProperty::clear:  name " << i->second->name, MAIA_DEBUG_IO);
    DEBUG("MProperty::clear:  solverId " << i->second->solverId, MAIA_DEBUG_IO);
    i->second->clear(); // remove the pointers of property
  }
  m_propertyMapLowercase->clear();

  for(zoneIterator i = m_zoneMap->begin(); i != m_zoneMap->end(); i++) {
    DEBUG("MProperty::clear()  name : " << i->second->name, MAIA_DEBUG_IO);
    DEBUG("MProperty::clear()  noSolvers : " << i->second->noSolvers, MAIA_DEBUG_IO);
    delete[] i->second->solvers;
  }
#ifdef MAIA_PRINT_PROPERTIES
  propsFromFile.close();
  propsDefault.close();
#endif
}


//---------------------------------------------------------------------------
/** \brief This method adds properties
 *
 *
 *
 *
 */
void Context::addProperty(MProperty* p) {
  //  TRACE();

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

//---------------------------------------------------------------------------
/** \brief This intializes the property Map
 *
 * Use ts function if you want to create a new property list. If you want
 * to load properties from file use the function readPropertyFile.
 *
 */
void Context::init() {
  //  TRACE();

  m_propertyMap = new propertyMap;
  m_propertyMapLowercase = new propertyMap;
}

//--------------------------------------------------------------------------
/** \brief This function checks if a property exists in general.
 *
 *
 */
MBool Context::propertyExists(const MString& name, MInt NotUsed(solverId)) {
  const MBool exists = !(m_propertyMap->find(name) == m_propertyMap->end());
  if(!exists) {
    // produce lowercase name to check for case-insensitive matches
    MString nameL = name;
    std::transform(nameL.begin(), nameL.end(), nameL.begin(), [](unsigned char c) { return std::tolower(c); });
    if(nameL == "postprocessing") { // delete this, when the two different spellings of postprocessing are delted!
      return exists;
    }
    const MBool existsCaseInsensitive = !(m_propertyMapLowercase->find(nameL) == m_propertyMapLowercase->end());
    if(existsCaseInsensitive) {
      mTerm(1, AT_,
            "Property " + name
                + " not found, but there was a case-insensitive match. Are you sure you spelled it correctly?");
    }
  }
  return exists;
}

/** \brief Checks existence of a solver property
 * details This function returns true, if the solver property with the requested name
 * exists for the requested solver.
 * \author Moritz Waldmann
 * \date 12.2018
 */
MBool Context::solverPropertyExists(const MString& name, MInt solverId) {
  if(propertyExists(name)) {
    m_pair = m_propertyMap->equal_range(name);                        // lookup all properties with the requested key
    for(propertyIterator i = m_pair.first; i != m_pair.second; i++) { // search for the solverId
      if(i->second->solverId == solverId) {                           // look for the requested solverId
        return true;
      }
    }
  }
  return false;
}

/** \brief Returns the number of elements of a property
 * \details This function returns the number of elements of a property.
 * If the property exists, it returns the value given by the function count(), which
 * is implemented in property.cpp
 * \author Moritz Waldmann
 * \date 12.2018
 */
MInt Context::propertyLength(const MString& name, MInt solverId) {
  if(solverPropertyExists(name, solverId)) {
    m_pair = m_propertyMap->equal_range(name);                        // lookup all properties with the requested key
    for(propertyIterator i = m_pair.first; i != m_pair.second; i++) { // search for the solverId
      if(i->second->solverId == solverId) {                           // look for the requested solverId
        return i->second->count();
      }
    }
  } else if(propertyExists(name)) {
    m_pair = m_propertyMap->equal_range(name);                        // lookup all properties with the requested key
    for(propertyIterator i = m_pair.first; i != m_pair.second; i++) { // search for the solverId
      if(i->second->solverId == m_noSolvers) {                        // look for the requested solverId
        return i->second->count();
      }
    }
  }
  return 0;
}


/// \brief Check the length of a given property, terminate if it does not match the given length
void Context::assertPropertyLength(const MString& name, const MInt length, const MInt solverId) {
  const MInt propLength = propertyLength(name, solverId);
  if(propLength != length) {
    const MString isSolver = (solverId != m_noSolvers) ? " for solver #" + std::to_string(solverId) : "";
    TERMM(1, "Property '" + name + "'" + isSolver + " has the wrong length, it is of length "
                 + std::to_string(Context::propertyLength(name)) + " but there should be " + std::to_string(length)
                 + " values!");
  }
}


/** \brief Returns the number of elements of a property
 * \details This function returns the number of elements of a property.
 * If the property exists, it returns its solverId
 * \author Moritz Waldmann
 * \date 12.2018
 */
MInt Context::propertySolverId(const MString& name, MInt solverId) {
  if(solverPropertyExists(name, solverId)) {
    m_pair = m_propertyMap->equal_range(name);                        // lookup all properties with the requested key
    for(propertyIterator i = m_pair.first; i != m_pair.second; i++) { // search for the solverId
      if(i->second->solverId == solverId) {                           // look for the requested solverId
        return i->second->solverId;
      }
    }
  } else if(propertyExists(name)) {
    m_pair = m_propertyMap->equal_range(name);                        // lookup all properties with the requested key
    for(propertyIterator i = m_pair.first; i != m_pair.second; i++) { // search for the solverId
      if(i->second->solverId == m_noSolvers) {                        // look for the requested solverId
        return i->second->solverId;
      }
    }
  }
  return 0;
}

//---------------------------------------------------------------------------
/** \brief Write the properties into a text file
 * \details This function closes the properties stream and writes all properties in the propertyMap to the same file
 * including extra information like the number of accesses \author Christoph Siewert \date 10.2010
 */
void Context::writePropertiesHumanReadable() {
  TRACE();

  // All domains write their own property access file
  {
    std::cerr << "rank " << globalDomainId() << " writing property access file. " << m_propertyFileOutputName
              << std::endl;
#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
    maia_properties.close();
#endif
    FILE* datei;
    datei = fopen(m_propertyFileOutputName.c_str(), "a+");
    MString tester;
    MChar buffer[100];
    MInt maxNameLength, maxElementLength;
    maxNameLength = maxElementLength = 0;
    fprintf(datei, "\n\n--------------------------property table--------------------------\n");
    for(auto& i : *m_propertyMap) {
      if(maxNameLength < ((MInt)i.first.length())) {
        maxNameLength = i.first.length();
      }
      for(MInt j = 0; j < i.second->count(); j++) {
        switch(i.second->type()) {
          case MBOOL: {
            sprintf(buffer, "%i", *i.second->asBool(j));
            tester.append(buffer);
            break;
          }
          case MINT: {
            sprintf(buffer, "%i", *i.second->asInt(j));
            tester.append(buffer);
            break;
          }
          case MFLOAT: {
            sprintf(buffer, "%G", *i.second->asFloat(j));
            tester.append(buffer);
            break;
          }
          case MSTRING: {
            tester.append(*i.second->asString(j));
            tester.append("  ");
            break;
          }
          default: {
            mTerm(1, AT_, "writePropertiesHumanReadable(): switch variable 'i->second->type' not matching any case");
          }
        }
        tester.append("; ");
      }
      if(maxElementLength < ((MInt)tester.length())) {
        maxElementLength = tester.length();
      }
      tester.clear();
    }
    maxNameLength = maxNameLength + 1;

    for(auto& i : *m_propertyMap) {
      switch(i.second->type()) {
        case MBOOL: {
          fprintf(datei, "%s   ", "bool");
          fprintf(datei, "%s", i.first.c_str());
          for(MInt j = 0; j < ((MInt)(maxNameLength - i.first.length())); j++)
            fprintf(datei, " ");
          fprintf(datei, "= ");
          for(MInt j = 0; j < i.second->count(); j++) {
            sprintf(buffer, "%i; ", *i.second->asBool(j));
            tester.append(buffer);
            fprintf(datei, "%i; ", i.second->boolField[j]);
          }
          for(MInt k = 0; k < ((MInt)(maxElementLength - tester.length())); k++)
            fprintf(datei, " ");
          break;
        }
        case MINT: {
          fprintf(datei, "%s    ", "int");
          fprintf(datei, "%s", i.first.c_str());
          for(MInt j = 0; j < ((MInt)(maxNameLength - i.first.length())); j++)
            fprintf(datei, " ");
          fprintf(datei, "= ");
          for(MInt j = 0; j < i.second->count(); j++) {
            sprintf(buffer, "%i; ", *i.second->asInt(j));
            tester.append(buffer);
            fprintf(datei, "%i; ", i.second->intField[j]);
          }
          for(MInt k = 0; k < ((MInt)(maxElementLength - tester.length())); k++)
            fprintf(datei, " ");
          break;
        }
        case MFLOAT: {
          fprintf(datei, "%s ", "double");
          fprintf(datei, "%s", i.first.c_str());
          for(MInt j = 0; j < ((MInt)(maxNameLength - i.first.length())); j++)
            fprintf(datei, " ");
          fprintf(datei, "= ");
          for(MInt j = 0; j < i.second->count(); j++) {
            sprintf(buffer, "%G; ", *i.second->asFloat(j));
            tester.append(buffer);
            fprintf(datei, "%G; ", i.second->floatField[j]);
          }
          for(MInt k = 0; k < ((MInt)(maxElementLength - tester.length())); k++)
            fprintf(datei, " ");
          break;
        }
        case MSTRING: {
          fprintf(datei, "%s   ", "char");
          fprintf(datei, "%s", i.first.c_str());
          for(MInt j = 0; j < ((MInt)(maxNameLength - i.first.length())); j++)
            fprintf(datei, " ");
          fprintf(datei, "= ");
          for(MInt j = 0; j < i.second->count(); j++) {
            tester.append(*i.second->asString(j));
            fprintf(datei, "\"%s\"; ", i.second->asString(j)->c_str());
          }
          for(MInt k = 0; k < ((MInt)(maxElementLength - tester.length() - 4)); k++)
            fprintf(datei, " ");
          break;
        }

        default: {
          fprintf(datei, "--------------ERROR--------------");
          fprintf(datei, "\n");
          fprintf(datei, "Undefined type of Property: ");
          fprintf(datei, "%s =", i.first.c_str());
          fprintf(datei, "\n");
          break;
        }
      }
      const MInt solverId = (i.second->solverId < m_noSolvers) ? i.second->solverId : -1;
      fprintf(datei, "SolverId %2i", solverId);
      fprintf(datei, " Accesses %i", i.second->noAccesses);
      fprintf(datei, " OldAccesses %i",
              i.second->noOldAccesses); // TODO labels:IO,totest,toremove not useful (anymore)?
      fprintf(datei, "\n");
      tester.clear();
    }
    fclose(datei);
  }
}

//---------------------------------------------------------------------------
/** \brief Sets flag to forbide property access
 * \details not yet activate, remove the comment
 * \author Christoph Siewert
 * \date 10.2010
 */
void Context::initializationProcessFinished() {
  //  TRACE();

#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
  maia_properties << endl << "--------------------------initialization process finshed--------------------------";
  maia_properties << endl
                  << "--------------------------property access not allowed anymore--------------------------" << endl
                  << endl;
#endif
  // chrs should be commented in, if properties regrouping is finished
  // writePropertiesHumanReadable();
}

//---------------------------------------------------------------------------
/** \brief Communicates properties to check if default falues match
 * \author Andreas Lintermann
 * \date 06.12.2010
 *
 * This function communicates properties available in the property map
 * in a tree-like structure and is called in environment->run().
 *
 * 1.  A binary tree is built to accelerate communication
 *
 * 2.  Starting at the leafs, properties are merged at father-nodes
 * 2.1 Check if properties from leafs do not violate equality restrictions
 *
 * 3. Finish at root node (processor with id 0)
 */
void Context::communicateProperties() {
  TRACE();

  NEW_TIMER_GROUP_STATIC(m_timerGroupId, "Context");
  NEW_TIMER(timerpcheck, "Property Validity Check", m_timerGroupId);
  RECORD_TIMER_START(timerpcheck);

  m_log << "==============================================================" << endl;
  m_log << "  Checking property validity... ";
#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
  maia_properties << "==============================================================" << endl;
  maia_properties << "  Checking property validity... ";
#endif

  // TEST
  /*
  MInt rank = MPI_COMM_WORLD.Get_rank();
  if(rank==1)
    {
      MInt* b = new MInt();
      *b=1;
      MProperty* a = new MProperty(1,"AAA",1,b);
      addProperty(a);
    }
  else
    {
      MInt* b = new MInt();
      *b=2;
      MProperty* a = new MProperty(1,"AAA",1,b);
      addProperty(a);
    }
  */

  MInt noDomains = globalNoDomains();
  BinaryPropertyMPITree tree(noDomains);
  // if(rank==0)
  // tree.printTree();

  MInt comm_partner = 0;

  list<MProperty*>* left = nullptr;
  list<MProperty*>* right = nullptr;

  // We need to receive something from below
  if(tree.getMyLeftMPISender() != nullptr || tree.getMyRightMPISender() != nullptr) {
    if(tree.getMyLeftMPISender() != nullptr) {
      comm_partner = *(tree.getMyLeftMPISender()->getObject());
      left = receiveProperties(comm_partner);
      checkPropertyViolation(comm_partner, left);
    }
    if(tree.getMyRightMPISender() != nullptr) {
      comm_partner = *(tree.getMyRightMPISender()->getObject());
      right = receiveProperties(comm_partner);
      checkPropertyViolation(comm_partner, right);
    }
  }

  // Send if we are not the almighty father of the tree
  if(tree.getMyMPIReceiver() != nullptr) {
    comm_partner = *(tree.getMyMPIReceiver()->getObject());
    sendProperties(comm_partner);
  }

  m_log << " --->  EVERYTHING IS FINE!!!!" << endl;
  m_log << "==============================================================" << endl;
#ifdef MAIA_WRITE_ACCESS_PROPERTIES_FILE
  maia_properties << " --->  EVERYTHING IS FINE!!!!" << endl;
  maia_properties << "==============================================================" << endl;
#endif
  RECORD_TIMER_STOP(timerpcheck);
}

/** \brief Receives a property list from another cpu.
 * \author Andreas Lintermann
 * \date 9.12.2010
 * \param[in] the rank of the processor to commuicate with
 * \return returning a list containing all properties from the sending cpu
 *
 * This function receives the properties from a given cpu in the
 * following order
 *
 * 1. the sizes of all arrays
 * 2. all integer values of the properties
 * 3. all float (double) values of the properties
 * 4. all string values of the properties
 * 5. all names of the properties
 * 6. all solverIds of the properties
 * 7. all offsets of the properties, required to find the right
 *    values for the associated property
 *
 * At the end a list of of received properties is built and returned for
 * further investigations.
 */
list<MProperty*>* Context::receiveProperties(MInt rank) {
  TRACE();
  // cerr << "Processor " <<  MPI_COMM_WORLD.Get_rank() << " is receiving from: " << rank << endl;

  MInt tag = 0;
  MPI_Status status;

  MInt listSizes[10];

  // Receive list sizes
  MPI_Recv(listSizes, 10, MPI_INT, rank, tag, MPI_COMM_WORLD, &status, AT_, "listSizes");

  // Receive integer array
  auto* int_ar = new MInt[listSizes[0]];
  MPI_Recv(int_ar, listSizes[0], MPI_INT, rank, tag, MPI_COMM_WORLD, &status, AT_, "int_ar");

  // Receive MBool array
  auto* bool_ar = new MBool[listSizes[1]];
  MPI_Recv(bool_ar, listSizes[1], MPI_C_BOOL, rank, tag, MPI_COMM_WORLD, &status, AT_, "bool_ar");

  // Receive float array
  auto* float_ar = new MFloat[listSizes[2]];
  MPI_Recv(float_ar, listSizes[2], MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, &status, AT_, "float_ar");

  // Receive string array
  MChar* char_ar = new MChar[listSizes[3]];
  MPI_Recv(char_ar, listSizes[3], MPI_CHAR, rank, tag, MPI_COMM_WORLD, &status, AT_, "char_ar");

  // Receive name array
  MChar* char_ar2 = new MChar[listSizes[4]];
  MPI_Recv(char_ar2, listSizes[4], MPI_CHAR, rank, tag, MPI_COMM_WORLD, &status, AT_, "char_ar2");

  // Receive offset array
  MInt offset_sizes = listSizes[5] + listSizes[6] + listSizes[7] + listSizes[8] + listSizes[9];
  auto* offset = new MInt[offset_sizes];
  MPI_Recv(offset, offset_sizes, MPI_INT, rank, tag, MPI_COMM_WORLD, &status, AT_, "offset");

  // Receive offset array
  auto* solverIds = new MInt[listSizes[9]];
  MPI_Recv(solverIds, listSizes[9], MPI_INT, rank, tag, MPI_COMM_WORLD, &status, AT_, "solverIds");

  // Build something we can give back and handle in the comparison
  auto* proplist = new list<MProperty*>();

  MInt boolOffset = listSizes[5];
  MInt floatOffset = listSizes[5] + listSizes[6];
  MInt stringOffset = listSizes[5] + listSizes[6] + listSizes[7];
  MInt nameOffset = listSizes[5] + listSizes[6] + listSizes[7] + listSizes[8];

  // Run over all elements
  MInt i = 0;
  for(; i < boolOffset; i++) {
    MInt low_off = 0;
    MInt low_off_name = 0;
    if(i > 0) {
      low_off = offset[i - 1];
      low_off_name = offset[nameOffset + i - 1];
    }

    auto* a = new MInt[offset[i] - low_off];
    for(MInt j = low_off, l = 0; j < offset[i]; j++, l++) {
      a[l] = int_ar[j];
    }

    MString name;
    for(MInt j = low_off_name; j < offset[nameOffset + i]; j++) {
      name += char_ar2[j];
    }

    MProperty* prop = new MProperty(offset[i] - low_off, name, solverIds[i], a);
    proplist->push_back(prop);

    delete[] a;
  }

  for(; i < floatOffset; i++) {
    MInt low_off = 0;
    MInt low_off_name = offset[nameOffset + i - 1];
    if(i > boolOffset) {
      low_off = offset[i - 1];
    }

    // cerr << offset[i]<< " " << low_off << endl;
    auto* a = new MBool[offset[i] - low_off];
    for(MInt j = low_off, l = 0; j < offset[i]; j++, l++) {
      a[l] = bool_ar[j];
    }

    MString name;
    for(MInt j = low_off_name; j < offset[nameOffset + i]; j++) {
      name += char_ar2[j];
    }

    MProperty* prop = new MProperty(offset[i] - low_off, name, solverIds[i], a);
    proplist->push_back(prop);

    delete[] a;
  }

  for(; i < stringOffset; i++) {
    MInt low_off = 0;
    MInt low_off_name = offset[nameOffset + i - 1];
    if(i > floatOffset) low_off = offset[i - 1];

    // cerr << offset[i]<< " " << low_off << endl;
    auto* a = new MFloat[offset[i] - low_off];
    for(MInt j = low_off, l = 0; j < offset[i]; j++, l++) {
      a[l] = float_ar[j];
    }

    MString name;
    for(MInt j = low_off_name; j < offset[nameOffset + i]; j++)
      name += char_ar2[j];

    MProperty* prop = new MProperty(offset[i] - low_off, name, solverIds[i], a);
    proplist->push_back(prop);

    delete[] a;
  }

  for(; i < nameOffset; i++) {
    MInt low_off = 0;
    MInt low_off_name = offset[nameOffset + i - 1];
    if(i > stringOffset) low_off = offset[i - 1];

    MString a;
    for(MInt j = low_off, l = 0; j < offset[i]; j++, l++)
      a += char_ar[j];

    MString name;
    for(MInt j = low_off_name; j < offset[nameOffset + i]; j++)
      name += char_ar2[j];

    MProperty* prop = new MProperty(1, name, solverIds[i], &a);
    proplist->push_back(prop);
  }

  // Check all properties
  // FOR DEBUG
  /*
  MInt l=0;
  for (list<MProperty*>::iterator it = proplist->begin(); it!=proplist->end();it++,l++)
    {
      cerr << (*it)->name << " " << (*it)->solverId << ": " ;
      for (MInt j=0; j<(*it)->elements;j++)
  {
    if(l<floatOffset)
      cerr << *((*it)->asInt(j)) << " ";
    else if(l<stringOffset)
      cerr << *((*it)->asFloat(j)) << " ";
    else if(l<nameOffset)
      cerr << *((*it)->asString(j)) << " ";
  }
      cerr << endl;
    }
  */
  delete[] offset;
  delete[] solverIds;
  delete[] int_ar;
  delete[] float_ar;
  delete[] bool_ar;
  delete[] char_ar;
  delete[] char_ar2;


  return proplist;
}

/** \brief Sends the property list to another cpu.
 * \author Andreas Lintermann
 * \date 9.12.2010
 * \param[in] the rank of the processor to commuicate with
 *
 * This function sends the own properties to a cpu given by rank.
 * The properties are sent in the following order
 *
 * 1. the sizes of all arrays
 * 2. all integer values of the properties
 * 3. all float (double) values of the properties
 * 4. all string values of the properties
 * 5. all names of the properties
 * 6. all solverIds of the properties
 * 7. all offsets of the properties, required to find the right
 *    values for the associated property
 *
 */
void Context::sendProperties(MInt rank) {
  TRACE();

  MInt tag = 0;

  // cerr << "Processor " <<  MPI_COMM_WORLD.Get_rank() << " is sending to: " << rank << endl;
  MInt listSizes[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  list<MProperty*> intProperties;
  list<MProperty*> floatProperties;
  list<MProperty*> stringProperties;
  list<MProperty*> boolProperties;

  for(auto& i : *m_propertyMap) {
    switch(i.second->type()) {
      case MINT: {
        intProperties.push_back(i.second);
        listSizes[0] += (i.second->elements);
        listSizes[5]++;
        break;
      }
      case MBOOL: {
        boolProperties.push_back(i.second);
        listSizes[1] += (i.second->elements);
        listSizes[6]++;
        break;
      }
      case MFLOAT: {
        floatProperties.push_back(i.second);
        listSizes[2] += (i.second->elements);
        listSizes[7]++;
        break;
      }
      case MSTRING: {
        stringProperties.push_back(i.second);
        listSizes[3] += (i.second->stringField->size());
        listSizes[8]++;
        break;
      }
      default: {
        mTerm(1, AT_, "Context::sendProperties(): switch variable 'i->second->type' not matching any case");
      }
    }
    listSizes[4] += (i.second->name.size());
    listSizes[9]++;
  }

  MInt offset_sizes = listSizes[5] + listSizes[6] + listSizes[7] + listSizes[8] + listSizes[9];
  // cerr << offset_sizes << endl;
  auto* offset = new MInt[offset_sizes];
  auto* solverIds = new MInt[listSizes[9]];
  /*
  cerr << "Number of INT properties: " << intProperties.size() << ", size = " << listSizes[0] <<  endl;
  cerr << "Number of FLOAT properties: " << floatProperties.size() << ", size = " << listSizes[1] << endl;
  cerr << "Number of STRING properties: " << stringProperties.size() << ", size = " << listSizes[2] << endl;
  cerr << "Number of NAMES: " << listSizes[7] << ", size = " << listSizes[3] << endl;
  */

  // Build integer array and fill sizes for offsets
  auto* int_ar = new MInt[listSizes[0]];
  MInt k = 0, m = 0, n = 0;
  for(auto& intPropertie : intProperties) {
    for(MInt j = 0; j < intPropertie->elements; j++) {
      int_ar[k] = intPropertie->intField[j];
      k++;
    }
    offset[m] = k;
    solverIds[m] = intPropertie->solverId;
    m++;
  }

  // Build float array
  auto* bool_ar = new MBool[listSizes[1]];
  k = 0;
  for(auto& boolPropertie : boolProperties) {
    for(MInt j = 0; j < boolPropertie->elements; j++) {
      bool_ar[k] = boolPropertie->boolField[j];
      k++;
    }


    offset[m] = k;
    solverIds[m] = boolPropertie->solverId;
    m++;
  }

  // Build MBool array and fill sizes for offsets
  auto* float_ar = new MFloat[listSizes[2]];
  k = 0;
  for(auto& floatPropertie : floatProperties) {
    for(MInt j = 0; j < floatPropertie->elements; j++) {
      float_ar[k] = floatPropertie->floatField[j];
      k++;
    }

    offset[m] = k;
    solverIds[m] = floatPropertie->solverId;
    m++;
  }

  // Build string array
  MChar* char_ar = new MChar[listSizes[3]];
  k = 0;
  for(auto& stringPropertie : stringProperties) {
    MString a = *(stringPropertie->asString());
    for(MInt j = 0; j < (MInt)a.size(); j++) {
      char_ar[k] = a[j];
      k++;
    }
    offset[m] = k;
    solverIds[m] = stringPropertie->solverId;
    m++;
  }

  // Build name array
  MChar* char_ar2 = new MChar[listSizes[4]];
  k = 0;
  for(auto& intPropertie : intProperties) {
    MString a = intPropertie->name;
    for(MChar j : a) {
      char_ar2[k] = j;
      k++;
    }
    offset[m] = k;
    m++;
  }
  for(auto& boolPropertie : boolProperties) {
    MString a = boolPropertie->name;
    for(MChar j : a) {
      char_ar2[k] = j;
      k++;
    }
    offset[m] = k;
    m++;
  }
  for(auto& floatPropertie : floatProperties) {
    MString a = floatPropertie->name;
    for(MChar j : a) {
      char_ar2[k] = j;
      k++;
    }
    offset[m] = k;
    m++;
  }
  for(auto& stringPropertie : stringProperties) {
    MString a = stringPropertie->name;
    for(MChar j : a) {
      char_ar2[k] = j;
      k++;
    }
    offset[m] = k;
    m++;
    n++;
  }

  // Send information about arrays
  MPI_Send(listSizes, 10, MPI_INT, rank, tag, MPI_COMM_WORLD, AT_, "listSizes");

  // Send integer array
  MPI_Send(int_ar, listSizes[0], MPI_INT, rank, tag, MPI_COMM_WORLD, AT_, "int_ar");

  // Send float array
  MPI_Send(bool_ar, listSizes[1], MPI_C_BOOL, rank, tag, MPI_COMM_WORLD, AT_, "bool_ar");

  // Send MBool array
  MPI_Send(float_ar, listSizes[2], MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, AT_, "float_ar");

  // Send string array
  MPI_Send(char_ar, listSizes[3], MPI_CHAR, rank, tag, MPI_COMM_WORLD, AT_, "char_ar");

  // Send name array
  MPI_Send(char_ar2, listSizes[4], MPI_CHAR, rank, tag, MPI_COMM_WORLD, AT_, "char_ar2");

  // Send offset information
  MPI_Send(offset, offset_sizes, MPI_INT, rank, tag, MPI_COMM_WORLD, AT_, "offset");

  // Send solverIds
  MPI_Send(solverIds, listSizes[9], MPI_INT, rank, tag, MPI_COMM_WORLD, AT_, "solverIds");

  // Clean up
  delete[] offset;
  delete[] solverIds;
  delete[] int_ar;
  delete[] float_ar;
  delete[] bool_ar;
  delete[] char_ar;
  delete[] char_ar2;
  intProperties.clear();
  floatProperties.clear();
  stringProperties.clear();
  boolProperties.clear();
}


/** \brief Checks if local properties differ from received properties.
 * \author Andreas Lintermann
 * \date 9.12.2010
 * \param[in] the rank the partner cpu
 * \param[in] the property list of the other cpu
 *
 * This function checks if local properties differ from the properties
 * received from the children and stops the program if values differ.
 *
 */
void Context::checkPropertyViolation(MInt partnerrank, list<MProperty*>* prop) {
  TRACE();

  MInt defaultInt = -1;
  MFloat defaultFloat = -1.0;
  MBool defaultBool = false;
  MString defaultString = " ";

  MInt rank = globalDomainId();

  // Set flag to indicate that property accesses should not be counted here
  m_checkingProperties = true;

  // Check local with other
  for(auto& it : *prop) {
    if(propertyExists(it->name)) {
      switch(it->type()) {
        case MINT: {
          for(MInt j = 0; j < it->elements; j++) {
            if(*(it->asInt(j)) != getSolverProperty<MInt>(it->name, it->solverId, AT_, &defaultInt, j)
               || it->solverId != propertySolverId(it->name, it->solverId)) {
              stringstream errorMessage;
              //        auto p = getProperty((*it)->name,(*it)->solverId, AT_, (MInt*) nullptr);
              errorMessage << "Poperty violation for '" << it->name << "' on CPUs: " << rank << "," << partnerrank;
              //                     << endl
              //                     << "received: name='" << (*it)->name << "', value='" << *((*it)->asInt(j))
              //                     << "', solverId='" << (*it)->solverId << "'" << endl
              //                     << "local: name = '" << p->name << "', value='" << *(p->asInt(j))
              //                     << "', solverId='" << p->solverId << "'";
              mTerm(1, AT_, errorMessage.str());
            }
          }
          break;
        }
        case MBOOL: {
          for(MInt j = 0; j < it->elements; j++) {
            if(*(it->asBool(j)) != getSolverProperty<MBool>(it->name, it->solverId, AT_, &defaultBool, j)
               || it->solverId != propertySolverId(it->name, it->solverId)) {
              stringstream errorMessage;
              errorMessage << "Poperty violation for '" << it->name << "' on CPUs: " << rank << "," << partnerrank;
              mTerm(1, AT_, errorMessage.str());
            }
          }
          break;
        }
        case MFLOAT: {
          for(MInt j = 0; j < it->elements; j++) {
            if(!approx(*(it->asFloat(j)), getSolverProperty<MFloat>(it->name, it->solverId, AT_, &defaultFloat, j),
                       MFloatEps)
               || it->solverId != propertySolverId(it->name, it->solverId)) {
              stringstream errorMessage;
              errorMessage << "Poperty violation for '" << it->name << "' on CPUs: " << rank << "," << partnerrank;
              mTerm(1, AT_, errorMessage.str());
            }
          }
          break;
        }
        case MSTRING: {
          for(MInt j = 0; j < it->elements; j++) {
            if(*(it->asString(j)) != getSolverProperty<MString>(it->name, it->solverId, AT_, &defaultString, j)
               || it->solverId != propertySolverId(it->name, it->solverId)) {
              stringstream errorMessage;
              errorMessage << "Poperty violation for '" << it->name << "' on CPUs: " << rank << "," << partnerrank;
              mTerm(1, AT_, errorMessage.str());
            }
          }
          break;
        }
        default: {
          mTerm(1, AT_, "Context::CheckPropertyViolation(): switch variable '(*it)->type' not matching any case");
        }
      }
    } else {
      MProperty* a;
      switch(it->type()) {
        case MINT: {
          a = new MProperty(it->elements, it->name, it->solverId, it->intField);
          break;
        }
        case MBOOL: {
          a = new MProperty(it->elements, it->name, it->solverId, it->boolField);
          break;
        }
        case MFLOAT: {
          a = new MProperty(it->elements, it->name, it->solverId, it->floatField);
          break;
        }
        case MSTRING: {
          a = new MProperty(it->elements, it->name, it->solverId, it->stringField);
          break;
        }
        default:
          a = new MProperty();
      }
      addProperty(a);
    }
  }

  // Reset flag after checking properties
  m_checkingProperties = false;

  // Clean up
  prop->clear();
}

/// Dump all properties to text file
void Context::dump(const MString& fileName) {
  ofstream os(fileName);
  for(auto&& pm : *m_propertyMap) {
    const auto name = pm.first;
    const auto prop = pm.second;
    os << name << " (";
    switch(prop->propertyType) {
      case MINT:
        os << "MINT";
        break;
      case MFLOAT:
        os << "MFLOAT";
        break;
      case MSTRING:
        os << "MSTRING";
        break;
      case MBOOL:
        os << "MBOOL";
        break;
      default:
        TERMM(1, "Bad type");
    }
    os << "): ";
    for(MInt i = 0; i < prop->count(); i++) {
      if(i != 0) {
        os << ", ";
      }
      switch(prop->propertyType) {
        case MINT:
          os << prop->intField[i];
          break;
        case MFLOAT:
          os << prop->floatField[i];
          break;
        case MSTRING:
          os << prop->stringField[i];
          break;
        case MBOOL:
          os << prop->boolField[i];
          break;
        default:
          TERMM(1, "Bad type");
      }
    }
    os << std::endl;
  }
}
