// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometry3d.h"

#include <bitset>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_set>
#include "COMM/mpioverride.h"
#include "GRID/partition.h"
#include "IO/parallelio.h"
#include "MEMORY/scratch.h"
#include "geometryadt.h"
#include "geometrycontext.h"
#include "globals.h"

using namespace std;

Geometry3D::Geometry3D(const MInt solverId_, const MPI_Comm comm) : Geometry<3>(solverId_, comm) {
  TRACE();

  NEW_TIMER_GROUP(tg_geometry, "Geometry (solver #" + std::to_string(solverId_) + ")");
  m_tg_geometry = tg_geometry;
  NEW_TIMER(t_geometryAll, "complete geometry", m_tg_geometry);
  m_t_geometryAll = t_geometryAll;
  g_tc_geometry.emplace_back("complete geometry", m_t_geometryAll);
  RECORD_TIMER_START(m_t_geometryAll);

  NEW_SUB_TIMER(t_initGeometry, "init geometry", m_t_geometryAll);
  g_tc_geometry.emplace_back("init geometry", t_initGeometry);
  RECORD_TIMER_START(t_initGeometry);
  m_log << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl;
  m_log << "##                                             Initializing Geometry in 3D                               "
           "           ##"
        << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl
        << endl;

  otherCalls++;

  m_adt = 0;
  m_elements = 0;
  m_noElements = 0;
  m_parGeomMemFactor = 1.0;

  // new parallel geometry
  if(m_parallelGeometry) {
    if(Context::propertyExists("parallelGeomFileName", solverId()))
      /*! \page propertyPage1
       \section parallelGeomFileName
       <code>MString Geometry::m_parallelGeomFileName</code>\n
       default = <code> no default </code>\n \n
       defines new name for parallely generated 3D geometry \n
       possible values are:
       <ul>
       <li> any string</li>
       </ul>
       Keywords: PARALLEL, GEOMETRY
      */
      m_parallelGeomFileName = Context::getSolverProperty<MString>("parallelGeomFileName", solverId(), AT_);
    else {
      stringstream errorMsg;
      errorMsg
          << "ERROR: parallel geometry is activated but no parallelGeomFileName is specified in the properties file"
          << endl;
      m_log << errorMsg.str();
      mTerm(1, AT_, errorMsg.str());
    }
    if(Context::propertyExists("gridInputFileName", solverId()))
      m_gridFileName = Context::getSolverProperty<MString>("gridInputFileName", solverId(), AT_);
    else {
      stringstream errorMsg;
      errorMsg
          << "ERROR: parallel geometry is activated but no parallelGeomFileName is specified in the properties file"
          << endl;
      m_log << errorMsg.str();
      mTerm(1, AT_, errorMsg.str());
    }

    /*! \page propertyPage1
     \section parGeomMemFactor
     <code>MString Geometry::m_parGeomMemFactor</code>\n
     default = <code> 1.0 </code>\n \n
     Can be increased to increase the number of cells the collector is allocated with.\n
     possible values are:
     <ul>
     <li> any double > 1.0</li>
     </ul>
     Keywords: PARALLEL, GEOMETRY
    */
    m_parGeomMemFactor = Context::getSolverProperty<MFloat>("parGeomMemFactor", solverId(), AT_, &m_parGeomMemFactor);
  } else {
    /*! \page propertyPage1
     \section communicateSegmentsSerial
     <code>MBool Geometry3D::m_communicateSegmentsSerial </code>\n
     default = <code>true</code>\n \n
     If enabled the geometry files are read only on rank 0 and the geometry information is broadcasted to all other
     ranks, since on a large number of ranks the parallel read of the geometry files from all ranks at the same time can
     be very slow.\n
     Keywords: GEOMETRY, IO
    */
    m_communicateSegmentsSerial = true;
    m_communicateSegmentsSerial =
        Context::getSolverProperty<MBool>("communicateSegmentsSerial", solverId(), AT_, &m_communicateSegmentsSerial);
  }

  // moving boundary
  m_mbelements = 0;
  m_noMBElements = 0;
  m_GFieldInitFromSTL = false;
  m_forceBoundingBox = false;
  m_levelSetIntfBndId = 0;
  m_noLevelSetIntfBndIds = 0;

  m_GFieldInitFromSTL = Context::getSolverProperty<MBool>("GFieldInitFromSTL", solverId(), AT_, &m_GFieldInitFromSTL);

  if(m_GFieldInitFromSTL) {
    m_levelSetIntfBndId = Context::getSolverProperty<MInt>("levelSetIntfBndId", solverId(), AT_, &m_levelSetIntfBndId);

    m_forceBoundingBox = Context::getSolverProperty<MBool>("forcedBoundingBox", solverId(), AT_, &m_forceBoundingBox);


    if(Context::propertyExists("bodyBndryCndIds", solverId())
       || Context::propertyExists("GFieldInitFromSTLBndCndIds", solverId())) {
      const MInt noBodyBndIds = Context::propertyLength("bodyBndryCndIds", solverId());
      const MInt noBndIds = Context::propertyLength("GFieldInitFromSTLBndCndIds", solverId());
      m_noLevelSetIntfBndIds = noBodyBndIds + noBndIds;
      if(domainId() == 0) {
        cerr << " noLevelSetIntfBndIds: " << m_noLevelSetIntfBndIds << endl;
      }
      m_levelSetIntfBndIds = new MInt[m_noLevelSetIntfBndIds];

      for(MInt i = 0; i < noBodyBndIds; i++) {
        m_levelSetIntfBndIds[i] = Context::getSolverProperty<MFloat>("bodyBndryCndIds", solverId(), AT_, i);
        if(domainId() == 0) {
          cerr << " levelSetIntfBndIds: " << m_levelSetIntfBndIds[i] << endl;
        }
      }

      for(MInt i = noBodyBndIds; i < m_noLevelSetIntfBndIds; i++) {
        m_levelSetIntfBndIds[i] =
            Context::getSolverProperty<MInt>("GFieldInitFromSTLBndCndIds", solverId(), AT_, i - noBodyBndIds);
        if(domainId() == 0) {
          cerr << " levelSetIntfBndIds: " << m_levelSetIntfBndIds[i] << endl;
        }
      }

    } else {
      m_levelSetIntfBndIds = nullptr;
      //       for(MInt i=0; i < m_noLevelSetIntfBndIds; i++){
      //         m_levelSetIntfBndIds[i] = m_levelSetIntfBndId;
      //       }
    }
  }

  MString testcaseDir = "./";
  testcaseDir = Context::getSolverProperty<MString>("testcaseDir", solverId(), AT_, &testcaseDir);

  /*! \page propertyPage1
    \section inputDir
    <code>MString Geometry2/3D::Geometry2/3D()::inputDir</code>\n
    default = <code>"./"</code>\n \n
    Specify input path for geometry property file relative to testcaseDir. \n \n
    Possible values are:
    <ul>
      <li>relative path</li>
    </ul>
    Keywords: <i>GENERAL, IO, GEOMETRY, DIRECTORY</i>
  */
  MString inputDir = "./";
  inputDir = Context::getSolverProperty<MString>("inputDir", solverId(), AT_, &inputDir);
  /*! \page propertyPage1
    \section geometryInputFileName
    <code>MString geometry3d::tmpFileName </code>\n
    default = <code>no default value</code>\n \n
    Name of the geometry property file
    possible values are:
    <ul>
    <li> string </li>
    </ul>
    Keywords: <i> geometry </i>
  */
  MString tmpFileName =
      testcaseDir + inputDir + Context::getSolverProperty<MString>("geometryInputFileName", solverId(), AT_);

  RECORD_TIMER_STOP(t_initGeometry);

  NEW_SUB_TIMER(t_readGeomFile, "read geometry property file", m_t_geometryAll);
  g_t_readGeomFile = t_readGeomFile;
  g_tc_geometry.emplace_back("read geometry property file", g_t_readGeomFile);
  RECORD_TIMER_START(g_t_readGeomFile);
  m_log << "  + Reading geometry property file " << endl;
  if(tmpFileName.find(".toml") == MString::npos) {
    geometryContext().readPropertyFile(NETCDF, tmpFileName.c_str());
  } else {
    geometryContext().readPropertyFile(TOML, tmpFileName.c_str());
  }
  m_bodyMap = geometryContext().getBodies();
  RECORD_TIMER_STOP(g_t_readGeomFile);

  /*! \page propertyPage1
    \section gridCutTest
    <code>MString Geometry3D::m_gridCutTest </code>\n
    default = <code>SAT</code>\n \n
    Method to use for checking if a cell is intersected by a geometry element.\n
    Possible values are:
    <ul>
    <li> "SAT" </li>
    <li> something else </li>
    </ul>
    Keywords: <i>GEOMETRY, TRIANGLES, CUT CELLS</i>
  */
  m_gridCutTest = "SAT";
  m_gridCutTest = Context::getSolverProperty<MString>("gridCutTest", solverId(), AT_, &m_gridCutTest);

  m_bodyIt = m_bodyMap.begin();

  m_noSegments = geometryContext().getNoSegments();
  m_segmentOffsets.resize(m_noSegments + 1);
  m_segmentOffsetsWithoutMB.resize(m_noSegments + 1); // The offsets of MB segments are deleted.


  Geometry3D::readSegments();

  NEW_SUB_TIMER(t_correctVertices, "correct vertices", m_t_geometryAll);
  g_tc_geometry.emplace_back("correct vertices", t_correctVertices);
  RECORD_TIMER_START(t_correctVertices);
  Geometry3D::correctVertexCoordinates();
  RECORD_TIMER_STOP(t_correctVertices);


  NEW_SUB_TIMER(t_buildAdt, "building adt", m_t_geometryAll);
  g_tc_geometry.emplace_back("building adt", t_buildAdt);
  RECORD_TIMER_START(t_buildAdt);

  m_adt = new GeometryAdt<3>(this);
  if(m_noElements > 0) m_adt->buildTree();
  RECORD_TIMER_STOP(t_buildAdt);

  // moving boundary
  if(m_GFieldInitFromSTL && m_noElements > 0) m_adt->buildTreeMB();

  printMemoryUsage();
  RECORD_TIMER_STOP(m_t_geometryAll);

  if(m_debugParGeom) {
    DISPLAY_ALL_GROUP_TIMERS(m_tg_geometry);
    m_log << endl;
  }

  m_log << endl;

  Geometry3D::collectGlobalMemoryUsage();
}

/*
 * brief: generates an auxiliary geometry from an ASCII file that can be used for whatever
 *
 * Thomas Schilden, December 2016
 */
Geometry3D::Geometry3D(const MInt solverId_, const MString& filename, const MPI_Comm comm) : Geometry(solverId_, comm) {
  TRACE();

  m_log << "reading the auxiliary stl " << filename << endl;

  m_GFieldInitFromSTL = 0;

  m_noElements = 0;
  m_noMBElements = 0;

  Geometry3D::countSegmentTrianglesASCII(filename, &m_noElements);

  m_log << "  number of Elements: " << m_noElements << endl;

  m_elements = new Collector<element<nDim>>(m_noElements, nDim, 0);

  MInt bla = 0;
  Geometry3D::readSegmentTrianglesASCII(filename, m_elements, 0, 0, &bla);

  elements = &(m_elements->a[0]);

  Geometry3D::calculateBoundingBox();

  m_adt = new GeometryAdt<3>(this);
  m_adt->buildTree();
}

Geometry3D::~Geometry3D() {
  TRACE();

  delete m_adt;
  delete m_elements;

  // moving boundary
  if(m_mbelements != nullptr) delete m_mbelements;
}

/** \brief Rebuilds the ADT tree
 *
 * \author Andreas Lintermann
 * \date 25.09.2015
 *
 **/
void Geometry3D::rebuildAdtTree() {
  delete m_adt;

  m_adt = new GeometryAdt<3>(this);
  m_adt->buildTree();

  printMemoryUsage();
}

/** \brief prints the current memory usage of this object
 *
 * \author Andreas Lintermann
 * \date 05.10.2015
 *
 **/
void Geometry3D::printMemoryUsage() {
  TRACE();

  m_log << "  + Used memmory for collectors: "
        << ((MFloat)m_elements->memoryUseage() + (MFloat)m_adt->memoryUsage()) / 1024.0 / 1024.0
               + (m_GFieldInitFromSTL ? ((MFloat)m_elements->memoryUseage() / 1024.0 / 1024.0) : 0.0)
        << " MB" << endl;
  m_log << "    * triangles:    " << (MFloat)m_elements->memoryUseage() / 1024.0 / 1024.0 << " MB" << endl;
  if(m_GFieldInitFromSTL)
    m_log << "    * mb triangles: " << (MFloat)m_mbelements->memoryUseage() / 1024.0 / 1024.0 << " MB" << endl;
  m_log << "    * adt tree:     " << (MFloat)m_adt->memoryUsage() / 1024.0 / 1024.0 << " MB" << endl;
}

void Geometry3D::collectGlobalMemoryUsage() {
  TRACE();

  MInt count = 3;
  if(m_GFieldInitFromSTL) count = 4;

  MFloatScratchSpace mem(noDomains() * count, AT_, "mem");

  MFloatScratchSpace snd(count, AT_, "snd");
  snd[0] = m_noElements;
  snd[1] = (MFloat)m_elements->memoryUseage() / 1024.0 / 1024.0;
  snd[2] = (MFloat)m_adt->memoryUsage() / 1024.0 / 1024.0;
  if(m_GFieldInitFromSTL) snd[3] = (MFloat)m_mbelements->memoryUseage() / 1024.0 / 1024.0;

  MPI_Gather(snd.getPointer(), count, MPI_DOUBLE, mem.getPointer(), count, MPI_DOUBLE, 0, mpiComm(), AT_,
             "snd.getPointer()", "mem.getPointer()");

  if(domainId() == 0) {
    ofstream ofl;
    ofl.open("memrep_geom.txt", std::ios_base::out | std::ios_base::trunc);
    ofl.precision(10);
    ofl << "###########################################################################################################"
           "#"
        << endl
        << "#" << endl;
    ofl << "#   domain   no_triangles   triangles[MB]   tree[MB]   tmb_triagles[MP]" << endl << "#" << endl;

    MFloatScratchSpace sums(count, AT_, "sums");
    for(MInt i = 0; i < count; i++)
      sums[i] = 0.0;

    MFloatScratchSpace min(count, AT_, "min");
    MFloatScratchSpace max(count, AT_, "max");
    for(MInt i = 0; i < count; i++) {
      min[i] = 9999999999999999.9;
      max[i] = -1.0;
    }

    for(MInt d = 0; d < noDomains(); d++)
      for(MInt i = 0; i < count; i++) {
        sums[i] += mem[d * count + i];
        if(mem[d * count + i] > max[i]) max[i] = mem[d * count + i];
        if(mem[d * count + i] < min[i]) min[i] = mem[d * count + i];
      }

    MString names[4] = {"no_triangles", "triangles", "tree", "mb_triangles"};
    for(MInt i = 0; i < count; i++) {
      ofl << "#     SUM " << names[i] << " = " << sums[i] << endl;
      ofl << "#     AVG " << names[i] << " = " << sums[i] / noDomains() << endl;
      ofl << "#     MIN " << names[i] << " = " << min[i] << endl;
      ofl << "#     MAX " << names[i] << " = " << max[i] << endl << endl;
    }
    ofl << "#" << endl
        << "###########################################################################################################"
           "#"
        << endl;

    for(MInt d = 0; d < noDomains(); d++) {
      ofl << d << "\t\t";
      for(MInt i = 0; i < count; i++)
        ofl << mem[d * count + i] << "\t\t";
      ofl << endl;
    }
    ofl.close();
  }
}


/** \brief counts the number of triangles in an ASCII STL file
 *
 * \author Andreas Lintermann
 * \date 20.07.2015
 *
 * \param[in] fileName the name of the file to open
 * \param[in] the number of entries to be returned (gets nulled in function again)
 *
 **/
inline void Geometry3D::countSegmentTrianglesASCII(MString fileName, MInt* noElements) {
  TRACE();

  *noElements = 0;

  // open file in ascii format
  ifstream ifl(fileName);

  if(!ifl) {
    stringstream errorMessage;
    errorMessage << " ERROR in segment::readSegment, couldn't find file : " << fileName << "!";
    mTerm(1, AT_, errorMessage.str());
  }

  char tmp[256];
  string text;
  // If number of elements unknown, count them...
  while(text.find("endsolid") == std::string::npos) {
    ifl.getline(tmp, 256);
    text = tmp;
    if(text.find("facet") != std::string::npos && text.find("endface") == std::string::npos)
      *noElements = *noElements + 1;
  }
  ifl.close();
}

/** \brief counts the number of triangles in a BINARY STL file
 *
 * \author Andreas Lintermann
 * \date 20.07.2015
 *
 * This function detects the file size in byte and substracts the header information
 * (80 bytes) and the the information on the number of triangles (4 bytes). Each
 * triangle carries a total of 50 bytes.
 *
 * \param[in] fileName the name of the file to open
 * \param[in] the number of entries to be returned (gets nulled in function again)
 *
 **/
inline void Geometry3D::countSegmentTrianglesBINARY(MString fileName, MInt* noElements) {
  TRACE();

  *noElements = 0;

  struct stat stat_buf {};
  MInt rc = stat(fileName.c_str(), &stat_buf);
  if(rc < 0) {
    stringstream errorMessage;
    errorMessage << " ERROR in segment::readSegment, couldn't determine file size of: " << fileName << "!";
    mTerm(1, AT_, errorMessage.str());
  }

  *noElements = (stat_buf.st_size - (80 + 4)) / 50;
}


/** \brief counts the number of triangles in a NETCDF STL file
 *
 * \author Andreas Lintermann
 * \date 21.07.2015
 *
 * This is stored in the file in the variable noTriangles.
 *
 * \param[in] fileName the name of the file to open
 * \param[in] the number of entries to be returned (gets nulled in function again)
 *
 **/
inline void Geometry3D::countSegmentTrianglesNETCDF(MString fileName, MInt* noElements, const MPI_Comm comm) {
  TRACE();

  using namespace maia::parallel_io;

  *noElements = 0;

  ParallelIo surface(fileName, PIO_READ, comm);
  surface.readScalar(noElements, "noTriangles");
}


/** \brief reads triangles from an ASCII file
 *
 * \author Andreas Lintermann
 * \date 20.07.2015
 *
 * \param[in] fileName the name of the file to open
 * \param[in] elemCollector pointer to the collector holding the triangles
 * \param[in] bndCndId the boundary condition id associated with this segment
 * \param[in] segmenId the id of the segment
 * \param[in] the pointer to the segment offset for later computation
 *
 **/
inline void Geometry3D::readSegmentTrianglesASCII(MString fileName, Collector<element<3>>* elemCollector, MInt bndCndId,
                                                  MInt segmentId, MInt* offset) {
  TRACE();

  // open file in ascii format
  ifstream ifl(fileName);

  if(!ifl) {
    stringstream errorMessage;
    errorMessage << " ERROR in segment::readSegment, couldn't find file : " << fileName << "!";
    mTerm(1, AT_, errorMessage.str());
  }

  element<3>* allelements = elemCollector->a;
  string text;
  char tmp[1024]{};
  while(text.find("endsolid") == std::string::npos) {
    ifl.getline(tmp, 1024);
    text = tmp;
    if(text.find("normal") != std::string::npos) {
      elemCollector->append();
      trim(text);
      std::vector<std::string> tokens;
      tokenize(text, tokens, " ", true);

      // Read normal vector components
      for(MInt i = 0; i < 3; i++) {
        tokens[i + 2] = trim(tokens[i + 2]);
        if(tokens[i + 2].find_first_not_of("0123456789.eE-+") != std::string::npos) {
          cerr << "ERROR: normal component " << tokens[i + 2] << " is not a number " << endl;
          TERM(-1);
        }

        allelements[*offset].m_normal[i] = stod(tokens[i + 2]);
      }
      // Jump expression : outer loop
      ifl.getline(tmp, 256);

      // Read vertices

      for(MInt j = 0; j < 3; j++) {
        ifl.getline(tmp, 1024);

        text = tmp;
        trim(text);
        std::vector<std::string> vertex;
        tokenize(text, vertex, " ", true);

        for(MInt i = 0; i < 3; i++) {
          vertex[i + 1] = trim(vertex[i + 1]);
          if(vertex[i + 1].find_first_not_of("0123456789.eE-+") != std::string::npos) {
            cerr << "ERROR: vertex component " << vertex[i + 1] << " is not a number " << endl;
            TERM(-1);
          }


          allelements[*offset].m_vertices[j][i] = stod(vertex[i + 1]);
        }
      }

      allelements[*offset].boundingBox();
      allelements[*offset].m_bndCndId = bndCndId;
      allelements[*offset].m_segmentId = segmentId;
      allelements[*offset].m_originalId = -1;

      if(m_GFieldInitFromSTL)
        for(MInt id = 0; id < m_noLevelSetIntfBndIds; id++)
          if(bndCndId == m_levelSetIntfBndIds[id]) m_noMBElements++;

      *offset = *offset + 1;
    }
  }

  ifl.close();
}

/** \brief swaps 4 bytes from little to big endian
 *
 * \author Andreas Lintermann
 * \date 22.07.2015
 *
 * \param[in] buf a pointer to a char array of 4 bytes
 **/
inline void Geometry3D::swap4BytesToBE(char* buf) {
  char tmp_byte = buf[0];
  buf[0] = buf[3];
  buf[3] = tmp_byte;
  tmp_byte = buf[1];
  buf[1] = buf[2];
  buf[2] = tmp_byte;
}

/** \brief determines the CPU type big or little endian
 *
 * \author Andreas Lintermann
 * \date 22.07.2015
 *
 * \return if the system is big endian or not
 **/
MInt Geometry3D::is_big_endian(void) {
  union {
    uint32_t i;
    char c[4];
  } bint = {0x01020304};

  return static_cast<MInt>(bint.c[0] == 1);
}

/** \brief reads triangles from a BINARY file and swaps to big endian
 *
 * \author Andreas Lintermann
 * \date 22.07.2015
 *
 * This function is called on big endian machines like the Juqueen.
 *
 * \param[in] fileName the name of the file to open
 * \param[in] elemCollector pointer to the collector holding the triangles
 * \param[in] bndCndId the boundary condition id associated with this segment
 * \param[in] segmenId the id of the segment
 * \param[in] the pointer to the segment offset for later computation
 *
 **/
inline void Geometry3D::readSegmentTrianglesBINARY_BE(MString fileName, Collector<element<3>>* elemCollector,
                                                      MInt bndCndId, MInt segmentId, MInt* offset) {
  TRACE();

  using facet_t = struct { float n[3], v1[3], v2[3], v3[3]; };
  facet_t facet;

  // open file in binary format
  FILE* fp = nullptr;

  if((fp = fopen(fileName.c_str(), "rb")) == nullptr) {
    stringstream errorMessage;
    errorMessage << " ERROR in segment::readSegment, couldn't find file : " << fileName << "!";
    mTerm(1, AT_, errorMessage.str());
  }

  element<3>* allelements = elemCollector->a;
  char header[85];
  MUshort ibuff2 = 0;
  if(!fread(header, 1, 84, fp)) {
    TERMM(-1, "ERROR: Memory error!");
  }

  for(MInt i = 0; fread(&facet, 48, 1, fp) > 0; i++) {
    if(!fread(&ibuff2, 2, 1, fp)) {
      TERMM(-1, "ERROR: Memory error!");
    }

    elemCollector->append();

    // reorder
    for(float& k : facet.n) {
      char* buf = (char*)(&k);
      swap4BytesToBE(buf);
    }
    for(float& k : facet.v1) {
      char* buf = (char*)(&k);
      swap4BytesToBE(buf);
    }
    for(float& k : facet.v2) {
      char* buf = (char*)(&k);
      swap4BytesToBE(buf);
    }
    for(float& k : facet.v3) {
      char* buf = (char*)(&k);
      swap4BytesToBE(buf);
    }


    for(MInt j = 0; j < nDim; j++) {
      allelements[*offset].m_normal[j] = facet.n[j];
      allelements[*offset].m_vertices[0][j] = facet.v1[j];
      allelements[*offset].m_vertices[1][j] = facet.v2[j];
      allelements[*offset].m_vertices[2][j] = facet.v3[j];
    }

    allelements[*offset].boundingBox();
    allelements[*offset].m_bndCndId = bndCndId;
    allelements[*offset].m_segmentId = segmentId;
    allelements[*offset].m_originalId = -1;

    if(m_GFieldInitFromSTL)
      for(MInt id = 0; id < m_noLevelSetIntfBndIds; id++)
        if(bndCndId == m_levelSetIntfBndIds[id]) m_noMBElements++;

    *offset = *offset + 1;
  }

  fclose(fp);
}


/** \brief reads triangles from a BINARY file and swaps to little endian
 *
 * \author Andreas Lintermann
 * \date 22.07.2015
 *
 * This function is called on little endian machines.
 *
 * \param[in] fileName the name of the file to open
 * \param[in] elemCollector pointer to the collector holding the triangles
 * \param[in] bndCndId the boundary condition id associated with this segment
 * \param[in] segmenId the id of the segment
 * \param[in] the pointer to the segment offset for later computation
 *
 **/
inline void Geometry3D::readSegmentTrianglesBINARY_LE(MString fileName, Collector<element<3>>* elemCollector,
                                                      MInt bndCndId, MInt segmentId, MInt* offset) {
  TRACE();

  using facet_t = struct { float n[3], v1[3], v2[3], v3[3]; };
  facet_t facet;

  // open file in binary format
  FILE* fp = nullptr;

  if((fp = fopen(fileName.c_str(), "rb")) == nullptr) {
    stringstream errorMessage;
    errorMessage << " ERROR in segment::readSegment, couldn't find file : " << fileName << "!";
    mTerm(1, AT_, errorMessage.str());
  }

  element<3>* allelements = elemCollector->a;
  char header[85];
  MUshort ibuff2 = 0;

  if(fread(header, 1, 84, fp) == 0u) {
    TERMM(-1, "ERROR: Memory error!");
  }

  for(MInt i = 0; fread(&facet, 48, 1, fp) > 0; i++) {
    if(fread(&ibuff2, 2, 1, fp) == 0u) {
      TERMM(-1, "ERROR: Memory error!");
    }

    elemCollector->append();

    for(MInt j = 0; j < nDim; j++) {
      allelements[*offset].m_normal[j] = facet.n[j];
      allelements[*offset].m_vertices[0][j] = facet.v1[j];
      allelements[*offset].m_vertices[1][j] = facet.v2[j];
      allelements[*offset].m_vertices[2][j] = facet.v3[j];
    }

    allelements[*offset].boundingBox();
    allelements[*offset].m_bndCndId = bndCndId;
    allelements[*offset].m_segmentId = segmentId;
    allelements[*offset].m_originalId = -1;

    if(m_GFieldInitFromSTL)
      for(MInt id = 0; id < m_noLevelSetIntfBndIds; id++)
        if(bndCndId == m_levelSetIntfBndIds[id]) m_noMBElements++;

    *offset = *offset + 1;
  }

  fclose(fp);
}


/** \brief reads triangles from a BINARY file
 *
 * \author Andreas Lintermann
 * \date 20.07.2015
 *
 * \param[in] fileName the name of the file to open
 * \param[in] elemCollector pointer to the collector holding the triangles
 * \param[in] bndCndId the boundary condition id associated with this segment
 * \param[in] segmenId the id of the segment
 * \param[in] the pointer to the segment offset for later computation
 *
 **/
inline void Geometry3D::readSegmentTrianglesNETCDF(MString fileName, Collector<element<3>>* elemCollector,
                                                   MInt bndCndId, MInt segmentId, MInt* offset, const MPI_Comm comm) {
  TRACE();

  using namespace maia::parallel_io;

  MInt offsetstart = *offset;
  MInt noTriangles = -1;
  ParallelIo surface(fileName, PIO_READ, comm);
  surface.readScalar(&noTriangles, "noTriangles");

  MInt noEntries = noTriangles;

  MFloatScratchSpace tmp_array(noEntries, AT_, "tmp_array");

  element<3>* allelements = elemCollector->a;
  for(MInt i = 0; i < noTriangles; i++) {
    elemCollector->append();

    allelements[*offset].m_bndCndId = bndCndId;
    allelements[*offset].m_segmentId = segmentId;
    allelements[*offset].m_originalId = -1;

    if(m_GFieldInitFromSTL)
      for(MInt id = 0; id < m_noLevelSetIntfBndIds; id++)
        if(bndCndId == m_levelSetIntfBndIds[id]) m_noMBElements++;

    *offset = *offset + 1;
  }

  MString normalnames[3] = {"normals0", "normals1", "normals2"};
  MString vertexnames[3][3] = {{"vertices00", "vertices01", "vertices02"},
                               {"vertices10", "vertices11", "vertices12"},
                               {"vertices20", "vertices21", "vertices22"}};

  // normals
  for(MInt d = 0; d < nDim; d++) {
    surface.setOffset(noEntries, 0);
    surface.readArray(tmp_array.begin(), normalnames[d]);

    for(MInt i = offsetstart, j = 0; i < offsetstart + noTriangles; i++, j++)
      allelements[i].m_normal[d] = tmp_array[j];
  }

  // vertices
  for(MInt vtx = 0; vtx < nDim; vtx++) {
    for(MInt coord = 0; coord < nDim; coord++) {
      surface.setOffset(noEntries, 0);
      surface.readArray(tmp_array.begin(), vertexnames[vtx][coord]);
      for(MInt i = offsetstart, j = 0; i < offsetstart + noTriangles; i++, j++) {
        allelements[i].m_vertices[vtx][coord] = tmp_array[j];
      }
    }
  }

  for(MInt i = offsetstart; i < offsetstart + noTriangles; i++) {
    allelements[i].boundingBox();
  }
}

/** \brief reads the segments in serial
 *
 * \author Andreas Lintermann
 * \date 17.09.2015
 *
 * Iterates over all segments in the geometry file and first counts the number of
 * elements to be read. Then allocates the according memory and fills the
 * geometry from the infomation in file.
 * Depending on the file type of the geomtry a different method is used for
 * counting an reading.
 *
 **/
inline Collector<element<3>>* Geometry3D::readSegmentsSerial() {
  TRACE();

  NEW_SUB_TIMER(t_countTriangles, "count triangles", m_t_readGeometry);
  g_tc_geometry.emplace_back("count triangles", t_countTriangles);
  RECORD_TIMER_START(t_countTriangles);

  if(domainId() == 0) { // Count segments only on rank 0 and broadcast total count
    if(m_noAllTriangles == 0) {
      MInt segmentId = 0;
      for(m_bodyIt = m_bodyMap.begin(); m_bodyIt != m_bodyMap.end(); m_bodyIt++) {
        // Do not create default body!
        if(m_bodyIt->second->name != "default") {
          for(MInt i = 0; i < m_bodyIt->second->noSegments; i++) {
            stringstream fileName;
            fileName << *geometryContext().getProperty("filename", segmentId)->asString();

            if(geometryContext().noPropertySegments("filename") == 1 && geometryContext().getNoSegments() > 1
               && (m_bodyIt->second->noSegments > 1 || m_bodyMap.size() > 1)) {
              fileName << "." << segmentId;
            }

            m_log << "    - file: " << fileName.str() << endl;

            MInt noElements = 0;
            switch(string2enum(m_geomFileType)) {
              case ASCII:
                countSegmentTrianglesASCII(fileName.str(), &noElements);
                break;
              case BINARY:
                countSegmentTrianglesBINARY(fileName.str(), &noElements);
                break;
              case NETCDF:
                countSegmentTrianglesNETCDF(fileName.str(), &noElements, MPI_COMM_SELF);
                break;
              default:
                break;
            }

            m_log << "      * no. triangles: " << noElements << endl;

            m_noElements += noElements;
            segmentId++;
          }
        }
      }

    } else {
      m_noElements = m_noAllTriangles;
    }
  }
  MPI_Bcast(&m_noElements, 1, MPI_INT, 0, mpiComm(), AT_, "m_noElements");
  RECORD_TIMER_STOP(t_countTriangles);

  m_log << "    - total no. triangles: " << m_noElements << endl << endl;
  m_log << "  + Reading segments ..." << endl;

  NEW_SUB_TIMER(t_readTriangles, "read triangles", m_t_readGeometry);
  g_tc_geometry.emplace_back("read triangles", t_readTriangles);
  RECORD_TIMER_START(t_readTriangles);

  auto* m_allelements = new Collector<element<nDim>>(m_noElements, nDim, 0);
  // Read segments only on rank 0 and communicate (if mode is not disabled)
  if(domainId() == 0 || !m_communicateSegmentsSerial) {
    vector<MInt> allBCs;
    MInt segmentId = 0;
    MInt counter = 0;

    // This loop reads all elements from all segments
    for(m_bodyIt = m_bodyMap.begin(); m_bodyIt != m_bodyMap.end(); m_bodyIt++) {
      // Do not create default body!
      if(m_bodyIt->second->name != "default") {
        for(MInt index = 0; index < m_bodyIt->second->noSegments; index++) {
          stringstream fileName;
          fileName << *geometryContext().getProperty("filename", segmentId)->asString();

          if(geometryContext().noPropertySegments("filename") == 1 && geometryContext().getNoSegments() > 1
             && (m_bodyIt->second->noSegments > 1 || m_bodyMap.size() > 1)) {
            fileName << "." << segmentId;
          }
          MInt bndCndId = *geometryContext().getProperty("BC", segmentId)->asInt();
          allBCs.push_back(bndCndId);

          m_log << "    - file: " << fileName.str() << endl;

          switch(string2enum(m_geomFileType)) {
            case ASCII:
              readSegmentTrianglesASCII(fileName.str(), m_allelements, bndCndId, segmentId, &counter);
              break;
            case BINARY:
              if(is_big_endian() != 0)
                readSegmentTrianglesBINARY_BE(fileName.str(), m_allelements, bndCndId, segmentId, &counter);
              else
                readSegmentTrianglesBINARY_LE(fileName.str(), m_allelements, bndCndId, segmentId, &counter);
              break;
            case NETCDF: {
              const MPI_Comm commSelf = MPI_COMM_SELF;
              const MPI_Comm comm = (m_communicateSegmentsSerial) ? commSelf : mpiComm();
              readSegmentTrianglesNETCDF(fileName.str(), m_allelements, bndCndId, segmentId, &counter, comm);
              break;
            }
            default:
              break;
          }

          segmentId++;
          m_segmentOffsets[segmentId] = counter;
          m_segmentOffsetsWithoutMB[segmentId] = m_segmentOffsetsWithoutMB[segmentId - 1]
                                                 + (m_segmentOffsets[segmentId] - m_segmentOffsets[segmentId - 1]);
          if(m_noLevelSetIntfBndIds > 0) {
            for(MInt i = 0; i < m_noLevelSetIntfBndIds; i++) {
              if(bndCndId == m_levelSetIntfBndIds[i])
                m_segmentOffsetsWithoutMB[segmentId] = m_segmentOffsetsWithoutMB[segmentId - 1];
            }
          }
        }
      }
    }
    m_log << endl;

    m_noAllBCs = allBCs.size();
    mAlloc(m_allBCs, m_noAllBCs, "m_allBCs", 0, AT_);
    for(MInt i = 0; i < m_noAllBCs; i++) {
      m_allBCs[i] = allBCs[i];
    }
  }

  // Broadcast geometry information
  if(m_communicateSegmentsSerial) {
    MFloatScratchSpace normalsBuffer(3 * m_noElements, AT_, "normalsBuffer");
    MFloatScratchSpace verticesBuffer(9 * m_noElements, AT_, "verticesBuffer");
    MIntScratchSpace idsBuffer(3 * m_noElements, AT_, "idsBuffer");
    // Fill buffers on rank 0 and broadcast
    if(domainId() == 0) {
      for(MInt i = 0; i < m_noElements; i++) {
        idsBuffer[i] = m_allelements->a[i].m_bndCndId;
        idsBuffer[m_noElements + i] = m_allelements->a[i].m_segmentId;
        idsBuffer[(2 * m_noElements) + i] = m_allelements->a[i].m_originalId;
        for(MInt j = 0; j < nDim; j++) {
          normalsBuffer[(i * nDim) + j] = m_allelements->a[i].m_normal[j];
          for(MInt k = 0; k < nDim; k++) {
            verticesBuffer[(i * 9) + (j * nDim) + k] = m_allelements->a[i].m_vertices[j][k];
          }
        }
      }
    }
    MPI_Bcast(normalsBuffer.getPointer(), 3 * m_noElements, maia::type_traits<MFloat>::mpiType(), 0, mpiComm(), AT_,
              "elemNormals");
    MPI_Bcast(verticesBuffer.getPointer(), 9 * m_noElements, maia::type_traits<MFloat>::mpiType(), 0, mpiComm(), AT_,
              "elemVertices");
    MPI_Bcast(idsBuffer.getPointer(), 3 * m_noElements, maia::type_traits<MInt>::mpiType(), 0, mpiComm(), AT_,
              "elemIds");

    // Assemble geometry information from buffers
    if(domainId() != 0) {
      for(MInt i = 0; i < m_noElements; i++) {
        m_allelements->append();
        for(MInt j = 0; j < nDim; j++) {
          m_allelements->a[i].m_normal[j] = normalsBuffer[(i * nDim) + j];
          for(MInt k = 0; k < nDim; k++) {
            m_allelements->a[i].m_vertices[j][k] = verticesBuffer[(i * 9) + (j * nDim) + k];
          }
        }
        m_allelements->a[i].boundingBox();
        m_allelements->a[i].m_bndCndId = idsBuffer[i];
        m_allelements->a[i].m_segmentId = idsBuffer[m_noElements + i];
        m_allelements->a[i].m_originalId = idsBuffer[(2 * m_noElements) + i];
      }
    }

    // Communicate remaining info from rank 0
    MPI_Bcast(&m_noMBElements, 1, maia::type_traits<MInt>::mpiType(), 0, mpiComm(), AT_, "m_noMBElements");
    MPI_Bcast(&m_segmentOffsets[0], m_segmentOffsets.size(), maia::type_traits<MInt>::mpiType(), 0, mpiComm(), AT_,
              "m_segmentOffsets");
    MPI_Bcast(&m_segmentOffsetsWithoutMB[0], m_segmentOffsetsWithoutMB.size(), maia::type_traits<MInt>::mpiType(), 0,
              mpiComm(), AT_, "m_segmentOffsetsWithoutMB");
    MPI_Bcast(&m_noAllBCs, 1, maia::type_traits<MInt>::mpiType(), 0, mpiComm(), AT_, "m_noAllBCs");
    if(domainId() != 0) {
      mAlloc(m_allBCs, m_noAllBCs, "m_allBCs", 0, AT_);
    }
    MPI_Bcast(m_allBCs, m_noAllBCs, maia::type_traits<MInt>::mpiType(), 0, mpiComm(), AT_, "m_allBCs");
  }

  RECORD_TIMER_STOP(t_readTriangles);
  return m_allelements;
}

/** \brief reads the segments in parallel
 *
 * \author Andreas Lintermann
 * \date 17.09.2015
 *
 * This function only reads the triangles that really belong to this domain.
 * The algorithm does the folloeing:
 *
 *   1.  Read the workload from the grid file. This is necessary to do the
 *       same decompositioning on the triangls as on the mesh.
 *   2.  Read basic geometry information, i.e., this reads the array containing the
 *       number of triangles contained in the according partitionCell. The size is the same
 *       as the partitionCellList.
 *   2.1 Count all local elements and exchange. This is necessary to determine the offsets
 *       for reading the triangles.
 *   2.2 Read the original triangle ids based on the offsets from 2.1, which are required
 *       to detect double entries. A set is created holding the real number of triangles,
 *       excluding multiple entries. This is the real number of triangles a domain has to
 *       add.
 *   3.  Init the collector holding the triangles base on the number obtained in 2.2
 *   4.  Read the triangles and skip those that have already been read before. Use
 *       a temporary binary array to determine if an entry has already been read. Use
 *       the original triangle id to determine such cases.
 *   4.1 Calculate the segment offsets.
 *   4.2 Get the segment ids and write the to the geometry using the offsets from 4.1
 *   4.3 Get the normals.
 *   4.4 Get the vertices.
 *   5.  Calculate the bounding box for all triangles
 *
 **/
inline Collector<element<3>>* Geometry3D::readSegmentsParallel() {
  TRACE();

  using namespace maia::parallel_io;
  using namespace maia::grid;

  NEW_SUB_TIMER(t_readTriangles, "read triangles", m_t_readGeometry);
  g_tc_geometry.emplace_back("read triangles", t_readTriangles);
  RECORD_TIMER_START(t_readTriangles);

  m_log << "  + Reading triangles ..." << endl;

  // 1. Read the workload from the grid file
  stringstream filename;
  filename << "out/" << m_gridFileName;
  m_log << "    - reading workload from gridfile " << filename.str() << endl;

  MInt partitionCellsCnt = 0;
  MInt noLocalPartitionCells = 0;
  MIntScratchSpace offset(noDomains() + 1, AT_, "offset");

  if(domainId() == 0) {
    ParallelIo grid(filename.str(), PIO_READ, MPI_COMM_SELF);
    partitionCellsCnt = (MInt)grid.getArraySize("partitionCellsGlobalId");

    MFloatScratchSpace partitionCellsWorkLoad(partitionCellsCnt, AT_, "partitionCellsWorkLoad");
    grid.setOffset(partitionCellsCnt, 0);
    grid.readArray(partitionCellsWorkLoad.begin(), "partitionCellsWorkload");

    if(noDomains() > 1) {
      m_log << "    - partitioning workload" << endl;
      optimalPartitioningSerial(&partitionCellsWorkLoad[0], partitionCellsCnt, noDomains(), &offset[0]);

      MPI_Bcast(offset.getPointer(), noDomains(), MPI_INT, 0, mpiComm(), AT_, "offset.getPointer()");
      MPI_Bcast(&partitionCellsCnt, 1, MPI_INT, 0, mpiComm(), AT_, "partitionCellsCnt");
    } else
      offset[0] = 0;
  } else {
    MPI_Bcast(offset.getPointer(), noDomains(), MPI_INT, 0, mpiComm(), AT_, "offset.getPointer()");
    MPI_Bcast(&partitionCellsCnt, 1, MPI_INT, 0, mpiComm(), AT_, "partitionCellsCnt");
  }

  if(domainId() < noDomains() - 1)
    noLocalPartitionCells = offset[domainId() + 1] - offset[domainId()];
  else
    noLocalPartitionCells = partitionCellsCnt - offset[domainId()];

  m_log << "    - no. local partitionCells: " << noLocalPartitionCells << endl;

  // 2. Read basic geometry information
  m_log << "    - reading basic geometry information" << endl;
  ParallelIo geomIO(m_parallelGeomFileName, PIO_READ, mpiComm());

  MInt noTriangles = 0;
  geomIO.readScalar(&noTriangles, "noTriangles");
  geomIO.setOffset(noLocalPartitionCells, offset[domainId()]);

  MIntScratchSpace noTrisPerPartitionCell(noLocalPartitionCells, AT_, "noTrisPerPartitionCell");
  geomIO.readArray(noTrisPerPartitionCell.getPointer(), "noTrisPerPartitionCell");

  // 2.1 Count all local elements and exchange
  MIntScratchSpace noAllLocalTris(noDomains(), AT_, "noAllLocalTris");
  MInt noL = 0;
  for(MInt i = 0; i < noLocalPartitionCells; i++) {
    noL += noTrisPerPartitionCell[i];
  }

  MPI_Allgather(&noL, 1, MPI_INT, noAllLocalTris.getPointer(), 1, MPI_INT, mpiComm(), AT_, "noL",
                "noAllLocalTris.getPointer()");

  MInt originalTriIdOffset = 0;
  // MInt noDomainSubComm = 0;
  for(MInt d = 0; d < domainId(); d++) {
    originalTriIdOffset += noAllLocalTris[d];
    // if(noAllLocalTris[d] > 0)
    // noDomainSubComm++;
  }

  // create subcommunicator
  /*
  m_log << "    - creating subcommunicator" << endl;
  MIntScratchSpace subCommDomainList (noDomainSubComm, AT_, "subCommDomainList");
  for(MInt d = 0, l = 0; d < domainId() ; d++)
    if(noAllLocalTris[d] > 0)
      {
  subCommDomainList[l] = d;
  l++;
      }

  MPI_Group tmp_group, subCommGroup;
  MPI_Comm subComm;

  MPI_Comm_group(mpiComm(), &tmp_group);

  MPI_Group_incl(tmp_group, noDomainSubComm, subCommDomainList.getPointer(), &subCommGroup, AT_ );

  MPI_Comm_create(mpiComm(), subCommGroup, &subComm, AT_, "subComm");

  //no need to read anything
  if(noL == 0)
    {
      m_log << "    - no triangles found for this domain" << endl << endl;
      m_noElements = 0;
      //return new Collector < element<nDim> > (0, nDim, 0);;
    }
  */

  // dummies
  MIntScratchSpace dummyInt(1, AT_, "dummyInt");
  MFloatScratchSpace dummyFloat(1, AT_, "dummyFloat");

  // 2.2 Read the original triangle ids
  MIntScratchSpace originalTriId(noL, AT_, "originalTriId");
  geomIO.setOffset(noL, originalTriIdOffset);
  if(noL == 0)
    geomIO.readArray(dummyInt.getPointer(), "originalTriId");
  else
    geomIO.readArray(originalTriId.getPointer(), "originalTriId");

  MIntScratchSpace segEntry(noL, AT_, "segEntry");
  if(noL == 0)
    geomIO.readArray(dummyInt.getPointer(), "segmentId");
  else
    geomIO.readArray(segEntry.getPointer(), "segmentId");

  vector<MInt> keepOffset;
  for(MInt i = 0; i < noL; i++) {
    pair<set<MInt>::iterator, MBool> ret = m_uniqueOriginalTriId.insert(originalTriId[i]);
    if(ret.second) keepOffset.push_back(i);
  }

  if(noL == 0)
    m_noElements = 0;
  else
    m_noElements = m_uniqueOriginalTriId.size();

  // cout << domainId() << "  size " << m_noElements << endl;
  // 3. Init the collector holding the triangles.
  auto* m_allelements = new Collector<element<nDim>>(m_noElements, nDim, 0);
  element<nDim>* allelements = m_allelements->a;

  for(MInt i = 0; i < m_noElements; i++)
    m_allelements->append();

  // 4.  Read the triangles and skip those that have already been read before.
  m_log << "    - actiually reading triangles" << endl;
  MString normalnames[3] = {"normals0", "normals1", "normals2"};
  MString vertexnames[3][3] = {{"vertices00", "vertices01", "vertices02"},
                               {"vertices10", "vertices11", "vertices12"},
                               {"vertices20", "vertices21", "vertices22"}};

  // 4.1 segmentOffsets
  m_log << "      * segmentOffsets" << endl;
  for(MInt i = 0; i < m_noElements; i++)
    m_segmentOffsets[segEntry[keepOffset[i]] + 1]++;

  MIntScratchSpace segmentCnt(m_noSegments, AT_, "segmentCnt");

  for(MInt i = 0; i < m_noSegments; i++) {
    m_ownSegmentId[i] = false;
    m_segmentOffsets[i + 1] += m_segmentOffsets[i];
    if(m_segmentOffsets[i] != m_segmentOffsets[i + 1]) m_ownSegmentId[i] = true;
    segmentCnt[i] = 0;
  }

  unordered_set<MInt> allBCs;

  // 4.2 segmentIds
  m_log << "      * segmentIds, bndCndIds, originalIds" << endl;
  for(MInt i = 0; i < m_noElements; i++) {
    MInt segmentId = segEntry[keepOffset[i]];
    /*! \page propertyPage1
    \section startTime
    <code>MInt* Geometry::m_allBCs </code>\n
    default = <code>no default</code>\n \n
    Is used to set the Boundary condition for a segment. \n \n
    Possible values are:
    <ul>
    <li> Positive Integer values associated with a boundary condition ID.</li>
    </ul>
    Keywords: <i> GENERAL, BNDRY, GEOMETRY </i>
    */
    MInt bndCndId = *geometryContext().getProperty("BC", segmentId)->asInt();
    MInt off = m_segmentOffsets[segmentId];
    allelements[segmentCnt[segmentId] + off].m_segmentId = segmentId;
    allelements[segmentCnt[segmentId] + off].m_bndCndId = bndCndId;
    allelements[segmentCnt[segmentId] + off].m_originalId = originalTriId[keepOffset[i]];
    allBCs.insert(bndCndId);
    segmentCnt[segmentId]++;
  }

  m_noAllBCs = allBCs.size();
  mAlloc(m_allBCs, m_noAllBCs, "m_allBCs", 0, AT_);
  MInt iter = 0;
  for(auto it = allBCs.begin(); it != allBCs.end(); ++it, iter++)
    m_allBCs[iter] = *it;

  // 4.3 normals
  m_log << "      * normals" << endl;
  for(MInt d = 0; d < nDim; d++) {
    // reset
    for(MInt i = 0; i < m_noSegments; i++)
      segmentCnt[i] = 0;

    MFloatScratchSpace normalEntry(noL, AT_, "normalEntry");
    if(noL == 0)
      geomIO.readArray(dummyFloat.getPointer(), normalnames[d]);
    else
      geomIO.readArray(normalEntry.getPointer(), normalnames[d]);

    for(MInt i = 0; i < m_noElements; i++) {
      MInt segmentId = segEntry[keepOffset[i]];
      MInt off = m_segmentOffsets[segmentId];
      allelements[segmentCnt[segmentId] + off].m_normal[d] = normalEntry[keepOffset[i]];
      segmentCnt[segmentId]++;
    }
  }

  // 4.4. vertices
  m_log << "      * vertices" << endl;
  for(MInt v = 0; v < nDim; v++)
    for(MInt d = 0; d < nDim; d++) {
      for(MInt i = 0; i < m_noSegments; i++)
        segmentCnt[i] = 0;

      MFloatScratchSpace vertexEntry(noL, AT_, "vertexEntry");
      if(noL == 0)
        geomIO.readArray(dummyFloat.getPointer(), vertexnames[v][d]);
      else
        geomIO.readArray(vertexEntry.getPointer(), vertexnames[v][d]);

      for(MInt i = 0; i < m_noElements; i++) {
        MInt segmentId = segEntry[keepOffset[i]];
        MInt off = m_segmentOffsets[segmentId];
        allelements[segmentCnt[segmentId] + off].m_vertices[v][d] = vertexEntry[keepOffset[i]];
        segmentCnt[segmentId]++;
      }
    }

  // 5. Calculate bounding box
  m_log << "      * calculating bounding box for the triangles" << endl;
  for(MInt i = 0; i < m_noElements; i++)
    allelements[i].boundingBox();

  for(MInt i = 0; i < m_noSegments; i++) {
    if(m_segmentOffsets[i + 1] - m_segmentOffsets[i] != 0) {
      m_log << "    - segment id " << i << endl;
      m_log << "      * number of triangles: " << m_segmentOffsets[i + 1] - m_segmentOffsets[i] << endl;
      m_log << "      * offsets:            " << m_segmentOffsets[i] << " - " << m_segmentOffsets[i + 1] << endl;
    }
  }
  m_log << endl;

  RECORD_TIMER_STOP(t_readTriangles);

  return m_allelements;
}

/** \brief reads the STL segments from file
 *
 * \author Rainhill Freitas, Andreas Lintermann
 * \date 20.07.2015
 *
 * Calls the according function to read the geometry.
 *
 **/
void Geometry3D::readSegments() {
  TRACE();

  NEW_SUB_TIMER(t_readGeometry, "read geometry", m_t_geometryAll);
  m_t_readGeometry = t_readGeometry;
  g_tc_geometry.emplace_back("read geometry", m_t_readGeometry);
  RECORD_TIMER_START(t_readGeometry);

  otherCalls++;

  // check binary or ascii
  m_geomFileType = "ASCII";
  if(geometryContext().propertyExists("geomFileType", 0))
    m_geomFileType = *(geometryContext().getProperty("geomFileType", 0)->asString(0));

  // check if we need to count the triangles
  m_noAllTriangles = 0;
  if(geometryContext().propertyExists("noAllTriangles", 0))
    m_noAllTriangles = *(geometryContext().getProperty("noAllTriangles", 0)->asInt(0));

  string tmp;

  // This loop only counts all elements of all segments
  for(MInt i = 0; i <= m_noSegments; i++) {
    m_segmentOffsets[i] = 0;
    m_segmentOffsetsWithoutMB[i] = 0;
  }

  mAlloc(m_ownSegmentId, m_noSegments, "m_ownSegmentId", true, AT_);

  m_log << "  + Geometry configuration:" << endl;
  m_log << "    - File type is set to: " << (m_parallelGeometry ? "NETCDF" : m_geomFileType) << endl;
  m_log << "    - endianess:           " << (is_big_endian() ? "big endian" : "little endian") << endl;
  m_log << "    - reading method:      " << (m_parallelGeometry ? "parallel" : "serial") << endl << endl;

  // trust me, this is just the definition, init will be performed in the following funtions
  Collector<element<3>>* m_allelements = nullptr;

  if(m_flowSolver && m_parallelGeometry) {
    m_allelements = readSegmentsParallel();
  } else {
    m_allelements = readSegmentsSerial();
  }

  if(m_GFieldInitFromSTL) {
    m_mbelements = new Collector<element<nDim>>(m_noMBElements, nDim, 0);
    m_noElements -= m_noMBElements;
  }

  MInt initSize = m_noElements;

  // add 15% overhead
  /*
  if(m_parallelGeometry)
    initSize *= m_parGeomMemFactor;
  */

  m_elements = new Collector<element<nDim>>(initSize, nDim, 0);
  elements = m_elements->a;
  if(m_GFieldInitFromSTL) mbelements = m_mbelements->a;

  MInt mbelem_counter = 0, elem_counter = 0;
  MBool mbElem = false;
  element<3>* allelements = m_allelements->a;

  for(MInt allelem_counter = 0; allelem_counter < m_allelements->size(); allelem_counter++) {
    if(m_GFieldInitFromSTL && m_noLevelSetIntfBndIds > 0) {
      for(MInt id = 0; id < m_noLevelSetIntfBndIds; id++) {
        if(allelements[allelem_counter].m_bndCndId == m_levelSetIntfBndIds[id]) {
          mbElem = true;
          break;
        }
      }

      // LevelSet Moving boundary elements
      if(mbElem) {
        m_mbelements->append();

        for(MInt j = 0; j < 3; j++) {
          mbelements[mbelem_counter].m_normal[j] = allelements[allelem_counter].m_normal[j];
          for(MInt i = 0; i < 3; i++) {
            mbelements[mbelem_counter].m_vertices[j][i] = allelements[allelem_counter].m_vertices[j][i];
          }
        }

        mbelements[mbelem_counter].boundingBox();
        mbelements[mbelem_counter].m_bndCndId = allelements[allelem_counter].m_bndCndId;
        mbelements[mbelem_counter].m_segmentId = allelements[allelem_counter].m_segmentId;
        mbelements[mbelem_counter].m_originalId = allelements[allelem_counter].m_originalId;
        mbelem_counter++;
        mbElem = false;
      }

      // Non movable elements
      else {
        m_elements->append();

        for(MInt j = 0; j < 3; j++) {
          elements[elem_counter].m_normal[j] = allelements[allelem_counter].m_normal[j];
          for(MInt i = 0; i < 3; i++) {
            elements[elem_counter].m_vertices[j][i] = allelements[allelem_counter].m_vertices[j][i];
          }
        }
        elements[elem_counter].boundingBox();
        elements[elem_counter].m_bndCndId = allelements[allelem_counter].m_bndCndId;
        elements[elem_counter].m_segmentId = allelements[allelem_counter].m_segmentId;
        elements[elem_counter].m_originalId = allelements[allelem_counter].m_originalId;
        elem_counter++;
      }
    } else {
      // LevelSet Moving boundary elements
      if(m_GFieldInitFromSTL && (m_levelSetIntfBndId == allelements[allelem_counter].m_bndCndId)) {
        m_mbelements->append();

        for(MInt j = 0; j < 3; j++) {
          mbelements[mbelem_counter].m_normal[j] = allelements[allelem_counter].m_normal[j];
          for(MInt i = 0; i < 3; i++) {
            mbelements[mbelem_counter].m_vertices[j][i] = allelements[allelem_counter].m_vertices[j][i];
          }
        }
        mbelements[mbelem_counter].boundingBox();
        mbelements[mbelem_counter].m_bndCndId = allelements[allelem_counter].m_bndCndId;
        mbelements[mbelem_counter].m_segmentId = allelements[allelem_counter].m_segmentId;
        mbelements[mbelem_counter].m_originalId = allelements[allelem_counter].m_originalId;
        mbelem_counter++;
      } else { // Non movable elements

        m_elements->append();

        for(MInt j = 0; j < 3; j++) {
          elements[elem_counter].m_normal[j] = allelements[allelem_counter].m_normal[j];
          for(MInt i = 0; i < 3; i++)
            elements[elem_counter].m_vertices[j][i] = allelements[allelem_counter].m_vertices[j][i];
        }
        elements[elem_counter].boundingBox();
        elements[elem_counter].m_bndCndId = allelements[allelem_counter].m_bndCndId;
        elements[elem_counter].m_segmentId = allelements[allelem_counter].m_segmentId;
        elements[elem_counter].m_originalId = allelements[allelem_counter].m_originalId;

        elem_counter++;
      }
    }
  }

  delete m_allelements;

  RECORD_TIMER_STOP(t_readGeometry);

  NEW_SUB_TIMER(t_calcBndBox, "calculating bounding box", m_t_geometryAll);
  g_tc_geometry.emplace_back("calculating bounding box", t_calcBndBox);
  RECORD_TIMER_START(t_calcBndBox);

  calculateBoundingBox();
  RECORD_TIMER_STOP(t_calcBndBox);

  // update m_haloElementOffset
  setHaloElementOffset(m_noElements);

  if(m_flowSolver && m_debugParGeom && m_noElements > 0) {
    writeParallelGeometryVTK("readSeg");
  }
}

/** \brief Write the current geometry to a VTK file for parallel computation
 *
 * \author Andreas Lintermann
 * \date 25.09.2015
 *
 * \param[in] filename the name of the file without extension
 **/
void Geometry3D::writeParallelGeometryVTK(MString filename) {
  TRACE();

  const MBool allTimerRunning = timers().isRunning(m_t_geometryAll);
  if(!allTimerRunning) {
    RECORD_TIMER_START(m_t_geometryAll);
  }
  NEW_SUB_TIMER(t_wrtParGeom, "writing parallel geometry", m_t_geometryAll);
  g_tc_geometry.emplace_back("writing parallel geometry", t_wrtParGeom);
  RECORD_TIMER_START(t_wrtParGeom);

  stringstream name;
  name << "out/" << domainId() << "_" << filename << ".vtk";
  ofstream st;
  st.open(name.str());
  st << "# vtk DataFile Version 3.0" << endl;
  st << "vtk output" << endl;
  st << "ASCII" << endl;
  st << "DATASET POLYDATA" << endl;
  st << "POINTS " << m_noElements * nDim << " float" << endl;

  for(MInt i = 0; i < m_noElements; i++)
    for(MInt v = 0; v < 3; v++) {
      for(MInt d = 0; d < 3; d++)
        st << elements[i].m_vertices[v][d] << " ";
      st << endl;
    }

  st << endl;
  st << "POLYGONS " << m_noElements << " " << m_noElements * (nDim + 1) << endl;
  for(MInt i = 0, j = 0; i < m_noElements; i++) {
    st << nDim << " ";
    for(MInt v = 0; v < 3; v++, j++)
      st << j << " ";
    st << endl;
  }

  st << endl;
  st << "CELL_DATA " << m_noElements << endl;
  st << "SCALARS domainId int 1" << endl;
  st << "LOOKUP_TABLE default" << endl;
  for(MInt i = 0; i < m_noElements; i++)
    st << domainId() << endl;
  st << endl;
  st << "SCALARS segmentId int 1" << endl;
  st << "LOOKUP_TABLE default" << endl;
  for(MInt i = 0; i < m_noElements; i++)
    st << elements[i].m_segmentId << endl;
  st << endl;
  st << "SCALARS originalCellId int 1" << endl;
  st << "LOOKUP_TABLE default" << endl;
  for(MInt i = 0; i < m_noElements; i++)
    st << elements[i].m_originalId << endl;
  st << endl;
  st << "SCALARS BC int 1" << endl;
  st << "LOOKUP_TABLE default" << endl;
  for(MInt i = 0; i < m_noElements; i++)
    st << elements[i].m_bndCndId << endl;
  st << "SCALARS area double 1" << endl;
  st << "LOOKUP_TABLE default" << endl;
  for(MInt i = 0; i < m_noElements; i++) {
    MFloat edge1[3];
    MFloat edge2[3];
    MFloat cross[3];
    for(MInt j = 0; j < nDim; j++) {
      edge1[j] = elements[i].m_vertices[1][j] - elements[i].m_vertices[0][j];
      edge2[j] = elements[i].m_vertices[2][j] - elements[i].m_vertices[0][j];
    }
    cross[0] = edge1[1] * edge2[2] - edge2[1] * edge1[2];
    cross[1] = edge1[2] * edge2[0] - edge2[2] * edge1[0];
    cross[2] = edge1[0] * edge2[1] - edge2[0] * edge1[1];

    MFloat area = sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
    st << area << endl;
  }
  st.close();
  RECORD_TIMER_STOP(t_wrtParGeom);
  if(!allTimerRunning) {
    RECORD_TIMER_STOP(m_t_geometryAll);
  }
}

/** \brief Calculates the global bounding box of the geometry
 *
 * \author Andreas Lintermann
 * \date 25.09.2015
 *
 * Works for both, parallel and serial geometry.
 *
 **/
void Geometry3D::calculateBoundingBox() {
  TRACE();

  m_log << "  + Calculating bounding box ..." << endl;
  m_log << "    - number of elements to consider: " << m_noElements << endl;

  // Communicate a reference triangle for those that have no triangles available
  MIntScratchSpace noGlobalElem(noDomains(), AT_, "noGlobalElem");
  noGlobalElem[domainId()] = m_noElements;

  MPI_Allgather(&m_noElements, 1, MPI_INT, noGlobalElem.getPointer(), 1, MPI_INT, mpiComm(), AT_, "m_noElements",
                "noGlobalElem.getPointer()");
  MInt refDomain = 0;
  for(MInt d = 0; d < noDomains(); d++) {
    if(noGlobalElem[d] > 0) {
      refDomain = d;
      break;
    }
  }

  std::array<MFloat, 2 * nDim> triMinMax;
  triMinMax.fill(std::nanf(""));

  if(domainId() == refDomain) {
    for(MInt j = 0; j < 2 * nDim; j++) {
      triMinMax[j] = elements[0].m_minMax[j];
    }
  }

  MPI_Bcast(triMinMax.data(), 2 * nDim, MPI_DOUBLE, refDomain, mpiComm(), AT_, "triMinMax.getPointer()");

  // Calculate bounding box of segment
  for(MInt j = 0; j < nDim; j++) {
    m_minMax[j] = triMinMax[j];
    m_minMax[j + nDim] = triMinMax[j + nDim];
  }

  for(MInt i = 0; i < m_noElements; i++) {
    for(MInt j = 0; j < nDim; j++) {
      m_minMax[j + nDim] =
          (m_minMax[j + nDim] < elements[i].m_minMax[j + nDim]) ? elements[i].m_minMax[j + nDim] : m_minMax[j + nDim];
      m_minMax[j] = (m_minMax[j] > elements[i].m_minMax[j]) ? elements[i].m_minMax[j] : m_minMax[j];
    }
  }
  if(m_noElements == 0) ASSERT(m_GFieldInitFromSTL, "");

  if(m_noElements > 0) {
    m_log << "    - local max: ";
    for(MInt i = 0; i < nDim; i++) {
      m_log << m_minMax[i + nDim] << " ";
    }
    m_log << endl;

    m_log << "    - local min: ";
    for(MInt i = 0; i < nDim; i++) {
      m_log << m_minMax[i] << " ";
    }
    m_log << endl;
  }

  if(m_parallelGeometry) {
    MPI_Allreduce(MPI_IN_PLACE, &m_minMax[0], nDim, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE", "m_minMax[0]");
    MPI_Allreduce(MPI_IN_PLACE, &m_minMax[0] + nDim, nDim, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                  "m_minMax[0]+nDim");
  }

  if(m_noElements > 0) {
    m_log << "    - max: ";
    for(MInt i = 0; i < nDim; i++) {
      m_log << m_minMax[i + nDim] << " ";
    }
    m_log << endl;

    m_log << "    - min: ";
    for(MInt i = 0; i < nDim; i++) {
      m_log << m_minMax[i] << " ";
    }
    m_log << endl << endl;
  }

  // Calculate bounding box of Moving Boundary segment
  if(m_GFieldInitFromSTL) {
    UpdateMBBoundingBox();
    m_log << " Bounding box for Moving Boundary segment:" << endl;
    m_log << " max = ";
    for(MInt i = 0; i < nDim; i++)
      m_log << m_mbminMax[i + nDim] << " ";

    m_log << endl;
    m_log << " min = ";
    for(MInt i = 0; i < nDim; i++)
      m_log << m_mbminMax[i] << " ";

    m_log << endl;

    if(m_forceBoundingBox) {
      // use the STL2levelset geometry as bounding box to define the length0 and center of gravity
      // this can be useful, to force a length onto a setup!
      for(MInt i = 0; i < nDim; i++) {
        if(m_mbminMax[i] < m_minMax[i]) {
          m_minMax[i] = m_mbminMax[i];
        }
        if(m_mbminMax[i + nDim] > m_minMax[i + nDim]) {
          m_minMax[i + nDim] = m_mbminMax[i + nDim];
        }
      }
    }
  }
}


void Geometry3D::correctVertexCoordinates() {
  TRACE();

  otherCalls++;


  MInt noEl = m_elements->size();
  MFloat epsilon = 0.0000001;
  MFloat maxCoordinate = F0;
  //---

  // find the maximum coordinate
  for(MInt el = 0; el < noEl; el++)
    for(MInt i = 0; i < 3; i++)
      for(MInt j = 0; j < 3; j++)
        maxCoordinate = mMax(ABS(elements[el].m_vertices[i][j]), maxCoordinate);

  // compute epsilon
  epsilon *= maxCoordinate;

  // remove small values around zero
  for(MInt el = 0; el < noEl; el++) {
    for(MInt i = 0; i < 3; i++) {
      for(MInt j = 0; j < 3; j++) {
        if(ABS(elements[el].m_vertices[j][i]) < epsilon) elements[el].m_vertices[j][i] = F0;
      }
    }
  }
}


/** \brief Determines all elements that are inside or intersect the target region with the separating axis theorem
   (SAT). \author Andreas Lintermann \date 24.11.2009

    The separating axis theorem (SAT) states, that two convex polyhedra, that are to be tested against intersection are
    disjoint, if they can be separated along either an axis parallel to a normal of a face of the first or the second
   polyhedra, or along an axis formed from the cross product of an edge from them.

    To make things simpler, an axis aligned gridcell is considered centered in the origin. Therefore, the triangle is
   translated accordingly.

    <b>Algorithm:</b>
    <ul>
    <li><b>1)</b> Create for all grid directions the normal with the triangle edges
    Then project the triangle vertices \f$v0,v1,v2\f$ onto these edges and the box as well.
    For latter purpose create a radius defined by:

    \f[r = h_x|a_x| + h_y|a_y| + h_z|a_z|,\f]

    where \f$h_i\f$ is the halflength of the box and \f$a_i\f$ is the grid direction.
    One of the directions for \f$a_i\f$ is 0, which reduces the above to:

    \f{eqnarray*}{
    r_x = h_y|a_y| + h_z|a_z|\\
    r_y = h_x|a_x| + h_z|a_z|\			\
    r_z = h_x|a_x| + h_y|a_y|
    \f}

    Then check the projected triangle points against the radius
    </li>

    <li><b>2)</b> Test overlap in the \f${x,y,z}\f$-directions
    find min, max of the triangle each direction, and test for overlap in
    that direction -- this is equivalent to testing a minimal AABB around
    the triangle against the AABB
     </li>
    <li><b>3)</b> Test if the box intersects the plane of the triangle -> compute plane equation of triangle:
   \f$n*x+d=0\f$
     </li>
    </ul>
    If all tests pass, the the box and the triangle are intersecting!

    This function replaces Geometry3D::getIntersectionElements.

*/
MInt Geometry3D::getIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList, MFloat cellHalfLength,
                                         const MFloat* const cell_coords) {
  // TRACE();

  getIECallCounter++;

  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  MFloat vert_mov[MAX_SPACE_DIMENSIONS][MAX_SPACE_DIMENSIONS];
  MFloat edge_mov[MAX_SPACE_DIMENSIONS][MAX_SPACE_DIMENSIONS];
  MFloat edge_fabs[MAX_SPACE_DIMENSIONS][MAX_SPACE_DIMENSIONS];
  MFloat min;
  MFloat max;
  MFloat d;
  MFloat p0;
  MFloat p1;
  MFloat rad;
  MFloat normal[MAX_SPACE_DIMENSIONS];
  MFloat vmin[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                       numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
  MFloat vmax[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                       numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized

  MInt combo[3][6] = {{0, 2, 0, 2, 1, 2}, {0, 2, 0, 2, 0, 1}, {0, 1, 0, 1, 1, 2}};

  MInt noReallyIntersectingNodes = 0;
  MBool cut;

  for(MInt n = 0; n < noNodes; n++) {
    cut = true;

    // Fill arrays (translate triangle)
    for(MInt i = 0; i < nDim; i++) {
      for(MInt j = 0; j < nDim; j++) {
        vert_mov[j][i] = elements[nodeList[n]].m_vertices[j][i] - cell_coords[i];
      }
      for(MInt j = 0; j < nDim; j++) {
        edge_mov[j][i] = vert_mov[(j + 1) % nDim][i] - vert_mov[j][i];
        edge_fabs[j][i] = fabs(edge_mov[j][i]);
      }
    }

    MInt mul;

    // Step 1) of the algorithm
    for(MInt i = 0; i < nDim; i++) {
      mul = 1;
      for(MInt l = 0, j = 2, k = 1; l < nDim; l++) {
        p0 = mul * edge_mov[i][j] * vert_mov[combo[i][2 * l]][k] - mul * edge_mov[i][k] * vert_mov[combo[i][2 * l]][j];
        p1 = mul * edge_mov[i][j] * vert_mov[combo[i][2 * l + 1]][k]
             - mul * edge_mov[i][k] * vert_mov[combo[i][2 * l + 1]][j];

        (p0 < p1) ? (min = p0, max = p1) : (min = p1, max = p0);
        rad = cellHalfLength * (edge_fabs[i][j] + edge_fabs[i][k]);

        if(min > rad || max < -rad) {
          cut = false;
          break;
        }

        mul *= -1;
        j -= l;
        k = 0;
      }
      if(!cut) break;
    }
    if(!cut) continue;

    // Step 2) of the algorithm
    for(MInt i = 0; i < nDim; i++) {
      min = max = vert_mov[0][i];
      for(MInt j = 1; j < nDim; j++) {
        if(vert_mov[j][i] < min) min = vert_mov[j][i];
        if(vert_mov[j][i] > max) max = vert_mov[j][i];
      }

      if(min > cellHalfLength || max < -cellHalfLength) {
        cut = false;
        break;
      }
    }
    if(!cut) continue;

    // Step 3) of the algorithm
    normal[0] = edge_mov[0][1] * edge_mov[1][2] - edge_mov[0][2] * edge_mov[1][1];
    normal[1] = edge_mov[0][2] * edge_mov[1][0] - edge_mov[0][0] * edge_mov[1][2];
    normal[2] = edge_mov[0][0] * edge_mov[1][1] - edge_mov[0][1] * edge_mov[1][0];

    d = 0.0;
    for(MInt i = 0; i < nDim; i++)
      d += normal[i] * vert_mov[0][i];
    d *= -1;

    for(MInt i = 0; i < nDim; i++) {
      if(normal[i] > 0.0f) {
        vmin[i] = -cellHalfLength;
        vmax[i] = cellHalfLength;
      } else {
        vmin[i] = cellHalfLength;
        vmax[i] = -cellHalfLength;
      }
    }

    if((normal[0] * vmin[0] + normal[1] * vmin[1] + normal[2] * vmin[2]) + d > 0.0f) continue;
    if((normal[0] * vmax[0] + normal[1] * vmax[1] + normal[2] * vmax[2]) + d >= 0.0f) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
    }
  }

  nodeList.resize(noReallyIntersectingNodes);

  getIECommCounter += noReallyIntersectingNodes;

  return noReallyIntersectingNodes;
}

/** \brief Determines all elements that are inside or intersect the target region
    \author Rainhill Freitas
    \date unknown

    Target region in this case means the gridcell.
    \bug labels:GEOM This code does not run on IBM_BLUE_GENE. For plane/line intersection Cramers rule is used. The
   determinant in the denominator of Cramer's rule can become 0 and this is not acceptable for IBM_BLUE_GENE. Use the
   other Geometry3D::getIntersectionElements instead.
  */
MInt Geometry3D::getIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  //  TRACE();

  getIECallCounter++;

  MInt noReallyIntersectingNodes = 0;
  // get all candidates for an intersection
  // Check for intersection...
  bitset<6> points[3];
  bitset<6> faceCodes[6];
  bitset<6> pCode;
  MInt spaceId;
  MInt spaceId1;
  MInt spaceId2;
  MInt edges[3][2] = {{0, 1}, {1, 2}, {0, 2}};
  // Corner points of the target regione i.e. p0=(targetRegion[targetPoints[0][0],
  //                                              targetRegion[targetPoints[0][1],
  //                                              targetRegion[targetPoints[0][2])
  MInt targetPoints[8][3] = {{0, 1, 2}, {3, 1, 2}, {0, 4, 2}, {3, 4, 2}, {0, 1, 5}, {3, 1, 5}, {0, 4, 5}, {3, 4, 5}};

  // Edges of the targetRegion (using targetPoints)
  MInt targetEdges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7},  // x-edges (i.e. parallel to x-axis)
                             {0, 2}, {1, 3}, {4, 6}, {5, 7},  // y-edes
                             {0, 4}, {1, 5}, {2, 6}, {3, 7}}; // z-edges
  MFloat e1[3];
  MFloat e2[3];
  MInt rejection = 0;
  MBool piercePointInside = 0;
  MBool triviallyAccepted = 0;

  // Each of the following arrays holds one different point for
  // all of the 6 planes. Points are built with the targetRegion
  // i.e. a = targetRegiont[pointAInPlane[0]],
  //          targetRegiont[pointAInPlane[1]],
  //          targetRegiont[pointAInPlane[2]])
  // etc.
  MInt pointAInPlane[6][3] = {{0, 1, 2}, {3, 4, 5}, {0, 1, 2}, {3, 4, 5}, {0, 1, 2}, {3, 4, 5}};

  MInt pointBInPlane[6][3] = {{0, 4, 2}, {3, 1, 5}, {3, 1, 2}, {0, 4, 5}, {3, 1, 2}, {3, 1, 5}};

  MInt pointCInPlane[6][3] = {{0, 1, 5}, {3, 4, 2}, {3, 1, 5}, {0, 4, 2}, {3, 4, 2}, {0, 1, 5}};
  MFloat a[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat b[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat c[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat d[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // start of piercing edge
  MFloat e[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // end of piercing edge
  MFloat s;
  MFloat gamma; // For pierce point calculation
  MFloat p2;
  MFloat q;
  MFloat epsilon = 0.00000000001;
  MFloat pP[3]; // piercePoint
  faceCodes[0] = IPOW2(0);
  faceCodes[1] = IPOW2(1);
  faceCodes[2] = IPOW2(2);
  faceCodes[3] = IPOW2(3);
  faceCodes[4] = IPOW2(4);
  faceCodes[5] = IPOW2(5);
  //   edges[0][0] = 0;
  //   edges[0][1] = 1;
  //   edges[1][0] = 1;
  //   edges[1][1] = 2;
  //   edges[2][0] = 0;
  //   edges[2][1] = 2;
  bitset<6> result;

  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  //  return noNodes;
  noReallyIntersectingNodes = 0;

  for(MInt n = 0; n < noNodes; n++) {
    // create edges (in 3D) AB, AC, BC and point A B C
    // Determine outcode (see Aftosmis, Solution Adaptive Cartesian Grid Methods for Aerodynamic flows ...)

    // Loop over all points of an element<3>
    for(MInt p = 0; p < nDim; p++) {
      points[p] = 0;
      // Calculate outcode for point
      for(MInt j = 0; j < nDim; j++) {
        if(elements[nodeList[n]].m_vertices[p][j] < targetRegion[j]) {
          points[p] |= faceCodes[2 * j];
        }
        if(elements[nodeList[n]].m_vertices[p][j] > targetRegion[j + nDim]) {
          points[p] |= faceCodes[2 * j + 1];
        }
      }
      //      m_log << points[p] << endl;
    }
    rejection = 0;
    // check outcode combinations for edges for trivial rejection
    for(MInt i = 0; i < nDim; i++) {
      if((points[edges[i][0]] & points[edges[i][1]]) != 0) {
        rejection++;
      } else {
        // If one point is inside region the element<3> is trivially accepted
        triviallyAccepted = false;
        for(MInt k = 0; k < nDim; k++) {
          if(points[k] == 0) {
            triviallyAccepted = true;
            break;
          }
        }
        if(triviallyAccepted) {
          break;
        }
        // No trivial rejection, check for rejection of subsegment:
        // For all pierce points!
        // 1. Calculate pierce point:
        //    a - determine plane for pierce point calculation
        //    b - calculate pierce point
        // 2. Check for rejection of new segment
        //    a - calculate new outcode
        //    b - check for containment in pierce planes face
        // 3. If all(!) pierce points are rejected -> reject edge
        // TODO labels:GEOM This algorith might get a problem if a triangle lies completely on
        //    a face.

        // 1.a
        result = (points[edges[i][0]] | points[edges[i][1]]);
        piercePointInside = false;
        for(MInt j = 0; j < 2 * nDim; j++) {
          if(result[j] == 1) {
            // pierce plane found
            // Algorithm taken von AFTOSMIS.  There seems to be an
            // error in the calculation of s in AFTOSMIS formula, the one below
            // should be correct!
            for(MInt k = 0; k < nDim; k++) {
              a[k] = targetRegion[pointAInPlane[j][k]];
              b[k] = targetRegion[pointBInPlane[j][k]];
              c[k] = targetRegion[pointCInPlane[j][k]];
              d[k] = elements[nodeList[n]].m_vertices[edges[i][0]][k];
              e[k] = elements[nodeList[n]].m_vertices[edges[i][1]][k];
            }
            gamma = ((e[0] - d[0]) * ((a[2] - c[2]) * (c[1] - b[1]) - (a[1] - c[1]) * (c[2] - b[2]))
                     - (e[1] - d[1]) * ((a[2] - c[2]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[2] - b[2]))
                     + (e[2] - d[2]) * ((a[1] - c[1]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[1] - b[1])));

            s = ((c[0] - d[0]) * ((a[2] - c[2]) * (c[1] - b[1]) - (a[1] - c[1]) * (c[2] - b[2]))
                 - (c[1] - d[1]) * ((a[2] - c[2]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[2] - b[2]))
                 + (c[2] - d[2]) * ((a[1] - c[1]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[1] - b[1])))
                / gamma;

            // 1. b Pierce point pP in plane j:
            for(MInt k = 0; k < nDim; k++) {
              pP[k] = d[k] + s * (e[k] - d[k]);
            }
            pCode = 0;
            // 2. a Calculate outcode for pierce point
            for(MInt k = 0; k < nDim; k++) {
              if(pP[k] < targetRegion[k]) {
                pCode |= faceCodes[2 * k];
              }
              if(pP[k] > targetRegion[k + nDim]) {
                pCode |= faceCodes[2 * k + 1];
              }
            }
            //	    m_log << " outcode of pierce point :" << pCode << endl;

          } else {
            continue;
          }
          // 2. b
          result = faceCodes[j];
          result.flip();
          result = (result & pCode);
          if(result == 0) { // -> is contained
            piercePointInside = true;
          }
        }
        // reject if all pierce points are off coresponding face
        if(!piercePointInside) {
          rejection++;
        }
      }
    }
    // Check for internal element<3>, i.e. completely inside target
    result = (points[0] | points[1] | points[2]);
    if(result == 0) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
      continue;
    }
    // If not all edges are rejected a cutting element<3> has been found
    //    m_log << rejection << endl;
    if(rejection < 3) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
      continue;
    }

    // Even if trivially rejected, check if targetRegion is completely inside triangle
    // For that check all edges of the cell against the plane of the triangle,
    // calculate all pierce points and their outcodes and finally check for containment
    // in the targetRegion. Before the pierce point calculation check if the current
    // edge is parallel to the triangle (evaluate the distance of the edgepoints to the plane)
    if(rejection >= 3) {
      for(MInt t = 0; t < 12; t++) { // Loop over all 12 edges of the targetRegion (cell-cube)
        for(MInt p = 0; p < nDim; p++) {
          e1[p] = targetRegion[targetPoints[targetEdges[t][0]][p]];
          e2[p] = targetRegion[targetPoints[targetEdges[t][1]][p]];
        }

        if(t < 4) {
          spaceId = 1;
          spaceId1 = 2;
          spaceId2 = 0;
        } else {
          if(t > 3 && t < 8) {
            spaceId = 2;
            spaceId1 = 0;
            spaceId2 = 1;
          } else {
            spaceId = 0;
            spaceId1 = 1;
            spaceId2 = 2;
          }
        }

        for(MInt k = 0; k < nDim; k++) {
          a[k] = elements[nodeList[n]].m_vertices[0][k];
          b[k] = elements[nodeList[n]].m_vertices[1][k];
          c[k] = elements[nodeList[n]].m_vertices[2][k];
        }
        if(approx(a[spaceId1], b[spaceId1], MFloatEps) && !approx(a[spaceId1], c[spaceId1], MFloatEps)) {
          // aspaceId != bspaceId, otherwise a and b would be the same point
          q = (e1[spaceId1] - a[spaceId1]) / (c[spaceId1] - a[spaceId1]);
          p2 = (e1[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
        } else {
          if(!approx(a[spaceId1], b[spaceId1], MFloatEps) && approx(a[spaceId1], c[spaceId1], MFloatEps)) {
            // aspaceId != cspaceId, otherwise a and c would be the same point
            p2 = (e1[spaceId1] - a[spaceId1]) / (b[spaceId1] - a[spaceId1]);
            q = (e1[spaceId] - a[spaceId] - p2 * (b[spaceId] - a[spaceId])) / (c[spaceId] - a[spaceId]);
          } else {
            if(approx(a[spaceId], b[spaceId], MFloatEps) && !approx(a[spaceId], c[spaceId], MFloatEps)) {
              // aspaceId1 != bspaceId1, otherwise a and b would be the same point
              q = (e1[spaceId] - a[spaceId]) / (c[spaceId] - a[spaceId]);
              p2 = (e1[spaceId1] - a[spaceId1] - q * (c[spaceId1] - a[spaceId1])) / (b[spaceId1] - a[spaceId1]);
            } else {
              if(!approx(a[spaceId], b[spaceId], MFloatEps) && approx(a[spaceId], c[spaceId], MFloatEps)) {
                // aspaceId1 != cspaceId1, otherwise a and c would be the same point
                p2 = (e1[spaceId] - a[spaceId]) / (b[spaceId] - a[spaceId]);
                q = (e1[spaceId1] - a[spaceId1] - p2 * (b[spaceId1] - a[spaceId1])) / (c[spaceId1] - a[spaceId1]);
              } else {
                // aspaceId1 != bspaceId1 && aspaceId1 != cspaceId1 && aspaceId != bspaceId && aspaceId != cspaceId
                q = ((e1[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                     - (b[spaceId1] - a[spaceId1]) * (e1[spaceId] - a[spaceId]))
                    / ((c[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                       - (b[spaceId1] - a[spaceId1]) * (c[spaceId] - a[spaceId]));
                p2 = (e1[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
              }
            }
          }
        }

        if(p2 * q >= 0 || p2 * q < 0) {
          // compute s
          gamma = a[spaceId2] + p2 * (b[spaceId2] - a[spaceId2]) + q * (c[spaceId2] - a[spaceId2]);
          s = (gamma - e1[spaceId2]) / (e2[spaceId2] - e1[spaceId2]);

          if(s < -epsilon || s > F1 + epsilon || p2 < -epsilon || q < -epsilon || (p2 + q) > F1 + epsilon) {
          } else {
            ////
            /*	if( edgeTriangleIntersection(elements[nodeList[n]].m_vertices[0],
                   elements[nodeList[n]].m_vertices[1],
                   elements[nodeList[n]].m_vertices[2], e1, e2) ){*/
            nodeList[noReallyIntersectingNodes] = nodeList[n];
            noReallyIntersectingNodes++;
            break;
          }
        }
      }
    }

    // write nodes in reallyIntersectingNodes
  }
  //  if(noNodes)
  if((noNodes != 0) && (noReallyIntersectingNodes == 0)) {
    //    m_log << "Rejected ! "<< ++rejectionCounter << endl;
  }

  nodeList.resize(noReallyIntersectingNodes);

  getIECommCounter += noReallyIntersectingNodes;

  return noReallyIntersectingNodes;
}


/** \brief Determines all elements that are inside or intersect the target region
 *
 *
 */
MInt Geometry3D::getIntersectionElementsTetraeder(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  TRACE();

  getIETCallCounter++;

  MInt noReallyIntersectingNodes = 0;
  // get all candidates for an intersection
  // Check for intersection...
  bitset<6> points[3];
  bitset<6> faceCodes[6];
  bitset<6> pCode;
  MInt edges[3][2] = {{0, 1}, {1, 2}, {0, 2}};
  // Corner points of the target regione i.e. p0=(targetRegion[targetPoints[0][0],
  //                                              targetRegion[targetPoints[0][1],
  //                                              targetRegion[targetPoints[0][2])
  MInt targetPoints[8][3] = {{0, 1, 2}, {3, 1, 2}, {0, 4, 2}, {3, 4, 2}, {0, 1, 5}, {3, 1, 5}, {0, 4, 5}, {3, 4, 5}};

  // Edges of the targetRegion (using targetPoints)
  MInt targetEdges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7},  // x-edges (i.e. parallel to x-axis)
                             {0, 2}, {1, 3}, {4, 6}, {5, 7},  // y-edes
                             {0, 4}, {1, 5}, {2, 6}, {3, 7}}; // z-edges
  MFloat e1[3];
  MFloat e2[3];
  MInt rejection = 0;
  MBool piercePointInside = 0;
  MBool triviallyAccepted = 0;

  // Each of the following arrays holds one different point for
  // all of the 6 planes. Points are built with the targetRegion
  // i.e. a = targetRegiont[pointAInPlane[0]],
  //          targetRegiont[pointAInPlane[1]],
  //          targetRegiont[pointAInPlane[2]])
  // etc.
  MInt pointAInPlane[6][3] = {{0, 1, 2}, {3, 4, 5}, {0, 1, 2}, {3, 4, 5}, {0, 1, 2}, {3, 4, 5}};

  MInt pointBInPlane[6][3] = {{0, 4, 2}, {3, 1, 5}, {3, 1, 2}, {0, 4, 5}, {3, 1, 2}, {3, 1, 5}};

  MInt pointCInPlane[6][3] = {{0, 1, 5}, {3, 4, 2}, {3, 1, 5}, {0, 4, 2}, {3, 4, 2}, {0, 1, 5}};
  MFloat a[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat b[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat c[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat d[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // start of piercing edge
  MFloat e[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // end of piercing edge
  MFloat s, gamma;                               // For pierce point calculation
  MFloat pP[3];                                  // piercePoint
  faceCodes[0] = IPOW2(0);
  faceCodes[1] = IPOW2(1);
  faceCodes[2] = IPOW2(2);
  faceCodes[3] = IPOW2(3);
  faceCodes[4] = IPOW2(4);
  faceCodes[5] = IPOW2(5);
  //   edges[0][0] = 0;
  //   edges[0][1] = 1;
  //   edges[1][0] = 1;
  //   edges[1][1] = 2;
  //   edges[2][0] = 0;
  //   edges[2][1] = 2;
  bitset<6> result;

  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  //  return noNodes;
  noReallyIntersectingNodes = 0;

  for(MInt n = 0; n < noNodes; n++) {
    // create edges (in 3D) AB, AC, BC and point A B C
    // Determine outcode (see Aftosmis, Solution Adaptive Cartesian Grid Methods for Aerodynamic flows ...)

    // Loop over all points of an element<3>
    for(MInt p = 0; p < nDim; p++) {
      points[p] = 0;
      // Calculate outcode for point
      for(MInt j = 0; j < nDim; j++) {
        if(elements[nodeList[n]].m_vertices[p][j] < targetRegion[j]) {
          points[p] |= faceCodes[2 * j];
        }
        if(elements[nodeList[n]].m_vertices[p][j] > targetRegion[j + nDim]) {
          points[p] |= faceCodes[2 * j + 1];
        }
      }
      //      m_log << points[p] << endl;
    }
    rejection = 0;
    // check outcode combinations for edges for trivial rejection
    for(auto& edge : edges) {
      if((points[edge[0]] & points[edge[1]]) != 0) {
        rejection++;
      } else {
        // If one point is inside region the element<3> is trivially accepted
        triviallyAccepted = false;
        for(auto& point : points) {
          if(point == 0) {
            triviallyAccepted = true;
            break;
          }
        }
        if(triviallyAccepted) {
          break;
        }
        // No trivial rejection, check for rejection of subsegment:
        // For all pierce points!
        // 1. Calculate pierce point:
        //    a - determine plane for pierce point calculation
        //    b - calculate pierce point
        // 2. Check for rejection of new segment
        //    a - calculate new outcode
        //    b - check for containment in pierce planes face
        // 3. If all(!) pierce points are rejected -> reject edge
        // TODO labels:GEOM This algorith might get a problem if a triangle lies completely on
        //    a face.

        // 1.a
        result = (points[edge[0]] | points[edge[1]]);
        piercePointInside = false;
        for(MInt j = 0; j < 2 * nDim; j++) {
          if(result[j] == 1) {
            // pierce plane found
            // Algorithm taken von AFTOSMIS.  There seems to be an
            // error in the calculation of s in AFTOSMIS formula, the one below
            // should be correct!
            for(MInt k = 0; k < nDim; k++) {
              a[k] = targetRegion[pointAInPlane[j][k]];
              b[k] = targetRegion[pointBInPlane[j][k]];
              c[k] = targetRegion[pointCInPlane[j][k]];
              d[k] = elements[nodeList[n]].m_vertices[edge[0]][k];
              e[k] = elements[nodeList[n]].m_vertices[edge[1]][k];
            }
            gamma = ((e[0] - d[0]) * ((a[2] - c[2]) * (c[1] - b[1]) - (a[1] - c[1]) * (c[2] - b[2]))
                     - (e[1] - d[1]) * ((a[2] - c[2]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[2] - b[2]))
                     + (e[2] - d[2]) * ((a[1] - c[1]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[1] - b[1])));

            s = ((c[0] - d[0]) * ((a[2] - c[2]) * (c[1] - b[1]) - (a[1] - c[1]) * (c[2] - b[2]))
                 - (c[1] - d[1]) * ((a[2] - c[2]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[2] - b[2]))
                 + (c[2] - d[2]) * ((a[1] - c[1]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[1] - b[1])))
                / gamma;

            // 1. b Pierce point pP in plane j:
            for(MInt k = 0; k < nDim; k++) {
              pP[k] = d[k] + s * (e[k] - d[k]);
            }
            pCode = 0;
            // 2. a Calculate outcode for pierce point
            for(MInt k = 0; k < nDim; k++) {
              if(pP[k] < targetRegion[k]) {
                pCode |= faceCodes[2 * k];
              }
              if(pP[k] > targetRegion[k + nDim]) {
                pCode |= faceCodes[2 * k + 1];
              }
            }
            //	    m_log << " outcode of pierce point :" << pCode << endl;

          } else {
            continue;
          }
          // 2. b
          result = faceCodes[j];
          result.flip();
          result = (result & pCode);
          if(result == 0) { // -> is contained
            piercePointInside = true;
          }
        }
        // reject if all pierce points are off coresponding face
        if(!piercePointInside) {
          rejection++;
        }
      }
    }
    // Check for internal element<3>, i.e. completely inside target
    result = (points[0] | points[1] | points[2]);
    if(result == 0) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
      continue;
    }
    // If not all edges are rejected a cutting element<3> has been found
    //    m_log << rejection << endl;
    if(rejection < 3) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
      continue;
    }

    // Even if trivially rejected, check if targetRegion is completely inside triangle
    // For that check all edges of the cell against the plane of the triangle,
    // calculate all pierce points and their outcodes and finally check for containment
    // in the targetRegion. Before the pierce point calculation check if the current
    // edge is parallel to the triangle (evaluate the distance of the edgepoints to the plane)
    if(rejection >= 3) {
      for(auto& targetEdge : targetEdges) { // Loop over all 12 edges of the targetRegion (cell-cube)
        for(MInt p = 0; p < nDim; p++) {
          e1[p] = targetRegion[targetPoints[targetEdge[0]][p]];
          e2[p] = targetRegion[targetPoints[targetEdge[1]][p]];
        }
        edgeTICallCounter--;
        if(edgeTriangleIntersection(elements[nodeList[n]].m_vertices[0], elements[nodeList[n]].m_vertices[1],
                                    elements[nodeList[n]].m_vertices[2], e1, e2)) {
          nodeList[noReallyIntersectingNodes] = nodeList[n];
          noReallyIntersectingNodes++;
          break;
        }
      }
    }

    // write nodes in reallyIntersectingNodes
  }
  //  if(noNodes)
  if(noNodes && !noReallyIntersectingNodes) {
    //    m_log << "Rejected ! "<< ++rejectionCounter << endl;
  }

  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}

/** \brief Determine intersection between an edge and a triangle
 *
 *
 */
inline MBool Geometry3D::edgeTriangleIntersection(MFloat* trianglePoint1, MFloat* trianglePoint2,
                                                  MFloat* trianglePoint3, MFloat* edgePoint1, MFloat* edgePoint2) {
  TRACE();

  edgeTICallCounter++;

  MFloat volume1 = 0;
  MFloat volume2 = 0;
  MFloat volume3 = 0;
  MInt maxRfnLvl = (sizeof(MInt) * 8) - 6;
  MFloat eps = 1.0 / FPOW2(maxRfnLvl);
  MFloat a[3] = {trianglePoint1[0], trianglePoint1[1], trianglePoint1[2]};
  MFloat b[3] = {trianglePoint2[0], trianglePoint2[1], trianglePoint2[2]};
  MFloat c[3] = {trianglePoint3[0], trianglePoint3[1], trianglePoint3[2]};
  MFloat d[3];
  // To determine an intersection between edge and triangle we have to:
  //  1. Check if the edge pierces the plane of the triangle
  //     - compute the signs of the signed tetrahedra build with the triangle points
  //     - and first and second edge point. Equal signs means edge doesn't pierce.
  //  2. Check if the pierce point lies inside the triangles
  //     - compute the signs of the 3 tetrahedra consisting of always the edge points
  //       and one edge of the triangle
  //     - if all signs are equal the edge pierces inside the triangle


  // 1.
  // Calculate signed tetraeder Volume of triangle and 1st point
  for(MInt i = 0; i < nDim; i++) {
    d[i] = edgePoint1[i];
  }
  volume1 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);

  for(MInt i = 0; i < nDim; i++) {
    d[i] = edgePoint2[i];
  }

  // Calculate signed tetraeder Volume of triangle and 2nd point
  volume2 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);

  if(volume1 * volume2 >= eps) {
    // Equal signs! no pierce point with triangle plane
    return false;
  }
  // 2.
  for(MInt i = 0; i < nDim; i++) {
    a[i] = edgePoint1[i];
    b[i] = trianglePoint1[i];
    c[i] = trianglePoint2[i];
    d[i] = edgePoint2[i];
  }
  volume1 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);

  for(MInt i = 0; i < nDim; i++) {
    a[i] = edgePoint1[i];
    b[i] = trianglePoint3[i];
    c[i] = trianglePoint1[i];
    d[i] = edgePoint2[i];
  }
  volume2 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);

  for(MInt i = 0; i < nDim; i++) {
    a[i] = edgePoint1[i];
    b[i] = trianglePoint2[i];
    c[i] = trianglePoint3[i];
    d[i] = edgePoint2[i];
  }
  volume3 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);
  //   if(volume1 >= eps && volume2 >= eps && volume3 >= eps)
  //     return true;
  //   if(volume1 <= eps && volume2 <= eps && volume3 <= eps)
  //     return true;

  //   if(volume1*volume1 < eps*eps || volume2*volume2 < eps*eps || volume3*volume3 < eps*eps  ){
  if(((volume1 < eps) && (volume1 > -eps)) || ((volume2 < eps) && (volume2 > -eps))
     || ((volume3 < eps) && (volume3 > -eps))) {
    m_log << " ! SINGULARITY IN TRIANGLE INTERSECTION ! At: " << a[0] << ", " << a[1] << ", " << a[2] << ";  " << d[0]
          << ", " << d[1] << ", " << d[2] << endl;
    //     cerr << " ! SINGULARITY IN TRIANGLE INTERSECTION ! At: " << a[0] << ", " << a[1] << ", " << a[2] << ";  " <<
    //     d[0] << ", " << d[1] << ", " << d[2] << endl;
    //       return true;
    return false;
  }

  if(volume1 >= 0.0 && volume2 >= 0.0 && volume3 >= 0.0) {
    return true;
  }
  if(volume1 <= 0.0 && volume2 <= 0.0 && volume3 <= 0.0) {
    return true;
  }

  return false;
}


/** \brief Determine intersection between an edge and a triangle
 *
 *
 */
inline MBool Geometry3D::edgeTriangleIntersectionLB(MFloat* trianglePoint1, MFloat* trianglePoint2,
                                                    MFloat* trianglePoint3, MFloat* edgePoint1, MFloat* edgePoint2) {
  TRACE();

  otherCalls++;

  MFloat volume1 = 0;
  MFloat volume2 = 0;
  MFloat volume3 = 0;
  MInt maxRfnLvl = (sizeof(MInt) * 8) - 6;
  MFloat eps = 1.0 / FPOW2(maxRfnLvl);
  MFloat a[3] = {trianglePoint1[0], trianglePoint1[1], trianglePoint1[2]};
  MFloat b[3] = {trianglePoint2[0], trianglePoint2[1], trianglePoint2[2]};
  MFloat c[3] = {trianglePoint3[0], trianglePoint3[1], trianglePoint3[2]};
  MFloat d[3];
  // To determine an intersection between edge and triangle we have to:
  //  1. Check if the edge pierces the plane of the triangle
  //     - compute the signs of the signed tetrahedra build with the triangle points
  //     - and first and second edge point. Equal signs means edge doesn't pierce.
  //  2. Check if the pierce point lies inside the triangles
  //     - compute the signs of the 3 tetrahedra consisting of always the edge points
  //       and one edge of the triangle
  //     - if all signs are equal the edge pierces inside the triangle


  // 1.
  // Calculate signed tetraeder Volume of triangle and 1st point
  for(MInt i = 0; i < nDim; i++) {
    d[i] = edgePoint1[i];
  }
  volume1 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);

  for(MInt i = 0; i < nDim; i++) {
    d[i] = edgePoint2[i];
  }

  // Calculate signed tetraeder Volume of triangle and 2nd point
  volume2 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);

  if(volume1 * volume2 >= eps) {
    // Equal signs! no pierce point with triangle plane
    return false;
  }
  // 2.
  for(MInt i = 0; i < nDim; i++) {
    a[i] = edgePoint1[i];
    b[i] = trianglePoint1[i];
    c[i] = trianglePoint2[i];
    d[i] = edgePoint2[i];
  }
  volume1 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);

  for(MInt i = 0; i < nDim; i++) {
    a[i] = edgePoint1[i];
    b[i] = trianglePoint3[i];
    c[i] = trianglePoint1[i];
    d[i] = edgePoint2[i];
  }
  volume2 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);

  for(MInt i = 0; i < nDim; i++) {
    a[i] = edgePoint1[i];
    b[i] = trianglePoint2[i];
    c[i] = trianglePoint3[i];
    d[i] = edgePoint2[i];
  }
  volume3 = (b[1] * a[2] * d[0] - d[1] * a[0] * c[2] + d[1] * a[0] * b[2] + c[1] * a[0] * d[2] - b[1] * a[0] * d[2]
             - a[1] * c[0] * d[2] + a[1] * d[0] * c[2] - a[1] * d[0] * b[2] + d[1] * a[2] * c[0] + d[1] * c[2] * b[0]
             - d[1] * a[2] * b[0] - d[1] * b[2] * c[0] - c[1] * d[2] * b[0] - c[1] * a[2] * d[0] + c[1] * b[2] * d[0]
             - b[1] * c[2] * d[0] + b[1] * c[0] * d[2] + a[1] * b[0] * d[2] - a[1] * b[0] * c[2] - b[1] * a[2] * c[0]
             + a[1] * c[0] * b[2] + c[1] * a[2] * b[0] - c[1] * a[0] * b[2] + b[1] * a[0] * c[2]);
  //   if(volume1 >= eps && volume2 >= eps && volume3 >= eps)
  //     return true;
  //   if(volume1 <= eps && volume2 <= eps && volume3 <= eps)
  //     return true;

  //   if(volume1*volume1 < eps*eps || volume2*volume2 < eps*eps || volume3*volume3 < eps*eps  ){
  if(((volume1 < eps) && (volume1 > -eps)) || ((volume2 < eps) && (volume2 > -eps))
     || ((volume3 < eps) && (volume3 > -eps))) {
    m_log << " ! SINGULARITY IN TRIANGLE INTERSECTION ! At: " << a[0] << ", " << a[1] << ", " << a[2] << ";  " << d[0]
          << ", " << d[1] << ", " << d[2] << endl;
    //     cerr << " ! SINGULARITY IN TRIANGLE INTERSECTION ! At: " << a[0] << ", " << a[1] << ", " << a[2] << ";  " <<
    //     d[0] << ", " << d[1] << ", " << d[2] << endl;
    return true;
    //     return false;
  }

  if(volume1 >= 0.0 && volume2 >= 0.0 && volume3 >= 0.0) return true;
  if(volume1 <= 0.0 && volume2 <= 0.0 && volume3 <= 0.0) return true;

  return false;
}


/** \brief This function determines if a line crosses a triangle.
 *
 *  \author Andreas Lintermann
 *  \date 25.01.2013
 *
 *  This is the simple case of getLineTriangleIntersection(...).
 *
 *  \param[in] p1 point 1 of the line
 *  \param[in] p2 point 2 of the line
 *  \param[in] v1 vertex 1 of the triangle
 *  \param[in] v2 vertex 2 of the triangle
 *  \param[in] v3 vertex 3 of the triangle
 *
 *  \return true if intersection is valid, false, otherwise
 *
 */
MBool Geometry3D::getLineTriangleIntersectionSimple(MFloat* p1, MFloat* p2, MFloat* v1, MFloat* v2, MFloat* v3) {
  // TRACE();

  MFloat radius = 0.0;
  MFloat intersection[3] = {0.0, 0.0, 0.0};
  MFloat normal[3] = {0.0, 0.0, 0.0};
  MFloat lambda2 = NAN;
  MFloat dist = 0.0;

  return (getLineTriangleIntersection(p1, p2, radius, v1, v2, v3, intersection, normal, &lambda2, &dist));
}

/** \brief This function determines if a line crosses a triangle.
 *
 *  \author Andreas Lintermann
 *  \date 25.01.2013
 *
 *  This is the simple case of getLineTriangleIntersection(...),
 *  additionally returns the distance.
 *
 *  \param[in] p1 point 1 of the line
 *  \param[in] p2 point 2 of the line
 *  \param[in] v1 vertex 1 of the triangle
 *  \param[in] v2 vertex 2 of the triangle
 *  \param[in] v3 vertex 3 of the triangle
 *  \param[in] dist the distance to the triangle
 *
 *  \return true if intersection is valid, false, otherwise
 *
 */
MBool Geometry3D::getLineTriangleIntersectionSimpleDistance(MFloat* p1, MFloat* p2, MFloat* v1, MFloat* v2, MFloat* v3,
                                                            MFloat* dist) {
  // TRACE();

  MFloat radius = 0.0;
  MFloat intersection[3] = {0.0, 0.0, 0.0};
  MFloat normal[3] = {0.0, 0.0, 0.0};
  MFloat lambda2 = NAN;

  return (getLineTriangleIntersection(p1, p2, radius, v1, v2, v3, intersection, normal, &lambda2, dist));
}


/** \brief This function determines if a line crosses a triangle.
 *
 *  \author Andreas Lintermann
 *  \date 08.01.2013
 *
 *  The algorithm consists of the following steps:
 *
 *  1. Get the normal of the triangle:
 *     - this is done by taking the cross-product of the spanning edges
 *     - additionally normalize for the steps 2., 2.a, 2.b
 *       + non-normalizing caused a problem, for very small triangles
 *       + the denominator in the following case was very small
 *  2. Get the cutting-point by solving the plane-equation
 *     - solve \f$\vec{n}\cdot\vec{p} = 0\f$ with \f$\vec{p} = \vec{p}_t + r * \vec{l}\f$ for \f$r\f$
 *     - \f$\vec{n}\f$ is the normal
 *     - \f$\vec{p}\f$ is a point in the plane
 *     - \f$\vec{p}_t\f$ is a point on the line (the transposed first point of the line)
 *     - \f$\vec{l}\f$ is the line
 *  2.a Check if line and triangle plane is parallel
 *     - ckeck if \f$\vec{n}\cdot\vec{l}=0\f$
 *     - if yes, proceed with 2.b, else with 3.
 *  2.b Check in 2D if we have an intersection
 *     - ckeck if both points of the line are in the plane
 *       (the line can still be just parallel)
 *     - if not in plane then we are finished
 *     - else perform 2D-test:
 *       + for each triangle edge perform cut with line
 *       + return the distance, then we are finished
 *  3. Check if the intesection point is inside the triangle
 *     - do this by evaluating the barycentric coordinates
 *
 *  \param[in] p1 point 1 of the line
 *  \param[in] p2 point 2 of the line
 *  \param[in] radius the radius of a particle
 *  \param[in] v1 vertex 1 of the triangle
 *  \param[in] v2 vertex 2 of the triangle
 *  \param[in] v3 vertex 3 of the triangle
 *  \param[in] intersection the intersection point to be returned
 *  \param[in] normal the normal of the cutting triagle to be returned
 *  \param[in] relative distance to cutpoint
 *
 *  \return true if intersection is valid, false, otherwise
 *
 */
MBool Geometry3D::getLineTriangleIntersection(const MFloat* const p1, const MFloat* const p2, const MFloat radius,
                                              const MFloat* const v1, const MFloat* const v2, const MFloat* const v3,
                                              MFloat* intersection, MFloat* normal, MFloat* lambda2, MFloat* dist) {
  // TRACE();

  std::array<MFloat, nDim> edge1{};
  std::array<MFloat, nDim> edge2{};
  std::array<MFloat, nDim> edge3{};
  std::array<MFloat, nDim> line{};
  std::array<MFloat, nDim> line_norm{};
  std::array<MFloat, nDim> normal_norm{};
  std::array<MFloat, nDim> trans_p1{};
  std::array<MFloat, nDim> trans_p2{};
  constexpr MFloat eps = 0.0000000000000001;
  constexpr MFloat eps2 = 0.00000001;

  // 1. Get the normal on the triangle
  for(MInt d = 0; d < nDim; d++) {
    edge1[d] = v2[d] - v1[d];
    edge2[d] = v3[d] - v1[d];
    edge3[d] = v3[d] - v2[d];
    line[d] = p2[d] - p1[d];

    trans_p1[d] = p1[d] - v1[d];
    trans_p2[d] = p2[d] - v1[d];
  }


  MFloat line_n_length = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    line_n_length += line[d] * line[d];
  }

  line_n_length = sqrt(line_n_length);
  const MFloat F1BLine_n_length = 1.0 / line_n_length;

  if(radius > 0.0)
    for(MInt d = 0; d < nDim; d++) {
      line[d] = (line_n_length + radius) * line[d] * F1BLine_n_length;
    }

  normal[0] = edge1[1] * edge2[2] - edge2[1] * edge1[2];
  normal[1] = edge1[2] * edge2[0] - edge2[2] * edge1[0];
  normal[2] = edge1[0] * edge2[1] - edge2[0] * edge1[1];

  MFloat n_length_sq = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    n_length_sq += normal[d] * normal[d];
  }

  // Normalize the normal
  const MFloat n_length = sqrt(n_length_sq);
  const MFloat F1BN_length = 1.0 / n_length;
  for(MInt d = 0; d < nDim; d++) {
    normal_norm[d] = normal[d] * F1BN_length;
    line_norm[d] = line[d] * F1BLine_n_length;
  }

  // 2. Get the cutting-point by solving the plane-equation

  // 2.a be careful, the dot-product of the normal and the line can 0  (they are perpendicular)
  const MFloat denom = (normal_norm[0] * line_norm[0] + normal_norm[1] * line_norm[1] + normal_norm[2] * line_norm[2]);

  // 2.b Check in 2D if we have an intersection
  if(fabs(denom) < eps) {
    // Check if the points are in the plane
    MFloat dot1 = trans_p1[0] * normal_norm[0] + trans_p1[1] * normal_norm[1] + trans_p1[2] * normal_norm[2];
    MFloat dot2 = trans_p2[0] * normal_norm[0] + trans_p2[1] * normal_norm[1] + trans_p2[2] * normal_norm[2];
    if(fabs(dot1) < eps && fabs(dot2) < eps) {
      // further testing
      MInt l_dim = 0;
      MFloat max = 0.0;
      for(MInt d = 0; d < nDim; d++) {
        if(normal_norm[d] > max) l_dim = d;
      }

      const std::array<const MFloat*, 3> verts = {v1, v2, v3};

      MFloat line2d[2] = {0.0, 0.0};
      MFloat pl_1[2] = {0.0, 0.0};
      for(MInt i = 0, j = 0; i < nDim; i++) {
        if(i != l_dim) {
          line2d[j] = p2[i] - p1[i];
          pl_1[j] = p1[i];
          j++;
        }
      }

      MBool cut = false;
      MFloat shortest = 10000000000.0;

      // for each triangle edge perform cut with line
      for(MInt d = 0; d < nDim; d++) {
        MInt next = (d + 1) % nDim;
        MFloat tri_edge[2] = {0.0, 0.0};
        MFloat tri_v1[2] = {0.0, 0.0};

        for(MInt i = 0, j = 0; i < nDim; i++) {
          if(i != l_dim) {
            tri_edge[j] = verts[next][i] - verts[d][i];
            tri_v1[j] = verts[d][i];
            j++;
          }
        }

        MFloat denom2 = tri_edge[1] * line2d[0];
        if(fabs(denom2) < eps) continue;

        MFloat denom3 = 1.0 - (line2d[1] * tri_edge[0] / denom2);
        if(fabs(denom3) < eps) continue;

        MFloat alpha =
            denom3 * (1.0 / line2d[0]) * (tri_v1[0] - pl_1[0] + (tri_edge[0] / tri_edge[1]) * (pl_1[1] - tri_v1[1]));

        // we have a cut
        if(alpha <= 1.0 && alpha >= 0.0) {
          MFloat beta = pl_1[1] / tri_edge[1] + alpha * line2d[1] / tri_edge[1] - tri_v1[1] / tri_edge[1];
          if(beta <= 1.0 && beta >= 0.0) {
            MFloat l = sqrt(alpha * alpha * line2d[0] * line2d[0] + alpha * alpha * line2d[1] * line2d[1]);
            if(l < shortest) shortest = l;
            cut = true;
            continue;
          }
        }

        // no cut
        else
          continue;
      }

      if(cut) {
        *dist = shortest;
        return true;
      } else {
        *dist = 0.0;
        return false;
      }

    } else {
      *dist = 0.0;
      return false;
    }
  }

  const MFloat r =
      (-normal_norm[0] * trans_p1[0] - normal_norm[1] * trans_p1[1] - normal_norm[2] * trans_p1[2]) / denom;

  // 3. Check if cut-point is inside the triangle
  std::array<MFloat, nDim> i_v1{};
  std::array<MFloat, nDim> i_v2{};
  std::array<MFloat, nDim> i_v3{};

  for(MInt d = 0; d < nDim; d++) {
    intersection[d] = p1[d] + r * line_norm[d];
    i_v1[d] = intersection[d] - v2[d];
    i_v2[d] = intersection[d] - v3[d];
    i_v3[d] = intersection[d] - v1[d];
  }

  *dist = r;

  if(r >= 0.0)
    *lambda2 = r * F1BLine_n_length;
  else
    *lambda2 = 1.1;

  // we are already finished if the cut-point has a larger distance than the second point
  if(fabs(*lambda2) > 1.0) return false;

  std::array<MFloat, nDim> normal1;
  std::array<MFloat, nDim> normal2;
  std::array<MFloat, nDim> normal3;

  normal1[0] = edge3[1] * i_v1[2] - i_v1[1] * edge3[2];
  normal1[1] = edge3[2] * i_v1[0] - i_v1[2] * edge3[0];
  normal1[2] = edge3[0] * i_v1[1] - i_v1[0] * edge3[1];

  normal2[0] = -edge2[1] * i_v2[2] + i_v2[1] * edge2[2];
  normal2[1] = -edge2[2] * i_v2[0] + i_v2[2] * edge2[0];
  normal2[2] = -edge2[0] * i_v2[1] + i_v2[0] * edge2[1];

  normal3[0] = edge1[1] * i_v3[2] - i_v3[1] * edge1[2];
  normal3[1] = edge1[2] * i_v3[0] - i_v3[2] * edge1[0];
  normal3[2] = edge1[0] * i_v3[1] - i_v3[0] * edge1[1];

  MFloat alpha;
  MFloat beta;
  MFloat gamma;
  alpha = beta = gamma = 0.0;

  for(MInt d = 0; d < nDim; d++) {
    alpha += normal[d] * normal1[d];
    beta += normal[d] * normal2[d];
    gamma += normal[d] * normal3[d];
  }

  alpha /= n_length_sq;
  beta /= n_length_sq;
  gamma /= n_length_sq;

  // check for inside
  if((alpha + beta + gamma <= 1.0 + eps2) && (alpha + beta + gamma >= 1.0 - eps2))
    if(alpha >= 0 && beta >= 0 && gamma >= 0)
      return true;
    else
      return false;
  else
    return false;
}


//----------------------------------------------------------------------------------------


/** \brief Returns the ids of all elements, cut by a orthogonal line (or rectangular region)
    \author Daniel Hartmann
    \date unknown

     Test is evaluated by the barycentric coordinates of the intersection point. Therefore a 3 by
     3 linear equation system is solved with Cramers rule. This function replaces
   Geometry3D::getLineIntersectionElementsOld1. \bug labels:GEOM This method does not run on IBM_BLUE_GENE and is
   hence not used anymore. This is because the determinant in the denominator of Cramers rule can become 0, which can
   not be handled by IBM_BLUE_GENE. Please use Geometry3D::getLineIntersectionElements instead.
  */
MInt Geometry3D::getLineIntersectionElementsOld2(MFloat* targetRegion, MInt* spaceDirection,
                                                 std::vector<MInt>& nodeList) {
  //  TRACE();

  m_getLIE2CallCounter++;

  MBool singleIntersectionPoint = 0;
  MInt noReallyIntersectingNodes = 0;
  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  MInt maxNoNodes = 120;
  MInt spaceId;
  MInt spaceId1;
  MInt spaceId2;
  array<MFloat, nDim> e1{};
  array<MFloat, nDim> e2{};
  MFloat epsilon = 0.00000000001;
  array<MFloat, nDim> a{};
  array<MFloat, nDim> b{};
  array<MFloat, nDim> c{};
  MFloat gamma;
  MFloat s;
  MFloat p;
  MFloat q;
  array<MFloat, nDim> pP{};
  auto** temp = new MFloat*[maxNoNodes];
  for(MInt k = 0; k < maxNoNodes; k++) {
    temp[k] = new MFloat[nDim];
  }

  for(MInt i = 0; i < nDim; i++) {
    e1[i] = targetRegion[i];
    e2[i] = targetRegion[i + nDim];
  }
  spaceId = spaceDirection[0];
  spaceId1 = spaceDirection[1];
  spaceId2 = spaceDirection[2];

  for(MInt n = 0; n < noNodes; n++) {
    for(MInt k = 0; k < nDim; k++) {
      a[k] = elements[nodeList[n]].m_vertices[0][k];
      b[k] = elements[nodeList[n]].m_vertices[1][k];
      c[k] = elements[nodeList[n]].m_vertices[2][k];
    }
    if(approx(a[spaceId1], b[spaceId1], MFloatEps) && !approx(a[spaceId1], c[spaceId1], MFloatEps)) {
      // aspaceId != bspaceId, otherwise a and b would be the same point
      q = (e1[spaceId1] - a[spaceId1]) / (c[spaceId1] - a[spaceId1]);
      p = (e1[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
    } else {
      if(!approx(a[spaceId1], b[spaceId1], MFloatEps) && approx(a[spaceId1], c[spaceId1], MFloatEps)) {
        // aspaceId != cspaceId, otherwise a and c would be the same point
        p = (e1[spaceId1] - a[spaceId1]) / (b[spaceId1] - a[spaceId1]);
        q = (e1[spaceId] - a[spaceId] - p * (b[spaceId] - a[spaceId])) / (c[spaceId] - a[spaceId]);
      } else {
        if(approx(a[spaceId], b[spaceId], MFloatEps) && !approx(a[spaceId], c[spaceId], MFloatEps)) {
          // aspaceId1 != bspaceId1, otherwise a and b would be the same point
          q = (e1[spaceId] - a[spaceId]) / (c[spaceId] - a[spaceId]);
          p = (e1[spaceId1] - a[spaceId1] - q * (c[spaceId1] - a[spaceId1])) / (b[spaceId1] - a[spaceId1]);
        } else {
          if(!approx(a[spaceId], b[spaceId], MFloatEps) && approx(a[spaceId], c[spaceId], MFloatEps)) {
            // aspaceId1 != cspaceId1, otherwise a and c would be the same point
            p = (e1[spaceId] - a[spaceId]) / (b[spaceId] - a[spaceId]);
            q = (e1[spaceId1] - a[spaceId1] - p * (b[spaceId1] - a[spaceId1])) / (c[spaceId1] - a[spaceId1]);
          } else {
            // aspaceId1 != bspaceId1 && aspaceId1 != cspaceId1 && aspaceId != bspaceId && aspaceId != cspaceId
            q = ((e1[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                 - (b[spaceId1] - a[spaceId1]) * (e1[spaceId] - a[spaceId]))
                / ((c[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                   - (b[spaceId1] - a[spaceId1]) * (c[spaceId] - a[spaceId]));
            p = (e1[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
          }
        }
      }
    }

    if(p * q >= 0 || p * q < 0) {
      // compute s
      gamma = a[spaceId2] + p * (b[spaceId2] - a[spaceId2]) + q * (c[spaceId2] - a[spaceId2]);
      s = (gamma - e1[spaceId2]) / (e2[spaceId2] - e1[spaceId2]);

      if(s < -epsilon || s > F1 + epsilon || p < -epsilon || q < -epsilon || (p + q) > F1 + epsilon) {
      } else {
        singleIntersectionPoint = true;

        // cut point pP
        for(MInt g = 0; g < nDim; g++) {
          pP[g] = e1[g] + s * (e2[g] - e1[g]);
          temp[noReallyIntersectingNodes][g] = pP[g];
        }

        if(noReallyIntersectingNodes == 0) {
          nodeList[noReallyIntersectingNodes] = nodeList[n];
          noReallyIntersectingNodes++;
        } else {
          for(MInt f = 0; f < noReallyIntersectingNodes; f++) {
            if((fabs(temp[f][0] - pP[0]) < epsilon) && (fabs(temp[f][1] - pP[1]) < epsilon)
               && (fabs(temp[f][2] - pP[2]) < epsilon)) {
              singleIntersectionPoint = false;
              break;
            }
          }
          if(singleIntersectionPoint) {
            nodeList[noReallyIntersectingNodes] = nodeList[n];
            noReallyIntersectingNodes++;
          }
        }
      }
    }
  }
  nodeList.resize(noReallyIntersectingNodes);

  // free memory
  for(MInt k = 0; k < maxNoNodes; k++) {
    delete[] temp[k];
  }
  delete[] temp;

  getLIE2CommCounter += noReallyIntersectingNodes;

  return noReallyIntersectingNodes;
}

/** \brief Returns the ids of all elements, cut by a ray by using the perp-dot operator.
    \author Andreas Lintermann
    \date 24.11.2009

    Let the following be definded:

    <ul>
    <li>triangle = \f$v0,v1,v2\f$</li>
    <li>edge = \f$e1,e2\f$</li>
    <li>normal on triangle plane = \f$n\f$</li>
    <li>\f$v1-v0 = u\f$</li>
    <li>\f$v2-v0 = v\f$</li>
    <li>intersection point = \f$w = su + tv\f$</li>
    </ul>

    The perp-dot operator is defined by

    \f{eqnarray*}{
    v^\perp = n \times v = (u*v)v - (v*v)u\\
    u^\perp = n \times u = (u*u)v - (u*v)u
    \f}
    Then we get:
    \f{eqnarray*}{
    & w & = su + tv\\
    \Leftrightarrow & n * v^\perp & = su * v^\perp + tv * v^\perp = su * v^\perp\\
    \Leftrightarrow & s & = \frac{w\times v^\perp}{u\times v^\perp}\\
    & & =\frac{(u*v)(w*v) - (v*v)(w*u)}{(u*v)(u*v) - (u*u)(v*v)}
    \f}
    \f{eqnarray*}{
    & w & = su + tv\\
    \Leftrightarrow & n * u^\perp & = su * u^\perp + tv * u^\perp = su * u^\perp\\
    \Leftrightarrow & t & = \frac{w\times u^\perp}{u\times u^\perp}\\
    & & =\frac{(u*v)(w*u) - (u*u)(w*v)}{(u*v)(u*v) - (u*u)(v*v)}
    \f}
*/
// TODO labels:GEOM unify with version without nodeList parameter!
MInt Geometry3D::getLineIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  // TRACE();

  m_getLIE2CallCounter++;

  const MFloat eps = 0.00000000001;

  MFloat u[nDim];
  MFloat v[nDim];
  MFloat w[nDim];

  constexpr MInt maxNoNodes = 20; // Should prevent vector reallocation in most cases
  std::vector<MFloat> ips{};
  ips.reserve(nDim * maxNoNodes);

  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  MFloat e1[3];
  MFloat ray[nDim];
  for(MInt i = 0; i < nDim; i++) {
    e1[i] = targetRegion[i];
    ray[i] = targetRegion[i + nDim] - e1[i];
  }

  MInt noReallyIntersectingNodes = 0;
  for(MInt n = 0; n < noNodes; n++) {
    const MInt nodeId = nodeList[n];
    MFloat** const verts = elements[nodeId].m_vertices;
    const MFloat* const vert0 = verts[0];
    const MFloat* const vert1 = verts[1];
    const MFloat* const vert2 = verts[2];

    // Fill triangle edges
    for(MInt i = 0; i < nDim; i++) {
      const MFloat v0 = vert0[i];
      u[i] = vert1[i] - v0;
      v[i] = vert2[i] - v0;
      w[i] = e1[i] - v0;
    }

    // Get normal
    const MFloat normal0 = u[1] * v[2] - u[2] * v[1];
    const MFloat normal1 = u[2] * v[0] - u[0] * v[2];
    const MFloat normal2 = u[0] * v[1] - u[1] * v[0];

    // dot products
    const MFloat a = -(normal0 * w[0] + normal1 * w[1] + normal2 * w[2]);
    const MFloat b = normal0 * ray[0] + normal1 * ray[1] + normal2 * ray[2];

    if(fabs(b) < eps) {
      if(fabs(a) < eps) {
        // ray lies in triangle plane
        // TODO labels:GEOM commented the following two lines; check what should be done in such case as ips is not
        // filled for the check later; occurs in case FV/3D_MP_postprocessing_coarse - no difference in results
        /* nodeList[noReallyIntersectingNodes] = nodeList[n]; */
        /* noReallyIntersectingNodes++; */
        m_log << "Found ray in triangle plane" << endl;
        continue;
      } else {
        // ray disjoint from plane
        continue;
      }
    }
    const MFloat r = a / b;

    const MFloat F1eps = F1 + eps;
    // is the ray going the right direction and is the scaling factor smaller than 1.0+eps
    if(r < -eps || r > F1eps) continue;


    const MFloat uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    const MFloat uv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    const MFloat vv = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];

    MFloat wu = 0.0;
    MFloat wv = 0.0;
    MFloat ip[nDim];
    // dot products
    for(MInt i = 0; i < nDim; i++) {
      ip[i] = e1[i] + r * ray[i];
      const MFloat tmp = ip[i] - vert0[i];
      wu += tmp * u[i];
      wv += tmp * v[i];
    }

    const MFloat D = uv * uv - uu * vv;

    const MFloat s = (uv * wv - vv * wu) / D;
    if(s < -eps || s > F1eps)
      // Intersection point is outside triangle
      continue;

    const MFloat t = (uv * wu - uu * wv) / D;
    if(t < -eps || (s + t) > F1eps)
      // Intersection point is outside triangle
      continue;

    // Here we are checking if such an intersection point has already been found.
    // This can be the case if the line intersects a triange edge or a vertex.
    if(noReallyIntersectingNodes == 0) {
      ASSERT(ips.empty(), "ips not empty");
      for(MInt i = 0; i < nDim; i++) {
        ips.push_back(ip[i]);
      }
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
    } else {
      MBool singleIntersectionPoint = true;
      for(MInt f = 0; f < noReallyIntersectingNodes; f++) {
        if((fabs(ips[f * nDim] - ip[0]) < eps) && (fabs(ips[f * nDim + 1] - ip[1]) < eps)
           && (fabs(ips[f * nDim + 2] - ip[2]) < eps)) {
          singleIntersectionPoint = false;
          // m_log << "Found multiple points" << endl;
          break;
        }
      }
      if(singleIntersectionPoint) {
        for(MInt i = 0; i < nDim; i++) {
          ips.push_back(ip[i]);
        }
        nodeList[noReallyIntersectingNodes] = nodeList[n];
        noReallyIntersectingNodes++;
        ASSERT(ips.size() == (MUlong)(noReallyIntersectingNodes * nDim),
               "Error: wrong vector size of ips. " + std::to_string(ips.size())
                   + " != " + std::to_string(noReallyIntersectingNodes * nDim));
      }
    }
  }

  nodeList.resize(noReallyIntersectingNodes);

  getLIE2CommCounter += noReallyIntersectingNodes;

  return noReallyIntersectingNodes;
}


MInt Geometry3D::getLineIntersectionElements(MFloat* targetRegion) {
  //  TRACE();
  TERMM(1, "should not be used anymore");

  m_getLIE2CallCounter++;

  const MFloat eps = 0.00000000001;

  MFloat u[nDim];
  MFloat v[nDim];
  MFloat normal[nDim];
  MFloat w[nDim];
  MFloat ray[nDim];
  MFloat* ip = nullptr;

  for(MInt i = 0; i < nDim; i++) {
    u[i] = v[i] = w[i] = numeric_limits<MFloat>::max();
  }

  stack<MFloat*> allIPs;

  MFloat a;
  MFloat b;
  MFloat r;
  MFloat D;
  MFloat uu;
  MFloat uv;
  MFloat vv;
  MFloat wu;
  MFloat wv;
  MFloat s;
  MFloat t;

  MInt noReallyIntersectingNodes = 0;

  std::vector<MInt> nodeList;

  MInt maxNoNodes = 30;
  auto** tmp = new MFloat*[maxNoNodes];
  for(MInt k = 0; k < maxNoNodes; k++) {
    tmp[k] = new MFloat[nDim];
  }
  MBool singleIntersectionPoint;

  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();

  MFloat e1[3];
  MFloat e2[3];
  for(MInt i = 0; i < nDim; i++) {
    e1[i] = targetRegion[i];
    e2[i] = targetRegion[i + nDim];
    ray[i] = e2[i] - e1[i];
  }

  for(MInt n = 0; n < noNodes; n++) {
    // Init dot products
    a = b = uu = uv = vv = wu = wv = 0.0;

    // Fill triangle edges
    for(MInt i = 0; i < nDim; i++) {
      u[i] = elements[nodeList[n]].m_vertices[1][i] - elements[nodeList[n]].m_vertices[0][i];
      v[i] = elements[nodeList[n]].m_vertices[2][i] - elements[nodeList[n]].m_vertices[0][i];
      w[i] = e1[i] - elements[nodeList[n]].m_vertices[0][i];
    }

    // Get normal
    normal[0] = u[1] * v[2] - u[2] * v[1];
    normal[1] = u[2] * v[0] - u[0] * v[2];
    normal[2] = u[0] * v[1] - u[1] * v[0];

    // dot products
    for(MInt i = 0; i < nDim; i++) {
      a += normal[i] * w[i];
      b += normal[i] * ray[i];
    }
    a *= -1;

    if(fabs(b) < eps) {
      if(fabs(a) < eps) {
        // ray lies in triangle plane
        nodeList[noReallyIntersectingNodes] = nodeList[n];
        //	noReallyIntersectingNodes++;
        m_log << "Found ray in triangle plane" << endl;
        continue;
      } else
        // ray disjoint from plane
        continue;
    }
    r = a / b;

    // is the ray going the right direction and is the scaling factor smaller than 1.0+eps
    if(r < -eps || r > F1 + eps) continue;

    // dot products
    ip = new MFloat[3];

    for(MInt i = 0; i < nDim; i++) {
      ip[i] = e1[i] + r * ray[i];
      uu += u[i] * u[i];
      uv += u[i] * v[i];
      vv += v[i] * v[i];
      w[i] = ip[i] - elements[nodeList[n]].m_vertices[0][i];
      wu += w[i] * u[i];
      wv += w[i] * v[i];
    }

    D = uv * uv - uu * vv;

    s = (uv * wv - vv * wu) / D;
    if(s < -eps || s > F1 + eps)
      // Intersection point is outside triangle
      continue;
    t = (uv * wu - uu * wv) / D;
    if(t < -eps || (s + t) > F1 + eps)
      // Intersection point is outside triangle
      continue;

    // Here we are checking if there are multiple intersection points.
    // This can be the case if the line intersects a triangle edge or a vertex.

    singleIntersectionPoint = true;

    if(noReallyIntersectingNodes == 0) {
      for(MInt i = 0; i < nDim; i++)
        tmp[noReallyIntersectingNodes][i] = ip[i];

      nodeList[noReallyIntersectingNodes] = nodeList[n];
      allIPs.push(tmp[noReallyIntersectingNodes]);

      noReallyIntersectingNodes++;
    } else {
      for(MInt f = 0; f < noReallyIntersectingNodes; f++) {
        if((fabs(tmp[f][0] - ip[0]) < eps) && (fabs(tmp[f][1] - ip[1]) < eps) && (fabs(tmp[f][2] - ip[2]) < eps)) {
          singleIntersectionPoint = false;
          // m_log << "Found multiple points" << endl;
          break;
        }
      }

      if(singleIntersectionPoint) {
        for(MInt i = 0; i < nDim; i++)
          tmp[noReallyIntersectingNodes][i] = ip[i];

        nodeList[noReallyIntersectingNodes] = nodeList[n];
        noReallyIntersectingNodes++;

      } else {
        delete[] ip;
      }
    }
  }

  nodeList.resize(noReallyIntersectingNodes);

  // free the temporary space
  for(MInt k = 0; k < maxNoNodes; k++) {
    delete[] tmp[k];
    tmp[k] = 0;
  }
  delete[] tmp;
  tmp = 0;

  getLIE2CommCounter += noReallyIntersectingNodes;

  return noReallyIntersectingNodes;
}


//==================================================
//==================================================


MBool Geometry3D::getClosestLineIntersectionLength(MInt bndCndId, const std::vector<MInt>& nodeList,
                                                   MFloat* targetRegion, MFloat* dist) {
  TRACE();

  MBool cut = false;
  *dist = 100000000000.0;

  MFloat tmp_length;

  MFloat eps = 0.00000000001;

  MFloat u[nDim];
  MFloat v[nDim];
  MFloat normal[nDim];
  MFloat w[nDim];
  MFloat ray[nDim];
  MFloat ip[nDim];

  for(MInt i = 0; i < nDim; i++) {
    ip[i] = w[i] = u[i] = v[i] = numeric_limits<MFloat>::max();
  }

  MFloat a, b, r;
  MFloat D, uu, uv, vv, wu, wv;
  MFloat s, t;

  MInt noReallyIntersectingNodes = 0;


  MInt maxNoNodes = 30;
  MFloat** tmp = new MFloat*[maxNoNodes];
  for(MInt k = 0; k < maxNoNodes; k++) {
    tmp[k] = new MFloat[nDim];
  }
  MBool singleIntersectionPoint;

  MFloat e1[3], e2[3];
  for(MInt i = 0; i < nDim; i++) {
    e1[i] = targetRegion[i];
    e2[i] = targetRegion[i + nDim];
    ray[i] = e2[i] - e1[i];
  }

  for(MInt n = 0; n < (signed)nodeList.size(); n++) {
    if(elements[nodeList[n]].m_bndCndId != bndCndId) continue;

    tmp_length = 0.0;

    // Init dot products
    a = b = uu = uv = vv = wu = wv = 0.0;

    // Fill triangle edges
    for(MInt i = 0; i < nDim; i++) {
      u[i] = elements[nodeList[n]].m_vertices[1][i] - elements[nodeList[n]].m_vertices[0][i];
      v[i] = elements[nodeList[n]].m_vertices[2][i] - elements[nodeList[n]].m_vertices[0][i];
      w[i] = e1[i] - elements[nodeList[n]].m_vertices[0][i];
    }

    // Get normal
    normal[0] = u[1] * v[2] - u[2] * v[1];
    normal[1] = u[2] * v[0] - u[0] * v[2];
    normal[2] = u[0] * v[1] - u[1] * v[0];

    // dot products
    for(MInt i = 0; i < nDim; i++) {
      a += normal[i] * w[i];
      b += normal[i] * ray[i];
    }
    a *= -1;

    if(fabs(b) < eps) {
      if(fabs(a) < eps) {
        // TODO labels:GEOM
        // ray lies in triangle plane
        noReallyIntersectingNodes++;
        cut = true;
        continue;
      } else
        // ray disjoint from plane
        continue;
    }
    r = a / b;

    // is the ray going the right direction and is the scaling factor smaller than 1.0+eps
    if(r < -eps || r > F1 + eps) continue;

    // dot products
    for(MInt i = 0; i < nDim; i++) {
      ip[i] = e1[i] + r * ray[i];
      uu += u[i] * u[i];
      uv += u[i] * v[i];
      vv += v[i] * v[i];
      w[i] = ip[i] - elements[nodeList[n]].m_vertices[0][i];
      wu += w[i] * u[i];
      wv += w[i] * v[i];
    }

    D = uv * uv - uu * vv;

    s = (uv * wv - vv * wu) / D;
    if(s < -eps || s > F1 + eps)
      // Intersection point is outside triangle
      continue;
    t = (uv * wu - uu * wv) / D;
    if(t < -eps || (s + t) > F1 + eps)
      // Intersection point is outside triangle
      continue;

    // Here we are checking if such an intersection point has already been found.
    // This can be the case if the line intersects a triange edge or a vertex.
    singleIntersectionPoint = true;

    if(noReallyIntersectingNodes == 0) {
      for(MInt i = 0; i < nDim; i++)
        tmp[noReallyIntersectingNodes][i] = ip[i];
      noReallyIntersectingNodes++;
      /*
      for(MInt dim=0; dim<nDim; dim++)
    tmp_length+=(ip[dim]-e1[dim])*(ip[dim]-e1[dim]);
      */
      for(MInt dim = 0; dim < nDim; dim++)
        tmp_length += (r * ray[dim]) * (r * ray[dim]);

      tmp_length = sqrt(tmp_length);

      if(tmp_length < *dist) *dist = tmp_length;

      cut = true;
    } else {
      for(MInt f = 0; f < noReallyIntersectingNodes; f++) {
        if((fabs(tmp[f][0] - ip[0]) < eps) && (fabs(tmp[f][1] - ip[1]) < eps) && (fabs(tmp[f][2] - ip[2]) < eps)) {
          singleIntersectionPoint = false;
          break;
        }
      }
      if(singleIntersectionPoint) {
        for(MInt i = 0; i < nDim; i++)
          tmp[noReallyIntersectingNodes][i] = ip[i];
        noReallyIntersectingNodes++;

        /*
        for(MInt dim=0; dim<nDim; dim++)
        tmp_length+=(ip[dim]-e1[dim])*(ip[dim]-e1[dim]);
        */
        for(MInt dim = 0; dim < nDim; dim++)
          tmp_length += (r * ray[dim]) * (r * ray[dim]);

        tmp_length = sqrt(tmp_length);


        if(tmp_length < *dist) *dist = tmp_length;

        cut = true;
      }
    }
  }

  // free the temporary space
  for(MInt k = 0; k < maxNoNodes; k++) {
    delete[] tmp[k];
    tmp[k] = nullptr;
  }
  delete[] tmp;
  tmp = nullptr;

  return cut;
}


MInt Geometry3D::getIntersectionMBElements(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  TRACE();

  otherCalls++;


  MInt noReallyIntersectingNodes = 0;
  // get all candidates for an intersection
  // Check for intersection...
  bitset<6> points[3];
  bitset<6> faceCodes[6];
  bitset<6> pCode;
  MInt spaceId;
  MInt spaceId1;
  MInt spaceId2;
  MInt edges[3][2] = {{0, 1}, {1, 2}, {0, 2}};
  // Corner points of the target regione i.e. p0=(targetRegion[targetPoints[0][0],
  //                                              targetRegion[targetPoints[0][1],
  //                                              targetRegion[targetPoints[0][2])
  MInt targetPoints[8][3] = {{0, 1, 2}, {3, 1, 2}, {0, 4, 2}, {3, 4, 2}, {0, 1, 5}, {3, 1, 5}, {0, 4, 5}, {3, 4, 5}};

  // Edges of the targetRegion (using targetPoints)
  MInt targetEdges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7},  // x-edges (i.e. parallel to x-axis)
                             {0, 2}, {1, 3}, {4, 6}, {5, 7},  // y-edes
                             {0, 4}, {1, 5}, {2, 6}, {3, 7}}; // z-edges
  MFloat e1[3];
  MFloat e2[3];
  MInt rejection;
  MBool piercePointInside;
  MBool triviallyAccepted;

  // Each of the following arrays holds one different point for
  // all of the 6 planes. Points are built with the targetRegion
  // i.e. a = targetRegiont[pointAInPlane[0]],
  //          targetRegiont[pointAInPlane[1]],
  //          targetRegiont[pointAInPlane[2]])
  // etc.
  MInt pointAInPlane[6][3] = {{0, 1, 2}, {3, 4, 5}, {0, 1, 2}, {3, 4, 5}, {0, 1, 2}, {3, 4, 5}};

  MInt pointBInPlane[6][3] = {{0, 4, 2}, {3, 1, 5}, {3, 1, 2}, {0, 4, 5}, {3, 1, 2}, {3, 1, 5}};

  MInt pointCInPlane[6][3] = {{0, 1, 5}, {3, 4, 2}, {3, 1, 5}, {0, 4, 2}, {3, 4, 2}, {0, 1, 5}};
  MFloat a[3]; // point in plane
  MFloat b[3]; // point in plane
  MFloat c[3]; // point in plane
  MFloat d[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // start of piercing edge
  MFloat e[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // end of piercing edge
  MFloat s;
  MFloat gamma; // For pierce point calculation
  MFloat p2;
  MFloat q;
  MFloat epsilon = 0.00000000001;
  MFloat pP[3]; // piercePoint
  faceCodes[0] = IPOW2(0);
  faceCodes[1] = IPOW2(1);
  faceCodes[2] = IPOW2(2);
  faceCodes[3] = IPOW2(3);
  faceCodes[4] = IPOW2(4);
  faceCodes[5] = IPOW2(5);
  //   edges[0][0] = 0;
  //   edges[0][1] = 1;
  //   edges[1][0] = 1;
  //   edges[1][1] = 2;
  //   edges[2][0] = 0;
  //   edges[2][1] = 2;
  bitset<6> result;
  m_adt->retrieveNodesMBElements(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  //  return noNodes;
  noReallyIntersectingNodes = 0;

  for(MInt n = 0; n < noNodes; n++) {
    // create edges (in 3D) AB, AC, BC and point A B C
    // Determine outcode (see Aftosmis, Solution Adaptive Cartesian Grid Methods for Aerodynamic flows ...)

    // Loop over all points of an element<3>
    for(MInt p = 0; p < nDim; p++) {
      points[p] = 0;
      // Calculate outcode for point
      for(MInt j = 0; j < nDim; j++) {
        if(mbelements[nodeList[n]].m_vertices[p][j] < targetRegion[j]) {
          points[p] |= faceCodes[2 * j];
        }
        if(mbelements[nodeList[n]].m_vertices[p][j] > targetRegion[j + nDim]) {
          points[p] |= faceCodes[2 * j + 1];
        }
      }
      //      m_log << points[p] << endl;
    }
    rejection = 0;
    // check outcode combinations for edges for trivial rejection
    for(auto& edge : edges) {
      if((points[edge[0]] & points[edge[1]]) != 0) {
        rejection++;
      } else {
        // If one point is inside region the element<3> is trivially accepted
        triviallyAccepted = false;
        for(MInt k = 0; k < nDim; k++) {
          if(points[k] == 0) {
            triviallyAccepted = true;
            break;
          }
        }
        if(triviallyAccepted) {
          break;
        }
        // No trivial rejection, check for rejection of subsegment:
        // For all pierce points!
        // 1. Calculate pierce point:
        //    a - determine plane for pierce point calculation
        //    b - calculate pierce point
        // 2. Check for rejection of new segment
        //    a - calculate new outcode
        //    b - check for containment in pierce planes face
        // 3. If all(!) pierce points are rejected -> reject edge
        // TODO labels:GEOM This algorith might get a problem if a triangle lies completely on
        //    a face.

        // 1.a
        result = (points[edge[0]] | points[edge[1]]);
        piercePointInside = false;
        for(MInt j = 0; j < 2 * nDim; j++) {
          if(result[j] == 1) {
            // pierce plane found
            // Algorithm taken von AFTOSMIS.  There seems to be an
            // error in the calculation of s in AFTOSMIS formula, the one below
            // should be correct!
            for(MInt k = 0; k < nDim; k++) {
              a[k] = targetRegion[pointAInPlane[j][k]];
              b[k] = targetRegion[pointBInPlane[j][k]];
              c[k] = targetRegion[pointCInPlane[j][k]];
              d[k] = mbelements[nodeList[n]].m_vertices[edge[0]][k];
              e[k] = mbelements[nodeList[n]].m_vertices[edge[1]][k];
            }
            gamma = ((e[0] - d[0]) * ((a[2] - c[2]) * (c[1] - b[1]) - (a[1] - c[1]) * (c[2] - b[2]))
                     - (e[1] - d[1]) * ((a[2] - c[2]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[2] - b[2]))
                     + (e[2] - d[2]) * ((a[1] - c[1]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[1] - b[1])));

            s = ((c[0] - d[0]) * ((a[2] - c[2]) * (c[1] - b[1]) - (a[1] - c[1]) * (c[2] - b[2]))
                 - (c[1] - d[1]) * ((a[2] - c[2]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[2] - b[2]))
                 + (c[2] - d[2]) * ((a[1] - c[1]) * (c[0] - b[0]) - (a[0] - c[0]) * (c[1] - b[1])))
                / gamma;

            // 1. b Pierce point pP in plane j:
            for(MInt k = 0; k < nDim; k++) {
              pP[k] = d[k] + s * (e[k] - d[k]);
            }
            pCode = 0;
            // 2. a Calculate outcode for pierce point
            for(MInt k = 0; k < nDim; k++) {
              if(pP[k] < targetRegion[k]) {
                pCode |= faceCodes[2 * k];
              }
              if(pP[k] > targetRegion[k + nDim]) {
                pCode |= faceCodes[2 * k + 1];
              }
            }
            //	    m_log << " outcode of pierce point :" << pCode << endl;

          } else {
            continue;
          }
          // 2. b
          result = faceCodes[j];
          result.flip();
          result = (result & pCode);
          if(result == 0) { // -> is contained
            piercePointInside = true;
          }
        }
        // reject if all pierce points are off coresponding face
        if(!piercePointInside) {
          rejection++;
        }
      }
    }
    // Check for internal element<3>, i.e. completely inside target
    result = (points[0] | points[1] | points[2]);
    if(result == 0) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
      continue;
    }
    // If not all edges are rejected a cutting element<3> has been found
    //    m_log << rejection << endl;
    if(rejection < 3) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
      continue;
    }

    // Even if trivially rejected, check if targetRegion is completely inside triangle
    // For that check all edges of the cell against the plane of the triangle,
    // calculate all pierce points and their outcodes and finally check for containment
    // in the targetRegion. Before the pierce point calculation check if the current
    // edge is parallel to the triangle (evaluate the distance of the edgepoints to the plane)
    if(rejection >= 3) {
      for(MInt t = 0; t < 12; t++) { // Loop over all 12 edges of the targetRegion (cell-cube)
        for(MInt p = 0; p < nDim; p++) {
          e1[p] = targetRegion[targetPoints[targetEdges[t][0]][p]];
          e2[p] = targetRegion[targetPoints[targetEdges[t][1]][p]];
        }

        if(t < 4) {
          spaceId = 1;
          spaceId1 = 2;
          spaceId2 = 0;
        } else {
          if(t > 3 && t < 8) {
            spaceId = 2;
            spaceId1 = 0;
            spaceId2 = 1;
          } else {
            spaceId = 0;
            spaceId1 = 1;
            spaceId2 = 2;
          }
        }

        for(MInt k = 0; k < nDim; k++) {
          a[k] = mbelements[nodeList[n]].m_vertices[0][k];
          b[k] = mbelements[nodeList[n]].m_vertices[1][k];
          c[k] = mbelements[nodeList[n]].m_vertices[2][k];
        }
        if(approx(a[spaceId1], b[spaceId1], MFloatEps) && !approx(a[spaceId1], c[spaceId1], MFloatEps)) {
          // aspaceId != bspaceId, otherwise a and b would be the same point
          q = (e1[spaceId1] - a[spaceId1]) / (c[spaceId1] - a[spaceId1]);
          p2 = (e1[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
        } else {
          if(!approx(a[spaceId1], b[spaceId1], MFloatEps) && approx(a[spaceId1], c[spaceId1], MFloatEps)) {
            // aspaceId != cspaceId, otherwise a and c would be the same point
            p2 = (e1[spaceId1] - a[spaceId1]) / (b[spaceId1] - a[spaceId1]);
            q = (e1[spaceId] - a[spaceId] - p2 * (b[spaceId] - a[spaceId])) / (c[spaceId] - a[spaceId]);
          } else {
            if(approx(a[spaceId], b[spaceId], MFloatEps) && !approx(a[spaceId], c[spaceId], MFloatEps)) {
              // aspaceId1 != bspaceId1, otherwise a and b would be the same point
              q = (e1[spaceId] - a[spaceId]) / (c[spaceId] - a[spaceId]);
              p2 = (e1[spaceId1] - a[spaceId1] - q * (c[spaceId1] - a[spaceId1])) / (b[spaceId1] - a[spaceId1]);
            } else {
              if(!approx(a[spaceId], b[spaceId], MFloatEps) && approx(a[spaceId], c[spaceId], MFloatEps)) {
                // aspaceId1 != cspaceId1, otherwise a and c would be the same point
                p2 = (e1[spaceId] - a[spaceId]) / (b[spaceId] - a[spaceId]);
                q = (e1[spaceId1] - a[spaceId1] - p2 * (b[spaceId1] - a[spaceId1])) / (c[spaceId1] - a[spaceId1]);
              } else {
                // aspaceId1 != bspaceId1 && aspaceId1 != cspaceId1 && aspaceId != bspaceId && aspaceId != cspaceId
                q = ((e1[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                     - (b[spaceId1] - a[spaceId1]) * (e1[spaceId] - a[spaceId]))
                    / ((c[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                       - (b[spaceId1] - a[spaceId1]) * (c[spaceId] - a[spaceId]));
                p2 = (e1[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
              }
            }
          }
        }

        if(p2 * q >= 0 || p2 * q < 0) {
          // compute s
          gamma = a[spaceId2] + p2 * (b[spaceId2] - a[spaceId2]) + q * (c[spaceId2] - a[spaceId2]);
          s = (gamma - e1[spaceId2]) / (e2[spaceId2] - e1[spaceId2]);

          if(s < -epsilon || s > F1 + epsilon || p2 < -epsilon || q < -epsilon || (p2 + q) > F1 + epsilon) {
          } else {
            ////
            /*	if( edgeTriangleIntersection(elements[nodeList[n]].m_vertices[0],
                   elements[nodeList[n]].m_vertices[1],
                   elements[nodeList[n]].m_vertices[2], e1, e2) ){*/
            nodeList[noReallyIntersectingNodes] = nodeList[n];
            noReallyIntersectingNodes++;
            break;
          }
        }
      }
    }

    // write nodes in reallyIntersectingNodes
  }
  //  if(noNodes)
  if(noNodes && !noReallyIntersectingNodes) {
    //    m_log << "Rejected ! "<< ++rejectionCounter << endl;
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}


MInt Geometry3D::getSphereIntersectionMBElements(MFloat* P, MFloat radius, std::vector<MInt>& nodeList) {
  // TRACE();

  otherCalls++;

  MInt noReallyIntersectingNodes = 0;

  // compute minimum target region that contains the sphere:
  MFloat target[6] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                      numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                      numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
  MFloat enlargeFactor = 1.05;
  for(MInt i = 0; i < nDim; i++) {
    target[i] = P[i] - radius * enlargeFactor;
    target[i + nDim] = P[i] + radius * enlargeFactor;
  }

  // fetch all possibly intersecting triangles in this region:
  m_adt->retrieveNodesMBElements(target, nodeList);
  const MInt noNodes = nodeList.size();

  MFloat A[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
  MFloat B[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
  MFloat C[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
  MFloat AB[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                  numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
  MFloat BC[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                  numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
  MFloat CA[3] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                  numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
  MFloat V[3], Q1[3], Q2[3], Q3[3], QC[3], QA[3], QB[3];

  // loop over all elements
  for(MInt n = 0; n < noNodes; n++) {
    // store the vertices of the current element<3> in A, B, C:
    for(MInt i = 0; i < nDim; i++) {
      A[i] = mbelements[nodeList[n]].m_vertices[0][i];
      B[i] = mbelements[nodeList[n]].m_vertices[1][i];
      C[i] = mbelements[nodeList[n]].m_vertices[2][i];
    }

    // compute separation algorithm from http://realtimecollisiondetection.net/blog/?p=103:
    for(MInt i = 0; i < nDim; i++) {
      A[i] = A[i] - P[i];
      B[i] = B[i] - P[i];
      C[i] = C[i] - P[i];
      AB[i] = B[i] - A[i];
      BC[i] = C[i] - B[i];
      CA[i] = A[i] - C[i];
    }
    V[0] = -AB[1] * CA[2] + AB[2] * CA[1];
    V[1] = -AB[2] * CA[0] + AB[0] * CA[2];
    V[2] = -AB[0] * CA[1] + AB[1] * CA[0];
    MFloat rr = radius * radius;
    MFloat d = F0;
    MFloat e = F0;
    MFloat aa = F0;
    MFloat ab = F0;
    MFloat ac = F0;
    MFloat bb = F0;
    MFloat bc = F0;
    MFloat cc = F0;
    MFloat e1 = F0;
    MFloat e2 = F0;
    MFloat e3 = F0;
    for(MInt i = 0; i < nDim; i++) {
      d += A[i] * V[i];
      e += V[i] * V[i];
      aa += A[i] * A[i];
      ab += A[i] * B[i];
      ac += A[i] * C[i];
      bb += B[i] * B[i];
      bc += B[i] * C[i];
      cc += C[i] * C[i];
      e1 += AB[i] * AB[i];
      e2 += BC[i] * BC[i];
      e3 += CA[i] * CA[i];
    }
    MBool sep1 = d * d > rr * e;
    MBool sep2 = (aa > rr) && (ab > aa) && (ac > aa);
    MBool sep3 = (bb > rr) && (ab > bb) && (bc > bb);
    MBool sep4 = (cc > rr) && (ac > cc) && (bc > cc);
    MFloat d1 = ab - aa;
    MFloat d2 = bc - bb;
    MFloat d3 = ac - cc;
    MFloat qq1 = F0;
    MFloat qq2 = F0;
    MFloat qq3 = F0;
    MFloat qq1c = F0;
    MFloat qq2a = F0;
    MFloat qq3b = F0;
    for(MInt i = 0; i < nDim; i++) {
      Q1[i] = A[i] * e1 - d1 * AB[i];
      Q2[i] = B[i] * e2 - d2 * BC[i];
      Q3[i] = C[i] * e3 - d3 * CA[i];
      QC[i] = C[i] * e1 - Q1[i];
      QA[i] = A[i] * e2 - Q2[i];
      QB[i] = B[i] * e3 - Q3[i];
      qq1 += Q1[i] * Q1[i];
      qq2 += Q2[i] * Q2[i];
      qq3 += Q3[i] * Q3[i];
      qq1c += Q1[i] * QC[i];
      qq2a += Q2[i] * QA[i];
      qq3b += Q3[i] * QB[i];
    }
    MBool sep5 = (qq1 > rr * e1 * e1) && (qq1c > 0);
    MBool sep6 = (qq2 > rr * e2 * e2) && (qq2a > 0);
    MBool sep7 = (qq3 > rr * e3 * e3) && (qq3b > 0);

    MBool separated = sep1 || sep2 || sep3 || sep4 || sep5 || sep6 || sep7;

    if(!separated) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
    }
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}


MInt Geometry3D::getLineIntersectionMBElements(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  TRACE();

  otherCalls++;


  MInt noReallyIntersectingNodes = 0;
  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  MFloat e1[3], e2[3];
  for(MInt i = 0; i < nDim; i++) {
    e1[i] = targetRegion[i];
    e2[i] = targetRegion[i + nDim];
  }

  for(MInt n = 0; n < noNodes; n++) {
    edgeTICallCounter--;
    if(edgeTriangleIntersection(mbelements[nodeList[n]].m_vertices[0], mbelements[nodeList[n]].m_vertices[1],
                                mbelements[nodeList[n]].m_vertices[2], e1, e2)) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
    }
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}


MInt Geometry3D::getLineIntersectionMBElements2(MFloat* targetRegion, MInt* spaceDirection, std::vector<MInt>& nodeList,
                                                MInt bcId) {
  TRACE();

  otherCalls++;


  MBool singleIntersectionPoint;
  MInt noReallyIntersectingNodes = 0;
  m_adt->retrieveNodesMBElements(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  // TAN maxNoNodes is low
  MInt maxNoNodes = noNodes; // 50;
  MInt spaceId, spaceId1, spaceId2;
  MFloat* e1 = new MFloat[nDim];
  MFloat* e2 = new MFloat[nDim];
  MFloat epsilon = 0.00000000001;
  MFloat* a = new MFloat[nDim];
  MFloat* b = new MFloat[nDim];
  MFloat* c = new MFloat[nDim];
  MFloat gamma, s, p, q;
  MFloat* pP = new MFloat[nDim];
  MFloat** temp = new MFloat*[maxNoNodes];
  for(MInt k = 0; k < maxNoNodes; k++) {
    temp[k] = new MFloat[nDim];
  }

  for(MInt i = 0; i < nDim; i++) {
    e1[i] = targetRegion[i];
    e2[i] = targetRegion[i + nDim];
  }
  spaceId = spaceDirection[0];
  spaceId1 = spaceDirection[1];
  spaceId2 = spaceDirection[2];

  for(MInt n = 0; n < noNodes; n++) {
    if(mbelements[nodeList[n]].m_bndCndId != bcId) continue;
    for(MInt k = 0; k < nDim; k++) {
      a[k] = mbelements[nodeList[n]].m_vertices[0][k];
      b[k] = mbelements[nodeList[n]].m_vertices[1][k];
      c[k] = mbelements[nodeList[n]].m_vertices[2][k];
    }
    if(approx(a[spaceId1], b[spaceId1], MFloatEps) && !approx(a[spaceId1], c[spaceId1], MFloatEps)) {
      // aspaceId != bspaceId, otherwise a and b would be the same point
      q = (e1[spaceId1] - a[spaceId1]) / (c[spaceId1] - a[spaceId1]);
      p = (e1[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
    } else {
      if(!approx(a[spaceId1], b[spaceId1], MFloatEps) && approx(a[spaceId1], c[spaceId1], MFloatEps)) {
        // aspaceId != cspaceId, otherwise a and c would be the same point
        p = (e1[spaceId1] - a[spaceId1]) / (b[spaceId1] - a[spaceId1]);
        q = (e1[spaceId] - a[spaceId] - p * (b[spaceId] - a[spaceId])) / (c[spaceId] - a[spaceId]);
      } else {
        if(approx(a[spaceId], b[spaceId], MFloatEps) && !approx(a[spaceId], c[spaceId], MFloatEps)) {
          // aspaceId1 != bspaceId1, otherwise a and b would be the same point
          q = (e1[spaceId] - a[spaceId]) / (c[spaceId] - a[spaceId]);
          p = (e1[spaceId1] - a[spaceId1] - q * (c[spaceId1] - a[spaceId1])) / (b[spaceId1] - a[spaceId1]);
        } else {
          if(!approx(a[spaceId], b[spaceId], MFloatEps) && approx(a[spaceId], c[spaceId], MFloatEps)) {
            // aspaceId1 != cspaceId1, otherwise a and c would be the same point
            p = (e1[spaceId] - a[spaceId]) / (b[spaceId] - a[spaceId]);
            q = (e1[spaceId1] - a[spaceId1] - p * (b[spaceId1] - a[spaceId1])) / (c[spaceId1] - a[spaceId1]);
          } else {
            // aspaceId1 != bspaceId1 && aspaceId1 != cspaceId1 && aspaceId != bspaceId && aspaceId != cspaceId
            q = ((e1[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                 - (b[spaceId1] - a[spaceId1]) * (e1[spaceId] - a[spaceId]))
                / ((c[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                   - (b[spaceId1] - a[spaceId1]) * (c[spaceId] - a[spaceId]));
            p = (e1[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
          }
        }
      }
    }

    if(p * q >= 0 || p * q < 0) {
      // compute s
      gamma = a[spaceId2] + p * (b[spaceId2] - a[spaceId2]) + q * (c[spaceId2] - a[spaceId2]);
      s = (gamma - e1[spaceId2]) / (e2[spaceId2] - e1[spaceId2]);

      if(s < -epsilon || s > F1 + epsilon || p < -epsilon || q < -epsilon || (p + q) > F1 + epsilon) {
      } else {
        singleIntersectionPoint = true;

        // cut point pP
        for(MInt g = 0; g < nDim; g++) {
          pP[g] = e1[g] + s * (e2[g] - e1[g]);
          temp[noReallyIntersectingNodes][g] = pP[g];
        }

        if(noReallyIntersectingNodes == 0) {
          nodeList[noReallyIntersectingNodes] = nodeList[n];
          noReallyIntersectingNodes++;
        } else {
          for(MInt f = 0; f < noReallyIntersectingNodes; f++) {
            if((fabs(temp[f][0] - pP[0]) < epsilon) && (fabs(temp[f][1] - pP[1]) < epsilon)
               && (fabs(temp[f][2] - pP[2]) < epsilon)) {
              singleIntersectionPoint = false;
              break;
            }
          }
          if(singleIntersectionPoint) {
            nodeList[noReallyIntersectingNodes] = nodeList[n];
            noReallyIntersectingNodes++;
          }
        }
      }
    }
  }
  nodeList.resize(noReallyIntersectingNodes);

  // free memory
  for(MInt k = 0; k < maxNoNodes; k++) {
    delete[] temp[k];
    temp[k] = nullptr;
  }
  delete[] temp;
  temp = nullptr;

  delete[] e1;
  e1 = nullptr;
  delete[] e2;
  e2 = nullptr;
  delete[] a;
  a = nullptr;
  delete[] b;
  b = nullptr;
  delete[] c;
  c = nullptr;
  delete[] pP;
  pP = nullptr;


  return noReallyIntersectingNodes;
}

void Geometry3D::MoveAllMBElementVertex(MFloat* dx) {
  TRACE();

  otherCalls++;


  for(MInt e = 0; e < m_noMBElements; e++) {
    for(MInt j = 0; j < nDim; j++) {
      for(MInt i = 0; i < nDim; i++) {
        mbelements[e].m_vertices[j][i] += dx[i];
      }
    }
  }
}

void Geometry3D::MoveMBElementVertex(MInt e, MInt v, MFloat* dx) {
  TRACE();

  otherCalls++;


  for(MInt i = 0; i < nDim; i++) {
    mbelements[e].m_vertices[v][i] += dx[i];
  }
}

void Geometry3D::ReplaceMBElementVertex(MInt e, MInt v, MFloat* np) {
  TRACE();

  otherCalls++;


  for(MInt i = 0; i < nDim; i++) {
    mbelements[e].m_vertices[v][i] = np[i];
  }
}

void Geometry3D::UpdateMBNormalVector(MInt e) {
  TRACE();

  otherCalls++;


  MFloat v1[3], v2[3], n[3], norm_n = NAN;
  for(MInt i = 0; i < nDim; i++) {
    v1[i] = mbelements[e].m_vertices[1][i] - mbelements[e].m_vertices[0][i];
    v2[i] = mbelements[e].m_vertices[2][i] - mbelements[e].m_vertices[0][i];
  }
  n[0] = v1[1] * v2[2] - v2[1] * v1[2];
  n[1] = v1[2] * v2[0] - v2[2] * v1[0];
  n[2] = v1[0] * v2[1] - v2[0] * v1[1];

  norm_n = sqrt(POW2(n[0]) + POW2(n[1]) + POW2(n[2]));
  for(MInt i = 0; i < nDim; i++)
    mbelements[e].m_normal[i] = n[i] / norm_n;
}

void Geometry3D::UpdateMBBoundingBox() {
  TRACE();

  otherCalls++;


  for(MInt e = 0; e < m_mbelements->size(); e++) {
    mbelements[e].boundingBox();
  }
  for(MInt j = 0; j < nDim; j++) {
    m_mbminMax[j] = mbelements[0].m_minMax[j];
    m_mbminMax[j + nDim] = mbelements[0].m_minMax[j + nDim];
  }
  for(MInt i = 0; i < m_mbelements->size(); i++) {
    for(MInt j = 0; j < nDim; j++) {
      m_mbminMax[j + nDim] = (m_mbminMax[j + nDim] < mbelements[i].m_minMax[j + nDim])
                                 ? mbelements[i].m_minMax[j + nDim]
                                 : m_mbminMax[j + nDim];
      m_mbminMax[j] = (m_mbminMax[j] > mbelements[i].m_minMax[j]) ? mbelements[i].m_minMax[j] : m_mbminMax[j];
    }
  }
}

void Geometry3D::UpdateADT() {
  TRACE();

  otherCalls++;

  m_adt->buildTreeMB();
}

void Geometry3D::writeSTL(const char* fileName) {
  TRACE();

  ofstream ofl;

  ofl.open(fileName);
  ofl.precision(12);

  if(ofl) {
    ofl << "solid PRO2STL version 1.0 part " << domainId() << endl;

    for(MInt i = 0; i < m_noElements; i++) {
      ofl << "  facet normal";
      for(MInt j = 0; j < 3; j++) {
        ofl << " " << elements[i].m_normal[j];
      }
      ofl << endl << "  outer loop" << endl;

      for(MInt k = 0; k < 3; k++) {
        ofl << "   vertex";
        for(MInt l = 0; l < 3; l++) {
          ofl << " " << elements[i].m_vertices[k][l];
        }
        ofl << endl;
      }

      ofl << "  endloop" << endl << " endfacet" << endl;
    }
    ofl << "endsolid PRO2STL version 1.0 part " << domainId() << endl;
  }
}

void Geometry3D::writeSTLMB(const char* fileName, MInt& noNodes, MInt*& nodeList) {
  TRACE();

  ofstream ofl;

  ofl.open(fileName);
  ofl.precision(12);

  if(ofl) {
    ofl << "solid PRO2STL version 1.0 part " << domainId() << endl;

    for(MInt i = 0; i < noNodes; i++) {
      ofl << "  facet normal";
      for(MInt j = 0; j < 3; j++) {
        ofl << " " << mbelements[nodeList[i]].m_normal[j];
      }
      ofl << endl << "  outer loop" << endl;

      for(MInt k = 0; k < 3; k++) {
        ofl << "   vertex";
        for(MInt l = 0; l < 3; l++) {
          ofl << " " << mbelements[nodeList[i]].m_vertices[k][l];
        }
        ofl << endl;
      }

      ofl << "  endloop" << endl << " endfacet" << endl;
    }
    ofl << "endsolid PRO2STL version 1.0 part " << domainId() << endl;
  }
}

void Geometry3D::writeADTAndSTLToNetCDF(const char* fileName) {
  TRACE();
  using namespace maia::parallel_io;

  ParallelIo::size_type count = 0;

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
  ParallelIo::size_type tempDim[2];
  ParallelIo::size_type tempDim3D[3];

  tempDim[0] = m_noElements;
  tempDim[1] = nDim;

  tempDim3D[0] = m_noElements;
  tempDim3D[1] = 3;
  tempDim3D[2] = nDim;

  parallelIo.defineArray(PIO_FLOAT, "minMax", 2 * nDim);

  parallelIo.defineArray(PIO_FLOAT, "vertexNormal", 2, tempDim);

  parallelIo.defineArray(PIO_FLOAT, "vertexCoordinates", 3, tempDim3D);

  parallelIo.defineArray(PIO_INT, "boundaryConditionId", m_noElements);

  /*
   *
   * <--- ParallelIo define solver
   *
   */

  MFloatScratchSpace floatScratch(9 * m_noElements, AT_, "tmp_float_array");
  MFloat* tmp = floatScratch.getPointer();
  MIntScratchSpace intScratch(m_noElements, AT_, "tmp_int_array");
  MInt* tmpI = intScratch.getPointer();

  for(MInt i = 0; i < 2 * nDim; i++) {
    tmp[i] = m_minMax[i];
  }
  count = 2 * nDim;
  parallelIo.setOffset(count, 0);
  parallelIo.writeArray(&tmp[0], "minMax");

  for(MInt i = 0; i < m_noElements; i++) {
    for(MInt j = 0; j < 3; j++) {
      tmp[i * 3 + j] = elements[i].m_normal[j];
    }
  }
  count = m_noElements;
  parallelIo.setOffset(count, 0, 2);
  parallelIo.writeArray(&tmp[0], "vertexNormal");

  for(MInt i = 0; i < m_noElements; i++) {
    for(MInt j = 0; j < 3; j++) {
      for(MInt k = 0; k < 3; k++) {
        tmp[3 * (i * 3 + j) + k] = elements[i].m_vertices[j][k];
      }
    }
  }
  count = m_noElements;
  parallelIo.setOffset(count, 0, 3);
  parallelIo.writeArray(&tmp[0], "vertexCoordinates");

  for(MInt i = 0; i < m_noElements; i++) {
    tmpI[i] = elements[i].m_bndCndId;
  }
  count = m_noElements;
  parallelIo.setOffset(count, 0);
  parallelIo.writeArray(&tmpI[0], "boundaryConditionId");

  /* If the tree needs to be saved
    //  parIdVar = file.add_var("parentNodeId", ncInt, noNodeDim);
    //  for(MInt i = 0; i < m_noElements; i++){
    //        tmp[i] = m_adt->m_nodes[i].m_parent;
    //  }
    //  parIdVar->put(&tmp[0], m_noElements);

    //  childIdVar = file.add_var("childNodeIds", ncInt, noNodeDim, noChildDim);
    //    tmp[2*i]   = m_adt->m_nodes[i].m_leftSubtree;
    //    tmp[2*i+1] = m_adt->m_nodes[i].m_rightSubtree;
    //  }
    //  childIdVar->put(&tmp[0], m_noElements, 2);

    //  levelVar = file.add_var("level", ncInt, noNodeDim);
      //    level[i] = m_adt->m_nodes[i].m_depth;
    //  }
    //  levelVar->put(&tmp[0], m_noElements);

    //  vertsVar = file.add_var("vertex", ncInt, noNodeDim);
      //    tmp[i] = m_adt->m_nodes[i].m_element;
    //  }
    //  vertsVar->put(&tmp[0], m_noElements);
  */
  //  m_log << Scratch::printSelf();
}

void Geometry3D::readSTLNetCDF(const char* fileName) {
  TRACE();

  otherCalls++;

  m_log << "domainId() = " << domainId() << endl;
  m_log << "noDomains = " << noDomains() << endl;

  ParallelIo::size_type count = 0;

  // WARNING: untested switch from NetCDF/Parallel netCDF to ParallelIo
  // The method previously used direct I/O calls, which were replaced by
  // ParallelIo methods in summer 2015. However, since the method was not
  // used by any of the testcases, this code is still *untested*. Thus,
  // if your code uses this part of the code, please make sure that the
  // I/O still works as expected and then remove this warning as well as
  // the subsequent TERMM().
  TERMM(1, "untested I/O method, please see comment for how to proceed");
  ParallelIo parallelIo(fileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);

  ParallelIo::size_type noVerticesVarSize = 0;
  noVerticesVarSize = parallelIo.getArraySize("noVertices");

  MInt domainStartRead;
  MInt domainEndRead;
  MInt elementsPerDomain = noVerticesVarSize / noDomains() + 1;
  domainStartRead = domainId() * elementsPerDomain;

  if(domainId() == noDomains() - 1) {
    m_noElements = noVerticesVarSize - ((noDomains() - 1) * (elementsPerDomain));
    domainEndRead = domainStartRead + m_noElements - 1;
  } else {
    m_noElements = elementsPerDomain;
    domainEndRead = domainStartRead + m_noElements - 1;
  }

  m_log << "no elements in file: " << noVerticesVarSize << endl;
  m_log << "m_noElements = " << m_noElements << endl;
  m_log << "domainStartRead = " << domainStartRead << endl;
  m_log << "domainEndRead = " << domainEndRead << endl;

  if(m_elements != nullptr) delete m_elements;
  m_elements = new Collector<element<nDim>>(m_noElements, nDim, 0);
  elements = m_elements->a;

  MFloat* tmp = new MFloat[9 * m_noElements];
  MInt* tmpI = new MInt[m_noElements];

  count = m_noElements;
  parallelIo.setOffset(count, 0, 2);
  parallelIo.readArray(tmp, "vertexNormal");

  for(MInt i = 0; i < m_noElements; i++) {
    m_elements->append();
    for(MInt j = 0; j < 3; j++) {
      elements[i].m_normal[j] = tmp[i * 3 + j];
    }
  }

  count = m_noElements;
  parallelIo.setOffset(count, 0, 3);
  parallelIo.readArray(tmp, "vertexCoordinates");

  for(MInt i = 0; i < m_noElements; i++) {
    for(MInt j = 0; j < 3; j++) {
      for(MInt k = 0; k < 3; k++) {
        elements[i].m_vertices[j][k] = tmp[3 * (i * 3 + j) + k];
      }
    }
    elements[i].boundingBox();
  }

  count = m_noElements;
  parallelIo.setOffset(count, 0);
  parallelIo.readArray(tmpI, "boundaryConditionId");

  for(MInt i = 0; i < m_noElements; i++) {
    for(MInt j = 0; j < 3; j++) {
      elements[i].m_bndCndId = tmpI[i];
    }
  }

  count = 2 * nDim;
  parallelIo.setOffset(count, 0);
  parallelIo.readArray(tmp, "minMax");
  for(MInt i = 0; i < 2 * nDim; i++) {
    m_minMax[i] = tmp[i];
  }

  /*
    for(MInt j=0; j < nDim; j++){
      m_minMax[j] = elements[0].m_minMax[j];
      m_minMax[j + nDim] = elements[0].m_minMax[j + nDim];
    }
    for(MInt i=0; i < m_noElements; i++){
      for(MInt j=0; j < nDim; j++){
        m_minMax[j + nDim] = (m_minMax[j + nDim] <
              elements[i].m_minMax[j + nDim]) ? elements[i].m_minMax[j + nDim] : m_minMax[j + nDim];
        m_minMax[j] = (m_minMax[j] > elements[i].m_minMax[j]) ? elements[i].m_minMax[j] : m_minMax[j];
      }
    }
  */

  delete[] tmp;
  delete[] tmpI;
}

/** \brief Returns unique edges of a given set segment id
 *
 * \author Andreas Lintermann
 * \date 17.09.2015
 *
 * param[in] segmentId the id of the segment (only for serial geometry)
 * return the edges of the segment
 *
 **/
inline vector<pair<MFloat*, MFloat*>> Geometry3D::GetUniqueSegmentEdges(MInt segmentId) {
  TRACE();

  vector<pair<MFloat*, MFloat*>> edges;

  for(MInt i = m_segmentOffsetsWithoutMB[segmentId]; i < m_segmentOffsetsWithoutMB[segmentId + 1]; i++) {
    // this runs over all edges of the triangle
    for(MInt e = 0; e < nDim; e++) {
      MFloat* p1 = elements[i].m_vertices[e];
      MFloat* p2 = elements[i].m_vertices[(e + 1) % nDim];

      // if not in yet, then add
      MInt del;
      MBool in = isEdgeAlreadyInCollection(edges, p1, p2, &del);
      if(!in)
        edges.emplace_back(p1, p2);
      else {
        auto it = edges.begin() + del;
        edges.erase(it);
      }
    }
  }

  return edges;
}


/** \brief Returns unique edges of a given set segment id for parallel geometry
 *
 * \author Andreas Lintermann
 * \date 17.09.2015
 *
 * The vertices are proviced by an array (see parameters)
 *
 * param[in] tri_vx the triangles in one array (v00,v01,v02,v10,v11,v12,v20,v21,v22,...)
 * param[in] keepOffsets only consider the triangles at the given offsets
 * param[in] size the number of triangles in the input array
 * return the edges of the segment
 *
 **/
inline vector<pair<MFloat*, MFloat*>> Geometry3D::GetUniqueSegmentEdgesParGeom(MFloat* tri_vx, MInt* keepOffsets,
                                                                               MInt size) {
  TRACE();

  vector<pair<MFloat*, MFloat*>> edges;

  for(MInt i = 0; i < size; i++) {
    MInt off = keepOffsets[i];

    // this runs over all edges of the triangle
    for(MInt k = 0; k < nDim; k++) {
      MFloat* p1 = &tri_vx[off + k * nDim];
      MFloat* p2 = &tri_vx[off + (((k + 1) * nDim) % (nDim * nDim))];


      MInt del;
      MBool in = isEdgeAlreadyInCollection(edges, p1, p2, &del);
      if(!in)
        edges.emplace_back(p1, p2);
      else {
        auto it = edges.begin() + del;
        edges.erase(it);
      }
    }
  }

  return edges;
}

/** \brief Checks if an edge given by two points is in a vector
 *
 * \author Andreas Lintermann
 * \date 17.09.2015
 *
 * \param[in] edges the vector containing the collection of edges
 * \param[in] p1 point 1 of the edge
 * \param[in] p2 point 2 of the edge
 * return if the edge is already contained
 *
 **/
inline MBool Geometry3D::isEdgeAlreadyInCollection(vector<pair<MFloat*, MFloat*>> edges, MFloat* p1, MFloat* p2,
                                                   MInt* to_delete) {
  TRACE();

  MInt pos = 0;
  for(auto& edge : edges) {
    MBool pt1_in = true;
    MBool pt2_in = true;

    MFloat* pcol1 = edge.first;
    MFloat* pcol2 = edge.second;

    for(MInt d = 0; d < nDim; d++)
      pt1_in = pt1_in && approx(p1[d], pcol1[d], MFloatEps);

    if(!pt1_in) {
      pt1_in = true;
      for(MInt d = 0; d < nDim; d++)
        pt1_in = pt1_in && approx(p1[d], pcol2[d], MFloatEps);
    }

    for(MInt d = 0; d < nDim; d++)
      pt2_in = pt2_in && approx(p2[d], pcol1[d], MFloatEps);

    if(!pt2_in) {
      pt2_in = true;
      for(MInt d = 0; d < nDim; d++)
        pt2_in = pt2_in && approx(p2[d], pcol2[d], MFloatEps);
    }

    if(pt1_in && pt2_in) {
      *to_delete = pos;
      return true;
    }
    pos++;
  }
  return false;
}


/** \brief This function gets all boundary vertices of an element<3> in a circular order
 *
 * \author Andreas Lintermann
 * \date 02.02.2010
 *
 * The algorithm is devided into two parts. First, all border edges are extracted by running over
 * all triangles of the segment and checking which edge appears only in one triangle. Then, the
 * vertices are extracted in a cirular order by running over all extracted edges and doing a walkthrough
 * element<3> by element<3>.
 *
 * param[in] segmentId the id of the segment (only for serial geometry)
 * param[in] tri_vx the triangles in one array (v00,v01,v02,v10,v11,v12,v20,v21,v22,...) (only for parallel geometry)
 * param[in] keepOffsets only consider the triangles at the given offsets (only for parallel geometry)
 * param[in] num the size of the first dimension of the returned 2D array
 * return the circular ordered edges of the segment
 *
 **/
MFloat** Geometry3D::GetBoundaryVertices(MInt segmentId, MFloat* tri_vx, MInt* keepOffsets, MInt size, MInt* num) {
  TRACE();

  vector<pair<MFloat*, MFloat*>> tmp_edges;
  // this runs over all triangles that have a certain segment id

  if(m_parallelGeometry)
    tmp_edges = GetUniqueSegmentEdgesParGeom(tri_vx, keepOffsets, size);
  else
    tmp_edges = GetUniqueSegmentEdges(segmentId);

  if(tmp_edges.size() < 2) {
    stringstream errorMsg;
    errorMsg << "ERROR: No boundary edges found for segment " << segmentId << endl;
    m_log << errorMsg.str();
    mTerm(1, AT_, errorMsg.str());
  }

  // Get distinct elements in a circular ordering
  vector<MFloat*> tmp_points;
  MBoolScratchSpace visited(tmp_edges.size(), AT_, "visited");
  for(MInt i = 0; i < (signed)tmp_edges.size(); i++)
    visited[i] = false;

  // insert first edge
  tmp_points.push_back(tmp_edges[0].first);
  tmp_points.push_back(tmp_edges[0].second);
  visited[0] = true;

  MFloat* first = tmp_points[0];
  MInt cmpPos = 1;

  for(MInt i = 1; i < (signed)tmp_edges.size(); i++) {
    for(MInt j = 1; j < (signed)tmp_edges.size(); j++) {
      // skip if alrady inserted
      if(visited[j]) continue;

      pair<MFloat*, MFloat*> edge = tmp_edges[j];

      // insert in this order
      if(vectorsEqual(edge.first, tmp_points[cmpPos])) {
        tmp_points.push_back(edge.second);
        visited[j] = true;
        cmpPos++;
      } else if(vectorsEqual(edge.second, tmp_points[cmpPos])) {
        tmp_points.push_back(edge.first);
        visited[j] = true;
        cmpPos++;
      }
    }
  }

  if(!vectorsEqual(tmp_points[tmp_points.size() - 1], first)) {
    stringstream errorMsg;
    errorMsg << "ERROR: Inconsistency in segment " << segmentId << endl;
    m_log << errorMsg.str();
    mTerm(1, AT_, errorMsg.str());
  }

  // Fill the array
  MFloat** vertices = new MFloat*[tmp_points.size() - 1];
  for(MInt i = 0; i < (MInt)tmp_points.size() - 1; i++) {
    vertices[i] = new MFloat[nDim];
    for(MInt j = 0; j < nDim; j++)
      vertices[i][j] = tmp_points[i][j];
  }

  *num = tmp_points.size() - 1;

  return vertices;
}

/** \brief Returns the area of a segment
 *
 * \author Andreas Lintermann
 * \date 17.09.2015
 *
 * \param[in] segmentId the id of the segment
 * return the area of the segment
 *
 **/
MFloat Geometry3D::GetBoundarySize(MInt segmentId) {
  TRACE();

  MInt offsetStart = 0;
  MInt offsetEnd = 0;
  if(m_parallelGeometry) {
    offsetStart = m_segmentOffsets[segmentId];
    offsetEnd = m_segmentOffsets[segmentId + 1];
  } else {
    offsetStart = m_segmentOffsetsWithoutMB[segmentId];
    offsetEnd = m_segmentOffsetsWithoutMB[segmentId + 1];
  }

  if(offsetStart == offsetEnd && !m_parallelGeometry)
    mTerm(1, AT_, "ERROR: You cannot choose a MB segment to calculate the characteristic Length");


  MFloat area = 0;
  std::array<MFloat, nDim> cross;
  std::array<MFloat, nDim> edge1;
  std::array<MFloat, nDim> edge2;

  for(MInt i = m_segmentOffsets[segmentId]; i < m_segmentOffsets[segmentId + 1]; i++) {
    for(MInt j = 0; j < nDim; j++) {
      edge1[j] = elements[i].m_vertices[1][j] - elements[i].m_vertices[0][j];
      edge2[j] = elements[i].m_vertices[2][j] - elements[i].m_vertices[0][j];
    }

    cross[0] = edge1[1] * edge2[2] - edge2[1] * edge1[2];
    cross[1] = edge1[2] * edge2[0] - edge2[2] * edge1[0];
    cross[2] = edge1[0] * edge2[1] - edge2[0] * edge1[1];

    area += sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]) / 2;
  }


  return area;
}

/** \brief Returns the area of a segment
 *
 * \author Andreas Lintermann
 * \date 17.09.2015
 *
 * \param[in] vertices the vertices as an array (v00,v01,v02,v10,v11,v12,v20,v21,v22,...)
 * \param[in] size the size of the array
 * return the area of the segment
 *
 **/
MFloat Geometry3D::GetBoundarySize(MFloat* vertices, MInt* keepOffset, MInt size) {
  TRACE();

  MFloat area = 0;
  std::array<MFloat, nDim> cross;
  std::array<MFloat, nDim> edge1;
  std::array<MFloat, nDim> edge2;

  for(MInt i = 0; i < size; i++) {
    MInt off = i;

    if(keepOffset != nullptr) off = keepOffset[i];

    for(MInt j = 0; j < nDim; j++) {
      edge1[j] = vertices[off + nDim + j] - vertices[off + j];
      edge2[j] = vertices[off + 2 * nDim + j] - vertices[off + j];
    }

    cross[0] = edge1[1] * edge2[2] - edge2[1] * edge1[2];
    cross[1] = edge1[2] * edge2[0] - edge2[2] * edge1[0];
    cross[2] = edge1[0] * edge2[1] - edge2[0] * edge1[1];

    area += sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]) / 2;
  }


  return area;
}


void Geometry3D::logStatistics() {
  TRACE();

  m_log << "******************** Geometry3D statistics ********************" << endl;
  m_log << "  Calls to getIntersectionElements: " << getIECallCounter << endl;
  m_log << "  Comm from getIntersectionElements: " << getIECommCounter << endl;
  m_log << "  Calls to getIntersectionElementsTetraeder: " << getIETCallCounter << endl;
  m_log << "  Calls to getLineIntersectionElements2: " << m_getLIE2CallCounter << endl;
  m_log << "  Comm from getLineIntersectionElements2: " << getLIE2CommCounter << endl;
  m_log << "  Calls to edgeTriangleIntersection: " << edgeTICallCounter << endl;
  m_log << "  Other calls in geometry3d: " << otherCalls << endl;
  m_log << "***************************************************************" << endl;
}

/** \brief Determines the ownership of a segment
 *
 * \author Andreas Lintermann
 * \date 18.09.2015
 *
 * \param[in] segmentId the segment id
 * \param[in] own pointer to the result of this domain owns the segment
 * \param[in] sumowners pointer to the result of the sum of the owners
 * \param[in] firstOwner pointer to the result of which is the first owner in the communicator
 * \param[in] owners array holding information which domain holds the segment
 *
 **/
void Geometry3D::determineSegmentOwnership(MInt segmentId, MInt* own, MInt* sumowners, MInt* firstOwner, MInt* owners) {
  TRACE();

  *firstOwner = -1;
  *sumowners = 0;

  if(m_ownSegmentId[segmentId])
    *own = 1;
  else
    *own = 0;


  MPI_Allgather(own, 1, MPI_INT, owners, 1, MPI_INT, mpiComm(), AT_, "own", "owners");
  for(MInt d = 0; d < noDomains(); d++) {
    if(owners[d] > 0) {
      if(*firstOwner < 0) *firstOwner = d;
      *sumowners += owners[d];
    }
  }
}

/** \brief This function gets the maximal radius for a boundary segment
 *
 * \author Andreas Lintermann
 * \date 21.01.2013
 *
 * This function first determines the area for all projections XY, XZ, YZ according to
 *
 * \f[A_{xy} = \frac{1}{2}\sum_{i=0}^{n-1}\left(x_i y_{i+1} - x_{i+1} y_i\right)\f]
 * \f[A_{xz} = \frac{1}{2}\sum_{i=0}^{n-1}\left(x_i z_{i+1} - x_{i+1} z_i\right)\f]
 * \f[A_{xy} = \frac{1}{2}\sum_{i=0}^{n-1}\left(y_i z_{i+1} - y_{i+1} z_i\right)\f]
 *
 * Based on the obtained size of the areas with the consideration of the special
 * cases that some areas can be zero the center of gravity is calculated using:;;;
 *
 * \f[ c_{\alpha}=\frac{1}{6A_{\alpha\beta}}\sum_{i=0}^{n-1}\left(\alpha_i + \alpha_{i+1}\right)\left(\alpha_i
 *\beta_{i+1}-\alpha_{i+1} \beta_i\right)\f],
 *
 * where \f$A_{ab}\neq 0\f$ and \f$\alpha,\beta\in\{x,y,z\}\f$.
 *
 * \param[in] vertices the vertices as pointer to a pointer
 * \param[in] num the number of vertices
 *
 **/
MFloat Geometry3D::getBndMaxRadius(MFloat** vertices, MInt num) {
  TRACE();

  MFloat areaXY = 0.0;
  MFloat areaXZ = 0.0;
  MFloat areaYZ = 0.0;
  MFloat cog_tmp[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  MFloat cog[3] = {0.0, 0.0, 0.0};

  for(MInt i = 0; i < num; i++) {
    MInt next = (i + 1) % num;
    MFloat xy = vertices[i][0] * vertices[next][1] - vertices[next][0] * vertices[i][1];
    MFloat xz = vertices[i][0] * vertices[next][2] - vertices[next][0] * vertices[i][2];
    MFloat yz = vertices[i][1] * vertices[next][2] - vertices[next][1] * vertices[i][2];

    areaXY += xy;
    areaXZ += xz;
    areaYZ += yz;

    cog_tmp[0][0] += (vertices[i][0] + vertices[next][0]) * xy;
    cog_tmp[0][1] += (vertices[i][0] + vertices[next][0]) * xz;
    cog_tmp[1][0] += (vertices[i][1] + vertices[next][1]) * xy;
    cog_tmp[1][1] += (vertices[i][1] + vertices[next][1]) * yz;
    cog_tmp[2][0] += (vertices[i][2] + vertices[next][2]) * xz;
    cog_tmp[2][1] += (vertices[i][2] + vertices[next][2]) * yz;
  }

  areaXY = fabs(areaXY) * 0.5;
  areaXZ = fabs(areaXZ) * 0.5;
  areaYZ = fabs(areaYZ) * 0.5;

  MFloat eps = 0.00000001;
  if(areaXY < eps) {
    if(areaYZ < eps) {
      cog[0] = 1.0 / (6.0 * areaXZ) * cog_tmp[0][1];
      cog[1] = vertices[0][1];
      cog[2] = 1.0 / (6.0 * areaXZ) * cog_tmp[2][0];
    } else if(areaXZ < eps) {
      cog[0] = vertices[0][0];
      cog[1] = 1.0 / (6.0 * areaYZ) * cog_tmp[1][1];
      cog[2] = 1.0 / (6.0 * areaYZ) * cog_tmp[2][1];
    } else {
      cog[0] = 1.0 / (6.0 * areaXZ) * cog_tmp[0][1];
      cog[1] = 1.0 / (6.0 * areaYZ) * cog_tmp[1][1];
      cog[2] = 1.0 / (6.0 * areaYZ) * cog_tmp[2][1];
    }
  } else if(areaXZ < eps) {
    if(areaYZ < eps) {
      cog[0] = 1.0 / (6.0 * areaXY) * cog_tmp[0][0];
      cog[1] = 1.0 / (6.0 * areaXY) * cog_tmp[1][0];
      cog[2] = vertices[0][2];
    } else {
      cog[0] = 1.0 / (6.0 * areaXY) * cog_tmp[0][0];
      cog[1] = 1.0 / (6.0 * areaYZ) * cog_tmp[1][1];
      cog[2] = 1.0 / (6.0 * areaYZ) * cog_tmp[2][1];
    }
  } else if(areaYZ < eps) {
    cog[0] = 1.0 / (6.0 * areaXZ) * cog_tmp[0][1];
    cog[1] = 1.0 / (6.0 * areaXY) * cog_tmp[1][0];
    cog[2] = 1.0 / (6.0 * areaXZ) * cog_tmp[2][0];
  } else {
    cog[0] = 1.0 / (6.0 * areaXZ) * cog_tmp[0][1];
    cog[1] = 1.0 / (6.0 * areaYZ) * cog_tmp[1][1];
    cog[2] = 1.0 / (6.0 * areaYZ) * cog_tmp[2][1];
  }

  m_log << "      * cog:                         " << cog[0] << " " << cog[1] << " " << cog[2] << endl;

  MFloat maxRadius = 0.0;

  for(MInt i = 0; i < num; i++) {
    MFloat tmp = (vertices[i][0] - cog[0]) * (vertices[i][0] - cog[0])
                 + (vertices[i][1] - cog[1]) * (vertices[i][1] - cog[1])
                 + (vertices[i][2] - cog[2]) * (vertices[i][2] - cog[2]);

    if(tmp > maxRadius) maxRadius = tmp;
  }

  return sqrt(maxRadius);
}

/** \brief deletes the current element<3> collector and reinitializes
 *
 * \author Andreas Lintermann
 * \date 11.01.2016
 *
 * \param[in] new_size the new size of the collector
 *
 **/
void Geometry3D::resizeCollector(MInt new_size) {
  TRACE();

  m_log << "      * resizing collector: " << m_noElements << " - " << m_noElements + new_size << endl;

  auto* tmp = new Collector<element<nDim>>(m_noElements + new_size, nDim, 0);

  element<nDim>* fromPtr = m_elements->a;
  element<nDim>* toPtr = tmp->a;

  for(MInt e = 0; e < m_noElements; e++) {
    tmp->append();
    copyElement(e, e, fromPtr, toPtr);
  }

  delete m_elements;
  m_elements = tmp;
  elements = m_elements->a;
}

/** \brief Adds an element<3> to the collector
 *
 * \author Andreas Lintermann
 * \date 25.09.2015
 *
 * \param[in] tri conatins the following information in MFloat: originalId, segmentId, bndCndId, normal, vertices
 *
 **/
void Geometry3D::addElement(MFloat* tri) {
  // TRACE();

  m_elements->append();
  m_noElements++;

  element<3>* tris = m_elements->a;

  MInt segmentId = tri[1];

  if(segmentId + 1 != m_noSegments)
    for(MInt c = m_noSegments; c > segmentId; c--) {
      copyElement(m_segmentOffsets[c - 1], m_segmentOffsets[c]);
    }

  for(MInt c = m_noSegments; c > segmentId; c--) {
    m_segmentOffsets[c] += 1;
  }

  // insert
  MInt pos = m_segmentOffsets[segmentId + 1] - 1;
  tris[pos].m_originalId = (MInt)tri[0];
  tris[pos].m_segmentId = segmentId;
  tris[pos].m_bndCndId = (MInt)tri[2];

  MInt j = 3;
  for(MInt d = 0; d < nDim; d++)
    tris[pos].m_normal[d] = tri[j++];
  for(MInt d1 = 0; d1 < nDim; d1++)
    for(MInt d2 = 0; d2 < nDim; d2++)
      tris[pos].m_vertices[d1][d2] = tri[j++];

  tris[pos].boundingBox();
}

/** \brief Copies an element<3> from one poistion to another
 *
 * \author Andreas Lintermann
 * \date 25.09.2015
 *
 * \param[in] from the location to copy the element<3> from
 * \param[in] to the location to copy the element<3> to
 *
 **/
void Geometry3D::copyElement(MInt from, MInt to) {
  // TRACE();
  element<3>* tris = m_elements->a;

  for(MInt d = 0; d < nDim; d++)
    tris[to].m_normal[d] = tris[from].m_normal[d];

  for(MInt d1 = 0; d1 < nDim; d1++)
    for(MInt d2 = 0; d2 < nDim; d2++)
      tris[to].m_vertices[d1][d2] = tris[from].m_vertices[d1][d2];

  for(MInt d = 0; d < 2 * nDim; d++)
    tris[to].m_minMax[d] = tris[from].m_minMax[d];

  tris[to].m_segmentId = tris[from].m_segmentId;
  tris[to].m_bndCndId = tris[from].m_bndCndId;
  tris[to].m_originalId = tris[from].m_originalId;
}

/** \brief Copies an element<3> from one poistion to another
 *
 * \author Andreas Lintermann
 * \date 25.09.2015
 *
 * \param[in] from the location to copy the element<3> from
 * \param[in] to the location to copy the element<3> to
 * \param[in] fromPtr the pointer to the from array
 * \param[in] toPtr the pointer to the to array
 *
 **/
void Geometry3D::copyElement(MInt from, MInt to, element<3>* fromPtr, element<3>* toPtr) {
  // TRACE();
  for(MInt d = 0; d < nDim; d++)
    toPtr[to].m_normal[d] = fromPtr[from].m_normal[d];

  for(MInt d1 = 0; d1 < nDim; d1++)
    for(MInt d2 = 0; d2 < nDim; d2++)
      toPtr[to].m_vertices[d1][d2] = fromPtr[from].m_vertices[d1][d2];

  for(MInt d = 0; d < 2 * nDim; d++)
    toPtr[to].m_minMax[d] = fromPtr[from].m_minMax[d];

  toPtr[to].m_segmentId = fromPtr[from].m_segmentId;
  toPtr[to].m_bndCndId = fromPtr[from].m_bndCndId;
  toPtr[to].m_originalId = fromPtr[from].m_originalId;
}
