// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometryadt.h"

#include <algorithm>
#include "MEMORY/collector.h"
#include "UTIL/debug.h"
#include "UTIL/timer.h"
#include "geometry.h"
#include "geometryadtnode.h"

using namespace std;

using namespace sortFunctions;

template <MInt nDim>
GeometryAdt<nDim>::GeometryAdt(const Geometry<nDim>* geometry) {
  TRACE();

  m_geometry = geometry;
  m_nodes = nullptr;

  // moving boundary
  m_mbnodes = nullptr;
}

template <MInt nDim>
GeometryAdt<nDim>::~GeometryAdt() {
  TRACE();

  delete[] m_nodes;

  // moving boundary
  delete[] m_mbnodes;
}

/** \brief Create an ADT from the geometry
 *
 *
 */
template <MInt nDim>
void GeometryAdt<nDim>::buildTree() {
  buildTreeNew();
}

template <MInt nDim>
MLong GeometryAdt<nDim>::memoryUsage() {
  return m_geometry->m_elements->size() * (5 * sizeof(MInt) + 2 * sizeof(MFloat)); // m_nodes[0].getStaticElementSize();
}

template <MInt nDim>
void GeometryAdt<nDim>::buildTreeOld() {
  TRACE();

  MInt noElements = m_geometry->m_elements->size();
  // If m_nodes already allocted delete it
  if(m_nodes) {
    delete[] m_nodes;
  }
  // Create nodes
  m_nodes = new GeometryAdtNode[mMax(noElements, 0)];
  // Holds the start id for elements at node
  vector<MInt> from(noElements, -1);
  // Holds the end id for elements at node
  vector<MInt> to(noElements, -1);

  MInt currentDir = 0;
  MInt nodeCounter = 0;
  MInt currentNode = 0;
  MInt currentDepth = 0;
  // Stores the position in the idArray for the last node of the left subtree
  MInt partitionId;
  // Stores the current left and right subtree nodes
  MInt left, right;
  // Stores the id of the element connected to the current node
  MInt currentElement;
  // Stores the nodes that will be processed later

  stack<MInt> nodeStack;

  // Initialize elements member for sorting function objects
  element<nDim>* elements = m_geometry->m_elements->a;

  // Set minimum and maximum values to geometry values (bounding box)
  copy_n(&m_geometry->m_minMax[0], 2 * nDim, &m_minMax[0]);

  m_log << " Building tree the old way .." << endl << endl;
  // Create Array with element ids, (on which the actual sorting
  // will be performed)
  vector<MInt>::iterator it, itFrom, itTo;
  vector<MInt> idArray(std::max(noElements, 1));
  MInt i = 0;
  for(it = idArray.begin(); it != idArray.end(); it++) {
    *it = i++;
  }


  // Set root node
  m_root = 0;
  nodeCounter = 0;
  currentNode = m_root;

  from.at(m_root) = idArray.front();
  to.at(m_root) = idArray.back();
  m_nodes[m_root].m_parent = -1;

  // For all nodes...
  for(;;) {
    // Set subtree ids
    /// If not a leaf...
    m_nodes[currentNode].m_depth = currentDepth;
    currentDir = currentDepth % (nDim * 2);
    // If more than one element left...
    if(from[currentNode] != to[currentNode]) {
      //       m_log << " -- " << " -- " << endl;

      left = ++nodeCounter;
      m_nodes[currentNode].m_leftSubtree = left;
      m_nodes[left].m_parent = currentNode;

      // If only two elements left, there is no right subtree
      if(to[currentNode] - from[currentNode] == 1) {
        right = 0;
        m_nodes[currentNode].m_rightSubtree = right;
      } else {
        right = ++nodeCounter;
        m_nodes[currentNode].m_rightSubtree = right;
        m_nodes[right].m_parent = currentNode;
      }
      // Set from and to iterators for sorting
      itFrom = idArray.begin();
      itFrom += from[currentNode];
      itTo = idArray.begin();
      // to-iterator must point behind last element!
      itTo += to[currentNode] + 1;

      // (sort parent ids)
      sort(itFrom, itTo, lessMinMax<nDim>(currentDir, elements));

      // mTerm(1, AT_);
      // Connect id to current node and increase from[currentNode]
      // (for that the element won't be connected twice)
      currentElement = idArray[from[currentNode]++];
      m_nodes[currentNode].m_element = currentElement;


      // find partition
      partitionId = (MInt)((to[currentNode] - from[currentNode]) / 2 + from[currentNode]);
      m_nodes[currentNode].m_partition = elements[idArray[partitionId]].m_minMax[currentDir];

      // Set min and max values
      m_nodes[currentNode].m_a = elements[idArray[from[currentNode] - 1]].m_minMax[currentDir];
      m_nodes[currentNode].m_b = elements[idArray[to[currentNode]]].m_minMax[currentDir];

      // Set from and to for sub-trees

      // Left subtree (lower values)
      from[left] = from[currentNode];
      to[left] = partitionId;

      // Right subtree (higher values)
      if(right) {
        from[right] = to[left] + 1;
        to[right] = to[currentNode];
      }
      // if only one element left ...
    } else {
      m_nodes[currentNode].m_leftSubtree = 0;
      m_nodes[currentNode].m_rightSubtree = 0;
      // Connect id to leaf
      currentElement = idArray[to[currentNode]];
      m_nodes[currentNode].m_element = currentElement;
      m_nodes[currentNode].m_a = elements[from[currentNode] - 1].m_minMax[currentDir];
      m_nodes[currentNode].m_b = elements[to[currentNode]].m_minMax[currentDir];
    }
    // Set current dir

    //   writeNode(m_root);

    // If left subtree exists ...
    if(m_nodes[currentNode].m_leftSubtree) {
      // increase depth and put rightSubtree Id on stack
      currentDepth++;
      if(m_nodes[currentNode].m_rightSubtree) nodeStack.push(nodeCounter);
      // ...go left,
      currentNode = m_nodes[currentNode].m_leftSubtree;

    } else {
      // else pop id from stack and go to id.
      if(!nodeStack.empty()) {
        currentNode = nodeStack.top();
        nodeStack.pop();
        currentDepth = m_nodes[m_nodes[currentNode].m_parent].m_depth + 1;
      } else {
        // Finished when empty stack and no left subtree.
        m_log << "ok" << endl;
        return;
      }
    }
  }
}

template <MInt nDim>
void GeometryAdt<nDim>::buildTreeNew() {
  TRACE();

  //  NEW_TIMER(timer1, (MChar*)"buildTree");

  //  START_TIMER (timer1);

  MInt noElements = m_geometry->m_elements->size();

  if(noElements == 1)
    mTerm(1,
          AT_,
          "Apparently this method can not build an adt with just one element. I observed that in 2D but did not "
          "debug further.");

  // If m_nodes already allocted delete it
  if(m_nodes) {
    delete[] m_nodes;
  }
  // Create nodes

  m_nodes = new GeometryAdtNode[mMax(noElements, 0)];
  // Holds the start id for elements at node
  MInt* from;
  from = new MInt[mMax(noElements, 0)];

  // Holds the end id for elements at node
  MInt* to;
  to = new MInt[mMax(noElements, 0)];

  MInt currentDir = 0;
  MInt nodeCounter = 0;
  MInt currentNode = 0;
  MInt currentDepth = 0;
  // Stores the position in the idArray for the last node of the left subtree
  MInt partitionId;
  // Stores the current left and right subtree nodes
  MInt left, right;
  // Stores the id of the element connected to the current node
  MInt currentElement;
  // Stores the nodes that will be processed later

  stack<MInt> nodeStack;

  // Initialize elements member for sorting function objects
  element<nDim>* elements = m_geometry->m_elements->a;

  for(MInt i = 0; i < 2 * nDim; i++)
    m_minMax[i] = m_geometry->m_minMax[i];


  m_log << "  + Building normal tree with " << noElements << " elements..." << endl;
  // Create Array with element ids, (on which the actual sorting
  // will be performed)
  vector<MInt>::iterator it, itFrom, itTo;
  vector<MInt> idArray(noElements);
  MInt i = 0;
  for(it = idArray.begin(); it != idArray.end(); it++)
    *it = i++;


  m_log << "    - setting root node" << endl;
  // Set root node
  m_root = 0;
  currentDir = 0;
  nodeCounter = 0;
  currentNode = m_root;

  from[m_root] = *idArray.begin();
  to[m_root] = *(idArray.end() - 1);
  m_nodes[m_root].m_parent = -1;
  m_nodes[m_root].m_depth = 0;

  left = 0;
  right = 0;

  if(noElements == 1) {
    m_nodes[m_root].m_leftSubtree = left;
    m_nodes[m_root].m_rightSubtree = right;

    // Connect id to leaf
    currentElement = idArray[to[m_root]];
    m_nodes[m_root].m_element = currentElement;
    m_nodes[m_root].m_a = elements[0].m_minMax[currentDir];
    m_nodes[m_root].m_b = elements[0].m_minMax[currentDir];

    m_log << "    - done building the tree" << endl;
    return;
  }


  // For all nodes...
  m_log << "    - inserting other nodes" << endl;
  nodeStack.push(m_root);
  currentDepth = 0;
  while(!nodeStack.empty()) {
    currentNode = nodeStack.top();
    nodeStack.pop();
    currentDepth = m_nodes[currentNode].m_depth;
    currentDir = currentDepth % (nDim * 2);

    // If more than one element left...
    if(from[currentNode] != to[currentNode]) {
      left = ++nodeCounter;
      m_nodes[currentNode].m_leftSubtree = left;
      m_nodes[left].m_parent = currentNode;
      m_nodes[left].m_depth = currentDepth + 1;

      // If only two elements left, there is no right subtree
      if(to[currentNode] - from[currentNode] == 1) {
        right = 0;
        m_nodes[currentNode].m_rightSubtree = right;
      } else {
        right = ++nodeCounter;
        m_nodes[currentNode].m_rightSubtree = right;
        m_nodes[right].m_parent = currentNode;
        m_nodes[right].m_depth = currentDepth + 1;
        nodeStack.push(right);
      }
      nodeStack.push(left);

      // Set from and to iterators for sorting
      itFrom = idArray.begin();
      itFrom += from[currentNode];
      itTo = idArray.begin();
      // to-iterator must point behind last element!
      itTo += to[currentNode] + 1;

      // (sort parent ids)
      sort(itFrom, itTo, lessMinMax<nDim>(currentDir, elements));

      // Connect id to current node and increase from[currentNode]
      // (for that the element won't be connected twice)
      currentElement = idArray[from[currentNode]++];
      m_nodes[currentNode].m_element = currentElement;

      // find partition (This is actually a KDTree!!!)
      partitionId = (MInt)((to[currentNode] - from[currentNode]) / 2 + from[currentNode]);
      m_nodes[currentNode].m_partition = elements[idArray[partitionId]].m_minMax[currentDir];

      // Set min and max values
      m_nodes[currentNode].m_a = elements[idArray[from[currentNode] - 1]].m_minMax[currentDir];
      m_nodes[currentNode].m_b = elements[idArray[to[currentNode]]].m_minMax[currentDir];

      // Set from and to for sub-trees

      // Left subtree (lower values)
      from[left] = from[currentNode];
      to[left] = partitionId;

      // Right subtree (higher values)
      if(right) {
        from[right] = to[left] + 1;
        to[right] = to[currentNode];
      }
      // if only one element left ...
    } else {
      right = 0;
      left = 0;
      m_nodes[currentNode].m_leftSubtree = left;
      m_nodes[currentNode].m_rightSubtree = right;

      // Connect id to leaf
      currentElement = idArray[to[currentNode]];
      m_nodes[currentNode].m_element = currentElement;
      m_nodes[currentNode].m_a = elements[from[currentNode] - 1].m_minMax[currentDir];
      m_nodes[currentNode].m_b = elements[to[currentNode]].m_minMax[currentDir];
    }
    // Set current dir

    // writeNode(m_root);
  }

  m_log << "    - done building the tree" << endl;
  // for(MInt k = 0; k < 10; k++) writeNode(k);
  delete[] from;
  delete[] to;
}

template <MInt nDim>
void GeometryAdt<nDim>::writeNode(MInt node) {
  m_log << " - Node : " << node << endl;
  m_log << "parent      : " << m_nodes[node].m_parent << endl;
  m_log << "partition   : " << m_nodes[node].m_partition << endl;
  m_log << "leftSubtree : " << m_nodes[node].m_leftSubtree << endl;
  m_log << "rightSubtree: " << m_nodes[node].m_rightSubtree << endl;
  m_log << "depth       : " << m_nodes[node].m_depth << endl;
  m_log << "a           : " << m_nodes[node].m_a << endl;
  m_log << "b           : " << m_nodes[node].m_b << endl;
}

/** \brief retrieves the nodes that intersect with a given target bounding box
 *
 * \author ???, Andreas Lintermann
 * \date 04.05.2016
 *
 * Traverses the tree and checks for an overlap between of the bounding box
 * of the given tree element and a specified target region.
 *
 * \param[in] targetRegion the target region defined by a bounding box
 * \param[in] noNodes will contain the number of intersected elements after the function call
 * \param[in] nodeList will contain the list of elements that are intersected after the function call
 **/
// TODO labels:GEOM,DLB optimize further? is a major time consumer during FV initSolutionStep during DLB
template <MInt nDim>
void GeometryAdt<nDim>::retrieveNodes(MFloat* targetRegion, std::vector<MInt>& nodeList) const {
  // TRACE();

  MInt noElements = m_geometry->m_elements->size();
  if(noElements == 0) {
    return;
  }

  // 2-D Mapping
  ASSERT(nodeList.empty(), "NodeList is not empty!");
  MFloat t_min[2 * MAX_SPACE_DIMENSIONS];
  MFloat t_max[2 * MAX_SPACE_DIMENSIONS];

  // t_min: min(bounding box: x, y, z), min(target: x, y, z)
  // t_max: max(target: x, y, z), max(bounding box: x, y, z)
  for(MInt i = 0; i < nDim; i++) {
    t_min[i] = m_minMax[i];
    t_min[i + nDim] = targetRegion[i];
    t_max[i] = targetRegion[i + nDim];
    t_max[i + nDim] = m_minMax[i + nDim];
  }

  // Init empty stack and start at first node
  MInt root = 0;
  stack<MInt> subtreeStack;

  // Infinite loop until complete tree is traversed
  for(;;) {
    // Check if the current nodes element-bounding box intersects target region
    const MInt currentElementId = m_nodes[root].m_element;
    const MFloat* const minMax = m_geometry->m_elements->a[currentElementId].m_minMax;
    MBool doesOverlap = true;

    for(MInt i = 0; i < 2 * nDim; i++)
      if(minMax[i] > t_max[i] || minMax[i] < t_min[i]) {
        doesOverlap = false;
        break;
      }

    // Inside target domain => add current element to return list
    if(doesOverlap) {
      nodeList.push_back(currentElementId);
    }

    const MInt currentDir = m_nodes[root].m_depth % (2 * nDim);

    // if right subtree is inside target domain push on stack
    const MInt right = m_nodes[root].m_rightSubtree;
    if(right > 0 && m_nodes[root].m_b >= t_min[currentDir] && m_nodes[root].m_partition <= t_max[currentDir])
      subtreeStack.push(right);

    // if left subtree is inside target domain set it as root
    const MInt left = m_nodes[root].m_leftSubtree;
    if(left > 0 && m_nodes[root].m_partition >= t_min[currentDir] && m_nodes[root].m_a <= t_max[currentDir]) {
      root = left;
      continue;
    }

    if(subtreeStack.empty()) {
      return;
    } else {
      root = subtreeStack.top();
      subtreeStack.pop();
    }
  }
}


template <MInt nDim>
void GeometryAdt<nDim>::writeTreeToDx() {
  TRACE();

  // Write dx output (intersected triangles)
  stringstream filename;
  // Write dx output (all triangles)
  filename.str("");
  filename << "dxOut_all";
  filename << ".dx";
  ofstream ofl(filename.str().c_str());
  //  ofl.open(filename.str().c_str());
  ofl << "object 1 class array type float rank 1 shape 3 items " << m_geometry->m_noElements * nDim << " data follows"
      << endl;
  for(MInt i = 0; i < m_geometry->m_noElements; i++) {
    for(MInt j = 0; j < nDim; j++) {
      for(MInt k = 0; k < nDim; k++) {
        ofl << m_geometry->m_elements->a[i].m_vertices[j][k] << " ";
      }
      ofl << endl;
    }
  }
  ofl << "object 2 class array type int rank 1 shape 3 items " << m_geometry->m_noElements << " data follows" << endl;
  MInt c = 0;
  for(MInt i = 0; i < m_geometry->m_noElements; i++) {
    for(MInt j = 0; j < nDim; j++) {
      ofl << c++ << " ";
    }
    ofl << endl;
  }
  ofl << R"(attribute "element type" string "triangles")" << endl;
  ofl << R"(attribute "ref" string "positions")" << endl;

  ofl << "object 3 class array type float rank 0 items " << m_geometry->m_noElements * 3 << " data follows" << endl;
  c = 0;
  for(MInt i = 0; i < m_geometry->m_noElements; i++) {
    for(MInt j = 0; j < nDim; j++) {
      ofl << c << " ";
    }
    c++;
    ofl << endl;
  }
  ofl << R"(attribute "dep" string "positions")" << endl;

  ofl << "object \"irregular positions irregular connections\" class field" << endl;
  ofl << "component \"positions\" value 1" << endl;
  ofl << "component \"connections\" value 2" << endl;
  ofl << "component \"data\" value 3" << endl;
  ofl << "" << endl;
  ofl.close();
}


template <MInt nDim>
void GeometryAdt<nDim>::splitTree(MInt /*noSubTrees*/) {
  /*
  geometryPropertyMap* pMap;
  string fileName = "NCTEST";

  GeometryIONetcdf* geometryIOBase = new GeometryIONetcdf;
  //  geometryIOBase->writeProperties((MChar*)fileName.c_str(), pMap);




  NcFile file(fileName.c_str(), NcFile::Replace);
  NcVar* tmpVar;
  NcDim* tmpDim;
  MString dimName;

  dimName = "nDim";
  tmpDim = file.add_dim(dimName.c_str(), nDim);

  dimName = "noPointsPerVertex";
  tmpDim = file.add_dim(dimName.c_str(), 3);

  dimName = "noVertices";
  tmpDim = file.add_dim(dimName.c_str(), 3);

  dimName = "noNodes";
  tmpDim = file.add_dim(dimName.c_str(), 3);

  dimName = "noChildIds";
  tmpDim = file.add_dim(dimName.c_str(), 3);
  */

  /*
  for( propertyMap::iterator i = pMap->begin(); i!=pMap->end(); i++) {
    dimName = i->first;
    switch( i->second->type() )
      {
      case INT:
  {
    if ( i->second->count() > 1 ){
      dimName.append("Dim");
      tmpDim = file.add_dim( dimName.c_str(), i->second->count());
      tmpVar = file.add_var(i->first.c_str(), ncInt, tmpDim);
      tmpVar->put(i->second->intField, i->second->count());
    }
    else {
      tmpVar = file.add_var(i->first.c_str(), ncInt, 0);
      tmpVar->put(i->second->intField, 1);
    }
    break;
  }

      case FLOAT:
  {
    if ( i->second->count() > 1 ){
      dimName.append("Dim");
      tmpDim = file.add_dim( dimName.c_str(), i->second->count());
      tmpVar = file.add_var(i->first.c_str(), ncDouble, tmpDim);
      tmpVar->put(i->second->floatField, i->second->count());
    }
    else {
      tmpVar = file.add_var(i->first.c_str(), ncDouble, 0);
      tmpVar->put(i->second->floatField, 1);
    }
    break;
  }

      case STRING:
  {
    if ( i->second->stringField[0].length() > 1 ){
      dimName.append("Dim");
      tmpDim = file.add_dim( dimName.c_str(), i->second->stringField[0].length());
      tmpVar = file.add_var(i->first.c_str(), ncChar, tmpDim);
      tmpVar->put(i->second->stringField[0].c_str(), i->second->stringField[0].length());
    }
    else {
      tmpVar = file.add_var(i->first.c_str(), ncChar, 0);
      tmpVar->put(i->second->stringField[0].c_str(), 1);
    }
    break;
  }

      default:
  {
    break;
  }
      }
  }
  */
}


template <MInt nDim>
void GeometryAdt<nDim>::buildTreeMB() {
  TRACE();
  /*
    //TAN FIX
  //  NEW_TIMER(timer1, (MChar*)"buildTreeMB");
  //  START_TIMER (timer1);
    MInt  noMBElements = m_geometry->m_mbelements->size();
    // If m_mbnodes already allocted delete it
    if(m_mbnodes){
      delete  [] m_mbnodes;
    }
    // Create nodes
    m_mbnodes = new GeometryAdtNode[noMBElements];
    // Holds the start id for elements at node
    MInt *from;
    from = new MInt[noMBElements];
    // Holds the end id for elements at node
    MInt *to;
    to = new MInt[noMBElements];
    MInt currentDir = 0;
    MInt nodeCounter = 0;
    MInt currentNode = 0;
    MInt currentDepth = 0;
    // Stores the position in the idArray for the last node of the left subtree
    MInt partitionId;
    // Stores the current left and right subtree nodes
    MInt left, right;
    // Stores the id of the element connected to the current node
    MInt currentElement;
    // Stores the nodes that will be processed later

    stack <MInt> nodeStack;

    // Initialize elements member for sorting function objects
    element* mbelements = m_geometry->m_mbelements->a;
    lessMinMax<nDim>::m_mbelements = mbelements;

    if(m_mbminMax) {delete [] m_mbminMax;}
    m_mbminMax = new MFloat[2*nDim];
    for(MInt i=0; i < 2*nDim; i ++){
      m_mbminMax[i] = m_geometry->m_mbminMax[i];
    }

    m_log << " Building tree...";
    // Create Array with element ids, (on which the actual sorting
    // will be performed)
    vector <MInt>::iterator it, itFrom, itTo;
    vector <MInt> idArray (noMBElements);
    MInt i=0;
    for (it = idArray.begin(); it!=idArray.end(); it++){
      *it=i++;
    }

    // Set root node
    m_root = 0;
    currentDir = 0;
    nodeCounter = 0;
    currentNode = m_root;

    from[m_root] = *idArray.begin();
    to[m_root] = *(idArray.end()-1);
    m_mbnodes[m_root].m_parent=-1;
    left = 0;
    right = 0;
    // For all nodes...

    for(;;){
    // Set subtree ids
      /// If not a leaf...
      m_mbnodes[currentNode].m_depth = currentDepth;
      currentDir = currentDepth % (nDim * 2);
      lessMinMax<nDim>::m_dir = currentDir;
      // If more than one element left...
      if(from[currentNode] != to[currentNode]){
  //       m_log << " -- " << " -- " << endl;

        left=++nodeCounter;
        m_mbnodes[currentNode].m_leftSubtree=left;
        m_mbnodes[ left ].m_parent = currentNode;

        // If only two elements left, there is no right subtree
        if(to[currentNode] - from[currentNode] == 1 ){
    right=0;
    m_mbnodes[currentNode].m_rightSubtree = right;
        }else{
    right=++nodeCounter;
    m_mbnodes[currentNode].m_rightSubtree = right;
    m_mbnodes[right].m_parent = currentNode;
        }
        // Set from and to iterators for sorting
        itFrom = idArray.begin();
        itFrom += from[currentNode];
        itTo = idArray.begin();
        // to-iterator must point behind last element!
        itTo += to[currentNode] + 1;

        // (sort parent ids)
        sort(itFrom, itTo,lessMinMax<nDim>::mbcomp);

        // mTerm(1, AT_);
        // Connect id to current node and increase from[currentNode]
        // (for that the element won't be connected twice)
        currentElement = idArray[from[currentNode]++];
        m_mbnodes[currentNode].m_element = currentElement;

        // find partition
        partitionId = (MInt)((to[currentNode] - from[currentNode]) / 2 + from[currentNode] );
        m_mbnodes[currentNode].m_partition = mbelements[idArray[partitionId]].m_minMax[currentDir];

        // Set min and max values
        m_mbnodes[currentNode].m_a = mbelements[idArray[from[currentNode] - 1]].m_minMax[currentDir];
        m_mbnodes[currentNode].m_b = mbelements[idArray[to[currentNode]]].m_minMax[currentDir];

        // Set from and to for sub-trees

        // Left subtree (lower values)
        from[left] = from[currentNode];
        to[left] =  partitionId;

        // Right subtree (higher values)
        if(right){
    from[right] = to[left] + 1;
    to[right] = to[currentNode];
        }
      // if only one element left ...
      }else{
        m_mbnodes[currentNode].m_leftSubtree=0;
        m_mbnodes[currentNode].m_rightSubtree=0;
        // Connect id to leaf
        currentElement = idArray[to[currentNode]];
        m_mbnodes[currentNode].m_element = currentElement;
        m_mbnodes[currentNode].m_a = mbelements[from[currentNode] - 1].m_minMax[currentDir];
        m_mbnodes[currentNode].m_b = mbelements[to[currentNode]].m_minMax[currentDir];

      }
      // Set current dir

      //   writeNode(m_root);

      // If left subtree exists ...
      if( m_mbnodes[currentNode].m_leftSubtree ){
        // increase depth and put rightSubtree Id on stack
        currentDepth++;
        if(m_mbnodes[currentNode].m_rightSubtree)
    nodeStack.push(nodeCounter);
        // ...go left,
        currentNode = m_mbnodes[currentNode].m_leftSubtree;

      }else{
        // else pop id from stack and go to id.
        if(!nodeStack.empty()){
    currentNode = nodeStack.top();
    nodeStack.pop();
    currentDepth = m_mbnodes[m_mbnodes[currentNode].m_parent].m_depth + 1;
        }else{
    // Finished when empty stack and no left subtree.
    m_log << "ok" << endl;
    //TAN FIXED
    delete []from;
    delete []to;
    return;
        }
      }
    }*/

  m_log << "geometry m_elements size = ";
  m_log << m_geometry->m_elements->size() << endl;

  MInt noMBElements = m_geometry->m_mbelements->size();
  // If m_mbnodes already allocted delete it

  if(m_mbnodes) {
    delete[] m_mbnodes;
  }
  // Create nodes

  m_mbnodes = new GeometryAdtNode[mMax(noMBElements, 0)];
  // Holds the start id for elements at node
  vector<MInt> from(noMBElements, -1);
  // Holds the end id for elements at node
  vector<MInt> to(noMBElements, -1);

  MInt currentDir = 0;
  MInt nodeCounter = 0;
  MInt currentNode = 0;
  MInt currentDepth = 0;
  // Stores the position in the idArray for the last node of the left subtree
  MInt partitionId;
  // Stores the current left and right subtree nodes
  MInt left, right;
  // Stores the id of the element connected to the current node
  MInt currentElement;
  // Stores the nodes that will be processed later

  stack<MInt> nodeStack;

  // Initialize elements member for sorting function objects
  element<nDim>* elements = m_geometry->m_mbelements->a;

  for(MInt i = 0; i < 2 * nDim; i++) {
    m_mbminMax[i] = m_geometry->m_mbminMax[i];
  }

  m_log << "  + Building MB tree with " << noMBElements << " elements..." << endl << endl;
  // Create Array with element ids, (on which the actual sorting
  // will be performed)
  vector<MInt>::iterator it, itFrom, itTo;
  vector<MInt> idArray(std::max(noMBElements, 1));
  MInt i = 0;
  for(it = idArray.begin(); it != idArray.end(); it++) {
    *it = i++;
  }

  // Set root node
  m_root = 0;
  currentDir = 0;
  nodeCounter = 0;
  currentNode = m_root;

  from.at(m_root) = idArray.front();
  to.at(m_root) = idArray.back();
  m_mbnodes[m_root].m_parent = -1;
  m_mbnodes[m_root].m_depth = 0;

  left = 0;
  right = 0;
  // For all nodes...

  nodeStack.push(m_root);
  currentDepth = 0;

  while(!nodeStack.empty()) {
    currentNode = nodeStack.top();
    nodeStack.pop();
    currentDepth = m_mbnodes[currentNode].m_depth;
    currentDir = currentDepth % (nDim * 2);

    // If more than one element left...
    if(from[currentNode] != to[currentNode]) {
      left = ++nodeCounter;
      m_mbnodes[currentNode].m_leftSubtree = left;
      m_mbnodes[left].m_parent = currentNode;
      m_mbnodes[left].m_depth = currentDepth + 1;

      // If only two elements left, there is no right subtree
      if(to[currentNode] - from[currentNode] == 1) {
        right = 0;
        m_mbnodes[currentNode].m_rightSubtree = right;
      } else {
        right = ++nodeCounter;
        m_mbnodes[currentNode].m_rightSubtree = right;
        m_mbnodes[right].m_parent = currentNode;
        m_mbnodes[right].m_depth = currentDepth + 1;
        nodeStack.push(right);
      }
      nodeStack.push(left);

      // Set from and to iterators for sorting
      itFrom = idArray.begin();
      itFrom += from[currentNode];
      itTo = idArray.begin();
      // to-iterator must point behind last element!
      itTo += to[currentNode] + 1;

      // (sort parent ids)
      sort(itFrom, itTo, lessMinMax<nDim>(currentDir, elements));

      // Connect id to current node and increase from[currentNode]
      // (for that the element won't be connected twice)
      currentElement = idArray[from[currentNode]++];
      m_mbnodes[currentNode].m_element = currentElement;


      // find partition (This is actually a KDTree!!!)
      partitionId = (MInt)((to[currentNode] - from[currentNode]) / 2 + from[currentNode]);
      m_mbnodes[currentNode].m_partition = elements[idArray[partitionId]].m_minMax[currentDir];

      // Set min and max values
      m_mbnodes[currentNode].m_a = elements[idArray[from[currentNode] - 1]].m_minMax[currentDir];
      m_mbnodes[currentNode].m_b = elements[idArray[to[currentNode]]].m_minMax[currentDir];

      // Set from and to for sub-trees

      // Left subtree (lower values)
      from[left] = from[currentNode];
      to[left] = partitionId;

      // Right subtree (higher values)
      if(right) {
        from[right] = to[left] + 1;
        to[right] = to[currentNode];
      }
      // if only one element left ...
    } else {
      right = 0;
      left = 0;
      m_mbnodes[currentNode].m_leftSubtree = left;
      m_mbnodes[currentNode].m_rightSubtree = right;
      // Connect id to leaf
      currentElement = idArray[to[currentNode]];
      m_mbnodes[currentNode].m_element = currentElement;
      m_mbnodes[currentNode].m_a = elements[from[currentNode] - 1].m_minMax[currentDir];
      m_mbnodes[currentNode].m_b = elements[to[currentNode]].m_minMax[currentDir];
    }
  }

  m_log << "ok" << endl;
}

template <MInt nDim>
void GeometryAdt<nDim>::retrieveNodesMBElements(const MFloat* targetRegion, std::vector<MInt>& nodeList) const {
  //   TRACE();

  // 2-D Mapping
  MFloat* t_min;
  MFloat* t_max;
  MFloat* x;
  MBool doesOverlap;
  ASSERT(nodeList.empty(), "NodeList is not empty!");

  t_min = new MFloat[2 * nDim];
  t_max = new MFloat[2 * nDim];
  x = new MFloat[2 * nDim];
  for(MInt i = 0; i < nDim; i++) {
    t_min[i] = m_mbminMax[i];
    t_min[i + nDim] = targetRegion[i];
    t_max[i] = targetRegion[i + nDim];
    t_max[i + nDim] = m_mbminMax[i + nDim];
  }

  MInt currentDir;
  GeometryElement<nDim>* currentElement;
  GeometryElement<nDim>* mbelements = m_geometry->m_mbelements->a;
  MInt currentElementId;
  MInt root;
  MInt left, right;
  root = 0;
  MInt searchedNodes = 0;
  stack<MInt> subtreeStack;
  for(;;) {
    currentDir = m_mbnodes[root].m_depth % (2 * nDim);
    currentElementId = m_mbnodes[root].m_element;
    currentElement = &mbelements[currentElementId];

    //     m_log << " Searching node " << root << endl;
    searchedNodes++;
    // Check if the current nodes element-bounding box intersects the target region
    doesOverlap = true;
    for(MInt i = 0; i < 2 * nDim; i++) {
      if(currentElement->m_minMax[i] <= t_max[i] && currentElement->m_minMax[i] >= t_min[i]) {
      } else {
        doesOverlap = false;
        break;
      }
    }

    if(doesOverlap) {
      // Inside target domain => add current element to return list
      nodeList.push_back(currentElementId);
      //       m_log << " Found intersection candidate " << currentElementId << endl;
    }
    // if right subtree is inside target domain push on stack
    right = m_mbnodes[root].m_rightSubtree;
    if(right) {
      if(m_mbnodes[root].m_b >= t_min[currentDir] && m_mbnodes[root].m_partition <= t_max[currentDir])
        subtreeStack.push(right);
    }
    // if left subtree is inside target domain set it as root
    left = m_mbnodes[root].m_leftSubtree;
    if(left) {
      if(m_mbnodes[root].m_partition >= t_min[currentDir] && m_mbnodes[root].m_a <= t_max[currentDir]) {
        root = left;
        continue;
      }
    }
    if(subtreeStack.empty()) {
      delete[] x;
      delete[] t_min;
      delete[] t_max;

      return;
    } else {
      root = subtreeStack.top();
      subtreeStack.pop();
    }
  }
}

template <MInt nDim>
MInt GeometryAdt<nDim>::get_root() {
  return m_root;
}

// Explicit instantiations for 2D and 3D
template class GeometryAdt<2>;
template class GeometryAdt<3>;
