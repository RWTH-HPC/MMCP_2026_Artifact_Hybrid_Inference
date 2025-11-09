// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <queue>
#include "INCLUDE/maiatypes.h"
#include "IO/infoout.h"
#include "globalvariables.h"

template <typename T>
class BinaryTreeNode {
 public:
  BinaryTreeNode<T>() {
    m_object = nullptr;
    m_father = nullptr;
    m_leftChild = nullptr;
    m_rightChild = nullptr;
  }

  explicit BinaryTreeNode<T>(T* obj) {
    m_object = obj;
    m_father = nullptr;
    m_leftChild = nullptr;
    m_rightChild = nullptr;
  }

  BinaryTreeNode<T>(T* obj, BinaryTreeNode<T>* father, BinaryTreeNode<T>* leftChild, BinaryTreeNode<T>* rightChild) {
    m_object = obj;
    m_father = father;
    m_leftChild = leftChild;
    m_rightChild = rightChild;
  }


  T* getObject() { return m_object; }
  BinaryTreeNode<T>* getFather() { return m_father; }
  BinaryTreeNode<T>* getLeftChild() { return m_leftChild; }
  BinaryTreeNode<T>* getRightChild() { return m_rightChild; }

  void setObject(T* obj) { m_object = obj; }
  void setFather(BinaryTreeNode<T>* node) { m_father = node; }
  void setLeftChild(BinaryTreeNode<T>* node) { m_leftChild = node; }
  void setRightChild(BinaryTreeNode<T>* node) { m_rightChild = node; }

 protected:
  T* m_object;
  BinaryTreeNode<T>* m_father;
  BinaryTreeNode<T>* m_leftChild;
  BinaryTreeNode<T>* m_rightChild;
};

using BinaryTreeIntNode = BinaryTreeNode<MInt>;

template <typename T>
class BinaryTree {
 public:
  BinaryTree<T>() { m_root = nullptr; }
  explicit BinaryTree<T>(BinaryTreeNode<T>* root) { m_root = root; }

  BinaryTreeNode<T>* getRoot() { return m_root; }

  MInt getMaxDepth() {
    if(m_root == nullptr) {
      return 0;
    }

    MInt lDepth = getNodeDepth(m_root->getLeftChild());
    MInt rDepth = getNodeDepth(m_root->getRightChild());

    if(lDepth > rDepth) {
      return (lDepth + 1);
    }
    return (rDepth + 1);
  }

 protected:
  BinaryTreeNode<T>* m_root;

  MInt getNodeDepth(BinaryTreeNode<T>* node) {
    if(node == nullptr) {
      return (0);
    }

    MInt lDepth = getNodeDepth(m_root->getLeftChild());
    MInt rDepth = getNodeDepth(m_root->getRightChild());

    if(lDepth > rDepth) {
      return (lDepth + 1);
    }
    return (rDepth + 1);
  }
};

using BinaryIntTree = BinaryTree<MInt>;


class BinaryPropertyMPITree : private BinaryIntTree {
 public:
  explicit BinaryPropertyMPITree(MInt num_procs) : m_myMPILocation(nullptr) {
    MInt rank = globalDomainId();
    m_queue = new std::queue<BinaryTreeIntNode*>;
    for(MInt i = 0; i < num_procs; i++) {
      auto* ins_int = new MInt();
      *ins_int = i;
      auto* ins = new BinaryTreeIntNode(ins_int);
      if(i == 0) {
        m_root = ins;
        m_queue->push(m_root);
      } else {
        insertBalancedNode(ins);
      }

      if(i == rank) {
        m_myMPILocation = ins;
      }
    }
  }
  ~BinaryPropertyMPITree() {
    for(MUlong i = 0; i < m_queue->size(); i++) {
      delete m_queue->front();
      m_queue->pop();
    }
    delete m_queue;
    m_queue = nullptr;
  }

  void printTree() { printTreeRoot(m_root, 0); }

  BinaryTreeIntNode* getMyMPILocation() {
    if(m_myMPILocation == nullptr) {
      return nullptr;
    }
    return m_myMPILocation;
  }

  BinaryTreeIntNode* getMyMPIReceiver() {
    if(m_myMPILocation != nullptr && m_myMPILocation->getFather() != nullptr) {
      return m_myMPILocation->getFather();
    }
    return nullptr;
  }

  BinaryTreeIntNode* getMyLeftMPISender() {
    if(m_myMPILocation != nullptr && m_myMPILocation->getLeftChild() != nullptr) {
      return m_myMPILocation->getLeftChild();
    }
    return nullptr;
  }
  BinaryTreeIntNode* getMyRightMPISender() {
    if(m_myMPILocation != nullptr && m_myMPILocation->getRightChild() != nullptr) {
      return m_myMPILocation->getRightChild();
    }
    return nullptr;
  }

 protected:
  BinaryTreeIntNode* m_myMPILocation;
  std::queue<BinaryTreeIntNode*>* m_queue;

  static void printTreeRoot(BinaryTreeIntNode* node, MInt level) {
    m_log << "L " << level << ":" << *(node->getObject()) << std::endl;
    if(node->getLeftChild() != nullptr) {
      printTreeRoot(node->getLeftChild(), level + 1);
    }
    if(node->getRightChild() != nullptr) {
      printTreeRoot(node->getRightChild(), level + 1);
    }
  }

  void insertBalancedNode(BinaryTreeIntNode* node) {
    while(!m_queue->empty()) {
      // full, descend
      if(m_queue->front()->getLeftChild() != nullptr && m_queue->front()->getRightChild() != nullptr) {
        m_queue->push(m_queue->front()->getLeftChild());
        m_queue->push(m_queue->front()->getRightChild());
        m_queue->pop();
      }
      // not full, fill up
      else {
        node->setFather(m_queue->front());
        if(m_queue->front()->getLeftChild() == nullptr) {
          m_queue->front()->setLeftChild(node);
        } else {
          m_queue->front()->setRightChild(node);
        }

        break;
      }
    }
  }
};

#endif // BINARYTREE_H
