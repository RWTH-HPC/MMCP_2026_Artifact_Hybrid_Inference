// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "INCLUDE/maiatypes.h"

class Application;
class Scratch;

/*!\class Environment
   \date begin: 03.01.20 change 00.00.00

   \brief Environment for the program

*/
class Environment {
 public:
  Environment(int, char**);
  ~Environment();

  MInt run();
  MInt end();
  static MInt m_argc;
  static MChar** m_argv;


 private:
  Application* mApplication;
  Scratch* mScratch;
  void parseCommandline();
  static void printStartupInformation();

  MString m_propertyFileInput = "properties.toml";
};

#endif
