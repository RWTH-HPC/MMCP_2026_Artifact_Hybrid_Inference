// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LOGTABLE_H
#define LOGTABLE_H

#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "INCLUDE/maiatypes.h"

namespace maia::logtable {

const MInt width = 80;

struct FrameEntry {
  const MString title;

  FrameEntry(MString title_) : title(title_){};
  virtual ~FrameEntry() = default;

  virtual MString buildString(const int shift = 0) = 0;
};

struct Data final : FrameEntry {
  const MString data;

  const MInt alignment;

  enum { LEFT, RIGHT };

  std::vector<std::shared_ptr<FrameEntry>> entries;

  // String version
  Data(MString title_, MString data_) : FrameEntry(title_), data(data_), alignment(LEFT){};

  // Integer version
  Data(MString title_, MInt data_) : FrameEntry(title_), data(std::to_string(data_)), alignment(RIGHT){};

  // Long version
  Data(MString title_, MLong data_) : FrameEntry(title_), data(std::to_string(data_)), alignment(RIGHT){};

  // Bool version
  Data(MString title_, MBool data_) : FrameEntry(title_), data(data_ ? "yes" : "no"), alignment(RIGHT){};

  // Floating point version
  Data(MString title_, MFloat data_) : FrameEntry(title_), data(std::to_string(data_)), alignment(RIGHT){};

  MString buildString(const int shift = 0) final {
    std::stringstream ss;

    ss << "| " << std::setw(width / 2 + shift) << std::left << title;
    if(alignment == LEFT) {
      ss << "| " << std::setw(width / 2 - 6) << std::left << data;
    } else if(alignment == RIGHT) {
      ss << "| " << std::setw(width / 2 - 6) << std::right << data;
    }
    ss << " |" << '\n';

    for(auto&& entry : entries) {
      ss << "| " << entry->buildString(shift - 2);
    }

    return ss.str();
  }

  template <typename T>
  Data& addData(const std::string& title_, T data_) {
    std::shared_ptr<Data> d = std::make_shared<Data>(title_, data_);
    entries.push_back(d);
    return *d;
  }

  void addBlank() {
    std::shared_ptr<Data> d = std::make_shared<Data>("", std::string(""));
    entries.push_back(d);
  }

  ~Data() final {
    for(auto&& en : entries) {
      en.reset();
    }
  }
};

struct Group : FrameEntry {
  std::vector<std::shared_ptr<FrameEntry>> entries;

  Group(MString title_) : FrameEntry(title_){};

  MString buildString(const int shift = 0) final {
    std::stringstream ss;

    ss << title << '\n';

    for(auto&& entry : entries) {
      ss << entry->buildString(shift);
    }

    return ss.str();
  }

  template <typename T>
  Data& addData(const std::string& title_, T data_) {
    std::shared_ptr<Data> d = std::make_shared<Data>(title_, data_);
    entries.push_back(d);
    return *d;
  }

  void addBlank() {
    std::shared_ptr<Data> d = std::make_shared<Data>("", std::string(""));
    entries.push_back(d);
  }
};

struct Frame {
  const MString bar = [&]() {
    std::stringstream tmp;
    for(int i = 0; i < width; i++) {
      tmp << "-";
    }
    return tmp.str();
  }();

  const MString title;

  std::vector<std::shared_ptr<FrameEntry>> entries;

  Frame(MString title_) : title(title_){};

  MString buildString() {
    std::stringstream ss;

    ss << bar << '\n';
    ss << title << '\n';
    ss << bar << '\n';

    MBool first = true;
    for(auto&& entry : entries) {
      if(!first) {
        ss << '\n';
      }
      ss << entry->buildString();
      first = false;
    }

    ss << bar << '\n';

    return ss.str();
  }

  Group& addGroup(const std::string& title_) {
    std::shared_ptr<Group> g = std::make_shared<Group>(title_);
    entries.push_back(g);
    return *g;
  }
};

} // namespace maia::logtable

#endif
