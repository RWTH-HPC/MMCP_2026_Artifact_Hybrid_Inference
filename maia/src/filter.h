// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FILTER_H_
#define FILTER_H_

#include <cmath>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "UTIL/debug.h"

/// \file This file contains various filter functions.

namespace maia {
namespace filter {

namespace slope {

namespace detail_ {

/// Auxiliary function to create two-way slope filters.
template <typename T, typename F>
T f2way(const T a, const T b, const T c, const T d, const T x, F&& f) {
  if(x < c) {
    return f(a, b, x);
  } else {
    return 1.0 - f(c, d, x);
  }
}

/// Auxiliary function to create a box slope filter (constant slope width)
template <MInt nDim, typename T, typename F>
T fbox(const T* const boxmin, const T* const boxmax, const T width, const T* const x, F&& f) {
  T r = 0.0;
  for(MInt i = 0; i < nDim; i++) {
    if(boxmin[i] - x[i] > 0.0) {
      r += (boxmin[i] - x[i]) * (boxmin[i] - x[i]);
    } else if(x[i] - boxmax[i] > 0.0) {
      r += (x[i] - boxmax[i]) * (x[i] - boxmax[i]);
    }
  }
  return 1.0 - f(0.0, width, std::sqrt(r));
}

/// Auxiliary function to create a box slope filter (variable slope width)
template <MInt nDim, typename T, typename F>
T fmultibox(const T* const boxmin, const T* const boxmax, const T* const width, const T* const x, F&& f) {
  T filter = 1.0;
  for(MInt i = 0; i < nDim; i++) {
    filter *= f2way(boxmin[i] - width[2 * i], boxmin[i], boxmax[i], boxmax[i] + width[2 * i + 1], x[i], f);
  }
  return filter;
}

/// Auxiliary function to create a sphere slope filter
template <MInt nDim, typename T, typename F>
T fsphere(const T* const center, const T radius, const T width, const T* const x, F&& f) {
  T r = 0.0;
  for(MInt i = 0; i < nDim; i++) {
    r += (x[i] - center[i]) * (x[i] - center[i]);
  }
  return 1.0 - f(0.0, width, std::sqrt(r) - radius);
}

} // namespace detail_

/// \brief Linear slope filter.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-08-11
///
/// \tparam T Datatype that supports substraction and divison.
/// \param[in] a Beginning of slope.
/// \param[in] b End of slope.
/// \param[in] x Point to interpolate.
///
/// \return 0.0 if x < a, 1.0 if x > b, value between 0.0 and 1.0 otherwise.
template <typename T>
T linear(const T a, const T b, const T x) {
  if(x < a) {
    return 0.0;
  } else if(x > b) {
    return 1.0;
  } else {
    return (x - a) / (b - a);
  }
}

template <typename T>
T linear2way(const T a, const T b, const T c, const T d, const T x) {
  return detail_::f2way(a, b, c, d, x, &linear<T>);
}

template <MInt nDim, typename T>
T linearbox(const T* const boxmin, const T* const boxmax, const T width, const T* const x) {
  return detail_::fbox<nDim>(boxmin, boxmax, width, x, &linear<T>);
}

template <MInt nDim, typename T>
T linearmultibox(const T* const boxmin, const T* const boxmax, const T* const width, const T* const x) {
  return detail_::fmultibox<nDim>(boxmin, boxmax, width, x, &linear<T>);
}

template <MInt nDim, typename T>
T linearsphere(const T* const center, const T radius, const T width, const T* const x) {
  return detail_::fsphere<nDim>(center, radius, width, x, &linear<T>);
}

/// \brief Cosine slope filter.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-08-11
///
/// \tparam T Datatype that supports substraction and divison.
/// \param[in] a Beginning of slope.
/// \param[in] b End of slope.
/// \param[in] x Point to interpolate.
///
/// \return 0.0 if x < a, 1.0 if x > b, value between 0.0 and 1.0 otherwise.
template <typename T>
T cos(const T a, const T b, const T x) {
  if(x < a) {
    return 0.0;
  } else if(x > b) {
    return 1.0;
  } else {
    return 0.5 * std::cos(PI / (b - a) * (x - b)) + 0.5;
  }
}

template <typename T>
T cos2way(const T a, const T b, const T c, const T d, const T x) {
  return detail_::f2way(a, b, c, d, x, &cos<T>);
}

template <MInt nDim, typename T>
T cosbox(const T* const boxmin, const T* const boxmax, const T width, const T* const x) {
  return detail_::fbox<nDim>(boxmin, boxmax, width, x, &cos<T>);
}

template <MInt nDim, typename T>
T cosmultibox(const T* const boxmin, const T* const boxmax, const T* const width, const T* const x) {
  return detail_::fmultibox<nDim>(boxmin, boxmax, width, x, &cos<T>);
}

template <MInt nDim, typename T>
T cossphere(const T* const center, const T radius, const T width, const T* const x) {
  return detail_::fsphere<nDim>(center, radius, width, x, &cos<T>);
}

template <MInt nDim, typename T>
T coscylinderzaxis(const T* const center, const T radius, const T width, const T* const x) {
  return detail_::fsphere<nDim>(center, radius, width, x, &cos<T>);
}

} // namespace slope

} // namespace filter
} // namespace maia


/// \brief Filter object for source terms
///
/// \author Marcus Wiens (marcus) <m.wiens@aia.rwth-aachen.de>
/// \date 2016-07-16
template <MInt nDim>
class Filter {
 public:
  Filter(const MInt solverId) : m_solverId(solverId) {}
  void init();

  MFloat filter(const MFloat* const point) const;
  void filter(const MFloat* const points, const MInt count, MFloat* const values);
  MString name() const;

 private:
  // Solver id for property calls
  const MInt m_solverId = -1;

  // Enumaration of possible filter and slope combinations
  enum class FilterType { cossphere, linearsphere, coscylinderzaxis, cosbox, linearbox, cosmultibox, linearmultibox };

  // Filter identificication
  FilterType m_filterId;
  // Minimum coordinates for box-shaped filter
  std::array<MFloat, nDim> m_filterRegionMin{};
  // Maximum coordinates for box-shaped filter
  std::array<MFloat, nDim> m_filterRegionMax{};
  // Center point for sphere-shaped filter
  std::array<MFloat, nDim> m_filterCenter{};
  // Radius for sphere-shaped filter
  MFloat m_filterRadius = -1.0;
  // Width of slope at the edge of the source term
  MFloat m_filterSlopeWidth = -1.0;
  // Width of slope at the edge of the source term (multibox version)
  std::array<MFloat, 2 * nDim> m_filterSlopeWidthMultiBox{};
  // Check for Initialization by init() method
  MBool m_isInitialized = false;
};


/// \fn void Filter<nDim>::init()
/// \brief Initializes the filter. Read all properties and set the filter
///        function
///
/// \author Marcus Wiens (marcus) <m.wiens@aia.rwth-aachen.de>
/// \date 2016-07-16
template <MInt nDim>
void Filter<nDim>::init() {
  // sets the shape of the filter
  MString filterShape;
  // sets the calculating function of the slope
  MString filterSlopeType;

  // Map to compare input filtertype with implemented filters
  std::map<MString, FilterType> filterMap;

  // Define filterlist
  // Note: if you add a new filter, do not forget to also adapt the switch
  // statements in filter() and name()
  filterMap["cossphere"] = FilterType::cossphere;
  filterMap["linearsphere"] = FilterType::linearsphere;
  filterMap["coscylinderzaxis"] = FilterType::coscylinderzaxis;
  filterMap["cosbox"] = FilterType::cosbox;
  filterMap["linearbox"] = FilterType::linearbox;
  filterMap["cosmultibox"] = FilterType::cosmultibox;
  filterMap["linearmultibox"] = FilterType::linearmultibox;

  /*! \page propertyPage1
    \section filterShape
    <code>MString Filter::filterShape</code>\n
    default = <code>none</code>\n \n
    Specify shape of filter.\n
    Possible values are:
    <ul>
      <li>"sphere"</li>
      <li>"cylinderzaxis"</li>
      <li>"box"</li>
      <li>"multibox"</li>
    </ul>
    Keywords: <i>FILTER</i>
  */
  filterShape = Context::getSolverProperty<MString>("filterShape", m_solverId, AT_);

  /*! \page propertyPage1
    \section filterSlopeType
    <code>MString Filter::filterSlopeType</code>\n
    default = <code>none</code>\n \n
    Specify slope type of filter.\n
    Possible values are:
    <ul>
      <li>"linear"</li>
      <li>"cos"</li>
    </ul>
    Keywords: <i>FILTER</i>
  */
  filterSlopeType = Context::getSolverProperty<MString>("filterSlopeType", m_solverId, AT_);

  // Check if filter type exists and sets the filter id. If the filter type
  // does not exists it throws an error
  auto search = filterMap.find(filterSlopeType + filterShape);
  if(search != filterMap.end()) {
    m_filterId = filterMap[filterSlopeType + filterShape];
  } else {
    TERMM(1, "Unknown filter type! Please check your property file.");
  }

  // Width of filter slope
  if(filterShape == "multibox") {
    for(MInt i = 0; i < 2 * nDim; i++) {
      m_filterSlopeWidthMultiBox[i] = Context::getSolverProperty<MFloat>("filterSlopeWidth", m_solverId, AT_, i);
    }
  } else {
    m_filterSlopeWidth = Context::getSolverProperty<MFloat>("filterSlopeWidth", m_solverId, AT_);
  }

  // read filter specific properties
  if(filterShape == "sphere" || filterShape == "cylinderzaxis") {
    // read properties of sphere filter
    for(MInt i = 0; i < nDim; i++) {
      // /*! \page propertyPage1
      //   \section filterCenter
      //   <code>MFloat DgCcApeSourceFiles::m_filterCenter </code>\n
      //   default = <code>none</code>\n \n
      //   The centor of the sphere that specifies the source term
      //   region.\n
      //   Keywords: <i>COUPLING, SOURCE_TERM, FILTER</i>
      // */
      m_filterCenter[i] = Context::getSolverProperty<MFloat>("filterCenter", m_solverId, AT_, i);
    }
    // /*! \page propertyPage1
    //   \section filterRadius
    //   <code>MFloat DgCcApeSourceFiles::m_filterRadius </code>\n
    //   default = <code>none</code>\n \n
    //   The radius of the sphere that specifies the source term
    //   region.\n
    //   Keywords: <i>COUPLING, SOURCE_TERM, FILTER</i>
    // */
    m_filterRadius = Context::getSolverProperty<MFloat>("filterRadius", m_solverId, AT_);
  } else {
    // read properties of box filter
    for(MInt i = 0; i < nDim; i++) {
      // /*! \page propertyPage1
      //   \section filterRegionMin
      //   <code>MFloat DgCcApeSourceFiles::m_filterRegionMin </code>\n
      //   default = <code>none</code>\n \n
      //   The coordinates of the lower left corner of the source term
      //   region.\n
      //   Keywords: <i>COUPLING, SOURCE_TERM, FILTER</i>
      // */
      m_filterRegionMin[i] = Context::getSolverProperty<MFloat>("filterRegionMin", m_solverId, AT_, i);

      // /*! \page propertyPage1
      //   \section filterRegionMax
      //   <code>MFloat DgCcApeSourceFiles::m_filterRegionMax</code>\n
      //   default = <code>none</code>\n \n
      //   The coordinates of the upper right corner of the source term
      //   region.\n
      //   Keywords: <i>COUPLING, SOURCE_TERM, FILTER</i>
      // */
      m_filterRegionMax[i] = Context::getSolverProperty<MFloat>("filterRegionMax", m_solverId, AT_, i);
    }
  }

  // Initialization is done!
  m_isInitialized = true;
}


/// \brief Calculates filter vaule for a given coordinates and returns it
///
/// \author Marcus Wiens (marcus) <m.wiens@aia.rwth-aachen.de>
/// \date 2016-07-16
template <MInt nDim>
MFloat Filter<nDim>::filter(const MFloat* const coordinates) const {
  // Check for correkt initialization of the filter object
  if(!m_isInitialized) {
    TERMM(1, "Filter Object was not initialized correctly!");
  }

  // selection of the type of filter
  switch(m_filterId) {
    // cosphere
    case FilterType::cossphere:
      // calculate filter value and return it
      return maia::filter::slope::cossphere<nDim>(&m_filterCenter[0], m_filterRadius, m_filterSlopeWidth, coordinates);

    // linearsphere
    case FilterType::linearsphere:
      TERMM(1, "linearsphere type of filter has not been tested yet!");

      // calculate filter value and return it
      return maia::filter::slope::linearsphere<nDim>(&m_filterCenter[0], m_filterRadius, m_filterSlopeWidth,
                                                     coordinates);
    // coscylinderzaxis
    case FilterType::coscylinderzaxis:
      // calculate filter value and return it
      // Simply uses sphere filter
      return maia::filter::slope::cossphere<2>(&m_filterCenter[0], m_filterRadius, m_filterSlopeWidth, coordinates);

    // cosbox
    case FilterType::cosbox:
      TERMM(1, "Cosbox type of filter has not been tested yet!");

      // calculate filter value and return it
      return maia::filter::slope::cosbox<nDim>(&m_filterRegionMin[0], &m_filterRegionMax[0], m_filterSlopeWidth,
                                               coordinates);

    // linearbox
    case FilterType::linearbox:
      TERMM(1, "linearbox type of filter has not been tested yet!");

      // calculate filter value and return it
      return maia::filter::slope::linearbox<nDim>(&m_filterRegionMin[0], &m_filterRegionMax[0], m_filterSlopeWidth,
                                                  coordinates);

    // cosmultibox
    case FilterType::cosmultibox:
      // calculate filter value and return it
      return maia::filter::slope::cosmultibox<nDim>(&m_filterRegionMin[0], &m_filterRegionMax[0],
                                                    &m_filterSlopeWidthMultiBox[0], coordinates);

    // linearmultibox
    case FilterType::linearmultibox:

      // calculate filter value and return it
      return maia::filter::slope::linearmultibox<nDim>(&m_filterRegionMin[0], &m_filterRegionMax[0],
                                                       &m_filterSlopeWidthMultiBox[0], coordinates);

    default:
      TERMM(1, "Unknown filter type! Please check your property file.");
      break;
  }
}


/// \brief Calculates filter values of a coordinates array based on the type of
/// filter
///
/// \author Marcus Wiens (marcus) <m.wiens@aia.rwth-aachen.de>
/// \date 2016-07-16
template <MInt nDim>
void Filter<nDim>::filter(const MFloat* const coordinates, const MInt count, MFloat* const values) {
  // Check for correct initialization of the filter object
  if(!m_isInitialized) {
    TERMM(1, "Filter Object was not initialized correctly!");
  }

  // selection of the type of filter
  switch(m_filterId) {
    // cosphere
    case FilterType::cossphere:
      // calculate all filter values for each point
      for(MInt j = 0; j < count; j++) {
        values[j] = maia::filter::slope::cossphere<nDim>(&m_filterCenter[0], m_filterRadius, m_filterSlopeWidth,
                                                         coordinates + j * nDim);
      }
      break;

    // linearsphere
    case FilterType::linearsphere:
      TERMM(1, "linearsphere type of filter has not been tested yet!");

      // calculate all filter values for each point
      for(MInt j = 0; j < count; j++) {
        values[j] = maia::filter::slope::linearsphere<nDim>(&m_filterCenter[0], m_filterRadius, m_filterSlopeWidth,
                                                            coordinates + j * nDim);
      }
      break;

    case FilterType::coscylinderzaxis:
      // calculate all filter values for each point
      for(MInt j = 0; j < count; j++) {
        values[j] = maia::filter::slope::cossphere<2>(&m_filterCenter[0], m_filterRadius, m_filterSlopeWidth,
                                                      coordinates + j * nDim);
      }
      break;

    // cosbox
    case FilterType::cosbox:
      TERMM(1, "Cosbox type of filter has not been tested yet!");

      // calculate all filter values for each point
      for(MInt j = 0; j < count; j++) {
        values[j] = maia::filter::slope::cosbox<nDim>(&m_filterRegionMin[0], &m_filterRegionMax[0], m_filterSlopeWidth,
                                                      coordinates + j * nDim);
      }
      break;
    // linearbox
    case FilterType::linearbox:
      TERMM(1, "linearbox type of filter has not been tested yet!");

      // calculate all filter values for each point
      for(MInt j = 0; j < count; j++) {
        values[j] = maia::filter::slope::linearbox<nDim>(&m_filterRegionMin[0], &m_filterRegionMax[0],
                                                         m_filterSlopeWidth, coordinates + j * nDim);
      }
      break;

    // cosmultibox
    case FilterType::cosmultibox:
      // calculate all filter values for each point
      for(MInt j = 0; j < count; j++) {
        values[j] = maia::filter::slope::cosmultibox<nDim>(&m_filterRegionMin[0], &m_filterRegionMax[0],
                                                           &m_filterSlopeWidthMultiBox[0], coordinates + j * nDim);
      }
      break;
    // linearmultibox
    case FilterType::linearmultibox:

      // calculate all filter values for each point
      for(MInt j = 0; j < count; j++) {
        values[j] = maia::filter::slope::linearmultibox<nDim>(&m_filterRegionMin[0], &m_filterRegionMax[0],
                                                              &m_filterSlopeWidthMultiBox[0], coordinates + j * nDim);
      }
      break;

    default:
      TERMM(1, "Unknown filter type! Please check your property file.");
      break;
  }
}


/// \brief Return filter name as a string.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-12-06
///
/// \return The internally used filter name.
template <MInt nDim>
MString Filter<nDim>::name() const {
  TRACE();

  MString filterName;
  switch(m_filterId) {
    case FilterType::cossphere:
      filterName = "cossphere";
      break;
    case FilterType::linearsphere:
      filterName = "linearsphere";
      break;
    case FilterType::coscylinderzaxis:
      filterName = "coscylinderzaxis";
      break;
    case FilterType::cosbox:
      filterName = "cosbox";
      break;
    case FilterType::linearbox:
      filterName = "linearbox";
      break;
    case FilterType::cosmultibox:
      filterName = "cosmultibox";
      break;
    case FilterType::linearmultibox:
      filterName = "linearmultibox";
      break;
    default:
      TERMM(1, "Unknown filter type! Please check your property file.");
      break;
  }

  return filterName;
}

#endif // FILTER_H_
