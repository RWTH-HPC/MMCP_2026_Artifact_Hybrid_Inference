// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiangrid.cpp"

/// Explicit instantiation of CartesianGrid for 3D
template class CartesianGrid<3>;

// The following line is only required to prevent missing symbols when linking
// in intel/production mode
template void CartesianGrid<3>::calculateNoOffspringsAndWorkload(Collector<void>*, MInt);

template MInt
CartesianGrid<3>::findContainingPartitionCell<false>(const MFloat* const coord, const MInt solverId,
                                                     function<MFloat*(MInt, MFloat* const)> correctCellCoord);
template MInt
CartesianGrid<3>::findContainingPartitionCell<true>(const MFloat* const coord, const MInt solverId,
                                                    function<MFloat*(MInt, MFloat* const)> correctCellCoord);

template MInt CartesianGrid<3>::findContainingLeafCell<false>(
    const MFloat* coord, std::function<MFloat*(MInt, MFloat* const)> correctCellCoord, const MInt solverId);

template MInt CartesianGrid<3>::findContainingLeafCell<true>(
    const MFloat* coord, std::function<MFloat*(MInt, MFloat* const)> correctCellCoord, const MInt solverId);

template MInt
CartesianGrid<3>::findContainingLeafCell<false>(const MFloat* coord, const MInt startId,
                                                std::function<MFloat*(MInt, MFloat* const)> correctCellCoord,
                                                const MInt solverId, const MBool allowNonLeafHalo = false);

template MInt
CartesianGrid<3>::findContainingLeafCell<true>(const MFloat* coord, const MInt startId,
                                               std::function<MFloat*(MInt, MFloat* const)> correctCellCoord,
                                               const MInt solverId, const MBool allowNonLeafHalo = false);

template MInt CartesianGrid<3>::findContainingHaloCell<false>(const MFloat* const coord, const MInt solverId,
                                                              MInt domainId, MBool onlyPartitionCells,
                                                              function<MFloat*(MInt, MFloat* const)> correctCellCoord);

template MInt CartesianGrid<3>::findContainingHaloCell<true>(const MFloat* const coord, const MInt solverId,
                                                             MInt domainId, MBool onlyPartitionCells,
                                                             function<MFloat*(MInt, MFloat* const)> correctCellCoord);
