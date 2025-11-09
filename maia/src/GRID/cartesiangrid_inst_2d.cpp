// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiangrid.cpp"

/// Explicit instantiation of CartesianGrid for 2D
template class CartesianGrid<2>;

// The following line is only required to prevent missing symbols when linking
// in intel/production mode
template void CartesianGrid<2>::calculateNoOffspringsAndWorkload(Collector<void>*, MInt);

template MInt
CartesianGrid<2>::findContainingPartitionCell<false>(const MFloat* const coord, const MInt solverId,
                                                     function<MFloat*(MInt, MFloat* const)> correctCellCoord);
template MInt
CartesianGrid<2>::findContainingPartitionCell<true>(const MFloat* const coord, const MInt solverId,
                                                    function<MFloat*(MInt, MFloat* const)> correctCellCoord);

template MInt CartesianGrid<2>::findContainingLeafCell<false>(
    const MFloat* coord, std::function<MFloat*(MInt, MFloat* const)> correctCellCoord, const MInt solverId);

template MInt CartesianGrid<2>::findContainingLeafCell<true>(
    const MFloat* coord, std::function<MFloat*(MInt, MFloat* const)> correctCellCoord, const MInt solverId);

template MInt
CartesianGrid<2>::findContainingLeafCell<false>(const MFloat* coord, const MInt startId,
                                                std::function<MFloat*(MInt, MFloat* const)> correctCellCoord,
                                                const MInt solverId, const MBool allowNonLeafHalo);

template MInt
CartesianGrid<2>::findContainingLeafCell<true>(const MFloat* coord, const MInt startId,
                                               std::function<MFloat*(MInt, MFloat* const)> correctCellCoord,
                                               const MInt solverId, const MBool allowNonLeafHalo);

template MInt CartesianGrid<2>::findContainingHaloCell<false>(const MFloat* const coord, const MInt solverId,
                                                              MInt domainId, MBool onlyPartitionCells,
                                                              function<MFloat*(MInt, MFloat* const)> correctCellCoord);

template MInt CartesianGrid<2>::findContainingHaloCell<true>(const MFloat* const coord, const MInt solverId,
                                                             MInt domainId, MBool onlyPartitionCells,
                                                             function<MFloat*(MInt, MFloat* const)> correctCellCoord);
