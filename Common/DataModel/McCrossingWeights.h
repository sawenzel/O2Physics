// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file McCrossingWeights.h
/// \brief Per MC collision: number of MC collisions in its bunch crossing, TF-boundary copy flag and weight (O2-7330)
/// \author Sandro Wenzel, CERN

#ifndef COMMON_DATAMODEL_MCCROSSINGWEIGHTS_H_
#define COMMON_DATAMODEL_MCCROSSINGWEIGHTS_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace mccrossing
{
DECLARE_SOA_COLUMN(CrossingMult, crossingMult, uint8_t);   //! MC collisions in the same bunch crossing (TF-boundary copies not counted)
DECLARE_SOA_COLUMN(IsTFCopy, isTFCopy, bool);              //! re-simulated copy of a collision of the previous timeframe; skip it
DECLARE_SOA_COLUMN(CrossingWeight, crossingWeight, float); //! weight restoring Poisson crossing multiplicities; 1 if no known issue
} // namespace mccrossing

DECLARE_SOA_TABLE(McCrossingWeights, "AOD", "MCCROSSWEIGHT", //! joinable with McCollisions
                  mccrossing::CrossingMult, mccrossing::IsTFCopy, mccrossing::CrossingWeight);
using McCrossingWeight = McCrossingWeights::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_MCCROSSINGWEIGHTS_H_
