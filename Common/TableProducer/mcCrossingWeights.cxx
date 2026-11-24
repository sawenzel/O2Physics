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

/// \file mcCrossingWeights.cxx
/// \brief Per MC collision: bunch-crossing multiplicity, TF-boundary copy flag and a weight correcting O2-7330
/// \author Sandro Wenzel, CERN
///
/// Productions simulated with MC-prod-2026-v9 to v14 put too many collisions into one bunch crossing
/// (https://its.cern.ch/jira/browse/O2-7330): the sampler revisited an occupied crossing a geometric number of times.
/// Which crossings are occupied stays correct, so the mean number of collisions per bunch, mu, is measured here from
/// the occupancy of the crossings in each dataframe. Each group of k MC collisions sharing a crossing then gets the weight
///   w(k) = P_Poisson(k | mu_b) / P_bug(k | mu_b)      (both for an occupied crossing)
/// when the production's MCProdInfo lists O2-7330 in its KnownIssues, and w = 1 otherwise.
///
/// Collisions of the last orbit of a timeframe are re-simulated at the start of the next one (orbitsEarly) and stored
/// again; those copies get isTFCopy = true and must be skipped when counting generated collisions.
///
/// Use: join aod::McCollisions with aod::McCrossingWeights; skip rows with isTFCopy(); weight generated and
/// reconstructed collisions (via their MC label) with crossingWeight().
///
/// QA histograms per dataframe: mu, rho = observed / Poisson-expected MC collisions (about exp(mu) where the bug
/// acted, 1 elsewhere) and its significance z, crossing multiplicities, weights, and the number of TF copies.

#include "Common/DataModel/McCrossingWeights.h"

#include "Common/Core/MetadataHelper.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsFT0/EventsPerBc.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

o2::common::core::MetadataHelper metadataInfo;

namespace
{
constexpr int MaxK = 20;                      // crossings with more collisions get the weight of MaxK
constexpr int MaxVisits = 60;                 // visits of one crossing summed in the O2-7330 model
constexpr double VisitProbabilityCut = 1e-15; // stop summing visits below this probability
constexpr int NewtonIterations = 30;          // iterations solving the occupancy equation for mu
constexpr int MaxStoredMult = 255;            // crossing multiplicity is stored as uint8_t

/// P(k | occupied) for one bunch with mean mu: correct sampler (Poisson) or the O2-7330 sampler
std::vector<double> multiplicityDistribution(double mu, bool buggy)
{
  std::vector<double> q(MaxK + 1, 0.);
  double norm = 1. - std::exp(-mu);
  double term = std::exp(-mu);
  for (int k = 1; k <= MaxK; ++k) {
    term *= mu / k;
    q[k] = term / norm;
  }
  if (!buggy) {
    return q;
  }
  // the crossing is visited V times, V geometric with P(another visit) = 1 - exp(-mu); each visit adds a draw from q
  double p = norm;
  std::vector<double> out(MaxK + 1, 0.), conv(MaxK + 1, 0.);
  conv[0] = 1.;
  double pv = 1. - p;
  for (int v = 1; v < MaxVisits && pv > VisitProbabilityCut; ++v) {
    std::vector<double> next(MaxK + 1, 0.);
    for (int k = 1; k <= MaxK; ++k) {
      for (int j = 0; j < k; ++j) {
        next[k] += conv[j] * q[k - j];
      }
    }
    conv = next;
    for (int k = 1; k <= MaxK; ++k) {
      out[k] += pv * conv[k];
    }
    pv *= p;
  }
  return out;
}
} // namespace

struct McCrossingWeightsTask {
  Produces<aod::McCrossingWeights> crossingWeights;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB url"};
  Configurable<int> mode{"mode", -1, "-1: from MCProdInfo KnownIssues; 0: never weight; 1: always apply the O2-7330 weights"};
  Configurable<std::string> issueId{"issueId", "O2-7330", "known-issue id that switches the weights on"};
  Configurable<std::string> mcProdInfoUser{"mcProdInfoUser", "aliprod", "owner of the MCProdInfo objects"};
  Configurable<bool> useBunchShape{"useBunchShape", true, "shape mu per bunch with FT0/Calib/EventsPerBc"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // per run
  int mRun = -1;
  bool mApply = false;
  int64_t mOrbitsPerTF = 32;
  std::vector<int> mCollidingBCs; // bunch slots of the colliding bunches
  std::vector<double> mShape;     // EventsPerBc shape, mean 1 over the colliding bunches
  std::unordered_map<int, int> mSlotIndex;
  std::string mProductionTag;

  void init(InitContext&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    mProductionTag = metadataInfo.isFullyDefined() ? metadataInfo.get("LPMProductionTag") : "";
    histos.add("mu", "mean collisions per colliding bunch per DF;#mu;DFs", HistType::kTH1D, {{200, 0., 0.4}});
    histos.add("rho", "observed / Poisson-expected MC collisions per DF;#rho;DFs", HistType::kTH1D, {{200, 0.8, 1.4}});
    histos.add("z", "excess of collisions over Poisson per DF (significance);z;DFs", HistType::kTH1D, {{200, -10., 40.}});
    histos.add("k", "MC collisions per occupied crossing;k;crossings", HistType::kTH1D, {{10, 0.5, 10.5}});
    histos.add("weight", "weight per MC collision;w;MC collisions", HistType::kTH1D, {{240, 0., 1.2}});
    histos.add("applied", "weights applied (1) or not (0), per DF;applied;DFs", HistType::kTH1D, {{2, -0.5, 1.5}});
    histos.add("copies", "TF-boundary copies per DF;copies;DFs", HistType::kTH1D, {{200, -0.5, 199.5}});
  }

  void setupRun(int run)
  {
    mRun = run;
    // CCDB access as in Common/Tools/EventSelectionModule.h
    auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run, mProductionTag);
    mOrbitsPerTF = runInfo.orbitsPerTF > 0 ? runInfo.orbitsPerTF : 32;
    mCollidingBCs.clear();
    mSlotIndex.clear();
    if (runInfo.grpLHC) {
      const auto& filling = runInfo.grpLHC->getBunchFilling();
      for (int bc = 0; bc < o2::constants::lhc::LHCMaxBunches; ++bc) {
        if (filling.testBC(bc)) {
          mSlotIndex[bc] = mCollidingBCs.size();
          mCollidingBCs.push_back(bc);
        }
      }
    }
    mShape.assign(mCollidingBCs.size(), 1.);
    if (useBunchShape) {
      auto epb = ccdb->getSpecific<o2::ft0::EventsPerBc>("FT0/Calib/EventsPerBc", (runInfo.sor + runInfo.eor) / 2);
      if (epb) {
        double sum = 0.;
        for (size_t i = 0; i < mCollidingBCs.size(); ++i) {
          mShape[i] = std::abs(epb->histogram[mCollidingBCs[i]]);
          sum += mShape[i];
        }
        for (auto& s : mShape) {
          s = sum > 0 ? s * mShape.size() / sum : 1.;
        }
      }
    }
    // the decision: an explicit mode, or the production's KnownIssues
    if (mode >= 0) {
      mApply = mode == 1;
    } else {
      auto info = o2::parameters::AggregatedRunInfo::getMCProdInfo(o2::ccdb::BasicCCDBManager::instance(), run, mProductionTag, mcProdInfoUser.value);
      auto it = info.find("KnownIssues");
      mApply = false;
      if (it != info.end()) {
        std::stringstream ss(it->second);
        std::string id;
        while (std::getline(ss, id, ',')) {
          mApply = mApply || id == issueId.value;
        }
      }
    }
    LOGP(info, "McCrossingWeights: run {} production {} orbitsPerTF {} colliding bunches {} weights {}", run,
         mProductionTag, mOrbitsPerTF, mCollidingBCs.size(), mApply ? "applied" : "not applied");
  }

  void process(aod::McCollisions const& mcCollisions, aod::BCs const& bcs, aod::Origins const& origins,
               soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions)
  {
    if (mcCollisions.size() == 0 || bcs.size() == 0) {
      return;
    }
    int run = bcs.begin().runNumber();
    if (run != mRun) {
      setupRun(run);
    }
    // timeframes of this dataframe: first orbit of each from Origins (DF id = run << 32 | first orbit)
    uint64_t firstOrbit = UINT64_MAX;
    int nTF = 0;
    for (const auto& origin : origins) {
      firstOrbit = std::min<uint64_t>(firstOrbit, origin.dataframeID() & 0xffffffffULL);
      ++nTF;
    }
    // MC collisions that have a reconstructed collision
    std::set<int64_t> hasReco;
    for (const auto& col : collisions) {
      if (col.has_mcCollision()) {
        hasReco.insert(col.mcCollisionId());
      }
    }
    // TF-boundary copies: before the first TF of the DF, or a second entry with the same crossing and vertex
    const size_t n = mcCollisions.size();
    std::vector<uint64_t> gbc(n);
    std::vector<bool> isCopy(n, false);
    std::map<std::tuple<uint64_t, float, float, float, float>, int64_t> seen;
    std::vector<int64_t> order(n);
    for (const auto& mc : mcCollisions) {
      gbc[mc.globalIndex() - mcCollisions.begin().globalIndex()] = mc.bc_as<aod::BCs>().globalBC();
    }
    for (size_t i = 0; i < n; ++i) {
      order[i] = i;
    }
    // keep the entry that has a reconstructed collision, else the first
    std::stable_sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
      return hasReco.count(a + mcCollisions.begin().globalIndex()) > hasReco.count(b + mcCollisions.begin().globalIndex());
    });
    for (const auto& i : order) {
      auto mc = mcCollisions.rawIteratorAt(i);
      if (firstOrbit != UINT64_MAX && gbc[i] / o2::constants::lhc::LHCMaxBunches < firstOrbit) {
        isCopy[i] = true;
        continue;
      }
      auto key = std::make_tuple(gbc[i], mc.posX(), mc.posY(), mc.posZ(), mc.t());
      if (!seen.emplace(key, i).second) {
        isCopy[i] = true;
      }
    }
    // crossing multiplicities
    std::unordered_map<uint64_t, int> kOf;
    int nColl = 0;
    for (size_t i = 0; i < n; ++i) {
      if (!isCopy[i]) {
        ++kOf[gbc[i]];
        ++nColl;
      }
    }
    // mu from occupancy: sum_b (1 - exp(-mu s_b)) * orbits = occupied crossings
    const double orbits = double(std::max(nTF, 1)) * mOrbitsPerTF;
    const double nOcc = kOf.size();
    double mu = 0.;
    if (!mShape.empty() && nOcc > 0) {
      mu = -std::log(std::max(1e-9, 1. - nOcc / (orbits * mShape.size())));
      for (int it = 0; it < NewtonIterations; ++it) {
        double f = -nOcc, d = 0.;
        for (const auto& s : mShape) {
          f += (1. - std::exp(-mu * s)) * orbits;
          d += s * std::exp(-mu * s) * orbits;
        }
        mu -= f / d;
      }
    }
    double expected = 0.;
    for (const auto& s : mShape) {
      expected += mu * s * orbits;
    }
    if (expected > 0) {
      histos.fill(HIST("mu"), mu);
      histos.fill(HIST("rho"), nColl / expected);
      histos.fill(HIST("z"), (nColl - expected) / std::sqrt(std::max(expected - nOcc, 1.)));
    }
    histos.fill(HIST("applied"), mApply ? 1 : 0);
    histos.fill(HIST("copies"), std::count(isCopy.begin(), isCopy.end(), true));
    for (const auto& [bc, k] : kOf) {
      histos.fill(HIST("k"), k);
    }
    // weights per bunch slot, computed once per DF
    std::unordered_map<int, std::vector<double>> wOfSlot;
    auto weight = [&](uint64_t bc, int k) -> double {
      if (!mApply || mu <= 0.) {
        return 1.;
      }
      int slot = bc % o2::constants::lhc::LHCMaxBunches;
      auto it = wOfSlot.find(slot);
      if (it == wOfSlot.end()) {
        auto idx = mSlotIndex.find(slot);
        double mub = mu * (idx != mSlotIndex.end() ? mShape[idx->second] : 1.);
        auto ok = multiplicityDistribution(mub, false);
        auto bug = multiplicityDistribution(mub, true);
        std::vector<double> w(MaxK + 1, 1.);
        for (int j = 1; j <= MaxK; ++j) {
          w[j] = bug[j] > 0 ? ok[j] / bug[j] : 0.;
        }
        it = wOfSlot.emplace(slot, std::move(w)).first;
      }
      return it->second[std::min(k, MaxK)];
    };
    for (size_t i = 0; i < n; ++i) {
      int k = isCopy[i] ? 0 : kOf[gbc[i]];
      float w = isCopy[i] ? 0.f : weight(gbc[i], k);
      if (!isCopy[i]) {
        histos.fill(HIST("weight"), w);
      }
      crossingWeights(static_cast<uint8_t>(std::min(k, MaxStoredMult)), isCopy[i], w);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<McCrossingWeightsTask>(cfgc)};
}
