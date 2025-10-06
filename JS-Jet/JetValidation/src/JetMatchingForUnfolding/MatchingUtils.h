#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <jetbase/JetContainer.h>
#include <jetbase/Jet.h>

// exact same phi wrapping & ΔR style as in JetUnfolding
inline float calculateDeltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  if (dphi > M_PI)  dphi -= 2.f * M_PI;
  if (dphi < -M_PI) dphi += 2.f * M_PI;
  return std::sqrt(deta * deta + dphi * dphi);
}

struct MatchResult {
  std::vector<const Jet*> matchedTruthJets;   // aligned 1:1 with matchedRecoJets
  std::vector<const Jet*> matchedRecoJets;    // aligned 1:1 with matchedTruthJets
  std::vector<int>        recoToTruthIndex;   // size = reco->size(), -1 if unmatched
  std::vector<bool>       recoMatched;        // size = reco->size()
  std::vector<bool>       truthMatched;       // size = truth->size()
};

// Greedy truth->reco match in descending pT, same cuts as JetUnfolding
inline MatchResult MatchJetsLikeUnfolding(
    JetContainer* jetsReco,
    JetContainer* jetsTruth,
    float dRMax      = 0.1f,
    float recoPtMin  = 5.0f,
    float etaMaxAbs  = 0.7f,
    float truthPtMin = 52.0f)
{
  MatchResult out;
  if (!jetsReco || !jetsTruth) return out;

  const size_t nR = jetsReco->size();
  const size_t nT = jetsTruth->size();

  out.recoToTruthIndex.assign(nR, -1);
  out.recoMatched.assign(nR, false);
  out.truthMatched.assign(nT, false);

  // sort indices by descending pT (truth, then reco)
  std::vector<size_t> tIdx(nT), rIdx(nR);
  std::iota(tIdx.begin(), tIdx.end(), 0);
  std::iota(rIdx.begin(), rIdx.end(), 0);

  std::sort(tIdx.begin(), tIdx.end(), [&](size_t i, size_t j){
    return (*jetsTruth)[i]->get_pt() > (*jetsTruth)[j]->get_pt();
  });
  std::sort(rIdx.begin(), rIdx.end(), [&](size_t i, size_t j){
    return (*jetsReco)[i]->get_pt() > (*jetsReco)[j]->get_pt();
  });

  // loop truth jets (greedy): pick the closest *unmatched* reco within ΔR and kinematic cuts
  for (size_t k = 0; k < tIdx.size(); ++k) {
    const auto* tj = (*jetsTruth)[tIdx[k]];
    const float tpt  = tj->get_pt();
    const float teta = tj->get_eta();
    const float tphi = tj->get_phi();

    if (tpt < truthPtMin || std::abs(teta) >= etaMaxAbs) continue;

    float bestDR = dRMax;
    ssize_t bestReco = -1;

    for (size_t m = 0; m < rIdx.size(); ++m) {
      const size_t ri = rIdx[m];
      if (out.recoMatched[ri]) continue;

      const auto* rj = (*jetsReco)[ri];
      const float rpt  = rj->get_pt();
      const float reta = rj->get_eta();
      const float rphi = rj->get_phi();

      if (rpt <= recoPtMin || std::abs(reta) > etaMaxAbs) continue;

      const float dR = calculateDeltaR(teta, tphi, reta, rphi);
      if (dR < bestDR) {
        bestDR = dR;
        bestReco = static_cast<ssize_t>(ri);
      }
    }

    if (bestReco >= 0) {
      out.matchedTruthJets.push_back(tj);
      out.matchedRecoJets.push_back((*jetsReco)[bestReco]);
      out.recoMatched[bestReco] = true;
      out.truthMatched[tIdx[k]] = true;
      out.recoToTruthIndex[bestReco] = static_cast<int>(tIdx[k]);
    }
  }

  return out;
}
