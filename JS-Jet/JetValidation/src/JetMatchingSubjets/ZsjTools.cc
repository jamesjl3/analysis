#include "ZsjTools.h"

#include <jetbase/Jet.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerDefs.h>
#include <jetbackground/TowerBackground.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>

using namespace fastjet;

namespace {
  inline double wrapPhi(double a){
    while (a >  M_PI) a -= 2*M_PI;
    while (a < -M_PI) a += 2*M_PI;
    return a;
  }

  // Try to retrieve geometry for (ieta,iphi) with a primary CaloId, and
  // fall back to other ids if needed (helps when retowers use HCALIN grid).
  static const RawTowerGeom* get_geom_flexible(
      RawTowerGeomContainer* gc, int ieta, int iphi,
      RawTowerDefs::CalorimeterId primary)
  {
    if (!gc) return nullptr;

    RawTowerDefs::CalorimeterId try_ids[] = {
      primary, RawTowerDefs::CEMC, RawTowerDefs::HCALIN, RawTowerDefs::HCALOUT
    };

    for (auto id : try_ids) {
      if (id == primary || id != primary) {
        const auto key = RawTowerDefs::encode_towerid(id, ieta, iphi);
        if (auto* g = gc->get_tower_geometry(key)) return g;
      }
    }
    return nullptr;
  }

  // Build massless PseudoJets from *all three* calos.
  // UE modulation is applied when doUnsub=true (v2, Psi2).
  static void BuildPseudoJetsFromRecoJet_AllCalos(
    Jet* jet,
    TowerInfoContainer* em,  RawTowerGeomContainer* geomEM,
    TowerInfoContainer* ih,  RawTowerGeomContainer* geomIH,
    TowerInfoContainer* oh,  RawTowerGeomContainer* geomOH,
    TowerBackground* bg,
    double v2, double psi2, bool doUnsub,
    std::vector<PseudoJet>& out)
  {
    out.clear();
    if (!jet) return;
    /*
    const std::vector<float>* ueEM = nullptr;
    const std::vector<float>* ueIH = nullptr;
    const std::vector<float>* ueOH = nullptr;
    if (bg) {
      try { ueEM = &bg->get_UE(0); } catch(...) {}
      try { ueIH = &bg->get_UE(1); } catch(...) {}
      try { ueOH = &bg->get_UE(2); } catch(...) {}
    }
    */
    auto push_from = [&](TowerInfoContainer* cont,
                         RawTowerGeomContainer* geom,
			 //                   const std::vector<float>* ueArr,
                         RawTowerDefs::CalorimeterId caloId,
                         unsigned comp_code_a, unsigned comp_code_b)
    {
      if (!cont || !geom) return;

      for (const auto& comp : jet->get_comp_vec()) {
        if (comp.first != (int)comp_code_a && comp.first != (int)comp_code_b) continue;

        const unsigned ch = comp.second;
        auto* twr = cont->get_tower_at_channel(ch);
        if (!twr) continue;

        const unsigned calokey = cont->encode_key(ch);
        const int ieta = cont->getTowerEtaBin(calokey);
        const int iphi = cont->getTowerPhiBin(calokey);

        // Be robust to grids (retower cases)
        const RawTowerGeom* tg = get_geom_flexible(geom, ieta, iphi, caloId);
        if (!tg) continue;

        const double eta = tg->get_eta();
        const double phi = tg->get_phi();

        double E = twr->get_energy();
	/*
	if (doUnsub && ueArr && ieta >= 0 && ieta < (int)ueArr->size()) {
          double UE = (*ueArr)[ieta];
          UE = UE * (1.0 + 2.0 * v2 * std::cos(2.0 * (phi - psi2)));
          E += UE;
        }
	*/
        const double pt = E / std::cosh(eta);
        const double px = pt * std::cos(phi);
        const double py = pt * std::sin(phi);
        const double pz = pt * std::sinh(eta);

        out.emplace_back(px, py, pz, E);
      }
    };

    // CEMC (typical comp codes seen in your code: 14, 29)
    push_from(em, geomEM,/* ueEM,*/ RawTowerDefs::CEMC,   14, 29);
    // HCALIN (15, 30)
    push_from(ih, geomIH,/* ueIH,*/ RawTowerDefs::HCALIN, 15, 30);
    // HCALOUT (16, 31)
    push_from(oh, geomOH,/* ueOH,*/ RawTowerDefs::HCALOUT,16, 31);
  }

  // pick leading two subjets, require each above threshold, return (theta,z)
  static bool SubjetZTheta_AKT01(const std::vector<PseudoJet>& constituents,
                                 double ptMinSubjet,
                                 double& theta, double& z)
  {
    if (constituents.empty()) return false;

    JetDefinition jd_sub(antikt_algorithm, 0.1);
    ClusterSequence cs_sub(constituents, jd_sub);
    auto subjets = sorted_by_pt(cs_sub.inclusive_jets());

    if (subjets.size() < 2) return false;

    const PseudoJet& sj1 = subjets[0];
    const PseudoJet& sj2 = subjets[1];

    std::cout << "SJ1 pt, eta ,phi " << sj1.pt() << ", " << sj1.eta() << ", " << sj1.phi_std() << "\n";
    std::cout << "SJ2 pt, eta ,phi " << sj2.pt() << ", " << sj2.eta() << ", " << sj2.phi_std() << "\n";
    
    if (sj1.pt() < ptMinSubjet || sj2.pt() < ptMinSubjet) return false;

    theta = sj1.delta_R(sj2);
    const double pt1 = sj1.pt();
    const double pt2 = sj2.pt();
    z = std::min(pt1, pt2) / (pt1 + pt2);

    return (std::isfinite(theta) && std::isfinite(z));
  }
} // anon

// ------------------ public API ------------------
bool ComputeZsjForJet(
  Jet* jet,
  TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  TowerBackground* bg,
  float v2, float psi2, bool doUnsub,           // NOTE: float here
  double etaCalMax, double R_jet, double /*z_cut*/, double /*beta*/, double pt_min_subjet,
  double& theta_sj, double& z_sj)
{
  theta_sj = std::numeric_limits<double>::quiet_NaN();
  z_sj     = std::numeric_limits<double>::quiet_NaN();

  // Build constituents from all three calos with geometry
  std::vector<PseudoJet> parts;
  BuildPseudoJetsFromRecoJet_AllCalos(
      jet, em, geomEM, ih, geomIH, oh, geomOH, bg,
      v2, psi2, doUnsub, parts);

  if (parts.empty()) return false;

  JetDefinition jd(antikt_algorithm, R_jet);
  ClusterSequence cs(parts, jd);
  auto jets = sorted_by_pt(cs.inclusive_jets());
  if (jets.empty()) return false;

  const PseudoJet& j0 = jets.front();
  if (std::abs(j0.eta()) > etaCalMax) return false;

  // AKT R=0.1 subjets; default 3 GeV for both reco & truth (match EMJetVal)
  return SubjetZTheta_AKT01(j0.constituents(), pt_min_subjet, theta_sj, z_sj);
}

// --- Short wrapper (what JetMatchingSubjets.cc currently calls) ---
bool ComputeZsjForJet(
  Jet* jet,
  TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  TowerBackground* bg,
  double R_jet, double etaCalMax,
  double& z_sj, double& theta_sj)
{
  double theta_tmp = std::numeric_limits<double>::quiet_NaN();
  double z_tmp     = std::numeric_limits<double>::quiet_NaN();

  const bool ok = ComputeZsjForJet(
      jet, em, ih, oh, geomEM, geomIH, geomOH, bg,
      /*v2=*/0.f, /*psi2=*/0.f, /*doUnsub=*/false,
      /*etaCalMax=*/etaCalMax, /*R_jet=*/R_jet,
      /*z_cut=*/0.0, /*beta=*/0.0, /*pt_min_subjet=*/0.0,
      /*theta_out=*/theta_tmp, /*z_out=*/z_tmp);

  if (ok) { z_sj = z_tmp; theta_sj = theta_tmp; }
  return ok;
}

// --- Truth (use same 3 GeV subjet threshold to mirror reco/EMJetVal) ---
bool ComputeZsjForTruthJet(
  Jet* truthJet,
  PHG4TruthInfoContainer* truthInfo,
  double R_jet, double etaCalMax,
  double& z_sj, double& theta_sj)
{
  z_sj = std::numeric_limits<double>::quiet_NaN();
  theta_sj = std::numeric_limits<double>::quiet_NaN();
  if (!truthJet || !truthInfo) return false;

  std::vector<PseudoJet> parts;
  parts.reserve(truthJet->size_comp());
  for (const auto& c : truthJet->get_comp_vec()) {
    const int id = static_cast<int>(c.second);
    if (auto* p = truthInfo->GetParticle(id)) {
      parts.emplace_back(p->get_px(), p->get_py(), p->get_pz(), p->get_e());
    }
  }
  if (parts.empty()) return false;

  JetDefinition jd(antikt_algorithm, R_jet);
  ClusterSequence cs(parts, jd);
  auto jets = sorted_by_pt(cs.inclusive_jets());
  if (jets.empty()) return false;

  const PseudoJet& j0 = jets.front();
  if (std::abs(j0.eta()) > etaCalMax) return false;

  return SubjetZTheta_AKT01(j0.constituents(), /*ptMin=*/0.0, theta_sj, z_sj);
}

// keep linker happy if referenced elsewhere
extern "C" void ZSJ_DebugPrintAndReset() {}
