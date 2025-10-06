#include "SubjetTools.h"
#include "MatchingUtils.h"

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
  
  static inline const RawTowerGeom*
  get_geom_flexible(int calo_id, int ieta, int iphi,
		    RawTowerGeomContainer* geomEM,
		    RawTowerGeomContainer* geomIH,
		    RawTowerGeomContainer* geomOH)
  {
    // calo_id: 0=CEMC, 1=HCALIN, 2=HCALOUT
    RawTowerGeomContainer* g = (calo_id==0? geomEM : (calo_id==1? geomIH : geomOH));
    if (!g) return nullptr;
    const RawTowerDefs::CalorimeterId cid = (calo_id==0? RawTowerDefs::CEMC :
					     calo_id==1? RawTowerDefs::HCALIN : RawTowerDefs::HCALOUT);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(cid, ieta, iphi);
    return g->get_tower_geometry(key);
  }
 
  /////// RECO TOWER BUILDER //////////////////
  static std::vector<fastjet::PseudoJet>
  BuildRecoConstituentsFromJet(const Jet* rjet,
			       TowerInfoContainer* cemc, TowerInfoContainer* hcalin, TowerInfoContainer* hcalout,
			       RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
			       TowerBackground* bg,
			       bool doUnsub, float v2 = 0.f, float psi2 = 0.f)
  {
    std::vector<fastjet::PseudoJet> out;
    if (!rjet) return out;
    
    auto handle_calo = [&](TowerInfoContainer* cont, int calo_id){
      return [=](unsigned int channel)->fastjet::PseudoJet {
	TowerInfo* tw = cont->get_tower_at_channel(channel);
	if (!tw) return fastjet::PseudoJet();
	const unsigned int key = cont->encode_key(channel);
	const int ieta = cont->getTowerEtaBin(key);
	const int iphi = cont->getTowerPhiBin(key);
	const RawTowerGeom* geo = get_geom_flexible(calo_id, ieta, iphi, geomEM, geomIH, geomOH);
	if (!geo) return fastjet::PseudoJet();
	
	/*
	// UE (off for pp)
	float UE = 0.f;
	if (doUnsub && bg) {
        // 0=CEMC,1=HCALIN,2=HCALOUT
        UE = bg->get_UE(calo_id).at(ieta);
        const float phi = geo->get_phi();
        UE *= (1.f + 2.f * v2 * std::cos(2.f * (phi - psi2)));
	}
	*/
	
	const double E  = tw->get_energy();
	const double eta = geo->get_eta();
	//	std::cout << "eta in towers --------------------- " << eta << std::endl;
	const double phi = geo->get_phi();
	//	std::cout << "phi in towers --------------------- " << phi << std::endl;
	const double pt  = E / std::cosh(eta);
	const double px = pt * std::cos(phi);
	const double py = pt * std::sin(phi);
	const double pz = pt * std::sinh(eta);
	
	return fastjet::PseudoJet(px, py, pz, E);
      };
    };
    
    auto mk_em   = handle_calo(cemc,   1);  //calo_id 0 is EMCAL but setting it 1 gives it the ihcal geometry (aka retowering)
    auto mk_ih   = handle_calo(hcalin, 1);  //calo_id 1
    auto mk_oh   = handle_calo(hcalout,2);  //calo_id 2
    
    // Accept both unsub/sub component codes: CEMC RETOWER(13,28) CEMC(14,29) HCALIN(15,30) HCALOUT(16,31)
    auto* nonconst = const_cast<Jet*>(rjet);
    for (const auto& comp : nonconst->get_comp_vec()) {
      const int   code = comp.first;
      const auto  chan = comp.second;
      
      fastjet::PseudoJet pj;
      if      (code == 13 || code == 28) pj = mk_em(chan);
      else if (code == 15 || code == 30) pj = mk_ih(chan);
      else if (code == 16 || code == 31) pj = mk_oh(chan);
      if (pj.E() > 0.0 && std::isfinite(pj.E())) out.push_back(pj);
    }
    
    return out;
  }
  // loop all channels in a TowerInfoContainer, keep towers in ΔR to (η0,φ0)
  static void collect_in_cone(
			      TowerInfoContainer* cont, int calo_id,
			      double eta0, double phi0, double R_jet,
			      RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
			      std::vector<fastjet::PseudoJet>& out)
  {
    if (!cont) return;
    const size_t nch = cont->size();              // TowerInfoContainer supports size()
    for (size_t ch = 0; ch < nch; ++ch)
      {
	TowerInfo* tw = cont->get_tower_at_channel(ch);
	if (!tw) continue;
	
	const unsigned int key = cont->encode_key(ch);
	const int ieta = cont->getTowerEtaBin(key);
	const int iphi = cont->getTowerPhiBin(key);
	
	// geometry: for retowered CEMC use HCALIN geom (calo_id=1)
	const RawTowerGeom* geo = get_geom_flexible(calo_id, ieta, iphi, geomEM, geomIH, geomOH);
	if (!geo) continue;
	
	const double eta = geo->get_eta();
	const double phi = geo->get_phi();
	const double dR  = calculateDeltaR(eta0, phi0, eta, phi);
	if (dR >= R_jet) continue;
	
	const double E  = tw->get_energy();
	if (!(std::isfinite(E) && E>0)) continue;
	const double pt = E / std::cosh(eta);
	const double px = pt * std::cos(phi);
	const double py = pt * std::sin(phi);
	const double pz = pt * std::sinh(eta);
	out.emplace_back(px, py, pz, E);
      }
  }
  
  /////// TRUTH BUILDER //////////////////
  // Build truth four-vectors for *this matched truth jet* from its stored particle components
  static bool BuildTruthConstituentsForJet(
					   Jet* truthJet,
					   PHG4TruthInfoContainer* truthinfo,
					   std::vector<fastjet::PseudoJet>& out)
  {
    out.clear();
    if (!truthJet || !truthinfo) return false;
    
    // Each component should correspond to a truth particle track id.
    // We retrieve the particle and turn it into a PseudoJet.
    for (const auto& comp : truthJet->get_comp_vec()) {
      const unsigned int trackid = comp.second;
      
      // Common sPHENIX API: GetParticle(trackid)
      PHG4Particle* p = truthinfo->GetParticle(trackid);
      if (!p) continue;
      
      const double px = p->get_px();
      const double py = p->get_py();
      const double pz = p->get_pz();
      double       E  = 0.0;
      
      // If energy is stored, use it; otherwise compute a massless E
      // (PHG4Particle typically has get_e(); if not, fall back to |p|).
#if defined(HAVE_PHG4PARTICLE_GET_E) || 1
      // Most builds have get_e()
      E = p->get_e();
      if (!std::isfinite(E) || E <= 0.0) {
        E = std::sqrt(px*px + py*py + pz*pz); // massless fallback
      }
#else
      E = std::sqrt(px*px + py*py + pz*pz);   // massless fallback
#endif
      
      // Skip non-finite or zero 4-vectors
      if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
	  !std::isfinite(E)  || E <= 0.0) {
	continue;
      }
      
      out.emplace_back(px, py, pz, E);
    }
    
    return !out.empty();
  }
  
  /////////// SUBJET CALCULATION/DEFINITION ///////////////////
  static bool SubjetsFromConstituents_R01(const std::vector<fastjet::PseudoJet>& constits,
					  double ptMin,
					  double& z_sj, double& theta_sj)
  {
    z_sj = std::numeric_limits<double>::quiet_NaN();
    theta_sj = std::numeric_limits<double>::quiet_NaN();
    
    if (constits.empty()) return false;
    
    fastjet::JetDefinition def01(fastjet::antikt_algorithm, 0.1);
    fastjet::ClusterSequence seq01(constits, def01);
    auto sj = fastjet::sorted_by_pt(seq01.inclusive_jets());
    if (sj.size() < 2) return false;
    
    const auto& a = sj[0];
    const auto& b = sj[1];
    if (a.pt() < ptMin || b.pt() < ptMin) return false;
    
    theta_sj = a.delta_R(b);
    if (theta_sj > 0.4) return false;
    
    const double pt1 = a.pt(), pt2 = b.pt();
    const double sum = pt1 + pt2;
    if (sum <= 0) return false;
    z_sj = std::min(pt1, pt2) / sum;
    if (z_sj < 0.1) return false; //z_cut
    return true;
  }
  static bool LeadingPtsFromConstituents_R01(
					     const std::vector<fastjet::PseudoJet>& constits,
					     double ptMin, double& pt1, double& pt2)
  {
    pt1 = pt2 = std::numeric_limits<double>::quiet_NaN();
    if (constits.empty()) return false;
    
    fastjet::JetDefinition def01(fastjet::antikt_algorithm, 0.1);
    fastjet::ClusterSequence seq01(constits, def01);
    auto sj = fastjet::sorted_by_pt(seq01.inclusive_jets());
    if (sj.size() < 2) return false;
    
    const auto& a = sj[0];
    const auto& b = sj[1];
    if (a.pt() < ptMin || b.pt() < ptMin) return false;
    
    pt1 = a.pt(); pt2 = b.pt();
    return true;
  }
} // anon
//////////// COMPUTE SUBJET OBSERVABLES FOR RECO ////////////////////
bool ComputeSubjetForJet(
			 Jet* rraw,
			 TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
			 RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
			 TowerBackground* bg,
			 float v2, float psi2, bool doUnsub,
			 double etaCalMax, double R_jet,
			 double z_cut, double beta, double pt_min_subjet,
			 double& z_sj, double& theta_sj)
{
  // build constituents from *this* jet only (all calos)
  auto constits = BuildRecoConstituentsFromJet(rraw, em, ih, oh, geomEM, geomIH, geomOH, bg, doUnsub, v2, psi2);
  if (constits.empty()) return false;
  
  // ensure parent R=0.4 jet is fiduciary (optional: check rraw->get_eta() against etaCalMax)
  if (std::abs(rraw->get_eta()) > etaCalMax) return false;
  
  return SubjetsFromConstituents_R01(constits, pt_min_subjet, z_sj, theta_sj);
}

//////////// COMPUTE SUBJET OBSERVABLES FOR TRUTH ////////////////////
bool ComputeSubjetForTruthJet(
			   Jet* truthJet, PHG4TruthInfoContainer* truthInfo,
			   double R_jet, double etaCalMax,
			   double pt_min_subjet,
			   double& z_sj, double& theta_sj)
{
  z_sj = std::numeric_limits<double>::quiet_NaN();
  theta_sj = std::numeric_limits<double>::quiet_NaN();
  
  if (!truthJet || !truthInfo) return false;
  if (std::abs(truthJet->get_eta()) > etaCalMax) return false;
  
  std::vector<fastjet::PseudoJet> constits;
  if (!BuildTruthConstituentsForJet(truthJet, truthInfo, constits)) return false;
  
  // Small-R reclustering of THIS jet's constituents: anti-kT, R=0.2
  return SubjetsFromConstituents_R01(constits, pt_min_subjet, z_sj, theta_sj);
}

bool LeadingSubjetPts_RecoJet(
			      Jet* rraw,
			      TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
			      RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
			      TowerBackground* bg,
			      float v2, float psi2, bool doUnsub,
			      double etaCalMax, double R_jet, double ptMinSubjet,
			      double& pt1, double& pt2)
{
  // Build reco constituents the same way as ComputeSubjetForJet
  auto constits = BuildRecoConstituentsFromJet(
					       rraw, em, ih, oh, geomEM, geomIH, geomOH, bg, doUnsub, v2, psi2);
  if (constits.empty()) return false;
  if (std::abs(rraw->get_eta()) > etaCalMax) return false;
  return LeadingPtsFromConstituents_R01(constits, ptMinSubjet, pt1, pt2);
}

bool LeadingSubjetPts_TruthJet(
			       Jet* truthJet, PHG4TruthInfoContainer* truthInfo,
			       double R_jet, double etaCalMax, double ptMinSubjet,
			       double& pt1, double& pt2)
{
  if (!truthJet || !truthInfo) return false;
  if (std::abs(truthJet->get_eta()) > etaCalMax) return false;
  
  std::vector<fastjet::PseudoJet> constits;
  if (!BuildTruthConstituentsForJet(truthJet, truthInfo, constits)) return false;
  return LeadingPtsFromConstituents_R01(constits, ptMinSubjet, pt1, pt2);
}

bool ComputeSubjetForTruthJet_FromTowersCone(
  Jet* truthJet,
  TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  double R_jet, double etaCalMax, double pt_min_subjet,
  double& z_sj, double& theta_sj)
{
  z_sj = theta_sj = std::numeric_limits<double>::quiet_NaN();
  if (!truthJet) return false;
  if (std::abs(truthJet->get_eta()) > etaCalMax) return false;

  std::vector<fastjet::PseudoJet> constits;
  // IMPORTANT: for CEMC we feed calo_id=1 to use HCALIN geometry (retower layout)
  collect_in_cone(em, /*calo_id=*/1, truthJet->get_eta(), truthJet->get_phi(), R_jet,
                  geomEM, geomIH, geomOH, constits);
  collect_in_cone(ih, /*calo_id=*/1, truthJet->get_eta(), truthJet->get_phi(), R_jet,
                  geomEM, geomIH, geomOH, constits);
  collect_in_cone(oh, /*calo_id=*/2, truthJet->get_eta(), truthJet->get_phi(), R_jet,
                  geomEM, geomIH, geomOH, constits);

  if (constits.empty()) return false;
  return SubjetsFromConstituents_R01(constits, pt_min_subjet, z_sj, theta_sj);
}

bool LeadingSubjetPts_TruthJet_FromTowersCone(
  Jet* truthJet,
  TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  double R_jet, double etaCalMax, double ptMinSubjet,
  double& pt1, double& pt2)
{
  pt1 = pt2 = std::numeric_limits<double>::quiet_NaN();
  if (!truthJet) return false;
  if (std::abs(truthJet->get_eta()) > etaCalMax) return false;

  std::vector<fastjet::PseudoJet> constits;
  collect_in_cone(em, /*calo_id=*/1, truthJet->get_eta(), truthJet->get_phi(), R_jet,
                  geomEM, geomIH, geomOH, constits);
  collect_in_cone(ih, /*calo_id=*/1, truthJet->get_eta(), truthJet->get_phi(), R_jet,
                  geomEM, geomIH, geomOH, constits);
  collect_in_cone(oh, /*calo_id=*/2, truthJet->get_eta(), truthJet->get_phi(), R_jet,
                  geomEM, geomIH, geomOH, constits);

  if (constits.empty()) return false;
  return LeadingPtsFromConstituents_R01(constits, ptMinSubjet, pt1, pt2);
}
