#include "JetMatchingSubjets.h"
#include "ZsjTools.h"         // <- Calculates zsj and theta_sj
#include "SubjetMatching.h"   // <— Matches subjets

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <jetbase/JetContainer.h>
#include <jetbase/Jet.h>

#include <centrality/CentralityInfo.h>
#include <jetbackground/TowerBackground.h>

#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>

#include <TH1D.h>
#include <TTree.h>
#include <RooUnfoldResponse.h>
#include <TParameter.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

#include <algorithm>  // std::upper_bound
#include <array>

extern "C" void ZSJ_DebugPrintAndReset();
namespace {
  template<class T> T* book1D(const char* name, const char* title, const std::vector<double>& edges) {
    return new T(name, title, static_cast<int>(edges.size()-1), edges.data());
  }
}

JetMatchingSubjets::JetMatchingSubjets(const std::string& recojetname,
                                       const std::string& truthjetname,
                                       const std::string& outputfilename)
  : SubsysReco(std::string("JetMatchingSubjets_")+recojetname+"_"+truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
{}

JetMatchingSubjets::~JetMatchingSubjets() {}

int JetMatchingSubjets::Init(PHCompositeNode*)
{
  
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  
  // Tree
  m_T = new TTree("T", "Jet matching for RooUnfold (pt and pt×z_sj)");
  m_T->Branch("event", &m_event, "event/I");
  m_T->Branch("cent",  &m_centrality, "cent/F");
  m_T->Branch("b",     &m_b, "b/F");
  
  m_T->Branch("reco_pt",    &v_reco_pt);
  m_T->Branch("reco_eta",   &v_reco_eta);
  m_T->Branch("reco_phi",   &v_reco_phi);
  m_T->Branch("reco_zsj",   &v_reco_zsj);
  m_T->Branch("reco_theta", &v_reco_theta);
  
  m_T->Branch("truth_pt",    &v_truth_pt);
  m_T->Branch("truth_eta",   &v_truth_eta);
  m_T->Branch("truth_phi",   &v_truth_phi);
  m_T->Branch("truth_zsj",   &v_truth_zsj);
  m_T->Branch("truth_theta", &v_truth_theta);

  // matched pairs
  m_T->Branch("match_reco_pt",    &v_match_reco_pt);
  m_T->Branch("match_reco_eta",   &v_match_reco_eta);
  m_T->Branch("match_reco_phi",   &v_match_reco_phi);
  m_T->Branch("match_reco_zsj",   &v_match_reco_zsj);
  m_T->Branch("match_reco_theta", &v_match_reco_theta);
  m_T->Branch("match_truth_pt",   &v_match_truth_pt);
  m_T->Branch("match_truth_eta",  &v_match_truth_eta);
  m_T->Branch("match_truth_phi",  &v_match_truth_phi);
  m_T->Branch("match_truth_zsj",  &v_match_truth_zsj);
  m_T->Branch("match_truth_theta",&v_match_truth_theta);
  m_T->Branch("match_dR",         &v_match_dR);
  /*
  // matched reco subjet observables
  m_T->Branch("match_reco_zsj_20_25", &v_match_reco_zsj_20_25);
  m_T->Branch("match_reco_theta_20_25", &v_match_reco_theta_20_25);
  m_T->Branch("match_reco_zsj_25_30", &v_match_reco_zsj_25_30);
  m_T->Branch("match_reco_theta_25_30", &v_match_reco_theta_25_30);

  // matched truth subjet observables
  m_T->Branch("match_truth_zsj_20_25", &v_match_truth_zsj_20_25);
  m_T->Branch("match_truth_theta_20_25", &v_match_truth_theta_20_25);
  m_T->Branch("match_truth_zsj_25_30", &v_match_truth_zsj_25_30);
  m_T->Branch("match_truth_theta_25_30", &v_match_truth_theta_25_30);
  */
  // fakes & misses
  m_T->Branch("fake_reco_pt",    &v_fake_reco_pt);
  m_T->Branch("fake_reco_eta",   &v_fake_reco_eta);
  m_T->Branch("fake_reco_phi",   &v_fake_reco_phi);
  m_T->Branch("fake_reco_zsj",   &v_fake_reco_zsj);
  m_T->Branch("fake_reco_theta", &v_fake_reco_theta);
  m_T->Branch("fake_truth_pt",    &v_fake_truth_pt);
  m_T->Branch("fake_truth_eta",   &v_fake_truth_eta);
  m_T->Branch("fake_truth_phi",   &v_fake_truth_phi);
  m_T->Branch("fake_truth_zsj",   &v_fake_truth_zsj);
  m_T->Branch("fake_truth_theta", &v_fake_truth_theta);

  // SubjetMatching.h
  m_T->Branch("match_sj", &v_match_sj);
  m_T->Branch("match_jet_and_sj", &v_match_jet_and_sj);
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingSubjets::process_event(PHCompositeNode* topNode)
{
  ++m_event;
  
  // clear vectors
  v_reco_pt.clear();  v_reco_eta.clear();  v_reco_phi.clear();  v_reco_zsj.clear();  v_reco_theta.clear();
  v_truth_pt.clear(); v_truth_eta.clear(); v_truth_phi.clear(); v_truth_zsj.clear(); v_truth_theta.clear();
  v_match_reco_pt.clear(); v_match_reco_eta.clear(); v_match_reco_phi.clear(); v_match_reco_zsj.clear(); v_match_reco_theta.clear();
  v_match_truth_pt.clear(); v_match_truth_eta.clear(); v_match_truth_phi.clear(); v_match_truth_zsj.clear(); v_match_truth_theta.clear();
  v_match_dR.clear();
  v_fake_reco_pt.clear(); v_fake_reco_eta.clear(); v_fake_reco_phi.clear(); v_fake_reco_zsj.clear(); v_fake_reco_theta.clear();
  v_fake_truth_pt.clear(); v_fake_truth_eta.clear(); v_fake_truth_phi.clear(); v_fake_truth_zsj.clear(); v_fake_truth_theta.clear();
  
  recoToTruth.clear();
  truthToReco.clear();
  v_match_sj.clear();
  v_match_jet_and_sj.clear();
  
  // Fetch nodes
  auto* jetsReco  = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  auto* jetsTruth = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  auto* cent      = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  
  TowerInfoContainer *em = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  TowerInfoContainer *ih = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *oh = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  
  // Background
  auto* bg = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");
  auto* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  
  if (!jetsReco || !jetsTruth || !cent || !em || !ih || !oh || !geomEM || !geomIH || !geomOH || !bg || !truthinfo) {
    std::cerr << "[JetMatchingSubjets] Missing node(s)." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  m_centrality = cent->get_centile(CentralityInfo::PROP::bimp);
  m_b          = cent->get_quantity(CentralityInfo::PROP::bimp);
  /*
    auto wrapPhi = [](double a){
    while (a >  M_PI) a -= 2*M_PI;
    while (a < -M_PI) a += 2*M_PI;
    return a;
    };
  */ 
  // Reco jets
  for (auto j: *jetsReco) {
    if (j->get_pt() < m_recoPtMin) continue; 
    if (std::abs(j->get_eta()) > m_etaRange.second) continue;
    
    v_reco_pt.push_back(j->get_pt());
    v_reco_eta.push_back(j->get_eta());
    v_reco_phi.push_back(j->get_phi());
    std::cout << __LINE__ << std::endl;
    auto& v = v_reco_phi;
    if (v.back() < -M_PI) v.back() += 2*M_PI;
    if (v.back() >  M_PI) v.back() -= 2*M_PI;
    std::cout << __LINE__ << std::endl;
    double z = std::numeric_limits<double>::quiet_NaN();
    double th = std::numeric_limits<double>::quiet_NaN();
    std::cout << __LINE__ << std::endl;
    const bool ok = ComputeZsjForJet(j,
				     em, ih, oh, geomEM, geomIH, geomOH, bg,
				     /*v2=*/0.f, /*psi2=*/0.f, /*doUnsub=*/false,
				     /*etaCalMax=*/1.1, /*R_jet=*/0.4,
				     /*z_cut=*/0.0, /*beta=*/0.0, /*pt_min_subjet=*/3.0,
				     th, z);
    std::cout << __LINE__ << std::endl;
    v_reco_zsj.push_back(ok ? z  : std::numeric_limits<double>::quiet_NaN());
    v_reco_theta.push_back(ok ? th : std::numeric_limits<double>::quiet_NaN());
    std::cout << __LINE__ << std::endl;
  }
  
  // Truth jets
  for (auto j: *jetsTruth) {
    if (j->get_pt() < m_truthPtMin) continue;
    if (std::abs(j->get_eta()) > m_etaRange.second) continue;
    std::cout << __LINE__ << std::endl;
    v_truth_pt.push_back(j->get_pt());
    v_truth_eta.push_back(j->get_eta());
    v_truth_phi.push_back(j->get_phi());
    std::cout << __LINE__ << std::endl;
    auto& v = v_truth_phi;
    if (v.back() < -M_PI) v.back() += 2*M_PI;
    if (v.back() >  M_PI) v.back() -= 2*M_PI;
    std::cout << __LINE__ << std::endl;
    double z, th;
    std::cout << __LINE__ << std::endl;
    if (ComputeZsjForTruthJet(j, truthinfo, 0.4, 1.1, z, th)) {
      v_truth_zsj.push_back(z);
      v_truth_theta.push_back(th);
      std::cout << __LINE__ << std::endl;
    } else {
      v_truth_zsj.push_back(NAN);
      v_truth_theta.push_back(NAN);
    }
  }
  std::cout << __LINE__ << std::endl;
  // Matching
  MatchJets1to1(jetsReco, jetsTruth);
  // ---- Fill matched (only keep zsj when BOTH reco & truth have valid subjets) ----
  for (auto& kv: truthToReco) {
    Jet* r = kv.second;
    Jet* t = kv.first;
    const float dr = DeltaR(r,t);
    std::cout << __LINE__ << std::endl;
    // Compute reco subjets (AKT R=0.1, pT>5), acceptance |eta|<0.6
    double zr = std::numeric_limits<double>::quiet_NaN();
    double tr = std::numeric_limits<double>::quiet_NaN();
    const bool okr = ComputeZsjForJet(
				      r,
				      em, ih, oh, geomEM, geomIH, geomOH, bg,
				      /*v2=*/0.f, /*psi2=*/0.f, /*doUnsub=*/false,
				      /*etaCalMax=*/1.1, /*R_jet=*/0.4,
				      /*z_cut=*/0.0, /*beta=*/0.0, /*pt_min_subjet=*/3.0,
				      tr, zr);
    std::cout << __LINE__ << std::endl;
    // Compute truth subjets (AKT R=0.1, pT>5), acceptance |eta|<0.6
    double zt = std::numeric_limits<double>::quiet_NaN();
    double tt = std::numeric_limits<double>::quiet_NaN();
    std::cout << __LINE__ << std::endl;
    const bool okt = ComputeZsjForTruthJet(
					   t, truthinfo,
					   /*R_jet=*/0.4, /*etaCalMax=*/1.1,
					   zt, tt);
    std::cout << __LINE__ << std::endl;
    /*
    std::vector<fastjet::PseudoJet> reco_sj, truth_sj;
    const double SJ_PT_MIN   = 3;   // analysis choice; keep angular-only matching policy
    const double SJ_DR_MAX   = 0.12;  //  subjet ΔR
    std::cout << __LINE__ << std::endl;
    bool haveRecoSJ  = BuildRecoSubjetsForJet(r, em, ih, oh, geomEM, geomIH, geomOH,
					      bg, 1.1, 0.4,
					      SJ_PT_MIN, reco_sj);
    bool haveTruthSJ = BuildTruthSubjetsForJet(t, truthinfo,
					       0.4, 1.1,
					       SJ_PT_MIN, truth_sj);
    std::cout << __LINE__ << std::endl;
    // 3) require BOTH sides to have ≥2 subjets and pass greedy matching
    bool subjet_gate = false;
    if (haveRecoSJ && haveTruthSJ) {
      auto sjres = MatchTruthTwoToRecoMany_Greedy(truth_sj, reco_sj, SJ_DR_MAX);
      subjet_gate = sjres.matched;
    }
    std::cout << __LINE__ << std::endl;
    v_match_sj.push_back(subjet_gate ? 1 : 0);
    v_match_jet_and_sj.push_back((subjet_gate) ? 1 : 0);
*/
    std::cout << __LINE__ << std::endl;
    // Only store as "matched" if BOTH are valid
    if (okr && okt /* && subjet_gate*/) {
      v_match_reco_pt.push_back(r->get_pt());
      v_match_reco_eta.push_back(r->get_eta());
      v_match_reco_phi.push_back(r->get_phi());
      std::cout << __LINE__ << std::endl;
      auto& v = v_match_reco_phi;
      if (v.back() < -M_PI) v.back() += 2*M_PI;
      if (v.back() >  M_PI) v.back() -= 2*M_PI;
      std::cout << __LINE__ << std::endl;
      v_match_reco_zsj.push_back(zr);
      v_match_reco_theta.push_back(tr);
      
      v_match_truth_pt.push_back(t->get_pt());
      v_match_truth_eta.push_back(t->get_eta());
      v_match_truth_phi.push_back(t->get_phi());
      std::cout << __LINE__ << std::endl;
      auto& b = v_match_truth_phi;
      if (b.back() < -M_PI) b.back() += 2*M_PI;
      if (b.back() >  M_PI) b.back() -= 2*M_PI;
      std::cout << __LINE__ << std::endl;
      v_match_truth_zsj.push_back(zt);
      v_match_truth_theta.push_back(tt);
      
      v_match_dR.push_back(dr);
      std::cout << __LINE__ << std::endl;
    } else {
      // Route partial/failed subjet cases to fakes
      if (okr) {
	v_fake_reco_pt.push_back(r->get_pt());
	v_fake_reco_eta.push_back(r->get_eta());
	v_fake_reco_phi.push_back(r->get_phi());
	std::cout << __LINE__ << std::endl;
	auto& v = v_fake_reco_phi;
	if (v.back() < -M_PI) v.back() += 2*M_PI;
	if (v.back() >  M_PI) v.back() -= 2*M_PI;
	std::cout << __LINE__ << std::endl;
	v_fake_reco_zsj.push_back(zr);
	v_fake_reco_theta.push_back(tr);
      }
      std::cout << __LINE__ << std::endl;
      if (okt) {
	v_fake_truth_pt.push_back(t->get_pt());
	v_fake_truth_eta.push_back(t->get_eta());
	v_fake_truth_phi.push_back(t->get_phi());
	std::cout << __LINE__ << std::endl;
	auto& v = v_fake_truth_phi;
	if (v.back() < -M_PI) v.back() += 2*M_PI;
	if (v.back() >  M_PI) v.back() -= 2*M_PI;
	std::cout << __LINE__ << std::endl;
	v_fake_truth_zsj.push_back(zt);
	v_fake_truth_theta.push_back(tt);
      }
      // If neither side had valid subjets, nothing goes to "matched".
      // Get a record of these via the fake_* entries above only if one side succeeded.
      // If neither succeeded here, they simply won't appear in match_* nor fake_* (both invalid).
    }
  }
  // Fakes (reco-only, outside matched loop)
  for (auto r: *jetsReco) {
    if (r->get_pt() < m_recoPtMin) continue;
    if (std::abs(r->get_eta()) > m_etaRange.second) continue;
    if (recoToTruth.find(r) != recoToTruth.end()) continue;
    std::cout << __LINE__ << std::endl;
    v_fake_reco_pt.push_back(r->get_pt());
    v_fake_reco_eta.push_back(r->get_eta());
    v_fake_reco_phi.push_back(r->get_phi());
    std::cout << __LINE__ << std::endl;
    auto& v = v_fake_reco_phi;
    if (v.back() < -M_PI) v.back() += 2*M_PI;
    if (v.back() >  M_PI) v.back() -= 2*M_PI;
    std::cout << __LINE__ << std::endl;
    double zr = std::numeric_limits<double>::quiet_NaN();
    double tr = std::numeric_limits<double>::quiet_NaN();
    std::cout << __LINE__ << std::endl;
    (void)ComputeZsjForJet(r,
			   em, ih, oh, geomEM, geomIH, geomOH, bg,
			   /*v2=*/0.f, /*psi2=*/0.f, /*doUnsub=*/false,
			   /*etaCalMax=*/1.1, /*R_jet=*/0.4,
			   /*z_cut=*/0.0, /*beta=*/0.0, /*pt_min_subjet=*/3.0,
			   tr, zr);
    
    v_fake_reco_zsj.push_back(zr);
    v_fake_reco_theta.push_back(tr);
  }
  std::cout << __LINE__ << std::endl;
  // Misses (truth-only)
  for (auto t: *jetsTruth) {
    if (t->get_pt() < m_truthPtMin) continue;
    if (std::abs(t->get_eta()) > m_etaRange.second) continue;
    if (truthToReco.find(t) != truthToReco.end()) continue;
    std::cout << __LINE__ << std::endl;
    double zt, tt;
    if (!ComputeZsjForTruthJet(t, truthinfo,
			       /*R_jet=*/0.4,
			       /*etaCalMax=*/1.1,
			       zt, tt)) {
      zt=NAN; tt=NAN;
    }
    std::cout << __LINE__ << std::endl;
    v_fake_truth_pt.push_back(t->get_pt());
    v_fake_truth_eta.push_back(t->get_eta());
    v_fake_truth_phi.push_back(t->get_phi());
    std::cout << __LINE__ << std::endl;

    auto& v = v_truth_phi;
    if (v.back() < -M_PI) v.back() += 2*M_PI;
    if (v.back() >  M_PI) v.back() -= 2*M_PI;
    std::cout << __LINE__ << std::endl;
    v_fake_truth_zsj.push_back(zt);
    v_fake_truth_theta.push_back(tt);
  }
  
  m_T->Fill();
  
  ZSJ_DebugPrintAndReset();

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingSubjets::End(PHCompositeNode*)
{
  PHTFileServer::get().cd(m_outputFileName);
  if (m_T) { m_T->Write(); m_T = nullptr; }
  // Store N(z-win) and sigma for this file
  TParameter<double>    pS("sigma_pb",  m_sigma_pb);
  TNamed                pName("sample_name", m_sample_name.c_str());
  TParameter<long long>("N_in_zwin", m_N_in_zwin).Write(); // number of entries written (|z|<=30)                                                                                                
  pS.Write();
  pName.Write();
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingSubjets::Reset(PHCompositeNode*) { return Fun4AllReturnCodes::EVENT_OK; }
void JetMatchingSubjets::Print(const std::string& what) const { }
  
// ---- helper definitions ----
static inline float wrapPhi(float dphi){
  while (dphi >  M_PI) dphi -= 2*M_PI;
  while (dphi < -M_PI) dphi += 2*M_PI;
  return dphi; 
}

float JetMatchingSubjets::DeltaR(Jet* a, Jet* b)
{
  const float dEta = a->get_eta() - b->get_eta(); 
  float dPhi = wrapPhi(a->get_phi() - b->get_phi());
  return std::sqrt(dEta*dEta + dPhi*dPhi);
}

void JetMatchingSubjets::MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets)
{
  struct Pair { float dr; Jet* r; Jet* t; };  
  std::vector<Pair> cand;
  for (auto r: *recoJets) {
    if (r->get_pt() < m_recoPtMin) continue;
    if (std::abs(r->get_eta()) > m_etaRange.second) continue;
    for (auto t: *truthJets) {
      if (t->get_pt() < m_truthPtMin) continue;
      if (std::abs(t->get_eta()) > m_etaRange.second) continue;
      float dr = DeltaR(r,t);
      if (dr < m_matchDRMax) cand.push_back({dr, r, t});
    }
  }
  std::sort(cand.begin(), cand.end(), [](const Pair& a, const Pair& b){ return a.dr < b.dr; });
  std::set<Jet*> usedR, usedT;
  recoToTruth.clear();
  truthToReco.clear();
  for (auto& p: cand) {
    if (usedR.count(p.r) || usedT.count(p.t)) continue;
    usedT.insert(p.t);
    usedR.insert(p.r);
    truthToReco[p.t] = p.r;
    recoToTruth[p.r] = p.t;     
  }
}

