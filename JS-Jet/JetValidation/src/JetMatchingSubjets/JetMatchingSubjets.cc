#include "JetMatchingSubjets.h"
#include "ZsjTools.h"

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

  // Fetch nodes
  auto* jetsReco  = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  auto* jetsTruth = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  auto* cent      = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  auto* em   = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
  auto* ih   = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
  auto* oh   = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
  auto* geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  auto* geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto* geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  auto* bg    = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub1");
  auto* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!jetsReco || !jetsTruth || !cent || !em || !ih || !oh || !geomEM || !geomIH || !geomOH || !bg || !truthinfo) {
    std::cerr << "[JetMatchingSubjets] Missing node(s)." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_centrality = cent->get_centile(CentralityInfo::PROP::bimp);
  m_b          = cent->get_quantity(CentralityInfo::PROP::bimp);

  // Reco jets
  for (auto j: *jetsReco) {
    if (j->get_pt() < m_recoPtMin) continue;
    if (std::abs(j->get_eta()) > m_etaRange.second) continue;

    v_reco_pt.push_back(j->get_pt());
    v_reco_eta.push_back(j->get_eta());
    v_reco_phi.push_back(j->get_phi());

    double z, th;
    if (ComputeZsjForJet(j, em, ih, oh, geomEM, geomIH, geomOH, bg, 0.0, 0.0,
                         /*doUnsub=*/false, 0.4, 0.1, 3.0, 1.1, z, th)) {
      v_reco_zsj.push_back(z);
      v_reco_theta.push_back(th);
    } else {
      v_reco_zsj.push_back(NAN);
      v_reco_theta.push_back(NAN);
    }
  }
  
  // Truth jets
  for (auto j: *jetsTruth) {
    if (j->get_pt() < m_truthPtMin) continue;
    if (std::abs(j->get_eta()) > m_etaRange.second) continue;
    
    v_truth_pt.push_back(j->get_pt());
    v_truth_eta.push_back(j->get_eta());
    v_truth_phi.push_back(j->get_phi());
    
    double z, th;
    if (ComputeZsjForTruthJet(j, truthinfo, 0.4, 0.1, 3.0, 1.1, z, th)) {
      v_truth_zsj.push_back(z);
      v_truth_theta.push_back(th);
    } else {
      v_truth_zsj.push_back(NAN);
      v_truth_theta.push_back(NAN);
    }
  }
  
  // Matching
  MatchJets1to1(jetsReco, jetsTruth);
  
  // Fill matched
  for (auto& kv: recoToTruth) {
    Jet* r = kv.first;
    Jet* t = kv.second;
    
    const float dr = DeltaR(r,t);
    
    v_match_reco_pt.push_back(r->get_pt());
    v_match_reco_eta.push_back(r->get_eta());
    v_match_reco_phi.push_back(r->get_phi());
    v_match_truth_pt.push_back(t->get_pt());
    v_match_truth_eta.push_back(t->get_eta());
    v_match_truth_phi.push_back(t->get_phi());
    v_match_dR.push_back(dr);
    
    double zr, tr, zt, tt;
    if (!ComputeZsjForJet(r, em, ih, oh, geomEM, geomIH, geomOH, bg, 0.0, 0.0, false, 0.4, 0.1, 3.0, 1.1, zr, tr)) {
      zr=NAN; tr=NAN;
    }
    if (!ComputeZsjForTruthJet(t, truthinfo, 0.4, 0.1, 3.0, 1.1, zt, tt)) {
      zt=NAN; tt=NAN;
    }
    v_match_reco_zsj.push_back(zr);
    v_match_reco_theta.push_back(tr);
    v_match_truth_zsj.push_back(zt);
    v_match_truth_theta.push_back(tt);

    // Fakes (reco-only)
    for (auto r: *jetsReco) {
      if (r->get_pt() < m_recoPtMin) continue;
      if (std::abs(r->get_eta()) > m_etaRange.second) continue;
      if (recoToTruth.find(r) != recoToTruth.end()) continue;
      
      double zr, tr;
      if (!ComputeZsjForJet(r, em, ih, oh, geomEM, geomIH, geomOH, bg, 0.0, 0.0, false, 0.4, 0.1, 3.0, 1.1, zr, tr)) {
	zr=NAN; tr=NAN;
      }
      v_fake_reco_pt.push_back(r->get_pt());
      v_fake_reco_eta.push_back(r->get_eta());
      v_fake_reco_phi.push_back(r->get_phi());
      v_fake_reco_zsj.push_back(zr);
      v_fake_reco_theta.push_back(tr);
    }
  } 
  // Misses (truth-only)
  for (auto t: *jetsTruth) {
    if (t->get_pt() < m_truthPtMin) continue;
    if (std::abs(t->get_eta()) > m_etaRange.second) continue;
    if (truthToReco.find(t) != truthToReco.end()) continue;
    
    double zt, tt;
    if (!ComputeZsjForTruthJet(t, truthinfo, 0.4, 0.1, 3.0, 1.1, zt, tt)) {
      zt=NAN; tt=NAN;
    }
    v_fake_truth_pt.push_back(t->get_pt());
    v_fake_truth_eta.push_back(t->get_eta());
    v_fake_truth_phi.push_back(t->get_phi());
    v_fake_truth_zsj.push_back(zt);
    v_fake_truth_theta.push_back(tt);
  }
  
  m_T->Fill();
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
    usedR.insert(p.r);
    usedT.insert(p.t);
    recoToTruth[p.r] = p.t;
    truthToReco[p.t] = p.r;
     
  }
}

