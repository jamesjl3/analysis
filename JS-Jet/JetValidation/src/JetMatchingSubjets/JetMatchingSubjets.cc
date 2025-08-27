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

#include <TH1F.h>
#include <TH1D.h>
#include <TTree.h>
#include <RooUnfoldResponse.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <unordered_map>

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
  if (m_ptEdges.empty()) m_ptEdges = {5,10,15,20,25,30,35,40,45,50,55};
  m_nPtBins = (int)m_ptEdges.size() - 1;
  
  if (m_zsjEdges.empty()) { const int n=10; for (int i=0;i<=n;++i) m_zsjEdges.push_back(0.5*i/n); }
  m_nZsjBins = (int)m_zsjEdges.size() - 1;
  
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  

  hRecoPt  = new TH1D("hRecoPt",  "Reco jet pT; p_{T} [GeV]; entries",
                      m_nPtBins, m_ptEdges.data());
  //std::cout << __LINE__ << std::endl;
  hTruthPt = new TH1D("hTruthPt", "Truth jet pT; p_{T} [GeV]; entries",
                      m_nPtBins, m_ptEdges.data());
  //std::cout << __LINE__ << std::endl;

  hRecoPt->SetDirectory(nullptr);
  //std::cout << __LINE__ << std::endl;
  hTruthPt->SetDirectory(nullptr);
  //std::cout << __LINE__ << std::endl;
  const double* xb = hRecoPt ->GetXaxis()->GetXbins()->GetArray();
  const double* yb = hTruthPt->GetXaxis()->GetXbins()->GetArray();
  hRespTempl = new TH2D("hRespTempl",";reco p_{T};truth p_{T}",
			hRecoPt->GetNbinsX(),  xb,
			hTruthPt->GetNbinsX(), yb);
  hRespTempl->SetDirectory(nullptr);
  
  respPt = new RooUnfoldResponse(hRecoPt, hTruthPt, hRespTempl, "respPt", "RooUnfold pt response");
  std::cout << __LINE__ << std::endl;

  //respPt = new RooUnfoldResponse(hRecoPt, hTruthPt, "respPt", "RooUnfold pt response");
  
  // Flattened (pT × z_sj)
  const int nFlat = m_nPtBins * m_nZsjBins;
  hRecoPtZ  = new TH1D("hRecoPtZ",  "Reco (p_{T}×z_{sj}) flattened;global bin index;entries",
                       nFlat, 0, nFlat);
  //std::cout << __LINE__ << std::endl;
  hTruthPtZ = new TH1D("hTruthPtZ", "Truth (p_{T}×z_{sj}) flattened;global bin index;entries",
                       nFlat, 0, nFlat);
  //std::cout << __LINE__ << std::endl;
  hRecoPtZ->SetDirectory(nullptr);
  //std::cout << __LINE__ << std::endl;
  hTruthPtZ->SetDirectory(nullptr);
  //std::cout << __LINE__ << std::endl;
  
  respPtZ = new RooUnfoldResponse(hRecoPtZ, hTruthPtZ, "respPtZ",
                                  "RooUnfold flattened p_{T}×z_{sj} response");
  //std::cout << __LINE__ << std::endl;
  // Optional: helper histograms with the edges (not filled)
  hPtEdges  = new TH1D("hPtEdges",  "pt edges helper",  m_nPtBins,  m_ptEdges.data());
  //std::cout << __LINE__ << std::endl;
  hZsjEdges = new TH1D("hZsjEdges", "z_{sj} edges helper", m_nZsjBins, m_zsjEdges.data());
  //std::cout << __LINE__ << std::endl;
  hPtEdges->SetDirectory(nullptr);
  //std::cout << __LINE__ << std::endl;
  hZsjEdges->SetDirectory(nullptr);
  //std::cout << __LINE__ << std::endl;
  // Tree
  m_T = new TTree("T", "Jet matching for RooUnfold (pt and pt×z_sj)");
  m_T->Branch("event", &m_event, "event/I");
  m_T->Branch("cent",  &m_centrality, "cent/F");
  m_T->Branch("b",     &m_b, "b/F");

  m_T->Branch("reco_pt",  &v_reco_pt);
  m_T->Branch("reco_eta", &v_reco_eta);
  m_T->Branch("reco_phi", &v_reco_phi);
  m_T->Branch("reco_zsj", &v_reco_zsj);

  m_T->Branch("truth_pt",  &v_truth_pt);
  m_T->Branch("truth_eta", &v_truth_eta);
  m_T->Branch("truth_phi", &v_truth_phi);
  m_T->Branch("truth_zsj", &v_truth_zsj);

  // matched pairs
  m_T->Branch("match_reco_pt",  &v_match_reco_pt);
  m_T->Branch("match_reco_eta", &v_match_reco_eta);
  m_T->Branch("match_reco_phi", &v_match_reco_phi);
  m_T->Branch("match_reco_zsj", &v_match_reco_zsj);
  m_T->Branch("match_truth_pt",  &v_match_truth_pt);
  m_T->Branch("match_truth_eta", &v_match_truth_eta);
  m_T->Branch("match_truth_phi", &v_match_truth_phi);
  m_T->Branch("match_truth_zsj", &v_match_truth_zsj);
  m_T->Branch("match_dR",        &v_match_dR);

  // fakes & misses
  m_T->Branch("fake_reco_pt",  &v_fake_reco_pt);
  m_T->Branch("fake_reco_eta", &v_fake_reco_eta);
  m_T->Branch("fake_reco_phi", &v_fake_reco_phi);
  m_T->Branch("fake_reco_zsj", &v_fake_reco_zsj);
  m_T->Branch("fake_truth_pt",  &v_fake_truth_pt);
  m_T->Branch("fake_truth_eta", &v_fake_truth_eta);
  m_T->Branch("fake_truth_phi", &v_fake_truth_phi);
  m_T->Branch("fake_truth_zsj", &v_fake_truth_zsj);

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingSubjets::process_event(PHCompositeNode* topNode)
{
  ++m_event;

  // clear vectors
  v_reco_pt.clear();  v_reco_eta.clear();  v_reco_phi.clear();  v_reco_zsj.clear();
  v_truth_pt.clear(); v_truth_eta.clear(); v_truth_phi.clear(); v_truth_zsj.clear();
  v_match_reco_pt.clear(); v_match_reco_eta.clear(); v_match_reco_phi.clear(); v_match_reco_zsj.clear();
  v_match_truth_pt.clear(); v_match_truth_eta.clear(); v_match_truth_phi.clear(); v_match_truth_zsj.clear();
  v_match_dR.clear();
  v_fake_reco_pt.clear(); v_fake_reco_eta.clear(); v_fake_reco_phi.clear(); v_fake_reco_zsj.clear();
  v_fake_truth_pt.clear(); v_fake_truth_eta.clear(); v_fake_truth_phi.clear(); v_fake_truth_zsj.clear();

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
  auto* geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  auto* bg    = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub1");

  if (!jetsReco || !jetsTruth || !cent || !em || !ih || !oh || !geomEM || !geomOH || !bg) {
    std::cerr << "[JetMatchingSubjets] Missing node(s)." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_centrality = cent->get_centile(CentralityInfo::PROP::bimp);
  m_b          = cent->get_quantity(CentralityInfo::PROP::bimp);

  // Precompute z_sj maps for all jets 
  std::unordered_map<const Jet*, double> zReco, zTruth;

  // Reco jets
  for (auto j: *jetsReco) {
    if (j->get_pt() < m_recoPtMin) continue;
    if (std::abs(j->get_eta()) > m_etaRange.second) continue;

    v_reco_pt.push_back(j->get_pt());
    v_reco_eta.push_back(j->get_eta());
    v_reco_phi.push_back(j->get_phi());

    double z, th;
    if (ComputeZsjForJet(j, em, ih, oh, geomEM, geomOH, bg, /*v2=*/0, /*psi2=*/0, /*doUnsub=*/true,
                         /*R_jet=*/0.4, /*R_sub=*/0.1, /*ptCutSubjet=*/3.0, z, th)) {
      zReco[j] = z;
      v_reco_zsj.push_back(z);
    } else {
      v_reco_zsj.push_back(NAN);
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
    if (ComputeZsjForJet(j, em, ih, oh, geomEM, geomOH, bg, /*v2=*/0, /*psi2=*/0, /*doUnsub=*/true,
                         /*R_jet=*/0.4, /*R_sub=*/0.1, /*ptCutSubjet=*/3.0, z, th)) {
      zTruth[j] = z;
      v_truth_zsj.push_back(z);
    } else {
      v_truth_zsj.push_back(NAN);
    }
  }
  // Matching
  MatchJets1to1(jetsReco, jetsTruth);

  auto findBin = [&](const std::vector<double>& edges, double x)->int {
    if (x < edges.front() || x >= edges.back()) return -1;
    auto it = std::upper_bound(edges.begin(), edges.end(), x);
    int idx = int(it - edges.begin()) - 1;
    return (idx>=0 && idx<int(edges.size())-1) ? idx : -1;
  };
  auto flatIndex = [&](int ptbin, int zbin)->int { return zbin + m_nZsjBins * ptbin; };

  // Fill matched
  for (auto& kv: recoToTruth) {
    Jet* r = kv.first;
    Jet* t = kv.second;

    const float dr = DeltaR(r,t);
    
    if (!respPt || !hRecoPt || !hTruthPt) {
      std::cerr << "RooUnfold objects not initialized!\n";
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    //    hRecoPt->Fill(r->get_pt());
    // hTruthPt->Fill(t->get_pt());
    // respPt->Fill(r->get_pt(), t->get_pt());

    const double xr = r->get_pt();
    const double yt = t->get_pt();
    const double w = 1.0;
    
    const int bx = hRecoPt->FindFixBin(xr);
    const int by = hTruthPt->FindFixBin(yt);
    const bool inx = (bx >= 1 && bx <= hRecoPt->GetNbinsX());
    const bool iny = (by >= 1 && by <= hTruthPt->GetNbinsX());

    if (inx && iny) {
      hRecoPt->Fill(xr);
      hTruthPt->Fill(yt);
      respPt->Fill(xr, yt, w);
    } else {
      //print of count how many are out-of-range
      std::cerr << "Out-of-range Fill: xr="<<xr<<" (bin "<<bx<<"/"<<hRecoPt->GetNbinsX()
		<< "), yt="<<yt<<" (bin "<<by<<"/"<<hTruthPt->GetNbinsX()<<")\n";
    }
    v_match_reco_pt.push_back(r->get_pt());
    v_match_reco_eta.push_back(r->get_eta());
    v_match_reco_phi.push_back(r->get_phi());
    v_match_truth_pt.push_back(t->get_pt());
    v_match_truth_eta.push_back(t->get_eta());
    v_match_truth_phi.push_back(t->get_phi());
    v_match_dR.push_back(dr);

    double zr = zReco.count(r) ? zReco[r] : NAN;
    double zt = zTruth.count(t) ? zTruth[t] : NAN;
    v_match_reco_zsj.push_back(zr);
    v_match_truth_zsj.push_back(zt);

    int r_ptbin = findBin(m_ptEdges, r->get_pt());
    int t_ptbin = findBin(m_ptEdges, t->get_pt());
    int r_zbin  = (std::isfinite(zr) ? findBin(m_zsjEdges, zr) : -1);
    int t_zbin  = (std::isfinite(zt) ? findBin(m_zsjEdges, zt) : -1);
    if (r_ptbin>=0 && t_ptbin>=0 && r_zbin>=0 && t_zbin>=0) {
      int r_flat = flatIndex(r_ptbin, r_zbin);
      int t_flat = flatIndex(t_ptbin, t_zbin);
      respPtZ->Fill(hRecoPtZ->GetBinCenter(r_flat+1), hTruthPtZ->GetBinCenter(t_flat+1));
      hRecoPtZ->Fill(hRecoPtZ->GetBinCenter(r_flat+1));
      hTruthPtZ->Fill(hTruthPtZ->GetBinCenter(t_flat+1));
    }
  }

  // Fakes (reco-only)
  for (auto r: *jetsReco) {
    if (r->get_pt() < m_recoPtMin) continue;
    if (std::abs(r->get_eta()) > m_etaRange.second) continue;
    if (recoToTruth.find(r) != recoToTruth.end()) continue;

    // Fakes (reco-only)
    const double xr = r->get_pt();
    const double w = 1.0; //arbitrary weight for now
    const int bx = hRecoPt->FindFixBin(xr);
    if (bx >= 1 && bx <= hRecoPt->GetNbinsX()) {
      hRecoPt->Fill(xr);
      respPt->Fake(xr, w);
    } else {
      std::cout << __LINE__ << std::endl; // optional debug log
    }
    
    double zr = zReco.count(r) ? zReco[r] : NAN;
    v_fake_reco_pt.push_back(r->get_pt());
    v_fake_reco_eta.push_back(r->get_eta());
    v_fake_reco_phi.push_back(r->get_phi());
    v_fake_reco_zsj.push_back(zr);

    int r_ptbin = findBin(m_ptEdges, r->get_pt());
    int r_zbin  = (std::isfinite(zr) ? findBin(m_zsjEdges, zr) : -1);
    if (r_ptbin>=0 && r_zbin>=0) {
      int r_flat = flatIndex(r_ptbin, r_zbin);
      respPtZ->Fake(hRecoPtZ->GetBinCenter(r_flat+1));
      hRecoPtZ->Fill(hRecoPtZ->GetBinCenter(r_flat+1));
    }
  }

  // Misses (truth-only)
  for (auto t: *jetsTruth) {
    if (t->get_pt() < m_truthPtMin) continue;
    if (std::abs(t->get_eta()) > m_etaRange.second) continue;
    if (truthToReco.find(t) != truthToReco.end()) continue;

    hTruthPt->Fill(t->get_pt());
    respPt->Miss(t->get_pt());

    double zt = zTruth.count(t) ? zTruth[t] : NAN;
    v_fake_truth_pt.push_back(t->get_pt());
    v_fake_truth_eta.push_back(t->get_eta());
    v_fake_truth_phi.push_back(t->get_phi());
    v_fake_truth_zsj.push_back(zt);

    int t_ptbin = findBin(m_ptEdges, t->get_pt());
    int t_zbin  = (std::isfinite(zt) ? findBin(m_zsjEdges, zt) : -1);
    if (t_ptbin>=0 && t_zbin>=0) {
      int t_flat = flatIndex(t_ptbin, t_zbin);
      respPtZ->Miss(hTruthPtZ->GetBinCenter(t_flat+1));
      hTruthPtZ->Fill(hTruthPtZ->GetBinCenter(t_flat+1));
    }
  }

  m_T->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}
int JetMatchingSubjets::End(PHCompositeNode*)
{
  PHTFileServer::get().cd(m_outputFileName);

  if (m_T) { m_T->Write(); m_T = nullptr; }
  //std::cout << __LINE__ << std::endl;
  
  if (respPt) {
    TH2* h2 = respPt->Hresponse();// uses template axes
    //std::cout << __LINE__ << std::endl;
    if (h2) {
      TH2* h2clone = static_cast<TH2*>(h2->Clone("hResponsePt"));
      //std::cout << __LINE__ << std::endl;
      h2clone->SetDirectory(nullptr);
      //std::cout << __LINE__ << std::endl;
      h2clone->Write();
      //std::cout << __LINE__ << std::endl;
      //delete h2clone;
      //std::cout << __LINE__ << std::endl;
    }
    //    delete respPt;
    //std::cout << __LINE__ << std::endl;
    if (!respPt) {
      //std::cout << __LINE__ << std::endl;
    }
  }

  // Write & delete the template spectra AFTER respPt is gone
  if (hRecoPt)   { hRecoPt->Write();   delete hRecoPt;   hRecoPt   = nullptr; }
  //std::cout << __LINE__ << std::endl;
  if (hTruthPt)  { hTruthPt->Write();  delete hTruthPt;  hTruthPt  = nullptr; }
  //std::cout << __LINE__ << std::endl;
  if (hRespTempl) { delete hRespTempl; hRespTempl=nullptr; }
  //std::cout << __LINE__ << std::endl;
  
  if (respPtZ) {
    TH2* h2 = respPtZ->Hresponse();
    //std::cout << __LINE__ << std::endl;
    if (h2) {
      TH2* h2clone = static_cast<TH2*>(h2->Clone("hResponsePtZ_flat"));
      //std::cout << __LINE__ << std::endl;
      h2clone->SetDirectory(nullptr);
      //std::cout << __LINE__ << std::endl;
      h2clone->Write();
      //std::cout << __LINE__ << std::endl;
      delete h2clone;
      //std::cout << __LINE__ << std::endl;
    }
    //    delete respPtZ; respPtZ = nullptr;
    //std::cout << __LINE__ << std::endl;
  }

  if (hRecoPtZ)  { hRecoPtZ->Write();  delete hRecoPtZ;  hRecoPtZ  = nullptr; }
  //std::cout << __LINE__ << std::endl;
  if (hTruthPtZ) { hTruthPtZ->Write(); delete hTruthPtZ; hTruthPtZ = nullptr; }
  //std::cout << __LINE__ << std::endl;
  if (hZsjEdges) { hZsjEdges->Write(); delete hZsjEdges; hZsjEdges = nullptr; }
  //std::cout << __LINE__ << std::endl;
  if (hPtEdges)  { hPtEdges->Write();  delete hPtEdges;  hPtEdges  = nullptr; }
  //std::cout << __LINE__ << std::endl;
  
  //  PHTFileServer::get().close();
  //std::cout << __LINE__ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingSubjets::Reset(PHCompositeNode*) { return Fun4AllReturnCodes::EVENT_OK; }
  
void JetMatchingSubjets::Print(const std::string& what) const { }
  
// ---- helper definitions ----
static inline float wrapPhi(float dphi){
  //std::cout << __LINE__ << std::endl;
  while (dphi >  M_PI) dphi -= 2*M_PI;
  //std::cout << __LINE__ << std::endl;
  while (dphi < -M_PI) dphi += 2*M_PI;
  //std::cout << __LINE__ << std::endl;
  return dphi; 
}

float JetMatchingSubjets::DeltaR(Jet* a, Jet* b)
{
  //std::cout << __LINE__ << std::endl;
  const float dEta = a->get_eta() - b->get_eta();
  //std::cout << __LINE__ << std::endl;
  float dPhi = wrapPhi(a->get_phi() - b->get_phi());
  //std::cout << __LINE__ << std::endl;
  return std::sqrt(dEta*dEta + dPhi*dPhi);
}

void JetMatchingSubjets::MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets)
{
  //std::cout << __LINE__ << std::endl;
  struct Pair { float dr; Jet* r; Jet* t; };
  //std::cout << __LINE__ << std::endl;
  std::vector<Pair> cand;
  //std::cout << __LINE__ << std::endl;
  for (auto r: *recoJets) {
    //std::cout << __LINE__ << std::endl;
    if (r->get_pt() < m_recoPtMin) continue;
    //std::cout << __LINE__ << std::endl;
    if (std::abs(r->get_eta()) > m_etaRange.second) continue;
    //std::cout << __LINE__ << std::endl;
    for (auto t: *truthJets) {
      //std::cout << __LINE__ << std::endl;
      if (t->get_pt() < m_truthPtMin) continue;
      //std::cout << __LINE__ << std::endl;
      if (std::abs(t->get_eta()) > m_etaRange.second) continue;
      //std::cout << __LINE__ << std::endl;
      float dr = DeltaR(r,t);
      //std::cout << __LINE__ << std::endl;
      if (dr < m_matchDRMax) cand.push_back({dr, r, t});
      //std::cout << __LINE__ << std::endl;
    }
    //std::cout << __LINE__ << std::endl;
  }
  std::sort(cand.begin(), cand.end(), [](const Pair& a, const Pair& b){ return a.dr < b.dr; });
  //std::cout << __LINE__ << std::endl;
  std::set<Jet*> usedR, usedT;
  //std::cout << __LINE__ << std::endl;
  recoToTruth.clear();
  //std::cout << __LINE__ << std::endl;
  truthToReco.clear();
  //std::cout << __LINE__ << std::endl;
  for (auto& p: cand) {
    //std::cout << __LINE__ << std::endl;
    if (usedR.count(p.r) || usedT.count(p.t)) continue;
    //std::cout << __LINE__ << std::endl;
    usedR.insert(p.r);
    //std::cout << __LINE__ << std::endl;
    usedT.insert(p.t);
    //std::cout << __LINE__ << std::endl;
    recoToTruth[p.r] = p.t;
    //std::cout << __LINE__ << std::endl;
    truthToReco[p.t] = p.r;
    //std::cout << __LINE__ << std::endl;
  }
}

