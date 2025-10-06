#include "JetMatchingForUnfolding.h"
#include "SubjetTools.h"
#include "MatchingUtils.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <jetbase/JetContainer.h>
#include <jetbase/JetMap.h>
#include <jetbase/Jet.h>

#include <centrality/CentralityInfo.h>
#include <jetbackground/TowerBackground.h>

#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>

#include <TDirectory.h>
#include <TH1.h>
#include <TH1D.h>
#include <TTree.h>
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

JetMatchingForUnfolding::JetMatchingForUnfolding(const std::string& recojetname,
                                       const std::string& truthjetname,
                                       const std::string& outputfilename)
  : SubsysReco(std::string("JetMatchingForUnfolding_")+recojetname+"_"+truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
{}

JetMatchingForUnfolding::~JetMatchingForUnfolding() {}

int JetMatchingForUnfolding::Init(PHCompositeNode*)
{
  
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  
  if (m_sjPtCuts.empty())
    {
      m_sjPtCuts = {0.0, 3.0, 6.0, 10.0};
    std::cout << "JetMatchingForUnfolding::Init - m_sjPtCuts was empty; using default {0,3,6,10} GeV\n";
    }
  
  // 3) book per-cut subjet spectra histograms
  h_sj1_reco.resize(m_sjPtCuts.size());
  h_sj2_reco.resize(m_sjPtCuts.size());
  h_sj1_truth.resize(m_sjPtCuts.size());
  h_sj2_truth.resize(m_sjPtCuts.size());
  
  for (size_t i = 0; i < m_sjPtCuts.size(); ++i)
    {
      const double cut = m_sjPtCuts[i];
      // Use integer-like label for name to keep it filesystem/ROOT friendly
      const auto tag = Form("ptMin%g", cut);
      
      h_sj1_reco[i]  = new TH1F(Form("h_sj1_reco_%s",  tag),
				Form("Reco leading subjet p_{T} (p_{T}^{subjet} > %.1f GeV);p_{T} [GeV];Counts", cut),
				11, 0, 60);
      h_sj2_reco[i]  = new TH1F(Form("h_sj2_reco_%s",  tag),
				Form("Reco subleading subjet p_{T} (p_{T}^{subjet} > %.1f GeV);p_{T} [GeV];Counts", cut),
				11, 0, 60);
      h_sj1_truth[i] = new TH1F(Form("h_sj1_truth_%s", tag),
				Form("Truth leading subjet p_{T} (p_{T}^{subjet} > %.1f GeV);p_{T} [GeV];Counts", cut),
				11, 0, 60);
      h_sj2_truth[i] = new TH1F(Form("h_sj2_truth_%s", tag),
				Form("Truth subleading subjet p_{T} (p_{T}^{subjet} > %.1f GeV);p_{T} [GeV];Counts", cut),
				11, 0, 60);
      
      h_sj1_reco[i]->Sumw2();  h_sj2_reco[i]->Sumw2();
      h_sj1_truth[i]->Sumw2(); h_sj2_truth[i]->Sumw2();
    }
  
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

int JetMatchingForUnfolding::process_event(PHCompositeNode* topNode)
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
  //JetContainer* jetsReco = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  auto* jetsTruth = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  //JetMap* jetsTruth = findNode::getClass<JetMap>(topNode, m_truthJetName);
  auto* cent      = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");

  m_centrality = (int)(100 * cent->get_centile(CentralityInfo::PROP::mbd_NS));
  m_b          = cent->get_quantity(CentralityInfo::PROP::bimp);

  // Towers (use SUB1 to mirror JetUnfolding.cc)
  auto* em =  findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");  //previously TOWERINFO_CALIB_CEMC_RETOWER
  auto* ih =  findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
  auto* oh =  findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
  
  // Geometry (unchanged)
  auto* geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto* geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  auto* geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  
  // Background
  auto* bg = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");
  auto* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  
  if (!jetsReco || !jetsTruth || !cent || !em || !ih || !oh || !geomEM || !geomIH || !geomOH || !bg || !truthinfo) {
    std::cerr << "[JetMatchingForUnfolding] Missing node(s)." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto wrapPhi = [](double a){
    while (a >  M_PI) a -= 2*M_PI;
    while (a < -M_PI) a += 2*M_PI;
    return a;
  };
  
  // Reco jets
  for (auto j: *jetsReco) {
    if (j->get_pt() < m_recoPtMin) continue; 
    if (std::abs(j->get_eta()) > m_etaRange.second) continue;
    
    v_reco_pt.push_back(j->get_pt());
    v_reco_eta.push_back(j->get_eta());
    v_reco_phi.push_back(wrapPhi(j->get_phi()));
    double z = std::numeric_limits<double>::quiet_NaN();
    double th = std::numeric_limits<double>::quiet_NaN();

    const bool ok = ComputeSubjetForJet(j,
			   em, ih, oh, geomEM, geomIH, geomOH, bg,
			   /*v2=*/0.f, /*psi2=*/0.f, /*doUnsub=*/false,
			   /*etaCalMax=*/0.7, /*R_jet=*/0.4,
			   /*z_cut=*/0.0, /*beta=*/0.0, /*pt_min_subjet=*/0.0,
			   z, th);
    
    v_reco_zsj.push_back(ok ? z  : std::numeric_limits<double>::quiet_NaN());
    v_reco_theta.push_back(ok ? th : std::numeric_limits<double>::quiet_NaN());
    
  }

  // Truth jets
  for (auto* j : *jetsTruth) {
    if (j->get_pt() < m_truthPtMin) continue;
    if (std::abs(j->get_eta()) > m_etaRange.second) continue;
    
    v_truth_pt.push_back(j->get_pt());
    v_truth_eta.push_back(j->get_eta());
    v_truth_phi.push_back(wrapPhi(j->get_phi()));
    
    double z  = std::numeric_limits<double>::quiet_NaN();
    double th = std::numeric_limits<double>::quiet_NaN();
    
    bool ok = false;
    if (m_useTruthFromSimTowers) {
      ok = ComputeSubjetForTruthJet_FromTowersCone(
						   j, em, ih, oh,
						   geomEM, geomIH, geomOH,
						   /*R_jet=*/0.4, /*etaCalMax=*/m_etaRange.second,
						   /*pt_min_subjet=*/0.0,
						   z, th);
    } else {
      ok = ComputeSubjetForTruthJet(
				    j, truthinfo,
				    /*R_jet=*/0.4, /*etaCalMax=*/m_etaRange.second,
				    /*pt_min_subjet=*/0.0,
				    z, th);
    }
    
    v_truth_zsj.push_back(ok ? z  : std::numeric_limits<double>::quiet_NaN());
    v_truth_theta.push_back(ok ? th : std::numeric_limits<double>::quiet_NaN());
  }

  // Matching
  // --- JetUnfolding-identical matching ---
  const float DR_MAX       = 0.1f;
  const float RECO_PT_MIN  = 5.0f;
  const float ETA_MAX_ABS  = 0.7f;
  const float TRUTH_PT_MIN = 52.0f;
  
  MatchResult mr = MatchJetsLikeUnfolding(
					  jetsReco, jetsTruth, DR_MAX, RECO_PT_MIN, ETA_MAX_ABS, TRUTH_PT_MIN);
  
  recoToTruth.clear();
  truthToReco.clear();
  for (size_t i = 0; i < mr.matchedRecoJets.size(); ++i) {
    Jet* r = const_cast<Jet*>(mr.matchedRecoJets[i]);
    Jet* t = const_cast<Jet*>(mr.matchedTruthJets[i]);
    recoToTruth[r] = t;
    truthToReco[t] = r;
  }

  // ---- Fill matched (only keep zsj when BOTH reco & truth have valid subjets) ----
  for (auto& kv: recoToTruth) {
    Jet* r = kv.first;
    Jet* t = kv.second;
    
    const float dr = calculateDeltaR(r->get_eta(), r->get_phi(),
				     t->get_eta(), t->get_phi());
    double zr = std::numeric_limits<double>::quiet_NaN();
    double tr = std::numeric_limits<double>::quiet_NaN();
    const bool okr = ComputeSubjetForJet(
				      r, em, ih, oh, geomEM, geomIH, geomOH, bg,
				      /*v2=*/0.f, /*psi2=*/0.f, /*doUnsub=*/false,
				      /*etaCalMax=*/0.7, /*R_jet=*/0.4,
				      /*z_cut=*/0.0, /*beta=*/0.0, /*pt_min_subjet=*/0.0,
				      zr, tr);
    // truth z, theta from the matched truth jet
    double zt = std::numeric_limits<double>::quiet_NaN();
    double tt = std::numeric_limits<double>::quiet_NaN();

    bool okt = false;
    if (m_useTruthFromSimTowers) {
      okt = ComputeSubjetForTruthJet_FromTowersCone(
						    t, em, ih, oh,
						    geomEM, geomIH, geomOH,
						    /*R_jet=*/0.4, /*etaCalMax=*/0.7, /*pt_min_subjet=*/0.0,
						    zt, tt);
    } else {
      okt = ComputeSubjetForTruthJet(
				     t, truthinfo,
				     /*R_jet=*/0.4, /*etaCalMax=*/0.7,
				     /*pt_min_subjet=*/0.0,
				     zt, tt);
    }
    
    // Only store as "matched" if BOTH are valid
    if (okr && okt) {
      v_match_reco_pt.push_back(r->get_pt());
      v_match_reco_eta.push_back(r->get_eta());
      v_match_reco_phi.push_back(wrapPhi(r->get_phi()));
      v_match_reco_zsj.push_back(zr);
      v_match_reco_theta.push_back(tr);
      
      v_match_truth_pt.push_back(t->get_pt());
      v_match_truth_eta.push_back(t->get_eta());
      v_match_truth_phi.push_back(wrapPhi(t->get_phi()));
      v_match_truth_zsj.push_back(zt);
      v_match_truth_theta.push_back(tt);
      
      v_match_dR.push_back(dr);

      for (size_t ic = 0; ic < m_sjPtCuts.size(); ++ic) {
	const double cut = m_sjPtCuts[ic];
	
	double r_pt1, r_pt2, t_pt1, t_pt2;
	const bool okRpts = LeadingSubjetPts_RecoJet(
						     r, em, ih, oh, geomEM, geomIH, geomOH, bg,
						     /*v2=*/0.f, /*psi2=*/0.f, /*doUnsub=*/false,
						     /*etaCalMax=*/0.7, /*R_jet=*/0.4,
						     /*ptMinSubjet=*/cut, r_pt1, r_pt2);
	
	bool okTpts = false;
	if (m_useTruthFromSimTowers) {
	  okTpts = LeadingSubjetPts_TruthJet_FromTowersCone(
							    t, em, ih, oh,
							    geomEM, geomIH, geomOH,
							    /*R_jet=*/0.4, /*etaCalMax=*/0.7,
							    /*ptMinSubjet=*/cut, t_pt1, t_pt2);
	} else {
	  okTpts = LeadingSubjetPts_TruthJet(
					     t, truthinfo, /*R_jet=*/0.4, /*etaCalMax=*/0.7,
					     /*ptMinSubjet=*/cut, t_pt1, t_pt2);
	}
	
	if (okRpts) { h_sj1_reco[ic]->Fill(r_pt1); h_sj2_reco[ic]->Fill(r_pt2); }
	if (okTpts) { h_sj1_truth[ic]->Fill(t_pt1); h_sj2_truth[ic]->Fill(t_pt2); }
      }
    } else {
      // Route partial/failed subjet cases to fakes
      if (okr) {
	v_fake_reco_pt.push_back(r->get_pt());
	v_fake_reco_eta.push_back(r->get_eta());
	v_fake_reco_phi.push_back(wrapPhi(r->get_phi()));
	v_fake_reco_zsj.push_back(zr);
	v_fake_reco_theta.push_back(tr);
      }
      if (okt) {
	v_fake_truth_pt.push_back(t->get_pt());
	v_fake_truth_eta.push_back(t->get_eta());
	v_fake_truth_phi.push_back(wrapPhi(t->get_phi()));
	v_fake_truth_zsj.push_back(zt);
	v_fake_truth_theta.push_back(tt);
      }
      // If neither side had valid subjets (case 3), nothing goes to "matched".
      // Still get a record of these via the fake_* entries above only if one side succeeded.
      // If neither succeeded here, they simply won't appear in match_* nor fake_* (both invalid).
    }
  }
  // Fakes (reco-only, outside matched loop)
  for (size_t ri = 0; ri < jetsReco->size(); ++ri) {
    Jet* r = (*jetsReco)[ri];
    if (r->get_pt() <= RECO_PT_MIN) continue;
    if (std::abs(r->get_eta()) > ETA_MAX_ABS) continue;
    if (ri < mr.recoMatched.size() && mr.recoMatched[ri]) continue;

    v_fake_reco_pt.push_back(r->get_pt());
    v_fake_reco_eta.push_back(r->get_eta());
    v_fake_reco_phi.push_back(wrapPhi(r->get_phi()));
    
    double zr = std::numeric_limits<double>::quiet_NaN();
    double tr = std::numeric_limits<double>::quiet_NaN();
    (void)ComputeSubjetForJet(r,
			   em, ih, oh, geomEM, geomIH, geomOH, bg,
			   /*v2=*/0.f, /*psi2=*/0.f, /*doUnsub=*/false,
			   /*etaCalMax=*/0.7, /*R_jet=*/0.4,
			   /*z_cut=*/0.0, /*beta=*/0.0, /*pt_min_subjet=*/0.0,
			   zr, tr);
    
    v_fake_reco_zsj.push_back(zr);
    v_fake_reco_theta.push_back(tr);
  }
  
  // Misses (truth-only)
  for (size_t ti = 0; ti < jetsTruth->size(); ++ti) {
    Jet* t = (*jetsTruth)[ti];
    if (t->get_pt() < TRUTH_PT_MIN) continue;
    if (std::abs(t->get_eta()) > ETA_MAX_ABS) continue;
    if (ti < mr.truthMatched.size() && mr.truthMatched[ti]) continue;
        
    double zt, tt;
    if (!ComputeSubjetForTruthJet(t, truthinfo,
			  /*R_jet=*/0.4,
			  /*etaCalMax=*/0.7,
			       0.0, zt, tt)) {
      zt=NAN; tt=NAN;
    }
    v_fake_truth_pt.push_back(t->get_pt());
    v_fake_truth_eta.push_back(t->get_eta());
    v_fake_truth_phi.push_back(wrapPhi(t->get_phi()));
    v_fake_truth_zsj.push_back(zt);
    v_fake_truth_theta.push_back(tt);
  }
  
  m_T->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1.h>

// ...

int JetMatchingForUnfolding::End(PHCompositeNode* /*topNode*/)
{
  std::cout << "JetMatchingForUnfolding::End - Writing output to "
            << m_outputFileName << std::endl;

  // Go to the output file managed by PHTFileServer
  PHTFileServer::get().cd(m_outputFileName);

  // 1) Write the TTree (if it exists)
  if (m_T)
  {
    m_T->Write();
  }

  // 2) Write histograms (group in a directory for tidiness)
  TDirectory* saveDir = gDirectory;
  TDirectory* dTop = gDirectory->GetDirectory("subjet_spectra");
  if (!dTop) dTop = gDirectory->mkdir("subjet_spectra");
  if (dTop) dTop->cd();

  // Make per-cut subdirs and write the 4 spectra for each cut
  for (size_t i = 0; i < m_sjPtCuts.size(); ++i)
  {
    const double cut = m_sjPtCuts[i];

    // subdir name safe for ROOT
    const TString subname = Form("ptMin_%.0f", cut);
    TDirectory* dCut = gDirectory->GetDirectory(subname);
    if (!dCut) dCut = gDirectory->mkdir(subname);
    if (dCut) dCut->cd();

    if (i < h_sj1_reco.size()  && h_sj1_reco[i])  h_sj1_reco[i]->Write();
    if (i < h_sj2_reco.size()  && h_sj2_reco[i])  h_sj2_reco[i]->Write();
    if (i < h_sj1_truth.size() && h_sj1_truth[i]) h_sj1_truth[i]->Write();
    if (i < h_sj2_truth.size() && h_sj2_truth[i]) h_sj2_truth[i]->Write();

    if (dTop) dTop->cd(); // step back up for the next subdir
  }

  if (saveDir) saveDir->cd();

  TParameter<double>    pS("sigma_pb",  m_sigma_pb);
  TNamed                pName("sample_name", m_sample_name.c_str());
  TParameter<long long>("N_in_zwin", m_N_in_zwin).Write(); // number of entries written (|z|<=30)                                                                                                
  pS.Write();
  pName.Write();
  
  std::cout << "JetMatchingForUnfolding::End - Done." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingForUnfolding::Reset(PHCompositeNode*) { return Fun4AllReturnCodes::EVENT_OK; }
void JetMatchingForUnfolding::Print(const std::string& what) const { }
  
