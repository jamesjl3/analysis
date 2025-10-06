#ifndef JETMATCHINGFORUNFOLDING_H
#define JETMATCHINGFORUNFOLDING_H

#include <fun4all/SubsysReco.h>
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <TH1D.h>
#include <TH2.h>
#include <TTree.h>
#include <RooUnfoldResponse.h>

#include <string>
#include <vector>
#include <map>
#include <limits>

class PHCompositeNode;

class JetMatchingForUnfolding : public SubsysReco
{
public:
  JetMatchingForUnfolding(const std::string& recojetname,
                     const std::string& truthjetname,
                     const std::string& outputfilename);
  ~JetMatchingForUnfolding() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;
  int Reset(PHCompositeNode*) override;
  void Print(const std::string& what = "") const override;

  // Matching helpers
  void MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets);

  // Config
  void setPtEdges(const std::vector<double>& e){ m_ptEdges = e; }
  void setZsjEdges(const std::vector<double>& e){ m_zsjEdges = e; }
  void setRecoPtMin(float x){ m_recoPtMin = x; }
  void setTruthPtMin(float x){ m_truthPtMin = x; }
  void setEtaRange(float lo, float hi){ m_etaRange = {lo,hi}; }
  void setMatchDRMax(float x){ m_matchDRMax = x; }
  void setPtBinning(const std::vector<double>& edges) { m_ptEdges = edges; }
  void setZsjBinning(const std::vector<double>& edges) { m_zsjEdges = edges; }
  void setZWindow(float zmin, float zmax) { m_zmin = zmin; m_zmax = zmax; }
  void setRecoConstituentJetNode(const std::string& n) { m_recoConstituentJetNode = n; }
  void setSubjetPtScan(const std::vector<double>& cuts) { m_sjPtCuts = cuts; }
  void useTruthFromSimTowers(bool v) { m_useTruthFromSimTowers = v; }
  
private:
  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;
  std::string m_recoConstituentJetNode;
  
  // analysis cuts
  float m_recoPtMin = 5.f;
  float m_truthPtMin = 52.f;
  std::pair<float,float> m_etaRange = {-0.7f, 0.7f}; // align with |eta|<1.1-R=0.7 for R=0.4
  float m_matchDRMax = 0.1f;

  // Matching maps (1↔1 greedy)
  std::map<Jet*, Jet*> recoToTruth;
  std::map<Jet*, Jet*> truthToReco;

  // binning
  std::vector<double> m_ptEdges;
  std::vector<double> m_zsjEdges;
  int m_nPtBins = 0;
  int m_nZsjBins = 0;

  TTree* m_T = nullptr;
  int   m_event = -1;
  float m_centrality = 0;
  float m_b = 0;
  // --- z-vertex and metadata for per-file bookkeeping ---
  float      m_vtx_z    = std::numeric_limits<float>::quiet_NaN();
  long long  m_N_in_zwin = 0;
  
  // configurable z-window (used only wanting to count inside the module)
  float      m_zmin = -30.f;
  float      m_zmax =  30.f;
  
  // optional: store per-file sigma and sample name (purely for provenance)
  double     m_sigma_pb = 0.0;
  std::string m_sample_name = "unknown";
  
  void setSigmaPb(double s) { m_sigma_pb = s; }
  void setSampleName(const std::string& n) { m_sample_name = n; }
  
  // reco jets
  std::vector<float> v_reco_pt, v_reco_eta, v_reco_phi, v_reco_zsj, v_reco_theta;
  // truth jets
  std::vector<float> v_truth_pt, v_truth_eta, v_truth_phi, v_truth_zsj, v_truth_theta;

  // matched
  std::vector<float> v_match_reco_pt, v_match_reco_eta, v_match_reco_phi, v_match_reco_zsj, v_match_reco_theta;
  std::vector<float> v_match_truth_pt, v_match_truth_eta, v_match_truth_phi, v_match_truth_zsj, v_match_truth_theta;
  std::vector<float> v_match_dR;

  // fakes (reco-only) & misses (truth-only)
  std::vector<float> v_fake_reco_pt, v_fake_reco_eta, v_fake_reco_phi, v_fake_reco_zsj, v_fake_reco_theta;
  std::vector<float> v_fake_truth_pt, v_fake_truth_eta, v_fake_truth_phi, v_fake_truth_zsj, v_fake_truth_theta;
  std::vector<double> m_sjPtCuts;  // e.g. {0,3,6,10}
  std::vector<TH1F*> h_sj1_reco, h_sj2_reco, h_sj1_truth, h_sj2_truth;  

  bool m_useTruthFromSimTowers{false};
};

#endif
