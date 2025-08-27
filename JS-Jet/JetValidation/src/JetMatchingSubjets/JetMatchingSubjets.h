#ifndef JETMATCHINGSUBJETS_H
#define JETMATCHINGSUBJETS_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <utility>
#include <set>
#include <map>
#include <functional>

class PHCompositeNode;
class TTree;
class TH1F;
class TH1D;
class TH2D;
class RooUnfoldResponse;
class Jet;
class JetContainer;
class CentralityInfo;

class JetMatchingSubjets : public SubsysReco
{
public:
  JetMatchingSubjets(const std::string& recojetname,
                     const std::string& truthjetname,
                     const std::string& outputfilename);
  ~JetMatchingSubjets() override;

  // Config
  void setEtaRange(double low, double high) { m_etaRange = {low, high}; }
  void setRecoPtMin(double x) { m_recoPtMin = x; }
  void setTruthPtMin(double x) { m_truthPtMin = x; }
  void setMatchDRMax(double x) { m_matchDRMax = x; }
  void setPtBinning(const std::vector<double>& edges) { m_ptEdges = edges; }
  void setZsjBinning(const std::vector<double>& edges) { m_zsjEdges = edges; }

  // Fun4All
  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;
  void Print(const std::string& what) const override;

private:
  // Helpers
  void MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets);
  static float DeltaR(Jet* a, Jet* b);

private:
  // Names / IO
  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  // Kinematic selections
  std::pair<double,double> m_etaRange { -0.7, 0.7 };
  double m_recoPtMin  = 5.0;
  double m_truthPtMin = 10.0;
  double m_matchDRMax = 0.3;

  // Binning
  std::vector<double> m_ptEdges;   // default set in Init()
  std::vector<double> m_zsjEdges;  // default set in Init()
  int m_nPtBins  = 0;
  int m_nZsjBins = 0;

  // Matching maps (1↔1 greedy)
  std::map<Jet*, Jet*> recoToTruth;
  std::map<Jet*, Jet*> truthToReco;

  // Tree
  TTree* m_T = nullptr;
  int   m_event = -1;
  float m_centrality = 0;
  float m_b = 0;

  // Per-event containers
  std::vector<float> v_reco_pt,  v_reco_eta,  v_reco_phi,  v_reco_zsj;
  std::vector<float> v_truth_pt, v_truth_eta, v_truth_phi, v_truth_zsj;

  // Matched pairs (parallel arrays, same length)
  std::vector<float> v_match_reco_pt,  v_match_reco_eta,  v_match_reco_phi,  v_match_reco_zsj;
  std::vector<float> v_match_truth_pt, v_match_truth_eta, v_match_truth_phi, v_match_truth_zsj;
  std::vector<float> v_match_dR;

  // Fakes / Misses
  std::vector<float> v_fake_reco_pt,  v_fake_reco_eta,  v_fake_reco_phi,  v_fake_reco_zsj;   // detector-only (no match)
  std::vector<float> v_fake_truth_pt, v_fake_truth_eta, v_fake_truth_phi, v_fake_truth_zsj;  // truth-only (missed)

  // JetMatchingSubjets.h (class members)
  TH1D* hRecoPt   = nullptr;
  TH1D* hTruthPt  = nullptr;
  RooUnfoldResponse* respPt = nullptr;
  
  TH1D* hRecoPtZ  = nullptr;
  TH1D* hTruthPtZ = nullptr;
  RooUnfoldResponse* respPtZ = nullptr;
  
  TH1D* hPtEdges  = nullptr;
  TH1D* hZsjEdges = nullptr;

  TH2D* hRespTempl = nullptr;
  
};

#endif
