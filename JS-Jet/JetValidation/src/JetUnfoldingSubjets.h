#ifndef JETUNFOLDINGSUBJETS_H
#define JETUNFOLDINGSUBJETS_H

#include <fun4all/SubsysReco.h>
#include <fastjet/PseudoJet.hh>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <TH2.h>

class Jet;
class JetContainer;
class TowerInfoContainer;
class RawTowerGeomContainer;
class TowerBackground;
class PHG4TruthInfoContainer;
class TH1F;
class TTree;
class RooUnfoldResponse;
class CentralityInfo;

class JetUnfoldingSubjets : public SubsysReco
{
public:
  JetUnfoldingSubjets(const std::string& recojetname = "AntiKt_Tower_r04",
                      const std::string& truthjetname = "AntiKt_Truth_r04",
                      const std::string& outputfilename = "JetUnfoldingSubjets.root");
  ~JetUnfoldingSubjets() override;

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;
  void Print(const std::string& what = "") const override;

private:
  void MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets, float dRMax);
  std::vector<fastjet::PseudoJet> BuildPseudoJets(Jet* jet, TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                                                  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
                                                  TowerBackground* bg, float v2, float psi2, bool doUnsub);
  std::vector<fastjet::PseudoJet> BuildTruthPseudoJets(Jet* truthJet, PHG4TruthInfoContainer* truthInfo);
  void AnalyzeMatchedJets(JetContainer* recoJets,
                          TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                          RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
                          TowerBackground* bg, float v2, float psi2,
                          PHG4TruthInfoContainer* truthInfo);
  void AnalyzeTruthJets(JetContainer* truthJets, PHG4TruthInfoContainer* truthInfo);

  // Input/output
  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  // Ranges
  std::pair<float, float> m_etaRange;
  std::pair<float, float> m_ptRange;
  int m_doTruthJets;
  
  // Tree and event info
  TTree* m_T;
  int m_event;
  float m_centrality;
  float m_impactparam;
  std::vector<float> m_pt, m_eta, m_phi;
  std::vector<float> m_pt_truth, m_eta_truth, m_phi_truth;

  // Matching
  std::map<Jet*, Jet*> recoToTruth;
  std::map<Jet*, Jet*> truthToReco;

  // 1D response and histograms
  //RooUnfoldResponse* m_response1D = nullptr;
  std::unique_ptr<RooUnfoldResponse> m_response1D;
  TH1F* hRecoJetPtMatched = nullptr;
  TH1F* hTruthJetPtMatched = nullptr;
  TH1F* hRecoJetPtUnfolded = nullptr;

  // z_sj histograms and response per pt bin
  std::vector<TH1F*> m_hRecoZsjMatched;
  std::vector<TH1F*> m_hTruthZsjMatched;
  //  std::vector<RooUnfoldResponse*> m_responseZsj;
  std::vector<std::unique_ptr<RooUnfoldResponse>> m_responseZsj;
  std::vector<TH1F*> m_hRecoZsjUnfolded;
};

#endif
