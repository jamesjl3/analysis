#pragma once

#include <fun4all/SubsysReco.h>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <fastjet/PseudoJet.hh>

class PHCompositeNode;
class JetContainer;
class Jet;
class TH1F;
class TTree;
class RooUnfoldResponse;
class TowerInfoContainer;
class RawTowerGeomContainer;
class TowerBackground;
class PHG4TruthInfoContainer;

class JetUnfoldingSubjets : public SubsysReco {
public:
    JetUnfoldingSubjets(const std::string& recojetname,
                        const std::string& truthjetname,
                        const std::string& outputfilename);
    ~JetUnfoldingSubjets() override;

    int Init(PHCompositeNode*) override;
    int process_event(PHCompositeNode*) override;
    int End(PHCompositeNode*) override;
    int Reset(PHCompositeNode*) override;
    void Print(const std::string& what="") const override;

private:
    using JetMap = std::map<Jet*, Jet*>;
    using JetMapPtr = std::unique_ptr<JetMap>;
  
    void MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets, float dRMax);
    
    std::vector<fastjet::PseudoJet> BuildPseudoJets(
        Jet* jet, TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
        RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
        TowerBackground* bg, float v2, float psi2, bool doUnsub) const;
    std::vector<fastjet::PseudoJet> BuildTruthPseudoJets(
        Jet* truthJet, PHG4TruthInfoContainer* truthInfo) const;
    void AnalyzeMatchedJets(JetContainer* recoJets,
        TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
        RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
        TowerBackground* bg, float v2, float psi2, PHG4TruthInfoContainer* truthInfo);
    void AnalyzeTruthJets(JetContainer* truthJets, PHG4TruthInfoContainer* truthInfo);

    std::string m_recoJetName;
    std::string m_truthJetName;
    std::string m_outputFileName;

    std::unique_ptr<TTree> m_T;
    int m_event;
    float m_centrality;
    float m_impactparam;
    std::vector<float> m_pt, m_eta, m_phi, m_pt_truth, m_eta_truth, m_phi_truth;

    std::unique_ptr<RooUnfoldResponse> m_response1D;
    std::unique_ptr<TH1F> hRecoJetPtMatched;    // Add this line
    std::unique_ptr<TH1F> hTruthJetPtMatched;   // Add this line
    std::unique_ptr<TH1F> hRecoJetPtUnfolded;   // Add this line

  JetMap recoToTruth; // Change to JetMap instead of unique_ptr<JetMap>
  JetMap truthToReco; // Change to JetMap instead of unique_ptr<JetMap>
  
  std::vector<std::unique_ptr<RooUnfoldResponse>> m_responseZsj;
  std::vector<std::unique_ptr<TH1F>> m_hRecoZsjMatched;
  std::vector<std::unique_ptr<TH1F>> m_hTruthZsjMatched;
  std::vector<std::unique_ptr<TH1F>> m_hRecoZsjUnfolded;
};

