#ifndef JETUNFOLDINGSUBJETS_H
#define JETUNFOLDINGSUBJETS_H

#include <fun4all/SubsysReco.h>
#include <vector>
#include <map>
#include <string>
#include <fastjet/PseudoJet.hh>

// Forward declarations for ROOT classes
class TTree;
class TH1F;
class RooUnfoldResponse;
class JetContainer;
class Jet;
class TowerInfoContainer;
class RawTowerGeomContainer;
class TowerBackground;
class PHG4TruthInfoContainer;
class PHCompositeNode;

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
    void Print(const std::string& what = "ALL") const override;
protected:
    // Matching and subjet analysis
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
        TowerBackground* bg, float v2, float psi2,
        PHG4TruthInfoContainer* truthInfo);
    void AnalyzeTruthJets(JetContainer* truthJets, PHG4TruthInfoContainer* truthInfo);

    // I/O and configuration
    std::string m_recoJetName;
    std::string m_truthJetName;
    std::string m_outputFileName;

    // Output tree and event-level info
    TTree* m_T = nullptr;
    int m_event = -1;
    float m_centrality = -1.0f;
    float m_impactparam = -1.0f;

    // Jet tree variables
    std::vector<float> m_pt;
    std::vector<float> m_eta;
    std::vector<float> m_phi;
    std::vector<float> m_pt_truth;
    std::vector<float> m_eta_truth;
    std::vector<float> m_phi_truth;

    // Jet matching maps
    std::map<Jet*, Jet*> recoToTruth;
    std::map<Jet*, Jet*> truthToReco;

    // 1D response and spectra
    RooUnfoldResponse* m_response1D = nullptr;
    TH1F* hRecoJetPtMatched = nullptr;
    TH1F* hTruthJetPtMatched = nullptr;
    TH1F* hRecoJetPtUnfolded = nullptr;

    // z_sj per-pt-bin histograms and responses
    std::vector<TH1F*> m_hRecoZsjMatched;
    std::vector<TH1F*> m_hTruthZsjMatched;
    std::vector<TH1F*> m_hRecoZsjUnfolded;
    std::vector<RooUnfoldResponse*> m_responseZsj;
};

#endif // JETUNFOLDINGSUBJETS_H
