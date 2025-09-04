#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>
#include <g4centrality/PHG4CentralityReco.h>
#include "../HIJetReco.C"
#include "../../src/JetMatchingSubjets/JetMatchingSubjets.h"
#include <jetbackground/TowerRhov1.h>
#include <jetbackground/DetermineTowerRho.h>
#include <towerrhoana/TowerRhoAna.h>
#include <caloembedding/caloTowerEmbed.h>
#include <jetbkgdsub/JetBkgdSub.h>
#include <TSystem.h>
#include <phool/recoConsts.h>
#include <jetbase/JetReco.h>
#include <G4_Global.C>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetCalib.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libJetValidation.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libTowerRhoAna.so)
R__LOAD_LIBRARY(libCaloEmbedding.so)
gSystem->Load("/sphenix/user/jamesj3j3/sPHENIX/install/lib/libcalo_reco.so");
R__LOAD_LIBRARY(libJetBkgdSub.so)

std::cout << "ROOT Dynamic Path: " << gSystem->GetDynamicPath() << std::endl;

#endif

//gSystem->mkdir("/sphenix/tg/tg01/jets/jamesj3j3/JetMatchingSubjets/", true);

void Fun4All_JetMatchingSubjets(const char *filelisttruth = "dst_truth_jet.list",
				const char *filelistcalo = "dst_calo_cluster.list",
				const char *filelistglobal = "dst_global.list",
				const char *filelisttruthg4hit = "g4hits.list",
				const char *outname = "outputest_JetMatchingSubjets.root",
				int start_event = 0, int end_event = 1000)
{

  std::string jetNode = "AntiKt_Tower_r04_Sub1";
  Fun4AllServer *se = Fun4AllServer::instance();
  int verbosity = 0;
  
  se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();

  PHG4CentralityReco *cent = new PHG4CentralityReco();
  cent->Verbosity(verbosity);
  cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
  se->registerSubsystem( cent );
  Enable::VERBOSITY = verbosity;
  
  HIJetReco();
  
  recoConsts::instance()->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  
  JetCalib* jetcalib = new JetCalib("JetCalib");
  jetcalib->set_InputNode("AntiKt_Tower_r04_Sub1");        // Uncalibrated jets                                                                                                                       
  jetcalib->set_OutputNode("AntiKt_Tower_r04_Sub1_Calib"); // Calibrated jets                                                                                                                         
  jetcalib->set_JetRadius(0.4);
  jetcalib->set_ZvrtxNode("GlobalVertexMap");
  jetcalib->set_ApplyEtaDependentCalib(true);
  jetcalib->set_ApplyZvrtxDependentCalib(true);
  jetcalib->Verbosity(0);

  

  //  myJetUnf->doUnsub(1);
  //  myJetUnf->doSeeds(0);
  //   myJetUnf->set_doSim(false);
  //  myJetUnf->doClusters(0);
  JetMatchingSubjets *myJetUnf = new JetMatchingSubjets("AntiKt_Tower_r04_Sub1_Calib", "AntiKt_Truth_r04", outname);  
  myJetUnf->setRecoConstituentJetNode("AntiKt_Tower_r04_Sub1");
  myJetUnf->setZWindow(-30.f, 30.f);
  myJetUnf->setEtaRange(-0.7, 0.7);
  myJetUnf->setRecoPtMin(5.0);
  myJetUnf->setTruthPtMin(10.0);
  myJetUnf->setMatchDRMax(0.2);
  myJetUnf->setPtBinning({5,10,15,20,25,30,35,40,45,50,55});
  myJetUnf->setZsjBinning({0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50});

  se->registerSubsystem(jetcalib);
  se->registerSubsystem(myJetUnf);

  Fun4AllInputManager *intrue = new Fun4AllDstInputManager("DSTtruth");
   intrue->AddListFile(filelisttruth); //for standalone Fun4All running
  //  intrue->AddFile(filelisttruth); //for batch running to condor
  se->registerInputManager(intrue);

  Fun4AllInputManager *in4 = new Fun4AllDstInputManager("DSTg4hit");
  in4->AddListFile(filelisttruthg4hit); //for standalone Fun4All running
  //in4->AddFile(filelisttruthg4hit); //for batch running to condor
  se->registerInputManager(in4);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTcalo");
  in2->AddListFile(filelistcalo); //for standalone Fun4All running
  //in2->AddFile(filelistcalo); //for batch running to condor
  se->registerInputManager(in2);

  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("DSTglobal");
  in3->AddListFile(filelistglobal); //for standalone Fun4All running
  //in3->AddFile(filelistglobal); //for batch running to condor
  se->registerInputManager(in3);
  //  se->run(1000);
  se->skip(start_event); // Skip events before start_event
  se->run(end_event - start_event);
  se->End();
  //  gSystem->Exit(0);
  return;
}
