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
#include "../../src/JetMatchingSubjets.h"
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


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libJetValidation.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libTowerRhoAna.so)
R__LOAD_LIBRARY(libCaloEmbedding.so)
gSystem->Load("/sphenix/user/jamesj3j3/sPHENIX/install/lib/libcalo_reco.so");
//R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libJetBkgdSub.so)

//R__LOAD_LIBRARY(libtrackbase_historic.so)
//gSystem->Load("/sphenix/user/jamesj3j3/sPHENIX/install/lib/libtrackbase_historic.so")
//gSystem->ProcessLine(".L /sphenix/user/jamesj3j3/sPHENIX/install/lib/libJetValidation.so")
std::cout << "ROOT Dynamic Path: " << gSystem->GetDynamicPath() << std::endl;

#endif
/*
void Fun4All_JetMatchingSubjets(const char *filelisttruth = "dst_truth_jet.list",
			   const char *filelistcalo = "dst_calo_cluster.list",
			   const char *filelistglobal = "dst_global.list",
			   const char *filelisttruthg4hit = "dst_truth.list",
			   const char *outname = "outputest_jetcorrection.root")
{*/
gSystem->mkdir("/sphenix/tg/tg01/jets/jamesj3j3/JetMatchingSubjets/", true);

void Fun4All_JetMatchingSubjets(const char *filelisttruth = "dst_truth_jet.list",
                           const char *filelistcalo = "dst_calo_cluster.list",
                           const char *filelistglobal = "dst_global.list",
                           const char *filelisttruthg4hit = "dst_truth.list",
                           const char *outname = "outputest_jetcorrection.root",
                           int start_event = 0, int end_event = 10000)
{
  //std::string updated_outname = std::string(outname) + "_event" + std::to_string(start_event) + ".root";
  std::string base = std::string(outname);
  if (base.find(".root") != std::string::npos)
    base = base.substr(0, base.find(".root"));
  std::string updated_outname = base + "_event" + std::to_string(start_event) + ".root";


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
  
  JetMatchingSubjets *myJetUnf = new JetMatchingSubjets("AntiKt_Tower_r04_Sub1", "AntiKt_Truth_r04", updated_outname.c_str());  
   HIJetReco();
  

  myJetUnf->setPtRange(1, 100);
  myJetUnf->setEtaRange(-1.1, 1.1);
  //  myJetUnf->doUnsub(1);
  myJetUnf->doTruth(1);
  //  myJetUnf->doSeeds(0);
  //   myJetUnf->set_doSim(false);
  //  myJetUnf->doClusters(0);
  se->registerSubsystem(myJetUnf);
  //We need to remove all DST inputs other than DSTcalo for data

  Fun4AllInputManager *intrue = new Fun4AllDstInputManager("DSTtruth");
  intrue->AddFile(filelisttruth);
  se->registerInputManager(intrue);

  Fun4AllInputManager *in4 = new Fun4AllDstInputManager("DSTg4hit");
  in4->AddFile(filelisttruthg4hit);
  se->registerInputManager(in4);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTcalo");
  in2->AddFile(filelistcalo);
  se->registerInputManager(in2);

  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("DSTglobal");
  in3->AddFile(filelistglobal);
  se->registerInputManager(in3);
  //  se->run(100);
  se->skip(start_event); // Skip events before start_event
  se->run(end_event - start_event);
  se->End();

  gSystem->Exit(0);
  return;

}
