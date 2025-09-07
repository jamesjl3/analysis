#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <GlobalVariables.C>
#include <G4_Input.C>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloWaveformProcessing.h>

#include <calowaveformsim/CaloWaveformSim.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/recoConsts.h>

#include <calotrigger/CaloTriggerEmulator.h>
#include <calowaveformsim/CaloEmulatorTreeMaker.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libtriggervalid.so)
R__LOAD_LIBRARY(libemulatortreemaker.so)
#endif

void Fun4All_Sim_Emulator(
    const int nEvents = 100000,
    const string &inputFile0 = "g4hits.list",
    const string &inputFile1 = "dst_calo_cluster.list",
    const string &inputFile2 = "pedestal.root",

    const string &outputFile = "DST_CALO_WAVEFORM_pp-0000000011-00000.root",
    const string &outdir = ".",
    const string &cdbtag = "ProdA_2023")
{

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  recoConsts *rc = recoConsts::instance();

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", cdbtag);
  rc->set_uint64Flag("TIMESTAMP", 7);
  CDBInterface::instance()->Verbosity(1);

  //===============
  // Input options
  //===============
  // verbosity setting (applies to all input managers)
  Input::VERBOSITY = 1;

  Input::READHITS = true;

  INPUTREADHITS::listfile[0] = inputFile0;
  INPUTREADHITS::listfile[1] = inputFile1;

  InputInit();

  InputRegister();

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  Enable::DSTOUT = false;
  Enable::DSTOUT_COMPRESS = false;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;

  CaloWaveformSim *caloWaveformSim;

  caloWaveformSim= new CaloWaveformSim();
  caloWaveformSim->set_detector_type(CaloTowerDefs::HCALIN);
  caloWaveformSim->set_detector("HCALIN");
  caloWaveformSim->set_nsamples(16);
  caloWaveformSim->set_highgain(false);
  caloWaveformSim->set_timewidth(0.2);
  caloWaveformSim->set_peakpos(6);
  caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  se->registerSubsystem(caloWaveformSim);

  caloWaveformSim= new CaloWaveformSim();
  caloWaveformSim->set_detector_type(CaloTowerDefs::HCALOUT);
  caloWaveformSim->set_detector("HCALOUT");
  caloWaveformSim->set_nsamples(16);
  caloWaveformSim->set_highgain(false);
  caloWaveformSim->set_timewidth(0.2);
  caloWaveformSim->set_peakpos(6);
  caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  se->registerSubsystem(caloWaveformSim);

  caloWaveformSim= new CaloWaveformSim();
  caloWaveformSim->set_detector_type(CaloTowerDefs::CEMC);
  caloWaveformSim->set_detector("CEMC");
  caloWaveformSim->set_nsamples(16);
  caloWaveformSim->set_highgain(false);
  caloWaveformSim->set_timewidth(0.2);
  caloWaveformSim->set_peakpos(6);
  caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  caloWaveformSim->set_calibName("cemc_pi0_twrSlope_v1_default");
  se->registerSubsystem(caloWaveformSim);
  
  CaloTowerBuilder *ca2 = new CaloTowerBuilder();
  ca2->set_detector_type(CaloTowerDefs::HCALOUT);
  ca2->set_nsamples(16);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::FAST);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
  se->registerSubsystem(ca2);

  ca2 = new CaloTowerBuilder();
  ca2->set_detector_type(CaloTowerDefs::HCALIN);
  ca2->set_nsamples(16);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::FAST);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
  se->registerSubsystem(ca2);

  ca2 = new CaloTowerBuilder();
  ca2->set_detector_type(CaloTowerDefs::CEMC);
  ca2->set_nsamples(16);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::FAST);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
  se->registerSubsystem(ca2);

  CaloTowerCalib *calib = new CaloTowerCalib();
  calib->set_detector_type(CaloTowerDefs::HCALOUT);
  calib->set_outputNodePrefix("TOWERSWAVEFORM_CALIB_");
  se->registerSubsystem(calib);

  calib = new CaloTowerCalib();
  calib->set_detector_type(CaloTowerDefs::HCALIN);
  calib->set_outputNodePrefix("TOWERSWAVEFORM_CALIB_");
  se->registerSubsystem(calib);

  calib = new CaloTowerCalib();
  calib->set_detector_type(CaloTowerDefs::CEMC);
  calib->set_outputNodePrefix("TOWERSWAVEFORM_CALIB_");
  se->registerSubsystem(calib);

  CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR");
  te->setTriggerType("JET");
  te->setNSamples(16);
  te->setTriggerSample(6);
  te->setTriggerDelay(4);
  te->Verbosity(0);
  te->setThreshold(1, 4, 6, 8);
//  te->setEmcalLUTFile("/sphenix/user/hanpuj/Jettrigger/calotriggeremulator/macros/emcal_lut.root");
//  te->setHcalinLUTFile("/sphenix/user/hanpuj/Jettrigger/calotriggeremulator/macros/ihcal_lut.root");
//  te->setHcaloutLUTFile("/sphenix/user/hanpuj/Jettrigger/calotriggeremulator/macros/ohcal_lut.root");
  te->setEmcalLUTFile("/sphenix/user/hanpuj/Jettrigger/lut/emcal_ll1_lut.root");
  te->setHcalinLUTFile("/sphenix/user/hanpuj/Jettrigger/lut/hcalin_ll1_lut.root");
  te->setHcaloutLUTFile("/sphenix/user/hanpuj/Jettrigger/lut/hcalout_ll1_lut.root");
  se->registerSubsystem(te);

  InputManagers();

  CaloEmulatorTreeMaker *tt1 = new CaloEmulatorTreeMaker("CaloEmulatorTreemaker","30GeVJetTestR04100k.root");
  tt1->UseCaloTowerBuilder(true);
  tt1->Verbosity(0);
  se->registerSubsystem(tt1);
  
  se->run(nEvents);
  CDBInterface::instance()->Print();  // print used DB files
  se->End();
  se->PrintTimer();
  gSystem->Exit(0);
}
