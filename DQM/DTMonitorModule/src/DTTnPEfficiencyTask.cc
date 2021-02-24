/*
 * \file DTTnPEfficiencyTask.cc
 *
 * \author L. Lunerti - INFN Bologna
 *
 */

#include "DQM/DTMonitorModule/src/DTTnPEfficiencyTask.h"

// Framework
#include "FWCore/Framework/interface/EventSetup.h"

// Geometry
#include "DataFormats/GeometryVector/interface/Pi.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"

//Math
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"

//Root
#include "TH1.h"
#include "TAxis.h"
#include "TRegexp.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <tuple>

DTTnPEfficiencyTask::DTTnPEfficiencyTask(const edm::ParameterSet& config) : 
  m_nEvents(0),
  m_muToken(consumes<reco::MuonCollection>(config.getUntrackedParameter<edm::InputTag>("inputTagMuons"))),
  m_primaryVerticesToken(consumes<std::vector<reco::Vertex>>(config.getUntrackedParameter<edm::InputTag>("inputTagPrimaryVertices"))),
  m_triggerResultsToken(consumes<edm::TriggerResults>(config.getUntrackedParameter<edm::InputTag>("trigResultsTag"))),
  m_triggerEventToken(consumes<trigger::TriggerEvent>(config.getUntrackedParameter<edm::InputTag>("trigEventTag"))),
  m_trigName(config.getUntrackedParameter<std::string>("trigName")),
  m_detailedAnalysis(config.getUntrackedParameter<bool>("detailedAnalysis")),
  m_probeSelector(config.getUntrackedParameter<std::string>("probeCut")),
  m_dxyCut(config.getUntrackedParameter<double>("probeDxyCut")),
  m_dzCut(config.getUntrackedParameter<double>("probeDzCut")),
  m_tagSelector(config.getUntrackedParameter<std::string>("tagCut")),
  m_borderCut(config.getUntrackedParameter<double>("borderCut")),
  m_lowPairMassCut(config.getUntrackedParameter<double>("lowPairMassCut")),
  m_highPairMassCut(config.getUntrackedParameter<double>("highPairMassCut")),
  m_dxCut(config.getUntrackedParameter<double>("dx_cut"))
{

  LogTrace("DTDQM|DTMonitorModule|DTTnPEfficiencyTask")
    << "[DTTnPEfficiencyTask]: Constructor" << std::endl;

}

DTTnPEfficiencyTask::~DTTnPEfficiencyTask() 
{

  LogTrace("DTDQM|DTMonitorModule|DTTnPEfficiencyTask")
    << "[DTTnPEfficiencyTask]: analyzed " << m_nEvents << " events" << std::endl;

}

void DTTnPEfficiencyTask::dqmBeginRun(const edm::Run& run, const edm::EventSetup& context) 
{
  bool changed = true;
  m_hltConfig.init(run,context,"HLT",changed);

  bool enableWildCard = true;

  TString tName = TString(m_trigName);
  TRegexp tNamePattern = TRegexp(tName,enableWildCard);

  for(unsigned iPath=0; iPath<m_hltConfig.size();++iPath)
  {
    TString pathName = TString(m_hltConfig.triggerName(iPath));
    if(pathName.Contains(tNamePattern)){
      m_trigIndices.push_back(static_cast<int>(iPath));
    }
  }
}

void DTTnPEfficiencyTask::bookHistograms(DQMStore::IBooker& iBooker, edm::Run const& run,
					 edm::EventSetup const& context) 
{

  LogTrace("DTDQM|DTMonitorModule|DTTnPEfficiencyTask") 
    << "[DTTnPEfficiencyTask]: bookHistograms" << std::endl;

  if (m_detailedAnalysis) 
    {
      std::string baseDir = topFolder() + "/detailed/";
      iBooker.setCurrentFolder(baseDir);
      
      LogTrace("DTDQM|DTMonitorModule|DTTnPEfficiencyTask")
	<< "[DTTnPEfficiencyTask]: booking histos in " << baseDir << std::endl;

      m_histos["probePt"] = iBooker.book1D("probePt", "probePt;probe p_{T} [GeV];Events", 125, 0., 250.);
      m_histos["probeEta"] = iBooker.book1D("probeEta", "probeEta;probe #eta;Events",24, -2.4, 2.4);
      m_histos["probePhi"] = iBooker.book1D("probePhi", "probePhi;probe #phi; Events",36, -TMath::Pi(), TMath::Pi());
      m_histos["probeNumberOfMatchedStations"] = iBooker.book1D("probeNumberOfMatchedStations", "probeNumberOfMatchedStations;Number of matched stations;Events",5, 0., 5);
      m_histos["pairMass"] = iBooker.book1D("pairMass", "pairMass", 25, 70., 120.);
    }
      
  for (int wheel = -2; wheel <= 2; ++wheel) 
    {
      bookWheelHistos(iBooker, wheel, "Task");
    }

  for (int station = -4; station <= 4; ++station) 
    {
      if (station == 0) continue;
      bookEndcapHistos(iBooker, station, "Task");
    }

  auto baseDir = topFolder() + "Task/";
  iBooker.setCurrentFolder(baseDir);

  MonitorElement* me_DT_pass_allCh = iBooker.book1D("DT_nPassingProbe_allCh", "DT_nPassingProbe_allCh", 20, 0.5, 20.5);
  MonitorElement* me_DT_fail_allCh = iBooker.book1D("DT_nFailingProbe_allCh", "DT_nFailingProbe_allCh", 20, 0.5, 20.5);

  MonitorElement* me_CSC_pass_allCh = iBooker.book2D("CSC_nPassingProbe_allCh", "CSC_nPassingProbe_allCh", 9, -4., 5., 4, 0.,4.5);
  MonitorElement* me_CSC_fail_allCh = iBooker.book2D("CSC_nFailingProbe_allCh", "CSC_nFailingProbe_allCh", 9, -4., 5., 4, 0.,4.5);

  MonitorElement* me_CSC_pass_allCh_1D = iBooker.book1D("CSC_nPassingProbe_allCh_1D", "CSC_nPassingProbe_allCh_1D", 9, -4., 5.);
  MonitorElement* me_CSC_fail_allCh_1D = iBooker.book1D("CSC_nFailingProbe_allCh_1D", "CSC_nFailingProbe_allCh_1D", 9, -4., 5.);

  MonitorElement* me_RPC_barrel_pass_allCh_1D = iBooker.book1D("RPC_nPassingProbe_Barrel_allCh_1D", "RPC_nPassingProbe_Barrel_allCh_1D", 20, 0.5, 20.5);
  MonitorElement* me_RPC_barrel_fail_allCh_1D = iBooker.book1D("RPC_nFailingProbe_Barrel_allCh_1D", "RPC_nFailingProbe_Barrel_allCh_1D", 20, 0.5, 20.5);

  MonitorElement* me_RPC_endcap_pass_allCh_1D = iBooker.book1D("RPC_nPassingProbe_Endcap_allCh_1D", "RPC_nPassingProbe_Endcap_allCh_1D", 9, -4., 5.);
  MonitorElement* me_RPC_endcap_fail_allCh_1D = iBooker.book1D("RPC_nFailingProbe_Endcap_allCh_1D", "RPC_nFailingProbe_Endcap_allCh_1D", 9, -4., 5.);

  me_DT_pass_allCh->setBinLabel(1 , "MB1/YB-2", 1);
  me_DT_pass_allCh->setBinLabel(2 , "MB2/YB-2", 1);
  me_DT_pass_allCh->setBinLabel(3 , "MB3/YB-2", 1);
  me_DT_pass_allCh->setBinLabel(4 , "MB4/YB-2", 1);
  me_DT_pass_allCh->setBinLabel(5 , "MB1/YB-1", 1);
  me_DT_pass_allCh->setBinLabel(6 , "MB2/YB-1", 1);
  me_DT_pass_allCh->setBinLabel(7 , "MB3/YB-1", 1);
  me_DT_pass_allCh->setBinLabel(8 , "MB4/YB-1", 1);
  me_DT_pass_allCh->setBinLabel(9 , "MB1/YB0", 1);
  me_DT_pass_allCh->setBinLabel(10, "MB2/YB0", 1);
  me_DT_pass_allCh->setBinLabel(11, "MB3/YB0", 1);
  me_DT_pass_allCh->setBinLabel(12, "MB4/YB0", 1);
  me_DT_pass_allCh->setBinLabel(13, "MB1/YB1", 1);
  me_DT_pass_allCh->setBinLabel(14, "MB2/YB1", 1);
  me_DT_pass_allCh->setBinLabel(15, "MB3/YB1", 1);
  me_DT_pass_allCh->setBinLabel(16, "MB4/YB1", 1);
  me_DT_pass_allCh->setBinLabel(17, "MB1/YB2", 1);
  me_DT_pass_allCh->setBinLabel(18, "MB2/YB2", 1);
  me_DT_pass_allCh->setBinLabel(19, "MB3/YB2", 1);
  me_DT_pass_allCh->setBinLabel(20, "MB4/YB2", 1);
  me_DT_pass_allCh->setAxisTitle("Number of passing probes", 2);

  me_DT_fail_allCh->setBinLabel(1 , "MB1/YB-2", 1);
  me_DT_fail_allCh->setBinLabel(2 , "MB2/YB-2", 1);
  me_DT_fail_allCh->setBinLabel(3 , "MB3/YB-2", 1);
  me_DT_fail_allCh->setBinLabel(4 , "MB4/YB-2", 1);
  me_DT_fail_allCh->setBinLabel(5 , "MB1/YB-1", 1);
  me_DT_fail_allCh->setBinLabel(6 , "MB2/YB-1", 1);
  me_DT_fail_allCh->setBinLabel(7 , "MB3/YB-1", 1);
  me_DT_fail_allCh->setBinLabel(8 , "MB4/YB-1", 1);
  me_DT_fail_allCh->setBinLabel(9 , "MB1/YB0", 1);
  me_DT_fail_allCh->setBinLabel(10, "MB2/YB0", 1);
  me_DT_fail_allCh->setBinLabel(11, "MB3/YB0", 1);
  me_DT_fail_allCh->setBinLabel(12, "MB4/YB0", 1);
  me_DT_fail_allCh->setBinLabel(13, "MB1/YB1", 1);
  me_DT_fail_allCh->setBinLabel(14, "MB2/YB1", 1);
  me_DT_fail_allCh->setBinLabel(15, "MB3/YB1", 1);
  me_DT_fail_allCh->setBinLabel(16, "MB4/YB1", 1);
  me_DT_fail_allCh->setBinLabel(17, "MB1/YB2", 1);
  me_DT_fail_allCh->setBinLabel(18, "MB2/YB2", 1);
  me_DT_fail_allCh->setBinLabel(19, "MB3/YB2", 1);
  me_DT_fail_allCh->setBinLabel(20, "MB4/YB2", 1);
  me_DT_fail_allCh->setAxisTitle("Number of failing probes", 2);

  me_CSC_pass_allCh->setBinLabel(1, "ME-4", 1);
  me_CSC_pass_allCh->setBinLabel(2, "ME-3", 1);
  me_CSC_pass_allCh->setBinLabel(3, "ME-2", 1);
  me_CSC_pass_allCh->setBinLabel(4, "ME-1", 1);
  me_CSC_pass_allCh->setBinLabel(6, "ME1", 1);
  me_CSC_pass_allCh->setBinLabel(7, "ME2", 1);
  me_CSC_pass_allCh->setBinLabel(8, "ME3", 1);
  me_CSC_pass_allCh->setBinLabel(9, "ME4", 1);
  for (int i=1; i<5; ++i){
    me_CSC_pass_allCh->setBinLabel(i, std::to_string(i), 2);
  }
  me_CSC_pass_allCh->setAxisTitle("Ring", 2);
  me_CSC_pass_allCh->setAxisTitle("Number of passing probes", 3);

  me_CSC_fail_allCh->setBinLabel(1, "ME-4", 1);
  me_CSC_fail_allCh->setBinLabel(2, "ME-3", 1);
  me_CSC_fail_allCh->setBinLabel(3, "ME-2", 1);
  me_CSC_fail_allCh->setBinLabel(4, "ME-1", 1);
  me_CSC_fail_allCh->setBinLabel(6, "ME1", 1);
  me_CSC_fail_allCh->setBinLabel(7, "ME2", 1);
  me_CSC_fail_allCh->setBinLabel(8, "ME3", 1);
  me_CSC_fail_allCh->setBinLabel(9, "ME4", 1);
  for (int i=1; i<5; ++i){
    me_CSC_fail_allCh->setBinLabel(i, std::to_string(i), 2);
  }
  me_CSC_fail_allCh->setAxisTitle("Ring", 2);
  me_CSC_fail_allCh->setAxisTitle("Number of failing probes", 3);

  me_CSC_pass_allCh_1D->setBinLabel(1, "ME-4", 1);
  me_CSC_pass_allCh_1D->setBinLabel(2, "ME-3", 1);
  me_CSC_pass_allCh_1D->setBinLabel(3, "ME-2", 1);
  me_CSC_pass_allCh_1D->setBinLabel(4, "ME-1", 1);
  me_CSC_pass_allCh_1D->setBinLabel(6, "ME1", 1);
  me_CSC_pass_allCh_1D->setBinLabel(7, "ME2", 1);
  me_CSC_pass_allCh_1D->setBinLabel(8, "ME3", 1);
  me_CSC_pass_allCh_1D->setBinLabel(9, "ME4", 1);
  me_CSC_pass_allCh_1D->setAxisTitle("Number of passing probes", 2);

  me_CSC_fail_allCh_1D->setBinLabel(1, "ME-4", 1);
  me_CSC_fail_allCh_1D->setBinLabel(2, "ME-3", 1);
  me_CSC_fail_allCh_1D->setBinLabel(3, "ME-2", 1);
  me_CSC_fail_allCh_1D->setBinLabel(4, "ME-1", 1);
  me_CSC_fail_allCh_1D->setBinLabel(6, "ME1", 1);
  me_CSC_fail_allCh_1D->setBinLabel(7, "ME2", 1);
  me_CSC_fail_allCh_1D->setBinLabel(8, "ME3", 1);
  me_CSC_fail_allCh_1D->setBinLabel(9, "ME4", 1);
  me_CSC_fail_allCh_1D->setAxisTitle("Number of failing probes", 2);

  me_RPC_barrel_pass_allCh_1D->setBinLabel(1 , "RB1/YB-2", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(2 , "RB2/YB-2", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(3 , "RB3/YB-2", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(4 , "RB4/YB-2", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(5 , "RB1/YB-1", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(6 , "RB2/YB-1", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(7 , "RB3/YB-1", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(8 , "RB4/YB-1", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(9 , "RB1/YB0", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(10, "RB2/YB0", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(11, "RB3/YB0", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(12, "RB4/YB0", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(13, "RB1/YB1", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(14, "RB2/YB1", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(15, "RB3/YB1", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(16, "RB4/YB1", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(17, "RB1/YB2", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(18, "RB2/YB2", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(19, "RB3/YB2", 1);
  me_RPC_barrel_pass_allCh_1D->setBinLabel(20, "RB4/YB2", 1);
  me_RPC_barrel_pass_allCh_1D->setAxisTitle("Number of passing probes", 2);

  me_RPC_barrel_fail_allCh_1D->setBinLabel(1 , "RB1/YB-2", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(2 , "RB2/YB-2", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(3 , "RB3/YB-2", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(4 , "RB4/YB-2", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(5 , "RB1/YB-1", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(6 , "RB2/YB-1", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(7 , "RB3/YB-1", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(8 , "RB4/YB-1", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(9 , "RB1/YB0", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(10, "RB2/YB0", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(11, "RB3/YB0", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(12, "RB4/YB0", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(13, "RB1/YB1", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(14, "RB2/YB1", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(15, "RB3/YB1", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(16, "RB4/YB1", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(17, "RB1/YB2", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(18, "RB2/YB2", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(19, "RB3/YB2", 1);
  me_RPC_barrel_fail_allCh_1D->setBinLabel(20, "RB4/YB2", 1);
  me_RPC_barrel_fail_allCh_1D->setAxisTitle("Number of failing probes", 2);

  me_RPC_endcap_pass_allCh_1D->setBinLabel(1, "RE-4", 1);
  me_RPC_endcap_pass_allCh_1D->setBinLabel(2, "RE-3", 1);
  me_RPC_endcap_pass_allCh_1D->setBinLabel(3, "RE-2", 1);
  me_RPC_endcap_pass_allCh_1D->setBinLabel(4, "RE-1", 1);
  me_RPC_endcap_pass_allCh_1D->setBinLabel(6, "RE1", 1);
  me_RPC_endcap_pass_allCh_1D->setBinLabel(7, "RE2", 1);
  me_RPC_endcap_pass_allCh_1D->setBinLabel(8, "RE3", 1);
  me_RPC_endcap_pass_allCh_1D->setBinLabel(9, "RE4", 1);
  me_RPC_endcap_pass_allCh_1D->setAxisTitle("Number of passing probes", 2);

  me_RPC_endcap_fail_allCh_1D->setBinLabel(1, "RE-4", 1);
  me_RPC_endcap_fail_allCh_1D->setBinLabel(2, "RE-3", 1);
  me_RPC_endcap_fail_allCh_1D->setBinLabel(3, "RE-2", 1);
  me_RPC_endcap_fail_allCh_1D->setBinLabel(4, "RE-1", 1);
  me_RPC_endcap_fail_allCh_1D->setBinLabel(6, "RE1", 1);
  me_RPC_endcap_fail_allCh_1D->setBinLabel(7, "RE2", 1);
  me_RPC_endcap_fail_allCh_1D->setBinLabel(8, "RE3", 1);
  me_RPC_endcap_fail_allCh_1D->setBinLabel(9, "RE4", 1);
  me_RPC_endcap_fail_allCh_1D->setAxisTitle("Number of failing probes", 2);

  m_histos["DT_nPassingProbe_allCh"] = me_DT_pass_allCh;
  m_histos["DT_nFailingProbe_allCh"] = me_DT_fail_allCh;

  m_histos["CSC_nPassingProbe_allCh"] = me_CSC_pass_allCh;
  m_histos["CSC_nFailingProbe_allCh"] = me_CSC_fail_allCh;

  m_histos["CSC_nPassingProbe_allCh_1D"] = me_CSC_pass_allCh_1D;
  m_histos["CSC_nFailingProbe_allCh_1D"] = me_CSC_fail_allCh_1D;

  m_histos["RPC_nPassingProbe_Barrel_allCh_1D"] = me_RPC_barrel_pass_allCh_1D;
  m_histos["RPC_nFailingProbe_Barrel_allCh_1D"] = me_RPC_barrel_fail_allCh_1D;

  m_histos["RPC_nPassingProbe_Endcap_allCh_1D"] = me_RPC_endcap_pass_allCh_1D;
  m_histos["RPC_nFailingProbe_Endcap_allCh_1D"] = me_RPC_endcap_fail_allCh_1D;

}

void DTTnPEfficiencyTask::analyze(const edm::Event& event, const edm::EventSetup& context) 
{
  ++m_nEvents;

  edm::Handle<reco::MuonCollection> muons;
  event.getByToken(m_muToken, muons);

  edm::Handle<std::vector<reco::Vertex>> vtxs;
  event.getByToken(m_primaryVerticesToken, vtxs);
  const reco::Vertex & vertex = vtxs->at(0);

  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByToken(m_triggerResultsToken, triggerResults);

  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByToken(m_triggerEventToken, triggerEvent);

  //DT variables
  std::vector<std::vector<int>> probe_coll_DT_wh;
  std::vector<std::vector<int>> probe_coll_DT_sec;
  std::vector<std::vector<int>> probe_coll_DT_sta;
  std::vector<std::vector<float>> probe_coll_DT_dx;
  std::vector<uint8_t> probe_coll_DT_staMatch;
  //CSC variables
  std::vector<std::vector<int>> probe_coll_CSC_zend;
  std::vector<std::vector<int>> probe_coll_CSC_ring;
  std::vector<std::vector<int>> probe_coll_CSC_sta;
  std::vector<std::vector<float>> probe_coll_CSC_dx;
  std::vector<uint8_t> probe_coll_CSC_staMatch;
  //RPC variables
  std::vector<std::vector<int>> probe_coll_RPC_region;
  std::vector<std::vector<int>> probe_coll_RPC_ring;
  std::vector<std::vector<int>> probe_coll_RPC_sta;
  std::vector<std::vector<int>> probe_coll_RPC_sec;
  std::vector<std::vector<int>> probe_coll_RPC_lay;
  std::vector<std::vector<int>> probe_coll_RPC_sub;
  std::vector<std::vector<int>> probe_coll_RPC_roll;
  std::vector<std::vector<float>> probe_coll_RPC_dx;
  std::vector<uint8_t> probe_coll_RPC_staMatch;
  //common tnp variables
  std::vector<unsigned> tag_indices;
  std::vector<unsigned> preSel_probe_indices;
  std::vector<unsigned> probe_indices;

  if (muons.isValid() && vtxs.isValid()) 
  {
    //Is there a better way to initialize two different type variables?
    for (auto[muon,muonColl_index] = std::tuple{std::vector<reco::Muon>::const_iterator{(*muons).begin()},0}; muon!=(*muons).end(); ++muon, ++muonColl_index){

      bool trigMatch = false;
      
      //Getting trigger infos for tag selection
      if(triggerResults.isValid() && triggerEvent.isValid())
      {
        const trigger::TriggerObjectCollection trigObjColl = triggerEvent->getObjects();
        trigMatch = hasTrigger(m_trigIndices,trigObjColl,triggerEvent,*muon);
      }
      
      //Probe selection
      if (m_probeSelector(*muon) &&
          (std::abs(muon->muonBestTrack()->dxy(vertex.position())) < m_dxyCut) &&
          (std::abs(muon->muonBestTrack()->dz(vertex.position()))  < m_dzCut))
      {
        preSel_probe_indices.push_back(muonColl_index);
      }
      //Tag selection
      if (m_tagSelector(*muon) &&
          trigMatch)
      {
        tag_indices.push_back(muonColl_index);
      }

    }//loop over muons
  }

    //Probe selection
    for(const auto i_tag : tag_indices){
      reco::Muon tag = (*muons).at(i_tag);
      float pt_max = 0.;
      int max_pt_idx;
      bool pair_found = false;

      for(const auto i_probe : preSel_probe_indices){
        //Prevent tag and probe to be the same object
        if(i_probe==i_tag) continue;

        reco::Muon preSel_probe = (*muons).at(i_probe);

	int pair_charge_product = tag.charge()*preSel_probe.charge();

	//check if tag+probe pair is compatible with Z decay
	if(pair_charge_product>0) continue;

        float pair_mass = (tag.polarP4()+preSel_probe.polarP4()).M();
	m_histos.find("pairMass")->second->Fill(pair_mass);

	if(pair_mass<m_lowPairMassCut || pair_mass>m_highPairMassCut) continue;

        float pair_pt = (tag.polarP4()+preSel_probe.polarP4()).Pt();
	if (pair_pt > pt_max){
	  pair_found = true;
	  pt_max = pair_pt;
	  max_pt_idx = i_probe;
	}
      }
      if (pair_found){
        probe_indices.push_back(max_pt_idx);
      }
    }

    //Fill probe dx + subdetector coordinates
    for (const auto i : probe_indices){
      //DT variables
      std::vector<int> probe_DT_wh;
      std::vector<int> probe_DT_sec;
      std::vector<int> probe_DT_sta;
      std::vector<float> probe_DT_dx;
      uint8_t DT_stationMatching = 0;
      //CSC variables
      std::vector<int> probe_CSC_zend;
      std::vector<int> probe_CSC_ring;
      std::vector<int> probe_CSC_sta;
      std::vector<float> probe_CSC_dx;
      uint8_t CSC_stationMatching = 0;
      //RPC variables
      std::vector<int> probe_RPC_region;
      std::vector<int> probe_RPC_ring;
      std::vector<int> probe_RPC_sta;
      std::vector<int> probe_RPC_sec;
      std::vector<int> probe_RPC_lay;
      std::vector<int> probe_RPC_sub;
      std::vector<int> probe_RPC_roll;
      std::vector<float> probe_RPC_dx;
      uint8_t RPC_stationMatching = 0;

      //Fill detailed plots
      if (m_detailedAnalysis)
      {
        m_histos.find("probeEta")->second->Fill((*muons).at(i).eta());
        m_histos.find("probePhi")->second->Fill((*muons).at(i).phi());
        m_histos.find("probeNumberOfMatchedStations")->second->Fill((*muons).at(i).numberOfMatchedStations());
        m_histos.find("probePt")->second->Fill((*muons).at(i).pt());
      }

      for (const auto chambMatch : (*muons).at(i).matches()){
        // look in DTs
        if (chambMatch.detector() == MuonSubdetId::DT){

          if (chambMatch.edgeX < m_borderCut && 
              chambMatch.edgeY < m_borderCut){

            DTChamberId chId(chambMatch.id.rawId());

            int wheel   = chId.wheel();
            int sector  = chId.sector();
            int station = chId.station();

            reco::MuonSegmentMatch closest_matchedSegment;
            double smallestDx = 999.;

            for (auto & seg : chambMatch.segmentMatches)
            {
              float dx = std::abs(chambMatch.x - seg.x);
              if (dx < smallestDx){
                smallestDx = dx;
                closest_matchedSegment = seg;
              }
            }

            DT_stationMatching = DT_stationMatching | (1 << (station-1));

            probe_DT_wh.push_back(wheel);
            probe_DT_sec.push_back(sector);
            probe_DT_sta.push_back(station);
            probe_DT_dx.push_back(smallestDx);
              
          }  
        }

        // look in CSCs
        else if (chambMatch.detector() == MuonSubdetId::CSC)
	{
          if (chambMatch.edgeX < m_borderCut && 
              chambMatch.edgeY < m_borderCut){

            CSCDetId chId(chambMatch.id.rawId());

            int zendcap = chId.zendcap();
            int ring    = chId.ring();
            int station = chId.station();

	    reco::MuonSegmentMatch closest_matchedSegment;
            double smallestDx = 999.;
	    for (auto & seg : chambMatch.segmentMatches)
	    {
	      float dx = std::abs(chambMatch.x - seg.x) ;
              if (dx < smallestDx){
	        smallestDx = dx;
	        closest_matchedSegment = seg;
	      }
	    }

	    CSC_stationMatching = CSC_stationMatching | (1 << (station-1));

            probe_CSC_zend.push_back(zendcap);
            probe_CSC_ring.push_back(ring);
            probe_CSC_sta.push_back(station);
            probe_CSC_dx.push_back(smallestDx);
          }  
	}

        // look in RPCs
        else if (chambMatch.detector() == MuonSubdetId::RPC)
	{
          if (chambMatch.edgeX < m_borderCut && 
              chambMatch.edgeY < m_borderCut){

            RPCDetId chId(chambMatch.id.rawId());

            int region    = chId.region(); // barrel if 0, endcap if -/+ 1
            int ring      = chId.ring(); // means wheel in the barrel and ring in the endcap
            int station   = chId.station();
            int sector    = chId.sector();
            int subsector = chId.subsector();
            int layer     = chId.layer();
            int roll      = chId.roll();

	    reco::MuonRPCHitMatch closest_matchedRPCHit;
            double smallestDx = 999.;
	    for (auto & seg : chambMatch.rpcMatches)
	    {
	      float dx = std::abs(chambMatch.x - seg.x);

              if (dx < smallestDx){
	        smallestDx = dx;
	        closest_matchedRPCHit = seg;
	      }
	    }

	    RPC_stationMatching = RPC_stationMatching | (1 << (station-1));

            probe_RPC_region.push_back(region);
            probe_RPC_ring.push_back(ring);
            probe_RPC_sta.push_back(station);
            probe_RPC_sec.push_back(sector);
            probe_RPC_lay.push_back(layer);
            probe_RPC_sub.push_back(subsector);
            probe_RPC_roll.push_back(roll);
            probe_RPC_dx.push_back(smallestDx);
          }  
	}
	else continue;
      }//loop over chamber matches

      //Fill DT variables
      probe_coll_DT_wh.push_back(probe_DT_wh);
      probe_coll_DT_sec.push_back(probe_DT_sec);
      probe_coll_DT_sta.push_back(probe_DT_sta);
      probe_coll_DT_dx.push_back(probe_DT_dx);
      probe_coll_DT_staMatch.push_back(DT_stationMatching);
      //Fill CSC variables
      probe_coll_CSC_zend.push_back(probe_CSC_zend);
      probe_coll_CSC_ring.push_back(probe_CSC_ring);
      probe_coll_CSC_sta.push_back(probe_CSC_sta);
      probe_coll_CSC_dx.push_back(probe_CSC_dx);
      probe_coll_CSC_staMatch.push_back(CSC_stationMatching);
      //Fill RPC variables
      probe_coll_RPC_region.push_back(probe_RPC_region);
      probe_coll_RPC_ring.push_back(probe_RPC_ring);
      probe_coll_RPC_sta.push_back(probe_RPC_sta);
      probe_coll_RPC_sec.push_back(probe_RPC_sec);
      probe_coll_RPC_lay.push_back(probe_RPC_lay);
      probe_coll_RPC_sub.push_back(probe_RPC_sub);
      probe_coll_RPC_roll.push_back(probe_RPC_roll);
      probe_coll_RPC_dx.push_back(probe_RPC_dx);
      probe_coll_RPC_staMatch.push_back(RPC_stationMatching);
    }//loop over probe collection

    //Loop over probes
    for(unsigned i=0; i<probe_indices.size(); ++i){
      uint8_t DT_matchPatt = probe_coll_DT_staMatch.at(i);
      uint8_t CSC_matchPatt = probe_coll_CSC_staMatch.at(i);
      uint8_t RPC_matchPatt = probe_coll_RPC_staMatch.at(i);

      //Loop over DT matches
      unsigned nDT_matches = probe_coll_DT_wh.at(i).size();
      for(unsigned j=0; j<nDT_matches; ++j){
        //DT variables
        int DT_wheel   = probe_coll_DT_wh.at(i).at(j);
        int DT_station = probe_coll_DT_sta.at(i).at(j);
        int DT_sector  = probe_coll_DT_sec.at(i).at(j);
        float DT_dx    = probe_coll_DT_dx.at(i).at(j);
        
        //Fill DT plots
        if ((DT_matchPatt & (1<<(DT_station-1))) != 0 && //avoids 0 station matching
            (DT_matchPatt & (1<<(DT_station-1))) !=DT_matchPatt) //avoids matching with the station under consideration only
        {
          if (DT_dx < m_dxCut)
          {
            std::string hName = std::string("DT_nPassingProbePerCh_W") + std::to_string(DT_wheel);  
            m_histos.find(hName)->second->Fill(DT_sector, DT_station);
            m_histos.find("DT_nPassingProbe_allCh")->second->Fill((DT_station) + 4*(DT_wheel+2));
          }
          else
          {
            std::string hName = std::string("DT_nFailingProbePerCh_W") + std::to_string(DT_wheel);  
            m_histos.find(hName)->second->Fill(DT_sector, DT_station);
            m_histos.find("DT_nFailingProbe_allCh")->second->Fill((DT_station) + 4*(DT_wheel+2));
          }
        }
      }

      //Loop over CSC matches
      unsigned nCSC_matches = probe_coll_CSC_zend.at(i).size();
      for(unsigned j=0; j<nCSC_matches; ++j){
        int CSC_zendcap = probe_coll_CSC_zend.at(i).at(j);
        int CSC_sta     = probe_coll_CSC_sta.at(i).at(j);
        int CSC_ring    = probe_coll_CSC_ring.at(i).at(j);
        float CSC_dx    = probe_coll_CSC_dx.at(i).at(j);

        //Fill CSC plots
        if ((CSC_matchPatt & (1<<(CSC_sta-1))) != 0 && //avoids 0 station matching
            (CSC_matchPatt & (1<<(CSC_sta-1))) !=CSC_matchPatt) //avoids matching with the station under consideration only
        {
          if (CSC_dx < m_dxCut)
          {
            m_histos.find("CSC_nPassingProbe_allCh")->second->Fill(CSC_zendcap*CSC_sta, CSC_ring);
            m_histos.find("CSC_nPassingProbe_allCh_1D")->second->Fill(CSC_zendcap*CSC_sta);
          }
          else
          {
            m_histos.find("CSC_nFailingProbe_allCh")->second->Fill(CSC_zendcap*CSC_sta, CSC_ring);
            m_histos.find("CSC_nFailingProbe_allCh_1D")->second->Fill(CSC_zendcap*CSC_sta);
          }

        }
      }

      //Loop over RPC matches
      unsigned nRPC_matches = probe_coll_RPC_region.at(i).size();
      for(unsigned j=0; j<nRPC_matches; ++j){
        //RPC variables
        int RPC_region = probe_coll_RPC_region.at(i).at(j);
        int RPC_ring   = probe_coll_RPC_ring.at(i).at(j);
        int RPC_sta    = probe_coll_RPC_sta.at(i).at(j);
        int RPC_sec    = probe_coll_RPC_sec.at(i).at(j);
        int RPC_lay    = probe_coll_RPC_lay.at(i).at(j);
        int RPC_subsec = probe_coll_RPC_sub.at(i).at(j);
        int RPC_roll   = probe_coll_RPC_roll.at(i).at(j);
        float RPC_dx   = probe_coll_RPC_dx.at(i).at(j);

        //Fill RPC plots
        if ((RPC_matchPatt & (1<<(RPC_sta-1))) != 0) //avoids 0 station matching
        {
          if (RPC_dx < m_dxCut)
          {
            //Barrel region
            if (RPC_region == 0){
              int barrel_histo_xcoord = RPC_sec;
              int barrel_histo_ycoord = get_barrel_histo_ycoord(RPC_ring,RPC_sta,RPC_sec,RPC_lay,RPC_subsec,RPC_roll);

              std::string hName = std::string("RPC_nPassingProbePerRoll_Barrel_W") + std::to_string(RPC_ring);  
              m_histos.find(hName)->second->Fill(barrel_histo_xcoord, barrel_histo_ycoord);

              std::string hName_1D = std::string("RPC_nPassingProbePerRoll_Barrel_1D_W") + std::to_string(RPC_ring);  
              m_histos.find(hName_1D)->second->Fill(barrel_histo_ycoord);

              m_histos.find("RPC_nPassingProbe_Barrel_allCh_1D")->second->Fill((RPC_sta) + 4*(RPC_ring+2));
            }
            //Endcap region
            else{
              int endcap_histo_xcoord = (6*(RPC_sec-1)) + RPC_subsec;
              int endcap_histo_ycoord = (3*(RPC_ring-2))+ RPC_roll;

              std::string hName = std::string("RPC_nPassingProbePerRoll_Endcap_Sta") + std::to_string(RPC_sta*RPC_region);  
              m_histos.find(hName)->second->Fill(endcap_histo_xcoord, endcap_histo_ycoord);

              std::string hName_1D = std::string("RPC_nPassingProbePerRoll_Endcap_1D_Sta") + std::to_string(RPC_sta*RPC_region);  
              m_histos.find(hName_1D)->second->Fill(endcap_histo_ycoord);

              m_histos.find("RPC_nPassingProbe_Endcap_allCh_1D")->second->Fill(RPC_region*RPC_sta);
            }
          }
          else
          {
            //Barrel region
            if (RPC_region == 0){
              int barrel_histo_xcoord = RPC_sec;
              int barrel_histo_ycoord = get_barrel_histo_ycoord(RPC_ring,RPC_sta,RPC_sec,RPC_lay,RPC_subsec,RPC_roll);

              std::string hName = std::string("RPC_nFailingProbePerRoll_Barrel_W") + std::to_string(RPC_ring);  
              m_histos.find(hName)->second->Fill(barrel_histo_xcoord, barrel_histo_ycoord);

              std::string hName_1D = std::string("RPC_nFailingProbePerRoll_Barrel_1D_W") + std::to_string(RPC_ring);  
              m_histos.find(hName_1D)->second->Fill(barrel_histo_ycoord);

              m_histos.find("RPC_nFailingProbe_Barrel_allCh_1D")->second->Fill((RPC_sta) + 4*(RPC_ring+2));
            }
            //Endcap region
            else{
              int endcap_histo_xcoord = (6*(RPC_sec-1)) + RPC_subsec;
              int endcap_histo_ycoord = (3*(RPC_ring-2))+ RPC_roll;

              std::string hName = std::string("RPC_nFailingProbePerRoll_Endcap_Sta") + std::to_string(RPC_sta*RPC_region);  
              m_histos.find(hName)->second->Fill(endcap_histo_xcoord, endcap_histo_ycoord);

              std::string hName_1D = std::string("RPC_nFailingProbePerRoll_Endcap_1D_Sta") + std::to_string(RPC_sta*RPC_region);  
              m_histos.find(hName_1D)->second->Fill(endcap_histo_ycoord);

              m_histos.find("RPC_nFailingProbe_Endcap_allCh_1D")->second->Fill(RPC_region*RPC_sta);
            }
          }
        }
      }
    }
}

void DTTnPEfficiencyTask::bookWheelHistos(DQMStore::IBooker& iBooker, 
					  int wheel, std::string folder) 
{

  auto baseDir = topFolder() + folder + "/";
  iBooker.setCurrentFolder(baseDir);

  LogTrace("DTDQM|DTMonitorModule|DTTnPEfficiencyTask")
    << "[DTTnPEfficiencyTask]: booking histos in " << baseDir << std::endl;

  auto hName_DT_pass = std::string("DT_nPassingProbePerCh_W") + std::to_string(wheel);    
  auto hName_DT_fail = std::string("DT_nFailingProbePerCh_W") + std::to_string(wheel);    

  auto hName_RPC_pass = std::string("RPC_nPassingProbePerRoll_Barrel_W") + std::to_string(wheel);    
  auto hName_RPC_fail = std::string("RPC_nFailingProbePerRoll_Barrel_W") + std::to_string(wheel);    

  auto hName_RPC_pass_1D = std::string("RPC_nPassingProbePerRoll_Barrel_1D_W") + std::to_string(wheel);    
  auto hName_RPC_fail_1D = std::string("RPC_nFailingProbePerRoll_Barrel_1D_W") + std::to_string(wheel);    

  MonitorElement* me_DT_pass = iBooker.book2D(hName_DT_pass.c_str(), hName_DT_pass.c_str(), 14, 0.5, 14.5, 4, 0., 4.5);
  MonitorElement* me_DT_fail = iBooker.book2D(hName_DT_fail.c_str(), hName_DT_fail.c_str(), 14, 0.5, 14.5, 4, 0., 4.5);

  MonitorElement* me_RPC_pass = iBooker.book2D(hName_RPC_pass.c_str(), hName_RPC_pass.c_str(), 12, 0.5, 12.5, 21, 0., 21.5);
  MonitorElement* me_RPC_fail = iBooker.book2D(hName_RPC_fail.c_str(), hName_RPC_fail.c_str(), 12, 0.5, 12.5, 21, 0., 21.5);

  MonitorElement* me_RPC_pass_1D = iBooker.book1D(hName_RPC_pass_1D.c_str(), hName_RPC_pass_1D.c_str(), 21, 0., 21.5);
  MonitorElement* me_RPC_fail_1D = iBooker.book1D(hName_RPC_fail_1D.c_str(), hName_RPC_fail_1D.c_str(), 21, 0., 21.5);

  me_DT_pass->setBinLabel(1, "MB1", 2);
  me_DT_pass->setBinLabel(2, "MB2", 2);
  me_DT_pass->setBinLabel(3, "MB3", 2);
  me_DT_pass->setBinLabel(4, "MB4", 2);
  for (int i=1; i<15; ++i){
    me_DT_pass->setBinLabel(i, std::to_string(i), 1);
  }
  me_DT_pass->setAxisTitle("Sector", 1);
  me_DT_pass->setAxisTitle("Number of passing probes", 3);

  me_DT_fail->setBinLabel(1, "MB1", 2);
  me_DT_fail->setBinLabel(2, "MB2", 2);
  me_DT_fail->setBinLabel(3, "MB3", 2);
  me_DT_fail->setBinLabel(4, "MB4", 2);
  for (int i=1; i<15; ++i){
    me_DT_fail->setBinLabel(i, std::to_string(i), 1);
  }
  me_DT_fail->setAxisTitle("Sector", 1);
  me_DT_fail->setAxisTitle("Number of failing probes", 3);

  me_RPC_pass->setBinLabel(1,  "RB1in B", 2);
  me_RPC_pass->setBinLabel(2,  "RB1in F", 2);
  me_RPC_pass->setBinLabel(3,  "RB1out B", 2);
  me_RPC_pass->setBinLabel(4,  "RB1out F", 2);
  if(std::abs(wheel)<2){
    me_RPC_pass->setBinLabel(5,  "RB2in B", 2);
    me_RPC_pass->setBinLabel(6,  "RB2in M", 2);
    me_RPC_pass->setBinLabel(7,  "RB2in F", 2);
    me_RPC_pass->setBinLabel(8,  "RB2out B", 2);
    me_RPC_pass->setBinLabel(9,  "RB2out F", 2);
  }
  else{
    me_RPC_pass->setBinLabel(5,  "RB2in B", 2);
    me_RPC_pass->setBinLabel(6,  "RB2in F", 2);
    me_RPC_pass->setBinLabel(7,  "RB2out B", 2);
    me_RPC_pass->setBinLabel(8,  "RB2out M", 2);
    me_RPC_pass->setBinLabel(9,  "RB2out F", 2);
  }
  me_RPC_pass->setBinLabel(10, "RB3- B", 2);
  me_RPC_pass->setBinLabel(11, "RB3- F", 2);
  me_RPC_pass->setBinLabel(12, "RB3+ B", 2);
  me_RPC_pass->setBinLabel(13, "RB3+ F", 2);
  me_RPC_pass->setBinLabel(14, "RB4- B", 2);
  me_RPC_pass->setBinLabel(15, "RB4- F", 2);
  me_RPC_pass->setBinLabel(16, "RB4+ B", 2);
  me_RPC_pass->setBinLabel(17, "RB4+ F", 2);
  me_RPC_pass->setBinLabel(18, "RB4-- B", 2);
  me_RPC_pass->setBinLabel(19, "RB4-- F", 2);
  me_RPC_pass->setBinLabel(20, "RB4++ B", 2);
  me_RPC_pass->setBinLabel(21, "RB4++ F", 2);
  for (int i=1; i<13; ++i){
    me_RPC_pass->setBinLabel(i, std::to_string(i), 1);
  }
  me_RPC_pass->setAxisTitle("Sector", 1);
  me_RPC_pass->setAxisTitle("Number of passing probes", 3);

  me_RPC_fail->setBinLabel(1,  "RB1in B", 2);
  me_RPC_fail->setBinLabel(2,  "RB1in F", 2);
  me_RPC_fail->setBinLabel(3,  "RB1out B", 2);
  me_RPC_fail->setBinLabel(4,  "RB1out F", 2);
  if(std::abs(wheel)<2){
    me_RPC_fail->setBinLabel(5,  "RB2in B", 2);
    me_RPC_fail->setBinLabel(6,  "RB2in M", 2);
    me_RPC_fail->setBinLabel(7,  "RB2in F", 2);
    me_RPC_fail->setBinLabel(8,  "RB2out B", 2);
    me_RPC_fail->setBinLabel(9,  "RB2out F", 2);
  }
  else{
    me_RPC_fail->setBinLabel(5,  "RB2in B", 2);
    me_RPC_fail->setBinLabel(6,  "RB2in F", 2);
    me_RPC_fail->setBinLabel(7,  "RB2out B", 2);
    me_RPC_fail->setBinLabel(8,  "RB2out M", 2);
    me_RPC_fail->setBinLabel(9,  "RB2out F", 2);
  }
  me_RPC_fail->setBinLabel(10, "RB3- B", 2);
  me_RPC_fail->setBinLabel(11, "RB3- F", 2);
  me_RPC_fail->setBinLabel(12, "RB3+ B", 2);
  me_RPC_fail->setBinLabel(13, "RB3+ F", 2);
  me_RPC_fail->setBinLabel(14, "RB4- B", 2);
  me_RPC_fail->setBinLabel(15, "RB4- F", 2);
  me_RPC_fail->setBinLabel(16, "RB4+ B", 2);
  me_RPC_fail->setBinLabel(17, "RB4+ F", 2);
  me_RPC_fail->setBinLabel(18, "RB4-- B", 2);
  me_RPC_fail->setBinLabel(19, "RB4-- F", 2);
  me_RPC_fail->setBinLabel(20, "RB4++ B", 2);
  me_RPC_fail->setBinLabel(21, "RB4++ F", 2);
  for (int i=1; i<13; ++i){
    me_RPC_fail->setBinLabel(i, std::to_string(i), 1);
  }
  me_RPC_fail->setAxisTitle("Sector", 1);
  me_RPC_fail->setAxisTitle("Number of failing probes", 3);

  me_RPC_pass_1D->setBinLabel(1,  "RB1in B", 1);
  me_RPC_pass_1D->setBinLabel(2,  "RB1in F", 1);
  me_RPC_pass_1D->setBinLabel(3,  "RB1out B", 1);
  me_RPC_pass_1D->setBinLabel(4,  "RB1out F", 1);
  if(std::abs(wheel)<2){
    me_RPC_pass_1D->setBinLabel(5,  "RB2in B", 1);
    me_RPC_pass_1D->setBinLabel(6,  "RB2in M", 1);
    me_RPC_pass_1D->setBinLabel(7,  "RB2in F", 1);
    me_RPC_pass_1D->setBinLabel(8,  "RB2out B", 1);
    me_RPC_pass_1D->setBinLabel(9,  "RB2out F", 1);
  }
  else{
    me_RPC_pass_1D->setBinLabel(5,  "RB2in B", 1);
    me_RPC_pass_1D->setBinLabel(6,  "RB2in F", 1);
    me_RPC_pass_1D->setBinLabel(7,  "RB2out B", 1);
    me_RPC_pass_1D->setBinLabel(8,  "RB2out M", 1);
    me_RPC_pass_1D->setBinLabel(9,  "RB2out F", 1);
  }
  me_RPC_pass_1D->setBinLabel(10, "RB3- B", 1);
  me_RPC_pass_1D->setBinLabel(11, "RB3- F", 1);
  me_RPC_pass_1D->setBinLabel(12, "RB3+ B", 1);
  me_RPC_pass_1D->setBinLabel(13, "RB3+ F", 1);
  me_RPC_pass_1D->setBinLabel(14, "RB4- B", 1);
  me_RPC_pass_1D->setBinLabel(15, "RB4- F", 1);
  me_RPC_pass_1D->setBinLabel(16, "RB4+ B", 1);
  me_RPC_pass_1D->setBinLabel(17, "RB4+ F", 1);
  me_RPC_pass_1D->setBinLabel(18, "RB4-- B", 1);
  me_RPC_pass_1D->setBinLabel(19, "RB4-- F", 1);
  me_RPC_pass_1D->setBinLabel(20, "RB4++ B", 1);
  me_RPC_pass_1D->setBinLabel(21, "RB4++ F", 1);
  me_RPC_pass->setAxisTitle("Number of passing probes", 2);

  me_RPC_fail_1D->setBinLabel(1,  "RB1in B", 1);
  me_RPC_fail_1D->setBinLabel(2,  "RB1in F", 1);
  me_RPC_fail_1D->setBinLabel(3,  "RB1out B", 1);
  me_RPC_fail_1D->setBinLabel(4,  "RB1out F", 1);
  if(std::abs(wheel)<2){
    me_RPC_fail_1D->setBinLabel(5,  "RB2in B", 1);
    me_RPC_fail_1D->setBinLabel(6,  "RB2in M", 1);
    me_RPC_fail_1D->setBinLabel(7,  "RB2in F", 1);
    me_RPC_fail_1D->setBinLabel(8,  "RB2out B", 1);
    me_RPC_fail_1D->setBinLabel(9,  "RB2out F", 1);
  }
  else{
    me_RPC_fail_1D->setBinLabel(5,  "RB2in B", 1);
    me_RPC_fail_1D->setBinLabel(6,  "RB2in F", 1);
    me_RPC_fail_1D->setBinLabel(7,  "RB2out B", 1);
    me_RPC_fail_1D->setBinLabel(8,  "RB2out M", 1);
    me_RPC_fail_1D->setBinLabel(9,  "RB2out F", 1);
  }
  me_RPC_fail_1D->setBinLabel(10, "RB3- B", 1);
  me_RPC_fail_1D->setBinLabel(11, "RB3- F", 1);
  me_RPC_fail_1D->setBinLabel(12, "RB3+ B", 1);
  me_RPC_fail_1D->setBinLabel(13, "RB3+ F", 1);
  me_RPC_fail_1D->setBinLabel(14, "RB4- B", 1);
  me_RPC_fail_1D->setBinLabel(15, "RB4- F", 1);
  me_RPC_fail_1D->setBinLabel(16, "RB4+ B", 1);
  me_RPC_fail_1D->setBinLabel(17, "RB4+ F", 1);
  me_RPC_fail_1D->setBinLabel(18, "RB4-- B", 1);
  me_RPC_fail_1D->setBinLabel(19, "RB4-- F", 1);
  me_RPC_fail_1D->setBinLabel(20, "RB4++ B", 1);
  me_RPC_fail_1D->setBinLabel(21, "RB4++ F", 1);
  me_RPC_fail_1D->setAxisTitle("Number of failing probes", 2);

  m_histos[hName_DT_pass] = me_DT_pass;
  m_histos[hName_DT_fail] = me_DT_fail;

  m_histos[hName_RPC_pass] = me_RPC_pass;
  m_histos[hName_RPC_fail] = me_RPC_fail;

  m_histos[hName_RPC_pass_1D] = me_RPC_pass_1D;
  m_histos[hName_RPC_fail_1D] = me_RPC_fail_1D;
}

void DTTnPEfficiencyTask::bookEndcapHistos(DQMStore::IBooker& iBooker, 
					  int station, std::string folder) 
{

  auto baseDir = topFolder() + folder + "/";
  iBooker.setCurrentFolder(baseDir);

  LogTrace("DTDQM|DTMonitorModule|DTTnPEfficiencyTask")
    << "[DTTnPEfficiencyTask]: booking histos in " << baseDir << std::endl;

  auto hName_RPC_pass = std::string("RPC_nPassingProbePerRoll_Endcap_Sta") + std::to_string(station);    
  auto hName_RPC_fail = std::string("RPC_nFailingProbePerRoll_Endcap_Sta") + std::to_string(station);    

  auto hName_RPC_pass_1D = std::string("RPC_nPassingProbePerRoll_Endcap_1D_Sta") + std::to_string(station);    
  auto hName_RPC_fail_1D = std::string("RPC_nFailingProbePerRoll_Endcap_1D_Sta") + std::to_string(station);    

  MonitorElement* me_RPC_pass = iBooker.book2D(hName_RPC_pass.c_str(), hName_RPC_pass.c_str(), 36, 0.5, 36.5, 6, 0.5, 6.5);
  MonitorElement* me_RPC_fail = iBooker.book2D(hName_RPC_fail.c_str(), hName_RPC_fail.c_str(), 36, 0.5, 36.5, 6, 0.5, 6.5);

  MonitorElement* me_RPC_pass_1D = iBooker.book1D(hName_RPC_pass_1D.c_str(), hName_RPC_pass_1D.c_str(), 6, 0.5, 6.5);
  MonitorElement* me_RPC_fail_1D = iBooker.book1D(hName_RPC_fail_1D.c_str(), hName_RPC_fail_1D.c_str(), 6, 0.5, 6.5);

  me_RPC_pass->setBinLabel(1, "R2_A", 2);
  me_RPC_pass->setBinLabel(2, "R2_B", 2);
  me_RPC_pass->setBinLabel(3, "R2_C", 2);
  me_RPC_pass->setBinLabel(4, "R3_A", 2);
  me_RPC_pass->setBinLabel(5, "R3_B", 2);
  me_RPC_pass->setBinLabel(6, "R3_C", 2);
  for (int i=1; i<37; ++i){
    me_RPC_pass->setBinLabel(i, std::to_string(i), 1);
  }
  me_RPC_pass->setAxisTitle("Sector", 1);
  me_RPC_pass->setAxisTitle("Number of passing probes", 3);
  me_RPC_pass->setTitle("RE" + std::to_string(station));

  me_RPC_fail->setBinLabel(1, "R2_A", 2);
  me_RPC_fail->setBinLabel(2, "R2_B", 2);
  me_RPC_fail->setBinLabel(3, "R2_C", 2);
  me_RPC_fail->setBinLabel(4, "R3_A", 2);
  me_RPC_fail->setBinLabel(5, "R3_B", 2);
  me_RPC_fail->setBinLabel(6, "R3_C", 2);
  for (int i=1; i<37; ++i){
    me_RPC_fail->setBinLabel(i, std::to_string(i), 1);
  }
  me_RPC_fail->setAxisTitle("Sector", 1);
  me_RPC_fail->setAxisTitle("Number of failing probes", 3);
  me_RPC_fail->setTitle("RE" + std::to_string(station));

  me_RPC_pass_1D->setBinLabel(1, "R2_A", 1);
  me_RPC_pass_1D->setBinLabel(2, "R2_B", 1);
  me_RPC_pass_1D->setBinLabel(3, "R2_C", 1);
  me_RPC_pass_1D->setBinLabel(4, "R3_A", 1);
  me_RPC_pass_1D->setBinLabel(5, "R3_B", 1);
  me_RPC_pass_1D->setBinLabel(6, "R3_C", 1);
  me_RPC_pass_1D->setAxisTitle("Number of passing probes", 2);
  me_RPC_pass_1D->setTitle("RE" + std::to_string(station));

  me_RPC_fail_1D->setBinLabel(1, "R2_A", 1);
  me_RPC_fail_1D->setBinLabel(2, "R2_B", 1);
  me_RPC_fail_1D->setBinLabel(3, "R2_C", 1);
  me_RPC_fail_1D->setBinLabel(4, "R3_A", 1);
  me_RPC_fail_1D->setBinLabel(5, "R3_B", 1);
  me_RPC_fail_1D->setBinLabel(6, "R3_C", 1);
  me_RPC_fail_1D->setAxisTitle("Number of failing probes", 2);
  me_RPC_fail_1D->setTitle("RE" + std::to_string(station));

  m_histos[hName_RPC_pass] = me_RPC_pass;
  m_histos[hName_RPC_fail] = me_RPC_fail;

  m_histos[hName_RPC_pass_1D] = me_RPC_pass_1D;
  m_histos[hName_RPC_fail_1D] = me_RPC_fail_1D;
}
bool DTTnPEfficiencyTask::hasTrigger(std::vector<int> & trigIndices,
	                             const trigger::TriggerObjectCollection & trigObjs,
	                             edm::Handle<trigger::TriggerEvent> & trigEvent,
	                             const reco::Muon & muon)
{
  float dRmatch = 999.;
  for (int trigIdx : trigIndices)
  {
    const std::vector<std::string> trigModuleLabels = m_hltConfig.moduleLabels(trigIdx); 
    const unsigned trigModuleIndex = trigModuleLabels.size() - 2;
    const unsigned hltFilterIndex = trigEvent->filterIndex(edm::InputTag(trigModuleLabels[trigModuleIndex],"","HLT"));
    if(hltFilterIndex < trigEvent->sizeFilters())
    {
      const trigger::Keys keys = trigEvent->filterKeys(hltFilterIndex);
      const trigger::Vids vids = trigEvent->filterIds(hltFilterIndex);
      const unsigned nTriggers = vids.size();

      for (unsigned iTrig=0; iTrig<nTriggers; ++iTrig)
      {
        trigger::TriggerObject trigObj = trigObjs[keys[iTrig]];
        float dR = deltaR(muon,trigObj);
        if(dR < dRmatch)
        dRmatch = dR;
      }
    }
  }
  return dRmatch < 0.1;
}

int DTTnPEfficiencyTask::get_barrel_histo_ycoord(int ring,
                                                 int station,
						 int sector,
						 int layer,
						 int subsector,
						 int roll)
{

  int ycoord;

  if (station < 3){
    //There are three rolls in RB2in for wheel=-1,0,+1 and in RB2out for wheel=-2,+2
    bool three_rolls = station == 2 && ((std::abs(ring)>1 && layer ==2) || (std::abs(ring)<2 && layer ==1));
    
    int layer_roll;
    
    if (!three_rolls){
      roll = roll>1 ? 2 : 1;
      int a = station==2 && std::abs(ring)<2 && layer==2 ? 3 : 2;
      layer_roll = (a*(layer-1)) + roll;
    }
    else{
      layer_roll = (2*(layer-1)) + roll;
    }
    
    ycoord = (4*(station-1)) + layer_roll;
  }
  else if (station == 3){
    roll = roll>1 ? 2 : 1;
    ycoord = 9 + (4*(station-3)) + (2*(subsector-1)) + roll;
  }
  else{
    int my_subsector = subsector;
    //The numbering scheme of subsector in sector 4 
    //of station 4 does not match the bins order in the plot:
    //_____SUBSECTOR_____|_______BIN_ORDERING_____
    // ++ --> subsector 4| RB4++ --> my_subsector 4                       
    //  + --> subsector 3| RB4-- --> my_subsector 3                   
    //  - --> subsector 2| RB4+  --> my_subsector 2                   
    // -- --> subsector 1| RB4-  --> my_subsector 1                   

    if (sector == 4){
      switch (subsector){
        case 1: my_subsector = 3;
	        break;
        case 2: my_subsector = 1;
	        break;
        case 3: my_subsector = 2;
	        break;
        case 4: my_subsector = 4;
	        break;
      }
    }
    roll = roll>1 ? 2 : 1;
    ycoord = 9 + (4*(station-3)) + (2*(my_subsector-1)) + roll;
  }

  return ycoord;
}
