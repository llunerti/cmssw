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

//Root
#include "TH1.h"
#include "TAxis.h"
#include "TRegexp.h"

#include <sstream>
#include <iostream>
#include <fstream>

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
  m_highPairMassCut(config.getUntrackedParameter<double>("highPairMassCut"))
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

      m_histos["probePt"] = iBooker.book1D("probePt", "probePt", 125, 0., 250.);
      m_histos["probeEta"] = iBooker.book1D("probeEta", "probeEta",24, -1.2, 1.2);
      m_histos["probePhi"] = iBooker.book1D("probePhi", "probePhi",36, -TMath::Pi(), TMath::Pi());
      m_histos["probeNumberOfMatchedStations"] = iBooker.book1D("probeNumberOfMatchedStations", "probeNumberOfMatchedStations",5, 0., 5);
      m_histos["pairMass"] = iBooker.book1D("pairMass", "pairMass", 25, 70., 120.);
    }
      
  for (int wheel = -2; wheel <= 2; ++wheel) 
    {
      bookWheelHistos(iBooker, wheel, "Task");
    }

  auto baseDir = topFolder() + "Task/";
  iBooker.setCurrentFolder(baseDir);

  MonitorElement* me_DT_pass_allCh = iBooker.book1D("DT_nPassingProbe_allCh", "DT_nPassingProbe_allCh", 20, 0.5, 20.5);
  MonitorElement* me_DT_fail_allCh = iBooker.book1D("DT_nFailingProbe_allCh", "DT_nFailingProbe_allCh", 20, 0.5, 20.5);

  MonitorElement* me_CSC_pass_allCh = iBooker.book2D("CSC_nPassingProbe_allCh", "CSC_nPassingProbe_allCh", 9, -4., 5., 4, 0.,4.5);
  MonitorElement* me_CSC_fail_allCh = iBooker.book2D("CSC_nFailingProbe_allCh", "CSC_nFailingProbe_allCh", 9, -4., 5., 4, 0.,4.5);

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
  me_CSC_fail_allCh->setAxisTitle("Ring", 2);
  me_CSC_fail_allCh->setAxisTitle("Number of failing probes", 3);

  m_histos["DT_nPassingProbe_allCh"] = me_DT_pass_allCh;
  m_histos["DT_nFailingProbe_allCh"] = me_DT_fail_allCh;

  m_histos["CSC_nPassingProbe_allCh"] = me_CSC_pass_allCh;
  m_histos["CSC_nFailingProbe_allCh"] = me_CSC_fail_allCh;

}

void DTTnPEfficiencyTask::analyze(const edm::Event& event, const edm::EventSetup& context) 
{
  ++m_nEvents;

  bool pairFound = false;

  edm::Handle<reco::MuonCollection> muons;
  event.getByToken(m_muToken, muons);

  edm::Handle<std::vector<reco::Vertex>> vtxs;
  event.getByToken(m_primaryVerticesToken, vtxs);
  const reco::Vertex & vertex = vtxs->at(0);

  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByToken(m_triggerResultsToken, triggerResults);

  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByToken(m_triggerEventToken, triggerEvent);

  reco::Muon dummy_muon {};

  std::vector<reco::Muon> tag_muons;
  std::vector<uint8_t> tag_muons_DT_stationMatching;
  std::vector<uint8_t> tag_muons_CSC_stationMatching;

  std::vector<reco::Muon> probe_muons;
  std::vector<uint8_t> probe_muons_DT_stationMatching;
  std::vector<uint8_t> probe_muons_CSC_stationMatching;
  std::vector<bool> probe_muons_firesIsoTrig;
  bool probeFiresTrig = false;


  std::vector<std::vector<std::array<int,3>>> probe_muons_DT_WhSecSta;
  std::vector<std::vector<std::array<int,3>>> tag_muons_DT_WhSecSta;
  std::vector<std::vector<float>> probe_muons_DT_dx;
  std::vector<std::vector<float>> tag_muons_DT_dx;

  std::vector<std::array<int,3>> dummy_DT_WhSecSta {{{-999,-999,-999}}};
  std::vector<std::array<int,3>> dummy_CSC_ZStaRing {{{-999,-999,-999}}};
  std::vector<float> dummy_dx {{999.}};

  std::vector<std::vector<std::array<int,3>>> probe_muons_CSC_ZStaRing;
  std::vector<std::vector<std::array<int,3>>> tag_muons_CSC_ZStaRing;
  std::vector<std::vector<float>> probe_muons_CSC_dx;
  std::vector<std::vector<float>> tag_muons_CSC_dx;

  //Selected tnp pair 
  std::pair<reco::Muon,reco::Muon> tnpPair_maxPt;
  std::pair<uint8_t,uint8_t> tnpPair_DT_stationMatching_maxPt;
  std::pair<uint8_t,uint8_t> tnpPair_CSC_stationMatching_maxPt;
  std::pair<std::vector<std::array<int,3>>,std::vector<std::array<int,3>>> tnpPair_DT_WhSecSta_maxPt;
  std::pair<std::vector<std::array<int,3>>,std::vector<std::array<int,3>>> tnpPair_CSC_ZStaRing_maxPt;
  std::pair<std::vector<float>,std::vector<float>> tnpPair_DT_dx_maxPt;
  std::pair<std::vector<float>,std::vector<float>> tnpPair_CSC_dx_maxPt;

  if (muons.isValid() && vtxs.isValid()) 
  {
    for (const auto & muon : (*muons)) 
      {
        bool isProbe = false;
        bool isTag = false;
	bool trigMatch = false;
	uint8_t DT_stationMatching = 0;
	uint8_t CSC_stationMatching = 0;

        std::vector<std::array<int,3>> probe_DT_WhSecSta;
        std::vector<std::array<int,3>> tag_DT_WhSecSta;
        std::vector<float> probe_DT_dx;
        std::vector<float> tag_DT_dx;

        std::vector<std::array<int,3>> probe_CSC_ZStaRing;
        std::vector<std::array<int,3>> tag_CSC_ZStaRing;
        std::vector<float> probe_CSC_dx;
        std::vector<float> tag_CSC_dx;

        //Getting trigger infos for tag selection
	if(triggerResults.isValid() && triggerEvent.isValid())
	{
	  const trigger::TriggerObjectCollection trigObjColl = triggerEvent->getObjects();
	  trigMatch = hasTrigger(m_trigIndices,trigObjColl,triggerEvent,muon);
	}

        
	//Probe selection
        if (m_probeSelector(muon) &&
	    (fabs(muon.muonBestTrack()->dxy(vertex.position())) < m_dxyCut) &&
	    (fabs(muon.muonBestTrack()->dz(vertex.position()))  < m_dzCut))
	{
	  isProbe = true;
          probe_muons.push_back(muon);
	  probe_muons_firesIsoTrig.push_back(trigMatch);
	}
        else
	{
	  probe_muons.push_back(dummy_muon);
	  probe_muons_DT_stationMatching.push_back(DT_stationMatching);
	  probe_muons_DT_WhSecSta.push_back(dummy_DT_WhSecSta);
          probe_muons_DT_dx.push_back(dummy_dx);

	  probe_muons_CSC_stationMatching.push_back(CSC_stationMatching);
	  probe_muons_CSC_ZStaRing.push_back(dummy_CSC_ZStaRing);
          probe_muons_CSC_dx.push_back(dummy_dx);

	  probe_muons_firesIsoTrig.push_back(false);
	}

        //Tag selection
        if (m_tagSelector(muon) &&
	    trigMatch)
	{
	  isTag = true;
          tag_muons.push_back(muon);
	}
	else
	{
	  tag_muons.push_back(dummy_muon);
	  tag_muons_DT_stationMatching.push_back(DT_stationMatching);
	  tag_muons_DT_WhSecSta.push_back(dummy_DT_WhSecSta);
          tag_muons_DT_dx.push_back(dummy_dx);

	  tag_muons_CSC_stationMatching.push_back(CSC_stationMatching);
	  tag_muons_CSC_ZStaRing.push_back(dummy_CSC_ZStaRing);
          tag_muons_CSC_dx.push_back(dummy_dx);
	}

	if (!isProbe && !isTag) continue;

        for (const auto chambMatch : muon.matches() ) 
          {
            // look in DTs
            if (chambMatch.detector() == MuonSubdetId::DT) 
	    {
              if (chambMatch.edgeX < m_borderCut && 
                  chambMatch.edgeY < m_borderCut)
                {
                  DTChamberId chId(chambMatch.id.rawId());

                  int wheel   = chId.wheel();
                  int sector  = chId.sector();
                  int station = chId.station();

	          reco::MuonSegmentMatch closest_matchedSegment;
                  double smallestDistance = 999.;
	          for (auto & seg : chambMatch.segmentMatches)
	          {
                    if ((chambMatch.x - seg.x) < smallestDistance){
	              smallestDistance = (chambMatch.x - seg.x);
	              closest_matchedSegment = seg;
	            }
	          }

	          DT_stationMatching = DT_stationMatching | (1 << (station-1));

	          if (isProbe){
                    probe_DT_WhSecSta.push_back({{wheel,sector,station}});
                    probe_DT_dx.push_back(std::abs(chambMatch.x - closest_matchedSegment.x));
	          }

	          if (isTag){
                    tag_DT_WhSecSta.push_back({{wheel,sector,station}});
                    tag_DT_dx.push_back(std::abs(chambMatch.x - closest_matchedSegment.x));
	          }

                }  
	      }

            // look in CSCs
            else if (chambMatch.detector() == MuonSubdetId::CSC)
	    {
              if (chambMatch.edgeX < m_borderCut && 
                  chambMatch.edgeY < m_borderCut)
                {
                  CSCDetId chId(chambMatch.id.rawId());

                  int zendcap = chId.zendcap();
                  int ring    = chId.ring();
                  int station = chId.station();

	          reco::MuonSegmentMatch closest_matchedSegment;
                  double smallestDistance = 999.;
	          for (auto & seg : chambMatch.segmentMatches)
	          {
                    if ((chambMatch.x - seg.x) < smallestDistance){
	              smallestDistance = (chambMatch.x - seg.x);
	              closest_matchedSegment = seg;
	            }
	          }

	          CSC_stationMatching = CSC_stationMatching | (1 << (station-1));

                
	          if (isProbe){
                    probe_CSC_ZStaRing.push_back({{zendcap,station,ring}});
                    probe_CSC_dx.push_back(std::abs(chambMatch.x - closest_matchedSegment.x));
	          }

	          if (isTag){
                    tag_CSC_ZStaRing.push_back({{zendcap,station,ring}});
                    tag_CSC_dx.push_back(std::abs(chambMatch.x - closest_matchedSegment.x));
	          }

              }  
	    }
	      else continue;

          }//loop over chamber matches

	if (isProbe){
	  probe_muons_DT_stationMatching.push_back(DT_stationMatching);
          probe_muons_DT_WhSecSta.push_back(probe_DT_WhSecSta);
          probe_muons_DT_dx.push_back(probe_DT_dx);

	  probe_muons_CSC_stationMatching.push_back(CSC_stationMatching);
          probe_muons_CSC_ZStaRing.push_back(probe_CSC_ZStaRing);
          probe_muons_CSC_dx.push_back(probe_CSC_dx);
	}
	if (isTag){
	  tag_muons_DT_stationMatching.push_back(DT_stationMatching);
          tag_muons_DT_WhSecSta.push_back(tag_DT_WhSecSta);
          tag_muons_DT_dx.push_back(tag_DT_dx);

	  tag_muons_CSC_stationMatching.push_back(CSC_stationMatching);
          tag_muons_CSC_ZStaRing.push_back(tag_CSC_ZStaRing);
          tag_muons_CSC_dx.push_back(tag_CSC_dx);
	}
      
      }//loop over muons
    }

    double pairPtMax = -999.;

    if (probe_muons.size()>0 && tag_muons.size()>0)
    {
      for (unsigned i=0; i<tag_muons.size(); ++i)
      {
        //Checking if tag muon is a dummy muon
        if (tag_muons.at(i).pt() < 0.001 && tag_muons.at(i).eta() < 0.001 && tag_muons.at(i).phi() < 0.001)
	  continue;

        for (unsigned j=0; j<probe_muons.size(); ++j)
        {
	  //Avoiding tag and probe to be the same object
          if (i==j) 
	    continue;

          //Checking if probe muon is a dummy muon
          if (probe_muons.at(j).pt() < 0.001 && probe_muons.at(j).eta() < 0.001 && probe_muons.at(j).phi() < 0.001)
	    continue;

	  int pair_charge_product = tag_muons.at(i).charge()*probe_muons.at(j).charge();
          math::PtEtaPhiMLorentzVector pairLorentzVector = tag_muons.at(i).polarP4() + probe_muons.at(j).polarP4();

	  //Fill the invariant mass plot for all pairs
	  if (pair_charge_product < 0) m_histos.find("pairMass")->second->Fill(pairLorentzVector.M());

          //If more than a tag+probe pair is found the highest pt one is selected
	  if (pair_charge_product < 0 && (pairLorentzVector.M() > m_lowPairMassCut && pairLorentzVector.M() < m_highPairMassCut) && pairLorentzVector.Pt()>pairPtMax)
	  {
            pairFound = true;

            pairPtMax = pairLorentzVector.Pt();
	    tnpPair_maxPt                     = std::make_pair(tag_muons.at(i),probe_muons.at(j));
            tnpPair_DT_stationMatching_maxPt  = std::make_pair(tag_muons_DT_stationMatching.at(i),probe_muons_DT_stationMatching.at(j));
            tnpPair_CSC_stationMatching_maxPt = std::make_pair(tag_muons_CSC_stationMatching.at(i),probe_muons_CSC_stationMatching.at(j));

            probeFiresTrig = probe_muons_firesIsoTrig.at(j);

	    tnpPair_DT_WhSecSta_maxPt    = std::make_pair(tag_muons_DT_WhSecSta.at(i),probe_muons_DT_WhSecSta.at(j));
	    tnpPair_DT_dx_maxPt          = std::make_pair(tag_muons_DT_dx.at(i),probe_muons_DT_dx.at(j));

	    tnpPair_CSC_ZStaRing_maxPt    = std::make_pair(tag_muons_CSC_ZStaRing.at(i),probe_muons_CSC_ZStaRing.at(j));
	    tnpPair_CSC_dx_maxPt          = std::make_pair(tag_muons_CSC_dx.at(i),probe_muons_CSC_dx.at(j));
	  }
        
        }
        
      }

    }

    if (pairFound)
    {
      //1) test the combination tag + probe
 
      //Fill detailed plots
      if (m_detailedAnalysis)
      {
        if(abs(tnpPair_maxPt.second.eta()) < 1.2){
          m_histos.find("probeEta")->second->Fill(tnpPair_maxPt.second.eta());
          m_histos.find("probePhi")->second->Fill(tnpPair_maxPt.second.phi());
          m_histos.find("probeNumberOfMatchedStations")->second->Fill(tnpPair_maxPt.second.numberOfMatchedStations());
        } 
        m_histos.find("probePt")->second->Fill(tnpPair_maxPt.second.pt());
      }
      
 
      //Fill DT numerator and denominator plots
      uint8_t DT_matchPatt = tnpPair_DT_stationMatching_maxPt.second;

      for(unsigned i = 0 ; i < tnpPair_DT_WhSecSta_maxPt.second.size(); ++i)
      {
        int wh  = tnpPair_DT_WhSecSta_maxPt.second.at(i)[0];
        int sec = tnpPair_DT_WhSecSta_maxPt.second.at(i)[1];
        int sta = tnpPair_DT_WhSecSta_maxPt.second.at(i)[2];

	float dx = tnpPair_DT_dx_maxPt.second.at(i);

        if ((DT_matchPatt & (1<<(sta-1))) != 0 && //avoids 0 station matching
            (DT_matchPatt & (1<<(sta-1))) !=DT_matchPatt && //avoids matching with the station under consideration only
            dx > 0.)
        {
          if (dx < 10.)
          {
            std::string hName = std::string("DT_nPassingProbePerCh_W") + std::to_string(wh);  
            m_histos.find(hName)->second->Fill(sec, sta);
            m_histos.find("DT_nPassingProbe_allCh")->second->Fill((sta) + 4*(wh+2));
          }
          else
          {
            std::string hName = std::string("DT_nFailingProbePerCh_W") + std::to_string(wh);  
            m_histos.find(hName)->second->Fill(sec, sta);
            m_histos.find("DT_nFailingProbe_allCh")->second->Fill((sta) + 4*(wh+2));
          }

        }

      }

      //Fill CSC numerator and denominator plots
      uint8_t CSC_matchPatt = tnpPair_CSC_stationMatching_maxPt.second;

      for(unsigned i = 0 ; i < tnpPair_CSC_ZStaRing_maxPt.second.size(); ++i)
      {
        int zendcap = tnpPair_CSC_ZStaRing_maxPt.second.at(i)[0];
        int sta     = tnpPair_CSC_ZStaRing_maxPt.second.at(i)[1];
        int ring    = tnpPair_CSC_ZStaRing_maxPt.second.at(i)[2];

	float dx = tnpPair_CSC_dx_maxPt.second.at(i);

        if ((CSC_matchPatt & (1<<(sta-1))) != 0 && //avoids 0 station matching
            (CSC_matchPatt & (1<<(sta-1))) !=CSC_matchPatt && //avoids matching with the station under consideration only
            dx > 0.)
        {
          if (dx < 10.)
          {
            m_histos.find("CSC_nPassingProbe_allCh")->second->Fill(zendcap*sta, ring);
          }
          else
          {
            m_histos.find("CSC_nFailingProbe_allCh")->second->Fill(zendcap*sta, ring);
          }

        }

      }

      //2) test the combination probe + tag
      DT_matchPatt = tnpPair_DT_stationMatching_maxPt.first;
      CSC_matchPatt = tnpPair_CSC_stationMatching_maxPt.first;

      if (m_tagSelector(tnpPair_maxPt.second) && probeFiresTrig) //consider probe as tag
      {
        //Fill detailed plots
        if (m_detailedAnalysis)
        {
          if(abs(tnpPair_maxPt.first.eta()) < 2.4){
            m_histos.find("probeEta")->second->Fill(tnpPair_maxPt.first.eta());
            m_histos.find("probePhi")->second->Fill(tnpPair_maxPt.first.phi());
            m_histos.find("probeNumberOfMatchedStations")->second->Fill(tnpPair_maxPt.first.numberOfMatchedStations());
          } 
          m_histos.find("probePt")->second->Fill(tnpPair_maxPt.first.pt());
        }

	//Fill DT numerator and denominator plots
        for(unsigned i = 0 ; i < tnpPair_DT_WhSecSta_maxPt.first.size(); ++i)
        {
          int wh  = tnpPair_DT_WhSecSta_maxPt.first.at(i)[0];
          int sec = tnpPair_DT_WhSecSta_maxPt.first.at(i)[1];
          int sta = tnpPair_DT_WhSecSta_maxPt.first.at(i)[2];

          float dx = tnpPair_DT_dx_maxPt.first.at(i);

          if ((DT_matchPatt & (1<<(sta-1))) != 0 && //avoids 0 station matching
              (DT_matchPatt & (1<<(sta-1))) !=DT_matchPatt && //avoids matching with the station under consideration only
              dx > 0.)
          {
            if (dx < 10.)
            {
              std::string hName = std::string("DT_nPassingProbePerCh_W") + std::to_string(wh);  
              m_histos.find(hName)->second->Fill(sec, sta);
              m_histos.find("DT_nPassingProbe_allCh")->second->Fill((sta) + 4*(wh+2));
            }
            else
            {
              std::string hName = std::string("DT_nFailingProbePerCh_W") + std::to_string(wh);  
              m_histos.find(hName)->second->Fill(sec, sta);
              m_histos.find("DT_nFailingProbe_allCh")->second->Fill((sta) + 4*(wh+2));
            }

          }

        } 

        //Fill CSC numerator and denominator plots
        for(unsigned i = 0 ; i < tnpPair_CSC_ZStaRing_maxPt.first.size(); ++i)
        {
          int zendcap = tnpPair_CSC_ZStaRing_maxPt.first.at(i)[0];
          int sta     = tnpPair_CSC_ZStaRing_maxPt.first.at(i)[1];
          int ring    = tnpPair_CSC_ZStaRing_maxPt.first.at(i)[2];

          float dx = tnpPair_CSC_dx_maxPt.first.at(i);

          if ((CSC_matchPatt & (1<<(sta-1))) != 0 && //avoids 0 station matching
              (CSC_matchPatt & (1<<(sta-1))) !=CSC_matchPatt && //avoids matching with the station under consideration only
              dx > 0.)
          {
            if (dx < 10.)
            {
              m_histos.find("CSC_nPassingProbe_allCh")->second->Fill(zendcap*sta, ring);
            }
            else
            {
              m_histos.find("CSC_nFailingProbe_allCh")->second->Fill(zendcap*sta, ring);
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

  MonitorElement* me_DT_pass = iBooker.book2D(hName_DT_pass.c_str(), hName_DT_pass.c_str(), 14, 0.5, 14.5, 4, 0., 4.5);
  MonitorElement* me_DT_fail = iBooker.book2D(hName_DT_fail.c_str(), hName_DT_fail.c_str(), 14, 0.5, 14.5, 4, 0., 4.5);

  me_DT_pass->setBinLabel(1, "MB1", 2);
  me_DT_pass->setBinLabel(2, "MB2", 2);
  me_DT_pass->setBinLabel(3, "MB3", 2);
  me_DT_pass->setBinLabel(4, "MB4", 2);
  me_DT_pass->setAxisTitle("Sector", 1);

  me_DT_fail->setBinLabel(1, "MB1", 2);
  me_DT_fail->setBinLabel(2, "MB2", 2);
  me_DT_fail->setBinLabel(3, "MB3", 2);
  me_DT_fail->setBinLabel(4, "MB4", 2);
  me_DT_fail->setAxisTitle("Sector", 1);

  m_histos[hName_DT_pass] = me_DT_pass;
  m_histos[hName_DT_fail] = me_DT_fail;
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
