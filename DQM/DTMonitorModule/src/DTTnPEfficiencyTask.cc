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

template <class T, size_t WH, size_t SEC, size_t STA> using WheelSectorStationMatrix = std::array<std::array<std::array<T, STA>, SEC>, WH>;

DTTnPEfficiencyTask::DTTnPEfficiencyTask(const edm::ParameterSet& config) : 
  m_nEvents(0),
  m_muToken(consumes<reco::MuonCollection>(config.getUntrackedParameter<edm::InputTag>("inputTagMuons"))),
  m_primaryVerticesToken(consumes<std::vector<reco::Vertex>>(config.getUntrackedParameter<edm::InputTag>("inputTagPrimaryVertices"))),
  m_triggerResultsToken(consumes<edm::TriggerResults>(config.getUntrackedParameter<edm::InputTag>("trigResultsTag"))),
  m_triggerEventToken(consumes<trigger::TriggerEvent>(config.getUntrackedParameter<edm::InputTag>("trigEventTag"))),
  m_trigName(config.getUntrackedParameter<std::string>("trigName")),
  m_isoTrigName(config.getUntrackedParameter<std::string>("isoTrigName")),
  m_detailedAnalysis(config.getUntrackedParameter<bool>("detailedAnalysis")),
  m_selector(config.getUntrackedParameter<std::string>("probeCut")),
  m_dxyCut(config.getUntrackedParameter<double>("probeDxyCut")),
  m_dzCut(config.getUntrackedParameter<double>("probeDzCut")),
  m_tagSelector(config.getUntrackedParameter<std::string>("tagCut")),
  m_borderCut(-10.)
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
    //std::cout<<"--------"<<pathName<<"----------"<<std::endl;
    if(pathName.Contains(tNamePattern))
      m_trigIndices.push_back(static_cast<int>(iPath));
  }

  tName = TString(m_isoTrigName);
  tNamePattern = TRegexp(tName,enableWildCard);

  for(unsigned iPath=0; iPath<m_hltConfig.size();++iPath)
  {
    TString pathName = TString(m_hltConfig.triggerName(iPath));
    //std::cout<<"--------"<<pathName<<"----------"<<std::endl;
    if(pathName.Contains(tNamePattern))
      m_isoTrigIndices.push_back(static_cast<int>(iPath));
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
      m_histos["probeEta"] = iBooker.book1D("probeEta", "probeEta",18, -1.2, 1.2);
      m_histos["probePhi"] = iBooker.book1D("probePhi", "probePhi",30, -TMath::Pi(), TMath::Pi());
      m_histos["probeNumberOfMatchedStations"] = iBooker.book1D("probeNumberOfMatchedStations", "probeNumberOfMatchedStations",5, 0., 5);
      m_histos["pairMass"] = iBooker.book1D("pairMass", "pairMass", 25, 70., 120.);
    }
      
  for (int wheel = -2; wheel <= 2; ++wheel) 
    {
      bookWheelHistos(iBooker, wheel, "Task");
    }

}

void DTTnPEfficiencyTask::analyze(const edm::Event& event, const edm::EventSetup& context) 
{
  ++m_nEvents;

  bool pairFound = true;

  edm::Handle<reco::MuonCollection> muons;
  event.getByToken(m_muToken, muons);

  edm::Handle<std::vector<reco::Vertex>> vtxs;
  event.getByToken(m_primaryVerticesToken, vtxs);
  const reco::Vertex & vertex = vtxs->at(0);

  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByToken(m_triggerResultsToken, triggerResults);

  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByToken(m_triggerEventToken, triggerEvent);

  std::vector<reco::Muon> probe_muons;
  std::vector<uint8_t> probe_muons_stationMatching;
  std::vector<WheelSectorStationMatrix <float, 5, 14, 4>> probe_muons_dx;
  std::vector<bool> probe_muons_firesTrig;
  std::vector<bool> probe_muons_firesIsoTrig;

  WheelSectorStationMatrix <float, 5, 14, 4> DummyWhSecSta_dx;
  for (auto& wh: DummyWhSecSta_dx)
  {
    for (auto& sec: wh)
    {
      for(auto& sta: sec)
      {
        sta = -999.;
      }
    }
  }
  
  std::map<std::array<unsigned,2>,std::vector<std::pair<reco::Muon,reco::Muon>>> tnp_pairs;
  std::map<std::array<unsigned,2>,std::vector<std::pair<uint8_t,uint8_t>>> tnp_pairs_stationMatching;
  std::map<std::array<unsigned,2>,std::vector<std::pair<WheelSectorStationMatrix<float, 5, 14, 4>,WheelSectorStationMatrix<float, 5, 14, 4>>>> tnp_pairs_dx;

  //Selected tnp pair 
  std::pair<reco::Muon,reco::Muon> tnpPair;
  std::pair<uint8_t,uint8_t> tnpPairStationMatching;
  std::pair<WheelSectorStationMatrix<float, 5, 14, 4>,WheelSectorStationMatrix<float, 5, 14, 4>> tnpPairHitSegDx;

  if (muons.isValid() && vtxs.isValid()) 
  {
    for (const auto & muon : (*muons)) 
      {
	uint8_t stationMatching = 0;
        WheelSectorStationMatrix <float, 5, 14, 4> WhSecSta_dx = DummyWhSecSta_dx;
        
        if (!m_selector(muon) ||
	    !(fabs(muon.muonBestTrack()->dxy(vertex.position())) < m_dxyCut) ||
	    !(fabs(muon.muonBestTrack()->dz(vertex.position()))  < m_dzCut))
          continue;

        probe_muons.push_back(muon);

	if(triggerResults.isValid() && triggerEvent.isValid())
	{
	  const trigger::TriggerObjectCollection trigObjColl = triggerEvent->getObjects();

	  probe_muons_firesTrig.push_back(hasTrigger(m_trigIndices,trigObjColl,triggerEvent,muon));
	  probe_muons_firesIsoTrig.push_back(hasTrigger(m_isoTrigIndices,trigObjColl,triggerEvent,muon));
	}
	else
	{
	  probe_muons_firesTrig.push_back(false);
	  probe_muons_firesIsoTrig.push_back(false);
	  
	}

        for (const auto chambMatch : muon.matches() ) 
          {
            // look only in DTs
            if (chambMatch.detector() != MuonSubdetId::DT)
              continue;
 
            if (chambMatch.edgeX < m_borderCut && 
                chambMatch.edgeY < m_borderCut)
              {
                DTChamberId chId(chambMatch.id.rawId());

                int wheel   = chId.wheel();
                int sector  = chId.sector();
                int station = chId.station();

		reco::MuonSegmentMatch matchedSegment;
		for (auto & seg : chambMatch.segmentMatches)
		{
		  if (seg.isMask())
		  {
		    matchedSegment = seg;
		    break;
		  }
		}

		stationMatching = stationMatching | (1 << (station-1));

		//WhSecSta_dx is filled reordering wheel labels
		// -2 ---> 0
		// -1 ---> 1
		//  0 ---> 2
		// +1 ---> 3
		// +2 ---> 4
		WhSecSta_dx[wheel+2][sector-1][station-1] = std::abs(chambMatch.x - matchedSegment.x);

              }  
          }//loop over chamber matches

	probe_muons_stationMatching.push_back(stationMatching);
	probe_muons_dx.push_back(WhSecSta_dx);
      
      }//loop over muons
    }


    //Selecting only oppositely charged muons pairs
    //with an invariant mass between 90 +/- 10 GeV
    if (probe_muons.size())
    {
      for(unsigned i=0; i<probe_muons.size(); i++)
      {
        if(!m_tagSelector(probe_muons.at(i)) //&&
	   /*!probe_muons_firesIsoTrig.at(i)*/) 
	  continue;
	for(unsigned j=0; j<probe_muons.size(); j++)
	{
          if (i==j) continue;
          std::array<unsigned,2> pairIndex {{i,j}};
	  std::sort(pairIndex.begin(),pairIndex.end());

	  int pair_charge_product = probe_muons.at(i).charge()*probe_muons.at(j).charge();
          math::PtEtaPhiMLorentzVector pairLorentzVector = probe_muons.at(i).polarP4() + probe_muons.at(j).polarP4();
	  if (pair_charge_product < 0 && (pairLorentzVector.M() > 80. && pairLorentzVector.M() < 100.))
	  {
	    std::pair tnp_pair = std::make_pair(probe_muons.at(i),probe_muons.at(j));
            std::pair tnp_pair_stationMatching = std::make_pair(probe_muons_stationMatching.at(i),probe_muons_stationMatching.at(j));
	    std::pair tnp_pair_dx = std::make_pair(probe_muons_dx.at(i),probe_muons_dx.at(j));

	    tnp_pairs[pairIndex].push_back(tnp_pair);
	    tnp_pairs_stationMatching[pairIndex].push_back(tnp_pair_stationMatching);
	    tnp_pairs_dx[pairIndex].push_back(tnp_pair_dx);
	  }
	}
      }
    }


    //If more than one tag+probe pair is
    //found, the highest pt pair is chosen
    std::array<unsigned,2> maxPtPairIdx;
    if(tnp_pairs.size()>0)
    {
      math::PtEtaPhiMLorentzVector maxPtTnpPairLorentzVector {};
      for (auto & p : tnp_pairs)
      {
        std::array<unsigned,2> pairIdx = p.first;
        math::PtEtaPhiMLorentzVector tagLorentzVector (probe_muons[pairIdx[0]].polarP4());
        math::PtEtaPhiMLorentzVector probeLorentzVector (probe_muons[pairIdx[1]].polarP4());
        math::PtEtaPhiMLorentzVector pair = tagLorentzVector + probeLorentzVector;
        if (pair.Pt() > maxPtTnpPairLorentzVector.Pt())
        {
	  maxPtPairIdx = pairIdx;
        }
      }
    }
    else
      pairFound = false;

    //If at least one tag + probe pair
    //is found then fill the histograms
    if (pairFound)
    {

      if (m_detailedAnalysis)
      {
        math::PtEtaPhiMLorentzVector tag (tnp_pairs[maxPtPairIdx][0].first.polarP4());
        math::PtEtaPhiMLorentzVector probe (tnp_pairs[maxPtPairIdx][0].second.polarP4());
        math::PtEtaPhiMLorentzVector TnpPair = tag + probe;
        m_histos.find("pairMass")->second->Fill(TnpPair.M());

	if(abs(tnp_pairs[maxPtPairIdx][0].first.eta()) < 1.2){
          m_histos.find("probeEta")->second->Fill(tnp_pairs[maxPtPairIdx][0].first.eta());
          m_histos.find("probePhi")->second->Fill(tnp_pairs[maxPtPairIdx][0].first.phi());
          m_histos.find("probeNumberOfMatchedStations")->second->Fill(tnp_pairs[maxPtPairIdx][0].first.numberOfMatchedStations());
	}

	if(abs(tnp_pairs[maxPtPairIdx][0].second.eta()) < 1.2){
          m_histos.find("probeEta")->second->Fill(tnp_pairs[maxPtPairIdx][0].second.eta());
          m_histos.find("probePhi")->second->Fill(tnp_pairs[maxPtPairIdx][0].second.phi());
          m_histos.find("probeNumberOfMatchedStations")->second->Fill(tnp_pairs[maxPtPairIdx][0].second.numberOfMatchedStations());
	}

        m_histos.find("probePt")->second->Fill(tnp_pairs[maxPtPairIdx][0].second.pt());
        m_histos.find("probePt")->second->Fill(tnp_pairs[maxPtPairIdx][0].first.pt());
      }

      for (unsigned wheelIdx=0; wheelIdx < 5; ++wheelIdx)
      {
        for (unsigned sectorIdx=0; sectorIdx < 14; ++sectorIdx)
        {
          for(unsigned stationIdx=0; stationIdx < 4; ++stationIdx)
          {
	    //Loop over the two possible combination (tag+probe, probe+tag)
	    //for the highest pt pair
            for(unsigned iPair=0; iPair<tnp_pairs[maxPtPairIdx].size(); ++iPair)
            {
	      uint8_t matchPatt = tnp_pairs_stationMatching[maxPtPairIdx][iPair].second;
              WheelSectorStationMatrix <float, 5, 14, 4> dx = tnp_pairs_dx[maxPtPairIdx][iPair].second;

              if ((matchPatt & (1<<stationIdx)) != 0 &&
                  (matchPatt & (1<<stationIdx)) !=matchPatt &&
                  dx[wheelIdx][sectorIdx][stationIdx] > 0.)
              {
                if (dx[wheelIdx][sectorIdx][stationIdx] < 10.)
                {
      	          std::string hName = std::string("nPassingProbePerCh_W") + std::to_string((int)wheelIdx-2);  
      	          m_histos.find(hName)->second->Fill(sectorIdx+1, stationIdx+1);
      	          m_histos.find("nPassingProbe_allCh")->second->Fill((stationIdx+1) + 4*(wheelIdx));
                }
                else
                {
      	          std::string hName = std::string("nFailingProbePerCh_W") + std::to_string((int)wheelIdx-2);  
      	          m_histos.find(hName)->second->Fill(sectorIdx+1, stationIdx+1);
      	          m_histos.find("nFailingProbe_allCh")->second->Fill((stationIdx+1) + 4*(wheelIdx));
                }
              }
            }

          }//loop over stations
        }//loop over sectors
      }//loop over wheels
 
    }
}

void DTTnPEfficiencyTask::bookWheelHistos(DQMStore::IBooker& iBooker, 
					  int wheel, std::string folder) 
{

  auto baseDir = topFolder() + folder + "/";
  iBooker.setCurrentFolder(baseDir);

  LogTrace("DTDQM|DTMonitorModule|DTTnPEfficiencyTask")
    << "[DTTnPEfficiencyTask]: booking histos in " << baseDir << std::endl;

  auto hName_pass = std::string("nPassingProbePerCh_W") + std::to_string(wheel);    
  auto hName_fail = std::string("nFailingProbePerCh_W") + std::to_string(wheel);    

  MonitorElement* me_pass = iBooker.book2D(hName_pass.c_str(), hName_pass.c_str(), 14, 0.5, 14.5, 4, 0., 4.5);
  MonitorElement* me_fail = iBooker.book2D(hName_fail.c_str(), hName_fail.c_str(), 14, 0.5, 14.5, 4, 0., 4.5);

  MonitorElement* me_pass_allCh = iBooker.book1D("nPassingProbe_allCh", "nPassingProbe_allCh", 20, 0.5, 20.5);
  MonitorElement* me_fail_allCh = iBooker.book1D("nFailingProbe_allCh", "nFailingProbe_allCh", 20, 0.5, 20.5);

  me_pass->setBinLabel(1, "MB1", 2);
  me_pass->setBinLabel(2, "MB2", 2);
  me_pass->setBinLabel(3, "MB3", 2);
  me_pass->setBinLabel(4, "MB4", 2);
  me_pass->setAxisTitle("Sector", 1);

  me_fail->setBinLabel(1, "MB1", 2);
  me_fail->setBinLabel(2, "MB2", 2);
  me_fail->setBinLabel(3, "MB3", 2);
  me_fail->setBinLabel(4, "MB4", 2);
  me_fail->setAxisTitle("Sector", 1);

  me_pass_allCh->setBinLabel(1 , "MB1/YB-2", 1);
  me_pass_allCh->setBinLabel(2 , "MB2/YB-2", 1);
  me_pass_allCh->setBinLabel(3 , "MB3/YB-2", 1);
  me_pass_allCh->setBinLabel(4 , "MB4/YB-2", 1);
  me_pass_allCh->setBinLabel(5 , "MB1/YB-1", 1);
  me_pass_allCh->setBinLabel(6 , "MB2/YB-1", 1);
  me_pass_allCh->setBinLabel(7 , "MB3/YB-1", 1);
  me_pass_allCh->setBinLabel(8 , "MB4/YB-1", 1);
  me_pass_allCh->setBinLabel(9 , "MB1/YB0", 1);
  me_pass_allCh->setBinLabel(10, "MB2/YB0", 1);
  me_pass_allCh->setBinLabel(11, "MB3/YB0", 1);
  me_pass_allCh->setBinLabel(12, "MB4/YB0", 1);
  me_pass_allCh->setBinLabel(13, "MB1/YB1", 1);
  me_pass_allCh->setBinLabel(14, "MB2/YB1", 1);
  me_pass_allCh->setBinLabel(15, "MB3/YB1", 1);
  me_pass_allCh->setBinLabel(16, "MB4/YB1", 1);
  me_pass_allCh->setBinLabel(17, "MB1/YB2", 1);
  me_pass_allCh->setBinLabel(18, "MB2/YB2", 1);
  me_pass_allCh->setBinLabel(19, "MB3/YB2", 1);
  me_pass_allCh->setBinLabel(20, "MB4/YB2", 1);
  me_pass_allCh->setAxisTitle("Number of passing probes", 2);

  me_fail_allCh->setBinLabel(1 , "MB1/YB-2", 1);
  me_fail_allCh->setBinLabel(2 , "MB2/YB-2", 1);
  me_fail_allCh->setBinLabel(3 , "MB3/YB-2", 1);
  me_fail_allCh->setBinLabel(4 , "MB4/YB-2", 1);
  me_fail_allCh->setBinLabel(5 , "MB1/YB-1", 1);
  me_fail_allCh->setBinLabel(6 , "MB2/YB-1", 1);
  me_fail_allCh->setBinLabel(7 , "MB3/YB-1", 1);
  me_fail_allCh->setBinLabel(8 , "MB4/YB-1", 1);
  me_fail_allCh->setBinLabel(9 , "MB1/YB0", 1);
  me_fail_allCh->setBinLabel(10, "MB2/YB0", 1);
  me_fail_allCh->setBinLabel(11, "MB3/YB0", 1);
  me_fail_allCh->setBinLabel(12, "MB4/YB0", 1);
  me_fail_allCh->setBinLabel(13, "MB1/YB1", 1);
  me_fail_allCh->setBinLabel(14, "MB2/YB1", 1);
  me_fail_allCh->setBinLabel(15, "MB3/YB1", 1);
  me_fail_allCh->setBinLabel(16, "MB4/YB1", 1);
  me_fail_allCh->setBinLabel(17, "MB1/YB2", 1);
  me_fail_allCh->setBinLabel(18, "MB2/YB2", 1);
  me_fail_allCh->setBinLabel(19, "MB3/YB2", 1);
  me_fail_allCh->setBinLabel(20, "MB4/YB2", 1);
  me_fail_allCh->setAxisTitle("Number of failing probes", 2);

  m_histos[hName_pass] = me_pass;
  m_histos[hName_fail] = me_fail;

  m_histos["nPassingProbe_allCh"] = me_pass_allCh;
  m_histos["nFailingProbe_allCh"] = me_fail_allCh;
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

    for (unsigned iTrig=0; nTriggers; ++iTrig)
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
