#ifndef DTTnPEfficiencyTask_H
#define DTTnPEfficiencyTask_H

/*
 * \file DTTnPEfficiencyTask.h
 *
 * \author L. Lunerti - INFN Bologna
 *
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/AnySelector.h"
 
#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include <vector>
#include <string>
#include <map>

class DTTnPEfficiencyTask : public DQMEDAnalyzer
{

public:
  /// Constructor
  DTTnPEfficiencyTask(const edm::ParameterSet& config);

  /// Destructor
  ~DTTnPEfficiencyTask() override;

protected:

  /// BeginRun
  void dqmBeginRun(const edm::Run& run, const edm::EventSetup& context) override;

  void bookHistograms(DQMStore::IBooker& iBooker, edm::Run const& run, edm::EventSetup const& context) override;

  /// Book wheel granularity histograms
  void bookWheelHistos(DQMStore::IBooker& iBooker, int wheel, std::string folder = "");

  /// Book endcap histograms
  void bookEndcapHistos(DQMStore::IBooker& iBooker, int stations, std::string folder = "");

  /// Return the top folder
  inline std::string topFolder() const { return "DT/10-Segment_TnP/"; };

  /// Analyze
  void analyze(const edm::Event& event, const edm::EventSetup& context) override;
  bool hasTrigger(std::vector<int> & trigIndices,const trigger::TriggerObjectCollection & trigObjs,edm::Handle<trigger::TriggerEvent> & trigEvent,const reco::Muon & muon);

  int get_barrel_histo_ycoord(int ring, int station, int sector, int layer, int subsector, int roll);

  /// To reset the MEs

private:

  int m_nEvents;

  edm::EDGetTokenT<reco::MuonCollection> m_muToken;
  edm::EDGetTokenT<std::vector<reco::Vertex>> m_primaryVerticesToken;
  edm::EDGetTokenT<edm::TriggerResults> m_triggerResultsToken;
  edm::EDGetTokenT<trigger::TriggerEvent> m_triggerEventToken;

  std::string m_trigName;
  HLTConfigProvider m_hltConfig;

  bool m_detailedAnalysis;

  //Probe selectors
  StringCutObjectSelector<reco::Candidate,true> m_probeSelector;
  double m_dxyCut;
  double m_dzCut;

  //Tag selectors
  StringCutObjectSelector<reco::Muon,true> m_tagSelector;

  //Trigger indices
  std::vector<int> m_trigIndices;

  std::map<std::string, MonitorElement*> m_histos;

  double m_borderCut;
  double m_lowPairMassCut;
  double m_highPairMassCut;
  double m_dxCut;

};

#endif
