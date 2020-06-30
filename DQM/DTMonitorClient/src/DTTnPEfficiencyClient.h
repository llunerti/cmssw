#ifndef DTTnPEfficiencyClient_H
#define DTTnPEfficiencyClient_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/DQMEDHarvester.h"

class DTTnPEfficiencyClient : public DQMEDHarvester {
public:
  DTTnPEfficiencyClient (const edm::ParameterSet & pSet);
  ~DTTnPEfficiencyClient () override;
protected:
  void beginRun(const edm::Run &, const edm::EventSetup &) override;
  void dqmEndJob(DQMStore::IBooker &, DQMStore::IGetter &) override;

  /// book the report summary

  void bookHistos(DQMStore::IBooker &, std::string folder);
  void dqmEndLuminosityBlock(DQMStore::IBooker &,
                             DQMStore::IGetter &,
                             edm::LuminosityBlock const &,
                             edm::EventSetup const &) override;

  /// Return the top folder
  inline std::string topFolder() const { return "DT/10-Segment_TnP/"; };
private:
  std::map<std::string,MonitorElement*> effHistos;
};

#endif
