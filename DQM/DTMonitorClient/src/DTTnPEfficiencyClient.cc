#include "DQM/DTMonitorClient/src/DTTnPEfficiencyClient.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

DTTnPEfficiencyClient::DTTnPEfficiencyClient (const edm::ParameterSet &pSet) :
  passNfailHistoNames(pSet.getUntrackedParameter<std::vector<std::string>>("histoNames"))
{
  edm::LogVerbatim("DQM|DTMonitorModule|DTTnPEfficiencyClient") << "DTTnPEfficiencyClient: Constructor called";
};

DTTnPEfficiencyClient::~DTTnPEfficiencyClient () {
  edm::LogVerbatim("DQM|DTMonitorModule|DTTnPEfficiencyClient") << "DTTnPEfficiencyClient: Constructor called";
};

void DTTnPEfficiencyClient::beginRun(const edm::Run& run, const edm::EventSetup& setup) {
}

void DTTnPEfficiencyClient::dqmEndLuminosityBlock(DQMStore::IBooker& ibooker,
                                                      DQMStore::IGetter& igetter,
                                                      edm::LuminosityBlock const& lumiSeg,
                                                      edm::EventSetup const& setup) {
  edm::LogVerbatim("DTDQM|DTMonitorClient|DTChamberEfficiencyClient") << "DTChamberEfficiencyClient: endluminosityBlock";
}

void DTTnPEfficiencyClient::dqmEndJob(DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter) {
  edm::LogVerbatim("DTDQM|DTMonitorClient|DTChamberEfficiencyClient") << "DTChamberEfficiencyClient: endRun";

  std::string outFolder = "Task";

  //bookHistos(ibooker,outFolder);
  
  ibooker.setCurrentFolder(topFolder() + outFolder + "/");
  std::string baseFolder = topFolder() + outFolder + "/";

  TH1::SetDefaultSumw2(kTRUE);

  for (auto s : passNfailHistoNames)
  {
    TH1::SetDefaultSumw2(kTRUE);
    
    std::string passHistoName = s.substr(0,s.find(":"));
    std::string failHistoName = s.substr(s.find(":")+1,s.length());

    std::cout<<"*******PASS HISTO NAME: "<<passHistoName<<"********"<<std::endl;
    std::cout<<"*******FAIL HISTO NAME: "<<failHistoName<<"********"<<std::endl;

    std::string histoName_pass = baseFolder + passHistoName;
    std::string histoName_fail = baseFolder + failHistoName;
    
    MonitorElement* me_pass = igetter.get(histoName_pass);
    MonitorElement* me_fail = igetter.get(histoName_fail);

    int fir = passHistoName.find("_");
    int sec = passHistoName.find("_",fir+1);
    std::string chName    = passHistoName.substr(0,fir);
    std::string specifier = passHistoName.substr(sec+1);
    std::string effHistoName = chName + "_chamberEff_" + specifier;

    if (!me_pass || !me_fail ) {
      edm::LogWarning("DTTnPEfficiencyClient") << "Monitor Element not available" << std::endl;
      return;
    }

    //1D histos
    if((int)me_pass->kind()==DQMNet::DQM_PROP_TYPE_TH1F && (int)me_fail->kind()==DQMNet::DQM_PROP_TYPE_TH1F){

      if (!(me_pass->getTH1F()) || !(me_fail->getTH1F())) {
        edm::LogWarning("DTTnPEfficiencyClient") << "Monitor Element not available" << std::endl;
        return;
      }

      TH1F* h1_pass = me_pass->getTH1F();
      TH1F* h1_fail = me_fail->getTH1F();

      const int nBinX_pass = h1_pass->GetNbinsX();
      const int nBinX_fail = h1_fail->GetNbinsX();

      if (nBinX_pass!=nBinX_fail){
        edm::LogWarning("DTTnPEfficiencyClient") << "Histograms with different number of bins: unable to compute the ratio" << std::endl;
        return;
      }

      TH1F* h1_den = (TH1F*)h1_pass->Clone();
      TH1F* h1_num = (TH1F*)h1_pass->Clone();
      h1_den->Sumw2();
      h1_num->Sumw2();
      h1_den->Add(h1_fail);

      h1_num->Divide(h1_den);
      TH1F* h1_ratio = (TH1F*)h1_num->Clone();
      
      effHistos[effHistoName] = ibooker.book1D(effHistoName,h1_ratio);
      effHistos[effHistoName]->setTitle(effHistoName);
      effHistos[effHistoName]->setAxisTitle("Efficiency", 2);
    }

    //2D histos
    if((int)me_pass->kind()==DQMNet::DQM_PROP_TYPE_TH2F && (int)me_fail->kind()==DQMNet::DQM_PROP_TYPE_TH2F){

      if (!(me_pass->getTH2F()) || !(me_fail->getTH2F())) {
        edm::LogWarning("DTTnPEfficiencyClient") << "Monitor Element not available: unable to compute the ratio" << std::endl;
        return;
      }

      TH2F* h2_pass = me_pass->getTH2F();
      TH2F* h2_fail = me_fail->getTH2F();

      const int nBinX_pass = h2_pass->GetNbinsX();
      const int nBinX_fail = h2_fail->GetNbinsX();
      const int nBinY_pass = h2_pass->GetNbinsY();
      const int nBinY_fail = h2_fail->GetNbinsY();

      if ((nBinX_pass!=nBinX_fail) || (nBinY_pass!=nBinY_fail)){
        edm::LogWarning("DTTnPEfficiencyClient") << "Histograms with different number of bins: unable to compute the ratio" << std::endl;
        return;
      }

      TH2F* h2_den = (TH2F*)h2_pass->Clone();
      TH2F* h2_num = (TH2F*)h2_pass->Clone();
      h2_den->Sumw2();
      h2_num->Sumw2();
      h2_den->Add(h2_fail);

      h2_num->Divide(h2_den);
      TH2F* h2_ratio = (TH2F*)h2_num->Clone();

      effHistos[effHistoName] = ibooker.book2D(effHistoName,h2_ratio);
      effHistos[effHistoName]->setTitle(effHistoName);
      effHistos[effHistoName]->setAxisTitle("Efficiency", 3);
    }
    
  }

  return;
}

void DTTnPEfficiencyClient::bookHistos(DQMStore::IBooker& ibooker, std::string folder) {
  ibooker.setCurrentFolder(topFolder() + folder + "/");
  
  effHistos["DT_chamberEff_allCh"] = ibooker.book1D("DT_chamberEff_allCh","DT_chamberEff_allCh", 20, 1., 15.);
  effHistos["CSC_chamberEff_allCh"] = ibooker.book2D("CSC_chamberEff_allCh","CSC_chamberEff_allCh", 9, -4., 5., 4, 0.,4.5);

  for (int wh = -2; wh <= 2; wh++) {
    std::string histoName = std::string("DT_chamberEff_W") + std::to_string(wh);

    effHistos[histoName] = ibooker.book2D(histoName.c_str(),histoName.c_str(), 14, 1., 15., 4, 0., 5.);
  }

  return;
}
