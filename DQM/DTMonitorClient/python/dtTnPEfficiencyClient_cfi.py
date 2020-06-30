import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

dtTnPEfficiencyClient = DQMEDHarvester("DTTnPEfficiencyClient",
                                       diagnosticPrescale = cms.untracked.int32(1))
