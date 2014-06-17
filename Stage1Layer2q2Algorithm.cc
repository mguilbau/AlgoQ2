///
/// \class l1t::Stage1Layer2q2Algorithm
///
/// \authors: Gian Michele Innocenti
///           R. Alex Barbieri
///
/// Description: q2 Algorithm HI

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2EtSumAlgorithmImp.h"
#include "L1Trigger/L1TCalorimeter/interface/PUSubtractionMethods.h"
#include "L1Trigger/L1TCalorimeter/interface/legacyGtHelper.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

l1t::Stage1Layer2q2Algorithm::Stage1Layer2q2Algorithm(CaloParamsStage1* params) : params_(params)
{

}


l1t::Stage1Layer2q2Algorithm::~Stage1Layer2q2Algorithm() {


}


void l1t::Stage1Layer2q2Algorithm::processEvent(const std::vector<l1t::CaloRegion> & regions,
							const std::vector<l1t::CaloEmCand> & EMCands,
							      std::vector<l1t::EtSum> * etsums) {

  double regionET=0.;  
  
  int counterregion=0;
  
  double sumET[L1CaloRegionDetId::N_PHI];
    
    for(std::vector<CaloRegion>::const_iterator region = regions.begin(); region != regions.end(); region++) {
   
    int ieta=region->hwEta();    
    if (ieta > 3 && ieta < 18) {
      continue;
    }

    int iphi=region->hwPhi();    
    regionET=region->hwPt();
    
    sumET[iphi] += regionET;
    counterregion++;
    
  }

  std::vector<double> q2(2,0.);
  double sumW=0.;
  for (unsigned int iphi=0; iphi<L1CaloRegionDetId::N_PHI; iphi++) {
// cos and sin must be come from in a LUT at the end...
      q2.at(0)+=sumET[iphi] * cos(2. * 3.1415927 * iphi * 1.0 / L1CaloRegionDetId::N_PHI);
      q2.at(1)+=sumET[iphi] * sin(2. * 3.1415927 * iphi * 1.0 / L1CaloRegionDetId::N_PHI);
  }
  
  double HFq2 = q2.at(0)*q2.at(0)+q2.at(1)*q2.at(1);
//  double psi2 = 0.5 * atan(q2y/q2x);
 
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > etLorentz(0,0,0,0);

  // convert back to hardware ET
  l1t::EtSum etTot (*&etLorentz,EtSum::EtSumType::kTotalEt,HFq2,0,0,0);

  std::vector<l1t::EtSum> *preGtEtSums = new std::vector<l1t::EtSum>();
  preGtEtSums->push_back(etTot);
        
  EtSumToGtScales(params_, preGtEtSums, etsums);

  delete preGtEtSums;

  // ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > etLorentz(0,0,0,0);

  // // convert back to hardware ET
  // l1t::EtSum etMiss(*&etLorentz,EtSum::EtSumType::kMissingEt,MET/jetLsb ,0,iPhiET,0);
  // l1t::EtSum htMiss(*&etLorentz,EtSum::EtSumType::kMissingHt,MHT/jetLsb ,0,iPhiHT,0);
  // l1t::EtSum etTot (*&etLorentz,EtSum::EtSumType::kTotalEt,sumET/jetLsb,0,0,0);
  // l1t::EtSum htTot (*&etLorentz,EtSum::EtSumType::kTotalHt,sumHT/jetLsb ,0,0,0);

  // std::vector<l1t::EtSum> *preGtEtSums = new std::vector<l1t::EtSum>();

  // preGtEtSums->push_back(etMiss);
  // preGtEtSums->push_back(htMiss);
  // preGtEtSums->push_back(etTot);
  // preGtEtSums->push_back(htTot);

  // // All algorithms
  // EtSumToGtScales(params_, preGtEtSums, etsums);

  // delete subRegions;
  // delete preGtEtSums;

}
