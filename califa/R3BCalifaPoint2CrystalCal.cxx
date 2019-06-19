// ------------------------------------------------------------
// -----      R3BPoint2CrystalCal header file           		-----
// -----      Created 1/03/2018 by E.Galiana          		-----
// -----                                                  -----
// ------------------------------------------------------------

#include "R3BCalifaPoint2CrystalCal.h"
#include "R3BCalifa.h"
#include "TClonesArray.h"
#include "FairRootManager.h"
#include "TRandom.h"
#include "TArrayD.h"
#include "TVector3.h"
#include "TMath.h"
#include <iostream>
#include <stdlib.h>
#include "R3BCalifaPoint.h"
#include "R3BCalifaCrystalCalData.h"

using std::cout;
using std::cerr;
using std::endl;


//R3BCalifaPoint2CrystalCal: Default Constructor
R3BCalifaPoint2CrystalCal::R3BCalifaPoint2CrystalCal() : FairTask("R3B CALIFA Point2CrystalCal") {
}

//R3BCalifaPoint2CrystalCal: Constructor
R3BCalifaPoint2CrystalCal::R3BCalifaPoint2CrystalCal(const TString& geoFile) : FairTask("R3B CALIFA Point2CrystalCal"),
								 fCalifaPointDataCA(NULL),
								 fCalifaCryCalDataCA(NULL),
	 							 fNonUniformity(0) {
}

//Virtual R3BCalifaPoint2CrystalCal: Destructor
R3BCalifaPoint2CrystalCal::~R3BCalifaPoint2CrystalCal() {  
  if (fCalifaPointDataCA) {
    fCalifaPointDataCA->Delete();
    delete fCalifaPointDataCA;
  }
  if (fCalifaCryCalDataCA) {
    fCalifaCryCalDataCA->Delete();
    delete fCalifaCryCalDataCA;
  }
  
  //Nota: delete TArrayF cryId[nHits], energyId[nHits], Sumaenergy[nHits];
  LOG(INFO) << "R3BCalifaPoint2CrystalCal: Delete instance" << FairLogger::endl;
}


// -----   Public method Init   --------------------------------------------
InitStatus R3BCalifaPoint2CrystalCal::Init() {
  //INPUT DATA
  //Point data
  FairRootManager* rootManager = FairRootManager::Instance();
  if (!rootManager) { return kFATAL;}
  
  fCalifaPointDataCA = (TClonesArray*)rootManager->GetObject("CrystalPoint");
  if (!fCalifaPointDataCA) { return kFATAL;}
  
  //OUTPUT DATA
  //CrystalCal data
  fCalifaCryCalDataCA = new TClonesArray("R3BCalifaCrystalCalData",10);
  rootManager->Register("CalifaCrystalCalData", "CALIFA Crystal Cal", fCalifaCryCalDataCA, kTRUE);
  
  return kSUCCESS;
}


// -----   Public method Execution   --------------------------------------------
void R3BCalifaPoint2CrystalCal::Exec(Option_t* option) {
  //Reading the Input -- Point data --
  Int_t nHits = fCalifaPointDataCA->GetEntries();
  if(!nHits) return;
  
  R3BCalifaPoint** pointData;
  pointData=new R3BCalifaPoint*[nHits];

  
  //Int_t crystalType;
  //Int_t crystalCopy;
  Int_t crystalId;
  //Double_t xIn;
  //Double_t yIn;
  //Double_t zIn;
  //Double_t xOut;
  //Double_t yOut;
  //Double_t zOut;
  //Double_t pxOut;
  //Double_t pyOut;
  //Double_t pzOut;
  Double_t Nf;
  Double_t Ns;
  Double_t time;
  Double_t energy;
  //TVector3 posIn;
  //TVector3 posOut;
  //TVector3 mom;  
  Double_t tot_energy=0.;
  Double_t Ary_cryId[nHits];
  Double_t Ary_energy[nHits];
  Double_t Ary_Nf[nHits];
  Double_t Ary_Ns[nHits];
  Double_t Sumaenergy[nHits];
  Int_t sorted[nHits];
  
  for(Int_t i = 0; i < nHits; i++) {
    pointData[i] = (R3BCalifaPoint*)(fCalifaPointDataCA->At(i));
    
    //crystalType = pointData[i]->GetCrystalType();
    //crystalCopy = pointData[i]->GetCrystalCopy();
    crystalId   = pointData[i]->GetCrystalId();
    //xIn = pointData[i]->GetXIn();
    //yIn = pointData[i]->GetYIn();
    //zIn	= pointData[i]->GetZIn();
    //xOut = pointData[i]->GetXOut();
    //yOut = pointData[i]->GetYOut();
    //zOut = pointData[i]->GetZOut();
    //pxOut = pointData[i]->GetPxOut();
    //pyOut = pointData[i]->GetPyOut();
    //pzOut = pointData[i]->GetPzOut();
    Nf = pointData[i]->GetNf();
    Ns = pointData[i]->GetNs();
    //pointData[i]->PositionIn(posIn);
    //pointData[i]->PositionOut(posOut);
    //pointData[i]->MomentumOut(mom);    
    time=pointData[i]->GetTime();
    energy=pointData[i]->GetEnergyLoss();
    
    Ary_cryId[i]=crystalId;
    Ary_energy[i]=energy;
    Ary_Nf[i]=Nf;
    Ary_Ns[i]=Ns;

  }
  
  TMath::Sort(nHits,Ary_cryId,sorted,kFALSE);//increasing order
  Int_t sameCrystal = 0;

  for(Int_t i=0; i<nHits; i++){
    if (i==0) AddCrystalCal(Ary_cryId[i], NUSmearing(Ary_energy[i]), Ary_Nf[i], Ary_Ns[i], time, tot_energy);
    else{
      if (Ary_cryId[sorted[i]==Ary_cryId[sorted[i-1]]]){
	sameCrystal++;
	((R3BCalifaCrystalCalData*)(fCalifaCryCalDataCA->At(i-sameCrystal)))->AddMoreEnergy(NUSmearing(Ary_energy[sorted[i]]));
	((R3BCalifaCrystalCalData*)(fCalifaCryCalDataCA->At(i-sameCrystal)))->AddMoreNf(Ary_Nf[sorted[i]]);
	((R3BCalifaCrystalCalData*)(fCalifaCryCalDataCA->At(i-sameCrystal)))->AddMoreNs(Ary_Ns[sorted[i]]);
      }
      else{ AddCrystalCal(Ary_cryId[i], NUSmearing(Ary_energy[i]), Ary_Nf[i], Ary_Ns[i], time, tot_energy);}
    }
    
  }  

}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BCalifaPoint2CrystalCal::EndOfEvent() {
  fCalifaCryCalDataCA->Clear();
  ResetParameters();
}

// -----   Public method Register   -------------------------------------------
void R3BCalifaPoint2CrystalCal::Register() {
  FairRootManager::Instance()->Register("CrystalCal", GetName(),
					fCalifaCryCalDataCA, kTRUE);
}

// -----   Public method Reset   ----------------------------------------------
void R3BCalifaPoint2CrystalCal::Reset() {
  fCalifaCryCalDataCA->Clear();
  ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method FinishEvent   ----------------------------------------------
void R3BCalifaPoint2CrystalCal::FinishEvent() {
}
// ----------------------------------------------------------------------------

// -----   Private method AddCrystalHit   --------------------------------------------
R3BCalifaCrystalCalData* R3BCalifaPoint2CrystalCal::AddCrystalCal(Int_t ident,
							   Double_t energy,
							   Double_t Nf,
							   Double_t Ns,
							   ULong64_t time,
							   Double_t tot_energy) {
  TClonesArray& clref = *fCalifaCryCalDataCA;
  Int_t size = clref.GetEntriesFast();
  LOG(INFO) << "-I- R3BCalifaPoint2CrystalCal: Adding CrystalCalData "
	    << " with unique identifier " << ident << " entering with "
	    << energy * 1e06 << " keV Nf=" << Nf << " Ns=" << Ns <<" Time="<<time
	    << " tot_energy=" <<tot_energy << FairLogger::endl;
  
  return new (clref[size]) R3BCalifaCrystalCalData(ident, energy, Nf, Ns, time, tot_energy);
  
}
// ----------------------------------------------------------------------------

// -----   Private method NUSmearing  --------------------------------------------
Double_t R3BCalifaPoint2CrystalCal::NUSmearing(Double_t inputEnergy) {
  // Very simple preliminary scheme where the NU is introduced as a flat random
  // distribution with limits fNonUniformity (%) of the energy value.
  //
  return gRandom->Uniform(inputEnergy - inputEnergy * fNonUniformity / 100,
			  inputEnergy + inputEnergy * fNonUniformity / 100);
}

// -----  Public method SetNonUniformity  ----------------------------------
void R3BCalifaPoint2CrystalCal::SetNonUniformity(Double_t nonU) {
  fNonUniformity = nonU;
  LOG(INFO) << "R3BCalifaPoint2CrystalCal::SetNonUniformity to " << fNonUniformity << " %" << FairLogger::endl;
}
