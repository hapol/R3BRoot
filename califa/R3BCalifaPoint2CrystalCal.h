// ------------------------------------------------------------
// -----      R3BCalifaPoint2CrystalCal header file           		-----
// -----      Created 1/03/2018 by E.Galiana          		-----
// -----                                                  -----
// ------------------------------------------------------------


#ifndef R3BCALIFADIGITIZER_H
#define R3BCALIFADIGITIZER_H

#include "R3BCalifa.h"
#include "FairTask.h"
#include "TClonesArray.h"
#include "R3BCalifaPoint.h"
#include "R3BCalifaCrystalCalData.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "string"

class R3BCalifaPoint2CrystalCal : public FairTask
{

 public:

  /** Default constructor **/  
  R3BCalifaPoint2CrystalCal();
  
 /** Standard contructor **/
 R3BCalifaPoint2CrystalCal(const TString& geoFile);

 /** Destructor **/
 ~R3BCalifaPoint2CrystalCal();

 /** Virtual method Init **/
 virtual InitStatus Init();

 /** Virtual method Exec **/
 virtual void Exec(Option_t* opt);
  
 /** Virtual method EndOffEvent **/ 
 virtual void EndOfEvent();
  
 /** Virtual method Register **/
 virtual void Register();
 
 /** Virtual method Reset **/
 virtual void Reset();
 
 /** Virtual method FinishEvent **/
 virtual void FinishEvent();
   
 //Accessors functions 
 void SetEnergyResolution(Double_t resolution){fResolution=resolution;}

  inline void ResetParameters()
	{
	};
 
	 /** Private method AddCrystalCal
   **
   ** Adds a CalifaCrystalCal data
   **/
  R3BCalifaCrystalCalData* AddCrystalCal(Int_t ident, 
  																			 Double_t energy,
  																			 Double_t Nf,
  																			 Double_t Ns,
  																			 ULong64_t time,
  																			 Double_t tot_energy=0.);
  																			   																			
  protected:
  
	  TClonesArray* fCalifaPointDataCA; //!  The crystal hit collection
	  TClonesArray* fCalifaCryCalDataCA; /**< Array with CALIFA Cal- output data. >*/
	  
  private:
      
    Int_t fResolution;  						//Resolution of the CALIFA calorimeter
		Double_t fNonUniformity;				//Adding some non-uniformity preliminary description
		
		
    /** Private method NUSmearing
     **
     ** Smears the energy according to some non-uniformity distribution
     ** Very simple preliminary scheme where the NU is introduced as a flat random
     ** distribution with limits fNonUniformity (%) of the energy value.
     **/
     Double_t NUSmearing(Double_t inputEnergy);
     
     /** Public method SetNonUniformity
     **
     ** Defines the fNonUniformity parameter in % deviation from the central value
     *@param nonU  Double parameter setting the maximum non-uniformity allowed
     **/
     void SetNonUniformity(Double_t nonU);
	
   ClassDef(R3BCalifaPoint2CrystalCal,1);
  
};

#endif
