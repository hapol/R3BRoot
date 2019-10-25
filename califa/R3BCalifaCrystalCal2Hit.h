#ifndef R3BCALIFACRYSTALCAL2HIT_H
#define R3BCALIFACRYSTALCAL2HIT_H

#include "FairTask.h"
#include "R3BCalifaGeometry.h"
#include "R3BCalifaHitData.h"

class TClonesArray;

class R3BCalifaCrystalCal2Hit : public FairTask
{

  public:
    /** Default constructor **/
    R3BCalifaCrystalCal2Hit();

    /** Destructor **/
    ~R3BCalifaCrystalCal2Hit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method Reset **/
    virtual void Reset();

    /** Public method SelectGeometryVersion
     **
     ** Defines the geometry
     *@param version  Integer parameter used to select the geometry:
     ** (see documentation /r3broot/cal/perlScripts/README))
     **/
    void SelectGeometryVersion(Int_t version);

    /** Public method SetAngularWindow
     **
     ** Sets the angular window open around the crystal with the largest energy
     ** to search for additional crystal hits and addback to the same cal hit
     **/
    void SetAngularWindow(Double_t deltaPolar, Double_t deltaAzimuthal, Double_t DeltaAngleClust = 0.);

    /** Public method SetClusteringAlgorithm
     **
     ** Select the clustering algorithm to be used and set the parameter to be used by one of them
     //  1  ->  square window
     //  2  ->  round window
     //  3  ->  advanced round window with opening proportional to the
     //         energy of the hit, need ParCluster1
     //  4  ->  advanced round window with opening proportional to the
     //         energy of the two hit, need ParCluster1 NOT IMPLEMENTED YET!
     **
     **/
    void SetClusteringAlgorithm(Int_t ClusteringAlgorithmSelector, Double_t ParCluster1);

    /** Public method SetCrystalThreshold
     **
     ** Defines the minimum energy requested in a crystal to be considered in a calorimeter Hit
     *@param thresholdEne  Double parameter used to set the threshold
     **/
    void SetCrystalThreshold(Double_t thresholdEne);

    /** Public method SetDRThreshold (for double reading)
     **
     ** Defines the minimum energy requested in a crystal to be considered in a calorimeter Hit
     *@param thresholdEne  Double parameter used to set the threshold
     **/
    void SetDRThreshold(Double_t DRthresholdEne);

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /** Virtual method Finish **/
    virtual void Finish();

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  protected:
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    TClonesArray* fCrystalHitCA;
    TClonesArray* fCalifaHitCA;

    Bool_t fOnline;                     // Selector for online data storage
    Int_t fGeometryVersion;             // Selecting the geometry of the CALIFA calorimeter
    Double_t fDeltaPolar;               // Angular window (polar angle)
    Double_t fDeltaAzimuthal;           // Angular window (azimuthal angle)
    Double_t fDeltaAngleClust;          // Angular opening used for the cluster condition
    Int_t fClusteringAlgorithmSelector; // Clustering algorithm selector
    Double_t fParCluster1;              // Clustering parameter 1
    Double_t fThreshold;                // Minimum energy requested in a crystal to be included in a Cal
    Double_t fDRThreshold;              // Threshold for selecting gamma or proton branch in double reading channels

    // Parameter class
    // R3BCalifaHitPar* fCalifaHitPar;

    R3BCalifaGeometry* fCalifaGeo;

  private:
    UInt_t nEvents;

    /** Private method AddHit
     **
     ** Adds a CalifaHit to the HitCollection
     **/
    R3BCalifaHitData* AddHit(UInt_t Nbcrystals,
                             Double_t ene,
                             Double_t Nf,
                             Double_t Ns,
                             Double_t pAngle,
                             Double_t aAngle,
                             ULong64_t time);

    ClassDef(R3BCalifaCrystalCal2Hit, 2);
};

#endif
