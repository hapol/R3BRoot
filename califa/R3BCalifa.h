#ifndef R3BCALIFA_H
#define R3BCALIFA_H

#include "R3BDetector.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include <map>

class TClonesArray;
class R3BCalifaPoint;
class FairVolume;
class TGeoRotation;
class R3BCalifaGeometry;

class R3BCalifa : public R3BDetector
{

  public:
    /** Default constructor */
    R3BCalifa();

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param trans   position
     *@param rot     rotation
     */
    R3BCalifa(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot = TGeoRotation());

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param combi   position + rotation
     */
    R3BCalifa(const TString& geoFile, const TGeoCombiTrans& combi = TGeoCombiTrans());

    /** Destructor **/
    ~R3BCalifa();

    /** Virtual method ProcessHits
     **
     ** Defines the action to be taken when a step is inside the
     ** active volume. Creates a R3BCalifaPoint or a
     ** R3BCalifaCrystalCalData and adds it to the collection.
     *@param vol  Pointer to the active volume
     **/
    virtual Bool_t ProcessHits(FairVolume* vol = 0);

    /** Virtual method BeginEvent
     **
     ** If verbosity level is set, print hit collection at the
     ** end of the event and resets it afterwards.
     **/
    virtual void BeginEvent();

    /** Virtual method EndOfEvent
     **
     ** Added support for R3BCalifaCrystalCalData
     ** If verbosity level is set, print hit collection at the
     ** end of the event and resets it afterwards.
     **/
    virtual void EndOfEvent();

    /** Virtual method Register
     **
     ** Registers the hit collection in the ROOT manager.
     **/
    virtual void Register();

    /** Accessor to the hit collection **/
    virtual TClonesArray* GetCollection(Int_t iColl) const;

    /** Virtual method Print
     **
     ** Screen output of hit collection.
     **/
    virtual void Print(Option_t* option = "") const;

    /** Virtual method Reset
     **
     ** Clears the hit collection
     **/
    virtual void Reset();

    /** Virtual method CopyClones
     **
     ** Copies the hit collection with a given track index offset
     *@param cl1     Origin
     *@param cl2     Target
     *@param offset  Index offset
     **/
    virtual void CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset);

    virtual Bool_t CheckIfSensitive(std::string name);

    /** Select the version of the CALIFA geometry
     **
     *@param version Geometry version (final BARREL+iPhos: 2020)
     **/
    void SelectGeometryVersion(Int_t version);

    virtual void Initialize();
    virtual void SetSpecialPhysicsCuts() {}

  private:
    /** Track information to be stored until the track leaves the
    active volume. **/
    Int_t fTrackID;                 //!  track index
    Int_t fTrackPID;                //!  particle identification
    Int_t fVolumeID;                //!  volume id
    Int_t fParentTrackID;           //!  parent track index
    Int_t fUniqueID;                //!  particle unique id (e.g. if Delta electron, fUniqueID=9)
    TLorentzVector fPosIn, fPosOut; //!  position
    TLorentzVector fMomIn, fMomOut; //!  momentum
    Double32_t fTime;               //!  time
    Double32_t fLength;             //!  length
    Double32_t fELoss;              //!  energy loss
    Double32_t fNf;                 //!  fast CsI(Tl) amplitude
    Double32_t fNs;                 //!  slow CsI(Tl) amplitude
    Int_t fPosIndex;                //!
    Int_t fNSteps;                  //!  Number of steps in the active volume
    Double32_t fEinc;               //!  Total incident energy
    TList* flGeoPar;                //!

    TClonesArray* fCalifaCollection;           //!  The point collection

    //! Defining functions for energy to light output calculation
    //    TF1 *tf_p_dNs;    //!
    //    TF1 *tf_p_dNf;    //!
    //    TF1 *tf_g_dNs;    //!
    //    TF1 *tf_g_dNf;    //!

    TF1* tf_dNf_dE; //!
    TF1* tf_dNs_dE; //!

    // Selecting the geometry of the CALIFA calorimeter (final BARREL+iPhos: 2020)
    Int_t fGeometryVersion;

    Double_t fCsIDensity;

    R3BCalifaGeometry* fCalifaGeo;

    /** Private method AddPoint
     **
     ** Adds a CalifaPoint to the HitCollection
     **/
    R3BCalifaPoint* AddPoint(Int_t trackID,
                             Int_t detID,
                             Int_t ident,
                             TVector3 posIn,
                             TVector3 pos_out,
                             TVector3 momIn,
                             TVector3 momOut,
                             Double_t time,
                             Double_t length,
                             Double_t eLoss,
                             Double_t Nf,
                             Double_t Ns);

    /** Private method AddCrystalCal
     **
     ** Adds a CalifaCrystalCalData to the HitCollection
     **/
    R3BCalifaCrystalCalData* AddCrystalHit(Int_t ident,
                                           Double_t energy,
                                           Double_t Nf,
                                           Double_t Ns,
                                           Double_t time,
                                           Double_t tot_energy);

    /** Private method ResetParameters
     **
     ** Resets the private members for the track parameters
     **/
    void ResetParameters();

    TGeoRotation* createMatrix(Double_t phi, Double_t theta, Double_t psi);

    ClassDef(R3BCalifa, 7);
};

inline void R3BCalifa::ResetParameters()
{
    fTrackID = fVolumeID = fParentTrackID = fTrackPID = fUniqueID = 0;
    fPosIn.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fPosOut.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fMomIn.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fMomOut.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fTime = fLength = fELoss = fNf = fNs = fEinc = 0;
    fPosIndex = 0;
    fNSteps = 0;
};

#endif
