#include "R3BCalifaCrystalCal2Hit.h"
#include "TGeoMatrix.h"
#include "TMath.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TRandom.h"

#include "TGeoManager.h"

#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaGeometry.h"

using std::cerr;
using std::cout;
using std::endl;

R3BCalifaCrystalCal2Hit::R3BCalifaCrystalCal2Hit()
    : FairTask("R3B CALIFA CrystalCal to Hit Finder")
    , fCrystalHitCA(NULL)
    , fCalifaHitCA(NULL)
    , fOnline(kFALSE)
{
    fGeometryVersion = 2020; // BARREL+iPhos
    fThreshold = 0.;         // no threshold
    fDRThreshold = 15000;    // in keV, for real data
    fDeltaPolar = 0.25;
    fDeltaAzimuthal = 0.25;
    fDeltaAngleClust = 0;
    fClusteringAlgorithmSelector = 1;
    fParCluster1 = 0;
    nEvents = 0;
    fCalifaGeo = 0;
    // fCalifaHitFinderPar=0;
}

R3BCalifaCrystalCal2Hit::~R3BCalifaCrystalCal2Hit()
{
    LOG(INFO) << "R3BCalifaCrystalCal2Hit: Delete instance";
    if (fCalifaHitCA)
        delete fCalifaHitCA;
    // do not delete fCalifaGeo. It's a pointer to a static instance
}

void R3BCalifaCrystalCal2Hit::SetParContainers()
{
    // // Get run and runtime database
    // FairRunAna* run = FairRunAna::Instance();
    // if (!run) LOG(fatal) << "R3BCalifaCrystalCal2Hit::SetParContainers: No analysis run";

    // FairRuntimeDb* rtdb = run->GetRuntimeDb();
    // if (!rtdb) LOG(fatal) << "R3BCalifaCrystalCal2Hit::SetParContainers: No runtime database";

    // fCalifaHitFinderPar = (R3BCalifaCrystalCal2HitPar*)(rtdb->getContainer("R3BCalifaCrystalCal2HitPar"));
    // if ( fVerbose && fCalifaHitFinderPar ) {
    //   LOG(INFO) << "R3BCalifaCrystalCal2Hit::SetParContainers() ";
    //   LOG(INFO) << "Container R3BCalifaCrystalCal2HitPar loaded ";
    // }
}

InitStatus R3BCalifaCrystalCal2Hit::Init()
{
    LOG(INFO) << "R3BCalifaCrystalCal2Hit::Init ";
    assert(!fCalifaHitCA); // in case someone calls Init() twice.
    FairRootManager* ioManager = FairRootManager::Instance();
    if (!ioManager)
        LOG(fatal) << "Init: No FairRootManager";
    fCrystalHitCA = (TClonesArray*)ioManager->GetObject("CalifaCrystalCalData");

    // Register output array
    fCalifaHitCA = new TClonesArray("R3BCalifaHitData", 50);
    if (!fOnline)
    {
        ioManager->Register("CalifaHitData", "CALIFA Hit", fCalifaHitCA, kTRUE);
    }
    else
    {
        ioManager->Register("CalifaHitData", "CALIFA Hit", fCalifaHitCA, kFALSE);
    }

    fCalifaGeo = R3BCalifaGeometry::Instance(fGeometryVersion);

    // Parameter retrieval from par container
    // ...

    return kSUCCESS;
}

InitStatus R3BCalifaCrystalCal2Hit::ReInit() { return kSUCCESS; }

void R3BCalifaCrystalCal2Hit::Exec(Option_t* opt)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // ALGORITHMS FOR HIT FINDING
    ULong64_t hitTime = 0;
    Double_t energy = 0.;         // caloHits energy
    Double_t Nf = 0.;             // caloHits Nf
    Double_t Ns = 0.;             // caloHits Ns
    Double_t polarAngle = 0.;     // caloHits reconstructed polar angle
    Double_t azimuthalAngle = 0.; // caloHits reconstructed azimuthal angle
    Double_t rhoAngle = 0.;       // caloHits reconstructed rho
    Double_t angle1, angle2;

    bool* usedCrystalHits = NULL; // array to control the CrystalHits
    Int_t crystalsInHit = 0;      // used crystals in each CalifaHitData
    Double_t testPolar = 0;
    Double_t testAzimuthal = 0;
    Double_t testRho = 0;
    bool takeCrystalInCluster;

    R3BCalifaCrystalCalData** crystalHit = NULL;

    Int_t crystalHits; // Nb of CrystalHits in current event
    crystalHits = fCrystalHitCA->GetEntries();

    if (crystalHits == 0)
        return;

    crystalHit = new R3BCalifaCrystalCalData*[crystalHits];
    usedCrystalHits = new bool[crystalHits];
    for (Int_t i = 0; i < crystalHits; i++)
    {
        crystalHit[i] = (R3BCalifaCrystalCalData*)fCrystalHitCA->At(i);
        usedCrystalHits[i] = 0;
    }

    // For the moment a simple analysis... more to come
    Int_t crystalWithHigherEnergy = 0;
    Int_t unusedCrystals = crystalHits;

    // removing those crystals with energy below the threshold
    for (Int_t i = 0; i < crystalHits; i++)
    {
        if (crystalHit[i]->GetEnergy() < fThreshold)
        {
            usedCrystalHits[i] = 1;
            unusedCrystals--;
        }
    }

    // S444 for crystals with double reading!!
    // taking only the proton branch in the cluster
    for (Int_t i = 0; i < crystalHits; i++)
    {
        for (Int_t j = i; j < crystalHits; j++)
        {
            if (abs(crystalHit[i]->GetCrystalId() - crystalHit[j]->GetCrystalId()) == 5000)
            {
                if (crystalHit[i]->GetCrystalId() > crystalHit[j]->GetCrystalId())
                { // i is protonbranch
                    if (crystalHit[i]->GetEnergy() < fDRThreshold)
                    { // take j, the gammabranch
                        usedCrystalHits[i] = 1;
                        unusedCrystals--;
                    }
                    else
                    { // take i, the protonbranch
                        usedCrystalHits[j] = 1;
                        unusedCrystals--;
                    }
                }
                else
                { // j is protonbranch
                    if (crystalHit[i]->GetEnergy() < fDRThreshold)
                    { // take i, the gammabranch
                        usedCrystalHits[j] = 1;
                        unusedCrystals--;
                    }
                    else
                    { // take j, the protonbranch
                        usedCrystalHits[i] = 1;
                        unusedCrystals--;
                    }
                }
            }
        }
    }

    //  int n_clusters = 0;

    while (unusedCrystals > 0)
    {
        // First, finding the crystal with higher energy from the unused crystalHits
        for (Int_t i = 1; i < crystalHits; i++)
        {
            if (!usedCrystalHits[i] && crystalHit[i]->GetEnergy() > crystalHit[crystalWithHigherEnergy]->GetEnergy())
                crystalWithHigherEnergy = i;
        }

        usedCrystalHits[crystalWithHigherEnergy] = 1;
        unusedCrystals--;
        crystalsInHit++;

        // Second, energy and angles come from the crystal with the higher energy
        hitTime = crystalHit[crystalWithHigherEnergy]->GetTime();
        energy = crystalHit[crystalWithHigherEnergy]->GetEnergy();
        Nf = crystalHit[crystalWithHigherEnergy]->GetNf();
        Ns = crystalHit[crystalWithHigherEnergy]->GetNs();
        fCalifaGeo->GetAngles(
            crystalHit[crystalWithHigherEnergy]->GetCrystalId(), polarAngle, azimuthalAngle, rhoAngle);

        // Third, finding closest hits and adding their energy
        // Clusterization: you want to put a condition on the angle between the highest
        // energy crystal and the others. This is done by using the TVector3 classes and
        // not with different DeltaAngle on theta and phi, to get a proper solid angle
        // and not a "square" one.                    Enrico Fiori
        TVector3 refAngle(1, 0, 0); // EF
        refAngle.SetTheta(polarAngle);
        refAngle.SetPhi(azimuthalAngle);
        for (Int_t i = 0; i < crystalHits; i++)
        {
            if (!usedCrystalHits[i])
            {
                fCalifaGeo->GetAngles(crystalHit[i]->GetCrystalId(), testPolar, testAzimuthal, testRho);

                takeCrystalInCluster = false;

                TVector3 testAngle(1, 0, 0); // EF
                testAngle.SetTheta(testPolar);
                testAngle.SetPhi(testAzimuthal);
                // Check if the angle between the two vectors is less than the reference angle.
                switch (fClusteringAlgorithmSelector)
                {
                    case 1:
                    { // square window
                        // Dealing with the particular definition of azimuthal
                        // angles (discontinuity in pi and -pi)
                        if (azimuthalAngle + fDeltaAzimuthal > TMath::Pi())
                        {
                            angle1 = azimuthalAngle - TMath::Pi();
                            angle2 = testAzimuthal - TMath::Pi();
                        }
                        else if (azimuthalAngle - fDeltaAzimuthal < -TMath::Pi())
                        {
                            angle1 = azimuthalAngle + TMath::Pi();
                            angle2 = testAzimuthal + TMath::Pi();
                        }
                        else
                        {
                            angle1 = azimuthalAngle;
                            angle2 = testAzimuthal;
                        }
                        if (TMath::Abs(polarAngle - testPolar) < fDeltaPolar &&
                            TMath::Abs(angle1 - angle2) < fDeltaAzimuthal)
                        {
                            takeCrystalInCluster = true;
                        }
                        break;
                    }
                    case 2: // round window
                        // The angle is scaled to a reference distance (e.g. here is
                        // set to 35 cm) to take into account Califa's non-spherical
                        // geometry. The reference angle will then have to be defined
                        // in relation to this reference distance: for example, 10° at
                        // 35 cm corresponds to ~6cm, setting a fDeltaAngleClust=10
                        // means that the gamma rays will be allowed to travel 6 cm in
                        // the CsI, no matter the position of the crystal they hit.
                        if (((refAngle.Angle(testAngle)) * ((testRho + rhoAngle) / (35. * 2.))) < fDeltaAngleClust)
                        {
                            takeCrystalInCluster = true;
                        }
                        break;
                    case 3:
                    { // round window scaled with energy
                        // The same as before but the angular window is scaled
                        // according to the energy of the hit in the higher energy
                        // crystal. It needs a parameter that should be calibrated.
                        Double_t fDeltaAngleClustScaled =
                            fDeltaAngleClust * (crystalHit[crystalWithHigherEnergy]->GetEnergy() * fParCluster1);
                        if (((refAngle.Angle(testAngle)) * ((testRho + rhoAngle) / (35. * 2.))) <
                            fDeltaAngleClustScaled)
                        {
                            takeCrystalInCluster = true;
                        }
                        break;
                    }
                    case 4: // round window scaled with the energy of the _two_ hits
                            //(to be tested and implemented!!)
                        // More advanced: the condition on the distance between the
                        // two hits is function of the energy of both hits
                        break;
                }

                if (takeCrystalInCluster)
                {
                    energy += crystalHit[i]->GetEnergy();
                    Nf += crystalHit[i]->GetNf();
                    Ns += crystalHit[i]->GetNs();
                    usedCrystalHits[i] = 1;
                    unusedCrystals--;
                    crystalsInHit++;
                }
            }
        }

        AddHit(crystalsInHit, energy, Nf, Ns, polarAngle, azimuthalAngle, hitTime);

        crystalsInHit = 0; // reset for next CalifaHitData

        // Finally, setting crystalWithHigherEnergy to the first unused
        // crystalHit (for the next iteration)
        for (Int_t i = 0; i < crystalHits; i++)
        {
            if (!usedCrystalHits[i])
            {
                crystalWithHigherEnergy = i;
                break;
            }
        }
    }

    //  std::cout << "# " << n_clusters << "\n--------\n";

    if (crystalHit)
        delete[] crystalHit;
    if (usedCrystalHits)
        delete[] usedCrystalHits;
}

void R3BCalifaCrystalCal2Hit::Reset()
{
    // Clear the CA structure
    LOG(DEBUG) << "Clearing CalifaHitData Structure";
    if (fCalifaHitCA)
        fCalifaHitCA->Clear();
}

void R3BCalifaCrystalCal2Hit::Finish() {}

void R3BCalifaCrystalCal2Hit::SelectGeometryVersion(Int_t version)
{
    fGeometryVersion = version;
    LOG(INFO) << "R3BCalifaCrystalCal2Hit::SelectGeometryVersion to " << fGeometryVersion;
}

void R3BCalifaCrystalCal2Hit::SetCrystalThreshold(Double_t thresholdEne)
{
    fThreshold = thresholdEne;

    LOG(INFO) << "R3BCalifaCrystalCal2Hit::SetCrystalThreshold to " << fThreshold << " keV.";
}

void R3BCalifaCrystalCal2Hit::SetDRThreshold(Double_t DRthresholdEne)
{
    fDRThreshold = DRthresholdEne;
    LOG(INFO) << "R3BCalifaCrystalCal2Hit::SetDRThreshold to " << fDRThreshold << " keV.";
}

void R3BCalifaCrystalCal2Hit::SetClusteringAlgorithm(Int_t ClusteringAlgorithmSelector, Double_t ParCluster1)
{
    // Select the clustering algorithm and the parameters of some of them
    // ClusteringAlgorithmSelector = 1  ->  square window
    // ClusteringAlgorithmSelector = 2  ->  round window
    // ClusteringAlgorithmSelector = 3  ->  advanced round window with opening proportional to the
    //                                     energy of the hit, need ParCluster1
    // ClusteringAlgorithmSelector = 4  ->  advanced round window with opening proportional to the
    //                                     energy of the two hit, need ParCluster1 NOT YET IMPLEMENTED!
    fClusteringAlgorithmSelector = ClusteringAlgorithmSelector;
    fParCluster1 = ParCluster1;
}

void R3BCalifaCrystalCal2Hit::SetAngularWindow(Double_t deltaPolar, Double_t deltaAzimuthal, Double_t DeltaAngleClust)
{
    //
    // Set the angular window open around the crystal with the largest energy
    // to search for additional crystal hits and addback to the same cal hit
    // [0.25 around 14.3 degrees, 3.2 for the complete calorimeter]
    fDeltaPolar = deltaPolar;
    fDeltaAzimuthal = deltaAzimuthal;
    fDeltaAngleClust = DeltaAngleClust;
}

R3BCalifaHitData* R3BCalifaCrystalCal2Hit::AddHit(UInt_t Nbcrystals,
                                                  Double_t ene,
                                                  Double_t Nf,
                                                  Double_t Ns,
                                                  Double_t pAngle,
                                                  Double_t aAngle,
                                                  ULong64_t time)
{
    TClonesArray& clref = *fCalifaHitCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BCalifaHitData(Nbcrystals, ene, Nf, Ns, pAngle, aAngle, time);
}

/*
 Double_t GetCMEnergy(Double_t theta, Double_t energy){
 //
 // Calculating the CM energy from the lab energy and the polar angle
 //
 Double_t beta = 0.8197505718204776;  //beta is 0.8197505718204776
 Double_t gamma = 1/sqrt(1-beta*beta);
 //Lorenzt boost correction
 //E' = gamma E + beta gamma P|| = gamma E + beta gamma P cos(theta)
 //In photons E=P
 Double_t energyCorrect = gamma*energy - beta*gamma*energy*cos(theta);

 return energyCorrect;
 }
*/

ClassImp(R3BCalifaCrystalCal2Hit)
