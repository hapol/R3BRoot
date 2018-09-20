// -------------------------------------------------------------------------
// -----                        R3BCalifa source file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -----       			Modified by 15/12/16 by P.Cabanelas         		   -----
// -----       			Modified by 1/03/18 by E.Galiana		         		   -----
// -------------------------------------------------------------------------


////////////////////////////////////////////////////////////////////////////
// Pablo Cabanelas - 15-12-2016:
// Changing classes name from Calo to Califa acording to name convention
//
// EliGaliana - 1/03/2018
// Removing old geometry versions and filling only PointData, not CrystalCal
////////////////////////////////////////////////////////////////////////////

#include "R3BCalifa.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "R3BCalifaPoint.h"
#include "R3BMCStack.h"
#include "TClonesArray.h"
#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoCone.h"
#include "TGeoMCGeometry.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPgon.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TMCProcess.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include "TVirtualMCStack.h"
#include "TGeoNode.h"
#include "R3BCalifaGeometry.h"

#include <iostream>
#include <stdlib.h>

using std::cout;
using std::cerr;
using std::endl;

#define U_MEV 931.4940954

inline double BETA(const double M, const double E_kin) { return sqrt(1. - M * M / ((M + E_kin) * (M + E_kin))); }

inline double GAMMA(const double M, const double E_kin) { return 1. + E_kin / M; }

R3BCalifa::R3BCalifa()
    : R3BCalifa("")
{
}

R3BCalifa::R3BCalifa(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BCalifa(geoFile, { trans, rot })
{
}

R3BCalifa::R3BCalifa(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BCalifa", kCALIFA, geoFile, combi)
{
    ResetParameters();
    fCrystal = NULL;
    fCaloCollection = new TClonesArray("R3BCalifaPoint");
    fPosIndex = 0;
    kGeoSaved = kFALSE;
    flGeoPar = new TList();
    flGeoPar->SetName(GetName());
    fNonUniformity = 0.;
    fGeometryVersion = 1;

    //  tf_p_dNs = new TF1("tf_p_dNs","-[0]*[1]*exp(-[1]*(x-[3]))+[2]",0,1000);
    //  tf_p_dNf = new TF1("tf_p_dNf","-[0]*[1]*exp(-[1]*(x-[3]))+[2]",0,1000);
    //  tf_g_dNs = new TF1("tf_g_dNs","[2]",0,1000);
    //  tf_g_dNf = new TF1("tf_g_dNf","[2]",0,1000);
    //
    //  tf_p_dNs->SetParameter(0,18.88 / 7.46);
    //  tf_p_dNs->SetParameter(1,0.0868);
    //  tf_p_dNs->SetParameter(2,4.228 / 7.46);
    //  tf_p_dNs->SetParameter(3,0);
    //
    //  tf_p_dNf->SetParameter(0,-18.88 / 7.46);
    //  tf_p_dNf->SetParameter(1,0.0868);
    //  tf_p_dNf->SetParameter(2,3.232 / 7.46);
    //  tf_p_dNf->SetParameter(3,0);
    //
    //  // tf_p_dNs->SetParameter(0,18.88);
    //  // tf_p_dNs->SetParameter(1,0.0868);
    //  // tf_p_dNs->SetParameter(2,4.228);
    //  // tf_p_dNs->SetParameter(3,4.117);
    //  // tf_p_dNs->SetParameter(4,4.259);
    //
    //  // tf_p_dNf->SetParameter(0,-32.66);
    //  // tf_p_dNf->SetParameter(1,0.07729);
    //  // tf_p_dNf->SetParameter(2,3.155);
    //  // tf_p_dNf->SetParameter(3,0);
    //  // tf_p_dNf->SetParameter(4,-3.947);
    //
    //  tf_g_dNs->SetParameter(2, tf_p_dNs->GetParameter(2));
    //  tf_g_dNf->SetParameter(2, tf_p_dNf->GetParameter(2));

    tf_dNf_dE = new TF1("tf_dNf_dE", "1./([0]+[1]*(x^[2])+[3]/(x^[4]))");
    tf_dNs_dE = new TF1("tf_dNs_dE", "1./([0]+[1]*(x^[2])+[3]/(x^[4]))");

    tf_dNf_dE->SetParameters(-1.79, 1.36e-2, 7.84e-1, 4.97, 1.75e-1);
    tf_dNs_dE->SetParameters(-1.24e2, 6.3e-3, 1.27, 1.262e2, 2.3e-3);
}

R3BCalifa::~R3BCalifa()
{

    if (flGeoPar)
    {
        delete flGeoPar;
    }
    if (fCaloCollection)
    {
        fCaloCollection->Delete();
        delete fCaloCollection;
    }

    //  delete tf_p_dNs;
    //  delete tf_p_dNf;
    //  delete tf_g_dNs;
    //  delete tf_g_dNf;
    delete tf_dNf_dE;
    delete tf_dNs_dE;
}

void R3BCalifa::Initialize()
{
    FairDetector::Initialize();

    LOG(INFO) << "R3BCalifa: initialisation" << FairLogger::endl;
    LOG(DEBUG) << "-I- R3BCalifa: Vol (McId) def" << FairLogger::endl;

    TGeoVolume* vol = gGeoManager->GetVolume("CalifaWorld");
    vol->SetVisibility(kFALSE);
}

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BCalifa::ProcessHits(FairVolume* vol)
{
    bool Print_Volinfo=kFALSE;//Info FairVolume* vol
    if(Print_Volinfo)
    {	
				Int_t fMCid = vol->getMCid();
				Int_t fCopyNo = vol->getCopyNo();
				Int_t fMotherId = vol->getMotherId();
				Int_t fMotherCopyNo = vol->getMotherCopyNo();		
				cout<<"-------   FairVolume info   ------"<<endl;
				cout<<"fMCid ="<<fMCid<<endl;
				cout<<"fCopyNo ="<<fCopyNo<<endl;
				cout<<"fMotherId ="<<fMotherId<<endl;
				cout<<"fMotherCopyNo ="<<fMotherCopyNo<<endl;			
		}
		
		
    // While tracking a single particle within a crystal (volume)
    // we can rely on the latest crystal information for each step
    if (gMC->IsTrackEntering() || fCrystal == NULL)
    {
        // Try to get crystal information from hash table
        // TODO: Still a performance benefit to use hash table?
				
				// Note: there are two crystals with the same NodeId,
				// which cause an error in the map filling 
				// check out this error in future geometry versions   
        gGeoManager->cd(gMC->CurrentVolPath());
        Int_t nodeId = gGeoManager->GetNodeId();
                
        std::map<Int_t, sCrystalInfo>::iterator it = fCrystalMap.find(nodeId);
        
        if (it == fCrystalMap.end())
        {
              
            // Not found in map => Create crystal information for crystal
            sCrystalInfo tmpInfo;
            memset(&tmpInfo, 0, sizeof(sCrystalInfo));//Fill by 0 sCrystalInfo 


						/*R3BCalifaGeometry* CalifaGeo;
						const char* path;
						int id;
						id=32;
 						Double_t polar;
						Double_t azimuthal;
						Double_t rho;
						
						path=CalifaGeo->GetCrystalVolumePath(id);
						cout<<"    HEREEEEE!  path="<<path<<endl;

						CalifaGeo->GetAngles(id, polar, azimuthal, rho);
						cout<<"id="<<id<<"  polar="<<polar<<"  azimuthal="<<azimuthal<<" rho="<<rho<<endl<<endl;*/

            if (GetCrystalInfo(tmpInfo))
            {
                fCrystal = &(fCrystalMap[nodeId] = tmpInfo);
                fCrystal->density = gGeoManager->GetCurrentVolume()->GetMaterial()->GetDensity();
            }
            else
                fCrystal = NULL;
        }
        else
            fCrystal = &(it->second);
            	
    }

    if (fCrystal == NULL)
    {
        // Still no crystal info at this point?
        // -> Something went wrong, but the user should
        //    already have been informed
        // -> Silently bail out
        return kFALSE;
    }

    if (fVerboseLevel > 1)
        LOG(INFO) << "R3BCalifa: Processing Points in Alveolus Nb " << fCrystal->volIdAlv << ", copy Nb "
                  << fCrystal->cpAlv << ", crystal copy Nb " << fCrystal->cpCry << " and unique crystal identifier "
                  << fCrystal->crystalId << FairLogger::endl;

    if (gMC->IsTrackEntering())
    {
        fELoss = 0.;
        fNf = 0.;
        fNs = 0.;
        fNSteps = 0; // FIXME
        fTime = gMC->TrackTime() * 1.0e09;
        fLength = gMC->TrackLength();
        gMC->TrackPosition(fPosIn);
        gMC->TrackMomentum(fMomIn);
        fEinc = gMC->Etot() - gMC->TrackMass(); // be aware!! Relativistic mass!
    }

    // Sum energy loss for all steps in the active volume
    Double_t dE = gMC->Edep() * 1000.;                          // in MeV
    Double_t post_E = (gMC->Etot() - gMC->TrackMass()) * 1000.; // in MeV
    TString ptype = gMC->GetStack()->GetCurrentTrack()->GetName();
    Double_t dx = gMC->TrackStep() * fCrystal->density;

    Double_t M_in = gMC->TrackMass() * 1000.;
    Double_t A_in = M_in / U_MEV;
    Double_t Z_in = gMC->TrackCharge();

    const double Z_CsI = 54.;
    const double A_CsI = 129.905;   // g/mol
    const double E_delta = 5.30227; // MeV
    const double m_e = .5109989461; // MeV
    const double slope_e = 1.33055;
    const double K = .307075; // MeV cm**2/mol
    // quenching
    const double q_1 = 0.0396113;
    const double q_2 = -0.0828619;
    const double q_3 = 0.780435;

    fELoss += dE / 1000.; // back to GeV

    if (dE > 0 && dx > 0)
    {

        //    cout << ptype << " E = " << post_E << " MeV, dE = " << dE << " MeV, dx = " << dx << " g/cm**2" << ", dE/dx
        //    = " << (dE/dx) << " MeV cm**2/g" << endl;

        if (fCrystal->fEndcapIdentifier == 1)
        {
            // CC Phoswich
            if (fCrystal->fPhoswichIdentifier == 1)
            {
                // LaBr
                fNf += dE;
            }
            else if (fCrystal->fPhoswichIdentifier == 2)
            {
                // LaCl
                fNs += dE;
            }
            else
            {
                LOG(ERROR) << "R3BCalifa: fPhoswichIdentifier not valid in R3BCalifa::ProcessHits(). "
                           << FairLogger::endl;
            }
        }
        else if (fCrystal->fEndcapIdentifier == 0)
        {
            //    if (ptype == "e-" || ptype == "e+" || ptype == "gamma") {
            //      fNs += tf_g_dNs->Integral(post_E, post_E + dE);
            //      fNf += tf_g_dNf->Integral(post_E, post_E + dE);
            //    } else if(ptype == "proton") {
            //      fNs += tf_p_dNs->Integral(post_E, post_E + dE);
            //      fNf += tf_p_dNf->Integral(post_E, post_E + dE);
            //    }

            if (ptype != "gamma" && post_E >= A_in * E_delta)
            {
                double beta_cut = BETA(M_in, A_in * E_delta);
                double gamma_cut = GAMMA(M_in, A_in * E_delta);
                double beta = BETA(M_in, post_E);
                double gamma = GAMMA(M_in, post_E);
                double T_cut = 2. * m_e * beta_cut * beta_cut * gamma_cut * gamma_cut /
                               (1. + 2. * gamma_cut * m_e / M_in + (m_e / M_in) * (m_e / M_in));
                double T_max = 2. * m_e * beta * beta * gamma * gamma /
                               (1. + 2. * gamma * m_e / M_in + (m_e / M_in) * (m_e / M_in));
                double C = 0.5 * K * Z_in * Z_in * Z_CsI / (A_CsI * beta * beta);

                // quenching
                double part1 =
                    q_1 / q_2 *
                    (1 / T_max - 1 / T_cut + (log(T_cut / T_max) + log((T_max - q_2) / (T_cut - q_2)) / q_2));
                double part2 =
                    q_1 * beta * beta / T_max * (log(T_cut / T_max) + log((T_max - q_2) / (T_cut - q_2)) / q_2);
                double N = 1 / T_cut - 1 / T_max - (beta * beta) / T_max * log(T_max / T_cut);
                double part3 = q_3 * N;
                double scaling = 1.;

                double dE_dxe = C * (log(T_max / T_cut) - beta * beta * (T_max - T_cut) / T_max);
                double dE_e = dE_dxe * dx;
                if (dE_e > dE)
                {
                    dE_e = dE;
                }
                if (T_max < 2.)
                {
                    scaling = (part1 + part2 + part3) / N;
                }
                else
                {
                    scaling = q_3;
                }

                //            cout << "  dE_e = " << dE_e << " MeV, N = " << N << ", scaling = " << scaling << endl;

                // std::cout << T_max << "  " << scaling << std::endl;
                fNf += (dE_e * scaling / (1 + slope_e)) / 1000.;
                fNs += (dE_e * scaling / (1. / slope_e + 1)) / 1000.;
                dE -= dE_e;
            }

            fNf += tf_dNf_dE->Eval(dE / dx) * dE / 1000.;
            fNs += tf_dNs_dE->Eval(dE / dx) * dE / 1000.;
        }
        else
        {
            LOG(ERROR) << "R3BCalifa: fEndcapIdentifier not valid in R3BCalifa::ProcessHits(). " << FairLogger::endl;
        }
    }

    fNSteps++;

    // Set additional parameters at exit of active volume. Create R3BCalifaPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {

        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fParentTrackID = gMC->GetStack()->GetCurrentParentTrackNumber();
        fVolumeID = vol->getMCid();
        fTrackPID = gMC->TrackPid();
        fUniqueID = gMC->GetStack()->GetCurrentTrack()->GetUniqueID();

        gMC->TrackPosition(fPosOut);
        gMC->TrackMomentum(fMomOut);

        if (fELoss == 0.)
            return kFALSE;

 
         if (gMC->IsTrackExiting()){
         
            const Double_t* oldpos;
            const Double_t* olddirection;
            Double_t newpos[3];
            Double_t newdirection[3];
            Double_t safety;
            gGeoManager->FindNode(fPosOut.X(),fPosOut.Y(),fPosOut.Z());
            oldpos = gGeoManager->GetCurrentPoint();
            olddirection = gGeoManager->GetCurrentDirection();
            
            for(Int_t i=0; i<3; i++){
                newdirection[i] = -1*olddirection[i];
            }
            gGeoManager->SetCurrentDirection(newdirection);
            safety = gGeoManager->GetSafeDistance();
            gGeoManager->SetCurrentDirection(-newdirection[0],
                               -newdirection[1],
                               -newdirection[2]);
            for (Int_t i=0; i<3; i++) {
                newpos[i] = oldpos[i] - (3*safety*olddirection[i]);
            }
            fPosOut.SetX(newpos[0]);
            fPosOut.SetY(newpos[1]);
            fPosOut.SetZ(newpos[2]);
           }
           AddPoint(fTrackID, fVolumeID, fCrystal->crystalType , fCrystal->crystalCopy , fCrystal->crystalId,
                   TVector3(fPosIn.X(),   fPosIn.Y(),   fPosIn.Z()),
                   TVector3(fPosOut.X(),  fPosOut.Y(),  fPosOut.Z()),
                   TVector3(fMomIn.Px(),  fMomIn.Py(),  fMomIn.Pz()),
                   TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
                   fTime, fLength, fELoss, fNf, fNs);
            

        // Increment number of CalifaPoints for this track
        R3BStack* stack = (R3BStack*)gMC->GetStack();
        stack->AddPoint(kCALIFA);

        ResetParameters();

    }

    return kTRUE;
}
// ----------------------------------------------------------------------------
// void R3BCalifa::SaveGeoParams(){
//
//  cout << " -I Save STS geo params " << endl;
//
//  TFolder *mf = (TFolder*) gDirectory->FindObjectAny("cbmroot");
//  cout << " mf: " << mf << endl;
//  TFolder *stsf = NULL;
//  if (mf ) stsf = (TFolder*) mf->FindObjectAny(GetName());
//  cout << " stsf: " << stsf << endl;
//  if (stsf) stsf->Add( flGeoPar0 ) ;
//  FairRootManager::Instance()->WriteFolder();
//  mf->Write("cbmroot",TObject::kWriteDelete);
//}




Bool_t R3BCalifa::GetCrystalInfo(sCrystalInfo& info)
{

    // Getting the Infos from Crystal Volumes
    Int_t cp1 = -1;
    Int_t volId1 = -1;
    Int_t cpAlv = -1;
    Int_t cpSupAlv = -1;
    Int_t volIdAlv = -1;
    Int_t volIdSupAlv = -1;
    Int_t cpCry = -1;
    Int_t volIdCry = -1;

    //Crystals Ids:	
    //Num.Id of Geant:volId1, volIdCry1, volIdAlv, volIdSupAlv
    //Copy num. of each volume:cp1, cpCry, cpAlv, cpSupAlv
    //																								BARREL info:								ENDCAP info:
    const char* bufferName = gMC->CurrentVolName();	//Current Volume name					Current Volume name
    volId1 = gMC->CurrentVolID(cp1);								//Crystal											Crystal	
    volIdCry = gMC->CurrentVolOffID(1, cpCry);			//Cry+Wrapping								Cry+Wrapping
    volIdAlv = gMC->CurrentVolOffID(2, cpAlv);			//Alveolus Inner							Alveolus_EC
    volIdSupAlv = gMC->CurrentVolOffID(3, cpSupAlv);//Alveolus Out								CalifaWorld
    // LOG(ERROR) << "TEST INITIAL. " <<  gMC->CurrentVolPath()<< FairLogger::endl;
		
    info.volIdAlv = volIdAlv;
    info.cpAlv = cpAlv;
    info.cpCry = cpCry;

    
    bool Print_Cryinfo=kFALSE;
    if(Print_Cryinfo)
    {
			cout<<">>> CALIFA main class <<<"<<endl;
		  cout<<"-------   Crystal info   ------"<<endl;	
		  cout<<"-- from gMC"<<endl;
		  cout<<"Path= "<<gMC->CurrentVolPath()<<endl;	
		  cout<<"1.Current Volume "<<endl;
		  cout<<"volId1_Name= "<<gMC->CurrentVolOffName(0)<<endl; //or bufferName
		  cout<<"volId1= "<<volId1<<endl;	  
		  cout<<"cp1= "<<cp1<<endl;
		  
		  cout<<"2."<<endl;
		  cout <<"volIdCry_Name= "<<gMC->CurrentVolOffName(1)<<endl;
		  cout<<"volIdCry= "<<volIdCry<<endl;
		  cout<<"cpCry= "<<cpCry<<endl;
		  
		  cout<<"3."<<endl;
		  cout <<"volIdAlv_Name= "<<gMC->CurrentVolOffName(2)<<endl;	
		  cout<<"volIdAlv= "<<volIdAlv<<endl;
		  cout<<"cpAlv= "<<cpAlv<<endl;
		  	 		  
		  cout<<"4."<<endl;
		  cout << "volIdSupAlv_Name= "<<gMC->CurrentVolOffName(3)<<endl;
		  cout<<"volIdSupAlv= "<<volIdSupAlv<<endl;
		  cout<<"cpSupAlv= "<<cpSupAlv<<endl;
			cout<<"----------"<<endl;
		}	

				/*NOTE: 1. gMC->VolName(id) is virtual (not defined) so not use it
								2. Use CurrentVolOffName without convert it into a const char, because it not works
						    			const char *Name = gMC->CurrentVolOffName(1) NO
						    			gMC->CurrentVolOffName(1) YES
				*/

        const char* alveolusECPrefix = "Alveolus_EC";
        const char* alveolusPrefix = "Alveolus_";
        const char* volumeNameCrystal = "";

        // Workaround to fix the hierarchy difference between Barrel and Endcap
        /*if (strncmp("CalifaWorld", gMC->CurrentVolOffName(3), 10) == 0)
        {
            volumeName = gMC->VolName(volIdAlv); //volumeName="Alveolus_Inner"   gMC->CurrentVolOffName(?);
            volumeNameCrystal = gMC->VolName(volId1); //volumeNameCrystal="Crystal" gMC->CurrentVolOffName(0);
        }*/        
        if (strncmp(alveolusECPrefix, gMC->CurrentVolOffName(2), 11) == 0) //ENDCAP GEOMETRY
        //NOTE:
        //  EndCap Geometry reproduces: Alveolus_EC/CrystalWithWrapping/Crystal
        //  if Alveolus_Inner is added to thE EndCap Geo change this gMC->CurrentVolOffName(3)
        { 
						//EndCap Volumes Info:		    
						//info.crystalType from 1 to 24
						//info.crystalCopy from 1 to 32
						//info.crystalId from 3000 to 3767
						
						/*
						 PATH:
							crystalType = ((iD - 3000) % 24) + 1;  from 1 to 24
							crystalCopy = (iD-3000 - crystalType + 1) / 24 + 1; from 1 to 32
							Int_t alveoliType[24]={1,1,2,2,3,3,4,4,5,6,7,8,9,9,10,10,11,11,12,12,13,13,14,14};
							
							Int_t wrappingType[24]={1,1,2,2,3,3,4,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
							

							nameVolume = TString::Format("/cave_1/CalifaWorld_0/Alveolus_EC_%i_%i/CrystalWithWrapping_%i_1/Crystal_%i_1",
								      alveoliType[crystalType-1], crystalCopy-1, wrappingType[crystalType-1], crystalType);
						*/
           
            info.crystalType = atoi(gMC->CurrentVolOffName(0) + 8); // converting to int the crystal index
            info.crystalCopy = cpAlv + 1;
            
            if (info.crystalType < 9)
            {
                // CC Phoswich crystals are combined in one crystal
                // Energies are contained in Nf -> LaBr and Ns -> LaCl

                // info.fEndcapIdentifier = 0 -> iPhos
                // info.fEndcapIdentifier = 1 -> Phoswich
                info.fEndcapIdentifier = 1;

                if (info.crystalType % 2 == 0)
                {
                    // info.fPhoswichIdentifier = 2 -> LaCl
                    info.fPhoswichIdentifier = 2;
                    //info.crystalType -= 1; 	//maybe this code line is wrong
                }
                else
                {
                    // info.fPhoswichIdentifier = 1 -> LaBr
                    info.fPhoswichIdentifier = 1;
                }
            }
                   		
            info.crystalId = 3000 + cpAlv * 24 + (info.crystalType - 1);
          
            if (info.crystalType > 24 || info.crystalType < 1 || info.crystalCopy > 32 || info.crystalCopy < 1 ||
                info.crystalId < 3000 || info.crystalId > 3767)
            {
                LOG(ERROR) << "R3BCalifa: Wrong crystal number in geometryVersion 16 (CC). " << FairLogger::endl;
                return kFALSE;
            }
        }
        else if (strncmp("Alveolus", gMC->CurrentVolOffName(3), 8) == 0)//BARREL GEOMETRY
        {
        
        	//Barrel Volumes Info:		
					//Crystal			    volId1=num.Id Geant    & cp1=1 copy num.
					//Cry+Wrapping    volIdCry1=num.Id Geant & cpCry=from 0 to 1 copy num.
					//Alveolus Inner  volIdAlv=num.Id Geant  & cpAlv=1 always copy num.
					//Alveolus Out    volIdSupAlv=num.Id Geant & cpSupAlv=from 0 to 31 copy num.								
        	//info.crystalType from 1 to 16
					//info.crystalCopy from 1 to 128
					//info.crystalId from 1 to 1952
					
					/*PATH:
									// First ring (single crystals)
									crystalType = 1;  //Alv type 1
									crystalCopy = iD;     //for Alv type 1 
									alveolusCopy = iD;    //Alv type 1 
									crystalInAlveolus =1;          //Alv type 1
								
									// Ring 2 - 16: 2x2 crystals
									crystalType = (Int_t)((iD-33)/128) + 2;  //Alv type (2 to 16)
									crystalCopy = ((iD-33)%128) + 1;         //CrystalCopy (1 to 128)
									alveolusCopy =(Int_t)(((iD-33)%128)/4) +1; //Alveolus copy (1 to 32)
									crystalInAlveolus = (iD-33)%4 + 1;//Crystal number in alveolus (1 to 4)
								
									Int_t alveoliType[16]={1,2,2,2,2,3,3,4,4,4,5,5,5,6,6,6};
				
									nameVolume = TString::Format( 
									"/cave_1/CalifaWorld_0/Alveolus_%i_%i/AlveolusInner_%i_1/CrystalWithWrapping_%i_%i_%i/Crystal_%i_%i_1",
									crystalType, alveolusCopy-1, 
									crystalType, alveoliType[crystalType-1], 
									crystalInAlveolus, crystalInAlveolus-1, 
									alveoliType[crystalType-1], crystalInAlveolus);
					*/
					
					
          info.crystalType = atoi(gMC->CurrentVolOffName(3) + 9); // converting to int the alveolus index
            
            if (info.crystalType == 1)
            {
                // only one crystal per alveoli in this ring, running from 1 to 32
                info.crystalCopy = cpSupAlv + 1;
                info.crystalId = cpSupAlv + 1;
            }
            else if (info.crystalType > 1 && info.crystalType < 17)
            {
                // running from 0*4+0+1=1 to 31*4+3+1=128
                info.crystalCopy = cpSupAlv * 4 + cpCry + 1;
                // running from 32+0*128+0*4+0+1=33 to 32+14*128+31*4+3+1=1952
                info.crystalId = 32 + (info.crystalType - 2) * 128 + cpSupAlv * 4 + cpCry + 1;
            }
            if (info.crystalType > 16 || info.crystalType < 1 || info.crystalCopy > 128 || info.crystalCopy < 1 ||
                info.crystalId > 1952 || info.crystalId < 1)
            {
                LOG(ERROR) << "R3BCalifa: Wrong crystal number in geometryVersion 16 (BARREL)." << FairLogger::endl;
                return kFALSE;
            }
        }
        else
        {
            LOG(ERROR) << "R3BCalifa: Impossible info.crystalType for geometryVersion 16." << FairLogger::endl;
            return kFALSE;
        }



						R3BCalifaGeometry* CalifaGeo;
						const char* path;					
 						Double_t polar;
						Double_t azimuthal;
						Double_t rho;
						Int_t id;

						//GetAngles ->done
						CalifaGeo->GetAngles(info.crystalId, polar, azimuthal, rho);
						//cout<<"CALIFA main: info.crystalId="<<info.crystalId<<"  polar="<<polar<<"  azimuthal="<<azimuthal<<" rho="<<rho<<endl<<endl;

						//GetPath ->done
						path=CalifaGeo->GetCrystalVolumePath(info.crystalId);
						//cout<<"CALIFA main: info.crystalId="<<info.crystalId<<endl;
						//cout<<"CALIFA main: path="<<path<<endl;

						//GetCrystalId ->not yet
						id=CalifaGeo->GetCrystalId(path);//NO va
						//cout<<">>> CALIFA main: path"<<gMC->CurrentVolPath()<<endl;
						//cout<<">>> CALIFA main: return id="<<id<<endl<<endl<<endl<<endl;
//double R3BCalifaGeometry::GetDistanceThroughCrystals(TVector3 &startVertex, TVector3 &direction, TVector3 *hitPos, int *numCrystals, int *crystalIds)

						//GetDistanceThroughCrystals
						TVector3 startVertex;
						TVector3 direction;
						TVector3 hitPos;
						Int_t *numCrystals;
						//*numCrystals=4; 
						Int_t *crystalIds;
						Double_t distance;
						startVertex={5,0,5};//need initialize these two vectors, otherwise not works at all
						direction={3,3,3};

						distance=CalifaGeo->GetDistanceThroughCrystals(startVertex,direction,&hitPos);
						//cout<<endl<<endl<<endl<<"distance="<<distance<<endl;

   return kTRUE;
 
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BCalifa::BeginEvent()
{

    //  if (! kGeoSaved ) {
    //      SaveGeoParams();
    //  cout << "-I STS geometry parameters saved " << endl;
    //  kGeoSaved = kTRUE;
    //  }
}
// -----   Public method EndOfEvent   -----------------------------------------
void R3BCalifa::EndOfEvent()
{
    if (fVerboseLevel)
        Print();

    fCaloCollection->Clear();

    ResetParameters();

    fCrystal = NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BCalifa::Register()
{
     FairRootManager::Instance()->Register("CrystalPoint", GetName(),
                                          fCaloCollection, kTRUE);
}
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BCalifa::GetCollection(Int_t iColl) const
{
    // HAPOL TODO -- DO I NEED TO RETURN A fCaloCrystalHitColletion????
    if (iColl == 0)
    {
        return fCaloCollection;
    }
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BCalifa::Print(Option_t* option) const
{
    Int_t nHits = fCaloCollection->GetEntriesFast();
    LOG(INFO) << "R3BCalifa: " << nHits << " points registered in this event" << FairLogger::endl;
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BCalifa::Reset()
{
    fCaloCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BCalifa::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(INFO) << "R3BCalifa: " << nEntries << " entries to add" << FairLogger::endl;
    TClonesArray& clref = *cl2;
    R3BCalifaPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = (R3BCalifaPoint*)cl1->At(i);
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BCalifaPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(INFO) << "R3BCalifa: " << cl2->GetEntriesFast() << " merged entries" << FairLogger::endl;
}

// -----   Private method AddPoint   --------------------------------------------
R3BCalifaPoint* R3BCalifa::AddPoint(Int_t trackID,
                                    Int_t detID,
                                    Int_t volid,
                                    Int_t copy,
                                    Int_t ident,
                                    TVector3 posIn,
                                    TVector3 posOut,
                                    TVector3 momIn,
                                    TVector3 momOut,
                                    Double_t time,
                                    Double_t length,
                                    Double_t eLoss,
                                    Double_t Nf,
                                    Double_t Ns)
{
    TClonesArray& clref = *fCaloCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
        LOG(INFO) << "R3BCalifa: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV"
                  << FairLogger::endl;
    return new (clref[size])
        R3BCalifaPoint(trackID, detID, volid, copy, ident, posIn, posOut, momIn, momOut, time, length, eLoss, Nf, Ns);
}

// -----  Public method SelectGeometryVersion  ----------------------------------
void R3BCalifa::SelectGeometryVersion(Int_t version) { fGeometryVersion = version; }

// -----  Public method SetNonUniformity  ----------------------------------
void R3BCalifa::SetNonUniformity(Double_t nonU)
{
    fNonUniformity = nonU;
    LOG(INFO) << "R3BCalifa::SetNonUniformity to " << fNonUniformity << " %" << FairLogger::endl;
}

Bool_t R3BCalifa::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("Crystal_"))
    {
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BCalifa)
