/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ------------------------------------------------------------
// -----             R3BCalifaJulichOnlineSpectra                 -----
// -----    Created 16/07/21  by J.L. Rodriguez-Sanchez   -----
// -----          Fill CalifaJulich online histograms             -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with CalifaJulich online data
 */

#include "R3BCalifaJulichOnlineSpectra.h"
#include "R3BAmsMappedData.h"
#include "R3BAmsStripCalData.h"
#include "R3BAmsHitData.h"
#include "R3BCalifaMappedData.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaHitData.h"

#include "R3BEventHeader.h"
#include "THttpServer.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TClonesArray.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

R3BCalifaJulichOnlineSpectra::R3BCalifaJulichOnlineSpectra()
    : FairTask("CalifaJulichOnlineSpectra", 1)
    , fMappedItemsCalifa(NULL)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
    , fMappedItemsSi(NULL)
    , fCalItemsSi(NULL)
    , fHitItemsSi(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(1)
    , fNbCrystals(16)
{
}

R3BCalifaJulichOnlineSpectra::R3BCalifaJulichOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsCalifa(NULL)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
    , fMappedItemsSi(NULL)
    , fCalItemsSi(NULL)
    , fHitItemsSi(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(1)
    , fNbCrystals(16)
{
}

R3BCalifaJulichOnlineSpectra::~R3BCalifaJulichOnlineSpectra()
{
    LOG(DEBUG) << "R3BCalifaJulichOnlineSpectra::Delete instance";
    if (fMappedItemsCalifa)
        delete fMappedItemsCalifa;
    if (fMappedItemsSi)
        delete fMappedItemsSi;
    if (fCalItemsCalifa)
        delete fCalItemsCalifa;
    if (fCalItemsSi)
        delete fCalItemsSi;
    if (fHitItemsSi)
        delete fHitItemsSi;
}

InitStatus R3BCalifaJulichOnlineSpectra::Init()
{
    LOG(INFO) << "R3BCalifaJulichOnlineSpectra::Init()";

    // Looking for FairRootManager
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectra::FairRootManager not found";

    // Create histograms for all the detectors

    // Get access to Mapped data
    fMappedItemsSi = (TClonesArray*)mgr->GetObject("AmsMappedData");
    if (!fMappedItemsSi)
    {
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectra::AmsMappedData not found";
        return kFATAL;
    }

    // Get access to Mapped data
    fMappedItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaMappedData");
    if (!fMappedItemsCalifa)
    {
        LOG(FATAL) << "R3BCalifaJulichOnlineSpectra::CalifaMappedData not found";
        return kFATAL;
    }

    // Get access to Cal data
    fCalItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaCrystalCalData");
    if (!fCalItemsCalifa)
        LOG(WARNING) << "R3BCalifaJulichOnlineSpectra::CalifaCrystalCalData not found";

    fCalItemsSi = (TClonesArray*)mgr->GetObject("AmsStripCalData");
    if (!fCalItemsSi)
        LOG(WARNING) << "R3BCalifaJulichOnlineSpectra::AmsStripCalData not found";

      fHitItemsSi = (TClonesArray*)mgr->GetObject("AmsHitData");
      if (!fHitItemsSi)
          LOG(WARNING) << "R3BCalifaJulichOnlineSpectra::AmsHitData not found";


    // Energy range for strips
    Double_t binsE = 200;
    Double_t minE = 0;
    Double_t maxE = 3500;
    char Name1[255];
    char Name2[255];
    char Name3[255];
    char Name4[255];

    // MAIN FOLDERs
    TFolder* mainfolCalifa = new TFolder("Califa", "Califa info");
    TFolder* mainfolSi = new TFolder("Si", "Si info");
    // Folders for mapped data
    TFolder* mapfolCalifa = new TFolder("Mapped", "Califa Mapped");
    TFolder* mapfolSi = new TFolder("MappedSi", "Si Mapped");
    // Folders for cal data
    TFolder* calfolCalifa = new TFolder("Cal", "Cal Califa info");
    TFolder* calfolSi = new TFolder("CalSi", "Cal Si info");
    // Folders for hit data and correlations
    TFolder* correlations = new TFolder("Correlations", "Hit Califa info");
    TFolder* hitfolSi = new TFolder("HitSi", "Hit Si info");

    // Mapped data Si
    fh2_EnergyVsStrip.resize(fNbDet);
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh2_energy_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Energy vs strip number for Si Det: %d", i + 1);
        fh2_EnergyVsStrip[i] = new TH2F(Name1, Name2, 64, 0, 64, binsE, minE, maxE);
        fh2_EnergyVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_EnergyVsStrip[i]->GetYaxis()->SetTitle("Energy [channels]");
        fh2_EnergyVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyVsStrip[i]->GetXaxis()->CenterTitle(true);
        fh2_EnergyVsStrip[i]->GetYaxis()->CenterTitle(true);
        fh2_EnergyVsStrip[i]->Draw("colz");
        mapfolSi->Add(fh2_EnergyVsStrip[i]);
    }
    mainfolSi->Add(mapfolSi);

    // Mapped data Califa
    fh1_MultiplicityGamma = new TH1F("MultGamma", "MultGamma", 40, 0, 20);
    fh1_MultiplicityGamma->GetXaxis()->SetTitle("Multiplicity Gamma");
    fh1_MultiplicityGamma->GetYaxis()->SetTitleOffset(1.4);
    fh1_MultiplicityGamma->GetXaxis()->CenterTitle(true);
    fh1_MultiplicityGamma->GetYaxis()->CenterTitle(true);
    fh1_MultiplicityGamma->Draw("col");
    mapfolCalifa->Add(fh1_MultiplicityGamma);

    fh1_MultiplicityProton = new TH1F("MultProton", "MultProton", 40, 0, 20);
    fh1_MultiplicityProton->GetXaxis()->SetTitle("Multiplicity Proton");
    fh1_MultiplicityProton->GetYaxis()->SetTitleOffset(1.4);
    fh1_MultiplicityProton->GetXaxis()->CenterTitle(true);
    fh1_MultiplicityProton->GetYaxis()->CenterTitle(true);
    fh1_MultiplicityProton->Draw("col");
    mapfolCalifa->Add(fh1_MultiplicityProton);

    fh1_EnergyCalifaCrystals.resize(fNbCrystals);
    fh2_Map_nf_ns.resize(fNbCrystals);
    for (Int_t i=0; i<fNbCrystals; i++)
    {
      if (i==0)  {sprintf(Name1, "fh1_energyBoxA_1_p"); sprintf(Name3, "fh2_Map_nf_ns_BoxA_1_p"); }
      if (i==1)  {sprintf(Name1, "fh1_energyBoxA_2_p"); sprintf(Name3, "fh2_Map_nf_ns_BoxA_2_p"); }
      if (i==2)  {sprintf(Name1, "fh1_energyBoxA_3_p"); sprintf(Name3, "fh2_Map_nf_ns_BoxA_3_p"); }
      if (i==3)  {sprintf(Name1, "fh1_energyBoxA_4_p"); sprintf(Name3, "fh2_Map_nf_ns_BoxA_4_p"); }
      if (i==4)  {sprintf(Name1, "fh1_energyBoxA_1_g"); sprintf(Name3, "fh2_Map_nf_ns_BoxA_1_g"); }
      if (i==5)  {sprintf(Name1, "fh1_energyBoxA_2_g"); sprintf(Name3, "fh2_Map_nf_ns_BoxA_2_g"); }
      if (i==6)  {sprintf(Name1, "fh1_energyBoxA_3_g"); sprintf(Name3, "fh2_Map_nf_ns_BoxA_3_g"); }
      if (i==7)  {sprintf(Name1, "fh1_energyBoxA_4_g"); sprintf(Name3, "fh2_Map_nf_ns_BoxA_4_g"); }
      if (i==8)  {sprintf(Name1, "fh1_energyBoxB_1_p"); sprintf(Name3, "fh2_Map_nf_ns_BoxB_1_p"); }
      if (i==9)  {sprintf(Name1, "fh1_energyBoxB_2_p"); sprintf(Name3, "fh2_Map_nf_ns_BoxB_2_p"); }
      if (i==10) {sprintf(Name1, "fh1_energyBoxB_3_p"); sprintf(Name3, "fh2_Map_nf_ns_BoxB_3_p"); }
      if (i==11) {sprintf(Name1, "fh1_energyBoxB_4_p"); sprintf(Name3, "fh2_Map_nf_ns_BoxB_4_p"); }
      if (i==12) {sprintf(Name1, "fh1_energyBoxB_1_g"); sprintf(Name3, "fh2_Map_nf_ns_BoxB_1_g"); }
      if (i==13) {sprintf(Name1, "fh1_energyBoxB_2_g"); sprintf(Name3, "fh2_Map_nf_ns_BoxB_2_g"); }
      if (i==14) {sprintf(Name1, "fh1_energyBoxB_3_g"); sprintf(Name3, "fh2_Map_nf_ns_BoxB_3_g"); }
      if (i==15) {sprintf(Name1, "fh1_energyBoxB_4_g"); sprintf(Name3, "fh2_Map_nf_ns_BoxB_4_g"); }

      sprintf(Name2, "Energy in Califa crystal: %d", i + 1);
      fh1_EnergyCalifaCrystals[i] = new TH1F(Name1, Name2, 3000, 0, 30000);
      fh1_EnergyCalifaCrystals[i]->GetXaxis()->SetTitle("Energy [channels]");
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->SetTitle("counts");
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh1_EnergyCalifaCrystals[i]->GetXaxis()->CenterTitle(true);
      fh1_EnergyCalifaCrystals[i]->GetYaxis()->CenterTitle(true);
      fh1_EnergyCalifaCrystals[i]->Draw("col");
      mapfolCalifa->Add(fh1_EnergyCalifaCrystals[i]);

      sprintf(Name4, "Ns vs Nf: %d", i + 1);
      fh2_Map_nf_ns[i] = new TH2F(Name3, Name4, 1000, -1000, 20000, 1000, -1000, 20000);
      fh2_Map_nf_ns[i]->GetXaxis()->SetTitle("nf");
      fh2_Map_nf_ns[i]->GetYaxis()->SetTitle("ns");
      fh2_Map_nf_ns[i]->GetYaxis()->SetTitleOffset(1.4);
      fh2_Map_nf_ns[i]->GetXaxis()->CenterTitle(true);
      fh2_Map_nf_ns[i]->GetYaxis()->CenterTitle(true);
      fh2_Map_nf_ns[i]->Draw("col");
      mapfolCalifa->Add(fh2_Map_nf_ns[i]);
    }
    mainfolCalifa->Add(mapfolCalifa);

    for (Int_t i=0;i<4;i++)
    {
      sprintf(Name1, "fh2_EnergyCrystal%d", i+1);
      sprintf(Name2, "Energy in Crystal %d", i+1);
      fh2_EnergyMapCalifa_SiStrip[i] = new TH2F(Name1, Name2, 64, 0, 64, 1000, 0, 300000);
      fh2_EnergyMapCalifa_SiStrip[i]->GetXaxis()->SetTitle("Strip Number x[0-31],y[32-64]");
      fh2_EnergyMapCalifa_SiStrip[i]->GetYaxis()->SetTitle("Energy");
      fh2_EnergyMapCalifa_SiStrip[i]->GetYaxis()->SetTitleOffset(1.4);
      fh2_EnergyMapCalifa_SiStrip[i]->GetXaxis()->CenterTitle(true);
      fh2_EnergyMapCalifa_SiStrip[i]->GetYaxis()->CenterTitle(true);
      fh2_EnergyMapCalifa_SiStrip[i]->Draw("");
      mapfolSi->Add(fh2_EnergyMapCalifa_SiStrip[i]);

      // sprintf(Name3, "fh2_PosX_PosY_crystal_%d", 4-i);
      // sprintf(Name4, "fh2_PosX_PosY crystal %d", 4-i);
      sprintf(Name3, "fh2_PosX_PosY_crystal_%d", i+1);
      sprintf(Name4, "fh2_PosX_PosY crystal %d", i+1);
      fh2_PosX_PosY_Califa[i] = new TH2F(Name3, Name4, 34, 0, 34, 34, 0, 34);
      fh2_PosX_PosY_Califa[i]->GetXaxis()->SetTitle("Strip x");
      fh2_PosX_PosY_Califa[i]->GetYaxis()->SetTitle("Strip Y");
      fh2_PosX_PosY_Califa[i]->GetYaxis()->SetTitleOffset(1.4);
      fh2_PosX_PosY_Califa[i]->GetXaxis()->CenterTitle(true);
      fh2_PosX_PosY_Califa[i]->GetYaxis()->CenterTitle(true);
      fh2_PosX_PosY_Califa[i]->Draw("");
      mapfolSi->Add(fh2_PosX_PosY_Califa[i]);
    }

    fh2_EnergyTotMapCalifa_SiStrip = new TH2F("fh2_EnergyTotMapCalifa_SiStrip", "fh2_EnergyTotMapCalifa_SiStrip", 64, 0, 64, 1000, 0, 300000);
    fh2_EnergyTotMapCalifa_SiStrip->GetXaxis()->SetTitle("Strip Number x[0-31],y[32-64]");
    fh2_EnergyTotMapCalifa_SiStrip->GetYaxis()->SetTitle("Energy");
    fh2_EnergyTotMapCalifa_SiStrip->GetYaxis()->SetTitleOffset(1.4);
    fh2_EnergyTotMapCalifa_SiStrip->GetXaxis()->CenterTitle(true);
    fh2_EnergyTotMapCalifa_SiStrip->GetYaxis()->CenterTitle(true);
    fh2_EnergyTotMapCalifa_SiStrip->Draw("");
    mapfolSi->Add(fh2_EnergyTotMapCalifa_SiStrip);

    // Cal data Si
    fh2_EnergyCalVsStrip.resize(fNbDet);
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh2_energyCal_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Energy calibrated vs strip number for Si Det: %d", i + 1);
        fh2_EnergyCalVsStrip[i] = new TH2F(Name1, Name2, 64, 0, 64, binsE, minE, maxE);
        fh2_EnergyCalVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_EnergyCalVsStrip[i]->GetYaxis()->SetTitle("Energy [keV]");
        fh2_EnergyCalVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyCalVsStrip[i]->GetXaxis()->CenterTitle(true);
        fh2_EnergyCalVsStrip[i]->GetYaxis()->CenterTitle(true);
        fh2_EnergyCalVsStrip[i]->Draw("col");
        mapfolSi->Add(fh2_EnergyCalVsStrip[i]);
    }
    mainfolSi->Add(calfolSi);

    // Cal data CALIFA
    fh1_EnergyCalCalifaCrystals.resize(fNbCrystals);
    for (Int_t i=0; i<fNbCrystals; i++)
    {
      sprintf(Name2, "Energy calibrated in Califa crystal: %d", i + 1);
      Int_t pad = 0;
      if (i<4)
      {
        if (i==0) {sprintf(Name1,  "fh1_energyCalBoxA_1_p"); pad=2;}
        if (i==1) {sprintf(Name1,  "fh1_energyCalBoxA_2_p"); pad=4;}
        if (i==2) {sprintf(Name1,  "fh1_energyCalBoxA_3_p"); pad=3;}
        if (i==3) {sprintf(Name1,  "fh1_energyCalBoxA_4_p"); pad=1;}
        fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 30000, 0, 300000);
      }

      if (i<8 && i>3)
      {
        if (i==4) {sprintf(Name1,  "fh1_energyCalBoxA_1_g"); pad=2; }
        if (i==5) {sprintf(Name1,  "fh1_energyCalBoxA_2_g"); pad=4; }
        if (i==6) {sprintf(Name1,  "fh1_energyCalBoxA_3_g"); pad=3; }
        if (i==7) {sprintf(Name1,  "fh1_energyCalBoxA_4_g"); pad=1; }
        fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 3000, 0, 30000);
      }

      if (i<12 && i>7)
      {
        if (i==8) {sprintf(Name1,  "fh1_energyCalBoxB_1_p"); pad=2; }
        if (i==9) {sprintf(Name1,  "fh1_energyCalBoxB_2_p"); pad=4; }
        if (i==10) {sprintf(Name1, "fh1_energyCalBoxB_3_p"); pad=3; }
        if (i==11) {sprintf(Name1, "fh1_energyCalBoxB_4_p"); pad=1; }
        fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 30000, 0, 300000);
      }

      if (i>11)
      {
        if (i==12) {sprintf(Name1, "fh1_energyCalBoxB_1_g"); pad=2;}
        if (i==13) {sprintf(Name1, "fh1_energyCalBoxB_2_g"); pad=4;}
        if (i==14) {sprintf(Name1, "fh1_energyCalBoxB_3_g"); pad=3;}
        if (i==15) {sprintf(Name1, "fh1_energyCalBoxB_4_g"); pad=1;}
       fh1_EnergyCalCalifaCrystals[i] = new TH1F(Name1, Name2, 3000, 0, 30000);
       }

      fh1_EnergyCalCalifaCrystals[i]->GetXaxis()->SetTitle("Energy [keV]");
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->SetTitle("counts");
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh1_EnergyCalCalifaCrystals[i]->GetXaxis()->CenterTitle(true);
      fh1_EnergyCalCalifaCrystals[i]->GetYaxis()->CenterTitle(true);
      fh1_EnergyCalCalifaCrystals[i]->Draw();
      calfolCalifa->Add(fh1_EnergyCalCalifaCrystals[i]);
    }
    mainfolCalifa->Add(calfolCalifa);


    // energy correlations CALIFA
    fh2_EnergyCorrelationsCrystals.resize(24);
    for (Int_t i=0; i<24; i++)
    {
      //BoxA
      if (i==0) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_2_p"); }
      if (i==1) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_3_p"); }
      if (i==2) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_4_p"); }
      if (i==3) {sprintf(Name1, "fh2_EnergyCorr_BoxA_2_3_p"); }
      if (i==4) {sprintf(Name1, "fh2_EnergyCorr_BoxA_2_4_p"); }
      if (i==5) {sprintf(Name1, "fh2_EnergyCorr_BoxA_3_4_p"); }
      if (i==6) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_2_g"); }
      if (i==7) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_3_g"); }
      if (i==8) {sprintf(Name1, "fh2_EnergyCorr_BoxA_1_4_g"); }
      if (i==9) {sprintf(Name1, "fh2_EnergyCorr_BoxA_2_4_g"); }
      if (i==10) {sprintf(Name1, "fh2_EnergyCorr_BoxA_3_4_g"); }
      if (i==11) {sprintf(Name1, "fh2_EnergyCorr_BoxA_2_3_g"); }
      //BoxB
      if (i==12) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_2_p"); }
      if (i==13) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_3_p"); }
      if (i==14) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_4_p"); }
      if (i==15) {sprintf(Name1, "fh2_EnergyCorr_BoxB_2_3_p"); }
      if (i==16) {sprintf(Name1, "fh2_EnergyCorr_BoxB_2_4_p"); }
      if (i==17) {sprintf(Name1, "fh2_EnergyCorr_BoxB_3_4_p"); }
      if (i==18) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_2_g"); }
      if (i==19) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_3_g"); }
      if (i==20) {sprintf(Name1, "fh2_EnergyCorr_BoxB_1_4_g"); }
      if (i==21) {sprintf(Name1, "fh2_EnergyCorr_BoxB_2_4_g"); }
      if (i==22) {sprintf(Name1, "fh2_EnergyCorr_BoxB_3_4_g"); }
      if (i==23) {sprintf(Name1, "fh2_EnergyCorr_BoxB_2_3_g"); }

      sprintf(Name2, "fh2_EnergyCorrelationsCrystals_%d", i);
      if (i<5 || (i>11 && i<18)) {fh2_EnergyCorrelationsCrystals[i] = new TH2F(Name1, Name2, 3000, 0, 300000, 3000,0,300000);}
      else {fh2_EnergyCorrelationsCrystals[i] = new TH2F(Name1, Name2, 3000, 0, 30000, 3000,0,30000);}

      fh2_EnergyCorrelationsCrystals[i]->GetYaxis()->SetTitleOffset(1.4);
      fh2_EnergyCorrelationsCrystals[i]->GetXaxis()->CenterTitle(true);
      fh2_EnergyCorrelationsCrystals[i]->GetYaxis()->CenterTitle(true);
      fh2_EnergyCorrelationsCrystals[i]->Draw("COLZ");
      correlations->Add(fh2_EnergyCorrelationsCrystals[i]);
    }

    fh2_EnergyCorrelationsAlvProton = new TH2F("fh2_EnergyCorrelationsAlvProton", "Energies Alv ProtonRange", 3000, 0, 300000, 3000,0,300000);
    fh2_EnergyCorrelationsAlvProton->GetXaxis()->SetTitle("Energy Alv1");
    fh2_EnergyCorrelationsAlvProton->GetYaxis()->SetTitle("Energy Alv2");
    fh2_EnergyCorrelationsAlvProton->GetYaxis()->SetTitleOffset(1.4);
    fh2_EnergyCorrelationsAlvProton->GetXaxis()->CenterTitle(true);
    fh2_EnergyCorrelationsAlvProton->GetYaxis()->CenterTitle(true);
    fh2_EnergyCorrelationsAlvProton->Draw("col");
    correlations->Add(fh2_EnergyCorrelationsAlvProton);

    fh2_EnergyCorrelationsAlvGamma = new TH2F("fh2_EnergyCorrelationsAlvGamma", "Energies Alv GammaRange", 300, 0, 30000, 300,0,30000);
    fh2_EnergyCorrelationsAlvGamma->GetXaxis()->SetTitle("Energy Alv1");
    fh2_EnergyCorrelationsAlvGamma->GetYaxis()->SetTitle("Energy Alv2");
    fh2_EnergyCorrelationsAlvGamma->GetYaxis()->SetTitleOffset(1.4);
    fh2_EnergyCorrelationsAlvGamma->GetXaxis()->CenterTitle(true);
    fh2_EnergyCorrelationsAlvGamma->GetYaxis()->CenterTitle(true);
    fh2_EnergyCorrelationsAlvGamma->Draw("COLZ");
    correlations->Add(fh2_EnergyCorrelationsAlvGamma);

    mainfolCalifa->Add(correlations);


    fh1_EnergyTotBoxA_g = new TH1F("fh1_EnergyTotBoxA_g", "fh1_EnergyTotBoxA_g",3000,0,30000);
    fh1_EnergyTotBoxA_g->Draw("col");
    fh1_EnergyTotBoxA_p = new TH1F("fh1_EnergyTotBoxA_p", "fh1_EnergyTotBoxA_p",3000,0,300000);
    fh1_EnergyTotBoxA_p->Draw("col");
    fh1_EnergyTotBoxB_g = new TH1F("fh1_EnergyTotBoxB_g", "fh1_EnergyTotBoxB_g",3000,0,30000);
    fh1_EnergyTotBoxB_g->Draw("col");
    fh1_EnergyTotBoxB_p = new TH1F("fh1_EnergyTotBoxB_p", "fh1_EnergyTotBoxB_p",3000,0,300000);
    fh1_EnergyTotBoxB_p->Draw("col");


    fh2_Energy_1A_1B = new TH2F("fh2_Energy_1A_1B", "fh2_Energy_1A_1B", 3000, 0, 300000, 3000, 0, 300000);
    fh2_Energy_1A_1B->GetXaxis()->SetTitle("Energy 1A");
    fh2_Energy_1A_1B->GetYaxis()->SetTitle("Energy 1B");
    fh2_Energy_1A_1B->Draw();
    correlations->Add(fh2_Energy_1A_1B);

    calfolCalifa->Add(fh1_EnergyTotBoxA_g);
    calfolCalifa->Add(fh1_EnergyTotBoxA_p);
    calfolCalifa->Add(fh1_EnergyTotBoxB_g);
    calfolCalifa->Add(fh1_EnergyTotBoxB_p);

    // Hit data Si
    fh2_PosX_PosY.resize(fNbDet);
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh2_PosX_PosY");
        sprintf(Name2, "fh2_PosX_PosY ");
        fh2_PosX_PosY[i] = new TH2F(Name1, Name2, 34, 0, 34, 34, 0, 34);
        fh2_PosX_PosY[i]->GetXaxis()->SetTitle("Strip X ");
        fh2_PosX_PosY[i]->GetYaxis()->SetTitle("Strip Y ");
        fh2_PosX_PosY[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_PosX_PosY[i]->GetXaxis()->CenterTitle(true);
        fh2_PosX_PosY[i]->GetYaxis()->CenterTitle(true);
        fh2_PosX_PosY[i]->Draw("col");
        mapfolSi->Add(fh2_PosX_PosY[i]);
    }
    mainfolSi->Add(hitfolSi);

    // Looking for FairRunOnline
    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);
    run->AddObject(mainfolSi);
    run->AddObject(mainfolCalifa);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_CalifaJulich", Form("/Objects/%s/->Reset_CalifaJulich_Histo()", GetName()));

    return kSUCCESS;
}

void R3BCalifaJulichOnlineSpectra::Reset_CalifaJulich_Histo()
{
    LOG(INFO) << "R3BCalifaJulichOnlineSpectra::Reset_CalifaJulich_Histo";

    // Mapped data
    fh1_MultiplicityGamma->Reset();
    fh1_MultiplicityProton->Reset();
    fh2_Energy_1A_1B->Reset();
    for (Int_t i=0;i<4;i++) {fh2_EnergyMapCalifa_SiStrip[i]->Reset(); fh2_PosX_PosY_Califa[i]->Reset();}

    for (Int_t i = 0; i < fNbDet; i++)
    {
        fh2_EnergyVsStrip[i]->Reset();
        fh2_EnergyCalVsStrip[i]->Reset();
        fh2_PosX_PosY[i]->Reset();
    }

    for (Int_t i = 0; i < fNbCrystals; i++)
    {
        fh1_EnergyCalifaCrystals[i]->Reset();
        fh1_EnergyCalCalifaCrystals[i]->Reset();
        fh1_EnergyHitCalifaCrystals[i]->Reset();
        fh2_Map_nf_ns[i]->Reset();
    }

    for (Int_t i=0;i<24;i++) { fh2_EnergyCorrelationsCrystals[i] -> Reset(); }
}

void R3BCalifaJulichOnlineSpectra::Exec(Option_t* option)
{
    // Fill mapped data

    if (fMappedItemsSi && fMappedItemsSi->GetEntriesFast() > 0)
    {
        auto nHits = fMappedItemsSi->GetEntriesFast();
        // auto nHitsSi2 = fMappedItemsSi->GetEntriesFast();
        // Float_t x=0.; Float_t y=0.;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            Int_t bin = 0;
            R3BAmsMappedData* hit = (R3BAmsMappedData*)fMappedItemsSi->At(ihit);
            if (!hit)
                continue;
            fh2_EnergyVsStrip[hit->GetDetectorId()]->Fill(hit->GetStripId(), hit->GetEnergy());

        }
    }
    if (fMappedItemsSi && fMappedItemsSi->GetEntriesFast() > 0 && fMappedItemsCalifa && fMappedItemsCalifa->GetEntriesFast() > 0)
    {
        auto nHitsSi = fMappedItemsSi->GetEntriesFast();
        auto nHitsCalifa = fMappedItemsCalifa->GetEntriesFast();
        auto nHitsSi2 = fMappedItemsSi->GetEntriesFast();
        Float_t x=0.; Float_t y=0.;
        Float_t energyTot = 0;
        Int_t crystals[4]; Float_t energies[4];
        for (Int_t i=0;i<4;i++) {crystals[i]=0; energies[i]=0;}

        for (Int_t ihit = 0; ihit < nHitsSi; ihit++)
        {
          Int_t bin=100; Int_t bin2=100;
            R3BAmsMappedData* hitSi = (R3BAmsMappedData*)fMappedItemsSi->At(ihit);
            if (!hitSi || hitSi->GetEnergy()<100)
                continue;

            if (hitSi->GetStripId()==2)  {bin=7 ;}
            if (hitSi->GetStripId()==3)  {bin=6 ;}
            if (hitSi->GetStripId()==4)  {bin=5 ;}
            if (hitSi->GetStripId()==5)  {bin=4 ;}
            if (hitSi->GetStripId()==6)  {bin=3 ;}
            if (hitSi->GetStripId()==7)  {bin=2;}
            if (hitSi->GetStripId()==8)  {bin=9 ;}
            if (hitSi->GetStripId()==9)  {bin=10;}
            if (hitSi->GetStripId()==10) {bin=16;}
            if (hitSi->GetStripId()==1)  {bin=8;}
            if (hitSi->GetStripId()==11) {bin=15;}
            if (hitSi->GetStripId()==12) {bin=14;}
            if (hitSi->GetStripId()==13) {bin=13;}
            if (hitSi->GetStripId()==14) {bin=12;}
            if (hitSi->GetStripId()==15) {bin=11;}
            if (hitSi->GetStripId()==16) {bin=1;}
            if (hitSi->GetStripId()==17) {bin=32;}
            if (hitSi->GetStripId()==18) {bin=22;}
            if (hitSi->GetStripId()==19) {bin=21;}
            if (hitSi->GetStripId()==20) {bin=20;}
            if (hitSi->GetStripId()==21) {bin=19;}
            if (hitSi->GetStripId()==22) {bin=18;}
            if (hitSi->GetStripId()==32) {bin=25;}
            if (hitSi->GetStripId()==23) {bin=17;}
            if (hitSi->GetStripId()==24) {bin=23;}
            if (hitSi->GetStripId()==25) {bin=24;}
            if (hitSi->GetStripId()==26) {bin=31;}
            if (hitSi->GetStripId()==27) {bin=30;}
            if (hitSi->GetStripId()==28) {bin=29;}
            if (hitSi->GetStripId()==29) {bin=28;}
            if (hitSi->GetStripId()==30) {bin=27;}
            if (hitSi->GetStripId()==31) {bin=26;}

            for (Int_t ihitCal = 0; ihitCal < nHitsCalifa; ihitCal++)
            {
              R3BCalifaMappedData* hit = (R3BCalifaMappedData*)fMappedItemsCalifa->At(ihitCal);
              if (!hit || hit->GetEnergy()<5000)
                  continue;
              Int_t index=1000;
              if (hit->GetCrystalId()==324) {index=3;}
              if (hit->GetCrystalId()==325) {index=2;}
              if (hit->GetCrystalId()==326) {index=1;}
              if (hit->GetCrystalId()==327) {index=0;}

              if(index==1000) {return;}

              if (hit->GetCrystalId()==324 || hit->GetCrystalId()==325 || hit->GetCrystalId()==326 || hit->GetCrystalId()==327 )
              {
                energyTot =+ hit->GetEnergy();
              }

              crystals[index]=1;
              energies[index]=hit->GetEnergy();

              fh2_EnergyMapCalifa_SiStrip[index]->Fill(hitSi->GetStripId(),hit->GetEnergy());
              fh2_EnergyTotMapCalifa_SiStrip->Fill(hitSi->GetStripId(), energyTot);

              if (bin<100)
              {
                for (Int_t ihit2 = 0; ihit2 < nHitsSi2; ihit2++)
                {
                  R3BAmsMappedData* hitSi2 = (R3BAmsMappedData*)fMappedItemsSi->At(ihit2);
                  if (!hitSi2)
                      continue;

                  if (hitSi2->GetStripId()==2 +32) {bin2=7 ;}
                  if (hitSi2->GetStripId()==3 +32) {bin2=6 ;}
                  if (hitSi2->GetStripId()==4 +32) {bin2=5 ;}
                  if (hitSi2->GetStripId()==5 +32) {bin2=4 ;}
                  if (hitSi2->GetStripId()==6 +32) {bin2=3 ;}
                  if (hitSi2->GetStripId()==7 +32) {bin2=2 ;}
                  if (hitSi2->GetStripId()==8 +32) {bin2=9 ;}
                  if (hitSi2->GetStripId()==9 +32) {bin2=10;}
                  if (hitSi2->GetStripId()==10+32) {bin2=16;}
                  if (hitSi2->GetStripId()==1 +32) {bin2=8 ;}
                  if (hitSi2->GetStripId()==11+32) {bin2=15;}
                  if (hitSi2->GetStripId()==12+32) {bin2=14;}
                  if (hitSi2->GetStripId()==13+32) {bin2=13;}
                  if (hitSi2->GetStripId()==14+32) {bin2=12;}
                  if (hitSi2->GetStripId()==15+32) {bin2=11;}
                  if (hitSi2->GetStripId()==16+32) {bin2=1 ;}
                  if (hitSi2->GetStripId()==17+32) {bin2=32;}
                  if (hitSi2->GetStripId()==18+32) {bin2=22;}
                  if (hitSi2->GetStripId()==19+32) {bin2=21;}
                  if (hitSi2->GetStripId()==20+32) {bin2=20;}
                  if (hitSi2->GetStripId()==21+32) {bin2=19;}
                  if (hitSi2->GetStripId()==22+32) {bin2=18;}
                  if (hitSi2->GetStripId()==32+32) {bin2=25;}
                  if (hitSi2->GetStripId()==23+32) {bin2=17;}
                  if (hitSi2->GetStripId()==24+32) {bin2=23;}
                  if (hitSi2->GetStripId()==25+32) {bin2=24;}
                  if (hitSi2->GetStripId()==26+32) {bin2=31;}
                  if (hitSi2->GetStripId()==27+32) {bin2=30;}
                  if (hitSi2->GetStripId()==28+32) {bin2=29;}
                  if (hitSi2->GetStripId()==29+32) {bin2=28;}
                  if (hitSi2->GetStripId()==30+32) {bin2=27;}
                  if (hitSi2->GetStripId()==31+32) {bin2=26;}

                  if (bin2<100)
                  {
                    x = bin; y = bin2;
                    fh2_PosX_PosY[hitSi2->GetDetectorId()]->Fill(x, y);

                    if (hit->GetCrystalId()==327 && energies[1]<4000 && energies[2]<4000 && energies[3]<4000 ) {fh2_PosX_PosY_Califa[0]->Fill(x,y);}
                    if (hit->GetCrystalId()==326 && energies[0]<4000 && energies[2]<4000 && energies[3]<4000 ) {fh2_PosX_PosY_Califa[1]->Fill(x,y);}
                    if (hit->GetCrystalId()==325 && energies[0]<4000 && energies[1]<4000 && energies[3]<4000 ) {fh2_PosX_PosY_Califa[2]->Fill(x,y);}
                    if (hit->GetCrystalId()==324 && energies[0]<4000 && energies[1]<4000 && energies[2]<4000 ) {fh2_PosX_PosY_Califa[3]->Fill(x,y);}
                  }
                }
              }
            }
        }
    }


    if (fMappedItemsCalifa && fMappedItemsCalifa->GetEntriesFast() > 0)
    {
        auto nHits = fMappedItemsCalifa->GetEntriesFast();
        Int_t multGamma=0; Int_t multProton=0;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaMappedData* hit = (R3BCalifaMappedData*)fMappedItemsCalifa->At(ihit);
            if (!hit)
                continue;
            Int_t index=1000;

            if (hit->GetCrystalId()==324) {index=3;}
            if (hit->GetCrystalId()==325) {index=2;}
            if (hit->GetCrystalId()==326) {index=1;}
            if (hit->GetCrystalId()==327) {index=0;}

            if (hit->GetCrystalId()==340) {index=7;}
            if (hit->GetCrystalId()==341) {index=6;}
            if (hit->GetCrystalId()==342) {index=5;}
            if (hit->GetCrystalId()==343) {index=4;}

            if (hit->GetCrystalId()==356) {index=11;}
            if (hit->GetCrystalId()==357) {index=10;}
            if (hit->GetCrystalId()==358) {index=9;}
            if (hit->GetCrystalId()==359) {index=8;}

            if (hit->GetCrystalId()==372) {index=15;}
            if (hit->GetCrystalId()==373) {index=14;}
            if (hit->GetCrystalId()==374) {index=13;}
            if (hit->GetCrystalId()==375) {index=12;}

            if(index==1000) {return;}

            if (index<4 || (index>7 && index<12)) {multProton++;}
            else {multGamma++;}

            fh1_EnergyCalifaCrystals[index]->Fill(hit->GetEnergy());
            fh2_Map_nf_ns[index]->Fill(hit->GetNf(),hit->GetNs());
        }

        fh1_MultiplicityGamma->Fill(multGamma);
        fh1_MultiplicityProton->Fill(multProton);
    }

    if (fCalItemsCalifa && fCalItemsCalifa->GetEntriesFast() > 0)
    {
        auto nHits = fCalItemsCalifa->GetEntriesFast();
        Int_t crystals[16];
        Float_t energies[16];
        Float_t energy_alv[4]; for (Int_t i=0;i<4;i++) {energy_alv[i]=0.0;}
        for (Int_t i=0;i<16;i++) {crystals[i]=0; energies[i]=0.0;}

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaCrystalCalData* hit = (R3BCalifaCrystalCalData*)fCalItemsCalifa->At(ihit);
            if (!hit)
                continue;
            Int_t index=1000;

            if (hit->GetCrystalId()==324) {index=3;}
            if (hit->GetCrystalId()==325) {index=2;}
            if (hit->GetCrystalId()==326) {index=1;}
            if (hit->GetCrystalId()==327) {index=0;}

            if (hit->GetCrystalId()==340) {index=7;}
            if (hit->GetCrystalId()==341) {index=6;}
            if (hit->GetCrystalId()==342) {index=5;}
            if (hit->GetCrystalId()==343) {index=4;}

            if (hit->GetCrystalId()==356) {index=11;}
            if (hit->GetCrystalId()==357) {index=10;}
            if (hit->GetCrystalId()==358) {index=9;}
            if (hit->GetCrystalId()==359) {index=8;}

            if (hit->GetCrystalId()==372) {index=15;}
            if (hit->GetCrystalId()==373) {index=14;}
            if (hit->GetCrystalId()==374) {index=13;}
            if (hit->GetCrystalId()==375) {index=12;}

            if (index==1000) {return;}

            fh1_EnergyCalCalifaCrystals[index]->Fill(hit->GetEnergy());

            crystals[index]=1;
            energies[index]=hit->GetEnergy();
        }

        for (Int_t i=0;i<4;i++)
        {
          energy_alv[0] += energies[i]; energy_alv[1] += energies[i+4];
          energy_alv[2] += energies[i+8]; energy_alv[3] += energies[i+12];
        }

        // BoxA ProtonRange
        if (crystals[0] && crystals[1])
        {
          fh2_EnergyCorrelationsCrystals[0]->Fill(energies[0],energies[1]);
          fh2_EnergyCorrelationsCrystals[0]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[0]->GetYaxis()->SetTitle("Energy [keV] crystal 2");
        }

        if (crystals[0] && crystals[2])
        {
          fh2_EnergyCorrelationsCrystals[1]->Fill(energies[0],energies[2]);
          fh2_EnergyCorrelationsCrystals[1]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[1]->GetYaxis()->SetTitle("Energy [keV] crystal 3");
        }

        if (crystals[0] && crystals[3])
        {
          fh2_EnergyCorrelationsCrystals[2]->Fill(energies[0],energies[3]);
          fh2_EnergyCorrelationsCrystals[2]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[2]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        if (crystals[1] && crystals[2])
        {
          fh2_EnergyCorrelationsCrystals[3]->Fill(energies[1],energies[2]);
          fh2_EnergyCorrelationsCrystals[3]->GetXaxis()->SetTitle("Energy [keV] crystal 2");
          fh2_EnergyCorrelationsCrystals[3]->GetYaxis()->SetTitle("Energy [keV] crystal 3");
        }

        if (crystals[1] && crystals[3])
        {
          fh2_EnergyCorrelationsCrystals[4]->Fill(energies[1],energies[3]);
          fh2_EnergyCorrelationsCrystals[4]->GetXaxis()->SetTitle("Energy [keV] crystal 2");
          fh2_EnergyCorrelationsCrystals[4]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        if (crystals[2] && crystals[3])
        {
          fh2_EnergyCorrelationsCrystals[5]->Fill(energies[2],energies[3]);
          fh2_EnergyCorrelationsCrystals[5]->GetXaxis()->SetTitle("Energy [keV] crystal 3");
          fh2_EnergyCorrelationsCrystals[5]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        //BoxA GammaRange
        if (crystals[4] && crystals[5])
          {
            fh2_EnergyCorrelationsCrystals[6]->Fill(energies[4],energies[5]);
            fh2_EnergyCorrelationsCrystals[6]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
            fh2_EnergyCorrelationsCrystals[6]->GetYaxis()->SetTitle("Energy [keV] crystal 6");
          }

        if (crystals[4] && crystals[6])
        {
          fh2_EnergyCorrelationsCrystals[7]->Fill(energies[4],energies[6]);
          fh2_EnergyCorrelationsCrystals[7]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[7]->GetYaxis()->SetTitle("Energy [keV] crystal 7");
        }

        if (crystals[4] && crystals[7])
        {
          fh2_EnergyCorrelationsCrystals[8]->Fill(energies[4],energies[7]);
          fh2_EnergyCorrelationsCrystals[8]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[8]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }

        if (crystals[5] && crystals[7])
        {
          fh2_EnergyCorrelationsCrystals[9]->Fill(energies[5],energies[7]);
          fh2_EnergyCorrelationsCrystals[9]->GetXaxis()->SetTitle("Energy [keV] crystal 6");
          fh2_EnergyCorrelationsCrystals[9]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }

        if (crystals[6] && crystals[7])
        {
          fh2_EnergyCorrelationsCrystals[10]->Fill(energies[6],energies[7]);
          fh2_EnergyCorrelationsCrystals[10]->GetXaxis()->SetTitle("Energy [keV] crystal 7");
          fh2_EnergyCorrelationsCrystals[10]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }
        if (crystals[5] && crystals[6])
        {
          fh2_EnergyCorrelationsCrystals[11]->Fill(energies[5],energies[6]);
          fh2_EnergyCorrelationsCrystals[11]->GetXaxis()->SetTitle("Energy [keV] crystal 6");
          fh2_EnergyCorrelationsCrystals[11]->GetYaxis()->SetTitle("Energy [keV] crystal 7");
        }

        //BoxB ProtonRange
        if (crystals[8] && crystals[9])
        {
          fh2_EnergyCorrelationsCrystals[12]->Fill(energies[8],energies[9]);
          fh2_EnergyCorrelationsCrystals[12]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[12]->GetYaxis()->SetTitle("Energy [keV] crystal 2");
        }

        if (crystals[8] && crystals[10])
        {
          fh2_EnergyCorrelationsCrystals[13]->Fill(energies[8],energies[10]);
          fh2_EnergyCorrelationsCrystals[13]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[13]->GetYaxis()->SetTitle("Energy [keV] crystal 3");
        }

        if (crystals[8] && crystals[11])
        {
          fh2_EnergyCorrelationsCrystals[14]->Fill(energies[8],energies[11]);
          fh2_EnergyCorrelationsCrystals[14]->GetXaxis()->SetTitle("Energy [keV] crystal 1");
          fh2_EnergyCorrelationsCrystals[14]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        if (crystals[9] && crystals[10])
        {
          fh2_EnergyCorrelationsCrystals[15]->Fill(energies[9],energies[10]);
          fh2_EnergyCorrelationsCrystals[15]->GetXaxis()->SetTitle("Energy [keV] crystal 2");
          fh2_EnergyCorrelationsCrystals[15]->GetYaxis()->SetTitle("Energy [keV] crystal 3");
        }

        if (crystals[9] && crystals[11])
        {
          fh2_EnergyCorrelationsCrystals[16]->Fill(energies[9],energies[11]);
          fh2_EnergyCorrelationsCrystals[16]->GetXaxis()->SetTitle("Energy [keV] crystal 2");
          fh2_EnergyCorrelationsCrystals[16]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        if (crystals[10] && crystals[11])
        {
          fh2_EnergyCorrelationsCrystals[17]->Fill(energies[10],energies[11]);
          fh2_EnergyCorrelationsCrystals[17]->GetXaxis()->SetTitle("Energy [keV] crystal 3");
          fh2_EnergyCorrelationsCrystals[17]->GetYaxis()->SetTitle("Energy [keV] crystal 4");
        }

        //BoxB GammaRange
        if (crystals[12] && crystals[13])
        {
          fh2_EnergyCorrelationsCrystals[18]->Fill(energies[12],energies[13]);
          fh2_EnergyCorrelationsCrystals[18]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[18]->GetYaxis()->SetTitle("Energy [keV] crystal 6");
        }

        if (crystals[12] && crystals[14])
        {
          fh2_EnergyCorrelationsCrystals[19]->Fill(energies[12],energies[14]);
          fh2_EnergyCorrelationsCrystals[19]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[19]->GetYaxis()->SetTitle("Energy [keV] crystal 7");
        }

        if (crystals[12] && crystals[15])
        {
          fh2_EnergyCorrelationsCrystals[20]->Fill(energies[12],energies[15]);
          fh2_EnergyCorrelationsCrystals[20]->GetXaxis()->SetTitle("Energy [keV] crystal 5");
          fh2_EnergyCorrelationsCrystals[20]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }

        if (crystals[13] && crystals[15])
        {
          fh2_EnergyCorrelationsCrystals[21]->Fill(energies[13],energies[15]);
          fh2_EnergyCorrelationsCrystals[21]->GetXaxis()->SetTitle("Energy [keV] crystal 6");
          fh2_EnergyCorrelationsCrystals[21]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }

        if (crystals[14] && crystals[15])
        {
          fh2_EnergyCorrelationsCrystals[22]->Fill(energies[14],energies[15]);
          fh2_EnergyCorrelationsCrystals[22]->GetXaxis()->SetTitle("Energy [keV] crystal 7");
          fh2_EnergyCorrelationsCrystals[22]->GetYaxis()->SetTitle("Energy [keV] crystal 8");
        }
        if (crystals[13] && crystals[14])
        {
          fh2_EnergyCorrelationsCrystals[23]->Fill(energies[13],energies[14]);
          fh2_EnergyCorrelationsCrystals[23]->GetXaxis()->SetTitle("Energy [keV] crystal 6");
          fh2_EnergyCorrelationsCrystals[23]->GetYaxis()->SetTitle("Energy [keV] crystal 7");
        }

        if (crystals[0] && crystals[8])
        {fh2_Energy_1A_1B->Fill(energies[0],energies[8]);}

      fh2_EnergyCorrelationsAlvProton->Fill(energy_alv[0],energy_alv[2]);
      fh2_EnergyCorrelationsAlvGamma->Fill(energy_alv[1],energy_alv[3]);

      fh1_EnergyTotBoxA_p->Fill(energy_alv[0]);
      fh1_EnergyTotBoxA_g->Fill(energy_alv[1]);
      fh1_EnergyTotBoxB_p->Fill(energy_alv[2]);
      fh1_EnergyTotBoxB_g->Fill(energy_alv[3]);
    }


    fNEvents += 1;
}

void R3BCalifaJulichOnlineSpectra::FinishEvent()
{
    if (fMappedItemsSi)
      {fMappedItemsSi->Clear();}

    if (fMappedItemsCalifa)
      {fMappedItemsCalifa->Clear();}

    if (fCalItemsCalifa)
      {fCalItemsCalifa->Clear();}

    if (fCalItemsSi)
      {fCalItemsSi->Clear();}

    if (fHitItemsCalifa)
      {fHitItemsCalifa->Clear();}

    if (fHitItemsSi)
      {fHitItemsSi->Clear();}
}

void R3BCalifaJulichOnlineSpectra::FinishTask()
{

   for (Int_t i = 0; i < fNbDet; i++)
      {
          if (fMappedItemsSi) { fh2_EnergyVsStrip[i]->Write(); }
          if (fCalItemsSi) { fh2_EnergyCalVsStrip[i]->Write(); }
          if (fHitItemsSi) { fh2_PosX_PosY[i]->Write(); }
      }

      if(fMappedItemsCalifa) fh1_MultiplicityGamma->Write();

      if(fMappedItemsCalifa && fMappedItemsSi)
      {
        for (Int_t i=0;i<4;i++) {fh2_EnergyMapCalifa_SiStrip[i]->Write(); fh2_PosX_PosY_Califa[i]->Write();}
        fh2_EnergyTotMapCalifa_SiStrip->Write();
      }
      for (Int_t i = 0; i < fNbCrystals; i++)
      {
          if (fMappedItemsCalifa){ fh1_EnergyCalifaCrystals[i]->Write(); fh2_Map_nf_ns[i]->Write();}
          if (fCalItemsCalifa){ fh1_EnergyCalCalifaCrystals[i]->Write(); }
          if (fHitItemsCalifa){ fh1_EnergyHitCalifaCrystals[i]->Write(); }
      }

      if (fCalItemsCalifa)
      {
        for (Int_t i=0;i<24;i++)
        {
          fh2_EnergyCorrelationsCrystals[i]->Write();
        }
        fh2_Energy_1A_1B->Write();
      }

      fh2_EnergyCorrelationsAlvProton->Write();
      fh2_EnergyCorrelationsAlvGamma->Write();

      fh1_EnergyTotBoxA_p->Write();
      fh1_EnergyTotBoxA_g->Write();
      fh1_EnergyTotBoxB_p->Write();
      fh1_EnergyTotBoxB_g->Write();

}

ClassImp(R3BCalifaJulichOnlineSpectra);
