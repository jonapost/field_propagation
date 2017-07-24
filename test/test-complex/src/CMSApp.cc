
#include "CMSApp.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4TransportationManager.hh"

#include "G4Material.hh"
#include "G4String.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"

#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"

#include "LogicalVolumeHandler.hh"

G4bool CMSApp::fgIsScoreActive = FALSE;

CMSApp::CMSApp()
: fInitialized(FALSE), 
  fECALMap(), fHCALMap(),
  fFluxElec(0), fFluxGamma(0), fFluxP(0), fFluxPi(0), fFluxK(0),
  fEdepElec(0), fEdepGamma(0), fEdepP(0), fEdepPi(0), fEdepK(0) {
  // Ctor..
  memset(fSensFlags, 0, kNvolumes*sizeof(G4bool));
  memset(fEdepECAL, 0, kNECALModules*sizeof(G4double));
  memset(fEdepHCAL, 0, kNHCALModules*sizeof(G4double));
  memset(fECALid, 0, kNECALModules*sizeof(G4int));
  memset(fHCALid, 0, kNHCALModules*sizeof(G4int));
  memset(fCubicVolumes, 0, kNvolumes*sizeof(G4double));

  if (fgIsScoreActive == FALSE)
    return;

//  TH1::AddDirectory(FALSE);
  fFluxElec = new TH1F("hFluxElec", "e+/e- flux/primary in ECAL", 50, 0., 2500.);
  fFluxElec->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxElec->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fFluxGamma = new TH1F("hFluxGamma", "Gamma flux/primary in ECAL", 50, 0., 2500.);
  fFluxGamma->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxGamma->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fFluxP = new TH1F("hFluxP", "Proton flux/primary in ECAL", 50, 0., 2500.);
  fFluxP->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxP->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fFluxPi = new TH1F("hFluxPi", "Pion flux/primary in ECAL", 50, 0., 2500.);
  fFluxPi->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxPi->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fFluxK = new TH1F("hFluxK", "Kaon flux/primary in ECAL", 50, 0., 2500.);
  fFluxK->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxK->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fEdepElec = new TH1F("hEdepElec", "Electron energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepElec->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepElec->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^3/primary]");
  fEdepGamma = new TH1F("hEdepGamma", "Gamma energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepGamma->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepGamma->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^3/primary]");
  fEdepP = new TH1F("hEdepP", "Proton energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepP->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepP->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^3/primary]");
  fEdepPi = new TH1F("hEdepPi", "Pion energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepPi->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepPi->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^3/primary]");
  fEdepK = new TH1F("hEdepK", "Kaon energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepK->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepK->GetYaxis()->SetTitle("Energy deposit density[MeV/cm^3/primary]");
//  TH1::AddDirectory(TRUE);
}

G4bool CMSApp::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return TRUE;

  // Loop unique volume id's
  const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  G4int necal = 0;
  G4int nhcal = 0;
  const G4String sEBRY = "EBRY";
  const G4String sEFRY = "EFRY";
  for (unsigned int ivol = 0; ivol < lvs->size(); ++ivol) {
  G4LogicalVolume *vol = (*lvs)[ivol];
    // ECAL cells
    G4String svol          = vol->GetName();
    ivol = LogicalVolumeHandler::Instance()->GetIndex(vol);
    G4bool isUsed = LogicalVolumeHandler::Instance()->IsUsed(ivol);
    if(!isUsed) 
      continue;

    if (!(svol.compare(0,sEBRY.length(),sEBRY)) || !(svol.compare(0,sEFRY.length(),sEFRY))  ) {
//      std::cout<< "---  indx = "<< ivol <<"  Name = "<<svol<<std::endl;
      fSensFlags[ivol] = TRUE;
      fECALMap[ivol] = necal;
      fECALid[necal] = ivol;
      ++necal;
    }
    // HCAL cells
    G4String smat = "";
    if(vol->GetMaterial()) smat = vol->GetMaterial()->GetName();
    if (smat == G4String("Scintillator")) {
      fSensFlags[ivol] = TRUE;
      fHCALMap[ivol] = nhcal;
      fHCALid[nhcal] = ivol;
      ++nhcal;
    }
  }

  // precompute cubic volumes 
  for (unsigned int ivol = 0; ivol < lvs->size(); ++ivol)
    if(fSensFlags[ivol])
      fCubicVolumes[ivol] =  (*lvs)[ivol]->GetSolid()->GetCubicVolume() / CLHEP::cm3;


  G4cout<< G4endl<<
  "=== CMSApp::Initialize: necal = "<< necal << "  nhcal = " << nhcal << G4endl << G4endl;
  fInitialized = TRUE;
  return TRUE;
}


void CMSApp::SteppingAction(const G4Step *step) {
  if (fgIsScoreActive == FALSE)
    return;

  if (!fInitialized)
    return; // FOR NOW

  // Check if the track is in the right volume and collect the
  // energy deposit and step length
  G4int ivol;
  G4int idtype;
  G4int mod;
  G4LogicalVolume *vol = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();

  ivol = LogicalVolumeHandler::Instance()->GetIndex(vol);
  idtype = 0;
  if(fSensFlags[ivol]) { 
    if (vol->GetName()[0] == 'E') idtype = 1;
    else idtype = 2;
    switch (idtype) {
      case 1:
        mod = fECALMap.find(ivol)->second;
        fEdepECAL[mod] += step->GetTotalEnergyDeposit()/CLHEP::MeV;
        break;
      case 2:
        mod = fHCALMap.find(ivol)->second;
        fEdepHCAL[mod]+= step->GetTotalEnergyDeposit()/CLHEP::MeV;
        break;
    } // end of idtype switch
  } // end of fSensFlags

  // Score in ECAL
  if (idtype==1) {
//      G4double capacity = fCubicVolumes[ivol];
      G4double capacity = 1.;
      const G4Track           *track = step->GetTrack();
      const G4DynamicParticle *dpart = track->GetDynamicParticle();
      G4int thePDG = track->GetParticleDefinition()->GetPDGEncoding();
      if( std::abs(thePDG)       == 11  ) {
          G4double totMomentum = dpart->GetTotalMomentum() / CLHEP::MeV; 
	  fFluxElec->Fill(totMomentum, step->GetStepLength() / CLHEP::cm /capacity);
	  fEdepElec->Fill(totMomentum, step->GetTotalEnergyDeposit() / CLHEP::MeV /capacity);
      } else if( thePDG          == 22 ||  thePDG == 0 ) {
          G4double totMomentum = dpart->GetTotalMomentum() / CLHEP::MeV; 
/*
            G4cout<< " Name Vol = "<< vol->GetName() << " position = ["
                   << step->GetPostStepPoint()->GetPosition().x() / CLHEP::cm << ", "
                   << step->GetPostStepPoint()->GetPosition().y() / CLHEP::cm << ", "
                   << step->GetPostStepPoint()->GetPosition().z() / CLHEP::cm << "] "
                   << " momentum = " << totMomentum
                   << " step length = " << step->GetStepLength() / CLHEP::cm
                   << " [cm]"<<G4endl;
*/
          fFluxGamma->Fill(totMomentum, step->GetStepLength() / CLHEP::cm /capacity);
          fEdepGamma->Fill(totMomentum, step->GetTotalEnergyDeposit() / CLHEP::MeV /capacity);
      } else if( thePDG          == 2212) {
          G4double totMomentum = dpart->GetTotalMomentum() / CLHEP::MeV; 
	  fFluxP->Fill(totMomentum, step->GetStepLength() / CLHEP::cm /capacity);
	  fEdepP->Fill(totMomentum, step->GetTotalEnergyDeposit() / CLHEP::MeV /capacity);
      } else if(std::abs(thePDG) == 211 ) {
          G4double totMomentum = dpart->GetTotalMomentum() / CLHEP::MeV; 
	  fFluxPi->Fill(totMomentum, step->GetStepLength() / CLHEP::cm /capacity);
	  fEdepPi->Fill(totMomentum, step->GetTotalEnergyDeposit() / CLHEP::MeV /capacity);
      } else if(std::abs(thePDG) == 321 ) {
          G4double totMomentum = dpart->GetTotalMomentum() / CLHEP::MeV; 
	  fFluxK->Fill(totMomentum, step->GetStepLength() / CLHEP::cm /capacity);
	  fEdepK->Fill(totMomentum, step->GetTotalEnergyDeposit() / CLHEP::MeV /capacity);
      }
  } 
}


void CMSApp::EndOfRunAction(G4int numprims) {
  G4double onePerPrims = 1.0/((G4double)numprims);
  const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();

  G4cout<< G4endl << "Energy deposit in ECAL [MeV/primary] " << G4endl;
  G4cout<<"================================================================================"<<G4endl;
  for (G4int i = 0; i < kNECALModules; ++i)
    G4cout<< "   volume " << (*lvs)[fECALid[i]]->GetName() << " edep =  "<< fEdepECAL[i] *onePerPrims << G4endl;
 
  G4cout<< "Energy deposit in HCAL [MeV/primary] "<< G4endl;
  G4cout<< "================================================================================"<<G4endl;
  for (G4int i = 0; i < kNHCALModules; ++i) 
    G4cout<< "   volume " << (*lvs)[fHCALid[i]]->GetName() << " edep =  "<< fEdepHCAL[i] *onePerPrims << G4endl;
  G4cout<< "================================================================================"<<G4endl;

  // if scoring is active 
  if (fgIsScoreActive == FALSE)
    return;

  TFile tf("hists.root","recreate");
   fFluxElec->Scale(onePerPrims);
   fFluxGamma->Scale(onePerPrims);
   fFluxP->Scale(onePerPrims);
   fFluxPi->Scale(onePerPrims);
   fFluxK->Scale(onePerPrims); 
   fEdepElec->Scale(onePerPrims);
   fEdepGamma->Scale(onePerPrims);
   fEdepP->Scale(onePerPrims);
   fEdepPi->Scale(onePerPrims);
   fEdepK->Scale(onePerPrims);   

   fFluxElec->Sumw2();
   fFluxGamma->Sumw2();
   fFluxP->Sumw2();
   fFluxPi->Sumw2();
   fFluxK->Sumw2(); 
   fEdepElec->Sumw2();
   fEdepGamma->Sumw2();
   fEdepP->Sumw2();
   fEdepPi->Sumw2();
   fEdepK->Sumw2();   

   fFluxElec->Write();
   fFluxGamma->Write();
   fFluxP->Write();
   fFluxPi->Write();
   fFluxK->Write(); 
   fEdepElec->Write();
   fEdepGamma->Write();
   fEdepP->Write();
   fEdepPi->Write();
   fEdepK->Write();   
}



