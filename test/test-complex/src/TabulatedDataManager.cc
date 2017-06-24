#ifdef USE_VECGEOM_NAVIGATOR
#undef USE_VECGEOM_NAVIGATOR
#define RESTORE_USE_VECGEOM_NAVIGATOR
#endif

#include "G4Track.hh"
#include "TabulatedDataManager.hh"

#include "TSystem.h"
#include "TFile.h"
#include "TError.h"
#include "TBits.h"
#include "TMath.h"

#include "TPartIndex.h"
#include "TEXsec.h"
#include "TMXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoExtension.h"

#include "G4ParticleChange.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

#include "G4ProcessManager.hh"

#include "G4Gamma.hh"

using CLHEP::GeV;
using CLHEP::cm;

#ifdef RESTORE_USE_VECGEOM_NAVIGATOR
#define USE_VECGEOM_NAVIGATOR
#endif

//#ifdef MAKESTAT
// unsigned long TabulatedDataManager::killedTracks =0;
//#endif

TabulatedDataManager *TabulatedDataManager::fgInstance = 0;
TGeoManager *TabulatedDataManager::fgGeom = 0; // Pointer to the geometry manager
G4bool TabulatedDataManager::fgIsUseRange = FALSE;

const char *TabulatedDataManager::tStatus[] = {"fAlive", "fStopButAlive", "fStopAndKill"};
G4int TabulatedDataManager::fgVerboseLevel = 0;

TabulatedDataManager *TabulatedDataManager::Instance() {
  if (fgInstance == 0) {
    // char* gdmlFileName = getenv("VP_GEOM_GDML");
    char *xsecFileName = getenv("VP_DATA_XSEC");
    char *fstaFileName = getenv("VP_DATA_FSTA");

    if (!fgGeom) {
      // if(gdmlFileName && xsecFileName && fstaFileName) {
      // fGeom = TGeoManager::Import(gdmlFileName);

      ::Error("TabulatedDataManager::Instance", "Missing pointer to TGeomManager");
      return 0;
    } else {
      if (xsecFileName && fstaFileName) {
        fgInstance = new TabulatedDataManager(fgGeom, xsecFileName, fstaFileName);
      } else {
        ::Error("TabulatedDataManager::Instance", "Missing VP_DATA_XSEC VP_DATA_FSTA");
        exit(1);
        return 0;
      }
    }
  }
  return fgInstance;
}

TabulatedDataManager::~TabulatedDataManager() {
  fgInstance = 0;
  delete[] fMatXsec;
  delete[] fElemXsec;
  delete[] fElemFstate;
  delete fDecay;
  delete fHasNCaptureAtRest;
}

TabulatedDataManager::TabulatedDataManager()
    : fNelements(0), fNmaterials(0), fElemXsec(0), fElemFstate(0), fMatXsec(0), fDecay(0), fHasNCaptureAtRest(0)
// fgGeom(0)
{}

TabulatedDataManager::TabulatedDataManager(TGeoManager *geom, const char *xsecfilename, const char *finalsfilename)
    : fNelements(0), fNmaterials(0), fElemXsec(0), fElemFstate(0), fMatXsec(0), fDecay(0), fHasNCaptureAtRest(0)
//, fgGeom(geom)
{
  // this is clone of TTabPhysMgr::TTabPhysMgr(TGeoManager* geom,
  //                 const char* xsecfilename, const char* finalsfilename):

  std::cout << "TabulatedDataManager - constructor called." << std::endl;

  if (fgGeom != geom) {
    Fatal("TabulateDataManager", "Conflicting pointers to TGeoManager");
  }

  // Open xsec_FTFP_BERT.root file and fstate_FTFP_BERT.root
  TFile *fxsec = TFile::Open(xsecfilename);
  if (!fxsec) {
    Fatal("TabulatedDataManager", "Cannot open %s", xsecfilename);
  }

  // Get the TPartIndex singleton object from the xsec file
  fxsec->Get("PartIndex");
  // check version of the data files
  if (fgVersion != TPartIndex::I()->Version()) {
    std::cerr << "\n\n*************************************************************\n"
              << "  ---------------------------ERROR-----------------------------\n"
              << "    Your xsec_*.root and fstate_*.root data files at           \n"
              << "    -> " << xsecfilename << "\n"
              << "    -> " << finalsfilename << "\n"
              << "    Version is       : " << TPartIndex::I()->VersionMajor() << "." << TPartIndex::I()->VersionMinor()
              << "." << TPartIndex::I()->VersionSub() << "\n"
              << "    Required version : " << GetVersion() << "\n"
              << "    Update your xsec_*.root and fstate_*.root data files !     "
              << "\n*************************************************************\n\n";
    exit(EXIT_FAILURE);
  }

  TFile *fstate = TFile::Open(finalsfilename);
  if (!fstate) {
    Fatal("TabulatedDataManager", "Cannot open %s", finalsfilename);
    exit(EXIT_FAILURE);
  }

  // Load elements from geometry
  TList *matlist = (TList *)geom->GetListOfMaterials();

  TIter next(matlist);
  TGeoMaterial *mat = 0;

  // get the decay table from the final state file
  fDecay = (TPDecay *)fstate->Get("DecayTable");

  // INFO: print number of materials in the current TGeoManager
  printf("#materials:= %d \n", matlist->GetSize());

  // First loop on all materials to mark used elements
  TBits elements(NELEM);
  while ((mat = (TGeoMaterial *)next())) {
    std::cout << "TabulatedDataManager> Checking material " << mat->GetName() << std::endl;
    if (!mat->IsUsed() || mat->GetZ() < 1.)
      continue;
    fNmaterials++;
    int nelem = mat->GetNelements();
    // Check if we are on the safe side; should exit otherwise
    if (nelem > MAXNELEMENTS) {
      Fatal("TabulatedDataManager", "Number of elements in %s is %d > MAXNELEMENTS=%d\n", mat->GetName(), nelem,
            MAXNELEMENTS);
    }
    for (int iel = 0; iel < nelem; ++iel) {
      double ad;
      double zd;
      double wd;
      mat->GetElementProp(ad, zd, wd, iel);
      if (zd < 1 || zd > NELEM) {
        Fatal("TabulatedDataManager", "In material %s found element with z=%d > NELEM=%d", mat->GetName(), (int)zd,
              NELEM);
      }
      elements.SetBitNumber(zd);
    }
  }
  fNelements = elements.CountBits();
  fElemXsec = new TEXsec *[NELEM];
  fElemFstate = new TEFstate *[NELEM];
  fMatXsec = new TMXsec *[fNmaterials];
  printf("Reading xsec and final states for %d elements in %d materials\n", fNelements, fNmaterials);

  // Loop elements and load corresponding xsec and final states
  int zel = elements.FirstSetBit();
  int nbits = elements.GetNbits();
  TEXsec *exsec;
  TEFstate *estate;

  // Load elements xsec data in memory
  ProcInfo_t procInfo1, procInfo2;
  gSystem->GetProcInfo(&procInfo1);

  while (zel < nbits) {

    exsec = TEXsec::GetElement(zel, 0, fxsec);
    fElemXsec[zel] = exsec;
    fElemXsec[zel]->SetIndex(zel); // quick access to the corresponding fstate
    estate = TEFstate::GetElement(zel, 0, fstate);
    fElemFstate[zel] = estate;
    printf("   loaded xsec data and states for: %s\n", TPartIndex::I()->EleSymb(zel));
    zel = elements.FirstSetBit(zel + 1);
    // init : does the particle have nuclear cpature at rest? array
    if (!fHasNCaptureAtRest) {
      G4int numParticles = TPartIndex::I()->NPart();
      fHasNCaptureAtRest = new G4bool[numParticles];
      for (G4int ip = 0; ip < numParticles; ++ip)
        fHasNCaptureAtRest[ip] = estate->HasRestCapture(ip);
    }
  }

  gSystem->GetProcInfo(&procInfo2);
  // long mem = (procInfo2.fMemResident - procInfo1.fMemResident)/1024;
  fxsec->Close();
  fstate->Close();

  // xsec and states now in memory
  // Go through all materials in the geometry and form the associated TMXsec
  // objects.

  int *z = new int[MAXNELEMENTS];
  int *a = new int[MAXNELEMENTS];
  float *w = new float[MAXNELEMENTS];
  fNmaterials = 0;
  next.Reset();

  while ((mat = (TGeoMaterial *)next())) {
    if (!mat->IsUsed())
      continue;
    int nelem = mat->GetNelements();
    // loop over the elements of the current material in order to obtain the
    // z, a, w, arrays of the elements of this material
    double ad;
    double zd;
    double wd;
    for (int iel = 0; iel < nelem; ++iel) {
      mat->GetElementProp(ad, zd, wd, iel);
      a[iel] = ad;
      z[iel] = zd;
      w[iel] = wd;
    }
    // Construct the TMXsec object that corresponds to the current material
    TMXsec *mxs = new TMXsec(mat->GetName(), mat->GetTitle(), z, a, w, nelem, mat->GetDensity(), true, fDecay);
    fMatXsec[fNmaterials++] = mxs;
    // Connect to TGeoMaterial
    mat->SetFWExtension(new TGeoRCExtension(new TOMXsec(mxs)));
  } // End of while

  delete[] z;
  delete[] a;
  delete[] w;

  int nelements = TEXsec::NLdElems();
  if (nelements != fNelements)
    Error("TabulatedDataManager", "Number of elements not matching");

  // INFO: print some info for checking
  printf("number of materials in fMatXsec[]:= %d\n", fNmaterials);
  for (int i = 0; i < fNmaterials; ++i)
    printf("   fMatXsec[%d]: %s\n", i, fMatXsec[i]->GetName());
}

//_____________________________________________________________________________
void TabulatedDataManager::EnergyLoss(G4int imat, const G4Track &atrack, const G4Step &astep,
                                      G4ParticleChange *particlechange, G4double energylimit) {
  G4int partIndex;    // GV particle index
  G4double kinEnergy; // kinetic energy of the particle in GeV

  G4double stepLength = astep.GetStepLength() / CLHEP::cm; // from mm->cm
  if (stepLength <= 0.0)
    return;

  partIndex = TPartIndex::I()->PartIndex(atrack.GetParticleDefinition()->GetPDGEncoding());
  if (partIndex >= TPartIndex::I()->NPartCharge())
    return;
  kinEnergy = astep.GetPreStepPoint()->GetKineticEnergy() / CLHEP::GeV; // from MeV->GeV

  if (imat < 0 || imat >= fNmaterials) {
    std::cout << "\n!!*******SAMPLING ENERY LOSS*****USING*GV-TABPHYS***************!!\n"
              << "***  Particle is       = " << TPartIndex::I()->PartName(partIndex) << " \n"
              << "***  Particle E kin.   = " << kinEnergy << " [GeV]"
              << " \n"
              << "***  Material index:   = " << imat << " OUT OF RANGE: [ 0 ," << fNmaterials << " ]!"
              << " \n"
              << "!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << std::endl;
    Fatal("TabulatedDataManager::EnergyLoss", "Material index is out of range!");
    exit(EXIT_FAILURE);
  }

#ifdef MAKESTAT
  // set the G4 'step defined by process' to the selected discreate process name
  G4VProcess *proc = (*(atrack.GetParticleDefinition()->GetProcessManager()->GetProcessList()))[1];
  G4String *strp = const_cast<G4String *>(&(proc->GetProcessName()));
  *strp = TPartIndex::I()->ProcName(2); // Ionization
#endif

  G4double range0 = 0.0;
  if (fgIsUseRange) {
    range0 = fMatXsec[imat]->Range(partIndex, kinEnergy);
    if (range0 < 0.)
      return; // there is no range -> there is no dE/dx

    if (stepLength >= range0) {
      particlechange->ProposeEnergy(0.0);
      particlechange->ProposeLocalEnergyDeposit(kinEnergy * GeV);
      if (HasRestProcess(partIndex))
        particlechange->ProposeTrackStatus(fStopButAlive);
      else
        particlechange->ProposeTrackStatus(fStopAndKill);
      return;
    }
  }

  // Get dE/dx and go for energy loss computation:
  G4double dedx = fMatXsec[imat]->DEdx(partIndex, kinEnergy);
  G4double edepo = dedx * stepLength;
  if (fgIsUseRange && edepo > 0.01 * kinEnergy) { // long step
    G4double invrange = fMatXsec[imat]->InvRange(partIndex, range0 - stepLength);
    edepo = kinEnergy - invrange;
  }

  G4double finalEkin = kinEnergy - edepo;

  // If Ekin-EnergyLoss is above the cut then update the current track.
  // Otherwise: particle is stopped, Ekin goes to energy deposit, status is set
  // to StopButAlive or StopAndKill depending on if the particle does/doesn't
  // have NuclearCaptureAtRest process.

  if (finalEkin > energylimit) {
    particlechange->ProposeEnergy(finalEkin * GeV); // from GeV->MeV
    particlechange->ProposeTrackStatus(fAlive);
    particlechange->ProposeLocalEnergyDeposit(edepo * GeV); // from GeV->MeV
  } else {
    particlechange->ProposeEnergy(0.0);
    particlechange->ProposeLocalEnergyDeposit(kinEnergy * GeV);
    if (finalEkin <= 0.0 && HasRestProcess(partIndex)) {
      particlechange->ProposeTrackStatus(fStopButAlive);
    } else {
      particlechange->ProposeTrackStatus(fStopAndKill);
      //#ifdef MAKESTAT
      //       ++killedTracks;
      //#endif
    }
  }

  if (fgVerboseLevel >= 2)
    std::cout << "\n=====COMPUTING ALONG-STEP ENERY LOSS==USING=GV-TABPHYS=(dE/dx)====\n"
              << "***  Particle is       = " << TPartIndex::I()->PartName(partIndex) << " \n"
              << "***  Particle E kin.   = " << kinEnergy << " [GeV]"
              << " \n"
              << "***  Selected reaction = index: " << TPartIndex::I()->ProcIndex("Ionisation")
              << " which is Ionisation. "
              << " \n"
              << "***  Material name     = " << fMatXsec[imat]->GetName() << " \n"
              << "***  Step lenght       = " << astep.GetStepLength() / cm << "[cm]"
              << " \n"
              << "***  dE/dx             = " << dedx << "[GeV/cm] Eloss:= " << edepo << " [GeV]"
              << "\n"
              << "\n==================================================================\n" << std::endl;
}

//_____________________________________________________________________________
G4double TabulatedDataManager::GetRange(const G4int imat, const G4Track &atrack) {
  G4double x = DBL_MAX;
  G4int partIndex;    // GV particle index
  G4double kinEnergy; // kinetic energy of the particle in GeV

  partIndex = TPartIndex::I()->PartIndex(atrack.GetParticleDefinition()->GetPDGEncoding());
  kinEnergy = atrack.GetKineticEnergy() / CLHEP::GeV; // from MeV->GeV

  if (imat < 0 || imat >= fNmaterials) {
    std::cout << "\n!!************GETTING-LINEAR-RANGE***USING*GV-TABPHYS************!!\n"
              << "***  Particle is       = " << TPartIndex::I()->PartName(partIndex) << " \n"
              << "***  Particle E kin.   = " << kinEnergy << " [GeV]"
              << " \n"
              << "***  Material index:   = " << imat << " OUT OF RANGE: [ 0 ," << fNmaterials << " ]!"
              << " \n"
              << "!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << std::endl;
    Fatal("TabulatedDataManager::GetInteractionLength", "Material index is out of range!");
    exit(EXIT_FAILURE);
  }

  // Get the linear range
  if (partIndex >= TPartIndex::I()->NPartCharge())
    return -1.0;

  x = fMatXsec[imat]->Range(partIndex, kinEnergy);
  return x; // from cm->mm
}

//_____________________________________________________________________________
G4double TabulatedDataManager::GetInteractionLength(G4int imat, const G4Track &atrack) {
  G4double x = DBL_MAX;
  G4int partIndex;    // GV particle index
  G4double kinEnergy; // kinetic energy of the particle in GeV
  G4double ptotal;    // total momentum

  partIndex = TPartIndex::I()->PartIndex(atrack.GetParticleDefinition()->GetPDGEncoding());
  kinEnergy = atrack.GetKineticEnergy() / CLHEP::GeV;                    // from MeV->GeV
  ptotal = atrack.GetDynamicParticle()->GetTotalMomentum() / CLHEP::GeV; // from MeV->GeV

  if (imat < 0 || imat >= fNmaterials) {
    std::cout << "\n!!********SAMPLING-INTERACTION-LENGTH***USING*GV-TABPHYS*********!!\n"
              << "***  Particle is       = " << TPartIndex::I()->PartName(partIndex) << " \n"
              << "***  Particle E kin.   = " << kinEnergy << " [GeV]"
              << " \n"
              << "***  Material index:   = " << imat << " OUT OF RANGE: [ 0 ," << fNmaterials << " ]!"
              << " \n"
              << "!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << std::endl;
    Fatal("TabulatedDataManager::GetInteractionLength", "Material index is out of range!");
    exit(EXIT_FAILURE);
  }

  // Get the total mean free path
  x = fMatXsec[imat]->Xlength(partIndex, kinEnergy, ptotal);
  return x * cm; // from cm->mm
}

//____________________________________________________________________________
// sampling element for interaction and type of interaction on that element
int TabulatedDataManager::SampleInteraction(const G4int imat, const G4Track &atrack, int &reactionid) {
  G4int partIndex;    // GV particle index
  G4double kinEnergy; // kinetic energy of the particle in GeV
  G4double ptotal;    // total momentum
  partIndex = TPartIndex::I()->PartIndex(atrack.GetParticleDefinition()->GetPDGEncoding());
  ptotal = atrack.GetDynamicParticle()->GetTotalMomentum() / CLHEP::GeV; // from MeV->GeV

  //  kinEnergy = atrack.GetKineticEnergy()/CLHEP::GeV; // from MeV->GeV
  // TO BE THE SAME AS G4 WITHOUT integral approach !!! IT IS DIFFERENT IT PROTO.
  // BUT NO EFFECT IF THE TRACKING LIMIT IS HIGHER THAN 3keV
  //  if(fgIsUseRange)
  //     kinEnergy = atrack.GetStep()->GetPreStepPoint()->GetKineticEnergy()/GeV;  // WHEN COMPARED TO G4
  //  else
  //     kinEnergy = atrack.GetStep()->GetPostStepPoint()->GetKineticEnergy()/GeV;

  kinEnergy = atrack.GetStep()->GetPreStepPoint()->GetKineticEnergy() / GeV; // WHEN COMPARED TO G4
  //    kinEnergy = atrack.GetStep()->GetPostStepPoint()->GetKineticEnergy()/GeV;   // WHEN COMPARED TO PROTO

  // sampling element for intercation based on element-wise relative tot-xsecs
  // and sampling the interaction itself based on the relative total xsections
  // on the selected element
  TEXsec *elemXsec = fMatXsec[imat]->SampleInt(partIndex, kinEnergy, reactionid, ptotal);

  // if deacy was selected -> elemXsec = NULL
  if (elemXsec || reactionid == 3) {
    if (fgVerboseLevel >= 2) {
      G4String elementName;
      if (reactionid == 3)
        elementName = "nothing";
      else
        elementName = elemXsec->GetName();
      std::cout << "\n=========SAMPLING INTERACTION====USING=GV-TABPHYS=================\n"
                << "***  Particle is       = " << TPartIndex::I()->PartName(partIndex) << " \n"
                << "***  Particle E kin.   = " << kinEnergy << " [GeV]"
                << " \n"
                << "***  Selected element  = " << elementName << " \n"
                << "***  Selected reaction = index: " << reactionid << " which is "
                << TPartIndex::I()->ProcName(reactionid) << " \n"
                << "==================================================================\n" << std::endl;
    }
#ifdef MAKESTAT
    // set the G4 'step defined by process' to the selected discreate process name
    G4VProcess *proc = (*(atrack.GetParticleDefinition()->GetProcessManager()->GetProcessList()))[1];
    G4String *strp = const_cast<G4String *>(&(proc->GetProcessName()));
    *strp = TPartIndex::I()->ProcName(reactionid);
#endif
    if (elemXsec)
      return elemXsec->Index();
    else
      return -1; // deacy -> no element selected
  } else {
    std::cout << "\n!!*******SAMPLING INTERACTION****USING*GV-TABPHYS***************!!\n"
              << "***  Particle is       = " << TPartIndex::I()->PartName(partIndex) << " \n"
              << "***  Particle E kin.   = " << kinEnergy << " [GeV]"
              << " \n"
              << "!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << std::endl;
    Fatal("TabulatedDataManager::SampleInteraction", "No element selected!");
    exit(EXIT_FAILURE);
  }
}

//_____________________________________________________________________________
// sampling final state, put 2ndaries into the particle change, update primary,
// do proper transformations if necessary
void TabulatedDataManager::SampleFinalState(const int elementindex, const int reactionid, const G4Track &atrack,
                                            G4ParticleChange *particlechange, double energylimit) {
  double totEdepo = 0.0;
  int nSecPart = 0;     // number of secondary particles per reaction
  const int *pid = 0;   // GeantV particle codes [nSecPart]
  const float *mom = 0; // momentum vectors the secondaries [3*nSecPart]
  float energyFst = 0;  // Ekin of primary after the interaction
  float kerma = 0;      // released energy
  float weightFst = 0;  // weight of the fstate (just a dummy parameter now)
  char isSurv = 0;      // is the primary survived the interaction
  int ebinindx = -1;    // energy bin index of the selected final state

  int partindex = TPartIndex::I()->PartIndex(atrack.GetParticleDefinition()->GetPDGEncoding());
  double kinEnergy = atrack.GetDynamicParticle()->GetKineticEnergy() / CLHEP::GeV;

  if (reactionid == 3) { // decay in fligth
    SampleDecayInFlight(partindex, atrack, particlechange, energylimit);
    particlechange->ProposeEnergy(0.0);
    particlechange->ProposeMomentumDirection(0.0, 0.0, 1.0); // not since <-Ekin =0.
    particlechange->ProposeTrackStatus(fStopAndKill);
    return;
  }

  isSurv = fElemFstate[elementindex]->SampleReac(partindex, reactionid, kinEnergy, nSecPart, weightFst, kerma,
                                                 energyFst, pid, mom, ebinindx);

  // it is the case of: pre-step energy sigma is not zero of this interaction
  //                    but post step is zero-> we don't have final state for
  //                    this interaction at the postStep energy bin-> do nothing
  // This can happen if interaction type is selected based on the pre-step energy
  // and there is some energy loss along the step. (i.e. when using ranges)
  if (isSurv && energyFst < 0) {
    particlechange->ProposeEnergy(kinEnergy * GeV);
    particlechange->SetNumberOfSecondaries(0);
    particlechange->ProposeTrackStatus(fAlive);
    return;
  }

  // setting the final state correction factor (we scale only the 3-momentums)
  //-get mass of the primary
  int primPDG = TPartIndex::I()->PDG(partindex); // GV part.code -> PGD code
  TParticlePDG *primPartPDG = TDatabasePDG::Instance()->GetParticle(primPDG);
  double primMass = primPartPDG->Mass(); // mass [GeV]
                                         //-compute corFactor = P_current/P_original = Pz_current/Pz_original
                                         // (normaly a check would be good but not necessary: if(ebinindx<0 -> ...)
  double orgPrimEkin = (TPartIndex::I()->EGrid())[ebinindx];
  double corFactor =
      std::sqrt(kinEnergy * (kinEnergy + 2.0 * primMass) / (orgPrimEkin * (orgPrimEkin + 2.0 * primMass)));

  //-if corFactor is set here to 1.0 --> no correction of the final states
  // corFactor = 1.0;

  // we should correct the kerma as well but we don't have enough information
  totEdepo = kerma;

  // store original (pre-interaction) direction of the primary particle
  double oldXdir = atrack.GetMomentumDirection().x();
  double oldYdir = atrack.GetMomentumDirection().y();
  double oldZdir = atrack.GetMomentumDirection().z();

  // set the non-corrected post-interaction Ekin. of the primary
  double postEkinOfParimary = energyFst;

  // check if we need to correct the post-interaction Ekin of the primary:
  // if the primary is survived and has non-zero Ekin --> compute its corrected Ekin
  if (isSurv && (postEkinOfParimary > 0.0)) { // survived
    // get corrected 3-momentum of the post-interaction primary
    double px = mom[0];
    double py = mom[1];
    double pz = mom[2];
    px *= corFactor;
    py *= corFactor;
    pz *= corFactor;
    // compute corrected P^2 in [GeV^2]
    double postPrimP2 = px * px + py * py + pz * pz;
    // recompute post-interaction Ekin of the primary with corrected 3-momentum
    postEkinOfParimary = std::sqrt(postPrimP2 + primMass * primMass) - primMass;
  }

  // now we can do the real job
  if (postEkinOfParimary > energylimit) { // survived even after the correction and the E-limit.
    double px = mom[0];
    double py = mom[1];
    double pz = mom[2];
    px *= corFactor;
    py *= corFactor;
    pz *= corFactor;
    // compute corrected P^2 in [GeV^2]
    double postPrimP2 = px * px + py * py + pz * pz;

    particlechange->ProposeTrackStatus(fAlive);
    particlechange->ProposeEnergy(postEkinOfParimary * GeV); // from GeV->MeV
    // rotate direction of primary particle
    double postPrimP = std::sqrt(postPrimP2); // total P [GeV]

    G4ThreeVector newDir(px / postPrimP, py / postPrimP, pz / postPrimP);
    RotateNewTrack(oldXdir, oldYdir, oldZdir, newDir);
    particlechange->ProposeMomentumDirection(newDir);
  } else {
    // Primary particle is stopped, Ekin goes to energy deposit, status is set
    // to StopButAlive or StopAndKill depending on if the particle does/doesn't
    // have NuclearCaptureAtRest or Deacy processes.
    totEdepo += postEkinOfParimary;
    particlechange->ProposeEnergy(0.0);
    particlechange->ProposeMomentumDirection(0.0, 0.0, 1.0); // not since <-Ekin =0.
    if (isSurv && postEkinOfParimary <= 0.0 && HasRestProcess(partindex)) {
      particlechange->ProposeTrackStatus(fStopButAlive);
    } else {
      particlechange->ProposeTrackStatus(fStopAndKill);
      //#ifdef MAKESTAT
      //        ++killedTracks;
      //#endif
    }
  }

  // go for the real secondary particles
  int isec = 0;
  int j = 0;
  if (isSurv)
    ++j; // skipp the first that is the post-interaction primary

  if (fgVerboseLevel >= 2)
    std::cout << "==========SECONDARIES=FROM=TAB.=PHYS===POST-STEP=============\n"
              << "***  Primary is a  : = " << TPartIndex::I()->PartName(partindex) << " \n"
              << "***  Kinetic Enery : = " << kinEnergy << " [GeV] "
              << " \n"
              << "***  Element is    : = " << fElemXsec[elementindex]->GetName() << " \n"
              << "***  Interaction is: = " << TPartIndex::I()->ProcName(reactionid) << " \n"
              << "***  LIST OF SECONDARIES: \n"
              << "***  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
              << " \n"
              << "***  Primary status: = " << tStatus[particlechange->GetTrackStatus()] << " \n"
              << "***  Kinetic Enery : = " << particlechange->GetEnergy() / GeV << " [GeV]"
              << " \n"
              << "***  Momentum dir. : = [ " << *(particlechange->GetMomentumDirection()) << " ]"
              << " \n"
              << "***  ______________________________________________________" << std::endl;

  std::vector<G4Track *> secTracks;
  G4int totalNumSec = 0;
  G4double time = atrack.GetGlobalTime();
  G4double weight = particlechange->GetParentWeight();

  for (isec = j; isec < nSecPart; ++isec) {
    if (pid[isec] >= TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
      int idummy = pid[isec] - 1000000000;
      int Z = idummy / 10000.;
      int A = (idummy - Z * 10000) / 10.;
      double secMass = TPartIndex::I()->GetAprxNuclearMass(Z, A);
      // get corrected 3-momentum of the post-interaction primary
      double px = mom[3 * isec];
      double py = mom[3 * isec + 1];
      double pz = mom[3 * isec + 2];
      px *= corFactor;
      py *= corFactor;
      pz *= corFactor;
      // compute corrected P^2 in [GeV^2]
      double secP2 = px * px + py * py + pz * pz; // total P^2 [GeV^2]
      // compute Ekin of the fragment and put to Edepo
      totEdepo += std::sqrt(secP2 + secMass * secMass) - secMass;
      //#ifdef MAKESTAT
      //       ++killedTracks;
      //#endif
      continue;
    }

    int secPDG = TPartIndex::I()->PDG(pid[isec]); // GV part.code -> PGD code
    TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
    double secMass = secPartPDG->Mass(); // mass [GeV]
    // get corrected 3-momentum of the post-interaction primary
    double px = mom[3 * isec];
    double py = mom[3 * isec + 1];
    double pz = mom[3 * isec + 2];
    px *= corFactor;
    py *= corFactor;
    pz *= corFactor;
    // compute corrected P^2 in [GeV^2]
    double secP2 = px * px + py * py + pz * pz; // total P^2 [GeV^2]
    // compute Ekin of the secondary
    double secEtot = std::sqrt(secP2 + secMass * secMass); // total E [GeV]
    double secEkin = secEtot - secMass;                    // kinetic energy in [GeV]

    // Check if Ekin of this secondary is above our energy limit and add it to
    // the list of secondary tracks.
    // Otherwise: Ekin of this secondary goes to energy deposit and check if
    // this secondary has NuclearCaptureAtRest process:
    // If it does has: insert into the list of secondary tracks with a status
    // of StopButAlive and Ekin = 0.0 (momentum direction is unimoprtant)).
    // If it doesn't have: there is nothing else to do! (drink a coffee!)
    if (secEkin > energylimit) {
      ++totalNumSec;

      G4ParticleDefinition *particleDef = G4ParticleTable::GetParticleTable()->FindParticle(secPDG);

      double secP = std::sqrt(secP2); // total P [GeV]
      G4ThreeVector newDir(px / secP, py / secP, pz / secP);
      RotateNewTrack(oldXdir, oldYdir, oldZdir, newDir);

      G4DynamicParticle *dynamicParticle = new G4DynamicParticle(particleDef, newDir, secEkin * GeV);

      G4Track *secTrack = new G4Track(dynamicParticle, time, atrack.GetPosition());
      secTrack->SetKineticEnergy(secEkin * GeV);
      secTrack->SetTrackStatus(fAlive);
      secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
      secTrack->SetWeight(weight);

      secTracks.push_back(secTrack);
      if (fgVerboseLevel >= 2)
        std::cout << "###  " << totalNumSec << "-th SECONDARY:"
                  << " \n"
                  << "***  Secondary is a : = " << TPartIndex::I()->PartName(pid[isec]) << " \n"
                  << "***  Momentum dir.  : = [ " << newDir[0] << ", " << newDir[1] << ", " << newDir[2] << " ]"
                  << " \n"
                  << "***  Kinetic Energy : = " << secEkin << " [GeV]"
                  << " \n"
                  << "***  ______________________________________________________" << std::endl;

    } else {
      totEdepo += secEkin;
      if (secEkin <= 0.0 && HasRestProcess(pid[isec])) {
        ++totalNumSec;

        G4ParticleDefinition *particleDef = G4ParticleTable::GetParticleTable()->FindParticle(secPDG);

        G4ThreeVector newDir(0, 0, 1); // not important since Ekin = 0.;
        G4DynamicParticle *dynamicParticle = new G4DynamicParticle(particleDef, newDir, 0.0);

        G4Track *secTrack = new G4Track(dynamicParticle, time, atrack.GetPosition());
        secTrack->SetKineticEnergy(0.0);
        secTrack->SetTrackStatus(fStopButAlive);
        secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
        secTrack->SetWeight(weight);

        secTracks.push_back(secTrack);
      } else {
        //#ifdef MAKESTAT
        //         ++killedTracks;
        //#endif
      } // else: do nothing just add secEkin to Edepo that is already done
    }
  }

  // all this verbosity story could be done in a nicer way using the partcile
  // change object after the for-loop below or our G4Track* vector here but I
  // don't care now (can be polished later)
  if (fgVerboseLevel >= 2)
    std::cout << "============================================================\n"
              << "---Total energy deposit := " << totEdepo << " [Gev] \n"
              << "============================================================\n" << std::endl;

  particlechange->SetNumberOfSecondaries(totalNumSec);
  for (G4int i = 0; i < totalNumSec; ++i)
    particlechange->AddSecondary(secTracks[i]);

  // Set the overall energy deposit
  particlechange->ProposeLocalEnergyDeposit(totEdepo * GeV); // from GeV->MeV
}

// Will be called only if the particle has decay or/and nuclear capture at rest!
// We have only nuclear capture at rest that includes final states from possible
// bound deacy as well (but we cannot see if it was nCapture or decay->always set
// to nCapture) So do nuclear capture at rest if the partcle has nCapture at rest
// and decay at rest otherwise.
//______________________________________________________________________________
void TabulatedDataManager::SampleFinalStateAtRest(const int imat, const G4Track &atrack,
                                                  G4ParticleChange *particlechange, double energylimit) {
  // fist we kill the current track; its energy ahs already been put into depo
  particlechange->ProposeTrackStatus(fStopAndKill);

  const double mecc = 0.00051099906; // e- mass c2 in [GeV]

  // sample one of the nuclear capture at rest final states for this particle
  // on the sampled element
  double totEdepo = 0.0;
  int nSecPart = 0;     // number of secondary particles per reaction
  const int *pid = 0;   // GeantV particle codes [nSecPart]
  const float *mom = 0; // momentum vectors the secondaries [3*nSecPart]
  float energyFst = 0;  // Ekin of primary after the interaction
  float kerma = 0;      // released energy
  float weightFst = 0;  // weight of the fstate (just a dummy parameter now)
  char isSurv = 0;      // is the primary survived the interaction
  G4int elementIndex = -1;
  TEFstate *elemfstate = 0;

  // get the GV particle index
  int partindex = TPartIndex::I()->PartIndex(atrack.GetParticleDefinition()->GetPDGEncoding());

  // check if particle is e+ : e+ annihilation at rest if $E_{limit}< m_{e}c^{2}$
  if (partindex == TPartIndex::I()->GetSpecGVIndex(1)) {
    if (energylimit < mecc) {
      double randDirZ = 1.0 - 2.0 * G4UniformRand();
      double randSinTheta = std::sqrt(1.0 - randDirZ * randDirZ);
      double randPhi = G4UniformRand() * CLHEP::twopi;
      double randDirX = randSinTheta * std::cos(randPhi);
      double randDirY = randSinTheta * std::sin(randPhi);

      particlechange->SetNumberOfSecondaries(2);

      //--------------------------------------------------------------
      G4ParticleDefinition *particleDef = G4Gamma::Definition();
      G4ThreeVector dir(randDirX, randDirY, randDirZ);

      G4DynamicParticle *dynamicParticle = new G4DynamicParticle(particleDef, dir, mecc * GeV);

      particlechange->AddSecondary(dynamicParticle);

      dynamicParticle = new G4DynamicParticle(particleDef, -1. * dir, mecc * GeV);

      particlechange->AddSecondary(dynamicParticle);
//---------------------------------------------------------------

/*
       std::vector<G4Track*> secTracks;
       G4double time = atrack.GetGlobalTime();
       G4double weight = particlechange->GetParentWeight();

       G4Track* secTrack =
                       new G4Track(dynamicParticle, time, atrack.GetPosition());
       secTrack->SetKineticEnergy(mecc*GeV);
       secTrack->SetTrackStatus(fAlive);
       secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
       secTrack->SetWeight(weight);

       particlechange->AddSecondary(secTrack);

       dynamicParticle = new G4DynamicParticle(particleDef, -1.0*dir, mecc*GeV);

       secTrack = new G4Track(dynamicParticle, time, atrack.GetPosition());
       secTrack->SetKineticEnergy(mecc*GeV);
       secTrack->SetTrackStatus(fAlive);
       secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
       secTrack->SetWeight(weight);

       particlechange->AddSecondary(secTrack);
*/

#ifdef MAKESTAT
      // if we need process statisctics
      G4VProcess *proc = (*(atrack.GetParticleDefinition()->GetProcessManager()->GetProcessList()))[1];
      G4String *strp = const_cast<G4String *>(&(proc->GetProcessName()));
      *strp = TPartIndex::I()->ProcName(9); // annihilation
#endif

      return;
    } else {
      return;
    }
  }

  // see if the particle has nuclear capture at rest and decay it if it doesn't.
  if (!fHasNCaptureAtRest[partindex]) { // decay
    isSurv = fDecay->SampleDecay(partindex, nSecPart, pid, mom);
#ifdef MAKESTAT
    // if we need process statisctics
    G4VProcess *proc = (*(atrack.GetParticleDefinition()->GetProcessManager()->GetProcessList()))[1];
    G4String *strp = const_cast<G4String *>(&(proc->GetProcessName()));
    *strp = TPartIndex::I()->ProcName(3); // decay
#endif
  } else { // nuclear capture at rest
    // sample one of the elements of the current material based on the relative
    // number atoms/volume
    elementIndex = fMatXsec[imat]->SampleElement();
    // get the corresponding element-wise final state pointer
    elemfstate = fElemFstate[elementIndex];
    double randn = G4UniformRand();
    isSurv = elemfstate->SampleRestCaptFstate(partindex, nSecPart, weightFst, kerma, energyFst, pid, mom, randn);
#ifdef MAKESTAT
    // if we need process statisctics
    G4VProcess *proc = (*(atrack.GetParticleDefinition()->GetProcessManager()->GetProcessList()))[1];
    G4String *strp = const_cast<G4String *>(&(proc->GetProcessName()));
    *strp = TPartIndex::I()->ProcName(6); // RestCapture
#endif
  }

  // Handle the final state !

  double randDirX = 0.;
  double randDirY = 0.;
  double randDirZ = 1.;
  double randSinTheta = 0.;
  double randPhi = 0.;

  if (nSecPart) {
    randDirZ = 1.0 - 2.0 * G4UniformRand();
    randSinTheta = std::sqrt(1.0 - randDirZ * randDirZ);
    randPhi = G4UniformRand() * CLHEP::twopi;
    randDirX = randSinTheta * std::cos(randPhi);
    randDirY = randSinTheta * std::sin(randPhi);
  }

  // isSurv should always be FALSE here because primary was stopped
  if (isSurv)
    std::cout << "\n---A stopped particle survived its rest process!!!---\n"
              << "---In TabulatedDataManager::SampleFinalStateAtRest---\n" << std::endl;

  // deposited energy goes to totEdepo
  totEdepo = kerma;

  // Go for the secondaries
  std::vector<G4Track *> secTracks;
  G4int totalNumSec = 0;
  G4double time = atrack.GetGlobalTime();
  G4double weight = particlechange->GetParentWeight();

  if (nSecPart && (fgVerboseLevel >= 2))
    std::cout << "============SECONDARIES=FROM=TAB.=PHYS===AT-REST============" << std::endl;

  for (int isec = 0; isec < nSecPart; ++isec) {
    if (pid[isec] >= TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
      int idummy = pid[isec] - 1000000000;
      int Z = idummy / 10000.;
      int A = (idummy - Z * 10000) / 10.;
      double secMass = TPartIndex::I()->GetAprxNuclearMass(Z, A);
      double px = mom[3 * isec];
      double py = mom[3 * isec + 1];
      double pz = mom[3 * isec + 2];
      double secPtot2 = px * px + py * py + pz * pz; // total P^2 [GeV^2]
      totEdepo += std::sqrt(secPtot2 + secMass * secMass) - secMass;
      //#ifdef MAKESTAT
      //       ++killedTracks;
      //#endif
      continue;
    }

    int secPDG = TPartIndex::I()->PDG(pid[isec]); // GV part.code -> PGD code
    TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
    double secMass = secPartPDG->Mass(); // mass [GeV]
    double px = mom[3 * isec];
    double py = mom[3 * isec + 1];
    double pz = mom[3 * isec + 2];
    double secPtot2 = px * px + py * py + pz * pz; // total P^2 [GeV^2]

    double secPtot = std::sqrt(secPtot2);                     // total P [GeV]
    double secEtot = std::sqrt(secPtot2 + secMass * secMass); // total E [GeV]
    double secEkin = secEtot - secMass;                       // kinetic energy in [GeV]

    // Check if Ekin of this secondary is above our energy limit and add it to
    // the list of secondary tracks.
    // Otherwise: Ekin of this secondary goes to energy deposit and check if
    // this secondary has NuclearCaptureAtRest process:
    // If it does has: insert into the list of secondary tracks with a status
    // of StopButAlive and Ekin = 0.0 (momentum direction is unimoprtant)).
    // If it doesn't have: there is nothing else to do! (drink a coffee!)
    if (secEkin > energylimit) {
      ++totalNumSec;

      G4ParticleDefinition *particleDef = G4ParticleTable::GetParticleTable()->FindParticle(secPDG);

      G4ThreeVector newDir(px / secPtot, py / secPtot, pz / secPtot);
      RotateNewTrack(randDirX, randDirY, randDirZ, newDir);

      G4DynamicParticle *dynamicParticle = new G4DynamicParticle(particleDef, newDir, secEkin * GeV);

      G4Track *secTrack = new G4Track(dynamicParticle, time, atrack.GetPosition());
      secTrack->SetKineticEnergy(secEkin * GeV);
      secTrack->SetTrackStatus(fAlive);
      secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
      secTrack->SetWeight(weight);

      secTracks.push_back(secTrack);
      if (fgVerboseLevel >= 2)
        std::cout << "###  " << totalNumSec << "-th SECONDARY:"
                  << " \n"
                  << "***  Secondary is a : = " << TPartIndex::I()->PartName(pid[isec]) << " \n"
                  << "***  Momentum dir.  : = [ " << newDir[0] << ", " << newDir[1] << ", " << newDir[2] << " ]"
                  << " \n"
                  << "***  Kinetic Energy : = " << secEkin << " [GeV]"
                  << " \n"
                  << "***  ______________________________________________________" << std::endl;

    } else {
      totEdepo += secEkin;
      // THIS CAN CAUSE PROBLEMS IF THE AT REST OR nCAPTURE FINAL STATE CONTAINS RECURSION !!!
      // BUT THAT IS A BUG IN PHYSCICS THEN !!!
      if (secEkin <= 0.0 && HasRestProcess(pid[isec])) {
        ++totalNumSec;

        G4ParticleDefinition *particleDef = G4ParticleTable::GetParticleTable()->FindParticle(secPDG);

        G4ThreeVector newDir(0, 0, 1); // not important since Ekin = 0.;
        G4DynamicParticle *dynamicParticle = new G4DynamicParticle(particleDef, newDir, 0.0);

        G4Track *secTrack = new G4Track(dynamicParticle, time, atrack.GetPosition());
        secTrack->SetKineticEnergy(0.0);
        secTrack->SetTrackStatus(fStopButAlive);
        secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
        secTrack->SetWeight(weight);

        secTracks.push_back(secTrack);
      } else {
        //#ifdef MAKESTAT
        //         ++killedTracks;
        //#endif
      } // else: do nothing just add secEkin to Edepo that is already done
    }
  } // end-loop over secondaries

  if (nSecPart && (fgVerboseLevel >= 2))
    std::cout << "============================================================\n"
              << "---Total energy deposit := " << totEdepo << " [Gev] \n"
              << "============================================================\n" << std::endl;

  particlechange->SetNumberOfSecondaries(totalNumSec);
  for (G4int i = 0; i < totalNumSec; ++i)
    particlechange->AddSecondary(secTracks[i]);

  // Set the overall energy deposit
  particlechange->ProposeLocalEnergyDeposit(totEdepo * GeV); // from GeV->MeV
}

void TabulatedDataManager::SampleDecayInFlight(const int partindex, const G4Track &atrack,
                                               G4ParticleChange *particlechange, double energylimit) {
  int nSecPart = 0;     // number of secondary particles per reaction
  const int *pid = 0;   // GeantV particle codes [nSecPart]
  const float *mom = 0; // momentum vectors the secondaries [3*nSecPart]
  char isSurv = 0;      // is the primary survived the interaction

  isSurv = fDecay->SampleDecay(partindex, nSecPart, pid, mom);
  // isSurv should always be FALSE here because primary was stopped
  if (isSurv)
    std::cout << "\n---       A particle survived its decay!!!       ---\n"
              << "---In TabulatedDataManager::SampleFinalStateAtRest---\n" << std::endl;

  if (nSecPart) {
    // Go for the secondaries
    std::vector<G4Track *> secTracks;
    G4int totalNumSec = 0;
    G4double time = atrack.GetGlobalTime();
    G4double weight = particlechange->GetParentWeight();
    G4double totEdepo = 0.0;

    G4double beta = atrack.GetDynamicParticle()->GetTotalMomentum() / atrack.GetTotalEnergy();
    G4double bx = atrack.GetMomentumDirection().x() * beta;
    G4double by = atrack.GetMomentumDirection().y() * beta;
    G4double bz = atrack.GetMomentumDirection().z() * beta;
    G4double b2 = bx * bx + by * by + bz * bz; // it is beta*beta
    G4double gam = 1.0 / std::sqrt(1.0 - b2);
    G4double gam2 = b2 > 0.0 ? (gam - 1.0) / b2 : 0.0;

    for (G4int isec = 0; isec < nSecPart; ++isec) {
      if (pid[isec] >= TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
        int idummy = pid[isec] - 1000000000;
        int Z = idummy / 10000.;
        int A = (idummy - Z * 10000) / 10.;
        double secMass = TPartIndex::I()->GetAprxNuclearMass(Z, A);
        double px = mom[3 * isec];
        double py = mom[3 * isec + 1];
        double pz = mom[3 * isec + 2];
        double secPtot2 = px * px + py * py + pz * pz; // total P^2 [GeV^2]
        totEdepo += std::sqrt(secPtot2 + secMass * secMass) - secMass;
        //#ifdef MAKESTAT
        //         ++killedTracks;
        //#endif
        continue;
      }

      int secPDG = TPartIndex::I()->PDG(pid[isec]); // GV part.code -> PGD code
      TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
      double secMass = secPartPDG->Mass(); // mass [GeV]
      double px = mom[3 * isec];
      double py = mom[3 * isec + 1];
      double pz = mom[3 * isec + 2];
      // compute corrected P^2 in [GeV^2]
      double secP2 = px * px + py * py + pz * pz;            // total P^2 [GeV^2]
      double secEtot = std::sqrt(secP2 + secMass * secMass); // total E [GeV]
      // double secEkin  = secEtot - secMass; //kinetic energy in [GeV]

      G4double bp = bx * px + by * py + bz * pz;
      px = px + gam2 * bp * bx + gam * bx * secEtot;
      py = py + gam2 * bp * by + gam * by * secEtot;
      pz = pz + gam2 * bp * bz + gam * bz * secEtot;
      secEtot = gam * (secEtot + bp);

      G4double secPtot = std::sqrt((secEtot - secMass) * (secEtot + secMass));
      G4double secEkin = secEtot - secMass;
      if (secEkin > energylimit) {
        G4ParticleDefinition *particleDef = G4ParticleTable::GetParticleTable()->FindParticle(secPDG);
        G4ThreeVector newDir(px / secPtot, py / secPtot, pz / secPtot);

        G4DynamicParticle *dynamicParticle = new G4DynamicParticle(particleDef, newDir, secEkin * GeV);

        G4Track *secTrack = new G4Track(dynamicParticle, time, atrack.GetPosition());
        secTrack->SetKineticEnergy(secEkin * GeV);
        secTrack->SetTrackStatus(fAlive);
        secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
        secTrack->SetWeight(weight);
        secTracks.push_back(secTrack);
        ++totalNumSec;
      } else {
        totEdepo += secEkin;
        if (secEkin <= 0.0 && HasRestProcess(pid[isec])) {
          ++totalNumSec;

          G4ParticleDefinition *particleDef = G4ParticleTable::GetParticleTable()->FindParticle(secPDG);

          G4ThreeVector newDir(0, 0, 1); // not important since Ekin = 0.;
          G4DynamicParticle *dynamicParticle = new G4DynamicParticle(particleDef, newDir, 0.0);

          G4Track *secTrack = new G4Track(dynamicParticle, time, atrack.GetPosition());
          secTrack->SetKineticEnergy(0.0);
          secTrack->SetTrackStatus(fStopButAlive);
          secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
          secTrack->SetWeight(weight);

          secTracks.push_back(secTrack);
        } else {
          //#ifdef MAKESTAT
          //           ++killedTracks;
          //#endif
        } // else: do nothing just add secEkin to Edepo that is already done
      }
    }

    particlechange->SetNumberOfSecondaries(totalNumSec);
    for (G4int i = 0; i < totalNumSec; ++i)
      particlechange->AddSecondary(secTracks[i]);

    // Set the overall energy deposit
    particlechange->ProposeLocalEnergyDeposit(totEdepo * GeV); // from GeV->MeV
  }
}

//______________________________________________________________________________
// Simple Multiple SCattering mode: change direction with a mean Theta (exctarct-
// ed by using tabxsec) and a random Phi in [0,2PI]
void TabulatedDataManager::ApplyMsc(G4int imat, const G4Track &atrack) {
  G4int partIndex;    // GV particle index
  G4double kinEnergy; // kinetic energy of the particle in GeV

  partIndex = TPartIndex::I()->PartIndex(atrack.GetParticleDefinition()->GetPDGEncoding());
  kinEnergy = atrack.GetKineticEnergy() / CLHEP::GeV; // from MeV->GeV

  if (imat < 0 || imat >= fNmaterials) {
    std::cout << "\n!!********APPLY-MSC*****USING*GV-TABPHYS***************!!\n"
              << "***  Particle is       = " << TPartIndex::I()->PartName(partIndex) << " \n"
              << "***  Particle E kin.   = " << kinEnergy << " [GeV]"
              << " \n"
              << "***  Material index:   = " << imat << " OUT OF RANGE: [ 0 ," << fNmaterials << " ]!"
              << " \n"
              << "!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << std::endl;
    Fatal("TabulatedDataManager::ApplyMsc", "Material index is out of range!");
  }

  double msTheta = fMatXsec[imat]->MS(partIndex, kinEnergy);
  double msPhi = G4UniformRand() * CLHEP::twopi;

  G4StepPoint *preStepPoint = atrack.GetStep()->GetPreStepPoint();

  // very hard, don't like to do something like this but now we need this
  RotateTrack(const_cast<G4ThreeVector &>(preStepPoint->GetMomentumDirection()), msTheta, msPhi);

  (const_cast<G4Track &>(atrack)).SetMomentumDirection(preStepPoint->GetMomentumDirection());
}

//_____________________________________________________________________________
// (oldxdir, oldydir, oldzdir) are the direction vector of parent track in lab.
// frame; direction vector of the current track, measured from local Z is in the
// G4ThreeVector &newDir; here we rotate it to lab. frame
void TabulatedDataManager::RotateNewTrack(double oldxdir, double oldydir, double oldzdir, G4ThreeVector &newdir) {
  const double one = 1.0;
  const double zero = 0.0;
  const double amin = 1.0e-10;
  const double one5 = 1.5;
  const double half = 0.5;

  double cosTheta0 = oldzdir;
  double sinTheta0 = std::sqrt(oldxdir * oldxdir + oldydir * oldydir);
  double cosPhi0;
  double sinPhi0;

  if (sinTheta0 > amin) {
    cosPhi0 = oldxdir / sinTheta0;
    sinPhi0 = oldydir / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = newdir.x();
  double h1 = sinTheta0 * newdir.z() + cosTheta0 * h0;
  double h2 = newdir.y();

  newdir.setX(h1 * cosPhi0 - h2 * sinPhi0);
  newdir.setY(h1 * sinPhi0 + h2 * cosPhi0);
  newdir.setZ(newdir.z() * cosTheta0 - h0 * sinTheta0);

  // renormalization: avoid 1/sqrt(x) computation by using the 1-th order
  //     Taylor aprx. around 1.0 that should be almost exact since the vector
  //     almost normalized!
  double delta = one5 - half * (newdir.x() * newdir.x() + newdir.y() * newdir.y() + +newdir.z() * newdir.z());
  newdir.setX(newdir.x() * delta);
  newdir.setY(newdir.y() * delta);
  newdir.setZ(newdir.z() * delta);
}

//______________________________________________________________________________
// G4ThreeVector newdir is the direction vector of the primary before scattering
// in lab frame measured from global Z; theta and phi are the scattering angles
// measured from the local Z frame.
void TabulatedDataManager::RotateTrack(G4ThreeVector &newdir, double theta, double phi) {
  const double one = 1.0;
  const double zero = 0.0;
  const double amin = 1.0e-10;
  const double one5 = 1.5;
  const double half = 0.5;

  double cosTheta0 = newdir.z(); // tracks.fZdirV[itrack];
  double sinTheta0 = std::sqrt(newdir.x() * newdir.x() + newdir.y() * newdir.y());
  double cosPhi0;
  double sinPhi0;
  double cosTheta = std::cos(theta);
  double sinTheta = std::sin(theta);

  if (sinTheta0 > amin) {
    cosPhi0 = newdir.x() / sinTheta0;
    sinPhi0 = newdir.y() / sinTheta0;
  } else {
    cosPhi0 = one;
    sinPhi0 = zero;
  }

  double h0 = sinTheta * std::cos(phi);
  double h1 = sinTheta0 * cosTheta + cosTheta0 * h0;
  double h2 = sinTheta * std::sin(phi);

  newdir.setX(h1 * cosPhi0 - h2 * sinPhi0);
  newdir.setY(h1 * sinPhi0 + h2 * cosPhi0);
  newdir.setZ(cosTheta * cosTheta0 - h0 * sinTheta0);

  // renormalization: -ensure normality to avoid accumulated numerical errors
  //    due to sequential calls of rotation; avoid 1/sqrt(x) computation by
  //    using the 1-th order Taylor aprx. around 1.0 that should be almost
  //    exact since the vector almost normalized!
  double delta = one5 - half * (newdir.x() * newdir.x() + newdir.y() * newdir.y() + newdir.z() * newdir.z());
  newdir.setX(newdir.x() * delta);
  newdir.setY(newdir.y() * delta);
  newdir.setZ(newdir.z() * delta);
}

//______________________________________________________________________________
char *TabulatedDataManager::GetVersion() {
  char *ver = new char[512];
  sprintf(ver, "%d.%d.%d", VersionMajor(), VersionMinor(), VersionSub());
  return ver;
}

G4bool TabulatedDataManager::HasRestProcess(int gvindex) {
  return fDecay->HasDecay(gvindex) || fHasNCaptureAtRest[gvindex] || (gvindex == TPartIndex::I()->GetSpecGVIndex(1));
}
