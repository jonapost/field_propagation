//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: NTSTDetectorMessenger.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 

#include "NTSTDetectorMessenger.hh"

#include "NTSTDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include <assert.h>

namespace {

const std::string CLASSICAL_RK4_STEPPER_NAME = "ClassicalRK4";
const std::string CASH_KARP_STEPPER_NAME = "CashKarp45";
const std::string DORMAND_PRINCE_STEPPER_NAME = "DormandPrince745";

std::map<std::string, NTSTDetectorMessenger::StepperType> STEPPERS =
{
    {CLASSICAL_RK4_STEPPER_NAME, NTSTDetectorMessenger::StepperType::ClassicalRK4},
    {CASH_KARP_STEPPER_NAME, NTSTDetectorMessenger::StepperType::CashKarp},
    {DORMAND_PRINCE_STEPPER_NAME, NTSTDetectorMessenger::StepperType::DormandPrince}
};

const std::string STANDARD_DRIVER_NAME = "G4MagInt_Driver";
const std::string NEW_STANRARD_DRIVER_NAME = "G4IntegrationDriver";
const std::string FSAL_DRIVER_NAME = "FSAL";

std::map<std::string, NTSTDetectorMessenger::DriverType> DRIVERS =
{
    {STANDARD_DRIVER_NAME, NTSTDetectorMessenger::DriverType::G4MagInt_Driver},
    {NEW_STANRARD_DRIVER_NAME, NTSTDetectorMessenger::DriverType::G4IntegrationDriver},
    {FSAL_DRIVER_NAME, NTSTDetectorMessenger::DriverType::G4FSALIntegrationDriver}
};

template <class KeyValueContainer>
const typename KeyValueContainer::mapped_type& safeAt(
    const KeyValueContainer& container,
    const typename KeyValueContainer::key_type& key)
{
    const auto it = container.find(key);
    assert(it != container.end());
    return it->second;
}

} // namespace


NTSTDetectorMessenger::NTSTDetectorMessenger(NTSTDetectorConstruction * NTSTDet):
    NTSTDetector(NTSTDet)
{ 
    NTSTdetDir = new G4UIdirectory("/NTST/");
    NTSTdetDir->SetGuidance("NTST detector control.");

    InputFileNameCmd = new G4UIcmdWithAString("/NTST/setInputFile", this);
    InputFileNameCmd->SetGuidance("Set input file name");
    InputFileNameCmd->SetParameterName("File",true);
    InputFileNameCmd->SetDefaultValue("SVT.dat");

    DisableDet = new G4UIcmdWithAString("/NTST/disable", this);
    DisableDet->SetGuidance("disable detetector");
    DisableDet->SetCandidates("none SVT DCH all");
    DisableDet->SetDefaultValue("none");

    DebugCmd = new G4UIcmdWithAnInteger("/NTST/setDebug",this);
    DebugCmd->SetGuidance("Set debug flag.");
    DebugCmd->SetParameterName("Debug",true);

    NSubLayer = new G4UIcmdWithAnInteger("/NTST/setNSubLayer",this);
    NSubLayer->SetGuidance("Set the number of SVT sublayers.");
    NSubLayer->SetParameterName("NSubLay",true);
    NSubLayer->SetDefaultValue(7);
    NSubLayer->SetRange("NSubLay<8");

#if 0
    MinimumDriverStep
        = new G4UIcmdWithADoubleAndUnit("/NTST/setOuterRadius",this);
    MinimumDriverStep->SetGuidance("Set Minimum Step for ");
    MinimumDriverStep->SetParameterName("MinimumStep",false,false);
    MinimumDriverStep->SetDefaultValue(0.1);
    MinimumDriverStep->SetDefaultUnit("mm");
    MinimumDriverStep->SetRange("Radius>0.");
#endif
  
    MotherOuterRadius = new G4UIcmdWithADoubleAndUnit("/NTST/setOuterRadius",this);
    MotherOuterRadius->SetGuidance("Set outer radius of the SVT mother volume");
    MotherOuterRadius->SetParameterName("Radius",false,false);
    MotherOuterRadius->SetDefaultValue(19.);
    MotherOuterRadius->SetDefaultUnit("cm");
    MotherOuterRadius->SetRange("Radius>0.");
  
    fieldStat = new G4UIcmdWithoutParameter("/NTST/getFieldStats",this);
    fieldStat->SetGuidance( "Return number calls to field routine" );

    StepperMethodCmd = new G4UIcmdWithAString("/NTST/stepper", this);
    StepperMethodCmd->SetGuidance("ClassicalRK4");

    DriverMethodCmd = new G4UIcmdWithAString("/NTST/driver" ,this);
    DriverMethodCmd->SetGuidance("G4MagInt_Driver, G4FSALIntegrationDriver");


}

NTSTDetectorMessenger::~NTSTDetectorMessenger()
{
    delete DebugCmd;
    delete MotherOuterRadius;
    delete StepperMethodCmd;
    delete DriverMethodCmd;
    delete NTSTdetDir;
}

void NTSTDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
    if (command == DebugCmd) {
        NTSTDetector->SetDebugCmd(DebugCmd->GetNewIntValue(newValue));
    }
   
    if (command == MotherOuterRadius) {
        NTSTDetector->SetOuterRadius(
                MotherOuterRadius->GetNewDoubleValue(newValue));
    }

    if (command == NSubLayer) {
        NTSTDetector->SetNSubLayer(
                NSubLayer->GetNewIntValue(newValue));
    }

    if (command == fieldStat) {
        NTSTDetector->GetFieldCallStats();
    }

    if (command == StepperMethodCmd) {
        NTSTDetector->SetStepperMethod(safeAt(STEPPERS, newValue));
    }

    if (command == DriverMethodCmd) {
        NTSTDetector->SetDriverMethod(safeAt(DRIVERS, newValue));
    }
}
