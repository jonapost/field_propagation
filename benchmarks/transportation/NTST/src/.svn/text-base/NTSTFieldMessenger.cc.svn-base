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
// $Id: NTSTFieldMessenger.cc,v 1.2 2007-10-26 09:51:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "NTSTFieldMessenger.hh"

#include "NTSTFieldSetup.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//////////////////////////////////////////////////////////////////////////////

NTSTFieldMessenger::NTSTFieldMessenger(NTSTFieldSetup* pEMfieldSetup)
  :fEMfieldSetup(pEMfieldSetup)
{ 
  NTSTfieldDir = new G4UIdirectory("/field/");
  NTSTfieldDir->SetGuidance("NTST field  control.");

  StepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  StepperCmd->SetGuidance("Select stepper type for magnetic field");
  StepperCmd->SetGuidance("StepperNumber = 0 for G4ExplicitEuler     :  StepperNumber = 1 for G4ImplicitEuler");
  StepperCmd->SetGuidance("StepperNumber = 2 for G4SimpleRunge       :  StepperNumber = 3 for G4SimpleHeum");
  StepperCmd->SetGuidance("StepperNumber = 4 for G4ClassicalRK4      :  StepperNumber = 5 for G4HelixExplicitEuler");
  StepperCmd->SetGuidance("StepperNumber = 6 for G4HelixImplicitEuler:  StepperNumber = 7 for G4HelixSimpleRunge");
  StepperCmd->SetGuidance("StepperNumber = 8 for G4CashKarpRKF45     :  StepperNumber = 9 for G4ExactHelixStepper");
  StepperCmd->SetGuidance("StepperNumber =10 for G4RKG3_Stepper      :  StepperNumber =11 for G4HelixMixedStepper");
  StepperCmd->SetParameterName("choice",true);
  StepperCmd->SetDefaultValue(4);
  StepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  FieldNameCmd = new G4UIcmdWithAnInteger("/field/setFieldType",this);
  FieldNameCmd->SetGuidance("Select type of Magnetic field");
  FieldNameCmd->SetGuidance("Number = 0 for UniformField(1.5 tesla,0,0)");
  FieldNameCmd->SetGuidance("Number = 1 for NTSTGradientField Bz=1.5 tesla+grad=0.0001 in xy plane");
  FieldNameCmd->SetGuidance("Number = 2 for TabulatedSolenoid Field");
  FieldNameCmd->SetGuidance("Number = 3 for Different UniformField, enter your value in src/NTSTFieldSetup.cc ");
  FieldNameCmd->SetParameterName("choice",true);
  FieldNameCmd->SetDefaultValue(4);
  FieldNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  UpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  UpdateCmd->SetGuidance("Update field settings");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
  UpdateStepperCmd = new G4UIcmdWithoutParameter("/field/updateStepper",this);
  UpdateStepperCmd->SetGuidance("Update  only Stepper");
  UpdateStepperCmd->AvailableForStates(G4State_Idle);
      
  GradFieldCmd = new G4UIcmdWithADouble("/field/setGradient",this);  
  GradFieldCmd->SetGuidance("Set Gradient for magnetic field.");
  GradFieldCmd->SetGuidance("This command has to be used only for GRADIENT FIELD");
  GradFieldCmd->SetParameterName("Grad",false);
  GradFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
 
  MinEpsilonCmd = new G4UIcmdWithADouble("/field/setMinEpsilon",this);  
  MinEpsilonCmd->SetGuidance("Set Minimum Epsilon for Propagation accuracy");
  MinEpsilonCmd->SetParameterName("MinEps",false);
  MinEpsilonCmd->AvailableForStates(G4State_Idle); 
 
  MaxEpsilonCmd = new G4UIcmdWithADouble("/field/setMaxEpsilon",this);  
  MaxEpsilonCmd->SetGuidance("Set Maximum Epsilon for Propagation accuracy");
  MaxEpsilonCmd->SetParameterName("MaxEps",false);
  MaxEpsilonCmd->AvailableForStates(G4State_Idle); 
 
  MinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);  
  MinStepCmd->SetGuidance("Define minimal step");
  MinStepCmd->SetGuidance("Magnetic field will be in Z direction.");
  MinStepCmd->SetParameterName("min step",false,false);
  MinStepCmd->SetDefaultUnit("mm");
  MinStepCmd->AvailableForStates(G4State_Idle);  
   
  fieldStat = new G4UIcmdWithoutParameter("/field/getFieldStats",this);
  fieldStat->SetGuidance( "Return number calls to field routine" );
  chordFinderStat = new G4UIcmdWithoutParameter("/field/getChordFinderStats",this);
  chordFinderStat->SetGuidance( "Return ChordFinderStatistics" );  
  
}

///////////////////////////////////////////////////////////////////////////////

NTSTFieldMessenger::~NTSTFieldMessenger()
{
  delete StepperCmd;
  delete FieldNameCmd;
  delete GradFieldCmd;
  delete MinStepCmd;
  delete MinEpsilonCmd;
  delete MaxEpsilonCmd;
  delete NTSTfieldDir;
  delete UpdateCmd;

   
}

////////////////////////////////////////////////////////////////////////////
//
//

void NTSTFieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{ 
  if( command == StepperCmd )
  { 
    fEMfieldSetup->SetStepperType(StepperCmd->GetNewIntValue(newValue));}
  if( command == FieldNameCmd )
  { 
    fEMfieldSetup->SetFieldName(FieldNameCmd->GetNewIntValue(newValue)); 
  }  
  if( command == UpdateCmd )
  { 
    fEMfieldSetup->CreateStepperAndChordFinder(); 
  }

 if( command == UpdateStepperCmd )
  { 
    fEMfieldSetup->SetStepper(); 
  }
  if( command == GradFieldCmd )
  { 
    fEMfieldSetup->SetGrad(GradFieldCmd->GetNewDoubleValue(newValue));
  }
  if( command == MinStepCmd )
  { 
    fEMfieldSetup->SetMinStep(MinStepCmd->GetNewDoubleValue(newValue));
  }

  if( command == MinEpsilonCmd )
  { 
    fEMfieldSetup->SetMinEpsilon(MinEpsilonCmd->GetNewDoubleValue(newValue));
  }
  if( command == MaxEpsilonCmd )
  { 
    fEMfieldSetup->SetMaxEpsilon(MaxEpsilonCmd->GetNewDoubleValue(newValue));
  }

  if (command == fieldStat) fEMfieldSetup->GetFieldCallStats();
  if (command == chordFinderStat) fEMfieldSetup->GetChordFinderStats(); 
}

//
//
/////////////////////////////////////////////////////////////////////////
