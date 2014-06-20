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
//
// $Id: NTSTFieldSetup.cc,v 1.2 2007-10-26 09:51:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   User Field setup class implementation.
//

#include "NTSTFieldSetup.hh"
#include "NTSTFieldMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ExactHelixStepper.hh"
#include "G4HelixMixedStepper.hh"
#include "G4NystromRK4.hh"
#include "NTSTField.hh"
//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

bool fieldFlag = true;

    NTSTFieldSetup::NTSTFieldSetup(G4MagneticField *pCommonField)
: fChordFinder(0),fEquation(0), fMagneticField(0),
    pAField1(0),pAField2(0),
    ffield(0),
    fStepper(0),fStepperType(4),
    fFieldName(0),fMinStep(0.01),fGradofField(0.000001)
{
    fMagneticField = pCommonField; 
    fFieldMessenger = new NTSTFieldMessenger(this) ; 

}

    NTSTFieldSetup::NTSTFieldSetup()
:  fChordFinder(0),fEquation(0), fMagneticField(0),
    pAField1(0),pAField2(0),
    ffield(0),
    fStepper(0),fStepperType(4),
    fFieldName(0),fMinStep(0.01),fGradofField(0.000001)
{
    fMagneticField = new G4UniformMagField( G4ThreeVector(0.0, 0.0, 0.0 ) );
    G4cout << " NTSTFieldSetup: magnetic field set to Uniform( 0.0, 0, 0 ) " << G4endl;
    InitialiseAll();
    fieldFlag = false;
}

    void
NTSTFieldSetup::InitialiseAll()
{


    fMinStep = 1.0*mm ; // minimal step of 1 mm is default
    fMaxEpsilon= 0.00001;
    fMinEpsilon= 0.001;
    fFieldManager = G4TransportationManager::GetTransportationManager()
        ->GetFieldManager();
    CreateStepperAndChordFinder();
}

////////////////////////////////////////////////////////////////////////////////

NTSTFieldSetup::~NTSTFieldSetup()
{
    // GetGlobalFieldManager()->SetDetectorField(0);
    //G4cout<<" Deleting NTSTFieldSetup"<<G4endl;
    //if(fMagneticField) delete fMagneticField;
    if(fChordFinder)   delete fChordFinder;
    if(fStepper)       delete fStepper;
    //G4cout<<"End of Deleting NTSTFieldSetup"<<G4endl;
}

/////////////////////////////////////////////////////////////////////////////
//
// Update field
//
#include "G4MagIntegratorDriver.hh"
void NTSTFieldSetup::CreateStepperAndChordFinder()
{
    SetField();  
    if(fEquation) delete fEquation;
    fEquation = new G4Mag_UsualEqRhs(fMagneticField);

    SetStepper();
    G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

    fFieldManager->SetDetectorField(fMagneticField );

    if(fChordFinder) delete fChordFinder;

    fChordFinder = new G4ChordFinder( fMagneticField, fMinStep,fStepper);

    fFieldManager->SetChordFinder( fChordFinder );

    fFieldManager->SetMinimumEpsilonStep( fMinEpsilon );    // Old value
    fFieldManager->SetMaximumEpsilonStep( fMaxEpsilon );      // FIX - old value
    fFieldManager->SetDeltaOneStep( 0.25 * mm );       // original value
    fFieldManager->SetDeltaIntersection( 0.10 * mm );  // original value

    G4cout << "Field Manager's parameters are " 
        << " minEpsilonStep= " << fFieldManager->GetMinimumEpsilonStep() << " "
        << " maxEpsilonStep= " << fFieldManager->GetMaximumEpsilonStep() << " " 
        << " deltaOneStep=   " << fFieldManager->GetDeltaOneStep() << " "
        << " deltaIntersection= " << fFieldManager->GetDeltaIntersection() 
        << G4endl;
    //G4MagInt_Driver *pDriver;
    //To have verbose from MagInt_Driver
    //fChordFinder->SetVerbose(1);
    //  pDriver=fpChordFinder->GetIntegratorDriver();
    //pDriver=new G4MagInt_Driver(fMinStep, 
    //                                   fStepper, 
    //                                   fStepper->GetNumberOfVariables(),2 );

    //fChordFinder->SetIntegrationDriver(pDriver);
    //fFieldManager->SetChordFinder( fpChordFinder );
    return;
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

//=============test template mode================

#include "TUniformMagField.hh"
#include "TMagFieldEquation.hh"
#include "TCashKarpRKF45.hh"
#include "TClassicalRK4.hh"

typedef NTSTField Field_t1;
typedef TUniformMagField Field_t2;
typedef TMagFieldEquation<Field_t1> Equation_t1;
typedef TMagFieldEquation<Field_t2> Equation_t2;
typedef TCashKarpRKF45<Equation_t1, 6> Stepper_t1;
typedef TCashKarpRKF45<Equation_t2, 6> Stepper_t2;
typedef TClassicalRK4<Equation_t1, 6> StepperRK4_t1;
typedef TClassicalRK4<Equation_t2, 6> StepperRK4_t2;
//===============================================

void NTSTFieldSetup::SetStepper()
{
    switch ( fStepperType ) 
    {
        case 0:  
            fStepper = new G4ExplicitEuler( fEquation ); 
            G4cout<<"G4ExplicitEuler is called"<<G4endl;     
            break;
        case 1:  
            fStepper = new G4ImplicitEuler( fEquation );      
            G4cout<<"G4ImplicitEuler is called"<<G4endl;     
            break;
        case 2:  
            fStepper = new G4SimpleRunge( fEquation );        
            G4cout<<"G4SimpleRunge is called"<<G4endl;     
            break;
        case 3:  
            fStepper = new G4SimpleHeum( fEquation );         
            G4cout<<"G4SimpleHeum is called"<<G4endl;     
            break;
        case 4:  
            fStepper = new G4ClassicalRK4( fEquation );       
            G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;     
            break;
        case 5:  
            fStepper = new G4HelixExplicitEuler( fEquation ); 
            G4cout<<"G4HelixExplicitEuler is called"<<G4endl;     
            break;
        case 6:  
            fStepper = new G4HelixImplicitEuler( fEquation ); 
            G4cout<<"G4HelixImplicitEuler is called"<<G4endl;     
            break;
        case 7:  
            fStepper = new G4HelixSimpleRunge( fEquation );   
            G4cout<<"G4HelixSimpleRunge is called"<<G4endl;     
            break;
        case 8:  
            fStepper = new G4CashKarpRKF45( fEquation );      
            G4cout<<"G4CashKarpRKF45 is called"<<G4endl;     
            break;
        case 9:  
            fStepper = new G4ExactHelixStepper( fEquation );       
            G4cout<<"G4ExactHelixStepper is called"<<G4endl;     
            break;
        case 10:  
            fStepper = new G4RKG3_Stepper( fEquation );       
            G4cout<<"G4RKG3_Stepper is called"<<G4endl;     
            break;
        case 11:  
            fStepper = new G4HelixMixedStepper( fEquation );

            G4cout<<"G4HelixMixedStepper is called"<<G4endl;     
            break;
        case 12:  
            fStepper = new G4HelixMixedStepper( fEquation,8 );

            G4cout<<"G4HelixMixedStepper is called"<<G4endl;     
            break;
        case 13:
            fStepper = new G4NystromRK4( fEquation);

            G4cout<<"G4NystromRK4 is called"<<G4endl;
            break;
            //=============test template mode================
        case 14:
            {
                if (fieldFlag)
                {
                    Equation_t1* pEquation = new Equation_t1( static_cast<Field_t1*>(fMagneticField) );
                    fStepper = new  Stepper_t1(pEquation);
                }
                else{

                    Equation_t2* pEquation = new Equation_t2( static_cast<Field_t2*>(fMagneticField) );
                    fStepper = new  Stepper_t2(pEquation);
                }
            }
            G4cout<<"Templated CashKarpRKF45 is called"<<G4endl;
            break;
        case 15:
            {
                if (fieldFlag)
                {
                    Equation_t1* pEquation = new Equation_t1( static_cast<Field_t1*>(fMagneticField) );
                    fStepper = new  StepperRK4_t1(pEquation);
                }
                else{

                    Equation_t2* pEquation = new Equation_t2( static_cast<Field_t2*>(fMagneticField) );
                    fStepper = new  StepperRK4_t2(pEquation);
                }
            }
            G4cout<<"Templated ClassicalRK4 is called"<<G4endl;
            break;
            //===============================================
        default: fStepper = 0;
    }
    return; 
}
#include "NTSTField.hh"
#include "NTSTGradientField.hh"
#include "NTSTTabulatedField3d.hh"
void NTSTFieldSetup::SetField()
{
    switch(fFieldName)
    { 
        case 0: G4cout<<"Field set to  UniformField(default)"<<G4endl;break;

        case 1: if(pAField1)delete pAField1;
                    pAField1= new NTSTGradientField(fGradofField);
                fMagneticField=pAField1;
                G4cout<<"Field set to Gradient Field"<<G4endl;break;
        case 2: pAField2= new NTSTTabulatedField3d("TableST5.dat", 0.0);
                fMagneticField=pAField2;
                G4cout<<"Field set to Tabulated Solenoid Field"<<G4endl;break;

        case 3: ffield= new NTSTField(1.5*tesla,0.,0.);
                fMagneticField=ffield; 
                G4cout<<"Field set to  UniformField with new value"<<G4endl;break;
        default: G4cout<<"Field set to  UniformField(default)"<<G4endl;break;
    };

}


void NTSTFieldSetup::GetChordFinderStats()
{
    if(fChordFinder){
        fChordFinder->PrintStatistics();
    }
    else{
        G4cout<<"Copy of ChordFinder doesn't exist"<<G4endl;
    }

}

void NTSTFieldSetup::GetFieldCallStats()
{
    if(fChordFinder){
        G4cout << "Number calls to field = " ;
        switch(fFieldName)
        { 
            case 0: G4cout<<"use NTST/getFieldStats command"<<G4endl;break;
            case 1: G4cout<< pAField1->GetCount()<<G4endl;pAField1->ClearCount();break;      
            case 2: G4cout<< pAField2->GetCount()<<G4endl;pAField2->ClearCount();break;
            case 3: G4cout<< ffield->GetCount()<<G4endl;ffield->ClearCount();break;
            default:G4cout<<"no field"<<G4endl;break;
        };
    } 
    else{
        G4cout<<"Field doesn't exist"<<G4endl;
    }  
}



////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  NTSTFieldSetup::GetGlobalFieldManager()
{
    return G4TransportationManager::GetTransportationManager()
        ->GetFieldManager();
}
