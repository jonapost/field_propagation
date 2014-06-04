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
// $Id: G4EquationOfMotion.hh 71664 2013-06-20 08:36:05Z gcosmo $
//
//
// class G4EquationOfMotion
//
// Class description:
//
// Abstract Base Class for the right hand size of the equation of
// motion of a particle in a field.

// History:
// - Created. J.Apostolakis
// -------------------------------------------------------------------

#ifndef G4_EquationOfMotion_DEF
#define G4_EquationOfMotion_DEF

#include "G4Types.hh"      // "globals.hh"
#include "G4Field.hh"   // required in inline method implementations

#include "G4ChargeState.hh"
template
<class State, class Field>
class G4EquationOfMotion 
{
  public:  // with description

     G4EquationOfMotion( Field *Field );
     virtual ~G4EquationOfMotion();
       // Constructor and virtual destructor. No operations.

     virtual void EvaluateRhsGivenB( State s, const  G4double B[3]) const = 0;
       // Given the value of the  field "B", this function 
       // calculates the value of the derivative dydx.
       // --------------------------------------------------------
       // This is the _only_ function a subclass must define.
       // The other two functions use Rhs_givenB.

     virtual void SetChargeMomentumMass(G4ChargeState particleCharge,
                                        G4double MomentumXc,
                                        G4double MassXc2) = 0;
       // Set the charge, momentum and mass of the current particle
       // --> used to set the equation's coefficients ...

     inline
     void RightHandSide(State s ) const;
       // This calculates the value of the derivative dydx at y.
       // It is the usual enquiry function.
       // ---------------------------
       // (It is not virtual, but calls the virtual function above.)

     void EvaluateRhsReturnB(State s, G4double Field[]  ) const;
       // Same as RHS above, but also returns the value of B.
       // Should be made the new default ? after putting dydx & B in a class.

     void GetFieldValue( const  G4double Point[4],
                                G4double Field[] )  const;
       // Obtain only the field - the stepper assumes it is pure Magnetic.
       // Not protected, because G4RKG3_Stepper uses it directly.

     const Field* GetFieldObj() const;
     void           SetFieldObj(Field* pField);

  private:
     // const int G4maximum_number_of_field_components = 24;
     enum { G4maximum_number_of_field_components = 24 } ;

     Field *itsField;

};

#include "G4EquationOfMotion.icc"

#endif /* G4_EquationOfMotion_DEF */
