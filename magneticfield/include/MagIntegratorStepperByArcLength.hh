// Nystrom stepper implemenation by Jason Suagee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 27 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
#ifndef MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYARCLENGTH_HH_
#define MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYARCLENGTH_HH_


#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorStepper.hh"

#include <iostream>
#include <assert.h>
using namespace std;

// G4MagIntegratorStepper_byArcLength, (a wrapper class so that CashKarp, ClassicalRK4 and others
// can interface with StepTracker when Stepper() is called, and to store total time).

template <class BaseStepper>

class MagIntegratorStepper_byArcLength : public BaseStepper {
public:

   MagIntegratorStepper_byArcLength(G4Mag_EqRhs *EquationRhs,
                                    G4int numberOfVariables = 6,
                                    G4int numStateVariables = 12);
   virtual ~MagIntegratorStepper_byArcLength();

   inline void ComputeRightHandSide(const G4double yInput[], G4double dydx[]);

   inline void Stepper(const G4double yInput[],
                       const G4double dydx[],
                             G4double hstep,
                             G4double yOutput[],
                             G4double yError [] );

   inline G4double DistChord() const;

private:
   //G4double yIn[10], dydx_copy[10], cached_dydx[10], last_function_evaluation[10];
   //G4double mass, inv_mass; // G4MagIntegratorStepper doesn't have these fields by default

   G4Mag_EqRhs  *m_fEq;

};


template <class BaseStepper>
inline
void MagIntegratorStepper_byArcLength<BaseStepper>::ComputeRightHandSide(const G4double yInput[], G4double dydx[]) {

   BaseStepper::ComputeRightHandSide(yInput, dydx);

   // Cache:
   // for (int i = 0; i < 6; i ++)
   //   last_function_evaluation[i] = dydx[i];

}

template <class BaseStepper>
inline
void MagIntegratorStepper_byArcLength<BaseStepper>::Stepper(const G4double yInput[],
      const G4double dydx[],
            G4double hstep,
            G4double yOutput[],
            G4double yError [] ) {


#ifdef TRACKING
   G4double diff[3];
   G4double last_postion_vector[3], second_to_last_postion_vector[3];

   while ( ! BaseStepper::mTracker -> check_that_wasnt_disgarded_by_Propagator( yInput, diff ) ) {
      cout << "Was disgarded by Propagator, " << diff[0] << ", " << diff[1] << ", " << diff[2] <<  endl;
      assert( BaseStepper::mTracker -> getBufferLength() >= 2);


      BaseStepper::mTracker -> set_last_time_val_was_accepted( true );


      second_to_last_postion_vector[0] = ( BaseStepper::mTracker -> get_buffer_ptr() -> at( BaseStepper::mTracker ->getBufferLength() - 2 ) )
            . at(2);
      second_to_last_postion_vector[1] = ( BaseStepper::mTracker -> get_buffer_ptr() -> at( BaseStepper::mTracker ->getBufferLength() - 2 ) )
            . at(3);
      second_to_last_postion_vector[2] = ( BaseStepper::mTracker -> get_buffer_ptr() -> at( BaseStepper::mTracker ->getBufferLength() - 2 ) )
            . at(4);

      last_postion_vector[0] = ( BaseStepper::mTracker -> get_buffer_ptr() -> at( BaseStepper::mTracker ->getBufferLength() - 1 ) )
            . at(2);
      last_postion_vector[1] = ( BaseStepper::mTracker -> get_buffer_ptr() -> at( BaseStepper::mTracker ->getBufferLength() - 1 ) )
            . at(3);
      last_postion_vector[2] = ( BaseStepper::mTracker -> get_buffer_ptr() -> at( BaseStepper::mTracker ->getBufferLength() - 1 ) )
            . at(4);

      if ( ( G4ThreeVector(yInput[0], yInput[1], yInput[2]) - G4ThreeVector(second_to_last_postion_vector[0],
            second_to_last_postion_vector[1], second_to_last_postion_vector[2]) ).mag() >=
            ( G4ThreeVector(yInput[0], yInput[1], yInput[2]) - G4ThreeVector(last_postion_vector[0],
                        last_postion_vector[1], last_postion_vector[2]) ).mag() ) {

         break;
      }

      BaseStepper::mTracker -> get_buffer_ptr() -> pop_back();
      BaseStepper::mTracker -> get_no_function_calls_buffer() -> pop_back();

   }

   //BaseStepper::mTracker -> check( yInput ); // Check to see last step' output is the same as this step's input. (Only 1st 3 coordinates)
#endif

   BaseStepper::Stepper( yInput, dydx, hstep, yOutput, yError );
   //( dynamic_cast<BaseStepper*>( this )) -> BaseStepper::Stepper( yIn, dydx_copy, hstep, yOutput, yError );

#ifdef TRACKING
   if ( BaseStepper::mTracker -> isArmed() ) {
      // BaseStepper::mTracker -> RecordResultOfStepper(yIn, dydx_copy);

      // nextFunctionEvaluation is done at the right endpoint of the integration step.
      // If this works should make nextFunctionEvaluation a member variable:
      G4double nextFunctionEvaluation[8]; // 8 for safety (only need 6).

      // Want to store velocity, not momentum, so wait to multiply yOutput[3..5] by FMass()
      // until after compute next function evaluation.
      BaseStepper::ComputeRightHandSide(yOutput, nextFunctionEvaluation);


      for (int i = 3; i < 6; i ++)
                  yOutput[i] /= m_fEq -> FMass();
      for (int i = 3; i < 6; i ++)
         nextFunctionEvaluation[i] /= m_fEq -> FMass();


      const G4CachedMagneticField *myField = (G4CachedMagneticField*)( BaseStepper::GetEquationOfMotion() -> GetFieldObj() );
      G4int no_function_calls = myField -> GetCountCalls();

      BaseStepper::mTracker -> RecordResultOfStepper(yOutput, nextFunctionEvaluation, no_function_calls); // Store as velocity (instead of mom.)

      for (int i = 3; i < 6; i ++)
            yOutput[i] *= m_fEq -> FMass();


      //BaseStepper::mTracker -> RecordResultOfStepper(yOutput, nextFunctionEvaluation);

      BaseStepper::mTracker -> UnArmTracker();
   }
#endif

}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byArcLength<BaseStepper>::DistChord() const{
   return BaseStepper::DistChord(); // If BaseStepper uses Richardson
                                    // extrapolation (auxStepper) then this might
                                    // not be producing the right results.
}

template <class BaseStepper>
inline MagIntegratorStepper_byArcLength<BaseStepper>::MagIntegratorStepper_byArcLength(
                                                                  G4Mag_EqRhs *EquationRhs,
                                                                  G4int numberOfVariables,
                                                                  G4int numStateVariables)
: BaseStepper( EquationRhs, numberOfVariables )//, numStateVariables)
{

   // Not really needed? (It is needed for MagIntegratorStepper_byTime because we need FMass(),
   // which is not in G4EquationOfMotion).
   m_fEq = EquationRhs;
}

template <class BaseStepper>
MagIntegratorStepper_byArcLength<BaseStepper>::~MagIntegratorStepper_byArcLength() { }


#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYARCLENGTH_HH_ */
