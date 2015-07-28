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

#define NO_STATE_VARIABLES 12


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

   G4double nextFunctionEvaluation[NO_STATE_VARIABLES]; // 8 for safety (only need 6).

   G4Mag_EqRhs  *m_fEq;
};


template <class BaseStepper>
inline
void MagIntegratorStepper_byArcLength<BaseStepper>::ComputeRightHandSide(const G4double yInput[], G4double dydx[]) {

   BaseStepper::ComputeRightHandSide(yInput, dydx);

}

template <class BaseStepper>
inline
void MagIntegratorStepper_byArcLength<BaseStepper>::Stepper(const G4double yInput[],
      const G4double dydx[],
            G4double hstep,
            G4double yOutput[],
            G4double yError [] ) {


#ifdef TRACKING
   std::vector<G4double> *last_pos_vector, *second_to_last_pos_vector;
   G4int buffer_length;

   while ( ! BaseStepper::mTracker -> check_that_wasnt_disgarded_by_Propagator( yInput ) ) {

      buffer_length = BaseStepper::mTracker -> getBufferLength();
      // Because of current bug with BogackiShampine45 (and possibly others):
      assert( buffer_length >= 2);

      BaseStepper::mTracker -> set_last_time_val_was_accepted( true );

      last_pos_vector =
                        &( BaseStepper::mTracker -> get_buffer_ptr() -> at( buffer_length - 1 ) );
      second_to_last_pos_vector =
                        &( BaseStepper::mTracker -> get_buffer_ptr() -> at( buffer_length - 2 ) );

      if (     ( G4ThreeVector(yInput[0], yInput[1], yInput[2])
                 - G4ThreeVector(second_to_last_pos_vector -> at(POSITION_SLOT + 0),
                                 second_to_last_pos_vector -> at(POSITION_SLOT + 1),
                                 second_to_last_pos_vector -> at(POSITION_SLOT + 2)) ).mag()
            >=
               ( G4ThreeVector(yInput[0], yInput[1], yInput[2])
                 - G4ThreeVector(last_pos_vector -> at(POSITION_SLOT + 0),
                                 last_pos_vector -> at(POSITION_SLOT + 1),
                                 last_pos_vector -> at(POSITION_SLOT + 2)) ).mag()     )
      {
         break; // We are done searching backwards for a good match.
      }

      BaseStepper::mTracker -> get_buffer_ptr() -> pop_back();
      BaseStepper::mTracker -> get_no_function_calls_buffer() -> pop_back();

   }
#endif

   BaseStepper::Stepper( yInput, dydx, hstep, yOutput, yError );
   //( dynamic_cast<BaseStepper*>( this )) -> BaseStepper::Stepper( yIn, dydx_copy, hstep, yOutput, yError );

#ifdef TRACKING
   if ( BaseStepper::mTracker -> isArmed() ) {
      // BaseStepper::mTracker -> RecordResultOfStepper(yIn, dydx_copy);

      // nextFunctionEvaluation is done at the right endpoint of the integration step.
      // If this works should make nextFunctionEvaluation a member variable:

      BaseStepper::ComputeRightHandSide(yOutput, nextFunctionEvaluation);

      // Want to store velocity, not momentum.
      for (int i = 3; i < 6; i ++)
                  yOutput[i] /= m_fEq -> FMass();
      // Divide by mass to get acceleration:
      for (int i = 3; i < 6; i ++)
         nextFunctionEvaluation[i] /= m_fEq -> FMass();


      const G4CachedMagneticField *myField = (G4CachedMagneticField*)( BaseStepper::GetEquationOfMotion() -> GetFieldObj() );
      G4int no_function_calls = myField -> GetCountCalls();

      BaseStepper::mTracker -> RecordResultOfStepper(yOutput, nextFunctionEvaluation, no_function_calls); // Store as velocity (instead of mom.)

      // Change back to momentum coordinates for next round of stepper:
      for (int i = 3; i < 6; i ++)
            yOutput[i] *= m_fEq -> FMass();

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
