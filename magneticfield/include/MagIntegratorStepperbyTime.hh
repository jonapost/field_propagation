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

#ifndef MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_
#define MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_


#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorStepper.hh"

#include "G4ThreeVector.hh"

#include <assert.h>
#include <iostream>
using namespace std;

//G4MagIntegratorStepper

template <class BaseStepper>

class MagIntegratorStepper_byTime : public BaseStepper {
public:
   inline
   MagIntegratorStepper_byTime(G4Mag_EqRhs *EquationRhs,
                                 G4int numberOfVariables = 6,
                                 G4int numStateVariables = 12);
   virtual ~MagIntegratorStepper_byTime();

   inline void ComputeRightHandSide(const G4double yInput[], G4double dydx[]);



   inline void Stepper(const G4double yInput[],
                       const G4double dydx[],
                             G4double hstep,
                             G4double yOutput[],
                             G4double yError [] );

   //inline void SetTrue_last_step_succeeded();

   //inline void Reset_last_step_succeeded();

   inline G4double DistChord() const;
   //inline G4int IntegratorOrder() const;


   //G4double mass, inv_mass;   // G4MagIntegratorStepper doesn't have these fields by default
                              // only public for quickness!
   //BaseStepper *baseStepper;
   //G4MagIntegratorStepper *baseStepper;

private:
   G4double yIn[10], dydx_copy[10], cached_dydx[10], last_function_evaluation[10];
   //G4double mass, inv_mass; // G4MagIntegratorStepper doesn't have these fields by default

   // G4bool last_step_succeeded;
   // last_step_succeeded now inherits from G4MagIntegratorStepper

   G4Mag_EqRhs  *m_fEq;

};


template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::ComputeRightHandSide(const G4double yInput[], G4double dydx[]) {
   /*if ( BaseStepper::get_last_step_succeeded() ) {
      dydx = cached_dydx;
      return;
   }
   */

   for (int i = 0; i < 10; i ++)
      yIn[i] = yInput[i];
   for (int i = 3; i < 6; i ++)
      yIn[i] *= 1. / m_fEq -> FMass();

   BaseStepper::ComputeRightHandSide(yIn, dydx);

   // Cache (for possible FSAL):
   //for (int i = 0; i < 6; i ++)
   //   last_function_evaluation[i] = dydx[i];

   for (int i = 3; i < 6; i ++)        // may want to make this from i = 0, (but really
                                       // doesn't matter since we're only using this
                                       // for Nystrom Steppers.
      dydx[i] *= m_fEq -> FMass();
   // Always feed this template class a Mag_UsualEqRhs_IntegrateByTime as EquationRhs
   //for (int i = 3; i < 6; i ++)
   //   dydx[i] *= inv_mass;
}

template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::Stepper(const G4double yInput[],
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

   /*
   while ( ! BaseStepper::mTracker -> check_that_wasnt_disgarded_by_Propagator( yInput, diff ) ) {
      cout << "Was disgarded by Propagator, " << diff[0] << ", " << diff[1] << ", " << diff[2] <<  endl;


      assert( BaseStepper::mTracker -> getBufferLength() >= 2);

      BaseStepper::mTracker -> set_last_time_val_was_accepted( true );
      BaseStepper::mTracker -> get_buffer_ptr() -> pop_back();
   }
   */
   //BaseStepper::mTracker -> check( yInput ); // Check to see last step' output is the same as this step's input. (Only 1st 3 coordinates)
#endif



   // Have to copy because yInput and dydx are constant in the function signature.
   for (int i = 0; i < 10; i ++){
      yIn[i] = yInput[i];
      dydx_copy[i] = dydx[i];
   }
   for (int i = 3; i < 6; i ++){
      yIn[i] *= 1. / m_fEq -> FMass();
      dydx_copy[i] *= 1. / m_fEq -> FMass();
   }

   /////// Getting velocity:
   G4double velocity = G4ThreeVector( yIn[3], yIn[4], yIn[5] ).mag();
   ///////

   // Convert hstep (which is in units of arclength) to hstep / velocity (which is units of time).
   BaseStepper::Stepper( yIn, dydx_copy, hstep / velocity, yOutput, yError );



   for (int i = 3; i < 6; i ++)
      yOutput[i] *= m_fEq -> FMass();

#ifdef TRACKING
   if ( BaseStepper::mTracker -> isArmed() ) {
      // BaseStepper::mTracker -> RecordResultOfStepper(yIn, dydx_copy);

      // nextFunctionEvaluation is done at the right endpoint of the integration step.
      // If this works should make nextFunctionEvaluation a member variable:
      G4double nextFunctionEvaluation[8]; // 8 for safety (only need 6).

      // Want to store velocity, not momentum, so wait to multiply yOutput[3..5] by FMass()
      // until after compute next function evaluation.
      ComputeRightHandSide(yOutput, nextFunctionEvaluation);

      for (int i = 3; i < 6; i ++)
            yOutput[i] /= m_fEq -> FMass();
      for (int i = 3; i < 6; i ++)
         nextFunctionEvaluation[i] /= m_fEq -> FMass();

      const G4CachedMagneticField *myField = (G4CachedMagneticField*)( BaseStepper::GetEquationOfMotion() -> GetFieldObj() );
      G4int no_function_calls = myField -> GetCountCalls();

      BaseStepper::mTracker -> RecordResultOfStepper(yOutput, nextFunctionEvaluation, no_function_calls);

      for (int i = 3; i < 6; i ++)
            yOutput[i] *= m_fEq -> FMass();

      BaseStepper::mTracker -> UnArmTracker();


   }
#endif

   // We are now recording position/momentum coordinates, and step units are in arclength

   //yOutput[7] = yInput[7] + hstep;

}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byTime<BaseStepper>::DistChord() const{
   return BaseStepper::DistChord(); // If BaseStepper uses Richardson
                                    // extrapolation (auxStepper) then this might
                                    // not be producing the right results.
}


template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::MagIntegratorStepper_byTime(
                                                                  G4Mag_EqRhs *EquationRhs,
                                                                  G4int numberOfVariables,
                                                                  G4int numStateVariables)
: BaseStepper( EquationRhs, numberOfVariables )//, numStateVariables)
{

   m_fEq = EquationRhs;
   // baseStepper = new BaseStepper(EquationRhs);  //, numberOfVariables, numStateVariables);
   for (int i = 0; i < 10; i ++)
      yIn[i] = 0.;

}
template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::~MagIntegratorStepper_byTime() {
   // delete baseStepper;

}


#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_ */

