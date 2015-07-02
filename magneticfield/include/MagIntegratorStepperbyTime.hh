/*
 * MagIntegratorStepperbyTime.hh
 *
 *  Created on: Jun 20, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_
#define MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_

#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorStepper.hh"

#include <iostream>
#include <assert.h>
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
   G4double mass, inv_mass; // G4MagIntegratorStepper doesn't have these fields by default

   // G4bool last_step_succeeded;
   // last_step_succeeded now inherits from G4MagIntegratorStepper

   G4Mag_EqRhs  *m_fEq;


};

#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_ */

/*
template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::ComputeRightHandSide2(const G4double yInput[], G4double dydx[]) {
   //cout  << "MagIntegratorStepper_byTime: "<< G4ThreeVector( yInput[3], yInput[4], yInput[5] ).mag() << endl;
   baseStepper->ComputeRightHandSide(yInput, dydx);

   // Always feed this template class a Mag_UsualEqRhs_IntegrateByTime as EquationRhs
   //for (int i = 3; i < 6; i ++)
   //   dydx[i] *= inv_mass;
}
*/

template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::ComputeRightHandSide(const G4double yInput[], G4double dydx[]) {
   //cout  << "MagIntegratorStepper_byTime: "<< G4ThreeVector( yInput[3], yInput[4], yInput[5] ).mag() << endl;
   //cout << "inside MagIntegratorStepper_byTime<BaseStepper>:: ComputeRightHandSide()" << endl;

   /*if ( BaseStepper::get_last_step_succeeded() ) {
      dydx = cached_dydx;
      return;
   }
   */

   for (int i = 0; i < 10; i ++)
      yIn[i] = yInput[i];
   for (int i = 3; i < 6; i ++)
      yIn[i] *= inv_mass;

   //baseStepper->ComputeRightHandSide(yIn, dydx);

   BaseStepper::ComputeRightHandSide(yIn, dydx);
   //BaseStepper::ComputeRightHandSide(yIn, dydx);
   //baseStepper -> ComputeRightHandSide(yIn, dydx);

   // Cache:
   for (int i = 0; i < 6; i ++)
      last_function_evaluation[i] = dydx[i];

   for (int i = 3; i < 6; i ++)
      dydx[i] *= mass;
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

   //cout << "inside MagIntegratorStepper_byTime<BaseStepper>:: Stepper()" << endl;

   for (int i = 0; i < 10; i ++){
      yIn[i] = yInput[i];
      dydx_copy[i] = dydx[i];
   }
   for (int i = 3; i < 6; i ++){
      yIn[i] *= inv_mass;
      dydx_copy[i] *= inv_mass;
   }
   //cout << "inside MagIntegratorStepper_byTime<BaseStepper>::Stepper()" << endl;
   //baseStepper->Stepper( yIn, dydx_copy, hstep, yOutput, yError );

   //baseStepper -> Stepper( yIn, dydx_copy, hstep, yOutput, yError );
   BaseStepper::Stepper( yIn, dydx_copy, hstep, yOutput, yError );

   //assert( yOutput[0] == yOutput[0] );

   for (int i = 3; i < 6; i ++)
      yOutput[i] *= mass;

   //BaseStepper::Reset_last_step_succeeded(); // Reset last_step_succeeded to false
}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byTime<BaseStepper>::DistChord() const{
   //return baseStepper->DistChord();
   //return BaseStepper::DistChord();
   //return baseStepper -> DistChord();
   return BaseStepper::DistChord();
}

/*
template <class BaseStepper>
inline
G4int MagIntegratorStepper_byTime<BaseStepper>::IntegratorOrder() const{

   cout << "Inside MagIntegratorStepper_byTime<BaseStepper>::IntegratorOrder()" << endl;

   if ( & BaseStepper::IntegratorOrder != 0 )
      return BaseStepper::IntegratorOrder();
   else
      G4Exception("MagIntegratorStepper_ByTime::IntegratorOrder", "GeomField0003",
                  FatalException, "BaseStepper does not implement IntegratorOrder()");
   //return baseStepper->IntegratorOrder();
   //return baseStepper -> IntergratorOrder();
}
*/

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
   mass = m_fEq -> FMass();
   inv_mass = 1. / mass;
   ///last_step_succeeded = false; // Might want to change??
   //BaseStepper::Reset_last_step_succeeded();

   //baseStepper = static_cast<BaseStepper*> ( static_cast<G4MagIntegratorStepper *>( this ) );

}
template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::~MagIntegratorStepper_byTime() {
   // delete baseStepper;
}

/*template <class BaseStepper>
inline void MagIntegratorStepper_byTime<BaseStepper>::SetTrue_last_step_succeeded() {
   //cached_dydx = &last_function_evaluation; // copy ptr
   for (int i = 0; i < 6; i ++)
         cached_dydx[i] = last_function_evaluation[i];
   BaseStepper::SetTrue_last_step_succeeded();

}
*/
//template <class BaseStepper>
//inline void MagIntegratorStepper_byTime<BaseStepper>::Reset_last_step_succeeded() {
//}


