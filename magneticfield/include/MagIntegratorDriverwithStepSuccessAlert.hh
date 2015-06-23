/*
 * MagIntegratorDriverwithStepSuccessAlert.hh
 *
 *  Created on: Jun 23, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_MAGINTEGRATORDRIVERWITHSTEPSUCCESSALERT_HH_
#define MAGNETICFIELD_INCLUDE_MAGINTEGRATORDRIVERWITHSTEPSUCCESSALERT_HH_

#include "G4MagIntegratorDriver.hh"
#include "MagIntegratorStepperbyTime.hh"

template <class BaseStepper>

class MagInt_Driver_withStepSuccessAlert<BaseStepper> : public G4MagInt_Driver {
public:
   // MagInt_Driver_withStepSuccessAlert();
   virtual ~MagInt_Driver_withStepSuccessAlert() {};

   inline void ReportStepSucceded();

   // No ReportStepFailed() -- For now MagIntegratorStepperByTime will assume a failed step unless told otherwise.
   // This is not the best way to do this since steps will succeed more often than they will fail
   // But it was the easiest to fit into the existing methods coming from G4MagInt_Driver

private:
   MagIntegratorStepperByTime<BaseStepper> *pIntStepper;



};

template <class BaseStepper>
inline void MagInt_Driver_withStepSuccessAlert<BaseStepper>::ReportStepSucceded() {
   pIntStepper -> StepSucceded();
}



#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORDRIVERWITHSTEPSUCCESSALERT_HH_ */
