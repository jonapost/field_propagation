/*
 * ChawlaSharmaWrapper.cc
 *
 *  Created on: Jun 17, 2015
 *      Author: jason
 */

#include "ChawlaSharmaWrapper.hh"
#include "G4SystemOfUnits.hh"


ChawlaSharmaWrapper::~ChawlaSharmaWrapper() {
	// TODO Auto-generated destructor stub
}

void ChawlaSharmaWrapper::Stepper(const G4double yInput[],
        const G4double dydx[],
              G4double hstep,
              G4double yOutput[],
              G4double yError [] ){

	G4double mass = m_fEq->FMass();
	G4double iMass = 1. / mass;

	   const G4int nvar = this->GetNumberOfVariables() ;
	   const G4int maxvar= GetNumberOfStateVariables();

	   G4int i;
	   // correction for Richardson Extrapolation.
	   G4double  correction = 1. / ( (1 << (IntegratorOrder()-1)) -1 ); // IntegratorOrder() - 1 because we
																		//falsely report to higher classes integrator order 4

	   //  Saving yInput because yInput and yOutput can be aliases for same array

	   for(i=0;i<nvar;i++) yInitial[i]=yInput[i];

	   for(i=3;i<nvar;i++) yInitial[i] *= iMass;

	   yInitial[7]= yInput[7];    // Copy the time in case ... even if not really needed
	   yMiddle[7] = yInput[7];  // Copy the time from initial value
	   yOneStep[7] = yInput[7]; // As it contributes to final value of yOutput ?
	   // yOutput[7] = yInput[7];  // -> dumb stepper does it too for RK4
	   for(i=nvar;i<maxvar;i++) yOutput[i]=yInput[i];
	   // yError[7] = 0.0;

	   G4double halfStep = hstep * 0.5;

	   // Comment (*):
	   // Do two half steps

	   // No renormalization yet:
	   //G4double m_ivelocity = 1. / std::sqrt(yInitial[3]*yInitial[3]+yInitial[4]*yInitial[4]+yInitial[5]*yInitial[5]);
	   //for(i = 3; i < 6; i++) yInitial[i] *= m_ivelocity;

	   DumbStepper  (yInitial, halfStep, yMiddle);



	   // Might put in for one of the set of coefficients:
	   // RightHandSide(yMiddle, dydxMid);

	   // Comment (**):
	   // Function evaluations f(x,y,y') are not at (x_k, y_k, y'_k)
	   // but rather at (x_k, y_k + alpha1*halfStep*y'_k, y'_k).
	   // There is no point in calculating dy/dx at yMiddle.
	   // See p. 376 of Chawla/Sharma.

	   DumbStepper  (yMiddle, halfStep, yOutput);
	   for(i=3;i<nvar;i++) yOutput[i] *= mass;

	   // Store midpoint, chord calculation

	   fMidPoint = G4ThreeVector( yMiddle[0],  yMiddle[1],  yMiddle[2]);

	   for(i=3;i<nvar;i++) yMiddle[i] *= mass;


	   // Do a full Step
	   DumbStepper(yInitial, hstep, yOneStep);
	   for(i=3;i<nvar;i++) yOneStep[i] *= mass;


	   for(i=0;i<nvar;i++) {
		  //std::cerr << "i = " << i << " ,yOutput = "  << yOutput[i] << " , yOneStep = " << yOneStep[i] << "\n";;
		  //if (yOutput[i] == yOneStep[i]){ std::cerr << "  were equal \n";} else std::cerr << "\n";
		  yError [i] = (yOutput[i] - yOneStep[i]) * correction ;
		  yOutput[i] += yError[i] ;  // Provides accuracy increased
												// by 1 order via the
												// Richardson Extrapolation.
												// Question: Is the order reported as 3 to MagIntegratorDriver
												// ( at MagIntegratorDriver.cc line 90 )
												// when it should be reported as 3+1 = 4 ?
	   }

	   //for(i=3;i<nvar;i++) yOutput[i] *= mass;


	   // Now renomalize velocityentum
	   //G4double velocityentum_sqrt =1. / std::sqrt( yOutput[3]*yOutput[3] + yOutput[4]*yOutput[4] + yOutput[5]*yOutput[5] );

	   //for(i = 3; i < 6; i ++){
	   //   yOutput[i] *= velocityentum_sqrt;
	   //}

	   fInitialPoint = G4ThreeVector( yInitial[0], yInitial[1], yInitial[2]);
	   fFinalPoint   = G4ThreeVector( yOutput[0],  yOutput[1],  yOutput[2]);

	   //cout << DistChord() << endl;

	   return ;
}

void ChawlaSharmaWrapper::DumbStepper(
										const G4double yIn[],
										// const G4double dydx[],
										G4double Step,
										G4double yOut[] ) {
	// dydx is not even used by DumbStepper. See Comment (**) above.

	G4double mass = m_fEq->FMass();
	G4double iMass = 1. / mass;
	G4int i;

	G4double K1[6], K2[6], K3[6];
	G4double pos[3], velocity[3], t;
	G4double temp_eval_pt[8];
	//G4double m_velocity_sqrd = 0., m_velocity;
	const G4double   a1 = 1./4. , a2 = 0., a3 = 1./4.,
 									b1 = 1./4., b2 = 0., b3 = 3./4.,
 									alpha1 = 0., alpha2 = 2./3., alpha3 = 2./3.,
 									beta21 = -1./9., beta31 = 2./9., beta32 = 0.,
 									gamma21 = 2./3., gamma31 = 1./3., gamma32 = 1./3.;

   /* Another set of coefficients to try:
    *
   const G4double   a1 = 0. , a2 = 1./2., a3 = 0.,
 									b1 = 0., b2 = 3./4., b3 = 1./4.,
 									alpha1 = 1./3., alpha2 = 1./3., alpha3 = 1.,
 									beta21 = 0., beta31 = 1., beta32 = -1./3.,
 									gamma21 = 1./3., gamma31 = -1., gamma32 = 2.;

   */

   // Initialise time to t0, needed when it is not updated by the integration.
   //        [ Note: Only for time dependent fields (usually electric)
   //                  is it neccessary to integrate the time.]
   yOut[7] = temp_eval_pt[7] = yIn[7]; // I suppose still do it this way??
   yOut[6] = temp_eval_pt[6] = yIn[6];
   // const G4int numberOfVariables= this->GetNumberOfVariables();
   // The number of variables to be integrated over
   // Does this actually vary??
   // Saving yInput because yInput and yOut can be aliases for same array

   for(i = 0; i < 3; i ++){
      pos[i] = yIn[i];
      velocity[i] = yIn[i+3] * iMass;
      //m_velocity_sqrd += velocity[i]*velocity[i];
   }

   G4double Delta_time = Step; // / G4ThreeVector( velocity[0], velocity[1], velocity[2] ).mag();
   G4double Delta_arc_len = Step * G4ThreeVector( velocity[0], velocity[1], velocity[2] ).mag();
   // redundant I know, but it's just for another test

   t = yIn[7];

	for(i = 0; i < 3; i ++){
		temp_eval_pt[i] = pos[i] + alpha1*Step*velocity[i];
	}
	for(i = 0; i < 3; i ++){
		temp_eval_pt[i+3] = velocity[i];
	}
	temp_eval_pt[6] = t + alpha1*Step;

	ComputeRightHandSide(temp_eval_pt, K1);

	for(i = 0; i < 3; i ++){
		temp_eval_pt[i] = pos[i] + alpha2*Step*velocity[i] + Step*Step*beta21*K1[i+3];
	}
	for(i = 0; i < 3; i ++){
		temp_eval_pt[i+3] = velocity[i] + Step*gamma21*K1[i+3];
	}
	temp_eval_pt[6] = t + alpha2*Step;

	ComputeRightHandSide(temp_eval_pt, K2);

	for(i = 0; i < 3; i ++){
		temp_eval_pt[i] = pos[i] + alpha3*Step*velocity[i] + Step*Step*
					(beta31*K1[i+3] + beta32*K2[i+3]);
	}
	for(i = 0; i < 3; i ++){
	temp_eval_pt[i+3] = velocity[i] + Step*(gamma31*K1[i+3] + gamma32*K2[i+3]);
	}
	temp_eval_pt[6] = t + alpha3*Step;

	ComputeRightHandSide(temp_eval_pt, K3);

	for(i = 0; i < 3; i ++)
	{
	// Accumulate increments with proper weights

		yOut[i] = pos[i] + Step*velocity[i] + Step*Step*
				(a1*K1[i+3] + a2*K2[i+3] + a3*K3[i+3]);
	}
	for(i = 0; i < 3; i ++){
		yOut[i+3] = velocity[i] + Step*( b1*K1[i+3] + b2*K2[i+3] + b3*K3[i+3] );
	}

	for (i = 3; i < 6; i ++){
				yOut[i] *= mass;
		}

	yOut[6] = yIn[6] + Delta_time;
	yOut[7] = yIn[7] + Delta_arc_len;
	return ;
}
