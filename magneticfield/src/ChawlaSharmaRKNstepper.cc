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


#include "ChawlaSharmaRKNstepper.hh"
using namespace CLHEP;

//#include <iomanip>
//using namespace std;

ChawlaSharmaRKNstepper::~ChawlaSharmaRKNstepper()
{
   delete[] yMiddle;
   delete[] dydxMid;
   delete[] yInitial;
   delete[] yOneStep;
   delete[] K1; delete[] K2; delete[] K3;
   delete[] pos; delete[] mom;
   delete[] temp_eval_pt;

}

void ChawlaSharmaRKNstepper::Stepper(  const G4double yInput[],
		            const G4double dydx[],
		                  G4double hstep,
		                  G4double yOutput[],
		                  G4double yError [] ){
		                  

   const G4int nvar = this->GetNumberOfVariables() ;
   const G4int maxvar= GetNumberOfStateVariables();
   
   G4int i;
   // correction for Richardson Extrapolation.
   G4double  correction = 1. / ( (1 << (IntegratorOrder()-1)) -1 ); // IntegratorOrder() - 1 because we 
                                                                    // report to higher classes integrator order 4
                                                                    // but this is only true with Richardson extrapolation
   
   //  Saving yInput because yInput and yOutput can be aliases for same array

   for(i=0;i<nvar;i++) yInitial[i]=yInput[i];

   yInitial[7]= yInput[7];    // Copy the time in case ... even if not really needed
   yMiddle[7] = yInput[7];  // Copy the time from initial value 
   yOneStep[7] = yInput[7]; // As it contributes to final value of yOutput ?
   // yOutput[7] = yInput[7];  // -> dumb stepper does it too for RK4
   for(i=nvar;i<maxvar;i++) yOutput[i]=yInput[i];
   // yError[7] = 0.0;         
   
   G4double halfStep = hstep * 0.5; 
   
   // Do two half steps 
   
   DumbStepper  (yInitial, halfStep, yMiddle);
   
   // Might put in for one of the set of coefficients:
   // RightHandSide(yMiddle, dydxMid);
   
   // Comment (**):
   // Function evaluations f(x,y,y') are not at (x_k, y_k, y'_k)
   // but rather at (x_k, y_k + alpha1*halfStep*y'_k, y'_k).
   // There is no point in calculating dy/dx at yMiddle.
   // See p. 376 of Chawla/Sharma.
   
   DumbStepper  (yMiddle, halfStep, yOutput); 

   // Store midpoint, chord calculation
   
   fMidPoint = G4ThreeVector( yMiddle[0],  yMiddle[1],  yMiddle[2]); 

   // Do a full Step
   DumbStepper(yInitial, hstep, yOneStep);
   for(i=0;i<nvar;i++) {
      yError [i] = (yOutput[i] - yOneStep[i]) * correction ;
      yOutput[i] += yError[i] ;  // Provides accuracy increased
                                            // by 1 order via the 
                                            // Richardson Extrapolation.
                                            // Question: Is the order reported as 3 to MagIntegratorDriver
                                            // ( at MagIntegratorDriver.cc line 90 )
                                            // when it should be reported as 3+1 = 4 ?
   }

   fInitialPoint = G4ThreeVector( yInitial[0], yInitial[1], yInitial[2]); 
   fFinalPoint   = G4ThreeVector( yOutput[0],  yOutput[1],  yOutput[2]); 

   return ;
}

void ChawlaSharmaRKNstepper::DumbStepper( 
										const G4double yIn[],
										// const G4double dydx[],
										G4double Step,
										G4double yOut[] ) {
	// dydx is not even used by DumbStepper. See Comment (**) above.
	
	G4int i;
	G4double t;
	const G4double   a1 = 1./4. , a2 = 0., a3 = 1./4.,
 									b1 = 1./4., b2 = 0., b3 = 3./4.,
 									alpha1 = 0., alpha2 = 2./3., alpha3 = 2./3.,
 									beta21 = -1./9., beta31 = 2./9., beta32 = 0.,
 									gamma21 = 2./3., gamma31 = 1./3., gamma32 = 1./3.;
 	
   /* Another set of coefficients to try (thi set of coefficients could make use of FSAL):
    *
   const G4double   a1 = 0. , a2 = 1./2., a3 = 0.,
 									b1 = 0., b2 = 3./4., b3 = 1./4.,
 									alpha1 = 1./3., alpha2 = 1./3., alpha3 = 1.,
 									beta21 = 0., beta31 = 1., beta32 = -1./3.,
 									gamma21 = 1./3., gamma31 = -1., gamma32 = 2.;
   
   */
   
   // Initialise time to t0, needed when it is not updated by the integration.
   //[ Note: Only for time dependent fields (usually electric)
   // is it neccessary to integrate the time.]
   yOut[7] = temp_eval_pt[7] = yIn[7]; // I suppose still do it this way??
   yOut[6] = temp_eval_pt[6] = yIn[6];
   // Saving yInput because yInput and yOut can be aliases for same array
   
   G4int index_last_mom_var = nposvars + nmomvars;

   for(i = 0; i < nposvars; i ++)
      pos[i] = yIn[i];
   for (i = nposvars; i < index_last_mom_var; i ++)
      mom[i] = yIn[i];

   t = yIn[7];
      
	for(i = 0; i < nposvars; i ++){
		temp_eval_pt[i] = pos[i] + alpha1*Step*mom[i];
	}
	for(i = nposvars; i < index_last_mom_var; i ++){
		temp_eval_pt[i] = mom[i - 3];
	}
	temp_eval_pt[7] = t + alpha1*Step;
	
	ComputeRightHandSide(temp_eval_pt, K1);
	
	for(i = 0; i < nposvars; i ++){
		temp_eval_pt[i] = pos[i] + alpha2*Step*mom[i] + Step*Step*beta21*K1[i+3];
	}
	for(i = nposvars; i < index_last_mom_var; i ++){
		temp_eval_pt[i] = mom[i - 3] + Step*gamma21*K1[i];
	}
	temp_eval_pt[7] = t + alpha2*Step;
   
	ComputeRightHandSide(temp_eval_pt, K2);

	for(i = 0; i < nposvars; i ++){
		temp_eval_pt[i] = pos[i] + alpha3*Step*mom[i] + Step*Step*
					(beta31*K1[i+3] + beta32*K2[i+3]);
	}
	for(i = nposvars; i < index_last_mom_var; i ++){
	temp_eval_pt[i] = mom[i - 3] + Step*(gamma31*K1[i] + gamma32*K2[i]);
	}
	temp_eval_pt[7] = t + alpha3*Step;
   
	ComputeRightHandSide(temp_eval_pt, K3);
   
	for(i = 0; i < nposvars; i ++)
	{
	// Accumulate increments with proper weights

		yOut[i] = pos[i] + Step*mom[i] + Step*Step*
				(a1*K1[i+3] + a2*K2[i+3] + a3*K3[i+3]);
	}
	for(i = nposvars; i < index_last_mom_var; i ++){
		yOut[i] = mom[i - 3] + Step*( b1*K1[i] + b2*K2[i] + b3*K3[i] );
	}
	
	return ;
}

G4double ChawlaSharmaRKNstepper::DistChord() const 
{
  // Estimate the maximum distance from the curve to the chord
  //
  //  We estimate this using the distance of the midpoint to 
  //  chord (the line between 
  // 
  //  Method below is good only for angle deviations < 2 pi, 
  //   This restriction should not a problem for the Runge cutta methods, 
  //   which generally cannot integrate accurately for large angle deviations.
  G4double distLine, distChord; 

  if (fInitialPoint != fFinalPoint) {
     distLine= G4LineSection::Distline( fMidPoint, fInitialPoint, fFinalPoint );
     // This is a class method that gives distance of Mid 
     //  from the Chord between the Initial and Final points.

     distChord = distLine;
  }else{
     distChord = (fMidPoint-fInitialPoint).mag();
  }
  return distChord;
}
