#include "ChawlaSharmaRKNstepper.hh"

ChawlaSharmaRKNstepper::~ChawlaSharmaRKNstepper()
{
   delete[] yMiddle;
   delete[] dydxMid;
   delete[] yInitial;
   delete[] yOneStep;
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
   G4double  correction = 1. / ( (1 << IntegratorOrder()) -1 );
   
   //  Saving yInput because yInput and yOutput can be aliases for same array

   for(i=0;i<nvar;i++) yInitial[i]=yInput[i];
   yInitial[7]= yInput[7];    // Copy the time in case ... even if not really needed
   yMiddle[7] = yInput[7];  // Copy the time from initial value 
   yOneStep[7] = yInput[7]; // As it contributes to final value of yOutput ?
   // yOutput[7] = yInput[7];  // -> dumb stepper does it too for RK4
   for(i=nvar;i<maxvar;i++) yOutput[i]=yInput[i];
   // yError[7] = 0.0;         
   
   G4double halfStep = hstep * 0.5; 
   
   // Comment (*):
   // Do two half steps 
   
   G4double m_mom = std::sqrt(yInitial[3]*yInitial[3]+yInitial[4]*yInitial[4]+yInitial[5]*yInitial[5]);
   for(i = 3; i < 6; i++) yInitial[i] /= m_mom;
   
   DumbStepper  (yInitial,  dydx,   halfStep, yMiddle);
   
   //RightHandSide(yMiddle, dydxMid);
   
   // Comment (**):
   // Function evaluations f(x,y,y') are not at (x_k, y_k, y'_k)
   // but rather at (x_k, y_k + alpha1*halfStep*y'_k, y'_k).
   // There is no point in calculating dy/dx at yMiddle.
   // See p. 376 of Chawla/Sharma.
   
   DumbStepper  (yMiddle, dydxMid, halfStep, yOutput); 

   // Store midpoint, chord calculation

   fMidPoint = G4ThreeVector( yMiddle[0],  yMiddle[1],  yMiddle[2]); 

   // Do a full Step
   DumbStepper(yInitial, dydx, hstep, yOneStep);
   for(i=0;i<nvar;i++) {
      yError [i] = yOutput[i] - yOneStep[i] ;
      yOutput[i] += yError[i]*correction ;  // Provides accuracy increased
                                            // by 1 order via the 
                                            // Richardson Extrapolation.
                                            // Question: Is the order reported as 3 to MagIntegratorDriver
                                            // ( at MagIntegratorDriver.cc line 90 )
                                            // when it should be reported as 3+1 = 4 ?
   }
   
   // Now renomalize momentum
   G4double momentum_sqrt = std::sqrt( yOutput[3]*yOutput[3] + yOutput[4]*yOutput[4] + yOutput[5]*yOutput[5] );
   
   for(i = 3; i < 6; i ++){
      yOutput[i] /= momentum_sqrt;
   }
   
   fInitialPoint = G4ThreeVector( yInitial[0], yInitial[1], yInitial[2]); 
   fFinalPoint   = G4ThreeVector( yOutput[0],  yOutput[1],  yOutput[2]); 

   return ;
}

void ChawlaSharmaRKNstepper::DumbStepper( 
										const G4double yIn[],
										const G4double dydx[],
										G4double Step,
										G4double yOut[] ) {
	// dydx is not even used by DumbStepper. See Comment (**) above.
	
	G4int i;
	G4double yTemp[8], K1[6], K2[6], K3[6];
	const G4double   a1 = 1./4. , a2 = 0., a3 = 1./4.,
 									b1 = 1./4., b2 = 0, b3 = 3./4.,
 									alpha1 = 0., alpha2 = 2./3., alpha3 = 2./3.,
 									beta21 = -1./9., beta31 = 2./9., beta32 = 0.,
 									gamma21 = 2./3., gamma31 = 1./3., gamma32 = 1./3.;
   
   // uncomment either of the following two lines (on lines 107 and 108)
   G4double c = 1. / m_fEq->FCof();
   //G4double c = 1.;
   // The first option undoes the scaling by m_fEq->FCof(), which is done in RightHandSide()
   // ( RightHandSide() in our case of our test program testPropagateMagField eventually
   // calls G4Mag_UsualEqRhs::EvaluateRhsGivenB(...) which does this scaling by 
   //  m_fEq->FCof() )
   // I have no idea why the first option works to produce better results?
   
   // Initialise time to t0, needed when it is not updated by the integration.
   //        [ Note: Only for time dependent fields (usually electric) 
   //                  is it neccessary to integrate the time.] 
   yOut[7] = yTemp[7] = yIn[7]; // I suppose still do it this way??
   yOut[6] = yTemp[6] = yIn[6];
   // const G4int numberOfVariables= this->GetNumberOfVariables(); 
   // The number of variables to be integrated over
   // Does this actually vary??
   //  Saving yInput because yInput and yOut can be aliases for same array
   
	for(i = 0; i < 3; i ++){
		yTemp[i] = yIn[i] + alpha1*Step*yIn[i+3];
	}
	for(i = 3; i < 6; i ++){
		yTemp[i] = yIn[i];
	}
	yTemp[7] = yIn[7] + alpha1*Step;
	
	G4double L1 = std::sqrt( yTemp[3]*yTemp[3] + yTemp[4]*yTemp[4] + yTemp[5]*yTemp[5] );
	// L1 is used to un-normalize momentum = (K1[3], K1[4], K1[5]) because 
	// RightHandSide(yTemp, K1) will normalize momentum (and multiply by a 
	// coefficient given by m_fEq->FCof() ).
	// We do this because we don't want to renormalize in between the 2 half steps in Comment (*)
	// Note: this is definitely not the best way to do this.
	// Should probably write our own RightHandSide() for this class.
	
	RightHandSide(yTemp, K1);

	for(i = 0; i < 3; i ++){
		yTemp[i] = yIn[i] + alpha2*Step*yIn[i+3] + Step*Step*beta21*c*L1*K1[3 + i];
	// Is "3 + i" because for this method we want y" which is yIn[3:5],
	// because yIn[3..5] are the components of the derivative of y'_k
	// (ComputeRightHandSide() treats the ODE system as a first order system
	// and returns a 6 vector. But we only care about the momentum components.)
	}
	for(i = 3; i < 6; i ++){
		yTemp[i] = yIn[i] + Step*gamma21*c*L1*K1[i];
	}
	yTemp[7] = yIn[7] + alpha2*Step;
   
   G4double L2 = std::sqrt( yTemp[3]*yTemp[3] + yTemp[4]*yTemp[4] + yTemp[5]*yTemp[5] );
   // Similar comment as for L1
   
	RightHandSide(yTemp, K2);

	for(i = 0; i < 3; i ++){
		yTemp[i] = yIn[i] + alpha3*Step*yIn[i+3] + Step*Step*
					(beta31*c*L1*K1[3 + i] + beta32*c*L2*K2[3 + i]);
	}
	for(i = 3; i < 6; i ++){
	yTemp[i] = yIn[i] + Step*(gamma31*c*L1*K1[i] + gamma32*c*L2*K2[i]);
	}
	yTemp[7] = yIn[7] + alpha3*Step;
   
   G4double L3 = std::sqrt( yTemp[3]*yTemp[3] + yTemp[4]*yTemp[4] + yTemp[5]*yTemp[5] );
   // Similar comment as for L1
   
	RightHandSide(yTemp, K3);

	for(i = 0; i < 3; i ++)
	{
	// Accumulate increments with proper weights

		yOut[i] = yIn[i] + Step*yIn[i+3] + Step*Step*
				(a1*c*L1*K1[3 + i] + a2*c*L2*K2[3 + i] + a3*c*L3*K3[3 + i]);
	}
	for(i = 3; i < 6; i ++){
		yOut[i] = yIn[i] + Step*(b1*c*L1*K1[i] + b2*c*L2*K2[i] + b3*c*L3*K3[i]);
	}
	
	//G4double normF = 1. / std::sqrt(yOut[3]*yOut[3] + yOut[4]*yOut[4] + yOut[5]*yOut[5]);
	//yOut[3] *= normF; yOut[4] *= normF; yOut[5] *= normF;
	
	// Nomalization above is commented out because we need to not 
	// normalize in between the half steps
	
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


