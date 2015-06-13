// Nystrom stepper implemenations and testing by Jason Suagee
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
   G4double  correction = 1. / ( (1 << (IntegratorOrder()-1)) -1 ); // IntegratorOrder() - 1 because we 
                                                                    //falsely report to higher classes integrator order 4
   
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
   
   G4double m_imom = 1. / std::sqrt(yInitial[3]*yInitial[3]+yInitial[4]*yInitial[4]+yInitial[5]*yInitial[5]);
   for(i = 3; i < 6; i++) yInitial[i] *= m_imom;
   
   DumbStepper  (yInitial, halfStep, yMiddle);
   
   /*G4double yMiddle_momentum_sqrt = std::sqrt( yMiddle[3]*yMiddle[3] + yMiddle[4]*yMiddle[4] + yMiddle[5]*yMiddle[5] );
   for(i = 3; i < 6; i ++){
      yMiddle[i] /= yMiddle_momentum_sqrt;
   }*/
   
   //RightHandSide(yMiddle, dydxMid);
   
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
										// const G4double dydx[],
										G4double Step,
										G4double yOut[] ) {
	// dydx is not even used by DumbStepper. See Comment (**) above.
	
	G4int i;
	G4double K1[6], K2[6], K3[6];
	G4double pos[3], mom[3], t;
	G4double temp_eval_pt[8];
	//G4double m_mom_sqrd = 0., m_mom;
	const G4double   a1 = 1./4. , a2 = 0., a3 = 1./4.,
 									b1 = 1./4., b2 = 0., b3 = 3./4.,
 									alpha1 = 0., alpha2 = 2./3., alpha3 = 2./3.,
 									beta21 = -1./9., beta31 = 2./9., beta32 = 0.,
 									gamma21 = 2./3., gamma31 = 1./3., gamma32 = 1./3.;
 	
   /*
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
   //  Saving yInput because yInput and yOut can be aliases for same array
   
   for(i = 0; i < 3; i ++){
      pos[i] = yIn[i];
      mom[i] = yIn[i+3];
      //m_mom_sqrd += mom[i]*mom[i];
   }
   /*m_mom = std::sqrt(m_mom_sqrd);
   for(i = 0; i < 3; i ++){
      mom[i] /= m_mom;
   }*/
   
   t = yIn[7];
      
	for(i = 0; i < 3; i ++){
		temp_eval_pt[i] = pos[i] + alpha1*Step*mom[i];
	}
	for(i = 0; i < 3; i ++){
		temp_eval_pt[i+3] = mom[i];
	}
	temp_eval_pt[7] = t + alpha1*Step;
	
	ComputeRightHandSide(temp_eval_pt, K1);
	
	for(i = 0; i < 3; i ++){
		temp_eval_pt[i] = pos[i] + alpha2*Step*mom[i] + Step*Step*beta21*K1[i+3];
	}
	for(i = 0; i < 3; i ++){
		temp_eval_pt[i+3] = mom[i] + Step*gamma21*K1[i+3];
	}
	temp_eval_pt[7] = t + alpha2*Step;
   
	ComputeRightHandSide(temp_eval_pt, K2);

	for(i = 0; i < 3; i ++){
		temp_eval_pt[i] = pos[i] + alpha3*Step*mom[i] + Step*Step*
					(beta31*K1[i+3] + beta32*K2[i+3]);
	}
	for(i = 0; i < 3; i ++){
	temp_eval_pt[i+3] = mom[i] + Step*(gamma31*K1[i+3] + gamma32*K2[i+3]);
	}
	temp_eval_pt[7] = t + alpha3*Step;
   
	ComputeRightHandSide(temp_eval_pt, K3);
   
	for(i = 0; i < 3; i ++)
	{
	// Accumulate increments with proper weights

		yOut[i] = pos[i] + Step*mom[i] + Step*Step*
				(a1*K1[i+3] + a2*K2[i+3] + a3*K3[i+3]);
	}
	for(i = 0; i < 3; i ++){
		yOut[i+3] = mom[i] + Step*( b1*K1[i+3] + b2*K2[i+3] + b3*K3[i+3] );
	}
	
	//G4double normF = 1. / std::sqrt(yOut[3]*yOut[3] + yOut[4]*yOut[4] + yOut[5]*yOut[5]);
	//yOut[3] *= normF; yOut[4] *= normF; yOut[5] *= normF;
	
	// Nomalization above is commented out because we need to not 
	// normalize in between the half steps
	
	return ;
}

/*
void
ChawlaSharmaRKNstepper::mEvaluateRhs( const G4double y[],
				           G4double dmom[] ) const
{ 
   G4double Point[4] = { y[0], y[1], y[2], y[7] };
   G4double B[3];
   m_fEq->GetFieldObj()->GetFieldValue( Point, B );
   
   G4double cof = m_fEq->FCof() ; // inv_momentum_magnitude;
   
   dmom[0] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dmom[1] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dmom[2] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)
   
   return ;
}
*/

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
