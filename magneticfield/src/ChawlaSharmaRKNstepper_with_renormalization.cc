#include "ChawlaSharmaRKNstepper_with_renormalization.hh"
using namespace CLHEP;


ChawlaSharmaRKNstepper_with_renormalization::~ChawlaSharmaRKNstepper_with_renormalization()
{
   delete[] yMiddle;
   delete[] dydxMid;
   delete[] yInitial;
   delete[] yOneStep;
}

void ChawlaSharmaRKNstepper_with_renormalization::Stepper(  const G4double yInput[],
		            const G4double dydx[],
		                  G4double hstep,
		                  G4double yOutput[],
		                  G4double yError [] ){
		                  
   const G4int nvar = this->GetNumberOfVariables() ;
   const G4int maxvar= GetNumberOfStateVariables();
   
   G4int i;
   
   //for (i = 0; i < 6; i ++)
   //         std::cerr << yInput[i] << ", ";
   //     std::cerr << std::endl;
   
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
   
   G4double m_mom = std::sqrt(yInitial[3]*yInitial[3]+yInitial[4]*yInitial[4]+yInitial[5]*yInitial[5]);
   for(i = 3; i < 6; i++) yInitial[i] /= m_mom;
   
   DumbStepper  (yInitial, halfStep, yMiddle);
   
   //std::cout << std::sqrt( yMiddle[3]*yMiddle[3] + yMiddle[4]*yMiddle[4] + yMiddle[5]*yMiddle[5] ) << std::endl;
   
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

void ChawlaSharmaRKNstepper_with_renormalization::DumbStepper( 
										const G4double yIn[],
										// const G4double dydx[],
										G4double Step,
										G4double yOut[] ) {
	// dydx is not even used by DumbStepper. See Comment (**) above.
	
	G4int i;
	G4double K1[3], K2[3], K3[3];
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
      //std::cerr << "yIn[i] -> mom[i] : " <<  yIn[i+3] << ", " << mom[i] << std::endl;
      //m_mom_sqrd += mom[i]*mom[i];
   }
   //for (i = 0; i < 3; i ++)
   //   std::cerr << "(Step, mom[i]) = " << Step << ", " << mom[i] << ";  ";
   //std::cerr << std::endl;
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
	
	mEvaluateRhs(temp_eval_pt, K1);
   
   //std::cerr << "K1 = { " << K1[0] << ", " << K1[1] << ", " << K1[2] << " }\n";
	
	for(i = 0; i < 3; i ++){
		temp_eval_pt[i] = pos[i] + alpha2*Step*mom[i] + Step*Step*beta21*K1[i];
	}
	for(i = 0; i < 3; i ++){
		temp_eval_pt[i+3] = mom[i] + Step*gamma21*K1[i];
	}
	temp_eval_pt[7] = t + alpha2*Step;
   
	mEvaluateRhs(temp_eval_pt, K2);
	
	//std::cerr << "K2 = { " << K2[0] << ", " << K2[1] << ", " << K2[2] << " }\n";

	for(i = 0; i < 3; i ++){
		temp_eval_pt[i] = pos[i] + alpha3*Step*mom[i] + Step*Step*
					(beta31*K1[i] + beta32*K2[i]);
	}
	for(i = 0; i < 3; i ++){
	temp_eval_pt[i+3] = mom[i] + Step*(gamma31*K1[i] + gamma32*K2[i]);
	}
	temp_eval_pt[7] = t + alpha3*Step;
   
	mEvaluateRhs(temp_eval_pt, K3);
   
   //std::cerr << "K3 = { " << K3[0] << ", " << K3[1] << ", " << K3[2] << " }\n";
   
	for(i = 0; i < 3; i ++)
	{
	// Accumulate increments with proper weights

		yOut[i] = pos[i] + Step*mom[i] + Step*Step*
				(a1*K1[i] + a2*K2[i] + a3*K3[i]);
	   //std::cerr << "pos = " <<  pos[i] << ", Step*mom[i] = " << Step*mom[i] << ", last part = " << Step*Step*(a1*K1[i] + a2*K2[i] + a3*K3[i]) << std::endl;
	}
	for(i = 0; i < 3; i ++){
		yOut[i+3] = mom[i] + Step*( b1*K1[i] + b2*K2[i] + b3*K3[i] );
	}
	
	//G4double normF = 1. / std::sqrt(yOut[3]*yOut[3] + yOut[4]*yOut[4] + yOut[5]*yOut[5]);
	//yOut[3] *= normF; yOut[4] *= normF; yOut[5] *= normF;
	
	// Nomalization above is commented out because we need to not 
	// normalize in between the half steps
	
	return ;
}


void
ChawlaSharmaRKNstepper_with_renormalization::mEvaluateRhs( const G4double y[],
				           G4double dmom[] ) const
{ 
   G4double Point[4] = { y[0], y[1], y[2], y[7] };
   G4double B[3];
   m_fEq->GetFieldObj()->GetFieldValue( Point, B );
   
   // Temp change:
   G4double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
   G4double inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );
   // Temp change.
   
   G4double cof = m_fEq->FCof() ; //* inv_momentum_magnitude;
   
   // Temp insertion:
   cof *= inv_momentum_magnitude;
   
   // try this:
   //G4double cof = inv_momentum_magnitude;
   //G4double cof = 1.;
   /*
   dmom[0] = (y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dmom[1] = (y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dmom[2] = (y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)
   */
   
   dmom[0] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dmom[1] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dmom[2] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)
   
   return ;
}

G4double ChawlaSharmaRKNstepper_with_renormalization::DistChord() const 
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
