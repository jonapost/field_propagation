#ifndef TCASHKARPRKF45
#define TCASHKARPRKF45

#include "G4LineSection.hh"
#include "G4MagIntegratorStepper.hh"

template
<class T_Equation, int N>
class TCashKarpRKF45 : public G4MagIntegratorStepper
{
	public:
		TCashKarpRKF45(T_Equation *EqRhs,
				G4int noIntegrationVariables=6, 
				G4bool primary=true)
			: G4MagIntegratorStepper(
					dynamic_cast<G4EquationOfMotion*>(EqRhs), 
					noIntegrationVariables),
			  fLastStepLength(0.), fAuxStepper(0),
			  fEquation_Rhs(EqRhs)
		{
			const G4int numberOfVariables = noIntegrationVariables;
			assert(numberOfVariables==N);
			
			fLastInitialVector = new G4double[N] ;
			fLastFinalVector = new G4double[N] ;
			fLastDyDx = new G4double[N];

			fMidVector = new G4double[N];
			fMidError =  new G4double[N];

			if( primary )
			{ 
				fAuxStepper = new TCashKarpRKF45(EqRhs, numberOfVariables, !primary);
			}
		}

		virtual ~TCashKarpRKF45()
		{
			delete[] fLastInitialVector;
			delete[] fLastFinalVector;
			delete[] fLastDyDx;
			delete[] fMidVector;
			delete[] fMidError;

			delete fAuxStepper;
		}

		inline void TRightHandSide(G4double y[], G4double dydx[]) 
		{fEquation_Rhs->T_Equation::TRightHandSide(y, dydx);}

		void Stepper(const G4double __restrict__ yInput[],
				const G4double __restrict__ dydx[],
				G4double Step,
				G4double __restrict__ yOut[],
				G4double __restrict__ yErr[])
		{
			// const G4int nvar = 6 ;
			// const G4double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
			G4int i;

			const G4double  b21 = 0.2 ,
				  b31 = 3.0/40.0 , b32 = 9.0/40.0 ,
				  b41 = 0.3 , b42 = -0.9 , b43 = 1.2 ,

				  b51 = -11.0/54.0 , b52 = 2.5 , b53 = -70.0/27.0 ,
				  b54 = 35.0/27.0 ,

				  b61 = 1631.0/55296.0 , b62 =   175.0/512.0 ,
				  b63 =  575.0/13824.0 , b64 = 44275.0/110592.0 ,
				  b65 =  253.0/4096.0 ,

				  c1 = 37.0/378.0 , c3 = 250.0/621.0 , c4 = 125.0/594.0 ,
				  c6 = 512.0/1771.0 ,
				  dc5 = -277.0/14336.0 ;

			const G4double dc1 = c1 - 2825.0/27648.0 ,  
				  dc3 = c3 - 18575.0/48384.0 ,
				  dc4 = c4 - 13525.0/55296.0 , 
				  dc6 = c6 - 0.25 ;

			// Initialise time to t0, needed when it is not updated by the integration.
			//       [ Note: Only for time dependent fields (usually electric) 
			//                 is it neccessary to integrate the time.] 
			//yOut[7] = yTemp[7]   = yIn[7]; 

			//  Saving yInput because yInput and yOut can be aliases for same array
			for(i=0;i<N;i++) 
			{
				yIn[i]=yInput[i];
			}
			// TRightHandSide(yIn, dydx) ;              // 1st Step

			for(i=0;i<N;i++) 
			{
				yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
			}
			TRightHandSide(yTemp, ak2) ;              // 2nd Step

			for(i=0;i<N;i++)
			{
				yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
			}
			TRightHandSide(yTemp, ak3) ;              // 3rd Step

			for(i=0;i<N;i++)
			{
				yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
			}
			TRightHandSide(yTemp, ak4) ;              // 4th Step

			for(i=0;i<N;i++)
			{
				yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
						b54*ak4[i]) ;
			}
			TRightHandSide(yTemp, ak5) ;              // 5th Step

			for(i=0;i<N;i++)
			{
				yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
						b64*ak4[i] + b65*ak5[i]) ;
			}
			TRightHandSide(yTemp, ak6) ;              // 6th Step

			for(i=0;i<N;i++)
			{
				// Accumulate increments with proper weights

				yOut[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]) ;
			}
			for(i=0;i<N;i++)
			{
				// Estimate error as difference between 4th and
				// 5th order methods

				yErr[i] = Step*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] +
						dc5*ak5[i] + dc6*ak6[i]) ;
			}
			for(i=0;i<N;i++)
			{
				// Store Input and Final values, for possible use in calculating chord
				fLastInitialVector[i] = yIn[i] ;
				fLastFinalVector[i]   = yOut[i];
				fLastDyDx[i]          = dydx[i];
			}
			// NormaliseTangentVector( yOut ); // Not wanted

			fLastStepLength =Step;

			return ;
		}

		G4double  DistChord()   const
		{
			G4double distLine, distChord; 
			G4ThreeVector initialPoint, finalPoint, midPoint;

			// Store last initial and final points (they will be overwritten in self-Stepper call!)
			initialPoint = G4ThreeVector( fLastInitialVector[0], 
					fLastInitialVector[1], fLastInitialVector[2]); 
			finalPoint   = G4ThreeVector( fLastFinalVector[0],  
					fLastFinalVector[1],  fLastFinalVector[2]); 

			// Do half a step using StepNoErr

			fAuxStepper->TCashKarpRKF45::Stepper( 
					fLastInitialVector, 
					fLastDyDx, 
					0.5 * fLastStepLength, 
					fMidVector,   
					fMidError );

			midPoint = G4ThreeVector( fMidVector[0], fMidVector[1], fMidVector[2]);       

			// Use stored values of Initial and Endpoint + new Midpoint to evaluate
			//  distance of Chord


			if (initialPoint != finalPoint) 
			{
				distLine  = G4LineSection::Distline( midPoint, initialPoint, finalPoint );
				distChord = distLine;
			}
			else
			{
				distChord = (midPoint-initialPoint).mag();
			}
			return distChord;
		}

		G4int IntegratorOrder() const { return 4; }

	private:
		TCashKarpRKF45(const TCashKarpRKF45&);
		TCashKarpRKF45& operator=(const G4CashKarpRKF45&);
		//private copy constructor and assignment operator.

	private:
		G4double ak2[N];
		G4double ak3[N];
		G4double ak4[N];
		G4double ak5[N];
		G4double ak6[N];
		G4double ak7[N];
		G4double yTemp[N];
		G4double yIn[N];
		// scratch space

		G4double  fLastStepLength;
		G4double* fLastInitialVector;
		G4double* fLastFinalVector;
		G4double* fLastDyDx; 
		G4double* fMidVector;
		G4double* fMidError;
		// for DistChord calculations

		TCashKarpRKF45* fAuxStepper; 
		T_Equation* fEquation_Rhs;
};

#endif /*TCashKARP_RKF45 */
