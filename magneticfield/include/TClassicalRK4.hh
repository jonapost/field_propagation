#include "TMagErrorStepper.hh"
#include "G4ThreeVector.hh"
#include "G4MagIntegratorStepper.hh"
template
<class T_Equation, int N>
class TClassicalRK4 : public  TMagErrorStepper
                      <TClassicalRK4<T_Equation, N>, T_Equation, N>
{
    public:  // with description

        static const double 
            IntegratorCorrection = 1./((1<<4)-1);

        TClassicalRK4(T_Equation *EqRhs, 
                      G4int numberOfVariables = 8)
            : TMagErrorStepper
              <TClassicalRK4<T_Equation, N>, T_Equation, N>
              (EqRhs, numberOfVariables),
              fEquation_Rhs(EqRhs)
        {
            unsigned int noVariables= 
                std::max(numberOfVariables,8); // For Time .. 7+1 
        }

        virtual ~TClassicalRK4(){;}

        __attribute__((always_inline)) 
        void 
        RightHandSideInl(G4double y[], G4double dydx[]) 
        { fEquation_Rhs->T_Equation::RightHandSideInl(y, dydx); }


        // A stepper that does not know about errors.
        // It is used by the MagErrorStepper stepper.

        __attribute__((always_inline)) 
        void  
        DumbStepper( const G4double  yIn[],
                     const G4double  dydx[],
                     G4double  h,
                     G4double  yOut[])
        // Given values for the variables y[0,..,n-1] and their derivatives
        // dydx[0,...,n-1] known at x, use the classical 4th Runge-Kutta
        // method to advance the solution over an interval h and return the
        // incremented variables as yout[0,...,n-1], which not be a distinct
        // array from y. The user supplies the routine RightHandSide(x,y,dydx),
        // which returns derivatives dydx at x. The source is routine rk4 from
        // NRC p. 712-713 .
        {
            G4int i;
            G4double  hh = h*0.5 , h6 = h/6.0  ;

            // Initialise time to t0, needed when it is not updated by the integration.
            //        [ Note: Only for time dependent fields (usually electric) 
            //                  is it neccessary to integrate the time.] 
            yt[7]   = yIn[7]; 
            yOut[7] = yIn[7];

            for(i=0;i<N;i++)
            {
                yt[i] = yIn[i] + hh*dydx[i] ;             // 1st Step K1=h*dydx
            }
            this->RightHandSideInl(yt,dydxt) ;                   // 2nd Step K2=h*dydxt

            for(i=0;i<N;i++)
            { 
                yt[i] = yIn[i] + hh*dydxt[i] ;
            }
            this->RightHandSideInl(yt,dydxm) ;                   // 3rd Step K3=h*dydxm

            for(i=0;i<N;i++)
            {
                yt[i]   = yIn[i] + h*dydxm[i] ;
                dydxm[i] += dydxt[i] ;                    // now dydxm=(K2+K3)/h
            }
            this->RightHandSideInl(yt,dydxt) ;                   // 4th Step K4=h*dydxt

            for(i=0;i<N;i++)    // Final RK4 output
            {
                yOut[i] = yIn[i]+h6*(dydx[i]+dydxt[i]+2.0*dydxm[i]); //+K1/6+K4/6+(K2+K3)/3
            }
            if ( N == 12 )  { this->NormalisePolarizationVector ( yOut ); }

        }  // end of DumbStepper ....................................................



    public:  // without description

        __attribute__((always_inline)) 
        G4int 
        IntegratorOrder() const { return 4; }

    private:
        TClassicalRK4(const TClassicalRK4&);
        TClassicalRK4& operator=(const TClassicalRK4&);
        // Private copy constructor and assignment operator.

    private:

        // G4int fNumberOfVariables ; // is set default to 6 in constructor
        G4double dydxm[N<8?8:N]; 
        G4double dydxt[N<8?8:N]; 
        G4double yt[N<8?8:N]; 
        // scratch space - not state 

        T_Equation *fEquation_Rhs;
};
