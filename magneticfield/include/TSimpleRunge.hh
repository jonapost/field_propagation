#ifndef TSimpleRunge_HH
#define TSimpleRunge_HH

#include "TMagErrorStepper.hh"
#include "G4ThreeVector.hh"

template
<class T_Equation, int N>
class TSimpleRunge : public TMagErrorStepper            
                     <TSimpleRunge<T_Equation, N>, T_Equation, N>
{

    public:  // with description

        static const double 
        IntegratorCorrection = 1./((1<<2)-1);

        TSimpleRunge(T_Equation* EqRhs, 
                     G4int numberOfVariables = 6)
            :  TMagErrorStepper
              <TSimpleRunge<T_Equation, N>, T_Equation, N>
              (EqRhs, numberOfVariables),
              fEquation_Rhs(EqRhs),
              fNumberOfVariables(numberOfVariables)
        {
            //default GetNumberOfStateVariables() == 12 
            assert (this->GetNumberOfStateVariables() <= 12);
        }


        ~TSimpleRunge(){;}


        __attribute__((always_inline)) 
        void 
        RightHandSideInl(G4double y[], G4double dydx[]) 
        { fEquation_Rhs->T_Equation::RightHandSideInl(y, dydx); }


        __attribute__((always_inline)) 
        void
        DumbStepper( const G4double  yIn[],
                     const G4double  dydx[],
                     G4double  h,
                     G4double  yOut[])
        {
            // Initialise time to t0, needed when it is not updated by the integration.
            yTemp[7] = yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO

            G4int i;

            for( i = 0; i < N; i++ ) 
            {
                yTemp[i] = yIn[i] + 0.5 * h*dydx[i] ;
            }

            this->RightHandSideInl(yTemp,dydxTemp);

            for( i = 0; i < N; i++ ) 
            {
                yOut[i] = yIn[i] + h * ( dydxTemp[i] );
            }
        } 

    public:  // without description

        __attribute__((always_inline)) 
        G4int  
        IntegratorOrder() const { return 2; }

    private:

        G4int fNumberOfVariables ;
        G4double dydxTemp[N>12?N:12];
        G4double yTemp[N>12?N:12]   ;

        T_Equation *fEquation_Rhs;
        // scratch space    
};

#endif /* TSimpleRunge_HH */
