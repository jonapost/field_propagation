
#ifndef TSIMPLEHEUM_HH
#define TSIMPLEHEUM_HH

#include "TMagErrorStepper.hh"
#include "G4ThreeVector.hh"

template
<class T_Equation, int N>
class TSimpleHeum : public  TMagErrorStepper
                    <TSimpleHeum<T_Equation, N>, T_Equation, N>
{

    public:  // with description
        static const double 
            IntegratorCorrection = 0.14285714285714285;

        TSimpleHeum(T_Equation *EqRhs, G4int num_variables=6)
            :  TMagErrorStepper
               <TClassicalRK4<T_Equation, N>, T_Equation, N>
               (EqRhs, numberOfVariables),
               fEquation_Rhs(EqRhs),
               fNumberOfVariables(num_variables)
        {
             assert( fNumberOfVariables == N );
        }

        virtual ~TSimpleHeum() {;}
        // Constructor and destructor.


        __attribute__((always_inline)) 
            void TRightHandSide(G4double y[], G4double dydx[]) 
            { fEquation_Rhs->T_Equation::TRightHandSide(y, dydx); }

        void DumbStepper( const G4double  yIn[],
                const G4double  dydx[],
                G4double  h,
                G4double  yOut[])
        {
            G4int i;
            for( i = 0; i < N; i++ ) 
            {
                yTemp[i] = yIn[i] + (1.0/3.0) * h *  dydx[i] ;
            }

            TRightHandSide(yTemp,dydxTemp);

            for( i = 0; i < N; i++ ) 
            {
                yTemp2[i] = yIn[i] + (2.0/3.0) * h * dydxTemp[i] ;
            }

            TRightHandSide(yTemp2,dydxTemp2);

            for( i = 0; i < N; i++ ) 
            {
                yOut[i] = yIn[i] + h * (0.25 * dydx[i] + 0.75 * dydxTemp2[i]);
            }

            if ( fNumberOfVariables == 12 ) {  this->NormalisePolarizationVector( yOut ); }
        }  


    public:  // without description

        __attribute__((always_inline)) 
            G4int IntegratorOrder() const { return 3; }

    private:

        G4int fNumberOfVariables;

        G4double dydxTemp[N]  ;
        G4double dydxTemp2[N] ;
        G4double yTemp[N]     ;
        G4double yTemp2[N]    ;

        T_Equation *fEquation_Rhs;
        // scratch space    
};

#endif /* TSIMPLEHEUM_HH */
