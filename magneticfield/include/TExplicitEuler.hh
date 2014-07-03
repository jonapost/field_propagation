#ifndef TExplicitEuler_HH
#define TExplicitEuler_HH

#include "TMagErrorStepper.hh"
#include "G4ThreeVector.hh"

template
<class Equation, int N>
class TExplicitEuler : public TMagErrorStepper            
                    <TExplicitEuler<T_Equation, N>, T_Equation, N>
{
    public:  // with description
    typedef Equation T_Equation;

    static const double 
    IntegratorCorrection = 1.;
    
    TExplicitEuler(T_Equation* EqRhs, 
                   G4int numberOfVariables = 6)
        :  TMagErrorStepper
           <TExplicitEuler<T_Equation, N>, T_Equation, N>
           (EqRhs, numberOfVariables),
           fEquation_Rhs(EqRhs){;}

    ~TExplicitEuler() {;}

    __attribute__((always_inline)) 
    void
    DumbStepper( const G4double  yIn[],
                 const G4double  dydx[],
                 G4double  h,
                 G4double  yOut[]        )
    {
        const G4int numberOfVariables = 
            this->GetNumberOfVariables();
        assert(numberOfVariables == N);
        // Initialise time to t0, needed when it is not updated by the integration.
        // yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO

        G4int i;

        for(i=0;i< N;i++)
        {
            yOut[i] = yIn[i] + h*dydx[i] ;             // 1st and only Step 
        }

        return ;
    }  

    public:  // without description

    __attribute__((always_inline)) 
    G4int  
    IntegratorOrder() const { return 1; }

    private:

    T_Equation* fEquation_Rhs;

};

#endif /* TExplicitEuler_HH */
