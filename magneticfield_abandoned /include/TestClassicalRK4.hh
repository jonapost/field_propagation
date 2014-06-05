//ignore this file. Not developed yet.

#include "G4MagErrorStepper.hh"

template
<class State, class Equation>
class TestClassicalRK4 : public G4MagErrorStepper 
{
  public:  
    TestClassicalRK4(Equation *EquationMotion, G4int numberOfVariables = 6) ;
    ~TestClassicalRK4() ;
    void DumbStepper(State s, G4double h) ;

  public:  
    G4int IntegratorOrder() const { return 4; }

  private:
    TestClassicalRK4(const G4ClassicalRK4&);
    TestClassicalRK4& operator=(const G4ClassicalRK4&);
};
