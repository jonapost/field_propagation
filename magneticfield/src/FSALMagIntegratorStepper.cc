

#include "FSALMagIntegratorStepper.hh"

// Constructor for stepper abstract base class. 
// 

FSALMagIntegratorStepper::FSALMagIntegratorStepper(G4EquationOfMotion* Equation,
					       G4int       num_integration_vars,
					       G4int       num_state_vars)
  : fEquation_Rhs(Equation),
    fNoIntegrationVariables(num_integration_vars),
    fNoStateVariables(num_state_vars)
    // fNumberOfVariables( std::max(num_var,fNoStateVariables) )
{
}

FSALMagIntegratorStepper::~FSALMagIntegratorStepper()
{
}

void FSALMagIntegratorStepper::ComputeRightHandSide( const G4double y[], G4double dydx[] ) 
{
  this->RightHandSide( y, dydx );
//	fEquation_Rhs->RightHandSide(y, dydx);
//        increasefNORHSCalls();


}

void FSALMagIntegratorStepper::increasefNORHSCalls(){
    //    std::cout<<"Yeah, I was called!";
    fNoRHSCalls++;
}


void FSALMagIntegratorStepper::RightHandSide( const  double y[], double dydx[] )
{
    fEquation_Rhs-> RightHandSide(y, dydx);
    increasefNORHSCalls();
}


