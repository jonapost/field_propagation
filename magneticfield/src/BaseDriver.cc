#include "BaseDriver.hh"

BaseDriver::BaseDriver(G4double hminimum,
                       G4EquationOfMotion* equation,
                       G4int numberOfComponents,
                       G4int  statisticsVerbosity):
    fMinimumStep(hminimum),
    fStepper(nullptr),
    fequation(equation),
    fnvar(numberOfComponents),
    fverb(statisticsVerbosity)
{
    system = [this](const state_type& y, state_type& dydx, double /*t*/){
        fequation->RightHandSide(y.data(),dydx.data());
    };
    System = [this](const state_type& y, state_type& dydx){
        fequation->RightHandSide(y.data(),dydx.data());
    };
}

BaseDriver::BaseDriver(G4double hminimum,
                       G4MagIntegratorStepper *pItsStepper,
                       G4int numberOfComponents,
                       G4int statisticsVerbosity):
    fMinimumStep(hminimum),
    fStepper(pItsStepper),
    fequation(pItsStepper->GetEquationOfMotion()),
    fnvar(numberOfComponents),
    fverb(statisticsVerbosity)
{

}

BaseDriver::~BaseDriver(){

}


G4double BaseDriver::ComputeNewStepSize(double /*errMaxNorm*/,
                                        double hstepCurrent ){
    return hstepCurrent;
}

void BaseDriver::SetVerboseLevel(G4int verb){
    fverb = verb;
}

