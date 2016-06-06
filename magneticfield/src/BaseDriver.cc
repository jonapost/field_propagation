#include "BaseDriver.hh"

BaseDriver::BaseDriver(G4EquationOfMotion* equation,
                       G4int integratedComponents,
                       G4int  verb):
    fequation(equation),fnvar(integratedComponents),fverb(verb)
{
    system = [this](const state_type& y, state_type& dydx, double /*t*/){
        fequation->RightHandSide(y.data(),dydx.data());
    };
}

BaseDriver::~BaseDriver(){

}

void BaseDriver::GetDerivatives(const G4FieldTrack& track, G4double dydx[] ){
    G4double y[G4FieldTrack::ncompSVEC];
    track.DumpToArray(y);
    fequation->RightHandSide(y,dydx);
}

G4double BaseDriver::ComputeNewStepSize(double errMaxNorm,
                                        double hstepCurrent ){
    return hstepCurrent;
}
