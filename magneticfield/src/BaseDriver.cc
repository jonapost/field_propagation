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


G4double BaseDriver::ComputeNewStepSize(double errMaxNorm,
                                        double hstepCurrent ){
    return hstepCurrent;
}

void BaseDriver::SetVerboseLevel(G4int verb){
    fverb = verb;
}

G4double BaseDriver::calcError(const state_type &y1, const state_type &y2, G4double hstep){

    G4double yErr[G4FieldTrack::ncompSVEC];
    for (int i = 0; i < fnvar; ++i){
        yErr[i] = y1[i] - y2[i];
    }

    G4double errPos2 = sqr(yErr[0]) + sqr(yErr[1]) + sqr(yErr[2]);
    G4double errMom2 = sqr(yErr[3]) + sqr(yErr[4]) + sqr(yErr[5]);
    G4double Mom2 = sqr(y1[3]) + sqr(y1[4]) + sqr(y1[5]);
    errMom2 /= Mom2;
    errPos2 /= (hstep*hstep);

    return std::max(sqrt(errPos2), sqrt(errMom2)*hstep);
}
