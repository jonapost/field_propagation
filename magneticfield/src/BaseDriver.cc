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
    G4double hnew;

    // Compute size of next Step for a failed step
    if(errMaxNorm > 1.0 )
    {
      // Step failed; compute the size of retrial Step.
      hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPshrnk()) ;
    } else if(errMaxNorm > 0.0 ) {
      // Compute size of next Step for a successful step
      hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPgrow()) ;
    } else {
      // if error estimate is zero (possible) or negative (dubious)
      hnew = max_stepping_increase * hstepCurrent;
    }

    return hnew;
}
