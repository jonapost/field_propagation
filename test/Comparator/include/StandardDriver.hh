#pragma once

#include "G4Types.hh"
#include "G4Mag_UsualEqRhs.hh"

#include <memory>

template<class Stepper>
class StandardDriver {
public:
    StandardDriver(G4double minimumStep, std::shared_ptr<G4Mag_UsualEqRhs> equation);
    void OneGoodStep(
        G4double y[], const G4double dydx[], G4double& trackLength, G4double htry,
            G4double epsRelMax, G4double& hdid, G4double& hnext);
    void RightHandSide(G4double y[], G4double dydx[]){
        fstepper.RightHandSide(y,dydx);
    }

    Stepper& stepper() {return fstepper;}
    void setSafety(G4double safety) {fsafety = safety;}
    void setTrials(G4int trials) {ftrials = trials;}
    void setMinimumStep(G4double minimumStep) {fminimumStep = minimumStep;}
    void setPShrink(G4double pShrink) {fpShrink = pShrink;}
    void setPGrow(G4double pGrow) {fpGrow = pGrow;}
    void setMaxStepIncrease(G4double maxStepIncrease)
        {fmaxStepIncrease = maxStepIncrease;}
private:

    std::shared_ptr<G4Mag_UsualEqRhs> fequation;

    Stepper fstepper;
    G4int ftrials;
    G4double fsafety;
    G4double fminimumStep;
    G4double fpShrink;
    G4double fpGrow;
    G4double fmaxStepIncrease;
    G4double fErrcon;
};

#include "StandardDriver.icc"
