#ifndef GUSTAFSSONDRIVER_HH
#define GUSTAFSSONDRIVER_HH

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"

#include <memory>

class GustafssonDriver {
public:
    GustafssonDriver(G4double minimumStep, G4MagIntegratorStepper* stepper);

    void OneGoodStep(G4double y[],
                     const G4double dydx[],
                     G4double& trackLength,
                     G4double htry,
                     G4double epsRelMax,
                     G4double& hdid,
                     G4double& hnext);

    const G4MagIntegratorStepper* GetStepper() const {return fstepper;}

    void setSafety(G4double safety) {fsafety = safety;}
    void setTrials(G4int trials) {ftrials = trials;}
    void setMinimumStep(G4double minimumStep) {fminimumStep = minimumStep;}
    void setMaxStepIncrease(G4double maxStepIncrease)
        {fmaxStepIncrease = maxStepIncrease;}

private:
    G4MagIntegratorStepper* fstepper;
    G4int ftrials;
    G4double fsafety;
    G4double fminimumStep;
    G4double KI;
    G4double KP;
    G4double fmaxStepIncrease;
};

#endif
