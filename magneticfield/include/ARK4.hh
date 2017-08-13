//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
// -------------------------------------------------------------------

#ifndef ARK4_HH
#define ARK4_HH

#include "G4MagErrorStepper.hh"

class ARK4 : public G4MagErrorStepper {
public:
    ARK4(G4EquationOfMotion* equation,
        G4int numberOfVariables = 6);

    ARK4(const ARK4&) = delete;
    ARK4& operator =(const ARK4&) = delete;

    virtual void DumbStepper(
        const G4double yIn[],
        const G4double dydx[],
        G4double hstep,
        G4double yOut[]) override;

    virtual G4int IntegratorOrder() const override
    {
        return 4;
    }

private:
    G4double particleCharge();

    // oscillation frequency * time step
    G4double integrationParameter(
        G4double charge,
        G4double momentum,
        G4double Bfield,
        G4double hstep);
};

#endif
