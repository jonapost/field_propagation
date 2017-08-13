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
//
// -------------------------------------------------------------------

#include "ARK4.hh"

#include "G4Mag_EqRhs.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4FieldTrack.hh"
#include "BettisFunctions.hh"
#include "Utils.hh"

namespace {

const G4double UNDEFINDED_CHARGE = 0;

G4double averageField(
    const G4double field2[],
    const G4double field3[],
    const G4double field4[])
{
    static const G4int NUMBER_OF_COMPONENTS = 3;

    G4double fieldSquared = 0;
    for (G4int i = 0; i < NUMBER_OF_COMPONENTS; ++i) {
        fieldSquared += sqr((field2[i] + field3[i] + field4[i]) / 3);
    }

    return sqrt(fieldSquared);
}

} // namespace

ARK4::ARK4(
    G4EquationOfMotion* equation,
    G4int numberOfVariables)
    : G4MagErrorStepper(equation, numberOfVariables)
{
}

void ARK4::DumbStepper(
    const G4double yIn[],
    const G4double dydx[],
    G4double hstep,
    G4double yOut[])
{
    const G4double halfStep = hstep * 0.5;

    G4double ytemp[G4FieldTrack::ncompSVEC];

    G4double ak2[G4FieldTrack::ncompSVEC];
    G4double ak3[G4FieldTrack::ncompSVEC];
    G4double ak4[G4FieldTrack::ncompSVEC];

    //TODO use G4EquationOfMotion::G4maximum_number_of_field_components = 24 ?
    G4double field2[3];
    G4double field3[3];
    G4double field4[3];

    // Initialise time to t0, needed when it is not updated by the integration.
    //        [ Note: Only for time dependent fields (usually electric)
    //                  is it neccessary to integrate the time.]
    magneticfield::copyValue(ytemp, yIn, magneticfield::Value1D::LabTime);
    magneticfield::copyValue(yOut, yIn, magneticfield::Value1D::LabTime);

    //RightHandSide(yIn, dydx);
    for (G4int i = 0; i < GetNumberOfVariables(); ++i) {
        ytemp[i] = yIn[i] + halfStep * dydx[i];
    }
    GetEquationOfMotion()->EvaluateRhsReturnB(ytemp, ak2, field2);

    for (G4int i = 0; i < GetNumberOfVariables(); ++i) {
        ytemp[i] = yIn[i] + halfStep * ak2[i];
    }
    GetEquationOfMotion()->EvaluateRhsReturnB(ytemp, ak3, field3);

    for (G4int i = 0; i < GetNumberOfVariables(); ++i) {
        ytemp[i] = yIn[i] + hstep * ak3[i];
    }
    GetEquationOfMotion()->EvaluateRhsReturnB(ytemp, ak4, field4);

    const G4double nu = integrationParameter(
            particleCharge(),
            magneticfield::extractValue(yIn, magneticfield::Value3D::Momentum),
            averageField(field2, field3, field4),
            hstep);

    const G4double phi3 = magneticfield::bettisFunction<3>(nu);
    const G4double phi4 = magneticfield::bettisFunction<4>(nu);

    const G4double b1 = 4 * phi4;
    const G4double b2 = 1 - 4 * phi3;
    const G4double b3 = 4 * phi3 - 8 * phi4;
    const G4double b4 = 4 * phi4;

    for (G4int i = 0; i < GetNumberOfVariables(); ++i) {
        yOut[i] = yIn[i] + hstep * (b1 * dydx[i] + b2 * ak2[i] +
                                    b3 * ak3[i] + b4 * ak4[i]);
    }
}

// TODO refactor G4EquationOfMotion::SetChargeMomentumMass
G4double ARK4::particleCharge()
{
    if (auto equation = dynamic_cast<G4Mag_EqRhs*>(GetEquationOfMotion())) {
        return equation->FCof() / (eplus * c_light);
    }
    return UNDEFINDED_CHARGE;
}

G4double ARK4::integrationParameter(
    G4double charge,
    G4double momentum,
    G4double Bfield,
    G4double hstep)
{
    return hstep * charge * Bfield * c_light / momentum;
}
