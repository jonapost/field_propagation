#include "G4ModifiedMidpoint.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldTrack.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4FieldUtils.hh"

#include "boost/numeric/odeint/stepper/modified_midpoint.hpp"

#include <memory>
#include <array>
#include <iostream>

const size_t INTEGRATED_COMPONENTS = 6;
const size_t TOTAL_COMPONENTS = G4FieldTrack::ncompSVEC;

const size_t NUMBER_OF_INTEGRATION_STEPS = 1000;
const size_t MAX_MIDPOINT_STEPS = 6;

using BoostState = std::array<double, INTEGRATED_COMPONENTS>;
using BoostMidpoint = typename boost::numeric::odeint::modified_midpoint<BoostState>;
using BoostMidpointDense = typename boost::numeric::odeint::modified_midpoint_dense_out<BoostState>;
using BoostEquation = std::function<void(const BoostState&, BoostState&, double)>;

using GeantState = double[TOTAL_COMPONENTS];
using GeantMidpoint = typename internal::G4ModifiedMidpoint;
using GeantEquation = std::shared_ptr<G4EquationOfMotion>;

void copy(GeantState& dst, const BoostState& src)
{
    memcpy(dst, src.data(), sizeof(double) * INTEGRATED_COMPONENTS);
}

void copy(BoostState& dst, const GeantState& src)
{
    memcpy(dst.data(), src, sizeof(double) * INTEGRATED_COMPONENTS);
}

void copy(GeantState& dst, const GeantState& src)
{
    memcpy(dst, src, sizeof(double) * TOTAL_COMPONENTS);
}

void copy(BoostState& dst, const BoostState& src)
{
    memcpy(dst.data(), src.data(), sizeof(double) * INTEGRATED_COMPONENTS);
}

template <typename T1, typename T2>
bool isEqual(const T1& container1, const T2& container2)
{
    for (size_t i = 0; i < INTEGRATED_COMPONENTS; ++i) {
        if (container1[i] != container2[i]) {
            return false;
        }
    }

    return true;
}

void reportDiverged(
    const GeantState& geant,
    const BoostState& boost,
    double curveLength)
{
    std::cerr << "geant midpoint and boost midpoint have diverged! \n";
    std::cout << "curveLength: " << curveLength << "\n";
    std::cout << "geant: " <<
        field_utils::makeVector(
            geant, field_utils::Value3D::Position) << "\n";
    std::cout << "boost: " <<
        field_utils::makeVector(
            boost, field_utils::Value3D::Position) << "\n";
}

void test(
    GeantMidpoint& geantMidpoint,
    const GeantEquation& geantEquation,
    const GeantState& geantState,
    BoostMidpoint& boostMidpoint,
    const BoostEquation& boostEquation,
    const BoostState& boostState,
    int nsteps)
{
    geantMidpoint.SetSteps(nsteps);
    boostMidpoint.set_steps(nsteps);

    GeantState geantY, geantDydx;
    copy(geantY, geantState);

    BoostState boostY, boostDydx;
    copy(boostY, boostState);

    double curveLength = 0;
    double stepLength = 1. * CLHEP::mm;

    for (size_t i = 0; i < NUMBER_OF_INTEGRATION_STEPS; ++i) {
        geantEquation->RightHandSide(geantY, geantDydx);
        boostEquation(boostY, boostDydx, curveLength);

        geantMidpoint.DoStep(geantY, geantDydx, geantY, stepLength);
        boostMidpoint.do_step_impl(
            boostEquation, boostY, boostDydx, curveLength, boostY, stepLength);
        
        curveLength += stepLength;

        if (!isEqual(geantY, boostY)) {
            reportDiverged(geantY, boostY, curveLength);
        }
    }
}

void testDense(
    GeantMidpoint& geantMidpoint,
    const GeantEquation& geantEquation,
    const GeantState& geantState,
    BoostMidpointDense& boostMidpointDense,
    const BoostEquation& boostEquation,
    const BoostState& boostState,
    int nsteps)
{
    geantMidpoint.SetSteps(nsteps);
    boostMidpointDense.set_steps(nsteps);
    
    GeantState geantY, geantDydx, geantYmid, geantDydxMid[MAX_MIDPOINT_STEPS];
    copy(geantY, geantState);
    
    BoostState boostY, boostDydx, boostYmid;
    BoostMidpointDense::deriv_table_type boostDydxMid(MAX_MIDPOINT_STEPS);
    copy(boostY, boostState);
    
    double curveLength = 0;
    double stepLength = 1. * CLHEP::mm;
    
    for (size_t i = 0; i < NUMBER_OF_INTEGRATION_STEPS; ++i) {
        geantEquation->RightHandSide(geantY, geantDydx);
        boostEquation(boostY, boostDydx, curveLength);
        
        geantMidpoint.DoStep(geantY, geantDydx, geantY, stepLength, geantYmid, geantDydxMid);
        boostMidpointDense.do_step(
            boostEquation, boostY, boostDydx, curveLength, boostY, stepLength, boostYmid, boostDydxMid);
        
        curveLength += stepLength;
        
        if (!isEqual(geantY, boostY)) {
            reportDiverged(geantY, boostY, curveLength);
        }
        if (!isEqual(geantYmid, boostYmid)) {
            reportDiverged(geantYmid, boostYmid, curveLength);
        }
    }
}

int main()
{
    G4DynamicParticle dynParticle(
            G4Proton::Definition(),
            G4ThreeVector(1, 0, 2).unit(),
            0.01 * CLHEP::GeV);

    auto track =
        std::make_shared<G4FieldTrack>(
            G4ThreeVector{0., 1., 0.}, // start position
            0, //LaboratoryTimeOfFlight
            dynParticle.GetMomentumDirection(),
            dynParticle.GetKineticEnergy(),
            dynParticle.GetMass(),
            dynParticle.GetCharge(),
            dynParticle.GetPolarization());

    auto field = std::make_shared<G4QuadrupoleMagField>(1 * CLHEP::tesla);
    auto equation = std::make_shared<G4Mag_UsualEqRhs>(field.get());

    equation->SetChargeMomentumMass(
            {
                dynParticle.GetCharge(),
                dynParticle.GetSpin(),
                dynParticle.GetMagneticMoment()
            },
            dynParticle.GetMomentum().mag(),
            dynParticle.GetMass());

    auto system = [equation, track](const BoostState& yIn, BoostState& dydxOut, double /*t*/)
    {
        GeantState y;
        GeantState dydx;

        //init y[0...11]
        track->DumpToArray(y);

        //copy yIn[0...5]
        copy(y, yIn);

        equation->RightHandSide(y, dydx);

        copy(dydxOut, dydx);
    };

    GeantMidpoint geantMidpoint(equation.get(), INTEGRATED_COMPONENTS);
    GeantState geantState;
    track->DumpToArray(geantState);

    BoostMidpoint boostMidpoint;
    BoostState boostState;
    copy(boostState, geantState);
    
    for (int steps : {2, 4, 6}) {
        test(geantMidpoint, equation, geantState,
             boostMidpoint, system, boostState,
             steps);
    }
    
    BoostMidpointDense boostMidpointDense;
    
    for (int steps : {2, 4, 6}) {
        testDense(geantMidpoint, equation, geantState,
                  boostMidpointDense, system, boostState,
                  steps);
    }

    return 0;
}

