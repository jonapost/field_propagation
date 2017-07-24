#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

using namespace CLHEP;


SteppingAction::SteppingAction(): G4UserSteppingAction(),
    fpField(nullptr),ncalls(0),fParticleGun(nullptr),time(0),
    out("/home/Dmitry/work/GSoC/field_propagation/test/SimplePropagation/log.txt")
{
    G4RunManager* runManager = G4RunManager::GetRunManager();

    fpField = static_cast<const DetectorConstruction*>
            (runManager->GetUserDetectorConstruction())->GetField();

    fParticleGun = static_cast<const  PrimaryGeneratorAction*>
            (runManager->GetUserPrimaryGeneratorAction())->GetParticleGun();

}



SteppingAction::~SteppingAction()
{
    out.close();
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    //B = 1 tesla, along Y axe.
    G4double BField = 1*tesla;

    const G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
    G4double mass = particle->GetPDGMass();
    G4double charge = particle->GetPDGCharge();

    G4double energy = fParticleGun->GetParticleEnergy();
    G4double mom2 = energy*(energy + 2*mass);
    G4ThreeVector momDirection = fParticleGun->GetParticleMomentumDirection();
    G4ThreeVector momentum = sqrt(mom2)*momDirection;

    G4ThreeVector position = fParticleGun->GetParticlePosition();

    G4double betta2 = mom2/(mass*mass + mom2);
    G4double gamma = 1./sqrt(1-betta2);

    G4double velocity = sqrt(betta2)*c_light;

    //G4cout<<"betta: "<<sqrt(betta2)<<G4endl;

    G4double stepLength = step->GetStepLength();

    G4double timeStep = stepLength/velocity;

    time += timeStep;

    G4double omega = charge*BField/(gamma*mass/c_squared);

    G4double momXZ = sqrt(sqr(momentum.x()) + sqr(momentum.z()))*c_light;

    G4double radius = momXZ/(charge*BField*c_squared);

    //G4cout<<"radius: "<<radius<<G4endl;

    G4ThreeVector post = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector pos_calc;


    //calc.setY(pre.y() + velocityVector.y()*timeStep);
    pos_calc.setY((position.getY() + velocity*momDirection.y()*time));
    G4double y_error2 = sqr(pos_calc.y() - post.y());


    //calc.setX(pre.x() + radius*sin(omega*timeStep));
    pos_calc.setX(radius*sin(omega*time));
    G4double x_error2 = sqr(pos_calc.x() - post.x());

    //calc.setZ(pre.z() - radius*cos(omega*timeStep));
    pos_calc.setZ(-radius*cos(omega*time));
    G4double z_error2 = sqr(pos_calc.z() - post.z());

    //G4cout<<2*pi/omega*velocityVector.y()<<G4endl;

    G4double pos_error = sqrt((x_error2 + y_error2 + z_error2)/*/sqr(stepLength)*/);

    /*
    if (error > 1e-5){
        G4cout<<"pre: "<<pre<<"\npost: "<<post<<"\ncalc: "<<calc<<G4endl;
    }*/


    if (/*error > perMillion*/ false){
        G4cout<<"ycalc: "<<pos_calc.y()<<" yget: "<<post.y()<<" y_error2: "<<y_error2<<G4endl;
        G4cout<<"xcalc: "<<pos_calc.x()<<" xget: "<<post.x()<<" x_error2: "<<x_error2<<G4endl;
        G4cout<<"zcalc: "<<pos_calc.z()<<" zget: "<<post.z()<<" z_error2: "<<z_error2<<G4endl;
    }



    //momentum error
    G4ThreeVector postMomentum = step->GetPostStepPoint()->GetMomentum();
    G4ThreeVector momCalc;
    momCalc.setY(velocity*momDirection.getY()*gamma*mass/c_light);
    G4double py_error2 = sqr(momCalc.getY() - postMomentum.getY());

    momCalc.setX(radius*omega*cos(omega*time)*gamma*mass/c_light);
    G4double px_error2 = sqr(momCalc.getX() - postMomentum.getX());

    momCalc.setZ(radius*omega*sin(omega*time)*gamma*mass/c_light);
    G4double pz_error2 = sqr(momCalc.getZ() - postMomentum.getZ());

    G4double mom_error = sqrt((px_error2 + py_error2 + pz_error2));


    G4cout<<"hstep: "<<stepLength<<" pos_error: "
          <<pos_error<<" mom_error: "<<mom_error<<G4endl;



    //print number of calls to field
    G4int TotCalls = fpField->GetCountCalls();
    G4int StepCalls = TotCalls - ncalls;
    ncalls = TotCalls;
    G4cout<<"calls to field: "<<StepCalls<<G4endl;


    out<<time*velocity<<"    "<<pos_error<<"    "<<mom_error<<"    "<<StepCalls<<G4endl;
  }



