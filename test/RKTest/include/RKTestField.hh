#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"


//G4QuadrupoleMagField SetupQMF(G4double gradient);
//G4UniformMagField SetupUMF( );
G4QuadrupoleMagField SetupQMF(G4double factor = 1., G4double gradient = 0.01*CLHEP::tesla);
G4UniformMagField SetupUMF(G4double factor = 1. ,
                           G4double x_field = 0.*CLHEP::tesla,
                           G4double y_field = 0.*CLHEP::tesla,
                           G4double z_field = 0.1*CLHEP::tesla);

G4double toStepLen(G4double theta, G4double R = 0);

G4double findRFrom2pt(G4double yOut2[], G4double yOut1[], G4double step_len);



