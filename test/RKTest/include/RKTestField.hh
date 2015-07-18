#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"


//G4QuadrupoleMagField SetupQMF(G4double gradient);
//G4UniformMagField SetupUMF( );
G4QuadrupoleMagField SetupQMF(G4double gradient = 0.01*CLHEP::tesla);
G4UniformMagField SetupUMF(G4double x_field = 0.*CLHEP::tesla,
                           G4double y_field = 0.*CLHEP::tesla,
                           G4double z_field = 0.1*CLHEP::tesla);
