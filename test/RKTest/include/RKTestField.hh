#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"

//G4QuadrupoleMagField SetupQMF(G4double gradient);
//G4UniformMagField SetupUMF( );
G4QuadrupoleMagField SetupQMF(G4double gradient = 0.01*CLHEP::tesla);
G4UniformMagField SetupUMF(G4double x_field = 0.*CLHEP::tesla,
                           G4double y_field = 0.*CLHEP::tesla,
                           G4double z_field = 0.1*CLHEP::tesla);
//G4double x_field ,
//G4double y_field ,
//G4double z_field

//#ifndef RKTEST_FIELD
//#define RKTEST_FIELD
//
//template<class MAG_FIELD>
//class RKTestField{
//private:
//    //members
//    MAG_FIELD theField;
//    G4double x_field, y_field, z_field; 	//for UMF
//    G4double gradient;					//For QMF
//public:
//    //member functions
//    RKTestField();
//    RKTestField(MAG_FIELD MField);
////    RKTestField(G4UniformMagField UMF);
//    ~RKTestField();
////    G4UniformMagField SetupUMF();
////    G4QuadrupoleMagField SetupQMF();
//};
//#endif //(Defined RKTestField