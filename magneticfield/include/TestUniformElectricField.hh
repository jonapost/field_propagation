#ifndef TESTUNIFORMELECTRICFIELD_HH
#define TESTUNIFORMELECTRICFIELD_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
//#include "G4ElectricField.hh"
#include "G4Field.hh"

class TestUniformElectricField : public G4Field
//: public G4ElectricField
{
  public: 
    TestUniformElectricField(const G4ThreeVector FieldVector);
      // A field with value equal to FieldVector.

    TestUniformElectricField(G4double vField,
                           G4double vTheta,
                           G4double vPhi     ) ;
       
    virtual ~TestUniformElectricField() ;

    TestUniformElectricField(const TestUniformElectricField &p);
    TestUniformElectricField& operator = (const TestUniformElectricField &p);
      // Copy constructor and assignment operator

    virtual void GetFieldValue(const G4double pos[4], G4double *field) const;

    virtual TestUniformElectricField* Clone() const;
  private:
  
    G4double fFieldComponents[6] ;
};

#endif
