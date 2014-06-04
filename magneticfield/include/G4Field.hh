/*
* this single file contains all the code for G4Field
* will add G4Field.cc later
* experiment only
*/

#ifndef G4FIELD_HH
#define G4FIELD_HH

#include "G4Types.hh"
#include "globals.hh"
template
<class Field>
class G4Field
{
  public:  

       void  GetFieldValue( const  double Point[4],
                                          double *fieldArr ) const
	{return fField->GetFieldValue(Point[4], *fieldArr);}

      G4Field( G4bool gravityOn= false);
      G4Field( const G4Field & );
      ~G4Field();
      inline G4Field& operator = (const G4Field &p); 
 
      G4bool   DoesFieldChangeEnergy() const
      {return fField->DoesFieldChangeEnergy();} 

      G4bool   IsGravityActive() const 
      {return fGravityActive;}
      inline void SetGravityActive( G4bool OnOffFlag );
      //G4Field* Clone() const
      //Implements cloning, needed by G4 MT
  private:
      G4bool  fGravityActive;
      Field*  fField;
};

inline void  G4Field::SetGravityActive( G4bool OnOffFlag )
{ 
  fGravityActive= OnOffFlag; 
} 
#endif /* G4FIELD_HH */
