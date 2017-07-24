#ifndef MaterialConverter_hh
#define MaterialConverter_hh 1

#include <vector>
class TGeoManager;

class MaterialConverter
{
 public:
   static MaterialConverter* Instance();
   static void SetTGeomManager(TGeoManager *tgeomgr){ fgTGeomMgr = tgeomgr;}
   int GetROOTMaterialID(int g4materialindx){return  g4materialindx;} // very simple now   

private:
   MaterialConverter(); 
   ~MaterialConverter();
   MaterialConverter(MaterialConverter const&);
   MaterialConverter& operator=(MaterialConverter const&);

   void Initialize();
    
 private:
   static MaterialConverter *fgInstance;
   static TGeoManager *fgTGeomMgr;
  
//   std::vector<int> fRootMatIndices; // [ key =  G4  Mat index ]
//   std::vector<int> fG4MatIndices;   // [ key = Root Mat index ]

};

#endif // MaterialConverter_hh 1
