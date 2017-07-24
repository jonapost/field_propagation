
#ifndef LOGICAL_VOLUME_HANDLER
#define LOGICAL_VOLUME_HANDLER

#include <vector>
#include <map>

#include "globals.hh"

class G4LogicalVolume;

class LogicalVolumeHandler{
 public:
   static LogicalVolumeHandler* Instance();
 
   G4int  GetIndex(G4LogicalVolume*);
   G4bool IsUsed(G4LogicalVolume*);
   G4bool IsUsed(G4int);

 private:
   LogicalVolumeHandler();
   void Initialize(); 
   LogicalVolumeHandler(const  LogicalVolumeHandler &);
   LogicalVolumeHandler &operator=(const  LogicalVolumeHandler &);


 private:
   static LogicalVolumeHandler* fgInstance;
   // temporary solution (till Geant4-10.2.beta) to get rid of the used logical // 
   // volume check that can take 8-10 minutes                                   //
   static const G4bool fgIsLogicalVolumeUsedMask[];         /** Mask for used logical volumes */
   std::map<G4LogicalVolume*,G4int> fLogicVolToIndex;       /** Map G4LogicalVolume* -> its index in the store */
 };
#endif
