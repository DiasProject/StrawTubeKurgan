#ifndef KurganRunAction_h
#define KurganRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;

class KurganRunAction : public G4UserRunAction
{
  public:
    KurganRunAction();
    virtual ~KurganRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    //void AddStrawNumber_Run(G4int StrawNumber){Layer_Number=StrawNumber;} 
    //void AddCounter1_Run(G4int Counter1){fCounter1 = Counter1;} 
    //void AddCounter2_Run(G4int Counter2){fCounter2 = Counter2;} 
    //void AddPositionDet_Run (G4double PosX, G4double PosY){fPosX1 = PosX; fPosY1 = PosY;} 
    //void AddPositionDet2_Run(G4double PosX, G4double PosY){fPosX2 = PosX; fPosY2 = PosY;} 

  private:
    //G4int Layer_Number = 0;
    //G4int fCounter1 = 0;
    //G4int fCounter2 = 0;
    //G4double fPosX1 = 0;
    //G4double fPosY1 = 0;
    //G4double fPosX2 = 0;
    //G4double fPosY2 = 0;
};
#endif
