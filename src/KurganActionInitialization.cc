#include "KurganActionInitialization.hh"
#include "KurganPrimaryGeneratorAction.hh"
#include "KurganRunAction.hh"
#include "KurganEventAction.hh"
#include "KurganSteppingAction.hh"

KurganActionInitialization::KurganActionInitialization(KurganDetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
  fDetConstruction(detConstruction)
{}

KurganActionInitialization::~KurganActionInitialization()
{}

void KurganActionInitialization::BuildForMaster() const
{
  KurganRunAction* runAction = new KurganRunAction;
  SetUserAction(runAction);
}

void KurganActionInitialization::Build() const
{
  SetUserAction(new KurganPrimaryGeneratorAction);

  KurganRunAction* runAction = new KurganRunAction;
  SetUserAction(runAction);
  
  KurganEventAction* eventAction = new KurganEventAction(runAction);
  SetUserAction(eventAction);
  
  SetUserAction(new KurganSteppingAction(fDetConstruction, eventAction));
}  
