#include "KurganAnalysis.hh"
#include "KurganRunAction.hh"
#include "KurganEventAction.hh"
#include "KurganPrimaryGeneratorAction.hh"
#include "KurganDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "comments.hh"

KurganRunAction::KurganRunAction()
: G4UserRunAction()
{ 
  // Register accumulable to the accumulable manager
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->RegisterAccumulable(fStrawNumber); 
  //accumulableManager->RegisterAccumulable(fPositionDet); 
  //accumulableManager->RegisterAccumulable(fPositionDet2); 

  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Creating histograms 
  analysisManager->CreateH1("Layer_Number","tube number", 15, 0., 15.); 
  //analysisManager->CreateH1("fCounter1","Counter 1", 100, 0., 100.); 
  //analysisManager->CreateH1("fCounter2","Counter 2", 100, 0., 100.); 
	
	//Creating 2D histograms 
  analysisManager->CreateH2("fPositionDet","fPositionDet", 130, -75.*cm, 75.*cm, 130, -75.*cm, 75.*cm); 
  analysisManager->CreateH2("fPositionDet2","fPositionDet2", 130, -75.*cm, 75.*cm, 130, -75.*cm, 75.*cm); 

  analysisManager->CreateNtuple("Kurgan", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("Layer_Number");
  //analysisManager->CreateNtupleDColumn("fCounter1");
  //analysisManager->CreateNtupleDColumn("fCounter2");
  analysisManager->FinishNtuple();

}

KurganRunAction::~KurganRunAction()
{
  delete G4AnalysisManager::Instance();
}

void KurganRunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  auto analysisManager = G4AnalysisManager::Instance();
  G4String fileName = "Kurgan";
  analysisManager->OpenFile(fileName);

}

void KurganRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->Merge();
  
  auto analysisManager = G4AnalysisManager::Instance();
  
  analysisManager->Write();
  analysisManager->CloseFile();

}
