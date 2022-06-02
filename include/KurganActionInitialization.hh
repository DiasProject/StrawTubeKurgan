#ifndef KurganActionInitialization_h
#define KurganActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.
class KurganDetectorConstruction;

class KurganActionInitialization : public G4VUserActionInitialization
{
  public:
    KurganActionInitialization(KurganDetectorConstruction*);
    virtual ~KurganActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    KurganDetectorConstruction* fDetConstruction;
};
#endif  
