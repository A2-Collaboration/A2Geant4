#include "G4RunManager.hh"
#include "A2SD.hh"
#include "A2Hit.hh"
#include "A2EventAction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include "stdio.h"

using namespace CLHEP;

A2SD::A2SD(G4String name,G4int Nelements):G4VSensitiveDetector(name)
{
  collectionName.insert(G4String("A2SDHits")+name);
  fCollection=NULL;

  fNelements=Nelements+1;//numbering starts from 1 not 0
  fhitID=new G4int[fNelements];
  for(G4int i=0;i<fNelements;i++)fhitID[i]=-1;
  fHits=new G4int[fNelements];
  for(G4int i=0;i<fNelements;i++)fHits[i]=0;
 
  fNhits=0;
  fHCID=-1;
}


A2SD::~A2SD()
{

}


void A2SD::Initialize(G4HCofThisEvent*)
{
  //Hits collections are deleted in G4RunManager, DoEventLoop
  //This calls StackPreviousEvent which deletes the G4Event
  //This deletes the hit collection for the event and thus A2Hit s
  //G4cout<<"void A2SD::Initialize "<<SensitiveDetectorName<<G4endl;
  fCollection = new A2HitsCollection(SensitiveDetectorName,collectionName[0]);


}


G4bool A2SD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{ 
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ((edep/keV == 0.)) return false;      
  // This TouchableHistory is used to obtain the physical volume
  // of the hit
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  
  G4VPhysicalVolume* volume=theTouchable->GetVolume();
  G4VPhysicalVolume* mothervolume=theTouchable->GetVolume(1);
  G4int id;
  //Get element copy number
  //TAPS volume  is contained in COVR which is the multiple placed volume!
  //For PbWO4 they have an additional Copy Number which should be added on to the COVR volume
  //if(volume->GetName().contains("TAPS")||volume->GetName().contains("PbWO"))id=mothervolume->GetCopyNo()+volume->GetCopyNo();
  if(mothervolume->GetName().contains("COVR"))id=mothervolume->GetCopyNo()+volume->GetCopyNo();
  else id = volume->GetCopyNo();
  //seperate ADC gates for TAPS
  if((mothervolume->GetName().contains("COVR"))&&(aStep->GetPreStepPoint()->GetGlobalTime()>2000*ns))return false;
  else if (aStep->GetPreStepPoint()->GetGlobalTime()>600*ns)return false; 

  // energy correction for non-linearity
  //if (aStep->GetTrack()->GetMaterial()->GetName().contains("G4_PLASTIC_SC_VINYLTOLUENE"))
  //  edep = GetEffectiveEnergyDeposit(aStep);

  //if(volume->GetName().contains("Pb")) G4cout<<volume->GetName()<<" id "<<id <<" "<<mothervolume->GetCopyNo()<<" "<<volume->GetCopyNo()<<" edep "<<edep/MeV<<G4endl;
  if (fhitID[id]==-1){
    //if this crystal has already had a hit
    //don't make a new one, add on to old one.   
    // G4cout<<"Make hit "<<fCollection<<G4endl;    
    A2Hit* Hit = new A2Hit;
    Hit->AddEnergy(edep);
    Hit->SetPos(aStep->GetPreStepPoint()->GetPosition());
    Hit->SetID(id);
    Hit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
    fhitID[id] = fCollection->insert(Hit) -1;
    fHits[fNhits++]=id;
  }
  else // This is not new
  {
    (*fCollection)[fhitID[id]]->AddEnergy(edep);
    // set more realistic hit times
    G4double time = aStep->GetPreStepPoint()->GetGlobalTime();
    if (volume->GetName().contains("TAPS"))
    {
      if (edep/MeV > 4. && time < (*fCollection)[fhitID[id]]->GetTime())
        (*fCollection)[fhitID[id]]->SetTime(time);
    }
    else if (volume->GetName().contains("CRYSTAL"))
    {
      if (edep/MeV > 2. && time < (*fCollection)[fhitID[id]]->GetTime())
        (*fCollection)[fhitID[id]]->SetTime(time);
    }
  }
  //G4cout<<"done "<<fNhits<<G4endl;
  return true;
}


void A2SD::EndOfEvent(G4HCofThisEvent* HCE)
{
  // G4cout<<"EndOfEvent( "<<fHCID<<" "<<fCollection<<" "<<fNhits<<G4endl;
  if(fHCID<0) fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  if(fNhits>0)  HCE->AddHitsCollection(fHCID,fCollection);
  //G4cout<<"EndOfEvent( "<<G4endl;
 
  //reset hit arrays
  for (G4int i=0;i<fNhits;i++) 
    {
      fhitID[fHits[i]]=-1;
      fHits[i]=0;
    }
  fNhits=0;
  //G4cout<<"EndOfEvent( done"<<G4endl;
}



void A2SD::clear()
{} 


void A2SD::DrawAll()
{} 


void A2SD::PrintAll()
{} 

G4double A2SD::GetEffectiveEnergyDeposit(const G4Step* aStep)
{
  // Weight energy deposit by scintillation-light quenching effects
  // Scintillation response is generally non-linear with energy
  // According to Birks it depends on de/dx
  //

  // density of EJ-204
  const Double_t den = 1.023; // [g/cm^3]

  // NIM B 170 (2000) 523
  G4double kB = 1.6e-2; // [g MeV^-1 cm^-2]
  G4double C  = 5.0e-5; // [(g MeV^-1 cm^-2)^2]
  kB /= den; // [cm MeV^-1]
  C /= den;  // [(cm MeV^-1)^2]

  G4double destep = aStep->GetTotalEnergyDeposit();
  if ((!destep)) return destep;
  G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  if (!charge) return destep;
  G4double stepl = aStep->GetStepLength()/cm;
  if (!stepl) return destep;

  return destep / (1. + kB*destep/stepl + C*destep*destep/stepl/stepl);
}

