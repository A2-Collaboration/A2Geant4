// He3 gas scintillator active target
// Original coding B.Strandberg 2013-2015
// Modified for A2Geant-Master
// J.R.M. Annand 18th June 2020
// J.R.M. Annand 29th June 2020 Add WLS plates option
//
#ifndef A2ActiveHe3_h
#define A2ActiveHe3_h 1

//includes
#include "A2Target.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4Region.hh"
#include "A2SD.hh"
#include "A2VisSD.hh"


class A2ActiveHe3: public A2Target
{
public:
  
  A2ActiveHe3();
  ~A2ActiveHe3();

  virtual G4VPhysicalVolume* Construct(G4LogicalVolume *MotherLogical, G4double=0);

  //functions used with macros
  void SetMakeMylarSections(G4int makesections) {fMakeMylarSections = makesections;}
  void SetOpticalSimulation(G4int optsim) {fOpticalSimulation = optsim; }
  G4int GetOpticalSimulation() { return fOpticalSimulation; }

  void SetScintillationYield(G4double scintyield);

  //used for communication with stepping action
  void SetMakeEpoxy(G4bool makeepoxy) {fMakeEpoxy = makeepoxy;}
  G4bool GetMakeEpoxy() {return fMakeEpoxy;}

  //functions that build parts of the detector
  void MakeVessel();
  void MakeWLS();
  void MakePCB();
  void MakeTeflonLayer();
  void MakeSiPMTs(){}
  void MakeMylarSections(G4int nrofsections);
  void DefineMaterials();
  void SetOpticalProperties();
  void SetIsOverlapVol(G4int isOv){ fIsOverlapVol = isOv; }
  void ReadParameters(const char*);
  void MakeSensitiveDetector();

  //this places all sub-parts into fMyLogic
  void PlaceParts();
  G4LogicalVolume* GetHeInsideTeflonLogic(){ return fHeInsideTeflonLogic; }

private:

  G4NistManager* fNistManager;
  G4int fIsOverlapVol;

  G4Region* fRegionAHe3;
  A2SD* fAHe3SD;
  A2VisSD* fAHe3VisSD;

  //------------------------------------------------------------------------
  //booleans controlling detector construction and simulation
  //------------------------------------------------------------------------
  G4int fMakeMylarSections;
  G4int fOpticalSimulation;
  G4bool fMakeEpoxy; //determine whether epoxy layer on top of pmt's is build or not
  G4double fScintYield;
  //
  G4int fIsWLS;

  //------------------------------------------------------------------------
  //volumes
  //------------------------------------------------------------------------
  //Logical and physical volumes that are part of every detector class
  G4LogicalVolume* fMotherLogic;        //Logical volume of the mother
  G4LogicalVolume* fMyLogic;            //Logical Volume for this detector
  G4VPhysicalVolume* fMyPhysi;          //Physical volume for this detector

  //Logical volumes of this specific detector class
  G4LogicalVolume* fVesselLogic;        //mother for Aluminum vessel and Be windows
  G4LogicalVolume* fHeOutsideTeflonLogic; //mother to fPCMLogic
  G4LogicalVolume* fWLSLogic;           //WLS plate logic
  G4LogicalVolume* fPCBLogic;           //mother for Printed circuit board
  G4LogicalVolume* fTeflonLogic;        //mother for Teflan cylinder and ends, mylar
  G4LogicalVolume* fHeInsideTeflonLogic; //mother for pmt's and epoxy layers

  G4LogicalVolume* fTeflonCylLogic;        //Teflan cylinder (need for optical properties)
  G4LogicalVolume* fTeflonCylEndLogic;     //Teflan cylinder end (need for optical properties)
  G4LogicalVolume* fPMTLogic;           //logic of single pmt cell (need for optical properties)
  G4LogicalVolume* fEpoxyLogic;         //logic of the epoxy layer (need for optical properties)
  G4LogicalVolume* fMylarLogic;         //logic of the Mylar window (need for optical properties)
  G4LogicalVolume* fMylarSectionLogic;  //logic for the mylar section window inside main cell

  G4VPhysicalVolume **fPMTPhysic;   //physical volumes to set LogicalBorderSurface
  G4VPhysicalVolume **fEpoxyPhysic;
  //------------------------------------------------------------------------
  //geometric parameters
  //------------------------------------------------------------------------

  //geometric parameters for MakeVessel()
  G4double fHeContainerR; //inner radius of the helium container
  G4double fHeContainerZ; //length of the helium container
  G4double fExtensionR; //inner radius of the extension tubes
  G4double fExtensionZU; //length upstream extension tube
  G4double fExtensionZD; //length downstream extension tube
  G4double fContainerThickness; //thickness of the metal vessel
  G4double fBeThickness; //thickness of the beryllium end windows

  // geometric parameters for MakeWLS()
  G4int fNwls;            // phi segmentation
  G4double fWLSthick;     // plate thickness
  G4double fRadClr;       // radial clearance to vessel
  G4double fLatClr;       // lateral clearance to vessel
  G4double fRwls1;        // radial pos. plate centre
  
  //geometric parameters for MakePCB()
  G4double fPCBThickness; //thickness of the printed circuit board
  G4double fPCBRadius;    //inner radius of the printed circuit board
  G4double fPCBZ;         //length of th pcb

  //geometric parameters for MakeTeflanLayer()
  G4double fTeflonThicknessEnd;  //thickness of teflan at end walls
  G4double fTeflonThicknessCyl;  //thickness of teflan covering the pcb cylinder
  G4double fTeflonR;       //inner radius of the teflan covering the pcb
  G4double fTeflonZ;
  G4double fMylarThickness;
  /*These parameters are used to create a flattening inside the teflon cylinder
   where the pmt's are attached. Without the flattening placing pmt's is a pain*/
  G4double fFlatteningWidth;
  G4double fFlatteningHeight;  

  //geometric parameters for MakeSiPMTs()
  G4double fPMTx;
  G4double fPMTy;
  G4double fPMTz;
  G4double fEpoxyz;
  G4double fPmtR; //radius from cylinder center to pmt center
  G4double fPMTEpoxyGap;
  G4double fEpoxyR; //radius from cylinder center to epoxy center

  G4int fNxplane; //nr of pmt's along x plane
  G4int fNyplane; //nr of pmt's along y plane
  G4int fNdiag1;
  G4int fNdiag2;
  G4int fNpmttotal;  //total nr of pmts

  G4double fOffsetx; //the center pos of first pmt from the teflon end wall x plane
  G4double fOffsety; //the center pos of first pmt from the teflon end wall y plane
  G4double fOffsetd1;
  G4double fOffsetd2;

  G4double fStepx; //step between pmt's in x plane
  G4double fStepy; //step between pmt's in y plane
  G4double fStepd1;
  G4double fStepd2;

  //geometric parameters for MakeMylarSections()
  G4double fMylarSecZ; //thickness of the sectioning windows
  G4double fMylarWinStep; //distance between windows
  G4double fMylarWinOffset; //the center pos of first window from the teflon end wall
  G4int fNMylarSecWins;  //number of sectioning windows  

} ;

#endif

