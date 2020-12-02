// He3 gas scintillator active target
// Original coding B.Strandberg 2013-2015
// Modified for A2Geant-Master
// J.R.M. Annand 18th June 2020
// J.R.M. Annand 29th June 2020 Add WLS plates option
//
#include "A2ActiveHe3.hh"

//G4 headers
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4NistManager.hh"    //to get elements
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "A2SD.hh" //NEW
#include "A2VisSD.hh" //NEW

//constructor
A2ActiveHe3::A2ActiveHe3() {

    //initiate all pointers to NULL

    //volumes
    fMotherLogic = NULL;
    fMyLogic = NULL;
    fMyPhysi = NULL;

    fVesselLogic = NULL;
    fHeOutsideTeflonLogic = NULL;
    fPCBLogic = NULL;
    fTeflonLogic = NULL;
    fHeInsideTeflonLogic = NULL;

    fTeflonCylLogic = NULL;
    fTeflonCylEndLogic = NULL;
    fPMTLogic = NULL;
    fEpoxyLogic = NULL;
    fMylarLogic = NULL;
    fMylarSectionLogic = NULL;

    fPMTPhysic = NULL;
    fEpoxyPhysic = NULL;

    fRegionAHe3 = new G4Region("AHe3");
    fAHe3SD = NULL;
    fAHe3VisSD = NULL;

    //initiate nist manager
    fNistManager=G4NistManager::Instance();

    fMakeMylarSections = 0;
    fOpticalSimulation = 0;
    fMakeEpoxy = 1;
    fScintYield = 65.;

    fIsOverlapVol = true;
    fIsWLS = 1;
    fNwls = 8;                       // phi segmentation
    fWLSthick = 3*mm;                // thickness WLS plate
    fRadClr = 2*mm;                  // radial clearance to vessel
    fLatClr = 5*mm;                  // lateral clearance to vessel
    ReadParameters("data/ActiveHe3.dat");
}


//destructor
A2ActiveHe3::~A2ActiveHe3() {
    //delete Rot;
    if(fRegionAHe3) delete fRegionAHe3; //remove the sensitive detector
    if(fAHe3SD) delete fAHe3SD; //remove the sensitive detector
    if(fAHe3VisSD) delete fAHe3VisSD; //remove the sensitive detector
}

/**********************************************************

This function is the main constructor function that is
called inside detector construction

**********************************************************/
G4VPhysicalVolume* A2ActiveHe3::Construct(G4LogicalVolume* MotherLogical, G4double Z0) {

    G4cout << "A2ActiveHe3::Construct() Building the active target." << G4endl;

    fMotherLogic = MotherLogical;
    DefineMaterials();
    //call functions that build parts of the detector
    MakeVessel();                    //builds and places parts inside fVesselLogic

    if(fIsWLS){
        MakeWLS();                       //builds wavelength shifting plates
    }
    else{
        MakePCB();                       //builds and places parts inside fPCBLogic
        MakeTeflonLayer();               //builds and places parts inside fTeflonLayer
        MakeSiPMTs();                    //places pmt's inside fTeflonLayer
        //
        if (fMakeMylarSections) {
            G4cout << "A2ActiveHe3::Construct() Sectioning the target." << G4endl;
            MakeMylarSections(fMakeMylarSections); //place sectioning mylar windows
        }
        else { G4cout << "A2ActiveHe3::Construct() Target not sectioned." << G4endl; }
    }
    //
    if (fOpticalSimulation) {
        G4cout << "A2ActiveHe3::Construct() Setting optical properties." << G4endl;
        SetOpticalProperties(); //set optical properties and surfaces
    }
    else { G4cout << "A2ActiveHe3::Construct() Optical properties not used." << G4endl; }
    //

    //place the separate parts into main logic of this detector
    PlaceParts();

    //construct the sensitive detector
    MakeSensitiveDetector();

    //place fMyLogic into MotherLogical
    fMyPhysi = new G4PVPlacement(0, G4ThreeVector(0,0,Z0) ,fMyLogic,"ActiveHe3",fMotherLogic,false,1,fIsOverlapVol);

    return fMyPhysi;
}


/**********************************************************

This function builds the aluminum vessel and the Be windows
and places them to member fVesselLogic

**********************************************************/
void A2ActiveHe3::MakeVessel() {
    //------------------------------------------------------------------------------
    //Define shapes (solids)
    //------------------------------------------------------------------------------

    G4Tubs *VesselCyl = new G4Tubs
            ("VesselCyl",
             0,                              //inner radius
             fHeContainerR+fContainerThickness,   //outer radius=biggest part
             fHeContainerZ/2 + fContainerThickness + fExtensionZU + fBeThickness,
             0.*deg,                         //start angle
             360.*deg);                      //final angle


    //main container cylinder
    G4Tubs *MainCell = new G4Tubs("MainCell",
                                  fHeContainerR,                     //inner radius
                                  fHeContainerR+fContainerThickness, //outer radius
                                  fHeContainerZ/2,                     //length
                                  0.*deg,                            //start angle
                                  360.*deg);                         //final angle

    //upstream extension tube cylinder
    G4Tubs *ExtCellU = new G4Tubs("ExtensionU",
                                  fExtensionR,                     //inner radius
                                  fExtensionR+fContainerThickness, //outer radius
                                  fExtensionZU/2,                     //length
                                  0.*deg,                          //start angle
                                  360.*deg);                       //final angle

    //downstream extension tube cylinder
    G4Tubs *ExtCellD = new G4Tubs("ExtensionD",
                                  fExtensionR,                     //inner radius
                                  fExtensionR+fContainerThickness, //outer radius
                                  fExtensionZD/2,                     //length
                                  0.*deg,                          //start angle
                                  360.*deg);                       //final angle

    //main cell end wall with a hole in the middle where extension
    //cylinder connects to main cell
    G4Tubs *MainCellEnd = new G4Tubs("MainCellEnd",
                                     fExtensionR,                       //inner radius
                                     fHeContainerR+fContainerThickness, //outer radius
                                     fContainerThickness/2,               //length
                                     0.*deg,                            //start angle
                                     360.*deg);                         //final angle
    //be windows
    G4Tubs *BerylliumWindow = new G4Tubs("BerylliumWindow",
                                         0,                     //inner radius
                                         fExtensionR,              //outer radius
                                         fBeThickness/2,             //length
                                         0.*deg,                   //start angle
                                         360.*deg);                //final angle

    //gas that fills the main cell outside teflon
    G4Tubs *HeMainCell = new G4Tubs("HeMainCell",
                                    0,                              //inner radius
                                    fHeContainerR,                  //outer radius
                                    fHeContainerZ/2,                //length
                                    0.*deg,                         //start angle
                                    360.*deg);                      //final angle

    //he inside extension tubes
    G4Tubs *HeExtCellU = new G4Tubs("HeExtCellU",
                                    0,                        //inner radius
                                    fExtensionR,              //outer radius
                                    (fExtensionZU+fContainerThickness)/2,     //length
                                    0.*deg,                   //start angle
                                    360.*deg);            //final angle
    G4Tubs *HeExtCellD = new G4Tubs("HeExtCellD",
                                    0,                        //inner radius
                                    fExtensionR,              //outer radius
                                    (fExtensionZD+fContainerThickness)/2,     //length
                                    0.*deg,                   //start angle
                                    360.*deg);            //final angle

    //create union solid for the He gas to create 1 logical volume
    //that can be set to sensitive detector
    G4UnionSolid* HeUnion1 = new G4UnionSolid
            ("HeUnion1",
             HeMainCell,
             HeExtCellU,
             0,  //no rotation
             G4ThreeVector(0, 0, -(fHeContainerZ+fExtensionZU+fContainerThickness)/2) ); //move extension cell

    G4UnionSolid* HeOutsideTeflon = new G4UnionSolid
            ("HeOutsideTeflon",
             HeUnion1,
             HeExtCellD,
             0,  //no rotation
             G4ThreeVector(0, 0, +(fHeContainerZ+fExtensionZD+fContainerThickness)/2) ); //move extension cell


    //------------------------------------------------------------------------------
    //Define logical volumes
    //------------------------------------------------------------------------------

    //invoke fVesselLogic
    fVesselLogic = new G4LogicalVolume
            (VesselCyl,   //solid
             fNistManager->FindOrBuildMaterial("G4_AIR"), //material
             "LfVesselLogic");


    G4LogicalVolume *LMainCell = new G4LogicalVolume
            (MainCell,   //solid
             fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"), //material
             "LogicMainCell");

    G4LogicalVolume *LExtCellU = new G4LogicalVolume
            (ExtCellU,   //solid
             fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"), //material
             "LogicExtCellU");

    G4LogicalVolume *LExtCellD = new G4LogicalVolume
            (ExtCellD,   //solid
             fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"), //material
             "LogicExtCellD");

    G4LogicalVolume *LMainCellEnd = new G4LogicalVolume
            (MainCellEnd,   //solid
             fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"), //material
             "LogicMainCellEnd");

    G4LogicalVolume *LBerylliumWindow = new G4LogicalVolume
            (BerylliumWindow,   //solid
             fNistManager->FindOrBuildMaterial("G4_Be"), //material
             "LogicBerylliumWindow");

    fHeOutsideTeflonLogic = new G4LogicalVolume
            (HeOutsideTeflon,
             fNistManager->FindOrBuildMaterial("ATGasMix"),
             "LogicHeOutsideTeflon"
             );

    //------------------------------------------------------------------------------
    //Set visual attributes
    //------------------------------------------------------------------------------

    G4VisAttributes* lblue  = new G4VisAttributes( G4Colour(0.0,0.0,0.75) );
    G4VisAttributes* grey   = new G4VisAttributes( G4Colour(0.5,0.5,0.5)  );
    G4VisAttributes* cyan   = new G4VisAttributes( G4Colour(0.0,1.0,1.0)  );

    fVesselLogic->SetVisAttributes(G4VisAttributes::Invisible);
    LMainCell->SetVisAttributes(grey);
    //LMainCell->SetVisAttributes(G4VisAttributes::Invisible);
    LExtCellU->SetVisAttributes(grey);
    //LExtCellU->SetVisAttributes(G4VisAttributes::Invisible);
    LExtCellD->SetVisAttributes(grey);
    //LExtCellD->SetVisAttributes(G4VisAttributes::Invisible);
    LMainCellEnd->SetVisAttributes(grey);
    //LMainCellEnd->SetVisAttributes(G4VisAttributes::Invisible);
    LBerylliumWindow->SetVisAttributes(lblue);
    //LBerylliumWindow->SetVisAttributes(G4VisAttributes::Invisible);
    fHeOutsideTeflonLogic->SetVisAttributes(cyan);
    //fHeOutsideTeflonLogic->SetVisAttributes(G4VisAttributes::Invisible);

    //------------------------------------------------------------------------------
    //Create placements of the logical volumes
    //------------------------------------------------------------------------------

    //place main vessel
    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0,0,0), //main vessel in the centre
                       LMainCell,            //main cell logic
                       "PMainCell",          //name
                       fVesselLogic,         //mother logic (vessel)
                       false,                //pMany = false, true is not implemented
                       11,fIsOverlapVol);                  //copy number


    //place main cell ends
    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , ( fHeContainerZ/2 + fContainerThickness/2 ) ),
                       LMainCellEnd,            //main cell end piece
                       "PMainCellEnd1",          //name
                       fVesselLogic,         //mother logic (vessel)
                       false,                //pMany = false, true is not implemented
                       12,fIsOverlapVol);                  //copy number

    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , -( fHeContainerZ/2 + fContainerThickness/2 ) ),
                       LMainCellEnd,            //main cell end piece
                       "PMainCellEnd2",          //name
                       fVesselLogic,         //mother logic (vessel)
                       false,                //pMany = false, true is not implemented
                       13,fIsOverlapVol);                  //copy number


    //place extension cells
    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , -( fHeContainerZ/2 + fContainerThickness + fExtensionZU/2 ) ),
                       LExtCellU,            //extension cell logic
                       "PExtCellU",          //name
                       fVesselLogic,         //mother logic (vessel)
                       false,                //pMany = false, true is not implemented
                       14,fIsOverlapVol);                  //copy number

    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , +( fHeContainerZ/2 + fContainerThickness + fExtensionZD/2 ) ),
                       LExtCellD,            //extension cell logic
                       "PExtCellD",          //name
                       fVesselLogic,         //mother logic (vessel)
                       false,                //pMany = false, true is not implemented
                       15,fIsOverlapVol);                  //copy number


    //place Be windows
    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , -( fHeContainerZ/2 + fContainerThickness + fExtensionZU + fBeThickness/2 ) ),
                       LBerylliumWindow,            //beryllium window logic
                       "PBeWindowU",          //name
                       fVesselLogic,         //mother logic (vessel)
                       false,                //pMany = false, true is not implemented
                       16,fIsOverlapVol);                  //copy number

    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , +( fHeContainerZ/2 + fContainerThickness + fExtensionZD + fBeThickness/2 ) ),
                       LBerylliumWindow,            //beryllium window logic
                       "PBeWindowD",          //name
                       fVesselLogic,         //mother logic (vessel)
                       false,                //pMany = false, true is not implemented
                       17,fIsOverlapVol);                  //copy number

} //end of MakeVessel()

//-----------------------------------------------------------------------------------
// WLS plates to fit inside vessel
// Use plastic scintillator material as its the same PVT base
// Each plate is trapesoid.
// fNwls phi segmentation...default 8-fold
// Looking down beam axis they form a fNwls-sided polygon. A small gap is maintained
// between each plate so that light will not readily pass from one plate to the next.
// WLS plates are positioned in PlaceParts()
// JRMA 29/06/2020
//-----------------------------------------------------------------------------------
void A2ActiveHe3::MakeWLS()
{
    G4double th = 360.0*deg/fNwls;
    G4double r0 = fHeContainerR - fRadClr;
    fRwls1 = r0*cos(th/2) - 0.5*fWLSthick;
    G4double z0 = fHeContainerZ/2 - fLatClr;
    G4double d0 = r0*sin(th/2) - 0.001*mm; // slightly smaller so plates dont touch
    G4double d1 = d0 - fWLSthick*tan(th/2);
    G4Trd* wls = new G4Trd("WLS-bar", z0,z0,d0,d1,fWLSthick/2);
    fWLSLogic = new G4LogicalVolume
            (wls,fNistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"),
             "LogicWLS");
    G4VisAttributes* blue   = new G4VisAttributes( G4Colour(0.0,0.0,1.0)  );
    fWLSLogic->SetVisAttributes(blue);
}


/**********************************************************

This function builds the printed circuit board and places
it to member fPCBLogic

**********************************************************/
void A2ActiveHe3::MakePCB() {

    //------------------------------------------------------------------------------
    //Define parameter values
    //------------------------------------------------------------------------------

    fPCBThickness = 1.*mm;
    fPCBRadius = fHeContainerR-4.5*mm; //define as a function of container R
    fTeflonThicknessEnd = 1.*mm;
    fPCBZ = fHeContainerZ - 20.*mm + 2*fTeflonThicknessEnd;
    //------------------------------------------------------------------------------
    //Define shapes (solids)
    //------------------------------------------------------------------------------

    //pcb is one cylinder
    G4Tubs *PCB = new G4Tubs("PCB",
                             //			   fPCBRadius,                     //inner radius
                             0,                     //inner radius
                             fPCBRadius+fPCBThickness,       //outer radius
                             fPCBZ/2,                  //length
                             0.*deg,                         //start angle
                             360.*deg);                      //final angle

    //------------------------------------------------------------------------------
    //Define logical volumes
    //------------------------------------------------------------------------------

    //invoke fPCBLogic
    fPCBLogic = new G4LogicalVolume
            (PCB,   //solid
             fNistManager->FindOrBuildMaterial("G4_C"), //material
             "LfPCBLogic");

    //------------------------------------------------------------------------------
    //Set visual attributes
    //------------------------------------------------------------------------------

    G4VisAttributes* lgreen = new G4VisAttributes( G4Colour(0.0,0.75,0.0) );

    fPCBLogic->SetVisAttributes(lgreen);

    //fPCBLogic is placed inside fVesselLogic in function PlaceParts()
} //end of MakePCB()


/**********************************************************

This function builds the teflan layer and Al-mylar windows
and places them to member fTeflanLogic

**********************************************************/
void A2ActiveHe3::MakeTeflonLayer() {

    //------------------------------------------------------------------------------
    //Define parameter values
    //------------------------------------------------------------------------------
    //fTeflonThicknessEnd = 1.*mm; //defined in MakePCB() to calculate pcb length
    fTeflonThicknessCyl = 0.5*mm;
    fTeflonR = fPCBRadius - fTeflonThicknessCyl; //reduced from original 95 to 70 mm, see doc
    fTeflonZ = fPCBZ - 2*fTeflonThicknessEnd; //such length that teflonz+2endcaps=pcbz
    fMylarThickness = 0.005*mm;

    /*parameters for the flat parts on teflon cyl where the pmt's go
    the minimum flattening height depending on R was calculated as
    h(fTeflonR) = fTeflonR - sqrt(fTeflonR^2 - (fPMTx/2)^2). The used
    number is more just in case. For h(R=95) = 0.047mm; h(R=46) = 0.098mm*/
    fFlatteningWidth = 6*mm;
    fFlatteningHeight = 0.12*mm;

    //------------------------------------------------------------------------------
    //Define shapes (solids)
    //------------------------------------------------------------------------------

    G4Tubs *TeflonLogicCyl = new G4Tubs
            ("TeflonLogicCyl",
             0,                              //inner radius
             fTeflonR+fTeflonThicknessCyl,   //outer radius
             fPCBZ/2,                        //same length as pcb
             0.*deg,                         //start angle
             360.*deg);                      //final angle


    //teflan cylinder covering pcb
    G4Tubs *TeflonCyl = new G4Tubs("TeflonCyl",
                                   fTeflonR,                       //inner radius
                                   fTeflonR+fTeflonThicknessCyl,   //outer radius
                                   fTeflonZ/2,                     //length
                                   0.*deg,                         //start angle
                                   360.*deg);                      //final angle

    G4Tubs *TeflonCylEnd = new G4Tubs("TeflonCylEnd",
                                      fExtensionR,                  //inner radius
                                      fTeflonR+fTeflonThicknessCyl, //outer radius
                                      fTeflonThicknessEnd/2,          //length
                                      0.*deg,                       //start angle
                                      360.*deg);                    //final angle

    G4Tubs *MylarWindow = new G4Tubs("MylarWindow",
                                     0,                       //inner radius
                                     fExtensionR,             //outer radius
                                     fMylarThickness/2,          //length
                                     0.*deg,                  //start angle
                                     360.*deg);               //final angle

    G4Box *Flattening = new G4Box("Flattening", fFlatteningWidth/2, fFlatteningHeight/2, fTeflonZ/2);

    //create rotations
    G4double theta = 90*deg;
    G4ThreeVector u = G4ThreeVector(std::cos(theta), -std::sin(theta), 0);
    G4ThreeVector v = G4ThreeVector(std::sin(theta), std::cos(theta), 0);
    G4ThreeVector w = G4ThreeVector(0, 0, 1);
    G4RotationMatrix *rotz90  = new G4RotationMatrix(u, v, w);

    theta = 45*deg;
    u = G4ThreeVector(std::cos(theta), -std::sin(theta), 0);
    v = G4ThreeVector(std::sin(theta), std::cos(theta), 0);
    w = G4ThreeVector(0, 0, 1);
    G4RotationMatrix *rotz45  = new G4RotationMatrix(u, v, w);

    theta = 135*deg;
    u = G4ThreeVector(std::cos(theta), -std::sin(theta), 0);
    v = G4ThreeVector(std::sin(theta), std::cos(theta), 0);
    w = G4ThreeVector(0, 0, 1);
    G4RotationMatrix *rotz135  = new G4RotationMatrix(u, v, w);

    G4UnionSolid* TefU1 = new G4UnionSolid
            ("TefU1",
             TeflonCyl,
             Flattening,
             rotz90,
             G4ThreeVector(fTeflonR-fFlatteningHeight/2, 0, 0) ); //move the flattening part

    G4UnionSolid* TefU2 = new G4UnionSolid
            ("TefU2",
             TefU1,
             Flattening,
             rotz90,
             G4ThreeVector(-fTeflonR+fFlatteningHeight/2, 0, 0) ); //move the flattening part

    G4UnionSolid* TefU3 = new G4UnionSolid
            ("TefU3",
             TefU2,
             Flattening,
             0,
             G4ThreeVector(0, fTeflonR-fFlatteningHeight/2, 0) ); //move the flattening part

    G4UnionSolid* TefU4 = new G4UnionSolid
            ("TefU4",
             TefU3,
             Flattening,
             0,
             G4ThreeVector(0, -fTeflonR+fFlatteningHeight/2, 0) ); //move the flattening part


    G4UnionSolid* TefU5 = new G4UnionSolid
            ("TefU5",
             TefU4,
             Flattening,
             rotz45,
             G4ThreeVector( -(fTeflonR-fFlatteningHeight/2)/sqrt(2), (fTeflonR-fFlatteningHeight/2)/sqrt(2), 0) );

    G4UnionSolid* TefU6 = new G4UnionSolid
            ("TefU6",
             TefU5,
             Flattening,
             rotz45,
             G4ThreeVector( (fTeflonR-fFlatteningHeight/2)/sqrt(2), -(fTeflonR-fFlatteningHeight/2)/sqrt(2), 0) );

    G4UnionSolid* TefU7 = new G4UnionSolid
            ("TefU7",
             TefU6,
             Flattening,
             rotz135,
             G4ThreeVector( -(fTeflonR-fFlatteningHeight/2)/sqrt(2), -(fTeflonR-fFlatteningHeight/2)/sqrt(2), 0) );

    G4UnionSolid* TefU8 = new G4UnionSolid
            ("TefU8",
             TefU7,
             Flattening,
             rotz135,
             G4ThreeVector( (fTeflonR-fFlatteningHeight/2)/sqrt(2), (fTeflonR-fFlatteningHeight/2)/sqrt(2), 0) );

    
    //he inside the teflon cylinder
    G4Tubs *HeTube = new G4Tubs("HeTube",
                                0,                       //inner radius
                                fTeflonR+fTeflonThicknessCyl,   //outer radius
                                fTeflonZ/2,               // same length as teflon main cyl
                                0.*deg,                         //start angle
                                360.*deg);                      //final angle

    //he gas = teflon main cyl (outer r) - teflon main cyl solid
    G4SubtractionSolid* HeInsideTeflon = new G4SubtractionSolid
            ("HeInsideTeflon",
             HeTube,
             TefU8 );

    //------------------------------------------------------------------------------
    //Define logical volumes
    //------------------------------------------------------------------------------

    //invoke fTeflonLogic
    fTeflonLogic = new G4LogicalVolume
            (TeflonLogicCyl,   //solid
             //fNistManager->FindOrBuildMaterial("G4_AIR"), //material
             fNistManager->FindOrBuildMaterial("ATGasMix"),
             "LfTeflonLogic");

    fTeflonCylLogic = new G4LogicalVolume
            (TefU8,   //solid
             fNistManager->FindOrBuildMaterial("ATTeflon"), //material
             "LTeflonCyl");

    fTeflonCylEndLogic = new G4LogicalVolume
            (TeflonCylEnd,   //solid
             fNistManager->FindOrBuildMaterial("ATTeflon"), //material
             "LTeflonCylEnd");

    fMylarLogic = new G4LogicalVolume
            (MylarWindow,   //solid
             fNistManager->FindOrBuildMaterial("ATMylarW"), //material
             "LMylarWindow");

    fHeInsideTeflonLogic = new G4LogicalVolume
            (HeInsideTeflon,   //solid
             fNistManager->FindOrBuildMaterial("ATGasMix"), //material
             "LHeInsideTeflon");

    //------------------------------------------------------------------------------
    //Set visual attributes
    //------------------------------------------------------------------------------

    G4VisAttributes* red    = new G4VisAttributes( G4Colour(0.3,0.0,0.0)  );
    G4VisAttributes* blue   = new G4VisAttributes( G4Colour(0.0,0.0,1.0)  );
    G4VisAttributes* cyan   = new G4VisAttributes( G4Colour(0.0,1.0,1.0)  );

    fTeflonLogic->SetVisAttributes(G4VisAttributes::Invisible);
    fTeflonCylLogic ->SetVisAttributes(red);
    fTeflonCylEndLogic->SetVisAttributes(red);
    fMylarLogic->SetVisAttributes(blue);
    fHeInsideTeflonLogic->SetVisAttributes(cyan);

    //------------------------------------------------------------------------------
    //Create placements of the logical volumes
    //------------------------------------------------------------------------------

    //place the teflon cylinder
    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0,0,0), //teflon cylinder in the centre
                       fTeflonCylLogic,                 //teflon cylinder logic
                       "PTeflonCyl",               //name
                       fTeflonLogic,            //mother logic (teflonLogic)
                       false,                //pMany = false, true is not implemented
                       31,fIsOverlapVol);                  //copy number

    //place cylinder end walls
    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , ( fTeflonZ/2 + fTeflonThicknessEnd/2 ) ),
                       fTeflonCylEndLogic,            //main cell end piece
                       "PTeflonCylEnd1",          //name
                       fTeflonLogic,         //mother logic (teflonLogic)
                       false,                //pMany = false, true is not implemented
                       32,fIsOverlapVol);                  //copy number

    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , -( fTeflonZ/2 + fTeflonThicknessEnd/2 ) ),
                       fTeflonCylEndLogic,            //main cell end piece
                       "PTeflonCylEnd2",          //name
                       fTeflonLogic,         //mother logic (teflonLogic)
                       false,                //pMany = false, true is not implemented
                       33,fIsOverlapVol);                  //copy number

    //place the mylar windows
    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , ( fTeflonZ/2 + fMylarThickness/2 ) ),
                       fMylarLogic,            //mylar window
                       "PMylarWindow1",          //name
                       fTeflonLogic,         //mother logic (teflonLogic)
                       false,                //pMany = false, true is not implemented
                       34,fIsOverlapVol);                  //copy number

    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0 , -( fTeflonZ/2 + fMylarThickness/2 ) ),
                       fMylarLogic,            //mylar window
                       "PMylarWindow2",          //name
                       fTeflonLogic,         //mother logic (teflonLogic)
                       false,                //pMany = false, true is not implemented
                       35,fIsOverlapVol);                  //copy number

} //end of MakeTeflonLayer()




/**********************************************************

This function builds mylar windows to section the main cell
The flattening for pmt's made it a bit tricky. Simply copied
over the logic that was used to create the teflon cylider,
flattening and helium inside it and scaled down the z dimension
to mylar thickness

**********************************************************/
void A2ActiveHe3::MakeMylarSections(G4int nrofsections) {

    //---------------------------------------------------------------------------
    //Define parameters
    //---------------------------------------------------------------------------
    fMylarSecZ = 0.005*mm;

    if (nrofsections == 8) {
        G4cout << "A2ActiveHe3::MakeMylarSections() Sectioning active target to 8." << G4endl;
        fMylarWinStep = 50.*mm;          //values with 8 sections
        fMylarWinOffset = 50.*mm;
        fNMylarSecWins = 7;
    }
    else if (nrofsections == 4) {
        G4cout << "A2ActiveHe3::MakeMylarSections() Sectioning active target to 4." << G4endl;
        fMylarWinStep = 100.*mm;         //values with 4 sections
        fMylarWinOffset = 100.*mm;
        fNMylarSecWins = 3;
    }
    else {
        G4cout << "A2ActiveHe3::MakeMylarSections() Unsupported number of active target sections defined, the target will not be sectioned." << G4endl;
        return;
    }

    //------------------------------------------------------------------------------
    //Create solid
    //------------------------------------------------------------------------------
    G4Tubs *MylarSecOuter = new G4Tubs("MylarSecOuter",
                                       fTeflonR,                       //inner radius
                                       fTeflonR+fTeflonThicknessCyl,   //outer radius
                                       fMylarSecZ/2,                     //length
                                       0.*deg,                         //start angle
                                       360.*deg);                      //final angle

    G4Box *MylarSecFlattening = new G4Box("MylarSecFlattening", fFlatteningWidth/2, fFlatteningHeight/2, fMylarSecZ/2);

    //create rotations
    G4double theta = 90*deg;
    G4ThreeVector u = G4ThreeVector(std::cos(theta), -std::sin(theta), 0);
    G4ThreeVector v = G4ThreeVector(std::sin(theta), std::cos(theta), 0);
    G4ThreeVector w = G4ThreeVector(0, 0, 1);
    G4RotationMatrix *rotz90  = new G4RotationMatrix(u, v, w);

    theta = 45*deg;
    u = G4ThreeVector(std::cos(theta), -std::sin(theta), 0);
    v = G4ThreeVector(std::sin(theta), std::cos(theta), 0);
    w = G4ThreeVector(0, 0, 1);
    G4RotationMatrix *rotz45  = new G4RotationMatrix(u, v, w);

    theta = 135*deg;
    u = G4ThreeVector(std::cos(theta), -std::sin(theta), 0);
    v = G4ThreeVector(std::sin(theta), std::cos(theta), 0);
    w = G4ThreeVector(0, 0, 1);
    G4RotationMatrix *rotz135  = new G4RotationMatrix(u, v, w);

    G4UnionSolid* MylU1 = new G4UnionSolid
            ("MylU1",
             MylarSecOuter,
             MylarSecFlattening,
             rotz90,
             G4ThreeVector(fTeflonR-fFlatteningHeight/2, 0, 0) ); //move the flattening part

    G4UnionSolid* MylU2 = new G4UnionSolid
            ("MylU2",
             MylU1,
             MylarSecFlattening,
             rotz90,
             G4ThreeVector(-fTeflonR+fFlatteningHeight/2, 0, 0) ); //move the flattening part

    G4UnionSolid* MylU3 = new G4UnionSolid
            ("MylU3",
             MylU2,
             MylarSecFlattening,
             0,
             G4ThreeVector(0, fTeflonR-fFlatteningHeight/2, 0) ); //move the flattening part

    G4UnionSolid* MylU4 = new G4UnionSolid
            ("MylU4",
             MylU3,
             MylarSecFlattening,
             0,
             G4ThreeVector(0, -fTeflonR+fFlatteningHeight/2, 0) ); //move the flattening part


    G4UnionSolid* MylU5 = new G4UnionSolid
            ("MylU5",
             MylU4,
             MylarSecFlattening,
             rotz45,
             G4ThreeVector( -(fTeflonR-fFlatteningHeight/2)/sqrt(2), (fTeflonR-fFlatteningHeight/2)/sqrt(2), 0) );

    G4UnionSolid* MylU6 = new G4UnionSolid
            ("MylU6",
             MylU5,
             MylarSecFlattening,
             rotz45,
             G4ThreeVector( (fTeflonR-fFlatteningHeight/2)/sqrt(2), -(fTeflonR-fFlatteningHeight/2)/sqrt(2), 0) );

    G4UnionSolid* MylU7 = new G4UnionSolid
            ("MylU7",
             MylU6,
             MylarSecFlattening,
             rotz135,
             G4ThreeVector( -(fTeflonR-fFlatteningHeight/2)/sqrt(2), -(fTeflonR-fFlatteningHeight/2)/sqrt(2), 0) );

    G4UnionSolid* MylU8 = new G4UnionSolid
            ("MylU8",
             MylU7,
             MylarSecFlattening,
             rotz135,
             G4ThreeVector( (fTeflonR-fFlatteningHeight/2)/sqrt(2), (fTeflonR-fFlatteningHeight/2)/sqrt(2), 0) );

    
    //mylar inside the teflon cylinder
    G4Tubs *MylarSec = new G4Tubs("MylarSec",
                                  0,                       //inner radius
                                  fTeflonR+fTeflonThicknessCyl,   //outer radius
                                  fMylarSecZ/2,
                                  0.*deg,                         //start angle
                                  360.*deg);                      //final angle

    //create the actual solid
    G4SubtractionSolid* MylarInsideTeflon = new G4SubtractionSolid
            ("MylarInsideTeflon",
             MylarSec,
             MylU8 );

    //------------------------------------------------------------------------------
    //Create logical volume
    //------------------------------------------------------------------------------
    fMylarSectionLogic = new G4LogicalVolume
            (MylarInsideTeflon,
             fNistManager->FindOrBuildMaterial("ATMylarW"),
             "LMylarSection");

    //------------------------------------------------------------------------------
    //Create placements of the logical volumes
    //------------------------------------------------------------------------------
    G4int nwindows;
    for (nwindows = 0; nwindows < fNMylarSecWins; nwindows++) {

        //place the mylar sectioning windows along z
        new G4PVPlacement
                (0,
                 G4ThreeVector(0, 0, -fTeflonZ/2 + fMylarWinOffset + nwindows*fMylarWinStep),
                 fMylarSectionLogic,
                 "PMylarSection",                  //name
                 fHeInsideTeflonLogic,            //mother logic
                 false,                   //pMany = false, true is not implemented
                 200+nwindows,fIsOverlapVol);                  //copy number

    } //end of for loop for mylar windows placement

} //end MakeMylarSections()




/**********************************************************

This function places the separate logical parts of the
detector inside fMyLogic

**********************************************************/
void A2ActiveHe3::PlaceParts() {

    //------------------------------------------------------------------------------
    //Define shape for fMyLogic and construct it (simple cylinder that
    //has room for everything)
    //------------------------------------------------------------------------------

    G4Tubs *MyLogicCyl = new G4Tubs
            ("MyLogicCyl",
             0,                              //inner radius
             fHeContainerR+fContainerThickness,   //outer radius
             fHeContainerZ/2 + fContainerThickness + fExtensionZU + fBeThickness,
             0.*deg,                         //start angle
             360.*deg);                      //final angle

    fMyLogic = new G4LogicalVolume
            (MyLogicCyl,                                  //solid
             fNistManager->FindOrBuildMaterial("G4_AIR"), //material
             "LfMyLogic");                                //name

    fMyLogic->SetVisAttributes(G4VisAttributes::Invisible);

    //------------------------------------------------------------------------------
    //Place parts of the detector
    //------------------------------------------------------------------------------
    if(!fIsWLS){
        //place the helium inside teflon cylinder
        new G4PVPlacement (0,  //no rotation
                           G4ThreeVector(0,0,0), //teflon cylinder in the centre
                           fHeInsideTeflonLogic,     //helium logic
                           "PHeInsideTeflon",        //name
                           fTeflonLogic,            //mother logic (teflonLogic)
                           false,                //pMany = false, true is not implemented
                           45,fIsOverlapVol);                  //copy number

        //place the teflon layer insice pcb
        new G4PVPlacement (0,                    //no rotation
                           G4ThreeVector(0,0,0), //center
                           fTeflonLogic,         //main cell logic
                           "PTeflonLogic",       //name
                           fPCBLogic,             //mother logic
                           false,                //pMany = false, true is not implemented
                           44,fIsOverlapVol);                  //copy number

        //place the pcb inside the outside helium
        new G4PVPlacement (0,                    //no rotation
                           G4ThreeVector(0,0,0), //center
                           fPCBLogic,            //main cell logic
                           "PPCBLogic",          //name
                           fHeOutsideTeflonLogic,//mother logic
                           false,                //pMany = false, true is not implemented
                           43,fIsOverlapVol);    //copy number
    }
    //place the outside helium inside vessel
    new G4PVPlacement (0,  //no rotation
                       G4ThreeVector(0, 0, 0),
                       fHeOutsideTeflonLogic,            //he gas logic
                       "PHeOutsideTeflon",          //name
                       fVesselLogic,         //mother logic (vessel)
                       false,                //pMany = false, true is not implemented
                       42,fIsOverlapVol);                  //copy number
    //
    // place WLS plates in He gas volume
    // loop round phi segments (sides of polygon)
    // NB copy numbers provisionally set 140 onwards...this may have to change
    //
    G4double th = 360*deg/fNwls;
    G4double th2 = 0.0;
    for(G4int iw=0; iw<fNwls; iw++){
        G4RotationMatrix* pp = new G4RotationMatrix();
        G4double th1 = iw*th;
        pp->rotateY(90*deg);            // align parallel beam axis
        pp->rotateX(th1);               // align for side of polygon
        G4double xw = fRwls1*cos(th2);  // coordinates for polygon side
        G4double yw = fRwls1*sin(th2);
        char nw[16];
        sprintf(nw,"PWLS_%d",iw);
        new G4PVPlacement (pp,
                           G4ThreeVector(xw, yw, 0),
                           fWLSLogic,
                           nw,
                           fHeOutsideTeflonLogic,//mother logic (He)
                           false,                //pMany = false
                           iw,fIsOverlapVol);    //copy number
        th2 -= th;
    }

    //place the vessel inside MyLogic
    new G4PVPlacement (0,                    //no rotation
                       G4ThreeVector(0,0,0), //center
                       fVesselLogic,         //main cell logic
                       "PVesselLogic",       //name
                       fMyLogic,             //mother logic
                       false,                //pMany = false, true is not implemented
                       41,fIsOverlapVol);    //copy number

}


/**********************************************************

This function sets the optical properties
for some materials and surfaces

**********************************************************/
void A2ActiveHe3::SetOpticalProperties() {

    //------------------------------------------------------------------------------
    //gas mixture
    //------------------------------------------------------------------------------
    const G4int nentries = 18;

    //Photon energies 1.88 - 3.87 eV --> 660 - 320 nm, even step of 20 nm
    G4double HeN_Energy[nentries] = {1.88*eV, 1.94*eV, 2.0*eV, 2.07*eV, 2.14*eV, 2.21*eV, 2.3*eV,
                                     2.38*eV, 2.48*eV, 2.58*eV, 2.7*eV, 2.82*eV, 2.95*eV, 3.1*eV,
                                     3.26*eV, 3.44*eV, 3.65*eV, 3.87*eV};

    //this gives a peak from 400 to 440 nm --> 3.1 to 2.95 eV
    G4double HeN_SCINT[nentries] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                    0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1};

    //See ref in Seian's paper; no need for higher accuracy in my opinion
    G4double HeN_RIND[nentries] =
    {1.00003, 1.00003, 1.00003, 1.00003, 1.00003, 1.00003, 1.00003, 1.00003, 1.00003,
     1.00003, 1.00003, 1.00003, 1.00003, 1.00003, 1.00003, 1.00003, 1.00003, 1.00003};

    //absorption length
    G4double HeN_ABSL[nentries] =
    {3500.*cm,3500.*cm,3500.*cm,3500.*cm,3500.*cm,3500.*cm,3500.*cm,3500.*cm, 3500.*cm,
     3500.*cm,3500.*cm,3500.*cm,3500.*cm,3500.*cm,3500.*cm,3500.*cm,3500.*cm, 3500.*cm};

    G4MaterialPropertiesTable* HeN_mt = new G4MaterialPropertiesTable();
    HeN_mt->AddProperty("FASTCOMPONENT", HeN_Energy, HeN_SCINT, nentries);
    HeN_mt->AddProperty("SLOWCOMPONENT", HeN_Energy, HeN_SCINT, nentries);
    HeN_mt->AddProperty("RINDEX",        HeN_Energy, HeN_RIND,  nentries);
    HeN_mt->AddProperty("ABSLENGTH",     HeN_Energy, HeN_ABSL,  nentries);

    HeN_mt->AddConstProperty("SCINTILLATIONYIELD", fScintYield/MeV   );
    HeN_mt->AddConstProperty("RESOLUTIONSCALE" ,   1.0        );

    //numbers from Jonh's optical simulation--------------
    // HeN_mt->AddConstProperty("FASTTIMECONSTANT",20.*ns);
    // HeN_mt->AddConstProperty("SLOWTIMECONSTANT",45.*ns);
    // HeN_mt->AddConstProperty("YIELDRATIO",1.0);

    //numbers from Active target simulation with gdml file
    // HeN_mt->AddConstProperty("FASTTIMECONSTANT",   0.5*ns     );// 8  ns  20
    // HeN_mt->AddConstProperty("SLOWTIMECONSTANT",   23.*ns     );// 23 ns  45
    // HeN_mt->AddConstProperty("YIELDRATIO",         0.5        );

    //numbers according to Seian's thesis-----------------
    HeN_mt->AddConstProperty("FASTTIMECONSTANT",   30.*ns     );
    HeN_mt->AddConstProperty("SLOWTIMECONSTANT",   60.*ns     );
    HeN_mt->AddConstProperty("YIELDRATIO",         1.0        );
    //----------------------------------------------------

    //add these properties to GasMix
    fNistManager->FindOrBuildMaterial("ATGasMix")->SetMaterialPropertiesTable(HeN_mt);

    //------------------------------------------------------------------------------
    //teflon properties
    //------------------------------------------------------------------------------

    //PTFE reflectivity from Janecek
    //"REFLECTIVITY SPECTRA FOR COMMONLY USED REFLECTORS"
    G4double PTFE_Reflectivity[nentries] =
    {0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,
     0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95};

    G4MaterialPropertiesTable* Teflon_mt = new G4MaterialPropertiesTable();
    Teflon_mt->AddProperty("REFLECTIVITY",HeN_Energy, PTFE_Reflectivity ,nentries);

    fNistManager->FindOrBuildMaterial("ATTeflon")->SetMaterialPropertiesTable(Teflon_mt);

    //------------------------------------------------------------------------------
    //mylar
    //------------------------------------------------------------------------------

    //flat 95% according to Seian's thesis
    G4double mylar_REFL[nentries] =
    {0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,
     0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95};

    //These seem to be the same as for quartz. In terms of the simulation
    //these numbers are not very relevant. Interpolated and extrapolated
    //results of Seian to match my energies and range
    // G4double mylar_RIND[nentries] =
    //   {1.51551, 1.5173, 1.51899, 1.52084, 1.52256, 1.52418, 1.52611,
    //    1.52752, 1.52931, 1.53097, 1.53321, 1.53548, 1.53828, 1.54184,
    //    1.54528, 1.55176, 1.55819, 1.56375};

    //http://www.filmetrics.com/refractive-index-database/PET/Estar-Melinex-Mylar
    G4double mylar_RIND[nentries] =
    {1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68,
     1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68};

    G4MaterialPropertiesTable* Mylar_mt = new G4MaterialPropertiesTable();
    Mylar_mt->AddProperty("REFLECTIVITY",HeN_Energy, mylar_REFL ,nentries);
    Mylar_mt->AddProperty("RINDEX",HeN_Energy,mylar_RIND, nentries);

    fNistManager->FindOrBuildMaterial("ATMylarW")->SetMaterialPropertiesTable(Mylar_mt);

    //------------------------------------------------------------------------------
    //Epoxy (best guesses I could do)
    //------------------------------------------------------------------------------

    //the film refractive index at 589 nm from Martin Sharratt's e-mail
    //no dependency known
    G4double Epoxy_RIND[nentries] =
    {1.54, 1.54, 1.54, 1.54, 1.54, 1.54, 1.54, 1.54, 1.54,
     1.54, 1.54, 1.54, 1.54, 1.54, 1.54, 1.54, 1.54, 1.54};

    //use 1, the pmt acceptance spectrum already includes the epoxy effect according
    //to e-mail from Sharratt
    G4double Epoxy_Transmittance[nentries] =
    {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    G4MaterialPropertiesTable *Epoxy_mt = new G4MaterialPropertiesTable();
    Epoxy_mt->AddProperty("RINDEX",HeN_Energy,Epoxy_RIND, nentries);
    Epoxy_mt->AddProperty("TRANSMITTANCE", HeN_Energy, Epoxy_Transmittance, nentries);

    fNistManager->FindOrBuildMaterial("ATEpoxy")->SetMaterialPropertiesTable(Epoxy_mt);

    //------------------------------------------------------------------------------
    //Silicon photomultipliers (currently set so that every hit is detected!)
    //------------------------------------------------------------------------------

    // detect every photon that hits the pmt. This can be used to correlate which
    // wavelengths are detected more efficiently
    G4double SiPMT_EFF[nentries] =
    {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // Photocathod reflectivity, currently 0, so detect all
    G4double SiPMT_REFL[nentries]   =
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    G4MaterialPropertiesTable *SiPMT_mt = new G4MaterialPropertiesTable();
    SiPMT_mt->AddProperty("EFFICIENCY",HeN_Energy,SiPMT_EFF, nentries);
    SiPMT_mt->AddProperty("REFLECTIVITY",HeN_Energy,SiPMT_REFL, nentries);

    fNistManager->FindOrBuildMaterial("ATSiPMT")->SetMaterialPropertiesTable(SiPMT_mt);

    //------------------------------------------------------------------------------
    //Create optical surfaces
    //------------------------------------------------------------------------------

    //teflon - only lambertian reflection, relatively OK simple approximation
    G4OpticalSurface* OptTeflonSurface =
            new G4OpticalSurface("OTeflonSurface",  unified, groundfrontpainted, dielectric_dielectric);

    //relate the materials table
    OptTeflonSurface->SetMaterialPropertiesTable(Teflon_mt);

    //set optical surfaces
    new G4LogicalSkinSurface("LSTeflonCylSurface",fTeflonCylLogic, OptTeflonSurface);
    new G4LogicalSkinSurface("LSTeflonEndSurface",fTeflonCylEndLogic, OptTeflonSurface);

    //------------------------------------------------------------------------------------

    //Mylar windows. model uses only reflection and absorption
    G4OpticalSurface* OptMylarSurface =
            new G4OpticalSurface("OMylarSurface",  unified, polishedfrontpainted, dielectric_dielectric);
    OptMylarSurface->SetMaterialPropertiesTable(Mylar_mt);

    new G4LogicalSkinSurface("LSMylarSurface", fMylarLogic, OptMylarSurface);
    if (fMakeMylarSections) {
        new G4LogicalSkinSurface("LSMylarSecSurface", fMylarSectionLogic, OptMylarSurface);
    }

    //------------------------------------------------------------------------------------

    //epoxy surface.  It should reflect, refract or absorb
    if (fMakeEpoxy) {
        G4OpticalSurface* OptEpoxySurface =
                new G4OpticalSurface("OEpoxySurface",  unified, ground, dielectric_dielectric);
        //janecek, moses 2010; ground sigmaalpha 12 deg = 0.21rad
        OptEpoxySurface->SetSigmaAlpha(0.21);
        OptEpoxySurface->SetMaterialPropertiesTable(Epoxy_mt);

        new G4LogicalSkinSurface("LSEpoxySurface", fEpoxyLogic, OptEpoxySurface);
    }

    //------------------------------------------------------------------------------------

    //pmt surface. The specified model uses ONLY reflection or absorption
    G4OpticalSurface* OptPmtSurface =
            new G4OpticalSurface("OPmtSurface",  unified, polishedfrontpainted, dielectric_dielectric);
    OptPmtSurface->SetMaterialPropertiesTable(SiPMT_mt);

    new G4LogicalSkinSurface("LSPmtSurfarce", fPMTLogic, OptPmtSurface);

} //end SetOpticalProperties()

void A2ActiveHe3::SetScintillationYield(G4double scintyield) {

    if (scintyield > 0.01) {
        fScintYield = scintyield;
        G4cout << "A2ActiveHe3::SetScintillationYield() Using user-set scintillation yield value " <<
                  fScintYield << " per MeV" << G4endl;
    }
    else {
        G4cout << "A2ActiveHe3::SetScintillationYield() Using default value for fScintYield " <<
                  fScintYield << " per MeV" << G4endl;
    }
}
void A2ActiveHe3::DefineMaterials()
{
    G4double density, fractionmass;
    G4int ncomponents;
    //G4double pressure, temperature, a, z;

    G4Material* BerylliumW = new G4Material("ATBerylliumW", density = 1.8480*g/cm3, ncomponents = 7);
    BerylliumW->AddElement(fNistManager->FindOrBuildElement(14), 0.06*perCent);      //Si
    BerylliumW->AddElement(fNistManager->FindOrBuildElement(4), 98.73*perCent);      //Be
    BerylliumW->AddElement(fNistManager->FindOrBuildElement(6), 0.15*perCent);       //C
    BerylliumW->AddElement(fNistManager->FindOrBuildElement(8), 0.75*perCent);       //O
    BerylliumW->AddElement(fNistManager->FindOrBuildElement(12), 0.08*perCent);      //Mg
    BerylliumW->AddElement(fNistManager->FindOrBuildElement(13), 0.1*perCent);       //Al
    BerylliumW->AddElement(fNistManager->FindOrBuildElement(26), 0.13*perCent);      //Fe

    G4Material* MylarW = new G4Material("ATMylarW", density = 1.4000*g/cm3, ncomponents = 3);
    MylarW->AddElement(fNistManager->FindOrBuildElement(1), 4.2*perCent);           //H
    MylarW->AddElement(fNistManager->FindOrBuildElement(6), 62.5*perCent);          //C
    MylarW->AddElement(fNistManager->FindOrBuildElement(8), 33.3*perCent);          //O

    G4Material* Teflon = new G4Material("ATTeflon", density = 2.2000*g/cm3, ncomponents = 2);
    Teflon->AddElement(fNistManager->FindOrBuildElement(6), 24*perCent);            //C
    Teflon->AddElement(fNistManager->FindOrBuildElement(9), 76*perCent);            //F


    //Gas mixture and He3 management----------------------------------------------------

    //defining He3 according to
    // http://hypernews.slac.stanford.edu/HyperNews/geant4/get/hadronprocess/731/1.html
    G4Element* ATHe3 = new G4Element("ATHe3", "ATHe3", ncomponents=1);
    ATHe3->AddIsotope((G4Isotope*)fNistManager->FindOrBuildElement(2)->GetIsotope(0), //isot. 0 = 3he
                      100.*perCent);

    //See HeGasDensity.nb. Might be an idea to include N2 effect to density
    //G4double he3density = 0.00247621*g/cm3; //20bar, calculated from ideal gas law
    G4double he3density = 0.0033*g/cm3; //20bar, calculated from ideal gas law

    G4Material* GasMix = new G4Material("ATGasMix", he3density, ncomponents = 2);
    GasMix->AddElement(ATHe3, 99.95*perCent);                                       //He3
    // GasMix->AddElement(fNistManager->FindOrBuildElement(2), 99.95*perCent);         //He4
    GasMix->AddElement(fNistManager->FindOrBuildElement(7), 0.05*perCent);           //N
    //Gas mixture and He3 end-----------------------------------------------------
    //----! IMPORTANT! Epoxy CURRENTLY TAKEN FROM A2 SIMULATION,------------------
    //NO IDEA WHETHER IT IS CORRECT OR NOT

    //Epoxy resin (C21H25Cl05) ***not certain if correct chemical formula or density
    //for colder temperature***:
    G4Material* EpoxyResin=new G4Material("EpoxyResin", density=1.15*g/cm3, ncomponents=4);
    EpoxyResin->AddElement(fNistManager->FindOrBuildElement(6), 21); //carbon
    EpoxyResin->AddElement(fNistManager->FindOrBuildElement(1), 25); //hydrogen
    EpoxyResin->AddElement(fNistManager->FindOrBuildElement(17), 1); //chlorine
    EpoxyResin->AddElement(fNistManager->FindOrBuildElement(8), 5);  //oxygen

    //Amine Hardener (C8H18N2) ***not certain if correct chemical formula or density
    //for colder temperature***:
    G4Material* Epoxy13BAC=new G4Material("Epoxy13BAC", density=0.94*g/cm3, ncomponents=3);
    Epoxy13BAC->AddElement(fNistManager->FindOrBuildElement(6), 8);
    Epoxy13BAC->AddElement(fNistManager->FindOrBuildElement(1), 18);
    Epoxy13BAC->AddElement(fNistManager->FindOrBuildElement(7), 2);

    //Epoxy adhesive with mix ratio 100:25 parts by weight resin/hardener
    // ***don't know density at cold temperature***:
    G4Material* Epoxy=new G4Material("ATEpoxy", density=1.2*g/cm3, ncomponents=2);
    Epoxy->AddMaterial(EpoxyResin, fractionmass=0.8);
    Epoxy->AddMaterial(Epoxy13BAC, fractionmass=0.2);
    //-----------------------------------------------------------------------------

    //just silicon for now
    G4Material* SiPMT = new G4Material("ATSiPMT", density = 2.329*g/cm3, ncomponents = 1);
    SiPMT->AddElement(fNistManager->FindOrBuildElement(14),100.*perCent);    //Si

}

//---------------------------------------------------------------------------
void A2ActiveHe3::ReadParameters(const char* file)
{
    //
    // Initialise detector parameters from file
    //
    char* keylist[] = { (char *)"AT-Dim:", (char *)"AT-WLS:",
                        (char *)"AT-Scint:", (char *)"Run-Mode:", NULL };
    enum { EAT_Dim, EAT_WLS, EAT_Scint, ERun_Mode, ENULL };
    G4int ikey, iread;
    G4int ierr = 0;
    char line[256];
    char delim[64];
    //char hname[32];
    //char fname[32];
    //G4double x,y,z,dx,dy,dz;
    FILE* pdata;
    if( (pdata = fopen(file,"r")) == NULL ){
        printf("Error opening detector parameter file: %s\n",file);
        return;
    }
    while( fgets(line,256,pdata) ){
        if( line[0] == '#' ) continue;       // comment #
        printf("%s\n",line);                 // record datum to log
        sscanf(line,"%s",delim);
        for(ikey=0; ikey<ENULL; ikey++)
            if(!strcmp(keylist[ikey],delim)) break;
        switch( ikey ){
        default:
            printf("Unrecognised delimiter: %s\n",delim);
            ierr++;
            break;
        case EAT_Dim:
            iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf",
                           &fHeContainerZ,&fHeContainerR,&fContainerThickness,
                           &fExtensionZU,&fExtensionZD,&fExtensionR,&fBeThickness);
            if( iread != 7 ) ierr++;
            break;
        case EAT_WLS:
            iread = sscanf(line,"%*s%d%lf%lf%lf",
                           &fNwls,&fWLSthick,&fRadClr,&fLatClr);
            if( iread != 4 ) ierr++;
            break;
        case EAT_Scint:
            iread =sscanf(line,"%*s%lf",
                          &fScintYield);
            if( iread != 1 ) ierr++;
            break;
        case ERun_Mode:
            iread = sscanf(line,"%*s%d%d%d%d",
                           &fIsOverlapVol,&fIsWLS,&fOpticalSimulation,&fMakeMylarSections);
            if( iread != 4 ) ierr++;
            break;
        }
        if( ierr ){
            printf("Fatal Error: invalid read of parameter line %s\n %s\n",
                   keylist[ikey],line);
            exit(-1);
        }
    }
}


/**********************************************************

This function builds the sensitive detector

**********************************************************/

void A2ActiveHe3::MakeSensitiveDetector(){
    if(!fAHe3SD){
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        fAHe3SD = new A2SD("AHe3SD",fNwls);
        SDman->AddNewDetector(fAHe3SD);
        fWLSLogic->SetSensitiveDetector(fAHe3SD);
        fRegionAHe3->AddRootLogicalVolume(fWLSLogic);
    }
}

