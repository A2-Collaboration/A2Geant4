// Abstract file-based event generator
// Author: Dominik Werthmueller, 2018

#include "G4ParticleDefinition.hh"

#include "A2FileGenerator.hh"

//______________________________________________________________________________
A2FileGenerator::A2FileGenerator(const char* filename)
{
    // Constructor.

    // init members
    fFileName = filename;
}

//______________________________________________________________________________
void A2FileGenerator::A2GenParticle_t::Print(const char* pre)
{
    // Print the content of this class.

    G4cout << pre << "Name                : " << (fDef ? fDef->GetParticleName() : "unknown") << G4endl
           << pre << "Momentum            : " << fP << G4endl
           << pre << "Energy              : " << fE << G4endl
           << pre << "Mass                : " << fM << G4endl
           << pre << "Vertex              : " << fX << G4endl
           << pre << "Vertex time         : " << fT << G4endl
           << pre << "Is tracked          : " << (fIsTrack ? "yes" : "no") << G4endl;
}

//______________________________________________________________________________
void A2FileGenerator::Print()
{
    // Print the content of this class.

    G4cout << "File name           : " << fFileName << G4endl
           << "Number of particles : " << fPart.size() << G4endl
           << "Primary vertex      : " << fVertex << G4endl
           << "Beam particle" << G4endl;
    fBeam.Print("  ");
    for (G4int i = 0; i < (G4int)fPart.size(); i++)
    {
        G4cout << "Particle " << i+1 << G4endl;
        fPart[i].Print("  ");
    }
}

