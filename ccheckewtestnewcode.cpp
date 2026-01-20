#include "lcio.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCEvent.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/LCTOOLS.h"
#include "EVENT/Vertex.h"
#include "IMPL/VertexImpl.h"
#include "UTIL/LCRelationNavigator.h"
// ROOT includes
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <array>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <cmath>

using namespace lcio;
using namespace std;

void ccheckewtestnewcode(const char* file) {
IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
lcReader->setReadCollectionNames({
"MCTruthRecoLink",
"RecoMCTruthLink",
"MCParticlesSkimmed",
"V0RecoParticles",
"PandoraPFOs",
"V0Vertices"
});

lcReader->open(file);
// Create canvases
TCanvas* c1 = new TCanvas("c1", "x MC vs RECO", 1200, 600);
TCanvas* c2 = new TCanvas("c2", "y MC vs RECO", 1200, 600);
TCanvas* c3 = new TCanvas("c3", "z MC vs RECO", 1200, 600);
TCanvas* c4 = new TCanvas("c4", "Differences", 1200, 600);
c1->Divide(2, 1);
c2->Divide(2, 1);
c3->Divide(2, 1);
c4->Divide(2, 1);

// Create histograms for MC positions
TH1F* histo = new TH1F("histo", "MC x position", 100, -5,5);
TH1F* histo1 = new TH1F("histo1", "MC y position", 100, -5,5);
TH1F* histo2 = new TH1F("histo2", "MC z position", 100, -5,5);
// Create histograms for V0 positions
TH1F* histo1vz = new TH1F("histo1vz", "V0 x position", 100,-700, 700);

TH1F* histo2vz = new TH1F("histo2vz", "V0 y position", 100, -700, 700);
TH1F* histo3vz = new TH1F("histo3vz", "V0 z position", 100, -700, 700);
// Create histograms for differences
TH1F* h_diff_x = new TH1F("h_diff_x", "x difference (MC-V0)", 100, -300, 300);
TH1F* h_diff_y = new TH1F("h_diff_y", "y difference (MC-V0)", 100, -300, 300);
TH1F* h_diff_z = new TH1F("h_diff_z", "z difference (MC-V0)", 100, -300, 300);
TH1F* histo_event = new TH1F("histo_event", "Events withoutV0", 100, -2000, 14800);
// Create ROOT file and tree
TFile* froot = new TFile("VZeroTree.root", "RECREATE");
TTree* t1 = new TTree("v0", "reconstructed vs true");

// Per-event variables
std::vector<double> xmc, ymc, zmc;
std::vector<double> xvzero, yvzero, zvzero;
std::vector<double> diff_x, diff_y, diff_z;

int has_v0_collection;
int event_number;
int num_matches;

// Branch for MC positions
t1->Branch("xmc", &xmc);
t1->Branch("ymc", &ymc);
t1->Branch("zmc", &zmc);
// Branch for V0 positions
t1->Branch("xvzero", &xvzero);
t1->Branch("yvzero", &yvzero);
t1->Branch("zvzero", &zvzero);
// Branch for differences
t1->Branch("diff_x", &diff_x);
t1->Branch("diff_y", &diff_y);
t1->Branch("diff_z", &diff_z);
// Branch for event info
t1->Branch("has_v0_collection", &has_v0_collection);
t1->Branch("event_number", &event_number);
t1->Branch("num_matches", &num_matches);
// Statistics
int total_matches = 0;
int events_without_v0 = 0;
int events_with_v0 = 0;
int nEventsInFile = lcReader->getNumberOfEvents();
cout << "File has " << nEventsInFile << " events" << endl;
// Main event loop
for (int i = 0; i < nEventsInFile; i++) {
if (i % 1000 == 0) cout << "Processing event " << i <<
endl;
// Clear per-event vectors
xmc.clear();
ymc.clear();
zmc.clear();
xvzero.clear();
yvzero.clear();
zvzero.clear();
diff_x.clear();
diff_y.clear();
diff_z.clear();
event_number = i;
num_matches = 0;
has_v0_collection = 0;
// Event-level vectors - now includes MC vertex info
paired with reco particles
std::vector<ReconstructedParticle*> reco_vec;
std::vector<Vertex*> ver_reco_vec;
std::vector<std::array<double, 3>> mc_vertex_vec; //
Store MC vertex for each reco particle
std::vector<ReconstructedParticle*> asos_vec;
std::vector<Vertex*> ver_asos_vec;
// Read event
LCEvent* evt = lcReader->readNextEvent();
LCCollection* linkCol = evt->getCollection("MCTruthRecoLink");
LCCollection* mcskimmed = evt->getCollection("MCParticlesSkimmed");
LCRelationNavigator* nav = new LCRelationNavigator(linkCol);
// Try to get V0Vertices collection
LCCollection* vero = nullptr;
try {
vero = evt->getCollection("V0Vertices");
has_v0_collection = 1;
events_with_v0++;
} catch (lcio::DataNotAvailableException&) {
has_v0_collection = 0;
events_without_v0++;
histo_event->Fill(i);
}
int mn = mcskimmed->getNumberOfElements();
// Loop through MC particles
for (int m = 0; m < mn; m++) {
MCParticle* mc = dynamic_cast<MCParticle*>
(mcskimmed->getElementAt(m));
if (!mc) continue;
int p = mc->getPDG();
// Look for V0 particles (K-short, Lambda, photon)
if (p == 310 || p == -310 || p == -5122 || p == 5122|| p == 22) {
std::vector<MCParticle*> daughters = mc->getDaughters();
if (daughters.size() == 0) continue;
// Get MC truth vertex position
double mc_x = mc->getVertex()[0];
double mc_y = mc->getVertex()[1];
double mc_z = mc->getVertex()[2];
// Fill MC histograms only once per V0 particle
histo->Fill(mc_x);
histo1->Fill(mc_y);
histo2->Fill(mc_z);
// Loop through daughters to find reconstructedparticles
for (size_t j = 0; j < daughters.size(); j++) {
std::vector<float> weight = nav->getRelatedToWeights(daughters.at(j));
if (weight.size() == 0) continue;
// Find best match by weight
int max1 = std::max_element(weight.begin(),
weight.end()) - weight.begin();
ReconstructedParticle* recopar_mass1 = dynamic_cast<ReconstructedParticle*>(
nav->getRelatedToObjects(daughters.at(j))[max1]);
if (!recopar_mass1) continue;
Vertex* reco_ver1 = recopar_mass1->getStartVertex();
if (reco_ver1 == nullptr) continue;
// FIXED: Store MC vertex position paired with each reconstructed particle
reco_vec.push_back(recopar_mass1);
ver_reco_vec.push_back(reco_ver1);
mc_vertex_vec.push_back({mc_x, mc_y, mc_z});
}
}
}
// Process V0Vertices collection if it exists
if (vero != nullptr) {
int nVertexrec = vero->getNumberOfElements();
for (int v = 0; v < nVertexrec; v++) {
Vertex* vertexzero = dynamic_cast<Vertex*>(vero->getElementAt(v));
if (!vertexzero) continue;
float ver_position[3] = {
vertexzero->getPosition()[0],
vertexzero->getPosition()[1],
vertexzero->getPosition()[2] };
	
histo1vz->Fill(ver_position[0]);
histo2vz->Fill(ver_position[1]);
histo3vz->Fill(ver_position[2]);
xvzero.push_back(ver_position[0]);
yvzero.push_back(ver_position[1]);
zvzero.push_back(ver_position[2]);

// Get associated particle
ReconstructedParticle* asos = vertexzero->getAssociatedParticle();
if (asos == nullptr) continue;
Vertex* ver_asos = asos->getStartVertex();
if (ver_asos == nullptr) continue;
asos_vec.push_back(asos);
ver_asos_vec.push_back(ver_asos);
}
}
// FIXED: Compare particles by vertex position, now with
correct MC vertex pairing
for (size_t j = 0; j < reco_vec.size(); j++) {
Vertex* v1 = ver_reco_vec[j];
if (v1 == nullptr) continue;
for (size_t l = 0; l < asos_vec.size(); l++) {
Vertex* v2 = ver_asos_vec[l];
if (v2 == nullptr) continue;
// Calculate 3D distance between vertices
float dx = v1->getPosition()[0] - v2->getPosition()[0];
float dy = v1->getPosition()[1] - v2->getPosition()[1];
float dz = v1->getPosition()[2] - v2->getPosition()[2];
float dist = sqrt(dx*dx + dy*dy + dz*dz);
// If vertices are close enough (within 0.1 mm),consider them matched
if (dist < 0.1) {
num_matches++;
total_matches++;
// Store MC vertex corresponding to this
match
xmc.push_back(mc_vertex_vec[j][0]);
ymc.push_back(mc_vertex_vec[j][1]);
zmc.push_back(mc_vertex_vec[j][2]);
// Store differences between MC and reconstructed
float diff_x_val = mc_vertex_vec[j][0] - v2->getPosition()[0];
float diff_y_val = mc_vertex_vec[j][1] - v2->getPosition()[1];
float diff_z_val = mc_vertex_vec[j][2] - v2->getPosition()[2];
diff_x.push_back(diff_x_val);
diff_y.push_back(diff_y_val);
diff_z.push_back(diff_z_val);
// Fill difference histograms
h_diff_x->Fill(diff_x_val);
h_diff_y->Fill(diff_y_val);
h_diff_z->Fill(diff_z_val);
}}}
// Fill tree for this event
t1->Fill();
delete nav;
}
// Print statistics
cout << "\n=== Analysis Summary ===" << endl;
cout << "Total events processed: " << nEventsInFile << endl;
cout << "Events with V0 collection: " << events_with_v0 <<
endl;
cout << "Events without V0 collection: " <<
events_without_v0 << endl;
cout << "Total matched particles: " << total_matches <<
endl;
// Draw histograms
c1->cd(1);
histo->Draw();
c1->cd(2);
histo1vz->Draw();
c1->Update();
c2->cd(1);
histo1->Draw();
c2->cd(2);
histo2vz->Draw();
c2->Update();
c3->cd(1);
histo2->Draw();
c3->cd(2);
histo3vz->Draw();
c3->Update();
c4->cd(1);
h_diff_x->Draw();
c4->cd(2);
h_diff_y->Draw();
c4->Update();
// Save canvases
c1->SaveAs("mc_vs_reco_x.png");
c2->SaveAs("mc_vs_reco_y.png");
c3->SaveAs("mc_vs_reco_z.png");
c4->SaveAs("differences.png");
// Write and close ROOT file
t1->Write();
froot->Write();
froot->Close();
// Close LCIO reader
lcReader->close();
// Clean up
delete c1;
delete c2;
delete c3;
delete c4;
delete froot;
delete lcReader;
cout << "\nAnalysis complete! Output saved to
VZeroTree.root" << endl;
	
}
