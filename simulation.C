#include "Particle.h"
#include "CYL.h"
#include "Point.h"
#include "VTX.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"

//generation of nevents vertexes with a multiplicity of particles each

void simulation(double nevents=100000, bool MS=true, bool assigned_distr=true){
    //MS: bool to turn on/off Multiple Scattering
    //if assigned_distr is true, extracts multiplicity and eta from a given distribution
    //if it's false, extracts eta from uniform (-2,2) and mult from a poissonian with mean value=mean_mult

    //monitoring CPU time
    TStopwatch timer;
    timer.Start();   // start timer

    //variables: data of the problem
    double mean_mult = 10; //mean multilplicity, to extract Mult from a Poissonian
    double pc = 1e3; // momentum in MeV, such that particle go straight

    double sigmaZ = 5.3; //cm
    double sigmaXY = 0.01; // cm
    double sigma_smearing_Z = 0.012; // cm (120 um)
    double sigma_smearing_AR = 0.003; // cm (30 um)

    double H = 27; //cm, detectors extension
    double H0 = 4*H; // beam pipe H
    //beam pipe: Berillium
    double R0 = 3; // cm
    double t0 = 0.08; // cm, thickness
    //Layer 1: Silicon
    double R1 = 4; //cm
    double t1 = 0.02; //cm
    //Layer 2: Silicon
    double R2 = 7; // cm
    double t2 = 0.02; // cm

    //set seed for gRandom
    int seed = 0;
    gRandom->SetSeed(seed);

    //generate cylinders
    CYL* beam_pipe = new CYL("beam_pipe", R0, H0, t0, "Be");
    CYL* layer_1 = new CYL("layer_1", R1, H, t1, "Si");
    CYL* layer_2 = new CYL("layer_2", R2, H, t2, "Si");

    //saving variables to save cpu processes in loops
    //getting RadLenghths
    double RadLen_bp = beam_pipe->GetX0();
    double RadLen_L1 = layer_1->GetX0();
    double RadLen_L2 = layer_2->GetX0();


    //Opening output files
    TFile hfile("hits_data.root","RECREATE");
    //declaring TTree
    TTree *tree = new TTree("tree", "TTree for hit data");

    //Declaring TClonesArrays
    //hits on layer 1
    TClonesArray *ptrhits1 = new TClonesArray("Point", 100);
    TClonesArray &hits1 = *ptrhits1; // hits1 is the object pointed by ptrhits1

    //hits on layer 2
    TClonesArray *ptrhits2 = new TClonesArray("Point", 100);
    TClonesArray &hits2 = *ptrhits2; //same for layer 2


    //declaring pointer of vertex and particle
    VTX* vertex = new VTX();
    Particle* part = new Particle();
    Point* intersection0 = new Point();
    Point* intersection1 = new Point();
    Point* intersection2 = new Point();


    //declaring branches
    tree->Branch("VertexMultiplicity", &vertex);
    tree->Branch("Points_L1", &ptrhits1);
    tree->Branch("Points_L2", &ptrhits2);


    //opening assigned distributions
    TFile *file_distributions = new TFile("kinem.root", "READ");
    TH1F* hMult = (TH1F*)file_distributions->Get("hm"); //multiplicity
    hMult->SetDirectory(0);
    TH1F* hEta = (TH1F*)file_distributions->Get("heta2"); //eta
    hEta->SetDirectory(0);



    //_____SIMULATION ON EVENTS_______
    // loop on events
    for(int ev=0; ev<nevents; ev++){
        if(assigned_distr) vertex->Initialize(mean_mult, hMult);
        else vertex->Initialize(mean_mult);

        if (ev%5000 == 0)cout << "Simulating vertex n. " << ev << endl;

        //counters for Hits
        int nHits1 = 0;
        int nHits2 = 0;

        for(int j=0; j<vertex->GetMult(); j++){

            //generate particle in vertex, with directions
            if(assigned_distr) part->Initialize(*vertex, hEta);
            else part->Initialize(*vertex);

            //transport to beam pipe and find intersections
            //then save them on *intersection0
            part->Transport(R0, H0, *intersection0);
            //calculate MS on beampipe
            part->MultipleScattering(pc, RadLen_bp, t0, MS);

            //transport to layer1, find intersection and MS
            part->Transport(R1, H, *intersection1);
            part->MultipleScattering(pc, RadLen_L1, t1, MS);
            //apply smearing on hit on layer1
            intersection1->SmearingZ(sigma_smearing_Z);
            intersection1->SmearingPhi(sigma_smearing_AR);
            //if there is intersection, save hit
            if (part->GetVisible()){
                new (hits1[nHits1]) Point(*intersection1);
                nHits1++;
            }
            
            //transport to layer2, find intersection and MS
            part->Transport(R2, H, *intersection2);
            part->MultipleScattering(pc, RadLen_L2, t1, MS);
            //smearing on layer2
            intersection2->SmearingZ(sigma_smearing_Z);
            intersection2->SmearingPhi(sigma_smearing_AR);
            //if there is intersection, save hit
            if (part->GetVisible()){
                new (hits2[nHits2]) Point(*intersection2);
                nHits2++;
            }
            
        }

        tree->Fill();
        ptrhits1->Clear();
        ptrhits2->Clear();

    }

    //write geometries in file
    hfile.cd();
    beam_pipe->Write("BeamPipeGeometry");
    layer_1->Write("Layer1Geometry");
    layer_2->Write("Layer2Geometry");

    //save all objects in file and close it
    hfile.Write();
    hfile.Close();

    timer.Stop();
    timer.Print();

}