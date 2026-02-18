#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "Particle.h"
#include "CYL.h"
#include "Point.h"
#include "VTX.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TStopwatch.h"
#include "TH1.h"
#include "TH3D.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include <cmath>
#include <vector>
#include <algorithm>

double get_eta(double theta);

void reconstruction(){

    //monitoring CPU time
    TStopwatch timer;
    timer.Start();   // start timer

    //constants
    const double kPi = 3.14159265358979323846;

    //creating directory for plots, if it doesnt exist
    if (gSystem->AccessPathName("plots")) gSystem->mkdir("plots",1);

    //setting minimum momentum for straight tracks
    double pc = 1e3; //MeV

    //conversion factor cm->um
    double cm_um = 10000;
    //setting bin width for vertex reconstruction
    double dz = 1.; // cm
    double acceptance = 30; // cm

    //reporting smearing parameters
    double sigma_smearing_Z = 0.012; // cm (120 um)
    double sigma_smearing_AR = 0.003; // cm (30 um)

    //reporting generation parameter (as known from beam)
    double sigmaXY = 0.01; // cm

    //open input file
    TFile hfile("hits_data.root");
    //reading TTree
    TTree *tree = (TTree*)hfile.Get("tree");

    //declaring vertex pointer and TClonesArrays 
    VTX *vertex = new VTX();
    TClonesArray *ptrhits1 = new TClonesArray("Point", 100);
    TClonesArray &hits1 = *ptrhits1; // hits1 is the object pointed by ptrhits1
    TClonesArray *ptrhits2 = new TClonesArray("Point", 100);
    TClonesArray &hits2 = *ptrhits2; //same for layer 2

    //declaring branches and setting addresses
    tree->SetBranchAddress("VertexMultiplicity", &vertex);
    tree->SetBranchAddress("Points_L1", &ptrhits1);
    tree->SetBranchAddress("Points_L2", &ptrhits2);

    //extracting cylinders geometry
    CYL *beam_pipe = (CYL*)hfile.Get("BeamPipeGeometry");
    CYL *layer_1 = (CYL*)hfile.Get("Layer1Geometry");
    CYL *layer_2 = (CYL*)hfile.Get("Layer2Geometry");

    //get radius of detectors
    double R1 = layer_1->GetR();
    double R2 = layer_2->GetR();
    //get Z extension of detector
    double H = layer_1->GetH();

    //declare histo to collect intersections (vertices)
    int nbins = std::ceil(2*acceptance/dz);
    TH1D* hIntersections = new TH1D("hIntersections", "Histogram of intersections", nbins, -acceptance, acceptance);

    //declare vector to collect intersections (vertices)
    std::vector<double> *vIntersections = new std::vector<double>();
    vIntersections->reserve(1000); //reserve 1000 slots in memory

    //TClonesArrays to collect reconstructed tracklets and particles
    TClonesArray *ptrtracklets = new TClonesArray("Particle", 300);
    TClonesArray &tracklets = *ptrtracklets;
    TClonesArray *ptrparticles = new TClonesArray("Particle", 150);
    TClonesArray &particles = *ptrparticles;

    //vectors to optimize CPU use in loops, and reserve them memory
    std::vector<double> h1_Phi, h1_Z; 
    std::vector<double> h2_Phi, h2_Z;
    h1_Phi.reserve(160), h1_Z.reserve(160);
    h2_Phi.reserve(160), h2_Z.reserve(160);

    //we use L1/L2 ThetaRMS and smearing in order to estimate DeltaPhi to reconstruct tracklets
    double ThetaRMS = layer_1->GetThetaRMS(pc);    
    //sigma Phi for smearing is the biggest one, from L1
    double sigma_smearing_Phi = sigma_smearing_AR/layer_1->GetR();
    //sigma Phi is also affected by gaussian distribution of x,y of vertex
    //considering, as the difference from the Phi of coord:
    //for L1, DeltaPhi1 = 3sigmaXY/R1
    //for L2, DeltaPhi2 = 3sigmaXY/R2
    //DeltaPhi = DeltaPhi1 - DeltaPhi2 because they are not indipendent

    //dPhi is 3*sigma, with sigma the square sum of errors
    double dPhi = std::sqrt(sigma_smearing_Phi*sigma_smearing_Phi +
                            ThetaRMS*ThetaRMS +
                            (sigmaXY*(1/R1 - 1/R2))*(sigmaXY*(1/R1 - 1/R2)));

    //cout << "Si ThetaRMS: " << ThetaRMS << "rad\n";
    //cout << "dPhi = " << dPhi << endl;

    //variabile size bins for multiplicity and Z, uniform for residuals
    std::vector<double> edges_multiplicity = {0, 2, 4, 6, 8, 10, 15, 20, 30, 40, 55, 70, 85, 100, 115, 140};
    int nbins_mult = edges_multiplicity.size() - 1;

    int nbins_residuals = 100;
    std::vector<double> edges_residuals(nbins_residuals + 1);
    for(int i=0; i <= nbins_residuals; i++) edges_residuals[i] = -1000. + i*(2*1000. / nbins_residuals);
  
    std::vector<double> edges_Ztrue = {-18,-15,-12,-9,-6,-4,-2,0,2,4,6,9,12,15,18};
    int nbins_Ztrue = edges_Ztrue.size() - 1;



    //___DECLARE HISTOS____
    gStyle->SetCanvasColor(10);
    gStyle->SetMarkerStyle(22);
    gStyle->SetMarkerColor(4);
    gStyle->SetMarkerSize(1);
    TH1D* hResiduals = new TH1D("hResiduals", "Residuals", nbins_residuals, edges_residuals.data()); //residuals z_rec-z_true
    TH1D* hPassedZ = new TH1D("hPassedZ", "", nbins_Ztrue, edges_Ztrue.data()); //z0
    TH1D* hPassedMult = new TH1D("hPassedMult", "", nbins_mult, edges_multiplicity.data()); //mult
    TH1D* hTotalZ = new TH1D("hTotalZ", "", nbins_Ztrue, edges_Ztrue.data()); //z0
    TH1D* hTotalMult = new TH1D("hTotalMult", "", nbins_mult, edges_multiplicity.data()); //mult

    //3D histo to get resolution
    TH3D* hResidualMultZtrue = new TH3D("hResidualMultZtrue", "", nbins_mult, edges_multiplicity.data(), nbins_Ztrue, edges_Ztrue.data(), nbins_residuals, edges_residuals.data());

    TH1D* hResMult = new TH1D("gResMultrec", "Resolution vs Multiplicity", nbins_mult, edges_multiplicity.data()); 
    TH1D* hResZ = new TH1D("gResZ", "Resolution vs Z_{true}", nbins_Ztrue, edges_Ztrue.data()); 
 

    //looking into TTree looping on events
    int nevents = tree->GetEntries();
    for(int ev=0; ev<nevents; ev++){ 
        if (ev%5000 == 0)cout << "Reconstructing vertex n. " << ev << endl;

        //select i-event
        tree->GetEvent(ev);
        //get number of hits on L1 and L2
        int nHits1=ptrhits1->GetEntries();
        int nHits2=ptrhits2->GetEntries();

        //counter of tracklets
        int n_track = 0;

        //faster loop to fill vectors
        for(int i=0; i<nHits1; i++){ //hits on L1
            Point* p1 = (Point*)hits1[i]; //pointer to i-th hit on L1
            h1_Phi.push_back(p1->GetPhi());
            h1_Z.push_back(p1->GetZ());
        }
        //loop on L2
        for(int j=0; j<nHits2; j++){ //hits on L1
            Point* p2 = (Point*)hits2[j]; //pointer to i-th hit on L1
            h2_Phi.push_back(p2->GetPhi());
            h2_Z.push_back(p2->GetZ());
        }

        //sort vectors

        //loop on hits
        for(int i=0; i<nHits1; i++){
            for(int j=0; j<nHits2; j++){

                //check if Hit1 and Hit2 have an azimuth difference less than dPhi
                if (std::abs(h1_Phi[i] - h2_Phi[j]) < dPhi){
                    double z1 = h1_Z[i];
                    double z2 = h2_Z[j];

                    //find intersections with z axis, formula for straight line, 2 points
                    double z0 = z1 - R1*(z2-z1)/(R2-R1); 

                    hIntersections->Fill(z0);  //fill histo
                    
                    vIntersections->push_back(z0); //fill vector with intersections

                    //calculate eta, phi (mean of phi1 and phi2) and create Particle
                    double tracklet_Phi = (h1_Phi[i] + h2_Phi[j])/2;
                    if (tracklet_Phi <0) tracklet_Phi += 2*kPi; //positive Phi
                    double tracklet_Theta = std::atan2((R2-R1), (z2-z1));
                    double tracklet_Eta = get_eta(tracklet_Theta);

                    //fill TClonesArray with tracklets
                    new(tracklets[n_track]) Particle(tracklet_Phi, tracklet_Eta, z0);
                    n_track++;

                }
            }
        }

        //get bin with maximum and bin limits
        int bin_max = hIntersections->GetMaximumBin();
        //double bin_width = hIntersections->GetBinWidth(bin_max);
        //double bin_center = hIntersections->GetBinCenter(bin_max);
        double z_low = hIntersections->GetXaxis()->GetBinLowEdge(bin_max);
        double z_up = hIntersections->GetXaxis()->GetBinUpEdge(bin_max);
        //int max_count_bin = hIntersections->GetMaximum(bin_max);

        //sort vector
        std::sort(vIntersections->begin(), vIntersections->end());
        //find lower and upper iterators, to limit z range in vector
        auto low_iter = std::lower_bound(vIntersections->begin(), vIntersections->end(), z_low);
        auto up_iter = std::upper_bound(vIntersections->begin(), vIntersections->end(), z_up);

        //calculate mean value of region between limits, in vector
        double z_mean = TMath::Mean(low_iter, up_iter);
        double z_rms = TMath::RMS(low_iter, up_iter); 
        
        //particles counter
        double n_particle = 0;

        //find multiplicity: particles have intersection within 1 sigma from z_mean
        for(int k=0; k<n_track; k++){
            if(std::abs(((Particle*)tracklets[k])->GetZ0() - z_mean) < z_rms){
                new(particles[n_particle]) Particle(*((Particle*)tracklets[k]));
                n_particle++;
            }
        }

        //DEBUG OUTPUT
        //cout << "nHits1 " << nHits1 << ",   nHits2 " << nHits2 << endl;
        //cout << "n_tracks = " << n_track <<",   bin_width = " << bin_width << ",   z_rms = " << z_rms << endl;
        //cout << "z_mean = " << z_mean << ",  z_rms = " << z_rms  << ",   z_center = " << bin_center << ",  max value = " << max_count_bin << endl;
        //cout << "particles " << n_particle<< endl;


        //getting z true and mult true value
        double z_true = vertex->GetZ();
        double mult_true = vertex->GetMult();

        //fill histo for generated vertices, to calculate efficiency
        hTotalZ->Fill(z_true);
        hTotalMult->Fill(mult_true);

        //fill histos if we accept the vertex (at least 2 tracklets)
        if(n_particle>=2){
            hResiduals->Fill((z_mean-z_true)*cm_um);  //residuals in um
            hPassedZ->Fill(z_true); //for efficiency; we use true values for this histo in order to not distort efficiency
            hPassedMult->Fill(mult_true);

            hResidualMultZtrue->Fill(mult_true, z_true, (z_mean-z_true)*cm_um); //fill histo for resolution

        }

        //clear our containers
        vIntersections->clear(); 
        hIntersections->Reset();
        ptrtracklets->Clear();
        ptrparticles->Clear();
        h1_Phi.clear();
        h1_Z.clear();
        h2_Phi.clear();
        h2_Z.clear();
        
    }

    //generate efficiency graphics
    //efficiency vs z0
    TEfficiency* gEffZtrue = new TEfficiency(*hPassedZ, *hTotalZ);
    //gEffZtrue->SetStatisticOption(TEfficiency::kBBayesian); //bayesian error
    //efficiency vs Mult
    TEfficiency* gEffMult = new TEfficiency(*hPassedMult, *hTotalMult);
    //gEffMult->SetStatisticOption(TEfficiency::kBBayesian); //bayesian error

    //__Calculate resolution_________________________________________
    double X_nbins = hResidualMultZtrue->GetXaxis()->GetNbins();
    double Y_nbins = hResidualMultZtrue->GetYaxis()->GetNbins();
    TH1D* proj_pointer = nullptr;

    //loop on multiplicity bins, to calculate rms of residuals and fill resolution graphic
    for(int i=1; i<=X_nbins; i++){
        proj_pointer = hResidualMultZtrue->ProjectionZ("proj_pointer", i, i, 1, Y_nbins); //take i-th X, all Y
        if (proj_pointer->GetEntries() > 20){        //only good statistics
            double rms = proj_pointer->GetRMS();        //get rms of residuals for i-th multiplicity -> resolution
            double rms_error = proj_pointer->GetRMSError(); //error of rms
            hResMult->SetBinContent(i, rms);     //fill TH1 with mult on X and rms of residuals on Y
            hResMult->SetBinError(i, rms_error); //set error 
        }
        delete proj_pointer;
    }

    //loop on z_true bins, to calculate rms of residuals and fill resolution graphic
    for(int j=1; j<=Y_nbins; j++){
        proj_pointer = hResidualMultZtrue->ProjectionZ("proj_pointer", 1, X_nbins, j, j); //take all X, j-th Y
        if (proj_pointer->GetEntries() > 20){        //only good statistics
            double rms = proj_pointer->GetRMS();        //get rms of residuals for i-th multiplicity -> resolution
            double rms_error = proj_pointer->GetRMSError(); //error of rms
            hResZ->SetBinContent(j, rms);
            hResZ->SetBinError(j, rms_error);
        }
        delete proj_pointer;
    }

    //_________DRAW HISTOS AND GRAPHICS_________________

    gStyle->SetOptStat(1111);
    //residuals
    TCanvas* cResiduals = new TCanvas("cResiduals","",1200,800);
    hResiduals->SetDirectory(0);
    hResiduals->SetTitle("Residuals");
    hResiduals->SetXTitle("(z_{rec}-z_{true}) [#mum]");
    hResiduals->SetYTitle(Form("Entries/(%.f #mum)", hResiduals->GetBinWidth(1)));
    hResiduals->Draw("PE1");
    cResiduals->Print("plots/residuals.png");

    gStyle->SetOptStat(0);

    //resolution vs generated multiplicity
    TCanvas* cResMult = new TCanvas("cResMult","",1200,800);
    hResMult->SetDirectory(0);
    hResMult->SetTitle("Resolution vs Multiplicity; Multiplicity; Resolution #sigma_{z} [#mum]");
    hResMult->Draw("PE1");
    cResMult->Print("plots/resolution_mult.png");

    //resolution vs ztrue
    TCanvas* cResZ = new TCanvas("cResZ","",1200,800);
    hResZ->SetDirectory(0);
    hResZ->SetTitle("Resolution vs z_{true}; z_{true} [cm]; Resolution #sigma_{z} [#mum]");
    hResZ->Draw("PE1");
    cResZ->Print("plots/resolution_ztrue.png");

    //efficiency vs ztrue
    TCanvas* cEffZtrue = new TCanvas("cEffZ0","",1200,800);
    gEffZtrue->SetDirectory(0);
    gEffZtrue->SetTitle("Efficiency vs Z_{true};z_{true} [cm];Efficiency");
    gEffZtrue->Draw("APE");
    cEffZtrue->Update();
    gEffZtrue->GetPaintedGraph()->SetMinimum(0);
    gEffZtrue->GetPaintedGraph()->SetMaximum(1.2);
    cEffZtrue->Update();
    cEffZtrue->Print("plots/efficiency_ztrue.png");

    //efficiency vs multiplicity
    TCanvas* cEffMult = new TCanvas("cEffMult","",1200,800);
    gEffMult->SetDirectory(0);
    gEffMult->SetTitle("Efficiency vs Multiplicity; Multiplicity; Efficiency");
    gEffMult->Draw("APE");
    cEffMult->Update();
    gEffMult->GetPaintedGraph()->SetMinimum(0);
    gEffMult->GetPaintedGraph()->SetMaximum(1.2);
    cEffMult->Update();
    cEffMult->Print("plots/efficiency_mult.png");

    timer.Stop();
    timer.Print();

}

double get_eta(double theta){
    return -std::log(std::tan(theta/2));
}

