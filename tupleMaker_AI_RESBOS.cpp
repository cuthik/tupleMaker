#ifndef makeAiProfile_cpp
#define makeAiProfile_cpp

/**
 * @file makeAiProfile.cpp
 * Description of this file
 *
 * @brief Creating Ai profiles.
 * The execution of AiMoments from Maarten.
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-02-18
 */
// package
#include "ReadResbosROOT.C"
#include "AiMoments.h"
#include "HistogramsVB.h"
// ROOT
#include "TStopwatch.h"
// C
#include <cstdlib>
#include <iomanip>

using TMath::Log2;

class TupleMaker_min{
    public:
        ~TupleMaker_min(){};
        TupleMaker_min(){};

        void MakeProfile(){
            ResbosRootNtuple ntup(input_fname);
            vector<AiMoments*> ai_files;
            vector<AiMoments*> ai_files_poswgt;
            vector<ResbosRootNtuple*> weights_trees;
            // init ai file, if there are some weights then init ai files and init weights
            size_t Nwgts = weights_fname.size() ? weights_fname.size() : 1;
            for (size_t i_wgts=0; i_wgts< Nwgts;i_wgts++){
                TString ai_outName="AiMoments";
                if (weights_fname.size()){
                    // non-central -- not using weight of 
                    weights_trees.push_back( new ResbosRootNtuple (weights_fname[i_wgts].Data(), "WEIGHT"));
                    ai_outName+="_"; ai_outName+= (int) i_wgts;
                }
                //ai_outName+=".root";
                ai_files.push_back(new AiMoments(ai_outName.Data()));
                ai_files[i_wgts]->Initialize();
                ai_outName.ReplaceAll("AiMoments","AiMoments_poswgt");
                ai_files_poswgt.push_back(new AiMoments(ai_outName.Data()));
                ai_files_poswgt[i_wgts]->Initialize();
            }
            if (weights_fname.size()){
                ntup.ConnectWeights(&weights_trees);
            }
            // loop tree
            Long64_t Nevents = ntup.GetEntries();
            for (int ievt=0; ievt< Nevents; ievt++){
                printEvent(ievt,Nevents);
                ntup.GetEntry(ievt);
                ntup.GetConnectedEntry(ievt);
                ntup.GetLorentzVector();
                for (size_t i_wgts=0; i_wgts<ai_files.size();i_wgts++){
                    double weight = ntup.wgt_corrected;
                    if (weights_fname.size()){
                        weight*=weights_trees[i_wgts]->wgt/ntup.wgt;
                    }
                    //if (fabs(weight)>10) printf(" weight % e   ratio % e", weight , weights_trees[i_wgts]->wgt/ntup.wgt);
                    //if ( fabs(fabs(ntup.v_VB.Rapidity())-3) > 0.2 ) continue;
                    ai_files[i_wgts]->Execute(
                            ntup.v_l1.Pt(), ntup.v_l1.Eta(), ntup.v_l1.Phi(), ntup.v_l1.M(), ntup.type1,
                            ntup.v_l2.Pt(), ntup.v_l2.Eta(), ntup.v_l2.Phi(), ntup.v_l2.M(), ntup.type2,
                            ntup.ebeam, weight
                            );
                    if (weight < 0) continue;
                    ai_files_poswgt[i_wgts]->Execute(
                            ntup.v_l1.Pt(), ntup.v_l1.Eta(), ntup.v_l1.Phi(), ntup.v_l1.M(), ntup.type1,
                            ntup.v_l2.Pt(), ntup.v_l2.Eta(), ntup.v_l2.Phi(), ntup.v_l2.M(), ntup.type2,
                            ntup.ebeam, weight
                            );
                }
            }
            printEvent(Nevents,Nevents);
            TString fname  = "";
            for (size_t i_wgts=0; i_wgts<ai_files.size();i_wgts++){
                fname  = ai_files[i_wgts]->GetFileName();
                ai_files[i_wgts]->Finalize();
                DumpInfoToFile( fname.Data() );
                fname  = ai_files_poswgt[i_wgts]->GetFileName();
                ai_files_poswgt[i_wgts]->Finalize();
                DumpInfoToFile( fname.Data() );
            }
            DumpInfoToFile(input_fname.Data());
        }

        void MakeNtupleWithWeights(bool makeProfiles=false){
            ResbosRootNtuple in_tree (input_fname);
            ResbosRootNtuple out_tree (output_fname,"RECREATE");
            in_tree.GetEntry(0);
            bool fillHistograms = true;
            // init hist
            HistogramsVB hists;
            if(fillHistograms){
                hists.output_name = TString(output_fname).ReplaceAll(".root",".hists.root");
                hists.currentUnitInMeV=1e3;
                hists.N_reweightings = weights_fname.size();
                hists.Init();
            }
            // init weight trees
            vector<ResbosRootNtuple*> weights_trees;
            for (size_t i = 0; i < weights_fname.size(); i++) {
                weights_trees.push_back( new ResbosRootNtuple (weights_fname[i], "WEIGHT"));
            }
            // connect variables orginal to new 
            out_tree.ConnectTree(&in_tree);
            out_tree.ConnectWeights(&weights_trees);
            // loop events
            Long64_t Nevents = out_tree.GetExpectedEntries();
            //Nevents = 1000000;
            for (int ievt=0; ievt<Nevents; ievt++){
                printEvent(ievt,Nevents);
                out_tree.GetConnectedEntry(ievt);
                out_tree.Fill();
                if (fillHistograms){
                    in_tree.GetLorentzVector();
                    hists.SetTruthKinematics(in_tree.v_VB, in_tree.typeVB, in_tree.v_l1, in_tree.type1, in_tree.v_l2, in_tree.type2);
                    hists.CalculateQuickReconstruction();
                    hists.SetWeight(in_tree.wgt_corrected, out_tree.vec_wgt);
                    hists.Fill();
                }
            }
            printEvent(Nevents,Nevents);
            // save outputs
            out_tree.Save();
            DumpInfoToFile(output_fname);
            if(fillHistograms){
                hists.Save();
                DumpInfoToFile(hists.output_name);
            }
            for (size_t i = 0; i < weights_fname.size(); i++) {
                delete weights_trees[i];
            }
            weights_trees.clear();
        }

        void DumpInfoToFile(TString file_name){
            // load tstring from config_dump
            TString s_infodump = "";
            if (!s_infodump.Length()) s_infodump = ResbosRootNtuple::GetStringFromFile("config_dump.txt");
            if (!s_infodump.Length()) s_infodump = ResbosRootNtuple::GetStringFromFile(input_fname.Data());
            // add weight list if needed
            for (size_t i = 0; i < weights_fname.size(); i++) s_infodump += TString::Format("weight_%d:   %s\n",i,weights_fname[i].Data());
            TFile *toFile = TFile::Open(file_name.Data(), "UPDATE");
            toFile->cd();
            TObjString o_infodump(s_infodump);
            o_infodump.Write("ConfigurationDump");
            toFile->Write();
            toFile->Close();
        }

        void printEvent(double evn, double total=0){
            if (evn==0){
                swatch.Start();
            }
            if( ! fmod( Log2(evn), 1 ) ) {
                cout << "Processed event: ";
                if (total !=0 ){
                    cout << " " ;
                    cout << fixed ;
                    cout << setfill('0') ;
                    cout << setw(5) ;
                    cout << setprecision(2) ;
                    cout << evn/total*100 ;
                    cout << "%  ";
                } 
                cout << setprecision(6) << setfill(' ') << setw(10) ;
                cout << (size_t) evn;
                // time
                double rtm = swatch.RealTime(); swatch.Continue();
                double Hz = evn/rtm;
                cout << "   time elapsed ";
                cout << rtm << "s";
                cout << " ( "<< Hz/1000 <<" kHz)";
                if (total !=0 ){
                    cout << "   remaining time ";
                    cout << (total-evn)/Hz;
                }
                cout << endl;
                //cout << "\r";
                //cout.flush();
            }
        }
        // data members
        TString output_fname;
        TString input_fname;
        vector<TString> weights_fname;
        TStopwatch swatch;
};

/**
 * Description of main program
 *
 */

void help(){
    cout << "USAGE1:" << endl;
    cout << " tupleMaker_AI_RESBOS input.root [weight.root ...]" << endl;
    cout << " Description : create AiMoments.root with angular profiles from RESBOS tree." << endl;
    cout << "               It will update the input.root with information about production (config_dump.txt)." << endl<<endl<<endl;

    cout << "USAGE2:" << endl;
    cout << " tupleMaker_AI_RESBOS -w output.root input.root weights_01.root [ weights_02.root ...]" << endl;
    cout << " Description : create output.root with  RESBOS tree and vector of weigths." << endl<<endl<<endl;

}

int main(int argc, const char * argv[]){

    TupleMaker_min tm; 
    if (argc < 2) {help(); return 1;}
    TString tmp = argv[1];
    if (!tmp.CompareTo("-h") || !tmp.CompareTo("--help")) {
        help();
        return 0;
    } else if (!tmp.CompareTo("-w")) {
        if (argc < 5 ) { help(); return 1;}
        //USAGE 2
        tm.output_fname = argv[2];
        tm.input_fname = argv[3];
        for (size_t i=4; i<argc; i++) tm.weights_fname.push_back(argv[i]);
        tm.MakeNtupleWithWeights();
        return 0;
    } else if(argc>=2) {
        // USAGE 1
        //tm.ebeam       = ResbosRootNtuple::GetConfValue<double>("config_dump.txt","ebeam");
        tm.input_fname = argv[1];
        if (argc > 2) for (int i_wt = 2; i_wt < argc; i_wt++){
            tm.weights_fname.push_back(argv[i_wt]);
        }
        tm.MakeProfile();
        return 0;
    } else {
        help(); return 1;
    }


    return 0;
}


#endif // makeAiProfile_cpp
