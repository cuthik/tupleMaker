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

#include "ReadResbosROOT.C"
#include "AiMoments.h"
#include <cstdlib>

class TupleMaker_min{
    public:
        ~TupleMaker_min(){};
        TupleMaker_min(){};

        void MakeProfile(){
            ResbosRootNtuple ntup(input_fname);
            AiMoments ai_file;
            ai_file.Initialize();
            for (int ievt=0; ievt<ntup.GetEntries(); ievt++){
                ntup.GetEntry(ievt);
                ntup.GetLorentzVector();
                ai_file.Execute(
                        ntup.v_l1.Pt(), ntup.v_l1.Eta(), ntup.v_l1.Phi(), ntup.v_l1.M(), ntup.type1,
                        ntup.v_l2.Pt(), ntup.v_l2.Eta(), ntup.v_l2.Phi(), ntup.v_l2.M(), ntup.type2,
                        ebeam, ntup.wgt_corrected
                        );
            }
            ai_file.Finalize();
        }

        void MakeNtupleWithWeights(){
            ResbosRootNtuple in_tree (input_fname);
            ResbosRootNtuple out_tree (output_fname,"RECREATE");
            in_tree.GetEntry(0);
            printf("  out first weight %f\n",in_tree.wgt);
            vector<ResbosRootNtuple*> weights_trees;
            for (size_t i = 0; i < weights_fname.size(); i++) {
                printf( "opening weight file: %s\n", weights_fname[i].Data());
                weights_trees.push_back( new ResbosRootNtuple (weights_fname[i], "WEIGHT"));
                weights_trees[i]->tree->Print();
                weights_trees[i]->tree->GetEntry(50);
                printf("  first weight %f\n",weights_trees[i]->wgt);
            }
            out_tree.ConnectTree(&in_tree);
            out_tree.ConnectWeights(&weights_trees);
            Long64_t Nevents = out_tree.GetExpectedEntries();
            Nevents = 10;
            for (int ievt=0; ievt<Nevents; ievt++){
                printf("EVENT %d\n",ievt);
                out_tree.GetConnectedEntry(ievt);
                out_tree.Fill();
            }
            out_tree.Save();
            for (size_t i = 0; i < weights_fname.size(); i++) {
                printf( " pointer to root tree %p\n", weights_trees[i]);
                delete weights_trees[i];
            }
            weights_trees.clear();
        }

        // data members
        TString output_fname;
        TString input_fname;
        vector<TString> weights_fname;
        float ebeam;
};

/**
 * Description of main program
 *
 */

void help(){
    cout << "USAGE1:" << endl;
    cout << " makeAiProfile ebeam input.root" << endl;
    cout << " Description : create AiMoments.root with angular profiles from RESBOS tree. Energy of beam is in GeV." << endl<<endl<<endl;

    cout << "USAGE2:" << endl;
    cout << " makeAiProfile -w output.root input.root [weights_01.root ...]" << endl;
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
    } else if(argc==3) {
        // USAGE 1
        tm.ebeam       = atoi(argv[1]);
        tm.input_fname = argv[2];
        tm.MakeProfile();
        return 0;
    } else {
        help(); return 1;
    }


    return 0;
}


#endif // makeAiProfile_cpp
