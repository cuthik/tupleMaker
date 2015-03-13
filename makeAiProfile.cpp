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

class AiMaker{
    public:
        ~AiMaker(){};
        AiMaker(){};

        void MakeProfile(){
            ResbosRootNtuple ntup(f_path_in);
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

        // data members
        TString f_path_in;
        float ebeam;
};

/**
 * Description of main program
 *
 */

void help(){
    cout << "USAGE:" << endl;
    cout << " makeAiProfile ebeam TODO_process input.root" << endl;
    cout << " Description : create AiMoments.root with angular profiles from RESBOS tree. Energy of beam is in GeV." << endl<<endl<<endl;
}

int main(int argc, const char * argv[]){

    if (argc < 3 ) { help(); return 1;}
    else {
        //HELP
        TString tmp = argv[1];
        if (!tmp.CompareTo("-h") || !tmp.CompareTo("--help")) {
            help();
            return 0;
        }

        AiMaker aim;
        aim.ebeam     = atoi(argv[1]);
        //TODO:aim.process = argv[2];
        aim.f_path_in = argv[2];
        aim.MakeProfile();
    }


    return 0;
}


#endif // makeAiProfile_cpp
