#ifndef tupleMaker3
#define tupleMaker3

/**
 * @file tupleMaker3.cxx
 * Description of this file
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2014-07-29
 */

#include "Output.hpp"

#include "TFile.h"
//#include "TTree.h"
#include "TMath.h"
//#include "TString.h"

#include <iostream>
#include <stdexcept>
//#include <vector>
//#include <cmath>

using TMath::Log2;
using std::runtime_error;
using std::vector;
using std::cout;
using std::endl;


class TupleMaker {
    public:
        TupleMaker():
            hepfilename       (""),
            //weightfileNames (),
            output            (0),
            hepfile           (0),
            f_out             (0),
            doConvertWeight   (false),
            doReweighting     (false),
            doSavePartonInfo  (false)
        {
            ClearEvent();
            ClearParticle();
        }
        ~TupleMaker(){}


        inline void SetKinematic(TString hepfile){ hepfilename = hepfile;}
        inline void SetConvertWeight (bool in=true) { doConvertWeight =in;}
        inline void SetDoReweighting (bool in=true) { doReweighting   =in;}
        inline void SetSavePartonInfo (bool in=true) { doSavePartonInfo   =in;}


        void SetOutName(TString outname){
            f_out = TFile::Open(outname, "RECREATE");
            if(!f_out){
                cout << "Can not open output file: " << outname.Data() << endl;
            }
        }


        inline void AddWeights(TString weightfile){ weightfileNames.push_back(weightfile); }

        inline void ClearEvent(){
            evn = 0;
            evt_wt = -1;
            vx = vy = vz = -999.9;
            Q2 = x1 = x2 = flav1 = flav2 = -999.9;
            weight_PDF.clear();
        }


        inline void ClearParticle(){
            id = origin = udk = 0;
            px = py = pz = E = -999.;
        }


        void MakeWeightRootfile(){
            runtime_error("Make Weight Rootfile not implemented yet.");
        }

        void OpenWeightFiles(){
            runtime_error("Open Weight files not implemented yet.");
        }

        void CloseWeightFiles(){
            runtime_error("Close Weight files not implemented yet.");
        }


        void Loop(){
            // CONVERT TXT TO ROOT FILE AND DIE
            if (doConvertWeight){ 
                MakeWeightRootfile();
                return;
            }


            // CREATE INPUT FOR PMCS
            // open input files
            if (doReweighting && weightfileNames.size()!=0 ) OpenWeightFiles();


            bool finished_file = false;
            // this loops over events
            while( ! finished_file ) {
                ClearEvent();
                int rc=0;
                // first line (event number and weight)
                rc = fscanf(hepfile,"%i %f", &evn, &evt_wt);
                if (rc!=1) runtime_error("Event line not found.");
                if( ! fmod( Log2(evn), 1 ) ) cout<<"Processed event: "<<evn<<endl;
                if( evn == 0 ) {
                    finished_file = true;
                    continue;
                } 

                // second line (beam position)
                rc=fscanf(hepfile,"%f %f %f", &vx, &vy, &vz);
                if (rc!=1) runtime_error("Vertex line not found.");

                // read weights

                if (doReweighting){
                    if(doSavePartonInfo){
                        // load parton info
                        runtime_error("Save Parton Info not implemented yet.");
                    } else {
                        // load weights
                        if (weightfileNames.size() == 0) runtime_error("Missing weight file names.");
                    }
                }

                output->NewEvent( evn, evt_wt, 0 ,
                                  vx, vy, vz,
                                  Q2, x1, x2, flav1, flav2,
                                  weight_PDF);

                // loop particles
                bool finished_particles=false;
                while (!finished_particles){
                    ClearParticle();

                    // read isajet id
                    rc=fscanf(hepfile,"%i",&id);
                    if (rc!=1) runtime_error("Missing isajet id.");

                    if(id==0){
                        finished_particles=true;
                        continue;
                    }

                    // read kinematics
                    rc=fscanf(hepfile,"%i",&px,&py,&pz,&E,&origin,&udk);
                    if (rc!=1) runtime_error("Missing particle kinematics.");

                    output->AddParticle( id, px,py,pz,E,origin);
                }

                output->Fill();
            }

            //close input files
            if (doReweighting && weightfileNames.size()!=0 ) CloseWeightFiles();
        }


        void Save(){
            f_out->Write();
            f_out->Close();
        }


    private:
        TString hepfilename;
        vector<TString> weightfileNames;
        Output * output;

        FILE * hepfile;
        TFile * f_out;

        // event related
        int evn;
        float evt_wt;
        float vx,vy,vz;
        float Q2, x1, x2, flav1, flav2;
        vector<float> weight_PDF;

        // particle related
        int id,origin,udk;
        float px,py,pz,E;


        bool doConvertWeight;
        bool doReweighting;
        bool doSavePartonInfo;

};


void help(){
    cout << "USAGE 1 :" << endl;
    cout << " tupleMaker3 output.root central.help [-f] [ weight_01.root weight_02.root ...]" << endl;
    cout << "    -f    save parton info (default off)" << endl << endl;
    cout << " Description : create input tree for pmcs from hep and weight files." << endl<<endl<<endl;


    cout << "USAGE 2 :" << endl;
    cout << " tupleMaker3 weitght.txt" << endl << endl;
    cout << " Description : create root weight file from weight text file." << endl<<endl<<endl;

    cout << "To only print this help use '-h' or '--help' as only argument." << endl;
}

/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){

    TupleMaker tm;

    // parse inputs : program output.root central.hep weight.txt ...
    if (argc < 2 ) { help(); return 1;}
    //if (argc == 2) { // create root file from weight text file (or help)
        //TString tmp = argv[1];
        //if (!tmp.CompareTo("-h") || !tmp.CompareTo("--help")) {
            //help();
            //return 0;
        //}
        //tm.SetConvertWeight(true);
        //tm.AddWeights(tmp);
        //tm.SetOutName(tmp.ReplaceAll(".txt",".root"));
    //}
    //if (argc > 2) { // minimum is output and input
        //tm.SetOutName(argv[1]);
        //tm.SetKinematic(argv[2]);
    //}
    //if (argc > 4 ){ // files and setting for weights
        //int start=3;
        //tm.SetDoReweighting(true);
        //// check for parton info saving
        //if ( !TString(argv[start]).CompareTo("-f") ){
            //tm.SetSavePartonInfo(true);
            //start++;
        //}
        //for (int i=start; i<=argc;i++){
            //tm.AddWeights(argv[i]);
        //}
    //}

    //tm.Loop();
    //tm.Save();

    return 0;
}


#endif // tupleMaker3
