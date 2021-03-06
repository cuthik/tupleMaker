#ifndef tupleMaker_DYRES
#define tupleMaker_DYRES

/**
 * @file tupleMaker_DYRES.cxx
 * Description of this file
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-01-20
 */

#include "Output.hpp"
#include "ReadGene.C"

#include "TFile.h"
//#include "TTree.h"
#include "TMath.h"
//#include "TString.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
//#include <vector>
#include <cstdio>
#include <cerrno>

using TMath::Log2;
using std::runtime_error;
using std::vector;
using std::ifstream;
using std::cout;
using std::endl;


class TupleMaker {
    public:
        TupleMaker():
            hepfilename       (""),
            //weightfileNames (),
            output            (0),
            //hepfile           (0),
            f_out             (0),
            debug             (false),
            doConvertWeight   (false),
            doReweighting     (false),
            doWeightRatio     (true),
            doSavePartonInfo  (false)
        {
            ClearEvent();
            ClearParticle();
        }
        ~TupleMaker(){
            if (f_out && f_out->IsOpen()) f_out->Close();
        }


        inline void SetKinematic(TString hepfile){ hepfilename = hepfile;}
        inline void SetConvertWeight (bool in=true) { doConvertWeight =in;}
        inline void SetDoReweighting (bool in=true) { doReweighting   =in;}
        inline void SetSavePartonInfo (bool in=true) { doSavePartonInfo   =in;}


        void SetOutName(TString outname){
            f_out = TFile::Open(outname, "RECREATE");
            if(!f_out){
                TString tmp =  "Can not open output file: "; tmp+=outname;
                throw runtime_error(tmp.Data());
            }
        }


        inline void AddWeights(TString weightfile){ weightfileNames.push_back(weightfile.Data()); }

        inline void ClearEvent(){
            evn = 0;
            evt_wt = -1;
            vx = vy = vz = -999.9;
            Q2 = x1 = x2 = flav1 = flav2 = -999.9;
            //weight_PDF.clear();
        }


        inline void ClearParticle(){
            id = origin = udk = 0;
            px = py = pz = E = -999.;
        }


        void GetEntries(TString filename, TString treename="Global"){
            f_out = TFile::Open(filename.Data(),"READ");
            TTree *tt = (TTree*) f_out->Get(treename);
            cout << tt->GetEntries() << endl;
            f_out->Close();
        }


        void MakeWeightRootfile(){
            TTree * tt = new TTree("weights","weights");
            tt->ReadFile(weightfileNames[0],"evt/I:weight/F");
        }

        void OpenWeightFiles(){
            const int N = weightfileNames.size();
            weight_evtn = new Int_t  [N];
            weight_val = new Float_t [N];

            for (size_t iw = 0; iw < weightfileNames.size(); iw++){
                cout << weightfileNames[iw].Data() << endl;
                weight_files.push_back( TFile::Open(weightfileNames[iw].Data(),"READ") );
                if (weight_files.at(iw) == 0) throw runtime_error(weightfileNames[iw].Prepend("Can not open weightfile ").Data());
                weight_trees.push_back( (TTree*)weight_files[iw]->Get("weights") ) ;

                weight_trees[iw]->SetBranchAddress("evt"    , & weight_evtn[iw] );
                weight_trees[iw]->SetBranchAddress("weight" , & weight_val[iw]  );
            }
        }

        void CloseWeightFiles(){
            for (vector<TFile *>::iterator it = weight_files.begin(); it!=weight_files.end(); ++it) (*it)->Close();
            delete[] weight_evtn ;
            delete[] weight_val ;
        }


        void Loop(){
            // CONVERT TXT TO ROOT FILE AND DIE
            if (doConvertWeight){ 
                MakeWeightRootfile();
                return;
            }


            // CREATE INPUT FOR PMCS
            // open input files
            cout << Form("Openning input file %s", hepfilename.Data()) << endl;
            if (doReweighting && weightfileNames.size()!=0 ) OpenWeightFiles();


            f_out->cd();
            output = new Output();
            output->Reset();


            bool finished_file = false;
            // this loops over events
            int ientry=0;
            DyresNtuple ntup( hepfilename.Data() );
            for (;ientry < ntup.GetEntries(); ientry++){
                if (debug && ientry>130) break;
                ClearEvent();
                int rc=0;
                ntup.GetEntry(ientry);
                evn = ntup.Evtnum;
                evt_wt = ntup.Evtwgt;
                if( ! fmod( Log2(evn), 1 ) ) cout<<"Processed event: "<<evn<<endl;
                vx = vy = vz = 0.0;
                // read weights
                if (doReweighting){
                    if(doSavePartonInfo){
                        // load parton info
                        throw runtime_error("Save Parton Info not implemented yet.");
                        if (weightfileNames.size() == 0) throw runtime_error("Missing weight file names.");
                    } else {
                        // load weights
                        weight_PDF.clear();
                        if (weightfileNames.size() == 0) throw runtime_error("Missing weight file names.");
                        for (size_t iw=0; iw<weight_trees.size(); iw++){
                            weight_trees[iw]->GetEntry(ientry);
                            if (weight_evtn[iw] != evn) throw runtime_error("Bad order of events");
                            double store_weight = weight_val[iw];
                            if(doWeightRatio) store_weight /= evt_wt;
                            weight_PDF.push_back(store_weight);
                        }
                        ientry++;
                    }
                }
                output->NewEvent( evn, evt_wt, 0 ,
                                  vx, vy, vz,
                                  Q2, x1, x2, flav1, flav2,
                                  weight_PDF);
                // loop particles
                bool finished_particles=false;
                //while (!finished_particles){
                for ( int ipart =0; ipart < ntup.Nump; ipart++){
                    ClearParticle();
                    id = ntup.Pid[ipart];
                    px = ntup.Ppx[ipart];
                    py = ntup.Ppy[ipart];
                    pz = ntup.Ppz[ipart];
                    E  = ntup.Ppe[ipart];
                    origin = ( ipart==0 ? 0 : 1); // only boson is unstable
                    output->AddParticlePDGID( id, px,py,pz,E,origin);
                }
                output->Fill();
                if (debug && ! fmod(Log2(evn),1) ){
                        cout<<"debug event dump: "                          ;
                        cout<< "evn("     << evn                    <<"), " ;
                        cout<< "evt_wt("  << evt_wt                 <<"), " ;
                        cout<< "v("       << vx<<", "<<vy<<", "<<vz <<"), " ;
                        cout<< "Q2("      << Q2                     <<"), " ;
                        cout<< "x12("     << x1<<", "<<x2           <<"), " ;
                        cout<< "flav12("  << flav1<<", "<<flav2     <<"), " ;
                        cout<< "pdf_wgt(" << weight_PDF[0] << ", "
                                          << weight_PDF[1] << ", "
                                          << weight_PDF[2] << ", "
                                          << weight_PDF[3] << ",...), ";
                        cout<< endl;

                }
            }
            if (doReweighting && weightfileNames.size()!=0 ) CloseWeightFiles();
        }


        void Save(){
            //if (output!=0) output->Write();
            if (output!=0) f_out=output->Tree()->GetCurrentFile();
            f_out->Write();
            f_out->Close();
        }


    private:
        TString hepfilename;
        vector<TString> weightfileNames;
        Output * output;

        //FILE * hepfile;
        TFile * f_out;

        vector<TFile *> weight_files;
        vector<TTree *> weight_trees;
        Int_t  * weight_evtn;
        Float_t * weight_val;

        // event related
        int evn;
        float evt_wt;
        float vx,vy,vz;
        float Q2, x1, x2, flav1, flav2;
        vector<float> weight_PDF;

        // particle related
        int id,origin,udk;
        float px,py,pz,E;

        bool debug;
        bool doConvertWeight;
        bool doReweighting;
        bool doWeightRatio;
        bool doSavePartonInfo;

};


void help(){
    cout << "USAGE 1 :" << endl;
    cout << " tupleMaker_DYRES output.root input.root [-f] [ weight_01.root weight_02.root ...]" << endl;
    cout << "    -f    save parton info (default off)" << endl << endl;
    cout << " Description : create input tree for pmcs from hep and weight files." << endl<<endl<<endl;


    cout << "USAGE 2 :" << endl;
    cout << " tupleMaker_DYRES -r weight.txt" << endl << endl;
    cout << " Description : create root weight file from weight text file." << endl<<endl<<endl;


    cout << "USAGE 3 :" << endl;
    cout << " tupleMaker_DYRES -c output.root" << endl << endl;
    cout << " Description : check number of entries in tree in file output.root." << endl<<endl<<endl;

    cout << "To only print this help use '-h' or '--help' as only argument." << endl;
}

/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){

    TupleMaker tm;

    //parse inputs : program output.root central.hep weight.txt ...
    if (argc < 3 ) { help(); return 1;}
    else {
        //HELP
        TString tmp = argv[1];
        if (!tmp.CompareTo("-h") || !tmp.CompareTo("--help")) {
            help();
            return 0;
        }

        if (!tmp.CompareTo("-r")){
            // USAGE 2
            tmp = argv[2];
            tm.SetConvertWeight(true);
            tm.AddWeights(tmp);
            tm.SetOutName(tmp.ReplaceAll(".txt",".root"));
        } else if (!tmp.CompareTo("-c")){
            // USAGE 3
            if (argc == 3){
                tm.GetEntries(argv[2]);
            } else {
                tm.GetEntries(argv[2],argv[3]);
            }
            return 0;
        } else {
            // USAGE 1
            tm.SetOutName(argv[1]);
            tm.SetKinematic(argv[2]);
            if (argc > 3 ){ // files and setting for weights
                int start=3;
                tm.SetDoReweighting(true);
                // check for parton info saving
                if ( !TString(argv[start]).CompareTo("-f") ){
                    tm.SetSavePartonInfo(true);
                    start++;
                }
                for (int i=start; i<argc;i++){
                    tm.AddWeights(argv[i]);
                }
            }
        }
    }

    try {
        tm.Loop();
    } catch (runtime_error& e){
        cout << "An error happend while looping:\n" << e.what() << endl;
        cout << "I will try to save it anyway, please check it carefully." << endl;
        tm.Save();
        return 2;
    }
    tm.Save();

    return 0;
}


#endif // tupleMaker_DYRES
