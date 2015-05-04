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
//#include "KinFile.C"
#include "KinFile.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TMath.h"
//#include "TString.h"

#include <iostream>
#include <iomanip>
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
            hepfile           (0),
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
            if ( f_out && f_out->IsOpen()) f_out->Close();
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

        void GetEntryWeights(int &ientry){
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
                    if (weight_evtn[iw] != evn) throw runtime_error( Form("Bad order of events: pmcs evt %d vs weight evt %i (of wgt %lu)",evn,weight_evtn[iw],iw));
                    double store_weight = weight_val[iw];
                    if(doWeightRatio) store_weight /= evt_wt;
                    weight_PDF.push_back(store_weight);
                }
                ientry++;
            }
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


        void Loop(){
            // CONVERT TXT TO ROOT FILE AND DIE
            if (doConvertWeight){ 
                MakeWeightRootfile();
                return;
            }


            // CREATE INPUT FOR PMCS
            // open input files
            cout << Form("Openning input file %s", hepfilename.Data()) << endl;
            FILE * hepfile;
            ifstream hepstream;
            bool useCFile = false;
            if (useCFile){
                hepfile = fopen(hepfilename.Data(),"r");
                if (hepfile==NULL){ perror ( "Opening failed: "); }
            } else { //use the ifstream
                hepstream.open(hepfilename.Data());
                if (!hepstream.is_open() || !hepstream.good()){ throw runtime_error(Form("Can not open file %s",hepfilename.Data()));}

            }
                //throw runtime_error(Form("Can not open input file `%s`, error number is %i", hepfilename.Data(),errno)); } 
            if (doReweighting && weightfileNames.size()!=0 ) OpenWeightFiles();


            f_out->cd();
            output = new Output();
            output->Reset();


            bool finished_file = false;
            // this loops over events
            int ientry=0;
            while( ! finished_file ) {
                if (debug && ientry>130) break;
                ClearEvent();
                int rc=0;
                // first line (event number and weight)
                if (useCFile){
                    rc = fscanf(hepfile,"%i %f", &evn, &evt_wt);
                    if (rc!=2) throw runtime_error("Event line not found.");
                } else { // use the ifstream
                    hepstream >> evn >> evt_wt;
                    if (!hepstream.good()) throw runtime_error("Event line not found.");
                }
                printEvent(evn);
                if( evn == 0 ) {
                    finished_file = true;
                    continue;
                } 

                // second line (beam position)
                if (useCFile){
                    rc=fscanf(hepfile,"%f %f %f", &vx, &vy, &vz);
                    if (rc!=3) throw runtime_error("Vertex line not found.");
                } else { // use the ifstream
                    hepstream >> vx >> vy >> vz;
                    if (!hepstream.good()) throw runtime_error("Vertex line not found.");
                }

                // read weights
                if (doReweighting) GetEntryWeights(ientry);

                output->NewEvent( evn, evt_wt, 0 ,
                                  vx, vy, vz,
                                  Q2, x1, x2, flav1, flav2,
                                  weight_PDF);


                // loop particles
                bool finished_particles=false;
                while (!finished_particles){
                    ClearParticle();

                    // read isajet id
                    if (useCFile){
                        rc=fscanf(hepfile,"%i",&id);
                        if (rc!=1) throw runtime_error("Missing isajet id.");
                    } else { // use the ifstream
                        hepstream >> id;
                        if (!hepstream.good()) throw runtime_error("Missing isajet id.");
                    }

                    if(id==0){
                        finished_particles=true;
                        continue;
                    }

                    // read kinematics
                    if (useCFile){
                        rc=fscanf(hepfile,"%f %f %f %f %i %i",&px,&py,&pz,&E,&origin,&udk);
                        if (rc!=6) throw runtime_error("Missing particle kinematics.");
                    } else { // use the ifstream
                        hepstream >> px >> py >> pz >> E >> origin >> udk;
                        if (!hepstream.good()) throw runtime_error("Missing particle && ematics.");
                    }

                    output->AddParticle( id, px,py,pz,E,origin);
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

            //close input files
                if (useCFile){
                    fclose(hepfile);
                } else { // use the ifstream
                    hepstream.close();
                }
            if (doReweighting && weightfileNames.size()!=0 ) CloseWeightFiles();
        }


        void Save(){
            //if (output!=0) output->Write();
            if (output!=0) f_out=output->Tree()->GetCurrentFile();
            f_out->Write();
            f_out->Close();
            f_out=0;
        }

        // KINEMATIC REWEIGHT PART;
        TString this_pmcs_path;
        TString this_kin_path;
        TString new_kin_path;
        TString new_pmcs_path;

        void load_VBleplep(Output *p, TLorentzVector &VB, TLorentzVector &l_part, TLorentzVector &l_anti){
            // find indeces
            size_t iVB     = 9999;
            size_t il_part = 9999;
            size_t il_anti = 9999;
            // get event
            int VB_PID = 0;
            for (size_t i = 0; i<p->_ana.npart; i++){ // find vector boson type
                int pid = p->_ana.ppid[i];
                if (pid ==  23 ){ VB_PID =  pid; iVB = i; break; }
                if (pid ==  24 ){ VB_PID =  pid; iVB = i; break; }
                if (pid == -24 ){ VB_PID =  pid; iVB = i; break; }
            }
            for (size_t i = 0; i<p->_ana.npart; i++){
                int pid = p->_ana.ppid[i];
                if (VB_PID ==  24){ // W+
                    if (il_part == 9999 && pid== 12) il_part=i; // nu
                    if (il_anti == 9999 && pid==-11) il_anti=i; // e+
                }
                if (VB_PID == -24){ // W-
                    if (il_part == 9999 && pid== 11) il_part=i; // e-
                    if (il_anti == 9999 && pid==-12) il_anti=i; // nubar
                }
                if (VB_PID ==  23){ // Z
                    if (il_part == 9999 && pid== 11) il_part=i; // e-
                    if (il_anti == 9999 && pid==-11) il_anti=i; // e+
                }
                if (il_part!=9999 && il_anti!=9999) break;
            }
            // set kinematics
            VB.SetPxPyPzE(
                    p->_ana.ppx[iVB],
                    p->_ana.ppy[iVB],
                    p->_ana.ppz[iVB],
                    p->_ana.pE [iVB]
                    );
            l_part.SetPxPyPzE(
                    p->_ana.ppx[il_part],
                    p->_ana.ppy[il_part],
                    p->_ana.ppz[il_part],
                    p->_ana.pE [il_part]
                    );
            l_anti.SetPxPyPzE(
                    p->_ana.ppx[il_anti],
                    p->_ana.ppy[il_anti],
                    p->_ana.ppz[il_anti],
                    p->_ana.pE [il_anti]
                    );
            return;
        }

        void CreateKinFile(){
            this_pmcs = new Output(this_pmcs_path);
            this_kin = new TKinFile(this_kin_path,"RECREATE");
            size_t Nevents = this_pmcs->GetEntries();
            //Nevents = int(1e5);
            for (size_t ievt = 0; ievt < Nevents; ievt++){
                printEvent(ievt,Nevents);
                this_pmcs -> GetEntry(ievt);
                TLorentzVector VB,l_part,l_anti;
                load_VBleplep(this_pmcs, VB, l_part, l_anti );
                this_kin -> Fill(VB, l_part, l_anti, this_pmcs->_ana.evwt[0] );
            }
            this_kin->Save();
            return;
        }

        void ReweightPMCSbyKin(){
            this_pmcs = new Output(this_pmcs_path);
            SetOutName(new_pmcs_path);
            new_pmcs  = new Output();
            new_pmcs ->Reset();

            this_kin = new TKinFile(this_kin_path ,"READ");
            new_kin  = new TKinFile(new_kin_path  ,"READ");

            //this_kin -> use4Dhist = false;
            //new_kin  -> use4Dhist = false;

            // loop old add weight
            int ientry=0;
            int Nevents = this_pmcs->GetEntries();
            //Nevents = int(1e5);
            for ( int i=0 ; i< Nevents; i++){
                this_pmcs->GetEntry(i);
                evn = i+1;
                printEvent(evn,Nevents);
                // load kinematics
                TLorentzVector VB,l_part,l_anti;
                load_VBleplep(this_pmcs, VB, l_part, l_anti );
                // get new weight
                double new_wt = this_kin->GetNewWeight(new_kin, VB, l_part, l_anti);
                // fill new weight
                new_pmcs->NewEventNewWeight( this_pmcs, new_wt * this_pmcs->_ana.evwt[0] );
                new_pmcs->AddParticles( this_pmcs );
                new_pmcs->Fill();
            }
            // save
            Save();
            delete this_kin;
            delete new_kin;
            return;
        }

        void AddWeightsToPMCS(){
            // open files / create new
            this_pmcs = new Output(this_pmcs_path);
            SetOutName(new_pmcs_path);
            new_pmcs  = new Output();
            new_pmcs->Reset();
            OpenWeightFiles();
            // loop old add weight
            int ientry=0;
            int Nevents = this_pmcs->GetEntries();
            for ( int i=0 ; i< Nevents; i++){
                this_pmcs->GetEntry(i);
                evn = i+1;
                printEvent(evn,Nevents);
                GetEntryWeights(ientry);
                new_pmcs->NewEvent     ( this_pmcs, weight_PDF );
                new_pmcs->AddParticles ( this_pmcs             );
                new_pmcs->Fill();
            }
            // save
            CloseWeightFiles();
            Save();
            return;
        }







    private:
        TString hepfilename;
        vector<TString> weightfileNames;
        Output * output;
        Output * this_pmcs;
        Output * new_pmcs;
        TKinFile * this_kin;
        TKinFile * new_kin;

        FILE * hepfile;
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


        TStopwatch swatch;
        bool debug;
        bool doConvertWeight;
        bool doReweighting;
        bool doWeightRatio;
        bool doSavePartonInfo;

};


void help(){
    cout << "USAGE 1 :" << endl;
    cout << " tupleMaker3 output.root central.help [-f] [ weight_01.root weight_02.root ...]" << endl;
    cout << "    -f    save parton info (default off)" << endl << endl;
    cout << " Description : create input tree for pmcs from hep and weight files." << endl<<endl<<endl;

    cout << "USAGE 2 :" << endl;
    cout << " tupleMaker3 -r weight.txt" << endl << endl;
    cout << " Description : create root weight file from weight text file." << endl<<endl<<endl;

    cout << "USAGE 3 :" << endl;
    cout << " tupleMaker3 -c output.root" << endl << endl;
    cout << " Description : check number of entries in tree in file output.root." << endl<<endl<<endl;

    cout << "USAGE 4 :" << endl;
    cout << " tupleMaker3 -k pmcs_in.root new_kin.root" << endl << endl;
    cout << " Description : will create new_kin.root containing TTree and nDim histogram for kinematic reweighting." << endl<<endl<<endl;

    cout << "USAGE 5 :" << endl;
    cout << " tupleMaker3 -k this_pmcs_in.root this_kin.root new_kin.root new_pmcs_in.root" << endl << endl;
    cout << " Description : create new_pmcs_in.root with weights calculated from this_kin.root and new_kin.root." << endl<<endl<<endl;

    cout << "USAGE 6 :" << endl;
    cout << " tupleMaker3 -a output.root input.root weight_01.root [...]" << endl << endl;
    cout << " Description : create output.root with additional weights from (at least one) weight file" << endl<<endl<<endl;

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
        } else if (!tmp.CompareTo("-k") ) {
            if ( argc == 4 ){
                // USAGE 4
                tm.this_pmcs_path = argv[2];
                tm.this_kin_path  = argv[3];
                tm.CreateKinFile();
            } else if ( argc == 6 ) {
                // USAGE 5
                tm.this_pmcs_path = argv[2];
                tm.this_kin_path  = argv[3];
                tm.new_kin_path   = argv[4];
                tm.new_pmcs_path  = argv[5];
                tm.ReweightPMCSbyKin();
            } else {
                help();
            }
            return 0;
        } else if (!tmp.CompareTo("-a") ) {
            // USAGE 6
            if ( argc > 4 ){
                tm.new_pmcs_path  = argv[2];
                tm.this_pmcs_path = argv[3];
                for (int i=4; i<argc;i++){
                    tm.AddWeights(argv[i]);
                }
                try {
                    tm.AddWeightsToPMCS();
                } catch (runtime_error& e){
                    cout << "An error happend while looping:\n" << e.what() << endl;
                    cout << "I will try to save it anyway, please check it carefully." << endl;
                    tm.Save();
                }
            } else {
                help();
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


#endif // tupleMaker3
