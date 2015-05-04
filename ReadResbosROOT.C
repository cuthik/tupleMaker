#ifndef ReadResbosROOT_C
#define ReadResbosROOT_C

/**
 * @file ReadResbosROOT.C
 * Description of this macro
 *
 * @brief Work with resbos root output files.
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-02-18
 */

#include "TString.h"
#include "TRegexp.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <vector>
#include <stdexcept>
#include <fstream>

using std::vector;
using std::runtime_error;


class ResbosRootNtuple{
    public :
        ~ResbosRootNtuple(){
            if(!hfile && hfile->IsOpen()) hfile->Close(); 
            if(!vec_wgt) delete vec_wgt;
        };
        ResbosRootNtuple(TString path="TestJob/resbos.root", TString options = "READ"){
            orig_tree = 0;
            weight_trees = 0;
            wgt_mean      = -999. ;
            if (options.CompareTo("WEIGHT",TString::kIgnoreCase)==0){
                hfile = TFile::Open(path.Data(),"READ");
                if(!hfile){
                    TString tmp =  "Can not open file: "; tmp+=path;
                    throw runtime_error(tmp.Data());
                }
                tree = (TTree *) hfile->Get("h10");
                tree -> SetBranchAddress("WTXX"  , & wgt   );
                tree -> SetBranchStatus("*",1);
            } else if (options.CompareTo("READ",TString::kIgnoreCase)==0){
                hfile = TFile::Open(path.Data(),"READ");
                if(!hfile){
                    TString tmp =  "Can not open file: "; tmp+=path;
                    throw runtime_error(tmp.Data());
                }
                tree = (TTree *) hfile->Get("h10");
                tree -> SetBranchAddress("Px_el" , & l1_Px );
                tree -> SetBranchAddress("Py_el" , & l1_Py );
                tree -> SetBranchAddress("Pz_el" , & l1_Pz );
                tree -> SetBranchAddress("E_el"  , & l1_E  );
                tree -> SetBranchAddress("Px_nu" , & l2_Px );
                tree -> SetBranchAddress("Py_nu" , & l2_Py );
                tree -> SetBranchAddress("Pz_nu" , & l2_Pz );
                tree -> SetBranchAddress("E_nu"  , & l2_E  );
                tree -> SetBranchAddress("Px_V"  , & VB_Px );
                tree -> SetBranchAddress("Py_V"  , & VB_Py );
                tree -> SetBranchAddress("Pz_V"  , & VB_Pz );
                tree -> SetBranchAddress("E_V"   , & VB_E  );
                tree -> SetBranchAddress("WT00"  , & wgt   );
                tree -> SetBranchStatus("*",1);
                // make correct weight => mean of weights should be equal to 1
                tree->Draw("WT00");
                //tree->GetHistogram()->SetDirectory(0);
                wgt_mean = tree->GetHistogram()->GetMean();
            } else if (options.CompareTo("RECREATE",TString::kIgnoreCase)==0){
                hfile = TFile::Open(path.Data(),"RECREATE");
                tree = (TTree*) new TTree("h10", "h10+weights");
            }
            ResetVars();
        }

        Long64_t GetEntries() { return tree->GetEntriesFast(); }

        Long64_t GetExpectedEntries() {
            Long64_t orig_entries = orig_tree->GetEntries();
            for (size_t i=0; i<weight_trees->size(); i++){
                if (orig_entries != weight_trees->at(i)->GetEntries()) return 0;
            }
            return orig_entries;
        }

        Long64_t GetEntry(Long64_t i=0) {
            Long64_t a = tree->GetEntry(i);
            wgt_corrected = wgt/wgt_mean;
            return  a;
        }

        Long64_t GetConnectedEntry(Long64_t ievt=0) {
            Long64_t a = orig_tree ? orig_tree->GetEntry(ievt) : 0;
            if(weight_trees) for (size_t i=0; i<weight_trees->size(); i++){
                Long64_t b = weight_trees->at(i)->GetEntry(ievt);
            }
            return  a;
        }

        void ResetVars(){
            l1_Px = -999. ;
            l1_Py = -999. ;
            l1_Pz = -999. ;
            l1_E  = -999. ;
            l2_Px = -999. ;
            l2_Py = -999. ;
            l2_Pz = -999. ;
            l2_E  = -999. ;
            VB_Px = -999. ;
            VB_Py = -999. ;
            VB_Pz = -999. ;
            VB_E  = -999. ;
            wgt   = -999. ;

            wgt_corrected = -999. ;

            TString s_infodump = "";
            if (!s_infodump.Length()) s_infodump = GetStringFromFile("config_dump.txt");
            if (!s_infodump.Length()) s_infodump =  hfile->Get("ConfigurationDump") ? ((TObjString*) hfile->Get("ConfigurationDump") )->GetString() : "";

            typeVB = GetValParameterInt(s_infodump,"V_pid");  //  24 -24  23
            type1  = GetValParameterInt(s_infodump,"el_pid"); // -13  13  13
            type2  = GetValParameterInt(s_infodump,"nu_pid"); //  14 -14 -13

            ebeam    = GetValParameterInt(s_infodump,"ebeam"); // 3.5e3;          // GeV
            kMu_mass = 105.6583715e-3; // GeV
            kEl_mass = 510.998928e-6;  // GeV
            kNu_mass = 0.;
        }

        double GetMassFromPDG(int pdg){
            switch(abs(pdg)){
                case 11: return kEl_mass; break;
                case 13: return kMu_mass; break;
                default: return kNu_mass;
            }
            return kNu_mass;
        }

        void GetLorentzVector(){
            v_VB.SetXYZT( VB_Px , VB_Py , VB_Pz , VB_E);
            v_l1.SetXYZM( l1_Px , l1_Py , l1_Pz , GetMassFromPDG(type1));
            v_l2.SetXYZM( l2_Px , l2_Py , l2_Pz , GetMassFromPDG(type2));

            // -- check if we take correct kinematics --
            // double diff1 = ( (v_l1+v_l2).Px() - v_VB.Px() ) / v_VB.Px() ;
            // double diff2 = ( (v_l1+v_l2).Py() - v_VB.Py() ) / v_VB.Py() ;
            // double diff3 = ( (v_l1+v_l2).Pz() - v_VB.Pz() ) / v_VB.Pz() ;
            // double diff4 = ( (v_l1+v_l2).M () - v_VB.M () ) / v_VB.M () ;
            // if ( fabs(diff1)>1e-2 || 
            //      fabs(diff2)>1e-2 ||
            //      fabs(diff3)>1e-2 ||
            //      fabs(diff4)>1e-2 
            //    )
            //     printf ( " kin rel diff:  momentum % e % e % e   Mass % e \n", 
            //             diff1,diff2,diff3,diff4
            //            );
        }

        void ConnectTree(ResbosRootNtuple * _orig_tree){
            // set pointer
            orig_tree=_orig_tree;
            // set brances
            if (hfile->IsWritable() ) {
                tree -> Branch("Px_el" , & orig_tree->l1_Px );
                tree -> Branch("Py_el" , & orig_tree->l1_Py );
                tree -> Branch("Pz_el" , & orig_tree->l1_Pz );
                tree -> Branch("E_el"  , & orig_tree->l1_E  );
                tree -> Branch("Px_nu" , & orig_tree->l2_Px );
                tree -> Branch("Py_nu" , & orig_tree->l2_Py );
                tree -> Branch("Pz_nu" , & orig_tree->l2_Pz );
                tree -> Branch("E_nu"  , & orig_tree->l2_E  );
                tree -> Branch("Px_V"  , & orig_tree->VB_Px );
                tree -> Branch("Py_V"  , & orig_tree->VB_Py );
                tree -> Branch("Pz_V"  , & orig_tree->VB_Pz );
                tree -> Branch("E_V"   , & orig_tree->VB_E  );
                tree -> Branch("WT00"  , & orig_tree->wgt   );
            }
        }

        void ConnectWeights(vector<ResbosRootNtuple*> * _wgt_trees){
            // set pointer
            weight_trees = _wgt_trees;
            vec_wgt = new vector<float>();
            // add weight branch
            if (hfile->IsWritable() ) {
                tree -> Branch("WTXX" , "vector<float>" , & vec_wgt   );
            }
        }

        void UpdateWeights(){
            vec_wgt->clear();
            for (size_t i=0; i<weight_trees->size(); i++){
                vec_wgt->push_back(weight_trees->at(i)->wgt/orig_tree->wgt);
            }
        }

        void Fill(){
            UpdateWeights();
            tree->Fill();
        }

        void Save(){
            hfile->cd();
            tree->Write();
            hfile->Write();
            hfile->Close();
        }

        // data members
        TTree * tree;
        TFile * hfile;

        Float_t l1_Px ;
        Float_t l1_Py ;
        Float_t l1_Pz ;
        Float_t l1_E  ;
        Float_t l2_Px ;
        Float_t l2_Py ;
        Float_t l2_Pz ;
        Float_t l2_E  ;
        Float_t VB_Px ;
        Float_t VB_Py ;
        Float_t VB_Pz ;
        Float_t VB_E  ;
        Float_t wgt   ;
        vector<float> *vec_wgt;

        int typeVB;
        int type1;
        int type2;
        Float_t ebeam;
        Float_t kMu_mass;
        Float_t kEl_mass;
        Float_t kNu_mass;
        Float_t wgt_mean   ;
        Float_t wgt_corrected   ;
        TLorentzVector v_VB;
        TLorentzVector v_l1;
        TLorentzVector v_l2;

        // reference pointers
        ResbosRootNtuple * orig_tree;
        vector<ResbosRootNtuple*> * weight_trees;

        // static functions
        static TString GetStringFromFile(TString file_path);
        static TString GetValParameterTString (TString config_dump, TString par_name);
        static int     GetValParameterInt     (TString config_dump, TString par_name);
};

TString ResbosRootNtuple::GetStringFromFile(TString file_path){
    TString output = "";
    std::ifstream ifs;
    ifs.open (file_path.Data(), std::ifstream::in);
    output.ReadFile(ifs);
    ifs.close();
    if(output.BeginsWith("root")){
        output="";
        TFile *f = TFile::Open(file_path,"READ"); // try to open it, maybe it's root file
        TObjString * sobj = (TObjString*) f->Get("ConfigurationDump");
        if(sobj) output = sobj->GetString();
    }
    return output;
}

// template<class T>
// T ResbosRootNtuple::GetValParameter<T>(TString config_dump, TString par_name){
//     return T();
// }

TString ResbosRootNtuple::GetValParameterTString(TString config_dump, TString par_name){
    TRegexp exp = par_name+":.*";
    TString line  =  config_dump(exp);
    exp = par_name+": *";
    TString remove = line(exp);
    return line.ReplaceAll(remove,"");
}

int ResbosRootNtuple::GetValParameterInt(TString config_dump, TString par_name){
    TString str  = GetValParameterTString(config_dump,par_name);
    return str.Atoi();
}





/**
 * The description of function.
 * 
 * @param _inArg Description of param
 * @return Description of return.
 */
void ReadResbosROOT(){

    return;
}

#endif // ReadResbosROOT_C
