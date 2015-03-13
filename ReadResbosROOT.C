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
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <vector>

using std::vector;


class ResbosRootNtuple{
    public :
        ~ResbosRootNtuple(){
            if(!hfile && hfile->IsOpen()) hfile->Close(); 
            if(!vec_wgt) delete vec_wgt;
        };
        ResbosRootNtuple(TString path="TestJob/resbos.root", TString options = "READ"){
            ResetVars();
            orig_tree = 0;
            weight_trees = 0;
            if (options.CompareTo("WEIGHT",TString::kIgnoreCase)==0){
                hfile = TFile::Open(path.Data(),"READ");
                printf(" weight constructor\n");
                tree = (TTree *) hfile->Get("h10");
                tree -> SetBranchAddress("WTXX"  , & wgt   );
                tree -> SetBranchStatus("*",1);
            } else if (options.CompareTo("READ",TString::kIgnoreCase)==0){
                hfile = TFile::Open(path.Data(),"READ");
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
                tree->GetHistogram()->SetDirectory(0);
                wgt_mean = tree->GetHistogram()->GetMean();
            } else if (options.CompareTo("RECREATE",TString::kIgnoreCase)==0){
                hfile = TFile::Open(path.Data(),"RECREATE");
                tree = (TTree*) new TTree("h10", "h10+weights");
            }
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
            Long64_t a = orig_tree->tree->GetEntry(ievt);
            for (size_t i=0; i<weight_trees->size(); i++){
                Long64_t b = weight_trees->at(i)->tree->GetEntry(ievt);
                printf("entry a %Lu b %Lu\n",a,b);
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

            wgt_mean      = -999. ;
            wgt_corrected = -999. ;

            type1 = -13;
            type2 =  13;

            ebeam    = 3.5e3;          // GeV
            kMu_mass = 105.6583715e-3; // GeV
            kEl_mass = 510.998928e-6;  // GeV
        }

        void GetLorentzVector(){
            v_l1.SetXYZM( l1_Px , l1_Py , l1_Pz , kMu_mass);
            v_l2.SetXYZM( l2_Px , l2_Py , l2_Pz , kMu_mass);
        }

        void ConnectTree(ResbosRootNtuple * _orig_tree){
            // set pointer
            orig_tree=_orig_tree;
            // set brances
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

        void ConnectWeights(vector<ResbosRootNtuple*> * _wgt_trees){
            // set pointer
            weight_trees = _wgt_trees;
            vec_wgt = new vector<float>();
            // add weight branch
            tree -> Branch("WTXX" , "vector<float>" , & vec_wgt   );
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
            printf(" execute save\n");
            tree->Write();
            printf(" end save\n");
            hfile->Write();
            hfile->Close();
        }

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

        int type1;
        int type2;
        Float_t ebeam;
        Float_t kMu_mass;
        Float_t kEl_mass;
        Float_t wgt_mean   ;
        Float_t wgt_corrected   ;
        TLorentzVector v_l1;
        TLorentzVector v_l2;

        // reference pointers
        ResbosRootNtuple * orig_tree;
        vector<ResbosRootNtuple*> * weight_trees;
};

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
