/**
 * @file ReadGeneNTUP.C
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-03-11
 */

#include "TFile.h"
#include "TChain.h"
//#include "THbookFile.h"
//#include "THbookTree.h"
#include "TLorentzVector.h"


// #define UseHBOOK 1

class GeneNtuple{
    public :
        ~GeneNtuple(){ if(hfile && hfile->IsOpen()) hfile->Close(); };
        GeneNtuple(TString path="TestJob/dyres.hbook", int lrecl=1024) {
            tree    = 0 ;
            hfile   = 0 ;
            Evtnum  = 0 ;
            Evtwgt  = 0 ;
            Nump    = 0 ;
            ntype = DyresRoot;
            reset_4v();
#ifdef UseHBOOK
            hfile = new THbookFile(path.Data(),lrecl);
            tree = (THbookTree*) hfile->Get(66);
            //tree = (TChain*) hfile->ConvertCWN(66);
#else
            hfile = TFile::Open(path.Data());
            switch(ntype){
                case DyresRoot :
                    tree = (TChain *) hfile->Get("h66");
                    break;
                case ResbosRoot :
                    tree = (TChain *) hfile->Get("h10");
                    break;
                default :
                    throw runtime_error("Uknown type generator ntuple.");
            }
#endif
            hfile->ls();
            tree -> Print();
            switch(ntype){
                case DyresRoot:
                    tree -> SetBranchAddress("Evtnum", &Evtnum , &b_Evtnum );
                    tree -> SetBranchAddress("Evtwgt", &Evtwgt , &b_Evtwgt );
                    tree -> SetBranchAddress("Nump"  , &Nump   , &b_Nump   );
                    tree -> SetBranchAddress("Pid"   ,  Pid    , &b_Pid    );
                    tree -> SetBranchAddress("Ppx"   ,  Ppx    , &b_Ppx    );
                    tree -> SetBranchAddress("Ppy"   ,  Ppy    , &b_Ppy    );
                    tree -> SetBranchAddress("Ppz"   ,  Ppz    , &b_Ppz    );
                    tree -> SetBranchAddress("Ppe"   ,  Ppe    , &b_Ppe    );
                    break;
                case ResboRoot:
                    tree -> SetBranchAddress("Px_el" , & l1_Ppx[1]  );
                    tree -> SetBranchAddress("Py_el" , & l1_Ppy[1]  );
                    tree -> SetBranchAddress("Pz_el" , & l1_Ppz[1]  );
                    tree -> SetBranchAddress("E_el"  , & l1_Ppe[1]   );
                    tree -> SetBranchAddress("Px_nu" , & l2_Ppx[2]  );
                    tree -> SetBranchAddress("Py_nu" , & l2_Ppy[2]  );
                    tree -> SetBranchAddress("Pz_nu" , & l2_Ppz[2]  );
                    tree -> SetBranchAddress("E_nu"  , & l2_Ppe[2]   );
                    tree -> SetBranchAddress("Px_V"  , & VB_Ppx[0]  );
                    tree -> SetBranchAddress("Py_V"  , & VB_Ppy[0]  );
                    tree -> SetBranchAddress("Pz_V"  , & VB_Ppz[0]  );
                    tree -> SetBranchAddress("E_V"   , & VB_Ppe[0]   );
                    tree -> SetBranchAddress("WT00"  , & Evtwgt );
                    tree -> SetBranchStatus("*",1);
                    break;
                default :
                    throw runtime_error("Uknown type generator ntuple.");
            }
            tree -> SetBranchStatus("*",1);
        }

        void reset_4v(){
            for (size_t ip = 0; ip < Nump; ip++){
                PP[ip].SetPxPyPzE( 0,0,0,0 );
            }
        }

        Long64_t GetEntries() { return tree->GetEntriesFast(); }
        Long64_t GetEntry(Long64_t i=0) {
            //tree->Print();
            Int_t rc = tree->GetEntry(i);
            if (ntype==ResboRoot) Evtnum = i;
            //Long64_t rc = tree->LoadTree(i);
            //Int_t rc = b_Nump->GetEntry(i);
            for (size_t ip = 0; ip < Nump; ip++){
                PP[ip].SetPxPyPzE( Ppx[ip], Ppy[ip], Ppz[ip], Ppe[ip] );
            }

            return rc;
        }


        enum GeneNtupType { DyresRoot=0, ResbosRoot };
        GeneNtupType ntype;

#ifdef UseHBOOK
        THbookTree * tree;
        THbookFile * hfile;
#else
        TChain * tree;
        TFile * hfile;
#endif
        // Declaration of leaf types
        // DYRES
        //struct {
        Int_t           Evtnum;
        Float_t         Evtwgt;
        Int_t           Nump;
        Int_t           Pid[10];   //[Nump]
        Float_t         Ppx[10];   //[Nump]
        Float_t         Ppy[10];   //[Nump]
        Float_t         Ppz[10];   //[Nump]
        Float_t         Ppe[10];   //[Nump]
        //} evt;
        // branches
        //TBranch        *b_Evtnum;   //!
        //TBranch        *b_Evtwgt;   //!
        //TBranch        *b_Nump;   //!
        //TBranch        *b_Pid;   //!
        //TBranch        *b_Ppx;   //!
        //TBranch        *b_Ppy;   //!
        //TBranch        *b_Ppz;   //!
        //TBranch        *b_Ppe;   //!
        // physical fourvectors
        Int_t           Npho;
        TLorentzVector  PP[10];
        //TLorentzVector & VB   = PP[0];
        //TLorentzVector & lep1 = PP[1];
        //TLorentzVector & lep2 = PP[2];
};


/**
 * The description of function.
 * 
 * @param _inArg Description of param
 * @return Description of return.
 */
void ReadGeneNTUP(){

    return;
}
