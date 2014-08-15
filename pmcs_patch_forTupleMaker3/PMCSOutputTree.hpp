#ifndef INC_PMCSOUTPUTTREE
#define INC_PMCSOUTPUTTREE


#include "TTree.h"
#include "TString.h"

#include "PMCSEMObj.hpp"
#include "PMCSEvent.hpp"
#include "PMCSWCand.hpp"

#include <vector>

using std::vector;


class PMCSOutputTree {
    public:


        PMCSOutputTree();
        ~PMCSOutputTree() {;}

        void Init();
        void Clear();
        void Fill();
        void Save(TString filename);

        void SetWenu(PMCSWCand * wcand, double scalarEt_All);
        void SetTruthWenu(PMCSEvent *pmcsevent, vector<PMCSEMObj> * emobj_gen, PMCSWCand * wcand);
        void SetPDFWeights(double * pdfarray, size_t n_pdf);
        void SetMassWeights(vector<double> * in);


    private:
        TTree * outTree;

        int    truth_evtn    ;

        double truth_VB_mass ;
        double truth_VB_pt   ;
        double truth_VB_phi  ;
        double truth_VB_eta  ;

        double truth_el_pt   ;
        double truth_el_phi  ;
        double truth_el_eta  ;

        double truth_nu_pt   ;
        double truth_nu_phi  ;
        double truth_nu_eta  ;

        double VB_mt         ;

        double MET_et        ;
        double MET_sumet     ;

        double el_pt         ;
        double el_phi        ;


        double hr_et         ;
        double hr_phi        ;
        double hr_prll       ;
        double hr_perp       ;


        double weight_evt    ;
        vector<double> *weights_mass ;
        vector<double> *weights_PDF  ;


};



#endif
