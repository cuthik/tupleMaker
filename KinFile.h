#ifndef _KinFile_H
#define _KinFile_H

/**
 * @file KinFile.h
 *
 * @brief A header file for class KinFile.
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-03-05
 */

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "THn.h"
#include "TH1D.h"

#include <vector>
#include <cstdio>

#include "TLVUtils.h"

using std::vector;

const double twopi = 2*TMath::Pi();
const double pi = TMath::Pi();

#define KDIM 4

class THnD_KIN : public THnD {
    public:
        ~THnD_KIN();
        THnD_KIN();
        THnD_KIN(TString name, TString title ="");

        void Fill();
        double GetBinContent();
        double ComputeIntegral();
        double Integral();
        void updateBin();
        bool isOutlayer();

        void CopyPointers(const THnD_KIN* from);

        // value helpers
        int      bin             ; // current bin index
        double   x[KDIM]         ; // array of current values (needed for passing to THn)
        Double_t *weight_ptr     ; //! pointer to event weight (need to copy it by hand)
        Double_t *val_ptrs[KDIM] ; //! array of pointers to variables (need to copy it by hand)


        // initial definitions
        static const int    kin_dim       ; ///< Number of variables
        static const int    kin_nbins[]   ; ///< Number of bins per each 
        static const char * var_name_1    ; ///< Title of axis of variable 1
        static const char * var_name_2    ; ///< Title of axis of variable 2
        static const char * var_name_3    ; ///< Title of axis of variable 3
        static const char * var_name_4    ; ///< Title of axis of variable 4
        static const double bin_edges_1[] ; ///< Bin edges of variable 1
        static const double bin_edges_2[] ; ///< Bin edges of variable 2
        static const double bin_edges_3[] ; ///< Bin edges of variable 3
        static const double bin_edges_4[] ; ///< Bin edges of variable 4

        ClassDef(THnD_KIN,1);
};


class TKinFile {
    public :
        TKinFile(TString name, TString opt="READ");
        ~TKinFile();;

        void setValues(const TLorentzVector &VB, const TLorentzVector &l1, const TLorentzVector &l2);

        void Fill(const TLorentzVector &VB, const TLorentzVector &l1, const TLorentzVector &l2, double wt);

        void normalizeKin();

        double getValue(const double qt, const double y, const double cos_theta_CS, const double phi_CS );

        double GetNewWeight(TKinFile * new_kin, const TLorentzVector &VB, const TLorentzVector &l_part, const TLorentzVector &l_anti );

        void Save();

        TFile * file;
        TTree * tree;
        THnD_KIN * histND;
        THnD_KIN * normND;

        double ebeam;
        bool use4Dhist;
        // tree var
        double weight;
        double qt;
        double y;
        double cos_theta_CS;
        double phi_CS;
        vector<double> aimom_CS;

    protected :
        TH1D * hw;
        TH1D * h1;
        TH1D * h2;
        TH1D * h3;
        TH1D * h4;
};

#endif //_KinFile_H
