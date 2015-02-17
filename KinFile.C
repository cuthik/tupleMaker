/**
 * @file KinFile.C
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-01-30
 */

const double twopi = 2*TMath:Pi();

class THnD_KIN : public THnD {
    public:
    ~THnD_KIN();
    THnD_KIN(TString name, TSting title =""):
        kin_dim     (4),
        kin_nbins ({2 ,2 ,2 ,2 }),
        THnD(name,title,kindim,dummy_nbin,0,0)
        {
            if(title.Length()==0) SetTitle(name.Data());
            // qt
            double qt_bins[] = {0.,100.,7000.};
            SetBinEdges(0,qt_bins);
            fAxes[0]->SetName("q_T");
            // y
            double y_bins[] = {-10.,0.,10.};
            SetBinEdges(1,y_bins);
            fAxes[1]->SetName("y");
            // theta_CS
            double th_bins[] = {-twopi,0.,twopi};
            SetBinEdges(2,th_bins);
            fAxes[2]->SetName("theta_CS");
            // phi_CS
            double ph_bins[] = {-twopi,0.,twopi};
            SetBinEdges(3,ph_bins);
            fAxes[3]->SetName("phi_CS");
        }
    protected:
        int    kin_dim      ; // (4),
        double kin_nbins[4] ; // ({2 ,2 ,2 ,2 }),

}


class TKinFile {
    TKinFile(TString name, TString opt="READ"){
        const char * c_opt=opt.Data();
        file = TFile::Open(name,opt);
        if ( !strcasecmp(c_opt,"READ") || !strcasecmp(c_opt,"CACHEREAD") ){
            histND = (THnD*)  file->Get("kin_hist");
            normND = (THnD*)  file->Get("norm_hist");
            tree   = (TTree*) file->Get("kin_tree");
            tree->SetBranchAddress("qt"       , &qt       );
            tree->SetBranchAddress("y"        , &y        );
            tree->SetBranchAddress("theta_CS" , &theta_CS );
            tree->SetBranchAddress("phi_CS"   , &phi_CS   );
        } else {
            histND = new THnD("kin_hist") ;
            normND = new THnD("norm_hist") ;
            tree   = new TTree("kin_tree","kin_tree");
            tree->SetBranchAddress("qt/D"       , &qt       );
            tree->SetBranchAddress("y/D"        , &y        );
            tree->SetBranchAddress("theta_CS/D" , &theta_CS );
            tree->SetBranchAddress("phi_CS/D"   , &phi_CS   );
        }
    }
    ~TKinFile();

    TFile * file;
    TTree * tree;
    THnD * histND;
    THnD * normND;

    // tree var
    double weight;
    double qt;
    double y;
    double theta_CS;
    double phi_CS;
}

/**
 * The description of function.
 * 
 * @param _inArg Description of param
 * @return Description of return.
 */
void KinFile(){

    return;
}
