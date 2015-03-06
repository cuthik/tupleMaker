#ifndef _KinFile_CXX
#define _KinFile_CXX

/**
 * @file KinFile.cxx
 *
 * @brief A source file for class KinFile.
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-03-05
 */

#include "KinFile.h"

ClassImp(THnD_KIN)

const int THnD_KIN::kin_dim = KDIM;
const int THnD_KIN::kin_nbins[] = {16,15,10,20} ; // {4,4,4,4};
const char * THnD_KIN::var_name_1 = "q_{T}"        ; const double THnD_KIN::bin_edges_1[] =  {0   , 3   , 6   , 9   , 12  , 15 , 18 , 21 , 24 , 27 , 30 , 35 , 40 , 50 , 100, 200, 1000 } ;// {0.,5.,10.,100.,7000.}   ;
const char * THnD_KIN::var_name_2 = "y"            ; const double THnD_KIN::bin_edges_2[] =  {0.  , 0.2 , 0.4 , 0.6 , 0.8 , 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 5.0 }       ;// {-10.,-5.,0.,5.,10.}    ;
const char * THnD_KIN::var_name_3 = "cos_theta_CS" ; const double THnD_KIN::bin_edges_3[] =  {-1.0, -0.8, -0.6, -0.4, -0.2, 0  , 0.2, 0.4, 0.6, 0.8, 1.0 }                                ;// {-1.,0.,1.}       ;
const char * THnD_KIN::var_name_4 = "phi_CS"       ; const double THnD_KIN::bin_edges_4[] =  { -pi , -2.7, -2.4, -2.1, -1.8, -1.5, -1.2 , -0.9 , -0.6 , -0.3 , 0   , 0.3 , 0.6 , 0.9 , 1.2 , 1.5, 1.8, 2.1, 2.4, 2.7, pi}                                 ;// {-twopi,-pi,0.,pi,twopi} ;
//const char * THnD_KIN::var_name_4 = "phi_CS"       ; const double THnD_KIN::bin_edges_4[] =  {0   , 0.3 , 0.6 , 0.9 , 1.2 , 1.5, 1.8, 2.1, 2.4, 2.7, pi }                                 ;// {-twopi,-pi,0.,pi,twopi} ; // 10

THnD_KIN::~THnD_KIN(){
};

THnD_KIN::THnD_KIN() : THnD() {
    bin = -1;
    for ( size_t i=0; i < KDIM; i++){
        val_ptrs[i] = 0;
        x[i] = 0;
    }
    weight_ptr = 0;
}


THnD_KIN::THnD_KIN(TString name, TString title) : 
    THnD(name.Data(),title.Data(),kin_dim,kin_nbins,  NULL , NULL ) {
    THnD::Sumw2();
    if(title.Length()==0) SetTitle(name.Data());
    ((TAxis *)fAxes[0]) ->SetTitle(var_name_1); SetBinEdges(0,bin_edges_1);
    ((TAxis *)fAxes[1]) ->SetTitle(var_name_2); SetBinEdges(1,bin_edges_2);
    ((TAxis *)fAxes[2]) ->SetTitle(var_name_3); SetBinEdges(2,bin_edges_3);
    ((TAxis *)fAxes[3]) ->SetTitle(var_name_4); SetBinEdges(3,bin_edges_4);

    bin = -1;
    for ( size_t i=0; i < KDIM; i++){
        val_ptrs[i] = 0;
        x[i] = 0;
    }
    weight_ptr = 0;
}

void THnD_KIN::Fill(){
    updateBin();
    THnD::Fill(x,*weight_ptr);
}

bool THnD_KIN::isOutlayer(){ // underflow/overflow
    for (size_t i=0; i< fNdimensions; i++){
        if ( x[i] < GetAxis(i)->GetXmin() )  return true;
        if ( x[i] > GetAxis(i)->GetXmax() )  return true;
    }
    return false;
}

double THnD_KIN::GetBinContent(){
    updateBin();
    if (isOutlayer()) return 0;
    return THnD::GetBinContent(bin);
}

double THnD_KIN::Integral(){
    // This is another very funny ROOT story:
    // The ComputeIntegral function will compute correctly everyting and than..
    // I don't know why is normalizing integral array. integral array is
    // defined
    //
    //     fIntegral[i]  = Sum_j=1^i bin_j
    //
    // so last item [GetBins()] contain integral and first [0] contains 0.
    // but after normalization the last bin contain 1 and all others the fraction
    //
    //     fIntegral[i]  = Sum_j=1^i bin_j / Sum_k=1^N bin_j
    //
    // so I decided to calculate back the integral from second item
    //
    //     fIntegral[1]  = bin_1 / Integral    =>     Integral = bin_1/fIntegral[1]
    //
    // this would be true if first bin has some value. The savest value is to
    // find first nonzero bin. Search code adapted from ComputeIntegral.
    //
    // compute integral normaly
    //if(fIntegralStatus == kValidInt && fIntegral && fIntegral[0]) return fIntegral[0];
    if(fIntegralStatus != kValidInt) ComputeIntegral();
    if(fIntegralStatus != kValidInt) return 0;
    Int_t* coord = new Int_t[fNdimensions];
    // setup coords in way that we skipp almost all underflows
    //coord[fNdimensions-1] = 0; // only last one
    Long64_t i = 0;
    Double_t VAL = 0;
    THnIter iter(this);
    while ((i = iter.Next(coord)) >= 0) {
        VAL = THnD::GetBinContent(i);
        // check whether the bin is regular
        bool regularBin = true;
        for (Int_t dim = 0; dim < fNdimensions; dim++) {
            if (coord[dim] < 1 || coord[dim] > GetAxis(dim)->GetNbins()) {
                regularBin = false;
                break;
            }
        }
        // if outlayer, skip
        if (!regularBin) continue;
        // as soon as we have non-zero value we can stop
        if (VAL!=0) break;
    }
    delete [] coord;
    fIntegral[0] = VAL/fIntegral[i+1];
    return fIntegral[0];
}

void THnD_KIN::updateBin(){
    // read out the values to one array
    for ( size_t i=0; i < KDIM; i++){
        x[i] = *(val_ptrs[i]);
    }
    bin = GetBin(x);
}





TKinFile::TKinFile(TString name, TString opt){
    const char * c_opt=opt.Data();
    file = TFile::Open(name,opt);
    file->cd();
    ebeam = 1960; // GeV -- tevatron is default
    normND = 0;
    tree = 0;
    if ( !strcasecmp(c_opt,"READ") || !strcasecmp(c_opt,"CACHEREAD") ){
        histND = (THnD_KIN*)  file->Get("kin_hist");
        //normND = (THnD_KIN*)  file->Get("norm_hist");
        //tree   = (TTree*) file->Get("kin_tree");
        //tree->SetBranchAddress("qt"       , &qt       );
        //tree->SetBranchAddress("y"        , &y        );
        //tree->SetBranchAddress("theta_CS" , &theta_CS );
        //tree->SetBranchAddress("phi_CS"   , &phi_CS   );
    } else {
        histND = new THnD_KIN("kin_hist")  ;
        //normND = new THnD_KIN("norm_hist") ;
        //tree   = new TTree("kin_tree","kin_tree");
        //tree->SetBranchAddress("qt/D"       , &qt       );
        //tree->SetBranchAddress("y/D"        , &y        );
        //tree->SetBranchAddress("theta_CS/D" , &theta_CS );
        //tree->SetBranchAddress("phi_CS/D"   , &phi_CS   );
        hw = new TH1D ("weights" , "weights" , 100         , -1           , 10);
        h1 = new TH1D (THnD_KIN::var_name_1, THnD_KIN::var_name_1, THnD_KIN::kin_nbins[0], THnD_KIN::bin_edges_1 );
        h2 = new TH1D (THnD_KIN::var_name_2, THnD_KIN::var_name_2, THnD_KIN::kin_nbins[1], THnD_KIN::bin_edges_2 );
        h3 = new TH1D (THnD_KIN::var_name_3, THnD_KIN::var_name_3, THnD_KIN::kin_nbins[2], THnD_KIN::bin_edges_3 );
        h4 = new TH1D (THnD_KIN::var_name_4, THnD_KIN::var_name_4, THnD_KIN::kin_nbins[3], THnD_KIN::bin_edges_4 );
    }
    qt = y = cos_theta_CS = phi_CS = 0.;
    histND->weight_ptr  = & weight       ;
    histND->val_ptrs[0] = & qt           ;
    histND->val_ptrs[1] = & y            ;
    histND->val_ptrs[2] = & cos_theta_CS ;
    histND->val_ptrs[3] = & phi_CS       ;
}

TKinFile::~TKinFile(){
    if (file && file->IsOpen()) file->Close();
    if (tree)   delete tree;
    if (histND) delete histND;
    if (normND) delete normND;
};

void TKinFile::setValues(const TLorentzVector &VB, const TLorentzVector &l1, const TLorentzVector &l2){
    qt = VB.Pt();
    y = VB.Rapidity();
    Utils::getCSFAngles(l1,-1,l2,ebeam,cos_theta_CS,phi_CS,&aimom_CS);
    return;
}

void TKinFile::Fill(const TLorentzVector &VB, const TLorentzVector &l1, const TLorentzVector &l2, double wt){
    weight = wt;
    setValues(VB, l1, l2);
    histND->Fill();
    hw->Fill( weight       , 1      );
    h1->Fill( qt           , weight );
    h2->Fill( y            , weight );
    h3->Fill( cos_theta_CS , weight );
    h4->Fill( phi_CS       , weight );
    return;
}

void TKinFile::normalizeKin(){
    normND = (THnD_KIN *) histND->Clone("norm");
    normND -> Scale(1./normND->Integral());
    return;
}

double TKinFile::getValue(const double qt, const double y, const double cos_theta_CS, const double phi_CS ){
    if (!normND) normalizeKin();
    this->qt           = qt;
    this->y            = y;
    this->cos_theta_CS = cos_theta_CS;
    this->phi_CS       = phi_CS;
    return normND->GetBinContent();
}

double TKinFile::GetNewWeight(TKinFile * new_kin, const TLorentzVector &VB, const TLorentzVector &l_part, const TLorentzVector &l_anti ){
    setValues(VB,l_part,l_anti);
    double this_wt = getValue(qt,y,cos_theta_CS,phi_CS);
    if (!this_wt) return 0;
    double new_wt = new_kin->getValue(qt,y,cos_theta_CS,phi_CS);
    return new_wt/this_wt;
}

void TKinFile::Save(){
    file->cd();
    histND->Write();
    file->Write();
    file->Close();
};



#endif //_KinFile_CXX
