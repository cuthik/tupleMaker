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

double THnD_KIN::GetBinContent(){
    updateBin();
    return THnD::GetBinContent(bin);
}

void THnD_KIN::updateBin(){
    // read out the values to one array
    Dump();

    for ( size_t i=0; i < KDIM; i++){
        printf("val_ptr %d %p %f", i, val_ptrs[i], *val_ptrs[i]);
        x[i] = *(val_ptrs[i]);
        printf(" x %p %f\n", &x[i], x[i]);
        printf("   -- axis %d ptr %p\n", i, fAxes[i]                 );
        //printf("   -- axis %d bin %d\n", i, GetAxis(i)->FindFixBin(x[i]));
    }
    bin = GetBin(x);
}





TKinFile::TKinFile(TString name, TString opt){
    const char * c_opt=opt.Data();
    file = TFile::Open(name,opt);
    printf("Tkinkonstruktor\n");
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
        h1 = new TH1D (THnD_KIN::var_name_1, THnD_KIN::var_name_1, THnD_KIN::kin_nbins[0], THnD_KIN::bin_edges_1      );
        h2 = new TH1D (THnD_KIN::var_name_2, THnD_KIN::var_name_2, THnD_KIN::kin_nbins[1], THnD_KIN::bin_edges_2      );
        h3 = new TH1D (THnD_KIN::var_name_3, THnD_KIN::var_name_3, THnD_KIN::kin_nbins[2], THnD_KIN::bin_edges_3      );
        h4 = new TH1D (THnD_KIN::var_name_4, THnD_KIN::var_name_4, THnD_KIN::kin_nbins[3], THnD_KIN::bin_edges_4      );
    }

    qt = y = cos_theta_CS = phi_CS = 0.;

    histND->weight_ptr  = & weight       ;
    histND->val_ptrs[0] = & qt           ;
    histND->val_ptrs[1] = & y            ;
    histND->val_ptrs[2] = & cos_theta_CS ;
    histND->val_ptrs[3] = & phi_CS       ;
}
TKinFile::~TKinFile(){
    printf("kinfile destruktor\n");
    file->Dump();
    if (file && file->IsOpen()) file->Close();
    if (tree)   delete tree;
    if (histND) delete histND;
    if (normND) delete normND;
};

void TKinFile::setValues(const TLorentzVector &VB, const TLorentzVector &l1, const TLorentzVector &l2){
    qt = VB.Pt();
    y = VB.Rapidity();
    Utils::getCSFAngles(l1,-1,l2,ebeam,cos_theta_CS,phi_CS,&aimom_CS);
    printf(" qt %p y %p costh %p phi %p\n", &qt, &y, &cos_theta_CS, &phi_CS);
    printf(" qt %f y %f costh %f phi %f\n", qt, y, cos_theta_CS, phi_CS);
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
    printf("integral old %f\n", histND->ComputeIntegral() );
    normND = (THnD_KIN *) histND->Clone("norm");
    normND -> Scale(1./normND->ComputeIntegral());
    printf("integral old %f new %f\n", histND->ComputeIntegral(), normND->ComputeIntegral() );
}

double TKinFile::getValue(const double qt, const double y, const double cos_theta_CS, const double phi_CS ){
    if (!normND) normalizeKin();
    this->qt           = qt;
    this->y            = y;
    this->cos_theta_CS = cos_theta_CS;
    this->phi_CS       = phi_CS;
    //return normND->GetBinContent();
    return histND->GetBinContent();
}

double TKinFile::GetNewWeight(TKinFile * new_kin, const TLorentzVector &VB, const TLorentzVector &l_part, const TLorentzVector &l_anti ){
    setValues(VB,l_part,l_anti);
    double new_wt = new_kin->getValue(qt,y,cos_theta_CS,phi_CS);
    double this_wt = getValue(qt,y,cos_theta_CS,phi_CS);
    return new_wt/this_wt;
}

void TKinFile::Save(){
    histND->Dump();
    file->cd();
    histND->Write();
    file->Write();
    file->Close();
};



#endif //_KinFile_CXX
