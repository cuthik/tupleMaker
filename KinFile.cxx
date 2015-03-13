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
const int THnD_KIN::kin_nbins[] = {43,15,10,20} ; // {4,4,4,4};
//const char * THnD_KIN::var_name_1 = "q_{T}"        ; const double THnD_KIN::bin_edges_1[] =  {0.5 , 3   , 6   , 9   , 12  , 15 , 18 , 21 , 24 , 27 , 30 , 35 , 40 , 50 , 100, 200, 1000 } ;// {0.,5.,10.,100.,7000.}   ;
const char * THnD_KIN::var_name_1 = "q_{T}"        ; const double THnD_KIN::bin_edges_1[] =  { 0.5    , 1.09   , 1.68   , 2.27   , 2.86   , 3.45   , 4.04   , 4.63   , 5.22   , 5.81   , 6.4    , 6.99   , 7.58   , 8.17   , 8.76   , 9.35   , 9.94   , 10.53  , 11.12  , 11.71  , 12.3   , 12.89  , 13.48  , 14.07  , 14.66  , 15.25  , 15.84  , 16.43  , 17.02  , 17.61  , 18.2   , 18.79  , 19.38  , 19.97  , 20.56  , 21.15  , 21.74  , 22.33  , 22.92  , 23.51  , 24.1   , 24.69  , 25.2   , 30.0 };

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
    //if (isOutlayer()) return 0; // not used any more, too many events were zero
    return THnD::GetBinContent(bin);
}


Double_t THnD_KIN::ComputeIntegral() {
   // Calculate the integral of the histogram
   // Taken from THnBase and corrected

   // delete old integral
   if (fIntegralStatus != kNoInt) {
      delete [] fIntegral;
      fIntegralStatus = kNoInt;
   }

   // check number of bins
   if (GetNbins() == 0) {
      Error("ComputeIntegral", "The histogram must have at least one bin.");
      return 0.;
   }

   // allocate integral array
   fIntegral = new Double_t [GetNbins() + 1];
   fIntegral[0] = 0.;

   // fill integral array with contents of regular bins (non over/underflow)
   Int_t* coord = new Int_t[fNdimensions];
   Long64_t i = 0;
   THnIter iter(this);
   while ((i = iter.Next(coord)) >= 0) {
      Double_t v = THnD::GetBinContent(i);
      fIntegral[i + 1] = fIntegral[i] + v;
   }
   delete [] coord;

   // check sum of weights
   if (fIntegral[GetNbins()] == 0.) {
      Error("ComputeIntegral", "No hits in bins (including over/underflow).");
      delete [] fIntegral;
      return 0.;
   }

   // normalize the integral array
   for (Long64_t j = 0; j < GetNbins(); ++j)
      fIntegral[j] = fIntegral[j] / fIntegral[GetNbins()];

   // set status to valid
   fIntegralStatus = kValidInt;
   return fIntegral[GetNbins()];
}



double THnD_KIN::Integral(){
    
    // After som time I realize that the function from THnBase sucks so I
    // copied the code a do my own integral.
    if(fIntegralStatus == kValidInt) return fIntegral[fNdimensions];
    return ComputeIntegral();
    // you can follow if you are interested..
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
    // 
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

void THnD_KIN::CopyPointers(const THnD_KIN * from){
    this->weight_ptr = from->weight_ptr;
    for (size_t i=0; i< fNdimensions;i++) this->val_ptrs[i] = from->val_ptrs[i];
}





TKinFile::TKinFile(TString name, TString opt){
    TH1::SetDefaultSumw2(true);
    const char * c_opt=opt.Data();
    file = TFile::Open(name,opt);
    file->cd();
    ebeam = 1960; // GeV -- tevatron is default
    normND = 0;
    tree = 0;
    use4Dhist = true;
    if ( !strcasecmp(c_opt,"READ") || !strcasecmp(c_opt,"CACHEREAD") ){
        histND = (THnD_KIN*)  file->Get("kin_hist");
        hw     = (TH1D*) file->Get("weights");
        h1     = (TH1D*) file->Get(THnD_KIN::var_name_1);
        h2     = (TH1D*) file->Get(THnD_KIN::var_name_2);
        h3     = (TH1D*) file->Get(THnD_KIN::var_name_3);
        h4     = (TH1D*) file->Get(THnD_KIN::var_name_4);
    } else {
        histND = new THnD_KIN("kin_hist")  ;
        hw = new TH1D ("weights" , "weights" , 100         , -1           , 10);
        h1 = new TH1D (THnD_KIN::var_name_1, THnD_KIN::var_name_1, THnD_KIN::kin_nbins[0], THnD_KIN::bin_edges_1 );
        h2 = new TH1D (THnD_KIN::var_name_2, THnD_KIN::var_name_2, THnD_KIN::kin_nbins[1], THnD_KIN::bin_edges_2 );
        h3 = new TH1D (THnD_KIN::var_name_3, THnD_KIN::var_name_3, THnD_KIN::kin_nbins[2], THnD_KIN::bin_edges_3 );
        h4 = new TH1D (THnD_KIN::var_name_4, THnD_KIN::var_name_4, THnD_KIN::kin_nbins[3], THnD_KIN::bin_edges_4 );

    }
    qt = y = cos_theta_CS = phi_CS = 0.;
    histND->weight_ptr  = & this->weight       ;
    histND->val_ptrs[0] = & this->qt           ;
    histND->val_ptrs[1] = & this->y            ;
    histND->val_ptrs[2] = & this->cos_theta_CS ;
    histND->val_ptrs[3] = & this->phi_CS       ;

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
    normND->CopyPointers(histND);
    normND -> Scale(1./normND->Integral());
    return;
}

double TKinFile::getValue(const double qt, const double y, const double cos_theta_CS, const double phi_CS ){
    if (!normND) normalizeKin();
    this->qt           = qt;
    this->y            = y;
    this->cos_theta_CS = cos_theta_CS;
    this->phi_CS       = phi_CS;
    if (use4Dhist) return normND->GetBinContent();
    else {
        double out = 1;
        out *= h1 ->GetBinContent( h1 ->FindFixBin( qt            ));
        out *= h2 ->GetBinContent( h2 ->FindFixBin( y             ));
        out *= h3 ->GetBinContent( h3 ->FindFixBin( cos_theta_CS  ));
        out *= h4 ->GetBinContent( h4 ->FindFixBin( phi_CS        ));
        out *= 1./h1->Integral();
        out *= 1./h2->Integral();
        out *= 1./h3->Integral();
        out *= 1./h4->Integral();
        return out;
    }
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
