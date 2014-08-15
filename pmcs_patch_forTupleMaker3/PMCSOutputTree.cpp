#include "PMCSOutputTree.hpp"

#include "TFile.h"

#include <iostream>
using std::cout;
using std::endl;


PMCSOutputTree::PMCSOutputTree() :
    outTree       ( 0  ),
    truth_VB_mass ( 0. ),
    truth_VB_pt   ( 0. ),
    truth_VB_phi  ( 0. ),
    truth_VB_eta  ( 0. ),
    truth_el_pt   ( 0. ),
    truth_el_phi  ( 0. ),
    truth_el_eta  ( 0. ),
    truth_nu_pt   ( 0. ),
    truth_nu_phi  ( 0. ),
    truth_nu_eta  ( 0. ),
    VB_mt         ( 0. ),
    MET_et        ( 0. ),
    MET_sumet     ( 0. ),
    el_pt         ( 0. ),
    el_phi        ( 0. ),
    hr_et         ( 0. ),
    hr_phi        ( 0. ),
    hr_prll       ( 0. ),
    hr_perp       ( 0. ),
    weight_evt    ( 0. ),
    weights_mass  ( new vector<double>()  ),
    weights_PDF   ( new vector<double>()  )
{

}


void PMCSOutputTree::Init(){
  outTree = new TTree("physics","physics of pmcs");

  outTree -> Branch( "truth_evtn" , &truth_evtn );

  outTree -> Branch( "truth_VB_mass" , &truth_VB_mass );
  outTree -> Branch( "truth_VB_pt"   , &truth_VB_pt   );
  outTree -> Branch( "truth_VB_phi"  , &truth_VB_phi  );
  outTree -> Branch( "truth_VB_eta"  , &truth_VB_eta  );

  outTree -> Branch( "truth_el_pt"   , &truth_el_pt   );
  outTree -> Branch( "truth_el_phi"  , &truth_el_phi  );
  outTree -> Branch( "truth_el_eta"  , &truth_el_eta  );

  outTree -> Branch( "truth_nu_pt"   , &truth_nu_pt   );
  outTree -> Branch( "truth_nu_phi"  , &truth_nu_phi  );
  outTree -> Branch( "truth_nu_eta"  , &truth_nu_eta  );

  outTree -> Branch( "VB_mt"         , &VB_mt         );

  outTree -> Branch( "MET_et"        , &MET_et        );
  outTree -> Branch( "MET_sumet"     , &MET_sumet     );

  outTree -> Branch( "el_pt"         , &el_pt         );
  outTree -> Branch( "el_phi"        , &el_phi        );


  outTree -> Branch( "hr_et"         , &hr_et         );
  outTree -> Branch( "hr_phi"        , &hr_phi        );
  outTree -> Branch( "hr_prll"       , &hr_prll       );
  outTree -> Branch( "hr_perp"       , &hr_perp       );


  outTree -> Branch( "weight_evt"    , &weight_evt    );
  outTree -> Branch( "weights_mass"  , weights_mass   );
  outTree -> Branch( "weights_PDF"   , weights_PDF    );

  Clear();
}


void PMCSOutputTree::Clear(){
    truth_VB_mass = -999. ;
    truth_VB_pt   = -999. ;
    truth_VB_phi  = -999. ;
    truth_VB_eta  = -999. ;
    truth_el_pt   = -999. ;
    truth_el_phi  = -999. ;
    truth_el_eta  = -999. ;
    truth_nu_pt   = -999. ;
    truth_nu_phi  = -999. ;
    truth_nu_eta  = -999. ;
    VB_mt         = -999. ;
    MET_et        = -999. ;
    MET_sumet     = -999. ;
    el_pt         = -999. ;
    el_phi        = -999. ;
    hr_et         = -999. ;
    hr_phi        = -999. ;
    hr_prll       = -999. ;
    hr_perp       = -999. ;
    weight_evt    = -999. ;
    weights_mass  ->clear();
    weights_PDF   ->clear();
}


void PMCSOutputTree::SetPDFWeights(double * pdfarray, size_t n_pdf){
    weights_PDF->clear();
    for(size_t i=0; i<n_pdf; i++){
        weights_PDF->push_back(pdfarray[i]);
    }
}


void PMCSOutputTree::SetMassWeights(vector<double> * in){
    weights_mass->clear();
    for(vector<double>::const_iterator it = in->begin(); it != in->end(); ++it) weights_mass->push_back((*it));
}


void PMCSOutputTree::Fill(){
    outTree->Fill();
}


void PMCSOutputTree::SetWenu(PMCSWCand * wcand, double scalarEt_All){
    VB_mt  = wcand-> TMass()                 ; // WMt       ;
    hr_et  = wcand-> MagRecoil()             ; // WRecoilPt ;
    MET_et = wcand-> GetMet().met()          ; // WMet      ;
    el_pt  = wcand-> GetElec().ppt()         ; // WElecPt   ;
    el_phi = wcand-> GetElec().pphi()        ;
    hr_phi = wcand-> GetRecoil().RecoilPhi() ;
    hr_prll   = wcand->UPara() ;
    hr_perp   = wcand->UPerp() ;
    MET_sumet = scalarEt_All   ;
}


void PMCSOutputTree::SetTruthWenu(PMCSEvent *pmcsevent, vector<PMCSEMObj> * emobj_gen, PMCSWCand * wcand){
    truth_evtn = pmcsevent->_truth_evtn;
    weight_evt = pmcsevent->GetEvtWeight();

    truth_VB_mass = pmcsevent->GetWZBoson() .Mass() ;
    truth_VB_pt   = pmcsevent->GetWZBoson() .ppt()  ;
    truth_VB_phi  = pmcsevent->GetWZBoson() .pphi() ;
    truth_VB_eta  = pmcsevent->GetWZBoson() .peta() ;


    for (vector<PMCSEMObj>::const_iterator it = emobj_gen->begin(); it != emobj_gen->end(); ++it){
        PMCSEMObj obj = (*it);
        if ( obj.getIndex() == wcand->GetElec().getIndex() ){
            truth_el_pt  = obj .ppt()  ;
            truth_el_eta = obj .peta() ;
            truth_el_phi = obj .pphi() ;
        }
        if (fabs(obj.ppid()) == 12.){
            truth_nu_pt  = obj .ppt()  ;
            truth_nu_eta = obj .peta() ;
            truth_nu_phi = obj .pphi() ;
        }
    }
}



void PMCSOutputTree::Save(TString filename){
    cout << "Writing ntuple file " << filename.Data() << ":/" << endl;
    TFile * file = TFile::Open(filename.Data(), "RECREATE");
    file->cd();
    outTree->Write();
    file->Close();
    cout << "Done." << endl;
}


