#ifndef _HistogramsVB_H
#define _HistogramsVB_H

/**
 * @file HistogramsVB.h
 *
 * @brief A header file for class HistogramVB.
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-03-19
 */

#include "MWQuickSim.h"

#include "THashList.h"
#include "TAxis.h"

#include <map>

using std::map;

/**
 * Description of class.
 *
 * Longer description of class
 */
class HistogramsVB {
    public:
         HistogramsVB(){
            output_name="hists.root";
            currentUnitInMeV=1.; // unit of input 1. == MeV; 1e-3 == GeV, 1e-6 == TeV
            pi = TMath::Pi();
            N_reweightings=0;
         };
        ~HistogramsVB(){};

        TH1D * h1(TString name){
            TH1D *h = (TH1D*) h_list(name.Data());
            if(!h) printf ("can not find object %s\n",name.Data());
            return h;
        }

        void Init(){
            GeV=1e3; // I don't like to write zeros
            // axis binning
            m_axis["VBpt"] = TAxis(120, 0, 120*GeV);
            m_axis["Lpt"]  = TAxis(120, 0, 120*GeV);
            m_axis["pt"]   = TAxis(120, 0, 120*GeV);
            m_axis["eta"]  = TAxis(60, -3, 3);
            m_axis["weta"]  = TAxis(60, -8, 8);
            m_axis["phi"]  = TAxis(50, -pi,pi);
            // truth
            book1d("truth_VB_m"  , "truth mass (VB) [MeV]" , m_axis["VBpt"] );
            book1d("truth_VB_mt" , "truth m_{T} (VB) [MeV]", m_axis["VBpt"] );
            book1d("truth_VB_pt" , "truth p_{T}(VB) [MeV]" , m_axis["VBpt"] );
            book1d("truth_VB_eta", "truth #eta(VB)"        , m_axis["weta"] );
            book1d("truth_l1_pt" , "truth p_{T}(l1) [MeV]" , m_axis["Lpt"]  );
            book1d("truth_l1_eta", "truth #eta(l1)"        , m_axis["weta"] );
            book1d("truth_l2_pt" , "truth p_{T}(l2) [MeV]" , m_axis["Lpt"]  );
            book1d("truth_l2_eta", "truth #eta(l2)"        , m_axis["weta"] );
            // reconstructed
            book1d("reco_VB_mt" , "m_T [MeV]"       , m_axis["VBpt"] );
            book1d("reco_VB_HR" , "HR [MeV]"        , m_axis["VBpt"] );
            book1d("reco_l_pt"  , "p_{T}(lep) [MeV]", m_axis["Lpt"]  );
            book1d("reco_l_eta" , "#eta(lep)"       , m_axis["eta"]  );
            book1d("reco_met_pt", "MET [MeV]"       , m_axis["Lpt"]  );
            bookReweighted();
        }

        void book1d(TString name, TString title, TAxis axis){
            TH1D *h =  new TH1D(name.Data(), Form("%s;%s;entries",name.Data(),title.Data()), axis.GetNbins(), axis.GetXmin(), axis.GetXmax());
            h->Sumw2();
            h->SetDirectory(0);
            h_list.Add(h);
        }

        TString getRewName(TString base, int index){
            if (index<0) return base;
            base += "_rew";
            base += index;
            return base;
        }

        void bookReweighted(){
            if(!N_reweightings) return;
            // create vector of all current names, duplicate with weight
            vector<TString> hnames;
            TIter iter(&h_list);
            for (TObject * obj = iter.Next(); obj!=0; obj=iter.Next()){
                hnames.push_back(obj->GetName());
            }
            // loop all names
            for (int i_rw=0; i_rw<N_reweightings; i_rw++){
                for (size_t i_h=0; i_h<hnames.size(); i_h++) {
                    TString hbase=hnames[i_h];
                    TH1* h = (TH1*) h1(hbase)->Clone(getRewName(hbase,i_rw).Data());
                    h->Reset();
                    h->SetDirectory(0);
                    h_list.Add( h );
                }
            }
        }

        void Save(){
            TFile* out = TFile::Open(output_name.Data(), "RECREATE");
            TIter iter(&h_list);
            for (TObject * obj = iter.Next(); obj!=0; obj=iter.Next()){
                obj->Write();
            }
            out->Write();
            out->Close();
        }

        double mT(const TLorentzVector &l1, const TLorentzVector &l2){
            return TMath::Sqrt(2.0*l1.Pt()*l2.Pt()*(1.0-TMath::Cos(l1.DeltaPhi(l2))));
        }

        void SetTruthKinematics(TLorentzVector _VB, int _VBpid, TLorentzVector _l1, int _l1pid, TLorentzVector _l2, int _l2pid){
            VB_pid = _VBpid;
            isW = isZ = false;
            if (abs(_VBpid)==24) isW=true;
            if (abs(_VBpid)==23) isZ=true;
            vVB_truth.SetPtEtaPhiM( _VB.Pt() * currentUnitInMeV, _VB.Eta(), _VB.Phi(), _VB.M()*currentUnitInMeV );

            if (isZ){ // first is lepton than anti-lepton
                if(_l1pid>_l2pid){
                    l1_pid = _l1pid;
                    l2_pid = _l2pid;
                    vl1_truth.SetPtEtaPhiM( _l1.Pt() * currentUnitInMeV, _l1.Eta(), _l1.Phi(), _l1.M()*currentUnitInMeV );
                    vl2_truth.SetPtEtaPhiM( _l2.Pt() * currentUnitInMeV, _l2.Eta(), _l2.Phi(), _l2.M()*currentUnitInMeV );
                } else {
                    l2_pid = _l1pid;
                    l1_pid = _l2pid;
                    vl2_truth.SetPtEtaPhiM( _l1.Pt() * currentUnitInMeV, _l1.Eta(), _l1.Phi(), _l1.M()*currentUnitInMeV );
                    vl1_truth.SetPtEtaPhiM( _l2.Pt() * currentUnitInMeV, _l2.Eta(), _l2.Phi(), _l2.M()*currentUnitInMeV );
                }
            }
            if (isW){ // first is lepton than neutrino
                if(abs(_l1pid)<abs(_l2pid)){
                    l1_pid = _l1pid;
                    l2_pid = _l2pid;
                    vl1_truth.SetPtEtaPhiM( _l1.Pt() * currentUnitInMeV, _l1.Eta(), _l1.Phi(), _l1.M()*currentUnitInMeV );
                    vl2_truth.SetPtEtaPhiM( _l2.Pt() * currentUnitInMeV, _l2.Eta(), _l2.Phi(), _l2.M()*currentUnitInMeV );
                } else {
                    l2_pid = _l1pid;
                    l1_pid = _l2pid;
                    vl2_truth.SetPtEtaPhiM( _l1.Pt() * currentUnitInMeV, _l1.Eta(), _l1.Phi(), _l1.M()*currentUnitInMeV );
                    vl1_truth.SetPtEtaPhiM( _l2.Pt() * currentUnitInMeV, _l2.Eta(), _l2.Phi(), _l2.M()*currentUnitInMeV );
                }

            }

            mT_truth = mT(vl1_truth,vl2_truth);

            isFiducial = true;
            if(vVB_truth.M() < 40*GeV) isFiducial=false;
            if(vVB_truth.M() > 120*GeV) isFiducial=false;
            if(TMath::Abs(vl1_truth.Eta()) > 2.5) isFiducial=false;
            if(TMath::Abs(vl2_truth.Eta()) > 2.5) isFiducial=false;
        }

        void CalculateQuickReconstruction(){
            if(isW) use_quick_sim = quick_sim.reconstructWEvent(vl1_truth, vl2_truth);
            if(isZ) use_quick_sim = quick_sim.reconstructZEvent(vl1_truth, vl2_truth);
        }

        void SetWeight(bool _wgt, vector<float>* _rew = 0){
            weight = _wgt;
            reweightings = _rew;
        }

        void fill1D(TString name, double value){
            for (int i = -1; i < N_reweightings; i++ ){
                double iw = weight;
                if(i>=0) iw *= reweightings->at(i);
                h1(getRewName(name,i))->Fill(value,iw);
            }
        }

        void Fill(){
            if(!isFiducial) return;
            fill1D("truth_VB_m"  , vVB_truth.M()  );
            fill1D("truth_VB_mt" , mT_truth       );
            fill1D("truth_VB_pt" , vVB_truth.Pt() );
            fill1D("truth_VB_eta" , vVB_truth.Eta() );

            fill1D("truth_l1_pt" , vl1_truth.Pt() );
            fill1D("truth_l2_pt" , vl2_truth.Pt() );
            fill1D("truth_l1_eta" , vl1_truth.Eta() );
            fill1D("truth_l2_eta" , vl2_truth.Eta() );


            if(use_quick_sim){
                fill1D("reco_VB_mt" , quick_sim.RecoTransverseMass      );
                fill1D("reco_VB_HR" , quick_sim.RecoHadronicRecoil.Pt() );
                if (isW){
                    fill1D("reco_l_pt"    , quick_sim.RecoMuon   .Pt() );
                    fill1D("reco_l_eta"    , quick_sim.RecoMuon   .Eta() );
                    fill1D("reco_met_pt"  , quick_sim.RecoETMiss .Pt() );
                }
                if (isZ){
                    fill1D("reco_l_pt"  , quick_sim.RecoMuon1 .Pt() );
                    fill1D("reco_l_pt"  , quick_sim.RecoMuon2 .Pt() );
                    fill1D("reco_l_eta"  , quick_sim.RecoMuon1 .Eta() );
                    fill1D("reco_l_eta"  , quick_sim.RecoMuon2 .Eta() );
                }

            }
        }


        // data members
        // configs
        TString output_name;
        double currentUnitInMeV;
        int N_reweightings;
        bool use_quick_sim;
        // kinematics and event info
        double mT_truth;
        TLorentzVector vVB_truth;
        TLorentzVector vl1_truth;
        TLorentzVector vl2_truth;
        double mT_reco;
        TLorentzVector vVB_reco;
        TLorentzVector vl1_reco;
        TLorentzVector vl2_reco;
        TLorentzVector vVB_transv_reco;
        TLorentzVector vHR_reco;
        double weight;
        vector<float>* reweightings;
        MWQuickSim quick_sim;
        bool isW;
        bool isZ;
        bool isFiducial;
        int VB_pid;
        int l1_pid;
        int l2_pid;
        // containers
        THashList h_list;
        map<const char*, TAxis> m_axis;

    private:
        double GeV;
        double pi;

};

#endif //_HistogramsVB_H
