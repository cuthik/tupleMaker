/**
 * @file plot_fast.C
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-04-14
 */


#include "ReadResbosROOT.C"
#include "TLVUtils.h"

#include "TCanvas.h"
#include "TStyle.h"
//#include "TFile.h"
#include "TStopwatch.h"
#include "TUnixSystem.h"

#include <map>
#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;
using TMath::Log2;

typedef map<TString,TH1D*> HMAP;

class FastPlot {
    public :
        ~FastPlot(){};
        FastPlot(TString _tname, TString _figname){
            filename = _tname;
            outdir   = "full";
            figname  = _figname;
            t_all.push_back( new ResbosRootNtuple (filename.Data()) );
            book(t_all.size()-1);
            tt = t_all.back();
            GeV = 1.;
        }

        void AddWeight(TString fname){
            t_all.push_back( new ResbosRootNtuple ( fname.Data(), "WEIGHT") );
            book(t_all.size()-1);
        }


        void MakePlots(){
            int N = tt->GetEntries();
            if (t_all.size() > 1){
                t_wgts.assign(t_all.begin()+1, t_all.end());
                tt->ConnectWeights(&t_wgts);
            }
            for (int ievt = 0; ievt < N; ievt++){
                /// get entry
                tt->GetEntry(ievt);
                tt->GetConnectedEntry(ievt);
                printEvent(ievt,N);
                //printf( "========= conf ebeam %f VB %d lep1 %d lep2 %d \n", tt->ebeam, tt->typeVB, tt->type1, tt->type2 );
                /// compute event
                tt->GetLorentzVector();
                Utils::getCSFAngles(tt->v_l1, ((tt->type1>0) ? -1 : 1), tt->v_l2, tt->ebeam, costhcs, phics, 0);
                l1pt=tt->v_l1.Pt();
                l2pt=tt->v_l2.Pt();
                l1aeta=fabs(tt->v_l1.Eta());
                l2aeta=fabs(tt->v_l2.Eta());
                Mll=(tt->v_l1+tt->v_l2).M();
                /// do cut
                outdir = "fiducial";
                if (
                     l1aeta > 2.4 || l1pt < 20*GeV ||
                     l2aeta > 2.4 || l2pt < 20*GeV ||
                     Mll < 66*GeV ||  Mll > 116*GeV
                    ) continue;
                //outdir = "y3";       if ( fabs(fabs(tt->v_VB.Rapidity())-3) > 0.2  ) continue;
                //outdir = "wg_neg";   if ( tt->wgt_corrected                 > 0    ) continue;
                //outdir = "wg_0";     if ( tt->wgt_corrected                 != 0    ) continue;
                //outdir = "wg_0band"; if ( fabs(tt->wgt_corrected)           > 0.01 ) continue;
                /// fill all
                for (size_t iwgt=0;iwgt<t_all.size();iwgt++) fill(iwgt);
            }
            outdir = "img_"+outdir;
            draw_hist();
            return;
        }

    private :


        void fill( int i = 0){
            double weight = tt->wgt_corrected;
            if (i!=0){
                weight*=t_all[i]->wgt/tt->wgt;
            }
            if (weight<0) return;
            get_hist("weight"  , i) -> Fill(weight   *1 );
            get_hist("qt"      , i) -> Fill(tt->v_VB.Pt()       *GeV  , weight);
            //get_hist("phi"     , i) -> Fill(tt->v_VB.Phi()      *1    , weight);
            //get_hist("Q"       , i) -> Fill(tt->v_VB.M()        *GeV  , weight);
            //get_hist("y"       , i) -> Fill(tt->v_VB.Rapidity() *1    , weight);
            //get_hist("lepPt"   , i) -> Fill(tt->v_l1.Pt()       *GeV  , weight);
            get_hist("phiCS"   , i) -> Fill(phics             *1    , weight);
            get_hist("costhCS" , i) -> Fill(costhcs           *1    , weight);
            if ( 0. <= fabs(tt->v_VB.Rapidity()) && fabs(tt->v_VB.Rapidity()) < 1.  ) get_hist("qt_0y1"    ,i) -> Fill(tt->v_VB.Pt()       *GeV  ,  weight);
            if ( 1. <= fabs(tt->v_VB.Rapidity()) && fabs(tt->v_VB.Rapidity()) < 2.  ) get_hist("qt_1y2"    ,i) -> Fill(tt->v_VB.Pt()       *GeV  ,  weight);
            if ( 2. <= fabs(tt->v_VB.Rapidity()) && fabs(tt->v_VB.Rapidity()) < 2.4 ) get_hist("qt_2y2.4"  ,i) -> Fill(tt->v_VB.Pt()       *GeV  ,  weight);
        }

        void book( int i = 0){
            book_hist( "weight"  , 100, -0.5, 6  , true  , i );
            book_hist( "qt"      , 500, 0.  , 100, false , i );
            book_hist( "qt_0y1"  , 500, 0.  , 100, false , i );
            book_hist( "qt_1y2"  , 500, 0.  , 100, false , i );
            book_hist( "qt_2y2.4", 500, 0.  , 100, false , i );
            //book_hist( "Q"       , 400, 40  , 120, false , i );
            //book_hist( "phi"     , 500, -4  , 4  , false , i );
            //book_hist( "y"       , 100, -5  , 5  , false , i );
            //book_hist( "lepPt"   , 500, 0.  , 100, false , i );
            book_hist( "phiCS"   , 500, -4  , 4  , false , i );
            book_hist( "costhCS" , 500, -1  , 1  , false , i );
        }


        void printEvent(double evn, double total=0){
            if (evn==0){
                swatch.Reset();
                swatch.Start();
            }
            if( ! fmod( Log2(evn), 1 ) ) {
                cout << "Processed event: ";
                if (total !=0 ){
                    cout << " " ;
                    cout << fixed ;
                    cout << setfill('0') ;
                    cout << setw(5) ;
                    cout << setprecision(2) ;
                    cout << evn/total*100 ;
                    cout << "%  ";
                } 
                cout << setprecision(6) << setfill(' ') << setw(10) ;
                cout << (size_t) evn;
                // time
                double rtm = swatch.RealTime(); swatch.Continue();
                double Hz = evn/rtm;
                cout << "   time elapsed ";
                cout << rtm << "s";
                cout << " ( "<< Hz/1000 <<" kHz)";
                if (total !=0 ){
                    cout << "   remaining time ";
                    cout << (total-evn)/Hz;
                }
                cout << endl;
                //cout << "\r";
                //cout.flush();
            }
        }

        void get_hist_name(TString &name, int i=0){
            if (i==0) return;
            name+="_ct";
            name+=Form("%02d",i);
        }

        TH1D * get_hist (TString name, int i=0){
            get_hist_name(name,i);
            return h[name.Data()];
        }

        void book_hist( TString name, int N, double min, double max, bool logY=false, int iwgt = 0 ){
            get_hist_name(name, iwgt);
            TString title = name + ";" + name + ";" + TString(logY ?  "log" : "entries");
            h[name.Data()] = new TH1D (name.Data(), title.Data(), N, min, max);
        }


        void draw_hist(){
            bool doDraw = false;
            sys.mkdir(outdir.Data(),true);
            outdir+="/h_";
            //
            TH1D *hi =0;
            TString filename = outdir+figname+".root";
            TFile *  outf = TFile::Open(filename.Data(),"RECREATE");
            TCanvas *c1 = new TCanvas("all", "all", 1600, 1200 );
            TCanvas *c2 = new TCanvas("one", "one", 800, 600 );
            c1->Divide(3,3);
            gStyle->SetOptStat(111110);
            int i = 1;
            for(HMAP::iterator it = h.begin(); it != h.end(); ++it) {
                hi = (*it).second;
                if(doDraw){
                    // separate canvas
                    c2->cd(0);
                    gPad->SetLogy(0);
                    hi->Draw();
                    filename = outdir+figname+"_"+hi->GetName()+TString(".pdf");
                    if (!strcmp(hi->GetYaxis()->GetTitle(),"log")) gPad->SetLogy(1); gPad->Update();
                    gPad->SaveAs(filename.Data());
                    // all in one canvas
                    c1->cd(i);
                    hi->Draw();
                    if (!strcmp(hi->GetYaxis()->GetTitle(),"log")) gPad->SetLogy(1); gPad->Update();
                }
                outf->cd();
                hi->Write();
                i++;
            }
            if(doDraw){
                c1->cd();
                filename = outdir+figname+"_all.pdf";
                c1->SaveAs(filename.Data());
            }
            delete c1;
            delete c2;
            outf->Write();
            outf->Close();
        }


        // DATA members
        TStopwatch swatch;
        TUnixSystem sys;
        HMAP h;
        TString outdir;
        TString infile;
        TString figname;
        TString filename;
        double GeV;
        vector<ResbosRootNtuple *> t_all;
        vector<ResbosRootNtuple *> t_wgts;
        ResbosRootNtuple * tt;
        double costhcs;
        double phics;
        double l1pt;
        double l2pt;
        double l1aeta;
        double l2aeta;
        double Mll;
};




int main(int argc, const char * argv[]){
    //plot_fast ("/home/cuth/workdir_etapfs/resbos/results/resbos_zmumu_lhc7_ct00_QTfull_100105.tree.root"  , "z0");
    //plot_fast ("/home/cuth/workdir_etapfs/resbos/results/resbos_wpmunu_lhc7_ct00_QTfull_100101.tree.root" , "wp");
    //plot_fast ("/home/cuth/workdir_etapfs/resbos/results/resbos_wmmunu_lhc7_ct00_QTfull_100101.tree.root" , "wm");
    //FastPlot a ("/home/cuth/workdir_etapfs/resbos/results/resbos_zmumu_lhc7_ct00_QTfull_100101.tree.root", "z0_100101");
    //a.AddWeight("/home/cuth/workdir_etapfs/resbos/results/resbos_zmumu_lhc7_ct01_QTfull_100101.weights.root");
    //a.AddWeight("/home/cuth/workdir_etapfs/resbos/results/resbos_zmumu_lhc7_ct02_QTfull_100101.weights.root");
    //a.MakePlots();
    
    if ( argc < 3 ) return 1;
    printf("out: %s  in: %s \n", argv[1],argv[2]);
    FastPlot a (argv[2],argv[1]);
    for (int i = 3; i<argc; i++){
        printf("wgt: %02d  %s \n", i-2,argv[i]);
        a.AddWeight(argv[i]);
    } 
    a.MakePlots();
    return 0;
}

