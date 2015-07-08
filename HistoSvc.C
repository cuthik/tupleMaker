#ifndef _HISTOSVC_H
#define _HISTOSVC_H

/**
 * @file HistoSvc.C
 *
 * @brief A header file for class HistoSvc.
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-06-18
 */

#include <TList.h>
#include <stdexcept>

/**
 * Description of class.
 *
 * Longer description of class
 */
class HistoSvc {
    public:
         HistoSvc(){};
        ~HistoSvc(){};
        //HistoSvc( const HistoSvc & other);


        void AddBinning(TString name, TString xtitle, int N, double xlo, double xhi ){
            TString title = ";";
            title+=xtitle;
            title+=";";
            TH1D * h = new TH1D(name.Data(),title.Data(),N,xlo,xhi);
            l_bins.Add(h);
        }
        //AddBinning(TString name, int N, double * xbins  ); // @TODO

        template <class T>
        void Book(TString name, TString v1, TString v2="", TString v3="");

        template <class T>
        void Fill(TString name, double v1, double v2=1, double v3=1, double v4=1);

        template<class T>
        T * Get(TString name);

        void Write(){
            for( TIter it(l_h.begin()); it!=l_h.end(); ++it){
                (*it)->Write();
            }
        }

        inline void SetWeight(double _in){ wgt = _in; }


    private:
        // data memebers
        TList l_bins;
        TList l_h;
        double wgt;

        // methods
        void getBinning(TString name, int &N, double &xlo, double &xhi, TString &tit){
            TH1D *h = (TH1D*) l_bins(name.Data());
            if (!h) throw invalid_argument(Form("There is no binning called `%s`", name.Data() ));
            N   = h->GetXaxis()->GetNbins();
            xlo = h->GetXaxis()->GetXmin();
            xhi = h->GetXaxis()->GetXmax();
            tit = h->GetXaxis()->GetTitle();
        }
};

template <> inline
void HistoSvc::Book<TH1D>(TString name, TString v1, TString v2, TString v3){
    int xN = 0;
    double xlo = 0;
    double xhi = 0;
    TString xtit = "";
    getBinning(v1,xN,xlo,xhi,xtit);
    TString title = name;
    title += ";";
    title += xtit;
    title += ";";
    TH1D* h = new TH1D(name,title,xN,xlo,xhi);
    h->Sumw2();
    l_h.Add(h);
}

template <> inline
void HistoSvc::Book<TH2D>(TString name, TString v1, TString v2, TString v3){
    int xN,yN,zN ;
    double xlo,ylo,zlo;
    double xhi,yhi,zhi;
    TString xtit,ytit,ztit;
    getBinning(v1,xN,xlo,xhi,xtit);
    getBinning(v2,yN,ylo,yhi,ytit);
    TString title = name;
    title += ";";
    title += xtit;
    title += ";";
    title += ytit;
    TH2D* h = new TH2D(name,title,xN,xlo,xhi,yN,ylo,yhi);
    h->Sumw2(); h->SetOption("COLZ");
    l_h.Add(h);
}

template <> inline
void HistoSvc::Book<TProfile>(TString name, TString v1, TString v2, TString v3){
    int xN,yN,zN ;
    double xlo,ylo,zlo;
    double xhi,yhi,zhi;
    TString xtit,ytit,ztit;
    getBinning(v1,xN,xlo,xhi,xtit);
    getBinning(v2,yN,ylo,yhi,ytit);
    TString title = name;
    title += ";";
    title += xtit;
    title += ";";
    title += ytit;
    TProfile* h = new TProfile(name,title,xN,xlo,xhi);
    l_h.Add(h);
}

template <> inline
void HistoSvc::Book<TProfile2D>(TString name, TString v1, TString v2, TString v3){
    int xN,yN,zN ;
    double xlo,ylo,zlo;
    double xhi,yhi,zhi;
    TString xtit,ytit,ztit;
    getBinning(v1,xN,xlo,xhi,xtit);
    getBinning(v2,yN,ylo,yhi,ytit);
    getBinning(v3,zN,zlo,zhi,ztit);
    TString title = name;
    title += ";";
    title += xtit;
    title += ";";
    title += ytit;
    title += ";";
    title += ztit;
    TProfile2D* h = new TProfile2D(name,title,xN,xlo,xhi,yN,ylo,yhi);
    l_h.Add(h);
}

template<class T> inline
T * HistoSvc::Get(TString name){
    T * h = (T*) l_h(name.Data());
    return h;
}

template <> inline
void HistoSvc::Fill<TH1D>(TString name, double v1, double v2, double v3, double v4){
    TH1D * h = (TH1D*)l_h(name.Data());
    if (!h) throw invalid_argument(Form("There is no histogram called `%s`", name.Data() ));
    h->Fill(v1,wgt);
}
template <> inline
void HistoSvc::Fill<TH2D>(TString name, double v1, double v2, double v3, double v4){
    TH2D * h = (TH2D*)l_h(name.Data());
    if (!h) throw invalid_argument(Form("There is no histogram called `%s`", name.Data() ));
    h->Fill(v1,v2,wgt);
}
template <> inline
void HistoSvc::Fill<TProfile>(TString name, double v1, double v2, double v3, double v4){
    TProfile * h = (TProfile*)l_h(name.Data());
    if (!h) throw invalid_argument(Form("There is no histogram called `%s`", name.Data() ));
    h->Fill(v1,v2,wgt);
}
template <> inline
void HistoSvc::Fill<TProfile2D>(TString name, double v1, double v2, double v3, double v4){
    TProfile2D * h = (TProfile2D*)l_h(name.Data());
    if (!h) throw invalid_argument(Form("There is no histogram called `%s`", name.Data() ));
    h->Fill(v1,v2,v3,wgt);
}

#endif //_HistoSvc_H
