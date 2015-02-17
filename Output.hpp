#ifndef Output_h
#define Output_h

#include "TTree.h"
#include "TDatabasePDG.h"
#include <string>
#include <vector>
#include "TBranch.h"
#include "TLeafI.h"
#include "math.h"

class Output {
    public:

        Output(TString filename);
        Output(TTree * t=0);

        ~Output() ;

        static const int nPDF=150;

        TTree * Tree() { return _tree;}
        void fixLeafOffsets( TBranch * b);

        const std::string ana_form;  
        struct anaBlock {
            Int_t           nevtp;
            Int_t           npart;
            Int_t           nvtx;
            Float_t         evnum[1000];   //[nevtp]
            Float_t         evrun[1000];   //[nevtp]
            Float_t         evwt[1000];   //[nevtp]
            Float_t         evxs[1000];   //[nevtp]
            Float_t         pE[1000];   //[npart]
            Float_t         pcid[1000];   //[npart]
            Float_t         pcnum[1000];   //[npart]
            Float_t         pdvtx[1000];   //[npart]
            Float_t         peta[1000];   //[npart]
            Float_t         pidx[1000];   //[npart]
            Float_t         pistable[1000];   //[npart]
            Float_t         pphi[1000];   //[npart]
            Float_t         ppid[1000];   //[npart]
            Float_t         ppt[1000];   //[npart]
            Float_t         ppvtx[1000];   //[npart]
            Float_t         ppx[1000];   //[npart]
            Float_t         ppy[1000];   //[npart]
            Float_t         ppz[1000];   //[npart]
            Float_t         vcid[1000];   //[nvtx]
            Float_t         vcnum[1000];   //[nvtx]
            Float_t         vct[1000];   //[nvtx]
            Float_t         vidx[1000];   //[nvtx]
            Float_t         visdisp[1000];   //[nvtx]
            Float_t         vpprt[1000];   //[nvtx]
            Float_t         vx[1000];   //[nvtx]
            Float_t         vy[1000];   //[nvtx]
            Float_t         vz[1000];   //[nvtx] 
#ifdef __USE_PDFS__
            Float_t        evflav1[1000]; 
            Float_t        evflav2[1000]; 
            Float_t        evqsq[1000];
            Float_t        evx1[1000];
            Float_t        evx2[1000];
#endif
#ifdef __USE_PDFS_RESBOS__
            Float_t        pdf_wgts[nPDF];
#endif
        };

        const std::string em_form;
        struct emBlock {
            Int_t           nelg;
            Int_t           nels;
            Int_t           nphg;
            Int_t           nphs;
            Float_t         elcalphis[120];   //[nels]
            Float_t         eleg[120];   //[nelg]
            Float_t         elelmergedEg[120];   //[nelg]
            Float_t         eles[120];   //[nels]
            Float_t         eletads[120];   //[nels]
            Float_t         eletag[120];   //[nelg]
            Float_t         eletas[120];   //[nels]
            Float_t         elfid[120];   //[nelg]
            Int_t           elhastrack[120];   //[nels]
            Float_t         eliso[120];   //[nels]
            Int_t           elmergedg[120];   //[nelg]
            Int_t           elpasshmtx[120];   //[nels]
            Int_t           elpassid1011[120];   //[nels]
            Float_t         elphig[120];   //[nelg]
            Float_t         elphis[120];   //[nels]
            Float_t         elphmergedEg[120];   //[nelg]
            Int_t           elpntg[120];   //[nels]
            Int_t           elpnts[120];   //[nelg]
            Float_t         elptg[120];   //[nelg]
            Float_t         elpts[120];   //[nels]
            Int_t           elpttr[120];   //[nels]
            Float_t         pheg[120];   //[nphg]
            Float_t         phes[120];   //[nphs]
            Float_t         phetads[120];   //[nphs]
            Float_t         phetag[120];   //[nphg]
            Float_t         phetas[120];   //[nphs]
            Float_t         phfid[120];   //[nphg]
            Int_t           phhastrack[120];   //[nels]
            Float_t         phiso[120];   //[nphs]
            Int_t           phpasshmtx[120];   //[nels]
            Float_t         phphig[120];   //[nphg]
            Float_t         phphis[120];   //[nphs]
            Int_t           phpntg[120];   //[nphs]
            Int_t           phpnts[120];   //[nphg]
            Float_t         phptg[120];   //[nphg]
            Float_t         phpts[120];   //[nphs]
        };

        const std::string met_form;    
        struct metBlock {
            Int_t           nmetg;
            Int_t           nmets;
            Float_t         metg[200];   //[nmetg]
            Float_t         metphig[200];   //[nmetg]
            Float_t         metphis[200];   //[nmets]
            Float_t         mets[200];   //[nmets]
            Float_t         metxg[200];   //[nmetg]
            Float_t         metxs[200];   //[nmets]
            Float_t         metyg[200];   //[nmetg]
            Float_t         metys[200];   //[nmets]
            Float_t         scalarg[200];   //[nmetg]
            Float_t         scalars[200];   //[nmets]
        };

        const std::string vtx_form;
        struct vtxBlock {
            Int_t           nvtxg;
            Int_t           nvtxs;
            Float_t         vtnds[150];   //[nvtxs]
            Float_t         vtxxs[150];   //[nvtxs]
            Float_t         vtxys[150];   //[nvtxs]
            Float_t         vtxzg[150];   //[nvtxg]
            Float_t         vtxzs[150];   //[nvtxs]
        };

        void Fill() ;
        void Write() ;

        void AddParticle( int id, float px, float py, float pz, float E, int origin);
        void AddParticlePDGID( int id, float px, float py, float pz, float E, int origin);


        void NewEvent( int evn, double evt_wt , int run , 
                float vx , float vy , float vz , 
                float Q2 , float x1 , float x2 ,
                float flav1 , float flav2 , std::vector<float> pdf_wgts );

        void NewEvent( int evn, double evt_wt , int run , float vx , float vy , float vz , 
                double evt_wt_PDF[], int length);

        void Reset();


        anaBlock _ana;
        emBlock _em;
        metBlock _met;
        vtxBlock _vtx;

        TTree * _tree;

        bool _cleanup;
        bool _written;
        TDatabasePDG * _pidDB;
};


#endif // Output_h

