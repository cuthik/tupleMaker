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
  
  Output(TTree * t=0);
  
  ~Output() {
    if(!_written) Write();
    if( _cleanup) delete _tree;
  }
  
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
    Float_t        pdf_wgts[51];
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
    
  void Fill() { _tree->Fill();}
  void Write() { 
    if( !_written) _tree->Write();
    _written = true;
  }
  
  void AddParticle( int id, float px, float py, float pz, float E, int origin)
  {
    double twopi = 2.*acos(-1.);
    int pdgid = _pidDB->ConvertIsajetToPdg(id);
    int npart = _ana.npart;
    _ana.pistable[npart] = origin;
    _ana.pE[npart] = E;
    _ana.ppx[npart] = px;
    _ana.ppy[npart] = py;
    _ana.ppz[npart] = pz;
    float phi, eta;
    float pt = px*px + py*py;
    pt = sqrt(pt);
    _ana.ppt[npart] = pt;
    float p = px*px + py*py + pz*pz;
    
    p = sqrt(p);
    if( fabs(px) > 0.000001 || fabs(py) > 0.0000001) {
      phi = atan2( py, px);
      if( phi < 0. ) phi += twopi;
      
      eta = 0.5 * log( (p+pz+0.0000001) / (p-pz+0.0000001) );
    } else {
      phi = -1.;
      if( pz > 0. )
	eta =  9999.;
      else
	eta = -9999.;
    }
    _ana.peta[npart] = eta;
    _ana.pphi[npart] = phi;
    _ana.ppid[npart] = pdgid;
    _ana.ppvtx[npart] = 0;
    ++_ana.npart;	
    if( _ana.nvtx == 0) {
      _ana.nvtx = 1;
      _ana.vx[0] = 0.;
      _ana.vy[0] = 0.;
      _ana.vz[0] = 0.;
    }
    if( std::abs(id) == 11 || std::abs(id) == 13 || std::abs(id) ==15 ) // neutrino
      {
	// update the met info 
	if( _met.nmetg == 0) {
	  _met.nmetg = 1;
	  _met.metg[0] = 0.;
	  _met.metxg[0] = 0.;
	  _met.metyg[0] = 0.;
	}
	_met.nmetg = 1;
	_met.metxg[0] += px;
	_met.metyg[0] += py;
	_met.metg[0] = sqrt( _met.metxg[0] * _met.metxg[0] + _met.metyg[0]*_met.metyg[0]);
      }
  }
  
 
  void NewEvent( int evn, double evt_wt , int run , 
                float vx , float vy , float vz , 
                float Q2 , float x1 , float x2 ,
                float flav1 , float flav2 , std::vector<float> pdf_wgts ) 
  {
    Reset();
    _ana.nvtx = 1;
    _ana.nevtp = 1;
    _ana.evnum[0] = evn;
    _ana.evwt[0] = evt_wt;
    _ana.evrun[0] = run;
    _ana.vx[0] = vx;
    _ana.vy[0] = vy;
    _ana.vz[0] = vz;
#ifdef __USE_PDFS__
    _ana.evflav1[0] = flav1;
    _ana.evflav2[0] = flav2;
    _ana.evqsq[0] = Q2;
    _ana.evx1[0] = x1;
    _ana.evx2[0] = x2;
#endif
#ifdef __USE_PDFS_RESBOS__
    if( pdf_wgts.size() >= 51 )
    {
      for( int i = 0 ; i < 51 ; i++ )
        _ana.pdf_wgts[i] = pdf_wgts[i];
    }
#endif
    }

  void NewEvent( int evn, double evt_wt , int run , float vx , float vy , float vz , 
		 double evt_wt_PDF[], int length) 
  {
    Reset();
    _ana.nvtx = 1;
    _ana.nevtp = 51;
    _ana.evnum[0] = evn;
    _ana.evwt[0] = evt_wt;

    // for different PDF sets
    for(int i=0; i<length; i++) {
#ifdef __USE_PDFS_RESBOS__
        _ana.pdf_wgts[i] = evt_wt_PDF[i];
#else
        _ana.evwt[i+1] = evt_wt_PDF[i];
#endif
    }

    _ana.evrun[0] = run;
    _ana.vx[0] = vx;
    _ana.vy[0] = vy;
    _ana.vz[0] = vz;
  }

  void Reset() 	
  {
    _written = false;
    _ana.nevtp = 0;
    _ana.npart=0;
    _ana.nvtx = 0;
    _em.nelg = 0;
    _em.nels=0;
    _em.nphg = 0;
    _em.nphs = 0;
    _met.nmetg = 0;
    _met.nmets = 0;
    _vtx.nvtxg = 0;
    _vtx.nvtxs = 0;
  }			
  
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

