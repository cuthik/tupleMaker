
#include "Output.hpp"
#include <iostream>
#include "TTree.h"
#include "TDatabasePDG.h"
#include <string>
#include "TBranch.h"
#include "TLeafI.h"





Output::Output(TTree * t) : _tree(t), _written(false),
#ifdef __USE_PDFS__
#ifdef __USE_PDFS_RESBOS__
            ana_form("nevtp/I:npart/I:nvtx/I:evnum[nevtp]/F:evrun[nevtp]/F:evwt[nevtp]/F:evxs[nevtp]/F:pE[npart]/F:pcid[npart]/F:pcnum[npart]/F:pdvtx[npart]/F:peta[npart]/F:pidx[npart]/F:pistable[npart]/F:pphi[npart]/F:ppid[npart]/F:ppt[npart]/F:ppvtx[npart]/F:ppx[npart]/F:ppy[npart]/F:ppz[npart]/F:vcid[nvtx]/F:vcnum[nvtx]/F:vct[nvtx]/F:vidx[nvtx]/F:visdisp[nvtx]/F:vpprt[nvtx]/F:vx[nvtx]/F:vy[nvtx]/F:vz[nvtx]/F:evflav1[nevtp]/F:evflav2[nevtp]/F:evqsq[nevtp]/F:evx1[nevtp]/F:evx2[nevtp]/F:pdf_wgts[150]/F"),
#else
            ana_form("nevtp/I:npart/I:nvtx/I:evnum[nevtp]/F:evrun[nevtp]/F:evwt[nevtp]/F:evxs[nevtp]/F:pE[npart]/F:pcid[npart]/F:pcnum[npart]/F:pdvtx[npart]/F:peta[npart]/F:pidx[npart]/F:pistable[npart]/F:pphi[npart]/F:ppid[npart]/F:ppt[npart]/F:ppvtx[npart]/F:ppx[npart]/F:ppy[npart]/F:ppz[npart]/F:vcid[nvtx]/F:vcnum[nvtx]/F:vct[nvtx]/F:vidx[nvtx]/F:visdisp[nvtx]/F:vpprt[nvtx]/F:vx[nvtx]/F:vy[nvtx]/F:vz[nvtx]/F:evflav1[nevtp]/F:evflav2[nevtp]/F:evqsq[nevtp]/F:evx1[nevtp]/F:evx2[nevtp]/F"),
#endif
#elif __USE_PDFS_RESBOS__
            ana_form("nevtp/I:npart/I:nvtx/I:evnum[nevtp]/F:evrun[nevtp]/F:evwt[nevtp]/F:evxs[nevtp]/F:pE[npart]/F:pcid[npart]/F:pcnum[npart]/F:pdvtx[npart]/F:peta[npart]/F:pidx[npart]/F:pistable[npart]/F:pphi[npart]/F:ppid[npart]/F:ppt[npart]/F:ppvtx[npart]/F:ppx[npart]/F:ppy[npart]/F:ppz[npart]/F:vcid[nvtx]/F:vcnum[nvtx]/F:vct[nvtx]/F:vidx[nvtx]/F:visdisp[nvtx]/F:vpprt[nvtx]/F:vx[nvtx]/F:vy[nvtx]/F:vz[nvtx]/F:pdf_wgts[45]/F"),
#else
			    ana_form("nevtp/I:npart/I:nvtx/I:evnum[nevtp]/F:evrun[nevtp]/F:evwt[nevtp]/F:evxs[nevtp]/F:pE[npart]/F:pcid[npart]/F:pcnum[npart]/F:pdvtx[npart]/F:peta[npart]/F:pidx[npart]/F:pistable[npart]/F:pphi[npart]/F:ppid[npart]/F:ppt[npart]/F:ppvtx[npart]/F:ppx[npart]/F:ppy[npart]/F:ppz[npart]/F:vcid[nvtx]/F:vcnum[nvtx]/F:vct[nvtx]/F:vidx[nvtx]/F:visdisp[nvtx]/F:vpprt[nvtx]/F:vx[nvtx]/F:vy[nvtx]/F:vz[nvtx]/F"),
#endif
			    em_form( "nelg/I:nels/I:nphg/I:nphs/I:elcalphis[nels]/F:eleg[nelg]/F:elelmergedEg[nelg]/F:eles[nelg]/F:eletads[nels]/F:eletag[nelg]/F:eletas[nels]/F:elfid[nelg]/F:elhastrack[nels]/I:eliso[nels]/F:elmergedg[nelg]/I:elpasshmtx[nels]/I:elpassid1011[nels]/I:elphig[nelg]/F:elphis[nels]/F:elphmergedEg[nelg]/F:elpntg[nels]/I:elpnts[nelg]/I:elptg[nelg]/F:elpts[nels]/F:elpttr[nels]/I:pheg[nphg]/F:phes[nphs]/F:phetads[nphs]/F:phetag[nphg]/F:phetas[nphs]/F:phfid[nphg]/F:phhastrack[nels]/I:phiso[nphs]/F:phpasshmtx[nels]/I:phphig[nphg]/F:phphis[nphs]/F:phpntg[nphs]/I:phpnts[nphg]/I:phptg[nphg]/F:phpts[nphs]/F"),			    
			    met_form("nmetg/I:nmets/I:metg[nmetg]/F:metphig[nmetg]/F:metphis[nmets]/F:mets[nmets]/F:metxg[nmetg]/F:metxs[nmets]/F:metyg[nmetg]/F:metys[nmets]/F:scalarg[nmetg]/F:scalars[nmets]/F"),
			    vtx_form("nvtxg/I:nvtxs/I:vtnds[nvtxs]/F:vtxxs[nvtxs]/F:vtxys[nvtxs]/F:vtxzg[nvtxg]/F:vtxzs[nvtxs]/F")
{ 
  if( _tree==0 ) {
    _tree = new TTree("Global","Fake pmcs output");
    _cleanup = true;
  } else {
    _cleanup = false;
  }
  
  TLeafI * index = 0;
  
  TBranch * b_ana =
    _tree->Branch("pmcs_ana", &_ana.nevtp, ana_form.c_str());
  
  index = (TLeafI*)b_ana->GetLeaf("nevtp");
  index->SetMaximum(1000);
  index = (TLeafI*)b_ana->GetLeaf("npart");
  index->SetMaximum(1000);
  index = (TLeafI*)b_ana->GetLeaf("nvtx");
  index->SetMaximum(1000);
  
  // root sucks. While it allows you to have a variable length array with the index in the same branch it doesn't
  // actually allow you to specify the span of the index variables at the appropriate moment thus all teh address 
  // offsets for the leaves are crap. Luckily in the HepTuple root code they have a workaround which I include a 
  // slightly modified version of here to get the leaf address offsets correct.
  
  fixLeafOffsets( b_ana);
  
  TBranch * b_em = 
    _tree->Branch("pmcs_em", &_em.nelg, em_form.c_str());
  
  index = (TLeafI*)b_em->GetLeaf("nelg");
  index->SetMaximum(120);
  index = (TLeafI*)b_em->GetLeaf("nels");
  index->SetMaximum(120);
  index = (TLeafI*)b_em->GetLeaf("nphg");
  index->SetMaximum(120);
  index = (TLeafI*)b_em->GetLeaf("nphs");
  index->SetMaximum(120);
  
  fixLeafOffsets( b_em );
  
  TBranch * b_met = 
    _tree->Branch("pmcs_met", &_met.nmetg, met_form.c_str());
  index = (TLeafI*)b_met->GetLeaf("nmetg");
  index->SetMaximum(200);
  index = (TLeafI*)b_met->GetLeaf("nmets");
  index->SetMaximum(200);
  
  fixLeafOffsets( b_met );
  
  TBranch * b_vtx =
    _tree->Branch("pmcs_vtx", &_vtx.nvtxg, vtx_form.c_str());
  
  index = (TLeafI*)b_vtx->GetLeaf("nvtxg");
  index->SetMaximum(150);
  index = (TLeafI*)b_vtx->GetLeaf("nvtxs");
  index->SetMaximum(150);

  fixLeafOffsets( b_vtx );
  
  _pidDB = TDatabasePDG::Instance();  
}

Output::~Output() {
    if(!_written) Write();
    if( _cleanup) delete _tree;
}

void Output::Fill() { _tree->Fill();}
void Write() { 
    if( !_written) _tree->Write();
    _written = true;
}

void Output::fixLeafOffsets( TBranch * b)
{
  // recalculates the addresses for all the leaves.
  // 
  // when constructing a branch with containing a variable length array with the index
  // variable in the same branch it is not possible to specify the span of the index variable
  // This span defaults to zero. When the addresses are asigned to the various leaves in the branch
  // it calculates the size of the particular leaf (variable length array) in the buffer by looking 
  // at the span of the index variable - 0 in this case! using the method TLeaf::GetLen(). 
  // The following code shoudl be applied to the branch after the spans of the index variables have been
  // specified manually using the TLeaf::SetMaximum method. This time the GetLen method calculates the correct offset.

  TObjArray * leaves = b->GetListOfLeaves();
  char * addr = b->GetAddress();
  int offset = 0;
  int nleaves = leaves->GetEntriesFast();
  
  // loop over the leaves:
  for( int i =0; i < nleaves; ++i) {
    TLeaf * leaf = (TLeaf *)leaves->UncheckedAt(i);
    leaf->SetAddress( addr + offset );
    int oldOffset = leaf->GetOffset();
    leaf->SetOffset( offset );
    //std::cout << " offset changed from : " << oldOffset << " to " << offset << std::endl;
    TLeaf * index = leaf->GetLeafCount();
    int nelements = 1;
    if( index ) {
      nelements = index->GetMaximum(); // deal with variable length arrays
    } else {
      nelements = leaf->GetLenStatic(); // deal with single variables and fixed length arrays
    }
    offset += leaf->GetLenType() * nelements;
  }
}
