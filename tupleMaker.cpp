#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "Output.hpp"
#include "TFile.h"
#include "TRandom2.h"

using namespace std;

int main( int argc, char * argv[])
{
  if( argc < 3) {
    std::cout << "usage : tupleMaker input-filenames output-filename <Unweight>" << std::endl;
    return 1;
  }
  std::string lastarg  = std::string( argv[argc-1] );
  
  int lastargisfilename = 0;
  if (lastarg.length() > 1 ){
    lastargisfilename =1;
  }
  if (!lastargisfilename){
    if( argc < 5) {
      std::cout << "usage : tupleMaker input-filenames output-filename <Unweight> <PDF reweight>" << std::endl;
      return 1;
    }
  }
  
  //
  // read in some control parameters
  //
  bool Unweight = 0;
  bool PDF_reweight = 0;
  if (!lastargisfilename){
    Unweight = atoi( argv[argc-2] );
    PDF_reweight = atoi( argv[argc-1] );
    std::cout<<"Unweight = "<<Unweight<<" PDF reweight = "<<PDF_reweight<<std::endl;
    if(PDF_reweight==1) {
      std::cout<<"========================================================================================================"<<std::endl;
      std::cout<<"Need to make sure you have right PDF text files resbos_P/resbos_weights_PDF_*.dat"<<std::endl;
      std::cout<<"========================================================================================================"<<std::endl;
    }

  }
  
  // get input and output file names
  std::vector<std::string> in_files;
  int nInputs = argc-(2-lastargisfilename);
  if(PDF_reweight) nInputs = argc-(3-lastargisfilename);

  for( int i = 1; i < nInputs;++i) {
    in_files.push_back( std::string( argv[i] ) );
    std::cout<<"Input file: "<<std::string(argv[i])<<std::endl;
  }

  int nOutput = argc-(2-lastargisfilename);
  if(PDF_reweight) nOutput = argc-(3-lastargisfilename);
  std::string out_filename = std::string( argv[nOutput] );
  std::cout<<"Output file: "<<out_filename<<std::endl;

  // setup outputfile
  TFile * of = new TFile(out_filename.c_str(),"RECREATE");
  if( of == 0) {
    std::cout << "Couldn't open output file : " << out_filename << std::endl;
    return 1;
  }
  
  Output * output = new Output();
  output->Reset();
  // initialise random number generator
  TRandom2 random;
  random.SetSeed(19742005);
  // loop over files
  Double_t StandardWeight =0;
  int nfiles = in_files.size();

  if (Unweight){
    std::cout << "All events will have weight one" << std::endl;
    std::cout << "Multiple files are merged using pass/fail based on the weight of events in the first file -- Standard Weight as shown below" << std::endl;
  }
  
  if (!Unweight)
    std::cout << "Unweighted events, Weights from generator are maintained"<< std::endl;

  //  
  // read event weights for different PDF files
  //
  std::ifstream f_weights[45];
  vector<float> pdfwgts;
  pdfwgts.resize(45);
  if(PDF_reweight) {
    char name[50];    
    for(int i=0; i<44; i++) {
      sprintf(name, "%s%d%s", "resbos_P/resbos_weights_PDF_", i+1, ".dat");
      f_weights[i].open(name, std::ios::in);
      if( !f_weights[i] ) std::cout<<"Could not find the weight file "<<name<<std::endl;
    }
  }

#ifdef __USE_PDFS_RESBOS__
  Unweight = true;
  ///
  {
    std::ifstream tmp_file;
    tmp_file.open( "weights.dat" , std::ios::in );
    if( tmp_file )
    {
      int num ; double wgt ;
      while( tmp_file >> num >> wgt )
      {
        if( wgt > StandardWeight )
          StandardWeight = wgt;
      }
      tmp_file.close();
    }
  }
  ///
  for( int i = 0 ; i < 45 ; i++ )
  {
    TString name;
    name.Form( "weights_%02i.hep" , i );
    f_weights[i].open( name.Data() , std::ios::in );
    if( !f_weights[i] ) std::cout<<"Could not find the weight file "<<name<<std::endl;
  }
#endif

  for(int i =0; i < nfiles; ++i) { 
    // open file
    std::cout << "Processing file: " << in_files[i].c_str() << std::endl;
    std::ifstream f((in_files[i]).c_str());
    if( !f ) {
      std::cout << "couldn't open file : " << in_files[i] << std::endl;
      return 1;
    }
    
    bool finished_file = false;
    // this loops over events
    while( ! finished_file ) {
      int evn;
      double evt_wt;
      double Q2,that,uhat,x1,x2,flav1,flav2 ;
#ifdef __USE_PDFS__
      f >> evn >> evt_wt >> Q2 >> that >> uhat >> x1 >> x2 >> flav1 >> flav2;
#else
      f >>  evn >> evt_wt;
#endif

      if(evn % 100000==0) std::cout<<"Processing event: "<<evn<<std::endl;

      if( evn == 0 ) {
	finished_file = true;
	continue;
      }

#ifdef __USE_PDFS_RESBOS__
      for( int j = 0 ; j<45 ; j++ )
      {
        double evn_tmp , wgt_tmp;
        f_weights[j] >> evn_tmp >> wgt_tmp;
        if( evn_tmp != evn )
        {
          cout << " WRONG EVENT NUMBER!!!! " << j << " " << evn_tmp << " " << evn << endl;
          return 1;
        }
        if( j == 0 ) evt_wt = wgt_tmp;
        if( evt_wt > 0 && pdfwgts[j] > 0 )
          pdfwgts[j] = wgt_tmp / evt_wt;
        else
          pdfwgts[j] = 1.0;
        if( pdfwgts[j] < 0 )
          pdfwgts[j] = 0.0;
      }
#endif
     
      float vx,vy,vz;
      f>> vx >> vy >>vz;
      
      bool finished_particles = false;
      if( !f.eof())
      {
	if (StandardWeight == 0. && Unweight){
	  StandardWeight = evt_wt;
	  std::cout << "Standard Weight  = " <<  StandardWeight << std::endl;
	}
      }
      else
      {
        finished_file = true;
      }

      Bool_t keeper = kTRUE;
      if (Unweight){
	Double_t weight_ratio=  TMath::Abs(evt_wt/StandardWeight);
	if (weight_ratio > 1.){
	  std::cout << "This event has a weight = " << weight_ratio <<" times the value of Standard Weight" << std::endl;
	}	  
	if (random.Rndm() > weight_ratio)
	  keeper = kFALSE;
      }
      if (Unweight)
      	evt_wt = TMath::Abs(evt_wt) / evt_wt ;
      
      //
      // read event weight for each PDF set
      //
      double weight_PDF[44];
      if(PDF_reweight) {
	for(int j=0; j<44; j++) (f_weights[j]) >> weight_PDF[i];
      }
      
      // save PDF information into the root file
      //      output->setPDFWeights(weight_PDF);
      if (keeper) {
	if(PDF_reweight) output->NewEvent( evn, evt_wt, 0 , vx,vy,vz, weight_PDF, 44);
	else 
   {
     output->NewEvent( evn, evt_wt, 0 , vx,vy,vz , Q2 , x1 , x2 , flav1 , flav2 , pdfwgts );
   }
      }
      
      //bool doFlip = false;
      //doFlip = ( random.Rndm() > 0.5); 
      
      // this loops over particles in an event
      while( ! finished_particles && !f.eof()) {
	int id;
	f >> id;
	if( id == 0 )
	  finished_particles = true;
	else {
	  float px,py,pz,E;
	  f>> px >> py >> pz >> E;
	  int origin, udk;
	  f >> origin >> udk;
	  //std::cout << id << " " << px << " " << py << " " << pz << " " << E << std::endl;
	  
	  // do the occaisional CP inversion
	  // flip W+ -> W-, e+->e-, nu -> nu-bar and invert all momenta
	  //if( doFlip ) {
	  // id is isajet ids
	  //  if( id != 10 ) id = -id;  // photon == anti-photon
	  //  px = -px;
	  //  py = -py;
	  //  pz = -pz;
	  //}
	  if(keeper)
	   output->AddParticle( id, px,py,pz,E,origin);
	}
      }
      if (keeper){
	output->Fill();	
      }
      if( f.eof()) finished_file = true;
    }
    f.close();
  }
  
  output->Write();
  of = output->Tree()->GetCurrentFile();
  of->Close();
  return 0;
};

