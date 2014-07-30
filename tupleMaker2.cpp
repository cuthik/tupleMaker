#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "Output.hpp"
#include "TFile.h"
#include "TRandom2.h"

#define ONE_WEIGHT_PER_FILE

using namespace std;

int main( int argc, char * argv[]){

    const int MAXPDF=51;

    // setup the seed and center
    TString center=argv[1];
    TString seed=argv[2];
    TString non_center= center.Data();
    if(argc>3) non_center=argv[3];

    TString in_file, out_file;
    in_file  .Form("resbos_%s_%s.hep"  ,center     .Data(),seed.Data());
    out_file .Form("resbos_%s_%s.root" ,non_center .Data(),seed.Data());

    cout << "input file: " << in_file << endl;
    cout << "output file: " << out_file << endl;

    TFile * of = TFile::Open(out_file.Data(),"RECREATE");
    if(!of){
        cout << "can not create outfile" << endl;
        return 1;
    }

#ifndef ONE_WEIGHT_PER_FILE
    // open weight files -- all at once
    cout << "openning reweighting files " << endl;
    std::fstream f_weights[MAXPDF];
    vector<float> pdfwgts;
    pdfwgts.resize(MAXPDF);
    TString tmp;
    if (center.Contains("central")){ // theo unc
        tmp.Form("weight_%s_%s.txt", "central" , seed.Data() ); f_weights[0].open( tmp.Data() , std::ios::in ); cout << tmp << endl;
        tmp.Form("weight_%s_%s.txt", "high_NS" , seed.Data() ); f_weights[1].open( tmp.Data() , std::ios::in ); cout << tmp << endl;
        tmp.Form("weight_%s_%s.txt", "low_NS"  , seed.Data() ); f_weights[2].open( tmp.Data() , std::ios::in ); cout << tmp << endl;
        tmp.Form("weight_%s_%s.txt", "high_S"  , seed.Data() ); f_weights[3].open( tmp.Data() , std::ios::in ); cout << tmp << endl;
        tmp.Form("weight_%s_%s.txt", "low_S"   , seed.Data() ); f_weights[4].open( tmp.Data() , std::ios::in ); cout << tmp << endl;
        for (int i=5; i<MAXPDF; i++){
            tmp.Form("weight_%s_%s.txt", "central" , seed.Data() ); f_weights[i].open( tmp.Data() , std::ios::in ); cout << tmp << endl;
        } 
    } else { //PDFs
        for (int i=0; i<MAXPDF; i++){
            tmp.Form("weight_%02i_%s.txt", i , seed.Data() ); f_weights[i].open( tmp.Data() , std::ios::in ); cout << tmp << endl;
        }
    }
#else
    // open weight files -- load only one
    cout << "openning reweighting file" << endl;
    std::fstream f_weights;
    TString tmp;
    tmp.Form( "weight_%s_%s.txt", non_center.Data() , seed.Data() );
    f_weights.open( tmp.Data(), std::ios::in );

    // get index of sample
    size_t i_non_center;
    if (center.Contains("central")){ // theo unc
        tmp.Form("%s", "central" ); if (tmp.CompareTo(non_center,TString::kIgnoreCase)==0) {i_non_center=0; } ;
        tmp.Form("%s", "high_NS" ); if (tmp.CompareTo(non_center,TString::kIgnoreCase)==0) {i_non_center=1; } ;
        tmp.Form("%s", "low_NS"  ); if (tmp.CompareTo(non_center,TString::kIgnoreCase)==0) {i_non_center=2; } ;
        tmp.Form("%s", "high_S"  ); if (tmp.CompareTo(non_center,TString::kIgnoreCase)==0) {i_non_center=3; } ;
        tmp.Form("%s", "low_S"   ); if (tmp.CompareTo(non_center,TString::kIgnoreCase)==0) {i_non_center=4; } ;
    } else { //PDFs
        for (int i=0; i<MAXPDF; i++){
            tmp.Form("%02i", i ); if (tmp.CompareTo(non_center,TString::kIgnoreCase)==0) {i_non_center=i; } ;
        }
    }
    cout << "name of sample: " << non_center << " " << i_non_center << endl;
#endif



    cout << "reading input tree" << endl;
    // process file
    Output * output = new Output();
    output->Reset();
    // initialise random number generator
    TRandom2 random;
    random.SetSeed(19742005);

    std::ifstream f ( in_file.Data());
    bool finished_file = false;
    // this loops over events
    while( ! finished_file ) {
        int evn;
        double evt_wt;
        double Q2,that,uhat,x1,x2,flav1,flav2 ;
        f >>  evn >> evt_wt;

        if(evn % 100000==0) std::cout<<"Processing event: "<<evn<<std::endl;

        if( evn == 0 ) {
            finished_file = true;
            continue;
        }

#ifndef ONE_WEIGHT_PER_FILE
        // load event weights -- all at once
        double weight_PDF[MAXPDF-1];
        for( int j = 0 ; j<MAXPDF ; j++ )
        {
            double evn_tmp , wgt_tmp;
            f_weights[j] >> evn_tmp >> wgt_tmp;
            if( evn_tmp != evn )
            {
                cout << " WRONG EVENT NUMBER!!!! " << j << " " << evn_tmp << " " << evn << endl;
                goto closeFile;
                return 1;
            }
            if( j == 0 ) evt_wt = wgt_tmp;
            if( evt_wt > 0 ) pdfwgts[j] = wgt_tmp / evt_wt;
            else pdfwgts[j] = 1.0;
            if( pdfwgts[j] < 0 ) pdfwgts[j] = 0.0;

            if(j>0)weight_PDF[j-1]=pdfwgts[j];
        }
#else
        // load event weights -- only one
        double weight_PDF[MAXPDF-1];
        for (int j=0; j < MAXPDF-1; j++){
            weight_PDF[j]=0;
        }
        if (i_non_center!=0){ // load non central weight
            double evn_tmp , wgt_tmp;
            f_weights >> evn_tmp >> wgt_tmp;
            if( evn_tmp != evn )
            {
                cout << " WRONG EVENT NUMBER!!!! " <<  evn_tmp << " " << evn << endl;
                goto closeFile;
                return 1;
            }
            if (wgt_tmp >0 ) weight_PDF[i_non_center-1]= (evt_wt>0) ? wgt_tmp/evt_wt : 1.;
        }
#endif


        // vertex 
        float vx,vy,vz;
        f>> vx >> vy >>vz;

        bool finished_particles = false;
        finished_file = f.eof();

        output->NewEvent( evn, evt_wt, 0 , vx,vy,vz, weight_PDF, MAXPDF-1);

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
                output->AddParticle( id, px,py,pz,E,origin);
            }
        }
        output->Fill();	
        if( f.eof()) {
            finished_file = true;
        }
    }

closeFile:
    f.close();
    output->Write();
    of = output->Tree()->GetCurrentFile();
    of->Close();
    return 0;
}
