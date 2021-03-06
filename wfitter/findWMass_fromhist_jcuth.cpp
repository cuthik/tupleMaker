#ifndef wfitter_jcuth
#define wfitter_jcuth
/**
 * @file findWMass_fromhist_jcuth.cpp
 * Template fit of the wmass histogram
 *
 * @brief Because I can not orient in code I reimplement it into object.
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2014-08-21
 */



#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <TRegexp.h>
#include <TApplication.h>
#include <vector>
#include <TRandom.h>
#include <TMinuit.h>
#include <TMath.h>

#include "wzfitter/genDist.hpp"
#include "wzfitter/genGaussMean.hpp" 
#include "wzfitter/genGauss.hpp"
#include "wzfitter/gen1Dspline.hpp"
#include "wzfitter/wzfitter.hpp"
#include "wzfitter/TNHist.hpp"
#include "wzfitter/config.hpp"

//#define ALLOW_BLINDING // blinding on local is not possible, it needs cafe Event
#ifdef ALLOW_BLINDING
#include "wmass_blinding_util/OffsetMass.hpp"
#include "wmass_blinding_util/OffsetWidth.hpp"
#include "wmass_blinding_util/BlindingAuthority.hpp"
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <stdexcept>
#include <ostream>
#include <sstream>
#include <cmath>

using namespace std;
using TMath::Log2;

// forward declaration of minimization function
void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


class FitHandler {
    public :
        FitHandler()
#ifdef ALLOW_BLINDING
            : BA(2) 
#endif
        {
            cout << "+--------------------------------------------------------------------------------+" << endl;
            cout << "|                   Running the FitHandler                                       |" << endl;
            cout << "|                   For more info run with --help                                |" << endl;
            cout << "+--------------------------------------------------------------------------------+" << endl;
            _lasthadler = this;
            myFitter = 0;
            outf     = 0;
            outTree  = 0;
            toyHist  = 0;
        };
        ~FitHandler(){
        };

        static FitHandler * _lasthadler;

        void SetupOptions(){
            // load input histogram from tree?
            doHistFromTree = false;
            if (mtpt.Contains("_tree")){
                doHistFromTree = true;
                mtpt.ReplaceAll("_tree","");
            }

            // use pdf or pt shape reweighting
            pdf_index=-1;
            TRegexp number="_[0-9]+";
            if (mtpt.Contains(number)){
                // find the number
                TSubString sub(mtpt(number));
                // get value
                TString val(sub.Data());
                val.ReplaceAll("_","");
                pdf_index = val.Atoi();
                // remove number
                mtpt.Remove(sub.Start(),sub.Length());
            }

            doToyMC = false;
            nToys = 1000;
            ToyEvents = 0;
            if (compare_option > 10){
                doToyMC = true;
                compare_option = compare_option-10;
                // /prj_root/7059/wmass2/stark/RunIIcAnal/PMCS14/cabout 
                //ToyEvents = 59709871; // result_wenD_IIb3.root
                //ToyEvents+= 55574965; // result_wenD_IIb4.root
                ToyEvents = 1612010;
            }

            dorebin = false;
            if ( (mtpt=="Mt" || mtpt=="MT" || mtpt=="mt" ) && (compare_option==2) ){
                configfile = "runWMassMT.config";
                histname = "hWcandMt_CC";
                tempname = "hWcandMt_CC_";
                distname = "tempWMThists.root";
                dirname = "default";
                dorebin=false;
            } else if ( (mtpt=="Pt" || mtpt=="PT" || mtpt=="pt" ) && (compare_option==2) ){ 
                configfile = "runWMassPT.config";
                histname = "hWcandElecPt_CC";
                tempname = "hWcandElecPt_CC_";
                distname = "tempWPThists.root";
                dirname = "default";
            } else if ( (mtpt=="MET" || mtpt=="Met" || mtpt=="met" ) && (compare_option==2) ){ 
                configfile = "runWMassMET.config";
                histname = "hWcandMet_CC";
                tempname = "hWcandMet_CC_";
                distname = "tempWMEThists.root";
                dirname = "default";
            } else if ((mtpt=="Width" || mtpt=="WIDTH" || mtpt=="width" ) && (compare_option==2) ){ 
                configfile = "runWWidth.config";
                histname = "hWcandMt_CC_Width";
                tempname = "hWcandMt_CC_Width_";
                distname = "tempWidthhists.root";
                ToyEvents = 100;
                dirname = "default";
            } else if ((mtpt=="Ptup1" || mtpt=="PTUP1" || mtpt=="ptup1" ) && (compare_option!=2) ) { 
                configfile = "runWMassPTUP1.config";   // UP1 = U||<0 
                histname = "WCandElecPtUpara1_Spatial_Match_0"; //0==CC, spatial means _just_ spatial match
                tempname = "hWcandElecPtUpara1_CC_";
                distname = "tempWPThists.root";
                dirname = "WCand_Hist";
            } else if ((mtpt=="Ptup2" || mtpt=="PTUP2" || mtpt=="ptup2" ) && (compare_option!=2) ) { 
                configfile = "runWMassPTUP2.config";    // UP2 = U||>0
                histname = "WCandElecPtUpara2_Spatial_Match_0"; //0==CC, spatial means _just_ spatial match
                tempname = "hWcandElecPtUpara2_CC_";
                distname = "tempWPThists.root";
                dirname = "WCand_Hist";
            } else if ((mtpt=="Pt" || mtpt=="PT" || mtpt=="pt" ) && (compare_option!=2) ) { 
                configfile = "runWMassPT.config";
                histname = "WCandElecPt_Spatial_Match_0"; //0==CC, spatial means _just_ spatial match
                tempname = "hWcandElecPt_CC_";
                distname = "tempWPThists.root";
                dirname = "WCand_Hist";
            } else if ((mtpt=="Mt" || mtpt=="MT" || mtpt=="mt" ) && (compare_option!=2) ){ 
                configfile = "runWMassMT.config";
                histname = "WCandMt_Spatial_Match_0";
                tempname = "hWcandMt_CC_";
                distname = "tempWMThists.root";
                dirname = "WCand_Hist";
            } else if ((mtpt=="Mtup1" || mtpt=="MTUP1" || mtpt=="mtup1" ) && (compare_option!=2) ){ 
                configfile = "runWMassMTUP1.config";               // UP1 = U||<0 
                histname = "WCandMtUpara1_Spatial_Match_0";
                tempname = "hWcandMtUpara1_CC_";
                distname = "tempWMThists.root";
                dirname = "WCand_Hist";
            } else if ((mtpt=="Mtup2" || mtpt=="MTUP2" || mtpt=="mtup2" ) && (compare_option!=2) ){ 
                configfile = "runWMassMTUP2.config";               // UP1 = U||>0 
                histname = "WCandMtUpara2_Spatial_Match_0";
                tempname = "hWcandMtUpara2_CC_";
                distname = "tempWMThists.root";
                dirname = "WCand_Hist";
            } else if ((mtpt=="Metup1" || mtpt=="METUP1" || mtpt=="metup1" ) && (compare_option!=2) ){ 
                configfile = "runWMassMETUP1.config";               // UP1 = U||<0 
                histname = "WCandMetUpara1_Spatial_Match_0";
                tempname = "hWcandMetUpara1_CC_";
                distname = "tempWMEThists.root";
                dirname = "WCand_Hist";
            } else if ((mtpt=="Metup2" || mtpt=="METUP2" || mtpt=="metup2" ) && (compare_option!=2) ){ 
                configfile = "runWMassMETUP2.config";               // UP1 = U||>0 
                histname = "WCandMetUpara2_Spatial_Match_0";
                tempname = "hWcandMetUpara2_CC_";
                distname = "tempWMEThists.root";
                dirname = "WCand_Hist";  
            }else if ((mtpt=="Met" || mtpt=="MET" || mtpt=="met" ) && (compare_option!=2) ){
                configfile = "runWMassMET.config";
                histname = "WCandMet_Spatial_Match_0";
                tempname = "hWcandMet_CC_";
                distname = "tempWMEThists.root";
                dirname = "WCand_Hist";
            } else if ((mtpt=="Width" || mtpt=="WIDTH" || mtpt=="width" ) && (compare_option!=2) ){ 
                configfile = "runWWidth.config";
                histname = "WCandMt_Spatial_Match_Width_0";
                tempname = "hWcandMt_CC_Width_";
                distname = "tempWWidthhists.root";
                dirname = "WCand_Hist";
            }else{
                throw runtime_error(Form("Didn't understand the options, or requested option not implimented (mtpt=%s,compare_option=%d)",mtpt.Data(),compare_option));
            }

            if (pdf_index!=-1){
                histname+="_PDF_";
                histname+=pdf_index;
            }


            doBackground = false;
            if (compare_option==0) doBackground=true;

            _iswidth=false;
            if ((mtpt=="Width" || mtpt=="WIDTH" || mtpt=="width" )){
                _iswidth=true; 
            }

        }

        void ReadConfigFile(){
            // Read the config file...Obtain global parameters....
            configfile.Prepend("./wzfitter/config/");
            config runconfig(configfile.Data());
            _npar       = runconfig.getInt("numparam");               // number of parameters in the fit (e.g M(W) is a single parameter)
            _ndim       = runconfig.getInt("ndimensions");            // dimensionality of the histograms being fit (e.g. 1D histogram of MT variable)
            _numbins    = runconfig.getIntArray("numbins");           // number of bins in the histogram mentioned on the previous line, this number is set in wz_epmcs when templates are made
            _lowb       = runconfig.getDoubleArray("lowrange");       // low edge of the histogram referred to on the previous line, this number is set in wz_epmcs when templates are made
            _highb      = runconfig.getDoubleArray("highrange");      // high edge of the histogram referred to on the previous line, this number is set in wz_epmcs when templates are made
            frange      = runconfig.getDoubleArray("fit_range");      // fitting range that is used when the histogram referred on the previous line is being fit, fit range cannot extend beyond lowrange nor highrange values
            vstart      = runconfig.getDoubleArray("start_value");    // MINUIT starting values for each parameter
            vstep       = runconfig.getDoubleArray("vstep");          // MINUIT step size in each parameter (in each dimension of parameter space), typically we're talking about one parameter only, M_W or Gamma_W
            maxiter     = runconfig.getInt("numiter");                // maximum number of iterations before MINUIT gives an answer
            sensitivity = (double)runconfig.getFloat("sensitivity");  // sensitivity of MINUIT fit (how close MINUIT needs to get to a minimum before it stops) 
                                                                      // at the minimum. The default tolerance is 0.1, and  the minimization will stop when the estimated vertical distance to the minimum (EDM) is less than 0.001*[tolerance]*UP (see SET ERR).
            _limits     = runconfig.getBool("limits");                // constraining the allowed parameter range during fitting
            vmin        = runconfig.getDoubleArray("vmin");           // if it is a constrained fit (_limits==true), this is lower limit, again this is an array, i.e. a limit for each parameter can be specified here, if there is more than one parameter
            vmax        = runconfig.getDoubleArray("vmax");           // if it is a constrained fit (_limits==true), this is upper limit, again this is an array, i.e. a limit for each parameter can be specified here, if there is more than one parameter
            _likely     = runconfig.getBool("likelihood");            // this flag is used for defining error in MINUIT, specifically "arglist[0]                                                                                                                                                                                                                        = .5;" (this is NOT a switch between loglikelihood fit and chi2 fit)
            _blinded    = runconfig.getBool("blinded");               // flag for blinding the output of the fit
            /*--------------------------------------------------------- 
              the remaining variables from config file are not read here 
              (wzfitter reads them) but they are relevant:
              _intover=runconfig.getBool("intover"); // include (when _intover==true) both under- and over-flow bins when normalizing templates to the data
              _fitwithover=runconfig.getBool("fitwithover"); include (when _fitwithover==true) both
              _intrange=runconfig.getBool("intrange"); // when _intrange==true normalize histograms (templates to data) using the same range as for the fitting, how is the normalization range specified in case of _intrange==false (Z mass, W width) ?
              _fitantirange=runconfig.getBool("fitantirange"); // if one wants to *exclude* the fit_range while performing the fit over the full range this flag must be set to true 
              */
        }

        TH1D * GetHistFromTree(TTree * tree, TString tree_var, TString tree_weight, TString hist_name){
            tree->Draw(tree_var.Data(),tree_weight.Data(),"goff");
            TH1D * h = (TH1D*)tree->GetHistogram()->Clone(hist_name.Data());
            h->SetDirectory(0);
            return h;
            //return GetHistFromTree(tt, tree_var.Data(), tree_weight.Data(), hist_name.Data());
        }

        // TH1D * GetHistFromTree(TTree * tt, string tree_var, string tree_weight, string hist_name){
        //     tt->Draw(tree_var.c_str(),tree_weight.c_str(),"goff");
        //     return (TH1D*)tt->GetHistogram()->Clone(hist_name.c_str());
        // }

        TString GetTreeVar(){
            TString tree_var = "";
            // setup the variable
            if (mtpt.Contains("mt",TString::kIgnoreCase)){
                tree_var="VB_mt";
                tree_var+=">>(300,50,200)";
            } else if (mtpt.Contains("met",TString::kIgnoreCase)){
                tree_var="MET_et";
                tree_var+=">>(200,0,200)";
            } else if (mtpt.Contains("pt",TString::kIgnoreCase)){
                tree_var="el_pt";
                tree_var+=">>(200,0,100)";
            } else throw runtime_error(Form("Load hist from tree not implemented for variable %s",mtpt.Data()));
            return tree_var;
        }

        void LoadInputHistogram(){
            if (doHistFromTree){
                ff = TFile::Open(infilename.Data(), "READ");
                tt = (TTree *) ff->Get("physics");
                TString tree_var = GetTreeVar();
                TString tree_weight = "weight_evt";

                // pdf reweighting
                if (pdf_index!=-1) {
                    tree_weight+="*weights_PDF[";
                    tree_weight+=pdf_index;
                    tree_weight+="]";
                }

                // get histogram from tree
                inHist = GetHistFromTree(tt, tree_var, tree_weight, "hdatain");
                ff->Close();
            } else {
                ff = TFile::Open(infilename.Data(), "READ");
                cout << "file and hist" << infilename << " " << histname << endl;
                if (!ff) throw runtime_error(infilename.Prepend("Can not open input file: ").Data());
                ff->cd(dirname);
                inHist = (TH1D *) gROOT->FindObject(histname.Data()); // ->Clone("hdatain");
                inHist->SetDirectory(0);
                ff->Close();
            }
            if (inHist->GetEntries() < 100000) throw runtime_error(Form("The histogram %s have insuficient statistics (entries=%f)'",histname.Data(),inHist->GetEntries()));
            if(dorebin){ //for MT 
                inHist->Rebin(2);
            }
        }

        void GetTemplates(){
            // Set up templates
            int dimension = 1;
            int option = 2; // directory option: 0 in file root; 1 in smeared dir; 2 in default dir;
            string spline_type = "mass";
            TString map_name = "histd1map_WMassTemplates";
            TString backround_filename = "fitter_backgrounds.root";
            if (_iswidth){
                spline_type = "width";
                map_name = "histd1map_WWidthTemplates";
            }
            // create template file
            //if (doHistFromTree){
            if (false){
                tempfilename = "templates_fromTree_last.root";
                option = 0; // everything is in file root
                // create map histogram TODO: read values from file
                int    Wmass_nbins   = 99;
                double Wmass_default = 80.419;
                double Wmass_step    = 0.002;
                double Wmass_min     = Wmass_default - 50 * Wmass_step ;
                double Wmass_max     = Wmass_default + 49 * Wmass_step ;
                // TODO: check if the template file is older than tree file
                // create file
                TFile * ftmp = new TFile (tempfilename,"RECREATE" );
                // put norm histogram (`tempname` ending on "0")
                TString tmp_histname = tempname ; tmp_histname+=0;
                TString tmp_weightname = "weight_evt" ;
                GetHistFromTree(tt, GetTreeVar(), tmp_weightname, tmp_histname )->Write();
                // put there histogram template histogram ( `tempname` + sample number)
                for (int i_w = 0; i_w < Wmass_nbins+1; i_w++){
                    tmp_histname = tempname; tmp_histname+=(i_w+1);
                    tmp_weightname = "weights_PDF*weights_mass["; 
                    tmp_weightname += (i_w);
                    tmp_weightname += "]";
                    GetHistFromTree(tt, GetTreeVar(), tmp_weightname, tmp_weightname)->Write();
                }
                // put there map histogream (`map_name`)
                TH1D * mapfile = new TH1D(map_name.Data(), "map template -- mass value", Wmass_nbins, Wmass_min, Wmass_max);
                mapfile -> Write();
                // save file
                ftmp->Write();
                ftmp->Close();
            }
            // load splines
            cout<<"creating spline function for "<< spline_type << ".."<<endl;
            templates = new gen1Dspline(const_cast<char*>(tempfilename.Data())
                    , const_cast<char*>(tempname.Data())
                    , const_cast<char*>(map_name.Data())
                    , dimension
                    , option
                    , backround_filename
                    , doBackground);
        }

        void Setup(int argc, const char * argv[]){

            // Parse command line
            if(argc!=6){
                help();
                throw runtime_error(Form("Bad number of input arguments: %d",argc));
            }

            outfilename    = argv[1] ;
            infilename     = argv[2] ;
            tempfilename   = argv[3] ;
            mtpt           = argv[4] ;
            compare_option = TString(argv[5]).Atoi() ;

            cout << "Input file: "     << infilename     << endl;
            cout << "Output file: "    << outfilename    << endl;
            cout << "Template file: "  << tempfilename   << endl;
            cout << "Compare option: " << compare_option << " using Method: " << mtpt << endl;

            SetupOptions();
            ReadConfigFile();
            LoadInputHistogram();
            GetTemplates();


            result = new TGraph(1);    result -> SetName("result"); result->SetMarkerStyle(20);
            scan   = new TGraph(200);  scan   -> SetName("scan");

            if (doToyMC) {
                TH1::AddDirectory(kFALSE);
                if (!outf) outf = TFile::Open(outfilename.Data(),"RECREATE");
                outTree = new TTree("ToyMC", "ToyMC");
                outTree->Branch("fit_val" , &paramval , "fit_val/D");
                outTree->Branch("fit_err" , &errval   , "fit_err/D");
                outTree->Branch("fit_est" , &estval   , "fit_est/D");
            }



            // Set fit ranges
            vector<double> d1range;
            d1range.push_back(frange[0]);  //32-48 or 65-95
            d1range.push_back(frange[1]);
            myRange.push_back(d1range);


            // Set up fitter
            //cout << " datahist1 : " << inHist->GetNbinsX() << " " << inHist->GetMean() << " "  << inHist->GetRMS() << " " << inHist->Integral() << endl;
            myFitter = new wzfitter(templates, inHist, const_cast<char*>(configfile.Data()), _likely );
            //cout << " datahist2 : " << inHist->GetNbinsX() << " " << inHist->GetMean() << " "  << inHist->GetRMS() << " " << inHist->Integral() << endl;
        }

        void SetupNewToy(){
            if (!myFitter) delete myFitter;
            if (!toyHist) delete toyHist;
            toyHist = (TH1D*) inHist->Clone("toyHist");
            toyHist->Reset();
            for (int i=0; i<ToyEvents; i++) toyHist->Fill(inHist->GetRandom());
            myFitter = new wzfitter(templates, toyHist, const_cast<char*>(configfile.Data()), _likely );
        }

        void Thefcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
#ifdef ALLOW_BLINDING
            if (_blinded) {
                BlindingAuthority BA(2);
                if(!_iswidth) par[0] = BA.Get_Offset_Mass()->CalculateTrueMassFromBlindedMass(par[0]);
                else par[0] = BA.Get_Offset_Width()->CalculateTrueWidthFromBlindedWidth(par[0]);

            }
#endif

            bool isBinnedLikelyHood = true;
            if (_likely) f = myFitter->doLogLikelihood(par,isBinnedLikelyHood, myRange);
            else         f = myFitter->doChiSquared   (par,                    myRange);

            //cout << "thefcn: " << par[0] << " : " << f << "  ("<< myRange[0][0] <<","<< myRange[0][1] <<")" <<  endl;
        }


        ////////////////////////////////////////////////////////////////////////////////////////////////
        // This function does the fit...regardless of what we're fitting 
        void doFit() {
            Double_t amin,edm,errdef;      // MINUIT variables...explained later
            Int_t nvpar,nparx,icstat;      // " " "

            // Initialize TMinuit
            TMinuit *gMinuit = new TMinuit(_npar);
            gMinuit->SetFCN(FCN);
            Double_t arglist[10];
            Int_t ierflg = 0;
            // mnexcm : minuit execute command ( command, parameter list, length of list, error flag)
            //  list of commands: http://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node18.html
            // SET PRINT : set verbosity
            arglist[0] = -1.; // -1: no output, ..., 3: maximum output
            gMinuit -> mnexcm("SET PRINT", arglist, 1, ierflg);
            // SET ERR initializes the error definition...
            // Sets the value of UP (default value= 1.), defining parameter errors. 
            // Minuit defines parameter errors as the change in parameter value 
            // required to change the function value by UP. Normally, for chisquared 
            // fits UP=1, and for negative log likelihood, UP=0.5.
            if(_likely) arglist[0] = .5;
            else arglist[0] = 1.;
            gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

            // Set starting values and step sizes for parameters
            // tells minuit the starting value and step size for varying the 
            // paramter(s) of the fit  
            // mnparm: implement one paramter ( param number, param name, stat value, step, low limit, upp limit, err flag)
            for(int i=0; i<_npar; i++)
            {
                if(_limits) gMinuit->mnparm(i, "a1", vstart[i], vstep[i], vmin[i], vmax[i],ierflg); // are there upper and lower fit limits?
                else        gMinuit->mnparm(i, "a1", vstart[i], vstep[i],       0,       0,ierflg); // or no limits?
                //      cout<<vstart[i]<<" "<<vstep[i]<<endl;
            }


            // Now ready for minimization step
            arglist[0] = maxiter;  // the maximum number of iterations
            arglist[1] = sensitivity;   // specifies required tolerance on the function value 
            // at the minimum. The default tolerance is 0.1, and 
            // the minimization will stop when the estimated vertical
            // distance to the minimum (EDM) is less than 
            // 0.001*[tolerance]*UP (see SET ERR).
            gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
            //Notice there are two parameters
            // MIGRAD is the command to commence minimization

            // Print results
            //                mnstat returns amin: the minimum
            //                edm: the estimated vertical distance remaining to minimum
            //                errdef: the value of UP defining parameter uncertainties
            //                nvpar: the number of currently variable parameters
            //                nparx: the highest (ext) parameter number defined by user
            //                icstat: status integer indicating  goodness of covariance
            //gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
            // Prints the values of the parameters at the time of the call
            //gMinuit->mnprin(_npar,amin);
            //Prints the covariance matrix
            //  gMinuit->mnmatu(1);

            // Get fit results: parameters and their errors:
            _gpar = new Double_t[_npar];
            _gerr = new Double_t[_npar];
            for(int i=0; i<_npar; i++)
            {
                gMinuit -> GetParameter(i, _gpar[i], _gerr[i]);
            }
            delete gMinuit; 
        }

        void ScanDistributions(){
            double *apar = new Double_t[_npar];
            double fixedLhood = myFitter->doLogLikelihood(_gpar,true,myRange);


            double mw = 80.300;
            if (_iswidth){
                mw = 1.971;
                if (_blinded) mw=paramval-0.1;
                for (int im = 0; im<201; im++){
                    apar[0]=mw + im/1000.;
                    //cout<<apar[0]<<" "<<(mw + im/1000.)<<endl;
                    double appar = apar[0];
#ifdef ALLOW_BLINDING
                    if (_blinded) apar[0] = BA.Get_Offset_Width()->CalculateTrueWidthFromBlindedWidth(apar[0]);
#endif
                    //double lval  = (myFitter->doLogLikelihood(apar,true,myRange))-fixedLhood;
                    double lval  = (myFitter->doChiSquared   (apar,     myRange));
                    scan->SetPoint(im,appar,lval); 
                }
            }
            else{
                mw = 80.300;
                mw = paramval - errval*3;
                if (_blinded) mw = paramval-0.1;
                for (int im = 0; im<201; im++){
                    //apar[0]=mw + im/1000.;
                    apar[0]=mw + im*errval*0.03;
                    //cout<<apar[0]<<" "<<(mw + im/1000.)<<endl;
                    double appar = apar[0];
#ifdef ALLOW_BLINDING
                    if (_blinded) apar[0] = BA.Get_Offset_Mass()->CalculateTrueMassFromBlindedMass(apar[0]);
#endif
                    double lval  = myFitter->doLogLikelihood(apar,true,myRange);
                    //double lval  = myFitter->doChiSquared   (apar,     myRange);
                    scan->SetPoint(im,appar,lval-fixedLhood); 
                }
            }
        }

        void FindBest(){
            TNHist* comphist;
            //ff = new TFile(infilename.Data(),"READ");
            //ff->cd(dirname);
            //inHist = (TH1D*)gROOT->FindObject(histname);

            if (ff){
                vector<int> fnumbins=_numbins;
                vector<double> flowb=_lowb;
                vector<double> fhighb=_highb;
                double *ppar;
                ppar=&paramval;

#ifdef ALLOW_BLINDING
                if (_blinded){
                    if(!_iswidth) ppar[0] = BA.Get_Offset_Mass()->CalculateTrueMassFromBlindedMass(ppar[0]);
                    else ppar[0] = BA.Get_Offset_Width()->CalculateTrueWidthFromBlindedWidth(ppar[0]);
                }
#endif

                comphist = templates->MakeDist(ppar,fnumbins,flowb,fhighb);
                best = new TH1D("best","bestfit",fnumbins[0],flowb[0],fhighb[0]);
                for(Int_t j=0;j<(comphist->GetSize());j++) {
                    best->SetBinContent(j,comphist->GetBinContent(j));
                }
                double rhi = best->Integral( best->FindBin(frange[0]),   best->FindBin(frange[1]));
                double dhi = inHist->Integral( best->FindBin(frange[0]),   best->FindBin(frange[1]));
                best->Scale(dhi/rhi);
            }

        }

        void CheckBlinding(){
#ifdef ALLOW_BLINDING
            if (BA.isConfigured() ) return;
            TFile * fpmcs = TFile::Open(tempfilename,"READ");
            ff = TFile::Open(infilename,"READ");
            //cout<<" the directory ?"<<fpmcs->FindObject("default")<<endl;
            BA.SetConfiguration(ff,fpmcs,compare_option);
            if ((compare_option == 0) && (_blinded == false)) {
                throw runtime_error ( "Blinding must be turned on to study data." );
            }
            else if ((_blinded==true) && (compare_option!=0) ){
                throw runtime_error ( "Blinding must be turned off to study fast/full MC." );
            }
#endif
        }

        void PerformFit(){

            //cout << "dumpstart" <<endl;
            //cout << histname     .Data() << endl;
            //cout << tempname     .Data() << endl;
            //cout << distname     .Data() << endl;
            //cout << configfile   .Data() << endl;
            //cout << dirname      .Data() << endl;
            //cout << outfilename  .Data() << endl;
            //cout << infilename   .Data() << endl;
            //cout << tempfilename .Data() << endl;
            //cout << mtpt         .Data() << endl;
            //cout << compare_option << endl;
            //cout << pdf_index      << endl;
            //cout << "dumpend" <<endl;

            int iToy = 0;
            for (; iToy <= nToys; iToy++){

                /////////////////////////////////////////
                // Before doFit always check with the BA
                //
                CheckBlinding();
                //
                /////////////////////////////////////////

                doFit();
                paramval = _gpar[0];
                errval   = _gerr[0];
                estval   = myFitter->doLogLikelihood(_gpar,true,myRange);
                if (iToy==0){
                    cout << "the value of the parameters are... "  << paramval << endl;
                    cout << "the errors of the parameters are... " << errval   << endl;
                    result->SetPoint(0,paramval,errval);
                    bool doScan=true;
                    if (doScan) ScanDistributions();
                    FindBest();
                }
                if (!doToyMC) break;
                if (iToy>0) 
                    outTree->Fill(); // first run is for full fit
                else {
                    //outTree->Fill(); // first run is for full fit
                    cout << "=Running ToyMC=" << endl;
                    cout << " N Toys: " << nToys 
                         << " events per toy: " << ToyEvents
                         << " Ntoys x events: " << ToyEvents*nToys
                         << " input hist entries: " << inHist->GetEntries()
                         << " fullMC/allToy ratio: " << inHist->GetEntries()/(ToyEvents*nToys)
                         << endl;
                }
                if( ! fmod( Log2(iToy), 1 ) ) cout<<"toys done: "<<iToy<<endl;
                SetupNewToy();
            }

        }

        void Save(){
            if (!outf) outf = TFile::Open(outfilename.Data(),"RECREATE");
            outf->cd();

            inHist->Write();
            best->Write();
            if(scan!=0) scan->Write();
            result->Write();

            outf->Write();
            outf->Close();
        }

        void help(){
            std::cout<< "Usage : ./findWMass_fromhist_x <output-filename> <input-filename> <template-filename> <Mt-or-Pt-or-MET-or-Width> <compare-option>" << std::endl;
            std::cout<<" Use <compare_option>=0  for comparing DATA with Fast-MC (PMCS),"<<endl
                <<"     <compare_option>=1  for comparing DATA with DATA,"<<endl
                <<"     <compare_option>=2  for comparing Fast-MC (PMCS) with Fast-MC (PMCS),"<<endl 
                <<"     <compare_option>=3  for comparing Full-MC (GEANT) with Fast-MC (PMCS)," << endl
                <<"     <compare_option>=4  for comparing Full-MC (GEANT) with Full-MC (GEANT)." << endl
                <<"     <compare_option>=5  for comparing DATA with Full-MC (GEANT)." << endl
                <<"     <compare_option>=6  for comparing fast-mc-reweighted with fast-mc." << endl
                <<"     <compare_option>+10  for runnig toy Monte-Carlos." << endl
                <<"     Only 0, 2 and 3 are useful." << endl;
        }


    private :
        TH1D              *inHist;      ///< Input histogram. The mass is measured from.
        TH1D              *toyHist;     ///< ToyMC histogram. The mass is measured from.
        TH1D              *best;        ///< Best template. The mass is measured from.
        TGraph            *scan;        ///< Chi2 or NLL scan. Value of chi2 or NLL per each template.
        TGraph            *result;      ///< W mass fit. x=fitted value, y=uncertainty.
        TTree             *outTree;     ///< Output tree for ToyHistograms.

#ifdef ALLOW_BLINDING
        BlindingAuthority  BA;          ///< Tool for blinding samples.
#endif

        wzfitter          *myFitter;    ///< Estimator for hist--template difference.
        gen1Dspline       *templates;   ///< Template holder. Handled as 1D splines.

        double            *_gpar;       ///< Pointer to parameter value array for minuit.
        double            *_gerr;       ///< Pointer to parameter error array for minuit.
        double             paramval;    ///< Final value.
        double             errval;      ///< Final error.
        double             estval;      ///< Final value of estimator.

        TFile             *ff;          ///< Pointer to currently opened file.
        TFile             *outf;        ///< Pointer to output file.
        TTree             *tt;          ///< Pointer to currently opened tree.

        vector<vector<double> > myRange; ///< Vector of ranges. Contain 2 item vectors: [ [min1,max1], [min2,max2], ...]

        // options
        TString        outfilename;
        TString        infilename;
        TString        tempfilename;
        TString        mtpt;
        int            compare_option;

        int            pdf_index;        ///< index of reweighted distribution

        // fitting options
        vector<double> frange;           ///< fitting range that is used when the histogram referred on the previous line is being fit, fit range cannot extend beyond lowrange nor highrange values
        vector<double> vstart;           ///< Start value of parameters.
        vector<double> vstep;            ///< Step value of parameters.
        vector<double> vmin;             ///< Limit minimum of parameters.
        vector<double> vmax;             ///< Limit maximum of parameters.
        Int_t          maxiter;          ///< Maximal number of iterations.
        Double_t       sensitivity;      ///< Required tolerance on the function value.
        bool           _limits;          ///< Apply limit to parameter?
        bool           _likely;          ///< Use NLL or chi2.

        Int_t          _npar;            ///< number of parameters in the fit (e.g M(W) is a single parameter)
        Int_t          _ndim;            ///< dimensionality of the histograms being fit (e.g. 1D histogram of MT variable)
        vector<int>    _numbins;         ///< number of bins in the histogram mentioned on the previous line, this number is set in wz_epmcs when templates are made
        vector<double> _lowb;            ///< low edge of the histogram referred to on the previous line, this number is set in wz_epmcs when templates are made
        vector<double> _highb;           ///< high edge of the histogram referred to on the previous line, this number is set in wz_epmcs when templates are made

        Int_t          nToys;            ///< number of ToyMC experiments.
        Int_t          ToyEvents;        ///< number of events per one ToyMC.

        TString        histname;         ///< Name of input histogram
        TString        tempname;         ///< Base name of template histograms
        TString        distname;         ///< External template histogram file. Not used!
        TString        configfile;       ///< Name of text file with configuration.
        TString        dirname;          ///< Name of directory insinde input root file.


        bool           doToyMC;          ///< Perform fitting of Monte-Carlo toys.
        bool           doHistFromTree;   ///< Load input hisogram from TTree.
        bool           doBackground;     ///< Switch to use background.
        bool           _iswidth;         ///< Switch whether we are fitting widht or mass.
        bool           dorebin;          ///< Switch to rebin input distribution.
        bool           _blinded;         ///< Switch to use blinding.

};


// Function for minization have to be global static  ...   WTF!??$^*#%)@#$!!!
// introduce static pointer to last created handler (similar to singleton, but with more than one instances allowed)
FitHandler * FitHandler::_lasthadler = 0;
// FCN ( number of parameters, gradient of function, function output, parameters, switch-case flag )
void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    FitHandler::_lasthadler->Thefcn(npar, gin, f, par, iflag);
}

/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){
    try{
        FitHandler fitter;
        fitter.Setup(argc,argv);
        fitter.PerformFit();
        fitter.Save();
    } catch (runtime_error e) {
        cout << Form("Terminating the fitting because of excetion:\n\t%s",e.what()) << endl;
        return 1;
    }

    return 0;
}


#endif // wfitter_jcuth
