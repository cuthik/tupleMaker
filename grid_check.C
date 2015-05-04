#ifndef GRID_CHECK_C
#define GRID_CHECK_C

/**
 * @file grid_check.C
 * Description of this file
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-04-14
 */

#include "TString.h"

#include <map>
#include <vector>

typedef vector<double> BinVec;
typedef vector<TString> StringVec;
typedef map<TString,BinVec> BinsMap;
typedef map<TString,double> ValMap;

class GridHandler{
    public:
        ~GridHandler(){};

        void ReadBinRange(){

        }

        void PrintBins(TString axis){
            BinVec *bb = & bins[axis.Data()];
            printf("binning of %s: [");
            int i = 0;
            for (BinVec::iterator it = bb->begin(); it != bb->end(); ++it ){
                if (!(i%10)) printf("\n    ");
                printf("%5f ,", (*it));
            }
            printf("\n  ]\n");
        }

        // data members
        StringVec names_var;
        StringVec names_value_sing;
        StringVec names_value_pert;
        BinsMap bins;
        ValMap values;


};

void grid_check(){
    GridHandler grid ("/etapfs03/atlashpc/cuth/resbos/grids/lhc7/wenu/ctnnlo-eigsets_W_kc1/y_w-_lhc7_ctnn00_gnw_kc1.out");
    grid.ReadBinRange();
    for (StringVec::iterator it=grid.names_var.begin(); it!=grid.names_var.end(); it++){
        grid.PrintBins((*it).Data())
    }

    return;
}



/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){

    return 0;
}


#endif // GRID_CHECK_C
