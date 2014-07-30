#ifndef test
#define test

/**
 * @file test.cxx
 * Description of this file
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2014-07-30
 */


#include "TString.h"
#include "Output.hpp"
#include <iostream>

using std::cout;
using std::endl;


/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){
    TString h="Bye...";
    Output out;
    cout << "Hello world! " << h << endl;

    return 0;
}


#endif // test
