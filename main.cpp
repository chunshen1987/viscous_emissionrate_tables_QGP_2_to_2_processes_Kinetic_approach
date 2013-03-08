#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Arsenal.h"
#include "Matrix_elements_sq.h"
#include "Phasespace_integrals.h"
#include "parameters.h"
#include "Stopwatch.h"

using namespace std;


int main()
{
   Stopwatch sw; 
   sw.tic();

   int channel = 0;
   string filename;
   
   //Compton scattering
   filename = "QGP_2to2_Compton";
   channel = 1;
   Calculate_emissionrates(channel, filename);
   
   //Pair Annihilation
   filename = "QGP_2to2_PairAnnihilation";
   channel = 2;
   Calculate_emissionrates(channel, filename);

/**********************************************************************/
// passed test
/**********************************************************************/

   sw.toc();
   cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
