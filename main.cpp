#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Arsenal.h"
#include "Stopwatch.h"
#include "QGP_2to2_Scattering_Kinetic.h"
#include "ParameterReader.h"

using namespace std;


int main(int argc, char** argv)
{
   Stopwatch sw; 
   sw.tic();

   int channel = 0;
   string filename;
   
   ParameterReader* paraRdr = new ParameterReader();
   paraRdr->readFromFile("parameters.dat");
   paraRdr->readFromArguments(argc, argv);
   
   QGP_2to2_Scattering_Kinetic test(paraRdr);

   //Compton scattering
   filename = "QGP_2to2_Compton";
   channel = 1;
   test.calculateEmissionrates(channel, filename);
   //Calculate_emissionrates(channel, filename);
   
   //Pair Annihilation
   filename = "QGP_2to2_PairAnnihilation";
   channel = 2;
   test.calculateEmissionrates(channel, filename);
   //Calculate_emissionrates(channel, filename);

/**********************************************************************/
// passed test
/**********************************************************************/

   sw.toc();
   cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
