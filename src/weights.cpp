#include "weights.h"
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

const float Weights::A =-2.01;
const float Weights::B =1.63; 
const float Weights::A2=-2.00;
const float Weights::B2=2.1535;
const float Weights::M2= 0.0783;
const float Weights::C2= 0.104;

long double* Weights::read_weight(const char* fileName, float& bias, int size)
   {
   long double* w= new long double[size];
   
   ifstream infile(fileName);
   if(!infile)
     {
     cerr << "\nCould not find weight file : \n" << fileName << endl;
     exit(1);
     }
   
   infile >> bias;
   for(int i=0; i<size; i++)
     infile >> w[i];
     
   infile.close();
   return w;
   }
   
long double* Weights::read_weight(const char* fileName, int size)
   {
   long double* w= new long double[size];
   
   ifstream infile(fileName);
   if(!infile)
     {
     cerr << "\nCould not find weight file : \n" << fileName << endl;
     exit(1);
     }
   
   for(int i=0; i<size; i++)
     infile >> w[i];
     
   infile.close();
   return w;
   }
 
   
   
      


   
     
