// DISOPRED V2.32
// By J. Ward and D. Jones
// Copyright (C) 2006 University College London
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <math.h>

#include "weights.h"

using namespace std;

enum amacid{ala, arg, asn, asp, cys,
            gln, glu, gly, his, ile,
	    leu, lys, met, phe, pro,
	    ser, thr, trp, tyr, val};


class protein 
{
private:
    static const int UCAP = 20, CAP= 26, WLEN= 15, spos=7;
    float* darr[UCAP], unarr[CAP];
    int pos;
    int plen, num, *conf;
    char* seq, *lab;
    long double* pred, *pred2;
    float thresh;
    Weights w;   /*ok for single prediction, to run multiple predictions from a single
		   executable this should be initialized and passed to protein object as an 
		   argument from main()*/
public:
    protein(const char*, const char*, int);
    void Npreds();             //predict N-termini
    void Mpreds(); 
    void Cpreds();             //predict C-termini
    void predict(string path); //obtain final prediction
};

protein::protein(const char* file, const char* path, int n) : pos(0), num(n), w(path) {
    float temp;
    int count, i, j;
    
    if(num < 1 || num > 10)
    {
	num = 5;
	cerr << "\nFalse positive rate: " << num << " not allowed\n"
	     << "Rate must be an integer in the range 1-10\n"
	     << "Using default: " << num << endl; 
    } 
    
    thresh= w.THolds[num-1];
    
    ifstream infile(file);
    if(!infile)  
    {
	cerr << "\nCould not open matrix file! : " << file << endl;
	exit(1);
    }
    
    infile >> plen;
    
    string tmpstr;
    conf = new int[plen+1];
    seq  = new char[plen+1];
    lab  = new char[plen+1];
    pred = new long double[plen+1];
    pred2= new long double[plen+1];
    
    
    
    for(i=0; i<UCAP; i++)
	darr[i]= new float[plen+1];
    
    for(i=0; i<plen; i++)
	infile >> seq[i];

    for(i=0; i<12; i++)
	getline(infile, tmpstr);
    
    for(i=0; i<plen; i++) {
	count =0;
	getline(infile,tmpstr);
	if (infile.bad()) {
	  cerr << "\nError while reading mtx file!";
	  exit(1);
	}

	stringstream s(tmpstr.c_str());
	for(j=0; j<26; j++) {
	    s >> temp;
	    if(!(j==0 || j==2 || j==21 || j>22))
		unarr[count++]=temp;
        }

	darr[ala][i]= unarr[0];
	darr[arg][i]= unarr[14];
	darr[asn][i]= unarr[11];
	darr[asp][i]= unarr[2];
	darr[cys][i]= unarr[1];
	darr[gln][i]= unarr[13];
	darr[glu][i]= unarr[3];
	darr[gly][i]= unarr[5];
	darr[his][i]= unarr[6];
	darr[ile][i]= unarr[7];
	darr[leu][i]= unarr[9];
	darr[lys][i]= unarr[8];
	darr[met][i]= unarr[10];
	darr[phe][i]= unarr[4];
	darr[pro][i]= unarr[12];
	darr[ser][i]= unarr[15];
	darr[thr][i]= unarr[16];
	darr[trp][i]= unarr[18];
	darr[tyr][i]= unarr[19];
	darr[val][i]= unarr[17];
	
	for(j=0; j<20; j++)
	    darr[j][i]=w.g[j]*darr[j][i]+w.c[j];
    }
    Npreds();
    Mpreds();   
    Cpreds();
}

void protein::Npreds() 
{
    for (int i=0; i<spos; i++)   
    {
	float fvec[315]={0};
	int feat=0;
	for(int j=(i-spos); j <(i+spos+1); j++)   
	{
	    if ( j<0) {
		feat+=UCAP;
		fvec[feat++]=1;
	    }
	    else   {
		for (int k=0; k < UCAP; k++)
		    fvec[feat++]=darr[k][j];
		fvec[feat++]=0;
	    }
	}
        pred[pos] = w.Nbias;
        
        for(int j=0; j<315; j++)
	    pred[pos]+=(fvec[j]*w.Nweights[j]);
	
        ++pos;
    }
}


void protein::Mpreds() 
{
    for (int i=spos; i<(plen-spos); i++)   {
	float fvec[300]={0};
	int feat=0;
	for(int j=(i-spos); j <(i+spos+1); j++)   
	{
	    for(int k=0; k <UCAP; k++)
		fvec[feat++]=darr[k][j];
	}
        
        pred[pos] = w.Mbias;
        
        for(int j=0; j<300; j++)
	    pred[pos]+=fvec[j]*w.Mweights[j];
        
        ++pos;
    }
}


void protein::Cpreds()   
{
    for (int i=(plen-spos); i<plen; i++)   {
	float fvec[315]={0};
	int feat=0;
        
	for(int j=(i-spos); j <(i+spos+1); j++)   {
	    if ( j<plen)   {
		for (int k=0; k < UCAP; k++)
	            fvec[feat++]=darr[k][j];
		fvec[feat++]=0;
	    }
	    else  {
                feat+=UCAP;
		fvec[feat++]=1;
	    }
	}
        
        pred[pos] = w.Cbias;
        
        for(int j=0; j<315; j++)
	    pred[pos]+=(fvec[j]*w.Cweights[j]);
        
        ++pos;
    }
}


void protein::predict(string path)  {
    int wid= 60;
    
    for(int i=0; i<plen; i++)
	pred[i]=1/(1+exp(w.A*pred[i]+w.B));
    
    ofstream outfileA;
    outfileA.setf(ios::fixed);
    
    string temp= path + ".horiz_d";
    string text= path + ".diso";
    
    for (int i=0; i<spos; i++)
	pred2[i]=pred[i];
    
    for (int i=spos; i<(plen-spos); i++)   
    {
        float fvec[15]={0};
        int feat=0;
        for(int k=(i-spos); k <(i+spos+1); k++)
	    fvec[feat++]=pred[k];
	
        pred2[i] = w.filtbias;
        
        for(int j=0; j<15; j++)
	    pred2[i]+=fvec[j]*w.filtweights[j];
        
        if(pred2[i]>0)
	    pred2[i]=1/(1+exp(w.A*pred2[i] + w.B2));
	else
	    pred2[i]=w.M2*pred2[i] + w.C2;
    }
    
    for(int i=(plen-spos); i<plen; i++)
	pred2[i]=pred[i];
    
    for(int i=0; i<plen; i++)  {     
	int j=0, c=0;
	while( j<10 )
	{
	    if(pred2[i] > w.THolds[j])
	    {c=9-j; break;}
	    j++;
	}
	conf[i]=c;
	
	if(pred2[i] >= thresh)
	    lab[i]= '*';
	else
	    lab[i]= '.';
    }
    
    outfileA.open(text.c_str());
    outfileA << "       -----DISOPRED version 2-----\n"
	     << "Disordered residues are marked with asterisks (*)\n"
	     << "   Ordered residues are marked with dots (.)\n"
	     << "Predictions at a false positive rate threshold of: " << num << "%\n\n"; 
    
    for(int i=0; i<plen; i++)
	outfileA << setw(5) << i+1 << ' ' << seq[i] << ' ' << lab[i] << ' ' 
		 << setw(8) << setprecision(3) << pred[i] << setw(8) << pred2[i] << endl;
    
    outfileA.close();
    
    
    outfileA.open(temp.c_str());
    outfileA << "DISOPRED predictions for a false positive rate threshold of: " << num << "%\n" << endl;
    
    outfileA << "conf: ";
    
    for(int i=0; i<plen; i++) {
	outfileA << conf[i];
	
	if( (i+1)%wid==0)  {
	    outfileA << endl << "pred: ";
	    for(int j=i+1-wid; j<i+1; j++)
		outfileA << lab[j];
	    
	    outfileA << endl << "  AA: ";
	    for(int j=i+1-wid; j<i+1; j++)
		outfileA << seq[j];
	    
	    outfileA << endl << "      ";
	    for(int j=i+1-wid; j<i+1; j++) {
		if(  (j+1)%10 ==0)
		    outfileA << setw(10) << j+1;
	    }
	    outfileA << endl << endl << "conf: ";
	}
	
	if(i==plen-1)  
	{
	    outfileA << endl << "  AA: ";
	    for(int j=wid*(i/wid); j<plen; j++)
		outfileA << seq[j];
	    
	    outfileA << endl<< "pred: ";
	    for(int j=wid*(i/wid); j<plen; j++)
		outfileA << lab[j];
	    
	    outfileA << endl<< "      ";
	    for(int j=wid*(i/wid); j<plen; j++) {
		if( (j+1)%10 ==0)
		    outfileA << setw(10) << j+1;
	    }
	    outfileA << endl << endl;
	}
    }
    
    outfileA << "Asterisks (*) represent disorder predictions and dots (.) \n"
	     << "prediction of order. The confidence estimates give a rough\n"
	     << "indication of the probability that each residue is disordered.\n" << endl;
    outfileA.close();
}   


int main(int argc, char* argv[]) {
    int rate;

    if(argc<4) {
	printf("ERROR: Need to include input argument and data file\n");
	return 0;
    }
    
    if(argc==4)
	rate=5;
    else
	rate= atoi(argv[4]);
    
    protein p1(argv[2], argv[3], rate);
    p1.predict(argv[1]);

    return 0;
}   
