#include <string>

class Weights
   {
   private:
     std::string path;
   public: 
     long double* linScale;
     long double* g, *c;                              //linearly rescale PSI-BLAST profiles
     long double* Nweights, *Mweights, *Cweights;     //Weights for C and N-termini and middle of protein
     long double* filtweights;                        //Weights for smoothing classifier
     long double* THolds;
     static const float A, B, A2, B2, M2, C2;
     float Nbias, Mbias, Cbias, filtbias;
       
     Weights(std::string p) : path(p)
        {
        using std::string; 
        string N= path + "Nweights.dat";
        string M= path + "Mweights.dat";
        string C= path + "Cweights.dat";
        string F= path + "2weights.dat";
        string L= path + "linScale.dat";
        
        Nweights= read_weight(N.c_str(), Nbias, 315);
        Mweights= read_weight(M.c_str(), Mbias, 300);
        Cweights= read_weight(C.c_str(), Cbias, 315);
        filtweights= read_weight(F.c_str(), filtbias, 15);
        linScale= read_weight(L.c_str(), 50);
        
        g= linScale;
        c= linScale+20;
        THolds= linScale+40;
        }
     long double* read_weight(const char*, float&, int);
     long double* read_weight(const char*, int); 
   };



