#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <time.h>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib> // for exit function
#include <tr1/random>
#include <numeric>
			
/*-----------------------------------------------------Header definitions---------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/

#ifndef RANDOMC_H
#define RANDOMC_H

const int nx=28;//grid cells for rain input data
const int ny=28;//grid cells for rain input data
const int TL=10;//juvenile development time (egg to adult)

using namespace std;
//using std::ifstream;
using namespace std::tr1;

const double PI=3.14159265;
const double TWOPI=6.28318531;
const long long int LONG_MAX=3147483647000;


/*----------------------struct for keeping track of global numbers------------------------*/
struct totals
	{
	
	long long int Jww,Jwd,Jdd,Jwr,Jrr,Jdr;
	long long int Mww,Mwd,Mdd,Mwr,Mrr,Mdr;
	long long int Vww,Vwd,Vdd,Vwr,Vrr,Vdr;
	long long int Fww,Fwd,Fdd,Fwr,Frr,Fdr;
	long long int JTot,MTot,VTot,FTot;
	double distW,distD,distR;
	};	
/*----------------------------------------------------------------------------------------*/


/*----------------------struct combining initial condition parameters---------------------*/
struct initials
{
	int NumAdultsWM; int NumAdultsWV; int NumAdultsWF;
	int NumAdultsDM; int NumAdultsDV; int NumAdultsDF;
	int NumAdultsRM; int NumAdultsRV; int NumAdultsRF;
	int NumJW[TL]; int NumJD[TL]; int NumJR[TL];
	string inputfile;
	string rainfile;
	double driver_time;
	double r_time;
	int NumDriver;
	int NumDriverPat;
	int NumR;
	int NumRes;
	char dist;
	int recPatFreq;
	};	
/*----------------------------------------------------------------------------------------*/

/*---------------------------number of settlements in each grid cell----------------------*/
vector<int>SetsPerCell[nx][ny]; 
/*----------------------------------------------------------------------------------------*/
		
/*------------------------state of population in each settlement---------------------------*/
struct Patch
{
	double x;
	double y;
	int sqx,sqy;
	long long int Jww[TL];
	long long int fJww[TL];
	long long int mJww[TL];
	long long int bJww[TL];
	long long int fJwd[TL];
	long long int mJwd[TL];
	long long int bJwd[TL];
	long long int Jdd[TL];
	long long int Jwr[TL];
	long long int fJwr[TL];
	long long int mJwr[TL];
	long long int bJwr[TL];
	long long int Jrr[TL];
	long long int Jdr[TL];
	long long int mbJrr[TL];
	long long int mbJdr[TL];
	long long int JTot;
	long long int MTot;
	long long int Mww,Mwd,Mdd,Mwr,Mrr,Mdr;
	long long int Vww,fVww,mVww,bVww,fVwd,mVwd,bVwd,Vdd,Vwr,fVwr,mVwr,bVwr,Vrr,Vdr;
	//int Fww,Fwd,Fdd,Fwr,Frr,Fdr;
	long long int Fwwww,Fwwwd,Fwwdd,Fwwwr,Fwwrr,Fwwdr;
	long long int fFwwww,fFwwwd,fFwwdd,fFwwwr,fFwwrr,fFwwdr;
	long long int mFwwww,mFwwwd,mFwwdd,mFwwwr,mFwwrr,mFwwdr;
	long long int bFwwww,bFwwwd,bFwwdd,bFwwwr,bFwwrr,bFwwdr;
	long long int fFwdww,fFwdwd,fFwddd,fFwdwr,fFwdrr,fFwddr;
	long long int mFwdww,mFwdwd,mFwddd,mFwdwr,mFwdrr,mFwddr;
	long long int bFwdww,bFwdwd,bFwddd,bFwdwr,bFwdrr,bFwddr;
	long long int Fddww,Fddwd,Fdddd,Fddwr,Fddrr,Fdddr;
	long long int Fwrww,Fwrwd,Fwrdd,Fwrwr,Fwrrr,Fwrdr;
	long long int fFwrww,fFwrwd,fFwrdd,fFwrwr,fFwrrr,fFwrdr;
	long long int mFwrww,mFwrwd,mFwrdd,mFwrwr,mFwrrr,mFwrdr;
	long long int bFwrww,bFwrwd,bFwrdd,bFwrwr,bFwrrr,bFwrdr;
	long long int Frrww,Frrwd,Frrdd,Frrwr,Frrrr,Frrdr;
	long long int Fdrww,Fdrwd,Fdrdd,Fdrwr,Fdrrr,Fdrdr;
	long long int AesFwwww,AesFwwwd,AesFwwdd,AesFwwwr,AesFwwrr,AesFwwdr;
	long long int AesfFwwww,AesfFwwwd,AesfFwwdd,AesfFwwwr,AesfFwwrr,AesfFwwdr;
	long long int AesmFwwww,AesmFwwwd,AesmFwwdd,AesmFwwwr,AesmFwwrr,AesmFwwdr;
	long long int AesbFwwww,AesbFwwwd,AesbFwwdd,AesbFwwwr,AesbFwwrr,AesbFwwdr;
	long long int AesfFwdww,AesfFwdwd,AesfFwddd,AesfFwdwr,AesfFwdrr,AesfFwddr;
	long long int AesmFwdww,AesmFwdwd,AesmFwddd,AesmFwdwr,AesmFwdrr,AesmFwddr;
	long long int AesbFwdww,AesbFwdwd,AesbFwddd,AesbFwdwr,AesbFwdrr,AesbFwddr;
	long long int AesFddww,AesFddwd,AesFdddd,AesFddwr,AesFddrr,AesFdddr;
	long long int AesFwrww,AesFwrwd,AesFwrdd,AesFwrwr,AesFwrrr,AesFwrdr;
	long long int AesfFwrww,AesfFwrwd,AesfFwrdd,AesfFwrwr,AesfFwrrr,AesfFwrdr;
	long long int AesmFwrww,AesmFwrwd,AesmFwrdd,AesmFwrwr,AesmFwrrr,AesmFwrdr;
	long long int AesbFwrww,AesbFwrwd,AesbFwrdd,AesbFwrwr,AesbFwrrr,AesbFwrdr;
	long long int AesFrrww,AesFrrwd,AesFrrdd,AesFrrwr,AesFrrrr,AesFrrdr;
	long long int AesFdrww,AesFdrwd,AesFdrdd,AesFdrwr,AesFdrrr,AesFdrdr;
	long long int LDMFwwww,LDMFwwwd,LDMFwwdd,LDMFwwwr,LDMFwwrr,LDMFwwdr;
	long long int LDMfFwwww,LDMfFwwwd,LDMfFwwdd,LDMfFwwwr,LDMfFwwrr,LDMfFwwdr;
	long long int LDMmFwwww,LDMmFwwwd,LDMmFwwdd,LDMmFwwwr,LDMmFwwrr,LDMmFwwdr;
	long long int LDMbFwwww,LDMbFwwwd,LDMbFwwdd,LDMbFwwwr,LDMbFwwrr,LDMbFwwdr;
	long long int LDMfFwdww,LDMfFwdwd,LDMfFwddd,LDMfFwdwr,LDMfFwdrr,LDMfFwddr;
	long long int LDMmFwdww,LDMmFwdwd,LDMmFwddd,LDMmFwdwr,LDMmFwdrr,LDMmFwddr;
	long long int LDMbFwdww,LDMbFwdwd,LDMbFwddd,LDMbFwdwr,LDMbFwdrr,LDMbFwddr;
	long long int LDMFddww,LDMFddwd,LDMFdddd,LDMFddwr,LDMFddrr,LDMFdddr;
	long long int LDMFwrww,LDMFwrwd,LDMFwrdd,LDMFwrwr,LDMFwrrr,LDMFwrdr;
	long long int LDMfFwrww,LDMfFwrwd,LDMfFwrdd,LDMfFwrwr,LDMfFwrrr,LDMfFwrdr;
	long long int LDMmFwrww,LDMmFwrwd,LDMmFwrdd,LDMmFwrwr,LDMmFwrrr,LDMmFwrdr;
	long long int LDMbFwrww,LDMbFwrwd,LDMbFwrdd,LDMbFwrwr,LDMbFwrrr,LDMbFwrdr;
	long long int LDMFrrww,LDMFrrwd,LDMFrrdd,LDMFrrwr,LDMFrrrr,LDMFrrdr;
	long long int LDMFdrww,LDMFdrwd,LDMFdrdd,LDMFdrwr,LDMFdrrr,LDMFdrdr;
	long long int MoveFwwww,MoveFwwwd,MoveFwwdd,MoveFwwwr,MoveFwwrr,MoveFwwdr;
	long long int MovefFwwww,MovefFwwwd,MovefFwwdd,MovefFwwwr,MovefFwwrr,MovefFwwdr;
	long long int MovemFwwww,MovemFwwwd,MovemFwwdd,MovemFwwwr,MovemFwwrr,MovemFwwdr;
	long long int MovebFwwww,MovebFwwwd,MovebFwwdd,MovebFwwwr,MovebFwwrr,MovebFwwdr;
	long long int MovefFwdww,MovefFwdwd,MovefFwddd,MovefFwdwr,MovefFwdrr,MovefFwddr;
	long long int MovemFwdww,MovemFwdwd,MovemFwddd,MovemFwdwr,MovemFwdrr,MovemFwddr;
	long long int MovebFwdww,MovebFwdwd,MovebFwddd,MovebFwdwr,MovebFwdrr,MovebFwddr;
	long long int MoveFddww,MoveFddwd,MoveFdddd,MoveFddwr,MoveFddrr,MoveFdddr;
	long long int MoveFwrww,MoveFwrwd,MoveFwrdd,MoveFwrwr,MoveFwrrr,MoveFwrdr;
	long long int MovefFwrww,MovefFwrwd,MovefFwrdd,MovefFwrwr,MovefFwrrr,MovefFwrdr;
	long long int MovemFwrww,MovemFwrwd,MovemFwrdd,MovemFwrwr,MovemFwrrr,MovemFwrdr;
	long long int MovebFwrww,MovebFwrwd,MovebFwrdd,MovebFwrwr,MovebFwrrr,MovebFwrdr;
	long long int MoveFrrww,MoveFrrwd,MoveFrrdd,MoveFrrwr,MoveFrrrr,MoveFrrdr;
	long long int MoveFdrww,MoveFdrwd,MoveFdrdd,MoveFdrwr,MoveFdrrr,MoveFdrdr;
	long long int MoveMww,MoveMwd,MoveMdd,MoveMwr,MoveMrr,MoveMdr;
	long double comp;
	char arrive;
	char release;
	long double mate_rate;
	vector<int> connecIND;
	vector<double> connecW;
	double TotW;
	char type;
	char Fixed;
	double distOrigin;
	double WaterTemp,WaterPerm;
	double alpha0;
};
/*----------------------------------------------------------------------------------------*/



/*----------------struct containt simulation timekeeping parameters-------------------------*/
struct Times
{
	int maxT;
	int initiate;
	int rec;
	int NumRuns;
	int interval;
	int yearnow;
};
/*----------------------------------------------------------------------------------------*/

/*---------------------struct containing model parameters----------------------------------*/
struct Pars
	{
	double muJ,muA,d,Fgamma,Mgamma,beta,theta,Frho,Mrho,xi,em,ef,LD,dL,muLD,psi,muAES,bias;
	double alpha0,alpha1,alpha2,delta,phi,kappa,al0var,mu,sig;
	int set,offset;
	int t_disp1,t_disp2,t_disp3,t_disp4;
	int t_hide1,t_hide2,t_wake1,t_wake2;
	double fwwww[16]; double fwwwd[16]; double fwwdd[16]; double fwwwr[16]; double fwwrr[16]; double fwwdr[16];
	double ffwwww[16]; double ffwwwd[16]; double ffwwdd[16]; double ffwwwr[16]; double ffwwrr[16]; double ffwwdr[16];
	double mfwwww[16]; double mfwwwd[16]; double mfwwdd[16]; double mfwwwr[16]; double mfwwrr[16]; double mfwwdr[16];
	double bfwwww[16]; double bfwwwd[16]; double bfwwdd[16]; double bfwwwr[16]; double bfwwrr[16]; double bfwwdr[16];
	double ffwdww[16]; double ffwdwd[16]; double ffwddd[16]; double ffwdwr[16]; double ffwdrr[16]; double ffwddr[16];
	double mfwdww[16]; double mfwdwd[16]; double mfwddd[16]; double mfwdwr[16]; double mfwdrr[16]; double mfwddr[16];
	double bfwdww[16]; double bfwdwd[16]; double bfwddd[16]; double bfwdwr[16]; double bfwdrr[16]; double bfwddr[16];
	double fddww[16]; double fddwd[16]; double fdddd[16]; double fddwr[16]; double fddrr[16]; double fdddr[16];
	double fwrww[16]; double fwrwd[16]; double fwrdd[16]; double fwrwr[16]; double fwrrr[16]; double fwrdr[16];
	double ffwrww[16]; double ffwrwd[16]; double ffwrdd[16]; double ffwrwr[16]; double ffwrrr[16]; double ffwrdr[16];
	double mfwrww[16]; double mfwrwd[16]; double mfwrdd[16]; double mfwrwr[16]; double mfwrrr[16]; double mfwrdr[16];
	double bfwrww[16]; double bfwrwd[16]; double bfwrdd[16]; double bfwrwr[16]; double bfwrrr[16]; double bfwrdr[16];
	double frrww[16]; double frrwd[16]; double frrdd[16]; double frrwr[16]; double frrrr[16]; double frrdr[16];
	double fdrww[16]; double fdrwd[16]; double fdrdd[16]; double fdrwr[16]; double fdrrr[16]; double fdrdr[16];
	};
/*----------------------------------------------------------------------------------------*/

/*-------------------------------function declarations-------------------------------------*/
		void RunOnceInt(double);
		void record(int);
		void RunMaxT();
		void RunNReps(int);
		void initiate();
		void PatPopulate(int,char);
		void OneStep(int);
		void PutDriver();
		void PutDriverPat(int);
		void JuvGetOlder();
		void VirginsMate();
		void AdultsMove();
		void LDM(char);
		void Hide();
		void Wake(int);
		void LayEggs();
		void JuvEmerge();
		void AdultsDie();
		void UpdateComp(int);
		void UpdateMate();
		void UpdateConnec(int);
		void SetFertility();
		int random_poisson(double);
		int random_binomial(int,double);
		double random_normal(double,double);
		int* random_multinom(int,long long int[6]);
		int* random_multinomEqualProb(int,int);
		int* random_multinom_var(int,int,double*,double);
		double dist(double,double,double,double);

/*----------------------------------------------------------------------------------------*/

/*--------------------------code for random number generation------------------------------*/
// Define 32 bit signed and unsigned integers.
// GieIange these definitions, if necessary, to match a particular platform
#if defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS)
   // 16 bit systems use long int for 32 bit integer
   typedef long int           int32;   // 32 bit signed integer
   typedef unsigned long int  uint32;  // 32 bit unsigned integer
#else
   // Most other systems use int for 32 bit integer
   typedef int                int32;   // 32 bit signed integer
   typedef unsigned int       uint32;  // 32 bit unsigned integer
#endif

// Define 64 bit signed and unsigned integers, if possible
#if (defined(__WINDOWS__) || defined(_WIN32)) && (defined(_MSC_VER) || defined(__INTEL_COMPILER))
   // Microsoft and other compilers under Windows use __int64
   typedef __int64            int64;   // 64 bit signed integer
   typedef unsigned __int64   uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#elif defined(__unix__) && (defined(_M_IX86) || defined(_M_X64))
   // Gnu and other compilers under Linux etc. use long long
   typedef long long          int64;   // 64 bit signed integer
   typedef unsigned long long uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#else
   // 64 bit integers not defined
   // You may include definitions for other platforms here
#endif


/***********************************************************************
System-specific user interface functions
***********************************************************************/

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(char * ErrorText);     // System-specific error reporting (userintf.cpp)


/***********************************************************************
Define random number generator classes
***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
#if 0
   // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else
   // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
public:
   CRandomMersenne(uint32 seed) {      // Constructor
      RandomInit(seed); LastInterval = 0;}
   void RandomInit(uint32 seed);       // Re-seedSetDirectory["~/Dropbox/YDriveBurkina/codesB/DoubleSex/Mosaic"];
   void RandomInitByArray(uint32 seeds[], int length); // Seed by more than 32 bits
   int IRandom (int min, int max);     // Output random integer
   int IRandomX(int min, int max);     // Output random integer, exact
   double Random();                    // Output random float
   uint32 BRandom();                   // Output random bits
private:
   void Init0(uint32 seed);            // Basic initialization procedure
   uint32 mt[MERS_N];                  // State vector
   int mti;                            // Index into mt
   uint32 LastInterval;                // Last interval length for IRandomX
   uint32 RLimit;                      // Rejection limit used by IRandomX
   enum TArch {LITTLE_ENDIAN1, BIG_ENDIAN1, NONIEEE}; // Definition of architecture
   TArch Architecture;                 // Conversion to float depends on architecture
};


class CRandomMother {             // Encapsulate random number generator
public:
   void RandomInit(uint32 seed);       // Initialization
   int IRandom(int min, int max);      // Get integer random number in desired interval
   double Random();                    // Get floating point random number
   uint32 BRandom();                   // Output random bits
   CRandomMother(uint32 seed) {   // Constructor
      RandomInit(seed);}
protected:
   uint32 x[5];                        // History buffer
};

#endif

