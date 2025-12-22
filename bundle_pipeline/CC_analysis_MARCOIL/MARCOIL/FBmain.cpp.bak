/***********************************************
************************************************
**                                            **
** Marcoil                                    **
** a program for predicting coiled-coil       **
** domains in protein sequences               **     
** Copyright (C) 2001 Mauro C. Delorenzi      **
**                                            **
************************************************
************************************************/

#include "FBmain.h"

int  	 emProb[kStates][21], trProb[kStates][kStates] ;

extern double CutoffForWriting;
extern bool	   verbose, optCmatrix;	
extern  int	   Code[27], Decode[22];
extern bool FBDdetailsD;
extern bool FBDdetailsS, FBDdetailsC, FBDdetailsL;
extern bool bWriteSeq;
extern double thresholds[6];
extern int    thresholdNb;

extern  int   seqNb, seqLen;
extern  int   seq[kMaxSeqLen+1];
extern  bool  completedFile;
extern  char  SeqName[kMaxSeqName+2];
extern  FILE *cFpWriteDom; 
extern  FILE *cFpWritePP, *cFpWritePPD,*cFpWriteCP; 
extern  FILE *cFpReadSeq;

/***********************************
* function-declarations  */

void  WriteProbProfile ( double PosStateProb[][kStates] );
void  DoFB (int seqLen, const int Seq[] ) ;
/*********************************** */

// ---------------------------------------------------------------------------
//		 FBmain 					
// ---------------------------------------------------------------------------

void FBmain(const char *transProbFile, const char *emissProbFile, char *seqFile )
{
int i, j, n, l;
FILE  * fprtr, * fprem,  * cFpWrite;
FILE  * fpwtr, * fpwem;
char 	*ReadFileName[11], *WriteFileName[11];
bool checkinput = false;

if (verbose)  cout  <<  "\n starting FBmain" <<  endl;
if  ((fprtr = fopen( transProbFile,"r" ))  == NULL )
	{ cerr << "\n could not open the file for transition P. with adress  " << transProbFile << "\n"<< endl; exit(1); }
else {ReadTransProb(fprtr);}  	
if (verbose)  cout << "\n finished ReadTransProb";
fclose( fprtr );


if  ((fprem = fopen(emissProbFile,"r" ))  == NULL )
	{ cerr << "\n could not open the file for emission P. with adress  " << emissProbFile << "\n" << endl; exit(1); }
else {ReadEmissProb(fprem);}  	
if (verbose)  cout << "\n finished ReadEmissProb";
fclose(fprem );

if (verbose)  cout << "\n optCmatrix = " <<  optCmatrix << endl;  

if (optCmatrix) //use only state 0 emiss from previous file,than use propensities from this another file
	{cout << "\n optCmatrix => reading lupas matrix with ReadProp";   ReadProp();}

if (checkinput)
	{
	fpwtr = fopen( "Outputs/Checks/W3.TransProb","w" );
	WriteTransProb(fpwtr); 
	fclose( fpwtr );
	if (verbose)  cout << "\n finished WriteTransProb" <<  endl;

	fpwem = fopen( "Outputs/Checks/W2.EmissProb","w" );
	WriteEmissProb(fpwem); 
	fclose( fpwem );
	if (verbose)  cout << "\n finished WriteEmissProb"<<  endl;
	}

if (FBDdetailsL)
	{
	if  ( (cFpWritePP = fopen("Outputs/ProbList","w" ))  == NULL )
		{ cerr << "\n could not open the file Outputs/ProbProfile for writing \n";}
	else   fprintf(  cFpWritePP, " COILED-COIL PROBABILITY LIST PER RESIDUE\n\n");
	}

if (FBDdetailsS)
	 {
	 if  ( (cFpWritePPD = fopen("Outputs/ProbPerState","w" ))  == NULL )
		 { cerr << "\n could not open the file Outputs/ProbPerState for writing \n";}
	else   fprintf(  cFpWritePPD, " COILED-COIL PROBABILITY LIST PER RESIDUE, PER EACH HEPTAD POSITION\n\n");
	}

if (FBDdetailsC)
	{
	 if  ( (cFpWriteCP = fopen("Outputs/CompactProfile","w" ))  == NULL )
		 { cerr << "\n could not open the file Outputs/CompactProfile for writing \n";}
	else   fprintf(  cFpWriteCP, " COILED-COIL PROBABILITY PER RESIDUE, COMPACT REPRESENTATION\n\n");
	}

if (FBDdetailsD)
	{
	 if  ( (cFpWriteDom = fopen("Outputs/Domains","w" ))  == NULL )
		 { cerr << "\n could not open the file Outputs/Domains for writing \n";}
	else   fprintf(  cFpWriteDom, " PREDICTED COILED-COIL DOMAINS: OVERVIEW\n\n");
	}
	
if  ((cFpReadSeq = fopen( seqFile, "r" ))  == NULL )
	{ cerr << "\n could not open the file for the sequences with adress  " << seqFile << "\n"; exit(1); }

if  ( (cFpWrite = fopen("Outputs/Checks/FBmain.check","w" ))  == NULL )
	{ cerr << "\n could not open the file Outputs/Checks/FBmain.check for writing \n";}

if (verbose)  cout << "\n  opened all Files" <<  endl;

completedFile = false; seqNb = 0;
while (! completedFile)
	{
	seqNb++;seqLen = 0;
	if ( ! (ReadSeq (seqNb, seqLen, seq, completedFile, SeqName) )   )
		{
		if (verbose)  cout << "\n ReadSeq done" <<  endl;
		if (seqLen > 0)	{ WriteWarning(seqNb,1);}	
		else	{seqNb--;}
		} 
	if (seqLen > 0)
		{ 
		if (verbose)   cout << "\n read seq-Nb " << seqNb << "  name = " << SeqName << endl;
		if (bWriteSeq) WriteSeq (seqNb, seqLen, seq, SeqName);
		else WriteSeqId(seqNb, SeqName);
		if (verbose)   cout << "\n processing seq-Nb " << seqNb << "  name = " << SeqName << endl;
		DoFB (seqLen, seq);
		}
	if ( (seqNb % kperiodStdout1) == 0)	 
		{
		if (seqNb == kperiodStdout1) printf ("One point %d sequences, one line %d\n", kperiodStdout1,kperiodStdout2); 
		printf (".");  fflush (stdout);
		}
	if ( (seqNb % kperiodStdout2) == 0)  printf ("$\n");  	 
	if ( (seqNb % kperiodStdout3) == 0)  printf ( "\n processed seq-nb %d \n", seqNb);
	}

fclose( cFpWrite );
	if (FBDdetailsL)  fclose( cFpWritePP );
	if (FBDdetailsS)  fclose( cFpWritePPD);
	if (FBDdetailsC)  fclose( cFpWriteCP);
	if (FBDdetailsD)  fclose( cFpWriteDom);

cout  <<  "\n processed  " << seqNb << "sequences" <<  endl;; 
}
/* FBMAIN
************************************/
/* 
************************************/
// ---------------------------------------------------------------------------
//		 DoFB 					
// ---------------------------------------------------------------------------
void  DoFB (int seqLen, const int Seq[] ) 
{
int      i;
float    BgrProbProfile[seqLen+2];
double   PosStateProb[seqLen+2][kStates];

if (verbose)   cout << "\n starting FM-algo"  << endl;
FBalgorithm (seqNb, seqLen, Seq , PosStateProb ); 
if (verbose)   cout << "\n after FM-algo  " << seqNb << endl;
if (verbose)   cout << "\n WRITING PROB PROFILE   "  << endl;
if  ((FBDdetailsL) || (FBDdetailsS)  ||  (FBDdetailsC) )      WriteProbProfile(PosStateProb);
for (i=1; i<= seqLen; i++)	BgrProbProfile[i] = (float) PosStateProb[i][0];
if (FBDdetailsD)    ParseIntoDomains (seqLen, BgrProbProfile, thresholdNb, thresholds);

}

/* DoFB
************************************/
// ---------------------------------------------------------------------------
//		 WriteProbProfile 					
// ---------------------------------------------------------------------------
void  WriteProbProfile ( double PosStateProb[][kStates] )
{
int i, st, k, bestPathM[seqLen+2] , bestPathcc[seqLen+2];  char c, mode;
double stateProb[9], max, ProbForWritingFBDdetails;
bool writing;

ProbForWritingFBDdetails = CutoffForWriting; 

if  (FBDdetailsS)
	fprintf( cFpWritePPD, "\n coiled-coil probability in percent per each heptad position");

for (i=1; i<=seqLen; i++)
	{
	stateProb[0] = PosStateProb[i][0];
	if  ((stateProb[0] <= (1- ProbForWritingFBDdetails))   || (i==1) || (i==seqLen)) 
		{
		max = -1;  bestPathM[i] = 23;  bestPathcc[i] = 23; 
		for (st=1; st<=7; st++)  
			{
			stateProb[st] = 0;
			for (k=0; k<=8; k++)  {stateProb[st] += PosStateProb[i][st+7*k]; 
						    }
			if (stateProb[st] > max)  {max = stateProb[st]; bestPathM[i] = st;}
			}
		bestPathcc[i] = bestPathM[i]  ; 	
		if  ( stateProb[0] > max)   {bestPathM[i] = 0;}   
		if  (FBDdetailsS)			
			{
			fprintf( cFpWritePPD, "\n     state:   a      b       c       d       e       f       g      cc(sum) ");
			fprintf( cFpWritePPD, "\npos=%4d% c", i, Decode[seq[i]] );  
			for (st=1; st<=7; st++)	 fprintf (cFpWritePPD, "  %6.2f\%" , 100 * stateProb[st]);
			fprintf( cFpWritePPD, "  %6.2f\%"  , 100 * (1 - stateProb[0] )  );
			fprintf( cFpWritePPD, " CC-ph: %c ;",  (bestPathcc[i] +96)   );
			if (bestPathM[i] == 0)  fprintf( cFpWritePPD, " max-st: 0");
			else   fprintf( cFpWritePPD, " max-st: %c", (bestPathM[i] +96) );
			}
	}	}

if  (FBDdetailsS)
	{
      fprintf(cFpWritePPD,"\n\n******************************************************************************"); 
	fprintf(cFpWritePPD,"**********************************************************************************\n"); 
	}

if  (FBDdetailsC)
	{
	 mode = 0;
	fprintf( cFpWriteCP, "\n\n [coiled-coil probability in percent] and heptad position with highest probability\n");
	 for (i=1; i<=seqLen; i++)
		 {
		 stateProb[0] = PosStateProb[i][0];
		 if ((stateProb[0] <= (1- ProbForWritingFBDdetails) ) || (i==1) || (i==seqLen)) 
			 {
			 writing = true;
			 mode ++;
			 fprintf(cFpWriteCP, " %4d%c[",i,Decode[seq[i]]);
			 if (stateProb[0] >= 0.9005) fprintf(cFpWriteCP, "0");
			 if (stateProb[0] >= 0.0005)  fprintf(cFpWriteCP,"%3.1f]%c", 100*(1 - stateProb[0]),(bestPathcc[i]+96));
			 else fprintf(cFpWriteCP,"%3.0f.]%c", 100*(1 - stateProb[0]),(bestPathcc[i]+96));
			 if (i==seqLen) fprintf(cFpWriteCP, "  **");
			 else
				 {
				 if  (mode ==7)  {mode = 0; fprintf(cFpWriteCP, "\n");}
				 else fprintf(cFpWriteCP, "");
			 }	}
		 else  {
			 if (writing)  {writing = false; fprintf(cFpWriteCP, "\n[..]\n"); mode=0;}
		 }	}
	fprintf(cFpWriteCP,"\n\n******************************************************************************"); 
	fprintf(cFpWriteCP,"**********************************************************************************\n"); 
	}

if (FBDdetailsL)
	{ 
	fprintf(  cFpWritePP, "  cc-probability in percent and best heptad phase\n");
	for (i=1; i<= seqLen; i++)
		{
		if ((  PosStateProb[i][0] <= (1- ProbForWritingFBDdetails))  || (i==1) || (i==seqLen))
			{
			fprintf(  cFpWritePP, "%4d %c ",  i , Decode[seq[i]] );
			fprintf(  cFpWritePP, "%5.1f %c\n", 100.0 * ((double) 1.0 - PosStateProb[i][0]), bestPathcc[i]+96); 
		}	}
	fprintf(cFpWritePP,"\n\n******************************************************************************"); 
	fprintf(cFpWritePP,"**********************************************************************************\n"); 
	}
}
/*  WriteProbProfile
************************************/
