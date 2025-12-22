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

#include <cstdlib> // for access to the Unix-system
#include <cstdio>
#include <stdlib.h> 

#include <iostream>
using namespace std; // for  cin, cout, <<  etc
#include <fstream> // for file streams like fout
#include <iomanip>  // for  endl;

#include "globals.h"

#include <cmath>
#include <cfloat>
#include <cstring>

/***********************************
* function-declarations  */

bool	ReadTransProb (FILE *cFrTrans);
bool	ReadEmissProb(FILE *cFrPar);
bool	ReadProp ();
void	  WriteEmissProb(FILE *cFwPar); 
void	  WriteTransProb (FILE	*cFwTrans); 

/* function-declarations 
************************************/
