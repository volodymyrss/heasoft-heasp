// definitions for heasp library

// all the includes

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <stdexcept>
#include <ctime>
#include <valarray>
#include <vector>

#include <CCfits/CCfits>

// use the std and CCfits namespaces

using namespace std;
using namespace CCfits;

// define Integer and Real

typedef int Integer;
typedef float Real;


#define HAVE_HEASP 1

// set up error statuses

enum{OK, NoSuchFile, NoData, NoChannelData, NoStatError, CannotCreate,
     NoEnergLo, NoEnergHi, NoSpecresp, NoEboundsExt, NoEmin, NoEmax,
     NoMatrixExt, NoNgrp, NoFchan, NoNchan, NoMatrix, CannotCreateMatrixExt,
     CannotCreateEboundsExt, InconsistentGrouping, InconsistentEnergies,
     InconsistentChannels, InconsistentUnits, UnknownXUnits, UnknownYUnits,
     InconsistentNumelt, InconsistentNumgrp};

const string SPerrorNames[] = 
  {"OK", "NoSuchFile", "NoData", "NoChannelData", "NoStatError", "CannotCreate",
   "NoEnergLo", "NoEnergHi", "NoSpecresp", "NoEboundsExt", "NoEmin", "NoEmax",
   "NoMatrixExt", "NoNgrp", "NoFchan", "NoNchan", "NoMatrix", "CannotCreateMatrixExt",
   "CannotCreateEboundsExt", "InconsistentGrouping", "InconsistentEnergies",
   "InconsistentChannels", "InconsistentUnits", "UnknownXUnits", "UnknownYUnits",
   "InconsistentNumelt", "InconsistentNumgrp"};

// the error message stack

static vector<string> SPerrorStack;

// some useful conversion factors

#define KEVTOA 12.3984191
#define KEVTOHZ 2.4179884076620228e17
#define KEVTOERG 1.60217733e-9
#define KEVTOJY 1.60217733e14
