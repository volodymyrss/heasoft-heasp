// Some useful general utility routines

#ifndef HAVE_SPutils
#include "SPutils.h"
#endif

// Arrays to store unit conversion information

static size_t nValidXUnits = 8;
static string validXUnits[] = {"keV",      "MeV",        "GeV",        "Hz", 
			"angstrom", "cm",         "micron",     "nm"};
static Real xConvFactor[] =   {1.0,        1000.0,       1.0e6,        1.0/KEVTOHZ, 
			KEVTOA,     KEVTOA/1.0e8, KEVTOA/1.0e4, KEVTOA/10.0};
static bool xIsWave[] =       {false,      false,        false,        false, 
			true,       true,         true,         true};

static size_t nValidYUnits = 16;
static string validYUnits[] = {
  "ph/cm^2/s",      "ph/cm^2/s/keV",  "ph/cm^2/s/Mev",  "ph/cm^2/s/Gev", 
  "ph/cm^2/s/Hz",   "ph/cm^2/s/A",    "ph/cm^2/s/cm",   "ph/cm^2/s/um",   
  "ph/cm^2/s/nm",   "ergs/cm^2/s",    "ergs/cm^2/s/Hz", "ergs/cm^2/s/A",  
  "ergs/cm^2/s/cm", "ergs/cm^2/s/um", "ergs/cm^2/s/nm", "Jy"};
static Real yConvFactor[] = {
  1.0,              1.0,              1.0e-3,           1.0e-6,
  KEVTOHZ,          KEVTOA,           KEVTOA/1.0e8,     KEVTOA/1.0e4,
  KEVTOA/10.0,      KEVTOERG,         KEVTOHZ/KEVTOERG, KEVTOA/KEVTOERG,
  KEVTOA/1.0e8/KEVTOERG, KEVTOA/1.0e4/KEVTOERG, KEVTOA/10.0/KEVTOERG, KEVTOHZ/KEVTOJY};
static bool yIsEnergy[] = {
  false,            false,            false,            false, 
  false,            false,            false,            false, 
  false,            true,             true,             true, 
  true,             true,             true,             true};
static bool yIsPerEnergy[] = {
  false,            true,             true,             true, 
  true,             false,            false,            false, 
  false,            false,            true,             false, 
  false,            false,            false,            true};
static bool yIsPerWave[] = {
  false,            false,            false,            false, 
  false,            true,             true,             true, 
  true,             false,            false,            true, 
  true,             true,             true,             false};

// Read the units associated with a column

void SPreadColUnits(ExtHDU& ext, string ColName, string& Units)
{
  // get the column index

  try {
    Column& Col = ext.column(ColName);
    int index = Col.index();

    // construct the keyword name

    stringstream nameStream;
    nameStream << "TUNIT" << index;

    string DefString(" ");
    Units = SPreadKey(ext, nameStream.str(), DefString);

  } catch(Table::NoSuchColumn&){
    Units = " ";
  }

  return;
}

// Write the units associated with a column

void SPwriteColUnits(Table& table, string ColName, string Units)
{
  // get the column index

  try {
    Column& Col = table.column(ColName);
    int index = Col.index();

    // construct the keyword name

    stringstream nameStream;
    nameStream << "TUNIT" << index;
    SPwriteKey(table, nameStream.str(), Units, " ");

  } catch(Table::NoSuchColumn&){
  }

  return;
}

// Return a FITS string specification for the longest string in the input
// vector<string>

string SPstringTform(const vector<string>& Data)
{

  size_t length = Data[0].size();
  for (size_t i=1; i<Data.size(); i++) {
    if ( Data[i].size() > length ) length = Data[i].size();
  }

  stringstream tmpStream;
  tmpStream << length;

  return tmpStream.str();
}

// copy from infile to outfile all HDUs which are not manipulated by this library 

Integer SPcopyHDUs(string infile, string outfile)
{

  size_t NumIgnore(5);
  string IgnoreNames[5] = {"SPECTRUM", "SPECRESP", "EBOUNDS", 
			   "MATRIX", "SPECRESP MATRIX"};

  auto_ptr<FITS> pInfile(0);
  auto_ptr<FITS> pOutfile(0);

  // open the input file

  try {
    pInfile.reset(new FITS(infile));
  } catch(...) {
    string msg = "Failed to read "+infile;
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  FITS InFITS(infile);

  // open the output file. If it doesn't exist then create it by copying the
  // primary header from the input file

  try {
    pOutfile.reset(new FITS(outfile, Write));
  } catch(...) {
    try {
      pOutfile.reset(new FITS(outfile, infile));
    } catch(...) {
      string msg = "Failed to create "+outfile+" for copying";
      SPreportError(CannotCreate, msg);
      return(CannotCreate);
    }
  }

  // Loop through the extensions in the input file

  bool done(false);
  Integer extNum(1);

  while (!done) {

    try {
      ExtHDU& inExtension = pInfile->extension(extNum);
      string ExtName = inExtension.name();

      bool found(false);
      for (size_t i=0; i<NumIgnore; i++) {
  	if ( ExtName == IgnoreNames[i] ) found = true;
      }
      if (!found) {
  	pOutfile->copy(pInfile->extension(extNum));
      }
      extNum++;
    } catch(...) {
      done = true;
    }

  }

  return(0);
}

// copy non-critical keywords from infile to outfile for the HDUnumber instance
// of the HDUname HDU.

Integer SPcopyKeys(string infile, string outfile, string HDUname, Integer HDUnumber)
{

  const vector<string> hduKeys;
  const vector<string> primaryKey;

  auto_ptr<FITS> pInfile(0);
  auto_ptr<FITS> pOutfile(0);

  // open the input file to the requested extension

  try {
    pInfile.reset(new FITS(infile,Read,HDUname,false,hduKeys,primaryKey,(int)HDUnumber));
  } catch(...) {
    string msg = "Failed to read "+HDUname+" in "+infile;
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  HDU *inHDU = &pInfile->extension(HDUname);

  // read in the keywords

  inHDU->readAllKeys();

  // and the output file to the requested extension

  try {
    pOutfile.reset(new FITS(outfile,Write,HDUname,false,hduKeys,primaryKey,(int)HDUnumber));
  } catch(...) {
    string msg = "Failed to open "+HDUname+" in "+outfile+" for writing.";
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  ExtHDU& outHDU = pOutfile->extension(HDUname);

  // copy keywords in categories TYP_CMPRS_KEY (20), TYP_CKSUM_KEY (100), TYP_WCS_KEY (110), 
  // TYP_REFSYS_KEY (120), and TYP_USER_KEY (150). This choice is currently hardwired into
  // CCfits - it may be necessary to change this at some point.

  outHDU.copyAllKeys(inHDU);

  return(0);
}

// Check whether valid X units

bool isValidXUnits(string xUnits)
{
  bool found = false;
  for (size_t i=0; i<nValidXUnits; i++) {
    if ( xUnits == validXUnits[i] ) found = true;
  }
  return(found);
}

// Calculate the unit conversion factor for energy/wavelength

Integer calcXfactor(string xUnits, bool& isWave, Real& xFactor)
{

  if ( !isValidXUnits(xUnits) ) {
    string msg = xUnits+" is not a valid unit for conversion";
    SPreportError(UnknownXUnits, msg);
    return(UnknownXUnits);
  }    

  for (size_t i=0; i<nValidXUnits; i++) {
    if ( xUnits == validXUnits[i] ) {
      isWave = xIsWave[i];
      xFactor = xConvFactor[i];
    }
  }
  return(OK);
}

// Check whether valid Y units

bool isValidYUnits(string yUnits)
{
  bool found = false;
  for (size_t i=0; i<nValidYUnits; i++) {
    if ( yUnits == validYUnits[i] ) found = true;
  }
  return(found);
}

// Calculate the unit conversion factor for the flux

Integer calcYfactor(string yUnits, bool& isEnergy, bool& perWave, bool& perEnergy, Real& yFactor)
{
  if ( !isValidYUnits(yUnits) ) {
    string msg = yUnits+" is not a valid unit for conversion";
    SPreportError(UnknownYUnits, msg);
    return(UnknownYUnits);
  }    

  for (size_t i=0; i<nValidYUnits; i++) {
    if ( yUnits == validYUnits[i] ) {
      isEnergy = yIsEnergy[i];
      perWave = yIsPerWave[i];
      perEnergy = yIsPerEnergy[i];
      yFactor = yConvFactor[i];
    }
  }
  return(OK);
}

// Add to the error stack

void SPreportError(int errorNumber, string optionalString)
{
  string msg = SPerrorNames[errorNumber] + ": " + optionalString;
  SPerrorStack.push_back(msg);

  std::cout << "error: " << msg << std::endl;

  return;
}

// Output the error stack

string SPgetErrorStack()
{
  string msg;
  for (size_t i=0; i<SPerrorStack.size(); i++) {
    msg += SPerrorStack[i] + "\n";
  }

  return msg;
}

// Clear the error stack

void SPclearErrorStack()
{
  SPerrorStack.clear();
  return;
}
