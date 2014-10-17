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
  "ph/cm^2/s",      "ph/cm^2/s/keV",  "ph/cm^2/s/MeV",  "ph/cm^2/s/GeV", 
  "ph/cm^2/s/Hz",   "ph/cm^2/s/A",    "ph/cm^2/s/cm",   "ph/cm^2/s/um",   
  "ph/cm^2/s/nm",   "ergs/cm^2/s",    "ergs/cm^2/s/Hz", "ergs/cm^2/s/A",  
  "ergs/cm^2/s/cm", "ergs/cm^2/s/um", "ergs/cm^2/s/nm", "Jy"};
static Real yConvFactor[] = {
  1.0,              1.0,              1.0e-3,           1.0e-6,
  KEVTOHZ,          KEVTOA,           KEVTOA/1.0e8,     KEVTOA/1.0e4,
  KEVTOA/10.0,      1.0/KEVTOERG,     KEVTOHZ/KEVTOERG, KEVTOA/KEVTOERG,
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

// copy non-critical columns from infile to outfile for the HDUnumber instance
// of the HDUname HDU.

Integer SPcopyCols(string infile, string outfile, string HDUname, Integer HDUnumber)
{
  return SPcopyCols(infile, outfile, HDUname, HDUname, HDUnumber, HDUnumber);
}

// copy non-critical columns from infile to outfile from the HDUnumber instance
// of the HDUname HDU to the outHDUnumber instance of the outHDUname HDU.
// note that this routine does not do anything yet since it is awaiting a change
// to CCfits

Integer SPcopyCols(string infile, string outfile, string HDUname, string outHDUname, 
		   Integer HDUnumber, Integer outHDUnumber)
{

  // The list of column names that should not be copied because they are already
  // being handled by the appropriate class.

  size_t NumIgnore(18);
  string IgnoreNames[18] = {"ENERG_LO", "ENERG_HI", "SPECRESP", "CHANNEL", 
			    "COUNTS", "RATE", "STAT_ERR", "SYS_ERR", "QUALITY",
			    "GROUPING", "AREASCAL", "BACKSCAL", "E_MIN", "E_MAX",
			    "F_CHAN", "N_CHAN", "MATRIX", "ORDER"};

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

  ExtHDU *inHDU = &pInfile->extension(HDUname);

  // and the output file to the requested extension

  try {
    pOutfile.reset(new FITS(outfile,Write,outHDUname,false,hduKeys,primaryKey,(int)outHDUnumber));
  } catch(...) {
    string msg = "Failed to open "+outHDUname+" in "+outfile+" for writing.";
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  ExtHDU& outHDU = pOutfile->extension(outHDUname);

  // loop round columns in the input file

  for (size_t i=0; i<(size_t)inHDU->numCols(); i++) {

    Column& thisCol = inHDU->column(i);
    string colname = thisCol.name();

    bool reserved = false;
    for (size_t j=0; j<NumIgnore; j++) {
      if ( colname == IgnoreNames[j] ) reserved = true;
    }

    // this column is safe to copy so go ahead and do so

    if ( !reserved ) {

    }

  }

  return(0);
}



// copy non-critical keywords from infile to outfile for the HDUnumber instance
// of the HDUname HDU.

Integer SPcopyKeys(string infile, string outfile, string HDUname, Integer HDUnumber)
{
  return SPcopyKeys(infile, outfile, HDUname, HDUname, HDUnumber, HDUnumber);
}

// copy non-critical keywords from infile to outfile from the HDUnumber instance
// of the HDUname HDU to the outHDUnumber instance of the outHDUname HDU.

Integer SPcopyKeys(string infile, string outfile, string HDUname, string outHDUname, 
		   Integer HDUnumber, Integer outHDUnumber)
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
    pOutfile.reset(new FITS(outfile,Write,outHDUname,false,hduKeys,primaryKey,(int)outHDUnumber));
  } catch(...) {
    string msg = "Failed to open "+outHDUname+" in "+outfile+" for writing.";
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  ExtHDU& outHDU = pOutfile->extension(outHDUname);

  // copy keywords in categories TYP_CMPRS_KEY (20), TYP_CKSUM_KEY (100), TYP_WCS_KEY (110), 
  // TYP_REFSYS_KEY (120), and TYP_USER_KEY (150). This choice is currently hardwired into
  // CCfits - it may be necessary to change this at some point.

  outHDU.copyAllKeys(inHDU);

  return(0);
}

// write the creating program and version id string into the CREATOR keyword in the
// specified file and HDUname extension

Integer SPwriteCreator(string filename, string HDUname, string creator) {
  return(SPwriteCreator(filename, HDUname, 0, creator));
}

// write the creating program and version id string into the CREATOR keyword in the
// specified file and HDUnumber'th HDUname extension.

Integer SPwriteCreator(string filename, string HDUname, Integer HDUnumber, string creator) {

  const vector<string> hduKeys;
  const vector<string> primaryKey;

  auto_ptr<FITS> pOutfile(0);

  // open the output file to the requested extension

  try {
    pOutfile.reset(new FITS(filename,Write,HDUname,false,hduKeys,primaryKey,(int)HDUnumber));
  } catch(...) {
    string msg = "Failed to open "+HDUname+" in "+filename+" for writing.";
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  ExtHDU& outHDU = pOutfile->extension(HDUname);

  outHDU.addKey("CREATOR", creator, "Program which created this extension");

  return(0);
}

// find the numbers of any extensions containing keyword keyname=keyvalue

vector<Integer> SPfindExtensions(string filename, string keyname, string value, Integer& Status)
{
  vector<Integer> ExtNumbers;

  if ( Status != 0 ) return ExtNumbers;

  const vector<string> primaryKey;

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(filename,Read,false,primaryKey));
  } catch(...) {
    Status = NoSuchFile;
    string msg = "Failed to open "+filename;
    SPreportError(Status, msg);
    return(ExtNumbers);
  }

  // Not case sensitive so do everything in upper case
  string testValue(value);
  transform(testValue.begin(),testValue.end(),testValue.begin(),::toupper);

  // Loop through the extensions

  try {  
    int ExtNo=1;
    while (true) {
      ExtHDU& ext = pInfile->extension(ExtNo);

      string keyRead = SPreadKey(ext, keyname, (string)" ");
      transform(keyRead.begin(),keyRead.end(),keyRead.begin(),::toupper);

      if ( keyRead == testValue ) ExtNumbers.push_back(ExtNo);
      ExtNo++;
    }
  } catch(...) {
  }

  return ExtNumbers;
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

// handy function to read an ascii file and place each row into its own
// element of a vector<string>

vector<string> SPreadStrings(const string & filename)
{
  vector<string> strings;
  ifstream File;
  string input;

  File.open(filename.c_str());
  if ( File.is_open() ) {
    getline(File, input);
    while ( !File.eof() ) {
      // remove any \r or \n which is still tacked on the end of the string
      size_t p = input.find_first_of("\r\n");
      if ( p == string::npos ) {
	strings.push_back(input);
      } else {
	strings.push_back(input.substr(0,p));
      }
      getline(File, input);
    }
  }

  return strings;

}

// handy function to divide a string into substrings delimited using delim

vector<string> SPtokenize(const string & str, const string & delim)
{
  vector<string> tokens;

  size_t p0 = str.find_first_not_of(delim);
  size_t p1 = string::npos;
  while(p0 != string::npos)
  {
    p1 = str.find_first_of(delim, p0);
    if(p1 != p0)
    {
      string token = str.substr(p0, p1 - p0);
      tokens.push_back(token);
    }
    p0 = str.find_first_not_of(delim, p1);
  }

  return tokens;
}

// handy function to do a partial match from a vector of strings

string SPmatchString(const string& str, const vector<string>& strArray, int& nmatch)
{
  nmatch = 0;
  string outstr;

  for (size_t i=0; i<strArray.size(); i++) {
    if (str.size() <= strArray[i].size()) {
      if ( str == strArray[i].substr(0,str.size()) ) {
	if ( outstr.size() == 0 ) outstr = strArray[i];
	nmatch++;
      }
    }
  }

  return outstr;
}

// handy function to convert a vector of strings into a vector or Reals

bool SPstring2Real(const vector<string>& str, vector<Real>& value)
{
  value.resize(str.size());
  for (size_t i=0; i<str.size(); i++) {
    istringstream iss(str[i]);
    if ( (iss >> value[i]).fail() ) return false;
  }
  return true;
}

// handy function to convert a string into a Real

bool SPstring2Real(const string& str, Real& value)
{
  istringstream iss(str);
  return !(iss >> value).fail();
}

// handy function to convert a vector of strings into a vector of doubles

bool SPstring2double(const vector<string>& str, vector<double>& value)
{
  value.resize(str.size());
  for (size_t i=0; i<str.size(); i++) {
    istringstream iss(str[i]);
    if ( (iss >> value[i]).fail() ) return false;
  }
  return true;
}

// handy function to convert a string into a double

bool SPstring2double(const string& str, double& value)
{
  istringstream iss(str);
  return !(iss >> value).fail();
}

// handy function to convert a vector of strings into a vector of floats

bool SPstring2float(const vector<string>& str, vector<float>& value)
{
  value.resize(str.size());
  for (size_t i=0; i<str.size(); i++) {
    istringstream iss(str[i]);
    if ( (iss >> value[i]).fail() ) return false;
  }
  return true;
}

// handy function to convert a string into a float

bool SPstring2float(const string& str, double& value)
{
  istringstream iss(str);
  return !(iss >> value).fail();
}

// handy function to convert a vector of strings into a vector of Integers

bool SPstring2Integer(const vector<string>& str, vector<Integer>& value)
{
  value.resize(str.size());
  for (size_t i=0; i<str.size(); i++) {
    istringstream iss(str[i]);
    if ( (iss >> value[i]).fail() ) return false;
  }
  return true;
}

// handy function to convert a string into an Integer

bool SPstring2Integer(const string& str, Integer& value)
{
  istringstream iss(str);
  return !(iss >> value).fail();
}

// convert a string of delimited range specifications of form n1-n2 meaning n1 to
// n2 or n3 meaning just n3.

bool SPrangeString2IntegerList(const string& str, const string& delim1, 
			       const string& delim2, vector<Integer>& list)
{
  vector<string> strElements;
  strElements = SPtokenize(str, delim1);

  for (size_t i=0; i<strElements.size(); i++) {
    size_t p = strElements[i].find_first_of(delim2);
    if ( p == string::npos ) {
      Integer value;
      bool status = SPstring2Integer(strElements[i], value);
      if (!status ) return(status);
      list.push_back(value);
    } else {
      Integer value1, value2;
      bool status = SPstring2Integer(strElements[i].substr(0,p), value1);
      if (!status ) return(status);
      status = SPstring2Integer(strElements[i].substr(p+1), value2);
      if (!status ) return(status);
      for (Integer j=value1; j<=value2; j++) list.push_back(j);
    }
  }
  return(true);

}


// calculate factors for shifting an array

void SPcalcShift(const vector<Real>& Low, const vector<Real>& High, 
		 const vector<Integer>& vStart, const vector<Integer>& vEnd, 
		 const vector<Real>& vShift, const vector<Real>& vFactor, 
		 vector<vector<size_t> >& fromIndex, vector<vector<Real> >& Fraction)
{
  // Note that Start and End are zero-based and are used to determine which
  // element numbers of the Low and High arrays are modified.

  size_t Nelts(Low.size());
  fromIndex.resize(Nelts);
  Fraction.resize(Nelts);

  // toIndex[i] is the list of elements to which i contributes
  // LowShift[i] and HighShift[i] are the shifted ranges for row i.
  // shifted[i] keeps track of whether a given channel has been shifted

  vector<vector<size_t> > toIndex(Nelts);
  vector<Real> LowShift(Nelts);
  vector<Real> HighShift(Nelts);
  vector<bool> shifted(Nelts,false);

  // Loop over entries in the input vectors of shifts

  for (size_t iShift=0; iShift<vStart.size(); iShift++) {

    Integer Start = vStart[iShift];
    if ( Start < 0 ) Start = 0;
    Integer End = vEnd[iShift];
    if ( End > Nelts - 1 ) End = Nelts - 1;
    Real Shift = vShift[iShift];
    Real Factor = vFactor[iShift];

    // Check whether we actually need to do anything for this shift

    if ( End >= Low[0] && Start <= High[Nelts-1] && !(Shift==0.0 && Factor==1.0) ) {

      // Loop over elements for this shift

      size_t j1 = 0;

      for (size_t iElt=(size_t)Start; iElt<=(size_t)End; iElt++) {

	LowShift[iElt] = Low[iElt]*Factor + Shift;
	HighShift[iElt] = High[iElt]*Factor + Shift;
	shifted[iElt] = true;

	// Trap out the case of no overlap between the shifted range and the
	// input ranges

	if ( HighShift[iElt] >= Low[0] && LowShift[iElt] <= High[Nelts-1] ) {

	  // Find the first j such that Low[iElt]*Factor + Shift >= Low[j]
	  // and call this j1

	  if ( LowShift[iElt] >= Low[0] ) {
	    while (j1<Nelts && LowShift[iElt] > Low[j1]) j1++;
	    j1--;
	  }

	  // Find the last j such that High[iElt]*Factor + Shift < High[j]
	  // and call this j2

	  size_t j2 = j1;
	  if ( HighShift[iElt] < High[Nelts-1] ) {
	    while (j2<Nelts && HighShift[iElt] > High[j2]) j2++;
	  } else {
	    j2 = Nelts - 1;
	  }
	  
	  // Loop over j1 to j2 setting the toIndex array

	  for (size_t j=j1; j<=j2; j++) toIndex[iElt].push_back(j);

	}

	// end loop over elements for this shift

      }

    }

    // end loop over shifts

  }

  // Now go through elements one more time to check for any which were not included
  // in a shift and specify that they are copied over unchanged

  for (size_t iElts=0; iElts<Nelts; iElts++) {
    if ( !shifted[iElts] ) {
      toIndex[iElts].push_back(iElts);
      LowShift[iElts] = Low[iElts];
      HighShift[iElts] = High[iElts];
    }
  }

  // Now reshuffle so that we have a list of elements which contribute to a target element
  // Loop over elements

  for (size_t iElt=0; iElt<Nelts; iElt++) {
  
    // search through the toIndex array looking for instances of this element
    // if found put the information in fromIndex. Calculate Fraction as the
    // fraction of the target element which comes from the donor element

    Real Diff = High[iElt] - Low[iElt];

    for (size_t i=0; i<Nelts; i++) {
      for (size_t j=0; j<toIndex[i].size(); j++) {
	if ( toIndex[i][j] == iElt ) {
	  fromIndex[iElt].push_back(i);
	  Real LowEnd = Low[iElt];
	  if ( LowShift[i] > LowEnd ) LowEnd = LowShift[i];
	  Real HighEnd = High[iElt];
	  if ( HighShift[i] < HighEnd ) HighEnd = HighShift[i];
	  Fraction[iElt].push_back((HighEnd-LowEnd)/Diff);
	}
      }
    }

  }

  return;
}

// if array is increasing the index is the last element in array <= target; if 
// target < array[0] then index = -1; if target >= array[N-1] then index = N-1.
// if array is decreasing the index is the last element in array >= target; if
// target > array[0] then index = -1; if target <= array[N-1] then index = N-1.
// template class T can be either vector or valarray of Real or Integer

template <class T> void SPfind(const T& array, const Real& target, Integer& index)
{
  size_t N (array.size()); 
  if (N < 2) {
    index = 0;
    return;
  }

  bool increasing(array[0] < array[N-1]);        
  if ( increasing ) {
    if (target < array[0]) {
      index = -1;
    } else if (target >= array[N-1]) {
      index = N-1;
    } else {
      index = 0;
      Integer upper = N-1;
      SPbisect(index, upper, array, target, increasing);
    }           
  } else {
    if (target > array[0]) {
      index = -1;
    } else if (target <= array[N-1]) {
      index = N-1;
    } else {
      index = 0;
      Integer upper = N-1;
      SPbisect(index, upper, array, target, increasing);
    }
  }

  return;
}

// required to make the linker instantiate correctly - options for the template

template void SPfind(const vector<Real>&, const Real&, Integer&);
template void SPfind(const RealArray&, const Real&, Integer&);
template void SPfind(const IntegerArray&, const Real&, Integer&);
template void SPfind(const valarray<Integer>&, const Real&, Integer&);

// handy routine to do bisection search between indices lower and upper

template <class T> void SPbisect(Integer& lower, Integer& upper, const T& array,
				 const Real& target, bool increasing)
{
  while ( upper - lower > 1) {
    Integer mid   = ( lower + upper )/ 2;
    if (increasing) {
      if (target < array[mid] ) upper = mid;
      else lower = mid;
    } else {
      if ( target > array[mid]) upper = mid;
      else lower = mid;      
    }
  }
  return;
}

// required to make the linker instantiate correctly - options for the template

template void SPbisect(Integer&, Integer&, const vector<Real>&, const Real&, bool);
template void SPbisect(Integer&, Integer&, const RealArray&, const Real&, bool);
template void SPbisect(Integer&, Integer&, const IntegerArray&, const Real&, bool);
template void SPbisect(Integer&, Integer&, const valarray<Integer>&, const Real&, bool);
