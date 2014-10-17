// arfII object code. Definitions in arfII.h

#ifndef HAVE_arfII
#include "arfII.h"
#endif

#ifndef HAVE_SPio
#include "SPio.h"
#endif


// Class arfII

// default constructor

arfII::arfII()
{
}


// Destructor

arfII::~arfII()
{
}

// read all ARFs from the file

Integer arfII::read(string filename)
{
  return(this->read(filename, 1));
}

// read all ARFs from the SPECRESP extension with EXTVER=ARFnumber

Integer arfII::read(string filename, Integer ARFnumber)
{
  Integer Status(OK);

  // find the number of ARFs

  Integer nARFs = NumberofARFs(filename, ARFnumber, Status);
  if ( Status != OK ) return(Status);

  // set up the RowNumber array (which is 1-based)

  vector<Integer> RowNumber(nARFs);
  for (size_t i=0; i<(size_t)nARFs; i++) RowNumber[i] = i + 1;

  // read the ARFs

  return(this->read(filename, ARFnumber, RowNumber));
}

// reading from ARF file.Somewhat inefficient but just reads each ARF individually

Integer arfII::read(string filename, Integer ARFnumber, vector<Integer> RowNumber)
{
  Integer Status(OK);

  for (size_t i=0; i<RowNumber.size(); i++) {
    arf inarf;
    Status = inarf.read(filename, ARFnumber, RowNumber[i]);
    if ( Status != OK ) return(Status);
    arfs.push_back(inarf);
  }

  return(OK);

}

// Deep copy

arfII& arfII::operator=(const arfII& beta)
{
  // Copy each constituent ARF object

  for (size_t i=0; i<beta.arfs.size(); i++) {
    arf onearf;
    onearf = beta.arfs[i];
    arfs.push_back(onearf);
  }

  return *this;

}

// Get arf object (counts from zero).

arf arfII::get(Integer number)
{
  arf ea;
  if ( number >= 0 && number < (Integer)arfs.size() ) ea = arfs[number];
  return ea;
}

// Push arf object into arfII object

void arfII::push(arf ea)
{
  arfs.push_back(ea);
  return;
}

// Return the number of ARS in the object

Integer arfII::NumberARFs()
{
  return arfs.size();
}

// Display information about the spectrum. Just runs arf::disp on each arf
// return as a string

string arfII::disp()
{
  ostringstream outstr;

  for (size_t i=0; i<arfs.size(); i++) {
    outstr << "ARF # " << i+1 << endl;
    outstr << arfs[i].disp() << endl;
  }

  return outstr.str();

}

// Clear information from the ARFs

void arfII::clear()
{
  for (size_t i=0; i<arfs.size(); i++) arfs[i].clear();
  return;
}

// Check completeness and consistency of information in arfs
  // if there is a problem then return diagnostic in string

string arfII::check()
{
  ostringstream outstr;

  for (size_t i=0; i<arfs.size(); i++) outstr << arfs[i].check() << endl;
  return outstr.str();
}

// Write arfs as type II file

Integer arfII::write(string filename)
{

  string Blank = " ";
  
  vector<string> ttype;
  vector<string> tform;
  vector<string> tunit;

  // Set the Version and ExtensionName for the output. Assume for the 
  // moment that these are the same for all input ARFs

  string ExtensionName = arfs[0].ExtensionName;
  string Version = arfs[0].Version;

  // create a new FITS file
  
  std::auto_ptr<FITS> pFits(0);
      
  try {                
    pFits.reset( new FITS(filename,Write) );
  } catch (FITS::CantCreate) {
    string msg = "Failed to create "+filename+" for SPECRESP extension";
    SPreportError(CannotCreate, msg);
    return(CannotCreate);       
  }

  // Hopefully all the ARFs we are writing have the same number of energy bins
  // but just in case find the maximum

  size_t Narfs = arfs.size();
  size_t MaxBins = arfs[0].EffArea.size();
  for (size_t i=1; i<Narfs; i++) {
    if (arfs[i].EffArea.size() > MaxBins) MaxBins = arfs[i].EffArea.size();
  }

  stringstream RepeatStream;
  RepeatStream << MaxBins;
  string Repeat(RepeatStream.str());

  // set up the column descriptors for those attributes which need to be 
  // written as columns

  ttype.push_back("ARF_NUM");
  tform.push_back("J");
  tunit.push_back(" ");

  vector<vector<Real> > LowEnergy;
  ttype.push_back("ENERG_LO");
  tform.push_back(Repeat+"E");
  tunit.push_back("counts/s");
  for (size_t i=0; i<Narfs; i++) {
    LowEnergy.push_back(arfs[i].LowEnergy);
  }

  vector<vector<Real> > HighEnergy;
  ttype.push_back("ENERG_HI");
  tform.push_back(Repeat+"E");
  tunit.push_back("counts/s");
  for (size_t i=0; i<Narfs; i++) {
    HighEnergy.push_back(arfs[i].HighEnergy);
  }

  vector<vector<Real> > EffArea;
  ttype.push_back("SPECRESP");
  tform.push_back(Repeat+"E");
  tunit.push_back("counts/s");
  for (size_t i=0; i<Narfs; i++) {
    EffArea.push_back(arfs[i].EffArea);
  }

  // now we need to find if any of the standard keywords differ between the 
  // ARFs so must be placed in a column. Quantities to check are 
  // Telescope, Instrument, Detector, Filter, 
  
  vector<string> Telescope;
  for (size_t i=0; i<Narfs; i++) {
    Telescope.push_back(arfs[i].Telescope);
  }
  if ( SPneedCol(Telescope) ) {
    ttype.push_back("TELESCOP");
    tform.push_back(SPstringTform(Telescope));
    tunit.push_back(" ");
  }
  
  vector<string> Instrument;
  for (size_t i=0; i<Narfs; i++) {
    Instrument.push_back(arfs[i].Instrument);
  }
  if ( SPneedCol(Instrument) ) {
    ttype.push_back("INSTRUME");
    tform.push_back(SPstringTform(Instrument));
    tunit.push_back(" ");
  }
  
  vector<string> Detector;
  for (size_t i=0; i<Narfs; i++) {
    Detector.push_back(arfs[i].Detector);
  }
  if ( SPneedCol(Detector) ) {
    ttype.push_back("DETNAM");
    tform.push_back(SPstringTform(Detector));
    tunit.push_back(" ");
  }

  vector<string> Filter;
  for (size_t i=0; i<Narfs; i++) {
    Filter.push_back(arfs[i].Filter);
  }
  if ( SPneedCol(Filter) ) {
    ttype.push_back("FILTER");
    tform.push_back(SPstringTform(Filter));
    tunit.push_back(" ");
  }
  
  // Create the new extension. First check for existing extensions to see whether
  // we need to give a version number.

  ExtMapConstIt itLow = pFits->extension().lower_bound(ExtensionName);
  ExtMapConstIt itHigh = pFits->extension().upper_bound(ExtensionName);

  Integer version(0);
  while (itLow != itHigh) {
    if (itLow->second->version() > version) version = itLow->second->version();
    itLow++;
  }

  version++;
  Table* parf = pFits->addTable(ExtensionName,arfs.size(),ttype,tform,tunit,BinaryTbl,version);
  Table& arf = *parf;

  // write the standard keywords that should always be the same for all ARFs

  SPwriteKey(arf, "HDUCLASS", (string)"OGIP", " ");
    
  SPwriteKey(arf, "HDUCLAS1", (string)"RESPONSE", " ");

  SPwriteKey(arf, "HDUCLAS2", (string)"SPECRESP", " ");

  SPwriteKey(arf, "HDUVERS", Version, "OGIP version number");

  // Write the columns. Column data will be written either as keyword, scalar 
  // or vector column in that order of priority.

  vector<Integer> ARFNum(Narfs);
  for (size_t i=0; i<Narfs; i++) ARFNum[i] = i + 1;
  SPwriteCol(arf, "ARF_NUM", ARFNum);

  SPwriteVectorCol(arf, "ENERG_LO", LowEnergy);

  SPwriteVectorCol(arf, "ENERG_HI", HighEnergy);

  SPwriteVectorCol(arf, "SPECRESP", EffArea);

  SPwriteCol(arf, "TELESCOP", Telescope);

  SPwriteCol(arf, "INSTRUME", Instrument);

  SPwriteCol(arf, "DETNAM", Detector); 

  SPwriteCol(arf, "FILTER", Filter);

  return(OK);

}

// Write arfs as type II file copying keywords and extensions from another file

Integer arfII::write(string filename, string copyfilename)
{
  return(this->write(filename, copyfilename, 1));
}

// Write arf as type II file copying keywords and extra extensions from another arf file. Use the HDUnumber instance of the arf extension in both files.

Integer arfII::write(string filename, string copyfilename, Integer HDUnumber)
{
  Integer Status(OK);

  Status = this->write(filename);

  if ( Status != OK ) return(Status);

  Status = SPcopyKeys(copyfilename, filename, "SPECRESP", HDUnumber);

  if ( Status != OK ) return(Status);

  Status = SPcopyHDUs(copyfilename, filename);

  return(Status);
}

