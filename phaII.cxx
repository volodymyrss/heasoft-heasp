// phaII object code. Definitions in SpectrumII.h

#ifndef HAVE_phaII
#include "phaII.h"
#endif

#ifndef HAVE_SPio
#include "SPio.h"
#endif


// Class phaII

// default constructor

phaII::phaII()
{
}


// Destructor

phaII::~phaII()
{
}

// read all spectra from the PHA file

Integer phaII::read(string filename)
{
  return(this->read(filename, 1));
}

// read all spectra from the PHA file extension with EXTVER=PHAnumber

Integer phaII::read(string filename, Integer PHAnumber)
{
  Integer Status(OK);

  // find the number of spectra

  Integer nSpectra = NumberofSpectra(filename, PHAnumber, Status);
  if ( Status != OK ) return(Status);

  // set up the SpectrumNumber array (which is 1-based)

  vector<Integer> SpectrumNumber(nSpectra);
  for (size_t i=0; i<(size_t)nSpectra; i++) SpectrumNumber[i] = i + 1;

  // read the spectra

  return(this->read(filename, PHAnumber, SpectrumNumber));
}

// reading from PHA file.Somewhat inefficient but just reads each spectrum individually

Integer phaII::read(string filename, Integer PHAnumber, vector<Integer> SpectrumNumber)
{
  Integer Status(OK);

  for (size_t i=0; i<SpectrumNumber.size(); i++) {
    pha inpha;
    Status = inpha.read(filename, PHAnumber, SpectrumNumber[i]);
    if ( Status != OK ) return(Status);
    phas.push_back(inpha);
  }

  return(Status);

}

// Deep copy

phaII& phaII::operator=(const phaII& beta)
{
  // Copy each constituent PHA object

  for (size_t i=0; i<beta.phas.size(); i++) {
    pha spectrum;
    spectrum = beta.phas[i];
    phas.push_back(spectrum);
  }

  return *this;

}

// Get pha object (counts from zero).

pha phaII::get(Integer number)
{
  pha spectrum;
  if ( number >= 0 && number < (Integer)phas.size() ) spectrum = phas[number];
  return spectrum;
}

// Push pha object into phaII object

void phaII::push(pha spectrum)
{
  phas.push_back(spectrum);
  return;
}

// Return the number of spectra in the object

Integer phaII::NumberSpectra()
{
  return phas.size();
}

// Display information about the spectrum. Just runs Spectrum::disp on each spectrum

string phaII::disp()
{
  ostringstream outstr;

  for (size_t i=0; i<phas.size(); i++) {
    outstr << "Spectrum # " << i+1 << endl;
    outstr << phas[i].disp() << endl;
  }

  return outstr.str();

}

// Clear information about the spectra

void phaII::clear()
{
  for (size_t i=0; i<phas.size(); i++) phas[i].clear();
  return;
}

// Check completeness and consistency of information in spectra
  // if there is a problem then return diagnostic in string

string phaII::check()
{
  ostringstream outstr;
 
  for (size_t i=0; i<phas.size(); i++) outstr << phas[i].check() << endl;

  // check for consistency of Poisserr and Datatype between spectra.

  bool poiss(phas[0].Poisserr);
  string type(phas[0].Datatype);
  for (size_t i=1; i<phas.size(); i++) {
    if ( poiss != phas[i].Poisserr ) {
      outstr << "Poisserr is inconsistent between spectrum 0 and spectrum " << i << endl;
    }
    if ( type != phas[i].Datatype ) {
      outstr << "Datatype is inconsistent between spectrum 0 and spectrum " << i << endl;
    }
  }

  return outstr.str();
}

// Write spectra as type II file

Integer phaII::write(string filename)
{

  string Blank = " ";
  
  vector<string> ttype;
  vector<string> tform;
  vector<string> tunit;

  // Set the Poisserr, Datatype and PHAVersion for the output. Assume for the 
  // moment that these are the same for all input spectra

  bool Poisserr = phas[0].Poisserr;
  string Datatype = phas[0].Datatype;
  string PHAVersion = phas[0].PHAVersion;

  // create a new FITS file
  
  std::auto_ptr<FITS> pFits(0);
      
  try {                
    pFits.reset( new FITS(filename,Write) );
  } catch (FITS::CantCreate) {
    string msg = "Failed to create "+filename+" for SPECTRUM extension";
    SPreportError(CannotCreate, msg);
    return(CannotCreate);
  }

  // Hopefully all the spectra we are writing have the same number of channels 
  // but just in case find the maximum

  size_t Nspec = phas.size();
  size_t MaxElements = phas[0].Pha.size();
  for (size_t i=1; i<Nspec; i++) {
    if (phas[i].Pha.size() > MaxElements) MaxElements = phas[i].Pha.size();
  }

  stringstream RepeatStream;
  RepeatStream << MaxElements;
  string Repeat(RepeatStream.str());

  // set up the column descriptors for those attributes which need to be 
  // written as columns

  bool isvector;
    
  ttype.push_back("SPEC_NUM");
  tform.push_back("J");
  tunit.push_back(" ");

  vector<vector<Integer> > Channel;
  for (size_t i=0; i<Nspec; i++) {
    Channel.push_back(phas[i].Channel);
  }
  if ( SPneedCol(Channel, isvector) ) {
    ttype.push_back("CHANNEL");
    tunit.push_back(" ");
    if ( isvector ) {
      tform.push_back(Repeat+"J");
    } else {
      tform.push_back("J");
    }
  }

  // Note that we assume all spectra have the same Datatype

  vector<vector<Real> > PhaRate;
  vector<vector<Integer> > PhaCount;

  if ( Datatype == "RATE" ) {
    ttype.push_back("RATE");
    tform.push_back(Repeat+"E");
    tunit.push_back("counts/s");
    for (size_t i=0; i<Nspec; i++) {
      PhaRate.push_back(phas[i].Pha);
    }
  } else {
    ttype.push_back("COUNTS");
    tform.push_back(Repeat+"J");
    tunit.push_back("counts");
    for (size_t i=0; i<Nspec; i++) {
      vector<Integer> Counts(phas[i].Pha.size());
      for (size_t j=0; j<phas[i].Pha.size(); j++) {
      	Counts[j] = (Integer)phas[i].Pha[j];
      }
      PhaCount.push_back(Counts);
    }
  }

  // Find which quantities need to be in columns and which can stay in keywords.
  // To be in a keyword the quantity must be the same for all spectra and all 
  // channels. We also need to find which columns can be scalar and which need 
  // to be vector. This is ugly.so wrap up test in a utility routine.

  // Note that we assume all spectra have the same POISSERR.

  vector<vector<Real> > StatError;
  if ( !Poisserr ) {
    for (size_t i=0; i<Nspec; i++) {
      StatError.push_back(phas[i].StatError);
    }
    if ( SPneedCol(StatError, isvector) ) {
      ttype.push_back("STAT_ERR");
      tunit.push_back(" ");
      if ( isvector ) {
	tform.push_back(Repeat+"E");
      } else {
	tform.push_back("E");
      }
    }
  }
	
  vector<vector<Real> > SysError;
  for (size_t i=0; i<Nspec; i++) {
    SysError.push_back(phas[i].SysError);
  }
  if ( SPneedCol(SysError, isvector) ) {
    ttype.push_back("SYS_ERR");
    tunit.push_back(" ");
    if ( isvector ) {
      tform.push_back(Repeat+"E");
    } else {
      tform.push_back("E");
    }
  }

  vector<vector<Integer> > Quality;
  for (size_t i=0; i<Nspec; i++) {
    Quality.push_back(phas[i].Quality);
  }
  if ( SPneedCol(Quality, isvector) ) {
    ttype.push_back("QUALITY");
    tunit.push_back(" ");
    if ( isvector ) {
      tform.push_back(Repeat+"I");
    } else {
      tform.push_back("I");
    }
  }

  vector<vector<Integer> > Group;
  for (size_t i=0; i<Nspec; i++) {
    Group.push_back(phas[i].Group);
  }
  if ( SPneedCol(Group, isvector) ) {
    ttype.push_back("GROUPING");
    tunit.push_back(" ");
    if ( isvector ) {
      tform.push_back(Repeat+"I");
    } else {
      tform.push_back("I");
    }
  }

  vector<vector<Real> > AreaScaling;
  for (size_t i=0; i<Nspec; i++) {
    AreaScaling.push_back(phas[i].AreaScaling);
  }
  if ( SPneedCol(AreaScaling, isvector) ) {
    ttype.push_back("AREASCAL");
    tunit.push_back(" ");
    if ( isvector ) {
      tform.push_back(Repeat+"E");
    } else {
      tform.push_back("E");
    }
  }

  vector<vector<Real> > BackScaling;
  for (size_t i=0; i<Nspec; i++) {
    BackScaling.push_back(phas[i].BackScaling);
  }
  if ( SPneedCol(BackScaling, isvector) ) {
    ttype.push_back("BACKSCAL");
    tunit.push_back(" ");
    if ( isvector ) {
      tform.push_back(Repeat+"E");
    } else {
      tform.push_back("E");
    }
  }

  // now we need to find if any of the standard keywords differ between the 
  // spectra so must be placed in a column. Quantities to check are Exposure, 
  // CorrectionScaling, SpectrumType, ResponseFile, AncillaryFile, 
  //  BackgroundFile, CorrectionFile, ChannelType, Telescope, Instrument, 
  // Detector, Filter, Datamode.

  vector<Real> Exposure;
  for (size_t i=0; i<Nspec; i++) {
    Exposure.push_back(phas[i].Exposure);
  }
  if ( SPneedCol(Exposure) ) {
    ttype.push_back("EXPOSURE");
    tform.push_back("E");
    tunit.push_back(" ");
  }

  vector<Real> CorrectionScaling;
  for (size_t i=0; i<Nspec; i++) {
    CorrectionScaling.push_back(phas[i].CorrectionScaling);
  }
  if ( SPneedCol(CorrectionScaling) ) {
    ttype.push_back("CORRSCAL");
    tform.push_back("E");
    tunit.push_back(" ");
  }

  vector<Integer> DetChans;
  for (size_t i=0; i<Nspec; i++) {
    DetChans.push_back(phas[i].DetChans);
  }
  if ( SPneedCol(DetChans) ) {
    ttype.push_back("DETCHANS");
    tform.push_back("J");
    tunit.push_back(" ");
  }

  // now check all the string variables

  vector<string> Spectrumtype;
  for (size_t i=0; i<Nspec; i++) {
    Spectrumtype.push_back(phas[i].Spectrumtype);
  }
  if ( SPneedCol(Spectrumtype) ) {
    ttype.push_back("HDUCLAS2");
    tform.push_back(SPstringTform(Spectrumtype));
    tunit.push_back(" ");
  }
  
  vector<string> ResponseFile;
  for (size_t i=0; i<Nspec; i++) {
    ResponseFile.push_back(phas[i].ResponseFile);
  }
  if ( SPneedCol(ResponseFile) ) {
    ttype.push_back("RESPFILE");
    tform.push_back(SPstringTform(ResponseFile));
    tunit.push_back(" ");
  }
  
  vector<string> AncillaryFile;
  for (size_t i=0; i<Nspec; i++) {
    AncillaryFile.push_back(phas[i].AncillaryFile);
  }
  if ( SPneedCol(AncillaryFile) ) {
    ttype.push_back("ANCRFILE");
    tform.push_back(SPstringTform(AncillaryFile));
    tunit.push_back(" ");
  }
  
  vector<string> BackgroundFile;
  for (size_t i=0; i<Nspec; i++) {
    BackgroundFile.push_back(phas[i].BackgroundFile);
  }
  if ( SPneedCol(BackgroundFile) ) {
    ttype.push_back("BACKFILE");
    tform.push_back(SPstringTform(BackgroundFile));
    tunit.push_back(" ");
  }
  
  vector<string> CorrectionFile;
  for (size_t i=0; i<Nspec; i++) {
    CorrectionFile.push_back(phas[i].CorrectionFile);
  }
  if ( SPneedCol(CorrectionFile) ) {
    ttype.push_back("CORRFILE");
    tform.push_back(SPstringTform(CorrectionFile));
    tunit.push_back(" ");
  }
  
  vector<string> ChannelType;
  for (size_t i=0; i<Nspec; i++) {
    ChannelType.push_back(phas[i].ChannelType);
  }
  if ( SPneedCol(ChannelType) ) {
    ttype.push_back("CHANTYPE");
    tform.push_back(SPstringTform(ChannelType));
    tunit.push_back(" ");
  }
  
  vector<string> Telescope;
  for (size_t i=0; i<Nspec; i++) {
    Telescope.push_back(phas[i].Telescope);
  }
  if ( SPneedCol(Telescope) ) {
    ttype.push_back("TELESCOP");
    tform.push_back(SPstringTform(Telescope));
    tunit.push_back(" ");
  }
  
  vector<string> Instrument;
  for (size_t i=0; i<Nspec; i++) {
    Instrument.push_back(phas[i].Instrument);
  }
  if ( SPneedCol(Instrument) ) {
    ttype.push_back("INSTRUME");
    tform.push_back(SPstringTform(Instrument));
    tunit.push_back(" ");
  }
  
  vector<string> Detector;
  for (size_t i=0; i<Nspec; i++) {
    Detector.push_back(phas[i].Detector);
  }
  if ( SPneedCol(Detector) ) {
    ttype.push_back("DETNAM");
    tform.push_back(SPstringTform(Detector));
    tunit.push_back(" ");
  }

  vector<string> Filter;
  for (size_t i=0; i<Nspec; i++) {
    Filter.push_back(phas[i].Filter);
  }
  if ( SPneedCol(Filter) ) {
    ttype.push_back("FILTER");
    tform.push_back(SPstringTform(Filter));
    tunit.push_back(" ");
  }
  
  vector<string> Datamode;
  for (size_t i=0; i<Nspec; i++) {
    Datamode.push_back(phas[i].Datamode);
  }
  if ( SPneedCol(Datamode) ) {
    ttype.push_back("DATAMODE");
    tform.push_back(SPstringTform(Datamode));
    tunit.push_back(" ");
  }
  
  // Do XSPECFilter here. Assume that all spectra have the same number of
  // XSPEC filter keywords

  for (size_t j=0; j<phas[0].XSPECFilter.size(); j++) {
    vector<string> XSPECFilter;
    for (size_t i=0; i<Nspec; i++) {
      XSPECFilter.push_back(phas[i].XSPECFilter[j]);
    }
    if ( SPneedCol(XSPECFilter) ) {
      ostringstream ColNameStream;
      ColNameStream << "XFLT" << setfill('0') << setw(4) << j+1;
      string ColName(ColNameStream.str());
      ttype.push_back(ColName);
      tform.push_back(SPstringTform(XSPECFilter));
      tunit.push_back(" ");
    }
  }

  // Create the new extension. First check for existing extensions to see whether
  // we need to give a version number.

  ExtMapConstIt itLow = pFits->extension().lower_bound("SPECTRUM");
  ExtMapConstIt itHigh = pFits->extension().upper_bound("SPECTRUM");

  Integer version(0);
  while (itLow != itHigh) {
    if (itLow->second->version() > version) version = itLow->second->version();
    itLow++;
  }

  version++;
  Table* pspectrum = pFits->addTable("SPECTRUM",phas.size(),ttype,tform,tunit,BinaryTbl,version);
  Table& spectrum = *pspectrum;

  // write the standard keywords that should always be the same for all spectra

  SPwriteKey(spectrum, "HDUCLASS", (string)"OGIP", " ");
    
  SPwriteKey(spectrum, "HDUCLAS1", (string)"SPECTRUM", " ");

  SPwriteKey(spectrum, "HDUCLAS3", Datatype, " ");
    
  SPwriteKey(spectrum, "HDUVERS", PHAVersion, "OGIP version number");

  // Write the columns. Column data will be written either as keyword, scalar 
  // or vector column in that order of priority.

  vector<Integer> SpecNum(Nspec);
  for (size_t i=0; i<Nspec; i++) SpecNum[i] = i + 1;
  SPwriteCol(spectrum, "SPEC_NUM", SpecNum);

  SPwriteVectorCol(spectrum, "CHANNEL", Channel);

  if ( Datatype == "RATE" ) {
    SPwriteVectorCol(spectrum, "RATE", PhaRate);
  } else {
    SPwriteVectorCol(spectrum, "COUNTS", PhaCount);
  }

  if ( !Poisserr ) {
    SPwriteVectorCol(spectrum, "STAT_ERR", StatError);
  }

  SPwriteVectorCol(spectrum, "SYS_ERR", SysError);

  SPwriteVectorCol(spectrum, "QUALITY", Quality);

  SPwriteVectorCol(spectrum, "GROUPING", Group);

  SPwriteVectorCol(spectrum, "AREASCAL", AreaScaling);

  SPwriteVectorCol(spectrum, "BACKSCAL", BackScaling);

  SPwriteCol(spectrum, "EXPOSURE", Exposure);

  SPwriteCol(spectrum, "CORRSCAL", CorrectionScaling);

  SPwriteCol(spectrum, "DETCHANS", DetChans);

  SPwriteCol(spectrum, "HDUCLAS2", Spectrumtype);

  SPwriteCol(spectrum, "RESPFILE", ResponseFile);

  SPwriteCol(spectrum, "ANCRFILE", AncillaryFile);

  SPwriteCol(spectrum, "BACKFILE", BackgroundFile);

  SPwriteCol(spectrum, "CORRFILE", CorrectionFile);

  SPwriteCol(spectrum, "CHANTYPE", ChannelType);

  SPwriteCol(spectrum, "TELESCOP", Telescope);

  SPwriteCol(spectrum, "INSTRUME", Instrument);

  SPwriteCol(spectrum, "DETNAM", Detector); 

  SPwriteCol(spectrum, "FILTER", Filter);

  SPwriteCol(spectrum, "DATAMODE", Datamode);

  for (size_t j=0; j<phas[0].XSPECFilter.size(); j++) {
    vector<string> XSPECFilter;
    for (size_t i=0; i<Nspec; i++) {
      XSPECFilter.push_back(phas[i].XSPECFilter[j]);
    }
    ostringstream ColNameStream;
    ColNameStream << "XFLT" << setfill('0') << setw(4) << j+1;
    string ColName(ColNameStream.str());
    SPwriteCol(spectrum, ColName, XSPECFilter);
  }

  return(OK);

}

// Write spectra as type II file copying extra keywords and extensions from
// another file

Integer phaII::write(string filename, string copyfilename)
{
  return(this->write(filename, copyfilename, 1));
}

// Write spectrum as type II file copying extra keywords and extensions from
// another file. The required spectrum is the HDUnumber instance of a spectrum
// extension in the file.

Integer phaII::write(string filename, string copyfilename, Integer HDUnumber)
{
  Integer Status(OK);

  Status = this->write(filename);

  if ( Status != OK ) return(Status);

  Status = SPcopyKeys(copyfilename, filename, "SPECTRUM", HDUnumber);

  if ( Status != OK ) return(Status);

  Status = SPcopyHDUs(copyfilename, filename);

  return(Status);
}
