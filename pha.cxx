// Spectrum object code. Definitions in Spectrum.h

#ifndef HAVE_pha
#include "pha.h"
#endif

#ifndef HAVE_grouping
#include "grouping.h"
#endif

#ifndef HAVE_SPio
#include "SPio.h"
#endif


// Class pha

// default constructor

pha::pha()
  : FirstChannel(0),
    Pha(),
    StatError(),
    SysError(),
    Channel(),
    Quality(),
    Group(),
    AreaScaling(),
    BackScaling(),
    Exposure(0.0),
    CorrectionScaling(0.0),
    DetChans(0),
    Poisserr(false),
    Datatype("COUNT"),
    PHAVersion("1.2.1"),
    Spectrumtype("TOTAL"),
    ResponseFile("NONE"),
    AncillaryFile("NONE"),
    BackgroundFile("NONE"),
    CorrectionFile("NONE"),
    FluxUnits(" "),
    ChannelType("PI"),
    Telescope(" "),
    Instrument(" "),
    Detector(" "),
    Filter(" "),
    Datamode(" "),
    XSPECFilter()
{
}


// Destructor

pha::~pha()
{
  // clear vectors with guaranteed reallocation
  vector<Real>().swap(Pha);
  vector<Real>().swap(StatError);
  vector<Real>().swap(SysError);
  vector<Integer>().swap(Channel);
  vector<Integer>().swap(Quality);
  vector<Integer>().swap(Group);
  vector<Real>().swap(AreaScaling);
  vector<Real>().swap(BackScaling);
  vector<string>().swap(XSPECFilter);
}


// reading from pha file. 

Integer pha::read(string filename)
{
  return(this->read(filename, 1, 1));
}

// reading from PHA file. this option allows multiple extensions in the same file

Integer pha::read(string filename, Integer PHAnumber)
{
  return(this->read(filename, PHAnumber, 1));
}

// reading from PHA file. For a type I file SpectrumNumber should be set to 1

Integer pha::read(string filename, Integer PHAnumber, Integer SpectrumNumber)
{

  const string hduName("SPECTRUM");
  string DefString;
  bool verbosity = FITS::verboseMode();


  // Read in the SPECTRUM extension number PHAnumber
  // and set up an object called spectrum with the contents

  const vector<string> hduKeys;
  const vector<string> primaryKey;

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(filename,Read,hduName,false,hduKeys,primaryKey,(int)PHAnumber));
  } catch(...) {
    string msg = "Failed to read "+hduName+" in "+filename;
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  ExtHDU& spectrum = pInfile->extension(hduName, (int)PHAnumber);
  
  // read the standard keywords and store in the object

  DefString = "UNKNOWN";

  ChannelType = SPreadKey(spectrum, "CHANTYPE", SpectrumNumber, DefString);

  PHAVersion = SPreadKey(spectrum, "HDUVERS", SpectrumNumber, DefString);
  if ( PHAVersion == "UNKNOWN" ) {
    PHAVersion = SPreadKey(spectrum, "HDUVERS1", SpectrumNumber, DefString);
    if ( PHAVersion == "UNKNOWN" ) PHAVersion = "1.2.1";
  }

  Telescope = SPreadKey(spectrum, "TELESCOP", SpectrumNumber, DefString);
  
  Instrument = SPreadKey(spectrum, "INSTRUME", SpectrumNumber, DefString);

  Detector = SPreadKey(spectrum, "DETNAM", SpectrumNumber, DefString);

  Filter = SPreadKey(spectrum, "FILTER", SpectrumNumber, DefString);

  Datamode = SPreadKey(spectrum, "DATAMODE", SpectrumNumber, DefString);

  DefString = "TOTAL";
  Spectrumtype = SPreadKey(spectrum, "HDUCLAS2", SpectrumNumber, DefString);

  DefString = "NONE";
  ResponseFile = SPreadKey(spectrum, "RESPFILE", SpectrumNumber, DefString);
    
  AncillaryFile = SPreadKey(spectrum, "ANCRFILE", SpectrumNumber, DefString);

  BackgroundFile = SPreadKey(spectrum, "BACKFILE", SpectrumNumber, DefString);

  CorrectionFile = SPreadKey(spectrum, "CORRFILE", SpectrumNumber, DefString);

  CorrectionScaling = SPreadKey(spectrum, "CORRSCAL", SpectrumNumber, (Real)1.0);

  Exposure = SPreadKey(spectrum, "EXPOSURE", SpectrumNumber, (Real)0.0);

  Poisserr = SPreadKey(spectrum, "POISSERR", SpectrumNumber, false);

  // Read the XFLT keywords

  bool done = false;
  int i = 0;
  while ( (i++)<=9998 && !done ) {
    ostringstream KeyStream;
    KeyStream << "XFLT" << setfill('0') << setw(4) << i;
    string KeyName(KeyStream.str());
    string KeyValue;
    DefString = "NONE";
    KeyValue = SPreadKey(spectrum, KeyName, SpectrumNumber, DefString);
    if (KeyValue == "NONE") {
      Real KeyReal = SPreadKey(spectrum, KeyName, SpectrumNumber, (Real)-999.0);
      if (KeyReal != -999.0) {
	ostringstream RStream;
	RStream << KeyReal;
	KeyValue = RStream.str();
      }
    }
    if (KeyValue != "NONE") {
      XSPECFilter.push_back(KeyValue);
    } else {
      done = true;
    }
  }

  // Check for TLMIN set for the CHANNEL column

  FITS::setVerboseMode(false);
  try {
    int ChannelIndex = spectrum.column("CHANNEL").index();
    ostringstream KeyStream;
    KeyStream << "TLMIN" << ChannelIndex;
    spectrum.readKey(KeyStream.str(),FirstChannel);
    FirstChannel = SPreadKey(spectrum, KeyStream.str(), SpectrumNumber, (Integer)1);
  } catch(Table::NoSuchColumn&) {
    FirstChannel = 1;
  } catch(HDU::NoSuchKeyword&) {
    FirstChannel = 1;
  }
  FITS::clearErrors();
  FITS::setVerboseMode(verbosity);

  // Get the number of detector channels which may differ from the actual number
  // of rows

  DetChans = SPreadKey(spectrum,"DETCHANS",SpectrumNumber,(Integer)0);

  // Read the CHANNEL column

  SPreadCol(spectrum,"CHANNEL",SpectrumNumber,Channel);

  // Read the data and set the datatype column appropriately

  SPreadCol(spectrum,"COUNTS",SpectrumNumber,Pha);
  if ( Pha.size() != 0 ) {
    Datatype = "COUNT";
  } else {
    SPreadCol(spectrum,"RATE",SpectrumNumber,Pha);
    Datatype = "RATE";
  }

  if ( Pha.size() == 0 ) {
    string msg = "Input spectrum has no channels";
    SPreportError(NoData, msg);
    return(NoData);
  }

  // If no CHANNEL column was read and the Channel array can be constructed do so

  if ( Channel.size() == 0 ) {
    if ( Pha.size() == (size_t)DetChans ) {
      Channel.resize(Pha.size());
      for (size_t i=0; i<Pha.size(); i++) {
	Channel[i] = FirstChannel + i;
      }
    } else {
      string msg = "Cannot read or construct an array of Channel numbers";
      SPreportError(NoChannelData, msg);
      return(NoChannelData);
    }
  }

  // Read the statistical error if poisserr is false

  if (!Poisserr) {
    SPreadCol(spectrum,"STAT_ERR",SpectrumNumber,StatError);
    if ( StatError.size() == 0 ) {
      string msg = "POISSERR is false but there is no STAT_ERR column";
      SPreportError(NoStatError, msg);
      return(NoStatError);
    }
  }
  
  // Read the systematic error if it is there otherwise set the array to 0

  FITS::setVerboseMode(false);
  try {
    SPreadCol(spectrum,"SYS_ERR",SpectrumNumber,SysError);
  } catch (...) {
  }
  FITS::clearErrors();
  FITS::setVerboseMode(verbosity);

  if ( SysError.size() == 1 || SysError.size() == 0 ) {
    SysError.resize(Pha.size());
    for (size_t i=0; i<Pha.size(); i++) SysError[i] = 0.0;
  }

  // Read the QUALITY

  SPreadCol(spectrum,"QUALITY",SpectrumNumber,Quality);
  if ( Quality.size() == 1 || Quality.size() == 0 ) {
    Quality.resize(Pha.size());
    for (size_t i=0; i<Pha.size(); i++) Quality[i] = 0;
  }

  // Read the GROUPING
  
  SPreadCol(spectrum,"GROUPING",SpectrumNumber,Group);
  if ( Group.size() == 1 || Group.size() == 0) {
    Group.resize(Pha.size());
    for (size_t i=0; i<Pha.size(); i++) Group[i] = 1;
  }

  // Read the AREASCAL

  SPreadCol(spectrum,"AREASCAL",SpectrumNumber,AreaScaling);
  if ( AreaScaling.size() == 1 ) {
    AreaScaling.resize(Pha.size());
    for (size_t i=1; i<Pha.size(); i++) AreaScaling[i] = AreaScaling[0];
  }
  if ( AreaScaling.size() == 0 ) {
    AreaScaling.resize(Pha.size());
    for (size_t i=0; i<Pha.size(); i++) AreaScaling[i] = 1.0;
  }

  // Read the BACKSCAL

  SPreadCol(spectrum,"BACKSCAL",SpectrumNumber,BackScaling);
  if ( BackScaling.size() == 1 ) {
    BackScaling.resize(Pha.size());
    for (size_t i=1; i<Pha.size(); i++) BackScaling[i] = BackScaling[0];
  }
  if ( BackScaling.size() == 0 ) {
    BackScaling.resize(Pha.size());
    for (size_t i=0; i<Pha.size(); i++) BackScaling[i] = 1.0;
  }

  return(OK);

}

// Deep copy

pha& pha::operator=(const pha& beta)
{
  // Copy the scalars

  FirstChannel = beta.FirstChannel;
  Exposure = beta.Exposure;
  CorrectionScaling = beta.CorrectionScaling;
  DetChans = beta.DetChans;
  Poisserr = beta.Poisserr;
  Datatype = beta.Datatype;
  Spectrumtype = beta.Spectrumtype;
  ResponseFile = beta.ResponseFile;
  AncillaryFile = beta.AncillaryFile;
  BackgroundFile = beta.BackgroundFile;
  CorrectionFile = beta.CorrectionFile;
  FluxUnits = beta.FluxUnits;
  ChannelType = beta.ChannelType;
  PHAVersion = beta.PHAVersion;
  Telescope = beta.Telescope;  
  Instrument = beta.Instrument;
  Detector = beta.Detector;  
  Filter = beta.Filter;
  Datamode = beta.Datamode;

  // copy the XFLT keywords
  XSPECFilter.resize(beta.XSPECFilter.size());
  for (size_t i=0; i<beta.XSPECFilter.size(); i++) XSPECFilter[i] = beta.XSPECFilter[i];

  // now copy the arrays

  Channel.resize(beta.Channel.size());
  for (size_t i=0; i<beta.Channel.size(); i++) Channel[i] = beta.Channel[i];
  Pha.resize(beta.Pha.size());
  for (size_t i=0; i<beta.Pha.size(); i++) Pha[i] = beta.Pha[i];
  StatError.resize(beta.StatError.size());
  for (size_t i=0; i<beta.StatError.size(); i++) StatError[i] = beta.StatError[i];
  SysError.resize(beta.SysError.size());
  for (size_t i=0; i<beta.SysError.size(); i++) SysError[i] = beta.SysError[i];

  Quality.resize(beta.Quality.size());
  for (size_t i=0; i<beta.Quality.size(); i++) Quality[i] = beta.Quality[i];
  Group.resize(beta.Group.size());
  for (size_t i=0; i<beta.Group.size(); i++) Group[i] = beta.Group[i];

  AreaScaling.resize(beta.AreaScaling.size());
  for (size_t i=0; i<beta.AreaScaling.size(); i++) AreaScaling[i] = beta.AreaScaling[i];
  BackScaling.resize(beta.BackScaling.size());
  for (size_t i=0; i<beta.BackScaling.size(); i++) BackScaling[i] = beta.BackScaling[i];

  return *this;

}

// Return the number of channels in the arrays

Integer pha::NumberChannels()
{
  return Pha.size();
}

// Display information about the spectrum

string pha::disp()
{
  ostringstream outstr;

  outstr << "Spectrum information : \n";
  outstr << "   Number of channels   = " << NumberChannels()<< endl;
  outstr << "   Detchans             = " << DetChans << endl;
  outstr << "   Exposure             = " << Exposure<< endl;
  outstr << "   CorrectionScaling    = " << CorrectionScaling<< endl;
  if ( Poisserr ) {
    outstr << "   Poisserr             = " << "true"<< endl;
  } else {
    outstr << "   Poisserr             = " << "false"<< endl;
  }
  outstr << "   Datatype             = " << Datatype<< endl;
  outstr << "   Spectrumtype         = " << Spectrumtype<< endl;
  outstr << "   ResponseFile         = " << ResponseFile<< endl;
  outstr << "   AncillaryFile        = " << AncillaryFile<< endl;
  outstr << "   BackgroundFile       = " << BackgroundFile<< endl;
  outstr << "   CorrectionFile       = " << CorrectionFile<< endl;
  outstr << "   FluxUnits            = " << FluxUnits<< endl;
  outstr << "   ChannelType          = " << ChannelType<< endl;
  outstr << "   PHAVersion           = " << PHAVersion<< endl;
  outstr << "   Telescope            = " << Telescope<< endl;
  outstr << "   Instrument           = " << Instrument<< endl;  
  outstr << "   Detector             = " << Detector<< endl;  
  outstr << "   Filter               = " << Filter<< endl;  
  outstr << "   Datamode             = " << Datamode<< endl;
  outstr << endl;
  if ( AreaScaling.size() > 0 ) {
    outstr << "   AreaScaling[0]       = " << AreaScaling[0]<< endl;
  }
  if ( BackScaling.size() > 0 ) {
    outstr << "   BackScaling[0]       = " << BackScaling[0]<< endl;
  }
  if ( Quality.size() > 0 ) {
    outstr << "   Quality[0]           = " << Quality[0]<< endl;
  }
  if ( Group.size() > 0 ) {
    outstr << "   Group[0]          = " << Group[0]<< endl;
  }
  outstr << endl;

  if ( XSPECFilter.size() > 0 ) {
    outstr << "   XSPEC filter keywords : " << endl;
    for (size_t i=0; i<XSPECFilter.size(); i++) {
      outstr << "          " << XSPECFilter[i] << endl;
    }
  }

  if ( Channel.size() > 1 ) outstr << "  Channel array of size " << Channel.size() << endl;
  if ( Pha.size() > 1 ) outstr << "  Pha array of size " << Pha.size() << endl;
  if ( StatError.size() > 1 ) outstr << "  StatError array of size " << StatError.size() << endl;
  if ( SysError.size() > 1  ) outstr << "  SysError array of size " << SysError.size() << endl;
  if ( Quality.size() > 1   ) outstr << "  Quality array of size " << Quality.size() << endl;
  if ( Group.size() > 1  ) outstr << "  Grouping array of size " << Group.size() << endl;
  if ( AreaScaling.size() > 1 ) outstr << "  AreaScaling array of size " << AreaScaling.size() << endl;
  if ( BackScaling.size() > 1 ) outstr << "  BackScaling array of size " << BackScaling.size() << endl;
  outstr << endl;

  // for debugging
  //  for (size_t i=0; i<Pha.size(); i++) {
  //    outstr << i << "  ";
  //    if ( Channel.size() > 1 ) outstr << "  " << Channel[i];
  //    outstr << "  " << Pha[i];
  //    if ( StatError.size() > 1 ) outstr << "  " << StatError[i];
  //    if ( SysError.size() > 1  ) outstr << "  " << SysError[i];
  //    if ( Quality.size() > 1   ) outstr << "  " << Quality[i];
  //    if ( Group.size() > 1  ) outstr << "  " << Group[i];
  //    if ( AreaScaling.size() > 1 ) outstr << "  " << AreaScaling[i];
  //    if ( BackScaling.size() > 1 ) outstr << "  " << BackScaling[i];
  //    outstr << endl;
  //  }
  // end of debugging

  return outstr.str(); 

}

// Clear the information in the spectrum

void pha::clear()
{
  FirstChannel = 0;

  Pha.clear();
  StatError.clear();
  SysError.clear();
  Channel.clear();
  Quality.clear();
  Group.clear();
  AreaScaling.clear();
  BackScaling.clear();

  Exposure = 0.0;
  CorrectionScaling = 0.0;

  DetChans = 0;
  Poisserr = false;
  Datatype = " ";
  PHAVersion = "";
  Spectrumtype = " ";
  ResponseFile = " ";
  AncillaryFile = " ";
  BackgroundFile = " ";
  CorrectionFile = " ";
  FluxUnits = " ";
  ChannelType = " ";
  Telescope = " ";
  Instrument = " ";
  Detector = " ";
  Filter = " ";
  Datamode = " ";

  XSPECFilter.clear();

  return;

}

// Check completeness and consistency of information in spectrum

string pha::check()
{
  ostringstream outstr;

  // Check for an exposure time > 0

  if ( Exposure <= 0.0 ) {
    outstr << "Exposure time is negative (" << Exposure << ")" << endl;
  }

  // Check for data

  if ( Pha.size() == 0 ) {
    outstr << "No data present" << endl;
  }

  // Check for consistency between the size of the Pha array and DetChans

  if ( Pha.size() > (size_t)DetChans ) {
    outstr << "Size of Pha array (" << Pha.size() << ") is greater than DetChans ("
	 << DetChans << ")" << endl;
  }

  // Check for consistency between the size of the Pha and StatError arrays

  if ( Pha.size() != StatError.size() && StatError.size() != 0 && StatError.size() != 1 ) {
    outstr << "Size of Pha array (" << Pha.size() << ") differs from size of StatError array ("
	 << StatError.size() << ")" << endl;
  }

  // Check for consistency between the size of the Pha and SysError arrays

  if ( Pha.size() != SysError.size() && SysError.size() != 0 && SysError.size() != 1 ) {
    outstr << "Size of Pha array (" << Pha.size() << ") differs from size of SysError array ("
	 << SysError.size() << ")" << endl;
  }

  // Check for consistency between the size of the Pha and Channel arrays

  if ( Pha.size() != Channel.size() && Channel.size() != 0 ) {
    outstr << "Size of Pha array (" << Channel.size() << ") differs from size of Channel array ("
	 << Channel.size() << ")" << endl;
  }

  // Check for consistency between the size of the Pha and Quality arrays

  if ( Pha.size() != Quality.size() ) {
    outstr << "Size of Pha array (" << Quality.size() << ") differs from size of Quality array ("
	 << Quality.size() << ")" << endl;
  }

  // Check for consistency between the size of the Pha and Group arrays

  if ( Pha.size() != Group.size() ) {
    outstr << "Size of Pha array (" << Group.size() << ") differs from size of Group array ("
	 << Group.size() << ")" << endl;
  }

  // Check for consistency between the size of the Pha and AreaScaling arrays

  if ( Pha.size() != AreaScaling.size() ) {
    outstr << "Size of Pha array (" << AreaScaling.size() << ") differs from size of AreaScaling array ("
	 << AreaScaling.size() << ")" << endl;
  }

  // Check for consistency between the size of the Pha and BackScaling arrays

  if ( Pha.size() != BackScaling.size() ) {
    outstr << "Size of Pha array (" << BackScaling.size() << ") differs from size of BackScaling array ("
	 << BackScaling.size() << ")" << endl;
  }

  // Check for consistency of Poisserr and presence of StatError

  if ( Poisserr && StatError.size() > 0 ) {
    outstr << "Poisserr is true but StatError present" << endl;
  }
  if ( !Poisserr && StatError.size() == 0 ) {
    outstr << "Poisserr is false but no StatError present" << endl;
  }

  // Check Datatype is either COUNT or RATE

  if ( Datatype != "COUNT" && Datatype != "RATE" ) {
    outstr << "Datatype (" << Datatype << ") must be COUNT or RATE" << endl;
  }

  // Check Spectrumtype is one of TOTAL, NET or BKG

  if ( Spectrumtype != "TOTAL" && Spectrumtype != "NET" && Spectrumtype != "BKG" ) {
    outstr << "Spectrumtype (" << Spectrumtype << ") must be TOTAL, NET or BKG" << endl;
  }

  return outstr.str();
}


// Write spectrum as type I file

Integer pha::write(string filename)
{

  string Blank = " ";

  vector<string> ttype;
  vector<string> tform;
  vector<string> tunit;

  // Create a new FITS file instance  

  std::auto_ptr<FITS> pFits(0);
      
  try {                
    pFits.reset( new FITS(filename,Write) );
  } catch (FITS::CantCreate) {
    string msg = "Failed to create "+filename+" for SPECTRUM extension";
    SPreportError(CannotCreate, msg);
    return(CannotCreate);       
  }

  // set up the column descriptors for those attributes which need to be 
  // written as columns

  ttype.push_back("CHANNEL");
  tform.push_back("J");
  tunit.push_back(" ");

  if ( Datatype == "RATE" ) {
    ttype.push_back("RATE");
    tform.push_back("E");
    tunit.push_back("counts/s");
  } else {
    ttype.push_back("COUNTS");
    tform.push_back("J");
    tunit.push_back("counts");
  }
     
  if ( !Poisserr ) {
    ttype.push_back("STAT_ERR");
    tform.push_back("E");
    tunit.push_back(" ");
  }

  if ( SPneedCol(SysError) ) {
    ttype.push_back("SYS_ERR");
    tform.push_back("E");
    tunit.push_back(" ");
  }

  if ( SPneedCol(Quality) ) {
    ttype.push_back("QUALITY");
    tform.push_back("I");
    tunit.push_back(" ");
  }

  if ( SPneedCol(Group) ) {
    ttype.push_back("GROUPING");
    tform.push_back("I");
    tunit.push_back(" ");
  }

  if ( SPneedCol(AreaScaling) ) {
    ttype.push_back("AREASCAL");
    tform.push_back("E");
    tunit.push_back(" ");
  }

  if ( SPneedCol(BackScaling) ) {
    ttype.push_back("BACKSCAL");
    tform.push_back("E");
    tunit.push_back(" ");
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
  Table* pspectrum = pFits->addTable("SPECTRUM",Pha.size(),ttype,tform,tunit,BinaryTbl,version);
  Table& spectrum = *pspectrum;

  // Write the standard keywords
  
  SPwriteKey(spectrum, "HDUCLASS", (string)"OGIP", Blank);
    
  SPwriteKey(spectrum, "HDUCLAS1", (string)"SPECTRUM", Blank);

  SPwriteKey(spectrum, "HDUCLAS2", Spectrumtype, Blank);
    
  SPwriteKey(spectrum, "HDUCLAS3", Datatype, Blank);
    
  SPwriteKey(spectrum, "CHANTYPE", ChannelType, "Channel type");

  SPwriteKey(spectrum, "HDUVERS", PHAVersion, "OGIP version number");

  SPwriteKey(spectrum, "TELESCOP", Telescope, Blank);

  SPwriteKey(spectrum, "INSTRUME", Instrument, Blank);

  SPwriteKey(spectrum, "DETNAM", Detector, Blank);

  SPwriteKey(spectrum, "FILTER", Filter, Blank);

  SPwriteKey(spectrum, "DATAMODE", Datamode, Blank);

  SPwriteKey(spectrum, "DETCHANS", DetChans, "Number of channels in spectrum");

  SPwriteKey(spectrum, "TLMIN1", FirstChannel, "First channel number");

  SPwriteKey(spectrum, "EXPOSURE", Exposure, "Exposure time");

  SPwriteKey(spectrum, "CORRSCAL", CorrectionScaling, "Scaling for correction file");

  SPwriteKey(spectrum, "POISSERR", Poisserr, "Is error Poisson ?");

  SPwriteKey(spectrum, "RESPFILE", ResponseFile, Blank);

  SPwriteKey(spectrum, "ANCRFILE", AncillaryFile, Blank);

  SPwriteKey(spectrum, "BACKFILE", BackgroundFile, Blank);

  SPwriteKey(spectrum, "CORRFILE", CorrectionFile, Blank);

  for (size_t i=0; i<XSPECFilter.size(); i++) {
    ostringstream KeyStream;
    KeyStream << "XFLT" << setfill('0') << setw(4) << i+1;
    string KeyName(KeyStream.str());
    SPwriteKey(spectrum, KeyName, XSPECFilter[i], Blank);
  }

  // Write the arrays - if an array is of size 1 or all the same value 
  // it will be written as a keyword

  SPwriteCol(spectrum, "CHANNEL", Channel, true);

  if ( Datatype == "RATE") {
    SPwriteCol(spectrum, "RATE", Pha, true);
  } else {
    SPwriteCol(spectrum, "COUNTS", Pha, true);
  }

  if (!Poisserr) SPwriteCol(spectrum, "STAT_ERR", StatError, true);

  SPwriteCol(spectrum, "SYS_ERR", SysError);

  SPwriteCol(spectrum, "QUALITY", Quality);

  SPwriteCol(spectrum, "GROUPING", Group);

  SPwriteCol(spectrum, "AREASCAL", AreaScaling);

  SPwriteCol(spectrum, "BACKSCAL", BackScaling);

  return(OK);
}

// Write spectrum as type I file copying extra keywords and extensions from
// another file

Integer pha::write(string filename, string copyfilename)
{
  return(this->write(filename, copyfilename, 1));
}

// Write spectrum as type I file copying extra keywords and extensions from
// another file. The required spectrum is the HDUnumber instance of a spectrum
// extension in the file.

Integer pha::write(string filename, string copyfilename, Integer HDUnumber)
{
  Integer Status(OK);

  Status = this->write(filename);

  if ( Status != OK ) {
    string msg = "Failed to write pha to "+filename;
    SPreportError(Status, msg);
    return(Status);
  }

  Status = SPcopyKeys(copyfilename, filename, "SPECTRUM", HDUnumber);

  if ( Status != OK ) {
    string msg = "Failed to copy additional keywords from "+copyfilename+" to "+filename;
    SPreportError(Status, msg);
    return(Status);
  }

  Status = SPcopyHDUs(copyfilename, filename);
  if ( Status != OK ) {
    string msg = "Failed to copy additional HDUs from "+copyfilename+" to "+filename;
    SPreportError(Status, msg);
    return(Status);
  }

  return(Status);
}

// Multiply by a constant

pha& pha::operator*=(const Real a)
{

  // Loop round channels multiplying by a constant factor

  for (size_t i=0; i<Pha.size(); i++) Pha[i] *= a;
  if ( StatError.size() > 0 ) for (size_t i=0; i<StatError.size(); i++) StatError[i] *= a;

  return *this;

}

// Add to another pha

pha& pha::operator+=(const pha& a)
{
  // should not get this but if sizes of spectra are different then give up

  if ( Pha.size() != a.Pha.size() ) return *this;

  // set up Pha and Error arrays in counts. if Poisserr then explicitly calculate
  // error in case we need it later

  vector<Real> Pha1(Pha);
  vector<Real> Error1(StatError);
  vector<Real> Pha2(a.Pha);
  vector<Real> Error2(a.StatError);

  if ( Datatype != "COUNT" ) {
    for (size_t i=0; i<Pha1.size(); i++) Pha1[i] *= Exposure;
    for (size_t i=0; i<Error1.size(); i++) Error1[i] *= Exposure;
  }
  if ( a.Datatype != "COUNT" ) {
    for (size_t i=0; i<Pha2.size(); i++) Pha2[i] *= Exposure;
    for (size_t i=0; i<Error2.size(); i++) Error2[i] *= Exposure;
  }

  if ( Poisserr ) {
    Error1.resize(Pha1.size());
    for (size_t i=0; i<Error1.size(); i++) Error1[i] = sqrt(Pha1[i]);
  }
  if ( a.Poisserr ) {
    Error2.resize(Pha2.size());
    for (size_t i=0; i<Error2.size(); i++) Error2[i] = sqrt(Pha2[i]);
  }
    
  // Sum the counts

  for (size_t i=0; i<Pha1.size(); i++) Pha1[i] += Pha2[i];

  // If either one of the spectra are not Poisson errors then add in quadrature

  if ( !(Poisserr && a.Poisserr) ) {
    for (size_t i=0; i<Error1.size(); i++) Error1[i] = sqrt(Error1[i]*Error1[i] + Error2[i]*Error2[i]);
  }

  // For systematic errors take the maximum for each channel

  vector<Real> Sys1(Pha1.size());
  vector<Real> Sys2(Pha2.size());

  if ( SysError.size() == 0 ) {
      for (size_t i=0; i<Sys1.size(); i++) Sys1[i] = 0.0;
  } else if ( SysError.size() == 1 ) {
      for (size_t i=0; i<Sys1.size(); i++) Sys1[i] = SysError[0];
  } else {
      for (size_t i=0; i<Sys1.size(); i++) Sys1[i] = SysError[i];
  }

  if ( a.SysError.size() == 0 ) {
      for (size_t i=0; i<Sys2.size(); i++) Sys2[i] = 0.0;
  } else if ( SysError.size() == 1 ) {
      for (size_t i=0; i<Sys2.size(); i++) Sys2[i] = a.SysError[0];
  } else {
      for (size_t i=0; i<Sys2.size(); i++) Sys2[i] = a.SysError[i];
  }

  Real SysMax(0.0);
  bool SysVec(false);
  for (size_t i=0; i<Sys1.size(); i++) {
    if ( Sys1[i] < Sys2[i] ) Sys1[i] = Sys2[i];
    if ( Sys1[i] > SysMax ) {
      if ( SysMax > 0.0 ) SysVec = true;
      SysMax = Sys1[i];
    }
  }

  // Sum the exposures. This may not be the right thing to do but higher level routines
  // can correct this later if necessary

  Exposure += a.Exposure;

  // Now put back in the object arrays

  if ( Datatype == "COUNT" ) {
    for (size_t i=0; i<Pha.size(); i++) Pha[i] = Pha1[i];
  } else {
    for (size_t i=0; i<Pha.size(); i++) Pha[i] = Pha1[i]/Exposure;
  }

  Poisserr = Poisserr && a.Poisserr;
  if ( !Poisserr ) {
    StatError.resize(Error1.size());
    if ( Datatype == "COUNT" ) {
      for (size_t i=0; i<Pha.size(); i++) StatError[i] = Error1[i];
    } else {
      for (size_t i=0; i<Pha.size(); i++) StatError[i] = Error1[i]/Exposure;
    }
  }

  if ( SysMax > 0.0 ) {
    if ( !SysVec ) {
      SysError.resize(1);
      SysError[0] = Sys1[0];
    } else {
      SysError.resize(Sys1.size());
      for (size_t i=0; i<Sys1.size(); i++) SysError[i] = Sys1[i];
    }
  }
    
  return *this;
}

// Check compatibility with another pha

Integer pha::checkCompatibility(const pha& a)
{
  // check that the channels match up

  if ( Pha.size() != a.Pha.size() ) return InconsistentChannels;
  if ( FirstChannel != a.FirstChannel ) return InconsistentChannels;
  for (size_t i=0; i<Pha.size(); i++) {
    if ( Channel[i] != a.Channel[i] ) return InconsistentChannels;
  }

  if ( FluxUnits != a.FluxUnits ) return InconsistentUnits;

  return OK;
}

// Select a subset of the channels

Integer pha::selectChannels(vector<Integer>& Start, vector<Integer>& End)
{
  vector<Real> tPha, tStatError, tSysError, tAreaScaling, tBackScaling;
  vector<Integer> tChannel, tQuality, tGroup;

  size_t iset = 0;
  bool first = true;
  size_t ich = 0;
  while ( ich<Channel.size() && iset<Start.size() ) {
    if ( Channel[ich] >= Start[iset] && Channel[ich] <= End[iset] ) {
      tPha.push_back(Pha[ich]);
      tStatError.push_back(StatError[ich]);
      tSysError.push_back(SysError[ich]);
      tAreaScaling.push_back(AreaScaling[ich]);
      tBackScaling.push_back(BackScaling[ich]);
      tChannel.push_back(Channel[ich]);
      tQuality.push_back(Quality[ich]);
      if ( first ) {
	tGroup.push_back(1);
	first = false;
      } else {
	tGroup.push_back(Group[ich]);
      }
    }
    while ( ich > (size_t)End[iset] ) {
      iset++;
      first = true;
    }
    ich++;
  }

  // Reset the arrays

  size_t Nnew = tPha.size();
  Pha.resize(Nnew);
  StatError.resize(Nnew);
  SysError.resize(Nnew);
  AreaScaling.resize(Nnew);
  BackScaling.resize(Nnew);
  Channel.resize(Nnew);
  Quality.resize(Nnew);
  Group.resize(Nnew);
  for ( size_t i=0; i<Nnew; i++) {
    Pha[i] = tPha[i];
    StatError[i] = tStatError[i];
    SysError[i] = tSysError[i];
    AreaScaling[i] = tAreaScaling[i];
    BackScaling[i] = tBackScaling[i];
    Channel[i] = tChannel[i];
    Quality[i] = tQuality[i];
    Group[i] = tGroup[i];
  }

  return OK;
}

// Set grouping array from Grouping object

Integer pha::setGrouping(grouping& GroupInfo)
{

  // check for compatibility

  if ( NumberChannels() != GroupInfo.size() ) {
    stringstream msg;
    msg << "Grouping array is size " << GroupInfo.size() << " but number of channels is " << NumberChannels();
    SPreportError(InconsistentGrouping, msg.str());
    return(InconsistentGrouping);
  }

  // loop through channels setting grouping array

  for (size_t i=0; i<(size_t)NumberChannels(); i++) Group[i] = GroupInfo.flag[i];

  // check for the presence of any bad channels signalled by input grouping flag
  // being zero

  bool isBad(false);
  for (size_t i=0; i<(size_t)NumberChannels(); i++) {
    if ( GroupInfo.flag[i] == 0 ) isBad = true;
  }

  // if isBad is set then need to change the Quality array. First resize if
  // necessary

  if ( isBad ) {

    // Loop over the quality array setting to 2 any channels for which
    // Group is 0 - also reset Group to 1 (this is consistent with old grppha
    // behaviour)

    for (size_t i=0; i<(size_t)NumberChannels(); i++) {
      if ( Group[i] == 0 ) {
	Quality[i] = 2;
	Group[i] = 1;
      }
    }

  }

  return(OK);
}

// Get grouping object from the grouping array

grouping pha::getGrouping()
{
  grouping groupInfo;
  groupInfo.loadFromVector(Quality, Group);
  return groupInfo;
}

// Get grouping between channels StartChannel and EndChannel using a 
// minimum number of counts per bin

grouping pha::getMinCountsGrouping(const Integer MinCounts, const Integer StartChannel,
			     const Integer EndChannel)
{

  size_t Nchan = NumberChannels();
  // first need to convert Pha data to counts if the type is RATE
  vector<Integer> dataCounts(Nchan);
  if ( Datatype == "RATE" ) {
    for (size_t i=0; i<Nchan; i++) dataCounts[i] = (Integer)round(Pha[i] * Exposure);
  } else {
    for (size_t i=0; i<Nchan; i++) dataCounts[i] = (Integer)round(Pha[i]);
  }

  // convert start and end channel to zero-based
  Integer start = StartChannel - FirstChannel;
  Integer end = EndChannel - FirstChannel;
  // calculate the grouping
  grouping groupInfo;
  groupInfo.loadMin(MinCounts, start, end, dataCounts);

  // and return it

  return groupInfo;
}

grouping pha::getMinCountsGrouping(const Integer MinCounts)
{
  Integer StartChannel = FirstChannel;
  Integer EndChannel = FirstChannel + NumberChannels() - 1;
  return getMinCountsGrouping(MinCounts, StartChannel, EndChannel);
}

// Get grouping (optionally between channels StartChannel and EndChannel) 
// using a minimum S/N. Optionally includes a background file as well.

grouping pha::getMinSNGrouping(const Real SignalToNoise, const Integer StartChannel,
			       const Integer EndChannel, const pha& Background)
{
  size_t Nchan(NumberChannels());

  // make vectors containing signal and noise
  vector<Real> Signal(Nchan), Noise(Nchan);

  for (size_t i=0; i<Nchan; i++) Signal[i] = Pha[i];

  if ( Poisserr ) {
    if ( Datatype == "RATE" ) {
      for (size_t i=0; i<Nchan; i++) Noise[i] = sqrt(Pha[i]*Exposure)/Exposure;
    } else {
      for (size_t i=0; i<Nchan; i++) Noise[i] = sqrt(Pha[i]);
    }
  } else {
    for (size_t i=0; i<Nchan; i++) Noise[i] = StatError[i];
  }

  if ( SysError.size() != 0 ) {
    for (size_t i=0; i<Nchan; i++) Noise[i] = sqrt(Noise[i]*Noise[i] + SysError[i]*SysError[i]*Pha[i]*Pha[i]);
  }

  // including the background 

  if ( Background.Pha.size() != 0 ) {

    // useful to construct the scaling argument for the background here
    vector<Real> Scale(Nchan);
    for (size_t i=0; i<Nchan; i++) {
      Scale[i] = AreaScaling[i]*BackScaling[i]
	/(Background.AreaScaling[i]*Background.BackScaling[i]);
      if ( Datatype == "COUNTS" ) {
	Scale[i] *= Exposure/Background.Exposure;
      }
    }

    // construct the background noise
    
    vector<Real> BackgroundNoise(Nchan);
    if ( Background.Poisserr ) {
      if ( Background.Datatype == "RATE" ) {
	for (size_t i=0; i<Nchan; i++) BackgroundNoise[i] = sqrt(Background.Pha[i]*Background.Exposure)/Background.Exposure;
      } else {
	for (size_t i=0; i<Nchan; i++) BackgroundNoise[i] = sqrt(Background.Pha[i]);
      }
    } else {
      for (size_t i=0; i<Nchan; i++) BackgroundNoise[i] = Background.StatError[i];
    }

    if ( Background.SysError.size() != 0 ) {
      for (size_t i=0; i<Nchan; i++) BackgroundNoise[i] = sqrt(BackgroundNoise[i]*BackgroundNoise[i] + Background.SysError[i]*Background.SysError[i]*Background.Pha[i]*Background.Pha[i]);
    }
    
    // substract background from the signal

    for (size_t i=0; i<Nchan; i++) {
      Signal[i] -= Background.Pha[i]*Scale[i];
      Noise[i] = sqrt(Noise[i]*Noise[i] 
	  + BackgroundNoise[i]*Scale[i]*BackgroundNoise[i]*Scale[i]);
    }

  }

  // convert start and end channel to zero-based
  Integer First = StartChannel - FirstChannel;
  Integer Last = EndChannel - FirstChannel;

  // set up the groupInfo.flag array to get the minimum signal to noise requested

  grouping groupInfo;

  bool ingroup = false;
  size_t istart = 0;
  Real signalSum = 0.0;
  Real noiseSum = 0.0;
  Real ratio;
  groupInfo.flag.resize(Nchan);
  for (size_t i=0; i<First; i++) groupInfo.flag[i] = 1;
  for (size_t i=First; i<=Last; i++) {
    signalSum += Signal[i];
    noiseSum = sqrt(noiseSum*noiseSum + Noise[i]*Noise[i]);
    ratio = 0.0;
    if ( noiseSum > 0.0 ) ratio = signalSum / noiseSum;
    if ( !ingroup ) {
      groupInfo.flag[i] = 1;
      istart = i;
      signalSum = Signal[i];
      noiseSum = Noise[i];
      ratio = 0.0;
      if ( noiseSum > 0.0 ) ratio = signalSum / noiseSum;
      if ( ratio < SignalToNoise ) ingroup = true;
    } else if ( ratio >= SignalToNoise ) {
      groupInfo.flag[i] = -1;
      ingroup = false;
    } else if ( ratio < SignalToNoise ) {
      groupInfo.flag[i] = -1;
    }
  }
  if ( ratio < SignalToNoise && ingroup ) {
    for (size_t i=istart; i<=Last; i++) groupInfo.flag[i] = 0;
  }
  for (size_t i=Last+1; i<Nchan; i++) groupInfo.flag[i] = 1;

  
  return groupInfo;
}

grouping pha::getMinSNGrouping(const Real SignalToNoise, const pha& Background)
{
  Integer StartChannel = FirstChannel;
  Integer EndChannel = FirstChannel + NumberChannels() - 1;
  return getMinSNGrouping(SignalToNoise, StartChannel, EndChannel, Background);
}

grouping pha::getMinSNGrouping(const Real SignalToNoise, const Integer StartChannel,
			       const Integer EndChannel)
{
  pha Background;
  return getMinSNGrouping(SignalToNoise, StartChannel, EndChannel, Background);
}

grouping pha::getMinSNGrouping(const Real SignalToNoise)
{
  Integer StartChannel = FirstChannel;
  Integer EndChannel = FirstChannel + NumberChannels() - 1;
  pha Background;
  return getMinSNGrouping(SignalToNoise, StartChannel, EndChannel, Background);
}

// Rebin channels based on the Grouping object

Integer pha::rebinChannels(grouping& GroupInfo)
{
  return this->rebinChannels(GroupInfo, "PROPAGATE");
}

// Rebin channels based on the Grouping object

Integer pha::rebinChannels(grouping& GroupInfo, string errorType)
{
  // check for compatibility

  if ( NumberChannels() != GroupInfo.size() ) {
    stringstream msg;
    msg << "Grouping array is size " << GroupInfo.size() << " but number of channels is " << NumberChannels();
    SPreportError(InconsistentGrouping, msg.str());
    return(InconsistentGrouping);
  }

  vector<Real> temp;
  Integer NumberOrigChannels(NumberChannels());

  // bin up the data

  GroupBin(Pha, SumMode, GroupInfo, temp);
  Pha.resize(temp.size());
  Pha = temp;
  temp.clear();

  // if necessary sum in quadrature the error

  if ( errorType == "PROPAGATE" ) {

    if ( !Poisserr ) {

      GroupBin(StatError, SumQuadMode, GroupInfo, temp);
      StatError.resize(temp.size());
      StatError = temp;
      temp.clear();

    }

  } else {

    // otherwise calculate the error using the requested option

    Real factor(1.0);
    if ( Datatype == "RATE" ) factor = Exposure;
    StatError.resize(Pha.size());
    Poisserr = false;

    if ( errorType == "GAUSS" ) {

      for (size_t i=0; i<StatError.size(); i++) StatError[i] = 1.0/sqrt(Pha[i]);

    } else if ( errorType == "POISS-0" ) {

      StatError.resize(0);
      Poisserr = true;

    } else if ( errorType == "POISS-1" ) {

      for (size_t i=0; i<StatError.size(); i++) StatError[i] = 1.0 + sqrt(factor*Pha[i]+0.75);

    } else if ( errorType == "POISS-2" ) {

      for (size_t i=0; i<StatError.size(); i++) StatError[i] = sqrt(factor*Pha[i]-0.25);

    } else if ( errorType == "POISS-3" ) {

      for (size_t i=0; i<StatError.size(); i++) StatError[i] = 0.5*(1.0 + sqrt(factor*Pha[i]+0.75)+sqrt(factor*Pha[i]-0.25));

    }

  }

  // if necessary sum in quadrature the systematic error

  if ( SysError.size() > 0 ) {

    GroupBin(SysError, SumQuadMode, GroupInfo, temp);
    SysError.resize(temp.size());
    SysError = temp;
    temp.clear();

  }

  // reset Channel array - need to take into account that not all channels
  // may be given values in the original spectrum

  Integer FirstOrigChannel(Channel[0]);
  Channel.resize(NumberChannels());
  for (size_t i=0; i<(size_t)NumberChannels(); i++) {
    Channel[i] = FirstOrigChannel + i;
  }
  DetChans -= NumberOrigChannels - NumberChannels();

  // redo Quality array - if any element making up a bin has bad quality then
  // the whole bin does

  vector<Integer> newQual;
  for (size_t i=0; i<Quality.size(); i++) {
    if ( GroupInfo.newBin(i) ) {
      newQual.push_back(Quality[i]);
    } else {
      if ( Quality[i] != 0 ) newQual[newQual.size()-1] = Quality[i];
    }
  }
  Quality.resize(newQual.size());
  for (size_t i=0; i<Quality.size(); i++) Quality[i] = newQual[i];

  // it doesn't make sense to try to preserve the grouping so reset

  Group.resize(NumberChannels());
  for (size_t i=0; i<Group.size(); i++) Group[i] = 1;

  // if there are AreaScaling or BackScaling arrays then average them

  if ( AreaScaling.size() > 1 ) {

    GroupBin(AreaScaling, MeanMode, GroupInfo, temp);
    AreaScaling.resize(temp.size());
    AreaScaling = temp;
    temp.clear();

  }

  if ( BackScaling.size() > 1 ) {

    GroupBin(BackScaling, MeanMode, GroupInfo, temp);
    BackScaling.resize(temp.size());
    BackScaling = temp;
    temp.clear();

  }

  return(OK);
}

// VS:
Integer pha::initChannels(Integer nchan) {
    FirstChannel = 0;

    Pha.clear();
    StatError.clear();
    SysError.clear();
    Channel.clear();
    Quality.clear();
    Group.clear();
    AreaScaling.clear();
    BackScaling.clear();

    Exposure = 0.0;
    CorrectionScaling = 0.0;

    DetChans = nchan;
    Poisserr = false;
    Datatype = "COUNT";
    PHAVersion = "1.2.1";
    Spectrumtype = "TOTAL";
    ResponseFile = "NONE";
    AncillaryFile = "NONE";
    BackgroundFile = "NONE";
    CorrectionFile = "NONE";
    FluxUnits = " ";
    ChannelType = "PI";
    Telescope = "mock_telescope";
    Instrument = "mock_instrument";
    Detector = "mock_detector";
    Filter = " ";
    Datamode = " ";

    XSPECFilter.clear();

    for (Integer i=0;i<nchan;i++) {
        Pha.push_back(0);
        StatError.push_back(0);
        SysError.push_back(0);
        Channel.push_back(i);
        Quality.push_back(0);
        Group.push_back(1);
        AreaScaling.push_back(0);
        BackScaling.push_back(0);
    }
}

// :VS

// Remaps counts up or down in channels.

Integer pha::shiftChannels(Integer Start, Integer End, Real Shift)
{
  Real Factor(1.0);
  return this->shiftChannels(Start, End, Shift, Factor);
}

Integer pha::shiftChannels(Integer Start, Integer End, Real Shift, Real Factor)
{
  vector<Integer> vStart(1,Start);
  vector<Integer> vEnd(1,End);
  vector<Real> vShift(1,Shift);
  vector<Real> vFactor(1,Factor);
  return this->shiftChannels(vStart, vEnd, vShift, vFactor);
}

Integer pha::shiftChannels(vector<Integer>& vStart, vector<Integer>& vEnd, vector<Real>& vShift, vector<Real>& vFactor)
{
  size_t Nchan(Channel.size());
  vector<Real> Low(Nchan);
  vector<Real> High(Nchan);
  for (size_t i=0; i<Nchan; i++) {
    Low[i] = Channel[i] - 0.5;
    High[i] = Channel[i] + 0.5;
  }
  return this->shiftChannels(Low, High, vStart, vEnd, vShift, vFactor);
}

Integer pha::shiftChannels(vector<Real>& Low, vector<Real>& High,
			   Integer Start, Integer End, Real Shift, Real Factor)
{
  vector<Integer> vStart(1,Start);
  vector<Integer> vEnd(1,End);
  vector<Real> vShift(1,Shift);
  vector<Real> vFactor(1,Factor);
  return this->shiftChannels(Low, High, vStart, vEnd, vShift, vFactor);
}

Integer pha::shiftChannels(vector<Real>& Low, vector<Real>& High, vector<Integer>& vStart,
			   vector<Integer>& vEnd, vector<Real>& vShift, 
			   vector<Real>& vFactor)
{
  size_t Nchan(Channel.size());

  // SPcalcShift requires vStart and vEnd to refer to element numbers within the
  // Low and High arrays so need to correct for the starting channel number
  vector<Integer> vStart0(vStart);
  vector<Integer> vEnd0(vEnd);
  if ( Channel[0] != 0 ) {
    for (size_t i=0; i<vStart0.size(); i++) {
      vStart0[i] -= Channel[0];
      vEnd0[i] -= Channel[0];
    }
  }
  
  // First set up vectors describing how to make the new pha
  // fromChannelElt[i] is the (zero-based) list of channel elements from which
  // the new channel i is calculated fromFraction[i] is the list of fractions
  // corresponding to fromChannel.

  vector<vector<size_t> > fromChannelElt(Nchan);
  vector<vector<Real> > fromFraction(Nchan);

  SPcalcShift(Low, High, vStart0, vEnd0, vShift, vFactor, fromChannelElt, fromFraction);

  // booleans to determine whether other vector quantities exist to be
  // shifted.

  bool isError(StatError.size() > 0);
  bool isSysError(SysError.size() > 0);
  bool isQuality(Quality.size() > 0);
  bool isGroup(Group.size() > 0);
  bool isArea(AreaScaling.size() > 0);
  bool isBack(BackScaling.size() > 0);

  // set up temporary arrays to store the output data and other quantities

  vector<Real> TmpPha(Nchan,0.0), TmpError(Nchan,0.0), TmpSysError(Nchan,0.0);
  vector<Integer> TmpQuality(Nchan,0), TmpGroup(Nchan,0);
  vector<Real> TmpArea(Nchan,0.0), TmpBack(Nchan,0.0);

  // loop over output channels

  for (size_t iChan=0; iChan<Nchan; iChan++) {

    for (size_t j=0; j<fromChannelElt[iChan].size(); j++) {
      size_t fromChan = fromChannelElt[iChan][j];
      Real fromFrac = fromFraction[iChan][j];
      TmpPha[iChan] += fromFrac*Pha[fromChan];
      if ( isError ) {
	TmpError[iChan] += fromFrac*StatError[fromChan]*fromFrac*StatError[fromChan];
      }
      if ( isSysError ) TmpSysError[iChan] += fromFrac*SysError[fromChan];
      if ( isQuality ) {
	if ( Quality[fromChan] != 0 ) TmpQuality[iChan] = 0;
      }
      if ( isGroup ) { 
	if ( Group[fromChan] == 1 ) TmpGroup[iChan] = 1;
      }
      if ( isArea ) TmpArea[iChan] += fromFrac*AreaScaling[fromChan];
      if ( isBack ) TmpBack[iChan] += fromFrac*BackScaling[fromChan];
    }

  }

  // transfer to object arrays

  for (size_t iChan=0; iChan<Nchan; iChan++) Pha[iChan] = TmpPha[iChan];
  if ( isError ) {
    for (size_t iChan=0; iChan<Nchan; iChan++) StatError[iChan] = sqrt(TmpError[iChan]);
  }
  if ( isSysError ) {
    for (size_t iChan=0; iChan<Nchan; iChan++) SysError[iChan] = TmpSysError[iChan];
  }
  if ( isQuality ) {
    for (size_t iChan=0; iChan<Nchan; iChan++) Quality[iChan] = TmpQuality[iChan];
  }
  if ( isGroup ) {
    for (size_t iChan=0; iChan<Nchan; iChan++) Group[iChan] = TmpGroup[iChan];
  }
  if ( isArea ) {
    for (size_t iChan=0; iChan<Nchan; iChan++) AreaScaling[iChan] = TmpArea[iChan];
  }
  if ( isBack ) {
    for (size_t iChan=0; iChan<Nchan; iChan++) BackScaling[iChan] = TmpBack[iChan];
  }

  return(OK);

}

// Convert flux units from whatever they are currently to ph/cm^/s. No need to use this if data are just
// in counts or counts/s.

Integer pha::convertUnits(vector<Real>& ChannelLowEnergy, vector<Real>& ChannelHighEnergy, string EnergyUnits)
{

  if ( FluxUnits == " " || FluxUnits == "counts" || FluxUnits == "counts/s" ) return(OK);

  // Find the flux units set and check for validity

  bool energy;
  bool perwave;
  bool perenergy;
  Real yfactor;

  Integer status(OK);

  status = calcYfactor(FluxUnits, energy, perwave, perenergy, yfactor);
  if ( status != OK ) {
    string msg = "Failed to calculate conversion factor for "+FluxUnits;
    SPreportError(status, msg);
    return(status);
  }

  // set up energy/wave conversion factors and check for valid units

  bool xwave;
  Real xfactor;

  status = calcXfactor(EnergyUnits, xwave, xfactor);
  if ( status != OK ) {
    string msg = "Failed to calculate conversion factor for "+EnergyUnits;
    SPreportError(status, msg);
    return(status);
  }

  // if nothing to be done then return

  if ( xfactor == 1.0 && yfactor == 1.0 ) return(OK);

  // set up arrays of mean energies/wavelengths and bin sizes
  // both are in keV.

  vector<Real> xmean(ChannelLowEnergy.size());
  vector<Real> xgmult(ChannelLowEnergy.size());
  vector<Real> xbinsize(ChannelLowEnergy.size());
  Real x2factor = xfactor*xfactor;

  for (size_t i=0; i<xmean.size(); i++) {
    if ( xwave ) {
      xmean[i] = xfactor*(1.0/ChannelLowEnergy[i] + 1.0/ChannelHighEnergy[i])/2.0;
      xgmult[i] = x2factor*((1.0/ChannelLowEnergy[i]) * (1.0/ChannelHighEnergy[i]));
      xbinsize[i] = xfactor*abs(1.0/ChannelHighEnergy[i] - 1.0/ChannelLowEnergy[i]);
    } else {
      xmean[i] = xfactor * (ChannelLowEnergy[i]+ChannelHighEnergy[i])/2.0;
      xbinsize[i] = xfactor * abs(ChannelHighEnergy[i]-ChannelLowEnergy[i]);
    }
  }

  // now do the conversion of Pha and StatError. Note the six different possibilities based on
  // whether flux is in photons or energy, is per wavelength or not, is per energy or not (both
  // per wave and per energy can be true at the same time).

  if ( !energy && !perwave && !perenergy ) {
    // this is the one we have already done - flux must be ph/cm^2/s
  } else if ( !energy && !perwave && perenergy ) {
    for (size_t i=0; i<Pha.size(); i++) Pha[i] *= yfactor*xbinsize[i];
    for (size_t i=0; i<StatError.size(); i++) StatError[i] *= yfactor*xbinsize[i];
  } else if ( !energy && perwave && !perenergy ) {
    for (size_t i=0; i<Pha.size(); i++) Pha[i] *= yfactor*xbinsize[i]/xgmult[i];
    for (size_t i=0; i<StatError.size(); i++) StatError[i] *= yfactor*xbinsize[i]/xgmult[i];
  } else if ( energy && !perwave && !perenergy ) {
    for (size_t i=0; i<Pha.size(); i++) Pha[i] *= yfactor/xmean[i];
    for (size_t i=0; i<StatError.size(); i++) StatError[i] *= yfactor/xmean[i];
  } else if ( energy && !perwave && perenergy ) {
    for (size_t i=0; i<Pha.size(); i++) Pha[i] *= yfactor*xbinsize[i]/xmean[i];
    for (size_t i=0; i<StatError.size(); i++) StatError[i] *= yfactor*xbinsize[i]/xmean[i];
  } else if ( energy && perwave && !perenergy ) {
    for (size_t i=0; i<Pha.size(); i++) Pha[i] *= yfactor*xbinsize[i]/xmean[i]/xgmult[i];
    for (size_t i=0; i<StatError.size(); i++) StatError[i] *= yfactor*xbinsize[i]/xmean[i]/xgmult[i];
  }

  //debug
  //  for (size_t i=0; i<xmean.size(); i++) cout << xmean[i] << " " << xgmult[i] << " " << xbinsize[i] << " " << Pha[i] << std::endl;

  FluxUnits = "ph/cm^2/s";
  return(OK);
}

//*******************************************************************************
// binary operations

pha operator+ (const pha& a, const pha& b)
{
  pha c(a);
  return c += b;
}

//*******************************************************************************
// Some utility routines. Definitions in pha.h


// Return the type of a PHA extension

// Integer PHAtype(string filename, Integer PHAnumber)

Integer PHAtype(string filename, Integer PHAnumber)
{
  Integer Status(OK);
  return(PHAtype(filename, PHAnumber, Status));
}

Integer PHAtype(string filename, Integer PHAnumber, Integer& Status)
{
  if ( Status != OK ) return 0;

  const string hduName("SPECTRUM");
  string DefString;

  // Read in the SPECTRUM extension number PHAnumber
  // and set up an object called spectrum with the contents

  const vector<string> hduKeys;
  const vector<string> primaryKey;

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(filename,Read,hduName,false,hduKeys,primaryKey,(int)PHAnumber));
  } catch(...) {
    Status = NoSuchFile;
    string msg = "Failed to read "+hduName+" in "+filename;
    SPreportError(Status, msg);
    return(0);
  }

  ExtHDU& spectrum = pInfile->extension(hduName, (int)PHAnumber);
  
  // Check for a HDUCLAS4 keyword of "TYPEII"

  DefString = "TYPEI";
  try {
    string hduclas4 = SPreadKey(spectrum, "HDUCLAS4", 1, DefString);
    if ( hduclas4 == "TYPEII" ) return 2;
  } catch (Table::NoSuchKeyword&) {
  }

  // Check for a type II extension by looking for the SPEC_NUM column

  try {
    Column& Col = spectrum.column("SPEC_NUM");
    if ( Col.rows() > 0 ) return 2;
  } catch (Table::NoSuchColumn&) {
  }

  // If there is no SPEC_NUM column then look for a RATE or COUNTS column
  // and check the TFORM

  try {
    Column& Col = spectrum.column("COUNTS");
    if ( Col.varLength() || Col.repeat() > 1 ) return 2;
  } catch (Table::NoSuchColumn&) {
    try {
      Column& Col = spectrum.column("RATE");
      if ( Col.varLength() || Col.repeat() > 1 ) return 2;
    } catch (Table::NoSuchColumn&) {
    }
  }
       
  // None of the checks for a type II extension panned out so it must be a type I

  return 1;
} 

// return 0 if COUNTS column exists and is integer or COUNTS column does not exist

bool IsPHAcounts(string filename, Integer PHAnumber)
{
  Integer Status(OK);
  return(IsPHAcounts(filename, PHAnumber, Status));
}

bool IsPHAcounts(string filename, Integer PHAnumber, Integer& Status)
{

  if ( Status != OK ) return 0;

  const string hduName("SPECTRUM");

  // Read in the SPECTRUM extension number PHAnumber
  // and set up an object called spectrum with the contents

  const vector<string> hduKeys;
  const vector<string> primaryKey;

  auto_ptr<FITS> pInfile(0);

  try{
    pInfile.reset(new FITS(filename,Read,hduName,false,hduKeys,primaryKey,(int)PHAnumber));
  } catch(...) {
    Status = NoSuchFile;
    string msg = "Failed to read "+hduName+" in "+filename;
    SPreportError(Status, msg);
    return(false);
  }

  ExtHDU& spectrum = pInfile->extension(hduName, (int)PHAnumber);

  try {
    Column& Col = spectrum.column("COUNTS");
    ValueType Type = Col.type();
    if ( Type == TSHORT || Type == TLONG ) return true;
  } catch (Table::NoSuchColumn&) {
  }

  return false;

}

// return the number of spectra in a type II PHA extension

Integer NumberofSpectra(string filename, Integer PHAnumber)
{
  Integer Status(OK);
  return(NumberofSpectra(filename, PHAnumber, Status));
}

Integer NumberofSpectra(string filename, Integer PHAnumber, Integer& Status)
{

  if ( Status != 0 ) return 0;

  const string hduName("SPECTRUM");

  // Read in the SPECTRUM extension number PHAnumber
  // and set up an object called spectrum with the contents

  const vector<string> hduKeys;
  const vector<string> primaryKey;

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(filename,Read,hduName,false,hduKeys,primaryKey,(int)PHAnumber));
  } catch(...) {
    Status = NoSuchFile;
    string msg = "Failed to read "+hduName+" in "+filename;
    SPreportError(Status, msg);
    return(0);
  }

  ExtHDU& spectrum = pInfile->extension(hduName, (int)PHAnumber);

  return (Integer)spectrum.rows();

}

vector<Integer> SpectrumExtensions(string filename)
{
  Integer Status(0);
  return SpectrumExtensions(filename, Status);
}

vector<Integer> SpectrumExtensions(string filename, Integer& Status)
{
  string value("SPECTRUM");
  string keyName("EXTNAME");
  return SPfindExtensions(filename, keyName, value, Status);
}

