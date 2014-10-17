// ResponseMatrix object code. Definitions in ResponseMatrix.h

#ifndef HAVE_rmf
#include "rmf.h"
#endif

#ifndef HAVE_arf
#include "arf.h"
#endif

#ifndef HAVE_grouping
#include "grouping.h"
#endif

#ifndef HAVE_SPio
#include "SPio.h"
#endif

#ifndef HAVE_SPutils
#include "SPutils.h"
#endif

// Class rmf

// default constructor

rmf::rmf()
{
}


// Destructor

rmf::~rmf()
{
}

// reading Matrix and Channel bounds extensions from RMF file. 

Integer rmf::read(string filename)
{
  return(this->read(filename, 1));
}

// reading Matrix and Channel bounds extensions from RMF file. This option allows multiple 
// extensions in the same file

Integer rmf::read(string filename, Integer RMFnumber)
{
  Integer Status(OK);
  
  // read the MATRIX extension

  Status = this->readMatrix(filename, RMFnumber);
  if ( Status != OK ) {
    return(Status);
};

  // try to read the EBOUNDS extension. If this fails try to read with RMFnumber = 1
  // under assumption that there are multiple MATRIX extensions and only one EBOUNDS
  // extension

  Status = this->readChannelBounds(filename, RMFnumber);
  if ( Status != OK ) {
    Status = this->readChannelBounds(filename, 1);
  }

  return(Status);
}

// reading channel bounds extension from RMF file. 

Integer rmf::readChannelBounds(string filename)
{
  return(this->readChannelBounds(filename, 1));
}

// reading channel bounds extension from PHA file. this option allows multiple extensions in the same file

Integer rmf::readChannelBounds(string filename, Integer EBDnumber)
{
  string hduName("EBOUNDS");
  string DefString;

  // Read in the Channel bounds extension number EBDnumber
  // and set up an object called ebd with the contents

  const vector<string> hduKeys;
  const vector<string> primaryKey;

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(filename,Read,hduName,false,hduKeys,primaryKey,(int)EBDnumber));
  } catch(...) {
    string msg = "Failed to read "+hduName+" in "+filename;
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  ExtHDU& ebd = pInfile->extension(hduName);
  
  // read the standard keywords and store in the object

  DefString = "UNKNOWN";
  ChannelType = SPreadKey(ebd, "CHANTYPE", DefString);

  EBDVersion = SPreadKey(ebd, "HDUVERS", DefString);
  if ( EBDVersion == "UNKNOWN" ) {
    EBDVersion = SPreadKey(ebd, "HDUVERS2", DefString);
    if ( EBDVersion == "UNKNOWN" ) {
      EBDVersion = SPreadKey(ebd, "EBDVERSN", DefString);
      if ( EBDVersion == "UNKNOWN" ) EBDVersion = "1.1.0";
    }
  }

  EBDExtensionName = SPreadKey(ebd, "EXTNAME", DefString);
  
  Telescope = SPreadKey(ebd, "TELESCOP", DefString);
  
  Instrument = SPreadKey(ebd, "INSTRUME", DefString);

  Detector = SPreadKey(ebd, "DETNAM", DefString);

  Filter = SPreadKey(ebd, "FILTER", DefString);

  // Read the E_MIN and E_MAX columns

  SPreadCol(ebd,"E_MIN",ChannelLowEnergy);
  SPreadCol(ebd,"E_MAX",ChannelHighEnergy);
  if ( ChannelLowEnergy.size() == 0 ) {
    string msg = "Failed to read any entries from E_MIN";
    SPreportError(NoEmin, msg);
    return(NoEmin);
  }
  if ( ChannelHighEnergy.size() == 0 ) {
    string msg = "Failed to read any entries from E_MAX";
    SPreportError(NoEmax, msg);
    return(NoEmax);
  }

  SPreadColUnits(ebd, "E_MIN", EnergyUnits);

  return(OK);
}



// reading Matrix extension from RMF file. 

Integer rmf::readMatrix(string filename)
{
  return(this->readMatrix(filename, 1));
}

// reading Matrix extension from PHA file. this option allows multiple extensions in the same file

Integer rmf::readMatrix(string filename, Integer RMFnumber)
{
  string hduName("MATRIX");
  string DefString;
  bool verbosity = FITS::verboseMode();
  FITS::setVerboseMode(true);

  // Read in the MATRIX extension number RMFnumber
  // and set up an object called rmf with the contents

  const vector<string> hduKeys;
  const vector<string> primaryKey;

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(filename,Read,hduName,false,hduKeys,primaryKey,(int)RMFnumber));
  } catch (...) {
    hduName = "SPECRESP MATRIX";
    FITS::clearErrors();
    try {
      pInfile.reset(new FITS(filename,Read,hduName,false,hduKeys,primaryKey,(int)RMFnumber));
    } catch(...) {
      string msg = "Failed to read MATRIX or SPECRESP MATRIX in "+filename;
      SPreportError(NoSuchFile, msg);
      return(NoSuchFile);
    }
  }

  ExtHDU& rmf = pInfile->extension(hduName);
  
  // read the standard keywords and store in the object

  DefString = "UNKNOWN";
  ChannelType = SPreadKey(rmf, "CHANTYPE", DefString);

  RMFVersion = SPreadKey(rmf, "HDUVERS", DefString);
  if ( RMFVersion == "UNKNOWN" ) {
    RMFVersion = SPreadKey(rmf, "HDUVERS2", DefString);
    if ( RMFVersion == "UNKNOWN" ) {
      RMFVersion = SPreadKey(rmf, "RMFVERSN", DefString);
      if ( RMFVersion == "UNKNOWN" ) RMFVersion = "1.1.0";
    }
  }

  RMFExtensionName = SPreadKey(rmf, "EXTNAME", DefString);
  
  Telescope = SPreadKey(rmf, "TELESCOP", DefString);
  
  Instrument = SPreadKey(rmf, "INSTRUME", DefString);

  Detector = SPreadKey(rmf, "DETNAM", DefString);

  Filter = SPreadKey(rmf, "FILTER", DefString);

  RMFType = SPreadKey(rmf, "HDUCLAS3", DefString);

  AreaScaling = SPreadKey(rmf, "EFFAREA", (Real)1.0);

  ResponseThreshold = SPreadKey(rmf, "LO_THRES", (Real)1.0);

  Integer Nrows;
  Nrows = SPreadKey(rmf, "NAXIS2", (Integer)0);

  Integer NtotGroups;
  NtotGroups = SPreadKey(rmf, "NUMGRP", (Integer)0);
  
  Integer NtotElts;
  NtotElts = SPreadKey(rmf, "NUMELT", (Integer)0);
 
  // Read the ENERG_LO and ENERG_HI columns

  SPreadCol(rmf,"ENERG_LO",LowEnergy);
  SPreadCol(rmf,"ENERG_HI",HighEnergy);
  if ( LowEnergy.size() == 0 ) {
    string msg = "Failed to read any entries from the ENERG_LO column";
    SPreportError(NoEnergLo, msg);
    return(NoEnergLo);
  }
  if ( HighEnergy.size() == 0 ) {
    string msg = "Failed to read any entries from the ENERG_HI column";
    SPreportError(NoEnergHi, msg);
    return(NoEnergHi);
  }

  SPreadColUnits(rmf, "ENERG_LO", EnergyUnits);

  // Get the number of groups for each energy bin

  SPreadCol(rmf,"N_GRP",NumberGroups);
  if ( NumberGroups.size() == 0 ) {
    string msg = "Failed to read any entries from the N_GRP column";
    SPreportError(NoNgrp, msg);
    return(NoNgrp);
  }

  // Set up the FirstGroup array - note that it counts from 0

  FirstGroup.resize(Nrows);
  Integer igrp = 0;
  for (size_t i=0; i<(size_t)Nrows; i++) {
    FirstGroup[i] = igrp;
    igrp += NumberGroups[i];
  }

  // Test for consistency between any value read from the NUMGRP keyword and the
  // sum of the N_GRP column.

  Integer Ntest(0);
  for (size_t i=0; i<(size_t)Nrows; i++) Ntest += NumberGroups[i];

  if ( NtotGroups != 0  && Ntest != NtotGroups ) {
    stringstream msg;
    msg << "The value of NUMGRP (" << NtotGroups << ") does not equal the sum of the N_GRP column (" << Ntest << ")";
    SPreportError(InconsistentNumgrp, msg.str());
    return(InconsistentNumgrp);
  }
  NtotGroups = Ntest;

  // Read the first channel and number of channels for each group - grab the whole data 
  // into temporary arrays then place in object. Also set the FirstElement array.

  vector<vector<Integer> > fchan, nchan;

  SPreadVectorCol(rmf, "F_CHAN", fchan);
  SPreadVectorCol(rmf, "N_CHAN", nchan);
  if ( fchan.size() == 0 ) {
    string msg = "Failed to read any entries from the F_CHAN column";
    SPreportError(NoFchan, msg);
    return(NoFchan);
  }
  if ( nchan.size() == 0 ) {
    string msg = "Failed to read any entries from the N_CHAN column";
    SPreportError(NoNchan, msg);
    return(NoNchan);
  }


  FirstChannelGroup.resize(NtotGroups);
  NumberChannelsGroup.resize(NtotGroups);
  FirstElement.resize(NtotGroups);
  size_t ipt = 0;
  size_t ielt = 0;
  for (size_t i=0; i<(size_t)Nrows; i++) {
    for (size_t j=0; j<(size_t)fchan[i].size(); j++) {
      if ( nchan[i][j] > 0 ) {
	FirstChannelGroup[ipt] = fchan[i][j];
	NumberChannelsGroup[ipt] = nchan[i][j];
	FirstElement[ipt] = ielt;
	ielt += NumberChannelsGroup[ipt];
	ipt++;
      }
    }
  }
  
  // Check for a TLMIN for the F_CHAN column

  try {
    int ChannelIndex = rmf.column("F_CHAN").index();
    ostringstream KeyStream;
    KeyStream << "TLMIN" << ChannelIndex;
    rmf.readKey(KeyStream.str(),FirstChannel);
    FirstChannel = SPreadKey(rmf, KeyStream.str(), (Integer)1);
  } catch(Table::NoSuchColumn&) {
    FirstChannel = 1;
  } catch(HDU::NoSuchKeyword&) {
    FirstChannel = 1;
  }

  // Test for consistency between any value read from the NUMELT keyword and the
  // sum of the N_CHAN column.

  Ntest = 0;
  for (size_t i=0; i<(size_t)NtotGroups; i++) Ntest += NumberChannelsGroup[i];

  if ( NtotElts != 0 && Ntest != NtotElts ) {
    stringstream msg;
    msg << "The value of NUMELT (" << NtotElts << ") does not equal the sum of the N_CHAN column (" << Ntest << ")";
    SPreportError(InconsistentNumelt, msg.str());
    return(InconsistentNumelt);
  }
  NtotElts = Ntest;

  // Read the matrix column into a temporary array then set the Matrix array
  // We have to be a bit careful here in case the file was created using fixed
  // length vectors. In this case the total size of the elements matrix read
  // will be larger than NtotElts with zero padding to the right in each row.
  // Use nchan to ensure we only load the required elements into Matrix.

  vector<vector<Real> > elements;
  SPreadVectorCol(rmf, "MATRIX", elements);
  if ( elements.size() == 0 ) {
    string msg = "Failed to read any entries from the MATRIX column";
    SPreportError(NoMatrix, msg);
    return(NoMatrix);
  }
  
  Matrix.resize(NtotElts);
  ipt = 0;
  for (size_t i=0; i<(size_t)Nrows; i++) {
    Integer NtoRead(0);
    for (size_t j=0; j<nchan[i].size(); j++) NtoRead += nchan[i][j];
    for (size_t j=0; j<(size_t)NtoRead; j++) {
      Matrix[ipt++] = elements[i][j];
    }
  }

  SPreadColUnits(rmf, "MATRIX", RMFUnits);

  // Read the optional order information

  vector<vector<Integer> > order;
  try {
    FITS::setVerboseMode(false);
    SPreadVectorCol(rmf, "ORDER", order);
  } catch (...) {
  }

  std::cout << "error stack:" << SPgetErrorStack() << std::endl;

  FITS::clearErrors();
  FITS::setVerboseMode(verbosity);

  if ( order.size() > 0 ) {
    OrderGroup.resize(NtotGroups);
    ipt = 0;
    for (size_t i=0; i<(size_t)Nrows; i++) {
      for (size_t j=0; j<(size_t)NumberGroups[i]; j++) {
	OrderGroup[ipt++] = order[j][i];
      }
    }
  }

  return(OK);
}

// update the FirstGroup and FirstElement arrays from NumberGroups and
// NumberChannelsGroup, respectively.

void rmf::update()
{
  // initialize
  FirstGroup.resize(NumberGroups.size());
  FirstElement.resize(NumberChannelsGroup.size());

  FirstGroup[0] = 0;
  for (size_t i=1; i<NumberGroups.size(); i++) {
    FirstGroup[i] = FirstGroup[i-1] + NumberGroups[i-1];
  }

  FirstElement[0];
  for (size_t i=1; i<NumberChannelsGroup.size(); i++) {
    FirstElement[i] = FirstElement[i-1] + NumberChannelsGroup[i-1];
  }

  return;
}

// initialize from an arf object. Copies members in common between arfs and rmfs

void rmf::initialize(const arf& a)
{
  EnergyUnits = a.EnergyUnits;
  Telescope   = a.Telescope;
  Instrument  = a.Instrument;
  Detector    = a.Detector;
  Filter      = a.Filter;
  
  LowEnergy.resize(a.LowEnergy.size());
  for (size_t i=0; i<LowEnergy.size(); i++) LowEnergy[i] = a.LowEnergy[i];
  HighEnergy.resize(a.HighEnergy.size());
  for (size_t i=0; i<HighEnergy.size(); i++) HighEnergy[i] = a.HighEnergy[i];

  return;
}

// Deep copy

rmf& rmf::operator= (const rmf& beta)
{

  FirstChannel = beta.FirstChannel;
  AreaScaling = beta.AreaScaling;
  ResponseThreshold = beta.ResponseThreshold;
  ChannelType = beta.ChannelType;
  RMFVersion = beta.RMFVersion;
  EBDVersion = beta.EBDVersion;
  Telescope = beta.Telescope;
  Instrument = beta.Instrument;
  Detector = beta.Detector;
  Filter = beta.Filter;
  RMFType = beta.RMFType;
  RMFExtensionName = beta.RMFExtensionName;
  EBDExtensionName = beta.EBDExtensionName;
  EnergyUnits = beta.EnergyUnits;
  RMFUnits = beta.RMFUnits;

  NumberGroups.resize(beta.NumberGroups.size());
  for (size_t i=0; i<NumberGroups.size(); i++) NumberGroups[i] = beta.NumberGroups[i];
  FirstGroup.resize(beta.FirstGroup.size());
  for (size_t i=0; i<FirstGroup.size(); i++) FirstGroup[i] = beta.FirstGroup[i];

  FirstChannelGroup.resize(beta.FirstChannelGroup.size());
  for (size_t i=0; i<FirstChannelGroup.size(); i++) FirstChannelGroup[i] = beta.FirstChannelGroup[i];
  NumberChannelsGroup.resize(beta.NumberChannelsGroup.size());
  for (size_t i=0; i<NumberChannelsGroup.size(); i++) NumberChannelsGroup[i] = beta.NumberChannelsGroup[i];
  FirstElement.resize(beta.FirstElement.size());
  for (size_t i=0; i<FirstElement.size(); i++) FirstElement[i] = beta.FirstElement[i];
  OrderGroup.resize(beta.OrderGroup.size());
  for (size_t i=0; i<OrderGroup.size(); i++) OrderGroup[i] = beta.OrderGroup[i];
  

  LowEnergy.resize(beta.LowEnergy.size());
  for (size_t i=0; i<LowEnergy.size(); i++) LowEnergy[i] = beta.LowEnergy[i];
  HighEnergy.resize(beta.HighEnergy.size());
  for (size_t i=0; i<HighEnergy.size(); i++) HighEnergy[i] = beta.HighEnergy[i];
  

  Matrix.resize(beta.Matrix.size());
  for (size_t i=0; i<Matrix.size(); i++) Matrix[i] = beta.Matrix[i];
  

  ChannelLowEnergy.resize(beta.ChannelLowEnergy.size());
  for (size_t i=0; i<ChannelLowEnergy.size(); i++) ChannelLowEnergy[i] = beta.ChannelLowEnergy[i];
  ChannelHighEnergy.resize(beta.ChannelHighEnergy.size());
  for (size_t i=0; i<ChannelHighEnergy.size(); i++) ChannelHighEnergy[i] = beta.ChannelHighEnergy[i];
  
  return *this;
}

// Return information

Integer rmf::NumberChannels()               // Number of spectrum channels 
{
  return ChannelLowEnergy.size();
}

Integer rmf::NumberEnergyBins()             // Number of response energies 
{
  return LowEnergy.size();
}

Integer rmf::NumberTotalGroups()            // Total number of response groups 
{
  return FirstChannelGroup.size();
}

Integer rmf::NumberTotalElements()          // Total number of response elements 
{
  return Matrix.size();
}

// Return the value for a particular channel and energy

Real rmf::ElementValue(Integer ChannelNumber, Integer EnergyBin)    
{                                           

  Integer GratingOrder = -999;

  return this->ElementValue(ChannelNumber, EnergyBin, GratingOrder);

}

// Return the value for a particular channel, energy and order. Use special case of
// GratingOrder = -999 for ignore order.

Real rmf::ElementValue(Integer ChannelNumber, Integer EnergyBin, Integer GratingOrder)
{                                          

  if ( ChannelNumber < FirstChannel || ChannelNumber >= FirstChannel+NumberChannels() || 
       EnergyBin < 0 || EnergyBin >= NumberEnergyBins() ) return 0.0;

  // loop round groups for this energy bin

  for(size_t i=FirstGroup[EnergyBin];i<(size_t)(FirstGroup[EnergyBin]+NumberGroups[EnergyBin]);i++) {

    if ( ( OrderGroup.size() > 0 && OrderGroup[i] == GratingOrder ) 
	 || GratingOrder == -999 ) {

      if( ChannelNumber >= FirstChannelGroup[i] && 
	  ChannelNumber < FirstChannelGroup[i]+NumberChannelsGroup[i]) {

	return(Matrix[FirstElement[i]+ChannelNumber-FirstChannelGroup[i]]);

      }

    }

  }

  return(0.0);

}

// Return vector of matrix values for a particular energy  

vector<Real> rmf::RowValues(Integer EnergyBin)
{

  Integer GratingOrder = -999;

  return this->RowValues(EnergyBin, GratingOrder);

}

// Return vector of matrix values for a particular energy and grating order
// Use GratingOrder = -999 as a special case for no check.

vector<Real> rmf::RowValues(Integer EnergyBin, Integer GratingOrder)
{
  vector<Real> values(NumberChannels());

  for (size_t i=0; i<values.size(); i++) values[i] = 0.0;

  // Loop round response groups for this energy

  for (size_t i=0; i<(size_t)NumberGroups[EnergyBin]; i++) {

    size_t igroup = i + FirstGroup[EnergyBin];
    size_t ivec = FirstChannelGroup[igroup] - FirstChannel;
    size_t ielt = FirstElement[igroup];

    // loop round elements in this group - adding them to the output array

    if ( ( OrderGroup.size() > 0 && OrderGroup[igroup] == GratingOrder ) 
	 || GratingOrder == -999 ) {

      for (size_t j=0; j<(size_t)NumberChannelsGroup[igroup]; j++) values[ivec+j] += Matrix[ielt+j];

    }

  }

  return values;
}

// Return vector of randomly generated channel numbers for a particular energy  

vector<Integer> rmf::RandomChannels(const Real energy, const Integer NumberPhotons)
{
  Integer GratingOrder(-999);

  return this->RandomChannels(energy, NumberPhotons, GratingOrder);
}

// Return vector of randomly generated channel numbers for a particular energy and
// grating order. Use GratingOrder = -999 as special case to ignore grating information

vector<Integer> rmf::RandomChannels(const Real energy, const Integer NumberPhotons, const Integer GratingOrder)
{
  vector<Integer> channel(NumberPhotons);

  // initialize the output array to -1s in the event that either the input energy is
  // outside the response range or that the response does not sum to unity and events
  // can fall off the end of the channels.

  for (size_t i=0; i<(size_t)NumberPhotons; i++) channel[i] = -1;

  size_t lower = 0;
  size_t upper = LowEnergy.size() - 1;

  // trap the case of the energy being outside the response range

  if ( energy < LowEnergy[lower] || energy > HighEnergy[upper] ) return channel;

  // find the energy bin associated with the input energy - assumes the energies are in increasing order

  size_t middle, energybin;
  while ( upper - lower > 1 ) {
    middle = (upper + lower)/2;
    if ( energy < HighEnergy[middle] ) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
  if ( energy > HighEnergy[lower] ) {
    energybin = upper;
  } else {
    energybin = lower;
  }

  // generate an array of size channel each element of which is the integrated response up to and
  // including that channel

  vector<Real> sumresponse(ChannelHighEnergy.size());
  for (size_t j=0; j<ChannelHighEnergy.size(); j++) sumresponse[j] = 0.0;

  sumresponse = this->RowValues(energybin, GratingOrder);

  for (size_t i=1; i<sumresponse.size(); i++) sumresponse[i] += sumresponse[i-1];

  // generate random numbers between 0 and 1

  vector<Real> RandomNumber(NumberPhotons);
  for (size_t i=0; i<(size_t)NumberPhotons; i++) RandomNumber[i] = (Real) HDmtDrand();

  // loop round the photons

  for (size_t i=0; i<(size_t)NumberPhotons; i++) {

    // find the array element containing this random number. note that we do
    // not assume that the total response sums to 1 - if the random number
    // exceeds the total response then we assume that the event fell off the
    // end of the channel array and return a -1

    lower = 0;
    upper = ChannelHighEnergy.size() - 1;
    if ( RandomNumber[i] <= sumresponse[upper] ) {
      while ( upper - lower > 1 ) {
	middle = (upper + lower)/2;
	if ( RandomNumber[i] < sumresponse[middle] ) {
	  upper = middle;
	} else {
	  lower = middle;
	}
      }
      if ( RandomNumber[i] > sumresponse[lower] ) {
	channel[i] = upper;
      } else {
	channel[i] = lower;
      }

      // correct the channel number for the first channel number in use in the response matrix

      channel[i] += FirstChannel;

    }

  }

  return channel;
}

// Display information about the RMF - return as a string

string rmf::disp()
{
  ostringstream outstr;

  outstr << "Response information : " << endl;

  outstr << "   FirstChannel        = " << FirstChannel << endl;
  outstr << "   AreaScaling         = " << AreaScaling << endl;
  outstr << "   ResponseThreshold   = " << ResponseThreshold << endl;
  outstr << "   ChannelType         = " << ChannelType << endl;
  outstr << "   RMFVersion          = " << RMFVersion << endl;
  outstr << "   EBDVersion          = " << EBDVersion << endl;
  outstr << "   Telescope           = " << Telescope << endl;
  outstr << "   Instrument          = " << Instrument << endl;
  outstr << "   Detector            = " << Detector << endl;
  outstr << "   Filter              = " << Filter << endl;
  outstr << "   RMFType             = " << RMFType << endl;
  outstr << "   RMFExtensionName    = " << RMFExtensionName << endl;
  outstr << "   EBDExtensionName    = " << EBDExtensionName << endl;

  outstr << "   EnergyUnits         = " << EnergyUnits << endl;
  outstr << "   RMFUnits            = " << RMFUnits << endl;

  outstr << "   NumberChannels      = " << NumberChannels() << endl;
  outstr << "   NumberEnergyBins    = " << NumberEnergyBins() << endl;
  outstr << "   NumberTotalGroups   = " << NumberTotalGroups() << endl;
  outstr << "   NumberTotalElements = " << NumberTotalElements() << endl;


  if ( NumberGroups.size() > 1 ) outstr << "   NumberGroups array of size " << NumberGroups.size() << endl;
  if ( FirstGroup.size() > 1 ) outstr << "   FirstGroup array of size " << FirstGroup.size() << endl;

  if ( FirstChannelGroup.size() > 1 ) outstr << "   FirstChannelGroup array of size " << FirstChannelGroup.size() << endl;
  if ( NumberChannelsGroup.size() > 1 ) outstr << "   NumberChannelsGroup array of size " << NumberChannelsGroup.size() << endl;
  if ( FirstElement.size() > 1 ) outstr << "   FirstElement array of size " << FirstElement.size() << endl;

  if ( OrderGroup.size() > 1 ) outstr << "   OrderGroup array of size " << OrderGroup.size() << endl;

  if ( LowEnergy.size() > 1 ) outstr << "   LowEnergy array of size " << LowEnergy.size() << endl;
  if ( HighEnergy.size() > 1 ) outstr << "   HighEnergy array of size " << HighEnergy.size() << endl;

  if ( Matrix.size() > 1 ) outstr << "   Matrix array of size " << Matrix.size() << endl;

  if ( ChannelLowEnergy.size() > 1 ) outstr << "   ChannelLowEnergy array of size " << ChannelLowEnergy.size() << endl;
  if ( ChannelHighEnergy.size() > 1 ) outstr << "   ChannelHighEnergy array of size " << ChannelHighEnergy.size() << endl;

  // debug
  //
  //  for (size_t i=0; i<NumberGroups.size(); i++) {
  //    outstr << i << "  " << NumberGroups[i];
  //    for (size_t j=0; j<NumberGroups[i]; j++) {
  //      outstr << "  " << FirstChannelGroup[j+FirstGroup[i]] << "  " << NumberChannelsGroup[j+FirstGroup[i]];
  //    }
  //    outstr << endl;
  //  }

  return outstr.str();
}

// Clear information from the response

void rmf::clear()
{
  FirstChannel = 0;

  NumberGroups.clear();
  FirstGroup.clear();
  FirstChannelGroup.clear();
  NumberChannelsGroup.clear();
  FirstElement.clear();
  OrderGroup.clear();

  LowEnergy.clear();
  HighEnergy.clear();
  Matrix.clear();
  ChannelLowEnergy.clear();
  ChannelHighEnergy.clear();

  AreaScaling = 0.0;
  ResponseThreshold = 0.0;

  EnergyUnits = " ";
  RMFUnits = " ";

  ChannelType = " ";
  RMFVersion = " ";
  EBDVersion = " ";
  Telescope = " ";
  Instrument = " ";
  Detector = " ";
  Filter = " ";
  RMFType = " ";
  RMFExtensionName = " ";
  EBDExtensionName = " ";

  return;
}

// Clear only the matrix from the response

void rmf::clearMatrix()
{
  NumberGroups.clear();
  FirstGroup.clear();
  FirstChannelGroup.clear();
  NumberChannelsGroup.clear();
  FirstElement.clear();
  OrderGroup.clear();
  LowEnergy.clear();
  HighEnergy.clear();
  Matrix.clear();

  return;
}

// Check completeness and consistency of information in the rmf
  // if there is a problem then return diagnostic in string

string rmf::check()
{
  ostringstream outstr;

  // check for presence of any data

  if ( Matrix.size() == 0 ) {
    outstr << "Matrix has no data" << endl;
  }

  // check size consistency between arrays - channels

  if ( ChannelLowEnergy.size() != ChannelHighEnergy.size() ) {
    outstr << "ChannelLowEnergy size (" << ChannelLowEnergy.size() 
	 << ") differs from ChannelHighEnergy size (" << ChannelHighEnergy.size() 
	 << ")" << endl;
  }

  // energies

  if ( LowEnergy.size() != HighEnergy.size() ) {
    outstr << "LowEnergy size (" << LowEnergy.size() 
	 << ") differs from HighEnergy size (" << HighEnergy.size() 
	 << ")" << endl;
  }

  if ( LowEnergy.size() != NumberGroups.size() ) {
    outstr << "LowEnergy size (" << LowEnergy.size() 
	 << ") differs from NumberGroups size (" << NumberGroups.size() 
	 << ")" << endl;
  }

  if ( LowEnergy.size() != FirstGroup.size() ) {
    outstr << "LowEnergy size (" << LowEnergy.size() 
	 << ") differs from FirstGroup size (" << FirstGroup.size() 
	 << ")" << endl;
  }

  // groups

  if ( FirstChannelGroup.size() != NumberChannelsGroup.size() ) {
    outstr << "FirstChannelGroup size (" << FirstChannelGroup.size() 
	 << ") differs from NumberChannelsGroup size (" << NumberChannelsGroup.size() 
	 << ")" << endl;
  }

  if ( FirstChannelGroup.size() != FirstElement.size() ) {
    outstr << "FirstChannelGroup size (" << FirstChannelGroup.size() 
	 << ") differs from FirstElement size (" << FirstElement.size() 
	 << ")" << endl;
  }

  if ( FirstChannelGroup.size() != OrderGroup.size() && OrderGroup.size() > 0 ) {
    outstr << "FirstChannelGroup size (" << FirstChannelGroup.size() 
	 << ") differs from OrderGroup size (" << OrderGroup.size() 
	 << ")" << endl;
  }

  // check that arrays have sensible values

  for (size_t i=0; i<NumberGroups.size(); i++) {
    if ( NumberGroups[i] < 0 ) {
      outstr << "NumberGroups has invalid value (" << NumberGroups[i] 
	   << ") for element " << i << endl;
    }
    if ( FirstGroup[i] < 0 || FirstGroup[i] >= (Integer)NumberChannelsGroup.size() ) {
      outstr << "FirstGroup has invalid value (" << FirstGroup[i] 
	   << ") for element " << i << endl;
    }
  }

  for (size_t i=0; i<FirstChannelGroup.size(); i++) {
    if ( FirstChannelGroup[i] < FirstChannel || FirstChannelGroup[i] >= (Integer)ChannelLowEnergy.size() ) {
      outstr << "FirstChannelGroup has invalid value (" << FirstChannelGroup[i] 
	   << ") for element " << i << endl;
    }
    if ( NumberChannelsGroup[i] < 0 || NumberChannelsGroup[i] >= (Integer)ChannelLowEnergy.size() ) {
      outstr << "NumberChannelsGroup has invalid value (" << NumberChannelsGroup[i] 
	   << ") for element " << i << endl;
    }
    if ( FirstElement[i] < 0 ) {
      outstr << "FirstElement has invalid value (" << FirstElement[i] 
	   << ") for element " << i << endl;
    }
  }      

  return outstr.str();
}

// Normalize the rmf so it sums to 1.0 for each energy bin

void rmf::normalize()
{

  // Loop over energies

  for (size_t ie=0; ie<(size_t)NumberEnergyBins(); ie++) {

    // sum up the response in this energy

    Real sumresp = 0.0;

    for (size_t i=0; i<(size_t)NumberGroups[ie]; i++) {
       size_t igrp = i + FirstGroup[ie];
       for (size_t j=0; j<(size_t)NumberChannelsGroup[igrp]; j++) {
         sumresp += Matrix[j+FirstElement[igrp]];
       }
    }

    // divide through by the summed response

    for (size_t i=0; i<(size_t)NumberGroups[ie]; i++) {
       size_t igrp = i + FirstGroup[ie];
       for (size_t j=0; j<(size_t)NumberChannelsGroup[igrp]; j++) {
         Matrix[j+FirstElement[igrp]] /= sumresp;
       }
    }

  }

  return;
}

void rmf::compress(const Real threshold)
{

  // Set up temporary object to store the output RMF and set it's threshold and first channel

  rmf work;
  work.ResponseThreshold = threshold;
  work.FirstChannel = FirstChannel;

  // Temporary array for the response for a given energy

  vector<Real> Response(ChannelLowEnergy.size());

  // Loop over energies

  for (size_t i=0; i<(size_t)LowEnergy.size(); i++) {

    // expand response matrix into a channel array for this energy

    Response = RowValues(i);

    // and add into the temporary response

    work.addRow(Response, LowEnergy[i], HighEnergy[i]);

    // end loop over energies

  }

  // copy new response arrays into current response

  NumberGroups.resize(work.NumberGroups.size());
  FirstChannelGroup.resize(work.FirstChannelGroup.size());
  NumberChannelsGroup.resize(work.NumberChannelsGroup.size());
  Matrix.resize(work.Matrix.size());

  for (size_t i=0; i<NumberGroups.size(); i++) {
    NumberGroups[i] = work.NumberGroups[i];
  }

  for (size_t i=0; i<FirstChannelGroup.size(); i++) {
    FirstChannelGroup[i] = work.FirstChannelGroup[i];
    NumberChannelsGroup[i] = work.NumberChannelsGroup[i];
  }

  for (size_t i=0; i<Matrix.size(); i++) {
    Matrix[i] = work.Matrix[i];
  }

  // update the FirstGroup and FirstElement arrays

  this->update();

  // set RMF threshold to the new value

  ResponseThreshold = threshold;

  return;
}

// Compress in channel space

Integer rmf::rebinChannels(grouping& GroupInfo)
{

  // check for consistency between grouping and number of channels

  if ( GroupInfo.size() != NumberChannels() ) {
    stringstream msg;
    msg << "Number of channels (" << NumberChannels() << ") differs from size of grouping array (" << GroupInfo.size() << ").";
    SPreportError(InconsistentGrouping, msg.str());
    return(InconsistentGrouping);
  }

  // Set up temporary object to store the output rmf and set its threshold and first channel

  rmf work;
  work.ResponseThreshold = ResponseThreshold;
  work.FirstChannel = FirstChannel;

  // Temporary array for the response for a given energy

  vector<Real> Response(NumberChannels());

  // Loop over energies

  for (size_t i=0; i<(size_t)LowEnergy.size(); i++) {

    // expand response matrix into a channel array for this energy

    Response = this->RowValues(i);

    // bin up the response for the energy

    vector<Real> binResponse;
    GroupBin(Response, SumMode, GroupInfo, binResponse);

    // update the compressed response arrays for this energy bin

    work.addRow(binResponse, LowEnergy[i], HighEnergy[i]);

    // end loop over energies

  }

  // copy new response arrays into current response

  NumberGroups.resize(work.NumberGroups.size());
  FirstChannelGroup.resize(work.FirstChannelGroup.size());
  NumberChannelsGroup.resize(work.NumberChannelsGroup.size());
  Matrix.resize(work.Matrix.size());

  for (size_t i=0; i<NumberGroups.size(); i++) {
    NumberGroups[i] = work.NumberGroups[i];
  }
  for (size_t i=0; i<NumberChannelsGroup.size(); i++) {
    FirstChannelGroup[i] = work.FirstChannelGroup[i];
    NumberChannelsGroup[i] = work.NumberChannelsGroup[i];
  }
  for (size_t i=0; i<Matrix.size(); i++) {
    Matrix[i] = work.Matrix[i];
  }

  // reset the FirstGroup and FirstElement arrays

  this->update();

  // now rebin the channel boundary arrays

  vector<Real> temp;
  GroupBin(ChannelLowEnergy, FirstEltMode, GroupInfo, temp);
  ChannelLowEnergy.resize(temp.size());
  ChannelLowEnergy = temp;
  GroupBin(ChannelHighEnergy, LastEltMode, GroupInfo, temp);
  ChannelHighEnergy.resize(temp.size());
  ChannelHighEnergy = temp;

  return(OK);
}

// Compress in energy space

Integer rmf::rebinEnergies(grouping& GroupInfo)
{

  // check for consistency between grouping and number of energy bins

  if ( GroupInfo.size() != NumberEnergyBins() ) {
    stringstream msg;
    msg << "Number of energy bins (" << NumberEnergyBins() << ") differs from size of grouping array (" << GroupInfo.size() << ").";
    SPreportError(InconsistentGrouping, msg.str());
    return(InconsistentGrouping);
  }

  // Set up temporary object to store the output rmf and set its threshold and first channel

  rmf work;
  work.ResponseThreshold = ResponseThreshold;
  work.FirstChannel = FirstChannel;

  // Temporary array for the response for a given energy

  vector<Real> Response(NumberChannels());
  Real eLow(0.0);

  // Loop over energies

  for (size_t i=0; i<(size_t)NumberEnergyBins(); i++) {

    // reset the Response array if this is the start of a new energy bin

    if ( GroupInfo.newBin(i) ) {
      for (size_t j=0; j<Response.size(); j++) Response[j] = 0.0;
      eLow = LowEnergy[i];
    }

    // expand response matrix into a channel array for this energy and accumulate

    vector<Real> RespArray(NumberChannels());
    RespArray = this->RowValues(i);
    for (size_t j=0; j<Response.size(); j++) Response[j] += RespArray[i];

    // if this is the last energy bin or the next one starts a new grouped bin then
    // update the compressed response arrays

    if ( i==(size_t)(NumberEnergyBins()-1) || GroupInfo.newBin(i+1) ) {
      work.addRow(Response, eLow, HighEnergy[i]);
    }

    // end loop over energy bins

  }

  // copy new response arrays into current response

  NumberGroups.resize(work.NumberGroups.size());
  FirstChannelGroup.resize(work.FirstChannelGroup.size());
  NumberChannelsGroup.resize(work.NumberChannelsGroup.size());
  Matrix.resize(work.Matrix.size());
  LowEnergy.resize(work.LowEnergy.size());
  HighEnergy.resize(work.HighEnergy.size());


  for (size_t i=0; i<NumberGroups.size(); i++) {
    NumberGroups[i] = work.NumberGroups[i];
  }
  for (size_t i=0; i<FirstChannelGroup.size(); i++) {
    FirstChannelGroup[i] = work.FirstChannelGroup[i];
  }
  for (size_t i=0; i<NumberChannelsGroup.size(); i++) {
    NumberChannelsGroup[i] = NumberChannelsGroup[i];
  }
  for (size_t i=0; i<Matrix.size(); i++) {
    Matrix[i] = work.Matrix[i];
  }
  for (size_t i=0; i<LowEnergy.size(); i++) {
    LowEnergy[i] = work.LowEnergy[i];
  }
  for (size_t i=0; i<HighEnergy.size(); i++) {
    HighEnergy[i] = work.HighEnergy[i];
  }

  // reset the FirstGroup and FirstElement arrays

  this->update();

  return(OK);
}

// Shift channels up or down.

Integer rmf::shiftChannels(Integer Start, Integer End, Real Shift)
{

  // switch Start and End to zero-based

  Integer Start0 = Start - FirstChannel;
  Integer End0 = End - FirstChannel;

  Integer N(NumberChannels());

  if ( End0 < 0 || Start0 > N-1 ) return(OK);

  // Set up temporary object to store the output rmf and set its threshold and first channel

  rmf work;
  work.ResponseThreshold = ResponseThreshold;
  work.FirstChannel = FirstChannel;

  // Temporary array for the response for a given energy

  vector<Real> Response(NumberChannels());

  Integer IntShift = (Integer)(abs(Shift)+1);
  Real FracShift = abs(Shift) - (IntShift-1);


  // Loop over energies

  for (size_t i=0; i<(size_t)LowEnergy.size(); i++) {

    // expand response matrix into a channel array for this energy

    Response = this->RowValues(i);

    vector<Real> OutResponse(N, 0.0);

    // shift the response for the energy

    Integer Ch = Start0;

    Integer Step = 1;
    if ( Shift < 0.0 ) Step = -1;
      
    while ( Ch <= End0 ) {
      Integer NewCh = Ch + Step*IntShift;
      if ( NewCh >= 0 && NewCh < N ) {
	OutResponse[NewCh] += FracShift * Response[Ch];
      }
      if ( NewCh-Step >= 0 && NewCh-Step < N ) {
	OutResponse[NewCh-Step] += (1-FracShift) * Response[Ch];
      }
      Ch++;
    }

    // update the compressed response arrays for this energy bin

    work.addRow(OutResponse, LowEnergy[i], HighEnergy[i]);

    // end loop over energies

  }

  // copy new response arrays into current response

  NumberGroups.resize(work.NumberGroups.size());
  FirstChannelGroup.resize(work.FirstChannelGroup.size());
  NumberChannelsGroup.resize(work.NumberChannelsGroup.size());
  Matrix.resize(work.Matrix.size());

  for (size_t i=0; i<NumberGroups.size(); i++) {
    NumberGroups[i] = work.NumberGroups[i];
  }
  for (size_t i=0; i<NumberChannelsGroup.size(); i++) {
    FirstChannelGroup[i] = work.FirstChannelGroup[i];
    NumberChannelsGroup[i] = work.NumberChannelsGroup[i];
  }
  for (size_t i=0; i<Matrix.size(); i++) {
    Matrix[i] = work.Matrix[i];
  }

  // reset the FirstGroup and FirstElement arrays

  this->update();

  return(OK);
}

// Write response matrix and channel bounds extensions. If file already exists then append.

Integer rmf::write(string filename)
{
  Integer Status(OK);
  Status = this->writeChannelBounds(filename);
  if ( Status != OK ) return(Status);
  Status = this->writeMatrix(filename);
  return(Status);
}

// Write response matrix and channel bounds extensions. If file already exists then append. Copy keywords and extra extensions from another file.

Integer rmf::write(string filename, string copyfilename)
{
  return(this->write(filename, copyfilename, 1));
}

// Write response matrix and channel bounds extensions. If file already exists then 
// append. Copy keywords and extra extensions from another file. Uses HDUnumber 
// instances of matrix and channel bounds extensions in both files.

Integer rmf::write(string filename, string copyfilename, Integer HDUnumber)
{
  Integer Status(OK);

  // this routines write the Matrix and ChannelBounds extensions while copying
  // extra keywords in these extensions from copyfilename. They do not copy extra
  // extensions

  Status = this->writeMatrix(filename, copyfilename, HDUnumber);
  if ( Status != OK ) return(Status);
  Status = this->writeChannelBounds(filename, copyfilename, HDUnumber);
  if ( Status != OK ) return(Status);

  // now copy the extra extensions

  Status = SPcopyHDUs(copyfilename, filename);

  return(Status);
}

// Write response matrix extension. If file already exists appends.

Integer rmf::writeMatrix(string filename)
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
    string msg = "Failed to create "+filename+" for MATRIX extension";
    SPreportError(CannotCreateMatrixExt, msg);
    return(CannotCreateMatrixExt);
  }

  // calculate the maximum number of groups and elements per row

  Integer Nrows = NumberGroups.size();
  
  Integer MaxGroups=0;
  for (size_t i=0; i<NumberGroups.size(); i++) {
    if ( NumberGroups[i] > MaxGroups ) MaxGroups = NumberGroups[i];
  }

  Integer MaxElts=0;
  for (size_t i=0; i<(size_t)Nrows; i++) {
    Integer NumElts=0;
    for (size_t j=0; j<(size_t)NumberGroups[i]; j++) {
      NumElts += NumberChannelsGroup[j+FirstGroup[i]];
    }
    if ( NumElts > MaxElts ) MaxElts = NumElts;
  }

  // set up the fchan and nchan vector arrays

  vector<vector<Integer> > fchan(Nrows), nchan(Nrows);
  for (size_t i=0; i<(size_t)Nrows; i++) {
    fchan[i].resize(MaxGroups);
    nchan[i].resize(MaxGroups);
    for (size_t j=0; j<(size_t)MaxGroups; j++) {
      fchan[i][j] = 0;
      nchan[i][j] = 0;
    }
    for (size_t j=0; j<(size_t)NumberGroups[i]; j++) {
      fchan[i][j] = FirstChannelGroup[j+FirstGroup[i]];
      nchan[i][j] = NumberChannelsGroup[j+FirstGroup[i]];
    }
  }

  // set up the array of matrix elements

  vector<vector<Real> > elements(Nrows);
  for (size_t i=0; i<(size_t)Nrows; i++) {
    Integer NumElts=0;
    for (size_t j=0; j<(size_t)NumberGroups[i]; j++) {
      NumElts += NumberChannelsGroup[j+FirstGroup[i]];
    }
    elements[i].resize(NumElts);
    for (size_t j=0; j<(size_t)NumElts; j++) {
      elements[i][j] = Matrix[j+FirstElement[FirstGroup[i]]];
    }
  }

  // if required set up the order vector array

  vector<vector<Integer> > order(Nrows);
  if ( OrderGroup.size() > 0 ) {

    for (size_t i=0; i<(size_t)Nrows; i++) { 
      order[i].resize(MaxGroups);
      for (size_t j=0; j<(size_t)MaxGroups; j++) {
	order[i][j] = 0;
      }
      for (size_t j=0; j<(size_t)NumberGroups[i]; j++) {
	order[i][j] = OrderGroup[j+FirstGroup[i]];
      }
    }

  }

  // set up the column descriptors for those attributes which need to be 
  // written as columns

  ttype.push_back("ENERG_LO");
  tform.push_back("E");
  tunit.push_back(EnergyUnits);

  ttype.push_back("ENERG_HI");
  tform.push_back("E");
  tunit.push_back(EnergyUnits);

  ttype.push_back("N_GRP");
  tform.push_back("I");
  tunit.push_back(" ");

  stringstream RepeatStream;
  RepeatStream << MaxGroups;
  string Repeat(RepeatStream.str());

  ttype.push_back("F_CHAN");
  if ( MaxGroups > 1 ) {
    tform.push_back(Repeat+"J");
  } else {
    tform.push_back("J");
  }
  tunit.push_back(" ");

  ttype.push_back("N_CHAN");
  if ( MaxGroups > 1 ) {
    tform.push_back(Repeat+"J");
  } else {
    tform.push_back("J");
  }
  tunit.push_back(" ");
     
  if ( SPneedCol(OrderGroup) ) {
    ttype.push_back("ORDER");
    if ( MaxGroups > 1 ) {
      tform.push_back(Repeat+"J");
    } else {
      tform.push_back("J");
    }
    tunit.push_back(" ");
  }

  RepeatStream.str("");
  RepeatStream << MaxElts;
  Repeat = RepeatStream.str();

  ttype.push_back("MATRIX");
  if ( MaxElts > 1 ) {
    tform.push_back("PE("+Repeat+")");
  } else {
    tform.push_back("E");
  }
  tunit.push_back(RMFUnits);

  // Create the new extension. First check for existing extensions to see whether
  // we need to give a version number.

  ExtMapConstIt itLow = pFits->extension().lower_bound("MATRIX");
  ExtMapConstIt itHigh = pFits->extension().upper_bound("MATRIX");

  Integer version(0);
  while (itLow != itHigh) {
    if (itLow->second->version() > version) version = itLow->second->version();
    itLow++;
  }

  version++;
  Table* prmf = pFits->addTable("MATRIX",Nrows,ttype,tform,tunit,BinaryTbl,version);
  Table& rmf = *prmf;

  // Write the standard keywords
  
  SPwriteKey(rmf, "HDUCLASS", (string)"OGIP", Blank);
    
  SPwriteKey(rmf, "HDUCLAS1", (string)"RESPONSE", Blank);

  SPwriteKey(rmf, "HDUCLAS2", (string)"RSP_MATRIX", Blank);
    
  SPwriteKey(rmf, "HDUCLAS3", RMFType, Blank);
    
  SPwriteKey(rmf, "CHANTYPE", ChannelType, "Channel type");

  SPwriteKey(rmf, "HDUVERS", RMFVersion, "OGIP version number");

  SPwriteKey(rmf, "TELESCOP", Telescope, Blank);

  SPwriteKey(rmf, "INSTRUME", Instrument, Blank);

  SPwriteKey(rmf, "DETNAM", Detector, Blank);

  SPwriteKey(rmf, "FILTER", Filter, Blank);

  SPwriteKey(rmf, "EFFAREA", AreaScaling, Blank);

  SPwriteKey(rmf, "LO_THRES", ResponseThreshold, Blank);

  SPwriteKey(rmf, "DETCHANS", NumberChannels(), "Number of channels in rmf");

  SPwriteKey(rmf, "NUMGRP", NumberTotalGroups(), "Total number of response groups");

  SPwriteKey(rmf, "NUMELT", NumberTotalElements(), "Total number of response elements");

  SPwriteKey(rmf, "TLMIN4", FirstChannel, "First channel number");

  // Write the arrays - if an array is of size 1 or all the same value 
  // it will be written as a keyword

  SPwriteCol(rmf, "ENERG_LO", LowEnergy, true);
  SPwriteCol(rmf, "ENERG_HI", HighEnergy, true);

  SPwriteCol(rmf, "N_GRP", NumberGroups, true);

  SPwriteVectorCol(rmf, "F_CHAN", fchan, true);
  SPwriteVectorCol(rmf, "N_CHAN", nchan, true);

  SPwriteVectorCol(rmf, "MATRIX", elements, true);

  if ( OrderGroup.size() > 0 ) {
    SPwriteVectorCol(rmf, "ORDER", order);
  }
  
  return(OK);
}

// Write response matrix extension. If file already exists appends. Copy extra
// keywords from the matrix extension in copyfilename. Note that unlike this
// method for some other classes extra extensions are not copied.

Integer rmf::writeMatrix(string filename, string copyfilename)
{
  return(this->writeMatrix(filename, copyfilename, 1));
}

// Write response matrix extension. If file already exists appends. Copy extra
// keywords from the HDUnumber instance of matrix extension in copyfilename. Note 
// that unlike this method for some other classes extra extensions are not copied.

Integer rmf::writeMatrix(string filename, string copyfilename, Integer HDUnumber)
{
  Integer Status(OK);

  Status = this->writeMatrix(filename);

  if ( Status != OK ) return(Status);

  Status = SPcopyKeys(copyfilename, filename, "MATRIX", HDUnumber);

  return(Status);
}


// Write channel bounds extension. If file already exists appends.

Integer rmf::writeChannelBounds(string filename)
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
    string msg = "Failed to create "+filename+" for EBOUNDS extension";
    SPreportError(CannotCreateEboundsExt, msg);
    return(CannotCreateEboundsExt);
  }

  // set up the column descriptors for those attributes which need to be 
  // written as columns

  ttype.push_back("CHANNEL");
  tform.push_back("I");
  tunit.push_back(" ");

  ttype.push_back("E_MIN");
  tform.push_back("E");
  tunit.push_back(EnergyUnits);

  ttype.push_back("E_MAX");
  tform.push_back("E");
  tunit.push_back(EnergyUnits);

  // Create the new extension

  Table* pebd = pFits->addTable("EBOUNDS",ChannelLowEnergy.size(),ttype,tform,tunit);
  Table& ebd = *pebd;

  // Write the standard keywords
  
  SPwriteKey(ebd, "HDUCLASS", (string)"OGIP", Blank);
    
  SPwriteKey(ebd, "HDUCLAS1", (string)"RESPONSE", Blank);

  SPwriteKey(ebd, "HDUCLAS2", (string)"EBOUNDS", Blank);
    
  SPwriteKey(ebd, "CHANTYPE", ChannelType, "Channel type");

  SPwriteKey(ebd, "HDUVERS", EBDVersion, "OGIP version number");

  SPwriteKey(ebd, "TELESCOP", Telescope, Blank);

  SPwriteKey(ebd, "INSTRUME", Instrument, Blank);

  SPwriteKey(ebd, "DETNAM", Detector, Blank);

  SPwriteKey(ebd, "FILTER", Filter, Blank);

  SPwriteKey(ebd, "DETCHANS", NumberChannels(), "Number of channels in ebd");

  // Generate and write the COLUMN array

  vector<Real> Channel(NumberChannels());
  for (size_t i=0; i<(size_t)NumberChannels(); i++) Channel[i] = i + FirstChannel;
  SPwriteCol(ebd, "CHANNEL", Channel);

  // Write the E_MIN and E_MAX arrays - if an array is of size 1 or all the 
  // same value it will be written as a keyword

  SPwriteCol(ebd, "E_MIN", ChannelLowEnergy);
  SPwriteCol(ebd, "E_MAX", ChannelHighEnergy);

  return(OK);
}

// Write channel bounds extension. If file already exists appends. Copy extra
// keywords from the channel bounds in copyfilename. Note that unlike this
// method for some other classes extra extensions are not copied.

Integer rmf::writeChannelBounds(string filename, string copyfilename)
{
  return(this->writeChannelBounds(filename, copyfilename, 1));
}

// Write channel bounds matrix extension. If file already exists appends. Copy extra
// keywords from the HDUnumber instance of channel bounds extension in copyfilename. Note 
// that unlike this method for some other classes extra extensions are not copied.

Integer rmf::writeChannelBounds(string filename, string copyfilename, Integer HDUnumber)
{
  Integer Status(OK);

  Status = this->writeChannelBounds(filename);
  if ( Status != OK ) return(Status);

  Status = SPcopyKeys(copyfilename, filename, "EBOUNDS", HDUnumber);

  return(Status);
}

// Merge arf and rmf

rmf& rmf::operator*=(const arf& a)
{
  // check that the arf and rmf are compatible
  // if not just return the current rmf

  if ( !checkCompatibility(a) ) return *this;

  // loop round energy bins multiplying appropriate elements of the rmf by
  // the effective area for this energy from the ARF

  for ( size_t i=0; i < (size_t)LowEnergy.size(); i++ ) {

    Real effarea = a.EffArea[i];

    // loop over response groups for this energy

    for ( size_t j=FirstGroup[i]; j<(size_t)(FirstGroup[i]+NumberGroups[i]); j++ ) {

      // loop over matrix elements for this response group

      for ( size_t k=FirstElement[j]; k<(size_t)(FirstElement[j]+NumberChannelsGroup[j]); k++ ) {

	Matrix[k] *= effarea;

      }

    }

  }     

  return *this;
}

rmf& rmf::operator+=(const rmf& r)
{

 // check that the two rmfs are compatible
 // if not just return the current rmf

  if ( !checkCompatibility(r) ) return *this;

  // temporary arrays for the response for each energy

  vector<Real> Response1(ChannelLowEnergy.size());
  vector<Real> Response2(ChannelLowEnergy.size());

  // temporary object to store the output rmf - set its threshold and first channel

  rmf work;
  work.ResponseThreshold = ResponseThreshold;
  work.FirstChannel = FirstChannel;

  // copy input rmf to get round r being defined as constant which produces
  // and error when using RowValues.

  rmf In(r);

  // loop round energy bins summing appropriate elements of the rmf

  for ( size_t i=0; i < (size_t)LowEnergy.size(); i++ ) {

    // expand both response matrices into a channel array for this energy

    Response1 = this->RowValues(i);
    Response2 = In.RowValues(i);

    // sum the two responses for this energy bin

    for (size_t j=0; j<Response1.size(); j++) {
      Response1[j] += Response2[j];
    }

    // update the compressed response arrays for this energy bin

    work.addRow(Response1, LowEnergy[i], HighEnergy[i]);

    // end loop over energies

  }

  // copy new response arrays into current response

  NumberGroups.resize(work.NumberGroups.size());
  FirstChannelGroup.resize(work.FirstChannelGroup.size());
  NumberChannelsGroup.resize(work.NumberChannelsGroup.size());
  Matrix.resize(work.Matrix.size());

  for (size_t i=0; i<NumberGroups.size(); i++) {
    NumberGroups[i] = work.NumberGroups[i];
  }
  for (size_t i=0; i<NumberChannelsGroup.size(); i++) {
    FirstChannelGroup[i] = work.FirstChannelGroup[i];
    NumberChannelsGroup[i] = work.NumberChannelsGroup[i];
  }
  for (size_t i=0; i<Matrix.size(); i++) {
    Matrix[i] = work.Matrix[i];
  }

  // reset the FirstGroup and FirstElement arrays

  this->update();

  return *this;
}


Integer rmf::checkCompatibility(const rmf& r)
{

  // check that the two rmf energy binnings are compatible

  if ( LowEnergy.size() != r.LowEnergy.size() ) return InconsistentEnergies;
  for ( size_t i=0; i < (size_t)LowEnergy.size(); i++ ) {
    if ( LowEnergy[i] != r.LowEnergy[i] ) return InconsistentEnergies;
    if ( HighEnergy[i] != r.HighEnergy[i] ) return InconsistentEnergies;
  }

 // check that the two rmf channel binnings are compatible

  if ( ChannelLowEnergy.size() != r.ChannelLowEnergy.size() ) return InconsistentChannels;
  for ( size_t i=0; i < (size_t)ChannelLowEnergy.size(); i++ ) {
    if ( ChannelLowEnergy[i] != r.ChannelLowEnergy[i] ) return InconsistentChannels;
    if ( ChannelHighEnergy[i] != r.ChannelHighEnergy[i] ) return InconsistentChannels;
  }

  if ( EnergyUnits != r.EnergyUnits ) return InconsistentUnits;
  if ( RMFUnits != r.RMFUnits ) return InconsistentUnits;

  return OK;
}

Integer rmf::checkCompatibility(const arf& a)
{

  if ( LowEnergy.size() != a.LowEnergy.size() ) return InconsistentEnergies;
  for ( size_t i=0; i < (size_t)LowEnergy.size(); i++ ) {
    if ( LowEnergy[i] != a.LowEnergy[i] ) return InconsistentEnergies;
    if ( HighEnergy[i] != a.HighEnergy[i] ) return InconsistentEnergies;
  }

  return OK;

}

// convert the energy units to keV if they are something else

Integer rmf::convertUnits()
{

  // set up energy/wave conversion factors and check for valid units

  bool xwave;
  Real xfactor;

  Integer status(OK);

  status = calcXfactor(EnergyUnits, xwave, xfactor);
  if ( status != OK ) return(status);

  if ( xfactor == 1.0 ) return(OK);

  if ( xwave ) {
    for (size_t i=0; i<LowEnergy.size(); i++) {
      Real temp(HighEnergy[i]);
      HighEnergy[i] = xfactor/LowEnergy[i];
      LowEnergy[i] = xfactor/temp;
    }
    for (size_t i=0; i<ChannelLowEnergy.size(); i++) {
      ChannelLowEnergy[i] = xfactor/ChannelLowEnergy[i];
      ChannelHighEnergy[i] = xfactor/ChannelHighEnergy[i];
    }
  } else {
    for (size_t i=0; i<LowEnergy.size(); i++) {
      LowEnergy[i] = xfactor*LowEnergy[i];
      HighEnergy[i] = xfactor*HighEnergy[i];
    }
    for (size_t i=0; i<ChannelLowEnergy.size(); i++) {
      ChannelLowEnergy[i] = xfactor*ChannelLowEnergy[i];
      ChannelHighEnergy[i] = xfactor*ChannelHighEnergy[i];
    }
  }

  // if necessary reverse channel low and high energies arrays

  for (size_t i=0; i<ChannelLowEnergy.size(); i++) {
    if ( ChannelLowEnergy[i] > ChannelHighEnergy[i] ) {
      Real temp(ChannelHighEnergy[i]);
      ChannelHighEnergy[i] = ChannelLowEnergy[i];
      ChannelLowEnergy[i] = temp;
    }
  }

  EnergyUnits = "keV";
  return(OK);
}

// reverse the rows, required if they are not in increasing order of energy
// note that this does not reverse the channels.
// note also that this does not the OrderGroup array

void rmf::reverseRows()
{

  // easiest just to set up a work rmf. Need to set the threshold and first channel in work rmf.

  rmf work;
  work.ResponseThreshold = ResponseThreshold;
  work.FirstChannel = FirstChannel;

  // loop through the response in reverse order extracting the response vector
  // and constructing new response

  size_t Nbins(LowEnergy.size());
  vector<Real> ResponseValues(NumberChannels());

  for (size_t i=0; i<Nbins; i++) {

    ResponseValues = this->RowValues(Nbins-i-1);

    work.addRow(ResponseValues, LowEnergy[Nbins-i-1], HighEnergy[Nbins-i-1]);

  }

  // Reset the arrays based on those in work

  for (size_t i=0; i<NumberChannelsGroup.size(); i++) {
    FirstChannelGroup[i] = work.FirstChannelGroup[i];
    NumberChannelsGroup[i] = work.NumberChannelsGroup[i];
  }

  for (size_t i=0; i<Matrix.size(); i++) {
    Matrix[i] = work.Matrix[i];
  }

  for (size_t i=0; i<Nbins; i++) {
    NumberGroups[i] = work.NumberGroups[i];
    LowEnergy[i] = work.LowEnergy[i];
    HighEnergy[i] = work.HighEnergy[i];
  }

  // Reset the FirstGroup and FirstElement arrays

  this->update();

  return;

}


void rmf::addRow(const vector<Real> Response, const Real eLow, const Real eHigh) 
{

  Integer NGroups(0);
  bool inGroup(false);
  bool First(true);

  for ( size_t j=0; j<Response.size(); j++ ) {

    if ( Response[j] >= ResponseThreshold ) {

      // if not in a response group then start a new one

      if ( !inGroup ) {

	NGroups++;

	FirstChannelGroup.push_back(j+FirstChannel);
	NumberChannelsGroup.push_back(1);
	Matrix.push_back(Response[j]);
	FirstElement.push_back(Matrix.size()-1);

	inGroup = true;

	if ( First ) {
	  FirstGroup.push_back(NumberChannelsGroup.size()-1);
	  First = false;
	}

	// otherwise add next response to this group

      } else {

	NumberChannelsGroup[NumberChannelsGroup.size()-1]++;
	Matrix.push_back(Response[j]);

      }

      // if response below threshold then end group if it is open

    } else {

      if ( inGroup ) inGroup = false;

    }

    // end loop over channels

  }

  // add the number of groups and response energies

  NumberGroups.push_back(NGroups);
  LowEnergy.push_back(eLow);
  HighEnergy.push_back(eHigh);

  return;

}

// version of addRow for multiple grating orders

void rmf::addRow(const vector<vector<Real> > Response, const Real eLow, const Real eHigh,
		 const vector<Integer> GratingOrder) 
{

  Integer NGroups(0);
  bool First(true);

  // loop over grating orders

  for ( size_t i=0; i<Response.size(); i++ ) {

    bool inGroup(false);

    for ( size_t j=0; j<Response[i].size(); j++ ) {

      if ( Response[i][j] >= ResponseThreshold ) {

	// if not in a response group then start a new one

	if ( !inGroup ) {

	  NGroups++;

	  FirstChannelGroup.push_back(j+FirstChannel);
	  NumberChannelsGroup.push_back(1);
	  Matrix.push_back(Response[i][j]);
	  FirstElement.push_back(Matrix.size()-1);
	  OrderGroup.push_back(GratingOrder[i]);

	  inGroup = true;

	  if ( First ) {
	    FirstGroup.push_back(NumberChannelsGroup.size()-1);
	    First = false;
	  }

	  // otherwise add next response to this group

	} else {

	  NumberChannelsGroup[NumberChannelsGroup.size()-1]++;
	  Matrix.push_back(Response[i][j]);

	}

	// if response below threshold then end group if it is open

      } else {

	if ( inGroup ) inGroup = false;

      }

      // end loop over channels

    }

    // end loop over orders

  }

  // add the number of groups and response energies

  NumberGroups.push_back(NGroups);
  LowEnergy.push_back(eLow);
  HighEnergy.push_back(eHigh);

  return;

}

// utility routines that are not methods for the rmf object

rmf operator* (const rmf& r, const arf& a){
  rmf rr(r);
  return rr *= a;
}

rmf operator* (const arf& a, const rmf& r){
  rmf rr(r);
  return rr *= a;
}

rmf operator+ (const rmf& a, const rmf& b){
  rmf r(a);
  return r += b;
}

// calculate the response vector for some energy given a gaussian width
// the gaussian is assumed to be in the units of energy,
// ChannelLowEnergy and ChannelHighEnergy

void calcGaussResp(const Real sigma, const Real energyLow, const Real energyHigh, 
                   const Real threshold, const vector<Real>& ChannelLowEnergy, 
		   const vector<Real>& ChannelHighEnergy, vector<Real>& ResponseVector)
{

  Real energy = (energyLow + energyHigh) / 2.0;
  Real winv = 1.0/sigma/sqrt(2.0);
  size_t N = ChannelLowEnergy.size();
  ResponseVector.resize(N);

  // store the part of the gaussian below the response energy range in alow

  Real alow = 0.5 * erf(winv*(ChannelLowEnergy[0]-energy));
  Real ahi(0.0);

  // loop round channels
  // This can be sped up a lot by starting with the centroid bin and
  // working outwards. Simple method: loop over all bins.

  for (size_t i=0; i<N; i++) {

    ahi = 0.5 * erf(winv*(ChannelHighEnergy[i]-energy));
    ResponseVector[i] = ahi - alow;
    alow = ahi;

  }

  return;
}

