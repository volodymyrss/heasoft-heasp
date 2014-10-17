// grouping class. Useful for setting grouping arrays and binning
// both spectra and responses.

#ifndef HAVE_grouping
#include "grouping.h"
#endif

// default constructor

grouping::grouping()
{
}

// constructor from an input array

grouping::grouping(vector<Integer> beta)
{
  flag.resize(beta.size());
  for(size_t i=0; i<flag.size(); i++) flag[i] = beta[i];
}

// destructor

grouping::~grouping()
{
}

// display grouping information - return as a string

string grouping::disp()
{
  ostringstream outstr;

  outstr << "Grouping information : " << endl;
  size_t nblock = 25;
  for (size_t i=0; i<(flag.size()/nblock); i++) {
    outstr << setw(5) << i*nblock+1 << " - " << setw(5) << (i+1)*nblock << "    ";
    for (size_t j=0; j<nblock; j++) outstr << flag[j+i*nblock+1] << " ";
    outstr << std::endl;
  }
  size_t nleft = flag.size() - nblock*(flag.size()/nblock);
  outstr << setw(5) << flag.size()-nleft+1 << " - " << setw(5) << flag.size() << "    ";
  for (size_t j=0; j<nleft; j++) outstr << flag[flag.size()-nleft+1+j] << " ";
  outstr << std::endl;

  return outstr.str();
}

// clear grouping information

void grouping::clear()
{
  flag.clear();
  return;
}

// read from an ascii file of grouping factors
// each line of the file consists of three numbers, the start, end, and binning factor
// Note that bin numbers are assumed to start at the value given by First..

Integer grouping::read(string filename, const Integer Number, const Integer First)
{

  Integer Status(OK);
  vector<Integer> StartBin, EndBin, BinFactor;

  Status = ReadBinFactors(filename, StartBin, EndBin, BinFactor);
  if ( Status != OK ) {
    string msg = "Failed to read "+filename;
    SPreportError(Status, msg);
    return(Status);
  }

  return(this->load(StartBin, EndBin, BinFactor, Number, First));
}

// set from a single binning factor

void grouping::load(const Integer BinSize, const Integer Number)
{

  vector<Integer> StartBin, EndBin, BinFactor;

  StartBin.push_back(1);
  EndBin.push_back(Number-1);
  BinFactor.push_back(BinSize);

  this->load(StartBin, EndBin, BinFactor, Number, 0);

  return;
}

// set from an array of binning factors

Integer grouping::load(const vector<Integer>& StartBin, const vector<Integer>& EndBin, 
		       const vector<Integer>& BinFactor, const Integer Number, const Integer First)
{

  // check that the input grouping won't run off the front or back of the flag array.

  for (size_t i=0; i<StartBin.size(); i++) {
    if ( StartBin[i] < First || EndBin[i] > First+Number-1 || BinFactor[i] <= 0 ) {
      stringstream msg;
      msg << "Start, End and binning factor = " << StartBin[i] << " " << EndBin[i] 
	  << " " << BinFactor[i] << ", input first channel = " << First
	  << " and number of channels = " << Number;
      SPreportError(InconsistentGrouping, msg.str());
      return(InconsistentGrouping);
    }
  }

  // create new StartBin and EndBin arrays which count from 0.

  vector<Integer> StartBin0(StartBin);
  vector<Integer> EndBin0(EndBin);

  for (size_t i=0; i<StartBin.size(); i++) {
    StartBin0[i] = StartBin[i] - First;
    EndBin0[i] = EndBin[i] - First;
  }
  
  // convert to the grouping flag array which is assumed to be Number long

  flag.resize(Number);

  // first set grouping for any elements before the start of binning

  for (size_t ich=0; ich<(size_t)StartBin0[0]; ich++) flag[ich] = 1;

  // loop through ranges for binning

  for (size_t ir=0; ir<BinFactor.size(); ir++) {

    if ( BinFactor[ir] == -1 ) {

      for (size_t ich=(size_t)StartBin0[ir]; ich<=(size_t)EndBin0[ir]; ich++) flag[ich] = -1;

    } else {

      for (size_t ich=(size_t)StartBin0[ir]; ich<=(size_t)EndBin0[ir]; ich+=(size_t)BinFactor[ir]) {

	flag[ich] = 1;
	size_t nbins = BinFactor[ir];
	if ( nbins > flag.size()-ich ) nbins = flag.size() - ich;
	for (size_t ibin=1; ibin<nbins; ibin++) {
	  flag[ich+ibin] = 0;
	}

      }

    }

  }

  /* set grouping for any elements after the end of binning */

  for (size_t ich=EndBin0[EndBin0.size()-1]+1; ich<=flag.size(); ich++) flag[ich] = 1;
  return(OK);
}

// return whether the current element is that start of a bin

bool grouping::newBin(const Integer i)
{
  if ( flag[i] == 1 ) return true;
  return false;
}

// return the number of elements in the grouping object

Integer grouping::size()
{
  return flag.size();
}

// bin an array based on the grouping factors
//   mode = SumMode       Sum
//   mode = SumQuadMode   Sum in quadrature
//   mode = MeanMode      Mean
//   mode = FirstEltMode  First ie value is that of first channel in bin
//   mode = LastEltMode   Last  ie value is that of last channel in bin */


template <class T> void GroupBin(const valarray<T>& inArray, const Integer mode, const grouping& GroupInfo, valarray<T>& outArray)
{
  vector<T> inTemp, outTemp;

  inTemp.resize(inArray.size());
  for (size_t ich=0; ich<inTemp.size(); ich++) inTemp[ich] = inArray[ich];

  GroupBin(inTemp, mode, GroupInfo, outTemp);

  outArray.resize(outTemp.size());
  for (size_t ich=0; ich<outArray.size(); ich++) outArray[ich] = outTemp[ich];
  
  return;
}

template <class T> void GroupBin(const vector<T>& inArray, const Integer mode, const grouping& GroupInfo, vector<T>& outArray)
{
  Integer ncount(1), outpt(0);
  bool first(true);

  for (size_t ich=0; ich<inArray.size(); ich++) {

    // if the input bin is the start of an output bin

    if (GroupInfo.flag[ich] == 1) {

      // if necessary modify current output bin
      if ( !first && mode == MeanMode ) outArray[outpt] /= ncount;
      if ( !first && mode == SumQuadMode ) outArray[outpt] = (T)sqrt((Real)outArray[outpt]);

      // start a new output bin
      first = false;
      ncount = 1;
      if ( mode == SumQuadMode ) {
	outArray.push_back(inArray[ich]*inArray[ich]);
      } else {
	outArray.push_back(inArray[ich]);
      }
      outpt = outArray.size() - 1;

    // else if input bin is not the start of an output bin

    } else if (GroupInfo.flag[ich] == 0) {
 
      if ( mode == SumMode || mode == MeanMode ) outArray[outpt] += inArray[ich];
      if ( mode == SumQuadMode ) outArray[outpt] += inArray[ich]*inArray[ich];
      if ( mode == LastEltMode ) outArray[outpt] = inArray[ich];
      ncount++;

    }

  }

  if ( ncount > 1 ) {
    if ( mode == MeanMode ) outArray[outpt] /= ncount;
    if ( mode == SumQuadMode ) outArray[outpt] = (T)sqrt((Real)outArray[outpt]);
  }

  return;

}

// required to make the linker instantiate correctly

template void GroupBin(const vector<Real>&, const Integer, const grouping&, vector<Real>&);
template void GroupBin(const vector<Integer>&, const Integer, const grouping&, vector<Integer>&);
template void GroupBin(const valarray<Real>&, const Integer, const grouping&, valarray<Real>&);
template void GroupBin(const valarray<Integer>&, const Integer, const grouping&, valarray<Integer>&);

// read a file with binning factors

Integer ReadBinFactors(string filename, vector<Integer>& StartBin, vector<Integer>& EndBin, vector<Integer>& BinFactor)
{

  // read the file setting up arrays of start, end, and binning factors. 

  ifstream infile;
  try {
    infile.open(filename.c_str(), ifstream::in);
  } catch(...) {
    string msg = "Failed to open "+filename;
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  Integer Start, End, Factor;

  if ( infile.is_open() ) {
    string instring;
    getline(infile, instring);
    while ( !infile.eof() ) {
      stringstream instream;
      instream << instring;
      instream >> Start >> End >> Factor;
      StartBin.push_back(Start);
      EndBin.push_back(End);
      BinFactor.push_back(Factor);
      getline(infile, instring);
    }
    infile.close();
  } else {
    string msg = "Failed to open "+filename;
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  return(OK);
}

