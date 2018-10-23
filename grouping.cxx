// grouping class. Useful for setting grouping arrays and binning
// both spectra and responses.

#ifndef HAVE_grouping
#include "grouping.h"
#endif

// default constructor

grouping::grouping()
  : flag()
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
  // clear vector with guaranteed reallocation
  vector<Integer>().swap(flag);
}

// display grouping information - return as a string

string grouping::disp()
{
  ostringstream outstr;

  if ( flag.size() > 0 ) {
    outstr << "Grouping information : " << endl;
    size_t nblock = 25;
    for (size_t i=0; i<(flag.size()/nblock); i++) {
      outstr << setw(5) << i*nblock << " - " << setw(5) << (i+1)*nblock-1 << "    ";
      for (size_t j=0; j<nblock; j++) outstr << flag[j+i*nblock] << " ";
      outstr << std::endl;
    }
    size_t nleft = flag.size() - nblock*(flag.size()/nblock);
    outstr << setw(5) << flag.size()-nleft << " - " << setw(5) << flag.size()-1 << "    ";
    for (size_t j=0; j<nleft; j++) outstr << flag[flag.size()-nleft+j] << " ";
    outstr << std::endl;
  } else {
    outstr << "No grouping information set yet" << endl;
  }

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

  StartBin.push_back(0);
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

      for (size_t ich=(size_t)StartBin0[ir]; ich<=(size_t)EndBin0[ir]; ich++) flag[ich] = 0;

    } else {

      for (size_t ich=(size_t)StartBin0[ir]; ich<=(size_t)EndBin0[ir]; ich+=(size_t)BinFactor[ir]) {

	flag[ich] = 1;
	size_t nbins = BinFactor[ir];
	if ( nbins > flag.size()-ich ) nbins = flag.size() - ich;
	for (size_t ibin=1; ibin<nbins; ibin++) {
	  flag[ich+ibin] = -1;
	}

      }

    }

  }

  /* set grouping for any elements after the end of binning */

  for (size_t ich=EndBin0[EndBin0.size()-1]+1; ich<flag.size(); ich++) flag[ich] = 1;
  return(OK);
}

// set from quality and grouping vectors from a pha object

Integer grouping::loadFromVector(const vector<Integer>& QualVector, const vector<Integer>& GroupVector)
{
  size_t Nchan(QualVector.size());
  if ( GroupVector.size() != Nchan ) {
    stringstream msg;
    msg << "Quality vector length (" << Nchan 
	<< ") not equal to Grouping vector length (" << GroupVector.size() 
	<< ").";
    SPreportError(InconsistentGrouping, msg.str());
    return(InconsistentGrouping);
  }

  // loop round channels setting the flag vector
  flag.resize(Nchan);
  for (size_t i=0; i<Nchan; i++) {
    if ( QualVector[i] == 0 ) {
      flag[i] = GroupVector[i];
    } else {
      flag[i] = 0;
    }
  }

  return(OK);
}


// set from a minimum based on a vector input

template <class T> Integer grouping::loadMin(const T Minimum, const vector<T>& Values)
{
  Integer StartChannel = 0;
  Integer EndChannel = Values.size()-1;
  return loadMin(Minimum, StartChannel, EndChannel, Values);
}

template <class T> Integer grouping::loadMin(const T Minimum, const Integer StartChannel,
			  const Integer EndChannel, const vector<T>& Values)
{
  // loop over channels accumulating up to Minimum. Any channels left over at
  // the end are ignored. Note that this implicitly assumes that StartChannel
  // and EndChannel are zero-based

  size_t Nchan = Values.size();
  size_t First(0), Last(Nchan-1);
  if ( StartChannel > 0 ) First = StartChannel;
  if ( EndChannel < Nchan-1 ) Last = EndChannel;
  bool ingroup = false;
  size_t istart = 0;
  T sum = 0;
  flag.resize(Nchan);
  for (size_t i=0; i<First; i++) flag[i] = 1;
  for (size_t i=First; i<=Last; i++) {
    sum += Values[i];
    if ( !ingroup ) {
      flag[i] = 1;
      istart = i;
      sum = Values[i];
      if ( sum < Minimum ) ingroup = true;
    } else if ( sum >= Minimum ) {
      flag[i] = -1;
      ingroup = false;
    } else if ( sum < Minimum ) {
      flag[i] = -1;
    }
  }
  if ( sum < Minimum && ingroup ) {
    for (size_t i=istart; i<=Last; i++) flag[i] = 0;
  }
  for (size_t i=Last+1; i<Nchan; i++) flag[i] = 1;

  return 0;
}

// required to make the linker instantiate correctly

template Integer grouping::loadMin(const Integer Minimum, const vector<Integer>& Values);
template Integer grouping::loadMin(const Real Minimum, const vector<Real>& Values);
template Integer grouping::loadMin(const Integer Minimum, const Integer StartChannel,
				   const Integer EndChannel, const vector<Integer>& Values);
template Integer grouping::loadMin(const Real Minimum, const Integer StartChannel,
				   const Integer EndChannel, const vector<Real>& Values);


// set grouping using optimal binning based on instrument FWHM and number of
// counts per resolution element.

Integer grouping::loadOptimal(const vector<Real>& FWHM, const vector<Integer>& Counts)
{

  // This performs optimal binning following Kaastra and Bleeker 2016, A&A 587, A151
  // The binning delta/FWHM is defined as
  //    1  if x <= 2.119 otherwise
  //    (0.08 + 7.0/x + 1.8/x^2)/(1 + 5.9/x)
  // where
  //    x = ln[N_r(1 + 0.2 ln R)]
  // and N_r is the number of counts per resolution element and R is the total
  // number of resolution elements.

  size_t Nchan = FWHM.size();
  if ( Nchan != Counts.size() ) return(InconsistentGrouping);

  // First estimate the total number of resolution elements. Since this enters
  // as ln(0.2ln R) it can be a rough estimate. See eqn B.1 of K&B.

  Real R(1);
  for (size_t i=0; i<Nchan; i++) R += 1.0/FWHM[i];
  Real logR = log(R);

  // Calculate the number of counts within the FWHM for each channel then
  // multiply by 1.314 to get the number of counts per resolution element.
  // This assumes the response is gaussian - I could improve this by also
  // including a vector with the fraction within the FWHM for each channel
  // but this is probably not going to make a significant difference

  vector<Real> Nr(Nchan);
  for (size_t i=0; i<Nchan; i++) {
    size_t start(0), end(Nchan-1);
    Real low = i - FWHM[i]/2.0;
    Real high = i + FWHM[i]/2.0;
    if ( round(low) >= 0 ) start = (size_t)round(low);
    if ( round(high) < Nchan ) end = (size_t)round(high);
    Nr[i] = 0.0;
    for (size_t k=start; k<=end; k++) Nr[i] += Counts[k];
    Nr[i] *= 1.314;
  }

  // Calculate the optimal bin size at each channel

  vector<Real> b(FWHM);
  for (size_t i=0; i<Nchan; i++) {
    Real x = log(Nr[i]*(1.0+0.2*logR));
    if ( x > 2.119 ) b[i] *= (0.08*x + 7.0 + 1.8/x) / (x + 5.9);
  }
  vector<Integer> bint(Nchan);
  for (size_t i=0; i<Nchan; i++) {
    bint[i] = (int)b[i];
    if ( bint[i] < 1 ) bint[i] = 1;
  }

  // Now set the grouping

  size_t ichan(0);
  flag.resize(Nchan);
  while ( ichan < Nchan ) {
    flag[ichan] = 1;
    size_t end = std::min(Nchan-1,ichan + bint[ichan] - 1);
    size_t a(end);
    for (size_t k=ichan+1; k<=end; k++) {
      if ( k + bint[k] -1 < a ) a = k + bint[k] - 1;
    }
    end = a;
    for (size_t k=ichan+1; k<=end; k++) flag[k] = -1;
    ichan = end + 1;
  }

  return(OK);
}

// set using optimal energy binning based on the instrument FWHM.
// note that this assumes that FWHM is that for each energy

Integer grouping::loadOptimalEnergy(const vector<Real>& FWHM, const vector<Integer>& Counts)
{
  // This performs optimal binning on the energy axis following Kaastra & Bleeker
  // 2016, A&A 587, A151. The binning delta/FWHM is defined as the minimum of 1 and
  // y where
  //     y = (1.404/x^{1/4})(1 + 18/x)
  //     x = N_r (1 + 0.1 ln R)
  // and N_r is the number of counts per resolution element and R is the total
  // number of resolution elements.

  size_t nE = FWHM.size();
  size_t nChan = Counts.size();
  
  // First estimate the total number of resolution elements. Since this enters
  // as 0.1lnR it can be a rough estimate

  Real R(1);
  for (size_t i=0; i<nChan; i++) R += 1.0/FWHM[i];
  Real logR = log(R);

  // Calculate the number of counts within the FWHM for each energy then
  // multiply by 1.314 to get the number of counts per resolution element.
  // This assumes the response is gaussian - I could improve this by also
  // including a vector with the fraction within the FWHM for each channel
  // but this is probably not going to make a significant difference

  vector<Real> Nr(nE);
  for (size_t i=0; i<nE; i++) {
    size_t start(0), end(nE-1);
    Real low = i - FWHM[i]/2.0;
    Real high = i + FWHM[i]/2.0;
    if ( round(low) >= 0 ) start = (size_t)round(low);
    if ( round(high) < nE ) end = (size_t)round(high);
    Nr[i] = 0.0;
    for (size_t k=start; k<=end; k++) Nr[i] += Counts[k];
    Nr[i] *= 1.314;
  }

  // Calculate the optimal bin size at each channel

  vector<Real> b(FWHM);
  for (size_t i=0; i<nE; i++) {
    Real x = Nr[i]*(1.0+0.1*logR);
    Real y = (1.404/(sqrt(sqrt(x))))*(1.0+18.0/x);
    if ( y > 1.0 ) b[i] *= y;
  }
  vector<Integer> bint(nE);
  for (size_t i=0; i<nE; i++) {
    bint[i] = (int)b[i];
    if ( bint[i] < 1 ) bint[i] = 1;
  }

  // Now set the grouping

  size_t iEn(0);
  flag.resize(nE);
  while ( iEn < nE ) {
    flag[iEn] = 1;
    size_t end = std::min(nE-1,iEn + bint[iEn] - 1);
    size_t a(end);
    for (size_t k=iEn+1; k<=end; k++) {
      if ( k + bint[k] -1 < a ) a = k + bint[k] - 1;
    }
    end = a;
    for (size_t k=iEn+1; k<=end; k++) flag[k] = -1;
    iEn = end + 1;
  }

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

// a set of channels starting with grouping=1 for the first and grouping=-1
// for the next are binned together. Any channels with grouping=0 are jumped
// over since they are assumed to be bad.

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

  outArray.clear();

  for (size_t ich=0; ich<inArray.size(); ich++) {

    // if the input bin is the start of an output bin

    if (GroupInfo.flag[ich] == 1) {

      // if necessary modify current output bin
      if ( !first && mode == MeanMode ) outArray[outpt] /= ncount;
      if ( !first && mode == SumQuadMode ) outArray[outpt] = (T)sqrt((Real)outArray[outpt]);

      // start a new output bin
      // note that the FirstEltMode case is handled silently because
      // if !SumQuadMode then outArray is set to current inArray and
      // if FirstEltMode then outArray is not modified if the GroupInfo.flag
      // is zero.
      first = false;
      ncount = 1;
      if ( mode == SumQuadMode ) {
	outArray.push_back(inArray[ich]*inArray[ich]);
      } else {
	outArray.push_back(inArray[ich]);
      }
      outpt = outArray.size() - 1;

    // else if input bin is not the start of an output bin

    } else if (GroupInfo.flag[ich] == -1) {
 
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

