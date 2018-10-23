// class definition for grouping class. Useful for setting grouping arrays and binning
// both spectra and responses.

#ifndef HAVE_HEASP
#include "heasp.h"
#endif

#ifndef HAVE_SPutils
#include "SPutils.h"
#endif

#define HAVE_grouping 1

class grouping{
 public:

  vector<Integer> flag;              // Grouping flag: 1=start of bin, -1=continuation of bin, 0=ignore

  // constructor

  grouping();

  // constructor from an input array

  grouping(vector<Integer>);

  // destructor

  ~grouping();

  // display grouping information - return as a string

  string disp();

  // clear grouping information

  void clear();

  // read from an ascii file of grouping factors

  Integer read(string filename, const Integer Number, const Integer First);

  // set from a single binning factor

  void load(const Integer Binsize, const Integer Number);

  // set from an array of binning factors

  Integer load(const vector<Integer>& StartBin, const vector<Integer>& EndBin, const vector<Integer>& BinFactor, const Integer Number, const Integer First);

  // set from quality and grouping vectors from a pha object

  Integer loadFromVector(const vector<Integer>& QualVector, const vector<Integer>& GroupVector);

  // set from a minimum with optional start and stop channel for the grouping

  template <class T> Integer loadMin(const T Minimum, const vector<T>& Values); 
  template <class T> Integer loadMin(const T Minimum, const Integer StartChannel, 
		  const Integer EndChannel, const vector<T>& Values);

  // set using optimal channel binning based on the instrument FWHM.
  // note that this assumes that FWHM is that for each channel, not for each energy

  Integer loadOptimal(const vector<Real>& FWHM, const vector<Integer>& Counts);

  // set using optimal energy binning based on the instrument FWHM.
  // note that this assumes that FWHM is that for each energy

  Integer loadOptimalEnergy(const vector<Real>& FWHM, const vector<Integer>& Counts);

  // return whether current element is start of new bin

  bool newBin(const Integer i);

  // return number of elements in grouping object

  Integer size();

};

// definition of the binning modes

enum{SumMode,SumQuadMode,MeanMode,FirstEltMode,LastEltMode};

// bin an array based on the grouping factors

template <class T> void GroupBin(const vector<T>&, const Integer, const grouping&, vector<T>&);
template <class T> void GroupBin(const valarray<T>&, const Integer, const grouping&, valarray<T>&);

// read a file with binning factors

Integer ReadBinFactors(string filename, vector<Integer>& StartBin, vector<Integer>& EndBin, vector<Integer>& BinFactor); 
