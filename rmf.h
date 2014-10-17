// Class definitions for rmf object

#ifndef HAVE_HEASP
#include "heasp.h"
#endif

#ifndef HAVE_SPutils
#include "SPutils.h"
#endif

#ifndef HAVE_arf
#include "arf.h"
#endif

#ifndef HAVE_grouping
#include "grouping.h"
#endif

#include "headas_rand.h"

#define HAVE_rmf 1

class rmf{
 public:

  Integer FirstChannel;              // First channel number 

  vector<Integer> NumberGroups;        // Number of response groups for this energy bin 
  vector<Integer> FirstGroup;          // First response group for this energy bin (counts from 0)

  vector<Integer> FirstChannelGroup;   // First channel number in this group 
  vector<Integer> NumberChannelsGroup; // Number of channels in this group 
  vector<Integer> FirstElement;        // First response element for this group (counts from 0)
  vector<Integer> OrderGroup;          // The grating order of this group 

  vector<Real> LowEnergy;              // Start energy of bin 
  vector<Real> HighEnergy;             // End energy of bin 

  vector<Real> Matrix;                 // Matrix elements 

  vector<Real> ChannelLowEnergy;       // Start energy of channel 
  vector<Real> ChannelHighEnergy;      // End energy of channel 

  Real AreaScaling;                 // Value of EFFAREA keyword 
  Real ResponseThreshold;           // Minimum value in response 

  string EnergyUnits;               // Energy units used
  string RMFUnits;                  // Units for RMF values

  string ChannelType;               // Value of CHANTYPE keyword 
  string RMFVersion;                // MATRIX extension format version 
  string EBDVersion;                // EBOUNDS extension format version 
  string Telescope;                             
  string Instrument;
  string Detector;
  string Filter;
  string RMFType;                   // Value of HDUCLAS3 keyword in MATRIX extension 
  string RMFExtensionName;          // Value of EXTNAME keyword in MATRIX extension 
  string EBDExtensionName;          // Value of EXTNAME keyword in EBOUNDS extension 

  // constructor

  rmf();

  // destructor

  ~rmf();

  // read file into object. 

  Integer read(string filename);
  Integer read(string filename, Integer RMFnumber);
  Integer readMatrix(string filename);
  Integer readMatrix(string filename, Integer RMFnumber);
  Integer readChannelBounds(string filename);
  Integer readChannelBounds(string filename, Integer RMFnumber);

  // update the FirstGroup and FirstElement arrays from NumberGroups and
  // NumberChannelsGroup, respectively.

  void update();

  // initialize from an arf object. Copies members in common between arfs and rmfs

  void initialize(const arf&);

  // Deep copy

  rmf& operator= (const rmf&);

  // Return information

  Integer NumberChannels();               // Number of spectrum channels 
  Integer NumberEnergyBins();             // Number of response energies 
  Integer NumberTotalGroups();            // Total number of response groups 
  Integer NumberTotalElements();          // Total number of response elements 

  Real ElementValue(Integer, Integer);    // Return the value for a particular channel
                                          // and energy
  Real ElementValue(Integer, Integer, Integer);  // ... and grating order

  vector<Real> RowValues(Integer);        // Return the array for a particular energy
  vector<Real> RowValues(Integer, Integer); // ... and grating order

  // Use the response matrix to generate random channel numbers for a photon 
  // of given energy (and grating order).

  vector<Integer> RandomChannels(const Real energy, const Integer NumberPhotons);
  vector<Integer> RandomChannels(const Real energy, const Integer NumberPhotons, const Integer GratingOrder);

  // Display information about the response - return as a string

  string disp();

  // Clear information from the response

  void clear();

  // Clear only the matrix from the response

  void clearMatrix();

  // Check completeness and consistency of information in the rmf
  // if there is a problem then return diagnostic in string

  string check();

  // Normalize the rmf so it sums to 1.0 for each energy bin

  void normalize();

  // Compress the rmf to remove all elements below the threshold value

  void compress(const Real threshold);

  // Rebin in either channel or energy space

  Integer rebinChannels(grouping&);
  Integer rebinEnergies(grouping&);

  // Remaps response up or down in channels

  Integer shiftChannels(Integer Start, Integer Stop, Real Shift);

  // Write response

  Integer write(string filename);
  Integer write(string filename, string copyfilename);
  Integer write(string filename, string copyfilename, Integer HDUnumber);

  Integer writeMatrix(string filename);
  Integer writeMatrix(string filename, string copyfilename);
  Integer writeMatrix(string filename, string copyfilename, Integer HDUnumber);

  Integer writeChannelBounds(string filename);
  Integer writeChannelBounds(string filename, string copyfilename);
  Integer writeChannelBounds(string filename, string copyfilename, Integer HDUnumber);

  // Merge ARF and rmf

  rmf& operator*=(const arf&);

  // add rmf's

  rmf& operator+=(const rmf&);

  // check compatibility with another rmf or ARF

  Integer checkCompatibility(const rmf&);
  Integer checkCompatibility(const arf&);

  // convert units. mainly useful if input energy/wavelengths are not in keV

  Integer convertUnits();

  // reverse the rows, required if they are not in increasing order of energy

  void reverseRows();

  // add a row to the response using an input response vector and energy range.

  void addRow(const vector<Real> Response, const Real eLow, const Real eHigh);
  void addRow(const vector<vector<Real> > Response, const Real eLow, const Real eHigh, const vector<Integer> GratingOrder);

};

// define these outside the class

// binary operators 

rmf operator* (const rmf&, const arf&);
rmf operator* (const arf&, const rmf&);
rmf operator+ (const rmf&, const rmf&);

// calculate the response vector for some energy given a gaussian width
// the gaussian is assumed to be in the units of energyLow, energyHigh,
// ChannelLowEnergy and ChannelHighEnergy

void calcGaussResp(const Real width, const Real energyLow, const Real energyHigh, 
                   const Real threshold, const vector<Real>& ChannelLowEnergy, 
		   const vector<Real>& ChannelHighEnergy, vector<Real>& ResponseVector);
