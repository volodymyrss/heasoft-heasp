// SWIG interface for heasp library. Based on heasp.h, pha.h, phaII.h, rmf.h, rmft.h, 
// arf.h, arfII.h, grouping.h table.h

%module heasp

%{
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <stdexcept>
#include <ctime>
#include <valarray>
#include <vector>

#include <CCfits/CCfits>

#include "heasp.h"
#include "pha.h"
#include "phaII.h"
#include "rmf.h"
#include "rmft.h"
#ifndef HAVE_arf
  #include "arf.h"
#endif
#include "arfII.h"
#ifndef HAVE_grouping
  #include "grouping.h"
#endif
#include "table.h"
%}

// use the std and CCfits namespaces

using namespace std;

// define Integer and Real

typedef int Integer;
typedef float Real;

// set up error statuses

enum{OK, NoSuchFile, NoData, NoChannelData, NoStatError, CannotCreate,
     NoEnergLo, NoEnergHi, NoSpecresp, NoEboundsExt, NoEmin, NoEmax,
     NoMatrixExt, NoNgrp, NoFchan, NoNchan, NoMatrix, CannotCreateMatrixExt,
     CannotCreateEboundsExt, InconsistentGrouping, InconsistentEnergies,
     InconsistentChannels, InconsistentUnits, UnknownXUnits, UnknownYUnits};

// SWIG includes and template definitions

%include "typemaps.i"
%include "std_string.i"
%include "std_vector.i"

%template(IntVector) vector<Integer>;
%template(floatVector) vector<Real>;
//%template(doubleVector) vector<Real>;

// *************************************************************************
// Class definitions for pha object

// python renames
%rename(__iadd__) pha::operator+=;
%rename(__imul__) pha::operator*=;

class pha{
 public:

  Integer FirstChannel;                 // First legal channel number

  vector<Real> Pha;                        // PHA data
  vector<Real> StatError;                  // Statistical error 
  vector<Real> SysError;                   // Statistical error 

  vector<Integer> Channel;                 // Channel number
  vector<Integer> Quality;                 // Data quality 
  vector<Integer> Group;                   // Data grouping 

  vector<Real> AreaScaling;                // Area scaling factor 
  vector<Real> BackScaling;                // Background scaling factor 

  Real Exposure;                        // Exposure time 
  Real CorrectionScaling;               // Correction file scale factor 

  Integer DetChans;                     // Total legal number of channels
  bool Poisserr;                        // If true, errors are Poisson 
  string Datatype;                      // "COUNT" for count data and "RATE" for count/sec 
  string PHAVersion;                    // PHA extension format version 

  string Spectrumtype;                  // "TOTAL", "NET", or "BKG" 

  string ResponseFile;                  // Response filename 
  string AncillaryFile;                 // Ancillary filename 
  string BackgroundFile;                // Background filename 
  string CorrectionFile;                // Correction filename 

  string FluxUnits;                     // Units for Pha and StatError

  string ChannelType;                   // Value of CHANTYPE keyword 
  string Telescope;                                          
  string Instrument;
  string Detector;
  string Filter;
  string Datamode;

  vector<string> XSPECFilter;           // Filter keywords 

  // constructor

  pha();

  // destructor

  ~pha();

  // read file into object. the third option is to read a pha object from a single row
  // of a type I file SpectrumNumber

  Integer read(string filename);
  Integer read(string filename, Integer PHAnumber);
  Integer read(string filename, Integer PHAnumber, Integer SpectrumNumber);

  // Deep copy

//  pha& operator= (const pha&);

  // Return information

  Integer NumberChannels();            // size of internal Arrays

  // Display information about the spectrum

  string disp();

  // Clear the information in the spectrum

  void clear();

  // Check completeness and consistency of information in spectrum

  string check();

  // Write spectrum as type I file

  Integer write(string filename);
  Integer write(string filename, string copyfilename);
  Integer write(string filename, string copyfilename, Integer HDUnumber);

  // Multiply by a constant

  pha& operator*= (const Real);

  // Add to another pha

  pha& operator+= (const pha&);

  // Check compatibility with another pha

  Integer checkCompatibility(const pha&);

  // Set grouping array from grouping object

  Integer setGrouping(grouping&);

  // Rebin channels

  Integer rebinChannels(grouping&, string);
  Integer rebinChannels(grouping&);

  // Shift channels

  Integer shiftChannels(Integer Start, Integer End, Real Shift);

  // Convert flux units from whatever they are currently to ph/cm^2/s. 
  // This requires an input the channel energy arrays from the rmf object and
  // the string specifying their units.

  Integer convertUnits(vector<Real>& ChannelLowEnergy, vector<Real>& ChannelHighEnergy, string EnergyUnits);

};

// Binary operation

// pha operator+ (const pha& a, const pha& b);

// Utility routines

// return the type of a PHA extension

Integer PHAtype(string filename, Integer PHAnumber); 

// return true if COUNTS column exists and is integer

bool IsPHAcounts(string filename, Integer PHAnumber); 

// return the number of spectra in a type II PHA extension

Integer NumberofSpectra(string filename, Integer PHAnumber); 

// Definition of the SpectrumII object. Just a wrap-up for a vector array of Spectrum objects

// *************************************************************************
// Class definitions for phaII object

class phaII{
 public:

  vector<pha> phas;           // vector of pha objects

  // constructor

  phaII();

  // destructor

  ~phaII();

  // read file into object. 

  Integer read(string filename);
  Integer read(string filename, Integer PHAnumber);
  Integer read(string filename, Integer PHAnumber, vector<Integer> SpectrumNumber);

  // Deep copy

//  phaII& operator= (const phaII&);

  // Get pha object (counts from zero).

  pha get(Integer number);

  // Push pha object into phaII object

  void push(pha spectrum);

  // Return information

  Integer NumberSpectra();          // Number of Spectra in the object

  // Display information about the spectra

  string disp();

  // Clear information about the spectra

  void clear();

  // Check completeness and consistency of information in spectra

  string check();

  // Write spectra as type II file

  Integer write(string filename);
  Integer write(string filename, string copyfilename);
  Integer write(string filename, string copyfilename, Integer HDUnumber);

};

// *************************************************************************
// Class definitions for rmf object

// python renames
%rename(__iadd__) rmf::operator+=;
%rename(__imul__) rmf::operator*=;


class rmf{
 public:

  Integer FirstChannel;              // First channel number 

  vector<Integer> NumberGroups;         // Number of response groups for this energy bin 
  vector<Integer> FirstGroup;           // First response group for this energy bin (counts from 0)

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

  // rmf& operator= (const rmf&);

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
  // of given energy  (and grating order).

  vector<Integer> RandomChannels(const Real energy, const Integer NumberPhotons);
  vector<Integer> RandomChannels(const Real energy, const Integer NumberPhotons, const Integer GratingOrder);

  // Display information about the spectrum

  string disp();

  // Clear information from the response

  void clear();

  // Clear only the matrix from the response

  void clearMatrix();

  // Check completeness and consistency of information in the rmf

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

//rmf operator* (const rmf&, const arf&);
//rmf operator* (const arf&, const rmf&);
//rmf operator+ (const rmf&, const rmf&);

// calculate the response vector for some energy given a gaussian width
// the gaussian is assumed to be in the units of energyLow, energyHigh,
// ChannelLowEnergy and ChannelHighEnergy

void calcGaussResp(const Real width, const Real energyLow, const Real energyHigh, 
                   const Real threshold, const vector<Real>& ChannelLowEnergy, 
		   const vector<Real>& ChannelHighEnergy, vector<Real>& ResponseVector);

// *************************************************************************
// Class definitions for ResponseMatrixTranspose object

class rmft{
 public:

  Integer FirstChannel;              // First channel number 

  vector<Integer> NumberGroups;         // Number of response groups for this channel bin 
  vector<Integer> FirstGroup;           // First response group for this channel bin (counts from 0)

  vector<Integer> FirstEnergyGroup;    // First energy bin in this group 
  vector<Integer> NumberEnergiesGroup; // Number of energy bins in this group 
  vector<Integer> FirstElement;        // First response element for this group (counts from 0)
  vector<Integer> OrderGroup;          // The grating order of this group 

  vector<Real> LowEnergy;              // Start energy of bin 
  vector<Real> HighEnergy;             // End energy of bin 

  vector<Real> Matrix;                 // Matrix elements 

  vector<Real> ChannelLowEnergy;       // Start energy of channel 
  vector<Real> ChannelHighEnergy;      // End energy of channel 

  Real AreaScaling;                 // Value of EFFAREA keyword 
  Real ResponseThreshold;           // Minimum value in response 

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

  rmft();

  // destructor

  ~rmft();

  // load object from a standard rmf

  void load(rmf&);

  // update the FirstGroup and FirstElement arrays from NumberGroups and
  // NumberEnergiesGroup, respectively.

  void update();

  // Deep copy

//  rmft& operator= (const rmft&);

  // Return information

  Integer NumberChannels();               // Number of spectrum channels 
  Integer NumberEnergyBins();             // Number of response energies 
  Integer NumberTotalGroups();            // Total number of response groups 
  Integer NumberTotalElements();          // Total number of response elements 

  Real ElementValue(Integer, Integer);    // Return the value for a particular channel
                                          // and energy

  vector<Real> RowValues(Integer);          // Return the array for a particular channel

  // Display information about the object

  string disp();

  // Clear information from the object

  void clear();

};

// *************************************************************************
// Class definitions for arf object

// python renames
%rename(__iadd__) arf::operator+=;

class arf{
 public:

  vector<Real> LowEnergy;                   // Start energy of bin
  vector<Real> HighEnergy;                  // End energy of bin

  vector<Real> EffArea;                     // Effective areas

  string EnergyUnits;                     // Units for energies
  string arfUnits;                        // Units for effective areas

  string Version;                        // SPECRESP extension format version
  string Telescope;                             
  string Instrument;
  string Detector;
  string Filter;
  string ExtensionName;               // Value of EXTNAME keyword in SPECRESP extension

  // constructor

  arf();

  // destructor

  ~arf();

  // read file into object. Third option is to read from a row of a type II file

  Integer read(string filename);
  Integer read(string filename, Integer ARFnumber);
  Integer read(string filename, Integer ARFnumber, Integer RowNumber);

  // Deep copy

//  arf& operator= (const arf&);

  // Return information

  Integer NumberEnergyBins();            // size of vector<Real>s

  // Display information about the arf - return as a string

  string disp();

  // Clear information from the arf

  void clear();

  // Check completeness and consistency of information in the arf
  // if there is a problem then return diagnostic in string

  string check();

  // Write arf

  Integer write(string filename);
  Integer write(string filename, string copyfilename);
  Integer write(string filename, string copyfilename, Integer HDUnumber);

  // Add arfs

  arf& operator+=(const arf&);

  Integer checkCompatibility(const arf&);

  Integer convertUnits();

};

// define this outside the class

//arf operator+ (const arf&, const arf&);

// *************************************************************************
// Definition of the arfII object. Just a wrap-up for a vector array of arf objects

class arfII{
 public:

  vector<arf> arfs;           // vector of arf objects

  // constructor

  arfII();

  // destructor

  ~arfII();

  // read file into object. 

  Integer read(string filename);
  Integer read(string filename, Integer ARFnumber);
  Integer read(string filename, Integer ARFnumber, vector<Integer> RowNumber);

  // Deep copy

//  arfII& operator= (const arfII&);

  // Get arf object (counts from zero).

  arf get(Integer number);

  // Push arf object into arfII object

  void push(arf ea);

  // Return information

  Integer NumberARFs();          // Number of ARFs in the object

  // Display information about the ARFs

  string disp();

  // Clear information from the ARFs

  void clear();

  // Check completeness and consistency of information in the arfs

  string check();

  // Write ARFs as type II file

  Integer write(string filename);
  Integer write(string filename, string copyfilename);
  Integer write(string filename, string copyfilename, Integer HDUnumber);

};

// *************************************************************************
// class definition for grouping class. Useful for setting grouping arrays and binning
// both spectra and responses.

class grouping{
 public:

  vector<Integer> flag;              // Grouping flag: 1=start of bin, 0=continuation of bin

  // constructor

  grouping();

  // constructor from an input array

  grouping(vector<Integer>);

  // destructor

  ~grouping();

  // display grouping information

  string disp();

  // clear grouping information

  void clear();

  // read from an ascii file of grouping factors

  Integer read(string, const Integer, const Integer);

  // set from a single binning factor

  void load(const Integer, const Integer);

  // set from an array of binning factors

  Integer load(const vector<Integer>& StartBin, const vector<Integer>& EndBin, const vector<Integer>& BinFactor, const Integer Number, const Integer First); 

  // return whether current element is start of new bin

  bool newBin(const Integer);

  // return number of elements in grouping object

  Integer size();

};

// definition of the binning modes

enum{SumMode,SumQuadMode,MeanMode,FirstEltMode,LastEltMode};

// bin an array based on the grouping factors

// template <class T> void GroupBin(const valarray<T>&, const Integer, const grouping&, valarray<T>&);
template <class T> void GroupBin(const vector<T>&, const Integer, const grouping&, vector<T>&);

// read a file with binning factors

Integer ReadBinFactors(string filename, vector<Integer>& StartBin, vector<Integer>& EndBin, vector<Integer>& BinFactor); 

// *************************************************************************
// class definition for individual parameters within the table model

class tableParameter{
 public:

  string Name;                  // Parameter name
  int InterpolationMethod;      // 0==linear, 1==log, -1==additional (non-interp_
  Real InitialValue;            // Initial value for fit
  Real Delta;                   // Delta for fit
  Real Minimum;                 // Hard lower-limit (should correspond to first tabulated value)
  Real Bottom;                  // Soft lower-limit
  Real Top;                     // Soft upper-limit
  Real Maximum;                 // Hard upper-limit (should correspond to last tabulated value)
  vector<Real> TabulatedValues; // Tabulated parameter values

  //constructor

  tableParameter();

  // destructor

  ~tableParameter();

  // display information about the table parameter

  string disp();

  // clear contents of the table parameter (mainly useful for Python)

  void clear();
};

// class definition for individual spectra (and additional spectra) within the table model

class tableSpectrum{
 public:

  vector<Real> Flux;
  vector<Real> ParameterValues;
  vector<vector<Real> > addFlux;

  //constructor

  tableSpectrum();

  // destructor

  ~tableSpectrum();

  // push an additional parameter spectrum

  void pushaddFlux(vector<Real>);

  // get an additional parameter spectrum

  vector<Real> getaddFlux(Integer Number);

  // display information about the table spectrum

  string disp();

  // clear contents of the table parameter (mainly useful for Python)

  void clear();

};


// class definition for table

class table{
 public:

  vector <tableParameter> Parameters;
  vector <tableSpectrum> Spectra;
  string ModelName;
  string ModelUnits;
  int NumIntParams;
  int NumAddParams;
  bool isError;
  bool isRedshift;
  bool isAdditive;
  vector<Real> Energies;
  string EnergyUnits;

  // constructor

  table();

  // destructor

  ~table();

  // read the table from a FITS file

  Integer read(string filename);

  // Push table Parameter object

  void pushParameter(const tableParameter& paramObject);

  // Push spectrum Parameter object

  void pushSpectrum(const tableSpectrum& spectrumObject);

  // Get table Parameter object (counts from zero)

  tableParameter getParameter(Integer number);

  // extract table Spectrum object (counts from zero)

  tableSpectrum getSpectrum(Integer number);

  // display information about the table

  string disp();

  // clear contents of table object (mainly useful for Python)

  void clear();

  // Check completeness and consistency of information in the table

  string check();

  // convert to standard units (keV and ph/cm^2/s).

  Integer convertUnits();

  // reverse the rows if energies are not increasing (as can occur eg after
  // converting from wavelength)

  void reverseRows();

  // write to a FITS file

  Integer write(string filename);

};


