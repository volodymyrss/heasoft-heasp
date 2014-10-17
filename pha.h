// Class definitions for PHA object

#ifndef HAVE_HEASP
#include "heasp.h"
#endif

#ifndef HAVE_SPutils
#include "SPutils.h"
#endif

#ifndef HAVE_grouping
#include "grouping.h"
#endif

#define HAVE_pha 1

class pha{
 public:

  Integer FirstChannel;                 // First legal channel number

  vector<Real> Pha;                     // PHA data
  vector<Real> StatError;               // Statistical error 
  vector<Real> SysError;                // Statistical error 

  vector<Integer> Channel;              // Channel number
  vector<Integer> Quality;              // Data quality 
  vector<Integer> Group;                // Data grouping 

  vector<Real> AreaScaling;             // Area scaling factor 
  vector<Real> BackScaling;             // Background scaling factor 

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

  pha& operator= (const pha&);

  // Return information

  Integer NumberChannels();            // size of internal Arrays

  // Display information about the spectrum - return as a string

  string disp();

  // Clear the information in the spectrum

  void clear();

  // Check completeness and consistency of information in spectrum
  // if there is a problem then return diagnostic in string

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

  Integer rebinChannels(grouping&);
  Integer rebinChannels(grouping&, string);

  // Shift channels

  Integer shiftChannels(Integer Start, Integer End, Real Shift);

  // Convert flux units from whatever they are currently to ph/cm^2/s. 
  // This requires as input the channel energy arrays from the rmf object and
  // the string specifying their units.

  Integer convertUnits(vector<Real>& ChannelLowEnergy, vector<Real>& ChannelHighEnergy, string EnergyUnits);

};

// Binary operation

pha operator+ (const pha& a, const pha& b);

// Utility routines

// return the type of a PHA extension

Integer PHAtype(string filename, Integer PHAnumber); 
Integer PHAtype(string filename, Integer PHAnumber, Integer& Status); 

// return true if COUNTS column exists and is integer

bool IsPHAcounts(string filename, Integer PHAnumber); 
bool IsPHAcounts(string filename, Integer PHAnumber, Integer& Status); 

// return the number of spectra in a type II PHA extension

Integer NumberofSpectra(string filename, Integer PHAnumber); 
Integer NumberofSpectra(string filename, Integer PHAnumber, Integer& Status); 




