// rmft object code. Definitions in rmft.h

#ifndef HAVE_rmft
#include "rmft.h"
#endif

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

// Class rmft

// default constructor

rmft::rmft()
  : FirstChannel(0),
    NumberGroups(),
    FirstGroup(),
    FirstEnergyGroup(),
    NumberEnergiesGroup(),
    FirstElement(),
    OrderGroup(),
    LowEnergy(),
    HighEnergy(),
    Matrix(),
    ChannelLowEnergy(),
    ChannelHighEnergy(),
    AreaScaling(1.0),
    ResponseThreshold(0.0),
    EnergyUnits("keV"),
    RMFUnits(" "),
    ChannelType("PI"),
    RMFVersion("1.1.0"),
    EBDVersion("1.1.0"),
    Telescope(" "),
    Instrument(" "),
    Detector(" "),
    Filter(" "),
    RMFType(" "),
    RMFExtensionName("MATRIX"),
    EBDExtensionName("EBOUNDS")
{
}


// Destructor

rmft::~rmft()
{
  // clear vectors with guaranteed reallocation
  vector<Integer>().swap(NumberGroups);
  vector<Integer>().swap(FirstGroup);
  vector<Integer>().swap(FirstEnergyGroup);
  vector<Integer>().swap(NumberEnergiesGroup);
  vector<Integer>().swap(FirstElement);
  vector<Integer>().swap(OrderGroup);
  vector<Real>().swap(LowEnergy);
  vector<Real>().swap(HighEnergy);
  vector<Real>().swap(Matrix);
  vector<Real>().swap(ChannelLowEnergy);
  vector<Real>().swap(ChannelHighEnergy);
}

// load object from a ResponseMatrix

void rmft::load(rmf& rmf)
{

  // Copy all the standard information which is in common between the objects

  FirstChannel = rmf.FirstChannel;

  LowEnergy.resize(rmf.LowEnergy.size());
  HighEnergy.resize(rmf.HighEnergy.size());
  for (size_t i=0; i<LowEnergy.size(); i++) LowEnergy[i] = rmf.LowEnergy[i];
  for (size_t i=0; i<HighEnergy.size(); i++) HighEnergy[i] = rmf.HighEnergy[i];

  ChannelLowEnergy.resize(rmf.ChannelLowEnergy.size());
  ChannelHighEnergy.resize(rmf.ChannelHighEnergy.size());
  for (size_t i=0; i<ChannelLowEnergy.size(); i++) ChannelLowEnergy[i] = rmf.ChannelLowEnergy[i];
  for (size_t i=0; i<ChannelHighEnergy.size(); i++) ChannelHighEnergy[i] = rmf.ChannelHighEnergy[i];

  AreaScaling = rmf.AreaScaling;
  ResponseThreshold = rmf.ResponseThreshold;

  EnergyUnits = rmf.EnergyUnits;
  RMFUnits = rmf.RMFUnits;

  ChannelType = rmf.ChannelType;
  RMFVersion = rmf.RMFVersion;
  EBDVersion = rmf.EBDVersion;
  Telescope = rmf.Telescope;
  Instrument = rmf.Instrument;
  Detector = rmf.Detector;
  Filter = rmf.Filter;
  RMFType = rmf.RMFType;
  RMFExtensionName = rmf.RMFExtensionName;
  EBDExtensionName = rmf.EBDExtensionName;

  // Grab the memory for the arrays - create temporary arrays for those of
  // size NumberTotalGroups since we can't assume that this will be the same as for
  // the ResponseMatrix. Set size of these arrays to ten times that of the ResponseMatrix.

  size_t NumberChannels(ChannelLowEnergy.size());

  NumberGroups.resize(NumberChannels);
  FirstGroup.resize(NumberChannels);

  size_t TempSize(10*rmf.FirstChannelGroup.size());

  vector<Integer> TempFirstEnergyGroup(TempSize);
  vector<Integer> TempNumberEnergiesGroup(TempSize);
  vector<Integer> TempFirstElement(TempSize);
  vector<Integer> TempOrderGroup(TempSize);

  Matrix.resize(rmf.Matrix.size());

  // Grab memory for and initialize temporary arrays for the channel

  size_t NumberEnergyBins(LowEnergy.size());

  vector<Integer> CurrentGroup(NumberEnergyBins);
  vector<Integer> GroupOffset(NumberEnergyBins);
  vector<Real>    ChannelVector(NumberEnergyBins);

  // Initialize counters for groups and matrix elements

  Integer igroup = 0;
  Integer ielement = -1;

  // Loop round channels

  for (size_t i=0; i<NumberChannels; i++) {

    for (size_t j=0; j<NumberEnergyBins; j++ ) ChannelVector[j] = 0.0;

    // Special case for first channel to initialize the temporary arrays

    if ( i == 0 ) {

      for (size_t j=0; j<NumberEnergyBins; j++) {

	if (rmf.NumberGroups[j] > 0 && 
            rmf.FirstChannelGroup[rmf.FirstGroup[j]] == rmf.FirstChannel) {
	  CurrentGroup[j] = rmf.FirstGroup[j];
	  GroupOffset[j] = 0;
	} else {
	  CurrentGroup[j] = -1;
	  GroupOffset[j] = -1;
	}

      }

    } else {

      // Calculate the new values of the temporary arrays

      for (size_t j=0; j<NumberEnergyBins; j++) {

	// If we haven't started a group yet for this energy

	if ( CurrentGroup[j] < 0 ) {
	  if (rmf.NumberGroups[j] > 0 && 
              rmf.FirstChannelGroup[rmf.FirstGroup[j]] == ((Integer)i+rmf.FirstChannel)) {
	    CurrentGroup[j] = rmf.FirstGroup[j];
	    GroupOffset[j] = 0;
	  }
		  
	} else {

	  // We have started a group. There are two possibilities for the last channel at this
          // energy - either it was part of a group in which case GroupOffset >= 0 and CurrentGroup
          // gives the group number or it was not part of a group in which case GroupOffset < 0
	  // and CurrentGroup gives the next group number

	  if ( GroupOffset[j] >= 0 ) {

	    // Check whether this channel takes us off the end of the current response energy group -
	    //  if so increment CurrentGroup and set GroupOffset to -1, if not just increment GroupOffset

	    if ( GroupOffset[j] == rmf.NumberChannelsGroup[CurrentGroup[j]]-1 ) {
	      GroupOffset[j] = -1;
	      CurrentGroup[j]++;
	    } else {
	      GroupOffset[j]++;
	    }

	  } else {

	    // We are not in a group so check whether this channel is the start of the next response
	    // energy group - if so set GroupOffset to 0 otherwise do nothing

	    if ( (Integer)i+rmf.FirstChannel == rmf.FirstChannelGroup[CurrentGroup[j]] ) GroupOffset[j] = 0;

	  }

	}

	// end loop over energy bins

      }

    }

    // Generate the response vector for this channel

    for (size_t j=0; j<NumberEnergyBins; j++ ) {

      if ( GroupOffset[j] >= 0 ) {
	ChannelVector[j] = rmf.Matrix[rmf.FirstElement[CurrentGroup[j]]+GroupOffset[j]];
      }

    }

    // Convert to the compressed form of the response

    NumberGroups[i] = 0;
    Integer counting = 0;
    for (size_t j=0; j<(size_t)rmf.NumberEnergyBins(); j++) {

      if (ChannelVector[j] > rmf.ResponseThreshold) {

	ielement++;
	Matrix[ielement] = ChannelVector[j];
	if (counting == 0) {
	  TempFirstEnergyGroup[igroup] = j;
	  TempFirstElement[igroup] = ielement;
	  NumberGroups[i]++;
	  if (NumberGroups[i]==1) FirstGroup[i] = igroup;
	}
	counting++;

      } else {

	if (counting > 0) {
	  TempNumberEnergiesGroup[igroup] = counting;
	  igroup++;
	  counting = 0;
	}

      }
    
    }

    // if we are still in a group then close it out

    if (counting > 0) {
      TempNumberEnergiesGroup[igroup] = counting+1;
      igroup++;
      counting = 0;
    }

    // end of loop over channels

  }

  // The total number of groups created is now igroup + 1 so we can now create the permanent
  // arrays and get rid of the temporarys

  size_t NumberTotalGroups = igroup;

  FirstEnergyGroup.resize(NumberTotalGroups);
  NumberEnergiesGroup.resize(NumberTotalGroups);
  FirstElement.resize(NumberTotalGroups);
  OrderGroup.resize(NumberTotalGroups);

  for (size_t i=0; i<NumberTotalGroups; i++) {
    FirstEnergyGroup[i] = TempFirstEnergyGroup[i];
    NumberEnergiesGroup[i] = TempNumberEnergiesGroup[i];
    FirstElement[i] = TempFirstElement[i];
    OrderGroup[i] = TempOrderGroup[i];
  }

  return;
}


// update the FirstGroup and FirstElement arrays from NumberGroups and
// NumberChannelsGroup, respectively.

void rmft::update()
{
  // initialize
  FirstGroup.resize(NumberGroups.size());
  FirstElement.resize(NumberEnergiesGroup.size());

  FirstGroup[0] = 0;
  for (size_t i=1; i<NumberGroups.size(); i++) {
    FirstGroup[i] = FirstGroup[i-1] + NumberGroups[i-1];
  }

  FirstElement[0];
  for (size_t i=1; i<NumberEnergiesGroup.size(); i++) {
    FirstElement[i] = FirstElement[i-1] + NumberEnergiesGroup[i-1];
  }

  return;
}

// Deep copy

rmft& rmft::operator= (const rmft& beta)
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
  for (size_t i=0; NumberGroups.size(); i++) NumberGroups[i] = beta.NumberGroups[i];
  FirstGroup.resize(beta.FirstGroup.size());
  for (size_t i=0; FirstGroup.size(); i++) FirstGroup[i] = beta.FirstGroup[i];
  

  FirstEnergyGroup.resize(beta.FirstEnergyGroup.size());
  for (size_t i=0; FirstEnergyGroup.size(); i++) FirstEnergyGroup[i] = beta.FirstEnergyGroup[i];
  NumberEnergiesGroup.resize(beta.NumberEnergiesGroup.size());
  for (size_t i=0; NumberEnergiesGroup.size(); i++) NumberEnergiesGroup[i] = beta.NumberEnergiesGroup[i];
  FirstElement.resize(beta.FirstElement.size());
  for (size_t i=0; FirstElement.size(); i++) FirstElement[i] = beta.FirstElement[i];
  OrderGroup.resize(beta.OrderGroup.size());
  for (size_t i=0; OrderGroup.size(); i++) OrderGroup[i] = beta.OrderGroup[i];
  

  LowEnergy.resize(beta.LowEnergy.size());
  for (size_t i=0; LowEnergy.size(); i++) LowEnergy[i] = beta.LowEnergy[i];
  HighEnergy.resize(beta.HighEnergy.size());
  for (size_t i=0; HighEnergy.size(); i++) HighEnergy[i] = beta.HighEnergy[i];
  

  Matrix.resize(beta.Matrix.size());
  for (size_t i=0; Matrix.size(); i++) Matrix[i] = beta.Matrix[i];
  

  ChannelLowEnergy.resize(beta.ChannelLowEnergy.size());
  for (size_t i=0; ChannelLowEnergy.size(); i++) ChannelLowEnergy[i] = beta.ChannelLowEnergy[i];
  ChannelHighEnergy.resize(beta.ChannelHighEnergy.size());
  for (size_t i=0; ChannelHighEnergy.size(); i++) ChannelHighEnergy[i] = beta.ChannelHighEnergy[i];
  
  return *this;
}

// Return information

Integer rmft::NumberChannels()               // Number of spectrum channels 
{
  return ChannelLowEnergy.size();
}

Integer rmft::NumberEnergyBins()             // Number of response energies 
{
  return LowEnergy.size();
}

Integer rmft::NumberTotalGroups()            // Total number of response groups 
{
  return FirstEnergyGroup.size();
}

Integer rmft::NumberTotalElements()          // Total number of response elements 
{
  return Matrix.size();
}

// Return the value for a particular channel and energy

Real rmft::ElementValue(Integer ChannelNumber, Integer EnergyBin)
{

  if ( ChannelNumber < FirstChannel || ChannelNumber >= FirstChannel+NumberChannels() || 
       EnergyBin < 0 || EnergyBin >= NumberEnergyBins() ) return 0.0;

  // loop round groups for this channel bin

  for(size_t i=FirstGroup[ChannelNumber];i<(size_t)(FirstGroup[ChannelNumber]+NumberGroups[ChannelNumber]);i++) {

    if( EnergyBin >= FirstEnergyGroup[i] && 
        EnergyBin < FirstEnergyGroup[i]+NumberEnergiesGroup[i]) {

      return(Matrix[FirstElement[i]+EnergyBin-FirstEnergyGroup[i]]);

    }

  }

  return 0.0;

}

// Return real vector for a particular channel  

vector<Real> rmft::RowValues(Integer Channel)
{
  vector<Real> values(NumberEnergyBins());

  for (size_t i=0; i<(size_t)NumberEnergyBins(); i++) values[i] = 0.0;

  // Loop round response groups for this channel

  for (size_t i=0; i<(size_t)NumberGroups[Channel]; i++) {

    size_t igroup = i + FirstGroup[Channel];
    size_t ivec = FirstEnergyGroup[igroup];
    size_t ielt = FirstElement[igroup];

    // loop round elements in this group - adding them to the output array

    for (size_t j=0; j<(size_t)NumberEnergiesGroup[igroup]; j++) values[ivec+j] += Matrix[ielt+j];

  }

  return values;
}

// Display information about the RMF - return as a string

string rmft::disp()
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

  if ( FirstEnergyGroup.size() > 1 ) outstr << "   FirstEnergyGroup array of size " << FirstEnergyGroup.size() << endl;
  if ( NumberEnergiesGroup.size() > 1 ) outstr << "   NumberEnergiesGroup array of size " << NumberEnergiesGroup.size() << endl;
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
  //      outstr << "  " << FirstEnergyGroup[j+FirstGroup[i]] << "  " << NumberEnergiesGroup[j+FirstGroup[i]];
  //    }
  //    outstr << endl;
  //  }

  return outstr.str();
}

void rmft::clear()
{
  FirstChannel = 0;

  NumberGroups.clear();
  FirstGroup.clear();
  FirstEnergyGroup.clear();
  NumberEnergiesGroup.clear();
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
