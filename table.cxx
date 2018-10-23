
#ifndef HAVE_table
#include "table.h"
#endif

#ifndef HAVE_SPio
#include "SPio.h"
#endif

#ifndef HAVE_SPutils
#include "SPutils.h"
#endif

//-------------------------------------------------------------------------------
// Class table

// default constructor

table::table()
  : Parameters(),
    Spectra(),
    ModelName(" "),
    ModelUnits(" "),
    NumIntParams(0),
    NumAddParams(0),
    isError(false),
    isRedshift(false),
    isAdditive(true),
    Energies(),
    EnergyUnits("keV"),
    LowEnergyLimit(0.0),
    HighEnergyLimit(0.0),
    Filename(" ")
{
}

// destructor

table::~table()
{
  // clear vectors with guaranteed reallocation
  vector<tableParameter>().swap(Parameters);
  vector<tableSpectrum>().swap(Spectra);
  vector<Real>().swap(Energies);
}

// Read the table from a FITS file

Integer table::read(string infilename)
{
  return this->read(infilename, true);
}

// Read the table from a FITS file. if loadAll is false then don't actualy read
// in the Spectra but set up the objects

Integer table::read(string infilename, bool loadAll)
{

  string hduName;
  string DefString;
  bool DefBool;
  Real DefReal;

  // open the FITS file at the primary extension

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(infilename, Read, false));
  } catch(...) {
    string msg = "Failed to read "+infilename;
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  Filename = infilename;

  PHDU& primary = pInfile->pHDU();

  // read the keywords we need from this extension

  DefString = " ";
  ModelName = SPreadKey(primary, "MODLNAME", DefString);
  DefString = "ph/cm^2/s";
  ModelUnits = SPreadKey(primary, "MODLUNIT", DefString);
  DefBool = true;
  isAdditive = SPreadKey(primary, "ADDMODEL", DefBool);
  DefBool = false;
  isRedshift = SPreadKey(primary, "REDSHIFT", DefBool);
  DefReal = -1.0;
  LowEnergyLimit = SPreadKey(primary, "LOELIMIT", DefReal);
  DefReal = -1.0;
  HighEnergyLimit = SPreadKey(primary, "HIELIMIT", DefReal);

  // at the moment no errors associated with model spectra

  isError = false;

  // Move to the parameters extension

  hduName = "PARAMETERS";
  ExtHDU& param = pInfile->extension(hduName);

  // set number of parameters

  NumIntParams = SPreadKey(param, "NINTPARM", (Integer)0);
  NumAddParams = SPreadKey(param, "NADDPARM", (Integer)0);

  size_t Nparams = NumIntParams + NumAddParams;

  // read the columns

  vector<string> name;
  vector<Integer> method, numbvals;
  vector<Real> initial, delta, minimum, bottom, top, maximum;
  vector<vector<Real> > value;

  SPreadCol(param, "NAME", name);
  SPreadCol(param, "METHOD", method);
  SPreadCol(param, "INITIAL", initial);
  SPreadCol(param, "DELTA", delta);
  SPreadCol(param, "MINIMUM", minimum);
  SPreadCol(param, "BOTTOM", bottom);
  SPreadCol(param, "TOP", top);
  SPreadCol(param, "MAXIMUM", maximum);
  SPreadCol(param, "NUMBVALS", numbvals);
  SPreadVectorCol(param, "VALUE", value);

  // loop round parameters setting up the objects and loading them

  for (size_t iPar=0; iPar<Nparams; iPar++) {
    tableParameter in;
    in.Name = name[iPar];
    in.InterpolationMethod = method[iPar];
    in.InitialValue = initial[iPar];
    in.Delta = delta[iPar];
    in.Minimum = minimum[iPar];
    in.Bottom = bottom[iPar];
    in.Top = top[iPar];
    in.Maximum = maximum[iPar];
    in.TabulatedValues.resize(numbvals[iPar]);
    for (size_t i=0; i<in.TabulatedValues.size(); i++) {
      in.TabulatedValues[i] = value[iPar][i];
    }

    // Force InterpolationMethod=-1 if there are no tabulated parameters
    if ( in.TabulatedValues.size() == 0 ) in.InterpolationMethod = -1;

    Parameters.push_back(in);
  }

  // Now read the energies

  hduName = "ENERGIES";
  ExtHDU& ener = pInfile->extension(hduName);

  vector<Real> eLow, eHigh;

  SPreadCol(ener, "ENERG_LO", eLow);
  SPreadCol(ener, "ENERG_HI", eHigh);

  size_t Nbins = eLow.size();

  Energies.resize(Nbins+1);
  for (size_t i=0; i<Nbins; i++) {
    Energies[i] = eLow[i];
  }
  Energies[Nbins] = eHigh[Nbins-1];
  EnergyUnits = "keV";

  // move to the SPECTRA extension

  hduName = "SPECTRA";
  ExtHDU& spec = pInfile->extension(hduName);

  // get the number of rows

  size_t Nspec = (size_t)SPreadKey(spec, "NAXIS2", (Integer)0);

  // construct the names of any columns with additive parameter spectra

  vector<string> addColName(NumAddParams);
  for (size_t iaF=0; iaF<(size_t)NumAddParams; iaF++) {
    ostringstream sname;
    sname << "ADDSP" << setfill('0') << setw(3) << iaF+1;
    addColName[iaF] = sname.str();
  }

  // If loading the spectra then Loop round taking each row in the extension at a time

  if ( loadAll ) {

    for (size_t iSpec=0; iSpec<Nspec; iSpec++) {

      tableSpectrum in;

      vector<Real> paramval;
      SPreadVectorColRow(spec, "PARAMVAL", iSpec+1, paramval);

      in.ParameterValues.resize(NumIntParams);
      for (size_t i=0; i<(size_t)NumIntParams; i++) {
	in.ParameterValues[i] = paramval[i];
      }

      vector<Real> intpspec;
      SPreadVectorColRow(spec, "INTPSPEC", iSpec+1, intpspec);
      in.Flux.resize(Nbins);
      for (size_t i=0; i<Nbins; i++) {
	in.Flux[i] = intpspec[i];
      }
      // not reading flux errors yet
      in.FluxError.resize(0);

      for (size_t iaF=0; iaF<(size_t)NumAddParams; iaF++) {

	vector<Real> addpspec;
	SPreadVectorColRow(spec, addColName[iaF], iSpec+1, addpspec);
	in.addFlux.push_back(addpspec);
	// not reading flux errors yet
	in.addFluxError.resize(0);

      }

      Spectra.push_back(in);

    }

  } else {

    // we are not loading the spectra so just set up dummy objects

    for (size_t iSpec=0; iSpec<Nspec; iSpec++) {

      tableSpectrum in;
      in.Flux.resize(0);
      in.FluxError.resize(0);
      in.ParameterValues.resize(0);
      in.addFlux.resize(0);
      in.addFluxError.resize(0);
      Spectra.push_back(in);
    
    }

  }

  return(OK);
}

// read the listed Spectra

template <class T> Integer table::readSpectra(T& spectrumList)
{

  // first check whether we have not already loaded the required spectra

  bool done = true;
  for (size_t iList=0; iList<spectrumList.size(); iList++) {
    if ( Spectra[spectrumList[iList]].Flux.size() == 0 ) {
      done = false;
      break;
    }
  }
  if ( done ) return(OK);

  // open the FITS file using the name which has been recorded in the table object 

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(Filename, Read, false));
  } catch(...) {
    string msg = "Failed to read "+Filename;
    SPreportError(NoSuchFile, msg);
    return(NoSuchFile);
  }

  // move to the SPECTRA extension

  string hduName = "SPECTRA";
  ExtHDU& spec = pInfile->extension(hduName);

  // construct the names of any columns with additive parameter spectra

  vector<string> addColName(NumAddParams);
  for (size_t iaF=0; iaF<(size_t)NumAddParams; iaF++) {
    ostringstream sname;
    sname << "ADDSP" << setfill('0') << setw(3) << iaF+1;
    addColName[iaF] = sname.str();
  }

  // Loop round reading each row required

  for (size_t iList=0; iList<spectrumList.size(); iList++) {

    size_t iSpec = spectrumList[iList];

    tableSpectrum& in = Spectra[iSpec];

    if (in.Flux.size() == 0 ) {

      vector<Real> paramval;
      SPreadVectorColRow(spec, "PARAMVAL", iSpec+1, paramval);

      in.ParameterValues.resize(NumIntParams);
      for (size_t i=0; i<(size_t)NumIntParams; i++) {
	in.ParameterValues[i] = paramval[i];
      }

      size_t Nbins = Energies.size()-1;

      vector<Real> intpspec;
      SPreadVectorColRow(spec, "INTPSPEC", iSpec+1, intpspec);
      in.Flux.resize(Nbins);
      for (size_t i=0; i<Nbins; i++) {
	in.Flux[i] = intpspec[i];
      }
      // not reading flux errors yet
      in.FluxError.resize(0);

      for (size_t iaF=0; iaF<(size_t)NumAddParams; iaF++) {

	vector<Real> addpspec;
	SPreadVectorColRow(spec, addColName[iaF], iSpec+1, addpspec);
	in.addFlux.push_back(addpspec);
	// not reading flux errors yet
	in.addFluxError.resize(0);

      }

    }

  }

  return(OK);
}

// required for linker instantiation
template Integer table::readSpectra(IntegerArray& spectrumList);
template Integer table::readSpectra(valarray<Integer>& spectrumList);
template Integer table::readSpectra(valarray<size_t>& spectrumList);
template Integer table::readSpectra(vector<size_t>& spectrumList);

// Push table Parameter object

void table::pushParameter(const tableParameter& paramObject)
{
  Parameters.push_back(paramObject);
  return;
}

// Push spectrum Parameter object

void table::pushSpectrum(const tableSpectrum& spectrumObject)
{
  Spectra.push_back(spectrumObject);
  return;
}

// Get table Parameter object (counts from zero)

tableParameter table::getParameter(Integer number)
{
  tableParameter paramObject;
  if ( number >= 0 && number < (Integer)Parameters.size() ) paramObject = Parameters[number];
  return paramObject;
}

// Get table Spectrum object (counts from zero)

tableSpectrum table::getSpectrum(Integer number)
{
  tableSpectrum spectrumObject;
  if ( number >= 0 && number < (Integer)Spectra.size() ) spectrumObject = Spectra[number];
  return spectrumObject;
}

// display information about the table - return as a string

string table::disp(const bool headerOnly)
{
  ostringstream outstr;

  outstr << "Table information : " <<endl;
  outstr << "Read from                          = " << Filename << endl;
  outstr << "Model name                         = " << ModelName << endl;
  outstr << "Model units                        = " << ModelUnits<< endl;
  outstr << "Number of interpolation parameters = " << NumIntParams << endl;
  outstr << "Number of additional parameters    = " << NumAddParams << endl;
  if ( isError ) outstr << "Model contains errors" << endl;
  if ( isRedshift ) outstr << "Model includes redshift" << endl;
  if ( isAdditive ) {
    outstr << "Model is additive" << endl;
  } else {
    outstr << "Model is multiplicative" << endl;
  }
  outstr << "Number of model energies           = " << Energies.size() << endl;
  outstr << "Energy units                       = " << EnergyUnits << endl;
  if ( LowEnergyLimit >= 0.0 )
    outstr << "Model below tabulated energies     = " << LowEnergyLimit << endl;
  if ( HighEnergyLimit >= 0.0 )
    outstr << "Model above tabulated energies     = " << HighEnergyLimit << endl;
  outstr << " " << endl;

  if ( headerOnly ) return outstr.str();

  for (size_t i=0; i<Parameters.size(); i++) {
    outstr << "Parameter " << i+1 << " : " << endl;
    outstr << Parameters[i].disp() << endl;
  }
  for (size_t i=0; i<Spectra.size(); i++) {
    outstr << "Spectrum " << i+1 << " : " << endl;
    outstr << Spectra[i].disp() << endl;
  }

  return outstr.str();
}

// clear contents of table object (mainly useful for Python)

void table::clear()
{
  Parameters.clear();
  Spectra.clear();
  ModelName = " ";
  ModelUnits = " ";
  NumIntParams = 0;
  NumAddParams = 0;
  isError = false;
  isRedshift = false;
  isAdditive = false;
  Energies.clear();
  EnergyUnits = " ";
  LowEnergyLimit = -1.0;
  HighEnergyLimit = -1.0;
  Filename = " ";
  return;
}

string table::check()
{
  ostringstream outstr;

  // check that energies are in increasing order

  bool isIncreasing(true);
  for (size_t i=0; i<Energies.size()-1; i++) {
    if ( Energies[i] > Energies[i+1] ) isIncreasing = false;
  }
  if ( !isIncreasing ) {
    outstr << "The Energies array is not in increasing order" << endl;
  }

  // check consistency of energy arrays. note that Flux may validly have size
  // zero if we have not loaded the complete Spectra objects

  for (size_t i=0; i<Spectra.size(); i++) {
    if ( Energies.size()-1 != Spectra[i].Flux.size() && Spectra[i].Flux.size() != 0 ) {
      outstr << "The size of the Energies array (" << Energies.size() 
	   << ") is not one more than that for the Flux array (" 
	   << Spectra[i].Flux.size() << ") for spectrum " << i << endl;
    }
    for (size_t j=0; j<Spectra[i].addFlux.size(); j++ ) {
      if ( Spectra[i].addFlux[j].size() != Spectra[i].Flux.size() ) {
	outstr << "The size of the addFlux array (" << Spectra[i].addFlux[j].size() 
	     << ") for additional parameter " << j+1 
	     << " is not equal to that for the Flux array (" 
	     << Spectra[i].Flux.size() << ") for spectrum in row " << i+1 << endl;
      }
    }
  }

  // check consistency of parameter arrays

  if ( NumIntParams + NumAddParams != (Integer)Parameters.size() ) {
    outstr << "The number of Parameters objects (" << Parameters.size() 
	 << ") does not match the sum of interpolated and additional parameters ("
	 << NumIntParams+NumAddParams << ")" << endl;
  }

  // check that there are the correct number of spectrum objects

  Integer Nspec(1);
  for (size_t i=0; i<(size_t)(NumIntParams+NumAddParams); i++) {
    if ( Parameters[i].InterpolationMethod >= 0 && Parameters[i].TabulatedValues.size() > 0 ) {
      Nspec *= Parameters[i].TabulatedValues.size();
    }
  }
  if ( Nspec != (Integer)Spectra.size() ) {
    outstr << "The number of Spectra objects (" << Spectra.size()
	 << ") does not match that expected from the Parameters objects ("
	 << Nspec << ")" << endl;
  }

  // check that the spectra all have the same size. Some spectra may validly have
  // zero size if they have not been read in yet.

  size_t Nbins(0);
  for (size_t iSpec=0; iSpec<(size_t)Nspec; iSpec++) {
    if ( Spectra[iSpec].Flux.size() != 0 ) {
      Nbins = Spectra[iSpec].Flux.size();
      break;
    }
  }

  for (size_t iSpec=0; iSpec<(size_t)Nspec; iSpec++) {
    if ( Spectra[iSpec].Flux.size() != Nbins && Spectra[iSpec].Flux.size() != 0 ) {
      outstr << "The spectrum for row " << iSpec+1 << " is not the same size as the first" << endl;
    }
    for (size_t iaF=0; iaF<Spectra[iSpec].addFlux.size(); iaF++) {
      if ( Spectra[iSpec].addFlux[iaF].size() != Nbins && Spectra[iSpec].addFlux[iaF].size() != 0 ) {
	outstr << "The " << iaF+1 << "th additional spectrum for row " 
	       << iSpec+1 << " is not the same size as the first spectrum" << endl;
      }
    }
  }


  return outstr.str();
}

// convert to standard units (keV and ph/cm^2/s). Note that this does not attempt to
// convert LowEnergyLimit or HighEnergyLimit because this cannot be done in all cases.

Integer table::convertUnits()
{
  Real xfactor(1.0);
  bool xwave(false);

  // set up energy/wave conversion factors and check for valid units

  Integer status(OK);

  status = calcXfactor(EnergyUnits, xwave, xfactor);
  if ( status != OK ) return(status);

  if ( xfactor != 1.0 ) {
    if ( xwave ) {
      for (size_t i=0; i<Energies.size(); i++) {
	Energies[i] = xfactor/Energies[i];
      }
    } else {
      for (size_t i=0; i<Energies.size(); i++) {
	Energies[i] *= xfactor;
      }
    }
  }

  EnergyUnits = "keV";

  // if the model is not multiplicative then we do not want to do any model
  // flux conversions because the units are just multiplicative factors

  if ( !isAdditive ) return(OK);

  bool energy(true);
  bool perwave(false);
  bool perenergy(false);
  Real yfactor(1.0);

  // Find the model units set and check for validity

  status = calcYfactor(ModelUnits, energy, perwave, perenergy, yfactor);
  if ( status != OK ) return(status);

  // if nothing more to be done leave now

  if ( yfactor == 1.0 ) return(OK);

  // set up arrays of mean energies/wavelengths and bin sizes
  // both are in keV.

  vector<Real> xmean(Energies.size()-1);
  vector<Real> xgmult(Energies.size()-1);
  vector<Real> xbinsize(Energies.size()-1);

  for (size_t i=0; i<xmean.size(); i++) {
    xmean[i] = (Energies[i]+Energies[i+1])/2.0;
    xgmult[i] = (Energies[i] * Energies[i+1]);
    xbinsize[i] = abs(Energies[i+1]-Energies[i]);
  }

  // Now do the conversions of Spectra.Flux and Spectra.addFlux. Note the six 
  // different possibilities based on whether flux is in photons or energy, 
  // is per wavelength or not, is per energy or not (both per wave and per 
  // energy can be true at the same time).

  for (size_t ispec=0; ispec<Spectra.size(); ispec++) {
    for (size_t i=0; i<Spectra[ispec].Flux.size(); i++) {
      if ( !energy && !perwave && !perenergy ) {
	// this is the one we have already done - flux must be ph/cm^2/s
      } else if ( !energy && !perwave && perenergy ) {
	Spectra[ispec].Flux[i] *= yfactor*xbinsize[i];
      } else if ( !energy && perwave && !perenergy ) {
	Spectra[ispec].Flux[i] *= yfactor*xbinsize[i]/xgmult[i];
      } else if ( energy && !perwave && !perenergy ) {
	Spectra[ispec].Flux[i] *= yfactor/xmean[i];
      } else if ( energy && !perwave && perenergy ) {
	Spectra[ispec].Flux[i] *= yfactor*xbinsize[i]/xmean[i];
      } else if ( energy && perwave && !perenergy ) {
	Spectra[ispec].Flux[i] *= yfactor*xbinsize[i]/xmean[i]/xgmult[i];
      }
    }
    for (size_t iadd=0; iadd<Spectra[ispec].addFlux.size(); iadd++) {
      for (size_t i=0; i<Spectra[ispec].addFlux[iadd].size(); i++) {
	if ( !energy && !perwave && !perenergy ) {
	  // this is the one we have already done - flux must be ph/cm^2/s
	} else if ( !energy && !perwave && perenergy ) {
	  Spectra[ispec].addFlux[iadd][i] *= yfactor*xbinsize[i];
	} else if ( !energy && perwave && !perenergy ) {
	  Spectra[ispec].addFlux[iadd][i] *= yfactor*xbinsize[i]/xgmult[i];
	} else if ( energy && !perwave && !perenergy ) {
	  Spectra[ispec].addFlux[iadd][i] *= yfactor/xmean[i];
	} else if ( energy && !perwave && perenergy ) {
	  Spectra[ispec].addFlux[iadd][i] *= yfactor*xbinsize[i]/xmean[i];
	} else if ( energy && perwave && !perenergy ) {
	  Spectra[ispec].addFlux[iadd][i] *= yfactor*xbinsize[i]/xmean[i]/xgmult[i];
	}
      }
    }
  }

  return(OK);
}

// reverse the rows. useful in the case that energies are not in increasing order

void table::reverseRows()
{

  // reverse the energies

  size_t Ne(Energies.size());
  vector<Real> TempE(Energies);
  for (size_t i=0; i<Ne; i++) Energies[i] = TempE[Ne-i-1];

  // loop round the tableSpectrum objects

  for (size_t iSpec=0; iSpec<Spectra.size(); iSpec++) {

    // reverse the Flux array

    size_t N(Spectra[iSpec].Flux.size());
    vector<Real> TempF(Spectra[iSpec].Flux);
    for (size_t i=0; i<N; i++) Spectra[iSpec].Flux[i] = TempF[N-i-1];

    // loop over any addFlux vectors

    for (size_t iaF=0; iaF<Spectra[iSpec].addFlux.size(); iaF++) {

      // reverse this addFlux array
      vector<Real> TempaF(Spectra[iSpec].addFlux[iaF]);
      for (size_t i=0; i<N; i++) Spectra[iSpec].addFlux[iaF][i] = TempaF[N-i-1];

    }

  }

  return;
}

// write to a FITS file

Integer table::write(string outfilename)
{

  vector<string> ttype;
  vector<string> tform;
  vector<string> tunit;

  // Create a new FITS file instance

  std::auto_ptr<FITS> pFits(0);

  try {                
    pFits.reset( new FITS(outfilename,Write) );
  } catch (FITS::CantCreate) {
    string msg = "Failed to create "+outfilename+" for table model file";
    SPreportError(CannotCreate, msg);
    return(CannotCreate);       
  }

  // Write the keywords to the primary header

  pFits->pHDU().addKey("HDUCLASS", "OGIP"," ");
  pFits->pHDU().addKey("HDUCLAS1", "XSPEC TABLE MODEL"," ");
  pFits->pHDU().addKey("HDUVERS", "1.1.0"," ");
  pFits->pHDU().addKey("MODLNAME", ModelName,"Table model name");
  pFits->pHDU().addKey("MODLUNIT", ModelUnits,"Table model units");
  pFits->pHDU().addKey("REDSHIFT", isRedshift,"Add redshift parameter?");
  pFits->pHDU().addKey("ADDMODEL", isAdditive,"Is model additive?");
  if ( LowEnergyLimit >= 0.0 )
    pFits->pHDU().addKey("LOELIMIT", LowEnergyLimit,
			 "Value of model below tabulated energies");
  if ( HighEnergyLimit >= 0.0 )
    pFits->pHDU().addKey("HIELIMIT", HighEnergyLimit,
			 "Value of model above tabulated energies");


  // Set up and create the PARAMETERS extension

  ttype.resize(10);
  tform.resize(10);
  tunit.resize(10);

  ttype[0] = "NAME";
  ttype[1] = "METHOD";
  ttype[2] = "INITIAL";
  ttype[3] = "DELTA";
  ttype[4] = "MINIMUM";
  ttype[5] = "BOTTOM";
  ttype[6] = "TOP";
  ttype[7] = "MAXIMUM";
  ttype[8] = "NUMBVALS";
  ttype[9] = "VALUE";

  tform[0] = "12A";
  tform[1] = "J";
  for (size_t i=2; i<8; i++) tform[i] = "E";
  tform[8] = "J";
  tform[9] = "PE";

  for (size_t i=0; i<10; i++) tunit[i] = " ";

  Table* pparam = pFits->addTable("PARAMETERS",Parameters.size(),ttype,tform,tunit);
  Table& param = *pparam;

  param.addKey("HDUCLASS", "OGIP"," ");
  param.addKey("HDUCLAS1", "XSPEC TABLE MODEL"," ");
  param.addKey("HDUCLAS2", "PARAMETERS"," ");
  param.addKey("HDUVERS", "1.0.0"," ");
  param.addKey("NINTPARM", NumIntParams,"Number of interpolation parameters ");
  param.addKey("NADDPARM", NumAddParams,"Number of additional parameters ");

  // write the parameter info. Note use of valarray because CCfits doesn't
  // support Column::write(vector<T>, Integer).

  size_t Nparams(Parameters.size());

  vector<string> names;
  for (size_t ipar=0; ipar<Nparams; ipar++) names.push_back(Parameters[ipar].Name);
  param.column("NAME").write(names,1);

  valarray<Integer> ivalues(Nparams);
  for (size_t ipar=0; ipar<Nparams; ipar++) ivalues[ipar] = Parameters[ipar].InterpolationMethod;
  param.column("METHOD").write(ivalues,1);

  valarray<Real> rvalues(Nparams);
  for (size_t ipar=0; ipar<Nparams; ipar++) rvalues[ipar] = Parameters[ipar].InitialValue;
  param.column("INITIAL").write(rvalues,1);

  for (size_t ipar=0; ipar<Nparams; ipar++) rvalues[ipar] = Parameters[ipar].Delta;
  param.column("DELTA").write(rvalues,1);

  for (size_t ipar=0; ipar<Nparams; ipar++) rvalues[ipar] = Parameters[ipar].Minimum;
  param.column("MINIMUM").write(rvalues,1);

  for (size_t ipar=0; ipar<Nparams; ipar++) rvalues[ipar] = Parameters[ipar].Bottom;
  param.column("BOTTOM").write(rvalues,1);

  for (size_t ipar=0; ipar<Nparams; ipar++) rvalues[ipar] = Parameters[ipar].Top;
  param.column("TOP").write(rvalues,1);

  for (size_t ipar=0; ipar<Nparams; ipar++) rvalues[ipar] = Parameters[ipar].Maximum;
  param.column("MAXIMUM").write(rvalues,1);

  for (size_t ipar=0; ipar<Nparams; ipar++) ivalues[ipar] = Parameters[ipar].TabulatedValues.size();
  param.column("NUMBVALS").write(ivalues,1);

  // workaround required here because Column::writeArrays will not take
  // vector<vector<T> > only vector<valarray<T> >.

  vector<valarray<Real> > pvalues(Nparams);
  for (size_t ipar=0; ipar<Nparams; ipar++) {
    pvalues[ipar].resize(Parameters[ipar].TabulatedValues.size());
    for (size_t i=0; i<Parameters[ipar].TabulatedValues.size(); i++) {
      pvalues[ipar][i] = Parameters[ipar].TabulatedValues[i];
    }
  }
  param.column("VALUE").writeArrays(pvalues,1);

  // Create the ENERGIES extension

  ttype.resize(2);
  tform.resize(2);
  tunit.resize(2);

  ttype[0] = "ENERG_LO";
  tform[0] = "E";
  tunit[0] = " ";

  ttype[1] = "ENERG_HI";
  tform[1] = "E";
  tunit[1] = " ";

  Table* penergies = pFits->addTable("ENERGIES",Energies.size()-1,ttype,tform,tunit);
  Table& energies = *penergies;

  energies.addKey("HDUCLASS", "OGIP"," ");
  energies.addKey("HDUCLAS1", "XSPEC TABLE MODEL"," ");
  energies.addKey("HDUCLAS2", "ENERGIES"," ");
  energies.addKey("HDUVERS", "1.0.0"," ");

  // write the energies

  size_t Nenergies(Energies.size()-1);
  rvalues.resize(Nenergies);

  for (size_t ien=0; ien<Nenergies; ien++) rvalues[ien] = Energies[ien];
  energies.column("ENERG_LO").write(rvalues,1);

  for (size_t ien=0; ien<Nenergies; ien++) rvalues[ien] = Energies[ien+1];
  energies.column("ENERG_HI").write(rvalues,1);

  // Create the SPECTRA extension
  // not writing flux errors yet

  ttype.resize(2+NumAddParams);
  tform.resize(2+NumAddParams);
  tunit.resize(2+NumAddParams);

  stringstream RepeatStream;
  RepeatStream << NumIntParams;
  string Repeat(RepeatStream.str());

  ttype[0] = "PARAMVAL";
  tform[0] = Repeat+"E";
  tunit[0] = " ";

  RepeatStream.str("");
  RepeatStream << Nenergies;
  Repeat = RepeatStream.str();

  ttype[1] = "INTPSPEC";
  tform[1] = Repeat+"E";
  tunit[1] = " ";

  for (size_t iadd=1; iadd<=(size_t)NumAddParams; iadd++) {
    RepeatStream.str("");
    RepeatStream << iadd;
    Repeat = RepeatStream.str();
    if ( iadd < 10 ) {
      ttype[1+iadd] = "ADDSP00"+Repeat;
    } else if ( iadd < 100 ) {
      ttype[1+iadd] = "ADDSP0"+Repeat;
    } else if ( iadd < 1000 ) {
      ttype[1+iadd] = "ADDSP"+Repeat;
    }
    RepeatStream.str("");
    RepeatStream << Nenergies;
    Repeat = RepeatStream.str();
    tform[1+iadd] = Repeat+"E";
    tunit[1+iadd] = " ";
  }

  Table* pspectra = pFits->addTable("SPECTRA",Spectra.size(),ttype,tform,tunit);
  Table& spectra = *pspectra;

  spectra.addKey("HDUCLASS", "OGIP"," ");
  spectra.addKey("HDUCLAS1", "XSPEC TABLE MODEL"," ");
  spectra.addKey("HDUCLAS2", "MODEL SPECTRA"," ");
  spectra.addKey("HDUVERS", "1.0.0"," ");

  // Write the spectra

  size_t Nspectra(Spectra.size());
  vector<valarray<Real> > rarray(Nspectra);

  if ( NumIntParams > 1 ) {
    for (size_t isp=0; isp<Nspectra; isp++) {
      rarray[isp].resize(NumIntParams);
      for (size_t j=0; j<(size_t)NumIntParams; j++) {
	rarray[isp][j] = Spectra[isp].ParameterValues[j];
      }
    }
    spectra.column("PARAMVAL").writeArrays(rarray,1);
  } else {
    rvalues.resize(Nspectra);
    for (size_t isp=0; isp<Nspectra; isp++) rvalues[isp] = Spectra[isp].ParameterValues[0];
    spectra.column("PARAMVAL").write(rvalues,1);
  }

  if ( Nenergies > 1 ) {
    for (size_t isp=0; isp<Nspectra; isp++) {
      rarray[isp].resize(Nenergies);
      for (size_t j=0; j<(size_t)Nenergies; j++) {
	rarray[isp][j] = Spectra[isp].Flux[j];
      }
    }
    spectra.column("INTPSPEC").writeArrays(rarray,1);

    for (size_t iadd=1; iadd<=(size_t)NumAddParams; iadd++) {
      for (size_t isp=0; isp<Nspectra; isp++) {
	rarray[isp].resize(Nenergies);
	for (size_t j=0; j<(size_t)Nenergies; j++) {
	  rarray[isp][j] = Spectra[isp].addFlux[iadd-1][j];
	}
      }
      spectra.column(ttype[iadd+1]).writeArrays(rarray,1);
    }
  } else {
    rvalues.resize(Nspectra);
    for (size_t isp=0; isp<Nspectra; isp++) rvalues[isp] = Spectra[isp].Flux[0];
    spectra.column("INTPSPEC").write(rvalues,1);

    for (size_t iadd=1; iadd<=(size_t)NumAddParams; iadd++) {
      for (size_t isp=0; isp<Nspectra; isp++) rvalues[isp] = Spectra[isp].addFlux[iadd-1][0];
      spectra.column(ttype[iadd+1]).write(rvalues,1);
    }

  }  


  return(OK);

}

// get values from table for input parameters using interpolation.

template <class T> Integer table::getValues(const T& parameterValues, const Real minEnergy, 
					    const Real maxEnergy, T& tableEnergyBins, 
					    T& tableValues, T& tableErrors)
{
  // if redshift is included then it will be the last entry in parameterValues so
  // ignore that for the moment but store the redshift factor

  size_t Nparin(parameterValues.size());
  Real zfact(1.0);
  if ( isRedshift) {
    zfact = 1.0 + parameterValues[Nparin-1];
    Nparin--;
  }

  if ( Nparin != (size_t)(NumIntParams+NumAddParams) ) {
    return(InconsistentNumTableParams);
  }
  
  // now split out interpolation from additional parameters
  vector<Real> interParamValues, addParamValues;
  vector<Integer> interParamIndex, addParamIndex;
  for (size_t i=0; i<Nparin; i++) {
    if ( Parameters[i].InterpolationMethod != -1 ) {
      interParamValues.push_back(parameterValues[i]);
      interParamIndex.push_back(i);
    } else {
      addParamValues.push_back(parameterValues[i]);
      addParamIndex.push_back(i);
    }
  }

  size_t Ninter = interParamValues.size();
  if ( Ninter != (size_t)NumIntParams) {
    return(InconsistentNumTableParams);
  }

  // now, test for limits.
  size_t totalBlockSize(1);
  vector<IntegerArray> prepRecordNumbers(Ninter);
  vector<bool> exactMatch(Ninter);
  for (size_t j= 0; j < Ninter; ++j) {
    Real value = interParamValues[j];
    const tableParameter& tabParam = Parameters[interParamIndex[j]];
    const vector<Real>& tabValue = tabParam.TabulatedValues;
    size_t N = tabValue.size();
    // accumulate product of numbers of parameter values,
    // will be used to calculate table offsets below.
    totalBlockSize *= N;
    size_t kExact (0);
    if ( N > 1 ) {
      // Condition for out-of-bounds test must be consistent with
      // the test for exactness below.  Note that tabValue's original
      // input comes from floats in a FITS file, not doubles.
      const Real fuzz = (value == 0.0) ? 0.0 : FUZZY;
      const Real magnitude = (value == 0.0) ? 1.0 : std::abs(value);
      if ((tabValue[0] - value)/magnitude > fuzz || 
	  (value - tabValue[N-1])/magnitude > fuzz ) {
 	return(TableParamValueOutsideRange);
      } else {
	for (; kExact < N-1; ++kExact) {
	  if ( (exactMatch[j] 
		= (std::abs((value - tabValue[kExact])/magnitude) <= fuzz)) )
	    break;
	}       
      }
    } else {
      exactMatch[j] = true;
    }

    if (exactMatch[j]) {
      // just store this here for now, will resize
      // correctly later.
      prepRecordNumbers[j].resize(1,kExact);
    } else {
      size_t n (1);
      // 2^(j+1) records needed.
      for (size_t l = 0; l <= j ; ++l ) if (!exactMatch[l]) n *= 2; 
      prepRecordNumbers[j].resize(n,0);
    }
  }
  
  IntegerArray blockOffset(Ninter,1);
  blockOffset[0] = totalBlockSize;

  size_t jj(0);
  do {
    totalBlockSize /= Parameters[interParamIndex[jj]].TabulatedValues.size(); 
    blockOffset[jj] = totalBlockSize;    
    ++jj;  
  } while ( jj < Ninter - 1);
  
  IntegerArray bracket(Ninter);
  for (size_t j = 0; j < Ninter; ++ j) {
    Real value = interParamValues[j];
    const tableParameter& tabParam = Parameters[interParamIndex[j]];
    const vector<Real>& tabValue = tabParam.TabulatedValues;
    if ( !exactMatch[j] ) {
      // bracket[j] is the arraypoint in TabulatedValues below the target
      // value, which is therefore straddled by (bracket[j],bracket[j]+1).
      // if the range is exceed, perform constant extrapolation.
      SPfind(tabValue,value,bracket[j]);
      if ( bracket[j] >= static_cast<int>(tabValue.size() - 1)) {
	exactMatch[j] = true;
	prepRecordNumbers[j][0] = tabValue.size() - 1;   
	for ( size_t l  = j ;  l < Ninter; ++l) {
	  prepRecordNumbers[l].resize(prepRecordNumbers[l].size()/2,0);   
	}
      } else if ( bracket[j] < 0) {
	exactMatch[j] = true;
	prepRecordNumbers[j][0] = 0;   
	for ( size_t l  = j ;  l < Ninter; ++l) {
	  prepRecordNumbers[l].resize(prepRecordNumbers[l].size()/2,0);   
	}      
      } else {
	if ( j == 0 ) {
	  prepRecordNumbers[0][0] = blockOffset[0]*(bracket[0]);   
	  prepRecordNumbers[0][1] = blockOffset[0]*(bracket[0] + 1);  
	} else {
	  IntegerArray& previous = prepRecordNumbers[j-1];
	  size_t MP = previous.size();
	  for (size_t k = 0; k < MP; ++k) {
	    prepRecordNumbers[j][2*k] = previous[k]  + blockOffset[j]*(bracket[j]);
	    prepRecordNumbers[j][2*k + 1] 
	      = previous[k]  + blockOffset[j]*(bracket[j] + 1);
	  }
	}
	
      }
    }  
    
    if (exactMatch[j]) {
      int kExact = prepRecordNumbers[j][0];
      if ( j == 0 ) {
	// this will evaluate correctly to the #of blocks
	// before the one containing the record of interest
	// because kExact is 0 based. 
	prepRecordNumbers[j][0] = blockOffset[0]*kExact;
      } else {
	IntegerArray& previous = prepRecordNumbers[j-1];
	prepRecordNumbers[j].resize(previous.size());
	size_t MP = previous.size();
	for (size_t k = 0; k < MP; ++k) {
	  prepRecordNumbers[j][k] = previous[k] + 
	    blockOffset[j]*kExact;
	}
      }
    }     
  }
  IntegerArray recordNumbers(prepRecordNumbers[Ninter-1]); 

  // now for each interpolation parameter i the input values lies between the
  // bracket[i] and bracket[i+1] entries in TabulatedValues. recordNumbers points
  // to the entries in Spectra which will be required.

  // calculate the fractions for each interpolation
  RealArray fraction(0.0, interParamValues.size());
  for (size_t j=0; j<fraction.size(); ++j) {
    if (!exactMatch[j]) {
      Real parValue = interParamValues[j];
      const tableParameter& tabParam = Parameters[interParamIndex[j]];
      const vector<Real>& tParValues = tabParam.TabulatedValues;
      Real x1 = tParValues[bracket[j]];
      Real x2 = tParValues[bracket[j]+1];
      if (tabParam.InterpolationMethod == 0) {
	fraction[j] = (parValue - x1)/(x2 - x1);
      } else {
	// we know x1 < parVal < x2.
	// now, we ought to check that x1,x2 > 0 earlier
	// than this point!
	fraction[j] = log(parValue/x1)/log(x2/x1);           
      }
    }
  }
  
  // find the range of tabulated energies between minEnergy and maxEnergy.
  // minEindex should be the last entry in Energies below minEnergy and
  // maxEindex should be the first entry in Energies above maxEnergy.
  // minEnergy and maxEnergy are assumed to be in observed frame while Energies
  // are in source frame so need to take into account (1+z) factor.

  size_t minEindex = 0;
  size_t maxEindex = Energies.size()-1;
  while ( Energies[minEindex] < minEnergy*zfact ) minEindex++;
  if ( minEindex > 0 ) minEindex--;
  while ( Energies[maxEindex] > maxEnergy*zfact ) maxEindex--;
  if ( maxEindex < Energies.size()-1 ) maxEindex++;

  const size_t NE(maxEindex-minEindex+1);

  // copy the table energies into the output array if the redshift is non-zero
  // then the output tableEnergyBins will be in the observed frame

  tableEnergyBins.resize(NE);
  if ( zfact == 1.0 ) {
    for (size_t ie=0; ie<NE; ie++) tableEnergyBins[ie] = Energies[ie+minEindex];
  } else {
    for (size_t ie=0; ie<NE; ie++) tableEnergyBins[ie] = Energies[ie+minEindex]/zfact;
  }

  // now set up to do the interpolation

  const size_t NP(interParamValues.size());
  const size_t NR(recordNumbers.size());
  const size_t NA(NumAddParams);

  // make sure that we have read in the Spectra objects which we will be using
  readSpectra(recordNumbers);

  // load all the spectra and variances into work arrays
  vector<RealArray> spectrumEntries(NR);
  vector<RealArray> varianceEntries(NR);

  // this should have the option of loading individual tableSpectrum elements
  // as necessary rather than holding the entire table model file in memory
  // ie the xspec readStrategy

  for (size_t i=0; i<NR; i++) {
    const tableSpectrum& tabSpec = Spectra[recordNumbers[i]];
    spectrumEntries[i].resize(NE);
    for (size_t ie=0; ie<NE; ie++) spectrumEntries[i][ie] = tabSpec.Flux[ie+minEindex];
    if ( isError ) {
      varianceEntries[i].resize(NE);
      for (size_t ie=0; ie<NE; ie++) {
	varianceEntries[i][ie] = tabSpec.FluxError[ie+minEindex]*tabSpec.FluxError[ie+minEindex];
      }
    }
  }

  // accumulate additional parameter spectra for each recordNumber
  for (size_t j=0; j<NR; j++) {
    const tableSpectrum& tabSpec = Spectra[recordNumbers[j]];
    for (size_t k=0; k<NA; k++) {
      Real parValue = addParamValues[k];
      for (size_t ie=0; ie<NE; ie++) spectrumEntries[j][ie] += tabSpec.addFlux[k][ie+minEindex] * parValue;
      if ( isError ) {
	for (size_t ie=0; ie<NE; ie++) varianceEntries[j][ie] += pow(2.0,tabSpec.addFluxError[k][ie+minEindex] * parValue);
      }
    }
  } 

  // and the interpolation

  if ( NR > 1 ) {

    size_t nd = NR;
    vector<RealArray> reducedSp(nd);
    vector<RealArray> reducedVar(nd);
    for (size_t i=0; i<nd; i++) reducedSp[i].resize(NE);
    if ( isError ) {
      for (size_t i=0; i<nd; i++) reducedVar[i].resize(NE);
    }
    for (int j=NP-1; j>=0; --j ) {
      Real factor = fraction[j];
      Real complfactor = 1.0-fraction[j];
      if ( !exactMatch[j] ) {
	nd /= 2;
	for (size_t k=0; k<nd; ++k) {
	  const size_t k2 = 2*k;
	  reducedSp[k] = complfactor*spectrumEntries[k2] + factor*spectrumEntries[k2+1];
	  if ( isError ) {
	    reducedVar[k] = complfactor*varianceEntries[k2]
		   + factor*varianceEntries[k2+1];       
	  }  
	}
      } else {
	for ( size_t k = 0; k < nd; ++k ) reducedSp[k] = spectrumEntries[k];
	if ( isError ) {
	  for (size_t k=0; k<nd; ++k) {
	    reducedVar[k] = varianceEntries[k];
	  }
	}
      }
      for (size_t i=0; i<nd; ++i) spectrumEntries[i] = reducedSp[i];
      if ( isError ) {
	for (size_t i=0; i<nd; ++i) varianceEntries[i] = reducedVar[i];
      }

    }

    if (nd != 1 ) {
      // we should never end up here
    }

  }

  // set the output arrays including time dilation factor if redshift is non-zero

  tableValues.resize(NE);
  if ( zfact == 1.0 ) {
    for (size_t ie=0; ie<NE; ie++) tableValues[ie] = spectrumEntries[0][ie];
  } else {
    for (size_t ie=0; ie<NE; ie++) tableValues[ie] = spectrumEntries[0][ie]/zfact;
  }
  if ( isError ) {
    tableErrors.resize(NE);
    if ( zfact == 1.0 ) {
      for (size_t ie=0; ie<NE; ie++) tableErrors[ie] = sqrt(varianceEntries[0][ie]);
    } else {
      for (size_t ie=0; ie<NE; ie++) tableErrors[ie] = sqrt(varianceEntries[0][ie])/zfact;
    }
  }

  return(OK);
}

// required for linker instantiation
template Integer table::getValues(const RealArray&, const Real, const Real,
				  RealArray&, RealArray&, RealArray&);
template Integer table::getValues(const vector<Real>&, const Real, const Real,
				  vector<Real>&, vector<Real>&, vector<Real>&);


//-------------------------------------------------------------------------------
// Class tableParameter

// default constructor

tableParameter::tableParameter()
  : Name(" "),
    InterpolationMethod(0),
    InitialValue(0.0),
    Delta(0.01),
    Minimum(0.0),
    Bottom(0.0),
    Top(1.0e6),
    Maximum(1.0e6),
    TabulatedValues()
{
}

// destructor

tableParameter::~tableParameter()
{
  // clear vector with guaranteed reallocation
  vector<Real>().swap(TabulatedValues);
}

// display information about the table parameter - return as a string

string tableParameter::disp()
{
  ostringstream outstr;

  outstr << "Parameter information : " << endl;
  outstr << "Parameter name       = " << Name << endl;
  outstr << "Interpolation method = ";
  if ( InterpolationMethod == 0 ) {
    outstr << "Linear interpolation" << endl;
  } else if ( InterpolationMethod == 1 ) {
    outstr << "Logarithmic interpolation" << endl;
  } else if ( InterpolationMethod == -1 ) {
    outstr << "Additional (non-interpolated)" << endl;
  } else {
    outstr << "Unrecognized interpolation method" << endl;
  }

  outstr << "Initial value        = " << InitialValue << endl;
  outstr << "Delta                = " << Delta << endl;
  outstr << "Minimum              = " << Minimum << endl;
  outstr << "Bottom               = " << Bottom << endl;
  outstr << "Top                  = " << Top << endl;
  outstr << "Maximum              = " << Maximum << endl;
  if ( InterpolationMethod != -1 ) {
    outstr << "Tabulated values     = ";
    for (size_t i=0; i<TabulatedValues.size(); i++) outstr << TabulatedValues[i] << "  ";
  }
  outstr << endl;

  return outstr.str();
}

// clear contents of the table parameter (mainly useful for Python)

void tableParameter::clear()
{
  Name = " ";
  InterpolationMethod = 0;
  InitialValue = 0.0;
  Delta = 0.0;
  Minimum = 0.0;
  Bottom = 0.0;
  Top = 0.0;
  Maximum = 0.0;
  TabulatedValues.clear();
  return;
}


//-------------------------------------------------------------------------------
// Class tableSpectrum

// default constructor

tableSpectrum::tableSpectrum()
  : Flux(),
    FluxError(),
    ParameterValues(),
    addFlux(),
    addFluxError()
{
}

// destructor

tableSpectrum::~tableSpectrum()
{
  // clear vectors with guaranteed reallocation
  vector<Real>().swap(Flux);
  vector<Real>().swap(ParameterValues);
  for (size_t i=0; i<addFlux.size(); i++) vector<Real>().swap(addFlux[i]);
  vector<vector<Real> >().swap(addFlux);

}

// push an additional parameter spectrum

void tableSpectrum::pushaddFlux(vector<Real> input)
{
  addFlux.push_back(input);
  return;
}

// get an additional parameter spectrum

vector<Real> tableSpectrum::getaddFlux(Integer Number)
{
  vector<Real> values;
  if ( Number >=0 && Number < (Integer)addFlux.size() ) {
    for (size_t i=0; i<addFlux[Number].size(); i++) values.push_back(addFlux[Number][i]);
  }
  return values;
}

// Display information about the table spectrum - return as a string

string tableSpectrum::disp()
{
  ostringstream outstr;

  outstr << "Spectrum information : " << endl;
  outstr << "Number of model flux bins        = " << Flux.size() << endl;
  outstr << "Parameter values                 = ";
  for (size_t i=0; i<ParameterValues.size(); i++) outstr << ParameterValues[i] << "  ";
  outstr << endl;
  outstr << "Number of additional flux arrays = " << addFlux.size();
  return outstr.str();
}

// clear contents of the table parameter (mainly useful for Python)

void tableSpectrum::clear()
{
  Flux.clear();
  ParameterValues.clear();
  addFlux.clear();
  return;
}
