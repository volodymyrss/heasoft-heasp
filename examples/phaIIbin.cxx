
#include "grouping.h"
#include "phaII.h"

using namespace std;

int main(int argc, char* argv[])
{
  string infile("testin.pha");
  string outfile("testout.pha");

  phaII inputSpectra;

  Integer Status(0);


  // read in all the spectra

  Status = inputSpectra.read(infile, 1);

  Integer Nspectra = inputSpectra.NumberSpectra();

  // loop round the spectra

  for (size_t i=0; i<(size_t)Nspectra; i++) {

    // set up the grouping object to rebin by a factor of 2

    grouping groupInfo;
    groupInfo.load(2, inputSpectra.phas[i].NumberChannels());

    // rebin this spectrum

    Status = inputSpectra.phas[i].rebinChannels(groupInfo);

  }

  // write the new spectra out copying extra keywords and extensions from
  // the input file

  inputSpectra.disp();

  Status = inputSpectra.write(outfile, infile);

  exit(0);
}
