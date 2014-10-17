
#include "rmf.h"
#ifndef HAVE_arf
#include "arf.h"
#endif

using namespace std;

int main(int argc, char* argv[])
{
  string rmffile("testin.rmf");
  string arffile("testin.arf");
  string outfile("testout.rmf");

  rmf inputRMF, outputRMF;
  arf inputARF;

  Integer Status(0);

  // read in the RMF and the ARF

  Status = inputRMF.read(rmffile);
  Status = inputARF.read(arffile);

  // remove elements from the RMF with values < 1.0e-6

  Real threshold(1.0e-6);
  inputRMF.compress(threshold);

  // multiply the compressed RMF and the ARF to make an output RMF

  if ( inputRMF.checkCompatibility(inputARF) ) {

    outputRMF = inputRMF * inputARF;

    // and write out the result copying any extra HDUs and keywords from
    // the input RMF
    
    Status = outputRMF.write(outfile, rmffile);

  }

  exit(0);
}
