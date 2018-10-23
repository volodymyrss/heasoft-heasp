// little program to run routines in the heasp library

#ifndef HAVE_HEASP
#include "heasp.h"
#endif
#ifndef HAVE_pha
#include "pha.h"
#endif
#ifndef HAVE_phaII
#include "phaII.h"
#endif
#ifndef HAVE_SPio
#include "SPio.h"
#endif
#ifndef HAVE_grouping
#include "grouping.h"
#endif
#ifndef HAVE_arf
#include "arf.h"
#endif
#ifndef HAVE_rmf
#include "rmf.h"
#endif
#ifndef HAVE_table
#include "table.h"
#endif

using namespace std;

// the program

int main(int argc, char* argv[])
{

  pha inpha;
  phaII inSpectra;
  grouping inGrouping;
  arf inARF;
  rmf inRMF;
  table inTable;

  Integer Status = 0;

  bool havePHA = false;
  bool havePHAII = false;
  bool haveARF = false;
  bool haveRMF = false;
  bool haveTable = false;
  bool verbosity = false;

  size_t istack=0;
  int iopt=0, nopt;
  string input, command, option;
  vector<string> words;
  vector<string> comstack;

  // store the commands

  string comstrings[] = {"quit", "exit", "arf", "debug", "response", "spectrum", "table"};
  vector<string> comvec;
  for (size_t i=0; i<7; i++) comvec.push_back(comstrings[i]);

  // check for any command line argument

  if ( argc > 1 ) {
    input = string(argv[1]);
    if ( input.find("@") != string::npos ) {
      comstack = SPreadStrings(input.substr(1,input.size()-1));
    } else if ( input.find(";") != string::npos ) {
      comstack = SPtokenize(input, ";");
    } else {
      comstack.push_back(input);
    }
  }


  // loop round commands seeking user input

  while (true) {

    // if necessary prompt the user for input

    if ( istack == comstack.size() ) {
      comstack.clear();
      istack = 0;
      cout << "SPMAN> ";
      getline(cin, input);
      if ( cin.eof() ) return(0);
      if ( input.find("@") != string::npos ) {
	comstack = SPreadStrings(input.substr(1,input.size()-1));
      } else if ( input.find(";") != string::npos ) {
	comstack = SPtokenize(input, ";");
      } else {
	comstack.push_back(input);
      }
    }

    input = comstack[istack++];

    // parse the next input line

    words = SPtokenize(input, " ");
    command = words[0];
    transform(command.begin(),command.end(),command.begin(),::tolower);
    nopt = (int) words.size();
    if ( nopt > 1 ) {
      iopt = 1;
      option = words[iopt];
    }

    int nmatch;
    string comfound = SPmatchString(command, comvec, nmatch);
    if ( nmatch > 1 ) {
      cout << "Multiple commands match: please be more specific." << std::endl;
      command = "?";
    } else if ( nmatch == 1 ) {
      command = comfound;
    }

    if ( command.substr(0,1) == "$" ) {

      system((input.substr(1,input.size()-1)).c_str());

    } else if ( command == "quit" || command == "exit" ) {

      return(0);

    } else if ( command == "?" || command == "help" ) {

      cout << "Possible commands are :" << endl;
      cout << " quit, exit: exit the program" << endl;
      cout << " arf       : operate on ARF" << endl;
      cout << " debug     : switch debug status (on/off)" << endl;
      cout << " response  : operate on RMF" << endl;
      cout << " spectrum  : operate on spectrum" << endl;
      cout << " table     : operate on table" << endl;
      cout << "To see additional help on a command type the command name followed by ?" << endl;
      
//--------------------------------arf command---------------------------------

    } else if ( command == "arf" ) {

      if ( option == "?" || option == "help" ) {
	cout << "Allowed options are :" << endl;
	cout << "    add filename    : add in the ARF read from filename" << endl;
	cout << "    disp            : write summary information about the ARF" << endl;
	cout << "    read filename   : read ARF from filename" << endl;
	cout << "    write filename  : write ARF to filename" << endl;

      } else if ( option.substr(0,1) != "r" && !haveARF ) {
	cout << "There is no ARF data available - use 'arf read filename'" << endl;

      } else if ( option.substr(0,1) == "a" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  arf inARF2;
	  Status = inARF2.read(filename);
	  if ( Status != 0 ) {
	    cout << "Error " << Status << " reading " << filename << endl;
	    Status = 0;
	  } else {
	    inARF += inARF2;
	  }
	}
	
      } else if ( option.substr(0,1) == "d" ) {
	cout << inARF.disp();
	
      } else if ( option.substr(0,1) == "r" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  Status = inARF.read(filename);
	  if ( Status != 0 ) {
	    cout << "Error " << Status << " reading " << filename << endl;
	    Status = 0;
	  } else {
	    if ( inARF.NumberEnergyBins() > 0 ) haveARF = true;
	  }
	}

      } else if ( option.substr(0,1) == "w" ) {
	if ( !haveARF ) {
	  cout << "There is no ARF data available" << endl;
	} else if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  Status = inARF.write(filename);
	  if ( Status != 0 ) {
	    cout << "Failed to write " << filename << endl;
	    Status = 0;
	  }
	}

      }

//--------------------------------debug command---------------------------------

    } else if ( command == "debug" ) {

      if ( option == "?" || option == "help" ) {
	cout << "The command debug toggles the CCfits verbose option" << endl;
      } else {
	verbosity = !verbosity;
	FITS::setVerboseMode(verbosity);
      }

//-----------------------------response command-------------------------------
// options related to the response (RMF)

    } else if ( command == "response" ) {

      if ( option == "?" || option == "help" ) {
	cout << "Allowed options are :" << endl;
	cout << "    add filename        : add RMF to another read from filename" << endl;
	cout << "    compress threshold  : remove RMF elements below threshold" << endl;
	cout << "    disp                : write summary information about RMF" << endl;
	cout << "    merge <filename>    : merge RMF and ARF (optionally from filename)"<< endl;
	cout << "    normalize           : normalize response to unity at each energy" << endl;
	cout << "    random energy #     : return random channels for input energy" << endl;
	cout << "    read filename #     : read RMF from filename" << endl;
	cout << "    rebin channel       : rebin in channel space" << endl;
	cout << "    write filename      : write RMF to filename" << endl;

      } else if ( option.substr(0,3) != "rea" && !haveRMF ) {
	cout << "There is no RMF data available - use 'response read filename'" << endl;

      } else if ( option.substr(0,1) == "a" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  rmf inRMF2;
	  Status = inRMF2.read(filename);
	  if ( Status != 0 ) {
	    cout << "Error " << Status << " reading " << filename << endl;
	    Status = 0;
	  } else {
	    inRMF += inRMF2;
	  }
	}
      
      } else if ( option.substr(0,1) == "c" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify threshold value" << endl;
	} else {
	  Real thresh;
	  if (SPstring2Real(words[++iopt], thresh)) {
	    inRMF.compress(thresh);
	  } else {
	    cout << "Failed to read value from " << words[iopt] << endl;
	  }
	}

      } else if ( option.substr(0,1) == "d" ) {
	cout << inRMF.disp();

      } else if ( option.substr(0,1) == "m" ) {
	if ( !haveARF && nopt <= 2 ) {
	  cout << "There is no ARF data available and no filename is given from which to read it" << endl;
	} else {
	  if ( !haveARF ) {
	    string filename(words[++iopt]);
	    Status = inARF.read(filename);
	    if ( Status != 0 ) {
	      cout << "Failed to read ARF from " << filename << endl;
	      Status = 0;
	    } else {
	      haveARF = true;
	    }
	  }
	  if ( haveARF ) {
	    if ( !inRMF.checkCompatibility(inARF) ) {
	      cout << "ARF and RMF are incompatible" << endl;
	    } else
	      inRMF *= inARF;
	  }
	}

      } else if ( option.substr(0,2) == "ra" ) {
	if ( nopt < 4 ) {
	  cout << "Please specify energy and number of photons" << endl;
	} else {
	  string senergy(words[++iopt]);
	  Real energy;
	  if ( !SPstring2Real(senergy,energy) ) {
	    cout << "Argument for energy must be a real number" << endl;
	  } else {
	    string snum(words[++iopt]);
	    Integer number;
	    if ( !SPstring2Integer(snum, number) ) {
	      cout << "Argument for number of photons must be an integer" << endl;
	    } else {
	      vector<Integer> channel(number);
	      channel = inRMF.RandomChannels(energy, number);
	      cout << "Channels: ";
	      for (size_t i=0; i<channel.size(); i++) cout << channel[i] << " ";
	      cout << endl;
	    }
	  }
	}

      } else if ( option.substr(0,3) == "rea" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  if ( nopt == 4 ) {
	    string snum(words[++iopt]);
	    Integer extvers;
	    if ( SPstring2Integer(snum, extvers) ) {
	      Status = inRMF.read(filename, extvers);
	      if ( Status != 0 ) {
		cout << "Error " << Status << " reading response from " << filename << " extension version " << extvers << endl;
		Status = 0;
	      }
	    } else {
	      cout << "Argument after filename must be an integer" << endl;
	    }
	  } else {
	    Status = inRMF.read(filename);
	    if ( Status != 0 ) {
	      cout << "Error " << Status << " reading response from " << filename << endl;
	      Status = 0;
	    }
	  }
	  if ( inRMF.NumberEnergyBins() > 0 ) haveRMF = true;
	}

      } else if ( option.substr(0,3) == "reb" ) {
	string suboption(words[++iopt]);
	if ( suboption == "?" || suboption == "help" ) {
	  cout << "Allowed options are :" << endl;
	  cout << "    channel filename    : rebin in channel space using grouping info from filename" << endl;
	  cout << "    energy  filename    : rebin in energy space using grouping info from filename" << endl;
	  cout << "    channel #           : rebin in channel space by a factor of #" << endl;
	  cout << "    energy  #           : rebin in energy space by a factor of #" << endl;
	} else if ( suboption.substr(0,1) == "c" || suboption.substr(0,1) == "e" ) {
	  if ( nopt <= 3 ) {
	    cout << "Please specify filename or integer factor" << endl;
	  } else {
	    string filename(words[++iopt]);
	    Integer binFactor;
	    if ( SPstring2Integer(filename, binFactor) ) {
	      inGrouping.load(binFactor, inRMF.NumberChannels());
	    } else {
	      Status = inGrouping.read(filename, inRMF.NumberChannels(), inRMF.FirstChannel);
	      if ( Status != 0 ) cout << "Error " << Status << " reading " << filename << endl;
	    }
	    if ( suboption.substr(0,1) == "c" ) {
	      Status = inRMF.rebinChannels(inGrouping);
	    } else {
	      Status = inRMF.rebinEnergies(inGrouping);
	    }
	    cout << "Done rebinning: Status = " << Status << endl;
	    if ( Status != 0 ) {
	      cout << "Error " << Status << " rebinning RMF" << endl;
	      Status = 0;
	    }
	  }	    
	} else {
	  cout << "Unrecognised option for 'response rebin'" << endl;
	}
	  
      } else if ( option.substr(0,1) == "w" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  Status = inRMF.writeChannelBounds(filename);
	  Status = inRMF.writeMatrix(filename);
	  if ( Status != 0 ) {
	    cout << "Failed to write " << filename << endl;
	    Status = 0;
	  }
	}

      } else if ( option.substr(0,1) == "n" ) {
	inRMF.normalize();

      } else {
	cout << option << " is not a recognized argument for response" << endl;
      }

//--------------------------------spectrum command---------------------------------

    } else if ( command == "spectrum" ) {

      if ( option == "?" || option == "help" ) {
	cout << "Allowed options are :" << endl;
	cout << "    disp                     : write summary information on spectrum" << endl;
	cout << "    read filename [#] [#] ...: read spectrum from file" << endl;
	cout << "    rebin filename           : rebin spectrum based on grouping info in filename" << endl;
	cout << "    rebin min #              : rebin spectrum to a minimum number of counts" << endl;
	cout << "    rebin sn #               : rebin spectrum to a minimum signal to noise" << endl;
	cout << "    rebin #                  : rebin spectrum by a factor #" << endl;
	cout << "    select # #               : select channels from # to #" << endl;
	cout << "    write filename [copyfn]  : write spectrum to file" << endl;

      } else if ( option.substr(0,3) != "rea" && !havePHA && !havePHAII ) {
	cout << "There is no spectrum data available" << endl;

      } else if ( option.substr(0,1) == "d" ) {
	if ( havePHA ) cout << inpha.disp();
	if ( havePHAII ) cout << inSpectra.disp();

      } else if ( option.substr(0,3) == "rea" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename [#] [#] ..." << endl;
	} else {
	  string filename(words[++iopt]);
	  FITS::setVerboseMode(false);
	  Integer type = PHAtype(filename, 1, Status);
	  vector<Integer> extNums = SpectrumExtensions(filename, Status);
	  cout << SPgetErrorStack() << endl;
	  FITS::clearErrors();
	  FITS::setVerboseMode(verbosity);
	  if ( Status == 0 ) {
	    cout << "File type is " << type << endl;
	    cout << "Number of spectrum extensions is " << extNums.size() << endl;
	    if ( type == 1 ) {
	      if ( extNums.size() > 1 ) {
		Integer MaxSpectrumNumber = extNums.size();
		vector<Integer> SpectrumNumber;
		if ( nopt > 3 ) {
		  SpectrumNumber.resize(nopt-3);
		  for (size_t i=0; i<(size_t)(nopt-3); i++) {
		    stringstream specnum;
		    specnum << words[++iopt];
		    specnum >> SpectrumNumber[i];
		  }
		} else {
		  SpectrumNumber.resize(MaxSpectrumNumber);
		  for (size_t i=0; i<(size_t)MaxSpectrumNumber; i++) SpectrumNumber[i] = i+1;
		}
		for (size_t i=0; i<SpectrumNumber.size(); i++) {
		  pha tmpPHA;
		  Status = tmpPHA.read(filename, SpectrumNumber[i]);
		  if ( Status == 0 ) {
		    inSpectra.push(tmpPHA);
		  } else {
		    cout << "Error reading " << SpectrumNumber[i] << " from "
			 << filename << ": Status = " << Status << endl;
		  }
		}
		if ( Status == 0 ) {
		  cout << "PHA type I file " << filename << " loaded" << endl;
		  if ( extNums.size() > 1 ) {
		    havePHAII = true;
		    cout << "into type II structure." << endl;
		  } else {
		    havePHA = true;
		    cout << "into type I structure." << endl;
		  }
		} else {
		  cout << "PHA error = " << Status << endl;
		}

	      } else {
		Status = inpha.read(filename);
	      }
	      if ( Status == 0 ) {
		cout << "PHA type I file " << filename << " loaded" << endl;
		if ( extNums.size() > 1 ) {
		  havePHAII = true;
		  cout << "into type II structure." << endl;
		} else {
		  havePHA = true;
		  cout << "into type I structure." << endl;
		}
	      } else {
		cout << "PHA error = " << Status << endl;
	      }
	    } else if ( type == 2 ) {
	      Integer MaxSpectrumNumber = NumberofSpectra(filename, 1, Status);
	      vector<Integer> SpectrumNumber;
	      if ( nopt > 3 ) {
		SpectrumNumber.resize(nopt-3);
		for (size_t i=0; i<(size_t)(nopt-3); i++) {
		  stringstream specnum;
		  specnum << words[++iopt];
		  specnum >> SpectrumNumber[i];
		}
	      } else {
		SpectrumNumber.resize(MaxSpectrumNumber);
		for (size_t i=0; i<(size_t)MaxSpectrumNumber; i++) SpectrumNumber[i] = i+1;
	      }
	      Status = inSpectra.read(filename, 1, SpectrumNumber);
	      if ( Status == 0 ) {
		havePHAII = true;
		cout << "PHA type II file " << filename << " loaded" << endl;
	      } else {
		cout << "PHA error = " << Status << endl;
	      }
	    }
	  } else {
	    cout << "Error " << Status << " reading " << filename << endl;
	    Status = 0;
	  }
	}

      } else if ( option.substr(0,3) == "reb" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename, min, sn or integer factor" << endl;
	} else {
	  string filename(words[++iopt]);
	  Integer binFactor;
	  if ( filename == "min" ) {
	    string smin(words[++iopt]);
	    Integer minCounts;
	    if ( SPstring2Integer(smin, minCounts) ) {
	      inGrouping = inpha.getMinCountsGrouping(minCounts);
	    }
	  } else if ( filename == "sn" ) {
	    string smin(words[++iopt]);
	    Real minSN;
	    if ( SPstring2Real(smin, minSN) ) {
	      inGrouping = inpha.getMinSNGrouping(minSN);
	    }
	  } else if ( SPstring2Integer(filename, binFactor) ) {
	    inGrouping.load(binFactor, inpha.NumberChannels());
	  } else {
	    Status = inGrouping.read(filename, inpha.NumberChannels(), inpha.FirstChannel);
	    if ( Status != 0 ) cout << "Error " << Status << " reading " << filename << endl;
	  }
	  Status = inpha.rebinChannels(inGrouping);
	  cout << "Done rebinning: Status = " << Status << endl;
	  if ( Status != 0 ) {
	    cout << "Error " << Status << " rebinning spectrum" << endl;
	    Status = 0;
	  }
	}	    

      } else if ( option.substr(0,1) == "s" ) {
	if ( nopt <= 3 ) {
	  cout << "Please specify start and end channels" << endl;
	} else {
	  string schan(words[++iopt]);
	  Integer start, end;
	  if ( !SPstring2Integer(schan, start) ) {
	    cout << "Failed to read start channel" << std::endl;
	  } else {
	    schan = words[++iopt];
	    if ( !SPstring2Integer(schan, end) ) {
	      cout << "Failed to read end channel" << std::endl;
	    } else {
	      vector<Integer> st(1);
	      vector<Integer> en(1);
	      st[0] = start;
	      en[0] = end;
	      inpha.selectChannels(st, en);
	    }
	  }
	}	    
	
      } else if ( option.substr(0,1) == "w" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  if ( havePHA ) {
	    if ( nopt <= 3 ) {
	      Status = inpha.write(filename);
	    } else {
	      string copyfn(words[++iopt]);
	      Status = inpha.write(filename,copyfn);
	    }
	  } else if ( havePHAII ) {
	    if ( nopt <= 3 ) {
	      Status = inSpectra.write(filename);
	    } else {
	      string copyfn(words[++iopt]);
	      Status = inSpectra.write(filename,copyfn);
	    }
	  }
	  if ( Status != 0 ) {
	    cout << "Failed to write " << filename << endl;
	    Status = 0;
	  }
	}
      }

//--------------------------------table command---------------------------------

    } else if ( command == "table" ) {

      if ( option == "?" || option == "help" ) {
	cout << "Allowed options are :" << endl;
	cout << "    disp             : write information about the table" << endl;
	cout << "    disp header      : write summary information about the table" << endl;
	cout << "    read filename    : read skeleton table from filename" << endl;
	cout << "    readall filename : read total table from filename" << endl;
	cout << "    write filename   : write table to filename" << endl;

      } else if ( option.substr(0,1) != "r" && !haveTable ) {
	cout << "There is no table data available - use 'table read filename'" << endl;

      } else if ( option.substr(0,1) == "d" ) {
	if ( nopt <= 2 ) {
	  cout << inTable.disp(false);
	} else {
	  string subOption(words[++iopt]);
	  if ( subOption.substr(0,1) == "h" ) cout << inTable.disp(true);
	}
      } else if ( option.substr(0,1) == "r" ) {
	if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  if ( option.length() > 4 ) {
	    Status = inTable.read(filename);
	  } else {
	    Status = inTable.read(filename, false);
	  }
	  if ( Status != 0 ) {
	    cout << "Error " << Status << " reading " << filename << endl;
	    Status = 0;
	  } else {
	    string msg = inTable.check();
	    if ( msg.size() == 0 ) {
	      haveTable = true;
	    } else {
	      cout << msg;
	    }
	  }
	}

      } else if ( option.substr(0,1) == "w" ) {
	if ( !haveTable ) {
	  cout << "There is no table data available" << endl;
	} else if ( nopt <= 2 ) {
	  cout << "Please specify filename" << endl;
	} else {
	  string filename(words[++iopt]);
	  Status = inTable.write(filename);
	  if ( Status != 0 ) {
	    cout << "Failed to write " << filename << endl;
	    Status = 0;
	  }
	}

      }

//--------------------------------unknown command---------------------------------

    } else {
      cout << "\"" << words[0] << "\"" << " is not recognised command" << endl;
      cout << "type ? or help to get a list of valid commands" << endl;
    }


    // clear parsed input

    words.clear();

  }

  return(0);
}

// **************************************************************************

