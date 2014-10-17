/* C version of little program to run routines in the heasp library via the C wrappers */

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "Cheasp.h"

#define TRUE 1
#define FALSE 0

#define LINESIZE 100

struct stack{
  char** strings;
  int    size;
};

/* function definitions */

void tokenize(char* str, char* delim, struct stack* tokens);
struct stack readstrings(char* filename);
void addstring(struct stack* stack, char* input);
void clearstack(struct stack*);
int getinputline(char *str, int limit);

/* the program */

int main(int argc, char* argv[])
{

  struct PHA inSpectrum;
  struct PHA *inSpectra;
  struct ARF inARF, inARF2;
  struct RMF inRMF, inRMF2;
  struct BinFactors inBinning;

  long NumberSpectra = 0;
  long *SpectrumNumber;
  int type;

  int Status = 0;

  int havePHA = FALSE;
  int havePHAII = FALSE;
  int haveARF = FALSE;
  int haveRMF = FALSE;
  int verbosity = FALSE;

  float thresh;

  int istack=0;
  int i, iopt, nopt;
  int binFactor;

  char singlechar;
  char *input, *command, *option, *suboption, *end, *filename;
  struct stack words;
  struct stack comstack;

  input   = (char *) malloc(LINESIZE);
  command = (char *) malloc(LINESIZE);
  option  = (char *) malloc(LINESIZE);
  filename = (char *) malloc(LINESIZE);
  suboption  = (char *) malloc(LINESIZE);

  /* check for any command line argument */

  if ( argc > 1 ) {
    strcpy(input, argv[1]);
    if ( strchr(input,'@') != NULL ) {
      comstack = readstrings(input++);
    } else if ( strchr(input,';') != NULL ) {
      tokenize(input, ";", &comstack);
    } else {
      addstring(&comstack,input);
    }
  }


  /* loop round commands seeking user input */

  while (TRUE) {

    /* if necessary prompt the user for input */

    if ( istack == comstack.size ) {
      clearstack(&comstack);
      istack = 0;
      printf("CSPMAN> ");
      getinputline(input, LINESIZE);
      if ( strchr(input,'@') != NULL ) {
	comstack = readstrings(++input);
      } else if ( strchr(input,';') != NULL ) {
	tokenize(input, ";", &comstack);
      } else {
	addstring(&comstack,input);
      }
    }

    strcpy(input, comstack.strings[istack++]);

    /* parse the next input line */

    iopt = 0;
    while (input[iopt]) {
      singlechar=input[iopt++];
      tolower(singlechar);
    }


    clearstack(&words);
    tokenize(input, " ", &words);
    if ( words.size == 0 ) {
      strcpy(command,"?");
    } else {
      strcpy(command, words.strings[0]);
    }

    nopt = words.size;
    if ( nopt > 1 ) {
      iopt = 1;
      strcpy(option, words.strings[iopt]);
    }


    if ( !strcmp(command,"quit") || !strcmp(command,"exit") ) {

      return(0);

    } else if ( command[0] == '$' ) {

      system(++input);

    } else if ( !strcmp(command,"?") || !strcmp(command,"help") ) {

      printf("Possible commands are :\n");
      printf( " quit, exit: exit the program\n");
      printf( " arf       : operate on ARF\n");
      printf( " debug     : switch debug status (on/off)\n");
      printf( " response  : operate on RMF\n");
      printf( " spectrum  : operate on spectrum\n");
      printf( "To see additional help on a command type the command name followed by ?\n");
      
/*--------------------------------arf command---------------------------------*/

    } else if ( !strcmp(command,"arf") ) {

      if ( !strcmp(option,"?") || !strcmp(option,"help") ) {
	printf( "Allowed options are :\n");
	printf( "    add filename    : add in the ARF read from filename\n");
	printf( "    disp            : write summary information about the ARF\n");
	printf( "    read filename   : read ARF from filename\n");
	printf( "    write filename  : write ARF to filename\n");

      } else if ( strcmp(option,"read") && !haveARF ) {
	printf( "There is no ARF data available - use 'arf read filename'\n");

      } else if ( !strcmp(option,"add") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify filename\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  Status = ReadARF(filename, 1, &inARF2);
	  if ( Status != 0 ) {
	    printf( "Error %d reading %s\n",Status,filename);
	    Status = 0;
	  } else {
	    Status = AddARF(&inARF, &inARF2);
	    if ( Status !=0 ) {
	      printf("Trying to add incompatible ARFs\n");
	      Status = 0;
	    }
	  }
	}
	
      } else if ( !strcmp(option,"disp") ) {
	DisplayARF(&inARF);
	
      } else if ( !strcmp(option,"read") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify filename\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  Status = ReadARF(filename, 1, &inARF);
	  if ( Status != 0 ) {
	    printf( "Error %d reading %s\n",Status,filename);
	    Status = 0;
	  } else {
	    if ( inARF.NumberEnergyBins > 0 ) haveARF = TRUE;
	  }
	}

      } else if ( !strcmp(option,"write") ) {
	if ( !haveARF ) {
	  printf( "There is no ARF data available\n");
	} else if ( nopt <= 2 ) {
	  printf( "Please specify filename\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  Status = WriteARF(filename, &inARF);
	  if ( Status != 0 ) {
	    printf( "Error %d writing %s\n",Status,filename);
	    Status = 0;
	  }
	}

      }

/*--------------------------------debug command---------------------------------*/

    } else if ( !strcmp(command,"debug") ) {

      if ( !strcmp(option,"?") || !strcmp(option,"help") ) {
	printf( "The command debug toggles the CCfits verbose option\n");
      } else {
	verbosity = !verbosity;
	if ( verbosity == TRUE ) {
	  SPsetCCfitsVerbose(0);
	} else {
	  SPsetCCfitsVerbose(1);
	}

      }

/*-----------------------------response command-------------------------------*/
/* options related to the response (RMF) */

    } else if ( !strcmp(command,"response") ) {

      if ( !strcmp(option,"?") || !strcmp(option,"help") ) {
	printf( "Allowed options are :\n");
	printf( "    add filename        : add RMF to another read from filename\n");
	printf( "    compress threshold  : remove RMF elements below threshold\n");
	printf( "    disp                : write summary information about RMF\n");
	printf( "    merge <filename>    : merge RMF and ARF (optionally from filename)\n");
	printf( "    normalize           : normalize response to unity at each energy\n");
	printf( "    read filename       : read RMF from filename\n");
	printf( "    rebin channel       : rebin in channel space\n");
	printf( "    write filename      : write RMF to filename\n");

      } else if ( strcmp(option,"read") && !haveRMF ) {
	printf( "There is no RMF data available - use 'response read filename'\n");

      } else if ( !strcmp(option,"add") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify filename\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  Status = ReadRMFMatrix(filename, 1, &inRMF2);
	  Status = ReadRMFEbounds(filename, 1, &inRMF2);
	  if ( Status != 0 ) {
	    printf( "Error %d writing %s\n",Status,filename);
	    Status = 0;
	  } else {
	    Status = AddRMF(&inRMF, &inRMF2);
	    if ( Status != 0 ) {
	      printf("Trying to add incompatible RMFs\n");
	      Status = 0;
	    }
	  }
	}
      
      } else if ( !strcmp(option,"compress") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify threshold value\n");
	} else {
	  thresh = strtof(words.strings[++iopt], &end);
	  if ( strcmp(words.strings[iopt],end) ) {
	    CompressRMF(&inRMF, thresh);
	  } else {
	    printf("Failed to read float from %s\n", words.strings[iopt]);
	  }
	}

      } else if ( !strcmp(option,"disp") ) {
	DisplayRMF(&inRMF);

      } else if ( !strcmp(option,"merge") ) {
	if ( !haveARF && nopt <= 2 ) {
	  printf( "There is no ARF data available and no filename is given from which to read it\n");
	} else {
	  if ( !haveARF ) {
	    strcpy(filename, words.strings[++iopt]);
	    Status = ReadARF(filename, 1, &inARF);
	    if ( Status != 0 ) {
	      printf( "Error %d reading %s\n",Status,filename);
	      Status = 0;
	    } else {
	      haveARF = TRUE;
	    }
	  }
	  if ( haveARF ) {
	    Status = MergeARFRMF(&inARF, &inRMF);
	    if ( Status != 0 ) {
	      printf( "ARF and RMF are incompatible\n");
	      Status = 0;
	    }
	  }
	}

      } else if ( !strcmp(option,"read") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify filename\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  Status = ReadRMFMatrix(filename, 1, &inRMF);
	  Status = ReadRMFEbounds(filename, 1, &inRMF);
	  if ( Status != 0 ) {
	    printf( "Error %d reading %s\n",Status,filename);
	    Status = 0;
	  } else {
	    if ( inRMF.NumberEnergyBins > 0 ) haveRMF = TRUE;
	  }
	}

      } else if ( !strcmp(option,"rebin") ) {
	strcpy(suboption, words.strings[++iopt]);
	if ( !strcmp(suboption,"?") || !strcmp(suboption,"help") ) {
	  printf( "Allowed options are :\n");
	  printf( "    channel filename    : rebin in channel space using grouping info from filename\n");
	  printf( "    energy  filename    : rebin in energy space using grouping info from filename\n");
	  printf( "    channel #           : rebin in channel space by a factor of #\n");
	  printf( "    energy  #           : rebin in energy space by a factor of #\n");
	} else if ( !strcmp(suboption,"channel") || !strcmp(suboption,"energy") ) {
	  if ( nopt <= 3 ) {
	    printf( "Please specify filename or integer factor\n");
	  } else {
	    strcpy(filename, words.strings[++iopt]);
	    binFactor = (int)strtol(filename, &end, 0);
	    if  ( strcmp(filename,end) ) {
	      inBinning.StartBin = (long *) malloc(sizeof(long));
	      inBinning.EndBin = (long *) malloc(sizeof(long));
	      inBinning.Binning = (long *) malloc(sizeof(long));
	      inBinning.Binning[0] = binFactor;
	      inBinning.NumberBinFactors = 1;
	      if ( !strcmp(suboption,"channel") ) {
		inBinning.StartBin[0] = inRMF.FirstChannel;
		inBinning.EndBin[0] = inRMF.NumberChannels + inRMF.FirstChannel - 1;
	      } else {
		inBinning.StartBin[0] = 0;
		inBinning.EndBin[0] = inRMF.NumberEnergyBins-1;
	      }
	    } else {
	      Status = SPReadBinningFile(filename, &inBinning);
	      if ( Status != 0 ) printf( "Error %d reading %s\n", Status, filename);
	      if ( !strcmp(suboption,"channel") ) {
		for (i=0; i<inBinning.NumberBinFactors; i++) {
		  inBinning.StartBin[i] += inRMF.FirstChannel;
		  inBinning.EndBin[i] += inRMF.FirstChannel;
		}
	      }
	    }
	    if ( !strcmp(suboption,"channel") ) {
	      Status = RebinRMFChannel(&inRMF, &inBinning);
	    } else {
	      Status = RebinRMFEnergy(&inRMF, &inBinning);
	    }
	    printf("Done rebinning: Status = %d\n", Status);
	    if ( Status != 0 ) {
	      printf( "Error %d rebinning RMF\n", Status);
	      Status = 0;
	    }
	  }
	} else {
	  printf( "Unrecognised option for 'response rebin'\n");
	}
	
      } else if ( !strcmp(option,"write") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify filename\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  Status = WriteRMFMatrix(filename, &inRMF);
	  Status = WriteRMFEbounds(filename, &inRMF);
	  if ( Status != 0 ) {
	    printf( "Error %d writing %s\n",Status,filename);
	    Status = 0;
	  }
	}

      } else if ( !strcmp(option,"normalize") ) {
	NormalizeRMF(&inRMF);

      } else {
	printf("%s is not a recognized argument for response\n", option);
      }

/*--------------------------------spectrum command---------------------------------*/

    } else if ( !strcmp(command,"spectrum") ) {

      if ( !strcmp(option,"?") || !strcmp(option,"help") ) {
	printf( "Allowed options are :\n");
	printf( "    disp                     : write summary information on spectrum\n");
	printf( "    read filename [#] [#] ...: read spectrum from file\n");
	printf( "    rebin filename           : rebin spectrum based on grouping info in filename\n");
	printf( "    rebin #                  : rebin spectrum by a factor #\n");
	printf( "    write filename           : write spectrum to file\n");

      } else if ( strcmp(option,"read") && !havePHA && !havePHAII ) {
	printf( "There is no spectrum data available\n");

      } else if ( !strcmp(option,"disp") ) {
	if ( havePHA ) DisplayPHAtypeI(&inSpectrum);
	if ( havePHAII ) DisplayPHAtypeII(NumberSpectra, &inSpectra);

      } else if ( !strcmp(option,"read") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify filename [#] [#] ...\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  type = ReturnPHAtype(filename, 1);
	  if ( type == 1 || type == 2 ) {
	    printf( "File type is %d\n", type);
	    if ( type == 1 ) {
	      Status = ReadPHAtypeI(filename, 1, &inSpectrum);
	      if ( Status == 0 ) {
		havePHA = TRUE;
		printf( "PHA type I file %s loaded\n",filename);
	      } else {
		printf( "PHA error = %d\n", Status);
	      }
	    } else if ( type == 2 ) {
	      NumberSpectra = ReturnNumberofSpectra(filename, 1);
	      if ( nopt > 3 ) {
		SpectrumNumber = (long *) malloc((nopt-3)*sizeof(long));
		for (i=0; i<(nopt-3); i++) {
		  SpectrumNumber[i] = strtol(words.strings[++iopt], &end, 0);
		}
	      } else {
		SpectrumNumber = (long *) malloc(NumberSpectra*sizeof(long));
		for (i=0; i<NumberSpectra; i++) SpectrumNumber[i] = i+1;
	      }
	      Status = ReadPHAtypeII(filename, 1, NumberSpectra, SpectrumNumber, &inSpectra);
	      if ( Status == 0 ) {
		havePHAII = TRUE;
		printf( "PHA type II file %s loaded\n", filename);
	      } else {
		printf( "PHA error = %d\n", Status);
	      }
	    } else {
	      printf( "Error finding type of %s\n", filename);
	      Status = 0;
	    }
	  } else {
	    printf( "Error opening %s with error = %d\n", filename, Status);
	    Status = 0;
	  }
	}

      } else if ( !strcmp(option,"rebin") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify filename or integer factor\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  binFactor = (int)strtol(filename, &end, 0);
	  if ( strcmp(filename, end) ) {
	    inBinning.StartBin = (long *) malloc(sizeof(long));
	    inBinning.EndBin = (long *) malloc(sizeof(long));
	    inBinning.Binning = (long *) malloc(sizeof(long));
	    inBinning.Binning[0] = binFactor;
	    inBinning.NumberBinFactors = 1;
	    inBinning.StartBin[0] = inSpectrum.FirstChannel;
	    inBinning.EndBin[0] = inSpectrum.NumberChannels + inSpectrum.FirstChannel - 1;
	  } else {
	    Status = SPReadBinningFile(filename, &inBinning);
	    if ( Status != 0 ) printf( "Error %d reading %s\n", Status, filename);
	    for (i=0; i<inBinning.NumberBinFactors; i++) {
	      inBinning.StartBin[i] += inSpectrum.FirstChannel;
	      inBinning.EndBin[i] += inSpectrum.FirstChannel;
	    }
	  }
	  Status = RebinPHA(&inSpectrum, &inBinning);
	  printf("Done rebinning: Status = %d\n", Status);
	  if ( Status != 0 ) {
	    printf("Error %d rebinning spectrum\n", Status);
	    Status = 0;
	  }
	}	    

      } else if ( !strcmp(option,"write") ) {
	if ( nopt <= 2 ) {
	  printf( "Please specify filename\n");
	} else {
	  strcpy(filename, words.strings[++iopt]);
	  Status = 0;
	  if ( havePHA ) {
	    Status = WritePHAtypeI(filename, &inSpectrum);
	  } else if ( havePHAII ) {
	    Status = WritePHAtypeII(filename, NumberSpectra, &inSpectra);
	  }
	  if ( Status != 0 ) {
	    printf( "Error %d writing %s\n",Status,filename);
	    Status = 0;
	  }
	}
      }
      

/*--------------------------------unknown command---------------------------------*/

    } else {
      printf( "\" %s \" is not a recognised command\n", command);
      printf( "type ? or help to get a list of valid commands\n");
    }


    /* clear parsed input */

    clearstack(&words);

  }

  return(0);
}


/* ************************************************************************** */
/* handy function to divide a string into substrings delimited using delim    */

void tokenize(char* str, char* delim, struct stack *tokens)
{
  char *word = NULL;

  word = strtok(str, delim);
  
  while ( word != NULL ) {

    addstring(tokens, word);
    word = strtok(NULL, delim);

  }

  return;
}

/* ************************************************************************** */
/* handy function to read an ascii file and place each row into its own       */
/* element of a stack                                                         */

struct stack readstrings(char *filename)
{
  struct stack input={0};

  FILE *fptr;
  char line[512];

  if ( (fptr = fopen(filename, "r")) == NULL ) {
    printf("Failed to open %s\n", filename);
    return input;
  }

  while ( fgets(line, 512, fptr) != NULL ) {
    input.size++;
    input.strings[input.size] = (char *) malloc(strlen(line));
    strcpy(input.strings[input.size], line);
  }

  return input;

}

/* ************************************************************************** */
/* add a string to a stack                                                    */

void addstring(struct stack* stack, char* input)
{
  struct stack tempstack;
  int i;

  tempstack.size = stack->size;
  tempstack.strings = (char **) malloc (tempstack.size*sizeof(char *));
  for (i=0; i<tempstack.size; i++) {
    tempstack.strings[i] = (char *) malloc(strlen(stack->strings[i])*sizeof(char));
    strcpy(tempstack.strings[i],stack->strings[i]);
    free(stack->strings[i]);
  }



  stack->size++;
  /*  if ( stack->strings == NULL ) { */
    stack->strings = (char **) malloc(stack->size*sizeof(char *));
    /*  } else { */
    /*   stack->strings = realloc(stack->strings,stack->size*sizeof(char *)); */
    /*  } */
  for (i=0; i<stack->size-1; i++) {
    stack->strings[i] = (char *) malloc(strlen(tempstack.strings[i])*sizeof(char));
    strcpy(stack->strings[i],tempstack.strings[i]);
  }
  i = stack->size - 1;
  stack->strings[i] = (char *) malloc(strlen(input)*sizeof(char));
  strcpy(stack->strings[i],input);

  clearstack(&tempstack);

  return;
}

/* ************************************************************************** */
/* clear stack                                                                */

void clearstack(struct stack* stack)
{
  int i;
  for (i=0; i<stack->size; i++) {
    free(stack->strings[i]);
  }
  stack->size = 0;

  return;
}

/* ************************************************************************** */
/* get an input line without the newline on the end                           */

int getinputline(char *str, int limit)
{
  int i;
  char c;

  i = 0;
  while (--limit > 0 && (c=getchar()) != EOF && c != '\n') str[i++] = c;
  str[i] = '\0';
  return(i);
}
