// Utility routines

#ifndef HAVE_SPio
#include "SPio.h"
#endif
#include "SPutils.h"

// Read a keyword from the primary header.

template<class T> T SPreadKey(PHDU& ext, string KeyName, T DefValue)
{

  T keyValue = DefValue;
  bool verbosity = FITS::verboseMode();
  FITS::setVerboseMode(false);

  // check whether there is a keyword called KeyName. If KeyName does not exist this
  // method fails silently and returns the input default value.

  try {
    ext.readKey(KeyName, keyValue);
  } catch (CCfits::HDU::NoSuchKeyword&) {
  } catch (CCfits::Keyword::WrongKeywordValueType&) {
  }

  FITS::clearErrors();
  FITS::setVerboseMode(verbosity);

  return keyValue;

}

// Read a keyword from an extension.

template<class T> T SPreadKey(ExtHDU& ext, string KeyName, T DefValue)
{

  T keyValue = DefValue;
  bool verbosity = FITS::verboseMode();
  FITS::setVerboseMode(false);

  // check whether there is a keyword called KeyName. If KeyName does not exist this
  // method fails silently and returns the input default value.

  try {
    ext.readKey(KeyName, keyValue);
  } catch (CCfits::HDU::NoSuchKeyword&) {
  } catch (CCfits::Keyword::WrongKeywordValueType&) {
  }

  FITS::clearErrors();
  FITS::setVerboseMode(verbosity);

  return keyValue;

}

// Read a keyword. Note for type I files RowNumber will necessarily be 1.

template<class T> T SPreadKey(ExtHDU& ext, string KeyName, Integer RowNumber, T DefValue)
{

  T keyValue = DefValue;
  bool verbosity = FITS::verboseMode();
  FITS::setVerboseMode(false);

  // First check whether there is a column with this name with a different value for
  // for each row. If there is no column then fall back to the keyword. Note that this
  // method fails silently and returns the input default value.

  try {
    vector<T> Data;
    ext.column(KeyName).read(Data,RowNumber,RowNumber);
    keyValue = Data[0];
  } catch (CCfits::Table::NoSuchColumn&) {
    // No column so fall back to looking for a keyword called KeyName
    try {
      ext.readKey(KeyName, keyValue);
    } catch (CCfits::HDU::NoSuchKeyword&) {
    } catch (CCfits::Keyword::WrongKeywordValueType&) {
    }
  }

  FITS::clearErrors();
  FITS::setVerboseMode(verbosity);

  return keyValue;

}


// Read a column into a valarray

template<class T>void SPreadCol(ExtHDU& ext, string ColName, valarray<T>& Data)
{
  return SPreadCol(ext, ColName, 1, Data);
}

// Read a column into a vector

template<class T>void SPreadCol(ExtHDU& ext, string ColName, vector<T>& Data)
{
  return SPreadCol(ext, ColName, 1, Data);
}

// Read a column into a valarray - for type I data the RowNumber is necessarily 1.

template<class T>void SPreadCol(ExtHDU& ext, string ColName, Integer RowNumber, valarray<T>& Data)
{

  T keyValue;
  bool verbosity = FITS::verboseMode();

  // first check whether there is a keyword called ColName

  try {
    FITS::setVerboseMode(false);
    ext.readKey(ColName, keyValue);
    Data.resize(1);
    Data = keyValue;
    return;
  } catch (HDU::NoSuchKeyword&) {
    FITS::clearErrors();
    FITS::setVerboseMode(verbosity);

    // no keyword so try to read the column

    try {
      Column& Col = ext.column(ColName);
      try {
	FITS::setVerboseMode(false);
	Col.read(Data,RowNumber);
      } catch(Column::WrongColumnType) {
	FITS::clearErrors();
	FITS::setVerboseMode(verbosity);
	if ( RowNumber == 1 ) {
	  Col.read(Data,1,Col.rows());
	}
      }
    } catch(Table::NoSuchColumn&){
      Data.resize(0);
      return;
    }
  }

  return;

}

// Read a column into a vector - for type I data the RowNumber is necessarily 1.

template<class T>void SPreadCol(ExtHDU& ext, string ColName, Integer RowNumber, vector<T>& Data)
{

  T keyValue;
  bool verbosity = FITS::verboseMode();

  // first check whether there is a keyword called ColName

  try {
    FITS::setVerboseMode(false);
    ext.readKey(ColName, keyValue);
    Data.resize(1);
    Data[0] = keyValue;
    return;
  } catch (HDU::NoSuchKeyword&) {
    FITS::clearErrors();
    FITS::setVerboseMode(verbosity);

    // no keyword so try to read the column

    try {
      Column& Col = ext.column(ColName);
      try {
	FITS::setVerboseMode(false);

	// workaround because Column::read(Data,RowNumber) has no option with 
	// Data as a vector
	valarray<T> tData;
	Col.read(tData,RowNumber);
	Data.resize(tData.size());
	for (size_t i=0; i<tData.size(); i++) Data[i] = tData[i];
      } catch(Column::WrongColumnType) {
	FITS::clearErrors();
	FITS::setVerboseMode(verbosity);
	if ( RowNumber == 1 ) {
	  Col.read(Data,1,Col.rows());
	}
      }
    } catch(Table::NoSuchColumn&){
      Data.resize(0);
      return;
    }
  }

  return;
}


// Read a string column into a vector. This is special case workaround because
// of the CCfits missing routine

void SPreadCol(ExtHDU& ext, string ColName, vector<string>& Data)
{

  string keyValue;
  bool verbosity = FITS::verboseMode();

  // first check whether there is a keyword called ColName

  try {
    FITS::setVerboseMode(false);
    ext.readKey(ColName, keyValue);
    Data.resize(1);
    Data[0] = keyValue;
    return;
  } catch (HDU::NoSuchKeyword&) {
    FITS::clearErrors();
    FITS::setVerboseMode(verbosity);

    // no keyword so try to read the column

    try {
      Column& Col = ext.column(ColName);
      Col.read(Data,1,Col.rows());
    } catch(Table::NoSuchColumn&){
      Data.resize(0);
      return;
    }
  }

  return;
}


// Read a vector column int a vector of valarrays.

template<class T>void SPreadVectorCol(ExtHDU& ext, string ColName, vector<valarray<T> >& Data)
{

  bool verbosity = FITS::verboseMode();

  // first try to read the column as a scalar in case it is a vector of length 1

  try {
    Column& Col = ext.column(ColName);
    try {
      FITS::setVerboseMode(false);
      valarray<T> tempData(Data.size());
      Col.read(tempData,1,Col.rows());
      Data.resize(tempData.size());
      for (size_t i=0; i<tempData.size(); i++) Data[i].resize(1);
      for (size_t i=0; i<tempData.size(); i++) Data[i][0] = tempData[i];
    } catch(...) {
      FITS::clearErrors();
      FITS::setVerboseMode(verbosity);
      Col.readArrays(Data,1,Col.rows());
    }
  } catch(Table::NoSuchColumn&){
    Data.resize(0);
  }
  return;

}

// Read a vector column int a vector of vectors.

template<class T>void SPreadVectorCol(ExtHDU& ext, string ColName, vector<vector<T> >& Data)
{
  vector<valarray<T> > Temp;
  SPreadVectorCol(ext, ColName, Temp);

  Data.clear();

  for (size_t i=0; i<Temp.size(); i++) {
    vector<T> V(Temp[i].size());
    for (size_t j=0; j<V.size(); j++) V[j] = Temp[i][j];
    Data.push_back(V);
  }

  return;
}

// Read a single row of a vector column into a valarray.

template<class T>void SPreadVectorColRow(ExtHDU& ext, string ColName, Integer RowNumber, valarray<T>& Data)
{

  bool verbosity = FITS::verboseMode();

  // first try to read the column as a scalar in case it is a vector of length 1

  try {
    Column& Col = ext.column(ColName);
    try {
      FITS::setVerboseMode(false);
      Col.read(Data,RowNumber,RowNumber);
    } catch(...) {
      FITS::clearErrors();
      FITS::setVerboseMode(verbosity);
      Col.read(Data,RowNumber);
    }
  } catch(Table::NoSuchColumn&){
    Data.resize(0);
  }
  return;

}

// Read a single row of a vector column into a vector.

template<class T>void SPreadVectorColRow(ExtHDU& ext, string ColName, Integer RowNumber, vector<T>& Data)
{

  // need to read into a valarray first since CCfits does not support reading
  // a single vector column row into a vector.

  valarray<T> Temp;
  SPreadVectorColRow(ext, ColName, RowNumber, Temp);

  Data.resize(Temp.size());
  for (size_t i=0; i<Temp.size(); i++) {
    Data[i] = Temp[i];  
  }

  return;
}

// Write a keyword - at present just a wrap-up of addKey

template <class T> void SPwriteKey(Table& table, string KeyName, T KeyValue, string Comment)

{
  table.addKey(KeyName, KeyValue, Comment);

  return;
}

// Write a column from a valarray. If the data size is 1 or all 
// values are the same then just write a keyword

template <class T> void SPwriteCol(Table& table, string ColName, valarray<T>& Data)
{

  if ( Data.size() == 0 ) {
    return;
  } else if ( !SPneedCol(Data) ) {
    table.addKey(ColName, Data[0], " ");
  } else {
    table.column(ColName).write(Data, 1);
  }

  return;
}

// Write a column from a vector. If the data size is 1 or all 
// values are the same then just write a keyword

template <class T> void SPwriteCol(Table& table, string ColName, vector<T>& Data)
{

  if ( Data.size() == 0 ) {
    return;
  } else if ( !SPneedCol(Data) ) {
    table.addKey(ColName, Data[0], " ");
  } else {
    table.column(ColName).write(Data, 1);
  }

  return;
}

// Write a column from a valarray. If the data size is 1 or all 
// values are the same then just write a keyword if the forceCol bool is false.

template <class T> void SPwriteCol(Table& table, string ColName, valarray<T>& Data, bool forceCol)
{

  if ( Data.size() == 0 ) {
    return;
  } else if ( !SPneedCol(Data) && !forceCol ) {
    table.addKey(ColName, Data[0], " ");
  } else {
    table.column(ColName).write(Data, 1);
  }

  return;
}

// Write a column from a vector. If the data size is 1 or all 
// values are the same then just write a keyword if the forceCol bool is false.

template <class T> void SPwriteCol(Table& table, string ColName, vector<T>& Data, bool forceCol)
{
  if ( Data.size() == 0 ) {
    return;
  } else if ( !SPneedCol(Data) && !forceCol ) {
    table.addKey(ColName, Data[0], " ");
  } else {
    table.column(ColName).write(Data, 1);
  }

  return;
}

// Write a column from a vector of valarrays. If the data size is 1 or all values 
// are the same then just write a keyword. If all values are the same within all
// valarrays then write a scalar column

template <class T> void SPwriteVectorCol(Table& table, string ColName, vector<valarray<T> >& Data)
{
  bool isvector;
  bool needcol;

  if ( Data.size() == 0 ) {
    return;
  } else {
    needcol = SPneedCol(Data, isvector);
    if ( !needcol && Data[0].size() > 0 ) {
      table.addKey(ColName, Data[0][0], " ");
    } else if ( needcol && !isvector) {
      valarray<T> tempData(Data.size());
      if ( Data[0].size() > 0 ) {
	for (size_t i=0; i<Data.size(); i++) tempData[i] = Data[i][0];
	table.column(ColName).write(tempData, 1);
      }
    } else {
      table.column(ColName).writeArrays(Data,1);
    }
  }

  return;
}

// Write a column from a vector of vectors. If the data size is 1 or all values 
// are the same then just write a keyword. If all values are the same within all
// valarrays then write a scalar column

template <class T> void SPwriteVectorCol(Table& table, string ColName, vector<vector<T> >& Data)
{
  vector<valarray<T> > Temp;

  for (size_t i=0; i<Data.size(); i++) {
    valarray<T> V(Data[i].size());
    for (size_t j=0; j<V.size(); j++) V[j] = Data[i][j];
    Temp.push_back(V);
  }

  SPwriteVectorCol(table, ColName, Temp);

  return;
}

// Write a column from a vector of valarrays. If the data size is 1 or all values 
// are the same then just write a keyword if the forceCol bool is false. If all 
// values are the same within all valarrays then write a scalar column

template <class T> void SPwriteVectorCol(Table& table, string ColName, vector<valarray<T> >& Data, bool forceCol)
{
  bool isvector;
  bool needcol;

  if ( Data.size() == 0 ) {
    return;
  } else {
    needcol = SPneedCol(Data, isvector);
    if ( !needcol && !forceCol && Data[0].size() > 0 ) {
      table.addKey(ColName, Data[0][0], " ");
    } else {
      if ( !isvector) {
	valarray<T> tempData(Data.size());
	if ( Data[0].size() > 0 ) {
	  for (size_t i=0; i<Data.size(); i++) tempData[i] = Data[i][0];
	  table.column(ColName).write(tempData, 1);
	}
      } else {
	table.column(ColName).writeArrays(Data,1);
      }
    }
  }

  return;
}

// Write a column from a vector of vectors. If the data size is 1 or all values 
// are the same then just write a keyword if the forceCol bool is false. If all 
// values are the same within all valarrays then write a scalar column

template <class T> void SPwriteVectorCol(Table& table, string ColName, vector<vector<T> >& Data, bool forceCol)
{
  vector<valarray<T> > Temp;

  for (size_t i=0; i<Data.size(); i++) {
    valarray<T> V(Data[i].size());
    for (size_t j=0; j<V.size(); j++) V[j] = Data[i][j];
    Temp.push_back(V);
  }

  SPwriteVectorCol(table, ColName, Temp, forceCol);

  return;
}

// check whether a given column is required and if it needs to be a vector column
// note that T is assumed to be a valarray or vector of a valarray of some type.

template <class T> bool SPneedCol(const T& Data, bool& isvector)
{

  isvector = false;
  for (size_t i=0; i<Data.size(); i++) {
    isvector = isvector || (Data[i].size() != 1);
  }

  bool needcol = isvector;
  if (!needcol) {
    for (size_t i=1; i<Data.size(); i++) {
      if (Data[i][0] != Data[0][0]) needcol = true;
    }
  }

  return needcol;
}

// check whether a given column is required. For this overloaded version
// T is assumed to be a valarray or a vector of a scalar of some type.

template <class T> bool SPneedCol(const T& Data)
{

  bool needcol = false;
  if ( Data.size() > 1 ) {
    for (size_t i=1; i<Data.size(); i++) {
      if (Data[i] != Data[0]) needcol = true;
    }
  }

  return needcol;
}


// required to make the linker instantiate correctly
// vector templates which are commented out require
// a Column::read (std::vector< S > &vals, long rows) method
// which is currently missing from CCfits.

template string SPreadKey(PHDU&, string, string);
template bool SPreadKey(PHDU&, string, bool);
template float SPreadKey(PHDU&, string, float);
template double SPreadKey(PHDU&, string, double);
template Integer SPreadKey(PHDU&, string, Integer);

template string SPreadKey(ExtHDU&, string, string);
template bool SPreadKey(ExtHDU&, string, bool);
template float SPreadKey(ExtHDU&, string, float);
template double SPreadKey(ExtHDU&, string, double);
template Integer SPreadKey(ExtHDU&, string, Integer);

template string SPreadKey(ExtHDU&, string, Integer, string);
template bool SPreadKey(ExtHDU&, string, Integer, bool);
template float SPreadKey(ExtHDU&, string, Integer, float);
template double SPreadKey(ExtHDU&, string, Integer, double);
template Integer SPreadKey(ExtHDU&, string, Integer, Integer);

template void SPreadCol(ExtHDU&, string, valarray<float>&);
template void SPreadCol(ExtHDU&, string, valarray<double>&);
template void SPreadCol(ExtHDU&, string, valarray<Integer>&);
template void SPreadCol(ExtHDU&, string, vector<float>&);
template void SPreadCol(ExtHDU&, string, vector<double>&);
template void SPreadCol(ExtHDU&, string, vector<Integer>&);
template void SPreadCol(ExtHDU&, string, vector<bool>&);
//template void SPreadCol(ExtHDU&, string, vector<string>&);

template void SPreadCol(ExtHDU&, string, Integer, valarray<float>&);
template void SPreadCol(ExtHDU&, string, Integer, valarray<double>&);
template void SPreadCol(ExtHDU&, string, Integer, valarray<Integer>&);
template void SPreadCol(ExtHDU&, string, Integer, vector<float>&);
template void SPreadCol(ExtHDU&, string, Integer, vector<double>&);
template void SPreadCol(ExtHDU&, string, Integer, vector<Integer>&);
template void SPreadCol(ExtHDU&, string, Integer, vector<bool>&);
//template void SPreadCol(ExtHDU&, string, Integer, vector<string>&);

template void SPreadVectorCol(ExtHDU&, string, vector<valarray<float> >& Data);
template void SPreadVectorCol(ExtHDU&, string, vector<valarray<double> >& Data);
template void SPreadVectorCol(ExtHDU&, string, vector<valarray<Integer> >& Data);
template void SPreadVectorCol(ExtHDU&, string, vector<vector<float> >& Data);
template void SPreadVectorCol(ExtHDU&, string, vector<vector<double> >& Data);
template void SPreadVectorCol(ExtHDU&, string, vector<vector<Integer> >& Data);

template void SPreadVectorColRow(ExtHDU&, string, Integer, valarray<float>& Data);
template void SPreadVectorColRow(ExtHDU&, string, Integer, valarray<double>& Data);
template void SPreadVectorColRow(ExtHDU&, string, Integer, valarray<Integer>& Data);
template void SPreadVectorColRow(ExtHDU&, string, Integer, vector<float>& Data);
template void SPreadVectorColRow(ExtHDU&, string, Integer, vector<double>& Data);
template void SPreadVectorColRow(ExtHDU&, string, Integer, vector<Integer>& Data);

template void SPwriteKey(Table&, string, Integer, string);
template void SPwriteKey(Table&, string, bool, string);
template void SPwriteKey(Table&, string, float, string);
template void SPwriteKey(Table&, string, double, string);
template void SPwriteKey(Table&, string, string, string);

template void SPwriteCol(Table&, string, valarray<Integer>&);
template void SPwriteCol(Table&, string, valarray<float>&);
template void SPwriteCol(Table&, string, valarray<double>&);
template void SPwriteCol(Table&, string, vector<Integer>&);
template void SPwriteCol(Table&, string, vector<float>&);
template void SPwriteCol(Table&, string, vector<double>&);
template void SPwriteCol(Table&, string, vector<string>&);

template void SPwriteCol(Table&, string, valarray<Integer>&, bool);
template void SPwriteCol(Table&, string, valarray<float>&, bool);
template void SPwriteCol(Table&, string, valarray<double>&, bool);
template void SPwriteCol(Table&, string, vector<Integer>&, bool);
template void SPwriteCol(Table&, string, vector<float>&, bool);
template void SPwriteCol(Table&, string, vector<double>&, bool);
template void SPwriteCol(Table&, string, vector<string>&, bool);

template void SPwriteVectorCol(Table&, string, vector<valarray<Integer> >&);
template void SPwriteVectorCol(Table&, string, vector<valarray<float> >&);
template void SPwriteVectorCol(Table&, string, vector<valarray<double> >&);
template void SPwriteVectorCol(Table&, string, vector<vector<Integer> >&);
template void SPwriteVectorCol(Table&, string, vector<vector<float> >&);
template void SPwriteVectorCol(Table&, string, vector<vector<double> >&);

template void SPwriteVectorCol(Table&, string, vector<valarray<Integer> >&, bool);
template void SPwriteVectorCol(Table&, string, vector<valarray<float> >&, bool);
template void SPwriteVectorCol(Table&, string, vector<valarray<double> >&, bool);
template void SPwriteVectorCol(Table&, string, vector<vector<Integer> >&, bool);
template void SPwriteVectorCol(Table&, string, vector<vector<float> >&, bool);
template void SPwriteVectorCol(Table&, string, vector<vector<double> >&, bool);

template bool SPneedCol(const vector<valarray<float> >&, bool&);
template bool SPneedCol(const vector<valarray<double> >&, bool&);
template bool SPneedCol(const vector<valarray<Integer> >&, bool&);
template bool SPneedCol(const vector<vector<float> >&, bool&);
template bool SPneedCol(const vector<vector<double> >&, bool&);
template bool SPneedCol(const vector<vector<Integer> >&, bool&);

template bool SPneedCol(const vector<Integer>&);
template bool SPneedCol(const vector<float>&);
template bool SPneedCol(const vector<double>&);
template bool SPneedCol(const vector<string>&);
template bool SPneedCol(const vector<bool>&);

template bool SPneedCol(const valarray<Integer>&);
template bool SPneedCol(const valarray<float>&);
template bool SPneedCol(const valarray<double>&);


// Read all keywords from the primary header into a vector of strings

vector<string> SPreadAllPrimaryKeywords(const string& filename)
{
  return SPreadAllKeywords(filename, string(""), (int)0);
}

// Read all keywords from an extension into a vector of strings. Use hduName="" as
// a special case of the primary header

vector<string> SPreadAllKeywords(const string& filename, const string& hduName, 
				 const int& hduNumber)
{
  vector<string> allKeywords;

  auto_ptr<FITS> pInfile(0);

  try {
    pInfile.reset(new FITS(filename, Read, false));
  } catch(...) {
    return allKeywords;
  }

  // read the keys and load into an associate array. Also get the comments and
  // history

  map<string,Keyword*> keywords;
  string comments;
  string history;
  
  if ( hduName.size() == 0 ) {
    PHDU& primary = pInfile->pHDU();
    primary.readAllKeys();
    keywords = primary.keyWord();
    comments = primary.getComments();
    history = primary.getHistory();
  } else {
    try {
      ExtHDU& ext = pInfile->extension(hduName, hduNumber);
      ext.readAllKeys();
      keywords = ext.keyWord();
      comments = ext.getComments();
      history = ext.getHistory();
    } catch(...) {
      return allKeywords;
    }
  }

  // loop round constructing the string to push into allKeywords

  allKeywords.resize(keywords.size());
  size_t i=0;
  for (map<string,Keyword*>::iterator it=keywords.begin();  it!=keywords.end(); it++) {

    stringstream keystream;
    keystream << setw(8) << left << it->first << setw(3) << " = " << setw(19);

    bool isString = ( it->second->keytype() == CCfits::Tstring );
    bool isLogical = ( it->second->keytype() == CCfits::Tlogical );
    if ( isLogical ) {
      bool tval;
      it->second->value(tval);
      if ( tval ) {
	keystream << right << "T";
      } else {
	keystream << right << "F";
      }
    } else {
      string tval;
      it->second->value(tval);
      if ( isString ) {
	tval = "'"+tval+"'";
	keystream << tval;
      } else {
	keystream << right << tval;
      }
    }

    keystream << setw(3) << " / ";
    size_t comWidth = 80 - keystream.str().size();
    keystream << setw(comWidth) << left << it->second->comment();

    allKeywords[i] = keystream.str();
    i++;
  }

  // get the COMMENT and HISTORY strings, split up and add to the allKeywords vector

  vector<string> allComments = SPtokenize(comments, "\n");
  for (size_t i=0; i<allComments.size(); i++)
    allKeywords.push_back(("COMMENT   "+allComments[i]));

  vector<string> allHistory = SPtokenize(history, "\n");
  for (size_t i=0; i<allHistory.size(); i++)
    allKeywords.push_back(("HISTORY   "+allHistory[i]));

  return allKeywords;
}
