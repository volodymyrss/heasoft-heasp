// Utility routines

#ifndef HAVE_SPio
#include "SPio.h"
#endif

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
    ext.column(KeyName).read(Data,RowNumber,1);
    keyValue = Data[0];
  } catch (CCfits::Table::NoSuchColumn&) {
    // No column so fall back to looking for a keyword called KeyName
    try {
      ext.readKey(KeyName, keyValue);
    } catch (CCfits::HDU::NoSuchKeyword&) {
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
    if ( !needcol ) {
      table.addKey(ColName, Data[0][0], " ");
    } else if ( needcol && !isvector) {
      valarray<T> tempData(Data.size());
      for (size_t i=0; i<Data.size(); i++) tempData[i] = Data[i][0];
      table.column(ColName).write(tempData, 1);
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
    if ( !needcol && !forceCol ) {
      table.addKey(ColName, Data[0][0], " ");
    } else {
      if ( !isvector) {
	valarray<T> tempData(Data.size());
	for (size_t i=0; i<Data.size(); i++) tempData[i] = Data[i][0];
	table.column(ColName).write(tempData, 1);
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
    isvector = isvector || (Data[i].size() > 1);
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
template Real SPreadKey(PHDU&, string, Real);
template Integer SPreadKey(PHDU&, string, Integer);

template string SPreadKey(ExtHDU&, string, string);
template bool SPreadKey(ExtHDU&, string, bool);
template Real SPreadKey(ExtHDU&, string, Real);
template Integer SPreadKey(ExtHDU&, string, Integer);

template string SPreadKey(ExtHDU&, string, Integer, string);
template bool SPreadKey(ExtHDU&, string, Integer, bool);
template Real SPreadKey(ExtHDU&, string, Integer, Real);
template Integer SPreadKey(ExtHDU&, string, Integer, Integer);

template void SPreadCol(ExtHDU&, string, valarray<Real>&);
template void SPreadCol(ExtHDU&, string, valarray<Integer>&);
template void SPreadCol(ExtHDU&, string, vector<Real>&);
template void SPreadCol(ExtHDU&, string, vector<Integer>&);
//template void SPreadCol(ExtHDU&, string, vector<string>&);

template void SPreadCol(ExtHDU&, string, Integer, valarray<Real>&);
template void SPreadCol(ExtHDU&, string, Integer, valarray<Integer>&);
template void SPreadCol(ExtHDU&, string, Integer, vector<Real>&);
template void SPreadCol(ExtHDU&, string, Integer, vector<Integer>&);
//template void SPreadCol(ExtHDU&, string, Integer, vector<string>&);

template void SPreadVectorCol(ExtHDU&, string, vector<valarray<Real> >& Data);
template void SPreadVectorCol(ExtHDU&, string, vector<valarray<Integer> >& Data);
template void SPreadVectorCol(ExtHDU&, string, vector<vector<Real> >& Data);
template void SPreadVectorCol(ExtHDU&, string, vector<vector<Integer> >& Data);

template void SPreadVectorColRow(ExtHDU&, string, Integer, valarray<Real>& Data);
template void SPreadVectorColRow(ExtHDU&, string, Integer, valarray<Integer>& Data);
template void SPreadVectorColRow(ExtHDU&, string, Integer, vector<Real>& Data);
template void SPreadVectorColRow(ExtHDU&, string, Integer, vector<Integer>& Data);

template void SPwriteKey(Table&, string, Integer, string);
template void SPwriteKey(Table&, string, bool, string);
template void SPwriteKey(Table&, string, Real, string);
template void SPwriteKey(Table&, string, string, string);

template void SPwriteCol(Table&, string, valarray<Integer>&);
template void SPwriteCol(Table&, string, valarray<Real>&);
template void SPwriteCol(Table&, string, vector<Integer>&);
template void SPwriteCol(Table&, string, vector<Real>&);
template void SPwriteCol(Table&, string, vector<string>&);

template void SPwriteCol(Table&, string, valarray<Integer>&, bool);
template void SPwriteCol(Table&, string, valarray<Real>&, bool);
template void SPwriteCol(Table&, string, vector<Integer>&, bool);
template void SPwriteCol(Table&, string, vector<Real>&, bool);
template void SPwriteCol(Table&, string, vector<string>&, bool);

template void SPwriteVectorCol(Table&, string, vector<valarray<Integer> >&);
template void SPwriteVectorCol(Table&, string, vector<valarray<Real> >&);
template void SPwriteVectorCol(Table&, string, vector<vector<Integer> >&);
template void SPwriteVectorCol(Table&, string, vector<vector<Real> >&);

template void SPwriteVectorCol(Table&, string, vector<valarray<Integer> >&, bool);
template void SPwriteVectorCol(Table&, string, vector<valarray<Real> >&, bool);
template void SPwriteVectorCol(Table&, string, vector<vector<Integer> >&, bool);
template void SPwriteVectorCol(Table&, string, vector<vector<Real> >&, bool);

template bool SPneedCol(const vector<valarray<Real> >&, bool&);
template bool SPneedCol(const vector<valarray<Integer> >&, bool&);
template bool SPneedCol(const vector<vector<Real> >&, bool&);
template bool SPneedCol(const vector<vector<Integer> >&, bool&);

template bool SPneedCol(const vector<Integer>&);
template bool SPneedCol(const vector<Real>&);
template bool SPneedCol(const vector<string>&);
template bool SPneedCol(const vector<bool>&);

template bool SPneedCol(const valarray<Integer>&);
template bool SPneedCol(const valarray<Real>&);


