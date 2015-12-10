////////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <algorithm>

#include "ReportTableBase.h"

////////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// ReportIterator methods
////////////////////////////////////////////////////////////////////////////////
TableIterator::TableIterator(ReportTableBase &r)
{
  // keep iterated class pointer
  m_reportTableBase = &r;
  // add iterator pointer to iterated object iterators list
  m_reportTableBase->m_allIterators.push_back(this);
  // define iterator to use, and initialize it
  if(m_reportTableBase->m_filteredTable.size()  == 0) {
    index = 0;
    m_tableIterator == m_reportTableBase->m_cells.begin();
  } else {
    index = -1;
    m_rowIterator = m_reportTableBase->m_filteredTable.begin();
  }
}
////////////////////////////////////////////////////////////////////////////////
TableIterator::~TableIterator()
{
  // remove iterator from the iterated object iterators list
  m_reportTableBase->m_allIterators.remove(this);
}
////////////////////////////////////////////////////////////////////////////////
TableIterator &TableIterator::begin(void)
{
  if(m_reportTableBase->m_filteredTable.size()  == 0)
    m_tableIterator = m_reportTableBase->m_cells.begin();
  else
    m_rowIterator = m_reportTableBase->m_filteredTable.begin();
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
TableIterator &TableIterator::end(void)
{
  if(m_reportTableBase->m_filteredTable.size()  == 0)
    m_tableIterator = m_reportTableBase->m_cells.end();
  else
    m_rowIterator = m_reportTableBase->m_filteredTable.end();
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
void TableIterator::moveTo(int count)
{
  begin();

  for(int i = 0 ; i < count ; i++)
    (*this)++;
}
////////////////////////////////////////////////////////////////////////////////
TableIterator TableIterator::operator++(int i)
{
  if(m_reportTableBase->m_filteredTable.size()  == 0)
    m_tableIterator++;
  else
    m_rowIterator++;
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
TableIterator TableIterator::operator--(int i)
{
  if(m_reportTableBase->m_filteredTable.size()  == 0)
    m_tableIterator--;
  else
    m_rowIterator--;
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
TableIterator TableIterator::operator+=(int i)
{
  if(m_reportTableBase->m_filteredTable.size()  == 0)
    m_tableIterator += i;
  else {
    for(int j = 0 ; j < i ; j++)
      m_rowIterator++;
  }
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
bool TableIterator::operator==(const TableIterator &o) const
{
  // different fltered list sizes means the iterators are different
  if(o.m_reportTableBase->m_filteredTable.size() != m_reportTableBase->m_filteredTable.size())
    return false;

  // if list has size greater than 0, then compare row iterators
  if(m_reportTableBase->m_filteredTable.size())
    return m_rowIterator == o.m_rowIterator;

  // else compare table iterators
  return m_tableIterator == o.m_tableIterator;
}
////////////////////////////////////////////////////////////////////////////////
bool TableIterator::operator!=(const TableIterator &o) const
{
  // different fltered list sizes means the iterators are different
  if(o.m_reportTableBase->m_filteredTable.size() != m_reportTableBase->m_filteredTable.size())
    return true;

  // if list has size greater than 0, then compare row iterators
  if(m_reportTableBase->m_filteredTable.size())
    return m_rowIterator != o.m_rowIterator;

  // else compare table iterators
  return m_tableIterator != o.m_tableIterator;
}
////////////////////////////////////////////////////////////////////////////////
vector<string> *TableIterator::operator*()
{
  if(m_reportTableBase->m_filteredTable.size()  == 0)
    return &(*m_tableIterator);
  return *m_rowIterator;
}
////////////////////////////////////////////////////////////////////////////////
vector<string> *TableIterator::operator->()
{
  if(m_reportTableBase->m_filteredTable.size()  == 0)
    return &(*m_tableIterator);
  return *m_rowIterator;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// table sorter methods
////////////////////////////////////////////////////////////////////////////////
bool ReportTableSortHandler::operator()(const vector<string> &a, const vector<string> &b)
{
  //if(a.size() <= m_sortColumn || b.size() <= m_sortColumn) return false;

  //string sa = a[m_sortColumn];
  //string sb = b[m_sortColumn];

  int ai = getInt(a[m_sortColumn].c_str());
  int bi = getInt(b[m_sortColumn].c_str());

  return ai < bi;
  //if(sa.length() < sb.length())
  //  return true;
  //return (sa.compare(sb) < 0);
}



////////////////////////////////////////////////////////////////////////////////
// Report table base methods
////////////////////////////////////////////////////////////////////////////////
ReportTableBase::ReportTableBase(const string &projectPath, const string &tableFilename)
: m_renderHeader(true),
  m_renderBorders(true),
  m_tableFilename(tableFilename),
  m_projectDir(projectPath),
  m_direction(0),
  m_idColumn(-1),
  m_renderingException(RENDERING_EXCEPTION_NONE)
{
  m_colTypes.resize(0);
}
////////////////////////////////////////////////////////////////////////////////
ReportTableBase::~ReportTableBase()
{
  clearView();
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableBase::clearView(void)
{
  // Delete all colType objects
  for(int i = 0; i < m_colTypes.size() ; i++)
    if(m_colTypes[i])
      delete m_colTypes[i];
  // set size to zero
  m_colTypes.clear();
}
////////////////////////////////////////////////////////////////////////////////
ReportColumnTypeBase *ReportTableBase::getColType(unsigned x) const
{
  // Check the column. If out-of-bounds, return NULL
  if(x >= m_colTypes.size())
    return NULL;
  // return a pointer to the column type object
  return m_colTypes[x];
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::setValue(unsigned row, unsigned col, string &value)
{
  // check bonds
  if(row >= m_cells.size())
    return ERROR;
  if(col >= m_cells[row].size())
    return ERROR;

  m_cells[row][col] = value;
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
string *ReportTableBase::getCell(unsigned x, unsigned y)
{
  // check row (Y) location. If out-of-bounds, return NULL
  if(y >= m_cells.size())
    return NULL;
  // get the row int an auxiliary variable
  vector<string> &aux = m_cells[y];
  // check the column within the row. If out-of-bounds, return NULL
  if(x >= aux.size())
    return NULL;
  // return a pointer to the string
  return &(aux[x]);
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableBase::getColumnId(int index)
{
  if(index < m_cells.size()) {
    if((m_cells[index].size() < m_idColumn) && (m_idColumn >= 0))
      return m_cells[index][m_idColumn];
  }
  string ret;
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
vector<string> ReportTableBase::getColumn(unsigned column)
{
  // define return structure
  vector<string> ret;

  // cycle thru all rows
  for(int i = 0 ; i < m_cells.size() ; i++) {
    // check row (column) location. If out-of-bounds, ignore
    if(column < m_cells[i].size())
      ret.push_back(m_cells[i][column]);
  }
  // return the structure
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
// Sort the table
////////////////////////////////////////////////////////////////////////////////
void ReportTableBase::sortTable(void)
{
  ReportTableSortHandler sorter(m_idColumn);
  std::sort(m_cells.begin(), m_cells.end(), sorter);
}
////////////////////////////////////////////////////////////////////////////////
// operations that return an iterator position
TableIterator ReportTableBase::begin(void)
{
  TableIterator tableIterator(*this);
  return tableIterator.begin();
}
////////////////////////////////////////////////////////////////////////////////
TableIterator ReportTableBase::end(void)
{
  TableIterator tableIterator(*this);
  return tableIterator.end();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::applyFilter(const unsigned column, const string &s)
{
  // clear the current filtered rows list
  m_filteredTable.clear();
  // check for undefined filter column, which means no filtering
  if(column == TABLE_FILTER_COL_NONE)
    return OK;
  // cycle thru all rows. If strings are equal, add to filtered rows list
  for(int i = 0 ; i < m_cells.size() ; i++) {
    // test if 'column' is out of bonds for this row
    if(column >= m_cells[i].size())
      return ERROR;
    // check if row satisfyes the condition. If so, add it to the list
    if(m_cells[i][column].compare(s) == 0)
      m_filteredTable.push_back(&(m_cells[i]));
  }

  // all existing iterators should point to "end" if filter changed
  list<TableIterator *>::iterator it;
  for(it = m_allIterators.begin() ; it != m_allIterators.end() ; it++)
    *(*it) = end();

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table element count (when filtered)
int ReportTableBase::getElementCount(void)
{
  if(m_filteredTable.size())
    return m_filteredTable.size();
  else
    return m_cells.size();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::getPageSplitStats(int itemsPerPage, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames, string &fnamePrefix, string &fnameSuffix)
{
  // calculate # of pages
  int elems = getElementCount();
  int nPages = elems / itemsPerPage;
  int elemsLeft = elems % itemsPerPage;
  if(elemsLeft > 0)
    nPages++;

  int i = 0;
  for(TableIterator it = begin() ; it != end() ; i++) {

    string fname = fnamePrefix;
    fname += parseInt(i);
    fname += fnameSuffix;
    fNames.push_back(fname);

    string start;
    if(m_idColumn < 0)
      start = parseInt(i);
    else {
      const vector<string> &aux = *(*it);
      start = parseInt(getInt(aux[m_idColumn].c_str()));
    }

    IDs.push_back(start);

    for(int j = 0 ; it != end() && j < itemsPerPage ; j++ , it++);

    it--;
    string end;
    if(m_idColumn < 0)
      end = parseInt(i);
    else {
      const vector<string> &aux = *(*it);
      end = parseInt(getInt(aux[m_idColumn].c_str()));
    }
    IDsEnd.push_back(end);
    it++;

  }

  return nPages;
}
////////////////////////////////////////////////////////////////////////////////
// Table loading and saving methods
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::loadTable(void)
{
  // make sure the table is empty;
  m_cells.empty();

  // get filename with full path
  string fn = composeFileName(m_projectDir, m_tableFilename);
  // open file to read
  ifstream file(fn.c_str(), ios::in);
    // if error, say so
  if(!file.is_open()) {
    //cerr << "ERROR: Report table: could not open file " << filename << endl;
    return ERROR;
  }

  // delimiter used to parse the file
  string delim = FILE_SEPARATOR;
  // parse the file
  string aux;
  // get one line from file
  while(getline(file, aux)) {
    // vector to store split cells
    vector<string> parts;
    // split the line
    stringSplit2(aux, parts, delim);
    // add the strings to the table
    m_cells.push_back(parts);
  }

  // close the file
  file.close();

  // return OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Saves the table in text format
int ReportTableBase::saveTable(void)
{
  // get filename with full path
  string fn = composeFileName(m_projectDir, m_tableFilename);
  // open file to write to
  ofstream file(fn.c_str(), ios::out);
    // if error, say so
  if(!file.is_open()) {
    //cerr << "ERROR: Report table: could not open file for writing: " << filename << endl;
    return ERROR;
  }

  // output file

  // output the table, seprated by FILE_SEPARATOR, one line per row
  for(int i = 0 ; i < m_cells.size() ; i++) {
    string line;
    for(int j = 0 ; j < m_cells[i].size() ; j++) {
      if(j) line += FILE_SEPARATOR;
      line += m_cells[i][j];
    }
    file << line << endl;
  }

  // close the file
  file.close();
  // return OK status
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::update(const ReportTableData &reportTableData)
{
  // load table
  if(loadTable() == ERROR)
    return ERROR;

  // check boundaries
//  if(reportTableData.updateRow < 0 || reportTableData.updateRow >= m_cells.size())
//    return ERROR;
//  if(reportTableData.updateCol < 0 || reportTableData.updateCol >= m_cells[reportTableData.updateRow].size())
//    return ERROR;

  // change data
//  m_cells[reportTableData.updateRow][reportTableData.updateCol] = reportTableData.updateText;

  // filter table, for a subset of rows
  applyFilter(reportTableData.filterCol, reportTableData.filterText);

  // iterate data
  TableIterator ti = begin();
  // cycle thru all rows
  for(; ti != end() ; ti++) {
    // get the row
    vector<string> &row = **ti;

    //cout << "Old value: " <<  row[reportTableData.updateCol] << endl;

    // process the row
    if(row.size() > reportTableData.updateCol)
      row[reportTableData.updateCol] = reportTableData.updateText;

    //cout << "New value: " <<  row[reportTableData.updateCol] << endl;
  }

  // save the table for future use
  return saveTable();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableBase::getData(ostream &ss, const ReportTableData &reportTableData, char fieldSep, char rowSep)
{
  // load table
  if(loadTable() == ERROR)
    return ERROR;

  // filter table, for a subset of rows
  applyFilter(reportTableData.filterCol, reportTableData.filterText);

  // iterate data
  bool first = true;
  TableIterator ti = begin();
  // cycle thru all rows
  for(; ti != end() ; ti++) {
    // put row separator if not the first row
    if(!first)
      ss << rowSep;
    // get the row
    vector<string> &row = **ti;

    // process the row
    for(int i = 0 ; i < row.size() ; i++) {
      // output col separator if not the first item
      if(i) ss << fieldSep;
      // output the cell contents
      ss << row[i];
    }
    // tell it's not the first row
    first = false;
  }
}
////////////////////////////////////////////////////////////////////////////////
// File name composition
string ReportTableBase::composeFileName(const string &projectDir, const string &fileName)
{
  // Compose output path
  string aux = projectDir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  // add filename
  aux += fileName;
  // return composed filename
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
