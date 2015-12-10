////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_BASE_H__
#define __REPORT_TABLE_BASE_H__
////////////////////////////////////////////////////////////////////////////////
#include <iterator>
#include <vector>

#include "ReportColumnTypes.h"
#include "spectrum.h"

////////////////////////////////////////////////////////////////////////////////
#define OK      1
#define ERROR  -1

#define FILE_SEPARATOR      ';'
#define FILE_SEPARATOR_TSV  '\t'


#define DEFAULT_ROOT_DIRECTORY  "."

#define TABLE_FILTER_COL_NONE     -1


#define RENDERING_EXCEPTION_NONE                        0
#define RENDERING_EXCEPTION_PROTEIN_HEADER              1
#define RENDERING_EXCEPTION_MAIN_PAGE                   2
#define RENDERING_EXCEPTION_PROTEIN_COVERAGE_PAGE       3
#define RENDERING_EXCEPTION_PROTEIN_COVERAGE_CSV_PAGE   4

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
class ReportTableData {

 public:

  // Update data -- used for updating a table
  int     filterCol;
  string  filterText;

  int     updateCol;
  string  updateText;

  int     sortColumn;
  int     direction;

  int     startRow;
  int     rows;

  // constructor
  ReportTableData() :
      filterCol(-1),
      updateCol(-1),
      sortColumn(-1),
      direction(0),
      startRow(0),
      rows(-1)
  {};

};

////////////////////////////////////////////////////////////////////////////////
// We need to declare class ReportTableBase before declaring TableIterator
class ReportTableBase;
////////////////////////////////////////////////////////////////////////////////
typedef vector<ReportColumnTypeBase *>::iterator   TableCtIterator;
////////////////////////////////////////////////////////////////////////////////
/*
 *  The iterator class knows the internals of the ReportTableBase, so that it
 *  may move from one element to the next.
 */
class TableIterator {

  // if we are iterating the list, we have the pointer to the current element of the list
  list<vector<string> *>::iterator m_rowIterator;
  // If not, we use the table iterator
  vector<vector<string> >::iterator m_tableIterator;
  // if the list is empty, we are iterating the table. we keep the index of the current position
  int index;

  // pointer to itererated class
  ReportTableBase *m_reportTableBase;

 public:

  // constructor and destructor
  TableIterator(ReportTableBase &x);        // --> adds iterator to iterator list in base class
  ~TableIterator();                         // --> removes iterator from iterator list in ReportTable class

  // Move the iterator to specific locations
  TableIterator &begin(void);
  TableIterator &end(void);

  void moveTo(int);

  // increment and decrement operators
  TableIterator operator++(int);
  TableIterator operator--(int);

  TableIterator operator+=(int);
  TableIterator operator-=(int);

  // Comparison operators
  bool operator==(const TableIterator&) const;
  bool operator!=(const TableIterator&) const;

  // Data access operators
  vector<string> *operator*();
  vector<string> *operator->();

};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class ReportTableSortHandler {

  int m_sortColumn;

  int m_sortDirection;


 public:

  ReportTableSortHandler(int col, int dir) : m_sortColumn(col), m_sortDirection(dir) {};
  ~ReportTableSortHandler() {};

  bool operator()(const vector<string> &a, const vector<string> &b);

};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class ReportTableBase {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // Table and column type holders, and iterator control

  // colHeadings contains the descriptive text that appears in the table files.
  vector<string>                  m_colHeadings;

  // Table and column type holders
  vector<ReportColumnTypeBase *>  m_colTypes;
  vector<vector<string> >         m_cells;

  // Filtered list, containing a pointer to the filtered table
  list<vector<string> *>          m_filteredTable;

  // List of iterators. Used to keep track of interators, in case of a table filter or order change
  list<TableIterator *>           m_allIterators;

  // table filename
  string  m_tableFilename;

  // Directories and data file names
  string  m_projectDir;

  // specifies if header should be rendered. Default is true.
  bool m_renderHeader;

  // specifies if borders should be rendered. Default is true.
  bool m_renderBorders;

  // filter column
  unsigned m_filterColumn;

  // ID column - used to to get IDs
  int m_idColumn;

  // Sort column - used to to sort the table
  int m_sortColumn;

  // table direction
  int m_direction;

  // table rendering exception. This is used to tell the renderer to use a specific method to render this table
  int m_renderingException;

  // specifies if the table is filtered
  bool m_filtered;


  //////////////////////////////////////////////////////////////////////////////
  // Method to change iterators, based on the iterator list, when a change occurs
  int invalidateIterators(void) {};

  //////////////////////////////////////////////////////////////////////////////
  // General filtering method. Filter a column by a string
  virtual int applyFilter(const unsigned column, const string &);


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build and edit the table

  // filename composition helper method
  virtual string composeFileName(const string &projectDir, const string &fileName);
  // edit the table
  virtual int setValue(unsigned row, unsigned col, string &value);


 public:


  // Constructors and destructor
  ReportTableBase();
  ReportTableBase(const string &projectPath, const string &tableFilename);
  // Destructor deletes colTypes vector
  virtual ~ReportTableBase();

  // add a data row to the data table
  virtual void addDataRow(vector<string> &row) {m_cells.push_back(row);};

  // find a data row in the table. Use 'column' as match column and data as comparison term. Returns row index or -1
  virtual int find(int column, string &data);

  //////////////////////////////////////////////////////////////////////////////
  // Methods to load & save data files

  // Load a table from file, given a filename
  virtual int loadTable(char separator = FILE_SEPARATOR);
  // Save a table to file, given a filename
  virtual int saveTable(char separator = FILE_SEPARATOR);
  // write the table to a stream
  virtual int writeTable(ostream &outstream, char separator = FILE_SEPARATOR, char endline = '\n');

  virtual void setTableFilename(const char *fn) {m_tableFilename = fn;};
  virtual void setTableFilename(string &fn)     {m_tableFilename = fn;};

  virtual void setProjectDir(const char *fn) {m_projectDir = fn;};
  virtual void setProjectDir(string &fn)     {m_projectDir = fn;};


  //////////////////////////////////////////////////////////////////////////////
  // Input/output methods

  // update table
  virtual int update(const ReportTableData &reportTableData);

  // get data
  virtual int getData(ostream &ss, const ReportTableData &reportTableData, char fieldSep, char rowSep);


  //////////////////////////////////////////////////////////////////////////////
  // Data comparison methods

  // compare the contents of two tables
  virtual int compare(ReportTableBase &o, vector<pair<int, int> > &res);


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build views

  // clear view
  virtual void clearView(void);
  // default view
  virtual void defineView(void)  {};
  // alternative view
  virtual void defineView2(void) {};
  // for generating images list
  virtual void defineViewImages(void) {};
  // specific view for no contigs mapping
  virtual void defineViewNoClusters(void) {};

  // add a heading
  virtual void addHeading(char *h)    {string a(h);m_colHeadings.push_back(a);};
  virtual void addHeading(string &h)  {m_colHeadings.push_back(h);};

  // method to switch off table header rendering
  virtual void noHeaders(void)        {m_renderHeader = false;};
  // query about headers rendering
  virtual bool doHeaders(void)        {return m_renderHeader;};
  // method to switch off table border rendering
  virtual void noBorders(void)        {m_renderBorders = false;};
  // query about border rendering
  virtual bool doBorders(void)        {return m_renderBorders;};

  virtual void setDirection(int st)   {m_direction = st;};
  // query about border rendering
  virtual int  getDirection(void)     {return m_direction;};

  // returns the # of elements (filtered)
  virtual int  getElementCount(void);
  virtual int  getPageSplitStats(int itemsPerPage, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames, string &fnamePrefix, string &fnameSuffix);


  // filter related methods. (Re)builds filtered table structure
  virtual int  applyFilter(const string &value)   {return applyFilter(m_filterColumn, value);};
  // define column to apply filter to, by index
  virtual void setFilterColumn(unsigned c)  {m_filterColumn = c;};

  // set the column ID (used for table sorting and to get column id, given it's index)
  virtual void setIdColumn(int id)          {m_idColumn = id;   };
  virtual string getColumnId(int index);
  // gets ID
  virtual void getId(vector<string> &IDs) {getFilteredId(IDs);};

  virtual void getFilteredId(vector<string> &IDs) {getFilteredId(IDs, m_idColumn);};
  virtual void getFilteredId(vector<string> &IDs, int column);


  // set/get the table rendering exception flag
  virtual void setRenderingException(int v) {m_renderingException = v;    };
  virtual int  getRenderingException(void)  {return m_renderingException; };

  // sort the table
  virtual void sortTable(void);
  // sort the table
  virtual void sortTable(int col)           {m_sortColumn=col; sortTable();};

  // Retrieve a columType from vector
  virtual ReportColumnTypeBase *getColType(unsigned x) const;
  // Retrieve a cell from table
  virtual string *getCell(unsigned x, unsigned y);
  // get a list of elements of a column for the entire table. This is usefull for iteration when the table is also being filtered for other purposes
  virtual vector<string> getColumn(unsigned column);

  // operations that return an iterator position for the table contents
  virtual TableIterator begin(void);
  virtual TableIterator end(void);
  // operations that return an iterator position for the column types vector
  virtual TableCtIterator ctBegin() {return m_colTypes.begin();};
  virtual TableCtIterator ctEnd()   {return m_colTypes.end();  };

  // must make iterator class friend
  friend class TableIterator;

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
