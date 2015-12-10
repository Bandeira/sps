////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_COLUMN_TYPES_H__
#define __REPORT_COLUMN_TYPES_H__
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
#include <list>

using namespace std;

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
class ReportParamsOption {

 public:

  string param;
  string option;
  string validator;

  ReportParamsOption(string p, string o, string v) :
    param(p), option(o), validator(v)
    {};

};
////////////////////////////////////////////////////////////////////////////////
class ReportColumnTypeBase {

 public:

  string  cssClass;       // CSS class name for HTML formating
  bool    leftDivider;    //
  bool    dynamic;        // specifies a cell with content dynamically used (content sent to server)
  string  columnLabel;    // label on the table column
  string  link;           // URL template for the link
  string  onClick;        // template for onClick HTML method
  string  id;             // template for field ID, needed to read or write data to
  string  validator;      // validator must be not null in order for the cell to be displayed
  int     displayLevel;   // display level. If input level < display level, then the item is not displayed

  // constructor used to initialize methods
  ReportColumnTypeBase() : leftDivider(false), dynamic(false), displayLevel(0) {};

  // virtual destructor to make class polymorfic
  virtual ~ReportColumnTypeBase() {};

};
////////////////////////////////////////////////////////////////////////////////
class ReportColumnTypeImageOnDemand : public ReportColumnTypeBase {

 public:

  string                      iconRenderer;       // renderer used to generate the icon. If empty, iconParams treated as image/URL
  //string                      iconParams;         // Icon path/image/URL
  vector<ReportParamsOption>  iconParams;         // Icon path/image/URL
  int                         iconDisplayLevel;   // display level for the icon. If input level < display level, then the item is not displayed
  string                      alt;              // alternative to display in place of icon


  string                      label;            // label to show for the link (defined by renderer and params)
  string                      renderer;         // Object name used for rendering the Image On Demand
  vector<ReportParamsOption>  params;           // parameters passed to the renderer object for  the Image On Demand
  int                         linkDisplayLevel; // display level for the link. If input level < display level, then the item is not displayed

  // When using a CGI call, the command is constructed in the following way:
  // /cgi-bin/<renderer> <params>
  //
  // when rendering local static pages, the renderer name is used to generate/request a render object by name (using a object factory model)
  // and <params> are passed to build the image

  ReportColumnTypeImageOnDemand() : iconDisplayLevel(0), linkDisplayLevel(0) {}

};
////////////////////////////////////////////////////////////////////////////////
class ReportColumnTypeString : public ReportColumnTypeBase {

 public:

  string  text;     // Text template for cell contents, button, input box.
  bool    isButton; // If true, a button is drawn with the text in the "text" field
  bool    isInput;  // if True, an input box is drawn


  // constructor used to initialize methods
  ReportColumnTypeString() : isButton(false), isInput(false) {};

};
////////////////////////////////////////////////////////////////////////////////
class ReportColumnTypeStringMultiple : public ReportColumnTypeBase {

 public:

  string  linkPrefix; // link filename prefix
  string  linkSuffix; // link filename suffix
  string  text;       // Text template for cell contents.

  // constructor used to initialize methods
  ReportColumnTypeStringMultiple() {};

};
////////////////////////////////////////////////////////////////////////////////
class ReportColumnTypeBox : public ReportColumnTypeBase {

 public:

  vector<ReportColumnTypeBase *> sequences; // Vector of several column types.

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
