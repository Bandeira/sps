////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererHtml.h"
#include "Tokenizer.h"
#include "copyright.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// Page component generation methods
////////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::buildNavigationBar(string &navBar, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames)
{
  navBar = "<div>";

  for(int i = 0 ; i < IDs.size() ; i++, navBar += "   ") {

    navBar += "<a href='";
    navBar += fNames[i];
    navBar += "'>[";
    navBar += IDs[i];
    navBar += '-';
    navBar += IDsEnd[i];
    navBar += "]</a>";
  }
  navBar += "</div>";
}
////////////////////////////////////////////////////////////////////////////////
// Page prolog and epilogue
////////////////////////////////////////////////////////////////////////////////
// prolog
int ReportRendererHtml::renderProlog(ReportBase *table, ostream &outstream)
{
  outstream << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"\n";
  outstream << "  \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n";
  outstream << "<HTML xmlns=\"http://www.w3.org/1999/xhtml\">\n";

  outstream << "<head>\n";
  outstream << "  <meta http-equiv=\"Content-Type\" content=\"text/shtml; charset=ISO-8859-1\" />\n";
//  outstream <<  "  <meta http-equiv=\"PRAGMA\" content=\"NO-CACHE\" />\n";
//  outstream <<  "  <meta http-equiv=\"CACHE-CONTROL\" content=\"NO-STORE, NO-CACHE, MUST-REVALIDATE, POST-CHECK=0, PRE-CHECK=0\" />\n";
//  outstream <<  "  <meta http-equiv=\"EXPIRES\" content=\"01 Jan 1970 00:00:00 GMT\" />\n";
//  outstream <<  "  <meta http-equiv=\"Last-Modified\" content=\"01 Jan 1970 00:00:00 GMT\" />\n";
//  outstream <<  "  <meta http-equiv=\"If-Modified-Since\" content=\"01 Jan 1970 00:00:00 GMT\" />\n";

  outstream << "  <title>UCSD Computational Mass Spectrometry Website</title>\n";
  outstream << "  <link href=\"styles/main.css\" rel=\"stylesheet\" type=\"text/css\" />\n";
  outstream << "  <link rel=\"shortcut icon\" href=\"images/favicon.ico\" type=\"image/icon\" />\n";
  outstream << "  <script src=\"scripts/util.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  outstream << "  <script src=\"scripts/download.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  outstream << "  <script src=\"scripts/render.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  outstream << "  <script src=\"scripts/inspect.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";

  outstream << "  <script type='text/javascript' src='js/prototype.js'></script>\n";
  outstream << "  <script type='text/javascript' src='js/scriptaculous.js?load=effects,builder'></script>\n";
  outstream << "  <script type='text/javascript' src='js/lightbox.js'></script>\n";

  outstream << "  <link rel='stylesheet' href='css/lightbox.css' type='text/css' media='screen' />\n";

  //outstream << "  <script src=\"js/mootools.js\" type=\"text/javascript\"></script>\n";
  //outstream << "  <script src=\"js/slimbox.js\" type=\"text/javascript\"></script>\n";
  //outstream << "  <link href=\"css/slimbox.css\" rel=\"stylesheet\" type=\"text/css\" media=\"screen\" />\n";

  outstream << "  <script src=\"js/sorttable.js\"></script>\n";



/*
  outstream << "  <script type=\"JavaScript\">\n";
  outstream << "  <!--\n";
  outstream << "  function myModalWin(filename, name)\n";
  outstream << "  {\n";
  outstream << "    if (! name)\n";
  outstream << "      name = \"name\";\n";

  outstream << "    if (window.showModalDialog)\n";
  outstream << "    {\n";
  outstream << "      alert('#1');\n";
  outstream << "      window.showModalDialog(filename, name, \"dialogWidth:800px; dialogHeight:600px\");\n";
  outstream << "    }\n";
  outstream << "    else\n";
  outstream << "    {\n";
  outstream << "      alert('#2');\n";
  outstream << "      window.open(filename, name, 'width=800, height=600, toolbar=no, directories=no, status=no, menubar=no, scrollbars=no, resizable=no, modal=yes');\n";
  outstream << "    }\n";
  outstream << "  }\n";

  outstream << "  function modalForm(form, name)\n";
  outstream << "  {\n";
  outstream << "    if (! name)\n";
  outstream << "      name = \"name\";\n";

  outstream << "    if (window.focus)\n";
  outstream << "    {\n";
  outstream << "      myModalWin(\"\", name);\n";
  outstream << "      form.target = name;\n";
  outstream << "    }\n";
  outstream << "  }\n";
  outstream << "  //-->\n";
  outstream << "  </script>\n";

  outstream << "  <script src=\"js/pager.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  outstream << "  <noscript>\n";
  outstream << "    <meta http-equiv=\"refresh\" content=\"10\">\n";
  outstream << "  </noscript>\n";

  outstream << "  <script language=\"JavaScript\">\n";
  outstream << "  <!--\n";
  outstream << "  var sURL = unescape(window.location.pathname);\n";
  outstream << "  function doLoad()\n";
  outstream << "  {\n";
  outstream << "      setTimeout( \"refresh()\", 10*1000 );\n";
  outstream << "  }\n";

  outstream << "  function refresh()\n";
  outstream << "  {\n";
  outstream << "      window.location.href = sURL;\n";
  outstream << "  }\n";
  outstream << "  //-->\n";
  outstream << "  </script>\n";

  outstream << "  <script language=\"JavaScript1.1\">\n";
  outstream << "  <!--\n";
  outstream << "  function refresh()\n";
  outstream << "  {\n";
  outstream << "      window.location.replace( sURL );\n";
  outstream << "  }\n";
  outstream << "  //-->\n";
  outstream << "  </script>\n";

  outstream << "  <script language=\"JavaScript1.2\">\n";
  outstream << "  <!--\n";
  outstream << "  function refresh()\n";
  outstream << "  {\n";
  outstream << "      window.location.reload(true);\n";
  outstream << "  }\n";
  outstream << "  //-->\n";
  outstream << "  </script>\n";

  outstream << "  <script language=\"JavaScript\">\n";
  outstream << "  <!--\n";
  outstream << "  function showImage(img)\n";
  outstream << "  {\n";
  outstream << "    var aux = 'data:image/png;base64,';\n";
  outstream << "    aux += img;\n";
  outstream << "      window.showModalDialog(aux);\n";
  outstream << "      myModalWin(aux, \"modalWindow\");\n";
  outstream << "  }\n";
  outstream << "  //-->\n";
  outstream << "  </script>\n";
*/
//   outstream << "  <script language=\"JavaScript\">\n";
//   outstream << "  <!--\n";
//   outstream << "  function needReload()\n";
//   outstream << "  {\n";
//   outstream << "      var loc=window.location.href.toString()\n";
//   outstream << "      if(loc.indexOf(\"?\") == -1)\n";
//   outstream << "      {\n";
//   outstream << "          location.replace(loc+\"?1\");\n";
//   outstream << "      }\n";
//   outstream << "      else\n"
//   outstream << "      {\n";
//   outstream << "          return false;\n";
//   outstream << "      }\n";
//   outstream << "  }\n";
//   outstream << "  //-->\n";
//   outstream << "  </script>\n";

//   outstream << "  <script language=\"JavaScript1.2\">\n";
//   outstream << "  <!--\n";
//   outstream << "  function refreshScroll()\n";
//   outstream << "  {\n";
//   outstream << "    sp = document.body.scrollTop;\n";
//   outstream << "    refresh();\n";
//   outstream << "    window.scrollTo(0, \"+sp+\");\n";
//   outstream << "  }\n";
//
//   outstream << "  function doLoadScroll()\n";
//   outstream << "  {\n";
//   outstream << "      setTimeout( \"refreshScroll()\", 10*1000 );\n";
//   outstream << "  }\n";
//
//   // function saves scroll position
//   outstream << "  function fScroll(val)\n";
//   outstream << "  {\n";
//   outstream << "      var hidScroll = document.getElementById('hidScroll');\n";
//   outstream << "      hidScroll.value = val.scrollTop;\n";
//   outstream << "  }\n";
//
//   // function moves scroll position to saved value
//   outstream << "  function fScrollMove(what)\n";
//   outstream << "  {\n";
//   outstream << "      var hidScroll = document.getElementById('hidScroll');\n";
//   outstream << "      document.getElementById(what).scrollTop = hidScroll.value;\n";
//   outstream << "  }\n";
//   outstream << "  //-->\n";
//   outstream << "  </script>\n";


  outstream << "<style>\n";

  outstream << "* {margin-top: 0; margin-bottom: 0;}\n";
  outstream << "html, body {height: 100%;}\n";
  outstream << ".wrapper {min-height: 100%; height: auto !important; height: 100%;margin: 0 auto -4.5em;}\n";
  outstream << ".footer, .push {height: 4.5em;}\n";

  outstream << "TD.ln {	height: 1px; background-color: #0055ff; PADDING: 0pt; MARGIN:	0pt; line-height: 1px;}\n";
  outstream << "TD.HSep 	{width: 	20pt;}\n";
  outstream << "TD.VHSep 	{height: 	20pt;}\n";


  outstream << "</style>\n";


  outstream << "</head>\n";
  outstream << "<body>\n";

  outstream << "<div class='wrapper' align='center'>";

  // Logo image
  outstream << "<div><h3 align=center><img src='data:image/jpg;base64,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' /></h3></div><br /><br />";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// epilog
int ReportRendererHtml::renderEpilog(ReportBase *table, ostream &outstream)
{
  outstream << "<div class='push'></div>";
  outstream << "</div>";

  outstream << "<div class='footer'><table width='100%'><tr><td class='VHSep'></td></tr><tr><td class='HSep'></td><td class='ln'></td><td class='HSep'></td></tr><tr><td class='HSep'></td><td class='Footer'>";
  outstream << REPORT_FOOTER_TEXT;
  outstream << "</td></tr></table></div>";

  outstream << "</body>\n";
  outstream << "</html>\n";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// between tables
int ReportRendererHtml::renderInterTable(ReportBase *table, ostream &outstream)
{
  // output a line break in a table row
  outstream << "<tr><td><br></td></tr>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering exceptions
///////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::breakProteinIntoChunks(vector<string> &in, int &count)
{
  count = 0;
  for(int i = 0 ; i < in.size() ; i++) {
    // string container
    string aux;
    for(int j = 0 ; j < in[i].length() ; j++) {
      // line break;
      if(count % 50 == 0) {
        aux += "<br>";
      // spacer
      } else if(count % 10 == 0) {
        aux += " ";
      }
      // copy element
      aux += in[i][j];
      count++;
    }
    in[i] = aux;
  }
}
///////////////////////////////////////////////////////////////////////////////
void ReportRendererHtml::colorProteinString(vector<string> &in, string &out)
{
  for(int i = 0 ; i < in.size() ; i++) {
    if(in[i].size()) {
      out += "<font color='";
      if(i % 2) {
        out += "0";
      } else {
        out += "#AAAAAA";
      }
      out += "'>";
      out += in[i];
      out += "</font>";
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream)
{
  // declare and initialize table row interator, to allow direct access to the table
  TableIterator ti = table->begin();

  // get table data
  string job          = (*(*ti))[0];
  string user         = (*(*ti))[1];
  string status       = (*(*ti))[2];
  string elapsed      = (*(*ti))[3];
  string log          = (*(*ti))[4];
  string contig       = (*(*ti))[5];
  string protein      = (*(*ti))[6];
  string cluster      = (*(*ti))[7];
  string fileStr      = (*(*ti))[8];
  string indicesStr   = (*(*ti))[9];

  vector<string> files, indices;
  // split indices and filenames into strings
  stringSplit2(fileStr, files, "|");
  stringSplit2(indicesStr, indices, "|");

  // output tables

  outstream << "<div id='bodyWrapper'>";
  outstream << "<div id='textWrapper'>";
  outstream << "<table class='mainform'>";
  outstream << "<tr> ";
  outstream << "    <th colspan='0'>Job Status</th>";
  outstream << "  </tr>";
  outstream << "  <tr><td>";
  outstream << "<table class='sched' width='100%'>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Job</span></th>";
  outstream << "    <td>" << job <<"</td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>User</span></th>";
  outstream << "    <td>" << user << "</td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Status</span></th>";
  outstream << "    <td style='background-color:lightgreen;'>" << status << " %</td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Elapsed</span></th>";
  outstream << "    <td>" << elapsed << "</td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Log</span></th>";
  outstream << "    <td><a href='spsplot.log'>" << log << "</a></td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399' rowspan='2'><span style='color:white'>Data</span></th>";
  outstream << "    <td><a href='" << contig << ".0.html'>Group by Contig</a></td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <td><a href='" << protein << ".html'>Group by Protein</a></td>";
  outstream << "  </tr>";
  outstream << "  <tr>";
  outstream << "    <th width='25%' bgcolor='#003399'><span style='color:white'>Cluster Data</span></th>";
  outstream << "    <td>";
  outstream << "      <a href='" << cluster << ".txt'>All Clusters (txt)</a><br>";

  for(int i = 0 ; i < files.size() && i < indices.size() ; i++)
    outstream << "      <a href='spectra." << indices[i] << ".0.html'>Group by <i>" << files[i] << "</i></a><br>";

  outstream << "    </td>";
  outstream << "  </tr>";
  outstream << "  </td></tr></table></td></tr>";
  outstream << "  <tr>";
  outstream << "    <td colspan='0' class='bottomline'>&nbsp;</td>";
  outstream << "  </tr>";
  outstream << "</table>";
  outstream << "</div></div>";
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::renderTableExceptionProteinHeader(ReportTableBase *table, ostream &outstream)
{
  // declare and initialize table row interator, to allow direct access to the table
  TableIterator ti = table->begin();

  // get the protein name, at indice 1
  string protienId    = (*(*ti))[0];
  string proteinName  = (*(*ti))[1];
  string contigs      = (*(*ti))[3];
  string spectra      = (*(*ti))[4];
  string aas          = (*(*ti))[5];
  string coverage     = (*(*ti))[6];
  string sequence     = (*(*ti))[7];

  // remove | from protein name
  for(int i = 0 ; i < proteinName.size() ; i++)
    if(proteinName[i] == '|')
      proteinName[i] = ' ';

  // format protein sequence
  int count;
  vector<string> breaked;
  string coloredProtein;
  stringSplit2(sequence, breaked, "|");
  breakProteinIntoChunks(breaked, count);
  colorProteinString(breaked, coloredProtein);

  // build legend for protein sequence
  string legend;
  int current = 1;
  while(current <= count) {
    legend += "<br>";
    legend += parseInt(current);
    current += 50;
  }

  // protein header HTML code
  outstream << "<table class='result' width='100%' style='border-spacing: 0px;'>";
  outstream << "<tr>";
  outstream << "<td colspan='0'><h2><i>" << proteinName << "</i></h2>";
  outstream << "<hr><b>" << contigs << " contigs, " << spectra << " spectra, " << aas << " amino acids, " << coverage << " coverage" << "</b></td>";
  outstream << "<td></td>";
  outstream << "</tr>";
  outstream << "<tr> ";
  outstream << "<td align='right'><tt>" << legend << "</tt></td>";
  outstream << "<td><tt>" << coloredProtein << "</tt></td>";
  outstream << "</tr>";
  outstream << "</table>";

  // link to protein detains page
  outstream << "<div>";
  outstream << "<a href='protein_details." << protienId << ".html'>Protein coverage</a>";
  outstream << "</div>";
  outstream << "<br>";

    // return status OK
  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Table rendering methods
///////////////////////////////////////////////////////////////////////////////
  // render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
int ReportRendererHtml::renderTable(ReportTableBase *table, ostream &outstream)
{
  switch(table->getRenderingException()) {

  case RENDERING_EXCEPTION_NONE:
    break;

  case RENDERING_EXCEPTION_PROTEIN_HEADER:
    return renderTableExceptionProteinHeader(table, outstream);
    break;

  case RENDERING_EXCEPTION_MAIN_PAGE:
    return renderTableExceptionMainPage(table, outstream);

  default:
    break;
  }


  // HTML class definition
  string aux = "class='result sortable'";
  // boder definition
  if(!table->doBorders())
    aux = "border='0'";
  // HTML table start
  outstream << "<table " << aux << ">";
  // render contents
  ReportRendererBase::renderTable(table, outstream);
  // HTML table terminator
  outstream << "</table>";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// render a header row method. Renders a table header row.
int ReportRendererHtml::renderTableHeaderRow(vector<string> &row, ostream &outstream)
{
  // render the header
  outstream << "<tr>";
  ReportRendererBase::renderTableHeaderRow(row, outstream);
  outstream << "</tr>";

  // return status OK
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
int ReportRendererHtml::renderTableRow(vector<string> &row, ostream &outstream)
{
  outstream << "<tr align='center' valign='middle'>";
  ReportRendererBase::renderTableRow(row, outstream);
  outstream << "</tr>";

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Building methods
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// build header cell comon method. Renders all ColumnTypes
int ReportRendererHtml::buildTableHeaderCell(ReportColumnTypeBase *base, stringstream &ss)
{
  // cell begin HTML tag with class
  ss << "<th><span style='color:white'>";
  // cell content
  ss << base->columnLabel;
  // cell end HTML tag
  ss << "</span></th>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// render cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
int ReportRendererHtml::buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &ss)
{
  // don't render cells with dynamic content
  if(base->dynamic)
    return OK;

  // check if this cell is a StringMultiple. If it is, we shouldn't issue <a> tags here.
  ReportColumnTypeStringMultiple *auxM = dynamic_cast<ReportColumnTypeStringMultiple*>(base);

  // auxiliary variables
  stringstream cls, link;

  // gather needed attributes
  if(base->link.size() && (auxM == NULL))
    link << parseTemplates(base->link, row); // parseTemplatesAll

  if(base->cssClass.size())
    cls << " class='" << base->cssClass << "'";


  // cell begin HTML tag with class
  ss << "<td align='center' valign='middle'" << cls.str() << ">";
  // Link section
  if(base->link.size()  && (auxM == NULL))
    ss << "<a href='" << link.str() << "'>";
  // process base class cell renderer
  ReportRendererBase::buildTableCell(base, row, ss);
  // Link section
  if(base->link.size()  && (auxM == NULL))
    ss << "</a>";
  // cell end HTML tag
  ss << "</td>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Renderers for the diferent cell types
int ReportRendererHtml::buildCellImageOnDemand(ReportColumnTypeImageOnDemand *ct, vector<string> *row, stringstream &ss)
{
  // auxiliary variables
  stringstream icon, label, url;

  //  cout << ct->iconDisplayLevel << " < " << m_displayLevel << " -- " << ct->iconRenderer << endl;

  // Icon path/image
  if(ct->iconParams.size()) {

    // result holder
    string image64;
    // get parsed parameters / URL
    vector<string> pars; // = parseTemplatesAll(ct->iconParams, row);
    parseParamsVector(pars, ct->iconParams, row);
    // if there is a renderer, use it to render the image
    if(ct->iconRenderer.size()  && ct->iconDisplayLevel < m_displayLevel) {
      //renderImage(ct->iconRenderer, parseTemplatesAll(pars, row));
      renderImage(ct->iconRenderer, pars);
      image64 = "data:image/png;base64,";
      image64 += m_image;
      ss << "<img src='" << image64 << "' />";
    } else
      ss << "<p>" << ct->alt <<"</p>";
  }


  // label to show for the link
  if(ct->label.size())
    label << parseTemplates(ct->label, row); // parseTemplatesAll
  // URL template to be used to get the image
  if(ct->renderer.size() && ct->linkDisplayLevel < m_displayLevel) {
    //url << ct->renderer << ' ' <<  parseTemplatesAll(ct->params, row);
    //renderImage(ct->renderer, parseTemplatesAll(ct->params, row));
    //outstream << "<img src='data:image/png;base64," << m_image << "' />";
    vector<string> pars;
    parseParamsVector(pars, ct->params, row);
    renderImage(ct->renderer, pars);

    //outstream << "<a href='" << url.str() << "'>";
    string image64 = "data:image/png;base64,";
    image64 += m_image;
    //outstream << "<a href='#' onClick='JavaScript: var a=\"" << m_image << "\";showImage(a);'>";
    ss << "<a href=\"" << image64 << "\" rel=\"lightbox\">";
  }

  // button
  if(ct->label.size())
    ss << label.str();
    //outstream << "<div onClick='javascript:alert(\"#1\");var a='" << m_image << "';showImage(a);'>" << label.str() << "</div>";

  // url for IOD - terminator
  if(ct->renderer.size())
    ss << "<a>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::buildCellString(ReportColumnTypeString *ct, vector<string> *row, stringstream &ss)
{
  // auxiliary variables
  stringstream onclick, id, text;

  // gather needed attributes
  if(ct->text.size())
    text << parseTemplates(ct->text, row); // parseTemplatesAll

  if(ct->onClick.size())
    onclick << " onclick='" << parseTemplates(ct->onClick, row) << "'"; // parseTemplatesAll

  if(ct->id.size())
    id << " ID='" << parseTemplates(ct->id, row) << "'"; // parseTemplatesAll


  // edit box
  if(ct->isInput) {
    ss << "<input type='text' style='text-transform: uppercase; width:100%'" << id.str() << onclick.str() << " />";
  }

  // button
  else if(ct->isButton) {
    ss << "<input type='button' value='" << text.str() <<"'" << onclick.str() << id.str() << "' />";

  // regular text
  } else {
    ss << "<p" << id.str() << onclick.str() << ">" << text.str() << "</p>";
  }

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Method to process 'multiple' table cells, to have several lines of data in the same table cell
int ReportRendererHtml::buildCellStringMultiple(ReportColumnTypeStringMultiple *ct, vector<string> *row, stringstream &ss)
{
  // auxiliary variables
  string text, link;
  vector<string> links, texts;

  // gather needed attributes
  if(ct->link.size())
    link = parseTemplates(ct->link, row); // parseTemplatesAll

  if(ct->text.size())
    text = parseTemplates(ct->text, row); // parseTemplatesAll

  StringExplode(link, links, true, "|");
  StringExplode(text, texts, true, "|");

  for(int i = 0 ; i < texts.size() ; i++) {

    // the ith link
    if(links.size() > i)
      ss << "<a href='" << ct->linkPrefix << links[i] << ct->linkSuffix << "'>";

    // regular text
    ss << "<p>" << texts[i] << "</p>";

    // the ith link
    if(links.size() > i)
      ss << "</a>";

    // line break
    if(i < texts.size()-1)
      ss << "<br>";
  }

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::buildCellBox(ReportColumnTypeBox *ct, vector<string> *row, stringstream &ss)
{
  // sequences box begin sequence
  ss << "<table align='center' border='0'><tr align='center'><td>";

  // sequences box cycle thru all cells
  vector<ReportColumnTypeBase *>::iterator it = ct->sequences.begin();
  for(bool begin = true; it != ct->sequences.end() ; it++) {

    // skip dynamic
    if( (*it)->dynamic )
      continue;

    // validate cell
    if( (*it)->validator.length()) {
      string aux = parseTemplates((*it)->validator, row); // parseTemplatesAll
      if(!aux.length()) {
        continue;
      }
    }

    // if there is a column label
    if((*it)->columnLabel.size()) {
      //new line, if it is not the first one
      if(!begin)
        ss << "<br><br>";
      // the column label, in bold
      ss << "<b>" << (*it)->columnLabel << "</b>";
      //new line
      ss << "<br>";
      // subsequent new lines between items will be inserted
      begin = false;
    }

    stringstream link;

    // gather needed attributes
    if((*it)->link.size())
      link << parseTemplates((*it)->link, row); // parseTemplatesAll

    // Link section
    if((*it)->link.size())
      ss << "<a href='" << link.str() << "'>";

    // call base class to render cell
    ReportRendererBase::buildTableCell(*it, row, ss);

    // Link section
    if((*it)->link.size())
      ss << "</a>";
  }

  // sequences box end sequence
  ss << "</td></tr></table>";

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportRendererHtml::renderImage(const string &object, vector<string> & params)
{
  //cout << object << " : " << paramList << endl;

  // clear image container
  m_image.clear();

  // Get the exemplar from the factory
  ReportModuleBase * module = ReportModuleFactory::getModule(object);
  if (module == 0)
    return ERROR;

  // get a clean module
  ReportModuleBase *moduleReal = module->clone();

  // get params as a vector
  //vector<string> params;
  //StringExplode(paramList, params, true);
  //stringSplit(paramList, params);

  // build parameter list for module
  int count = params.size();
  char *aux[count+2];

  // Fill parameters structure
  for(int i = 0 ; i < params.size() ; i++)
    aux[i+1] = (char *)params[i].c_str();

  // add executable location/name
  string exeAux = m_exeDir; exeAux += '/' ; exeAux += object;
  aux[0] = (char *)exeAux.c_str();
  // set terminator
  aux[count+1] = 0;

  // invoke
  try {
    moduleReal->invoke(count+1, aux);
  }
  catch(...) {
    return ERROR;
  }

  // get the generated image
  moduleReal->getData(m_image);

  // delete the object
  delete moduleReal;

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
