///////////////////////////////////////////////////////////////////////////////
#include "DataComparer.h"
///////////////////////////////////////////////////////////////////////////////
using namespace std;
namespace specnets {
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void DataCompare::setType(string &t)
{
  if(!t.compare("pklbin"))
    setType(0);
  else if(!t.compare("abruijn"))
    setType(1);
  else if(!t.compare("bin"))
    setType(2);
  else if(!t.compare("txt"))
    setType(3);
}
///////////////////////////////////////////////////////////////////////////////
// add an object
int DataCompare::addObject(string filename, int index)
{
  switch(m_type) {

    // pklbin
    case 0:
      return m_specset[index].loadPklBin(filename.c_str());
      break;

    // abruijn graph
    case 1:
      return Load_abinfo(filename.c_str(), m_abinfo[index]);
      break;

    // bin array
    case 2:
      return Load_binArray(filename.c_str(), m_binArray[index]);
      break;

    // report tables
    case 3:
      int found = filename.find_last_of("/\\");
      string path = filename.substr(0,found);
      string fn   = filename.substr(found+1);

      m_reportTable[index].setProjectDir(path);
      m_reportTable[index].setTableFilename(fn);
      m_reportTable[index].loadTable();
      break;

  }
}
///////////////////////////////////////////////////////////////////////////////
// Compare specsets
int DataCompare::compareSpecsets(void)
{
  // data holder for comparison
  SpecCompareData cd;

  // do the comparison
  int ret = m_specset[0].compare(m_specset[1], cd);

  // return 0 if 0 specs found different
  if(!ret) return 0;

  if(ret == -1) {
    cout << "Different number of spectra in both files." << endl;
    cout << "File 1 has " << cd[0].second[0] << endl;
    cout << "File 2 has " << cd[0].second[1] << endl;
    return 1;
  }

  if(m_level == 0)
    cout << "Files differ" << endl;
  else {
    cout << "Total spectra: " << m_specset[0].size() << endl;
    cout << "Different spectra: " << ret << endl;
  }

  if(m_level == 2) {
    cout << "Different spectra indexes: ";
    for(int i = 0 ; i < cd.size() ; i++)
      if(cd[i].second.size())
        cout << i << " ";
    cout << endl;
  }

  if(m_level == 3)
    for(int i = 0 ; i < cd.size() ; i++)
      if(cd[i].second.size())
        if(cd[i].second[0] == -1)
          cout << "Spectrum " << cd[i].first << " size is different: " << cd[i].second[1] << " vs " << cd[i].second[2] << endl;
        else
          cout << "Different peaks for spectrum " << cd[i].first << ": " << cd[i].second.size() << endl;

  if(m_level == 4)
    for(int i = 0 ; i < cd.size() ; i++)
      if(cd[i].second.size())
        if(cd[i].second[0] == -1)
          cout << "Spectrum " << cd[i].first << " size is different: " << cd[i].second[1] << " vs " << cd[i].second[1] << endl;
        else {
          cout << "Different peak indices for spectrum " << cd[i].first << ": ";
          for(int j = 0 ; j < cd[i].second.size() ; j++)
            cout << cd[i].second[j] << " [(" << m_specset[0][cd[i].first][cd[i].second[j]][0] <<
                                       " - " << m_specset[0][cd[i].first][cd[i].second[j]][1] <<
                                     ") ; (" << m_specset[1][cd[i].first][cd[i].second[j]][0] <<
                                       " - " << m_specset[1][cd[i].first][cd[i].second[j]][1] <<
                                       ")] ";
          cout << endl;
        }

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
// Compare abruijns
int DataCompare::compareAbruijn(void)
{
  // data holder for comparison
  abDiff_t cd;

  // do the comparison
  int ret = compareAbInfo(m_abinfo[0], m_abinfo[1], cd);

  // return 0 if 0 specs found different
  if(!ret) return 0;

  // check total number of contigs
  if(cd.size1 != cd.size2) {
    cout << "Different number of contigs in both abruijn graphs." << endl;
    cout << "File 1 has " << cd.size1 << endl;
    cout << "File 2 has " << cd.size2 << endl;
    return 1;
  } else
    cout << "Files have " << cd.size1 << " contigs." << endl;


  for(int i = 0 ; i < cd.contigs.size() ; i++) {

    // check contig indices 1-by-1
    if(cd.contigs[i].id1 != cd.contigs[i].id2) {
      cout << "Contigs at index " << i << "are different." << endl;
      cout << "For file 1, contig " << cd.contigs[i].id1 << endl;
      cout << "For file 2, contig " << cd.contigs[i].id2 << endl;
      continue;
    }

    // check number of stars
    if(cd.contigs[i].size1 != cd.contigs[i].size2) {
      cout << "Contig " << cd.contigs[i].id1 << " has a different number of stars." << endl;
      cout << "For file 1, contig " << cd.contigs[i].size1 << endl;
      cout << "For file 2, contig " << cd.contigs[i].size2 << endl;
      continue;
    }

    cout << "Contigs " << cd.contigs[i].id1 << " have " << cd.contigs[i].size1 << " stars." << endl;

    // check reversed flags
    if(cd.contigs[i].reversedList.size()) {
      cout << "Stars with different reverse flag: ";
      for(int j = 0 ; j < cd.contigs[i].reversedList.size() ; j++)
        cout << cd.contigs[i].reversedList[j] << " ";
      cout << endl;
    }

    // check for contigs with different number of nodes
    if(cd.contigs[i].nodesList.size1 != cd.contigs[i].nodesList.size2) {
      cout << "Contigs " << cd.contigs[i].id1 << " have different number of nodes." << endl;
      cout << "For file 1, nodes: " << cd.contigs[i].nodesList.size1 << endl;
      cout << "For file 2, nodes: " << cd.contigs[i].nodesList.size2 << endl;
      continue;
    }

    // check inside the nodes
    for(int j = 0 ; j < cd.contigs[i].nodesList.nodeData.size() ; j++) {

      // check vertex sizes
      if(cd.contigs[i].nodesList.nodeData[j].size1i != cd.contigs[i].nodesList.nodeData[j].size2i) {
        cout << "Nodes idx " << cd.contigs[i].nodesList.nodeData[j].idx << " have different sizes." << endl;
        cout << "For file 1, size: " << cd.contigs[i].nodesList.nodeData[j].size1i << endl;
        cout << "For file 2, size: " << cd.contigs[i].nodesList.nodeData[j].size2i << endl;
        continue;
      }

      // check peaks
      for(int k = 0 ; k < cd.contigs[i].nodesList.nodeData[j].peaks.size() ; k++) {

        // check peak stars Idxs
        if(cd.contigs[i].nodesList.nodeData[j].peaks[k].starId1 == cd.contigs[i].nodesList.nodeData[j].peaks[k].starId2) {
          cout << "Star idx in node " << cd.contigs[i].nodesList.nodeData[j].idx << " have different indexes." << endl;
          cout << "For file 1, size: " << cd.contigs[i].nodesList.nodeData[j].peaks[k].starId1 << endl;
          cout << "For file 2, size: " << cd.contigs[i].nodesList.nodeData[j].peaks[k].starId2 << endl;
          continue;
        }

        // check peak values
        if(cd.contigs[i].nodesList.nodeData[j].peaks[k].val1 == cd.contigs[i].nodesList.nodeData[j].peaks[k].val2) {
          cout << "Star idx in node " << cd.contigs[i].nodesList.nodeData[j].idx << " have different values." << endl;
          cout << "For file 1, size: " << cd.contigs[i].nodesList.nodeData[j].peaks[k].val1 << endl;
          cout << "For file 2, size: " << cd.contigs[i].nodesList.nodeData[j].peaks[k].val2 << endl;
          continue;
        }

      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
// Compare bin arrays
int DataCompare::compareBinArray(void)
{
  // data holder for comparison
  vector<pair<int,int> > cd;

  // do the comparison
  int ret = Compare_binArray(m_binArray[0], m_binArray[1], cd);

  // return 0 if 0 specs found different
  if(!ret) return 0;

  // check bin array length
  if(ret == -1) {
    cout << "Bin arrays have different number of rows." << endl;
    cout << "Bin array 1 has " << cd[0].first << endl;
    cout << "Bin array 2 has " << cd[0].second << endl;
    return -1;
  }

  // check for empty arrays.
  if(ret == -2) {
    cout << "Bin arrays are empty." << endl;
    return -1;
  }

  // check bin array width
  if(ret == -3) {
    cout << "Bin arrays have different number of columns." << endl;
    cout << "Bin array 1 has " << cd[0].first << endl;
    cout << "Bin array 2 has " << cd[0].second << endl;
    return -1;
  }

  // check for empty rows.
  if(ret == -4) {
    cout << "Bin arrays have empty rows." << endl;
    return -1;
  }

  cout << ret << " differences found." << cd.size() << endl;

  // output differences
  for(int i = 0 ; i < cd.size() ; i++) {
    cout << "row " << cd[i].first << ", column " << cd[i].second << endl;
    cout << "  < : " << m_binArray[0][cd[i].first][cd[i].second] << endl;
    cout << "  > : " << m_binArray[1][cd[i].first][cd[i].second] << endl;    
  }


  return ret;
}
///////////////////////////////////////////////////////////////////////////////
// Compare report tables
int DataCompare::compareReportTables(void)
{
  // data holder for comparison
  vector<pair<int,int> > cd;

  // do the comparison
  int ret = m_reportTable[0].compare(m_reportTable[1], cd);

  // return 0 if 0 specs found different
  if(!ret) return 0;

  // check bin array length
  if(ret == -1) {
    cout << "Report tables have different number of rows." << endl;
    cout << "Report table 1 has " << cd[0].first << endl;
    cout << "Report table 2 has " << cd[0].second << endl;
    return -1;
  }

  // check for empty arrays.
  if(ret == -2) {
    cout << "Report tables are empty." << endl;
    return -1;
  }

  // check bin array width
  if(ret == -3) {
    cout << "Report tables have different number of columns." << endl;
    cout << "Report table 1 has " << cd[0].first << endl;
    cout << "Report table 2 has " << cd[0].second << endl;
    return -1;
  }

  // check for empty rows.
  if(ret == -4) {
    cout << "Report tables have empty rows." << endl;
    return -1;
  }

  cout << ret << " differences found." << cd.size() << endl;

  // output differences
  for(int i = 0 ; i < cd.size() ; i++) {
    cout << "row " << cd[i].first << ", column " << cd[i].second << endl;
    cout << "  < : " << m_reportTable[0].getCell(cd[i].first, cd[i].second) << endl;
    cout << "  > : " << m_reportTable[1].getCell(cd[i].first, cd[i].second) << endl;
  }

  return ret;
}
///////////////////////////////////////////////////////////////////////////////
// execute the comparison
int DataCompare::exec(void)
{
  switch(m_type) {

    // pklbin
    case 0:
      return compareSpecsets();
      break;

    // abruijn
    case 1:
      return compareAbruijn();
      break;

    // binArray
    case 2:
      return compareBinArray();
      break;

    // report tables
    case 3:
      return compareReportTables();
      break;

  }
}
///////////////////////////////////////////////////////////////////////////////
};
///////////////////////////////////////////////////////////////////////////////
