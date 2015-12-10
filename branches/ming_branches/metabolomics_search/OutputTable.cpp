/*
 * Table.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: aguthals
 */

#include "OutputTable.h"

namespace specnets
{
  OutputTable::OutputTable(void)
  {

  }

  void OutputTable::setValue(int row, int col, string value, bool header)
  {

    if (row < 0 || row > values.size()) {
      DEBUG_MSG("Setting (" << row << ", " << col << ") = " << value);
      ERROR_MSG("Invalid assignment of row " << row << " to table with size " << values.size());
      return;
    }
    if (col < 0 || col > values[row].size()) {
      DEBUG_MSG("Setting (" << row << ", " << col << ") = " << value);
      ERROR_MSG("Invalid assignment of col " << col << " to table row " << row << " with size " << values[row].size());
      return;
    }

    values[row][col].first = value;
    values[row][col].second = header;
  }

  /**
   * Prints contents of values to file
   * @param filename
   * @return true of file was written successfully, false if not
   */
  bool OutputTable::printToCSV(const char* filename)
  {
    string delim = CSV_SEP;
    FILE* output = fopen(filename, "w");

    if (output == NULL) {
      ERROR_MSG("Could not write to file " << filename);
      return false;
    }
    for (int i = 0; i < values.size(); i++) {
      for (int j = 0; j < values[i].size(); j++) {
        fprintf(output, "%s%s", values[i][j].first.c_str(), delim.c_str());
      }
      fprintf(output, "\n");
    }
    fclose(output);
    return true;
  }
}
