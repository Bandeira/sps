// Header Include
#include "Logger.h"
#include "SpectrumPairSet.h"
#include <algorithm>

using namespace std;

namespace specnets
{


bool SpectrumPairComparator (SpectrumPair sp1, SpectrumPair sp2)
{
	float val1 = sp1.score1 + sp1.score2;
	float val2 = sp2.score1 + sp2.score2;
	return (val1 >= val2);
}

	//---------------------------------------------------------------------------
  SpectrumPairSet::SpectrumPairSet(void)
  {
    resize(0);
  }

  //---------------------------------------------------------------------------
  SpectrumPairSet::SpectrumPairSet(unsigned int nsize)
  {
    thePairs.resize(nsize);
  }

  //---------------------------------------------------------------------------
  unsigned int SpectrumPairSet::size(void) const
  {
    return thePairs.size();
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::resize(unsigned int newSize,
                               SpectrumPair newPairExemplar)
  {
    thePairs.resize(newSize, newPairExemplar);
  }

  //---------------------------------------------------------------------------
  SpectrumPair const & SpectrumPairSet::operator[](unsigned int index) const
  {
    return thePairs[index];
  }

  //---------------------------------------------------------------------------
  SpectrumPair & SpectrumPairSet::operator[](unsigned int index)
  {
    return thePairs[index];
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::sort_pairs()
  {
	  // using function as comp
	    sort (thePairs.begin(), thePairs.end(), SpectrumPairComparator);
  }
  //---------------------------------------------------------------------------
  int SpectrumPairSet::loadFromBinaryFile(const std::string & filename)
  {
    FILE *fp;
    unsigned int numEntries, vCount, resIdx, vIdx;
    fp = fopen(filename.c_str(), "r");
    if (fp == 0)
    {
      ERROR_MSG("Can not open: " << filename);
      return -1;
    }

    fread(&numEntries, sizeof(unsigned int), 1, fp); // Number of entries in the file
    resize(numEntries);
    if (numEntries == 0)
    {
      fclose(fp);
      ERROR_MSG("No entries in file: " << filename);
      return 0;
    }
    vCount = operator[](0).loadSz();
    float *data = (float *)new float[vCount];
    int read_result;
    vector<float> dataV(vCount);
    for (unsigned int resIdx = 0; resIdx < numEntries; resIdx++)
    {
      read_result = fread(data, sizeof(float), vCount, fp);
      if (read_result != vCount)
      {
        resize(0);
        ERROR_MSG("Can not read " << filename);
        return -2;
      }
      for (vIdx = 0; vIdx < 2; vIdx++)
      {
        dataV[vIdx] = data[vIdx] - 1; // Spectrum indices are 1-based
      }
      for (vIdx = 2; vIdx < vCount; vIdx++)
      {
        dataV[vIdx] = data[vIdx];
      }
      operator[](resIdx).load(dataV);
    }
    delete[] data;
    fclose(fp);
    return numEntries;
  }

  //---------------------------------------------------------------------------
  bool SpectrumPairSet::saveToBinaryFile(const std::string & filename)
  {
    FILE *fp;
    vector<float> dataV;
    float data[6];

    fp = fopen(filename.c_str(), "w");
    if (fp == 0)
    {
      ERROR_MSG("Can not open: " << filename);
      return false;
    }
    unsigned int numEntries = thePairs.size();
    fwrite(&numEntries, sizeof(unsigned int), 1, fp); // Number of entries in the file
    for (unsigned int i = 0; i < numEntries; i++)
    {
      thePairs[i].serialize(dataV);
      for (unsigned int j = 0; j < 2; j++)
      {
        data[j] = dataV[j] + 1; // Spectrum indices are 1-based
      }
      for (unsigned int j = 2; j < dataV.size(); j++)
      {
        data[j] = dataV[j];
      }
      fwrite(data, sizeof(float), dataV.size(), fp);
    }
    fclose(fp);
    return true;
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::push_back(const SpectrumPair & newPair)
  {
    thePairs.push_back(newPair);
  }

} // namespace specnets
