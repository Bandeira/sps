/*
 * SpecSet.h
 *
 *  Created on: Aug 30, 2011
 *      Author: jsnedecor
 */

#ifndef SPECSET_H_
#define SPECSET_H_

#include "spectrum.h"

namespace specnets
{
  class SpecSet
  {
    friend class SpecPlot;
    friend class ContPlot;
    friend class SpsPlot;
    friend class ReportTableGenerator;

  public:
    SpecSet(unsigned int sz = 0);

    /**
     * TODO: add description
     *
     *@param filename
     */
    SpecSet(char *filename);

    /**
     * TODO: add description
     *
     *@param i
     *@return
     */
    Spectrum &operator[](unsigned int i);

    const Spectrum & operator[](unsigned int i) const;

    /**
     * Copy contents from another SpecSet
     *
     *@param other
     *@return
     */
    SpecSet &operator=(const SpecSet &other);

    /**
     * Add Spectrum to back of SpecSet
     *
     *@param other
     *@return
     */
    void push_back(const Spectrum & x);

    /**
     * Insert Spectrum into SpecSet
     */
    void insert(vector<Spectrum>::iterator position,
                vector<Spectrum>::iterator start,
                vector<Spectrum>::iterator end);

    vector<Spectrum>::iterator begin(void);

    vector<Spectrum>::iterator end(void);

    /**
     * Get scan at scan number
     * @param scan_num - scan number (indexed from 1!)
     */
    Spectrum* getScan(int scan_num);

    /**
     * TODO: add description
     *
     *@return
     */
    unsigned int size() const;

    /**
     * TODO: add description
     */
    unsigned int resize(unsigned int newSize);

    /**
     * Sets the peak tolerance in classic or ppm form
     */
    void setPeakTolerance(float tolerance, bool applyPPM = false);

    /**
     * Sets the parent mass tolerance in classic or ppm form
     */
    void setParentMassTolerance(float tolerance, bool applyPPM = false);

    /**
     * Appends the set of spectra from another spectset to this one.
     */
    void appendSpecSet(SpecSet other);

    /**
     * Clear all associated psmList values in SpecSet
     *
     */
    void clearPsms(void);

    /**
     * TODO: add description
     *
     *@param results
     *@param resolution
     *@return
     */
    vector<list<int> > &getMassesHistogram(vector<list<int> > &results,
                                           float resolution = 1.0);

    /**
     * TODO: add description
     *
     *@param newResolution
     *@param enforceRounding
     */
    void setResolution(float newResolution, bool enforceRounding);

    /**
     * TODO: add description
     *
     *@param tolerance
     *@param ionOffset
     *@param includeY0k
     */
    void addZPMpeaks(float tolerance, float ionOffset, bool includeY0k);

    /**
     * TODO: add description
     *
     */
    bool maximumParsimony(void);

    /**
     * TODO: add description
     *
     *@param newTotalInt
     */
    void normalize(float newTotalInt = 1);

    /**
     * TODO: add description
     *
     *@param features
     *@param featureValue
     *@param output
     *@return
     */
    template<class T> unsigned int extract(vector<T> &features,
                                           T featureValue,
                                           SpecSet &output);

    bool LoadSpecSet_pkl_mic(const char* filename);
    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_pkl(const char *filename);

    /**
     *  LoadSpecSet_ms2: Loads a set of spectra from a .ms2 file. Spectra must be separated
     *  by at least one empty line.
     *
     *   Note: If a spectrum has more than 1 header then this function only uses the last header.
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_ms2(const char *filename);

    /**
     * LoadSpecSet_mgf: Loads a set of spectra from a .mgf file. Recognizes and processes
     * these MGF tags: BEGIN IONS, END IONS, CHARGE, PEPMASS, PEPTIDE
     *
     * Note: If a spectrum has more than 1 header then this function only uses the last header.
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_mgf(const char *filename);

    /**
     * Loads a specset with the accompanying peptide spectrum matches with input being pklbin
     *
     *@param spectra_filename
     *@param annotation_filename
     *@return
     */
    unsigned int
    LoadSpecSet_pklbin_with_annotation(const char * spectra_filename,
                                       const char * annotation_filename);

    /**
     * Loads a specset with the accompanying peptide spectrum matches with input being mgf
     *
     *@param spectra_filename
     *@param annotation_filename
     *@return
     */
    unsigned int
    LoadSpecSet_mgf_with_annotation(const char * spectra_filename,
                                    const char * annotation_filename);

    /**
     * LoadSpecSet_prms: Loads a set of spectra from a .prms file, as output but the
     * current version of pepnovo_prms (2007/04/09). Text file with multiple spectra,
     * each delimited in the following format:
     *
     * - ">> <original msms file index> <scan/index # in msms file> <single-spectrum file name>"
     *    -- scan/index # in msms file is stored in the Spectrum.scan field.
     *    -- The whole info line is stored in the Spectrum.info field.
     *  - Set of peaks in "<peak mass> <peak intensity/score>" format
     *  - Spectrum terminated by empty line
     *  - Last pair of mass/intensity corresponds to the corrected sum-of-amino-acid-masses parent mass
     *
     *  Note: If a spectrum has more than 1 header then this function only uses the last header.
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_prms(const char *filename);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    unsigned int LoadSpecSet_prmsv3(const char *filename);

    /**
     * LoadSpecSet_pklbin - loads a set of spectra in binary format. File format
     * 1 int - number of spectra in the file
     * numSpecs shorts - number of peaks per spectrum in the file.
     * arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
     *
     *@param filename Name of the pklbin file to load
     *@param countSpectraOnly If then only counts number of spectra in the file without loading the spectra (defaults to false)
     *@return
     */
    unsigned int LoadSpecSet_pklbin(const char *filename,
                                    bool countSpectraOnly = false);

    /**
     * VERSION 1: Loads all data fields of each spectrum from an open pklbin file
     * @param fp Open file pointer where the number of spectra and version of
     *    already been read from. This will be closed upon return
     * @param numSpecs Number of spectra to read. This has already been read from the file
     * @param oldVersion true if reading in the original pklbin format (pre version 1).
     * @param subversion Sub-version of this format. Currently not used (at version 1, subversion 0)
     * @return true if the SpecSet was successfully loaded, false if not. In either case the file will
     *    be closed
     */
    bool loadPklBin_1(FILE* fp, int numSpecs, bool oldVersion, char subversion);

    /**
     * VERSION 2: Loads all data fields of each spectrum from an open pklbin file
     * @param fp Open file pointer where the number of spectra and version of
     *    already been read from. This will be closed upon return
     * @param numSpecs Number of spectra to read. This has already been read from the file
     * @param subversion Sub-version of this format. Currently not used (at version 2, subversion 0)
     * @return true if the SpecSet was successfully loaded, false if not. In either case the file will
     *    be closed
     */
    bool loadPklBin_2(FILE* fp, int numSpecs, char subversion);

    /**
     * Loads all data fields of each spectrum in a pklbin file. Backwards compatible with all old versions.
     *
     *@param filename
     *@param psmFileame
     *@param peaksFilename
     *@return
     */
    int loadPklBin(const char * filename,
                   const char * psmFileame = 0x0,
                   const char * peaksFilename = 0x0);
    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_info(const char *filename);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_mgf(const char *filename, const char* activation = 0);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_pkl(const char *filename);

    /**
     * Saves the annotation_peptide field of each spectrum to the indexed line in output file
     *@param filename output filename, if a spectrum has an annotation_peptide it will be written to the same line # as its index
     *@return true if file was written successfully, false otherwise
     */
    bool SaveAnnotations(const char* filename);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_ms2(const char* filename);

    /**
     * saves the .bin component of a set of spectra in .bin format containing (per row)
     *
     * 1 int  - scan #
     * 1 int - ms level
     *
     *@param filename
     *@return
     */
    short SaveSpecSet_bin(const char *filename);

    /**
     * TODO: add description
     *
     *@param filename
     *@param filename
     *@return
     */
    short SaveSpecSet_pklbin(const char *filename, const char *binFilename =
        NULL);

    /**
     * TODO: add description
     *
     *@param filename
     *@return
     */
    void SaveScanNums(const char *filename);

    /**
     * SaveSpecSet_pklbin - saves a set of spectra in binary format. File format
     * 1 int - number of spectra in the file
     * numSpecs shorts - number of peaks per spectrum in the file
     * arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
     *
     *@param filename
     *@param specs
     *@return
     */
    short SaveSpecSet_pklbin(char *filename, vector<Spectrum *> specs);

    /**
     * TODO: add description
     *
     *@param filename
     *@param specs
     *@return
     */
    int saveMatchedProts(const char *filename);

    /**
     * TODO: Current pklbin format. Includes scan # and msLevel.
     *
     *@param filename
     *@return
     */
    int savePklBin(const char *filename,
                   const char * psmFileame = 0x0,
                   const char * peaksFilename = 0x0);

  protected:
    sps::vector<Spectrum> specs;
  };
}
#endif /* SPECSET_H_ */
