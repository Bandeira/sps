#ifndef PRUNEENVELOPES_H
/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 8/05/09
 */
#define PRUNEENVELOPES_H

#include <ctype.h>
#include <cstdio>
#include <cmath>
#include <list>
#include <string>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <set>
#include <map>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include "../ion.h"
#include "../label.h"
#include "../range.h"
#include "../spectrum.h"
#include "../aminoacid.h"
#include "../dbg_print.h"
#include "../twovalues.h"
#include "../db_fasta.h"
#include "../utils.h"
#include "../batch.h"
#include "../inputParams.h"
#include "../setmerger.h"
#include "../ms1.h"

using namespace std;

/**
 * Determines the charge of masses in MS/MS spectra as explained by their isotopic
 * envelopes and converts them to monoisotopic masses with the summed intensity
 * of their isotopic envelopes.
 */
class PruneEnvelopes {
public:
	/**
	 * Destination of pruned spectra.
	 */
	SpecSet prunedSpecs;
	
	SpecSet pklSpecs;
	
	/**
	 * Universal TwoValues container.
	 */
	TwoValues<float> twovalues;
	
	/**
	 * Holds charge and score of all peaks in all pruned spectrum.
	 */
	vector<vector<TwoValues<float> > > specChargeScore;
	
	/**
	 * IsoEnvelope instance that extracts and scores envelopes. Contains expected distributions isotopic envelopes.
	 */
	IsoEnvelope env;
	
	/**
	 * Default constructor
	 *
	 *@return
	 */
	PruneEnvelopes();
	~PruneEnvelopes();
	
	/**
	 * Constructor that loads a IsoEnvelope from file. An error message is displayed if the file cannot be loaded.
	 *
	 *@param filename IsoEnvelope file
	 *@return
	 */
	PruneEnvelopes(const char* filename);
	
	/**
	 * Constructor that assigns env an IsoEnvelope instance.
	 *
	 *@param env2 IsoEnvelope instance
	 *@return
	 */
	PruneEnvelopes(IsoEnvelope& env2);
	
	/**
	 * Prunes a pkl spectrum by extracting all isotopic envelopes and converting them to single monoisotopic masses with summed intensities.
	 *
	 *@param spectrum pkl input spectrum
	 *@param putHere output spectrum with all peaks converted to charge 1
	 *@return putChargeScore output container of a charge and score for each monoisotopic peak
	 *@param peakTol peak tolerance of pkl spectrum
	 *@param threshold minimum score of each mass/charge combination. Typically around 0.5.
	 *@param strictMode passed to IsoEnvelope's extract and score functions.
	 *@return
	 */
	void getPrunedSpectrum(Spectrum& spectrum, Spectrum* putHere, vector<TwoValues<float> >* putChargeScore, float peakTol, float threshold, bool strictMode, bool usePPM = false);
	bool uploadPklBin(const char* filename, float peakTol, float threshold, bool strictMode, bool usePPM = false);
	bool uploadMGF(const char* filename, float peakTol, float threshold, bool strictMode, bool usePPM);
	void transferAllSpecs(float peakTol, float threshold, bool strictMode, bool usePPM = false);
	
	float getMonoisotopicmass(float mass, unsigned short charge);
	//bool outputSpectrumMicCompare(char* outfile, char* micfile, char* pklfile, float peakTol);
};
#endif
