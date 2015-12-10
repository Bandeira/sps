#ifndef MS1_H
#define MS1_H

#include "spectrum.h"
#include <set>

/**
 * TODO: add description
 *
 *@param spec1
 *@param spec2
 *@return
 */
float BinnedMatchProduct(Spectrum &spec1, Spectrum &spec2);

/**
 * TODO: add description
 *
 *@param run1
 *@param run2
 *@param peakTol
 *@param matchedScans
 *@param binSpectra
 */
void AlignChromatography(SpecSet &run1, SpecSet &run2, float peakTol, vector<
		TwoValues<unsigned int> > &matchedScans, bool binSpectra = true);

/**
 *  Isotopic envelope classes & functions.
 */
class IsoEnvelope {

	/**
	 * [i][0] - peptide mass.
	 */
	vector<vector<float> > envelopes;

	/**
	 * Number of isotopes considered in each envelope.
	 * [i][1+j:end] - Propensities of each isotopic peak of mass monoisotopic+j.
	 */
	unsigned int envelopeSize;

	unsigned int findEnvelope(float mass);
	void makeStrict(unsigned int idxEnv, vector<float> &out);

public:

	IsoEnvelope() {
		envelopeSize = 0;
	}

	/**
	 * TODO: add description
	 *
	 *@param filename
	 *@return
	 */
	bool LoadModel(const char *filename);

	/**
	 * TODO: add description
	 *
	 *@param monoisotopicMass
	 *@param massEnvelope
	 *@param strictMode
	 *@return
	 */
	float ScoreEnvelope(float monoisotopicMass, vector<float> &massEnvelope,
			bool strictMode = false);

	/**
	 * TODO: add description
	 *
	 *@param monoisotopicMass
	 *@param charge
	 *@param spec
	 *@param peakTol
	 *@param massEnvelope
	 *@param strictMode
	 */
	void ExtractEnvelope(float monoisotopicMass, unsigned short charge,
			Spectrum &spec, float peakTol, vector<float> &massEnvelope,
			bool strictMode = false);
	
	void ExtractEnvelope(float monoisotopicMass, unsigned short charge,
			Spectrum &spec, float peakTol, vector<float> &massEnvelope,
			vector<list<unsigned int> >& indicesUsed, bool strictMode = false);

	/**
	 * TODO: add description
	 *
	 *@param values
	 *@param avoidZeros true if zeros are to be avoided and false otherwise.
	 *@param startIdx
	 */
	void normalize(vector<float> &values, bool avoidZeros = true,
			unsigned int startIdx = 0);
};

#endif
