#ifndef BATCH_H
#define BATCH_H

#include "spectrum.h"
#include "clusters.h"
#include "utils.h"
#include "PepnovoTags.h"
#include <vector>
#include <list>
#include <fstream>
#include <map>
#include <set>
#include <time.h>

using namespace std;

/**
 * Stores the results of aligning two spectra.
 */
class PairAlignRes {
public:

	/**
	 * Total summed score of the alignment.
	 */
	float score;

	/**
	 * Indices of the matched PRMs in spectrum 1.
	 */
	vector<int> idxMatched1;

	/**
	 * Indices of the matched PRMs in spectrum 2.
	 */
	vector<int> idxMatched2;
};


/**
 * TODO: add description
 *
 */
class PairAlign {
public:

	/**
	 * Index of the first spectrum to align.
	 */
	int spec1;

	/**
	 * Index of the second spectrum to align.
	 */
	int spec2;

	/**
	 * Shifts to compute.
	 */
	vector<float> shifts;

	/**
	 * Storage for the computed results.
	 */
	vector<PairAlignRes> results;

	PairAlign() {
	}
};


/**
 * TODO: add description
 *
 */
class Results_SP {
public:

	/**
	 * Index of the first spectrum to align.
	 */
	int spec1;

	/**
	 * Index of the second spectrum to align.
	 */
	int spec2;

	/**
	 * Match scores in spectrum1.
	 */
	float score1;

	/**
	 * Match scores in spectrum2.
	 */
	float score2;

	/**
	 * TODO: add description
	 *
	 *@param other
	 *@return
	 */
	virtual Results_SP &operator=(const Results_SP &other) {
		spec1 = other.spec1;
		spec2 = other.spec2;
		score1 = other.score1;
		score2 = other.score2;
		return *this;
	}

	/**
	 * TODO: add description
	 *
	 *@param out
	 *@param separator
	 */
	void output(ostream &out, const char separator) {
		out << spec1 + 1 << separator << spec2 + 1 << separator << score1
				<< separator << score2 << endl;
	}

	/**
	 * TODO: add description
	 *
	 *@param v
	 */
	void serialize(vector<float> &v) {
		v.resize(4);
		v[0] = spec1;
		v[1] = spec2;
		v[2] = score1;
		v[3] = score2;
	}

	/**
	 * TODO: add description
	 *
	 *@return
	 */
	unsigned int loadSz() {
		return 4;
	}

	/**
	 * TODO: add description
	 *
	 *@param v
	 */
	void load(vector<float> &v) {
		if (v.size() < 4)
			return;
		spec1 = (unsigned int) v[0];
		spec2 = (unsigned int) v[1];
		score1 = v[2];
		score2 = v[3];
	}
};

class Results_OCC {
public:

	/**
	 * Index of the first spectrum to align.
	 */
	int spec1;

	/**
	 * Index of the second spectrum to align.
	 */
	int spec2;
	
	float shift1;
	float shift2;
	
	/**
	 * Match scores in spectrum1.
	 */
	float score1;

	/**
	 * Match scores in spectrum2.
	 */
	float score2;
	
	int matchedPeaks;
	/**
	 * TODO: add description
	 *
	 *@param other
	 *@return
	 */
	virtual Results_OCC &operator=(const Results_OCC &other) {
		spec1 = other.spec1;
		spec2 = other.spec2;
		shift1 = other.shift1;
		shift2 = other.shift2;
		score1 = other.score1;
		score2 = other.score2;
		matchedPeaks = other.matchedPeaks;
		return *this;
	}
	
	/**
	 * TODO: add description
	 *
	 *@param out
	 *@param separator
	 */
	void output(FILE* out, const char* separator = " ") {
		
		fprintf(out, "%d%s%d%s%.5f%s%.5f%s%.5f%s%.5f%s%d\n", spec1, separator, spec2, separator, shift1, separator, shift2, separator, score1, separator, score2, separator, matchedPeaks);
	}

	/**
	 * TODO: add description
	 *
	 *@param out
	 *@param separator
	 */
	void output(ostream &out, const char separator = ' ') {
		out << spec1 << separator << spec2 << separator << shift1 << separator << shift2 << separator << score1
				<< separator << score2 << separator << matchedPeaks << endl;
	}

	/**
	 * TODO: add description
	 *
	 *@param v
	 */
	void serialize(vector<float> &v) {
		v.resize(7);
		v[0] = spec1;
		v[1] = spec2;
		v[2] = shift1;
		v[3] = shift2;
		v[4] = score1;
		v[5] = score2;
		v[6] = matchedPeaks;
	}

	/**
	 * TODO: add description
	 *
	 *@return
	 */
	unsigned int loadSz() {
		return 7;
	}

	/**
	 * TODO: add description
	 *
	 *@param v
	 */
	void load(vector<float> &v) {
		if (v.size() < 7)
			return;
		spec1 = (unsigned int) v[0];
		spec2 = (unsigned int) v[1];
		shift1 = v[2];
		shift2 = v[3];
		score1 = v[4];
		score2 = v[5];
		matchedPeaks = (unsigned int) v[6];
	}
};

/**
 * TODO: add description
 *
 */
class Results_PA;

/**
 * TODO: add description
 *
 */
class Results_ASP: public Results_SP {
public:
    Results_ASP() : shift1() {}

	/**
	 * Value of the symmetric shift in the Almost Same Peptide
	 * problem (first shift is always zero).
	 */
	float shift1;


	/**
	 * TODO: add description
	 *
	 *@param other
	 *@return
	 */
	virtual Results_ASP &operator=(const Results_ASP &other) {
		spec1 = other.spec1;
		spec2 = other.spec2;
		score1 = other.score1;
		score2 = other.score2;
		shift1 = other.shift1;
		return *this;
	}

	/**
	 * TODO: add description
	 *
	 *@param out
	 *@param separator
	 */
	void output(ostream &out, const char separator) {
		out << spec1 + 1 << separator << spec2 + 1 << separator << shift1
				<< separator << score1 << separator << score2 << endl;
	}

	/**
	 * TODO: add description
	 *
	 *@param o
	 *@return
	 */
	Results_ASP &operator=(Results_PA &o);

	/**
	 * TODO: add description
	 *
	 *@param v
	 */
	void serialize(vector<float> &v) {
		v.resize(5);
		v[0] = spec1;
		v[1] = spec2;
		v[2] = shift1;
		v[3] = score1;
		v[4] = score2;
	}

	/**
	 * TODO: add description
	 *
	 *@return
	 */
	unsigned int loadSz() {
		return 5;
	}

	/**
	 * TODO: add description
	 *
	 *@param v
	 */
	void load(vector<float> &v) {
		if (v.size() < 5)
			return;
		spec1 = (unsigned int) v[0];
		spec2 = (unsigned int) v[1];
		shift1 = v[2];
		score1 = v[3];
		score2 = v[4];
	}
};

/**
 * TODO: add description
 *
 */
class Results_PA: public Results_ASP {
public:
    Results_PA() : shift2() {}

	/**
	 * Value of the second shift in the Pairwise Alignment problem
	 * (first shift is stored in shift1).
	 */
	float shift2;

	/**
	 * TODO: add description
	 *
	 *@param other
	 *@return
	 */
	virtual Results_PA &operator=(const Results_PA &other) {
		spec1 = other.spec1;
		spec2 = other.spec2;
		score1 = other.score1;
		score2 = other.score2;
		shift1 = other.shift1;
		shift2 = other.shift2;
		return *this;
	}

	/**
	 * TODO: add description
	 *
	 *@param other
	 *@return
	 */
	virtual Results_PA &operator=(const Results_ASP &other) {
		spec1 = other.spec1;
		spec2 = other.spec2;
		score1 = other.score1;
		score2 = other.score2;
		shift1 = other.shift1;
		shift2 = 0;
		return *this;
	}

	/**
	 * TODO: add description
	 *
	 *@param out
	 *@param separator
	 */
	void output(ostream &out, const char separator) {
		out << spec1 + 1 << separator << spec2 + 1 << separator << shift1
				<< separator << shift2 << separator << score1 << separator
				<< score2 << endl;
	}

	/**
	 * TODO: add description
	 *
	 *@param v
	 */
	void serialize(vector<float> &v) {
		v.resize(6);
		v[0] = spec1;
		v[1] = spec2;
		v[2] = shift1;
		v[3] = shift2;
		v[4] = score1;
		v[5] = score2;
	}

	/**
	 * TODO: add description
	 *
	 *@return
	 */
	unsigned int loadSz() {
		return 6;
	}

	/**
	 * TODO: add description
	 *
	 *@param v
	 */
	void load(vector<float> &v) {
		if (v.size() < 6)
			return;
		spec1 = (unsigned int) v[0];
		spec2 = (unsigned int) v[1];
		shift1 = v[2];
		shift2 = v[3];
		score1 = v[4];
		score2 = v[5];
	}
};

/**
 * TODO: add description
 *
 *@param o
 *@return
 */
inline Results_ASP &Results_ASP::operator=(Results_PA &o) {
	spec1 = o.spec1;
	spec2 = o.spec2;
	shift1 = o.shift1;
	score1 = o.score1;
	score2 = o.score2;
	return *this;
}

// Results_CS - alignment results for connector spectra (spectra connecting pairs of contigs)
//   spec1, spec2 - indices of paired contigs
//   specC - index of connector spectrum
//   shift1, shift2 - shift from spec1 to spec2, shift from spec1 to specC (respectively)


/**
 * TODO: add description
 *
 */
class Results_CS: public Results_PA {
public:

	/**
	 *
	 * Value of the second shift in the Pairwise Alignment problem
	 * (first shift is stored in shift1).
	 */
	int specC;

	/**
	 * Boolean value indicating whether spec2 was matched as-is (==false) or
	 * reversed (==true).
	 */
	bool spec2rev;

	/**
	 *
	 *@param out
	 *@param separator
	 */
	void output(ostream &out, const char separator) {
		out << spec1 + 1 << separator << spec2 + 1 << separator << specC + 1
				<< separator << shift1 << separator << shift2 << separator
				<< score1 << separator << score2 << endl;
	}
};

/**
 * TODO: add description
 *
 *@param filename
 *@param aligns
 *@return
 */
bool Load_pabatch(const char *filename, vector<PairAlign> &aligns);

/**
 * TODO: add description
 *
 *@param filename
 *@param aligns
 *@return
 */
bool Save_results(const char *filename, vector<PairAlign> &aligns);

bool combineResultsOCC(const char* filename, const char* outfile, SpecSet& spectra, vector<int>& charges, vector<int>& matchedPeaks, vector<float>& minRatio);

bool Save_resultsOCC(const char *filename, list<Results_OCC> &aligns);

bool Load_resultsOCC(const char *filename, list<Results_OCC> &aligns);

/**
 * TODO: add description
 *
 *@param filename
 *@param aligns
 *@param sep
 *@return
 */
bool Save_resultsCS(const char *filename, list<Results_CS> &aligns, const char sep);
/*
 bool Load_resultsASP(char *filename, vector<Results_ASP> &results);
 bool Load_resultsASPbin(char *filename, vector<Results_ASP> &results);
 bool Load_resultsASPbin_multiple(char *filename, vector<Results_ASP> &results);
 bool Load_resultsPA(char *filename, vector<Results_PA> &results);
 bool Save_resultsASPbin(char *filename, list<Results_ASP> &results);
 bool Save_resultsASPbin(char *filename, vector<Results_ASP> &results);
 bool Save_resultsPAbin(char *filename, list<Results_PA> &results);
 bool Load_results_bin(char *filename, vector<Results_ASP> &results);
 bool Load_results_bin(char *filename, vector<Results_PA> &results);
 bool Load_resultsASPbin(char *filename, vector<Results_ASP> &results);
 bool Save_resultsASPbin(char *filename, vector<Results_ASP> &results);
 bool Load_resultsPAbin(char *filename, vector<Results_PA> &results);
 */

/**
 * TODO: add description
 *
 *@param filename
 *@param results
 *@return
 */
bool Load_resultsPA(const char *filename, vector<Results_PA> &results);

/**
 * TODO: add description
 *
 *@param filename
 *@param results
 *@return
 */
bool Load_results(const char *filename, vector<Results_ASP> &results);

/**
 * TODO: add description
 *
 *@param filename
 *@param results
 *@return
 */
bool Load_results(const char *filename, vector<Results_PA> &results);

/**
 * TODO: add description
 *
 *@param filename
 *@param results
 *@return
 */
template<class T> bool Load_results_bin(const char *filename, vector<T> &results);

/**
 * TODO: add description
 *
 *@param filename
 *@param results
 *@return
 */
template<class T> bool Load_results_bin_multiple(const char *filename,
		vector<T> &results);

/**
 * TODO: add description
 *
 *@param filename
 *@param count
 *@param it
 *@return
 */
template<class TIter> bool Save_results_bin(const char *filename, unsigned int count,
		TIter it) {
	FILE *fp; vector<float> dataV; float data[6];
	unsigned int i,j;

	fp = fopen(filename,"w");
	if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return false; }
	fwrite(&count,sizeof(unsigned int),1,fp);  // Number of entries in the file
	for(i=0;i<count;i++) {
		it->serialize(dataV);   it++;
		for(j=0;j<2;j++) data[j]=dataV[j]+1;    // Spectrum indices are 1-based
		for(j=2;j<dataV.size();j++) data[j]=dataV[j];
		fwrite(data,sizeof(float),dataV.size(),fp);
	}
	fclose(fp);
	return true;
}


/**
 * Computes Almost Same Peptide matches between spectra.
 *
 *@param specSet
 *@param baseSpectraIdx
 *@param aaDiff
 *@param minShift
 *@param maxShift
 *@param pmTol
 *@param peakTol
 *@param minMatchRatio
 *@param results
 *@param ratios
 *@param means
 *@param tagsMatchFlank
 *@param tagsMatchCount
 *@param resolution
 *@param symmetryOffset
 */
void getPairAlignsASP(SpecSet &specSet, vector<int> &baseSpectraIdx,
		short aaDiff, float minShift, float maxShift, float pmTol,
		float peakTol, float minMatchRatio, list<Results_ASP> &results, list<
				TwoValues<float> > &ratios, vector<TwoValues<float> > &means,
		vector<float> &varTerms, vector<vector<TTag> > &tags,
		unsigned int tagsMatchFlank, unsigned int tagsMatchCount,
		float resolution = 0.1, float symmetryOffset = 0);


/**
 *
 * Computes Pairwise Alignments between spectra (TODO: NEEDS DEBUGGING).
 *
 *@param specSet
 *@param pmTol
 *@param peakTol
 *@param results
 */
void getPairAlignsPA(SpecSet &specSet, float pmTol, float peakTol, vector<
		Results_PA> &results);

/**
 * Paraneters for contig/mate-pair alignments
 *@connectors mate-pair spectra
 *@contigs contig spectra
 *@mpsToAlign all matepairs to align (only checked by getMatepairContigAlignment_thread)
 *@outputFile temporary output file of alignments
 *@minRatio minimum overlapping intensity of all alignments to report
 *@minNumMatchedPeaks minimum number of matching peaks of all alignments to report
 *@minAbsShift minimum absolute shift to consider
 *@allowedJumps allowed amino acid jumps (not utilized)
 *@precCharge minimum precursor charge of all mate-pairs to align
 *@checkConsPeaks whether or not to report minNumMatchedPeaks as the minimum number of consecutive matching peaks on the contig
 *@results output alignment results
 */
struct MatepairContigAlignParams {
	SpecSet* connectors;
	SpecSet* contigs;
	list<int> mpsToAlign;
	string outputFile;
	
	float minRatio;
	int minNumMatchedPeaks;
	float minAbsShift;
	AAJumps* allowedJumps;
	short precCharge;
	bool checkConsPeaks;
	
	list<Results_OCC> results;
};

//Copies all fields in MatepairContigAlignParams except results
void copy(MatepairContigAlignParams& from, MatepairContigAlignParams& to);

/**
 * Computes the best alignment between a matepair and a contig
 *@idx1 index of matepair
 *@idx2 index of contig
 *@params alignment parameters (has pointer to contig  and matepair spectra)
 *@result best alignment result (if score is above threshold)
 *@return true if the best alignment is above thresholds in params, false if not
 */
bool alignMatepairContig(int mpIdx, int contigIdx, MatepairContigAlignParams* params, Results_OCC& result);

/**
 * Computes all pair-wise alignments between matepairs and contigs with parallel threads
 *@inparams alignment parameters
 *@num_threads number of threads to spawn with pthread to do the heavy computation
 *@return
 */
void getMatepairContigAlignment(MatepairContigAlignParams* inparams, unsigned int num_threads);

/**
 * Computes best alignments for each matepair specified in vparams->pairsToCompute against all contigs. Can be called by multiple pthreads. 
 *@vparams alignment parameters
 *@return
 */
void* getMatepairContigAlignment_thread(void* vparams);

/**
 * Reads the temporary file  written by getMatepairContigAlignment_thread and deletes the file
 *@buf_filename temporary buffer of Results_OCC. Methid will attempt to delete it after reading
 *@results list of Results_OCC to dump alignment results in file
 *@return
 */
bool transferTempBuffer( const char* buf_filename, list<Results_OCC>& results);

/**
 * Parameters for contig/contig alignments
 *@contigs contig spectra
 *@pairsToCompute all contig/contig pairs to align (only checked by getContigContigAlignment_thread)
 *@minScore minimum overlapping score of all alignments to compute (overlapping peaks * overlapping intensity ratio)
 *@minNumMatchedPeaks minimum number of overlapping peaks
 *@minAbsShift minimum shift to consider
 *@allowedJumps allowed amino acid jumps (not utilized)
 *@results output alignment results
 */
struct ContigContigAlignParams {
	SpecSet* contigs;
	list<TwoValues<int> > pairsToCompute;
	
	float minScore;
	int minNumMatchedPeaks;
	float minAbsShift;
	AAJumps* allowedJumps;
	
	list<Results_OCC> results;
};

//Copies all fields in ContigContigAlignParams except results
void copy(ContigContigAlignParams& from, ContigContigAlignParams& to);

/**
 * Computes the best alignment of a pair of contigs
 *@idx1 index of 1st contig in pair
 *@idx2 index of 2nd contig in pair
 *@params alignment parameters (has pointer to contig spectra)
 *@result best alignment result (if score is above threshold)
 *@return true if the best alignment is above thresholds in params, false if not
 */
bool alignContigs(int idx1, int idx2, ContigContigAlignParams* params, Results_OCC& result);

/**
 * Computes all contig pair-wise alignments with parallel threads
 *@inparams alignment parameters
 *@num_threads number of threads to spawn with pthread to do the heavy computation
 *@return
 */
void getContigContigAlignment(ContigContigAlignParams* inparams, unsigned int num_threads = 1);

/**
 * Computes best alignments for each pair specified in vparams->pairsToCompute. Can be called by multiple pthreads. 
 *@vparams alignment parameters
 *@return
 */
void* getContigContigAlignment_thread(void* vparams);

/**
 * Computes Pairwise Alignments between spectra (TODO: NEEDS DEBUGGING).
 *
 *@param specSet
 *@param startIdx
 *@param endIdx
 *@param pmTol
 *@param minRatio
 *@param minPeakAreaOvlp
 *@param minNumMatchedPeaks
 *@param allowedJumps
 *@param minAbsShift
 *@param results
 *@param ratios
 *@param numMatchedPeaks
 *@param means
 *@param varTerms
 *@param alignStats
 *@param specStats
 */
void getPairAlignsPA2(SpecSet &specSet, unsigned int startIdx,
		unsigned int endIdx, float peakTol, float pmTol, float minRatio,
		float minPeakAreaOvlp, short minNumMatchedPeaks, AAJumps &allowedJumps,
		float minAbsShift, list<Results_PA> &results,
		list<TwoValues<float> > &ratios,
		list<TwoValues<int> > &numMatchedPeaks,
		vector<TwoValues<float> > &means, vector<float> &varTerms, list<vector<
				float> > &alignStats, vector<vector<float> > &specStats);


/**
 * Extends consensus multiple alignments by aligning to single PRM spectra.
 *
 *@param clusters
 *@param specSet
 *@param startIdx
 *@param endIdx
 *@param aaDiff
 *@param pmTol
 *@param peakTol
 *@param minMatchRatio
 *@param results
 *@param ratios
 *@param meansClst
 *@param meansSpecs
 *
 */
void getPairAlignsPAext(Clusters &clusters, SpecSet &specSet, int startIdx,
		int endIdx, short aaDiff, float pmTol, float peakTol,
		float minMatchRatio, list<Results_PA> &results,
		list<TwoValues<float> > &ratios, vector<TwoValues<float> > &meansClst,
		vector<TwoValues<float> > &meansSpecs);

/**
 * TODO: add description
 *
 *@param prmSpec
 *@param cSpec
 *@param cSpecRev
 *@param aaDiff
 *@param pmTol
 *@param peakTol
 *@param result
 */
void getPairAlignsPAext_aux(Spectrum &prmSpec, Spectrum &cSpec,
		Spectrum &cSpecRev, short aaDiff, float pmTol, float peakTol,
		Results_PA &result);

/**
 * Connect consensus multiple alignments by finding single PRM spectra to
 * connect pairs of contigs (multiple alignments).
 *
 *@param contigs
 *@param specSet
 *@param peakTol
 *@param minRatio
 *@param minPeakDist
 *@param mergedPeakCount
 *@param results
 *@param ratios
 */
void findConnectorSpectra(Clusters &contigs, SpecSet &specSet, float peakTol,
		float minRatio, float minPeakDist, int mergedPeakCount,
		list<Results_CS> &results, list<vector<float> > &ratios);

#endif
