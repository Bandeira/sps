#ifndef __SPECTRUM_LIBRARY_H__
#define __SPECTRUM_LIBRARY_H__

/*! @file MsLibrary.h
	\brief Holds the spectrum library class.
*/

#include "../PepNovo/SpectraAggregator.h"
#include "../PepNovo/AnnotationFile.h"
#include "../MsCluster/MsClusterIncludes.h"
#include "../MsCluster/Cluster.h"



class Config; // fwd dclr
struct MsParameterStruct; // fwd dclr
class  AllScoreModels; // fwd dclr

/// the increment in the similarity values (for the pvalue computations)
const float CDF_BIN_SIZE = 0.005;

/// the number of bins for similarity values
const size_t NUM_CDF_BINS = 201;

inline
size_t computeCdfBin(float x)
{
	if (x<=0.0)
		return 0;
	if (x>=1.0)
		return 201;
	return (static_cast<size_t>(x*200+0.499));
}



/*! @struct LibraryEntry
	\brief holds the information necessary for spectral comparisons.
*/
struct LibraryEntry {
	bool operator< (const LibraryEntry& rhs) const
	{
		return (header.getMOverZ() < rhs.header.getMOverZ());
	}

	size_t				 position; /// 
	SingleSpectrumHeader header;
	DistancePeakList	 distancePeaks;
};


/*! @class SpectrumLibrary
	\brief Holds the spectrum library (for rapid peptide identification)
*/
class SpectrumLibrary {
public:

	/*! \brief Loads a set of previously learned cdfs for similarity of random spectra pairs.

		The cdfs should be created per model since they assume certain settings such as
		peak density and fragment/precursor mass tolerance. The actual values that are loaded
		are the cdfs of the probability of observing a random pair of spectra with similarity>X
		(these get converted to cdf values for matches <= X).
	*/
	void loadCdfs(const Config* config);


	/*! \brief Writes a set of previously learned cdfs for similarity of random spectra pairs.

	The actual values that are written are the cdfs of the probability of observing a random 
	pair of spectra with similarity>X.
	*/
	void writeCdfs(const Config* config) const;


	/*! \brief Learns the cdf parameters for having a match with similarity<=s when comparing two random spectra.

		This function uses two lists of annotated spectra (preferably from different speicies).
		Each test spectrum is compared to various numbers of background spectra to determine the probability
		of a random match.

		@param backgroundSpectra The "library" of spectra (should be large).
		@param testSpectraList The spectra that are compared to the library.
		@param config The models Config class.
	*/
	void learnCdfsForPValues(const string& backgroundSpectra, const string& testSpectraList,
					  const Config* config);

	/*! \brief reads spectra from files and writes them as a library (in dat format)

		This function reads a set of spectra and writes them into a dat file as LibraryEntry.

		@param model The scoring models.
		@param params the parameter sctructure; should inclde:
		--list - The "library" of spectra, preferable mgf files with sequence annotations (SEQ=...)
		--output-name - the name of the dat file prefixes
		--out-dir     - the directory where the files should be written
	*/
	size_t createLibrary(AllScoreModels* model, const MsParameterStruct* params);

	/*! \brief reads spectra and reports statistics
	*/
	void libraryStats(const MsParameterStruct* params, const Config* config);

	/*! \brief reads a set of spectra into memory as a library.

	The function also loads the models precomputed cdfs for pvalue computations (MODEL_lib_cdfs.txt)

	The params should include:
	--dat-list - the library dat files
	--min-mz (optional) - the minimal m/z of a library entry
	--max-mz (optional) - the maximal m/z of a library entry
	*/
	void   readLibrary(const MsParameterStruct* params, const Config* config);



	/*! \brief Checks how well the similarity functions can be used to perform spectral library identification.
	
	The function expects params to include two list of mgf files having spectra of the same peptidesL: \c spectraListToLoad (the library of identified peptides),
	and \c spectraListToTest (the list of spectra that we would like to test). If the lists share certain scans (idnetified
	by the "TITLE=" field), the benchmark aborts.

	@param config The score models Config class.
	@param params The structure with the run's parameters.
	*/
	void benchmarkLibrary(const Config* config,
						  const MsParameterStruct* params);



	/*! \breif identify test files using the library.

	For each spectrum file being tested, a file with the suffix "_libres.txt" is created. It lists
	for each scan the best library match found (title, peptide, charge) also states match similarity
	and p-value and quality score (unless --no-sqs is used).

	@param params should include:
	--library-test (spectra files to test),
	--window parameter (m/z size to consider),
	--out-dir where to write results
	*/
	void identifyUsingLibrary(const MsParameterStruct* params, const Config* config, float maxPValue = 1.0) const;



	/*! \breif computes the probability that a random match between n pairs can have a similarity >=X

		Assumes each pair has prob cdf[X] of having a similarity <= X. Pvalue is computed as 1.0 - (cdf[X])^n.
	*/
	float computeMatchPValue(float similarity, int numPeptidesCompared) const;

private:

	vector<LibraryEntry> libraryEntries_;
	vector<double>  cdf_;
};







#endif



