#ifndef OVERLAPOUT_H
/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 6/29/09
 */
#define OVERLAPOUT_H

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

using namespace std;

struct intCmp {
	bool operator()( int s1, int s2 ) const {
		return s1 < s2;
	}
};

/**
 * Outputs a semicolon-delimeted overlap of small contig spectra, multiple target proteins, and/or larger orbcal spectra.
 * Also computes statistics of aligning all contigs with all orbcals.
 */
class OverlapOut {
public:
	
	/** 
	 * Maps each contig index to a spectral alignment of it and the target peptide
	 */
	SpecSet overlaps;
	
	/**
	 * Small contig spectra
	 */
	SpecSet contigs;
	
	/**
	 * Larger orbcal spectra
	 */
	SpecSet orbcals;
	
	/**
	 * Fasta file of all target protein sequences
	 */
	DB_fasta fasta;
	
	/**
	 * Maps each contig index to its matched target protein and reverse flag
	 */
	vector<vector<int> > prot_match;
	
	/**
	* Maps each orbcal index to its database matched protein, if it exists
	*/
	BufferedLineReader orbcal_pep;
	
	const char* DELIM;

int Totpeaks;
int Bpeaks;
int Ypeaks;
vector<string> contigSeq;
	vector<string> orbcalSeq;
	
	OverlapOut();
	~OverlapOut();
	
	/**
	 * Constructor for an OverlapOut instance that can output a overlap between small contig spectra and a target protein
	 *
	 *@param over_file file to load into overlaps
	 *@param contigs_file file to load into contigs
	 *@param protein_match file to load into prot_match
	 *@param fasta_file file to load into fasta
	 *@return
	 */
	OverlapOut(const char* over_file, const char* contigs_file, const char* protein_match, const char* fasta_file, const char* delim = ",");
	
	/**
	 * Constructor for an OverlapOut instance that can output a overlap between small contig spectra, orbcal spectra, and a target protein
	 *
	 *@param over_file file to load into overlaps
	 *@param contigs_file file to load into contigs
	 *@param protein_match file to load into prot_match
	 *@param fasta_file file to load into fasta
	 *@param orbcal_seq file to load into orbcal_pep
	 *@param orbcal_spec file to load into orbcals
	 *@return
	 */
	OverlapOut(const char* over_file, const char* contigs_file, const char* protein_match, const char* fasta_file, const char* orbcal_seq, const char* orbcal_spec, const char* delim = ",");
	
	OverlapOut &operator=(const OverlapOut &other);
	
	/**
	 * Outputs an overlap between small contig spectra and a specified protein
	 *
	 *@param filename output file
	 *@param proteinIndex index of target protein in fasta file
	 *@param range validate an AA mass in a contig when its mass is within +- range of the corresponding mass in the target spectrum
	 *@return true if overlap was written successfully, false otherwise
	 */
	bool outputOverlap(const char* filename, int proteinIndex, float range, int line_length = 25);
	
	bool outputOverlapComp(const char* filename, const char* componentsFile, int proteinIndex, float range, int line_length = 25);
	
	bool outputOrbcalOverlapComp(const char* filename, const char* componentFile, int proteinIndex, float range, short precCharge, float tolerance, int line_length = 25);
	
	/**
	 * Outputs an overlap between small contig spectra, larger orbcal spectra, and a specified protein
	 *
	 *@param filename output file
	 *@param proteinIndex index of target protein in fasta file
	 *@param precCharge display all orbcals with at least precCharge
	 *@param range validate an AA mass in a contig when its mass is within +- range of the corresponding mass in the target spectrum
	 *@param tolerance an orbcal peak will be displayed if it is within +/- tolerance of an actual peak in the orbcal spectrum
	 *@return true if overlap was written successfully, false otherwise
	 */
	bool outputOrbcalOverlap(const char* filename, int proteinIndex, short precCharge, float range, float tolerance, int line_length = 25);
	
	/**
	 * Outputs an overlap between small contig spectra, larger orbcal spectra, and a specified protein
	 * that only displays contigs/orbcals with a minimum number of matching peaks with other orbcals/contigs
	 *
	 *@param filename output file
	 *@param proteinIndex index of target protein in fasta file
	 *@param precCharge display all orbcals with at least precCharge
	 *@param range validate an AA mass in a contig when its mass is within +- range of the corresponding mass in the target spectrum
	 *@param tolerance an orbcal peak will be displayed if it is within +/- tolerance of an actual peak in the orbcal spectrum
	 *@param match only display contigs/orbcals that have at least match matching peaks with each other
	 *@param checkConsPeaks true if match applies to consecutive matching peaks, false otherwise
	 *@return true if overlap was written successfully, false otherwise
	 */
	bool outputOrbcalContigFiltOverlap(const char* filename, int proteinIndex, short precCharge, float range, float tolerance, int match, bool checkConsPeaks, int line_length = 25);
	
	/**
	 * Gets all contig and/or orbcal ids that match the specified protein and orders them by the first overlapping index in the protein.
	 *
	 *@param protein_index index of target protein in fasta file
	 *@param precCharge order orbcals with at least precCharge
	 *@param range validate an AA mass in a contig when its mass is within +- range of the corresponding mass in the target spectrum
	 *@param tolerance validate an orbcal peak if it is within +/- tolerance of an actual peak in the orbcal spectrum
	 *@param index_order address of resulting contig order. For each val in index_order, val[0] == contig id, val[1] == start overlap index
	 *@param orbcal_order address of resulting orbcal order. For each val in orbcal_order, val[0] == orbcal id, val[1] == start overlap index
	 *@param curPair empty TwoValues
	 *@param orbcal if true, prune orbcals as well as contigs
	 *@return
	 */
	void pruneContigs(int proteinIndex, float range, list<TwoValues<int> >& index_order, TwoValues<int>& curPair);
	
	void pruneOrbcals(int proteinIndex, short precCharge, float tolerance, list<TwoValues<int> >& orbcal_order, TwoValues<int>& curPair);
	
	/**
	 * Outputs the statistics of an alignment of all contigs against all orbcals when compared to actual contig/overlaps. The output format is ";"-delimeted.
	 *
	 *@param outfile file name for statistics
	 *@param alignment alignment results from getOrbcalContigAlignment in batch.cpp
	 *@param precCharge use all orbcals with at least precCharge
	 *@param range validate an AA mass in a contig when its mass is within +- range of the corresponding mass in the target spectrum
	 *@param tolerance validate an orbcal peak if it is within +/- tolerance of an actual peak in the orbcal spectrum
	 *@param startMinNumMatchedPeaks beginning number of minimum matched peaks for stats
	 *@param endMinNumMatchedPeaks ending number of minimum matched peaks for stats
	 *@param stepMinNumMatchedPeaks step size of minimum matched peaks for stats
	 *@param startMinRatio beginning minRatio for stats
	 *@param endMinRatio end minRatio for stats
	 *@param stepMinRatio step size minRatio for stats
	 *@return
	 */
	void outputStats(char* outfile, list<Results_OCC>& alignment, short precCharge, float range, float tolerance, float parentMassTol, int startMinNumMatchedPeaks, int endMinNumMatchedPeaks, int stepMinNumMatchedPeaks, float startMinRatio, float endMinRatio, float stepMinRatio);
	
	
	/**
	 * Outputs the statistics of an alignment of all contigs against all orbcals when compared to actual contig/overlaps. The output format is ";"-delimeted.
	 *
	 *@param outfile file name for statistics
	 *@param inAlignmentFile filename of alignment results
	 *@param precCharge use all orbcals with at least precCharge
	 *@param range validate an AA mass in a contig when its mass is within +- range of the corresponding mass in the target spectrum
	 *@param tolerance validate an orbcal peak if it is within +/- tolerance of an actual peak in the orbcal spectrum
	 *@param startMinNumMatchedPeaks beginning number of minimum matched peaks for stats
	 *@param endMinNumMatchedPeaks ending number of minimum matched peaks for stats
	 *@param stepMinNumMatchedPeaks step size of minimum matched peaks for stats
	 *@param startMinRatio beginning minRatio for stats
	 *@param endMinRatio end minRatio for stats
	 *@param stepMinRatio step size minRatio for stats
	 *@return
	 */
	void outputStats(char* outfile, char* inAlignmentFile, short precCharge, float range, float tolerance, float parentMassTol, int startMinNumMatchedPeaks, int endMinNumMatchedPeaks, int stepMinNumMatchedPeaks, float startMinRatio, float endMinRatio, float stepMinRatio);
	
protected:
	
vector<int> orbprot;
vector<int> contprot;

map<int, int> contigComponents;
bool checkComps;
	/**
	 * Identifies which constructor was used
	 */
	int usesOrbcal;
	
	/**
	 * Identifies which contigs have already been reversed
	 */
	set<int> reversed;
	
	/**
	 * Maps a protein index to list of contig ids that have a peak at index. Populated after calling pruneContigsOrbcals
	 */
	vector<set<int> > proteinToContig;
	
	/**
	 * Maps a protein index to list of orbcal ids that have a peak at index. Populated after calling pruneContigsOrbcals
	 */
	vector<set<int> > proteinToOrbcal;
	
	vector<list<int> > proteinToOrbSt;
	vector<list<int> > proteinToConSt;
	vector<bool> orbcalsFiltered;
	vector<bool> contigsFiltered;
	vector<map<int, float> > orbcalContShifts;
	vector<map<int, float> > contContShifts;
	vector<map<int, float> > orbcalOrbcalShifts;
	
	set<int> mySet;
	list<unsigned int> newL;
	map<unsigned int, list<unsigned int> > newM;
	map<unsigned int, Results_OCC> newMRes;
	
	list<float> myListFloat;
	
	/** Writes a 2 or 3-way overlap of the specified protein.
	 *
	 *@param filiname file to be written
	 *@param proteinIndex index of target protein
	 *@param showOrbcal 1 if orbcal overlap is to be shown, 0 otherwise
	 *@param precCharge display all orbcals with at least precCharge
	 *@param range a contig AA will be matched to a peptide AA if its mass is within +- range of the peptide AA
	 *@param tolerance a peak will be displayed if it is within +- tolerance of an actual peak in orbcal
	 *@param filterOrbcal only display contigs/orbcals that have at least match matching peaks with each other
	 *@param checkConsPeaks true if match applies to consecutive matching peaks, false otherwise
	 *@return true if overlap was written successfully, false otherwise
	 */
	bool outputOverlap(const char* filename, int proteinIndex, short showOrbcal, float range, float tolerance, short precCharge, int filterOrbcal, bool checkConsPeaks, int line_length = 25);
	
	/**
	 * Writes a contig overlapped with the protein to file and returns the index
	 * of the last AA written before end_pep_index is reached (-1 if the end of the contig
	 * was written).
	 *
	 *@param output pointer to file output stream
	 *@param peptide char pointer to target peptide
	 *@param PepSpec pointer to Spectrum of target peptide
	 *@param contig_ident index of contig to be written
	 *@param start_pep_index index of first character of the target peptide on the current line
	 *@param prev_c index of last AA matched from the target contig ( 0 if just beginning contig overlap)
	 *@param end_pep_index index of last character of the target peptide on the current line + 1
	 *@param range a contig AA will be matched to a peptide AA if its mass is within +- range of the peptide AA
	 *@return	last written index of the contig (-1 if the last contig character was written)
	 */
	int outputContig(FILE* output, char* peptide, Spectrum* PepSpec, int contig_ident, int start_pep_index, int prev_c, int end_pep_index, float range);
	
	/**
	 * Writes an orbcal overlapped with the protein to file and returns the index
	 * of the last AA written before end_pep_index is reached (-1 if the end of the contig
	 * was written).
	 *
	 *@param output pointer to file output stream
	 *@param peptide char pointer to target peptide
	 *@param orbcal_ident index of orbcal to be written
	 *@param start_pep_index index of first character of the target peptide on the current line
	 *@param start_overlap index of frist characer of peptide to be aligned with orbcal
	 *@param prev_c index of last AA matched from the target orbcal ( 0 if just beginning contig overlap)
	 *@param end_pep_index index of last character of the target peptide on the current line + 1
	 *@param tolerance an orbcal break will be displayed if its mass is within +- tolerance of an actual orbcal break
	 *@return	last written index of the orbcal (-1 if the last orbcal character was written)
	 */
	int outputOrbcalContig(FILE* output, char* peptide, int orbcal_ident, int start_pep_index, int start_overlap, int prev_c, int end_pep_index, float tolerance);
	
	void getCorrectShifts(int index, int start, int pepIndex, vector<list<int> >& start_order_other, vector<set<int> >& proteinToOrbcal, float tolrange, float pmTol, bool contig);
	
	TwoValues<float> getContigOverlapScore(int idx1, int idx2, float shiftFor, float shiftRev, float pmTol, bool debug = false);
	
	void getCorrectConnShifts(int index, int start, int pepIndex, vector<list<int> >& start_order_orbcal, vector<set<int> >& proteinToOrbcal, float tolerance, float pmTol);
	
	void getContigDistance(int index, vector<list<int> >& start_order_other, float pmTol);
	
	/**
	 * Gets the mass of an AA
	 *@param aa pointer to character that corresponds to an AA mass
	 *@return	mass of AA, 0.0 if the AA is not recognized
	 */
	float getMass(char *aa);
	
	/**
	 * Outputs the proper break between 2 orbcal AAs
	 *
	 *@param output pointer to file output stream
	 *@param orbcal pointer to spectrum of target orbcal
	 *@param bpeak_ possible prefix mass
	 *@param ypeak_ possible suffix mass
	 *@param tolerance a peak will be displayed if it is within +- tolerance of an actual peak in orbcal
	 *@return
	 */
	void outputOrbcalBreak(FILE* output, Spectrum* orbcal, float bpeak_, float ypeak_, float tolerance);
	
	/**
	 * Determines if an orbcal has at least filterOrbcal matching peaks with 2 or more contigs
	 * 
	 *@param index orbcal index
	 *@param start starting protein index of overlap
	 *@param tolerance
	 *@param filterOrbcal minimum number of matching peaks with at least 2 other contigs
	 *@return	true if orbcal is filtered, false if not
	 */
	bool isFilteredOrbcal(int index, int start, float tolerance, int filterOrbcal, float peakTol, bool countConsec);
	
	/**
	 * Determines if a contig has at least filterOrbcal matching peaks with at least 1 other orbcal
	 * 
	 *@param index contig index
	 *@param filterOrbcal minimum number of matching peaks with at least 1 other orbcal
	 *@return	true if contig is filtered, false if not
	 */
	bool isFilteredContig(int index, int filterOrbcal, float peakTol, bool countConsec);
	
	/**
	 * Maps all matching peak indicies in proteinToOrbcal to specified orbcal
	 * 
	 *@param index orbcal index
	 *@param start starting protein index of overlap
	 *@param tolerance
	 *@return
	 */
	void rememberOrbcal(int index, int start, float tolerance);
	
	/**
	 * Maps all matching peak indicies in contigToOrbcal to specified contig
	 * 
	 *@param index contig index
	 *@param range
	 *@param PepSpec pointer to spectrum of target protein
	 *@return
	 */
	void rememberContig(int index, float range, Spectrum* PepSpec);
	
	unsigned int numMatchingPeaks(Spectrum& contig, Spectrum& orbcal, float contigOffset, float peakTol, bool countConsec);
	
	/**
	 * Determines if an orbcal spectrum contains one of the two peaks within tolerance
	 * 
	 *@param orbcal
	 *@param bpeak_
	 *@param ypeak_
	 *@param tolerance
	 *@return	true if orbcal contains one of the two peaks within tolerance, false otherwise
	 */
	bool hasPeak(Spectrum* orbcal, float bpeak_, float ypeak_, float tolerance);
};
#endif
