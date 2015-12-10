#ifndef DB_FASTA_H

/** */
#define DB_FASTA_H

#include "aminoacid.h"
#include "spectrum.h"
#include <vector>
#include <cstdio>
#include <iostream>

using namespace std;

/**
 * TODO: add description
 *
 */
class DB_fasta {
protected:

	/**
	 * TODO: add description
	 *
	 */
	static char emptyStr[];

	/**
	 * Protein descriptions stored as \0-terminated strings
	 *
	 */
	vector<char *> desc;

	/**
	 * Protein sequences stored as \0-terminated strings
	 *
	 */
	vector<char *> sequences;

public:

    /**
     * Protein ID descriptors stored as \0-terminated strings
     *
     */
    vector<char *> IDs;

	/**
	 * Contains each protein's cumulative masses as a Spectrum.
	 */
	vector<Spectrum> masses;

	/**
	 * TODO: add description
	 *
	 */
	~DB_fasta();

	/**
	 * TODO: add description
	 *
	 *@param id
	 *@return
	 */
	char *operator[](char *id) const;

	/**
	 * TODO: add description
	 *
	 *@param index
	 *@return
	 */
	char *operator[](int index) const {
		return sequences[index];
	}
	
	void setLine(int index, char* newLine) {
		sequences[index] = newLine;
	}

	/**
	 * TODO: add description
	 *
	 *@param index
	 *@param putHere
	 */
	void getMassesIdx(int index, vector<float> &putHere) {
		getMasses(sequences[index], putHere);
	}

	/**
	 * TODO: add description
	 *
	 *@param index
	 *@return
	 */
	Spectrum &getMassesSpec(int index);

	/**
	 * Replaces all amino acids prevAA with repAA (e.g., all 'I' become 'L')
	 * on all protein sequences.
	 *
	 *@param prevAA Previous amino acid character
	 *@param repAA Replacement amino acid character
	 */
	void replaceAA(char prevAA, char repAA);

	/**
	 * Get protein ID text from protein index.
	 *
	 *@param index Protein index
	 *@return Protein ID string (char *)
	 */
	char *getID(unsigned int index) const {
		if (index<0 or index>=IDs.size()) return (char *)0;
		else return IDs[index];
	}

	/**
	 * Provide access to protein sequences
	 *
	 *@param index - Protein index
	 *@return (char *) to \0-terminated string (protein sequence)
	 */
	char *getSequence(unsigned int index) const {
		if (index<0 or index>=sequences.size()) return (char *)0;
		else return sequences[index];
	}

	/**
	 * TODO: add description
	 *
	 *@param index
	 *@return
	 */
	char *getDesc(unsigned int index) const {
		if (index<0 or index>=desc.size()) return (char *)0;
		else return desc[index];
	}

	/**
	 * TODO: add description
	 *
	 *@return
	 */
	unsigned int size() const {
		return IDs.size();
	}

	/**
	 * TODO: add description
	 *
	 */
	void reset();

	/**
	 * TODO: add description
	 *
	 *@param filename
	 *@return
	 */
	unsigned int Load(const char *filename);

	/**
	 * TODO: add description
	 *
	 *@param out
	 */
	void output(ostream &out);

	/**
	 * Searches all protein sequences for the string peptide.
	 *
	 *@param peptide
	 *@param matches List of (protein index, peptide) for matched strings
	 *@return Number of matches
	 */
	unsigned int find(const char *peptide, list<pair<int,string> > &matches);
};

/**
 * TODO: add description
 *
 *@param filename
 *@param db
 *@param protID1
 *@param protID2
 *@param matchedIndices
 *@return
 */
bool Load_clustalw(const char *filename, DB_fasta &db, int protID1, int protID2,
		vector<TwoValues<int> > &matchedIndices);


/**
 * TODO: add description
 *
 *@param filename
 *@param db
 *@param matchedProts
 *@param matchedAAs
 *@return
 */
bool Load_clustalw_multiple(const char *filename, DB_fasta &db,
		vector<TwoValues<int> > &matchedProts,
		vector<vector<TwoValues<int> > > &matchedAAs);

/**
 * TODO: add description
 *
 */
class DB_index {


	/**
	 * TODO: add description
	 *
	 */
	class Tag {
	public:


		/**
		 * Tag text.
		 */
		char *text;

		/**
		 * List of tag instances as (protein ID, starting amino acid position).
		 */
		list<pair<int, short> > insts;

		/**
		 * TODO: add description
		 *
		 */
		Tag() {
			text = (char *) 0;
		}

		/**
		 * TODO: add description
		 *
		 *@param other
		 */
		Tag(const Tag &other) {
			if (other.text == (char *) 0) {
				text = (char *) 0;
				insts.clear();
				return;
			}
			text = (char *) malloc(strlen(other.text) + 1);
			strcpy(text, other.text);
			insts = other.insts;
		}

		/**
		 * Destructor.
		 */
		~Tag() {
			if (text != (char *) 0)
				free(text);
		}

		/**
		 * TODO: add description
		 *
		 *@param newText
		 *@param textLen
		 */
		void setTag(char *newText, short textLen) {
			if (text != (char *) 0)
				free(text);
			text = (char*) malloc(textLen + 1);
			strcpy(text, newText);
		}

		/**
		 * TODO: add description
		 *
		 *@param protId
		 *@param tagPos
		 */
		void addInst(int protId, short tagPos) {
			pair<int, short> tmp(protId, tagPos);
			insts.push_back(tmp);
		}

		/**
		 * TODO: add description
		 *
		 *@param newInstance
		 */
		void addInst(pair<int, short> &newInstance) {
			insts.push_back(newInstance);
		}

		/**
		 * TODO: add description
		 *
		 *@param other
		 *@return
		 */
		bool operator==(const Tag &other) {
			return strcmp(text, other.text) == 0;
		}

		/**
		 * TODO: add description
		 *
		 *@param other
		 *@return
		 */
		Tag &operator=(const Tag &other) {
			if (text != (char *) 0)
				free(text);
			text = (char *) malloc(strlen(other.text) + 1);
			strcpy(text, other.text);
			insts = other.insts;
			return *this;
		}
	};

	/**
	 * TODO: add description
	 *
	 */
	vector<vector<vector<Tag> > > index;

	/**
	 * Coefficients for hash1.
	 */
	vector<short> coeffs1;

	/**
	 *  Coefficients for hash2.
	 */
	vector<short> coeffs2;


	/**
	 * TODO: add description
	 *
	 */
	short tagLength;

	/**
	 * TODO: add description
	 *
	 */
	short indexDim1;

	/**
	 * TODO: add description
	 *
	 */
	short indexDim2;

	/**
	 * TODO: add description
	 *
	 *@param tagLength
	 */
	void hash_init(short tagLength);

	/**
	 * TODO: add description
	 *
	 *@param hashIdx
	 *@param tag
	 *@return
	 */
	short hash(short hashIdx, char *tag);

	/**
	 * TODO: add description
	 *
	 *@param tag1
	 *@param tag2
	 *@return
	 */
	bool compareTags(char *tag1, char *tag2) {
		for (short tagIndex = 0; tagIndex < tagLength; tagIndex++)
			if (tag1[tagIndex] != tag2[tagIndex])
				return false;
		return true;
	}

	/**
	 * TODO: add description
	 *
	 *@param t
	 */
	void doNothing(char *t) {}

public:

	/**
	 * TODO: add description
	 *
	 *@param db
	 *@param newIndexDim1
	 *@param newIndexDim2
	 *@param newTagLength
	 */
	DB_index(DB_fasta &db, short newIndexDim1, short newIndexDim2,
			short newTagLength);


	/**
	 * TODO: add description
	 *
	 *@param db
	 *@param newIndexDim1
	 *@param newIndexDim2
	 *@param newTagLength
	 */
	void buildIndex(DB_fasta &db, short newIndexDim1, short newIndexDim2,
			short newTagLength);

	/**
	 * TODO: add description
	 *
	 *@param tag
	 *@param proteinID
	 *@param tagPos
	 */
	void insert(char *tag, int proteinID, short tagPos);

	/**
	 * TODO: add description
	 *
	 *@param tag
	 *@param location
	 *@return
	 */
	bool find(char *tag, list<pair<int, short> > **location);

	/**
	 * TODO: add description
	 *
	 *@param tag
	 *@param db
	 *@param peptides List of (protein index, peptide) for matched tags
	 *@param minMatchFlanking
	 *@param flankPref
	 *@param tolPref
	 *@param flankSuff
	 *@param tolSuff
	 *@return
	 */
	bool find(char *tag, DB_fasta &db, list<pair<int, string> > &peptides,
			short minMatchFlanking, float flankPref, float tolPref,
			float flankSuff, float tolSuff);
};


/**
 * TODO: add description
 *
 *@param proteinHits
 */
void MaximumParsimony(vector<list<int> > &proteinHits);

/**
 * MatchSpecToPeptide - scores the direct/reversed match of the spectrum spec
 * against a given peptide sequence; the match scores are defined as the summed
 * intensity of matched peaks in spec.
 *
 *@param spec PRM or MS/MS spectrum
 *@param peptide Peptide string
 *@param peakTol Peak mass tolerance
 *@param ionOffset Ion offset added to theoretical peptide masses to match masses in spec
 *@param filterMatched Select matched peaks in spec
 *@param pepMatches Indices of matched peptide masses (1-to-1 correspondence to spec if filterMatched set to true)
 *@return Summed score of matched peaks (in spec)
 */
float MatchSpecToPeptide(Spectrum &spec, const char *peptide, float peakTol,
		float ionOffset, bool filterMatched, vector<int> *pepMatches = 0);

/**
 * MatchSpecsToPeps - scores the match of each spectrum in specs against all entries
 *  in peptides (peptide substrings are also matched); each spectrum is assigned to a
 *  peptide with maximal match score.
 *
 *@param specs
 *@param peptides
 *@param peakTol
 *@param ionOffset
 *@param ctermMass
 *@param noisePenalty
 *@param scores
 *@param ionOffsets
 *@param pepMatchScores
 *@param pepMatchedSpecs
 *@param matchSuffix
 *@return
 */
bool MatchSpecsToPeps(SpecSet &specs, SpecSet &peptides, float peakTol,
		float ionOffset, float ctermMass, float noisePenalty,
		vector<TwoValues<float> > &scores, vector<float> &ionOffsets,
		vector<TwoValues<float> > &pepMatchScores,
		vector<list<unsigned int> > &pepMatchedSpecs, bool matchSuffix = true);

/**
 * TODO: add description
 *
 *@param specs
 *@param peptides
 *@param peakTol
 *@param ionOffset
 *@param ctermMass
 *@param scores
 *@param ionOffsets
 *@param pepMatchScores
 *@param pepMatchedSpecs
 *@param matchSuffix
 *@
 *@return
 */
bool MatchSpecsToPeps_old(SpecSet &specs, SpecSet &peptides, float peakTol,
		float ionOffset, float ctermMass, vector<float> &scores,
		vector<float> &ionOffsets, vector<TwoValues<float> > &pepMatchScores,
		vector<list<unsigned int> > &pepMatchedSpecs, bool matchSuffix = true);

/**
 * TODO: add description
 *
 */
void unit_MatchSpecsToPeps();

#endif
