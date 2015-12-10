/*
 * ExecReportSPSStats.h
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

#ifndef EXECREPORTSPSSTATS_H_
#define EXECREPORTSPSSTATS_H_

// Module Includes
#include "ExecBase.h"
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"

// External Includes
#include "db_fasta.h"
#include "PeptideSpectrumMatchSet.h"
#include "MappedSPSStatTable.h"
#include "MappedContigStatTable.h"
#include "utils.h"

// System Includes
#include <string>
#include <vector>

namespace specnets
{
  class ExecReportSPSStats: public ExecBase
  {
  public:
    ExecReportSPSStats(void);

    ExecReportSPSStats(const ParameterList & inputParams);

    virtual ~ExecReportSPSStats(void);

    virtual ExecBase * clone(const ParameterList & input_params) const;

    virtual bool invoke(void);

    virtual bool loadInputData(void);

    virtual bool saveOutputData(void);

    virtual bool saveInputData(std::vector<std::string> & filenames);

    virtual bool loadOutputData(void);

    virtual std::vector<ExecBase *> const & split(int numSplit);

    virtual bool merge(void);

    virtual bool validateParams(std::string & error);

    void FilterFastaProteins(set<int>& target_proteins,
                             vector<string>& put_proteins);

    void FilterSpecIds(vector<string>& put_peptides);

    void ChopContigEnds(int endsChop);

  private:
    SpecSet* m_contigs; // input contigs
    SpecSet* m_stars; // star spectra assmebled into contigs
    abinfo_t* m_abinfo; // abinfo detailing what was assembled in each contig
    SpecSet* m_overlaps; // matchma results overlapping contigs with target proteins
    vector<vector<int> >* m_protMatch; // matchma results assigning contigs to target proteins
    DB_fasta* m_fasta; // target proteins
    MS2ScoringModel* m_model; // scoring model for annotating spectra with b/y ions
    PeptideSpectrumMatchSet* m_specIDs; // database search results for assembled spectra
    vector<string>* m_targetProts;
    SpecSet* m_proteinSpectra;
    vector<string>* m_peptides;
    bool ownInput; //! Does this object "own" the input data structures (and hence have to free them)

    MappedSpecnets* m_mappedProj; // SPS project mapped to target proteins
    bool ownOutput;

  };

}

#endif /* EXECREPORTSPSSTATS_H_ */
