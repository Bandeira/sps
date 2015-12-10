// Header Includes
#include "ExecSpecNetsPropagation.h"

// Module Includes
#include "ExecTagSearch.h"
#include "Logger.h"
//#include "FileUtils.h"

// SpecNets Includes
#include "SpectralPairs.h"
#include "tags.h"
#include "tuple.h"
#include "utils.h"
#include "db_fasta.h"

using namespace std;
using namespace specnets;

namespace specnets
{

	ExecSpecNetsPropagation::ExecSpecNetsPropagation(void) :
        m_msms_spectra(0x0), m_prm_spectra(0x0), m_star_spectra(0x0),
        m_svm(0x0), m_svm_model(0x0), m_pairs(0x0), m_db(0x0), ownInput(true),
        m_propagation_pairs(0x0), m_propagation_info(0x0),
        m_psms_spectra(0x0), ownOutput(true)
  {
    m_name = "ExecSpecNetsPropagation";
    m_type = "ExecSpecNetsPropagation";
  }

  // -------------------------------------------------------------------------

	ExecSpecNetsPropagation::ExecSpecNetsPropagation(const ParameterList & inputParams) :
        ExecBase(inputParams), m_msms_spectra(0x0), m_prm_spectra(0x0),
        m_star_spectra(0x0), m_svm(0x0), m_svm_model(0x0), m_pairs(0x0), m_db(0x0), ownInput(true),
        m_propagation_pairs(0x0), m_propagation_info(0x0),
        m_psms_spectra(0x0), ownOutput(true)
  {
    m_name = "ExecSpecNetsPropagation";
    m_type = "ExecSpecNetsPropagation";
  }

  // -------------------------------------------------------------------------

	ExecSpecNetsPropagation::ExecSpecNetsPropagation(const ParameterList & inputParams,
                                        SpecSet * msms_spectra,
                                        SpecSet * prm_spectra,
                                        SpecSet * star_spectra,
                                        ExecSvmStatistics * svm,
                                        MS2ScoringModel * svm_model,
                                        SpectrumPairSet * pairs,
                                        DB_fasta * db,
                                        PeptideSpectrumMatchSet * psms,
                                        SpectrumPairSet * propagation_pairs,
                                        vector<vector<int> > * propagation_info,
                                        SpecSet * psms_spectra) :
    ExecBase(inputParams), m_msms_spectra(msms_spectra), m_prm_spectra(prm_spectra),
        m_star_spectra(star_spectra), m_svm(svm), m_svm_model(svm_model),
        m_pairs(pairs), m_db(db), m_psms(psms), ownInput(false),
        m_propagation_pairs(propagation_pairs), m_propagation_info(propagation_info),
        m_psms_spectra(psms_spectra), ownOutput(false)
  {
    m_name = "ExecSpecNetsPropagation";
    m_type = "ExecSpecNetsPropagation";
  }

  // -------------------------------------------------------------------------

	ExecSpecNetsPropagation::~ExecSpecNetsPropagation(void)
  {
    if (ownInput) {
      if (m_msms_spectra)
        delete m_msms_spectra;
      if (m_prm_spectra)
        delete m_prm_spectra;
      if (m_star_spectra)
        delete m_star_spectra;
      if (m_pairs)
        delete m_pairs;
      if (m_db)
        delete m_db;
    }
    if (ownOutput) {
      if (m_propagation_pairs)
        delete m_propagation_pairs;
      if (m_propagation_info)
        delete m_propagation_info;
      if(m_psms_spectra)
        delete m_psms_spectra;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecSpecNetsPropagation::clone(const ParameterList & inputParams) const
  {
    return new ExecSpecNetsPropagation(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecSpecNetsPropagation::invoke(void)
  {
//    return invoke_deprecated();
    if(!m_msms_spectra or !m_prm_spectra or !m_star_spectra or !m_pairs) {
    	ERROR_MSG("ERROR: Input data not available [in ExecSpecNetsPropagation::invoke()]");
    	return false;
    }

    if(!m_psms_spectra or !m_svm or !m_svm_model) {
    	ERROR_MSG("ERROR: Output container objects not available [in ExecSpecNetsPropagation::invoke()]");
    	return false;
    }

    unsigned int numSpecs = m_star_spectra->size();
    SpecSet *specSetRev = new SpecSet;  // Reversed versions of every star spectrum
    DEBUG_TRACE;
    m_psms_spectra->resize(numSpecs);
    DEBUG_TRACE;
    for(unsigned int idxSpec = 0; idxSpec < numSpecs ; idxSpec++)
      (*m_psms_spectra)[idxSpec] = (*m_star_spectra)[idxSpec];
    DEBUG_TRACE;

    /*
     * Load InsPecT annotations if available
     */
    if(m_psms) {
      	PeptideSpectrumMatchSet * psms = new PeptideSpectrumMatchSet(*m_psms); // InsPecT PSMs used to initialize the annotations
//      psms->loadInspectResultsFile(m_params.getValue("INSPECT_PSMS").c_str(), true);
//      DEBUG_TRACE;
        psms->addSpectra(m_psms_spectra);
        DEBUG_TRACE;

/*
        // Use InsPecT IDs to get initial mp/midx
        ParameterList tagsearchParams(m_params);
        tagsearchParams.setValue("MAX_NUM_TAGS","-1");
        ExecTagSearch inspectPSMprocessor( tagsearchParams,
                                           m_psms_spectra,
                                           m_db,
                                           (vector<unsigned int> *) 0);
        DEBUG_TRACE;
        if( not inspectPSMprocessor.invoke() ) return false;
        DEBUG_TRACE;
*/

cerr<<"ExecSpecNetsPropagation::invoke() -- >> Outputting m_psms_mp/midx:\n";
for(unsigned int i=0; i<m_psms_spectra->size(); i++) {
	cerr<<" -- -- spectrum "<<i<<", scan = "<<(*m_psms_spectra)[i].scan<<", seq = ";
	if(!(*m_psms_spectra)[i].psmList.empty()) {
		psmPtr psm = (*m_psms_spectra)[i].psmList.front();
		cerr<<psm->m_annotation<<", midx:\n";
		for(unsigned int j=0; j < psm->m_matchedPeaks.size(); j++)
			cerr<<psm->m_matchedPeaks[j][0]<<" "<<psm->m_matchedPeaks[j][1]<<"\n";
	} else cerr<<"\n";
}
cerr.flush();


        delete psms;
    }

    /*
     * Initialize data structures for projection of spectrum identifications
     */

    specSetRev->resize(numSpecs);
    vector<float> svmScores(numSpecs);
    map<pair<float, int>, pair<int, bool> > rankedPSMs;  // Ranked set of PSMs
                                      // Used to keep track of best projection and select the next PSM to project from
                                      //   Key pair: SVM score, spectrum index
                                      //   Value pair: index of projection source spectrum (-1 if none), true = reversed spectrum
    map<pair<float, int>, pair<int, bool> >::iterator rankedPSMs_iter;
    map<pair<float, int>, pair<int, bool> >::reverse_iterator rankedPSMs_riter;
    pair<pair<float,int>, pair<int,bool> > rankedPSMelmt;

    float minExpInt = (float) m_params.getValueDouble("MIN_PERC_EXPINT");
    float minTP = (float) m_params.getValueDouble("MIN_PERC_TP");
    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");
    unsigned int minMatchedPeaks = m_params.exists("MIN_MATCHED_PEAKS")?(int) m_params.getValueInt("MIN_MATCHED_PEAKS"):0;
    bool projectionKeepAllPRMs = m_params.getValue("SPECNETS_PROJ_TYPE","all") == "all";  // Possible values are "all" / "matched"

    AAJumps aminoacids(1); // Used to convert projected spectra to annotation strings
    unsigned int numMatchedPairs = 0;       // Number of pairs used for propagation
    vector<bool> processed(numSpecs);

    if(m_propagation_pairs) m_propagation_pairs->resize(m_pairs->size());
    if(m_propagation_info) {
    	m_propagation_info->resize(numSpecs);
    	for (unsigned int i=0; i<numSpecs; i++) {
    		(*m_propagation_info)[i].resize(5);
    		for (unsigned int j=0; j<5; j++)
    			(*m_propagation_info)[i][j]=0;
    	}
    }

    DEBUG_VAR(numSpecs);
    int annotatedCount=0;
    for(unsigned int i=0; i<numSpecs; i++) {
      DEBUG_VAR(i);

      (*m_star_spectra)[i].reverse(0, &(*specSetRev)[i]);
      processed[i] = false;
      svmScores[i] = -1000000;

      bool annotated = not (*m_psms_spectra)[i].psmList.empty();
      DEBUG_VAR(annotated);
      if( annotated ) {
        if(m_propagation_info) {
          (*m_propagation_info)[i][0] = 0;   // Level of propagation/propagation (none yet)
          (*m_propagation_info)[i][1] = i;   // Spectrum where annotation comes from (self)
          (*m_propagation_info)[i][2] = -1;  // N/A
          (*m_propagation_info)[i][3] = -1;  // N/A
          (*m_propagation_info)[i][4] = -1;  // N/A
        }

        DEBUG_VAR((*m_psms_spectra)[i].psmList.front()->m_score);
// TODO: Once it's available
        //   svmScores[i] = (*m_psms_spectra)[i].psmList.front()->m_score;
// Meanwhile just add up matched-PRM scores in star spectrum - doesn't work because PRMs are not filtered yet
//        svmScores[i] = 0;
//        for(unsigned int p=0; p<(*m_psms_spectra)[i].size(); p++)
//          svmScores[i] += (*m_psms_spectra)[i][p][1];
        svmScores[i] = MatchSpecToPeptide((*m_psms_spectra)[i],
        		                          (*m_psms_spectra)[i].psmList.front()->m_annotation.c_str(),
        		                          peakTol, 0, false);

        rankedPSMelmt.first.first = svmScores[i];
        rankedPSMelmt.first.second = i;       // spectrum index
        rankedPSMelmt.second.first = -1;      // source of annotation is self
        rankedPSMelmt.second.second = false;  // not reversed
        rankedPSMs.insert(rankedPSMelmt);

        annotatedCount++;
      }
      DEBUG_VAR(svmScores[i]);
    }
    DEBUG_MSG("Number of annotated spectra: "<<annotatedCount);

    vector<list<int> > neighs(numSpecs);  // Neighbors per spectrum
    for(unsigned int i=0; i<m_pairs->size(); i++) {
    	neighs[(*m_pairs)[i].spec1].push_back((*m_pairs)[i].spec2);
    	neighs[(*m_pairs)[i].spec2].push_back((*m_pairs)[i].spec1);
    }

    /*
     * Propagate spectrum identifications
     */

    Spectrum curProj, bestProj;
    float curScore, curScoreRev;
    list<int> neighToProcess;     // Index of the non-annotated node during propagation (data structure is a side effect of using ProjectSpectrum: here always contains a single element)
    list<int> curDeltas, curBestDeltas;
    vector<TwoValues<int> > matchesSpec;  // Indices of projection-matched peaks in annotated / non-annotated spectra

    while(not rankedPSMs.empty()) {
cerr<<"rankedPSMs.size() = "<<rankedPSMs.size()<<endl;
      // Process the highest scoring PSM and remove it from rankedPSMs
      rankedPSMs_riter = rankedPSMs.rbegin();
      rankedPSMelmt = (*rankedPSMs_riter);
      rankedPSMs.erase(rankedPSMelmt.first);
      unsigned int idxCurSpec     = rankedPSMelmt.first.second;

cerr<<"Processing rankedPSMelmt = ("<<rankedPSMelmt.first.first<<","<<rankedPSMelmt.first.second<<") / ("<<rankedPSMelmt.second.first<<","<<rankedPSMelmt.second.second<<")\n";

      if(processed[idxCurSpec] and rankedPSMelmt.second.first>=0) {
        cerr<<"ERROR: attempting to project onto a spectrum (index "<<idxCurSpec<<") that is already annotated!\n";
        exit(-1);
      }

//TODO: Keep PRMs with score zero at true prefix masses with no peaks in PRM spectra (helps quality of annotations over poor star spectra)

      //
      // ---------------------------------------------------------------------------------------------
      //   Stage 1: Finalize projection (if there was any projection better than initial annotation)
      // ---------------------------------------------------------------------------------------------
      //
      if(rankedPSMelmt.second.first >= 0) {
        SpecSet *curSpecSet;
        unsigned int idxPropgtnSpec = rankedPSMelmt.second.first;
        bool reversedProjection     = rankedPSMelmt.second.second;
        if( reversedProjection ) curSpecSet = specSetRev; else curSpecSet = m_star_spectra;

        // Reproduce best projection
        neighToProcess.clear();   neighToProcess.push_front(idxCurSpec);
        curDeltas.clear();
        curScore=0;
        curBestDeltas.clear();
        curProj.resize(0);
        matchesSpec.resize(0);

cerr<<"ExecSpecNetsPropagation.invoke(), before ProjectSpectrum:\n";
cerr<<"ExecSpecNetsPropagation.invoke(), projecting from "<<idxPropgtnSpec<<" in m_psms_spectra:\n";
(*m_psms_spectra)[idxPropgtnSpec].output(cerr);
cerr<<"ExecSpecNetsPropagation.invoke(), projecting to "<<idxCurSpec<<" in curSpecSet:\n";
(*curSpecSet)[idxCurSpec].output(cerr);

        ProjectSpectrum((*curSpecSet),
        		        (*m_psms_spectra)[idxPropgtnSpec],
        		        neighToProcess,
        		        curDeltas,
        		        curScore,
        		        curBestDeltas,
        		        peakTol,
        		        &curProj,
        		        &matchesSpec,
        		        minMatchedPeaks);

cerr<<"ExecSpecNetsPropagation.invoke().matchesSpec:\n";
for(unsigned int p=0; p<matchesSpec.size(); p++)
    cerr<<matchesSpec[p][0]<<"\t "<<matchesSpec[p][1]<<endl;

        (*curSpecSet)[idxCurSpec].copyNP(curProj);
        if(projectionKeepAllPRMs)
        { // Propagate the annotation, keep all annotation peaks even if not from this spectrum (accumulated evidence)
          for(unsigned int p=0; p<curProj.size(); p++)
            curProj[p][1]=0;
          for(unsigned int p=0; p<matchesSpec.size(); p++)
              curProj[matchesSpec[p][0]] = (*curSpecSet)[idxCurSpec][matchesSpec[p][1]];
        }
        else
        { // Propagate the annotation, keep only annotated peaks
          curProj.resize(matchesSpec.size());
          for(unsigned int p=0; p<matchesSpec.size(); p++)
              curProj[p] = (*curSpecSet)[idxCurSpec][matchesSpec[p][1]];
        }

cerr << "-- old peaks for spectrum " << idxCurSpec <<" with parent mass = "<<(*m_psms_spectra)[idxCurSpec].parentMass<<" and scan # = "<<(*m_psms_spectra)[idxCurSpec].scan<<":\n";
(*m_psms_spectra)[idxCurSpec].output(cerr);
cerr.flush();
//        (*m_spectra)[idxCurSpec] = curProj;
        (*m_psms_spectra)[idxCurSpec] = curProj;
cerr << "-- propagated peaks for spectrum " << idxCurSpec <<" with parent mass = "<<(*m_psms_spectra)[idxCurSpec].parentMass<<" and scan # = "<<(*m_psms_spectra)[idxCurSpec].scan<< ":\n";
(*m_psms_spectra)[idxCurSpec].output(cerr);
cerr.flush();

        // Update output data structures

		DEBUG_TRACE;

        // Statistics for best propagation
        float maxSpecScore = 0;
        for(unsigned int p=0; p<(*curSpecSet)[idxCurSpec].size(); p++) maxSpecScore += (*curSpecSet)[idxCurSpec][p][1];
        if(m_propagation_info) {
          (*m_propagation_info)[idxCurSpec][0] = (*m_propagation_info)[idxPropgtnSpec][0] + 1;
          (*m_propagation_info)[idxCurSpec][1] = idxPropgtnSpec;
          (*m_propagation_info)[idxCurSpec][2] = curDeltas.empty() ? 0 : curDeltas.front();
          (*m_propagation_info)[idxCurSpec][3] = maxSpecScore>0 ? (int)round(10000*( svmScores[idxCurSpec] / maxSpecScore )) : 0;
          (*m_propagation_info)[idxCurSpec][4] = (int)round(10000*(((float)matchesSpec.size())/((float)(*m_psms_spectra)[idxPropgtnSpec].size())));
        }

        DEBUG_TRACE;

        // Add propagated annotation as spectrum PSM
        (*m_psms_spectra)[idxCurSpec].psmList.clear();
        psmPtr propagatedPSM(new PeptideSpectrumMatch),
        	   propgtnPSM = (*m_psms_spectra)[idxPropgtnSpec].psmList.front();
        aminoacids.getPeptideFromSpectrum((*m_psms_spectra)[idxCurSpec], propagatedPSM->m_annotation, peakTol);
        propagatedPSM->m_spectrum = &((*m_psms_spectra)[idxCurSpec]);
        DEBUG_VAR(propagatedPSM->m_spectrum);
        propagatedPSM->m_scanNum = (*m_psms_spectra)[idxCurSpec].scan;
        DEBUG_VAR(propagatedPSM->m_scanNum);
        propagatedPSM->m_dbIndex = propgtnPSM->m_dbIndex;
        DEBUG_VAR(propagatedPSM->m_dbIndex);
        DEBUG_TRACE;
// TODO: Handle sequence extensions/deletions appropriately
        if(projectionKeepAllPRMs) {
        	// Keeping all PRMs so protein matches are exactly the same
        	propagatedPSM->m_matchedPeaks.resize(propgtnPSM->m_matchedPeaks.size());
        	for(unsigned int p=0; p < propgtnPSM->m_matchedPeaks.size(); p++)
        		propagatedPSM->m_matchedPeaks[p] = propgtnPSM->m_matchedPeaks[p];
        } else {
        	propagatedPSM->m_matchedPeaks.resize(matchesSpec.size());
        	for(unsigned int p=0; p < matchesSpec.size(); p++)
        		propagatedPSM->m_matchedPeaks[p].set(matchesSpec[p][0], propgtnPSM->m_matchedPeaks[matchesSpec[p][0]][1]);
        }
//(*m_psms_midx)[idxCurSpec].output(cerr);
//cerr.flush();
        DEBUG_TRACE;
        (*m_psms_spectra)[idxCurSpec].psmList.push_back(propagatedPSM);

        DEBUG_TRACE;

        // Update set of pairs used for propagations
        if(m_propagation_pairs) {
          (*m_propagation_pairs)[numMatchedPairs].spec1 = idxCurSpec;
          (*m_propagation_pairs)[numMatchedPairs].spec2 = idxPropgtnSpec;
          (*m_propagation_pairs)[numMatchedPairs].shift1 = curDeltas.empty() ? 0 : curDeltas.front();
          (*m_propagation_pairs)[numMatchedPairs].score1 = svmScores[idxCurSpec];
          (*m_propagation_pairs)[numMatchedPairs].score2 = -1;   // missing value
        }

        DEBUG_TRACE;

        // Set PSM midx/mp
/* --->>> FIX to midx/mp in PSM
        if(m_psms_mp and m_psms_midx) {
          (*m_psms_mp)[idxCurSpec].resize(3);
          (*m_psms_mp)[idxCurSpec][0] = (*m_psms_mp)[idxPropgtnSpec][0];      // Same protein match as source of projection
          (*m_psms_mp)[idxCurSpec][1] = max(0,(*m_psms_mp)[idxPropgtnSpec][1]) + 1;  // Number of mods set to number of projection iterations
          (*m_psms_mp)[idxCurSpec][2] = (*m_psms_mp)[idxPropgtnSpec][2] ^ (reversedProjection ? 1 : 0);  // Reversed with relation to protein match (yes=1/no=0)
//cerr << "-- midx for spectrum " << idxCurSpec << " (mp = [" << (*m_psms_mp)[idxCurSpec][0] << "," << (*m_psms_mp)[idxCurSpec][1] << "," << (*m_psms_mp)[idxCurSpec][2] << "]:\n";

          DEBUG_TRACE;

          if(projectionKeepAllPRMs) {
            (*m_psms_midx)[idxCurSpec] = (*m_psms_midx)[idxPropgtnSpec];  // Keeping all PRMs so protein matches are exactly the same
          } else {
            (*m_psms_midx)[idxCurSpec].resize( (*m_propagation_peaks)[numMatchedPairs].size() );
            for(unsigned int p = 0; p < (*m_propagation_peaks)[numMatchedPairs].size(); p++) {
              unsigned int idx_source_peak      = (*m_propagation_peaks)[numMatchedPairs][p][1];
              unsigned int idx_protein_peak     = (*m_psms_midx)[idxPropgtnSpec][idx_source_peak][1];  // Propagated peaks are always matched to the protein (by definition)
              (*m_psms_midx)[idxCurSpec][p].set( p , idx_protein_peak );
            }
          }
//(*m_psms_midx)[idxCurSpec].output(cerr);
//cerr.flush();
          DEBUG_TRACE;

        }
*/

        numMatchedPairs++;

      }  // Finalize projection
      else  // Enforce original PSM PRMs
      {
    	  psmPtr psm = (*m_psms_spectra)[idxCurSpec].psmList.front();
    	  aminoacids.getPRMMasses( psm->m_annotation,
    			                   (*m_psms_spectra)[idxCurSpec]);

    	  // TEMPORARY DEBUG code to test propagation while ExecTagSearch is crashing
    	  int startProtPos;
    	  psm->m_matchedPeaks.resize((*m_psms_spectra)[idxCurSpec].size());
    	  if(idxCurSpec==1) startProtPos = 6;
    	  if(idxCurSpec==2) startProtPos = 4;
    	  if(idxCurSpec==3) startProtPos = 6;
    	  if(idxCurSpec==4) startProtPos = 4;
    	  for(unsigned int k=0; k<psm->m_matchedPeaks.size(); k++)
    		  psm->m_matchedPeaks[k].set(k,k+startProtPos);
      }

      //
      // ---------------------------------------------------------------------------------
      //   Stage 2: attempt to propagate final annotation to all non-annotated neighbors
      // ---------------------------------------------------------------------------------
      //
      PeptideSpectrumMatch svmPSM;
      for(list<int>::iterator curNeigh=neighs[idxCurSpec].begin(); curNeigh!=neighs[idxCurSpec].end(); curNeigh++) {
cerr<<"Projecting to neighbor "<<*curNeigh<<" of idxCurSpec = "<<idxCurSpec<<", processed["<<*curNeigh<<"] = "<<processed[*curNeigh]<<endl;
        if(processed[*curNeigh]) {
cerr<<"Skipping *curNeigh = "<<*curNeigh<<" from idxCurSpec = "<<idxCurSpec<<" (processed[*curNeigh] = "<<processed[*curNeigh]<<")\n";
          continue;  // Do not re-process spectra with finalized annotations
        }

        // Project using non-reversed spectrum
        neighToProcess.clear();   neighToProcess.push_front(*curNeigh);
        curDeltas.clear();        curBestDeltas.clear();
        curProj.resize(0);        matchesSpec.resize(0);
        curScore=0;               curScoreRev = 0;
cerr<<"Projecting *curNeigh = "<<*curNeigh<<" using annotated base idxCurSpec = "<<idxCurSpec<<"\n";
        ProjectSpectrum(*m_star_spectra, (*m_psms_spectra)[idxCurSpec], neighToProcess, curDeltas, curScore, curBestDeltas, peakTol, &curProj, &matchesSpec, minMatchedPeaks);

        // Compute SVM match score

        svmPSM.m_charge = (*m_msms_spectra)[*curNeigh].parentCharge;
        svmPSM.m_scanNum = (*curNeigh) + 1;
        DEBUG_VAR(svmPSM.m_scanNum);
        aminoacids.getPeptideFromSpectrum(curProj, svmPSM.m_annotation, peakTol);
        DEBUG_VAR(svmPSM.m_annotation);

        if(*curNeigh==0) svmPSM.m_annotation="PEPSIDE";

#if 0  // Need to fix SVM scores/models - current scores are nonsensical
        DEBUG_VAR(svmPSM.m_annotation);
        if (!m_svm->getSvmScore( &(*m_msms_spectra)[*curNeigh],
                                 &(*m_prm_spectra)[*curNeigh],
                                 &(*m_star_spectra)[*curNeigh],
                                 svmPSM,
                                 m_svm_model))
        {
          DEBUG_MSG("Unable to generate svm score in ExecSpecNetsPropagation::invoke()");
          DEBUG_VAR(idxCurSpec);
          DEBUG_VAR(*curNeigh);
          return false;
        }
        DEBUG_VAR(svmPSM.m_scanNum);
        DEBUG_VAR(svmPSM.m_annotation);
        DEBUG_VAR(svmPSM.m_score);
        curScore = svmPSM.m_score;
#else
        // Compute score as sum of matched PRM scores
        curScore=0; for(unsigned int p=0; p<matchesSpec.size(); p++) {
          curScore+=(*m_star_spectra)[*curNeigh][matchesSpec[p][1]][1];
        DEBUG_VAR(curScore);
//cerr<<"--- curScore += ("<<(*specSetRev)[*curNeigh][matchesSpec[p][1]][0]<<","<<(*specSetRev)[*curNeigh][matchesSpec[p][1]][1]<<") = "<<curScore<<endl;
        }

#endif

cerr<<"--- got scores = "<<curScore<<" / "<<svmScores[*curNeigh]<<" :"; curProj.output(cerr);
//cerr<<"--- --- ---\n";

        // Project using reversed spectrum
        curDeltas.clear();        curBestDeltas.clear();
        curProj.resize(0);        matchesSpec.resize(0);
//cerr<<"Projecting *curNeigh = "<<*curNeigh<<" using annotated base idxCurSpec = "<<idxCurSpec<<" (specSetRev)\n";
        ProjectSpectrum(*specSetRev, (*m_psms_spectra)[idxCurSpec], neighToProcess, curDeltas, curScoreRev, curBestDeltas, peakTol, &curProj, &matchesSpec, minMatchedPeaks);

#if 0  // Need to fix SVM scores/models - current scores are nonsensical
        aminoacids.getPeptideFromSpectrum(curProj, svmPSM.m_annotation, peakTol);
        DEBUG_VAR(svmPSM.m_annotation);
        if (!m_svm->getSvmScore( &(*m_msms_spectra)[*curNeigh],
                                 &(*m_prm_spectra)[*curNeigh],
                                 &(*m_star_spectra)[*curNeigh],
                                 svmPSM,
                                 m_svm_model))
        {
          DEBUG_MSG("Unable to generate svm score in ExecSpecNetsPropagation::invoke()");
          DEBUG_VAR(idxCurSpec);
          DEBUG_VAR(*curNeigh);
          return false;
        }
        DEBUG_VAR(svmPSM.m_scanNum);
        DEBUG_VAR(svmPSM.m_annotation);
        DEBUG_VAR(svmPSM.m_score);
        curScoreRev = svmPSM.m_score;

#else
        curScoreRev=0; for(unsigned int p=0; p<matchesSpec.size(); p++) {
          curScoreRev+=(*specSetRev)[*curNeigh][matchesSpec[p][1]][1];
//cerr<<"--- curScoreRev += ("<<(*specSetRev)[*curNeigh][matchesSpec[p][1]][0]<<","<<(*specSetRev)[*curNeigh][matchesSpec[p][1]][1]<<") = "<<curScoreRev<<endl;
        }
#endif

cerr<<"--- got scoresRev = "<<curScoreRev<<" / "<<svmScores[*curNeigh]<<" :"; curProj.output(cerr);
//cerr<<"--- --- ---\n";

        // Update rankedPSMs if this projection improves the neighbor's annotation score
        if( curScore > svmScores[*curNeigh] or curScoreRev > svmScores[*curNeigh] ) {
          // Remove current entry from rankedPSMs (if any)
//cerr<<"rankedPSMs.size() = "<<rankedPSMs.size()<<" before erasing ("<<svmScores[*curNeigh]<<","<<*curNeigh<<")\n";
          rankedPSMs.erase( make_pair<float,int> ( svmScores[*curNeigh] , *curNeigh ) );
//cerr<<"rankedPSMs.size() = "<<rankedPSMs.size()<<" after erasing ("<<svmScores[*curNeigh]<<","<<*curNeigh<<")\n";

          // Update rankedPSMs with new projected-annotation source for *curNeigh
          rankedPSMelmt.first.first   = max( curScore, curScoreRev );
          rankedPSMelmt.first.second  = *curNeigh;
          rankedPSMelmt.second.first  = idxCurSpec;
          rankedPSMelmt.second.second = curScoreRev > curScore;
          rankedPSMs.insert(rankedPSMelmt);

          svmScores[*curNeigh] = rankedPSMelmt.first.first;
//cerr<<"Added rankedPSMelmt = ("<<rankedPSMelmt.first.first<<","<<rankedPSMelmt.first.second<<") / ("<<rankedPSMelmt.second.first<<","<<rankedPSMelmt.second.second<<")\n";
        }
      }  // Propagation to non-annotated neighbors

      processed[idxCurSpec] = true;
//cerr<<"Done processing idxCurSpec = "<<idxCurSpec<<", processed["<<idxCurSpec<<"] = "<<processed[idxCurSpec]<<endl;

    } // while(not rankedPSMs.empty())


//    for(unsigned int i=0; i<numSpecs; i++) {
//    	if(m_propagation_info and (*m_propagation_info)[i][0]>0) (*m_propagation_info)[i][1]++; // Convert indices to 1-based instead of 0-based
//    	if(m_psms_spectra and (*m_psms_spectra)[i].size()==0)
//    		(*m_psms_spectra)[i]=(*m_star_spectra)[i];  // Recover original spectra for unsuccessful propagations
//    }

    if(m_propagation_pairs) m_propagation_pairs->resize(numMatchedPairs);

    delete specSetRev;
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecNetsPropagation::invoke_deprecated(void)
  {
	  return false;
#if 0
	  SpecSet *propagationSpecs = new SpecSet;

    if(!m_psms_spectra or !m_pairs or !propagationSpecs)
      return false;

    float minExpInt = (float) m_params.getValueDouble("MIN_PERC_EXPINT");
    float minTP = (float) m_params.getValueDouble("MIN_PERC_TP");
    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");
    unsigned int minMatchedPeaks = m_params.exists("MIN_MATCHED_PEAKS")?(int) m_params.getValueInt("MIN_MATCHED_PEAKS"):0;
    unsigned int numSpecs = m_psms_spectra->size();
    propagationSpecs->resize(numSpecs);
    SpecSet *specSetRev = new SpecSet;  // Reversed versions of every proximately-annotated spectrum
    specSetRev->resize(numSpecs);

    AAJumps aminoacids(1); // Used to convert projected spectra to annotation strings

    if(m_propagation_peaks) m_propagation_peaks->resize( min( (unsigned int)m_psms_spectra->size() , (unsigned int)m_pairs->size() ) );
    if(m_propagation_pairs) m_propagation_pairs->resize(m_pairs->size());
    if(m_propagation_info) {
      m_propagation_info->resize(numSpecs);
      for (unsigned int i=0; i<numSpecs; i++) {
        (*m_propagation_info)[i].resize(5);
        for (unsigned int j=0; j<5; j++)
          (*m_propagation_info)[i][j]=0;
      }
    }

    unsigned int numMatchedPairs = 0;       // Number of pairs used for propagation
    vector<bool> annotated(numSpecs);
    vector<vector<short> > annotatedFromFile;
    int annotatedCount=0;
/*
    if (m_params.exists("INPUT_ANNOTATED")) {
      Load_binArray(m_params.getValue("INPUT_ANNOTATED"), annotatedFromFile);
      if(annotatedFromFile.size()!=numSpecs or annotatedFromFile[0].size()!=2) { cerr<<"ERROR reading "<<m_params.getValue("INPUT_ANNOTATED")<<"!\n"; return -1; }
      for(unsigned int i=0; i<numSpecs; i++) {
        (*propagationSpecs)[i]=(*m_psms_spectra)[i];
        if(annotatedFromFile[i].size()!=2) { cerr<<"ERROR reading "<<m_params.getValue("INPUT_ANNOTATED")<<" at entry number "<<i<<"!\n"; return -1; }
        annotated[i]=(annotatedFromFile[i][0]>=minExpInt-0.0001 and annotatedFromFile[i][1]>=minTP-0.0001);
        if(annotated[i]) {
          if(m_propagation_info) {
            (*m_propagation_info)[i][3]=annotatedFromFile[i][0];
            (*m_propagation_info)[i][4]=annotatedFromFile[i][1];
          }
          annotatedCount++;
        }
      }
      cout<<"Number of annotated spectra: "<<annotatedCount<<endl;
    }
*/

    DEBUG_VAR(numSpecs);
    for(unsigned int i=0; i<numSpecs; i++) {
      DEBUG_VAR(i);
      (*propagationSpecs)[i]=(*m_psms_spectra)[i];

      annotated[i] = not (*m_psms_spectra)[i].psmList.empty();
      DEBUG_VAR(annotated[i]);
      if( annotated[i] ) {
        if(m_propagation_info) {
          (*m_propagation_info)[i][0] = 0;   // Level of propagation/propagation (none yet)
          (*m_propagation_info)[i][1] = i;   // Spectrum where annotation comes from (self)
          (*m_propagation_info)[i][2] = -1;  // N/A
          (*m_propagation_info)[i][3] = -1;  // N/A
          (*m_propagation_info)[i][4] = -1;  // N/A
        }
        annotatedCount++;
      }
    }
    cout<<"Number of annotated spectra: "<<annotatedCount<<endl;

    vector<list<int> > neighs(numSpecs);  // Neighbors per spectrum
    for(unsigned int i=0; i<m_pairs->size(); i++) {
      neighs[(*m_pairs)[i].spec1].push_back((*m_pairs)[i].spec2);
      neighs[(*m_pairs)[i].spec2].push_back((*m_pairs)[i].spec1);
    }

    vector<short> scheduled(numSpecs);
    for(unsigned int i=0; i<numSpecs; i++)
      scheduled[i]=0;  // Indicates the iteration on which a spectrum is to be processed
    vector<list<int> > annNeighs(numSpecs);  // Annotated neighbors to consider when annotating the current spectrum
    list<int> toProcess; toProcess.clear();  // Spectra eligible for annotation on the next iteration
    int annIter=1,              // Number of hops between current iteration and initially annotated spectra
    lastInIteration=0,      // Last spectrum index in toProcess with the current annIter
    lastInNextIteration=-1; // Last spectrum index in toProcess with the next annIter
    for(unsigned int specIdx=0; specIdx<numSpecs; specIdx++)
      if(annotated[specIdx])
        for(list<int>::iterator curNeigh=neighs[specIdx].begin(); curNeigh!=neighs[specIdx].end(); curNeigh++)
          if(not annotated[*curNeigh]) {
            if(scheduled[*curNeigh]==0) {
              scheduled[*curNeigh] = 1;
              annNeighs[*curNeigh].clear();
              toProcess.push_back(*curNeigh);
              lastInIteration = *curNeigh;
            }
            annNeighs[*curNeigh].push_back(specIdx);
          }

    /*
     * Project spectrum identifications
     */
    for(list<int>::iterator curSpec=toProcess.begin(); curSpec!=toProcess.end(); curSpec++) {
      Spectrum curProj, bestProj;
      float bestScore=0, curScore=0;
      list<int> neighToProcess;     // Index of the node being annotated (propagated onto); data structure is a side effect of using ProjectSpectrum: this always contains the single element *curSpec
      list<int> curDeltas, bestDeltas, curBestDeltas;
      vector<TwoValues<int> > matchesSpec;
      list<int>::iterator curNeigh;
      int bestNeigh, bestMatchedPeakCount=0;
      bool bestRev;  // Indicates whether the best propagation was to a reversed version of the spectrum
      float maxSpecScore=0; for(unsigned int i=0; i<(*m_psms_spectra)[*curSpec].size(); i++) maxSpecScore+=(*m_psms_spectra)[*curSpec][i][1];

      bestNeigh = annNeighs[*curSpec].empty()?-1:annNeighs[*curSpec].front();  // initialize with first neighbor

      // Find the annotated neighbor with the best propagation onto the current spectrum
      (*m_psms_spectra)[*curSpec].reverse(0, &(*specSetRev)[*curSpec]);
      curDeltas.clear();   bestDeltas.clear();   curBestDeltas.clear();
      neighToProcess.clear();   neighToProcess.push_front(*curSpec);
      if(m_propagation_peaks) (*m_propagation_peaks)[numMatchedPairs].resize(0);

      /* ProjectSpectrum is used to project the non-annotated spectrum (*curSpec) onto each
       * annotated neighbor (curNeigh). At the end of the loop bestProj contains the modified
       * spectrum of the annotated neighbor best-matched to *curSpec (masses modified by delta
       * and peak scores increased by matched peaks in *curSpec)
       */
      for(curNeigh=annNeighs[*curSpec].begin(); curNeigh!=annNeighs[*curSpec].end(); curNeigh++) {
        curDeltas.clear();   curScore=0;   matchesSpec.resize(0);
cerr<<"Projecting *curNeigh = "<<*curNeigh<<", *curSpec = "<<*curSpec<<"\n";
        ProjectSpectrum(*m_psms_spectra, (*m_psms_spectra)[*curNeigh], neighToProcess, curDeltas, curScore, curBestDeltas, peakTol, &curProj, &matchesSpec, minMatchedPeaks);
        curScore=0; for(unsigned int p=0; p<matchesSpec.size(); p++) curScore+=(*m_psms_spectra)[*curSpec][matchesSpec[p][1]][1];
cerr<<"--- got scores = "<<curScore<<" / "<<bestScore<<" :"; curProj.output(cerr);
cerr<<"--- --- ---\n";
        if(matchesSpec.size()>0 and curScore>bestScore) {
          bestScore = curScore;
          bestDeltas.clear();
          bestDeltas.splice(bestDeltas.end(),curBestDeltas);
          bestNeigh = *curNeigh;
          bestProj = curProj;
          bestRev  = false;
          bestMatchedPeakCount = matchesSpec.size();
          (*propagationSpecs)[*curSpec] = (*m_psms_spectra)[*curSpec];
          (*propagationSpecs)[*curSpec].resize(matchesSpec.size());
          if(m_propagation_peaks) (*m_propagation_peaks)[numMatchedPairs].resize(matchesSpec.size()); // Update sets of peaks matched by propagations
          for(unsigned int p=0; p<matchesSpec.size(); p++) {
            if(propagationSpecs) (*propagationSpecs)[*curSpec][p] = (*m_psms_spectra)[*curSpec][matchesSpec[p][1]];
            if(m_propagation_peaks) (*m_propagation_peaks)[numMatchedPairs][p].set(p,matchesSpec[p][0]);
          }
        }

        curDeltas.clear();   curScore=0;   matchesSpec.resize(0);
cerr<<"Projecting *curNeigh = "<<*curNeigh<<", *curSpec = "<<*curSpec<<" (specSetRev)\n";
        ProjectSpectrum((*specSetRev), (*m_psms_spectra)[*curNeigh], neighToProcess, curDeltas, curScore, curBestDeltas, peakTol, &curProj, &matchesSpec, minMatchedPeaks);
        curScore=0; for(unsigned int p=0; p<matchesSpec.size(); p++) curScore+=(*specSetRev)[*curSpec][matchesSpec[p][1]][1];
cerr<<"--- got scores = "<<curScore<<" / "<<bestScore<<" :"; curProj.output(cerr);
cerr<<"--- --- ---\n";
        if(matchesSpec.size()>0 and curScore>bestScore) {
          bestScore = curScore;
          bestDeltas.clear();
          bestDeltas.splice(bestDeltas.end(),curBestDeltas);
          bestNeigh = *curNeigh;
          bestProj = curProj;
          bestRev  = true;
          bestMatchedPeakCount = matchesSpec.size();
          (*propagationSpecs)[*curSpec] = (*specSetRev)[*curSpec];
          (*propagationSpecs)[*curSpec].resize(matchesSpec.size());
          if(m_propagation_peaks)
            (*m_propagation_peaks)[numMatchedPairs].resize(matchesSpec.size()); // Update sets of peaks matched by propagations
          for(unsigned int p=0; p<matchesSpec.size(); p++) {
            if(propagationSpecs) (*propagationSpecs)[*curSpec][p] = (*specSetRev)[*curSpec][matchesSpec[p][1]];
            if(m_propagation_peaks) (*m_propagation_peaks)[numMatchedPairs][p].set(p,matchesSpec[p][0]);
          }
        }
      }

      scheduled[*curSpec] = false;
      if(bestProj.size()>0) {
        annotated[*curSpec] = true;
        //      (*m_psms_spectra)[*curSpec] = bestProj;   // Propagate the annotation, keep all annotation peaks even if not from this spectrum (accumulated evidence)
cerr << "-- old peaks for spectrum " << *curSpec << ":\n";
(*m_psms_spectra)[*curSpec].output(cerr);
cerr.flush();
        (*m_psms_spectra)[*curSpec] = (*propagationSpecs)[*curSpec];   // Propagate the annotation, keep only annotated peaks
cerr << "-- propagated peaks for spectrum " << *curSpec << ":\n";
(*m_psms_spectra)[*curSpec].output(cerr);
cerr.flush();

        // Statistics for best propagation
        if(m_propagation_info) {
          (*m_propagation_info)[*curSpec][0] = annIter;              (*m_propagation_info)[*curSpec][1] = bestNeigh;
          (*m_propagation_info)[*curSpec][2] = bestDeltas.empty()?0:bestDeltas.front();
          (*m_propagation_info)[*curSpec][3] = maxSpecScore>0?(int)round(10000*(bestScore/maxSpecScore)):0;
          (*m_propagation_info)[*curSpec][4] = (int)round(10000*(((float)bestMatchedPeakCount)/((float)(*m_psms_spectra)[bestNeigh].size())));
        }

        // Add propagated annotation as spectrum PSM
        (*m_psms_spectra)[*curSpec].psmList.clear();
        psmPtr propagatedPSM(new PeptideSpectrumMatch);
        aminoacids.getPeptideFromSpectrum((*m_psms_spectra)[*curSpec], propagatedPSM->m_annotation, peakTol);
        propagatedPSM->m_spectrum = &((*m_psms_spectra)[*curSpec]);
        propagatedPSM->m_scanNum = (*m_psms_spectra)[*curSpec].scan;
        (*m_psms_spectra)[*curSpec].psmList.push_back(propagatedPSM);

        // Update set of pairs used for propagations
        if(m_propagation_pairs) {
          (*m_propagation_pairs)[numMatchedPairs].spec1 = *curSpec;
          (*m_propagation_pairs)[numMatchedPairs].spec2 = bestNeigh;
          (*m_propagation_pairs)[numMatchedPairs].shift1 = m_propagation_info ? (*m_propagation_info)[*curSpec][2] : 0;
          (*m_propagation_pairs)[numMatchedPairs].score1 = bestScore;
          (*m_propagation_pairs)[numMatchedPairs].score2 = -1;   // missing value
        }

        // Set PSM midx/mp
        if(m_psms_mp and m_psms_midx) {
          (*m_psms_mp)[*curSpec].resize(3);
          (*m_psms_mp)[*curSpec][0] = (*m_psms_mp)[bestNeigh][0];
          (*m_psms_mp)[*curSpec][1] = annIter;
          (*m_psms_mp)[*curSpec][2] = bestRev ? 1 : 0;
cerr << "-- midx for spectrum " << *curSpec << " (mp = [" << (*m_psms_mp)[*curSpec][0] << "," << (*m_psms_mp)[*curSpec][1] << "," << (*m_psms_mp)[*curSpec][2] << "]:\n";

          (*m_psms_midx)[*curSpec].resize( (*m_propagation_peaks)[numMatchedPairs].size() );
          for(unsigned int p = 0; p < (*m_propagation_peaks)[numMatchedPairs].size(); p++) {
            unsigned int idx_propagation_peak = (*m_propagation_peaks)[numMatchedPairs][p][0];
            unsigned int idx_source_peak = (*m_propagation_peaks)[numMatchedPairs][p][1];
            unsigned int idx_protein_peak = (*m_psms_midx)[bestNeigh][idx_source_peak][1];
            (*m_psms_midx)[*curSpec][p].set( idx_propagation_peak , idx_protein_peak );
          }
(*m_psms_midx)[*curSpec].output(cerr);
cerr.flush();
        }
        numMatchedPairs++;

        // Check if any of the current spectrum's neighbors needs processing
        for(curNeigh=neighs[*curSpec].begin(); curNeigh!=neighs[*curSpec].end(); curNeigh++) {
          if(not annotated[*curNeigh]) {
            if(scheduled[*curNeigh]==0) {
              scheduled[*curNeigh] = annIter+1;
              annNeighs[*curNeigh].clear();
              toProcess.push_back(*curNeigh);
              lastInNextIteration = *curNeigh;
            }
            if(scheduled[*curNeigh]>annIter)
              annNeighs[*curNeigh].push_back(*curSpec);
          }
        }
      }

      if((*curSpec)==lastInIteration) {
        annIter++;
        lastInIteration = lastInNextIteration;
        lastInNextIteration=-1;
      }
    }

//    for(unsigned int i=0; i<numSpecs; i++) {
//      if(m_propagation_info and (*m_propagation_info)[i][0]>0) (*m_propagation_info)[i][1]++; // Convert indices to 1-based instead of 0-based
//      if(propagationSpecs and (*propagationSpecs)[i].size()==0)
//        (*propagationSpecs)[i]=(*m_psms_spectra)[i];  // Recover original spectra for unsuccessful propagations
//    }
    if(m_propagation_peaks) m_propagation_peaks->resize(numMatchedPairs);
    if(m_propagation_pairs) m_propagation_pairs->resize(numMatchedPairs);

    delete specSetRev, propagationSpecs;
    return true;
#endif
  }

  // -------------------------------------------------------------------------

  bool ExecSpecNetsPropagation::loadInputData(void)
  {
  	return false;
/*
  	if (ownInput) {
      if (!m_spectra)
        m_spectra = new SpecSet;
      if (!m_pairs)
      	m_pairs = new SpectrumPairSet;
    }
    m_spectra->resize(0);
    m_pairs->resize(0);

    if (ownOutput) {
      if (!m_propagation_peaks)
      	m_propagation_peaks = new SpecSet;
      if (!m_propagation_pairs)
      	m_propagation_pairs = new SpectrumPairSet;
      if (!m_propagation_info)
        m_propagation_info = new vector<vector<int> > ;
    }
    m_propagation_peaks->resize(0);
    m_propagation_pairs->resize(0);
    m_propagation_info->resize(0);

    if (m_params.exists("AMINO_ACID_MASSES")) {
      AAJumps tmpJumps(-1);
      tmpJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true); // Set global defaults for amino acid masses
    }

    if (!m_params.exists("INPUT_SPECS_PKLBIN")) {
      ERROR_MSG("Parameters are incomplete. INPUT_SPECS_PKLBIN is missing.");
      return false;
    }
    else if (m_spectra->LoadSpecSet_pklbin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str())
        <= 0 or m_spectra->size() == 0) {
      ERROR_MSG("Error reading input spectra from " <<  m_params.getValue("INPUT_SPECS_PKLBIN"));
      return false;
    }

    if (!m_params.exists("INPUT_ALIGNS")) {
      ERROR_MSG("Parameters are incomplete. INPUT_ALIGNS is missing.");
      return false;
    }
    else if ( !m_pairs->loadFromBinaryFile(m_params.getValue("INPUT_ALIGNS").c_str()) ) {
      ERROR_MSG("Error reading set of pairs from " << m_params.getValue("INPUT_ALIGNS"));
      return false;
    }

    return true;
*/
  }

  // -------------------------------------------------------------------------

  bool ExecSpecNetsPropagation::saveOutputData(void)
  {
    if (m_psms_spectra and m_params.exists("OUTPUT_SPECS_PROJ"))
    	m_psms_spectra->SaveSpecSet_pklbin(m_params.getValue("OUTPUT_SPECS_PROJ").c_str());

//    if (m_propagation_peaks and m_params.exists("OUTPUT_SPECS_MATCHIDX"))
//    	m_propagation_peaks->SaveSpecSet_pklbin(m_params.getValue("OUTPUT_SPECS_MATCHIDX").c_str());

    if (m_propagation_pairs and m_params.exists("OUTPUT_ALIGNS"))
    	m_propagation_pairs->saveToBinaryFile(m_params.getValue("OUTPUT_ALIGNS").c_str());

    if (m_propagation_info and m_params.exists("OUTPUT_ANNOTINFO"))
    	Save_binArray(m_params.getValue("OUTPUT_ANNOTINFO").c_str(), (*m_propagation_info));

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecNetsPropagation::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecNetsPropagation::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecSpecNetsPropagation::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecNetsPropagation::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecNetsPropagation::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    VALIDATE_PARAM_EXIST("MIN_PERC_EXPINT");
    VALIDATE_PARAM_EXIST("MIN_PERC_TP");

    m_isValid = true;
    return true;
  }

} // namespace specnets

