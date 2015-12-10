#include "MSBlast.h"
#include "auxfun.h"
#include "DeNovoSolutions.h"



bool getScanStringFromTitle(const string& title, string& scanStr)
{
	const size_t len = title.length();
	if (len<=7)
		return false;

	if (title[len-1] == 'a' && title[len-2] == 't' && title[len-3]=='d' && 
		title[len-6] == '.' && title[len-4]== '.')
	{
		size_t pos = len-7;
		int numDots = 0;
		while (pos>0 && numDots<2)
		{
			--pos;
			if (title[pos] == '.')
				++numDots;
		}
		if (numDots == 2)
		{
			string scanStr = title.substr(pos+1,len-7-pos);
			return true;
		}
	}
	return false;
}


void MSBlastCollector::generateMsBlastSequences(const char* msb_name,
												const vector<string>& list_vector, 
												AllScoreModels *model,
												float min_filter_prob,
												size_t maxSequencesPerSet,
												bool indOutputCumulativeProbs)
{
	config_ = model->get_config();
	PeptideRankScorer *rank_model = model->get_rank_model_ptr(1);
	PrmGraph *prm_ptr = NULL;
	vector<PrmGraph *> prm_ptrs;

	int total_benchmark = 0;
	int correct_benchmark = 0;
	double totalScore = 0.0;

	const string dnv_name = msb_name + std::string("_dnv.txt");
	const string full_name = msb_name + std::string("_full.txt");
	ofstream dnv_stream(dnv_name.c_str());
	ofstream full_stream(full_name.c_str());

	if (! dnv_stream.good() || ! full_stream.good())
		error("problem creating msb files: ",msb_name);

	for (size_t fileIdx=0; fileIdx<list_vector.size(); fileIdx++) 
	{
		const char *spectra_file = list_vector[fileIdx].c_str();

		SpectraAggregator sa;
		sa.initializeFromSpectraFilePath(spectra_file, config_);

		SpectraList sl(sa);
		sl.selectAllAggregatorHeaders();


		cout << "Processing: " << fileIdx << "  " << spectra_file << " ("
			 << sl.getNumHeaders() << " spectra)" << endl;

		if (sl.getNumHeaders()<=0)
			continue;

		for (size_t sc=0; sc<sl.getNumHeaders(); sc++)
		{
			const SingleSpectrumHeader* header = sl.getSpectrumHeader(sc);

			AnnotatedSpectrum as;
			if (! as.readSpectrum(sa, header))
			{
				SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
				nonConstHeader->setSpectraFileIndexInList(fileIdx);
				header->printStats(dnv_stream);
				dnv_stream << "#could not read spectrum correctly..." << endl;
				cout << sc << "\t";
				header->printStats(cout, false);
				cout << " #could not read spectrum correctly..." << endl;
				continue;
			}

			SingleSpectrumHeader* nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
			nonConstHeader->setSpectraFileIndexInList(fileIdx);
			
			header->printStats(dnv_stream, false);
			cout << sc << "\t";
			header->printStats(cout, false);
			
			if (as.getNumPeaks()<5)
			{
				dnv_stream << endl << "# too few peaks..." << endl;
				cout << " # too few peaks..." << endl;
				continue;
			}

			vector<SeqPath> solutions;
			solutions.clear();
			if ( as.getCharge() > model->get_max_score_model_charge())
			{
				dnv_stream << endl << "# Charge " << as.getCharge() << " not supported yet..." << endl << endl;
				cout << " # Charge " << as.getCharge() << " not supported yet..." << endl;
				continue;
			}

			bool perform_rerank=false;
			int rerank_size_idx = NEG_INF;
			int num_sols_to_generate_before_ranking=300;
			float spectrum_quality = 0.0;
		
			if (1)
			{
				vector<mass_t> pms_with_19;
				vector<int>    charges;
				pms_with_19.clear();
				charges.clear();		
			
				// output m/z and prob values for the different charge states
				vector<PmcSqsChargeRes> pmcsqs_res;
				model->selectPrecursorMassesAndCharges(config_, as, pms_with_19, charges, &pmcsqs_res);
				
				if (pmcsqs_res.size()>charges[0])
				{
					const float sqs = pmcsqs_res[charges[0]].sqs_prob;
					if (sqs<min_filter_prob)
					{
						dnv_stream << endl << "# low quality, skipping: " << sqs << endl << endl;
						cout	   << " # low quality, skipping: " << sqs << endl;
						continue;
					}
					else
						dnv_stream << "\tSQS: " << sqs << endl;
				}

				if (pmcsqs_res.size()>charges[0])
					spectrum_quality = pmcsqs_res[charges[0]].sqs_prob;

				if (prm_ptrs.size()<pms_with_19.size())
					prm_ptrs.resize(pms_with_19.size(),NULL);

				if (pms_with_19[0]<100.0)
				{
					dnv_stream << endl << "# Could not process spectrum..." << endl << endl;
					cout << " # Could not process spectrum..." << endl;
					continue;
				}
				
				const int num_rerank_sols_per_charge[] = {0,300,300,300,100,100,100,100,100,100,100,100,100,100};
				const int num_rerank_per_charge =  num_rerank_sols_per_charge[charges[0]];
				if (rank_model && num_sols_to_generate_before_ranking<num_rerank_per_charge)
						num_sols_to_generate_before_ranking=num_rerank_per_charge;
				
				generate_denovo_solutions_from_several_pms(
					prm_ptrs,
					model,
					&as,
					true, 
					num_sols_to_generate_before_ranking,
					MIN_MSB_DENOVO_SEQ_LENGTH,
					MAX_MSB_DENOVO_SEQ_LENGTH,
					pms_with_19,
					charges,
					solutions,
					false);

				// use charge of top scoring solution
				if (solutions.size()>1)
				{
					const int sol_charge = solutions[0].charge;
					for (size_t j=1; j<solutions.size(); j++)
					{
						if (solutions[j].charge != sol_charge)
						{
							if (j==solutions.size()-1)
							{
								solutions.pop_back();
							}
							else
							{
								solutions[j]=solutions[solutions.size()-1];
								solutions.pop_back();
								j--;
							}
						}
					}
				}
			}

			if (rank_model && solutions.size()>0)
			{
				rerank_size_idx = config_->calc_size_idx(solutions[0].charge,solutions[0].pm_with_19);
				if (rank_model->get_ind_part_model_was_initialized(solutions[0].charge,rerank_size_idx))
					perform_rerank=true;
			}

			vector<score_pair> score_idx_pairs;
			if (perform_rerank)
			{
				rank_model->scoreDenovoSequences(solutions, as, score_idx_pairs, rerank_size_idx);
				for (size_t i=0; i<score_idx_pairs.size(); i++)
					solutions[score_idx_pairs[i].idx].rerank_score = score_idx_pairs[i].score;
				sort(score_idx_pairs.begin(),score_idx_pairs.end());
			}
			else
			{
				score_idx_pairs.resize(solutions.size());
				for (size_t i=0; i<solutions.size(); i++)
					score_idx_pairs[i].idx=i;
			}

			// for debug
			bool had_pep = false;
			bool had_correct = false;

			if (solutions.size() == 0)
			{
				dnv_stream << endl << "# No solutions found." << endl << endl;
				cout << " # No solutions found." << endl;
			}
			else 
			{
				dnv_stream << "#Index\t";
				dnv_stream << "RnkScr\t";
				if (indOutputCumulativeProbs)
					dnv_stream << "CumProb\t";

				dnv_stream << "PnvScr\tN-Gap\tC-Gap\t[M+H]\tCharge\tSequence" << endl;

				if ( indOutputCumulativeProbs)
				{
					for (size_t i=0; i<solutions.size() && i<maxSequencesPerSet; i++)
					{
						const int idx = (perform_rerank ? score_idx_pairs[i].idx : i);
						const vector<PathPos>& positions = solutions[idx].positions;
						solutions[idx].prm_ptr->calc_amino_acid_probs(solutions[idx],i);
					}
					const int first_sol_charge = solutions[0].charge;
					const int first_sol_size_idx = config_->calc_size_idx(first_sol_charge,solutions[0].pm_with_19);
					CumulativeSeqProbModel* csp_model = (CumulativeSeqProbModel* )model->get_cumulative_seq_prob_model_ptr(0);
			
					csp_model->calc_cumulative_seq_probs(first_sol_charge, first_sol_size_idx, 
						spectrum_quality, score_idx_pairs, solutions); 
				}

				for (size_t i=0; i<solutions.size() && i<maxSequencesPerSet; i++) 
				{
					const int idx = (perform_rerank ? score_idx_pairs[i].idx : i);
					mass_t c_gap=solutions[idx].pm_with_19 - solutions[idx].c_term_mass;
					if (c_gap<24.0)
						c_gap = 0;

					dnv_stream << setprecision(3) << fixed << i << "\t";
					if (perform_rerank)
					{
						dnv_stream << score_idx_pairs[i].score << "\t";
					}
					else
						dnv_stream << -999 << "\t";

					if (indOutputCumulativeProbs)
						dnv_stream << solutions[idx].cumulative_seq_prob << "\t";

					dnv_stream << solutions[idx].path_score << "\t";
					dnv_stream << solutions[idx].n_term_mass << "\t";
					dnv_stream << c_gap << "\t";
					dnv_stream << solutions[idx].pm_with_19 << "\t";
					dnv_stream << solutions[idx].charge << "\t";
					dnv_stream << solutions[idx].seq_str;	

					if (header->getPeptideStr().length() >2)
					{
						if (solutions[idx].check_if_correct(header->getPeptideStr(),config_))
						{
							dnv_stream << " *";
							if (! had_correct)
							{
								correct_benchmark++;
								had_correct=true;
							}
						}
						had_pep=true;
					}
					dnv_stream << endl;
				}
			}
			if (had_pep) // for annotated spectra (benchmark)
				total_benchmark++;

			dnv_stream << endl;

			// create MSBlast seqeunces
			MSBSequenceSet seqSet;
			seqSet.convertSeqPathsToMSBSequences(config_, solutions, maxSequencesPerSet);

			string msbLine;
			seqSet.createMSBFullLine(msbLine);
			full_stream << msbLine << endl;
			if (seqSet.getMSBSequences().size()>0)
			{
				cout << "\t" << seqSet.getMSBSequences().size() <<" \t" << seqSet.getMSBSequences()[0].msbScore;
			}
			else
				cout << "\t" << 0;
			
			addToExistingSets(seqSet, maxSequencesPerSet*3);

			if (had_pep)
			{
				const string& pep = header->getPeptideStr();
				const vector<MSBSequence>& msbs = seqSet.getMSBSequences();
				float maxScore=0.0;
				for (size_t i=0; i<msbs.size(); i++)
				{
					float score = msbs[i].calcMatchScore(config_, pep);
					if (score > maxScore)
						maxScore = score;
				}
				totalScore += maxScore;
				cout << " s: " << maxScore;
			}
			if (had_correct)
				cout << " *";

			cout << endl;
		}
	}

	dnv_stream.close();
	full_stream.close();

	if (totalScore != 0.0)
		cout << "TotalScore = " << totalScore << endl;

	cout << endl << "Done..." << endl;
	cout << "Created" << dnv_name << endl;
	cout << "Created " << full_name << endl;

	/////////////////////////////////////////////////////////////////
	// this part works only if the spectra are annotated (benchmark)
	/////////////////////////////////////////////////////////////////
	if (total_benchmark>0)
	{
		cout << "#Correct spectra " << correct_benchmark << "/" << total_benchmark << " (" <<
			fixed << setprecision(3) << (double)correct_benchmark/(double)total_benchmark << ")" << endl;
	}

}


bool compMSBStringScore(const MSBString& lhs, const MSBString& rhs)
{
	return (lhs.msbScore > rhs.msbScore);
}


bool compMsbStringVector(const vector<MSBString>& lhs, const vector<MSBString>& rhs)
{
	return (lhs.size()>0 &&
		    (rhs.size() == 0 || lhs[0].msbScore > rhs[0].msbScore));
}

void MSBlastCollector::generateMsBlastQuery(const char* name, size_t maxQuerySize, float minMsbScore)
{
	for (size_t i=0; i<msbStringSets.size(); i++)
		sort(msbStringSets[i].begin(), msbStringSets[i].end(), compMSBStringScore);

	sort(msbStringSets.begin(), msbStringSets.end(), compMsbStringVector);

	string query_file = name + std::string("_query.txt");
	ofstream query_stream(query_file.c_str());

	if (! query_stream.good())
		error("could not open query file for writing: ",query_file.c_str());

	cout << "Creating " << query_file;
	cout.flush();

	size_t usedSize =0;
	size_t idx=0;
	while (idx<msbStringSets.size() && usedSize + 256 <maxQuerySize)
	{
		vector<MSBString>& vec = msbStringSets[idx];

		if (vec[0].msbScore<minMsbScore)
			break;

		ostringstream oss;

		oss << vec[0].fileIndex << "\t" << vec[0].scanNumber << "\t" << setprecision(3) << fixed << vec[0].mz << "\t" << vec.size() << "\t" << vec[0].msbScore << "\t";
		oss << vec[0].seqStr;
		size_t i;
		for (i=1; i<vec.size(); i++)
			oss << "-" << vec[i].seqStr;

		const string line = oss.str();
		size_t numChars=0;
		for (i=0; i<line.length(); i++)
			if (line[i]>='A' && line[i]<='Z')
				numChars++;

		usedSize += numChars;
		query_stream << line << endl;
		idx++;
	}
	query_stream.close();
	cout << " ...Done" << endl;
	cout << "Wrote " << idx << " sets of sequences (" << usedSize << " chars)" << endl;
}



void findMinScoreAndPos(const vector<MSBString>& vec, float& score, size_t& pos)
{
	if (vec.size() == 0)
	{
		score = -9999.0;
		pos   = 0;
		return;
	}

	score = vec[0].msbScore;
	pos	  = 0;

	for (size_t i=1; i<vec.size(); i++)
	{
		if (vec[i].msbScore<score)
		{
			score = vec[i].msbScore;
			pos = i;
		}
	}
}


void MSBlastCollector::addToExistingSets(MSBSequenceSet& sequenceSet, size_t maxSequencesPerSet)
{
	const vector<MSBSequence>& sequences = sequenceSet.getMSBSequences();
	vector<MSBString> msbStrings(sequences.size());

	for (size_t i=0; i<sequences.size(); i++)
	{
		msbStrings[i].seqStr      = sequences[i].makeSeqString(config_);
		msbStrings[i].msbScore    = sequences[i].msbScore;
		msbStrings[i].fileIndex   = sequences[i].fileIndex;
		msbStrings[i].mz	      = sequences[i].mz;
		msbStrings[i].scanNumber  = sequences[i].scanNumber;
	}

	size_t setIdx = MAX_SIZE_T;
	for (size_t i=0; i<msbStrings.size(); i++)
	{
		map<MSBString,size_t>::iterator it = coveredStrings.find(msbStrings[i]);
		if (it != coveredStrings.end())
		{
			setIdx = it->second;
			break;
		}
	}


	if (setIdx == MAX_SIZE_T)
	{
		setIdx = msbStringSets.size();
		msbStringSets.resize(setIdx+1);
	}

	vector<MSBString>& vec = msbStringSets[setIdx];

	// add to vector
	for (size_t i=0; i<msbStrings.size(); i++)
	{
		float  minScore=0.0;
		size_t minPos=0;

		if (vec.size() >= maxSequencesPerSet)
			findMinScoreAndPos(vec, minScore, minPos);

		for (size_t i=0; i<msbStrings.size(); i++)
		{
			const MSBString& newStr = msbStrings[i];
			if (vec.size()>=maxSequencesPerSet && newStr.msbScore <= minScore)
				continue;
			size_t j;
			for (j=0; j<vec.size(); j++)
				if (vec[j] == newStr)
					break;

			if (j<vec.size())
			{
				if (vec[j].msbScore < newStr.msbScore)
					vec[j] = newStr; // replace with higher score source
				continue;
			}

			coveredStrings[newStr]=setIdx;

			if (vec.size()<maxSequencesPerSet)
			{
				vec.push_back(newStr);
			}
			else
				vec[minPos] = newStr;

			if (vec.size() >= maxSequencesPerSet)
				findMinScoreAndPos(vec, minScore, minPos);
		}
	}
}



