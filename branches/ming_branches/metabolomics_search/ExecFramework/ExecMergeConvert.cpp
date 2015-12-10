/*
 * ExecMergeConvert.cpp
 *
 *  Created on: Oct 5, 2011
 *      Author: aguthals
 */

#include <algorithm>

#include "ExecMergeConvert.h"

using namespace std;

namespace specnets {

bool ExecMergeConvert::loadSpecsetMultiple(const string& specFilePaths,
		SpecSet* spectra, const char* separator) {
	list<string> fileList;
	int idxUse = 0;
	splitText(specFilePaths.c_str(), fileList, separator);
	SpecSet* tempSpecs = new SpecSet();
	for (list<string>::iterator fileIt = fileList.begin(); fileIt
			!= fileList.end(); fileIt++) {
		if (!loadSpecset(*fileIt, tempSpecs)) {
			goto load_fail;
		}
		spectra->resize(spectra->size() + tempSpecs->size());

		for (int i = 0; i < tempSpecs->size(); i++) {
			(*spectra)[idxUse] = (*tempSpecs)[i];
			++idxUse;
		}
	}

	delete tempSpecs;

	return true;

	load_fail:

	delete tempSpecs;

	return false;
}

bool ExecMergeConvert::loadSpecset(const string& specFilePath, SpecSet* spectra) {

	// extract the file extension to determine its format
	vector<string> pathParts;
	const char* filename = specFilePath.c_str();
	splitText(filename, pathParts, ".");
	string extension = pathParts[pathParts.size() - 1];

	// transform file extension to lower case for easy comparisons
	std::transform(extension.begin(), extension.end(), extension.begin(),
			::tolower);

	if (extension == "mgf") {
		if (!spectra->LoadSpecSet_mgf(filename)) {
			goto load_fail;
		}
	} else if (extension == "pklbin") {
		if (!spectra->LoadSpecSet_pklbin(filename)) {
			goto load_fail;
		}
	} else if (extension == "prms" or extension == "prmsv3") {
		if (!spectra->LoadSpecSet_prmsv3(filename)) {
			goto load_fail;
		}
	} else if (extension == "pkl") {
		if (!spectra->LoadSpecSet_pkl(filename)) {
			goto load_fail;
		}
	} else if (extension == "mzxml") {
		vector<short> msLevel;
		if (LoadMzxml(filename, *spectra, &msLevel, 2) == 0) {
			goto load_fail;
		}
	} else { // try loading a directory of pkl or mic files if all else fails
		vector<string> possibleDirFiles = directoryContents(specFilePath, "",
				"", true);

		if (possibleDirFiles.size() == 0) {
			ERROR_MSG("/'" << specFilePath << "\' is an empty directory, " <<
					"is an un-reachable directory, or is not named by a " <<
					"supported spectrum-set file format (.mgf, .pklbin, .prms, .pkl, or .mzxml)");
			goto load_fail;
		}

		int idxUse = 0;
		for (int i = 0; i < possibleDirFiles.size(); i++) {
			string fullPath = specFilePath + "/" + possibleDirFiles[i];
			if (possibleDirFiles[i].find(".pkl") == string::npos
					and possibleDirFiles[i].find(".mic") == string::npos) {
				WARN_MSG("Skipping file \'" << fullPath << "\'");
				continue;
			}

			SpecSet tempSpecs;
			if (!tempSpecs.LoadSpecSet_pkl_mic(fullPath.c_str())) {
				ERROR_MSG("Failed to load file \'" << fullPath << "\' in pkl or mic format!");
				goto load_fail;
			}

			spectra->resize(idxUse + tempSpecs.size());
			for (int j = 0; j < tempSpecs.size(); j++) {
				(*spectra)[idxUse] = tempSpecs[j];
				++idxUse;
			}
		}
	}

	DEBUG_MSG("Loaded " << spectra->size() << " spectra from \'" << filename << "\'");

	return true;

	load_fail:

	ERROR_MSG("Failed to load from \'" << filename << "\'!");

	return false;

}

bool ExecMergeConvert::saveSpecset(const string& specFilePath, SpecSet* spectra) {

	// extract the file extension to determine its format
	vector<string> pathParts;
	const char* filename = specFilePath.c_str();
	splitText(filename, pathParts, ".");
	string extension = pathParts[pathParts.size() - 1];

	// transform file extension to lower case for easy comparisons
	std::transform(extension.begin(), extension.end(), extension.begin(),
			::tolower);

	DEBUG_MSG("Saving " << spectra->size() << " spectra to " << filename << " in \'" << extension << "\' format ...");

	if (extension == "mgf") {
		if (!spectra->SaveSpecSet_mgf(filename)) {
			goto save_fail;
		}
	} else if (extension == "pklbin") {
		if (!spectra->SaveSpecSet_pklbin(filename)) {
			goto save_fail;
		}
	} else if (extension == "pkl") {
		if (!spectra->SaveSpecSet_pkl(filename)) {
			goto save_fail;
		}
	} else {
		ERROR_MSG("Unrecognized file extension \'" << extension << "\'");
		goto save_fail;
	}

	return true;

	save_fail:

	ERROR_MSG("Failed to save to \'" << filename << "\'!");
	return false;

}

pair<bool, list<string> > ExecMergeConvert::getFileList(
		const string& listFilePath) {
	BufferedLineReader blr;
	list<string> fileList;

	if (blr.Load(listFilePath.c_str()) <= 0) {
		ERROR_MSG("Could not load list of files from \'" << listFilePath << "\'");
		return pair<bool, list<string> > (false, fileList);
	}

	for (int lineIdx = 0; lineIdx < blr.size(); lineIdx++) {
		string file = blr.getline(lineIdx);
		file = replaceAll(file, "\n", "");
		if (file.length() > 0) {
			fileList.push_back(file);
		}
	}
	blr.reset();

	return pair<bool, list<string> > (true, fileList);
}

ExecMergeConvert::ExecMergeConvert(void) :
	m_inputSpecsetIdx(0x0), ownInput(true), m_outputSpecset(0x0), ownOutput(
			true) {
	m_name = "ExecMergeConvert";
	m_type = "ExecMergeConvert";
	m_recordedInput.resize(0);
	m_recordedOutput.resize(0);
	m_separator = ";";
}

ExecMergeConvert::ExecMergeConvert(const ParameterList & inputParams) :

	ExecBase(inputParams), m_inputSpecsetIdx(0x0), ownInput(true),
			m_outputSpecset(0x0), ownOutput(true) {

	m_name = "ExecMergeConvert";
	m_type = "ExecMergeConvert";
	m_recordedInput.resize(0);
	m_recordedOutput.resize(0);
	m_separator = ";";
}

ExecMergeConvert::ExecMergeConvert(const ParameterList & inputParams, vector<
		pair<int, int> >* inputSpecsetIdx) :

	ExecBase(inputParams), m_inputSpecsetIdx(inputSpecsetIdx), ownInput(false),
			m_outputSpecset(0x0), ownOutput(true) {

	m_name = "ExecMergeConvert";
	m_type = "ExecMergeConvert";
	m_recordedInput.resize(0);
	m_recordedOutput.resize(0);
	m_separator = ";";
}

ExecMergeConvert::ExecMergeConvert(const ParameterList & inputParams, vector<
		pair<int, int> >* inputSpecsetIdx, SpecSet * outputSpecset) :

	ExecBase(inputParams), m_inputSpecsetIdx(inputSpecsetIdx), ownInput(false),
			m_outputSpecset(outputSpecset), ownOutput(true) {

	m_name = "ExecMergeConvert";
	m_type = "ExecMergeConvert";
	m_recordedInput.resize(0);
	m_recordedOutput.resize(0);
	m_separator = ";";
}

ExecMergeConvert::~ExecMergeConvert(void) {
	if (ownInput) {
		delete m_inputSpecsetIdx;
	}
	if (ownOutput) {
		delete m_outputSpecset;
	}
}

ExecBase * ExecMergeConvert::clone(const ParameterList & inputParams) const {
	return new ExecMergeConvert(inputParams);
}

bool ExecMergeConvert::invoke(void) {
	const char* separator = m_separator.c_str();

	if (m_params.exists("TOLERANCE_PEAK_PPM")) {
		float pktol = m_params.getValueFloat("TOLERANCE_PEAK_PPM");
		m_outputSpecset->setPeakTolerance(pktol, true);
		m_recordedProcessing << "Setting peak tolerance to " << pktol
				<< " ppm\n";
		DEBUG_MSG("Setting peak tolerance to " << pktol
				<< " ppm");
	}

	if (m_params.exists("TOLERANCE_PEAK")) {
		float pktol = m_params.getValueFloat("TOLERANCE_PEAK");
		m_outputSpecset->setPeakTolerance(pktol, false);
		m_recordedProcessing << "Setting peak tolerance to " << pktol
				<< " Da\n";
		DEBUG_MSG("Setting peak tolerance to " << pktol
				<< " Da");
	}

	if (m_params.exists("TOLERANCE_PM_PPM")) {
		float pktol = m_params.getValueFloat("TOLERANCE_PM_PPM");
		m_outputSpecset->setParentMassTolerance(pktol, true);
		m_recordedProcessing << "Setting parent mass tolerance to " << pktol
				<< " ppm\n";
		DEBUG_MSG("Setting parent mass tolerance to " << pktol
				<< " ppm");
	}

	if (m_params.exists("TOLERANCE_PM")) {
		float pktol = m_params.getValueFloat("TOLERANCE_PM");
		m_outputSpecset->setParentMassTolerance(pktol, false);
		m_recordedProcessing << "Setting parent mass tolerance to " << pktol
				<< " Da\n";
		DEBUG_MSG("Setting parent mass tolerance to " << pktol
				<< " Da");
	}

	if (m_outputSpecset->size() == 0) {
		WARN_MSG("Detected 0 spectra, skipping ExecMergeConvert::invoke()");
		return true;
	}

	if (m_params.exists("FILTER_CHARGE")) {
		if (m_params.exists("FILTER_OUTPUT_CHARGES")) {
			WARN_MSG("Ignoring FILTER_CHARGE, using FILTER_OUTPUT_CHARGES instead");
		}
		m_params.addIfDoesntExist("FILTER_OUTPUT_CHARGES", m_params.getValue(
				"FILTER_CHARGE"));
	}

	if (m_params.exists("FILTER_INPUT_CHARGES")) {
		list<string> chargeRanges;
		splitText(m_params.getValue("FILTER_INPUT_CHARGES").c_str(),
				chargeRanges, separator);

		if (chargeRanges.size() != m_inputSpecsetIdx->size()) {
			ERROR_MSG("Found " << chargeRanges.size() << " input charge filters for " << m_inputSpecsetIdx->size() << " input files, must have same for both");
			return false;
		}

		int rangeIdx = 0;
		set<short> locRangeVals;
		for (list<string>::iterator rangeIt = chargeRanges.begin(); rangeIt
				!= chargeRanges.end(); rangeIt++) {

			locRangeVals.clear();
			if (!getRanges(rangeIt->c_str(), locRangeVals)) {
				ERROR_MSG("Failed to parse charge range \'" << *rangeIt << "\'");
				return false;
			}

			m_recordedProcessing << "Keeping spectra of precursor charge "
					<< *rangeIt << " from indices "
					<< (*m_inputSpecsetIdx)[rangeIdx].first << " to "
					<< (*m_inputSpecsetIdx)[rangeIdx].second << "\n";
			DEBUG_MSG("Keeping spectra of precursor charge "
					<< *rangeIt << " from indices "
					<< (*m_inputSpecsetIdx)[rangeIdx].first << " to "
					<< (*m_inputSpecsetIdx)[rangeIdx].second);

			for (int i = (*m_inputSpecsetIdx)[rangeIdx].first; i
					<= (*m_inputSpecsetIdx)[rangeIdx].second; i++) {
				if (locRangeVals.count((*m_outputSpecset)[i].parentCharge) == 0) {
					(*m_outputSpecset)[i].resize(0);
				}
			}
			++rangeIdx;
		}
	}

	if (m_params.getValueInt("FIX_CHARGE_ZEROS", 0) == 1) {
		DEBUG_MSG("Guessing the charge for any charge 0 spectrum where the parent mass is less than the last mass in the spectrum");
		for (int i = 0; i < m_outputSpecset->size(); i++) {
			if ((*m_outputSpecset)[i].parentCharge == 0
					&& (*m_outputSpecset)[i].parentMass + AAJumps::massMH
							< (*m_outputSpecset)[i][(*m_outputSpecset)[i].size()
									- 1][0]) {

				while ((*m_outputSpecset)[i].parentMass + AAJumps::massMH
						< (*m_outputSpecset)[i][(*m_outputSpecset)[i].size()
								- 1][0]) {
					(*m_outputSpecset)[i].parentCharge += 1;
					(*m_outputSpecset)[i].parentMass
							= ((*m_outputSpecset)[i].parentMass
									* (float) (*m_outputSpecset)[i].parentCharge)
									- (AAJumps::massHion
											* ((float) (*m_outputSpecset)[i].parentCharge
													- 1.0));
				}

				m_recordedProcessing << "Assigned charge "
						<< (*m_outputSpecset)[i].parentCharge
						<< " to spectrum " << i << "\n";

				DEBUG_MSG("Assigned charge "
						<< (*m_outputSpecset)[i].parentCharge << " to spectrum " << i);

			}
		}
	}

	if (m_params.exists("CONVERT_CHARGE")) {

		list<string> chargeRanges;
		splitText(m_params.getValue("CONVERT_CHARGE").c_str(), chargeRanges,
				separator);

		map<short, short> chargeRangeVals;
		const char* rangeSep = ">";

		for (list<string>::iterator rangeIt = chargeRanges.begin(); rangeIt
				!= chargeRanges.end(); rangeIt++) {
			vector<string> rangeConvert;
			splitText(rangeIt->c_str(), rangeConvert, rangeSep);

			if (rangeConvert.size() != 2) {
				ERROR_MSG("Failed to parse charge range conversion \'"
						<< rangeIt->c_str() << "\', must be something like \'4-7" << rangeSep << "3\'");
				return false;
			}

			set<short> rangeVals;
			if (!getRanges(rangeConvert[0].c_str(), rangeVals)) {
				ERROR_MSG("Failed to parse charge range \'" << rangeConvert[0] << "\'");
				return false;
			}
			short newCharge = (short) getInt(rangeConvert[1].c_str());
			m_recordedProcessing << "Converting all charge " << rangeConvert[0]
					<< " spectra to charge " << newCharge << "\n";
			DEBUG_MSG("Converting all charge " << rangeConvert[0] << " spectra to charge " << newCharge);

			for (set<short>::iterator cIt = rangeVals.begin(); cIt
					!= rangeVals.end(); cIt++) {
				chargeRangeVals[*cIt] = newCharge;
			}
		}

		for (int i = 0; i < m_outputSpecset->size(); i++) {
			if (chargeRangeVals.count((*m_outputSpecset)[i].parentCharge) > 0) {
				float prevCharegF = (float) (*m_outputSpecset)[i].parentCharge;
				short newCharge =
						chargeRangeVals[(*m_outputSpecset)[i].parentCharge];
				float newChargeF = (float) newCharge;
				(*m_outputSpecset)[i].parentMZ
						= ((*m_outputSpecset)[i].parentMass + ((newChargeF
								- 1.0) * AAJumps::massHion)) / newChargeF;
				(*m_outputSpecset)[i].parentCharge = newCharge;
			}
		}
	}

	//float samePrecTol = 0.001;

	if (m_params.getValueInt("MERGE_CONSECUTIVE", 0) > 1) {
		int allowedConsec = m_params.getValueInt("MERGE_CONSECUTIVE");
		SpecSet mergedSpecs(m_outputSpecset->size());

		m_recordedProcessing << "Merging every " << allowedConsec
				<< " consecutive spectra from the same precursor\n";
		DEBUG_MSG("Merging every " << allowedConsec << " consecutive spectra from the same precursor");

		int rootSpecIdx = 0;
		mergedSpecs[rootSpecIdx] = (*m_outputSpecset)[0];
		float lastPM = (*m_outputSpecset)[0].parentMass;
		int curConsec = 0;
		for (int i = 1; i < m_outputSpecset->size(); i++) {
			if (curConsec < allowedConsec) {

				mergedSpecs[rootSpecIdx].mergeClosestPeaks(
						(*m_outputSpecset)[i], 2);
				++curConsec;
			} else {
				/*
				 if (curConsec < allowedConsec) {
				 WARN_MSG("Only merged " << curConsec << " consecutive spectra from idx "
				 << rootSpecIdx << " to " << i - 1 << ". Parent masses not equal ("
				 << lastPM << " and " << (*m_outputSpecset)[i].parentMass << ")");
				 }
				 */
				curConsec = 0;
				++rootSpecIdx;
				mergedSpecs[rootSpecIdx] = (*m_outputSpecset)[i];
				lastPM = (*m_outputSpecset)[i].parentMass;
			}
		}
		mergedSpecs.resize(rootSpecIdx + 1);

		*m_outputSpecset = mergedSpecs;
	}

	if (m_params.getValueInt("MERGE_PARALLEL", 0) > 1) {
		int allowedMerged = m_params.getValueInt("MERGE_PARALLEL");

		m_recordedProcessing
				<< "Merging spectra at parallel indices in the first "
				<< allowedMerged << " files\n";
		DEBUG_MSG("Merging spectra at parallel indices in the first " << allowedMerged << " files");

		int expectedParallel = (*m_inputSpecsetIdx)[0].second + 1;

		for (int i = 1; i < allowedMerged; i++) {
			if (i >= m_inputSpecsetIdx->size()) {
				ERROR_MSG("Only loaded " << m_inputSpecsetIdx->size() <<
						" files, need at least " << allowedMerged << " to merge parallel indices!");
				return false;
			}

			int detectedParallel = ((*m_inputSpecsetIdx)[i].second
					- (*m_inputSpecsetIdx)[i].first) + 1;

			if (detectedParallel != expectedParallel) {
				ERROR_MSG("Detected " << expectedParallel << " spectra in the first file, need exactly "
						<< expectedParallel << " in the first " << allowedMerged <<
						" files to merge parallel indices (found " << detectedParallel << " spectra in file " << i+1 << ")");
				return false;
			}
		}

		SpecSet mergedSpecs(expectedParallel);

		for (int i = 0; i < expectedParallel; i++) {
			mergedSpecs[i] = (*m_outputSpecset)[i];
			int parallelIdx = i + expectedParallel;
			int numMerged = 1;
			while (numMerged < allowedMerged) {
				/*
				 if (!MZRange::EqualWithinRange(mergedSpecs[i].parentMass,
				 allSpecs[parallelIdx].parentMass,
				 0.001)) {
				 cerr << "ERROR: Spectra at parallel index " << i
				 << " have different parent masses ("
				 << mergedSpecs[i].parentMass << ") and ("
				 << allSpecs[parallelIdx].parentMass << ")\n";
				 return -1;
				 }
				 */

				mergedSpecs[i].mergeClosestPeaks(
						(*m_outputSpecset)[parallelIdx], 2);
				parallelIdx += expectedParallel;
				++numMerged;
			}
		}
		*m_outputSpecset = mergedSpecs;
	}

	if (m_params.exists("FILTER_OUTPUT_CHARGES")) {
		list<string> chargeRanges;
		splitText(m_params.getValue("FILTER_OUTPUT_CHARGES").c_str(),
				chargeRanges, separator);
		set<short> rangeVals;

		for (list<string>::iterator rangeIt = chargeRanges.begin(); rangeIt
				!= chargeRanges.end(); rangeIt++) {

			set<short> locRangeVals;
			if (!getRanges(rangeIt->c_str(), rangeVals)) {
				ERROR_MSG("Failed to parse charge range \'" << *rangeIt << "\'");
				return false;
			}
			rangeVals.insert(locRangeVals.begin(), locRangeVals.end());
		}

		m_recordedProcessing << "Keeping spectra of precursor charge "
				<< m_params.getValue("FILTER_OUTPUT_CHARGES") << "\n";
		DEBUG_MSG("Keeping spectra of precursor charge " << m_params.getValue("FILTER_OUTPUT_CHARGES"));

		for (int i = 0; i < m_outputSpecset->size(); i++) {
			if (rangeVals.count((*m_outputSpecset)[i].parentCharge) == 0) {
				(*m_outputSpecset)[i].resize(0);
			}
		}
	}

	if (m_params.exists("MASS_OFFSET")) {
		float massOffset = m_params.getValueFloat("MASS_OFFSET");

		DEBUG_MSG("Adding " << massOffset << " to every mass");
		m_recordedProcessing << "Adding " << massOffset << " to every mass\n";

		for (int i = 0; i < m_outputSpecset->size(); i++) {
			for (int j = 0; j < (*m_outputSpecset)[i].size(); j++) {
				(*m_outputSpecset)[i][j][0] += massOffset;
			}
		}
	}

	if (m_params.exists("ACTIVATION")) {
		string activStr = m_params.getValue("ACTIVATION");
		std::transform(activStr.begin(), activStr.end(), activStr.begin(),
				::toupper);
		Spectrum tempSpec;

		if (activStr == "CID") {
			tempSpec.msFragType = Spectrum::FragType_CID;
		} else if (activStr == "HCD") {
			tempSpec.msFragType = Spectrum::FragType_HCD;
		} else if (activStr == "ETD") {
			tempSpec.msFragType = Spectrum::FragType_ETD;
		} else {
			ERROR_MSG("Unrecognized activation method \'" << activStr << "\'");
			return false;
		}

		m_recordedProcessing << "Setting activation method to \'" << activStr
				<< "\' for all spectra\n";
		DEBUG_MSG("Setting activation method to \'" << activStr << "\' for all spectra");

		for (int i = 0; i < m_outputSpecset->size(); i++) {
			(*m_outputSpecset)[i].msFragType = tempSpec.msFragType;
		}
	}

	if (m_params.getValueInt("RANK_FILTER", 0) > 0) {
		int rankFiltK = m_params.getValueInt("RANK_FILTER");
		float windowRadius = m_params.getValueFloat("RANK_FILTER_RADIUS",
				AAJumps::minAAmass - 1.0);

		m_recordedProcessing
				<< "Applying rank filtering, only keeping the top "
				<< rankFiltK << " peaks(s) in every +/- " << windowRadius
				<< " Da window\n";

		DEBUG_MSG("Applying rank filtering, only keeping the top " << rankFiltK <<
				" peaks(s) in every +/- " << windowRadius << " Da window");

		for (int i = 0; i < m_outputSpecset->size(); i++) {
			(*m_outputSpecset)[i].rankFilterPeaks(rankFiltK, windowRadius);
		}

	}

	if (m_params.exists("REVERSE_OFFSET")) {
		float revOffset = m_params.getValueFloat("REVERSE_OFFSET");

		for (int i = 0; i < m_outputSpecset->size(); i++) {
			(*m_outputSpecset)[i].reverse(revOffset);
		}
	}

	if (m_params.getValueInt("SET_SCAN_NUMS", 0) > 0) {
		DEBUG_MSG("Setting the scan number of each spectrum to its one-based index");
		m_recordedProcessing
				<< "Setting the scan number of each spectrum to its index\n";
		for (int i = 0; i < m_outputSpecset->size(); i++) {
			(*m_outputSpecset)[i].scan = (unsigned int) (i + 1);
		}
	}

	return true;
}

bool ExecMergeConvert::loadInputData(void) {
	if (ownInput && (!m_inputSpecsetIdx)) {
		m_inputSpecsetIdx = new vector<pair<int, int> > ;
	}
	if (ownOutput && (!m_outputSpecset)) {
		m_outputSpecset = new SpecSet();
	}
	m_inputSpecsetIdx->resize(0);
	m_outputSpecset->resize(0);
	m_recordedInput.resize(0);

	const char* separator = m_separator.c_str();

	bool loadSpecResult;
	int idxUse = 0;

	SpecSet* locSpecs = new SpecSet();

	// Load spectra from text files of paths
	if (m_params.exists("INPUT_SPECTRA_LISTS")) {
		list<string> listFilesList;
		pair<bool, list<string> > fileListResult;
		splitText(m_params.getValue("INPUT_SPECTRA_LISTS").c_str(),
				listFilesList, separator);
		for (list<string>::iterator fileListIt = listFilesList.begin(); fileListIt
				!= listFilesList.end(); fileListIt++) {
			fileListResult = getFileList(*fileListIt);
			if (!fileListResult.first) {
				return false;
			}
			for (list<string>::iterator fileIt = fileListResult.second.begin(); fileIt
					!= fileListResult.second.end(); fileIt++) {

				loadSpecResult = loadSpecset(*fileIt, locSpecs);

				if (!loadSpecResult) {
					delete locSpecs;
					return false;
				}
				m_recordedInput.push_back(pair<string, int> (*fileIt,
						locSpecs->size()));
				m_outputSpecset->resize(idxUse + locSpecs->size());
				m_inputSpecsetIdx->push_back(pair<int, int> (idxUse, idxUse
						+ locSpecs->size() - 1));
				for (int i = 0; i < locSpecs->size(); i++) {
					(*m_outputSpecset)[idxUse] = (*locSpecs)[i];
					++idxUse;
				}
			}
		}
	}

	// Load spectra from direct paths
	if (m_params.exists("INPUT_SPECTRA")) {
		list<string> listFiles;
		splitText(m_params.getValue("INPUT_SPECTRA").c_str(), listFiles,
				separator);
		for (list<string>::iterator fileIt = listFiles.begin(); fileIt
				!= listFiles.end(); fileIt++) {
			loadSpecResult = loadSpecset(*fileIt, locSpecs);
			if (!loadSpecResult) {
				delete locSpecs;
				return false;
			}

			m_recordedInput.push_back(pair<string, int> (*fileIt,
					locSpecs->size()));
			m_outputSpecset->resize(idxUse + locSpecs->size());
			m_inputSpecsetIdx->push_back(pair<int, int> (idxUse, idxUse
					+ locSpecs->size() - 1));
			for (int i = 0; i < locSpecs->size(); i++) {
				(*m_outputSpecset)[idxUse] = (*locSpecs)[i];
				++idxUse;
			}
		}
	}

	delete locSpecs;

	return true;
}

bool ExecMergeConvert::saveInputData(std::vector<std::string> & filenames) {
	return true;
}

bool ExecMergeConvert::saveOutputData(void) {
	const char* separator = m_separator.c_str();

	bool splitOut = m_params.getValueInt("SPLIT_OUTPUT", 0) > 0;

	int curIdx = 0, specsPerFile = 0;

	// Save merged spectra as single SpecSet
	if (m_params.exists("OUTPUT_SPECTRA")) {
		list<string> listFiles;
		splitText(m_params.getValue("OUTPUT_SPECTRA").c_str(), listFiles,
				separator);

		if (splitOut) {
			int numOut = (int) listFiles.size();
			specsPerFile = ((int) m_outputSpecset->size()) / (numOut);
			int remainder = ((int) m_outputSpecset->size()) % (numOut);
			if (remainder > 0) {
				++specsPerFile;
			}
			DEBUG_MSG("Spltting output spectra into " << numOut << " sets of ~" << specsPerFile << " spectra");
			m_recordedProcessing << "Spltting output spectra into " << numOut
					<< " sets of ~" << specsPerFile << " spectra\n";
		}
		for (list<string>::iterator fileIt = listFiles.begin(); fileIt
				!= listFiles.end(); fileIt++) {

			if (splitOut) {
				int startIdx = curIdx;
				int endIdx = min(curIdx + specsPerFile - 1,
						((int) m_outputSpecset->size()) - 1);
				SpecSet tempSpecs(endIdx - startIdx + 1);
				int locIdx = 0;
				for (int i = startIdx; i <= endIdx; i++) {
					tempSpecs[locIdx] = (*m_outputSpecset)[i];
					++locIdx;
				}
				if (!saveSpecset(*fileIt, &tempSpecs)) {
					return false;
				}
				m_recordedOutput.push_back(pair<string, int> (*fileIt,
						tempSpecs.size()));
				curIdx = endIdx + 1;
			} else {
				if (!saveSpecset(*fileIt, m_outputSpecset)) {
					return false;
				}
				m_recordedOutput.push_back(pair<string, int> (*fileIt,
						m_outputSpecset->size()));
			}
		}
	}

	// Save a record of what was loaded, merged, filtered, and/or saved
	if (ownInput and ownOutput) {
		stringstream filePath;
		filePath << m_params.m_inputFilePath << ".info";

		DEBUG_MSG("Saving activity record to \'" << filePath.str() << "\' ...");
		bool success = saveActivity(filePath.str());
	}

	return true;
}

bool ExecMergeConvert::loadOutputData(void) {

	return true;
}

std::vector<ExecBase *> const & ExecMergeConvert::split(int numSplit) {

	m_subModules.resize(0);

	return m_subModules;
}

bool ExecMergeConvert::merge(void) {

	return true;
}

bool ExecMergeConvert::validateParams(std::string & error) {
	m_isValid = false;

	VALIDATE_PARAM_EXIST("OUTPUT_SPECTRA");

	if (!m_params.exists("INPUT_SPECTRA_LISTS")) {
		VALIDATE_PARAM_EXIST("INPUT_SPECTRA");
	} else {
		VALIDATE_PARAM_EXIST("INPUT_SPECTRA_LISTS");
	}

	m_isValid = true;
	return true;
}

bool ExecMergeConvert::saveActivity(const string& infoFilePath) {

	FILE* output = fopen(infoFilePath.c_str(), "w");
	if (output == NULL) {
		ERROR_MSG("Failed to save activity record to file \'" << infoFilePath << "\'");
		return false;
	}

	fprintf(output, "<Input File>\t<# Spectra Loaded>\n");
	for (int i = 0; i < m_recordedInput.size(); i++) {
		fprintf(output, "%s\t%d\n", m_recordedInput[i].first.c_str(),
				m_recordedInput[i].second);
	}

	fprintf(output, "\n%s\n", m_recordedProcessing.str().c_str());

	fprintf(output, "<Output File>\t<# Spectra Saved>\n");
	for (int i = 0; i < m_recordedOutput.size(); i++) {
		fprintf(output, "%s\t%d\n", m_recordedOutput[i].first.c_str(),
				m_recordedOutput[i].second);
	}

	fclose(output);

	return true;
}

}
