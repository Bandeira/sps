/*
 * ExecMergePaired.cpp
 *
 *  Created on: Dec 12, 2011
 *      Author: aguthals
 */

#include "ExecMergePaired.h"

#include "ExecMergeConvert.h"
#include "PairedSpecSet.h"

using namespace std;

namespace specnets {
ExecMergePaired::ExecMergePaired(void) :
	m_mergedSpecset(0x0), ownInput(true), m_inputSpecset(0x0), ownOutput(true) {
	m_name = "ExecMergePaired";
	m_type = "ExecMergePaired";
}

ExecMergePaired::ExecMergePaired(const ParameterList & inputParams) :
	ExecBase(inputParams), m_mergedSpecset(0x0), ownInput(true),
			m_inputSpecset(0x0), ownOutput(true) {
	m_name = "ExecMergePaired";
	m_type = "ExecMergePaired";
}

ExecMergePaired::ExecMergePaired(const ParameterList & inputParams,
		SpecSet * pairedSpectra) :
	ExecBase(inputParams), m_mergedSpecset(0x0), ownInput(false),
			m_inputSpecset(pairedSpectra), ownOutput(true) {
	m_name = "ExecMergePaired";
	m_type = "ExecMergePaired";
}

ExecMergePaired::ExecMergePaired(const ParameterList & inputParams,
		SpecSet * pairedSpectra, SpecSet * mergedSpectra) :
	ExecBase(inputParams), m_mergedSpecset(mergedSpectra), ownInput(false),
			m_inputSpecset(pairedSpectra), ownOutput(false) {
	m_name = "ExecMergePaired";
	m_type = "ExecMergePaired";
}

ExecMergePaired::~ExecMergePaired(void) {
	if (ownInput) {
		delete m_inputSpecset;
	}
	if (ownOutput) {
		delete m_mergedSpecset;
	}
}

ExecBase * ExecMergePaired::clone(const ParameterList & inputParams) const {
	return new ExecMergePaired(inputParams);
}

bool ExecMergePaired::invoke(void) {

	if (m_params.exists("TOLERANCE_PEAK_PPM")) {
		m_inputSpecset->setPeakTolerance(m_params.getValueFloat(
				"TOLERANCE_PEAK_PPM"), true);
	} else if (m_params.exists("TOLERANCE_PEAK")) {
		m_inputSpecset->setPeakTolerance(m_params.getValueFloat(
				"TOLERANCE_PEAK"), false);
	}

	if (m_params.exists("TOLERANCE_PM_PPM")) {
		m_inputSpecset->setParentMassTolerance(m_params.getValueFloat(
				"TOLERANCE_PM_PPM"), true);
	} else if (m_params.exists("TOLERANCE_PM")) {
		m_inputSpecset->setParentMassTolerance(m_params.getValueFloat(
				"TOLERANCE_PM"), false);
	}

	float maxPeakDensity = m_params.getValueFloat("MAX_PEAK_DENSITY", 100000.0);
	DEBUG_VAR(maxPeakDensity);
	int finalStage = m_params.getValueInt("MERGE_TO_STAGE", 5);
	DEBUG_VAR(finalStage);
	PairedSpecSet paired(m_inputSpecset);
	if (finalStage >= 1) {
		paired.mergePRMsStage1(maxPeakDensity);
	}
	if (finalStage >= 2) {
		paired.mergePRMsStage2(maxPeakDensity);
	}
	if (finalStage >= 3) {
		paired.mergePRMsStage3(maxPeakDensity);
	}
	if (finalStage >= 4) {
		paired.mergePRMsStage4(maxPeakDensity);
	}
	if (finalStage >= 5) {
		paired.mergePRMsStage5(maxPeakDensity);
	}

	(*m_mergedSpecset) = *(paired.mergedSet);

	int maxSize = 0, maxIdx = -1;
	for (int i = 0; i < m_mergedSpecset->size(); i++) {
		if ((*m_mergedSpecset)[i].size() > maxSize) {
			maxSize = (*m_mergedSpecset)[i].size();
			maxIdx = i;
		}
	}

	//(*m_mergedSpecset)[maxIdx].output(cout);
	return true;
}

bool ExecMergePaired::loadInputData(void) {
	if (ownInput && (!m_inputSpecset)) {
		m_inputSpecset = new SpecSet();
	}
	if (ownOutput && (!m_mergedSpecset)) {
		m_mergedSpecset = new SpecSet();
	}

	if (m_params.exists("INPUT_PAIRED_SPECTRA")) {
		if (!ExecMergeConvert::loadSpecset(m_params.getValue(
				"INPUT_PAIRED_SPECTRA"), m_inputSpecset)) {
			return false;
		}
	} else if (m_params.exists("INPUT_CID_SPECTRA") || m_params.exists(
			"INPUT_HCD_SPECTRA") || m_params.exists("INPUT_ETD_SPECTRA")) {
		SpecSet* tempCID = tempCID = new SpecSet();
		int numPaired = 0;
		if (m_params.exists("INPUT_CID_SPECTRA")) {
			if (!ExecMergeConvert::loadSpecsetMultiple(m_params.getValue(
					"INPUT_CID_SPECTRA"), tempCID)) {
				delete tempCID;
				return false;
			}
			for (int i = 0; i < tempCID->size(); i++) {
				(*tempCID)[i].msFragType = Spectrum::FragType_CID;
			}
			numPaired++;
		}

		SpecSet* tempHCD = new SpecSet();
		if (m_params.exists("INPUT_HCD_SPECTRA")) {
			if (!ExecMergeConvert::loadSpecsetMultiple(m_params.getValue(
					"INPUT_HCD_SPECTRA"), tempHCD)) {
				delete tempCID;
				delete tempHCD;
				return false;
			}
			for (int i = 0; i < tempHCD->size(); i++) {
				(*tempHCD)[i].msFragType = Spectrum::FragType_HCD;
			}
			numPaired++;
		}

		SpecSet* tempETD = new SpecSet();
		if (m_params.exists("INPUT_ETD_SPECTRA")) {
			if (!ExecMergeConvert::loadSpecsetMultiple(m_params.getValue(
					"INPUT_ETD_SPECTRA"), tempETD)) {
				delete tempCID;
				delete tempHCD;
				delete tempETD;
				return false;
			}
			for (int i = 0; i < tempETD->size(); i++) {
				(*tempETD)[i].msFragType = Spectrum::FragType_ETD;
			}
			numPaired++;
		}

		m_inputSpecset->resize(tempCID->size() + tempHCD->size()
				+ tempETD->size());

		for (int i = 0; i < tempCID->size(); i++) {
			(*m_inputSpecset)[i * numPaired] = (*tempCID)[i];
		}
		int offset = (tempCID->size() > 0) ? 1 : 0;
		for (int i = 0; i < tempHCD->size(); i++) {
			(*m_inputSpecset)[(i * numPaired) + offset] = (*tempHCD)[i];
		}
		offset = (tempHCD->size() > 0) ? offset + 1 : offset;
		for (int i = 0; i < tempETD->size(); i++) {
			(*m_inputSpecset)[(i * numPaired) + offset] = (*tempETD)[i];
		}
		DEBUG_VAR(m_inputSpecset->size());
		delete tempCID;
		delete tempHCD;
		delete tempETD;
	} else {
		WARN_MSG("No input spectra parameter!");
	}
	return true;
}

bool ExecMergePaired::saveOutputData(void) {
	if (m_params.exists("OUTPUT_MERGED_SPECTRA")) {
		if (!ExecMergeConvert::saveSpecset(m_params.getValue(
				"OUTPUT_MERGED_SPECTRA"), m_mergedSpecset)) {
			return false;
		}
	} else {
		WARN_MSG("No output spectra parameter!");
	}
	return true;
}

bool ExecMergePaired::saveInputData(std::vector<std::string> & filenames) {
	return true;
}

bool ExecMergePaired::loadOutputData(void) {

	return true;
}

std::vector<ExecBase *> const & ExecMergePaired::split(int numSplit) {

	m_subModules.resize(0);

	return m_subModules;
}

bool ExecMergePaired::merge(void) {

	return true;
}

bool ExecMergePaired::validateParams(std::string & error) {
	m_isValid = false;
	m_isValid = true;
	return true;
}
}
