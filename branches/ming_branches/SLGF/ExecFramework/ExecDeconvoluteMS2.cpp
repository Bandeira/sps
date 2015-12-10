/*
 * ExecDeconvoluteMS2.cpp
 *
 *  Created on: Feb 28, 2011
 *      Author: aguthals
 */

#include "ExecMergeConvert.h"

// Header includes
#include "ExecDeconvoluteMS2.h"

using namespace specnets;
using namespace std;

namespace specnets {

ExecDeconvoluteMS2::ExecDeconvoluteMS2(void) :
	m_inputSpectra(0x0), m_inputEnv(0x0), m_outputSpectra(0x0), ownInput(true),
			ownOutput(true) {
	m_name = "ExecDeconvoluteMS2";
	m_type = "ExecDeconvoluteMS2";
}

ExecDeconvoluteMS2::ExecDeconvoluteMS2(const ParameterList & inputParams) :

	ExecBase(inputParams), m_inputSpectra(0x0), m_inputEnv(0x0),
			m_outputSpectra(0x0), ownInput(true), ownOutput(true) {

	m_name = "ExecDeconvoluteMS2";
	m_type = "ExecDeconvoluteMS2";
}

ExecDeconvoluteMS2::ExecDeconvoluteMS2(const ParameterList & inputParams,
		SpecSet * inputSpectra, IsoEnvelope * m_inputEnv) :
	ExecBase(inputParams), m_inputSpectra(inputSpectra),
			m_inputEnv(m_inputEnv), m_outputSpectra(0x0), ownInput(false),
			ownOutput(true) {

	m_name = "ExecDeconvoluteMS2";
	m_type = "ExecDeconvoluteMS2";
}

ExecDeconvoluteMS2::ExecDeconvoluteMS2(const ParameterList & inputParams,
		SpecSet * inputSpectra, IsoEnvelope * m_inputEnv,
		SpecSet * outputSpectra) :

	ExecBase(inputParams), m_inputSpectra(inputSpectra),
			m_inputEnv(m_inputEnv), m_outputSpectra(outputSpectra), ownInput(
					false), ownOutput(false) {

	m_name = "ExecDeconvoluteMS2";
	m_type = "ExecDeconvoluteMS2";
}

ExecDeconvoluteMS2::~ExecDeconvoluteMS2(void) {
	if (ownInput) {
		delete m_inputSpectra;
		delete m_inputEnv;
	}
	if (ownOutput) {
		delete m_outputSpectra;
	}
}

void ExecDeconvoluteMS2::setOwnInput(bool _ownInput) {
	ownInput = _ownInput;
}

ExecBase * ExecDeconvoluteMS2::clone(const ParameterList & inputParams) const {
	return new ExecDeconvoluteMS2(inputParams);
}

bool ExecDeconvoluteMS2::invoke(void) {

	if (m_params.exists("TOLERANCE_PPM")) {
		m_outputSpectra->setPeakTolerance(m_params.getValueFloat(
				"TOLERANCE_PPM"), true);
	} else if (m_params.exists("TOLERANCE_PEAK")) {
		m_outputSpectra->setPeakTolerance(m_params.getValueFloat(
				"TOLERANCE_PEAK"), false);
	}

	float threshold = m_params.getValueFloat("MAX_KLDiv", 0.50);

	int start = (m_params.exists("IDX_START")) ? m_params.getValueInt(
			"IDX_START") : 0;
	int end = (m_params.exists("IDX_END")) ? m_params.getValueInt("IDX_END")
			: m_inputSpectra->size() - 1;

	m_outputSpectra->resize(m_inputSpectra->size());
	DeconvSpectrum workingSpec;
	for (int i = start; i <= end; i++) {
		workingSpec = (DeconvSpectrum &) (*m_inputSpectra)[i];
		workingSpec.AssignChargesKLDiv(m_inputEnv, threshold);
		workingSpec.ConvertPeaksChargeOne((*m_outputSpectra)[i]);
		//(*m_outputSpectra)[i].output(cout);
		//break;
	}
	return true;
}

bool ExecDeconvoluteMS2::loadInputData(void) {
	ownInput = true;
	m_inputSpectra = new SpecSet();
	m_outputSpectra = new SpecSet();
	m_inputEnv = new IsoEnvelope();

	//---------------------------------
	// Load spectrum data
	//---------------------------------
	if (m_params.exists("INPUT_SPECTRA")) {
		if (!ExecMergeConvert::loadSpecsetMultiple(m_params.getValue("INPUT_SPECTRA"),
				m_inputSpectra)) {
			return false;
		}
	}
	if (m_inputSpectra->size() == 0) {
		ERROR_MSG("Input spectra size is 0!, did you specify INPUT_SPECTRA or INPUT_SPECTRA_PKLBIN?");
		return false;
	}

	if (m_params.exists("INPUT_ISO_ENV")) {
		DEBUG_MSG("Loading input iso envelope from <" << m_params.getValue("INPUT_ISO_ENV") << "> ...");
		if (!m_inputEnv->LoadModel(m_params.getValue("INPUT_ISO_ENV").c_str())) {
			ERROR_MSG("Could not load " << m_params.getValue("INPUT_ISO_ENV"));
			return false;
		}
	} else {
		ERROR_MSG("Must specify INPUT_ISO_ENV");
		return false;
	}

	return true;
}

bool ExecDeconvoluteMS2::saveInputData(std::vector<std::string> & filenames) {
	//SpecSet m_inputContigs; // the input spectra
	std::string spectraFilename = getName() + "_spectra.pklbin";
	m_params.setValue("INPUT_SPECTRA_PKLBIN", spectraFilename);
	if (!fileExists(spectraFilename)) {
		m_inputSpectra->SaveSpecSet_pklbin(spectraFilename.c_str());
	}

	// No method for saving IsoEnvelope

	// Have to set up the output files also so the params will be correct on reload
	m_params.setValue("OUTPUT_SPECTRA_PKLBIN", getName() + "_spectra_z.pklbin");

	std::string paramFilename = getName() + ".params";
	m_params.writeToFile(paramFilename);

	filenames.push_back(paramFilename); // Parameter file MUST be first in vector
	filenames.push_back(spectraFilename);

	return true;
}

bool ExecDeconvoluteMS2::saveOutputData(void) {
	if (m_params.exists("OUTPUT_SPECTRA")) {
		if (!ExecMergeConvert::saveSpecset(m_params.getValue("OUTPUT_SPECTRA"),
				m_outputSpectra)) {
			return false;
		}
	}
	return true;
}

bool ExecDeconvoluteMS2::loadOutputData(void) {
	if (m_params.exists("OUTPUT_SPECTRA")) {
		if (!ExecMergeConvert::loadSpecset(m_params.getValue("OUTPUT_SPECTRA"),
				m_outputSpectra)) {
			return false;
		}
	}
	return true;
}

std::vector<ExecBase *> const & ExecDeconvoluteMS2::split(int numSplit) {

	m_subModules.resize(0);

	if (numSplit < 2) {
		DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
		return m_subModules;
	}

	if (m_inputSpectra->size() == 0) {
		DEBUG_MSG("Must have at least one spectrum");
		return m_subModules;
	}

	m_subModules.resize(numSplit);

	long num_ops = 0;
	for (int i = 0; i < m_inputSpectra->size(); i++) {
		num_ops += (*m_inputSpectra)[i].size();
	}

	long numOpsPerChild = num_ops / ((long) numSplit);
	int globalIdx = 0, childIdx = 0;
	int startIdx, endIdx;

	DEBUG_MSG("Splitting into " << numSplit << " children");
	for (int i = 0; i < numSplit; i++) {

		num_ops = 0;
		startIdx = globalIdx;
		endIdx = globalIdx;

		while (globalIdx < m_inputSpectra->size() && (num_ops <= numOpsPerChild
				|| i == numSplit - 1)) {

			num_ops += (*m_inputSpectra)[globalIdx].size();

			++globalIdx;
		}
		endIdx = globalIdx - 1;

		ParameterList childParams(m_params);
		childParams.setValue("IDX_START", parseInt(startIdx));
		childParams.setValue("IDX_END", parseInt(endIdx));

		ExecBase * theClone = new ExecDeconvoluteMS2(childParams,
				m_inputSpectra, m_inputEnv);

		theClone ->setName(makeName(m_name, i));
		m_subModules[i] = theClone;
	}

	DEBUG_MSG("Splitting success");
	DEBUG_TRACE;
	return m_subModules;
}

bool ExecDeconvoluteMS2::merge(void) {
	if (m_subModules.size() == 0) {
		DEBUG_MSG("No children found when merging");
		return false;
	}

	int num_specs = m_inputSpectra->size();

	m_outputSpectra->resize(num_specs);
	int specIdx = 0;

	DEBUG_MSG("Merging " << m_subModules.size() << " children");
	for (int child = 0; child < m_subModules.size(); child++) {
		ExecDeconvoluteMS2* theChild =
				(ExecDeconvoluteMS2*) m_subModules[child];
		SpecSet* childSpecs = theChild->m_outputSpectra;

		for (int i = 0; i < childSpecs->size(); i++) {
			(*m_outputSpectra)[specIdx] = (*childSpecs)[i];
			++specIdx;
		}
	}
	m_outputSpectra->resize(specIdx);

	DEBUG_MSG("Merging success");
	return true;
}

bool ExecDeconvoluteMS2::validateParams(std::string & error) {
	VALIDATE_PARAM_EXIST("INPUT_SPECTRA");
	VALIDATE_PARAM_EXIST("OUTPUT_SPECTRA");
	return true;
}
}

