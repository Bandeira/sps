/**
 * Annotation of spectra.
 */

#ifndef SPECTRUM_ANNOT_STATISTICS_H_
#define SPECTRUM_ANNOT_STATISTICS_H_

#include "spectrum.h"

#include <string>
#include <cmath>
#include <set>

#include "spectrum_scoring.h"
#include "utils.h"
#include "aminoacid.h"


class SpectrumAnnotStatistics {
public:

		/**
		 * Calculate explained intensity for input spectrum
		 * @param spectrum - annotated Spectrum (@see spectrum)
		 * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
		 * information on ftIonFragment names)
		 */
		float percentExplainedIntensity(Spectrum &spectrum, string &ionNames);

		/**
		 * Calculate percentage of explained peaks for input spectrum
		 * @param spectrum - annotated Spectrum (@see spectrum)
		 * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
		 * information on ftIonFragment names)
		 */
		float percentExplainedPeaks(Spectrum &spectrum, string &ionNames);

		/**
		 * Calculate percentage of observed ions vs. possible observed ions for input spectrum (N-1 per ion type)
		 * @param spectrum - annotated Spectrum (@see spectrum)
		 * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
		 * information on ftIonFragment names)
		 */
		float percentObservedIons(Spectrum &spectrum, string &ionNames);

		/**
		 * Count total peaks
		 * @param spectrum - spectrum to be considered
		 */
		int totalPeaks(Spectrum &spectrum);

		/**
		 * PPM precursor mass error, defined as 10^6*
		 * (Experimental PrecursorMass-TheoreticalPrecursorMass)/TheoreticalPrecursorMass.
		 * Choose the smallest PPM error for TheoreticalPrecursorMass=mono isotopic mass
		 * and TheoreticalPrecursorMass=one-C13 mass;
		 * Note that ParentMass = Sum(amino acid masses)+mass(H2O)+charge*mass(Hion) and TheoreticalPrecursorMass=ParentMass/charge.
		 * @param spectrum - spectrum to be considered.
		 * @param charge - use this as the charge instead of charge defined on spectrum
		 */
		float parentMassErrorPPM(Spectrum &spectrum, int charge = -1);

		/**
		 * Parent mass error in Da, defined as Exeprimental_ParentMass-ParentMass,
		 * choose the smallest parent mass error for ParentMass=mono isotopic mass and
		 * ParentMass=one-C13 mass.
		 * @param spectrum - spectrum to be considered
		 * @param jumps - amino acid masses to use for getting peptide mass (@see aminoacid.h)
		 * @param charge - use this as the charge instead of charge defined on spectrum
		 */
		float parentMassErrorDa(Spectrum &spectrum, int charge = -1);

		/**
		 * a break is defined as a peptide bond for which we observe at least one b-ion or y-ion with charge 1 or 2. A peptide of length N has N-1 possible breaks,
		 *  percentage of observed breaks is defined as number_of_observed_breaks / (N-1).
		 * 	@param spectrum - spectrum to be considered
		 */
		float observedBreaks(Spectrum &spectrum, string &ionNames);
};

#endif /* SPECTRUM_ANNOT_STATISTICS_H_ */
