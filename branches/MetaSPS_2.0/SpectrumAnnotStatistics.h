/**
 * Annotation of spectra.
 */

#ifndef SPECTRUM_ANNOT_STATISTICS_H_
#define SPECTRUM_ANNOT_STATISTICS_H_

#include <string>
#include <cmath>
#include <set>

#include "utils.h"
#include "aminoacid.h"
#include "PeptideSpectrumMatchSet.h"
#include "spectrum.h"

namespace specnets
{
  class SpectrumAnnotStatistics
  {
  public:

    /**
     * Calculate explained intensity for input spectrum
     * @param psm - PeptideSpectrumMatch with pointer to spectrum (@see PeptideSpectrumMatch)
     * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
     * information on ftIonFragment names)
     */
    static float percentExplainedIntensity(PeptideSpectrumMatch &psm, string &ionNames);

    /**
	 * Calculate explained intensity for input psmSet
	 * @param psmSet - PeptideSpectrumMatchSet with pointer to spectr (@see PeptideSpectrumMatchSet)
	 * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
	 * information on ftIonFragment names)
	 */
	static float percentExplainedIntensity(PeptideSpectrumMatchSet &psmSet,
			string &ionNames);

    /**
     * Calculate percentage of explained peaks for input spectrum
     * @param psm - PeptideSpectrumMatch with pointer to spectrum (@see PeptideSpectrumMatch)
     * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
     * information on ftIonFragment names)
     */
    static float percentExplainedPeaks(PeptideSpectrumMatch &psm, string &ionNames);

    /**
	 * Calculate percentage of explained peaks for input psmSet
	 * @param psmSet - PeptideSpectrumMatchSet with pointer to spectr (@see PeptideSpectrumMatchSet)
	 * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
	 * information on ftIonFragment names)
	 */
	static float percentExplainedPeaks(PeptideSpectrumMatchSet &psmSet,
			string &ionNames);

    /**
     * Calculate percentage of observed ions vs. possible observed ions for input spectrum (N-1 per ion type)
     * @param psm - PeptideSpectrumMatch with pointer to spectrum (@see PeptideSpectrumMatch)
     * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
     * information on ftIonFragment names)
     */
    static float percentObservedIons(PeptideSpectrumMatch &psm, string &ionNames);

    /**
	 * Calculate percentage of observed ions vs. possible observed ions for input psmSet
	 * @param psmSet - PeptideSpectrumMatchSet with pointer to spectr (@see PeptideSpectrumMatchSet)
	 * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
	 * information on ftIonFragment names)
	 */
	static float percentObservedIons(PeptideSpectrumMatchSet &psmSet,
			string &ionNames);

    /**
     * Count total peaks
     * @param spectrum - spectrum to be considered
     */
    static int totalPeaks(PeptideSpectrumMatch &psm);

    /**
     * PPM precursor mass error, defined as 10^6*
     * (Experimental PrecursorMass-TheoreticalPrecursorMass)/TheoreticalPrecursorMass.
     * Choose the smallest PPM error for TheoreticalPrecursorMass=mono isotopic mass
     * and TheoreticalPrecursorMass=one-C13 mass;
     * Note that ParentMass = Sum(amino acid masses)+mass(H2O)+charge*mass(Hion) and TheoreticalPrecursorMass=ParentMass/charge.
     * @param psm - PeptideSpectrumMatch with pointer to spectrum (@see PeptideSpectrumMatch)
     * @param charge - use this as the charge instead of charge defined on spectrum
     */
    static float parentMassErrorPPM(PeptideSpectrumMatch &psm, int charge = -1);

    /**
     * Parent mass error in Da, defined as Exeprimental_ParentMass-ParentMass,
     * choose the smallest parent mass error for ParentMass=mono isotopic mass and
     * ParentMass=one-C13 mass.
     * @param psm - PeptideSpectrumMatch with pointer to spectrum (@see PeptideSpectrumMatch)
     * @param charge - use this as the charge instead of charge defined on spectrum
     */
    static float parentMassErrorDa(PeptideSpectrumMatch &psm, int charge = -1);

    /**
     * a break is defined as a peptide bond for which we observe at least one b-ion or y-ion with charge 1 or 2. A peptide of length N has N-1 possible breaks,
     *  percentage of observed breaks is defined as number_of_observed_breaks / (N-1).
     * 	@param spectrum - spectrum to be considered
     */
    static float observedBreaks(PeptideSpectrumMatch &psm, string &ionNames);

    /**
	 * Calculate percentage of observed peptide breaks input psmSet
	 * @param psmSet - PeptideSpectrumMatchSet with pointer to spectr (@see PeptideSpectrumMatchSet)
	 * @param ionNames - string indicating which ion names to consider (@see spectrum_scoring for
	 * information on ftIonFragment names)
	 */
	static float observedBreaks(PeptideSpectrumMatchSet &psmSet,
			string &ionNames);
  };
}
#endif /* SPECTRUM_ANNOT_STATISTICS_H_ */
