/**
 * Annotation of ion statistics
 */

/**
 * Calculate explained intensity for input spectrum
 * @param spectrum - annotated Spectrum (@see spectrum)
 * @param ionNames - comma delimited string indicating which ion names to consider (@see spectrum_scoring for
 * information on ftIonFragment names)
 * @param jumps - amino acid masses to use. @see aminoacid.h.
 */
#include "spectrum_annot_statistics.h"

using namespace std;

float SpectrumAnnotStatistics::percentExplainedIntensity(Spectrum &spectrum,
                                                         string &ionNames)
{
  if (spectrum.annotation_peptide.empty())
  {
    cerr << "Warning: Spectrum " << spectrum.scan << " is not annotated!"
        << endl;
    return -1.0;
  }

  AAJumps jumps(1);

  vector < string > ion_types;

  splitText(ionNames.c_str(), ion_types, (const char*)",");

  vector<float> masses;

  double peptide_mass = jumps.getPeptideMass(spectrum.annotation_peptide);

  float total_intensity = 0;
  float identified_intensity = 0;

  for (unsigned i = 0; i < spectrum.peakList.size(); i++)
  {
    if (spectrum.peakList[i][0] >= AAJumps::minAAmass
        && spectrum.peakList[i][0] <= (peptide_mass - 37))
    { // ignore small masses or masses larger than we expect.
      total_intensity += spectrum.peakList[i][1]; // add to total intensity regardless of annotation
      const ftIonFragment* curr_frag = spectrum.annotation[i].first;
      if (curr_frag != NULL)
      { // make sure peak is annotated
        if (ionNames.compare("all") == 0)
        { // include all identified peaks
          identified_intensity += spectrum.peakList[i][1];
        }
        else
        {
          for (int j = 0; j < ion_types.size(); j++)
          {
            if (ion_types[j].compare(curr_frag->name) == 0)
            { // include only listed peaks
              identified_intensity += spectrum.peakList[i][1];
            }
          }
        }
      }
    }
  }
  return (identified_intensity / total_intensity) * 100;
}

/**
 * Calculate percentage of explained peaks for input spectrum
 * @param spectrum - annotated Spectrum (@see spectrum)
 * @param ionNames - comma delimited string  indicating which ion names to consider (@see spectrum_scoring for
 * information on ftIonFragment names)
 */
float SpectrumAnnotStatistics::percentExplainedPeaks(Spectrum &spectrum,
                                                     string &ionNames)
{
  if (spectrum.annotation_peptide.empty())
  {
    cerr << "Warning: Spectrum " << spectrum.scan << " is not annotated!"
        << endl;
    return -1.0;
  }
  AAJumps jumps(1);

  vector < string > ion_types;

  splitText(ionNames.c_str(), ion_types, (const char*)",");

  double peptide_mass = jumps.getPeptideMass(spectrum.annotation_peptide);

  float total_peaks = 0;
  float identified_peaks = 0;

  for (unsigned i = 0; i < spectrum.peakList.size(); i++)
  {
    if (spectrum.peakList[i][0] >= AAJumps::minAAmass
        && spectrum.peakList[i][0] <= (peptide_mass - 37))
    { // ignore small masses or masses larger than we expect.
      total_peaks++;// add to total regardless of annotation
      const ftIonFragment* curr_frag = spectrum.annotation[i].first;
      if (curr_frag != NULL)
      { // make sure peak is annotated
        if (ionNames.compare("all") == 0)
        { // include all identified peaks
          identified_peaks++;
        }
        else
        {
          for (int j = 0; j < ion_types.size(); j++)
          {
            if (ion_types[j].compare(curr_frag->name) == 0)
            { // include only listed peaks
              identified_peaks++;
            }
          }
        }
      }
    }
  }
  return ((float)identified_peaks / (float)total_peaks) * 100;
}

/**
 * Calculate percentage of observed ions vs. possible observed ions for input spectrum (N-1 per ion type)
 * @param spectrum - annotated Spectrum (@see spectrum)
 * @param ionNames - comma delimited string  indicating which ion names to consider (@see spectrum_scoring for
 * information on ftIonFragment names)
 */
float SpectrumAnnotStatistics::percentObservedIons(Spectrum &spectrum,
                                                   string &ionNames)
{
  if (spectrum.annotation_peptide.empty())
  {
    cerr << "Warning: Spectrum " << spectrum.scan << " is not annotated!"
        << endl;
    return -1.0;
  }
  AAJumps jumps(1);

  vector < string > ion_types;

  splitText(ionNames.c_str(), ion_types, (const char*)",");

  double peptide_mass = jumps.getPeptideMass(spectrum.annotation_peptide);

  int peptide_length = jumps.getPeptideLength(spectrum.annotation_peptide);

  int num_ions; //number of different ion types
  int identified_ions = 0; //number of identified ion peaks

  if (ionNames.compare("all") == 0)
  {
    num_ions = spectrum.ionTypes.size();
  }
  else
  {
    num_ions = ion_types.size();
  }

  for (unsigned i = 0; i < spectrum.peakList.size(); i++)
  {
    if (spectrum.peakList[i][0] >= AAJumps::minAAmass
        && spectrum.peakList[i][0] <= (peptide_mass - 37))
    { // ignore small masses or masses larger than we expect.
      const ftIonFragment* curr_frag = spectrum.annotation[i].first;
      if (curr_frag != NULL)
      { // make sure peak is annotated
        if (ionNames.compare("all") == 0)
        { // include all identified peaks
          identified_ions++;
        }
        else
        {
          for (int j = 0; j < ion_types.size(); j++)
          {
            if (ion_types[j].compare(curr_frag->name) == 0)
            { // include only listed peaks
              identified_ions++;
            }
          }
        }
      }
    }
  }
  return ((float)identified_ions / ((float)num_ions * ((float)peptide_length
      - 1))) * 100;
}

/**
 * Count total peaks in range we're looking at
 * @param spectrum - spectrum to be considered
 */
int SpectrumAnnotStatistics::totalPeaks(Spectrum &spectrum)
{
  if (spectrum.annotation_peptide.empty())
  {
    cerr << "Warning: Spectrum " << spectrum.scan << " is not annotated!"
        << endl;
    return -1;
  }

  int total_peaks = 0;

  AAJumps jumps(1);

  double peptide_mass = jumps.getPeptideMass(spectrum.annotation_peptide);

  for (unsigned i = 0; i < spectrum.peakList.size(); i++)
  {
    if (spectrum.peakList[i][0] >= AAJumps::minAAmass
        && spectrum.peakList[i][0] <= (peptide_mass - 37))
    { // ignore small masses or masses larger than we expect.
      total_peaks++;// add to total regardless of annotation
    }
  }
  return total_peaks;
}

/**
 * PPM precursor mass error, defined as 10^6*
 * (Experimental PrecursorMass-TheoreticalPrecursorMass)/TheoreticalPrecursorMass.
 * Choose the smallest PPM error for TheoreticalPrecursorMass=mono isotopic mass
 * and TheoreticalPrecursorMass=one-C13 mass;
 * Note that ParentMass = Sum(amino acid masses)+mass(H2O)+charge*mass(Hion) and TheoreticalPrecursorMass=ParentMass/charge.
 * @param spectrum - spectrum to be considered.
 */
float SpectrumAnnotStatistics::parentMassErrorPPM(Spectrum &spectrum,
                                                  int charge)
{
  if (spectrum.annotation_peptide.empty())
  {
    cerr << "Warning: Spectrum " << spectrum.scan << " is not annotated!"
        << endl;
    return -1.0;
  }

  AAJumps jumps(1);

  if (charge <= 0)
  { //charge is not defined
    charge = spectrum.parentCharge;
  }

  double peptide_mass = jumps.getPeptideMass(spectrum.annotation_peptide);
  double monoisotopic_mass = peptide_mass;
  double oneisotopic_mass = peptide_mass + AAJumps::massHion; /* one C13 value */

  /* set to parent mass/charge */
  monoisotopic_mass = (monoisotopic_mass + ((double)charge * AAJumps::massHion)
      + AAJumps::massH2O) / (double)charge;
  oneisotopic_mass = (oneisotopic_mass + ((double)charge * AAJumps::massHion)
      + AAJumps::massH2O) / (double)charge;

#ifdef DEBUG
  std::cout << "monoisotopic_mass " << monoisotopic_mass << endl;
  std::cout << "oneisotopic_mass " << oneisotopic_mass << endl;
  std::cout << "charge " << charge << endl;
#endif

#ifdef DEBUG
  std::cout << "ParentMZ " << spectrum.parentMZ << endl;
#endif

  /* get PPM */
  float monoisotopic_ppm = 1000000.0 * ((spectrum.parentMZ - monoisotopic_mass)
      / monoisotopic_mass);
  float oneisotopic_ppm = 1000000.0 * ((spectrum.parentMZ - oneisotopic_mass)
      / monoisotopic_mass);

#ifdef DEBUG
  std::cout << "monoisotopic_ppm " << monoisotopic_ppm << endl;
  std::cout << "oneisotopic_ppm " << oneisotopic_ppm << endl;
#endif

  if (abs(monoisotopic_ppm) > abs(oneisotopic_ppm))
  {
    return oneisotopic_ppm;
  }
  else
  {
    return monoisotopic_ppm;
  }
}

/**
 * Parent mass error in Da, defined as Exeprimental_ParentMass-ParentMass,
 * choose the smallest parent mass error for ParentMass=mono isotopic mass and
 * ParentMass=one-C13 mass.
 * @param spectrum - spectrum to be considered
 * @param charge - use this instead of defined charge on spectrum
 */
float SpectrumAnnotStatistics::parentMassErrorDa(Spectrum &spectrum, int charge)
{
  if (spectrum.annotation_peptide.empty())
  {
    cerr << "Warning: Spectrum " << spectrum.scan << " is not annotated!"
        << endl;
    return -1.0;
  }

  if (charge <= 0)
  { //charge is not defined
    charge = spectrum.parentCharge;
  }

  AAJumps jumps(1);

  double peptide_mass = jumps.getPeptideMass(spectrum.annotation_peptide);
  double monoisotopic_mass = peptide_mass;
  double oneisotopic_mass = peptide_mass + AAJumps::massHion; /* one C13 value */

  /* set to parent mass/charge */
  monoisotopic_mass = (monoisotopic_mass + (charge * AAJumps::massHion)
      + AAJumps::massH2O) / charge;
  oneisotopic_mass = (oneisotopic_mass + (charge * AAJumps::massHion)
      + AAJumps::massH2O) / charge;

  /* get error in Da */
  float monoisotopic_da = spectrum.parentMZ - monoisotopic_mass;
  float oneisotopic_da = spectrum.parentMZ - oneisotopic_mass;

  if (abs(monoisotopic_da) > abs(oneisotopic_da))
  {
    return oneisotopic_da;
  }
  else
  {
    return monoisotopic_da;
  }
}

/**
 * a break is defined as a peptide bond for which we observe at least one b-ion or y-ion with charge 1 or 2. A peptide of length N has N-1 possible breaks,
 *  percentage of observed breaks is defined as number_of_observed_breaks / (N-1).
 * 	@param spectrum - spectrum to be considered.
 *  @param ionNames - comma delimited string  indicating which ion names to consider (@see spectrum_scoring for
 * information on ftIonFragment names)
 */
float SpectrumAnnotStatistics::observedBreaks(Spectrum &spectrum,
                                              string &ionNames)
{
  if (spectrum.annotation_peptide.empty())
  {
    cerr << "Warning: Spectrum " << spectrum.scan << " is not annotated!"
        << endl;
    return -1.0;
  }

  AAJumps jumps(1);

  int peptide_length = jumps.getPeptideLength(spectrum.annotation_peptide);
  float peptide_mass = jumps.getPeptideMass(spectrum.annotation_peptide);

  vector < string > ion_types;

  splitText(ionNames.c_str(), ion_types, (const char*)",");

  set<short> breaks; //keep track of breaks for which we've seen an ion
  set<short>::iterator breaks_it; //breaks iterator

  breaks_it = breaks.begin();

  for (unsigned i = 0; i < spectrum.peakList.size(); i++)
  {
    if (spectrum.peakList[i][0] >= AAJumps::minAAmass
        && spectrum.peakList[i][0] <= (peptide_mass - 37))
    { // ignore small masses or masses larger than we expect.
      const ftIonFragment* curr_frag = spectrum.annotation[i].first;
      if (curr_frag != NULL)
      {
        if (ionNames.compare("all") == 0)
        { // include all identified peaks
          breaks_it = breaks.insert(breaks_it,
                                    spectrum.annotation[i].second);
        }
        else
        {
          for (int j = 0; j < ion_types.size(); j++)
          {
            if (ion_types[j].compare(curr_frag->name) == 0)
            {
              breaks_it = breaks.insert(breaks_it,
                                        spectrum.annotation[i].second);
            }
          }
        }
      }
    }
  }
  return ((float)breaks.size() / ((float)peptide_length - 1)) * 100;
}

