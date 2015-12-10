#include "mzrange.h"

/**
 * Checks if two 1D values are equal within a range
 * @param center1
 * @param center2
 * @param range
 * @return true if both centers are within radius of each other, false otherwise
 */
bool MZRange::EqualWithinRange(float center1, float center2, float radius)
{
  return (center1 >= center2 - radius) && (center1 <= center2 + radius);
}

/**
 * Computes the radius of a mass with PPM tolerance
 * @param in_center mass
 * @param PPM_tol PPM tolerance
 * @return radius of mass error
 */
float MZRange::GetPPMRadius(float in_center, float PPM_tol)
{
  return in_center * PPM_tol * PPM_FACTOR;
}

/**
 * Computes the weighted average of a sequence of values
 * @param values
 * @param weights all should be >= 0. If one score is < 0, an offset is added to all scores
 * @return weighted average [0], summed score [1]
 */
TwoValues<float> MZRange::WeightedAverage(vector<float>& values,
                                          vector<float>& scores)
{
  if (values.size() != scores.size())
  {
    cerr << "ERROR: Values (size=" << values.size() << ") and scores (size="
        << scores.size() << ") vectors must be same size\n";
    return TwoValues<float> (0, 0);
  }

  float total_score = 0;
  float new_average = 0;
  float min_score = 1.0;
  float offset = 0;
  for (int i = 0; i < scores.size(); i++)
  {
    total_score += scores[i];
    min_score = min(min_score, scores[i]);
  }

  //check for negative scores
  if (min_score < 0)
  {
    //make all scores positive
    offset = 0 - min_score;
    total_score += offset * ((float)scores.size());
  }

  //add up weighted average
  for (int i = 0; i < values.size(); i++)
  {
    new_average += values[i] * ((offset + scores[i]) / total_score);
  }

  return TwoValues<float> (new_average, total_score);
}

/**
 * Inserts a mzrange into existing mzrange_bins
 * @param mzrange_bins each average mzrange mapped to a sequence of mzranges it is derived from
 * @param new_mzrange mzrange to insert
 * @return average of overlapping mzranges
 */
MZRange MZRange::InsertMZRange(map<MZRange, vector<MZRange> >& mzrange_bins,
                               MZRange& new_mzrange)
{
  vector<MZRange> overlapping;

  if (mzrange_bins.count(new_mzrange) > 0)
  {
    //just add the mzrange if it's unique
    overlapping.push_back(new_mzrange);
    mzrange_bins[new_mzrange] = overlapping;
    return MZRange(new_mzrange);
  }

  //get all overlapping mzrange for this range
  overlapping = mzrange_bins[new_mzrange];
  overlapping.push_back(new_mzrange);

  mzrange_bins.erase(new_mzrange);

  //compute average mzrange
  MZRange average_mzrange;
  average_mzrange.MergeMZRanges(overlapping);
  mzrange_bins[average_mzrange] = overlapping;

  return MZRange(average_mzrange);
}

/**
 * Sets the center and radius of this MZRange to the weighted average of
 *   a sequence of MZRanges. applyPPM is set to false since the new
 *   weighted tolerance will equal the radius
 * @param peaks sequence of MZRanges
 * @param scores sequence of scores parallel to peaks
 * @return summed score of the new center. This is the size of centers if
 *   scores are not specified.
 */
float MZRange::MergeMZRanges(vector<MZRange>& peaks, vector<float>* scores)
{
  vector<float> use_scores(peaks.size(), 1.0);
  if (scores != 0)
  {
    use_scores = *scores;
  }

  vector<float> centers(peaks.size());
  vector<float> radii(peaks.size());

  for (int i = 0; i < peaks.size(); i++)
  {
    centers[i] = peaks[i].getCenter();
    radii[i] = peaks[i].getRadius();
  }

  TwoValues<float> weightedCenters;
  TwoValues<float> weightedRadii;

  weightedCenters = WeightedAverage(centers, use_scores);
  weightedRadii = WeightedAverage(radii, use_scores);

  set(weightedCenters[0], weightedRadii[0], false);

  return weightedCenters[1];
}
