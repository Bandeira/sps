/**
 @file mzrange.h
 */

#ifndef MZRANGE_H
#define MZRANGE_H

#include <vector>

#include "twovalues.h"
#include "utils.h"

using namespace std;

/**
 * 1-dimensional range for use with mass spectrum m/z values.
 * Supports mass tolerance in both classic and PPM form.
 *
 * Handy for storing (peak masses +/- tolerance) in a set, map, or list.
 * This could also be handy if specnets ever utilizes PPM tolerance, as
 * the tolerance of every contig peak would be different.
 */
class MZRange
{
public:

  //STATIC CLASS METHODS

  /**
   * Checks if two 1D values are equal within a range
   * @param center1
   * @param center2
   * @param range
   * @return true if both centers are within radius of each other, false otherwise
   */
  static bool EqualWithinRange(float center1, float center2, float radius);

  /**
   * Computes the radius of a mass with PPM tolerance
   * @param in_center mass
   * @param PPM_tol PPM tolerance
   * @return radius of mass error
   */
  static float GetPPMRadius(float in_center, float PPM_tol);

  /**
   * Computes the weighted average of a sequence of values
   * @param values
   * @param weights all should be >= 0. If one score is < 0, an offset is added to all scores
   * @return weighted average [0], summed score [1]
   */
  static TwoValues<float> WeightedAverage(vector<float>& values,
                                          vector<float>& scores);

  /**
   * Inserts a mzrange into existing mzrange_bins
   * @param mzrange_bins each average mzrange mapped to a sequence of mzranges it is derived from
   * @param new_mzrange mzrange to insert
   * @return average of overlapping mzranges
   */
  static MZRange InsertMZRange(map<MZRange, vector<MZRange> >& mzrange_bins,
                               MZRange& new_mzrange);

  //NON-STATIC CLASS METHODS

  /**
   * Default constructor. Sets center and tolerance to 0 without PPM
   */
  MZRange()
  {
    set(0, 0, false);
  }

  /**
   * MZRange constructor
   * @param in_center the center of this range
   * @param in_tolerance the +/- tolerance of the center
   * @param in_applyPPM if true, the tolerance is taken as the PPM tolerance
   *   (radius = center * tolerance * PPM_FACTOR). If false, radius == tolerance
   */
  MZRange(float in_center, float in_tolerance, bool in_applyPPM = false)
  {
    set(in_center, in_tolerance, in_applyPPM);
  }

  /**
   * MZRange copy constructor
   * @param other MZRange instance to copy
   */
  MZRange(const MZRange& other)
  {
    set(other.getCenter(), other.getTolerance(), other.getPPM());
  }

  /**
   * Assignment operator
   */
  MZRange& operator=(const MZRange& other)
  {
    set(other.getCenter(), other.getTolerance(), other.getPPM());
  }

  /**
   * Assignment operator for float
   */
  MZRange& operator=(const float other)
  {
    setCenter(other);
  }

  /**
   * Less than operator.
   */
  bool operator<(const MZRange& other) const
  {
    return upperBound < other.getLowerBound();
  }

  /**
   * Greater than operator.
   */
  bool operator>(const MZRange& other) const
  {
    return lowerBound > other.getUpperBound();
  }

  /**
   * Equals operator. 2 ranges are equal if they overlap.
   */
  bool operator==(const MZRange& other) const
  {
    float other_center = other.getCenter();
    float other_rad = other.getRadius();

    return (EqualWithinRange(getLowerBound(), other_center, other_rad)
        || EqualWithinRange(getUpperBound(), other_center, other_rad)
        || EqualWithinRange(other.getLowerBound(), center, radius)
        || EqualWithinRange(other.getUpperBound(), center, radius));
  }

  /**
   * Resets all fields of this MZRange. Called by constructors.
   * @param in_center
   * @param in_tolerance
   * @param in_applyPPM
   */
  void set(float in_center, float in_tolerance, bool in_applyPPM = false)
  {
    center = in_center;
    tolerance = in_tolerance;
    applyPPM = in_applyPPM;
    computeRadius();
  }

  /**
   * Sets the center if this range.
   * @param in_center
   */
  void setCenter(float in_center)
  {
    center = in_center;
    computeRadius();
  }

  /**
   * Sets the tolerance of this range, recomputes radius.
   * @param in_tolerance
   */
  void setTolerance(float in_tolerance)
  {
    tolerance = in_tolerance;
    computeRadius();
  }

  /**
   * Sets the PPM application of this MZRange
   * @param in_applyPPM
   */
  void setPPM(bool in_applyPPM)
  {
    applyPPM = in_applyPPM;
    computeRadius();
  }

  /**
   * Returns the center of this MZRange
   * @return center
   */
  float getCenter() const
  {
    return center;
  }

  /**
   * Returns the tolerance of this MZRange
   * @return tolerance
   */
  float getTolerance() const
  {
    return tolerance;
  }

  /**
   * Returns the radius of this range. If applyPPM is false,
   *   this is equivalent to the tolerance
   * @return radius
   */
  float getRadius() const
  {
    return radius;
  }

  /**
   * Returns the lower bound of this range
   * @return lowerBound
   */
  float getLowerBound() const
  {
    return lowerBound;
  }

  /**
   * Returns the upper bound of this range
   * @return upperBound
   */
  float getUpperBound() const
  {
    return upperBound;
  }

  /**
   * Returns whether or not PPM is applied to this range
   * @return applyPPM
   */
  float getPPM() const
  {
    return applyPPM;
  }

  /**
   * Sets the center and radius of this MZRange to the weighted average of
   *   a sequence of MZRanges. applyPPM is set to false since the new
   *   weighted tolerance will equal the radius
   * @param peaks sequence of MZRanges
   * @param scores sequence of scores parallel to peaks. If not specified,
   *   all centers are scored equally.
   * @return summed score of the new center. This is the size of centers if
   *   scores are not specified.
   */
  float MergeMZRanges(vector<MZRange>& peaks, vector<float>* scores = 0);

private:
  //private data fields, these are updated automatically

  //1-D center of the range, equidistant from upper and lower bounds
  float center;

  //Tolerance applied in classic form or in PPM
  float tolerance;

  //Whether or not to apply tolerance in PPM
  bool applyPPM;

  //Applied radius of range
  float radius;

  //center - radius. On hand for faster comparison operations
  float lowerBound;

  //center + radius. On hand for faster comparison operations
  float upperBound;

  /*
   * Called internally to recompute radius and lower/upper bound
   *   if values change
   */
  void computeRadius()
  {
    if (applyPPM)
      radius = GetPPMRadius(center, tolerance);
    else
      radius = tolerance;

    lowerBound = center - radius;
    upperBound = center + radius;
  }
};

#endif // MZRANGE_H
