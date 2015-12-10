/**
 @file mzrange.h
 */

#ifndef MZRANGE_H
#define MZRANGE_H

#include <vector>
#include <list>
#include <map>
#include <iostream>

using namespace std;

namespace specnets
{

  //local ppm peak tolerance = mass * ppm * PPM_FACTOR
  extern const float PPM_FACTOR;

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
     * @param in_mass mass
     * @param PPM_tol PPM tolerance
     * @return radius of mass error
     */
    static float GetPPMRadius(float in_mass, float PPM_tol);

    /**
     * Computes the weighted average of a sequence of MZRanges. Weight of each
     *   MZRange is the intensity.
     * @param ranges sequence of MZRanges.
     * @return a MZRange with weighted average center, weighted average tolerance,
     *   and summed intensity over all in ranges
     */
    static MZRange WeightedAverage(list<MZRange>* ranges);

    /**
     * Inserts a mzrange into existing mzrange_bins
     * @param mzrange_bins each coalesced mzrange mapped to a sequence of mzranges it is derived from
     * @param new_mzrange mzrange to insert
     * @param mergeType same as in MergeMZRanges
     * @param ranges_to_insert specifies a bin of MZRanges to associate with the merged. Used during
     *   recursive calls.
     * @return coalesced MZRange of overlapping ranges
     */
    //static MZRange InsertMZRange(map<MZRange, list<
    //		MZRange>*>* mzrange_bins, MZRange& new_mzrange, short mergeType,
    //		list<MZRange>* ranges_to_insert = 0);

    //NON-STATIC CLASS METHODS

    /**
     * Default constructor. Sets center and tolerance to 0 without PPM
     */
    MZRange()
    {
      set(0, 0, 0);
    }

    /**
     * MZRange constructor
     * @param in_center the center of this range
     * @param in_tolerance the +/- tolerance of the center
     * @param in_applyPPM if true, the tolerance is taken as the PPM tolerance
     *   (radius = center * tolerance * PPM_FACTOR). If false, radius == tolerance
     */
    MZRange(const float in_mass,
            const float in_intensity,
            const float in_tolerance)
    {
      set(in_mass, in_intensity, in_tolerance);
    }

    /**
     * MZRange copy constructor
     * @param other MZRange instance to copy
     */
    MZRange(const MZRange& other)
    {
      set(other.getMass(), other.getIntensity(), other.getTolerance());
    }

    /**
     * Assignment operator
     */
    MZRange& operator=(const MZRange& other)
    {
      if (this == &other) {
        return *this;
      }
      set(other.getMass(), other.getIntensity(), other.getTolerance());
      return *this;
    }

    /**
     * Less than operator.
     */
    bool operator<(const MZRange& other) const
    {
      return upperBound < other.getLowerBound();
    }

    /**
     * Less than operator for float.
     */
    bool operator<(const float other) const
    {
      return upperBound < other;
    }

    /**
     * Greater than operator.
     */
    bool operator>(const MZRange& other) const
    {
      return lowerBound > other.getUpperBound();
    }

    /**
     * Greater than operator for float.
     */
    bool operator>(float other) const
    {
      return lowerBound > other;
    }

    /**
     * Equals operator. 2 ranges are equal if they overlap.
     */
    bool operator==(const MZRange& other) const
    {
      float other_center = other.getMass();
      float other_tol = other.getTolerance();

      return (EqualWithinRange(getLowerBound(), other_center, other_tol)
          || EqualWithinRange(getUpperBound(), other_center, other_tol)
          || EqualWithinRange(other.getLowerBound(), mass, tolerance)
          || EqualWithinRange(other.getUpperBound(), mass, tolerance));
    }

    /**
     * Not equals operator for float
     */
    bool operator!=(float other) const
    {
      return !((*this) == other);
    }

    /**
     * Not equals operator. 2 ranges are equal if they overlap.
     */
    bool operator!=(const MZRange& other) const
    {
      return !((*this) == other);
    }

    /**
     * Equals operator for float
     */
    bool operator==(float other) const
    {
      return EqualWithinRange(mass, other, tolerance);
    }

    /**
     * <= operator
     */
    bool operator<=(const MZRange& other) const
    {
      return (*this < other) || (*this == other);
    }

    /**
     * <= operator for floats
     */
    bool operator<=(const float other) const
    {
      return (*this < other) || (*this == other);
    }

    /**
     * >= operator
     */
    bool operator>=(const MZRange& other) const
    {
      return (*this > other) || (*this == other);
    }

    /**
     * >= operator for floats
     */
    bool operator>=(const float other) const
    {
      return (*this > other) || (*this == other);
    }

    /**
     * += operator
     */
    MZRange& operator+=(const MZRange& other)
    {
      set(getMass() + other.getMass(),
          getIntensity() + other.getIntensity(),
          getTolerance() + other.getTolerance());
      return *this;
    }

    /**
     * += operator for floats
     */
    MZRange& operator+=(const float other)
    {
      setMass(getMass() + other);
      return *this;
    }

    /**
     * -= operator
     */
    MZRange& operator-=(const MZRange& other)
    {
      set(getMass() - other.getMass(),
          getIntensity() + other.getIntensity(),
          getTolerance() + other.getTolerance());
      return *this;
    }

    /**
     * -= operator floats
     */
    MZRange& operator-=(const float other)
    {
      setMass(getMass() - other);
      return *this;
    }

    /**
     * *= operator floats
     */
    MZRange& operator*=(const float other)
    {
      setMass(getMass() * other);
      return *this;
    }

    /**
     * /= operator floats
     */
    MZRange& operator/=(const float other)
    {
      setMass(getMass() / other);
      return *this;
    }

    /**
     * Addition operator
     */
    MZRange operator+(const MZRange& other) const
    {
      return MZRange(*this) += other;
    }

    /**
     * Addition operator for floats
     */
    MZRange operator+(const float other) const
    {
      return MZRange(*this) += other;
    }

    /**
     * Subtraction operator
     */
    MZRange operator-(const MZRange& other) const
    {
      return MZRange(*this) -= other;
    }

    /**
     * Subtraction operator for floats
     */
    MZRange operator-(const float other) const
    {
      return MZRange(*this) -= other;
    }

    /**
     * * operator for floats
     */
    MZRange operator*(const float other) const
    {
      return MZRange(*this) *= other;
    }

    /**
     * / operator for floats
     */
    MZRange operator/(const float other) const
    {
      return MZRange(*this) /= other;
    }

    /**
     * Resets all fields of this MZRange. Called by constructors.
     * @param in_center
     * @param in_tolerance
     * @param in_applyPPM
     */
    void set(const float in_mass,
             const float in_intensity,
             const float in_tolerance)
    {
      mass = in_mass;
      intensity = in_intensity;
      tolerance = in_tolerance;
      computeRadius();
    }

    void set(MZRange& other)
    {
      mass = other.getMass();
      intensity = other.getIntensity();
      tolerance = other.getTolerance();
      computeRadius();
    }

    void setMass(const float _mass)
    {
      mass = _mass;
      computeRadius();
    }

    void setIntensity(const float _intensity)
    {
      intensity = _intensity;
    }

    /**
     * Sets the tolerance of this range, recomputes radius.
     * @param in_tolerance
     */
    void setTolerance(const float in_tolerance)
    {
      tolerance = in_tolerance;
      computeRadius();
    }

    float getMass() const
    {
      return mass;
    }

    float getIntensity() const
    {
      return intensity;
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
     * Sets the center and radius of this MZRange to the center MZRange of
     *   a sequence of MZRanges.
     * @param peaks sequence of MZRanges
     * @param mergeType specifies how to compute center MZRange
     *   0: MZRange with lowest tolerance is taken as center
     *   1: Take weighted average of centers and tolerances (intensity is the "weight")
     *   2: MZRange with the highest intensity is taken as the center
     *   3: the first MZRange in the list is taken as the center
     * @return
     */
    void MergeMZRanges(list<MZRange>* peaks, short mergeType);

    void printToStream(ostream &output) const
    {
      output << "mass=" << getMass() << ", intensity=" << getIntensity()
          << ", tolerance=" << getTolerance() << endl;
    }

  protected:
    //private data fields, these are updated automatically

    //1-D center of the range, equidistant from upper and lower bounds
    float mass;

    float intensity;

    //Tolerance applied in classic form or in PPM
    float tolerance;

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
      lowerBound = mass - tolerance;
      upperBound = mass + tolerance;
    }
  };

}

#endif // MZRANGE_H
