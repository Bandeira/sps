/*
 * SpectrumLibrary.h
 *
 *  Created on: Apr 21, 2011
 *      Author: jsnedecor
 */

#ifndef SPECTRUMLIBRARY_H_
#define SPECTRUMLIBRARY_H_

// Module Includes
#include "PeptideSpectrumMatch.h"
#include "PeptideSpectrumMatchSet.h"
#include "spectrum.h"
#include "Logger.h"

// System Includes
#include <vector>

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#endif

using namespace std;

namespace specnets
{
  /*! \brief Class to build spectrum library

   This class takes a single peptide spectrum match set and
   spectrum set and takes the top scoring PSM for each spectrum.
   */
  class SpectrumLibrary
  {
  public:
    /*! \brief Return associated  SpectrumLibraryMatch for
     that vector position.
     */
    psmPtr & operator[](unsigned int i);

    /*! \brief Return associated  SpectrumLibraryMatch
     for that vector position.
     */
    const psmPtr & operator[](unsigned int i) const;

    /*! \brief Set value of one SpectrumLibraryMatchSet to another

     */
    SpectrumLibrary & operator=(SpectrumLibrary &other);

    /*! \brief Returns size of m_library vector

     */
    unsigned int size();

    /*! \brief Resizes m_library vector

     @param newSize the new size of the parameter vector
     */
    unsigned int resize(unsigned int newSize);

    /*! \brief adds new values to peptide spectrum library.

     @param spectrum candidate spectrum for current psm
     @param psm candidate spectrum psm
     @param score used to compare to existing spectrum library.
     @param scoreAscending indicates whether we want a higher or lower score.
     True = p-value (higher values worse), false = MQScore (higher values better)
     @return true if this spectrum has been added to the library, false if not.
     */
    bool addToLibrary(Spectrum &spectrum,
                      PeptideSpectrumMatch &psm,
                      float score,
                      bool scoreAscending = 0);

    /*! \brief adds new values to peptide spectrum library with key as parameter.

      @param spectrum candidate spectrum for current psm
      @param psm candidate spectrum psm
      @param score used to compare to existing spectrum library.
      @param key indicates the key we wish to store this peptide under
      @param scoreAscending indicates whether we want a higher or lower score.
      True = p-value (higher values worse), false = MQScore (higher values better)
      @return true if this spectrum has been added to the library, false if not.
      */
    bool addToLibrary(Spectrum &spectrum,
                      PeptideSpectrumMatch &psm,
                      float score,
                      string &key,
                      bool scoreAscending = 0);

    /*! \brief returns spectrum library match for annotation and charge
     *
     * Returns matching spectrum if current annotation and charge is defined.
     * Returns false if annotation and charge do not have candidate spectrum.
     * @param annotation - Input peptide annotation
     * @param charge - Input charge specification.
     */
    bool getSpectrumLibraryMatch(string &annotation,
                                 int charge,
                                 psmPtr &outputPsm);

    PeptideSpectrumMatchSet m_library; //! vector to contain top match spectra for library.
  protected:
    /*! \brief build m_map from m_library
     *
     */
    void buildKeyMap(void);

    /*! \brief helper function for addToLibrary
     *
     */
    void addSpectrumMatch(Spectrum &spectrum,
                          PeptideSpectrumMatch &psm,
                          string &key,
                          float score);

    vector<float> m_libraryScores; //! vector of score associated with m_library
    std::tr1::unordered_map<string, unsigned int> m_map; //! associates annotation_charge with the index in m_library.
  };
}

#endif /* SPECTRUMLIBRARY_H_ */
