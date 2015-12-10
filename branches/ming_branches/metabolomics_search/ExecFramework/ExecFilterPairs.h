#ifndef __ExecFilterPairs_H__
#define __ExecFilterPairs_H__

// Module Includes
#include "ExecBase.h"
#include "SpectrumPairSet.h"

// External Includes
#include "spectrum.h"

// System Includes
#include <string>

namespace specnets
{
  /*! \brief Execution class for filtering pairs of spectra

   */
  class ExecFilterPairs : public ExecBase
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief The exemplar constructor.

     Generally this constructor should not be used. It is only used by the
     module execution factory in order to create an exemplar object (without
     valid parameters) which is then used to create a real (valid) object
     using the clone() method.
     @sa clone()
     */
    ExecFilterPairs(void);

    /*! \brief The default constructor

     This is the default constructor. A valid set of parameters must be
     supplied in the input parameter. The parameters can then be verified
     using the validateParams() method.
     @sa validateParams()
     @param ParameterList structure containing all input parameters necessary for execution
     */
    ExecFilterPairs(const ParameterList & inputParams);

    /*! \brief Constructor which takes both the parameters and the input data

     Constructor which takes both the parameters and the input data structure
     containing the spectra to be filtered.
     The parameters can be verified using the validateParams() method.
     @sa validateParams()
     @param inputSpectra Set of input spectra to be aligned
     */
    ExecFilterPairs(const ParameterList & inputParams,
                    SpecSet * inputSpectra);

    /*! \brief Constructor which takes both the parameters and the data structures

     Constructor which takes both the parameters and the data structures necessary for
     executing the alignment. When using this constructor no external data need be read in.
     The parameters can be verified using the validateParams() method.
     @sa validateParams()
     @param inputSpectra Set of input spectra to be aligned
     @param inputSpectraMS2 Set of input spectra to be aligned (non-scored MS/MS spectra, used for cosine matching)
     @param filteredPairs Output set of filtered spectrum pairs
     @param ratios
     @param means
     @param varTerms
     @param alignStats
     @param specStats
     */
    ExecFilterPairs(const ParameterList & inputParams,
                    SpecSet * inputSpectra,
                    SpecSet * inputSpectraMS2,
                    SpectrumPairSet * filteredPairs,
                    std::vector<TwoValues<float> > * ratios,
                    std::vector<TwoValues<float> > * means,
                    std::vector<float> * varTerms,
                    std::list<std::vector<float> > * alignStats,
                    std::vector<std::vector<float> > * specStats,
                    std::vector<unsigned int> * idxKept,
                    std::vector<TwoValues<float> > * pvalues);

    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ExecFilterPairs(void);
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Creates a new module of the virtual class with the given params.

     The parameters should be sufficient for the derived class to be invoked properly.

     @return A pointer to the newly created object
     */
    virtual ExecBase * clone(const ParameterList & inputParams) const;

    //@}

    //! \name MODIFIERS
    //@{
    /*! \brief Executes the module.

     In order to call this method succesfully all the necessary data for
     execution must already be loaded into memory (data members). This can
     be accomplished using the loadInputData() method.

     @return True if execution finished successfully, false otherwise.
     @sa loadInputData()
     */
    virtual bool invoke(void);

    /*! \brief Loads the input data from the files specified in the params.

     Loads all the data from files into the data members. This method is
     primarily used by the execution module to load necessary data when
     executing in a separate process..

     @return True if data was loaded successfully, false otherwise.
     @sa ExecBase(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadInputData(void);

    /*! \brief Saves all the result data to files specified in the params.

     Saves all the result data generated by the module to files specified
     by values in the params. This is used to either save the data permanantly,
     or to be loaded back in after remote execution is finished and results
     need to be merged by the merge() method.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecBase(const ParameterList & input_params), loadInputData(), merge()
     */
    virtual bool saveOutputData(void);

    /*! \brief Saves all the internal data into files specified in the params.

     Saves all the data required for an external process to execute the
     module. The external process would call loadInputData() to reload the
     data into the members before calling invoke().

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecBase(const ParameterList & input_params), loadInputData()
     */
    virtual bool saveInputData(std::vector<std::string> & filenames);

    /*! \brief Loads the output data from the files specified in the params.

     Loads all the data from output files into the data members. The purpose
     of this is to ready "child" modules to be merged back together after
     being executed separately.

     @return True if data was loaded successfully, false otherwise.
     @sa ExecBase(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadOutputData(void);

    /*! \brief Splits the module into multiple "children" for parallel execution.

     Divides the work required by the module into a vector of sub-modules
     that can be executed in parallel (by the the ParallelExecution() class.

     @param numSplit Number of separate modules to split into
     @return The set of sub-modules (children) that the original module has been split into.
     An empty vector implies an error.
     @sa merge()
     */
    virtual std::vector<ExecBase *> const & split(int numSplit);

    /*! \brief Merges the child modules back into a complete result

     This method is only called when split is used for parallel execution
     and after the execution of all children has been completed. This method
     will merge the results generated by each of the child modules into
     one cohesive result as if the module had been run as a single entity.

     @return True if merge could be performed successfully, false otherwise.
     @sa split(int numSplit)
     */
    virtual bool merge(void);

    /*! \brief Performs validation of the input parameters.

     Checks the parameters structure provided in the constructor to see if they
     are sufficient and correct to invoke the module. Also sets the internal
     validity flag so that isValid() will return the correct result.

     @param error A description of the error (if any occurs)
     @return True if the parameters for the module are valid, false otherwise.
     @sa isValid()
     */
    virtual bool validateParams(std::string & error);
    //@}

  private:
    SpecSet * m_inputSpectra; //! The set of input spectra to be aligned
    SpecSet * m_inputSpectraMS2; //! The set of input spectra to be aligned (non-scored MS/MS spectra, used for cosine matching)
    bool ownInput; //! Does this object "own" the input data pointers (and hence have to free them)

    SpectrumPairSet * m_filteredPairs;

    std::vector<TwoValues<float> > * m_ratios;
    std::vector<TwoValues<float> > * m_means;
    std::vector<float> * m_varTerms;
    list<vector<float> > * m_alignStats;
    std::vector<vector<float> > * m_specStats;

    bool ownOutput; //! Does this object "own" the output data pointers (and hence have to free them)
  };

} // namespace specnets

#endif // __ExecFilterPairs_H__
