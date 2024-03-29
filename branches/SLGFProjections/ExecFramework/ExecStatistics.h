#ifndef __ExecStatistics_H__
#define __ExecStatistics_H__

// External Includes
#include "ExecBase.h"
#include "ExecMergeConvert.h"
#include "ParameterList.h"
#include "SpectrumAnnotStatistics.h"
#include "SpectrumAnnotParameterList.h"
#include "spectrum.h"
#include "PeptideSpectrumMatch.h"

// System Includes
#include <string>
#include <vector>

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_set>
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_set>
#  include <unordered_map>
#endif

namespace specnets
{
  /*! \brief Class to generate statistics on a per spectrum basis

   This class is used to generate the statistical input to SVM models
   */
  class ExecStatistics : public ExecBase
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
    ExecStatistics(void);

    /*! \brief The default constructor

     This is the default constructor. A valid set of parameters must be
     supplied in the input parameter. The parameters can then be verified
     using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     */
    ExecStatistics(const ParameterList & inputParams);

    /*! \brief Constructor which takes both the parameters and the data structures

     Constructor which takes both the parameters and the data structures necessary for
     generating statistics. When using this constructor no external data need be read in.
     The parameters can be verified using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     @param spectra Input spectra to be considered
     @param model Fragmentation model for MS2 spectra
     @param statsParams Parameters for which statistics to generate
     @param peptideResults Peptide results from Inspect or Specnets.
     @param spectraStats Output spectral statistics
     @param spectraHeader Output header for spectral statistics
     */
    ExecStatistics(const ParameterList & inputParams,
                   SpecSet * spectra,
                   MS2ScoringModel * model,
                   SpectrumAnnotParameterList * statsParams,
                   PeptideSpectrumMatchSet * peptideResults,
                   vector<vector<float> > * spectraStats,
                   vector<string> * spectraHeader);
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ExecStatistics(void);
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Creates a new module of the virtual class with the given params.

     The parameters should be sufficient for the derived class to be invoked properly.

     @return A pointer to the newly created object
     */
    virtual ExecBase * clone(const ParameterList & inputParams) const;
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
     @sa ExecStatistics(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadInputData(void);

    /*! \brief Saves all the result data to files specified in the params.

     Saves all the result data generated by the module to files specified
     by values in the params. This is used to either save the data permanantly,
     or to be loaded back in after remote execution is finished and results
     need to be merged by the merge() method.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecStatistics(const ParameterList & input_params), loadInputData(), merge()
     */
    virtual bool saveOutputData(void);

    /*! \brief Saves all the internal data into files specified in the params.

     Saves all the data required for an external process to execute the
     module. The external process would call loadInputData() to reload the
     data into the members before calling invoke(). The user passes a vector
     that will contain the names of all files necessary to run the module as
     a separate process. The first file in this list will always be the main
     parameter file.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecStatistics(const ParameterList & input_params), loadInputData()
     */
    virtual bool saveInputData(std::vector<std::string> & filenames);

    /*! \brief Loads the output data from the files specified in the params.

     Loads all the data from output files into the data members. The purpose
     of this is to ready "child" modules to be merged back together after
     being executed separately.

     @return True if data was loaded successfully, false otherwise.
     @sa ExecStatistics(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadOutputData(void);

    /*! \brief Splits the module into multiple "children" for parallel execution.

     Divides the work required by the module into a vector of sub-modules
     that can be executed in parallel (by the the ParallelExecution() class.
     The split method should divide the work into "nodes" and sub-divide into
     "cpus" for each node. Because the precise implementation of parallel execution
     is not known apriori, nodes may be thought of as large batches, while the
     number of CPUs can be thought of as modules which will be run simultaneously
     within those batches. For example, split may be called with numNodes = 1
     and numCpus = 4 on a quad-processor machine to prep for a threaded implementation.
     Or it may be called with numNodes = 4 and numCpus = 1 to prepare for
     distributed execution on 4 separate machines on a network.

     @param numNodes Number of separate nodes
     @param numCpus Number of CPUs per node
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
     @sa split(int numNodes, int numCpus)
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
    SpecSet * m_spectra; //! The set of spectra for which we will be generating statistics
    bool ownInput; //! Does this object "own" the input data pointers (and hence have to free them)
    MS2ScoringModel * m_model; //! The scoring model used for ion statistics
    SpectrumAnnotParameterList * m_statsParams; //! The statistics that we are considering
    PeptideSpectrumMatchSet * m_peptideResults; //! The peptide results we're considering.
    std::vector<vector<float> > * m_spectraStats; //! Vector of vectors containing peptide statistics
    std::vector<string> * m_spectraStatsHeader; //!Vector containing m_spectraStats header information

    bool ownOutput; //! Does this object "own" the output data pointers (and hence have to free them)
  };

} // namespace specnets

#endif // __ExecStatistics_H__
