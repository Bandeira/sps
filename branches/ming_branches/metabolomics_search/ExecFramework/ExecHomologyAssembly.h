#ifndef __ExecHomologyAssembly_H__
#define __ExecHomologyAssembly_H__

// Module Includes
#include "ExecBase.h"

// SpecNets Includes
#include "spectrum.h"
#include "db_fasta.h"
#include "SetMerger.h"

// System Includes
#include <string>
#include <vector>

namespace specnets
{
	using namespace std;

	/*! \brief Execution class for sequence-tag-based matching of spectra to database sequences

   Matches spectra to database peptide sequences by deriving tags from each spectrum
   and matching the resulting tags to a protein sequence database.

   */
  class ExecHomologyAssembly : public ExecBase
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
  	ExecHomologyAssembly(void);

    /*! \brief Constructor which takes only the set of parameters

     Constructor which takes only the set of parameters required to execute the
     module. The parameters specify all data files for input and output, and requires
     that the user call loadInputData() in order to read all data from files into
     the members of the object. The parameters can be verified using the
     validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     */
  	ExecHomologyAssembly(const ParameterList & inputParams);

    /*! \brief Constructor which takes both the parameters and the data structures

     Constructor which takes both the parameters and the data structures necessary for
     executing the alignment. When using this constructor no external data need be read in.
     The parameters can be verified using the validateParams() method.
     If matchedSpectra is not NULL then invoke() returns only info for matched spectra
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     @param contigs Input spectra as contigs (used instead of spectra if both are set)
     @param spectra Input spectra as a set of spectra
     @param db Protein sequences
     @param matchedProts Per-spectrum: best-matched protein (-1 if none), # mods, prefix/suffix match (cols 0-2, respectively)
     @param matchedSpectra Matched versions of spectra significantly matched to a protein
     @param matchedSpectraIndices Indices of spectra significantly matched to a protein
     */
	  ExecHomologyAssembly(const ParameterList & inputParams,
				                 SpecSet * spectra,
				                 DB_fasta * db,
				                 SpecSet * contigSpectra,
				                 SpecSet * cspsMatchedIndices,
				                 vector<vector<int> > * cspsMatchedProts);

    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ExecHomologyAssembly(void);
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Creates a new module of the virtual class with the given params.

     The parameters should be sufficient for the derived class to be invoked properly.

     @return A pointer to the newly created object
     */
    virtual ExecBase * clone(const ParameterList & input_params) const;

    //@}

    //! \name MODIFIERS
    //@{
    /*! \brief Executes the module.

     In order to call this method succesfully all the necessary data for
     execution must already be loaded into memory (data members). This can
     be accomplished using the loadInputData() method or by calling the
     appropriate constructor.

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
     @sa ExecBase(const ParameterList & input_params), saveInputData()
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
  	SpecSet * m_spectra; //! Input spectra
	  DB_fasta * m_db; //! Input database of protein sequences
	  bool ownInput; //! Does this object "own" the input data pointers (and hence has to free them)

  	SpecSet * m_contigs; //! Assembled contig spectra
	  SpecSet * m_cspsMatchedIndices; //! Per-contig indices of matched spectrum/protein masses
	  vector<vector<int> > * m_cspsMatchedProts; //! Per-contig: best-matched protein (-1 if none), # mods, prefix/suffix match (cols 0-2, respectively)
	  bool ownOutput; //! Does this object "own" the output data pointers (and hence has to free them)

	  // Internal
	  SetMerger m_components; //! Sets of spectra assembled in the same contigs/spectral-networks
	  vector<vector<list<TwoValues<int> > > > m_abVertices;  //! ABruijn vertices for the selected contigs/spectral-networks
	  vector<bool> m_specFlipped; //! Indicates whether spectra are/are-not reversed in their output orientation
  };

} // namespace specnets

#endif // __ExecHomologyAssembly_H__