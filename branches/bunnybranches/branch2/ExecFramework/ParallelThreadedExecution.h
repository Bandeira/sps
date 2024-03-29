#ifndef __ParallelThreadedExecution_H__
#define __ParallelThreadedExecution_H__

// Module Includes
#include "ExecBase.h"
#include "ParallelExecution.h"

// System Includes
#include <pthread.h>
#include <vector>


namespace specnets
{
  /*! \brief Executes a specnets module in parallel on multiple threads

   */
  class ParallelThreadedExecution : public ParallelExecution
  {
  public:

    //! \name CONSTRUCTORS
    //@{
    /*! \brief The constructor

     An execution module is passed to the constructor and it is this module that
     will be executed when invoke is called.
     @sa invoke()
     @param moduleExec the module that will be executed
     */
    ParallelThreadedExecution(ExecBase * moduleExec);
    //@}

    //! \name DESTRUCTOR
    //@{
    ~ParallelThreadedExecution(void);
    //@}

    //! \name MODIFIERS
    //@{
    /*! \brief Executes the module in parallel on separate threads

     Calls methods on the original module to run portions of that module
     in parallel on separate threads. The number of threads used will be
     equal to the numSplits param. The split() method of the original module
     will be used to obtain the sub-modules for execution passing numSplits
     directly to the module. Since all operations are performed in memory
     no calls to the original module to save or load data to files are
     necessary in this implementation. When all threads have completed, the
     merge() method of the original module is called to merge all the sub-module
     back into a choesive result set.

     @param numSplits The number of sub-modules that will be used to execute the module
     @return True if execution of all sub-modules finished without error, false otherwise.
     */
    virtual bool invoke(int numSplits);

    /*! \brief Executes a vector of modules in parallel on separate threads

     Calls the invoke() method on the modules to run each module in parallel on
     separate threads. The number of threads used will be equal to the number
     of modules in the vector.

     @param vecModules The vector of modules to run in parallel
     @return True if execution of all modules finished without error, false otherwise.
     */
    virtual bool invoke(std::vector<ExecBase *> & vecModules);
    //@}

  protected:

  private:
    //! \name NOT IMPLEMENTED
    //@{
    ParallelThreadedExecution(void);
    ParallelThreadedExecution(const ParallelThreadedExecution & that);
    //@}

  };

} // namespace specnets

#endif // __ParallelThreadedExecution_H__
