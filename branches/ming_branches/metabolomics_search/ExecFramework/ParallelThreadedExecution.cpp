// Header Include
#include "ExecBase.h"
#include "Logger.h"
#include "ParallelThreadedExecution.h"

// System Includes
#include <vector>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
// This is the function executed by the individual threads
// -------------------------------------------------------------------------
static void * threadFunc(void * userData)
{
  DEBUG_TRACE;
  ExecBase * module = (ExecBase *)userData;
  DEBUG_TRACE;
  module->invoke();
  DEBUG_TRACE;
}

// -------------------------------------------------------------------------
ParallelThreadedExecution::ParallelThreadedExecution(ExecBase * moduleExec) :
  ParallelExecution(moduleExec)
{
  // EMPTY
}

// -------------------------------------------------------------------------
ParallelThreadedExecution::~ParallelThreadedExecution(void)
{
  // EMPTY
}

// -------------------------------------------------------------------------
bool ParallelThreadedExecution::invoke(int numSplits)
{
  DEBUG_VAR(numSplits);

  vector<ExecBase *> const & subModules = m_moduleExec->split(numSplits);

  DEBUG_VAR(subModules.size());

  vector<pthread_t> vecThread(numSplits); 
  for (int i = 0; i < numSplits; i++)
  {
    pthread_t thread;

    DEBUG_MSG("Creating thread " << i);

    int returnValue = pthread_create(&thread,
                                     NULL,
                                     &threadFunc,
                                     (void *)subModules[i]);

    DEBUG_MSG("Created thread " << i);

    if (returnValue)
    {
      ERROR_MSG("Return code from pthread_create() is " << returnValue );
      return false;
    }

    vecThread[i] = thread;
    DEBUG_MSG("Joining thread " << i);
    //pthread_join(vecThread[i], NULL);

    //subModules[i]->invoke();
  }

  DEBUG_MSG("Waiting for threads");
  for (int i = 0; i < numSplits; i++)
  {
    pthread_join(vecThread[i], NULL);
  }

  DEBUG_MSG("All threads completed");

  // Merge back the results
  m_moduleExec->merge();

  DEBUG_TRACE;

  return true;
}

// -------------------------------------------------------------------------
bool ParallelThreadedExecution::invoke(vector<ExecBase *> & vecModules)
{
  vector<pthread_t> vecThread(vecModules.size()); 
  for (int i = 0; i < vecModules.size(); i++)
  {
    pthread_t thread;

    DEBUG_MSG("Creating thread " << i);

    int returnValue = pthread_create(&thread,
                                     NULL,
                                     &threadFunc,
                                     (void *)vecModules[i]);

    DEBUG_MSG("Created thread " << i);

    if (returnValue)
    {
      ERROR_MSG("Return code from pthread_create() is " << returnValue );
      return false;
    }

    vecThread[i] = thread;
    DEBUG_MSG("Joining thread " << i);
    //pthread_join(vecThread[i], NULL);

    //subModules[i]->invoke();
  }

  DEBUG_MSG("Waiting for threads");
  for (int i = 0; i < vecModules.size(); i++)
  {
    pthread_join(vecThread[i], NULL);
  }

  DEBUG_MSG("All threads completed");

  return true;
}

