#ifndef _SgeGridMonitor_H_
#define _SgeGridMonitor_H_

// Module Includes

// External Module Includes

// System Includes
#include <map>
#include <string>

namespace specnets
{
  class SgeGridMonitor
  {
  public:
    enum JobStatus
    {
      UNKNOWN = 0,
      DELETION = 1,
      ERROR = 2,
      HOLD = 3,
      RUNNING = 4,
      RESTARTED = 5,
      SUSPENDED = 6,
      WAITING = 7,
      DONE = 8
    };
    
    SgeGridMonitor(void);
    ~SgeGridMonitor(void);

    JobStatus getJobStatus(const std::string & jobId);
    JobStatus getJobStatus(int jobId);
    
    bool refreshInfo(const std::string & sgePath);

    std::string submitJob(const std::string & sgePath, 
                          const std::string & jobFilename, 
                          const std::string & params,
                          const std::string & gridType);
    
  private:
    std::map<std::string, JobStatus> jobs;
  };    

} //namespace specnets

#endif // _GridUtils_H_

