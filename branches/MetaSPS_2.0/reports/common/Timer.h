#ifndef __SPS_TIMER_H__
#define __SPS_TIMER_H__

#include <ctime>


class Timer_c {

  clock_t t1, t2;
  bool running;
  
 public:
  
 Timer_c() : running(true) {start();};
 ~Timer_c() {};
 
 void   start(void)   {t1 = clock(); running = true;};
 double restart(void) {double ret = get(); t1 = clock(); return ret;};
 double stop(void)    {double ret = get(); running = false; return ret;};
 
 double get(void)     {t2 = clock(); if(running) return ((double)(t2 - t1) / CLOCKS_PER_SEC);return 0.0;};

};

#endif


