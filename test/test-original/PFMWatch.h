// A "timer" measuring hardware performance counter between start and stop
// this is based directly on the libpfm4.4 (perf) library
// this code is adapted from the self-basic.c example of libpfm 
// @author: Sandro Wenzel (sandro.wenzel@cern.ch)

#ifndef PFMWATCH_H
#define PFMWATCH_H


#include <cstdlib>

#include <sys/types.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <locale.h>
#include <sys/ioctl.h>
#include <err.h>
#include <vector>
#include <perfmon/pfmlib_perf_event.h>
#include <iomanip>

#define GROUPSIZE 8 // 8 is the maximum number of counted events ( without multiplexing on my Ivybridge -- HT switched off )


struct PFMWatch 
{
  struct perf_event_attr attr[GROUPSIZE];

  int fd[GROUPSIZE], ret;

  std::vector<std::string> EVENTSTRING; 

  uint64_t count[GROUPSIZE], startvalues[3*GROUPSIZE], stopvalues[3*GROUPSIZE];
  double countoverhead[GROUPSIZE];

  void printSummary(int reps=1)
  {
    std::cout << "# PFM Summary (for measured code section) " << std::endl;
    std::cout << "# --------------------------------------- " << std::endl;
    for(unsigned int i=0;i<GROUPSIZE;i++)
      {
	std::cout << std::left << "# " <<std::setw(25)<<EVENTSTRING[i] << "\t: " << std::fixed << count[i]/(1.*reps) << std::endl;
      }
  }

  uint64_t getCounter(unsigned int event)
  {
    return count[event];
  }

  std::string getEventName(unsigned int event)
  {
    return EVENTSTRING[event];
  }

  uint64_t getNumberOfEvents()
  {
    return GROUPSIZE;
  }



  // this should probably be a singleton class
  PFMWatch() 
  { 
    EVENTSTRING.resize(GROUPSIZE);

    //    EVENTSTRING[0]="cs"; 
    //    EVENTSTRING[1]="migrations"; 
    EVENTSTRING[2]="cycles"; 
    EVENTSTRING[3]="instructions"; 
    EVENTSTRING[4]="branch-misses"; 
    EVENTSTRING[5]="L1-dcache-load-misses"; 
    EVENTSTRING[6]="L1-icache-load-misses"; 
    EVENTSTRING[7]="branches"; 
    // EVENTSTRING[2]="cache-misses"; 
    
    // EV/ENTSTRING[4]="cycles"; 
    // EVENTSTRING[5]="instructions"; 
    EVENTSTRING[0]="idle-cycles-frontend"; 
    EVENTSTRING[1]="cs"; //idle-cycles-backend"; 
 

    ret = pfm_initialize();
    if (ret != PFM_SUCCESS)
      errx(1, "cannot initialize library: %s", pfm_strerror(ret));

    for(unsigned int i=0;i<GROUPSIZE;++i)
      {
	memset(&attr[i], 0, sizeof(perf_event_attr)); // initialize attr;
	count[i]=0;

	/*
	 * 1st argument: event string
	 * 2nd argument: default privilege level (used if not specified in the event string)
	 * 3rd argument: the perf_event_attr to initialize
	 */
	ret = pfm_get_perf_event_encoding( EVENTSTRING[i].c_str(), PFM_PLM0|PFM_PLM3, &attr[i], NULL, NULL);
	if (ret != PFM_SUCCESS)
	  errx(1, "cannot find encoding: %s", pfm_strerror(ret));

	/*
	 * request timing information because event may be multiplexed
	 * and thus it may not count all the time. The scaling information
	 * will be used to scale the raw count as if the event had run all
	 * along
	 */
	attr[i].read_format = PERF_FORMAT_TOTAL_TIME_ENABLED|PERF_FORMAT_TOTAL_TIME_RUNNING;
	attr[i].disabled=1;

	/*
	 * create the event and attach to self
	 * Note that it attaches only to the main thread, there is no inheritance
	 * to threads that may be created subsequently.
	 *
	 * if mulithreaded, then getpid() must be replaced by gettid()
	 */
	if( i == 0 )
	  {
	    fd[i] = perf_event_open(&attr[i], 0, -1, -1, 0);
	    if (fd[i] < 0) 
	      err(1, "cannot create main event");
	  }
	else
	  {
	    // create other events in same event group
	    fd[i]= perf_event_open(&attr[i], 0,-1, fd[0], 0); // this is now an event group with fd
	    if (fd[i] < 0) 
	      err(1, "cannot create child event");
	  }
      }


    for(unsigned int i=0;i<GROUPSIZE;++i)
      {
	// start the counter for each event
	ret = ioctl(fd[i], PERF_EVENT_IOC_ENABLE, 0);
	if (ret)
	  err(1, "ioctl(enable) failed");
      }
  } 


  ~PFMWatch()
  {
    /*
     * stop counting
     */
    for(unsigned int i=0;i<GROUPSIZE;++i)
      {
	ret = ioctl(fd[i], PERF_EVENT_IOC_DISABLE, 0);
	if (ret)
	  err(1, "ioctl(disable) failed");

	close(fd[i]);
      }
    
    /* free libpfm resources cleanly */
    pfm_terminate();
  }


  void Start(){ 
    /*
     * start counting now
     */
    for(unsigned int i=0;i<GROUPSIZE;++i)
      {
	ret = read(fd[i], startvalues + i*3,  sizeof(uint64_t)*3);
	if (ret != sizeof(uint64_t)*3)
	  err(1, "cannot read results: %s", strerror(errno));
      }
  }

  void Stop()
  {   
    /*
     * read out and calculate counters
     */
    for(unsigned int i=0;i<GROUPSIZE;++i)
      {
	ret = read(fd[i], stopvalues + i*3, sizeof(uint64_t)*3);
	if (ret != sizeof(uint64_t)*3)
	  err(1, "cannot read results: %s", strerror(errno));

	if (stopvalues[3*i+2])
	  count[i] = (uint64_t)((double) ( stopvalues[3*i+0] - startvalues[3*i+0] ) * ( stopvalues[3*i+1] - startvalues[3*i+1] )/ ( 1.*(stopvalues[3*i+2] - startvalues[3*i+2] )));
      }
  }


  void HeatUp(){ for(unsigned int i=0;i<5;++i){ Start();Stop(); } }

  double getDeltaSecs() { return (count[1]); }
  //  double getDeltaSecs() { return (count[0]-countoverhead[0])/(count[1]-countoverhead[1]); }

  // find out about the background overhead
  double getOverhead(int N)
  {
    HeatUp();
    unsigned long long Taccum=0L;
    for(unsigned int i=0;i<GROUPSIZE;++i)
      {
	countoverhead[i] = 0;
      }

    for(unsigned int i=0; i < N; i++)
      {
	Start();
	Stop();
	for(unsigned int i=0;i<GROUPSIZE;++i)
	  {
	    //  fprintf(stderr,"#count of event %d is %ld\n", i, count[i]);
	    countoverhead[i] += count[i];
	  }
      }
    for(unsigned int i=0;i<GROUPSIZE;++i)
      {
	countoverhead[i] /= (1.*N);
	fprintf(stderr,"#OVERHEAD of %s = %lf\n", EVENTSTRING[i].c_str(),countoverhead[i]);
      }

    return countoverhead[1];
    //    return Taccum/(1.*N);
  }
};

#endif
