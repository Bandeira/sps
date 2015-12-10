///////////////////////////////////////////////////////////////////////////////
#ifndef __CONDITIONAL_H__
#define __CONDITIONAL_H__
///////////////////////////////////////////////////////////////////////////////
// Conditional compilation repository
//
// Used to set differeces when compiling for Linux, Cygwin or Unix
///////////////////////////////////////////////////////////////////////////////
// Linux section

#if defined(__linux__)


// Includes for Linux
#include <cxxabi.h>
#include <execinfo.h>
#include <sys/wait.h>
#include <sys/sysinfo.h>


// CygWin Section
#elif defined(__MINGW32__) || defined(__CYGWIN__)

#include <windows.h>

#endif


#if defined(__MINGW32__)
#define abort() return -9

#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif // defined(__MINGW32__) || defined(__CYGWIN__)

#if !defined(ACCESSPERMS)
#define ACCESSPERMS (S_IRWXU|S_IRWXG|S_IRWXO)
#endif


///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
