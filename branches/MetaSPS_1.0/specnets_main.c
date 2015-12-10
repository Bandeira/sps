/*
 * MATLAB Compiler: 4.13 (R2010a)
 * Date: Tue Dec  7 15:52:56 2010
 * Arguments: "-B" "macro_default" "-R" "-nojvm" "-R" "-nojit" "-m" "-W" "main"
 * "-T" "link:exe" "specnets.m" "-o" "specnets" 
 */
#include <stdio.h>
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

extern mclComponentData __MCC_specnets_component_data;

#ifdef __cplusplus
}
#endif

static HMCRINSTANCE _mcr_inst = NULL;

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultPrintHandler(const char *s)
{
  return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultErrorHandler(const char *s)
{
  int written = 0;
  size_t len = 0;
  len = strlen(s);
  written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
  if (len > 0 && s[ len-1 ] != '\n')
    written += mclWrite(2 /* stderr */, "\n", sizeof(char));
  return written;
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifndef LIB_specnets_C_API
#define LIB_specnets_C_API /* No special import/export declaration */
#endif

LIB_specnets_C_API 
bool MW_CALL_CONV specnetsInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
  if (_mcr_inst != NULL)
    return true;
  if (!mclmcrInitialize())
    return false;
    {
        mclCtfStream ctfStream = 
            mclGetEmbeddedCtfStream((void *)(specnetsInitializeWithHandlers), 
                                    79340);
        if (ctfStream) {
            bResult = mclInitializeComponentInstanceEmbedded(   &_mcr_inst,
                                                                
                                                     &__MCC_specnets_component_data,
                                                                true, 
                                                                NoObjectType, 
                                                                ExeTarget,
                                                                error_handler, 
                                                                print_handler,
                                                                ctfStream, 
                                                                79340);
            mclDestroyStream(ctfStream);
        } else {
            bResult = 0;
        }
    }  
    if (!bResult)
    return false;
  return true;
}

LIB_specnets_C_API 
bool MW_CALL_CONV specnetsInitialize(void)
{
  return specnetsInitializeWithHandlers(mclDefaultErrorHandler, mclDefaultPrintHandler);
}
LIB_specnets_C_API 
void MW_CALL_CONV specnetsTerminate(void)
{
  if (_mcr_inst != NULL)
    mclTerminateInstance(&_mcr_inst);
}

int run_main(int argc, const char **argv)
{
  int _retval;
  /* Generate and populate the path_to_component. */
  char path_to_component[(PATH_MAX*2)+1];
  separatePathName(argv[0], path_to_component, (PATH_MAX*2)+1);
  __MCC_specnets_component_data.path_to_component = path_to_component; 
  if (!specnetsInitialize()) {
    return -1;
  }
  argc = mclSetCmdLineUserData(mclGetID(_mcr_inst), argc, argv);
  _retval = mclMain(_mcr_inst, argc, argv, "specnets", 1);
  if (_retval == 0 /* no error */) mclWaitForFiguresToDie(NULL);
  specnetsTerminate();
  mclTerminateApplication();
  return _retval;
}

int main(int argc, const char **argv)
{
  if (!mclInitializeApplication(
    __MCC_specnets_component_data.runtime_options, 
    __MCC_specnets_component_data.runtime_option_count))
    return 0;

  return mclRunMain(run_main, argc, argv);
}