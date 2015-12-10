/*
 * MATLAB Compiler: 4.13 (R2010a)
 * Date: Tue Dec  7 15:52:42 2010
 * Arguments: "-B" "macro_default" "-R" "-nojvm" "-R" "-nojit" "-m" "-W" "main"
 * "-T" "link:exe" "csps.m" "-o" "csps" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_csps_session_key[] = {
    '6', '2', '7', '0', '3', 'B', 'C', '6', '2', '7', '9', 'F', 'A', '7', 'A',
    '5', 'B', 'E', '2', 'C', '9', 'B', '1', 'F', 'E', 'E', '0', '7', '3', '3',
    'C', '8', '8', '6', 'A', '5', '5', 'E', 'A', '2', 'F', '4', 'B', '0', '5',
    'C', 'E', '8', '1', '3', '3', '5', '9', '4', 'F', '9', '2', '0', 'A', '9',
    '8', '3', '4', '8', '9', '7', '8', 'F', '4', '6', '5', '6', 'C', 'B', '6',
    'F', '5', '5', '1', '1', '6', '0', 'D', '5', '8', 'F', '0', '8', 'D', 'C',
    'A', '0', 'D', '1', '5', '2', 'B', '0', '8', 'D', '0', 'F', '3', '3', '1',
    'C', 'F', '6', '4', '1', '3', 'F', '8', '6', '9', 'A', 'E', 'C', 'E', '2',
    '0', '3', '3', '4', '4', '7', 'E', '0', '0', 'D', '7', '8', 'E', 'C', '2',
    '7', 'D', 'E', 'D', 'E', '1', '6', 'A', 'B', 'C', '3', '0', '9', 'A', 'E',
    '2', 'A', 'C', 'F', 'E', 'B', 'D', '7', '1', '0', 'F', 'B', 'B', 'D', 'F',
    '4', '2', '5', '1', 'B', '5', '3', 'C', 'E', '1', 'D', '9', '5', '8', '4',
    '5', 'B', '7', '0', 'F', '9', '5', '7', '1', '6', '3', '3', '3', 'E', 'A',
    '8', 'B', '9', 'E', '7', '6', '5', '2', '3', 'E', '7', 'A', '1', 'F', '6',
    '4', 'D', 'C', '8', '5', '9', 'F', '5', 'D', 'B', 'B', 'A', 'E', '7', '8',
    'E', 'B', 'C', '0', 'C', '6', 'F', 'C', 'C', '1', '2', 'F', 'C', '8', '1',
    'B', 'D', 'B', 'E', '4', 'F', '4', 'C', 'C', '4', '4', '0', 'E', 'A', '9',
    '1', '\0'};

const unsigned char __MCC_csps_public_key[] = {
    '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9', '2',
    'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1', '0', '1',
    '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B', '0', '0', '3',
    '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1', '0', '0', 'C', '4',
    '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3', 'A', '5', '2', '0', '6',
    '5', '8', 'F', '6', 'F', '8', 'E', '0', '1', '3', '8', 'C', '4', '3', '1',
    '5', 'B', '4', '3', '1', '5', '2', '7', '7', 'E', 'D', '3', 'F', '7', 'D',
    'A', 'E', '5', '3', '0', '9', '9', 'D', 'B', '0', '8', 'E', 'E', '5', '8',
    '9', 'F', '8', '0', '4', 'D', '4', 'B', '9', '8', '1', '3', '2', '6', 'A',
    '5', '2', 'C', 'C', 'E', '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4',
    'D', '0', '8', '5', 'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2',
    'E', 'D', 'E', '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6',
    '3', '7', '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E',
    '6', '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
    '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1', 'B',
    'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9', '9', '0',
    '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0', 'B', '6', '1',
    'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B', '5', '8', 'F', 'C',
    '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6', 'E', 'B', '7', 'E', 'C',
    'D', '3', '1', '7', '8', 'B', '5', '6', 'A', 'B', '0', 'F', 'A', '0', '6',
    'D', 'D', '6', '4', '9', '6', '7', 'C', 'B', '1', '4', '9', 'E', '5', '0',
    '2', '0', '1', '1', '1', '\0'};

static const char * MCC_csps_matlabpath_data[] = 
  { "csps/", "$TOOLBOXDEPLOYDIR/", "$TOOLBOXMATLABDIR/general/",
    "$TOOLBOXMATLABDIR/ops/", "$TOOLBOXMATLABDIR/lang/",
    "$TOOLBOXMATLABDIR/elmat/", "$TOOLBOXMATLABDIR/randfun/",
    "$TOOLBOXMATLABDIR/elfun/", "$TOOLBOXMATLABDIR/specfun/",
    "$TOOLBOXMATLABDIR/matfun/", "$TOOLBOXMATLABDIR/datafun/",
    "$TOOLBOXMATLABDIR/polyfun/", "$TOOLBOXMATLABDIR/funfun/",
    "$TOOLBOXMATLABDIR/sparfun/", "$TOOLBOXMATLABDIR/scribe/",
    "$TOOLBOXMATLABDIR/graph2d/", "$TOOLBOXMATLABDIR/graph3d/",
    "$TOOLBOXMATLABDIR/specgraph/", "$TOOLBOXMATLABDIR/graphics/",
    "$TOOLBOXMATLABDIR/uitools/", "$TOOLBOXMATLABDIR/strfun/",
    "$TOOLBOXMATLABDIR/imagesci/", "$TOOLBOXMATLABDIR/iofun/",
    "$TOOLBOXMATLABDIR/audiovideo/", "$TOOLBOXMATLABDIR/timefun/",
    "$TOOLBOXMATLABDIR/datatypes/", "$TOOLBOXMATLABDIR/verctrl/",
    "$TOOLBOXMATLABDIR/codetools/", "$TOOLBOXMATLABDIR/helptools/",
    "$TOOLBOXMATLABDIR/demos/", "$TOOLBOXMATLABDIR/timeseries/",
    "$TOOLBOXMATLABDIR/hds/", "$TOOLBOXMATLABDIR/guide/",
    "$TOOLBOXMATLABDIR/plottools/", "toolbox/local/",
    "$TOOLBOXMATLABDIR/datamanager/", "toolbox/compiler/" };

static const char * MCC_csps_classpath_data[] = 
  { "" };

static const char * MCC_csps_libpath_data[] = 
  { "" };

static const char * MCC_csps_app_opts_data[] = 
  { "" };

static const char * MCC_csps_run_opts_data[] = 
  { "-nojvm", "-nojit" };

static const char * MCC_csps_warning_state_data[] = 
  { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_csps_component_data = { 

  /* Public key data */
  __MCC_csps_public_key,

  /* Component name */
  "csps",

  /* Component Root */
  "",

  /* Application key data */
  __MCC_csps_session_key,

  /* Component's MATLAB Path */
  MCC_csps_matlabpath_data,

  /* Number of directories in the MATLAB Path */
  37,

  /* Component's Java class path */
  MCC_csps_classpath_data,
  /* Number of directories in the Java class path */
  0,

  /* Component's load library path (for extra shared libraries) */
  MCC_csps_libpath_data,
  /* Number of directories in the load library path */
  0,

  /* MCR instance-specific runtime options */
  MCC_csps_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_csps_run_opts_data,
  /* Number of MCR global runtime options */
  2,
  
  /* Component preferences directory */
  "csps_A189C2E1C12AC89AC392B6F9E822BD87",

  /* MCR warning status data */
  MCC_csps_warning_state_data,
  /* Number of MCR warning status modifiers */
  1,

  /* Path to component - evaluated at runtime */
  NULL

};

#ifdef __cplusplus
}
#endif


