#!/bin/bash


################################################################################
# List of routine tests
################################################################################

TestsSpsSmall=('test_sps_small/test_sps_small' 'test_sps_small/test_sps_small_MetaSPS' 'test_sps_small/test_sps_small_genoMS' 'test_sps_small/test_sps_small_genoMS_MetaSPS' 'test_sps_small/test_sps_small_noclusters' 'test_sps_small/test_sps_small_noclusters_MetaSPS' 'test_sps_small/test_sps_small_noclusters_genoMS' 'test_sps_small/test_sps_small_noclusters_genoMS_MetaSPS' 'test_sps_small/test_small_genoMS' 'test_sps_small/test_small_noclusters_genoMS');

TestsShort=('test_snets_tiny' 'test_genoms_small' 'test_sps_ox40L_LC');

TestsLong=('test_genoms_big' 'test_sps_aBTLA_LC' 'test_genoms_mzXML');

TestsSilac=('test_sps_silac' 'test_sps_silac_no_boost' 'test_snets_silac')

TestsAll=(`echo ${TestsShort[*]} ${TestsSpsSmall[*]} ${TestsSilac[*]} ${TestsLong[*]}`)



RUN_SCRIPT=./run_nogrid.sh
SILENT=0

MOVE_SOURCE=report
MOVE_TARGET=/data/sites/projects/sps/jcanhita/Reference

################################################################################
# Get a directory listing (only params files)
################################################################################
function fileList()
{
  local __resultvar=$1

  local LIST2=()
  for f in *.params; do
    LIST2+=($f)
  done
  
  if [[ "$__resultvar" ]]; then
      eval $__resultvar="'${LIST2[@]}'"
  else
      echo "${LIST2[@]}"
  fi
}
################################################################################
#
################################################################################
function delFiles()
{
  declare -a argAry1=("$@")

  for i in ${argAry1[@]}
  do
    rm -f $i
  done
}
################################################################################
#
################################################################################
function appentToParams()
{
  declare -a argAry1=("$@")

  for i in ${argAry1[@]}
  do
    #echo "Appeding to $i"
    echo -e "\n\nREPORT_DYNAMIC=0\n" >> $i
  done

}
################################################################################
#
################################################################################
function testDirectory()
{
  if [ ! -d "$1" ]; then
    echo "$1 not found."
    return 1
  fi
  return 0
}
################################################################################
# Restore
################################################################################
function performRestore()
{
  ProcessMultiple performRestore $1 "Restoring params files in"
  # Get return value
  return_val=$?
  # Test return value. if 2, means 
  if [ "$return_val" -ne 2 ]; then
    return $return_val
  fi

  # checking for test directory existence
  testDirectory $1
  return_val=$?
  if [ "$return_val" -eq 1 ]
  then
    return 1
  fi

  echo "Restoring params files in $1"

  # Move to directory
  pushd "$1" > /dev/null

  DIR=()
  # Get the params files' names
  fileList DIR
  # delete the params files
  delFiles ${DIR[@]}
  # restore from svn
  svn up

  # go back
  popd > /dev/null
  
  return 0
}
################################################################################
# Add
################################################################################
function performAdd()
{
  ProcessMultiple performAdd $1 "Processing params files in"
  # Get return value
  return_val=$?
  # Test return value. if 2, means 
  if [ "$return_val" -ne 2 ]; then
    return $return_val
  fi

  # checking for test directory existence
  testDirectory $1
  return_val=$?
  if [ "$return_val" -eq 1 ]
  then
    return 1
  fi

  echo "Processing params files in $1"

  # Move to directory
  pushd "$1" > /dev/null

  DIR=()
  # Get the params files' names
  fileList DIR
  # delete the params files
  appentToParams ${DIR[@]}

  # go back
  popd > /dev/null
  
  return 0
}
################################################################################
# Params:
# $1 - function to call
# $2 - test(s)
# $3 - message
################################################################################
function ProcessMultiple()
{
  local PROC_LIST=()
  local RET=0
  
  # Process silac tests
  if [ $2 == "silac" ]; then
    echo "$3 silac tests."
    PROC_LIST=${TestsSilac[@]}

    # Process long tests
  elif [ $2 == "long" ]; then
    echo "$3 long tests."
    PROC_LIST=${TestsLong[@]}

  # Process short tests
  elif [ $2 == "short" ]; then
    echo "$3 short tests."
    PROC_LIST=${TestsShort[@]}

  # Process small tests
  elif [ $2 == "small" ]; then
    echo "$3 small tests."
    PROC_LIST=${TestsSpsSmall[@]}

  # Process all tests
  elif [ $2 == "all" ]; then
    echo "$3 all tests."
    PROC_LIST=${TestsAll[@]}

  # If no group test specified, exit
  else
    return 2
  fi
  
  
  # Process the tests
  for i in ${PROC_LIST[@]}
  do
    # Call the function 
    $1 $i
    # Get return value
    return_val=$?
    # Test return value
    if [ "$return_val" -eq 1 ]; then
      RET=1
    fi
  done
  
  # Exit
  return $RET
}
################################################################################
#
################################################################################
function showHelp()
{
  echo
  echo "AddStatic: Spectral Networks QC script helper, version 1.0"
  echo
  echo "Options:"
  echo
  echo "-a --add     <DIR>           Add REPORT_DYNAMIC=0 to specified reports' params files"
  echo "-r --restore <DIR>           Restore original params files in given directories"
  echo
  echo "Available test groups are: small, short, long, all"
  echo
  echo
}
################################################################################
# Main program
################################################################################
  TMP=`getopt --name=$0 -a --longoptions=add:,restore:,help -o a:,r:,h -- $@`
  eval set -- $TMP

  until [ $1 == -- ]; do
      case $1 in
          -a|--add)
              option="ADD"
              target=$2 # option param is in position 2
              shift;; # remove the option's param
          -r|--restore)
              option="RESTORE"
              target=$2 # option param is in position 2
              shift;; # remove the option's param
          -h|--help)
              option="HELP"
              ;; # no parameter taken by the bar option
      esac
      shift # move the arg list to the next option or '--'
  done
  shift # remove the '--', now $1 positioned at first argument if any


# set appropriate directory
currDir=`pwd`


case $option in
"ADD")
  performAdd $target
  ;;
"RESTORE")
  performRestore $target
  ;;
"HELP")
  showHelp
  ;;
*)
  echo "No command option specified. Type 'addStatic.sh --help' for more information."
  $(exit 1)
  ;;
esac
