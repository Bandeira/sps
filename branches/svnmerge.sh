#!/bin/bash

BRANCH=svn://usa.ucsd.edu/sps/branches/$1
#BRANCH=svn://usa.ucsd.edu/sps/branches/ming_branches/trunk_remove_mzxml_4

export TEMP_BRANCH_PREFIX=`echo ${BRANCH} | sed 's/[0-9]*$//'`
export TEMP_BRANCH_SUFFIX=`echo ${BRANCH} | sed 's/[^0-9]*\([0-9]*$\)/\1/'`

#echo $TEMP_BRANCH_SUFFIX
#echo $TEMP_BRANCH_PREFIX

TEMP_BRANCH_SUFFIX=`expr $TEMP_BRANCH_SUFFIX + 1`
#echo $TEMP_BRANCH_SUFFIX

TEMP_BRANCH=$TEMP_BRANCH_PREFIX$TEMP_BRANCH_SUFFIX
#echo $TEMP_BRANCH

#return 0

#TEMP_BRANCH=${BRANCH}_1
#TEMP_BRANCH=svn://usa.ucsd.edu/sps/branches/ming_branches/trunk_remove_mzxml_5
TRUNK=svn://usa.ucsd.edu/sps/trunk

svn copy $TRUNK $TEMP_BRANCH -m "merging from trunk on branch $BRANCH"
svn co $TEMP_BRANCH

export local_new_branch_dir=`echo $TEMP_BRANCH | awk -F/ '{print $NF}'`
export local_old_branch_dir=`echo $BRANCH | awk -F/ '{print $NF}'`
cd $local_old_branch_dir
#cd $local_branch_dir

export REVISION_NUMBER=`svn log  --stop-on-copy | grep  line | grep r | grep 2013 | grep '|' | tail -n 1 | awk '{print $1}'`

cd ..
cd $local_new_branch_dir

echo $local_branch_dir
echo $REVISION_NUMBER

echo svn merge -$REVISION_NUMBER:HEAD $BRANCH .
