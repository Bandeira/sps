#!/bin/bash


SOURCE=../../../../bin
COMMON=../common
NOTICES=../../../notices
VERSION=$(svnversion ${SOURCE})

GNUPLOT=../../../gnuplot
GNUPLOT_BIN=${GNUPLOT}/linux-x64
GNUPLOT_EXE=${GNUPLOT_BIN}/gnuplot.dynamic
GNUPLOT_LIB=${GNUPLOT_BIN}/lib
GNUPLOT_COMP=${GNUPLOT}/components
GNUPLOT_NOTICES=${GNUPLOT}/notices

DEPENDENCIES_SCRIPT=../../../copyDeps.sh

TARGET=sps
TARGET_BIN=${TARGET}/bin
TARGET_RESOURCE_DIR=${TARGET_BIN}/resources
TARGET_RESOURCE_WEB=${TARGET_RESOURCE_DIR}/css
TARGET_MODEL_DIR=${TARGET_BIN}
TARGET_LIBS=${TARGET_BIN}/libs
TARGET_LIB_GP=${TARGET_BIN}/lib


rm -rf ${TARGET}



mkdir ${TARGET}
mkdir ${TARGET_BIN}
mkdir ${TARGET_RESOURCE_DIR}
mkdir ${TARGET_RESOURCE_WEB}
mkdir ${TARGET_MODEL_DIR}
mkdir ${TARGET_LIBS}
mkdir ${TARGET_LIB_GP}


cp ${NOTICES}/*         ${TARGET}
cp ${GNUPLOT_NOTICES}/* ${TARGET}

cp -r other/*     ${TARGET}
cp -r ${COMMON}/* ${TARGET}

cp -r ${SOURCE}/resources/css/*               ${TARGET_RESOURCE_WEB}
cp -r ${SOURCE}/resources/Models_mscluster    ${TARGET_RESOURCE_DIR}
cp -r ${SOURCE}/resources/Models_pepnovo      ${TARGET_RESOURCE_DIR}
cp -r ${SOURCE}/DBs_GenoMS                    ${TARGET_BIN}

chmod 0777 ${TARGET_BIN}/DBs_GenoMS
chmod 0777 ${TARGET_BIN}/DBs_GenoMS/*

cp ${SOURCE}/clustalw           ${TARGET_BIN}
cp ${SOURCE}/convert            ${TARGET_BIN}
cp ${SOURCE}/MsCluster_bin      ${TARGET_BIN}
cp ${SOURCE}/PepNovo_bin        ${TARGET_BIN}
cp ${SOURCE}/MSGFDB.jar         ${TARGET_BIN}
cp ${SOURCE}/GenoMS.jar         ${TARGET_BIN}
cp ${SOURCE}/AA*                ${TARGET_BIN}

cp ${SOURCE}/model_cid.txt      ${TARGET_MODEL_DIR}
cp ${SOURCE}/model_etd.txt      ${TARGET_MODEL_DIR}
cp ${SOURCE}/model_prm.txt      ${TARGET_MODEL_DIR}

cp ${SOURCE}/main_execmodule  ${TARGET_BIN}
cp ${SOURCE}/main_specnets    ${TARGET_BIN}

cp ${SOURCE}/specplot         ${TARGET_BIN}
cp ${SOURCE}/contplot         ${TARGET_BIN}
cp ${SOURCE}/spsReport        ${TARGET_BIN}
cp ${SOURCE}/spsReportServer  ${TARGET_BIN}
cp ${SOURCE}/clusterinfo      ${TARGET_BIN}

cp -r ${GNUPLOT_EXE}    ${TARGET_BIN}/gnuplot
cp -r ${GNUPLOT_LIB}    ${TARGET_BIN}
cp -r ${GNUPLOT_COMP}/* ${TARGET_BIN}

#${DEPENDENCIES_SCRIPT} ${TARGET_BIN}/gnuplot        ${TARGET_LIB_GP}  &> /dev/null
${DEPENDENCIES_SCRIPT} ${TARGET_BIN}/clustalw       ${TARGET_LIBS}    &> /dev/null
${DEPENDENCIES_SCRIPT} ${TARGET_BIN}/convert        ${TARGET_LIBS}    &> /dev/null
${DEPENDENCIES_SCRIPT} ${TARGET_BIN}/MsCluster_bin  ${TARGET_LIBS}    &> /dev/null
${DEPENDENCIES_SCRIPT} ${TARGET_BIN}/PepNovo_bin    ${TARGET_LIBS}    &> /dev/null
rm ${TARGET_LIB_GP}/libc.so.6
rm ${TARGET_LIBS}/libc.so.6


rm -rf ${TARGET}/.svn
rm -rf ${TARGET}/*/.svn
rm -rf ${TARGET}/*/*/.svn
rm -rf ${TARGET}/*/*/*/.svn
rm -rf ${TARGET}/*/*/*/*/.svn

              
zip -r -9 sps-linux-x64-dynamic-${VERSION}.zip ${TARGET}       