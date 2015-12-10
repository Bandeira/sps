#!/bin/bash


SOURCE=../../../../bin
COMMON=../common
NOTICES=../../../notices
VERSION=$(svnversion ${SOURCE})

GNUPLOT=../../../gnuplot
GNUPLOT_BIN=${GNUPLOT}/linux-x64
GNUPLOT_EXE=${GNUPLOT_BIN}/gnuplot.static
GNUPLOT_COMP=${GNUPLOT}/components
GNUPLOT_NOTICES=${GNUPLOT}/notices


TARGET=sps
TARGET_BIN=${TARGET}/bin
TARGET_RESOURCE_DIR=${TARGET_BIN}/resources
TARGET_RESOURCE_WEB=${TARGET_RESOURCE_DIR}/css
TARGET_MODEL_DIR=${TARGET_BIN}


rm -rf ${TARGET}


mkdir ${TARGET}
mkdir ${TARGET_BIN}
mkdir ${TARGET_RESOURCE_DIR}
mkdir ${TARGET_RESOURCE_WEB}
mkdir ${TARGET_MODEL_DIR}


cp ${NOTICES}/*         ${TARGET}
cp ${GNUPLOT_NOTICES}/* ${TARGET}

cp -r other/*     ${TARGET}
cp -r ${COMMON}/* ${TARGET}

cp -r ${SOURCE}/resources/css/*               ${TARGET_RESOURCE_WEB}
cp -r ${SOURCE}/resources/Models_mscluster    ${TARGET_RESOURCE_DIR}
cp -r ${SOURCE}/resources/Models_pepnovo      ${TARGET_RESOURCE_DIR}
cp -r ${SOURCE}/DBs_GenoMS                    ${TARGET_BIN}

chmod 0777 ${BIN}/DBs_GenoMS
chmod 0777 ${BIN}/DBs_GenoMS/*

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
cp -r ${GNUPLOT_COMP}/* ${TARGET_BIN}


rm -rf ${TARGET}/.svn
rm -rf ${TARGET}/*/.svn
rm -rf ${TARGET}/*/*/.svn
rm -rf ${TARGET}/*/*/*/.svn
rm -rf ${TARGET}/*/*/*/*/.svn

              
zip -r -9 sps-linux-x64-static-${VERSION}.zip ${TARGET}       