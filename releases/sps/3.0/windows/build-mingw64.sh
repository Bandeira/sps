#!/bin/bash


SOURCE=../../../../bin
COMMON=../common
NOTICES=../../../notices
VERSION=$(svnversion ${SOURCE})

GNUPLOT=../../../gnuplot
GNUPLOT_BIN=${GNUPLOT}/windows
GNUPLOT_COMP=${GNUPLOT}/components
GNUPLOT_NOTICES=${GNUPLOT}/notices

TARGET=sps
BIN=${TARGET}/bin
RESOURCE_DIR=${BIN}/resources
RESOURCE_WEB=${RESOURCE_DIR}/css
MODEL_DIR=${BIN}


rm -rf ${TARGET}


mkdir ${TARGET}
mkdir ${BIN}
mkdir ${RESOURCE_DIR}
mkdir ${RESOURCE_WEB}
mkdir ${MODEL_DIR}


cp ${NOTICES}/*         ${TARGET}
cp ${GNUPLOT_NOTICES}/* ${TARGET}

cp -r other/*           ${TARGET}
cp -r mingw64/*         ${TARGET}
cp -r ${COMMON}/*       ${TARGET}


cp /bin/sh.exe              ${BIN}
cp /bin/cygbz2-1.dll        ${BIN}
cp /bin/cygX11-6.dll        ${BIN}
cp /bin/cygXau-6.dll        ${BIN}
cp /bin/cygXaw-7.dll        ${BIN}
cp /bin/cygXdmcp-6.dll      ${BIN}
cp /bin/cygXpm-4.dll        ${BIN}
cp /bin/cygcairo-2.dll      ${BIN}
cp /bin/cygcrypto-0.9.8.dll ${BIN}
cp /bin/cygcurl-4.dll       ${BIN}
cp /bin/cygexpat-1.dll      ${BIN}
cp /bin/cygfontconfig-1.dll ${BIN}
cp /bin/cygfreetype-6.dll   ${BIN}
cp /bin/cyggcc_s-1.dll      ${BIN}
cp /bin/cyggd-2.dll         ${BIN}
cp /bin/cygiconv-2.dll      ${BIN}
cp /bin/cygicudata38.dll    ${BIN}
cp /bin/cygicuuc38.dll      ${BIN}
cp /bin/cygidn-11.dll       ${BIN}
cp /bin/cygintl-8.dll       ${BIN}
cp /bin/cygjpeg-7.dll       ${BIN}
cp /bin/cygjpeg-8.dll       ${BIN}
cp /bin/cygncursesw-10.dll  ${BIN}
cp /bin/cygpng12.dll        ${BIN}
cp /bin/cygpng14-14.dll     ${BIN}
cp /bin/cygreadline7.dll    ${BIN}
cp /bin/cygssh2-1.dll       ${BIN}
cp /bin/cygssl-0.9.8.dll    ${BIN}
cp /bin/cygstdc++-6.dll     ${BIN}
cp /bin/cyguuid-1.dll       ${BIN}
cp /bin/cygwin1.dll         ${BIN}
cp /bin/cygxcb-1.dll        ${BIN}
cp /bin/cygz.dll            ${BIN}
cp /bin/gnuplot.dll         ${BIN}
cp /bin/libcairo-2.dll      ${BIN}
cp /bin/libfontconfig-1.dll ${BIN}
cp /bin/libstdc++-6.dll     ${BIN}
cp /bin/cp.exe              ${BIN}
cp /bin/cygattr-1.dll       ${BIN}
cp /bin/cygxerces-c-3-0.dll ${BIN}
cp /bin/cygxerces-c-3-1.dll ${BIN}
cp /bin/mkdir.exe           ${BIN}
cp /bin/rm.exe              ${BIN}
cp /bin/rmdir.exe           ${BIN}
cp /bin/diff.exe            ${BIN}

cp -r ${SOURCE}/resources/css/*            ${RESOURCE_WEB}
cp -r ${SOURCE}/resources/Models_mscluster ${RESOURCE_DIR}
cp -r ${SOURCE}/resources/Models_pepnovo   ${RESOURCE_DIR}
cp -r ${SOURCE}/DBs_GenoMS                 ${BIN}

chmod 0777 ${BIN}/DBs_GenoMS
chmod 0777 ${BIN}/DBs_GenoMS/*

cp ${SOURCE}/clustalw.exe        ${BIN}
cp ${SOURCE}/convert.exe         ${BIN}
cp ${SOURCE}/MsCluster_bin.exe   ${BIN}
cp ${SOURCE}/PepNovo_bin.exe     ${BIN}
cp ${SOURCE}/MSGFDB.jar          ${BIN}
cp ${SOURCE}/GenoMS.jar          ${BIN}
cp ${SOURCE}/AA*                 ${BIN}

cp ${SOURCE}/main_execmodule.exe   ${BIN}
cp ${SOURCE}/main_specnets.exe     ${BIN}

cp ${SOURCE}/specplot.exe          ${BIN}
cp ${SOURCE}/contplot.exe          ${BIN}
cp ${SOURCE}/spsReport.exe         ${BIN}
cp ${SOURCE}/spsReportServer.exe   ${BIN}
cp ${SOURCE}/clusterinfo.exe       ${BIN}

cp -r ${GNUPLOT_BIN}/*      ${BIN}
cp -r ${GNUPLOT_COMP}/*     ${BIN}

cp ${SOURCE}/model_cid.txt  ${MODEL_DIR}
cp ${SOURCE}/model_etd.txt  ${MODEL_DIR}
cp ${SOURCE}/model_prm.txt  ${MODEL_DIR}



rm -rf ${TARGET}/.svn
rm -rf ${TARGET}/*/.svn
rm -rf ${TARGET}/*/*/.svn
rm -rf ${TARGET}/*/*/*/.svn
rm -rf ${TARGET}/*/*/*/*/.svn


zip -r -9 sps-mingw64-${VERSION}.zip ${TARGET}

