#!/bin/bash


SOURCE=../../../../bin
NOTICES=../../../notices
VERSION=$(svnversion ${SOURCE})

GNUPLOT=../../../gnuplot
GNUPLOT_BIN=${GNUPLOT}/linux-x86
GNUPLOT_EXE=${GNUPLOT_BIN}/gnuplot.static
GNUPLOT_COMP=${GNUPLOT}/components
GNUPLOT_NOTICES=${GNUPLOT}/notices


TARGET=specplot
TARGET_BIN=${TARGET}/bin
TARGET_RESOURCE_DIR=${TARGET_BIN}/resources


rm -rf ${TARGET}


mkdir ${TARGET}
mkdir ${TARGET_RESOURCE_DIR}

cp ${NOTICES}/*         ${TARGET}
cp ${GNUPLOT_NOTICES}/* ${TARGET}

cp -r other/* ${TARGET}

cp ${SOURCE}/AA*          ${TARGET_BIN}
cp ${SOURCE}/model_*.txt  ${TARGET_BIN}

cp ${SOURCE}/specplot ${TARGET_BIN}

cp -r ${GNUPLOT_EXE}    ${TARGET_BIN}/gnuplot
cp -r ${GNUPLOT_COMP}/* ${TARGET_BIN}


rm -rf ${TARGET}/*/.svn
rm -rf ${TARGET}/*/*/.svn
rm -rf ${TARGET}/*/*/*/.svn
rm -rf ${TARGET}/*/*/*/*/.svn

              
zip -r -9 specplot-linux-x86-static-${VERSION}.zip ${TARGET}
                                  
                   