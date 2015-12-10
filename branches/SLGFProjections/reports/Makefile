################################################################################
# Definitions
################################################################################
BASE_DIR := ../..
#-------------------------------------------------------------------------------
include $(BASE_DIR)/Makefiles/Makefile.defs

################################################################################
# Conditions
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.cond

################################################################################
# Compile rules - general
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.general

################################################################################
# Child directories
################################################################################
SUBDIRS = \
  common \
  spsReports.dir \
  contplot.dir \
  specplot.dir \
  uu64Encode.dir \
  clusterInfo.dir \
  spsReportServer.dir


.PHONY: subdirs $(SUBDIRS) clean all

ifeq ($(type),debug)
all: subdirs
	cp spsReports.dir/spsReport spsReport
	cp spsReportServer.dir/spsReportServer spsReportServer
	cp specplot.dir/specplot specplot
	cp contplot.dir/contplot contplot
	cp uu64Encode.dir/uu64Encode uu64Encode
	cp clusterInfo.dir/clusterinfo
else
all: subdirs
	cp spsReports.dir/spsReport spsReport
	cp spsReportServer.dir/spsReportServer spsReportServer
	cp specplot.dir/specplot specplot
	cp contplot.dir/contplot contplot
	cp uu64Encode.dir/uu64Encode uu64Encode
	cp clusterInfo.dir/clusterinfo clusterinfo
	$(STRIP) spsReport spsReportServer specplot contplot uu64Encode clusterinfo
endif


subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	cd spsReports.dir ; make clean ; cd ..
	cd spsReportServer.dir ; make clean ; cd ..
	cd specplot.dir ; make clean ; cd ..
	cd contplot.dir ; make clean ; cd ..
	cd uu64Encode.dir ; make clean ; cd ..
	cd clusterInfo.dir ; make clean ; cd ..
	cd common ; make clean ; cd ..
	rm spsReport spsReportServer specplot contplot uu64Encode clusterinfo depends

################################################################################

