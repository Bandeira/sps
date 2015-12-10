SUBDIRS = Common MLlib PepNovo MsCluster
	
.PHONY: subdirs $(SUBDIRS) clean all
	
subdirs: $(SUBDIRS)
	
$(SUBDIRS):
	$(MAKE) -C $@ LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"

all: subdirs

Models:
	cd ./Models ; ./dos2unix *.* ; dos2unix */*.*

clean:
	cd ./Common ; $(MAKE) clean
	cd ./MLlib ; $(MAKE) clean
	cd ./PepNovo ; $(MAKE) clean
	cd ./MsCluster ; $(MAKE) clean
