bin_PROGRAMS=gtfcuff

gtfcuff_CPPFLAGS = -I$(top_srcdir)/lib/gtf
gtfcuff_LDFLAGS = -L$(top_builddir)/lib/gtf
gtfcuff_LDADD = -lgtf

.PHONY: pwd
pwd:
	@echo -n "pwd "; pwd
	@echo "LDFLAGS: $(gtfcuff_LDFLAGS)"

gtfcuff_SOURCES = gtfcuff.h gtfcuff.cc \
				  cuffitem.h cuffitem.cc \
				  quantitem.h quantitem.cc \
				  main.cc
