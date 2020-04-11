bin_PROGRAMS=gtfsimulator

gtfsimulator_CPPFLAGS = -I$(top_srcdir)/lib/gtf -I$(top_srcdir)/lib/util
gtfsimulator_LDFLAGS = -L$(top_srcdir)/lib/gtf -L$(top_srcdir)/lib/util
gtfsimulator_LDADD = -lgtf -lutil

gtfsimulator_SOURCES = main.cc \
					   simulator.h simulator.cc 
