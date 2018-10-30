bin_PROGRAMS=gtfmerge

gtfmerge_CPPFLAGS = -pthread -I../lib/gtf -std=c++11
gtfmerge_LDFLAGS = -L../lib/gtf -pthread
gtfmerge_LDADD = -lgtf

gtfmerge_SOURCES = genome1.h genome1.cc \
				   gtfmerge.h gtfmerge.cc \
				   config.h config.cc \
				   main.cc
