bin_PROGRAMS = bamkit

bamkit_CPPFLAGS = -std=c++11 -I$(HTSLIB)

# static link
bamkit_LDFLAGS = -pthread
bamkit_LDADD = $(HTSLIB)/libhts.a -lbz2 -lz

bamkit_SOURCES = hit.h hit.cc \
				 bamkit.h bamkit.cc \
				 config.h config.cc \
				 util.h util.cc \
				 main.cc
