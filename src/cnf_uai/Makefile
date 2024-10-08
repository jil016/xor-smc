

# Remove the first character in the following line if compile with BOOST Graph library
BOOST = -DLINUX -DLONGLONG_COST -DWIDE_STRING -DLONGDOUBLE_PROB

# Same thing for XML
#XML = XMLFLAG
#XMLI = -I/usr/include/libxml2 -DXMLFLAG
# -DXMLFLAG -DMAXCSP

ifeq ($(BOOST),-DBOOST)
	SRCS = $(SRCSFULL)
else
	SRCS = $(SRCSFULL:tb2boostgraph.cpp=)
endif

ifeq ($(XML),XMLFLAG)
	OBJS = $(SRCS:.cpp=.o) -lxml2 -lgmp
#	-lxml2 -lz -lm -lpthread
else
	OBJS = $(SRCS:.cpp=.o) -lgmp
endif


CC = g++
ICC = icc


ifeq ($(findstring aiopt, $(HOSTNAME)), aiopt)
	ILOGBASE  = /Applications/CPLEX_Studio2211
else
	ILOGBASE  = /Applications/CPLEX_Studio2211
endif

# CPLEX PATH
ILOGLIBFORMAT = arm64_osx/static_pic	# MACOS
ILOGFLAGS = -DILOGLUE -DIL_STD -I$(ILOGBASE)/cpoptimizer/include/ -I$(ILOGBASE)/concert/include/ -I$(ILOGBASE)/cplex/include/
ILOGLIBS  = -L$(ILOGBASE)/cpoptimizer/lib/$(ILOGLIBFORMAT) -L$(ILOGBASE)/concert/lib/$(ILOGLIBFORMAT) -L$(ILOGBASE)/cplex/lib/$(ILOGLIBFORMAT) -lcp -lilocplex -lcplex -lconcert -lm -lpthread

CFLAGS = -g -std=c++17 -Wall $(ILOGFLAGS) $(BOOST) $(XMLI) -DNARYCHAR -DWCSPFORMATONLY

# GMP PATH -- MACOS
GMPFLAGS = -I/opt/homebrew/Cellar/gmp/6.3.0/include/
GMPLIBS = -L/opt/homebrew/Cellar/gmp/6.3.0/lib/

ICCFLAGS = 
#-O3 -ipo -static -DNDEBUG

PFLAGS =
# -pg -g3 # use these options to catch numerical overflows

BUILD = $(CC) -shared

default: all

%.o : %.cpp %.hpp
	$(CC) $(OFLAGS) $(CFLAGS) $(GMPFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

%.o : %.cpp %.h
	$(CC) $(OFLAGS) $(CFLAGS) $(GMPFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

%.o : %.cpp
	$(CC) $(OFLAGS) $(CFLAGS) $(GMPFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

all: SMC


SMC: main_smc.cpp utils.o utils.h smc.o smc.h
	$(CC) $(OFLAGS) $(CFLAGS) $(GMPFLAGS) -o $@ $< utils.o smc.o $(ILOGLIBS) $(GMPLIBS) -L. -ldl #-no-pie -lstdc++fs

clean:
	rm SMC *.o

