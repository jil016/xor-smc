

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

# ILOGBASE  = /home/li4255/.local/ibm/ILOG/CPLEX_Studio201
ifeq ($(findstring aiopt, $(HOSTNAME)), aiopt)
	ILOGBASE  = /home/jinzhao/.local/ibm/ILOG/CPLEX_Studio2211
	# /home/jinzhao/jinzhao/CPLEX/CPLEX_Studio201
else
	ILOGBASE  = /home/jinzhao/.local/ibm/ILOG/CPLEX_Studio2211
	# /home/jinzhao/jinzhao/CPLEX/CPLEX_Studio201
#/opt/ibm/ILOG/CPLEX_Studio201
endif
#ILOGLIBFORMAT = x86-64_sles10_4.1/static_pic
ILOGLIBFORMAT = x86-64_linux/static_pic
ILOGFLAGS = -DILOGLUE -DIL_STD -I$(ILOGBASE)/cpoptimizer/include/ -I$(ILOGBASE)/concert/include/ -I$(ILOGBASE)/cplex/include/
ILOGLIBS  = -L$(ILOGBASE)/cpoptimizer/lib/$(ILOGLIBFORMAT) -L$(ILOGBASE)/concert/lib/$(ILOGLIBFORMAT) -L$(ILOGBASE)/cplex/lib/$(ILOGLIBFORMAT) -lcp -lilocplex -lcplex -lconcert -lm -lpthread

CFLAGS = -g -std=c++17 -Wall $(ILOGFLAGS) $(BOOST) $(XMLI) -DNARYCHAR -DWCSPFORMATONLY


ICCFLAGS = 
#-O3 -ipo -static -DNDEBUG

PFLAGS = 

# -pg -g3 # use these options to catch numerical overflows

BUILD = $(CC) -shared

default: all

%.o : %.cpp %.hpp
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

%.o : %.cpp %.h
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

%.o : %.cpp
	$(CC) $(OFLAGS) $(CFLAGS) $(LIBFLAGS) -c -o $@  $< $(PFLAGS)

all: SL SPC

SL: main_sl.cpp graph.o graph.h utils.o utils.h shelter_location.o shelter_location.h
	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< graph.o utils.o shelter_location.o $(ILOGLIBS) -L. -lgmp -ldl -no-pie -lstdc++fs

SPC: main_spc.cpp utils.o utils.h supply_net.o supply_net.h supply_chain.o supply_chain.h
	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< utils.o supply_net.o supply_chain.o $(ILOGLIBS) -L. -lgmp -ldl -no-pie -lstdc++fs

# SMC: main.cpp graph.o graph.h utils.o utils.h shelter_location.o shelter_location.h
# 	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< graph.o utils.o shelter_location.o $(ILOGLIBS) -L. -lgmp -ldl -no-pie


clean:
	rm SL SPC *.o

