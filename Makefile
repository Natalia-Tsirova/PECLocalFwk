# Check if an installation path of PECFwk is provided
ifeq ($(PEC_FWK_INSTALL), )
  $(error Mandatory environment variable PEC_FWK_INSTALL is not set)
endif


# Set Boost variables
# A bit of special setup is required for lxplus
ifeq ($(shell hostname | grep -q lxplus; echo $$?), 0)
  ifeq ($(BOOST_ROOT), )
    # Set the path to Boost appropriate for lxplus
    BOOST_ROOT = /afs/cern.ch/sw/lcg/external/Boost/1.50.0_python2.7/x86_64-slc5-gcc47-opt
  endif
  
  BOOST_INCLUDE = $(BOOST_ROOT)/include/boost-1_50
  BOOST_LIB_POSTFIX = -gcc47-mt-1_50
else
  ifeq ($(BOOST_ROOT), )
    $(error Mandatory environment variable BOOST_ROOT is not set)
  endif
  
  BOOST_INCLUDE = $(BOOST_ROOT)/include
  BOOST_LIB_POSTFIX = 
endif

BOOST_LIB = $(BOOST_ROOT)/lib


# Define the flags to control make
CC = g++
INCLUDE = -Icore/include -Iextensions/include -I$(shell root-config --incdir) -I$(BOOST_INCLUDE)
OPFLAGS = -O2
CFLAGS = -Wall -Wextra -fPIC -std=c++11 $(INCLUDE) $(OPFLAGS)
#LDFLAGS = $(shell root-config --libs) -lTreePlayer -lHistPainter \
# -L$(BOOST_LIB) -lboost_filesystem$(BOOST_LIB_POSTFIX) $(PEC_FWK_INSTALL)/lib/libpecfwk.a \
# -Wl,-rpath=$(BOOST_LIB)
SOURCES = $(shell ls core/src/ | grep .cpp) $(shell ls extensions/src/ | grep .cpp)
OBJECTS = $(SOURCES:.cpp=.o)
# See http://www.gnu.org/software/make/manual/make.html#Substitution-Refs

# Define the search paths for the targets and dependencies. See
# http://www.gnu.org/software/make/manual/make.html#Directory-Search
vpath %.hpp core/include:extensions/include
vpath %.cpp core/src:extensions/src

# Define the phony targets
.PHONY: clean

# The default rule
all: libpecfwk.a libpecfwk.so

libpecfwk.a: $(OBJECTS)
	@ mkdir -p lib
	@ rm -f lib/$@
	@ ar -cq lib/$@ $+
# '$@' is expanded to the target, '$+' expanded to all the dependencies. See
# http://www.gnu.org/savannah-checkouts/gnu/make/manual/html_node/Automatic-Variables.html

libpecfwk.so: $(OBJECTS)
	@ mkdir -p lib
	@ rm -f lib/$@
	@ $(CC) -shared -Wl,-soname,$@.1 -o $@.1.0 $+
	@ mv $@.1.0 lib/
	@ ln -sf $@.1.0 lib/$@
	
%.o: %.cpp
	@ $(CC) $(CFLAGS) -c $< -o $@
# '$<' is expanded to the first dependency

clean:
	@ rm -f *.o
# '@' prevents the command from being printed to stdout
