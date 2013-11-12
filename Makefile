# Define some macros
#MACROS = -D B_EFF_PATH=\"$(shell cd ..; pwd)/BTagEff/\"

# Define the flags to control make
CC = g++
INCLUDE = -Icore/include -Iextensions/include -I$(shell root-config --incdir)
OPFLAGS = -O2
CFLAGS = -Wall -Wextra -fPIC -std=c++11 $(INCLUDE) $(OPFLAGS) $(MACROS)
#LDFLAGS = $(shell root-config --libs) -lTreePlayer -lHistPainter
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
