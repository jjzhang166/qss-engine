#Target OS variables
OS ?= unix
DEBUG ?= False

#Compiler and Linker
CXX         := gcc

#The Directories, Source, Includes, Objects, Binary 
SRCDIR      := .
BUILDDIR    := $(SRCDIR)/obj
SRCEXT      := c
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CXXFLAGS 	:= -Wall
ifeq ($(DEBUG),True)
CXXFLAGS 	+= -g  
endif
INC         := -I. -I./usr/include
RMS 		:= rm -rf

# Source files.
GRAPHSRC = $(SRCDIR)/graph_profile.c

# Objects
GRAPHOBJ=$(addprefix $(BUILDDIR)/grp_, $(notdir $(GRAPHSRC:.c=.o)))

# Make dependencies
DEPS = $(GRAPHOBJ:.o=.d)

all: graph 

$(BUILDDIR)/grp_%.o : $(SRCDIR)/%.c
	$(CXX) $(INC) $(CXXFLAGS) -MM -MT $@ -MF $(patsubst %.o,%.d,$@) $<
	$(CXX) $(INC) -c $< -o $@ $(CXXFLAGS) 

graph: $(GRAPHOBJ)
	 
$(GRAPHOBJ): | $(BUILDDIR)

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

-include $(DEPS)

.PHONY: clean

clean:
ifeq ($(OS),win32)
	$(RMS) $(DEPS) $(GRAPHOBJ)
else
	$(RMS) $(DEPS) $(GRAPHOBJ)
endif

help:
	@echo "make DEBUG=<True|False> OS=<unix|win32|mac>"
	@echo "Default values:"
	@echo ""
	@echo "DEBUG=False"
	@echo "OS=unix"
