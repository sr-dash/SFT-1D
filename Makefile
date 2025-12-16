# Set FORTRAN90 compiler
# On my Mac the pdflatex location is
# /usr/local/texlive/2023/bin/universal-darwin/pdflatex
FC = gfortran

# Location of files for netcdf library
# Replace it according to your system configuration. 

NETCDF = -I/data/sdash/Softwares/netcdf/include
NETCDFLIB = -L/data/sdash/Softwares/netcdf/lib -lnetcdff

# Set compiler flags
FFLAGS = -O3 -Wall -Wno-unused-variable -fcheck=all -Wtabs -Wextra -ffast-math -funroll-loops --param max-unroll-times=5
parallelflags = -fopenmp

# Set build parameters
TARGET = SFT_1D

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

OBJDIR = obj
BINDIR = bin
DOCDIR = doc
SRCDIR = src

OBJFILES = main.o variables.o write_data.o output.o grid_SFT.o flows.o evolSFT.o init_condition.o
FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(OBJDIR)

# Rule to build the fortran files

%.o: $(SRCDIR)/%.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	@cd $(SRCDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) $(parallelflags) -J $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: $(SRCDIR)/%.F90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	@cd $(SRCDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) $(parallelflags) -J $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(parallelflags) $(FFLAGS) -J $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(NETCDFLIB)

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.sh.*  *.log $(DOCDIR)/user*.pdf

PDF:	
	@cd doc; make pdf; make clean

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

help:
	@printf "\nusage : make <commands> \n\nthe following commands are available : \n\n"
	@printf "<> : Generates the SFT executable in the bin directory.\n\n"
	@printf "PDF : Creates the PDF documentation files in the doc directory.\n\n"
	@printf "clean : Removes the executable and supporting compilation files.\n\n"
	@printf "help : Shows the available commands.\n\n"
	@printf "\n"

# All the dependencies
variables.o:variables.f90
output.o:output.f90 variables.o
grid_SFT.o:grid_SFT.f90 flows.o variables.o
write_data.o:write_data.f90 variables.o grid_SFT.o
flows.o:flows.f90 variables.o
evolSFT.o:evolSFT.f90 variables.o grid_SFT.o flows.o
init_condition.o:init_condition.f90 variables.o grid_SFT.o

main.o: main.f90 variables.o write_data.o output.o grid_SFT.o flows.o evolSFT.o init_condition.o
