# Makefile for FEM1D
# REVISION
#   HNG, Feb 07,2013
# source and binary directory
SRC = .
BIN = ../bin

# serial fortran compiler
FC = gfortran
FFLAGS = -std=f2003 -g -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wall -Wline-truncation -fno-trapping-math -fno-backslash -Wunused
#FFLAGS = -std=f2008 -O5
#FFLAGS = std03 -O3
FC_COMPILE = $(FC) $(FFLAGS)

# object tree
OBJ = $(BIN)/set_precision.o   \
			$(BIN)/global.o          \
			$(BIN)/read_input.o          \
			$(BIN)/mesh.o            \
			$(BIN)/interpolate.o     \
			$(BIN)/form.o            \
			$(BIN)/bc.o              \
			$(BIN)/solve.o           \
			$(BIN)/plot.o            \
			$(BIN)/fem1d.o

# main executable for 1D FEM
fem1d = $(BIN)/fem1d

default: fem1d

clean:
	rm -f *.mod $(BIN)/*.o core

all: clean fem1d

fem1d: $(OBJ)
	$(FC_COMPILE) -o $(fem1d) $(OBJ)

$(BIN)/%.o: $(SRC)/%.f90
	${FC_COMPILE} ${FFLAGS} -c -o $@ $<
