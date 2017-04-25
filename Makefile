# Makefile for FEM1D
# DEVELOPER
#   Hom Nath Gharti
#   hngharti_AT_gmail_DOT_com 
#
# HISTORY
#   Feb 07,2013, HNG
#
# source and binary directory

default: fem1d

all: fem1d

clean:
	(cd src; make clean)

fem1d:
	(cd src; make $@)
