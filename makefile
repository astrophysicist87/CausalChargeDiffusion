SHELL=/bin/sh

SFLDIR=./special_function_library/cmlib

SRCS= \
CausalChargeDiffusion.cpp \
gauss_quadrature.cpp \
lib.cpp \
$(SFLDIR)/chf.cpp \
$(SFLDIR)/chfs.cpp \
$(SFLDIR)/complex.cpp \
$(SFLDIR)/g.cpp \
$(SFLDIR)/f211.cpp \
$(SFLDIR)/f212.cpp \
$(SFLDIR)/f213.cpp \
$(SFLDIR)/polyrt.cpp \
$(SFLDIR)/digam.cpp \
$(SFLDIR)/ebznew.cpp

HDRS= \
gauss_quadrature.h \
defs.h \
lib.h \
$(SFLDIR)/cmlib.h \
$(SFLDIR)/complex.h \
$(SFLDIR)/protom.h

MAKEFILE=makefile

COMMAND=run_CCD

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS=  -pg -g -O3
#WARNFLAGS= -ansi -pedantic -Wall -W \
#   -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
#   -Wwrite-strings -fshort-enums -fno-common 
WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include
 
$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

CausalChargeDiffusion.o : CausalChargeDiffusion.cpp defs.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c CausalChargeDiffusion.cpp -o CausalChargeDiffusion.o

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
