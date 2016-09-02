######################################################################
#			Makefile for DISOPRED3                       #
######################################################################

#### CXX	Command to execute C++ compiler
#### CC		Command to execute C compiler
#### CXXFLAGS	Flags to pass to C++ compiler
#### CCFLAGS	Flags to pass to C compiler
#### LIBS	Libraries

INCDIR		= .
CXX      	= g++
CXXFLAGS 	= -Wall -O
CC       	= gcc
CCFLAGS  	= -Wall -O
LIBS		= -lm

all:		chkparse disopred2 diso_neu_net diso_neighb combine svm-predict

install:
		mkdir ../bin/
		mv -t ../bin chkparse disopred2 diso_neu_net diso_neighb combine svm-predict

clean:
		rm -f chkparse disopred2 diso_neu_net diso_neighb combine svm-predict *~
		rm -rf ../bin/

disopred2:	pred_dis.cpp weights.cpp
		$(CXX) $(CXXFLAGS) pred_dis.cpp weights.cpp $(LIBS) -o disopred2

diso_neu_net:	disord_pred.c disord_net.h
		$(CC) $(CCFLAGS) disord_pred.c $(LIBS) -o diso_neu_net

diso_neighb:	diso_neighb.c
		$(CC) $(CCFLAGS) diso_neighb.c -o diso_neighb

combine:	disordcomb_pred.c
		$(CC) $(CCFLAGS) disordcomb_pred.c $(LIBS) -o combine

svm-predict:	svm-predict.c svm.cpp svm.h
		$(CXX) $(CXXFLAGS) svm-predict.c svm.cpp $(LIBS) -o svm-predict

chkparse:	chkparse.c
		$(CC) $(CFLAGS) chkparse.c $(LIBS) -o chkparse
