GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
CFLAGS = -Wall -O2 -std=gnu99 $(shell pkg-config --cflags  cairo)
CFLAGS += -DVERSION=\"$(GIT_VERSION)\"
CFLAGS += -DGRAPHICS
GRAPHLIBS = $(shell pkg-config --libs cairo)
LIBS = -lm -lgsl -lgslcblas
CC = gcc

HDRS = structures.h
OBJS = main.o collision_rules.o front.o graphics.o

default:$(OBJS)
	$(CC)  $(CFLAGS) $(OBJS) $(GRAPHLIBS) $(LIBS) -o ../disks-event

clean:
	-/bin/rm	$(OBJS) ../disks-event


main.o:			Makefile $(HDRS)
front.o:		Makefile $(HDRS)
graphics.o:		Makefile $(HDRS)
collision_rules.o:	Makefile $(HDRS)
