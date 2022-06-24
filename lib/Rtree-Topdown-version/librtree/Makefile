CC=gcc
RM=rm -f
AR=ar
CPPFLAGS=-Iinclude
LDFLAGS=-L.
LDLIBS=-lrtree

SRCS=src/card.c \
			src/index.c \
			src/node.c \
			src/rect.c \
			src/split_l.c

OBJS=$(subst .c,.o,$(SRCS))

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	OFLAGS = -dynamiclib
endif
ifeq ($(UNAME_S),Darwin)
	CCFLAGS += -fPIC
	OFLAGS = -shared
endif

all: lib

lib: $(OBJS)
	$(AR) -cvq libtree.a $(OBJS)

depend: .depend

.depend: $(SRCS)
	$(RM) -f ./.depend
	$(CC) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)
	$(RM) librtree.a

dist-clean: clean
	$(RM) *~ .depend

include .depend
