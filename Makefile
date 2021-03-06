# /****************************************************************
# Copyright (C) Lucent Technologies 1997
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that both that the copyright notice and this
# permission notice and warranty disclaimer appear in supporting
# documentation, and that the name Lucent Technologies or any of
# its entities not be used in advertising or publicity pertaining
# to distribution of the software without specific, written prior
# permission.
#
# LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
# IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
# IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
# ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
# THIS SOFTWARE.
# ****************************************************************/

CFLAGS = -g -Wall -O2

CC = gcc
CPP = g++

YACC = bison -y
YACC = yacc
YFLAGS = -d

OFILES = b.o main.o parse.o proctab.o tran.o lib.o run.o lex.o addon.o edlib.o md5.o

SOURCE = awk.h ytab.c ytab.h proto.h awkgram.y end_adapter.h lex.c b.c main.c \
	maketab.c parse.c lib.c run.c tran.c proctab.c addon.c md5.c

LISTING = awk.h proto.h awkgram.y lex.c b.c main.c maketab.c parse.c \
	lib.c run.c tran.c addon.c md5.c

SHIP = README FIXES $(SOURCE) ytab[ch].bak makefile  \
	 awk.1

UNAME = $(shell uname -s)

bioawk:ytab.o $(OFILES)
	$(CPP) $(CFLAGS) ytab.o $(OFILES) $(ALLOC) -o $@ -lm -lz
	cp bioawk bioawk_cas

$(OFILES):	awk.h ytab.h proto.h addon.h end_adapter.h

ytab.o:	awk.h proto.h awkgram.y
	$(YACC) $(YFLAGS) awkgram.y
	mv y.tab.c ytab.c
	mv y.tab.h ytab.h
	$(CC) $(CFLAGS) -c ytab.c

proctab.c:	maketab
	./maketab >proctab.c

edlib.o:	edlib.h
	$(info OS $(UNAME))
ifeq ($(UNAME),Darwin)
	cp edlib.o.mac edlib.o
else ifeq ($(UNAME),Linux)
	cp edlib.o.linux edlib.o
else
	cp edlib.o.orig edlib.o
endif

ytab.h:

maketab:	ytab.h maketab.c
	$(CC) $(CFLAGS) maketab.c -o maketab

names:
	@echo $(LISTING)

clean:
	rm -fr a.out *.o *.obj maketab maketab.exe *.bb *.bbg *.da *.gcov *.gcno *.gcda awk ytab.o proctab.c *.dSYM

debug: CFLAGS = -g3 -Wall
debug: clean
debug: bioawk

all: bioawk