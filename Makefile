LIBS_PATH=libs

DEBUG=1

ifdef DEBUG
	OPT = -O0 -g -ggdb
else
	OPT = -O3
endif

CFLAGS = -Wall -Wextra -std=c99 $(OPT)
OBJFLAGS = -fPIC
LINKFLAGS = -lalign -lstrbuf -lpthread -lz

INCS=-I $(LIBS_PATH) -I src
LIBS=-L $(LIBS_PATH)/string_buffer -L src
LINK=-lalign -lstrbuf -lpthread -lz

# Compile and bundle all non-main files into library
SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)

all: bin/smith_waterman src/libalign.a

# Build libraries only if they're downloaded
src/libalign.a: $(OBJS)
	[ -d libs/string_buffer ] && cd libs && $(MAKE)
	ar -csru src/libalign.a $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(OBJFLAGS) $(INCS) -c $< -o $@

bin/smith_waterman: src/tools/sw_cmdline.c src/libalign.a | bin
	$(CC) -o bin/smith_waterman $(SRCS) $(CFLAGS) $(TGTFLAGS) $(INCS) $(LIBS) src/tools/sw_cmdline.c $(LINKFLAGS)

examples: src/libalign.a
	cd examples; $(MAKE) LIBS_PATH=$(abspath $(LIBS_PATH))

prepare_libs:
	sudo apt-get install zlib1g-dev

bin:
	mkdir -p bin

clean:
	rm -rf bin src/*.o src/libalign.a
	cd examples && $(MAKE) clean

.PHONY: all clean examples
