# Makefile template for shared library

CC = clang # C compiler
CFLAGS = -fPIC -Wall -Wextra -O2 -g # C flags
LDFLAGS = -shared -Wl,-whole-archive -L . -lwsmp64 -Wl,-no-whole-archive -L /home/bates/git/julia/usr/lib -lopenblas -lpthread -lgfortran -lm # linking flags
RM = rm -f  # rm command
TARGET_LIB = libwsmp.so # target lib

.PHONY: all
all: ${TARGET_LIB}

$(TARGET_LIB): 
	$(CC) ${LDFLAGS} -o $@ $^

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB}
