# Makefile template for shared library

CC = clang # C compiler
CFLAGS = -fPIC -Wall -Wextra -O2 -g # C flags
LDFLAGS = -shared -Wl,-whole-archive -L ../lib -lwsmp64 -Wl,-no-whole-archive -L /home/bates/git/julia/usr/lib -lopenblas -lpthread -lgfortran -lm # linking flags
RM = rm -f  # rm command
TARGET_LIB = libwsmp.so # target lib

SRCS = wsmp.c  # source files
OBJS = $(SRCS:.c=.o)

.PHONY: all
all: ${TARGET_LIB}

$(TARGET_LIB): $(OBJS)
	$(CC) ${LDFLAGS} -o $@ $^

$(SRCS:.c=.d):%.d:%.c
	$(CC) $(CFLAGS) -MM $< >$@

include $(SRCS:.c=.d)

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB} ${OBJS} $(SRCS:.c=.d)
