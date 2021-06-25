CC=		gcc
CFLAGS= -Wall -g -O2 -Wextra
LIBS=	-lm -lz
SRC_DIR=	./src
SRC=	./src
OBJ_DIR=	./obj
BIN_DIR=	./bin
BIN=	./bin
MKDIR_P = mkdir -p
PROG=	bam-amplicon

SRCS=	$(wildcard $(SRC_DIR)/*.c)
OBJS_=	$(notdir $(SRCS))
OBJS=	$(OBJS_:%.c=%.o)
AOBJS=	$(OBJS:%=$(OBJ_DIR)/%)



all: $(BIN_DIR) $(OBJ_DIR) $(PROG)
#all: $(BIN_DIR) $(OBJ_DIR) ls

bam-amplicon: $(OBJS)
	$(CC) $(CFLAGS) $(LIBS) $(AOBJS) -o $(BIN_DIR)/$@

%.o: $(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $(LIBS) $^ -o $(OBJ_DIR)/$@	

ls: $(OBJS) 
	ls -l

$(BIN_DIR): 
	mkdir -p $(BIN_DIR)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR)
