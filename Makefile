##########################
# Makefile for mmdecoder #
##########################

LIB = -lm

OBJECTS = \
	toolbox.o\

APP1 = mmdecoder
SRC1 = mmdecoder.c
OBJ1 = mmdecoder.o

DATE = $(shell date +\%Y-\%m-\%d)

###########
# Targets #
###########

default:
	make gcc

$(APP1): $(OBJ1) $(OBJECTS)
	$(CC) -o $(APP1) $(CFLAGS) $(OBJ1) $(OBJECTS) $(LIB)

clean:
	rm -f *.o $(APP1)

#################
# Architectures #
#################

gcc:
	make $(APP1) CC="gcc" CFLAGS="-O2 -Wall -Werror"	

###################
# Inference Rules #
###################

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

