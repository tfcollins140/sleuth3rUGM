#
# Makefile by Tommy E. Cathey of NESC (919)541-1500
#
#   Modified by D. Donato: Tuesday, Feb. 15, 2005
#

VPATH=./GD
GD_LIB = ./GD

#
# choose a compiler
#
#CC = cc
#CC = gcc
CC = mpicc

#
# select flags
#   -UNDEBUG (turns asserts on; for development code only)
#   -DNDEBUG (turns asserts off; for production code; faster execution)
#   -DMPI (if running on an MPI machine else -UMPI)
#
CFLAGS= -g -Og -UNDEBUG -DMPI -I$(GD_LIB) -I/usr/local/mpich2-1.0/include

CLIBS = -L./ -L./GD/ -L/usr/local/mpich2-1.0/lib -lgd -lm -lc
#CLIBS = -L./ -L./GD/ -lgd -lm -lc

SRCS_W_HDRS   = stats_obj.c timer_obj.c proc_obj.c transition_obj.c coeff_obj.c landclass_obj.c deltatron.c growth.c output.c utilities.c spread.c random.c scenario_obj.c igrid_obj.c gdif_obj.c pgrid_obj.c memory_obj.c wgrid_obj.c grid_obj.c color_obj.c driver.c input.c

SRCS_WO_HDRS  = main.c

SRCS = ${SRCS_W_HDRS} ${SRCS_WO_HDRS}

HDRS  = ${SRCS_W_HDRS:.c=.h} globals.h ugm_typedefs.h ugm_defines.h ugm_macros.h
OBJS = ${SRCS:.c=.o}

grow : $(OBJS)
	$(CC) -I$(GD_LIB) $(OBJS) -o grow $(CLIBS)
clean :
	rm $(OBJS)
clean_all :
	rm $(OBJS) grow
