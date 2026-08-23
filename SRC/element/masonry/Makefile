include ../../../Makefile.def

OBJS       = MasonPan12.o \
	MasonPan3D.o \
	BeamGT.o \
	MasonryMacro2D.o \
	MasonryMacro3D.o \
	OPS_MasonryMacro.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
