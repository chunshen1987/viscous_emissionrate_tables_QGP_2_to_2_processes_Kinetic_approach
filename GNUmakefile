# ===========================================================================
#  Makefile photon_emission                           Chun Shen May 5, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS = -O3

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	QGPphotonEmissionrates.e
endif

SRC		=	main.cpp Arsenal.cpp QGP_2to2_Scattering_Kinetic.cpp ParameterReader.cpp \
                  Physicalconstants.cpp gauss_quadrature.cpp \

INC		= 	Arsenal.h Stopwatch.h QGP_2to2_Scattering_Kinetic.h ParameterReader.h\
                  Physicalconstants.h gauss_quadrature.h


# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(HDF5LD) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(LDFLAGS) $(HDF5FLAGS) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
main.cpp : Arsenal.h Stopwatch.h QGP_2to2_Scattering_Kinetic.h ParameterReader.h
QGP_2to2_Scattering_Kinetic.cpp : QGP_2to2_Scattering_Kinetic.h Physicalconstants.h ParameterReader.h gauss_quadrature.h Arsenal.h
ParameterReader.cpp : ParameterReader.h Arsenal.h
Arsenal.cpp : Arsenal.h
Physicalconstants.cpp : Physicalconstants.h
gauss_quadrature.cpp : gauss_quadrature.h
