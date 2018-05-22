#############################################################################
#   
#	Makefile for DISSECT
#   
#############################################################################

#############################################################################
# Leave blank after "=" to disable using BOOST or ZLIB libraries; put "= 1" to enable
# ZLIB can only be used with BOOST. Otherwise it will be ignored.
BOOSTLIB = 1
ZLIB = 1
BGEN = 1

#############################################################################
# Select MPI implementation
MPI = MPICH2
# MPI = OPENMPI

#############################################################################
# Put C++ compiler here. It is platform specific. 
CXX = mpic++
# CXX = CC
# CXX = icc

#############################################################################
# scaLAPACK implementation
SCALAPACK = MKL
# SCALAPACK = LIBSCI

#############################################################################
# PATHs
MKL_PATH = /opt/intel/mkl
BOOST_PATH = 
BGENFOLDER = /home/ocanela/software/bgen/bgen

#############################################################################
# Dynamic linking (1 for dynamic linking, empty otherwise)
DYNAMIC = 1


#############################################################################
# Print traceback (1 for creating traceback, empty otherwise)
TRACEBACK = 1


#############################################################################
# Additional flags for the compiler
CXXFLAGS = -g -O3 #-Wall -pg
#-O3

CRM = rm

##############################################################################
#
#       Build script
#
##############################################################################

#----------------------------------------------------------------------

ifeq ($(MYHOSTNAME),eddie)
  MKL_PATH=/exports/applications/apps/SL6/intel/mkl
  CXX = mpic++
endif

ifeq ($(MYHOSTNAME),archer)
  MKL_PATH=/opt/intel/composerxe/mkl
  CXX = CC
# CXXFLAGS = -I $(MKL_PATH)/include -openmp 
  DYNAMIC = 
endif

#----------------------------------------------------------------------

CXXFLAGS += -fopenmp -m64
LIB = 

# Dynamic linking?
ifneq ($(DYNAMIC), 1)
  CXXFLAGS += -static
endif

ifeq ($(BOOSTLIB), 1)
  ifneq ($(BOOST_PATH),)
    CXXFLAGS += -I $(BOOST_PATH)
  else
    ifdef BOOST_DIR
      CXXFLAGS += -I $(BOOST_DIR)
    endif
  endif
endif

ifeq ($(BGEN),1)
  LIB += -L$(BGENFOLDER)/build/ -L$(BGENFOLDER)/build/3rd_party/zstd-1.1.0
endif


# scaLAPACK implementation?
ifeq ($(SCALAPACK),MKL)
  CXXFLAGS += -I $(MKL_PATH)/include
  
  ifeq ($(MPI), OPENMPI)
    MKLMPILIB = mkl_blacs_openmpi_lp64
  else
    MKLMPILIB = mkl_blacs_intelmpi_lp64
  endif
  
  ifeq ($(DYNAMIC), 1)
    LIB += -L$(MKL_PATH)/lib/intel64/ -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -l$(MKLMPILIB) -Wl,-rpath=$(MKL_PATH)/lib/intel64/
  else
    LIB += $(MKL_PATH)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKL_PATH)/lib/intel64/libmkl_intel_lp64.a $(MKL_PATH)/lib/intel64/libmkl_core.a $(MKL_PATH)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group $(MKL_PATH)/lib/intel64/lib$(MKLMPILIB).a  
  endif
endif

LIB += -ldl -lpthread -lm

# CXXFLAGS += -I $(MKL_PATH)/include -fopenmp -m64 #-Wall 

OUTPUT = dissect

# Define some specific flags

ifeq ($(BOOSTLIB),1)
  CXXFLAGS += -DBOOSTLIB
  ifeq ($(ZLIB),1)
    LIB += -lboost_iostreams
    CXXFLAGS += -DZLIB
  endif
endif

ifeq ($(TRACEBACK),1)
  CXXFLAGS += -DCREATETRACEBACK
endif

ifeq ($(BGEN),1)
  CXXFLAGS += -DBGEN
  CXXBGENFLAGS += -std=c++11 -I$(BGENFOLDER)/genfile/include/ -I$(BGENFOLDER)/3rd_party/zstd-1.1.0/lib 
  LIB += -lbgen -lzstd -lz
endif

SRC = options.cpp communicator.cpp matrix.cpp misc.cpp genotype.cpp reml.cpp gwas.cpp covariancematrix.cpp phenotype.cpp \
	covariate.cpp auxiliar.cpp simulatephenotype.cpp range.cpp message.cpp pca.cpp results.cpp analysis.cpp predictphenotype.cpp \
	test.cpp main.cpp singlereml.cpp multireml.cpp kernel.cpp blockmatrix.cpp accuracybysnp.cpp labeledmatrix.cpp groupeffects.cpp \
	glm.cpp glmm.cpp pcagentemp.cpp igwas.cpp mpresiduals.cpp gwasmp.cpp #finetest.cpp
SRC2 = genotypebgen.cpp
HDR = global.h communicator.h main.h options.h matrix.h misc.h genotype.h reml.h gwas.h covariancematrix.h phenotype.h \
	covariate.h auxiliar.h simulatephenotype.h range.h message.h pca.h results.h analysis.h predictphenotype.h \
	test.h singlereml.h multireml.h kernel.h blockmatrix.h accuracybysnp.h labeledmatrix.h groupeffects.h \
	glm.h glmm.h pcagentemp.h igwas.h mpresiduals.h #finetest.h
OBJ = $(SRC:.cpp=.o)
OBJ2 = $(SRC2:.cpp=.o)

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(OBJ2) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OBJ2) : $(HDR)
	$(CXX) $(CXXFLAGS) $(CXXBGENFLAGS) -c $*.cpp


$(OUTPUT) : $(OBJ)
$(OUTPUT) : $(OBJ2)

FORCE:

clean:
	$(CRM) *.o *~ dissect -f

