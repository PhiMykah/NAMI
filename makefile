# ---------
# PATHS
# ---------
DT = Datatypes
BTS_PATH = Tools/BTS
MOD = $(BTS_PATH)/$(MMOD)
MMOD = Modules
IO = FileIO
# ---------
# VARIABLES
# ---------

CXX=g++
CXXFLAGS= -g -Wall -std=c++17 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp
DC = DataContainers
ES = EsimModules
# IS = IsimModules
INCLUDES = Datatypes/$(DC).o $(MOD)/$(ES).o
BTS = $(BTS_PATH)/BTS.h

# Module variables
KMN = kmeansNANI
NN = nani

# Algorithm Variables
MSD = MeanSquareDeviation
EC = ExtendedComparison
READ = ReadNPY
CS = ComplementarySimilarity
MED = Medoid
OUTL = Outlier
DS = DiversitySelection
NI = NewIndex

BTS = $(DC).o $(ES).o $(READ).o $(MSD).o $(EC).o $(CS).o $(MED).o $(OUTL).o $(DS).o $(NI).o $(NN).o #$(IS).o 

OBJ_FILES = $(DT)/$(DC).o \
            $(MOD)/$(ES).o \
			$(IO)/$(READ).o \
            $(BTS_PATH)/$(MSD).o \
            $(BTS_PATH)/$(EC).o \
            $(BTS_PATH)/$(CS).o \
            $(BTS_PATH)/$(MED).o \
            $(BTS_PATH)/$(OUTL).o \
            $(BTS_PATH)/$(DS).o \
            $(BTS_PATH)/$(NI).o \
			$(MMOD)/$(KMN)/$(NN).o
			# $(MOD)/$(IS).o

# ----------------
# MAKEFILE SCRIPTS
# ----------------

all: main clean

main: $(BTS)
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) Tests/main.cpp -o main

clean:
	cd $(DT) && rm -f *.o
	cd $(BTS_PATH) && rm -f *.o
	cd $(MOD) && rm -f *.o
	cd $(IO) && rm -f *.o
	cd $(MMOD) && rm -rf */*.o

# datatest: $(DC).o
# 	$(CXX) $(CXXFLAGS) $(DT)/$(DC).o data_containers_test.cpp -o datatest

logictest: $(BTS)
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) Tests/logictest.cpp -o logictest

iotest: $(DC).o $(READ).o
	$(CXX) $(CXXFLAGS) $(DT)/$(DC).o $(IO)/$(READ).o Tests/io_test.cpp -o io_test
	make clean

kmeanstest: $(BTS)
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) Tests/kmeans_test.cpp -o kmeans_test
	
alatest: $(BTS)
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) Tests/ala10_test.cpp -o ala10_test

# Data Types and Containers Object
$(DC).o:
	$(CXX) $(CXXFLAGS) -c Datatypes/$(DC).cpp -o Datatypes/$(DC).o

# Esim Modules Object
$(ES).o: $(DT)/$(DC).o
	$(CXX) $(CXXFLAGS) -c $(MOD)/$(ES).cpp -o $(MOD)/$(ES).o

$(IS).o: $(DT)/$(DC).o
	$(CXX) $(CXXFLAGS) -c $(MOD)/$(IS).cpp -o $(MOD)/$(IS).o

# Read NPY Object
$(READ).o: $(DT)/$(DC).o
	$(CXX) $(CXXFLAGS) -c $(IO)/$(READ).cpp -o $(IO)/$(READ).o

# Mean Squared Deviation Object
# Requires:
#	- Default includes
$(MSD).o: $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(MSD).cpp -o $(BTS_PATH)/$(MSD).o

# Extended Comparison Object
# Requires:
# 	- Mean Square Deviation
#	- Default includes
$(EC).o: $(MSD).o $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(EC).cpp -o $(BTS_PATH)/$(EC).o

# Complimentary Similarity Object
# Requires:
# 	- Extended Comparison
#	- Default includes
$(CS).o: $(EC).o $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(CS).cpp -o $(BTS_PATH)/$(CS).o

# Medoid Calculations Object
# Requires:
# 	- Extended Comparison
#	- Complimentary Similarities
#	- Default includes
$(MED).o:  $(EC).o $(CS).o $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(MED).cpp -o $(BTS_PATH)/$(MED).o

# Outlier Trimming and Calculations Object
# Requires:
# 	- Extended Comparison
#	- Complimentary Similarities
#	- Default includes
$(OUTL).o: $(EC).o $(CS).o $(INCLUDES) 
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(OUTL).cpp -o $(BTS_PATH)/$(OUTL).o

# Get New Index N Object
# Requires:
# 	- Extended comparison
#	- Default includes
$(NI).o: $(EC).o $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(NI).cpp -o $(BTS_PATH)/$(NI).o

# Diversity Selection Object
# Requires:
# 	- Outlier
#	- Medoid
#	- Get New Index N
#	- Default includes
$(DS).o: $(OUTL).o $(MED).o $(NI).o $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(DS).cpp -o $(BTS_PATH)/$(DS).o

# kmeansNANI Nani Object
# Requires:
#	- Default includes
#	- Complimentary Similarities
#	- Diversity Selection
$(NN).o: $(DS).o $(CS).o $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(MMOD)/$(KMN)/$(NN).cpp -o $(MMOD)/$(KMN)/$(NN).o