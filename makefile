# ---------
# PATHS
# ---------
DT = Datatypes
BTS_PATH = Tools/BTS
MOD = $(BTS_PATH)/Modules

# ---------
# VARIABLES
# ---------

CXX=g++
CXXFLAGS= -g -Wall $(INC_PARAMS)
DC = DataContainers
ES = EsimModules
INCLUDES = Datatypes/$(DC).o $(MOD)/$(ES).o
BTS = $(BTS_PATH)/BTS.h

# Algorithm Variables
MSD = MeanSquareDeviation
EC = ExtendedComparison
CS = ComplementarySimilarity
MED = Medoid
OUTL = Outlier
DS = DiversitySelection
NI = NewIndex
ALI = Align

BTS = $(DC).o $(ES).o $(MSD).o $(EC).o $(CS).o $(MED).o $(OUTL).o $(DS).o $(NI).o $(ALI).o

OBJ_FILES = $(DT)/$(DC).o \
            $(MOD)/$(ES).o \
            $(BTS_PATH)/$(MSD).o \
            $(BTS_PATH)/$(EC).o \
            $(BTS_PATH)/$(CS).o \
            $(BTS_PATH)/$(MED).o \
            $(BTS_PATH)/$(OUTL).o \
            $(BTS_PATH)/$(DS).o \
            $(BTS_PATH)/$(NI).o \
            $(BTS_PATH)/$(ALI).o

# ----------------
# MAKEFILE SCRIPTS
# ----------------

all: main clean

main: $(BTS)
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) main.cpp -o main

clean:
	cd $(DT) && rm -f *.o
	cd $(BTS_PATH) && rm -f *.o
	cd $(MOD) && rm -f *.o

datatest: $(DC)
	$(CXX) $(CXXFLAGS) data_containers_test.cpp -o datatest

# Data Types and Containers Object
$(DC).o:
	$(CXX) $(CXXFLAGS) -c Datatypes/$(DC).cpp -o Datatypes/$(DC).o

# Esim Modules Object
$(ES).o: Datatypes/$(DC).o
	$(CXX) $(CXXFLAGS) -c $(MOD)/$(ES).cpp -o $(MOD)/$(ES).o

# Mean Squared Deviation Object
# Requires:
#	- Default includes
$(MSD).o: $(INCLUDES) # $(MSD).cpp $(MSD).h
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(MSD).cpp -o $(BTS_PATH)/$(MSD).o

# Extended Comparison Object
# Requires:
# 	- Mean Square Deviation
#	- Default includes
$(EC).o: $(MSD).o $(INCLUDES) # $(EC).cpp $(EC).h 
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(EC).cpp -o $(BTS_PATH)/$(EC).o

# Complimentary Similarity Object
# Requires:
# 	- Extended Comparison
#	- Default includes
$(CS).o: $(EC).o $(INCLUDES) # $(CS).cpp $(CS).h
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(CS).cpp -o $(BTS_PATH)/$(CS).o

# Medoid Calculations Object
# Requires:
# 	- Extended Comparison
#	- Complimentary Similarities
#	- Default includes
$(MED).o:  $(EC).o $(CS).o $(INCLUDES) # $(MED).cpp $(MED).h 
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(MED).cpp -o $(BTS_PATH)/$(MED).o

# Outlier Trimming and Calculations Object
# Requires:
# 	- Extended Comparison
#	- Complimentary Similarities
#	- Default includes
$(OUTL).o: $(EC).o $(CS).o $(INCLUDES) # $(OUTL).cpp $(OUTL).h 
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(OUTL).cpp -o $(BTS_PATH)/$(OUTL).o

# Diversity Selection Object
# Requires:
# 	- Outlier
#	- Medoid
#	- Default includes
$(DS).o: $(OUTL).o $(MED).o $(INCLUDES) # $(DS).cpp $(DS).h
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(DS).cpp -o $(BTS_PATH)/$(DS).o

# Get New Index N Object
# Requires:
#	- Default includes
$(NI).o: $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(NI).cpp -o $(BTS_PATH)/$(NI).o

# Alignment Operations Object
# Requires:
#	- Default includes
$(ALI).o: $(INCLUDES)
	$(CXX) $(CXXFLAGS) -c $(BTS_PATH)/$(ALI).cpp -o $(BTS_PATH)/$(ALI).o

