all: test clean

test:
	g++ -g -c Tools/Data_Containers.cpp -o Tools/Data_Containers.o
	g++ -g -c Tools/Esim_Modules.cpp -o Tools/Esim_Modules.o
	g++ -g -c Tools/BTS.cpp -o Tools/BTS.o
	g++ -g Tools/BTS.o Tools/Data_Containers.o Tools/Esim_Modules.o main.cpp -o main

datatest:
	g++ -c Tools/Data_Containers.cpp -o Tools/Data_Containers.o
	g++ Tools/Data_Containers.o data_containers_test.cpp -o datatest

clean:
	rm Tools/*.o