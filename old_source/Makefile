#!/bin/bash

gsllib=-L/Users/xywu/opt/gsl/lib
gslinclude=-I/Users/xywu/opt/gsl/include

INIReader.o:./inih/INIReader.h ./inih/INIReader.cpp
	g++ ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -static -std=c++14 -c ./inih/INIReader.cpp
eos_base.o:./eos/eos_base.h ./eos/eos_base.cpp
	g++ ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -static -std=c++14 -c ./eos/eos_base.cpp
eos.o:./eos/eos.h ./eos/eos.cpp
	g++ ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -static -std=c++14 -c ./eos/eos.cpp
eos_neos.o:./eos/eos_hotQCD.h ./eos/eos_neos.cpp
	g++ ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -static -std=c++14 -c ./eos/eos_neos.cpp
ini.o:./inih/ini.h ./inih/ini.c
	g++  ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -static -std=c++14 -c ./inih/ini.c
parameter.o:parameter.h parameter.cpp
	g++  ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -static -std=c++14 -c parameter.cpp
reader.o:reader.h reader.cpp
	g++  ${gsllib} ${gslinclude} -lgsl -lgslcblas  -O3 -static -std=c++14 -c reader.cpp
grid.o:grid.h grid.cpp
	g++  ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -static -std=c++14 -c grid.cpp
main.o:main.cpp
	g++  ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -static -std=c++14 -c main.cpp
main:main.o INIReader.o eos.o eos_base.o eos_neos.o  ini.o parameter.o reader.o grid.o 
	g++ ${gsllib} ${gslinclude} -lgsl -lgslcblas -O3 -std=c++14  ini.o  INIReader.o eos.o eos_base.o eos_neos.o  parameter.o reader.o grid.o main.o -o main

clean:
	rm -f *.o main
