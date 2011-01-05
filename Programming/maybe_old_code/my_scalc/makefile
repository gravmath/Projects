
scalc: scalc.o param.o array.o tensor.o msfield.o util.o
	g++ scalc.o array.o tensor.o msfield.o param.o util.o -O -o scalc 
scalc.o: scalc.cpp 
	g++ -c scalc.cpp
param.o: param.cpp
	g++ -c param.cpp
array.o: array.cpp
	g++ -c array.cpp
tensor.o: tensor.cpp
	g++ -c tensor.cpp
msfield.o: msfield.cpp
	g++ -c msfield.cpp
util.o: util.cpp
	g++ -c util.cpp


