
testCode: TestCode.o
	g++ -g -o TestCode TestCode.o -std=c++17

TestCode.o: TestCode.cpp qbMatrix.h
	g++ -o TestCode.o -c TestCode.cpp -std=c++17

.PHONY: clean

clean:
	rm *.o TestCode
