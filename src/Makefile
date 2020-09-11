run: main.o SG.o getSG.o extendAnchor.o mergeSG.o unionFind.o
	g++ -std=c++11 -o ../bin/dragom.exe main.o SG.o getSG.o extendAnchor.o mergeSG.o unionFind.o
main.o: main.cpp
	g++ -std=c++11 -c main.cpp
SG.o: SG.cpp
	g++ -std=c++11 -c SG.cpp
getSG.o: getSG.cpp
	g++ -std=c++11 -c getSG.cpp
mergeSG.o: mergeSG.cpp
	g++ -std=c++11 -c mergeSG.cpp
extendAnchor.o: extendAnchor.cpp
	g++ -std=c++11 -c extendAnchor.cpp
unionFind.o: unionFind.cpp
	g++ -std=c++11 -c unionFind.cpp
clean:
	rm -rf *.o *.exe
