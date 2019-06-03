ga_exe : main.o gh.o
	g++ -o ga_exe main.o gh.o

main.o: main.cpp
	g++ -c -o main.o main.cpp

gh.o: GraphHandler.cpp
	g++ -c -o gh.o GraphHandler.cpp
clean:
	rm *.o ga_exe
