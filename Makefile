#this target will compile all the files
parprog: placementRouting.cpp
	g++ -o parprog placementRouting.cpp
#here placementRouting.cpp is the dependency for the target main.o
clean:
	rm -rf parprog
